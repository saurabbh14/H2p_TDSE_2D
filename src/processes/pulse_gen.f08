!> This module contains the subroutine for generating the laser pulse
!> and the envelope functions.  Uses an allocatable array of per-pulse
!> parameters so any number of lasers (1..N) can be defined.
module pulse_mod
    use global_vars, only: dp, Nt, output_data_dir, dt, time
    use InputVars,    only: LaserParams
    use data_au
    use FFTW3
    use omp_lib
    implicit none
    private
    public:: pulse_param

    !> Internal per-pulse data — one instance per laser
    type :: single_pulse_data
        character(150) :: envelope_shape = ""
        real(dp) :: lambda     = 0._dp     ! wavelength in a.u. after init
        real(dp) :: omega      = 0._dp     ! angular frequency in a.u.
        real(dp) :: tp         = 0._dp     ! pulse duration in a.u.
        real(dp) :: t_mid      = 0._dp     ! pulse midpoint in a.u.
        real(dp) :: alpha0     = 0._dp     ! quiver amplitude in a.u.
        real(dp) :: phi        = 0._dp     ! phase in rad
        real(dp) :: rise_time  = 0._dp     ! rise time in a.u. (trapezoidal)
        real(dp) :: pulse_offset = 0._dp   ! internal offset (always 0 for now)
        real(dp), allocatable :: alpha_t(:)   ! quiver field  alpha(t)
        real(dp), allocatable :: E_field(:)   ! electric field E(t)
        real(dp), allocatable :: A_field(:)   ! vector potential A(t)
        real(dp), allocatable :: env(:)       ! envelope g(t)
    end type single_pulse_data

    !> Main pulse object — holds an array of pulses plus the total fields.
    type :: pulse_param
        integer :: N_pulses = 0
        type(single_pulse_data), allocatable :: pulses(:)
        real(dp), allocatable :: alpha_t(:)   ! total quiver field
        real(dp), allocatable :: El(:)        ! total electric field
        real(dp), allocatable :: Al(:)        ! total vector potential
    contains
        procedure :: initialize    => initialize_from_lasers
        procedure :: param_print   => print_pulse_param
        procedure :: generate      => generate_pulse
        procedure :: write_to_file => write_pulse_to_file
        procedure :: deallocate_env => deallocate_envelope
        procedure :: deallocate_field => deallocate_field
        procedure :: spectra       => field_spectra
        procedure :: deallocate_all
    end type pulse_param

contains

    !> Initialise internal pulse data from the LaserParams array read by
    !> the input module.  Replaces the old read_pulse_params (no file I/O).
    subroutine initialize_from_lasers(this, input_lasers, N_lasers)
        class(pulse_param), intent(inout) :: this
        type(LaserParams), intent(in)     :: input_lasers(:)
        integer, intent(in)               :: N_lasers
        integer :: i

        this%N_pulses = N_lasers
        allocate(this%pulses(N_lasers))

        do i = 1, N_lasers
            associate(pl => this%pulses(i), inp => input_lasers(i))
                pl%envelope_shape = inp%envelope
                pl%lambda     = inp%lambda
                pl%tp         = inp%tp
                pl%t_mid      = inp%t_mid
                pl%alpha0     = inp%alpha0
                pl%phi        = inp%phi
                pl%rise_time  = inp%rise_time
                pl%pulse_offset = 0._dp
            end associate
        end do

        ! Convert units and compute derived quantities
        do i = 1, N_lasers
            associate(pl => this%pulses(i))
                pl%tp        = pl%tp        / au2fs
                pl%t_mid     = pl%t_mid     / au2fs
                pl%rise_time = pl%rise_time / au2fs
                pl%omega     = (1._dp / (pl%lambda * 1.e-7_dp)) * cm2au
                pl%phi       = pl%phi       * pi
            end associate
        end do
    end subroutine initialize_from_lasers

    !> Print all pulse parameters.
    subroutine print_pulse_param(this)
        class(pulse_param), intent(in) :: this
        integer :: i
        real(dp) :: field_strength, intensity

        print*, "Laser parameters:"
        print*, "Number of pulses:", this%N_pulses
        do i = 1, this%N_pulses
            associate(pl => this%pulses(i))
                field_strength = pl%alpha0 * pl%omega**2
                intensity      = field_strength**2 * 3.509e16_dp
                print*, "--- Laser #", i, "---"
                print*, "  Envelope shape:  ", trim(pl%envelope_shape)
                print*, "  Wavelength:      ", sngl(pl%lambda), "nm"
                print*, "  Quiver amplitude:", sngl(pl%alpha0), "a.u."
                print*, "  Field strength:  ", sngl(field_strength), "a.u."
                print*, "  Intensity:       ", sngl(intensity), "W/cm2"
                print*, "  Pulse width (tp):", sngl(pl%tp * au2fs), "fs"
                print*, "  Pulse midpoint:  ", sngl(pl%t_mid * au2fs), "fs"
                print*, "  Phase (phi):     ", sngl(pl%phi / pi), "pi"
                print*, "  Rise time:       ", sngl(pl%rise_time * au2fs), "fs"
            end associate
        end do
        print*, "------------------------------------------------------"
    end subroutine print_pulse_param

    !> Generate the complete pulse: loop over all defined pulses,
    !> compute individual envelopes/fields, then sum into total arrays.
    subroutine generate_pulse(this)
        use differentiation, only: central_diff_on_grid
        class(pulse_param), intent(inout) :: this
        integer :: i, k, n_cycles
        real(dp) :: ttime, TU_eff, tp_eff

        print*
        print*, "Pulse generation..."

        ! Allocate per-pulse arrays
        do i = 1, this%N_pulses
            associate(pl => this%pulses(i))
                allocate(pl%alpha_t(Nt), pl%E_field(Nt), pl%A_field(Nt), pl%env(Nt))
                pl%alpha_t  = 0._dp
                pl%E_field  = 0._dp
                pl%A_field  = 0._dp
                pl%env      = 0._dp
            end associate
        end do

        ! Allocate total-field arrays
        allocate(this%alpha_t(Nt), this%El(Nt), this%Al(Nt))
        this%alpha_t = 0._dp
        this%El      = 0._dp
        this%Al      = 0._dp

        ! Generate each pulse
        do i = 1, this%N_pulses
            call generate_single(this%pulses(i), i)
        end do

        ! Sum contributions into total fields
        do i = 1, this%N_pulses
            associate(pl => this%pulses(i))
                this%alpha_t = this%alpha_t + pl%alpha_t
                this%Al      = this%Al      + pl%A_field
                this%El      = this%El      + pl%E_field
            end associate
        end do

        print*, "Pulse generation complete."
        call this%write_to_file()
    end subroutine generate_pulse

    !> Generate a single pulse (internal helper).
    subroutine generate_single(pl, idx)
        use differentiation, only: central_diff_on_grid
        type(single_pulse_data), intent(inout) :: pl
        integer, intent(in) :: idx
        integer :: k, n_cycles
        real(dp) :: ttime, TU_eff, tp_eff
        character(20) :: label

        write(label, '(A,I0)') 'Laser', idx

        select case(trim(pl%envelope_shape))
        case("cos2")
            do k = 1, Nt
                pl%env(k)    = cos2(time(k), pl%tp, pl%t_mid, pl%pulse_offset)
                pl%alpha_t(k) = pl%alpha0 * pl%env(k) &
                    & * cos(pl%omega * (time(k) - pl%t_mid - pl%pulse_offset) + pl%phi)
            end do
            call central_diff_on_grid(pl%alpha_t, Nt, dt, pl%A_field)
            call central_diff_on_grid(pl%A_field, Nt, dt, pl%E_field)
            pl%E_field = -pl%E_field

        case("sin2")
            n_cycles = int(pl%tp * pl%omega / (2._dp * pi)) + 1
            n_cycles = max(n_cycles, 1)
            tp_eff = 2._dp * pi * real(n_cycles, dp) / pl%omega
            do k = 1, Nt
                pl%env(k)    = sin2(time(k), tp_eff, pl%t_mid, pl%pulse_offset)
                pl%alpha_t(k) = pl%alpha0 * pl%env(k) &
                    & * sin(pl%omega * (time(k) - pl%t_mid - pl%pulse_offset - tp_eff/2) + pl%phi)
                pl%A_field(k) = sin2_vector_pulse(time(k), tp_eff, pl%t_mid, pl%alpha0, &
                    & pl%omega, pl%phi, pl%pulse_offset)
                pl%E_field(k) = sin2_electric_pulse(time(k), tp_eff, pl%t_mid, pl%alpha0, &
                    & pl%omega, pl%phi, pl%pulse_offset)
            end do

        case("gaussian")
            do k = 1, Nt
                pl%env(k)    = gaussian(time(k), pl%tp, pl%t_mid)
                pl%alpha_t(k) = pl%alpha0 * pl%env(k) &
                    & * cos(pl%omega * (time(k) - pl%t_mid - pl%pulse_offset) + pl%phi)
                pl%A_field(k) = gaussian_vector_pulse(time(k), pl%tp, pl%t_mid, pl%alpha0, &
                    & pl%omega, pl%phi, pl%pulse_offset)
                pl%E_field(k) = gaussian_electric_pulse(time(k), pl%tp, pl%t_mid, pl%alpha0, &
                    & pl%omega, pl%phi, pl%pulse_offset)
            end do

        case("trapezoidal")
            n_cycles = int(pl%rise_time * pl%omega / (2._dp * pi)) + 1
            n_cycles = max(n_cycles, 1)
            TU_eff = 2._dp * pi * real(n_cycles, dp) / pl%omega
            print*, trim(label) // " trapezoidal boundary times (fs):"
            print'(a,f10.4)', "  Pulse start: ", (pl%t_mid - pl%tp/2 - TU_eff + pl%pulse_offset)*au2fs
            print'(a,f10.4)', "  Rise end:    ", (pl%t_mid - pl%tp/2 + pl%pulse_offset)*au2fs
            print'(a,f10.4)', "  Flat end:    ", (pl%t_mid + pl%tp/2 + pl%pulse_offset)*au2fs
            print'(a,f10.4)', "  Pulse end:   ", (pl%t_mid + pl%tp/2 + TU_eff + pl%pulse_offset)*au2fs
            do k = 1, Nt
                ttime = time(k) - pl%t_mid - pl%pulse_offset
                pl%env(k)    = trapezoidal(time(k), pl%tp, pl%t_mid, TU_eff, pl%pulse_offset)
                pl%alpha_t(k) = pl%alpha0 * pl%env(k) * cos(pl%omega * ttime + pl%phi)
                pl%A_field(k) = trapezoidal_vector_pulse(time(k), pl%omega, pl%phi, &
                    & pl%alpha0, pl%tp, pl%t_mid, pl%pulse_offset, TU_eff)
                pl%E_field(k) = trapezoidal_electric_pulse(time(k), pl%omega, pl%phi, &
                    & pl%alpha0, pl%tp, pl%t_mid, pl%pulse_offset, TU_eff)
            end do
        case default
            print*, trim(label) // ": Default pulse shape is CW."
        end select
    end subroutine generate_single

    !> Calculate field spectra using FFTW.
    subroutine field_spectra(this)
        use global_vars, only: prop_par_FFTW, pulse_data_dir
        class(pulse_param), intent(inout) :: this
        integer :: k
        character(150) :: filename
        type(C_PTR) :: planTF, planTB, p_in, p_out
        complex(C_DOUBLE_COMPLEX), pointer:: E_dum_in(:), E_dum_out(:)
        integer:: field_spec_tk

        p_in = fftw_alloc_complex(int(Nt, C_SiZE_T))
        call c_f_pointer(p_in,E_dum_in,[Nt])
        p_out = fftw_alloc_complex(int(Nt, C_SiZE_T))
        call c_f_pointer(p_out,E_dum_out,[Nt])

        call fftw_initialize_threads
        print*, "FFTW plan creation for pulse spectra ..."
        call fftw_create_c2c_plans(E_dum_in, E_dum_out, Nt, &
            & planTF, planTB, prop_par_FFTW)
        print*, "Done setting up FFTW."

        E_dum_in = this%El
        call fftw_execute_dft(planTF,E_dum_in, E_dum_out)
        E_dum_in = E_dum_out/sqrt(dble(Nt))

        write(filename,fmt='(a,a)') adjustl(trim(pulse_data_dir)), 'field_spectra.out'
        open(newunit=field_spec_tk, file=filename,status="unknown")
        do k = Nt/2+1, Nt
            write(field_spec_tk,*) -(Nt + 1 - k) * 2 *pi/(dt * Nt), &
                & real(E_dum_in(k)), imag(E_dum_in(k)), abs(E_dum_in(k))
        end do
        do k = 1, Nt/2
            write(field_spec_tk,*) (k-1)*2*pi/(dt*Nt), real(E_dum_in(k)), &
                & imag(E_dum_in(k)), abs(E_dum_in(k))
        end do
        close(field_spec_tk)

        call fftw_destroy_plan(planTF)
        call fftw_destroy_plan(planTB)
        call fftw_free(p_in)
        call fftw_free(p_out)
    end subroutine field_spectra

    !> Write individual pulse data and total fields to files.
    subroutine write_pulse_to_file(this)
        use global_vars, only: pulse_data_dir
        class(pulse_param), intent(in) :: this
        integer :: i, k, unit_tk, tk_ef, tk_vf, tk_kf
        real(dp) :: A_check(Nt), alpha_t_check(Nt)
        character(150) :: filename
        character(8)   :: idx_str

        ! Per-pulse files
        do i = 1, this%N_pulses
            write(idx_str, '(I0)') i
            associate(pl => this%pulses(i))
                write(filename,fmt='(a,a,a,a)') adjustl(trim(pulse_data_dir)), &
                    & 'envelope', trim(idx_str), '.out'
                open(newunit=unit_tk, file=filename,status="unknown")
                do k = 1, Nt
                    write(unit_tk,*) time(k)*au2fs, pl%env(k)
                end do
                close(unit_tk)

                write(filename,fmt='(a,a,a,a)') adjustl(trim(pulse_data_dir)), &
                    & 'field', trim(idx_str), '.out'
                open(newunit=unit_tk, file=filename,status="unknown")
                do k = 1, Nt
                    write(unit_tk,*) time(k)*au2fs, pl%E_field(k), pl%A_field(k), pl%alpha_t(k)
                end do
                close(unit_tk)
            end associate
        end do

        ! Total-field files
        write(filename,fmt='(a,a)') adjustl(trim(pulse_data_dir)), &
            & 'Total_electric_field.out'
        open(newunit=tk_ef, file=filename,status="unknown")
        write(filename,fmt='(a,a)') adjustl(trim(pulse_data_dir)), &
            & 'Total_vector_field.out'
        open(newunit=tk_vf, file=filename,status="unknown")
        write(filename,fmt='(a,a)') adjustl(trim(pulse_data_dir)), &
            & 'Total_KH_field.out'
        open(newunit=tk_kf, file=filename,status="unknown")

        do k = 1, Nt
            A_check(k)      = -sum(this%El(1:k))*dt
            alpha_t_check(k) =  sum(this%Al(1:k))*dt
            write(tk_ef, *) time(k)*au2fs, this%El(k)
            write(tk_vf, *) time(k)*au2fs, this%Al(k), A_check(k)
            write(tk_kf, *) time(k)*au2fs, this%alpha_t(k), alpha_t_check(k)
        end do
        close(tk_ef)
        close(tk_vf)
        close(tk_kf)

        print*, "Done writing field information in the files."
    end subroutine write_pulse_to_file

    ! ==================================================================
    !  Envelope functions (unchanged)
    ! ==================================================================
    function cos2(t, tp, t_mid, pulse_offset)
        real(dp), intent(in) :: t, tp, t_mid, pulse_offset
        real(dp) :: cos2
        if (t > (t_mid+pulse_offset-tp/2) .and. t < (t_mid+pulse_offset+tp/2)) then
            cos2 = cos((t - t_mid - pulse_offset)*pi/tp)**2
        else
            cos2 = 0._dp
        end if
    end function cos2

    function sin2(t, tp, t_mid, pulse_offset)
        real(dp), intent(in) :: t, tp, t_mid, pulse_offset
        real(dp) :: sin2
        if (t > (t_mid+pulse_offset-tp/2) .and. t < (t_mid+pulse_offset+tp/2)) then
            sin2 = sin((t - t_mid - pulse_offset + tp/2)*pi/tp)**2
        else
            sin2 = 0._dp
        end if
    end function sin2

    function trapezoidal(t, tp, t_mid, rise_time, pulse_offset)
        real(dp), intent(in) :: t, tp, t_mid, rise_time, pulse_offset
        real(dp) :: trapezoidal, slope, yc, teff
        teff = t - pulse_offset
        if (teff >= t_mid - (tp/2 + rise_time) .and. teff <= t_mid - tp/2) then
            slope = 1._dp/rise_time
            yc = (t_mid - (tp/2 + rise_time)) * slope
            trapezoidal = slope * teff - yc
        elseif (teff > t_mid - tp/2 .and. teff <= t_mid + tp/2) then
            trapezoidal = 1._dp
        elseif (teff > t_mid + tp/2 .and. teff <= t_mid + (tp/2 + rise_time)) then
            slope = -1._dp/rise_time
            yc = (t_mid + tp/2 + rise_time) * slope
            trapezoidal = slope * teff - yc
        else
            trapezoidal = 0._dp
        end if
    end function trapezoidal

    function gaussian(t, tp, t_mid)
        implicit none
        real(dp):: t, tp, t_mid
        real(dp):: gaussian, fwhm
        fwhm = (4._dp * log(2._dp)) / tp**2
        gaussian = exp(-fwhm * (t - t_mid)**2)
    end function gaussian

    ! ==================================================================
    !  Analytic vector and electric fields (unchanged signatures)
    ! ==================================================================
    function gaussian_vector_pulse(t, tp, t_mid, alpha0, omega, phase, pulse_offset)
        implicit none
        real(dp), intent(in) :: t, tp, t_mid, alpha0, omega, phase, pulse_offset
        real(dp) :: gaussian_env, fwhm, theta
        real(dp) :: gaussian_vector_pulse
        theta = omega * (t-t_mid-pulse_offset) + phase
        fwhm = (4._dp * log(2._dp)) / tp**2
        gaussian_env = exp(-fwhm * (t - t_mid)**2)
        gaussian_vector_pulse = -alpha0 * gaussian_env * ((t-t_mid) * 2._dp * fwhm * cos(theta) &
            & + omega * sin(theta))
    end function gaussian_vector_pulse

    function gaussian_electric_pulse(t, tp, t_mid, alpha0, omega, phase, pulse_offset)
        implicit none
        real(dp), intent(in) :: t, tp, t_mid, alpha0, omega, phase, pulse_offset
        real(dp) :: gaussian_env, fwhm, theta
        real(dp) :: gaussian_electric_pulse
        theta = omega * (t-t_mid-pulse_offset) + phase
        fwhm = (4._dp * log(2._dp)) / tp**2
        gaussian_env = exp(-fwhm * (t - t_mid)**2)
        gaussian_electric_pulse = alpha0 * gaussian_env * (t*t * 4._dp *fwhm * fwhm  &
            & + 2._dp * fwhm + omega * omega) * cos(theta)
    end function gaussian_electric_pulse

    function sin2_vector_pulse(t, tp, t_mid, alpha0, omega, phase, pulse_offset)
        implicit none
        real(dp), intent(in) :: t, tp, t_mid, alpha0, omega, phase, pulse_offset
        real(dp) :: sin2_env, theta, t_local, inv_tp
        real(dp) :: sin2_vector_pulse
        if (t > (t_mid+pulse_offset-tp/2) .and. t < (t_mid+pulse_offset+tp/2)) then
            t_local = t - t_mid - pulse_offset + tp/2
            inv_tp = 1._dp/tp
            theta = omega * t_local + phase
            sin2_env = sin(pi * t_local * inv_tp)**2
            sin2_vector_pulse = alpha0 * ( pi* inv_tp * sin(2._dp * pi * t_local * inv_tp) * sin(theta) &
                & + omega * sin2_env * cos(theta))
        else
            sin2_vector_pulse = 0._dp
        end if
    end function sin2_vector_pulse

    function sin2_electric_pulse(t, tp, t_mid, alpha0, omega, phase, pulse_offset)
        implicit none
        real(dp), intent(in) :: t, tp, t_mid, alpha0, omega, phase, pulse_offset
        real(dp) :: sin2_env, theta, t_local, inv_tp
        real(dp) :: sin2_electric_pulse
        if (t > (t_mid+pulse_offset-tp/2) .and. t < (t_mid+pulse_offset+tp/2)) then
            t_local = t - t_mid - pulse_offset + tp/2
            inv_tp = 1._dp/tp
            theta = omega * t_local + phase
            sin2_env = sin(pi * t_local * inv_tp)**2
            sin2_electric_pulse = -alpha0 * ( 2._dp * pi*pi * inv_tp*inv_tp * cos(2._dp * pi * t_local * inv_tp) &
                & * sin(theta) + 2._dp * omega * pi * inv_tp * sin(2._dp * pi * t_local * inv_tp) * cos(theta) &
                & - omega * omega * sin2_env * sin(theta))
        else
            sin2_electric_pulse = 0._dp
        end if
    end function sin2_electric_pulse

    function trapezoidal_vector_pulse(t, omega, phase, alpha0, tp, t_mid, pulse_offset, rise_time)
        real(dp), intent(in) :: t, omega, phase, alpha0, tp, t_mid, pulse_offset, rise_time
        real(dp) :: trapezoidal_vector_pulse
        real(dp) :: t_start, t_local, theta, phi0
        real(dp) :: TU, TF
        TU = rise_time
        TF = tp + 2._dp * TU
        t_start = t_mid - tp/2._dp - TU + pulse_offset
        t_local = t - t_start
        phi0 = -omega * (tp/2._dp + TU) + phase
        theta = omega * t_local + phi0
        if (t_local >= 0._dp .and. t_local <= TU) then
            trapezoidal_vector_pulse = -alpha0 / TU &
                & * ( omega * t_local * sin(theta) - cos(theta) + cos(phi0) )
        elseif (t_local > TU .and. t_local <= TU + tp) then
            trapezoidal_vector_pulse = -alpha0 * omega &
                & * ( sin(theta) + (cos(phi0) - cos(omega * TU + phi0)) / (omega * TU) )
        elseif (t_local > TU + tp .and. t_local <= TF) then
            trapezoidal_vector_pulse = -alpha0 / TU &
                & * ( (TF - t_local) * omega * sin(theta) + cos(theta) - cos(omega * TF + phi0) )
        else
            trapezoidal_vector_pulse = 0._dp
        end if
    end function trapezoidal_vector_pulse

    function trapezoidal_electric_pulse(t, omega, phase, alpha0, tp, t_mid, pulse_offset, rise_time)
        real(dp), intent(in) :: t, omega, phase, alpha0, tp, t_mid, pulse_offset, rise_time
        real(dp) :: trapezoidal_electric_pulse
        real(dp) :: t_start, t_local, theta, phi0
        real(dp) :: TU, TF
        TU = rise_time
        TF = tp + 2._dp * TU
        t_start = t_mid - tp/2._dp - TU + pulse_offset
        t_local = t - t_start
        phi0 = -omega * (tp/2._dp + TU) + phase
        theta = omega * t_local + phi0
        if (t_local >= 0._dp .and. t_local <= TU) then
            trapezoidal_electric_pulse = alpha0 / TU &
                & * ( omega**2 * t_local * cos(theta) + 2._dp * omega * sin(theta) - 2._dp * omega * sin(phi0) )
        elseif (t_local > TU .and. t_local <= TU + tp) then
            trapezoidal_electric_pulse = alpha0 * omega**2 &
                & * ( cos(theta) + 2._dp * (sin(omega * TU + phi0) - sin(phi0)) / (omega * TU) )
        elseif (t_local > TU + tp .and. t_local <= TF) then
            trapezoidal_electric_pulse = alpha0 / TU &
                & * ( omega**2 * (TF - t_local) * cos(theta) - 2._dp * omega * sin(theta) &
                &     + 2._dp * omega * sin(omega * TF + phi0) )
        else
            trapezoidal_electric_pulse = 0._dp
        end if
    end function trapezoidal_electric_pulse

    ! ==================================================================
    !  Deallocation routines
    ! ==================================================================
    subroutine deallocate_envelope(this)
        class(pulse_param), intent(inout) :: this
        integer :: i
        if (allocated(this%pulses)) then
            do i = 1, size(this%pulses)
                if (allocated(this%pulses(i)%env)) deallocate(this%pulses(i)%env)
            end do
        end if
    end subroutine deallocate_envelope

    subroutine deallocate_field(this)
        class(pulse_param), intent(inout) :: this
        integer :: i
        if (allocated(this%pulses)) then
            do i = 1, size(this%pulses)
                if (allocated(this%pulses(i)%E_field))  deallocate(this%pulses(i)%E_field)
                if (allocated(this%pulses(i)%A_field))  deallocate(this%pulses(i)%A_field)
                if (allocated(this%pulses(i)%alpha_t))  deallocate(this%pulses(i)%alpha_t)
            end do
        end if
    end subroutine deallocate_field

    subroutine deallocate_all(this)
        class(pulse_param), intent(inout) :: this
        print*
        print*, "Cleaning up pulse variables ..."
        call deallocate_envelope(this)
        call deallocate_field(this)
        if (allocated(this%El))      deallocate(this%El)
        if (allocated(this%Al))      deallocate(this%Al)
        if (allocated(this%alpha_t)) deallocate(this%alpha_t)
        if (allocated(this%pulses))  deallocate(this%pulses)
        print*,"Done"
    end subroutine deallocate_all

end module pulse_mod