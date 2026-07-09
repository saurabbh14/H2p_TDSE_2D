!> This module contains the subroutine for generating the laser pulse
!> and the envelope functions.
module pulse_mod
    use global_vars, only: dp, Nt, output_data_dir, dt, time
    use data_au
    use FFTW3
    use omp_lib
    implicit none
    private
    public:: pulse_param

    type :: pulse_param
        character(150) :: envelope_shape_laser1, envelope_shape_laser2
        real(dp) :: tp1, fwhm, t_mid1, rise_time1
        real(dp) :: tp2, t_mid2, rise_time2
        real(dp) :: alpha01, alpha02, phi1, phi2
        real(dp) :: lambda1, lambda2
        real(dp) :: omega1, omega2
        real(dp) :: pulse_offset1, pulse_offset2
        real(dp), allocatable :: alpha_t(:)
        real(dp), allocatable :: alpha_t1(:), alpha_t2(:)
        real(dp), allocatable :: El(:), Al(:)
        real(dp), allocatable :: E21(:), E22(:)
        real(dp), allocatable :: A21(:), A22(:)
        real(dp), allocatable :: g1(:), g2(:)
  
    contains

        procedure :: read => read_pulse_params
        procedure :: initialize => initialize_pulse_param
        procedure :: param_print => print_pulse_param
        procedure :: generate => generate_pulse
        procedure :: write_to_file => write_pulse_to_file
        procedure :: deallocate_env => deallocate_envelope
        procedure :: deallocate_field => deallocate_field
        procedure :: spectra => field_spectra
        procedure :: deallocate_all
    end type pulse_param

contains
  
    ! Read pulse parameters from the input file
    subroutine read_pulse_params(this, input_path)
        class(pulse_param), intent(inout) :: this
        character(2000), intent(in) :: input_path
        ! file tokens
        integer :: input_tk

        ! Intermediate variables for pulse_param components
        character(150) :: envelope_shape_laser1, envelope_shape_laser2
        real(dp) :: lambda1, lambda2, tp1, tp2, t_mid1, t_mid2
        real(dp) :: alpha01, alpha02, phi1, phi2, rise_time1, rise_time2

        namelist /laser_param/envelope_shape_laser1, envelope_shape_laser2, &
        & lambda1, lambda2, tp1, tp2, t_mid1, t_mid2, alpha01, alpha02, & 
        & phi1, phi2, rise_time1, rise_time2

        open(newunit=input_tk, file=adjustl(trim(input_path)), status='old')
        read(input_tk,nml=laser_param)
        close(input_tk)

        ! Assign values to the pulse_param components
        this%envelope_shape_laser1 = envelope_shape_laser1
        this%envelope_shape_laser2 = envelope_shape_laser2
        this%lambda1 = lambda1
        this%lambda2 = lambda2
        this%tp1 = tp1
        this%tp2 = tp2
        this%t_mid1 = t_mid1
        this%t_mid2 = t_mid2
        this%alpha01 = alpha01
        this%alpha02 = alpha02
        this%phi1 = phi1
        this%phi2 = phi2
        this%rise_time1 = rise_time1
        this%rise_time2 = rise_time2
    end subroutine read_pulse_params
    
    ! A subroutine for printing the pulse parameters
    subroutine print_pulse_param(this)
        class(pulse_param), intent(in) :: this
        print*, "Laser parameters:"
        print*, "Laser #1:"
        print*, "Envelope shape:", trim(this%envelope_shape_laser1)
        print*, "Lambda:", this%lambda1, "nm"
        print*, "Quiver amplitude (alpha0):", this%alpha01, "a.u."
        print*, "Electric field strength (derived):", this%alpha01 * this%omega1**2, "a.u."
        print*, "Pulse width (tp):", this%tp1 * au2fs, "fs"
        print*, "Pulse midpoint:", this%t_mid1 * au2fs, "fs"
        print*, "phi1:", this%phi1, "pi"
        print*, "Rise time:", this%rise_time1 * au2fs, "fs"
        print*, "Laser #2:"
        print*, "Envelope shape:", trim(this%envelope_shape_laser2)
        print*, "Lambda:", this%lambda2, "nm"
        print*, "Quiver amplitude (alpha0):", this%alpha02, "a.u."
        print*, "Electric field strength (derived):", this%alpha02 * this%omega2**2, "a.u."
        print*, "Pulse width (tp):", this%tp2 * au2fs, "fs"
        print*, "Pulse midpoint:", this%t_mid2 * au2fs, "fs"
        print*, "phi2:", this%phi2, "pi"
        print*, "Rise time:", this%rise_time2 * au2fs, "fs"
        print*, "------------------------------------------------------"
        print*, "Final pulse parameters:"
        print*, "Wavelength 1 =", sngl(this%lambda1), "nm"
        print*, "Phase 1 =", sngl(this%phi1)
        print*, "Quiver amplitude 1 =", sngl(this%alpha01), "a.u."
        print*, "Field strength 1 =", sngl(this%alpha01 * this%omega1**2), "a.u."
        print*, "Intensity 1 =", sngl((this%alpha01 * this%omega1**2)**2 * 3.509e16_dp), "W/cm2"
        print*, "Wavelength 2 =", sngl(this%lambda2), "nm"
        print*, "Phase 2 =", sngl(this%phi2)
        print*, "Quiver amplitude 2 =", sngl(this%alpha02), "a.u."
        print*, "Field strength 2 =", sngl(this%alpha02 * this%omega2**2), "a.u."
        print*, "Intensity 2 =", sngl((this%alpha02 * this%omega2**2)**2 * 3.509e16_dp), "W/cm2"
        print*, "Wave duration =", sngl(this%tp1*au2fs), "fs"
        print*, "------------------------------------------------------"
    end subroutine print_pulse_param

    ! A subroutine for initializing the pulse parameters
    subroutine initialize_pulse_param(this)
        class(pulse_param), intent(inout) :: this
        ! Initialize the pulse parameters
        this%tp1 = this%tp1 / au2fs  
        this%tp2 = this%tp2 / au2fs  
        this%t_mid1 = this%t_mid1 / au2fs   
        this%t_mid2 = this%t_mid2 / au2fs   
        this%rise_time1 = this%rise_time1 / au2fs
        this%rise_time2 = this%rise_time2 / au2fs
        this%omega1 = (1._dp / (this%lambda1 * 1.e-7_dp)) * cm2au
        this%omega2 = (1._dp / (this%lambda2 * 1.e-7_dp)) * cm2au
        this%phi1 = this%phi1 * pi
        this%phi2 = this%phi2 * pi
    end subroutine initialize_pulse_param


    ! A subroutine for defining the field 
    subroutine generate_pulse(this)
        use differentiation, only: central_diff_on_grid
        class(pulse_param), intent(inout) :: this
        integer :: k
        integer :: n_cycles
        real(dp) :: ttime, TU_eff, tp_eff

        print*
        print*, "Pulse generation..."
    
        ! Initialize arrays
        allocate(this%alpha_t(Nt))
        allocate(this%alpha_t1(Nt), this%alpha_t2(Nt))
        allocate(this%El(Nt), this%Al(Nt))
        allocate(this%E21(Nt), this%E22(Nt))
        allocate(this%A21(Nt), this%A22(Nt))
        allocate(this%g1(Nt), this%g2(Nt))
        this%alpha_t = 0.0_dp
        this%El = 0.0_dp
        this%Al = 0.0_dp
        this%alpha_t1 = 0.0_dp; this%alpha_t2 = 0._dp
        this%E21 = 0.0_dp; this%E22 = 0.0_dp
        this%A21 = 0.0_dp; this%A22 = 0.0_dp
        this%g1 = 0.0_dp; this%g2 = 0.0_dp
        this%pulse_offset1 = 0.0_dp
        this%pulse_offset2 = 0.0_dp

        ! Calculate the envelope shapes
        ! Envelope shape for laser 1
        select case(trim(this%envelope_shape_laser1))
        case("cos2")
            do k = 1, Nt
                this%g1(k) = cos2(time(k), this%tp1, this%t_mid1, this%pulse_offset1)
                this%alpha_t1(k) = this%alpha01 * this%g1(k) * cos(this%omega1 * (time(k) - this%t_mid1 &
                  & - this%pulse_offset1) + this%phi1)  
            enddo
            call central_diff_on_grid(this%alpha_t1, Nt, dt, this%A21)
            call central_diff_on_grid(this%A21, Nt, dt, this%E21)
            this%E21 = -this%E21
        
        case("sin2")
            n_cycles = int(this%tp1 * this%omega1 / (2._dp * pi)) + 1
            n_cycles = max(n_cycles, 1)
            tp_eff = 2._dp * pi * real(n_cycles, dp) / this%omega1
            do k = 1, Nt
                this%g1(k) = sin2(time(k), tp_eff, this%t_mid1, this%pulse_offset1)
                this%alpha_t1(k) = this%alpha01 * this%g1(k) * sin(this%omega1 * (time(k) - this%t_mid1 &
                    & - this%pulse_offset1 - tp_eff/2) + this%phi1)  
                this%A21(k) = sin2_vector_pulse(time(k), tp_eff, this%t_mid1, this%alpha01, this%omega1, &
                    & this%phi1, this%pulse_offset1)
                this%E21(k) = sin2_electric_pulse(time(k), tp_eff, this%t_mid1, this%alpha01, this%omega1, &
                    & this%phi1, this%pulse_offset1)
            enddo
            !call central_diff_on_grid(this%alpha_t1, Nt, dt, this%A21)
            !call central_diff_on_grid(this%A21, Nt, dt, this%E21)
            !this%E21 = -this%E21
        
        case("gaussian")
            do k = 1, Nt
                this%g1(k) = gaussian(time(k), this%tp1, this%t_mid1)
                this%alpha_t1(k) = this%alpha01 * this%g1(k) * cos(this%omega1 * (time(k) - this%t_mid1 &
                    & - this%pulse_offset1) + this%phi1)  
                this%A21(k) = gaussian_vector_pulse(time(k), this%tp1, this%t_mid1, this%alpha01, this%omega1, &
                    & this%phi1, this%pulse_offset1)
                this%E21(k) = gaussian_electric_pulse(time(k), this%tp1, this%t_mid1, this%alpha01, this%omega1, &
                    & this%phi1, this%pulse_offset1)
            enddo
            !call central_diff_on_grid(this%alpha_t1, Nt, dt, this%A21)
            !call central_diff_on_grid(this%A21, Nt, dt, this%E21)
            !this%E21 = -this%E21
        
        case("trapezoidal")
            n_cycles = int(this%rise_time1 * this%omega1 / (2._dp * pi)) + 1
            n_cycles = max(n_cycles, 1)
            TU_eff = 2._dp * pi * real(n_cycles, dp) / this%omega1
            print*, "Laser1 trapezoidal boundary times (fs):"
            print'(a,f10.4)', "  Pulse start: ", (this%t_mid1 - this%tp1/2 - TU_eff + this%pulse_offset1)*au2fs
            print'(a,f10.4)', "  Rise end:    ", (this%t_mid1 - this%tp1/2 + this%pulse_offset1)*au2fs
            print'(a,f10.4)', "  Flat end:    ", (this%t_mid1 + this%tp1/2 + this%pulse_offset1)*au2fs
            print'(a,f10.4)', "  Pulse end:   ", (this%t_mid1 + this%tp1/2 + TU_eff + this%pulse_offset1)*au2fs
            do k = 1, Nt
                ttime = time(k) - this%t_mid1 - this%pulse_offset1 
                this%g1(k) = trapezoidal(time(k), this%tp1, this%t_mid1, TU_eff, this%pulse_offset1)
                this%alpha_t1(k) = this%alpha01 * this%g1(k) * cos(this%omega1 * ttime + this%phi1)
                this%A21(k) = trapezoidal_vector_pulse(time(K), this%omega1, this%phi1, &
                    & this%alpha01, this%tp1, this%t_mid1, this%pulse_offset1, TU_eff) 
                this%E21(k) = trapezoidal_electric_pulse(time(K), this%omega1, this%phi1, &
                    & this%alpha01, this%tp1, this%t_mid1, this%pulse_offset1, TU_eff) 
            enddo
        case default
            print*, "Laser1: Default pulse shape is CW."
        end select
        ! Envelope shape for laser 2
        select case(trim(this%envelope_shape_laser2))
        case("cos2")
            !this%tp1 = this%tp1/(1-2/pi) ! check this
            do k = 1, Nt 
                this%g2(k) = cos2(time(k), this%tp2, this%t_mid2, this%pulse_offset2)
                this%alpha_t2(k) = this%alpha02 * this%g2(k) * cos(this%omega2 * (time(k) - this%t_mid2 &
                  & - this%pulse_offset2) + this%phi2)                
            enddo
            call central_diff_on_grid(this%alpha_t2, Nt, dt, this%A22)
            call central_diff_on_grid(this%A22, Nt, dt, this%E22)
            this%E22 = -this%E22
        
        case("sin2")
            n_cycles = int(this%tp2 * this%omega2 / (2._dp * pi)) + 1
            n_cycles = max(n_cycles, 1)
            tp_eff = 2._dp * pi * real(n_cycles, dp) / this%omega2
            do k = 1, Nt
                this%g2(k) = sin2(time(k), tp_eff, this%t_mid2, this%pulse_offset2)
                this%alpha_t2(k) = this%alpha02 * this%g2(k) * sin(this%omega2 * (time(k) - this%t_mid2 &
                  & - this%pulse_offset2 + tp_eff/2) + this%phi2)
                this%A22(k) = sin2_vector_pulse(time(k), tp_eff, this%t_mid2, this%alpha02, this%omega2, &
                    & this%phi2, this%pulse_offset2)
                this%E22(k) = sin2_electric_pulse(time(k), tp_eff, this%t_mid2, this%alpha02, this%omega2, &
                    & this%phi2, this%pulse_offset2)  
            enddo
            !call central_diff_on_grid(this%alpha_t2, Nt, dt, this%A22)
            !call central_diff_on_grid(this%A22, Nt, dt, this%E22)
            !this%E22 = -this%E22
        
        case("gaussian")
            do k = 1, Nt
                this%g2(k) = gaussian(time(k), this%tp2, this%t_mid2)
                this%alpha_t2(k) = this%alpha02 * this%g2(k) * cos(this%omega2 * (time(k) - this%t_mid2 &
                  & - this%pulse_offset2) + this%phi2)
                this%A22(k) = gaussian_vector_pulse(time(k), this%tp2, this%t_mid2, this%alpha02, this%omega2, &
                    & this%phi2, this%pulse_offset2)
                this%E22(k) = gaussian_electric_pulse(time(k), this%tp2, this%t_mid2, this%alpha02, this%omega2, &
                    & this%phi2, this%pulse_offset2)
            enddo
            !call central_diff_on_grid(this%alpha_t2, Nt, dt, this%A22)
            !call central_diff_on_grid(this%A22, Nt, dt, this%E22)
            !this%E22 = -this%E22
        
        case("trapezoidal")
            n_cycles = int(this%rise_time2 * this%omega2 / (2._dp * pi)) + 1
            n_cycles = max(n_cycles, 1)
            TU_eff = 2._dp * pi * real(n_cycles, dp) / this%omega2
            print*, "Laser2 trapezoidal boundary times (fs):"
            print'(a,f10.4)', "  Pulse start: ", (this%t_mid2 - this%tp2/2 - TU_eff + this%pulse_offset2)*au2fs
            print'(a,f10.4)', "  Rise end:    ", (this%t_mid2 - this%tp2/2 + this%pulse_offset2)*au2fs
            print'(a,f10.4)', "  Flat end:    ", (this%t_mid2 + this%tp2/2 + this%pulse_offset2)*au2fs
            print'(a,f10.4)', "  Pulse end:   ", (this%t_mid2 + this%tp2/2 + TU_eff + this%pulse_offset2)*au2fs
            do k = 1, Nt
                this%g2(k) = trapezoidal(time(k), this%tp2, this%t_mid2, TU_eff, this%pulse_offset2)
                this%alpha_t2(k) = this%alpha02 * this%g2(k) * cos(this%omega2 * (time(k) - this%t_mid2 &
                  & - this%pulse_offset2) + this%phi2)                
                this%A22(k) = trapezoidal_vector_pulse(time(K), this%omega2, this%phi2, &
                    & this%alpha02, this%tp2, this%t_mid2, this%pulse_offset2, TU_eff)
                this%E22(k) = trapezoidal_electric_pulse(time(K), this%omega2, this%phi2, &
                    & this%alpha02, this%tp2, this%t_mid2, this%pulse_offset2, TU_eff)
            enddo
        case default
            print*, "Laser2: Default pulse shape is CW."
        end select

        ! Generate the Quiver-field
        this%alpha_t = this%alpha_t1 + this%alpha_t2
        
        ! Generate Vector field A = d(alpha_t) / dt
        this%Al = this%A21 + this%A22
        
        !Generate Electric field E = - dA / dt
        this%El = this%E21 + this%E22

        print*, "Pulse generation complete."
        call this%write_to_file()
    end subroutine generate_pulse

    ! A subroutine to calculate the field spectra using FFTW
    subroutine field_spectra(this)
        use global_vars, only: prop_par_FFTW, pulse_data_dir
        class(pulse_param), intent(inout) :: this
        integer :: k
        character(150) :: filename
        type(C_PTR) :: planTF, planTB, p_in, p_out
        complex(C_DOUBLE_COMPLEX), pointer:: E_dum_in(:), E_dum_out(:)
        ! file tokens
        integer:: field_spec_tk

        ! Creating aligned memory for FFTW
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
        ! write the field spectra to file
        do k = Nt/2+1, Nt
            write(field_spec_tk,*) -(Nt + 1 - k) * 2 *pi/(dt * Nt), & 
                & real(E_dum_in(k)), imag(E_dum_in(k)), abs(E_dum_in(k))
        enddo
        do k = 1, Nt/2
            write(field_spec_tk,*) (k-1)*2*pi/(dt*Nt), real(E_dum_in(k)), & 
                & imag(E_dum_in(k)), abs(E_dum_in(k))
        enddo
        close(field_spec_tk)

        call fftw_destroy_plan(planTF)
        call fftw_destroy_plan(planTB)
        call fftw_free(p_in)
        call fftw_free(p_out)
    end subroutine field_spectra

    ! A subroutine for writing the pulse to files
    subroutine write_pulse_to_file(this)
        use global_vars, only: pulse_data_dir
        class(pulse_param), intent(in) :: this
        integer :: k
        real(dp) :: A_check(Nt), alpha_t_check(Nt)
        character(150) :: filename
 
        ! file tokens
        integer:: envelope1_tk, envelope2_tk
        integer:: field1_tk, field2_tk
        integer:: elec_field_tk, vec_field_tk, kh_field_tk
    
        write(filename,fmt='(a,a)') adjustl(trim(pulse_data_dir)), 'envelope1.out'
        open(newunit=envelope1_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a)') adjustl(trim(pulse_data_dir)), 'envelope2.out'
        open(newunit=envelope2_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a,f4.2,a,i0,a)') adjustl(trim(pulse_data_dir)), &
            & 'electric_field1_alpha', this%alpha01,'_width',Int(this%tp1*au2fs),'.out'
        open(newunit=field1_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a,f6.4,a,i0,a)') adjustl(trim(pulse_data_dir)), &
            & 'electric_field2_alpha', this%alpha02,'_width',Int(this%tp2*au2fs),'.out'
        open(newunit=field2_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a,f4.2,a)') adjustl(trim(pulse_data_dir)), &
            & 'Total_electric_field_phi', this%phi2/pi,'pi.out'
        open(newunit=elec_field_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a,f4.2,a)') adjustl(trim(pulse_data_dir)), &
            & 'Total_vector_field_phi', this%phi2/pi, 'pi.out'
        open(newunit=vec_field_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a,f4.2,a)') adjustl(trim(pulse_data_dir)), &
            & 'Total_KH_field_phi', this%phi2/pi, 'pi.out'
        open(newunit=kh_field_tk, file=filename,status="unknown")

        timeloop: do k = 1, Nt
            ! Check elctric field correctness by reverse generating the vector potential
            A_check(k) = -sum(this%El(1:k))*dt
            alpha_t_check(k) = sum(this%Al(1:k))*dt
            write(field1_tk,*) time(k)*au2fs, this%E21(k), this%A21(k), this%alpha_t1(k)
            write(field2_tk,*) time(k)*au2fs, this%E22(k), this%A22(k), this%alpha_t2(k)
            write(envelope1_tk,*) time(k)*au2fs, this%g1(k)
            write(envelope2_tk,*) time(k)*au2fs, this%g2(k)
            write(elec_field_tk,*) time(k)*au2fs, this%El(k)
            write(vec_field_tk,*) time(k)*au2fs, this%Al(k), A_check(k)
            write(kh_field_tk,*) time(k)*au2fs, this%alpha_t(k), alpha_t_check(k)
        enddo timeloop

        print*, "Done writing field information in the files."
        close(envelope1_tk)
        close(envelope2_tk)
        close(field1_tk)
        close(field2_tk)
        close(elec_field_tk)
        close(vec_field_tk)
        close(kh_field_tk)

    end subroutine write_pulse_to_file
  
    ! pulse envelope functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cos2(time, tp, t_mid, pulse_offset)
        real(dp), intent(in) :: time, tp, t_mid, pulse_offset
        real(dp) :: cos2
        if (time .gt. (t_mid+pulse_offset-tp/2) .and. time .lt. (t_mid+pulse_offset+tp/2)) then
            cos2 = cos((time - t_mid-pulse_offset)*pi/tp)**2      
        else
            cos2 = 0._dp
        endif
    end function cos2

    function sin2(time, tp, t_mid, pulse_offset)
        real(dp), intent(in) :: time, tp, t_mid, pulse_offset
        real(dp) :: sin2
        if (time .gt. (t_mid+pulse_offset-tp/2) .and. time .lt. (t_mid+pulse_offset+tp/2)) then
            sin2 = sin((time - t_mid-pulse_offset+tp/2)*pi/tp)**2      
        else
            sin2 = 0._dp
        endif
    end function sin2

    function trapezoidal(time, tp, t_mid, rise_time, pulse_offset)
        real(dp), intent(in) :: time, tp, t_mid, rise_time, pulse_offset
        real(dp) :: trapezoidal, slope, yc, teff
        teff = time - pulse_offset
        if (teff .ge. t_mid - (tp/2 + rise_time) .and. teff .le. t_mid - tp/2) then
            slope = 1._dp/rise_time
            yc = (t_mid - (tp/2 + rise_time)) * slope
            trapezoidal = slope * teff - yc
        elseif (teff .gt. t_mid - tp/2 .and. teff .le. t_mid + tp/2) then
            trapezoidal = 1._dp
        elseif (teff .gt. t_mid + tp/2 .and. teff .le. t_mid + (tp/2 + rise_time)) then 
            slope = -1._dp/rise_time
            yc = (t_mid + tp/2 + rise_time) * slope
            trapezoidal = slope * teff - yc
        else
            trapezoidal = 0._dp
        endif
    end function trapezoidal
 
    function gaussian(time, tp, t_mid)
        implicit none
        real(dp):: time, tp, t_mid
        real(dp):: gaussian, fwhm
 
        fwhm = (4._dp * log(2._dp)) / tp**2
        gaussian = exp(-fwhm * (time - t_mid)**2)
    end function gaussian

    ! Analytic vector and electric fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function gaussian_vector_pulse(time, tp, t_mid, alpha0, omega, phase, pulse_offset)
        implicit none
        real(dp), intent(in) :: time, tp, t_mid, alpha0, omega, phase, pulse_offset
        real(dp) :: gaussian, fwhm, theta
        real(dp) :: gaussian_vector_pulse
        theta = omega * (time-t_mid-pulse_offset) + phase
        fwhm = (4._dp * log(2._dp)) / tp**2
        gaussian = exp(-fwhm * (time - t_mid)**2)
        gaussian_vector_pulse = -alpha0 * gaussian * ((time-t_mid) * 2._dp * fwhm * cos(theta) &
            & + omega * sin(theta))
    end function gaussian_vector_pulse

    function gaussian_electric_pulse(time, tp, t_mid, alpha0, omega, phase, pulse_offset)
        implicit none
        real(dp), intent(in) :: time, tp, t_mid, alpha0, omega, phase, pulse_offset
        real(dp) :: gaussian, fwhm, theta, t
        real(dp) :: gaussian_electric_pulse
        t = time - t_mid
        theta = omega * (time-t_mid-pulse_offset) + phase
        fwhm = (4._dp * log(2._dp)) / tp**2
        gaussian = exp(-fwhm * t*t)
        gaussian_electric_pulse = alpha0 * gaussian * (t*t * 4._dp *fwhm * fwhm  &
            & + 2._dp * fwhm + omega * omega) * cos(theta)
    end function gaussian_electric_pulse

    function sin2_vector_pulse(time, tp, t_mid, alpha0, omega, phase, pulse_offset)
        implicit none
        real(dp), intent(in) :: time, tp, t_mid, alpha0, omega, phase, pulse_offset
        real(dp) :: sin2, theta, t_local, inv_tp
        real(dp) :: sin2_vector_pulse
        if (time .gt. (t_mid+pulse_offset-tp/2) .and. time .lt. (t_mid+pulse_offset+tp/2)) then
            t_local = time - t_mid - pulse_offset + tp/2
            inv_tp = 1._dp/tp
            theta = omega * t_local + phase
            sin2 = sin(pi * t_local * inv_tp)**2

            sin2_vector_pulse = alpha0 * ( pi* inv_tp * sin(2._dp * pi * t_local * inv_tp) * sin(theta) &
                & + omega * sin2 * cos(theta))
        else
            sin2_vector_pulse = 0._dp
        endif
    end function sin2_vector_pulse

    function sin2_electric_pulse(time, tp, t_mid, alpha0, omega, phase, pulse_offset)
        implicit none
        real(dp), intent(in) :: time, tp, t_mid, alpha0, omega, phase, pulse_offset
        real(dp) :: sin2, theta, t_local, inv_tp
        real(dp) :: sin2_electric_pulse
        if (time .gt. (t_mid+pulse_offset-tp/2) .and. time .lt. (t_mid+pulse_offset+tp/2)) then
            t_local = time - t_mid - pulse_offset + tp/2
            inv_tp = 1._dp/tp
            theta = omega * t_local + phase
            sin2 = sin(pi * t_local * inv_tp)**2

            sin2_electric_pulse = -alpha0 * ( 2._dp * pi*pi * inv_tp*inv_tp * cos(2._dp * pi * t_local * inv_tp) & 
                & * sin(theta) + 2._dp * omega * pi * inv_tp * sin(2._dp * pi * t_local * inv_tp) * cos(theta) &
                & - omega * omega * sin2 * sin(theta))
        else
            sin2_electric_pulse = 0._dp
        endif
    end function sin2_electric_pulse

    !> Analytic vector potential A(t) = d(alpha_t)/dt with boundary-matched
    !> integration constants so that A(t=0)=0 and A(t=TF)=0.
    !> The caller must pass an integer-cycle-adjusted rise_time (TU) for
    !> internal boundary matching between Case II and Case III.
    function trapezoidal_vector_pulse(time, omega, phase, alpha0, tp, t_mid, pulse_offset, rise_time)
        real(dp), intent(in) :: time, omega, phase, alpha0, tp, t_mid, pulse_offset, rise_time
        real(dp) :: trapezoidal_vector_pulse
        real(dp) :: t_start, t_local, theta, phi0
        real(dp) :: TU, TF

        TU = rise_time
        TF = tp + 2._dp * TU
        t_start = t_mid - tp/2._dp - TU + pulse_offset
        t_local = time - t_start

        ! Initial phase at pulse start
        phi0 = -omega * (tp/2._dp + TU) + phase
        theta = omega * t_local + phi0

        if (t_local .ge. 0._dp .and. t_local .le. TU) then
            ! Case I: Rise
            trapezoidal_vector_pulse = -alpha0 / TU &
                & * ( omega * t_local * sin(theta) - cos(theta) + cos(phi0) )
        elseif (t_local .gt. TU .and. t_local .le. TU + tp) then
            ! Case II: Flat top
            trapezoidal_vector_pulse = -alpha0 * omega &
                & * ( sin(theta) + (cos(phi0) - cos(omega * TU + phi0)) / (omega * TU) )
        elseif (t_local .gt. TU + tp .and. t_local .le. TF) then
            ! Case III: Fall
            trapezoidal_vector_pulse = -alpha0 / TU &
                & * ( (TF - t_local) * omega * sin(theta) + cos(theta) - cos(omega * TF + phi0) )
        else
            trapezoidal_vector_pulse = 0._dp
        endif
    end function trapezoidal_vector_pulse

    !> Analytic electric field E(t) = -dA/dt with boundary-matched
    !> integration constants so that E(t=0)=0 and E(t=TF)=0.
    !> The caller must pass an integer-cycle-adjusted rise_time (TU) for
    !> internal boundary matching between Case II and Case III.
    function trapezoidal_electric_pulse(time, omega, phase, alpha0, tp, t_mid, pulse_offset, rise_time)
        real(dp), intent(in) :: time, omega, phase, alpha0, tp, t_mid, pulse_offset, rise_time
        real(dp) :: trapezoidal_electric_pulse
        real(dp) :: t_start, t_local, theta, phi0
        real(dp) :: TU, TF

        TU = rise_time
        TF = tp + 2._dp * TU
        t_start = t_mid - tp/2._dp - TU + pulse_offset
        t_local = time - t_start

        ! Initial phase at pulse start
        phi0 = -omega * (tp/2._dp + TU) + phase
        theta = omega * t_local + phi0

        if (t_local .ge. 0._dp .and. t_local .le. TU) then
            ! Case I: Rise
            trapezoidal_electric_pulse = alpha0 / TU &
                & * ( omega**2 * t_local * cos(theta) + 2._dp * omega * sin(theta) - 2._dp * omega * sin(phi0) )
        elseif (t_local .gt. TU .and. t_local .le. TU + tp) then
            ! Case II: Flat top
            trapezoidal_electric_pulse = alpha0 * omega**2 &
                & * ( cos(theta) + 2._dp * (sin(omega * TU + phi0) - sin(phi0)) / (omega * TU) )
        elseif (t_local .gt. TU + tp .and. t_local .le. TF) then
            ! Case III: Fall
            trapezoidal_electric_pulse = alpha0 / TU &
                & * ( omega**2 * (TF - t_local) * cos(theta) - 2._dp * omega * sin(theta) &
                &     + 2._dp * omega * sin(omega * TF + phi0) )
        else
            trapezoidal_electric_pulse = 0._dp
        endif
    end function trapezoidal_electric_pulse

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! A subroutine for deallocating the envelope arrays
    subroutine deallocate_envelope(this)
        class(pulse_param), intent(inout) :: this
        if (allocated(this%g1)) then
            deallocate(this%g1)
        end if
        if (allocated(this%g2)) then
            deallocate(this%g2)
        end if
    end subroutine deallocate_envelope
    ! A subroutine for deallocating the field arrays
    subroutine deallocate_field(this)
        class(pulse_param), intent(inout) :: this
        if (allocated(this%E21)) then
            deallocate(this%E21)
        end if
        if (allocated(this%E22)) then
            deallocate(this%E22)
        end if
        if (allocated(this%A21)) then
            deallocate(this%A21)
        end if
        if (allocated(this%A22)) then
            deallocate(this%A22)
        end if
        if (allocated(this%alpha_t1)) then
            deallocate(this%alpha_t1)
        end if
        if (allocated(this%alpha_t2)) then
            deallocate(this%alpha_t2)
        end if
    end subroutine deallocate_field
    ! A subroutine for deallocating all arrays
    subroutine deallocate_all(this)
    !    type(pulse_param) :: this
        class(pulse_param), intent(inout) :: this

        print*
        print*, "Cleaning up pulse variables ..."
        call deallocate_envelope(this)
        call deallocate_field(this)
        if (allocated(this%El)) then
            deallocate(this%El)
        end if
        if (allocated(this%Al)) then
            deallocate(this%Al)
        end if
        if (allocated(this%alpha_t)) then
            deallocate(this%alpha_t)
        end if
        print*,"Done"
    end subroutine deallocate_all
  !------------------------------------------------------------------------------

end module pulse_mod
 

