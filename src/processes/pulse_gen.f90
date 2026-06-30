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
        real(dp) :: e01, e02, phi1, phi2
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
        real(dp) :: E01, E02, phi1, phi2, rise_time1, rise_time2

        namelist /laser_param/envelope_shape_laser1, envelope_shape_laser2, &
        & lambda1, lambda2, tp1, tp2, t_mid1, t_mid2, E01, E02, & 
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
        this%E01 = E01
        this%E02 = E02
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
        print*, "Electric field strength:", this%E01, "a.u."
        print*, "Pulse width (tp):", this%tp1 * au2fs, "fs"
        print*, "Pulse midpoint:", this%t_mid1 * au2fs, "fs"
        print*, "phi1:", this%phi1, "pi"
        print*, "Rise time:", this%rise_time1 * au2fs, "fs"
        print*, "Laser #2:"
        print*, "Envelope shape:", trim(this%envelope_shape_laser2)
        print*, "Lambda:", this%lambda2, "nm"
        print*, "Electric field strength:", this%E02, "a.u."
        print*, "Pulse width (tp):", this%tp2 * au2fs, "fs"
        print*, "Pulse midpoint:", this%t_mid2 * au2fs, "fs"
        print*, "phi2:", this%phi2, "pi"
        print*, "Rise time:", this%rise_time2 * au2fs, "fs"
        print*, "------------------------------------------------------"
        print*, "Final pulse parameters:"
        print*, "Wavelength 1 =", sngl(this%lambda1), "nm"
        print*, "Phase 1 =", sngl(this%phi1)
        print*, "Field strength =", sngl(this%e01), "a.u.", sngl(this%e01*e02au), "V/m"
        print*, "Intensity =", sngl(this%e01**2*3.509e16_dp), "W/cm2"
        print*, "Wavelength 2 =", sngl(this%lambda2), "nm"
        print*, "Phase 2 =", sngl(this%phi2)
        print*, "Field strength =", sngl(this%e02), "a.u.", sngl(this%e02*e02au), "V/m"
        print*, "Intensity =", sngl(this%e02**2*3.509e16_dp), "W/cm2"
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
        real(dp) :: A01, A02
        real(dp) :: ttime

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

        ! Calculate amplitudes
        A01 = this%e01 / this%omega1
        A02 = this%e02 / this%omega2
        ! Calculate the envelope shapes
        ! Envelope shape for laser 1
        select case(trim(this%envelope_shape_laser1))
        case("cos2")
            !this%tp1 = this%tp1/(1-2/pi) ! check this
            do k = 1, Nt
                this%g1(k) = cos2(time(k), this%tp1, this%t_mid1, this%pulse_offset1)
                this%alpha_t1(k) = this%E01 * this%g1(k) * cos(this%omega1 * (time(k) - this%t_mid1 &
                  & - this%pulse_offset1) + this%phi1)  
            enddo
            call central_diff_on_grid(this%alpha_t1, Nt, dt, this%A21)
            call central_diff_on_grid(this%A21, Nt, dt, this%E21)
            this%E21 = -this%E21
        case("gaussian")
            do k = 1, Nt
                this%g1(k) = gaussian(time(k), this%tp1, this%t_mid1)
                this%alpha_t1(k) = this%E01 * this%g1(k) * cos(this%omega1 * (time(k) - this%t_mid1 &
                  & - this%pulse_offset1) + this%phi1)  
            enddo
            call central_diff_on_grid(this%alpha_t1, Nt, dt, this%A21)
            call central_diff_on_grid(this%A21, Nt, dt, this%E21)
            this%E21 = -this%E21
        case("trapezoidal")
            do k = 1, Nt
                ttime = time(k) - this%t_mid1 - this%pulse_offset1 
                this%g1(k) = trapezoidal(time(k), this%tp1, this%t_mid1, this%rise_time1)
                this%alpha_t1(k) = this%E01 * this%g1(k) * cos(this%omega1 * ttime + this%phi1)
                this%A21(k) = trapezoidal_vector_pulse(time(K), this%omega1, this%phi1, &
                    & this%E01, this%tp1, this%t_mid1, this%pulse_offset1, this%rise_time1) 
                this%E21(k) = trapezoidal_electric_pulse(time(K), this%omega1, this%phi1, &
                    & this%E01, this%tp1, this%t_mid1, this%pulse_offset1, this%rise_time1) 
            enddo
        case default
            print*, "Laser1: Default pulse shape is CW."
        end select
        ! Envelope shape for laser 2
        select case(trim(this%envelope_shape_laser2))
        case("cos2")
            !this%tp1 = this%tp1/(1-2/pi) ! check this
            do k = 1, Nt 
                this%g2(k) = cos2(time(k), this%tp1, this%t_mid2, this%pulse_offset2)
                this%alpha_t2(k) = this%E02 * this%g2(k) * cos(this%omega2 * (time(k) - this%t_mid2 &
                  & - this%pulse_offset2) + this%phi2)                
            enddo
            call central_diff_on_grid(this%alpha_t2, Nt, dt, this%A22)
            call central_diff_on_grid(this%A22, Nt, dt, this%E22)
            this%E22 = -this%E22
        case("gaussian")
            do k = 1, Nt
                this%g2(k) = gaussian(time(k), this%tp2, this%t_mid2)
                this%alpha_t2(k) = this%E02 * this%g2(k) * cos(this%omega2 * (time(k) - this%t_mid2 &
                  & - this%pulse_offset2) + this%phi2)
            enddo
            call central_diff_on_grid(this%alpha_t2, Nt, dt, this%A22)
            call central_diff_on_grid(this%A22, Nt, dt, this%E22)
            this%E22 = -this%E22
        case("trapezoidal")
            do k = 1, Nt
                this%g2(k) = trapezoidal(time(k), this%tp2, this%t_mid2, this%rise_time2)
                this%alpha_t2(k) = this%E02 * this%g2(k) * cos(this%omega2 * (time(k) - this%t_mid2 &
                  & - this%pulse_offset2) + this%phi2)                
            enddo
            this%A22 = 0._dp
            this%E22 = 0._dp
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
            & 'electric_field1_E', this%E01,'_width',Int(this%tp1*au2fs),'.out'
        open(newunit=field1_tk, file=filename,status="unknown")
        write(filename,fmt='(a,a,f6.4,a,i0,a)') adjustl(trim(pulse_data_dir)), &
            & 'electric_field2_E', this%E02,'_width',Int(this%tp2*au2fs),'.out'
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
            write(field1_tk,*) time(k)*au2fs, this%E21(k), this%A21(k), this%alpha_t1(k)
            write(field2_tk,*) time(k)*au2fs, this%E22(k), this%A22(k), this%alpha_t2(k)
            write(envelope1_tk,*) time(k)*au2fs, this%g1(k)
            write(envelope2_tk,*) time(k)*au2fs, this%g2(k)
            write(elec_field_tk,*) time(k)*au2fs, this%El(k)
            write(vec_field_tk,*) time(k)*au2fs, this%Al(k)
            write(kh_field_tk,*) time(k)*au2fs, this%alpha_t(k)            
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

    function trapezoidal(time, tp, t_mid, rise_time)
        real(dp) :: time, tp, t_mid, rise_time
        real(dp) :: trapezoidal, slope, yc 
        if (time .ge. t_mid - (tp/2 + rise_time) .and. time .le. t_mid - tp/2) then
            slope = 1._dp/rise_time
            yc = (t_mid - (tp/2 + rise_time)) * slope
            trapezoidal = slope * time - yc
        elseif (time .gt. t_mid - tp/2 .and. time .le. t_mid + tp/2) then
            trapezoidal = 1._dp
        elseif (time .gt. t_mid + tp/2 .and. time .le. t_mid+(tp/2 + rise_time)) then 
            slope = -1._dp/rise_time
            yc = (t_mid + tp/2 + rise_time) * slope
            trapezoidal = slope * time - yc
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

    function trapezoidal_vector_pulse(time, omega, phase, E0, tp, t_mid, pulse_offset, rise_time)
        real(dp) :: time, ttime, tp, t_mid, rise_time, t_start, t_end, t_u
        real(dp) :: trapezoidal_vector_pulse, pulse_offset, omega, phase, E0 
        real(dp) :: slope, yc, rise_time_new
        integer :: rise_cycles

        rise_cycles = int(rise_time*omega/(2*pi)) + 1
        rise_time_new = 2*pi*rise_cycles/omega
        pulse_offset = 0._dp ! dummy (change this later)
        t_start = 2.5_dp / au2fs ! CW envelope start
        t_end = t_start + tp + 2*rise_time_new  ! CW envelope end
        t_u = t_start + rise_time_new
        ttime = time - t_start
        
        if (time .ge. t_start .and. time .le. t_u) then
            slope = 1._dp/rise_time_new
            yc = t_start 
            trapezoidal_vector_pulse = -slope * E0 * (omega * (time - yc) * sin(omega * ttime + phase) &
                & - cos(omega * ttime + phase ) + cos(phase) ) 
            if (time-yc .le. 1e-2) then
                print*, "Vector field:"
                print*, "rise cycles: ", rise_cycles
                print*, "t-start =", t_start*au2fs, ' fs'
                print*, "t-end =", t_end*au2fs, ' fs'
                print*, "New rise-time =", rise_time_new*au2fs, ' fs'
                print*, "slope =", slope, ' E0 =', E0
                print*, "t = 0, A(t) =", trapezoidal_vector_pulse,&
                    & ', [time =', time*au2fs, ' fs]'
            endif

        elseif (time .gt. t_u .and. time .le. t_u + tp) then
            trapezoidal_vector_pulse = -E0 * omega * (sin(omega * ttime + phase) &
                & + (cos(phase) - cos(omega*rise_time_new + phase))/(omega*rise_time_new))

        elseif (time .gt. t_start + rise_time_new + tp .and. time .le. t_end) then 
            slope = -1._dp/rise_time_new
            yc = t_end 
            trapezoidal_vector_pulse = -slope * E0 * (omega * (time - yc) * sin(omega * ttime + phase) &
                & + cos(omega * ttime + phase) - cos(omega * (t_end-t_start) + phase) )
        else
            trapezoidal_vector_pulse = 0._dp
        endif
    end function trapezoidal_vector_pulse

    function trapezoidal_electric_pulse(time, omega, phase, E0, tp, t_mid, pulse_offset, rise_time)
        real(dp) :: time, ttime, tp, t_mid, rise_time, t_start, t_end, t_u
        real(dp) :: trapezoidal_electric_pulse, pulse_offset, omega, phase, E0 
        real(dp) :: slope, yc, rise_time_new
        integer :: rise_cycles

        rise_cycles = int(rise_time*omega/(2*pi)) + 1
        rise_time_new = 2*pi*rise_cycles/omega
        pulse_offset = 0._dp ! dummy (change this later)
        t_start = 2.5_dp / au2fs ! CW envelope start
        t_end = t_start + tp + 2*rise_time_new  ! CW envelope end
        t_u = t_start + rise_time_new
        ttime = time - t_start
        
        if (time .ge. t_start .and. time .le. t_u) then
            slope = 1._dp/rise_time_new
            yc = t_start 
            trapezoidal_electric_pulse = -slope * E0 * (omega**2 * (time - yc) * cos(omega * ttime + phase) &
                & + 2*omega*sin(omega * ttime + phase ) - 2*omega*sin(phase)) 
            if (time-yc .le. 1e-2) then
                print*, "Electric field:"
                print*, "rise cycles: ", rise_cycles
                print*, "t-start =", t_start*au2fs, ' fs'
                print*, "t-end =", t_end*au2fs, ' fs'
                print*, "New rise-time =", rise_time_new*au2fs, ' fs'
                print*, "slope =", slope, ' E0 =', E0
                print*, "t = 0, A(t) =", trapezoidal_electric_pulse,&
                    & ', [time =', time*au2fs, ' fs]'
            endif

        elseif (time .gt. t_u .and. time .le. t_u + tp) then
            trapezoidal_electric_pulse = -E0 * omega**2 * (cos(omega * ttime + phase) &
                & + 2 * (sin(omega*rise_time_new + phase)- sin(phase))/(omega*rise_time_new))

        elseif (time .gt. t_start + rise_time_new + tp .and. time .le. t_end) then 
            slope = -1._dp/rise_time_new
            yc = t_end 
            trapezoidal_electric_pulse = -slope * E0 * (omega**2 * (time - yc) * cos(omega * ttime + phase) &
                & + 2*omega*sin(omega * ttime + phase) - 2*omega*sin(omega * (t_end-t_start) + phase) )
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
 

