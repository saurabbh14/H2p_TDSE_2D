module split_operator_mod
    use global_vars, only: dp
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: split_operator_type
    ! Type for split-operator propagation and FFTW management
    type :: split_operator_type
        ! Kinetic propagator for half and full time steps
        complex(dp), allocatable :: kprop_half(:), kprop_full(:)
        ! Potential propagator for each electronic state
        complex(dp), allocatable :: vprop(:,:)
        ! FFTW plan and memory pointers
        type(C_PTR) :: planF, planB, p_in, p_out
        ! FFTW input/output arrays
        complex(C_DOUBLE_COMPLEX), pointer:: psi_in(:), psi_out(:)
    contains
        procedure :: fft_initialize      ! Initialize FFTW plans and memory
        procedure :: kprop_gen           ! Generate kinetic propagators
        procedure :: vprop_gen           ! Generate potential propagators
        procedure :: split_operator_step      ! Apply split-operator step
    end type split_operator_type

contains
    !> Initialize FFTW plans and memory for split-operator propagation
    subroutine fft_initialize(this)
        use global_vars, only: NR, prop_par_FFTW
        use FFTW3
        class(split_operator_type), intent(inout) :: this

        print*
        print*, "FFTW intialization ..."
        print*

        ! Creating aligned memory for FFTW
        this%p_in = fftw_alloc_complex(int(NR, C_SiZE_T)) 
        call c_f_pointer(this%p_in,this%psi_in,[NR])
        this%p_out = fftw_alloc_complex(int(NR, C_SiZE_T)) 
        call c_f_pointer(this%p_out,this%psi_out,[NR])

        call fftw_initialize_threads
        print*, "FFTW plan creation ..."
        call fftw_create_c2c_plans(this%psi_in, this%psi_out, NR, & 
            & this%planF, this%planB, prop_par_FFTW)
        print*, "Done setting up FFTW."

    end subroutine fft_initialize

    !> Generate kinetic propagators for half and full time steps
    subroutine kprop_gen(this)
        use global_vars, only: NR, dt, m_red, PR
        use data_au, only: im
        class(split_operator_type), intent(inout) :: this

        allocate(this%kprop_half(NR), this%kprop_full(NR))
        
        this%kprop_half = exp(-im *dt * pR*pR /(4._dp*m_red))  ! pR**2 /2 * red_mass UND Half time step
        this%kprop_full = exp(-im *dt * pR*pR /(2._dp*m_red))  
         
    end subroutine kprop_gen

    !> Generate potential propagators for all electronic states
    subroutine vprop_gen(this)
        use global_vars, only: NR, Nstates, dt, adb, dp
        use data_au, only: im
        class(split_operator_type), intent(inout) :: this
        integer:: j
        
        allocate(this%vprop(NR,Nstates))
        do j = 1, Nstates           
            this%vprop(:,j) = exp(-im * 0.5_dp * dt * adb(:,j)) !+0.8d0*R(i)*E(K)))!+H_ac(i,j))) !         
        end do

    end subroutine vprop_gen

    !> Apply split-operator step to wavefunction psi_ges
    subroutine split_operator_step(this, psi_ges)
        use global_vars, only: NR, Nstates
        use FFTW3
        class(split_operator_type), intent(inout) :: this
        complex(dp), intent(inout):: psi_ges(NR, Nstates)
        integer:: j

        do j = 1, Nstates
            this%psi_in = (0._dp, 0._dp)
            this%psi_out = (0._dp, 0._dp)
            this%psi_in(:) = psi_ges(:,J)  ! Hilfsgroesse
            call fftw_execute_dft(this%planF, this%psi_in, this%psi_out)
            this%psi_in = this%psi_out * this%kprop_half
            call fftw_execute_dft(this%planB, this%psi_in, this%psi_out)
            this%psi_in = this%psi_out / dble(NR)
            psi_ges(:,J) = this%psi_in(:)      
        end do
    end subroutine split_operator_step

end module