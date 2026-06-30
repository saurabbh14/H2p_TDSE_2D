module split_operator_2d_mod
    use global_vars, only: dp
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: split_operator_2d_type
    ! Type for split-operator propagation and FFTW management
    type :: split_operator_2d_type
        character(20) :: gauge ! Guage Type: "length" or "velocity"
        ! Kinetic propagator for half and full time steps
        complex(dp), allocatable :: kprop_full(:,:)
        ! Potential propagator 
        complex(dp), allocatable :: vprop(:,:), vcol_prop(:)
        ! gauge transform factor
        complex(dp), allocatable :: gauge_transform(:,:)
        ! FFTW plan and memory pointers
        type(C_PTR) :: planF, planB, p_in, p_out
        ! FFTW input/output arrays
        complex(C_DOUBLE_COMPLEX), pointer:: psi_in(:,:), psi_out(:,:)
    contains
        procedure :: fft_initialize      ! Initialize FFTW plans and memory
        procedure :: split_operator_initialize ! Initialize memory and functions based on gauge choice
        procedure :: kprop_gen_len       ! Generate length-guage kinetic propagators 
        procedure :: vprop_gen_len       ! Generate length-guage potential propagators
        procedure :: kprop_gen_vel       ! Generate velocity-guage kinetic propagators 
        procedure :: vprop_gen_vel       ! Generate velocity-guage potential propagators
        procedure :: split_operator_step      ! Apply split-operator step
        procedure :: finalize             ! Clean up FFTW resources
    end type split_operator_2d_type

contains
    !> Initialize FFTW plans and memory for split-operator propagation
    subroutine fft_initialize(this)
        use global_vars, only: NR, Nx, prop_par_FFTW
        use FFTW3
        class(split_operator_2d_type), intent(inout) :: this

        print*
        print*, "FFTW intialization ..."
        print*

        ! Creating aligned memory for FFTW
        this%p_in = fftw_alloc_complex(int(NR * Nx, C_SiZE_T)) 
        call c_f_pointer(this%p_in,this%psi_in,[NR,Nx])
        this%p_out = fftw_alloc_complex(int(NR * Nx, C_SiZE_T)) 
        call c_f_pointer(this%p_out,this%psi_out,[NR,Nx])

        call fftw_initialize_threads
        print*, "FFTW plan creation ..."
        call fftw_create_c2c_2d_plans(this%psi_in, this%psi_out, NR, Nx, & 
            & this%planF, this%planB, prop_par_FFTW)
        print*, "Done setting up FFTW."

    end subroutine fft_initialize

    subroutine split_operator_initialize(this) 
        use global_vars, only: NR, Nx
        class(split_operator_2d_type), intent(inout) :: this

        allocate(this%kprop_full(NR,Nx))
        allocate(this%vprop(NR,Nx), this%vcol_prop(NR))
        allocate(this%gauge_transform(NR,Nx))

        select case(this%gauge)
        case("length")
            call this%kprop_gen_len()
            
        case("velocity")
            call this%vprop_gen_vel() 
        end select  
    end subroutine

    !> Generate kinetic propagators for half and full time steps
    subroutine kprop_gen_len(this)
        use global_vars, only: Nx, dt, m_red, PR, px
        use data_au, only: im
        class(split_operator_2d_type), intent(inout) :: this
        integer:: j
        
        do j = 1, Nx
            this%kprop_full(:,j) = exp(-im * dt * (pR(:) * pR(:) / (2._dp*m_red) + 0.5_dp * px(j) * px(j))) 
        end do
         
    end subroutine kprop_gen_len

    subroutine kprop_gen_vel(this, A)
        use global_vars, only: Nx, dt, m_red, PR, px, lam, kap, x, R
        use data_au, only: im
        class(split_operator_2d_type), intent(inout) :: this
        integer :: j
        real(dp) :: A
        
        do j = 1, Nx
            this%kprop_full(:,j) = exp(-im * dt * ((pR(:) + lam * A)**2  / (2._dp*m_red) &
                & + 0.5_dp * (px(j) + kap * A)**2)) 
            this%gauge_transform(:,j) = exp(im * A * (x(j) + R(:)))
        end do
         
    end subroutine kprop_gen_vel

    !> Generate potential propagators 
    subroutine vprop_gen_len(this, E, A)
        use global_vars, only: Nx, dt, pot, dp, kap, lam, R, x
        use data_au, only: im
        class(split_operator_2d_type), intent(inout) :: this
        integer :: j
        real(dp) :: E, A

        do j = 1, Nx          
            this%vprop(:,j) = exp(-im * 0.5_dp * dt * (pot(:,j) + (kap*x(j)*E + lam*R(:)*E)))
            this%gauge_transform(:,j) = exp(-im * A * (x(j) + R(:)))          
        end do
        this%vcol_prop = exp(-im * 0.5_dp * dt / R)

    end subroutine vprop_gen_len

    subroutine vprop_gen_vel(this)
        use global_vars, only: dt, pot, dp, R
        use data_au, only: im
        class(split_operator_2d_type), intent(inout) :: this
                
        this%vprop = exp(-im * 0.5_dp * dt * pot)  
        this%vcol_prop = exp(-im * 0.5_dp * dt / R)
        
    end subroutine vprop_gen_vel

    !> Apply split-operator step to wavefunction psi_ges
    subroutine split_operator_step(this, psi, region)
        use global_vars, only: NR, Nx
        use FFTW3
        class(split_operator_2d_type), intent(inout) :: this
        complex(dp), intent(inout):: psi(NR, Nx)
        character(20) :: region

        this%psi_in = (0._dp, 0._dp)
        this%psi_out = (0._dp, 0._dp)
        select case(adjustl(trim(region)))
        case("inner-xR")
            this%psi_in = psi * this%vprop  ! Hilfsgroesse
            call fftw_execute_dft(this%planF, this%psi_in, this%psi_out)
            this%psi_in = this%psi_out * this%kprop_full
            call fftw_execute_dft(this%planB, this%psi_in, this%psi_out)
            this%psi_in = this%psi_out * this%vprop 
            psi = this%psi_in / dble(NR*Nx)
        
        case("outer-x")
            print*, "to be implemented"

        case("outer-R")
            print*, "to be implemented"
        
        case("outer-xR")
            if (this%gauge == "velocity") then
                psi = psi * this%kprop_full
            else
                psi = psi * this%kprop_full * this%gauge_transform
            end if

        end select

    end subroutine split_operator_step

    !> Clean up FFTW plans and memory
    subroutine finalize(this)
        use FFTW3
        class(split_operator_2d_type), intent(inout) :: this

        call fftw_destroy_plan(this%planF)
        call fftw_destroy_plan(this%planB)
        call fftw_free(this%p_in)
        call fftw_free(this%p_out)

    end subroutine finalize

end module split_operator_2d_mod
