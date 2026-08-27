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
        procedure :: itp_initialize      ! Initialize memory and imaginary-time propagators (2D ITP)
        procedure :: kprop_gen_len       ! Generate length-guage kinetic propagators 
        procedure :: vprop_gen_len       ! Generate length-guage potential propagators
        procedure :: vprop_gen_kh        ! Generate KH-gauge (time-dep) potential propagators
        procedure :: kprop_gen_vel       ! Generate velocity-guage kinetic propagators 
        procedure :: vprop_gen_vel       ! Generate velocity-guage potential propagators
        procedure :: split_operator_step      ! Apply split-operator step
        procedure :: finalize             ! Clean up FFTW resources
    end type split_operator_2d_type

contains
    !> Initialize FFTW plans and memory for split-operator propagation.
    !! `parallel` optionally overrides the global `prop_par_FFTW` flag (used by
    !! the 2D ITP module, which follows the `itp_fftw` input setting instead).
    subroutine fft_initialize(this, parallel)
        use global_vars, only: NR, Nx, prop_par_FFTW
        use FFTW3
        class(split_operator_2d_type), intent(inout) :: this
        character(*), intent(in), optional :: parallel
        character(len=10) :: par_flag

        if (present(parallel)) then
            par_flag = parallel
        else
            par_flag = prop_par_FFTW
        end if

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
            & this%planF, this%planB, par_flag)
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

    !> Initialize memory and the *imaginary-time* propagators used by the 2D ITP.
    !! Instead of the real-time factors exp(-i dt H) this builds
    !!   kprop_full = exp(-dt_itp * (pR^2/(2 m_red) + px^2/(2 m_eff)))
    !!   vprop      = exp(-0.5 * dt_itp * pot)
    !! so that `split_operator_step(psi, "inner-xR")` performs one Strang step
    !! of imaginary-time propagation without any further changes.
    !! The field-free (Lab-frame) Hamiltonian is used: the ITP looks for the
    !! bound eigenstates of the unperturbed system, hence no gauge/laser terms.
    subroutine itp_initialize(this, dt_itp)
        use global_vars, only: NR, Nx, m_red, m_eff, pR, px, pot
        class(split_operator_2d_type), intent(inout) :: this
        real(dp), intent(in) :: dt_itp
        integer :: j

        if (.not. allocated(this%kprop_full)) allocate(this%kprop_full(NR,Nx))
        if (.not. allocated(this%vprop)) allocate(this%vprop(NR,Nx))
        if (.not. allocated(this%vcol_prop)) allocate(this%vcol_prop(NR))
        if (.not. allocated(this%gauge_transform)) allocate(this%gauge_transform(NR,Nx))

        do j = 1, Nx
            this%kprop_full(:,j) = cmplx(exp(-dt_itp * (pR(:)**2 / (2._dp*m_red) &
                & + px(j)**2 / (2._dp*m_eff))), 0._dp, dp)
        end do
        this%vprop = cmplx(exp(-0.5_dp * dt_itp * pot), 0._dp, dp)
        ! Not used in imaginary time (the 1/R repulsion is already part of pot),
        ! kept neutral so that any accidental use is a no-op.
        this%vcol_prop = (1._dp, 0._dp)
        this%gauge_transform = (1._dp, 0._dp)

    end subroutine itp_initialize

    !> Generate kinetic propagators for half and full time steps
    !! Kinetic energy: pR²/(2*m_red) + px²/(2*m_eff), with m_eff the effective
    !! electron mass (m1+m2)/(m1+m2+1) of the (R, x) coordinate system.
    subroutine kprop_gen_len(this)
        use global_vars, only: Nx, dt, m_red, m_eff, PR, px
        use data_au, only: im
        class(split_operator_2d_type), intent(inout) :: this
        integer:: j
        
        do j = 1, Nx
            this%kprop_full(:,j) = exp(-im * dt * (pR(:) * pR(:) / (2._dp*m_red) &
                & + px(j) * px(j) / (2._dp*m_eff))) 
        end do
         
    end subroutine kprop_gen_len

    subroutine kprop_gen_vel(this, A)
        use global_vars, only: Nx, dt, m_red, m_eff, PR, px, lam, kap, x, R
        use data_au, only: im
        class(split_operator_2d_type), intent(inout) :: this
        integer :: j
        real(dp) :: A
        
        do j = 1, Nx
            this%kprop_full(:,j) = exp(-im * dt * ((pR(:) + lam * A)**2  / (2._dp*m_red) &
                & + (px(j) + kap * A)**2 / (2._dp*m_eff))) 
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

    !> Generate KH-gauge potential propagator (no E-field term — laser coupling in potential)
    subroutine vprop_gen_kh(this, pot_kh)
        use global_vars, only: dt, dp, R
        use data_au, only: im
        class(split_operator_2d_type), intent(inout) :: this
        real(dp), intent(in) :: pot_kh(:,:)

        this%vprop = exp(-im * 0.5_dp * dt * pot_kh)
        this%vcol_prop = exp(-im * 0.5_dp * dt / R)  ! dummy, not used in KH

    end subroutine vprop_gen_kh

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
