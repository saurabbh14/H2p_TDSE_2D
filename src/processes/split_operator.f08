module split_operator_mod
    use global_vars, only: dp
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: split_operator_type
    ! Type for split-operator propagation and FFTW management
    type :: split_operator_type
        character(20) :: gauge ! Gauge Type: "length" or "velocity"
        ! Kinetic propagator for half and full time steps
        complex(dp), allocatable :: kprop_half(:), kprop_full(:)
        ! Potential propagator for each electronic state
        complex(dp), allocatable :: vprop(:,:)
        ! Gauge transform factor (used in velocity gauge for absorber region)
        complex(dp), allocatable :: gauge_transform(:)
        ! FFTW plan and memory pointers
        type(C_PTR) :: planF, planB, p_in, p_out
        ! FFTW input/output arrays
        complex(C_DOUBLE_COMPLEX), pointer:: psi_in(:), psi_out(:)
    contains
        procedure :: fft_initialize      ! Initialize FFTW plans and memory
        procedure :: kprop_gen_len       ! Generate length-gauge kinetic propagators
        procedure :: kprop_gen_vel       ! Generate velocity-gauge kinetic propagators
        procedure :: vprop_gen           ! Generate potential propagators
        procedure :: gauge_transform_gen ! Generate gauge transform factor
        procedure :: split_operator_step      ! Apply split-operator step
        procedure :: finalize             ! Clean up FFTW resources
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

    !> Generate length-gauge kinetic propagators for half and full time steps
    subroutine kprop_gen_len(this)
        use global_vars, only: NR, dt, m_red, PR
        use data_au, only: im
        class(split_operator_type), intent(inout) :: this

        allocate(this%kprop_half(NR), this%kprop_full(NR))
        
        this%kprop_half = exp(-im *dt * pR*pR /(4._dp*m_red))  ! pR**2 /2 * red_mass UND Half time step
        this%kprop_full = exp(-im *dt * pR*pR /(2._dp*m_red))  
         
    end subroutine kprop_gen_len

    !> Generate velocity-gauge kinetic propagators for half and full time steps
    !! Uses shifted momentum: p -> p + lam*A
    subroutine kprop_gen_vel(this, A)
        use global_vars, only: NR, dt, m_red, PR, lam
        use data_au, only: im
        class(split_operator_type), intent(inout) :: this
        real(dp), intent(in) :: A

        if (.not. allocated(this%kprop_half)) allocate(this%kprop_half(NR))
        if (.not. allocated(this%kprop_full)) allocate(this%kprop_full(NR))
        
        this%kprop_half = exp(-im * dt * (pR(:) + lam * A)**2 / (4._dp * m_red))
        this%kprop_full = exp(-im * dt * (pR(:) + lam * A)**2 / (2._dp * m_red))
         
    end subroutine kprop_gen_vel

    !> Generate gauge transform factor for velocity gauge
    !! U(R) = exp(-i * lam * A * R) used to convert wavefunction between gauges
    subroutine gauge_transform_gen(this, A)
        use global_vars, only: NR, lam, R
        use data_au, only: im
        class(split_operator_type), intent(inout) :: this
        real(dp), intent(in) :: A

        if (.not. allocated(this%gauge_transform)) allocate(this%gauge_transform(NR))
        
        this%gauge_transform(:) = exp(-im * lam * A * R(:))
         
    end subroutine gauge_transform_gen

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
    subroutine split_operator_step(this, psi_ges, apply_gauge_transform)
        use global_vars, only: NR, Nstates
        use FFTW3
        class(split_operator_type), intent(inout) :: this
        complex(dp), intent(inout):: psi_ges(NR, Nstates)
        logical, intent(in), optional :: apply_gauge_transform
        integer:: j
        logical :: apply_gt

        apply_gt = .false.
        if (present(apply_gauge_transform)) apply_gt = apply_gauge_transform

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

        ! Apply gauge transform after kinetic step if requested (for velocity gauge absorber output)
        if (apply_gt .and. allocated(this%gauge_transform)) then
            do j = 1, Nstates
                psi_ges(:,j) = psi_ges(:,j) * this%gauge_transform(:)
            end do
        end if

    end subroutine split_operator_step

    !> Clean up FFTW plans and memory
    subroutine finalize(this)
        use FFTW3
        class(split_operator_type), intent(inout) :: this

        call fftw_destroy_plan(this%planF)
        call fftw_destroy_plan(this%planB)
        call fftw_free(this%p_in)
        call fftw_free(this%p_out)

    end subroutine finalize

end module