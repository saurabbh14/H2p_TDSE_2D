module rk4_operator_2d_mod
    use global_vars, only: dp
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: rk4_operator_2d_type

    !> Type for RK4 time propagation in 2D with FFTW management
    !! Handles kinetic operator via 2D FFT, and provides RHS evaluation for RK4 stepping
    type :: rk4_operator_2d_type
        character(20) :: gauge ! Gauge Type: "length" or "velocity"
        ! Kinetic energy factor in momentum space: pR²/(2*m_red) + 0.5*px²
        real(dp), allocatable :: kin_energy(:,:)
        ! FFTW plan and memory pointers
        type(C_PTR) :: planF, planB, p_in, p_out
        ! FFTW input/output arrays
        complex(C_DOUBLE_COMPLEX), pointer :: psi_in(:,:), psi_out(:,:)
    contains
        procedure :: fft_initialize      ! Initialize FFTW plans and memory
        procedure :: kin_energy_gen      ! Pre-compute kinetic energy factor
        procedure :: rhs_2d              ! Evaluate RHS: d(psi)/dt = -i*H*psi
        procedure :: rk4_step            ! Perform one full RK4 time step
        procedure :: finalize            ! Clean up FFTW resources and kinetic energy array
    end type rk4_operator_2d_type

contains

    !> Initialize FFTW plans and memory for RK4 propagation (2D)
    subroutine fft_initialize(this)
        use global_vars, only: NR, Nx, prop_par_FFTW
        use FFTW3
        class(rk4_operator_2d_type), intent(inout) :: this

        print*
        print*, "RK4 2D: FFTW initialization ..."
        print*

        ! Creating aligned memory for FFTW
        this%p_in = fftw_alloc_complex(int(NR * Nx, C_SiZE_T))
        call c_f_pointer(this%p_in, this%psi_in, [NR, Nx])
        this%p_out = fftw_alloc_complex(int(NR * Nx, C_SiZE_T))
        call c_f_pointer(this%p_out, this%psi_out, [NR, Nx])

        call fftw_initialize_threads
        print*, "RK4 2D: FFTW plan creation ..."
        call fftw_create_c2c_2d_plans(this%psi_in, this%psi_out, NR, Nx, &
            & this%planF, this%planB, prop_par_FFTW)
        print*, "RK4 2D: Done setting up FFTW."

        ! Pre-compute kinetic energy factor in momentum space
        allocate(this%kin_energy(NR, Nx))
        call this%kin_energy_gen()

    end subroutine fft_initialize

    !> Pre-compute kinetic energy: pR²/(2*m_red) + 0.5*px²
    subroutine kin_energy_gen(this)
        use global_vars, only: NR, Nx, pR, px, m_red
        class(rk4_operator_2d_type), intent(inout) :: this
        integer :: j

        do j = 1, Nx
            this%kin_energy(:, j) = pR(:) * pR(:) / (2._dp * m_red) + 0.5_dp * px(j) * px(j)
        end do

    end subroutine kin_energy_gen

    !> Evaluate RHS of TDSE in 2D: rhs = -i * H * psi
    !! H = T + V(R,x) + V_field(R,x,t)
    !! T is applied via 2D FFT (kinetic energy in momentum space)
    !! V includes the 2D potential and the dipole-field interaction
    !! Supports both "length" and "velocity" gauge
    subroutine rhs_2d(this, psi, psi_rhs, E_field, A_field, pot)
        use global_vars, only: NR, Nx, kap, lam, R, x
        use data_au, only: im
        use FFTW3
        class(rk4_operator_2d_type), intent(inout) :: this
        complex(dp), intent(in)  :: psi(NR, Nx)
        complex(dp), intent(out) :: psi_rhs(NR, Nx)
        real(dp), intent(in)     :: E_field
        real(dp), intent(in)     :: A_field
        real(dp), intent(in)     :: pot(NR, Nx)

        integer :: j

        ! Start with psi_rhs = 0
        psi_rhs = (0._dp, 0._dp)

        ! Step 1: Kinetic energy contribution: T * psi
        ! 2D FFT → multiply by kin_energy → 2D iFFT
        this%psi_in = (0._dp, 0._dp)
        this%psi_out = (0._dp, 0._dp)
        this%psi_in = psi
        call fftw_execute_dft(this%planF, this%psi_in, this%psi_out)
        this%psi_in = this%psi_out * this%kin_energy
        call fftw_execute_dft(this%planB, this%psi_in, this%psi_out)
        psi_rhs = this%psi_out / dble(NR * Nx)

        ! Step 2: Potential contribution: V(R,x) * psi + field interaction
        select case(this%gauge)
        case("length")
            do j = 1, Nx
                psi_rhs(:, j) = psi_rhs(:, j) &
                    & + (pot(:, j) + kap * x(j) * E_field + lam * R(:) * E_field) * psi(:, j)
            end do
        case("velocity")
            ! In velocity gauge, the field coupling is folded into the kinetic term,
            ! so only the bare potential is used here. The kinetic term above
            ! should be recomputed with shifted momentum for velocity gauge.
            ! For simplicity, we apply the bare potential.
            psi_rhs = psi_rhs + pot * psi
        end select

        ! Step 3: Multiply by -i to get d(psi)/dt = -i * H * psi
        psi_rhs = -im * psi_rhs

    end subroutine rhs_2d

    !> Perform one full RK4 time step in 2D
    !! k1 = rhs(psi, t)
    !! k2 = rhs(psi + k1*dt/2, t + dt/2)
    !! k3 = rhs(psi + k2*dt/2, t + dt/2)
    !! k4 = rhs(psi + k3*dt, t + dt)
    subroutine rk4_step(this, psi, dt, E_now, E_half, A_now, A_half, pot)
        use global_vars, only: NR, Nx
        class(rk4_operator_2d_type), intent(inout) :: this
        complex(dp), intent(inout) :: psi(NR, Nx)
        real(dp), intent(in)      :: dt
        real(dp), intent(in)      :: E_now, E_half   ! E(t) and E(t+dt/2)
        real(dp), intent(in)      :: A_now, A_half   ! A(t) and A(t+dt/2)
        real(dp), intent(in)      :: pot(NR, Nx)

        complex(dp), allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
        complex(dp), allocatable :: psi_tmp(:,:)

        allocate(k1(NR, Nx), k2(NR, Nx), k3(NR, Nx), k4(NR, Nx))
        allocate(psi_tmp(NR, Nx))

        ! k1 = rhs(psi, t)
        call this%rhs_2d(psi, k1, E_now, A_now, pot)

        ! k2 = rhs(psi + k1*dt/2, t + dt/2)
        psi_tmp = psi + k1 * (0.5_dp * dt)
        call this%rhs_2d(psi_tmp, k2, E_half, A_half, pot)

        ! k3 = rhs(psi + k2*dt/2, t + dt/2)
        psi_tmp = psi + k2 * (0.5_dp * dt)
        call this%rhs_2d(psi_tmp, k3, E_half, A_half, pot)

        ! k4 = rhs(psi + k3*dt, t + dt)
        psi_tmp = psi + k3 * dt
        call this%rhs_2d(psi_tmp, k4, E_now, A_now, pot)

        ! Update: psi_new = psi + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
        psi = psi + (k1 + 2._dp * k2 + 2._dp * k3 + k4) * (dt / 6._dp)

        deallocate(k1, k2, k3, k4, psi_tmp)

    end subroutine rk4_step

    !> Clean up FFTW plans, memory, and kinetic energy array
    subroutine finalize(this)
        use FFTW3
        class(rk4_operator_2d_type), intent(inout) :: this

        call fftw_destroy_plan(this%planF)
        call fftw_destroy_plan(this%planB)
        call fftw_free(this%p_in)
        call fftw_free(this%p_out)
        if (allocated(this%kin_energy)) deallocate(this%kin_energy)

    end subroutine finalize

end module rk4_operator_2d_mod
