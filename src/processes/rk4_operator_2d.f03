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
        procedure :: rk4_step_kh         ! Perform one full RK4 time step in KH gauge
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
        use global_vars, only: NR, Nx, kap, lam, R, x, pR, px, m_red
        use data_au, only: im
        use FFTW3
        class(rk4_operator_2d_type), intent(inout) :: this
        complex(dp), intent(in)  :: psi(NR, Nx)
        complex(dp), intent(out) :: psi_rhs(NR, Nx)
        real(dp), intent(in)     :: E_field
        real(dp), intent(in)     :: A_field
        real(dp), intent(in)     :: pot(NR, Nx)

        integer :: j
        real(dp), allocatable :: kin_shifted(:,:)

        ! Start with psi_rhs = 0
        psi_rhs = (0._dp, 0._dp)

        ! Step 1: Kinetic energy contribution: T * psi
        ! 2D FFT → multiply by kin_energy → 2D iFFT
        this%psi_in = (0._dp, 0._dp)
        this%psi_out = (0._dp, 0._dp)
        this%psi_in = psi
        call fftw_execute_dft(this%planF, this%psi_in, this%psi_out)

        select case(this%gauge)
        case("length")
            this%psi_in = this%psi_out * this%kin_energy
        case("velocity")
            ! Compute shifted kinetic energy on the fly:
            ! (pR + lam*A)²/(2*m_red) + 0.5*(px + kap*A)²
            allocate(kin_shifted(NR, Nx))
            do j = 1, Nx
                kin_shifted(:, j) = (pR(:) + lam * A_field)**2 / (2._dp * m_red) &
                    & + 0.5_dp * (px(j) + kap * A_field)**2
            end do
            this%psi_in = this%psi_out * kin_shifted
            deallocate(kin_shifted)
        end select

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
            ! In velocity gauge, field coupling is in the shifted kinetic term.
            ! Apply only the bare potential (no E-field coupling).
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
    subroutine rk4_step(this, psi, dt, E_now, E_half, E_next, A_now, A_half, A_next, pot)
        use global_vars, only: NR, Nx
        class(rk4_operator_2d_type), intent(inout) :: this
        complex(dp), intent(inout) :: psi(NR, Nx)
        real(dp), intent(in)      :: dt
        real(dp), intent(in)      :: E_now, E_half, E_next  ! E(t), E(t+dt/2), E(t+dt)
        real(dp), intent(in)      :: A_now, A_half, A_next  ! A(t), A(t+dt/2), A(t+dt)
        real(dp), intent(in)      :: pot(NR, Nx)

        complex(dp), allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
        complex(dp), allocatable :: psi_tmp(:,:)

        allocate(k1(NR, Nx), k2(NR, Nx), k3(NR, Nx), k4(NR, Nx))
        allocate(psi_tmp(NR, Nx))

        ! k1 = rhs(psi, t)              using E(t), A(t)
        call this%rhs_2d(psi, k1, E_now, A_now, pot)

        ! k2 = rhs(psi + k1*dt/2, t+dt/2) using E(t+dt/2), A(t+dt/2)
        psi_tmp = psi + k1 * (0.5_dp * dt)
        call this%rhs_2d(psi_tmp, k2, E_half, A_half, pot)

        ! k3 = rhs(psi + k2*dt/2, t+dt/2) using E(t+dt/2), A(t+dt/2)
        psi_tmp = psi + k2 * (0.5_dp * dt)
        call this%rhs_2d(psi_tmp, k3, E_half, A_half, pot)

        ! k4 = rhs(psi + k3*dt,   t+dt)   using E(t+dt), A(t+dt)
        psi_tmp = psi + k3 * dt
        call this%rhs_2d(psi_tmp, k4, E_next, A_next, pot)

        ! Update: psi_new = psi + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
        psi = psi + (k1 + 2._dp * k2 + 2._dp * k3 + k4) * (dt / 6._dp)

        deallocate(k1, k2, k3, k4, psi_tmp)

    end subroutine rk4_step

    !> Perform one full RK4 time step in KH gauge
    !! Uses time-dependent KH potential (laser coupling already in potential, no E-field term)
    !! pot_now  = V_KH(R, x, t)
    !! pot_half = V_KH(R, x, t + dt/2)
    !! pot_next = V_KH(R, x, t + dt)
    subroutine rk4_step_kh(this, psi, dt, pot_now, pot_half, pot_next)
        use global_vars, only: NR, Nx
        class(rk4_operator_2d_type), intent(inout) :: this
        complex(dp), intent(inout) :: psi(NR, Nx)
        real(dp), intent(in)      :: dt
        real(dp), intent(in)      :: pot_now(NR, Nx), pot_half(NR, Nx), pot_next(NR, Nx)

        complex(dp), allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
        complex(dp), allocatable :: psi_tmp(:,:)

        ! E_field = 0 and A_field = 0 since laser coupling is in the potential
        real(dp), parameter :: E_zero = 0._dp, A_zero = 0._dp

        allocate(k1(NR, Nx), k2(NR, Nx), k3(NR, Nx), k4(NR, Nx))
        allocate(psi_tmp(NR, Nx))

        ! k1 = rhs(psi, t) with pot_now
        call this%rhs_2d(psi, k1, E_zero, A_zero, pot_now)

        ! k2 = rhs(psi + k1*dt/2, t + dt/2) with pot_half
        psi_tmp = psi + k1 * (0.5_dp * dt)
        call this%rhs_2d(psi_tmp, k2, E_zero, A_zero, pot_half)

        ! k3 = rhs(psi + k2*dt/2, t + dt/2) with pot_half
        psi_tmp = psi + k2 * (0.5_dp * dt)
        call this%rhs_2d(psi_tmp, k3, E_zero, A_zero, pot_half)

        ! k4 = rhs(psi + k3*dt, t + dt) with pot_next
        psi_tmp = psi + k3 * dt
        call this%rhs_2d(psi_tmp, k4, E_zero, A_zero, pot_next)

        ! Update: psi_new = psi + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
        psi = psi + (k1 + 2._dp * k2 + 2._dp * k3 + k4) * (dt / 6._dp)

        deallocate(k1, k2, k3, k4, psi_tmp)

    end subroutine rk4_step_kh

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
