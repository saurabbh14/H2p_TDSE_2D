module rk4_operator_mod
    use global_vars, only: dp
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: rk4_operator_type

    !> Type for RK4 time propagation with FFTW management
    !! Handles kinetic operator via FFT, and provides RHS evaluation for RK4 stepping
    type :: rk4_operator_type
        ! Kinetic energy factor in momentum space: p²/(2*m_red)
        real(dp), allocatable :: kin_energy(:)
        ! FFTW plan and memory pointers
        type(C_PTR) :: planF, planB, p_in, p_out
        ! FFTW input/output arrays
        complex(C_DOUBLE_COMPLEX), pointer :: psi_in(:), psi_out(:)
    contains
        procedure :: fft_initialize      ! Initialize FFTW plans and memory
        procedure :: kin_energy_gen      ! Pre-compute kinetic energy factor
        procedure :: rhs_1d              ! Evaluate RHS: d(psi)/dt = -i*H*psi
        procedure :: rk4_step            ! Perform one full RK4 time step
        procedure :: finalize            ! Clean up FFTW resources and kinetic energy array
    end type rk4_operator_type

contains

    !> Initialize FFTW plans and memory for RK4 propagation
    subroutine fft_initialize(this)
        use global_vars, only: NR, prop_par_FFTW
        use FFTW3
        class(rk4_operator_type), intent(inout) :: this

        print*
        print*, "RK4: FFTW initialization ..."
        print*

        ! Creating aligned memory for FFTW
        this%p_in = fftw_alloc_complex(int(NR, C_SiZE_T))
        call c_f_pointer(this%p_in, this%psi_in, [NR])
        this%p_out = fftw_alloc_complex(int(NR, C_SiZE_T))
        call c_f_pointer(this%p_out, this%psi_out, [NR])

        call fftw_initialize_threads
        print*, "RK4: FFTW plan creation ..."
        call fftw_create_c2c_plans(this%psi_in, this%psi_out, NR, &
            & this%planF, this%planB, prop_par_FFTW)
        print*, "RK4: Done setting up FFTW."

        ! Pre-compute kinetic energy factor in momentum space
        allocate(this%kin_energy(NR))
        call this%kin_energy_gen()

    end subroutine fft_initialize

    !> Pre-compute kinetic energy: p²/(2*m_red)
    subroutine kin_energy_gen(this)
        use global_vars, only: pR, m_red
        class(rk4_operator_type), intent(inout) :: this

        this%kin_energy(:) = pR(:) * pR(:) / (2._dp * m_red)

    end subroutine kin_energy_gen

    !> Evaluate RHS of TDSE: rhs = -i * H * psi
    !! H = T + V_diag + V_coupling(t)
    !! T is applied via FFT (kinetic energy in momentum space)
    !! V_diag is the adiabatic potential for each electronic state
    !! V_coupling is the off-diagonal dipole coupling matrix
    subroutine rhs_1d(this, psi, psi_rhs, E_field, mu_all, adb)
        use global_vars, only: NR, Nstates
        use data_au, only: im
        use blas_interfaces_module, only: zgemv
        use FFTW3
        class(rk4_operator_type), intent(inout) :: this
        complex(dp), intent(in)  :: psi(NR, Nstates)
        complex(dp), intent(out) :: psi_rhs(NR, Nstates)
        real(dp), intent(in)     :: E_field
        real(dp), intent(in)     :: mu_all(Nstates, Nstates, NR)
        real(dp), intent(in)     :: adb(NR, Nstates)

        integer :: i, j
        complex(dp) :: tout(Nstates, Nstates)
        complex(dp) :: psi_Nstates(Nstates), psi_Nstates1(Nstates)

        ! Start with psi_rhs = 0
        psi_rhs = (0._dp, 0._dp)

        ! Step 1: Kinetic energy contribution: T * psi
        ! For each electronic state: FFT → multiply by kin_energy → iFFT
        do j = 1, Nstates
            this%psi_in = (0._dp, 0._dp)
            this%psi_out = (0._dp, 0._dp)
            this%psi_in(:) = psi(:, j)
            call fftw_execute_dft(this%planF, this%psi_in, this%psi_out)
            ! Multiply by kinetic energy in momentum space
            this%psi_in = this%psi_out * this%kin_energy(:)
            call fftw_execute_dft(this%planB, this%psi_in, this%psi_out)
            ! Normalize and add to RHS
            psi_rhs(:, j) = this%psi_out(:) / dble(NR)
        end do

        ! Step 2: Diagonal potential contribution: V_diag * psi
        do j = 1, Nstates
            psi_rhs(:, j) = psi_rhs(:, j) + adb(:, j) * psi(:, j)
        end do

        ! Step 3: Off-diagonal dipole coupling: V_coupling(t) * psi
        ! at each grid point, apply the dipole matrix (no exponentiation)
        do i = 1, NR
            call pulse_direct(tout, mu_all(:, :, i), E_field)
            psi_Nstates(:) = psi(i, :)
            call zgemv('N', int(Nstates), int(Nstates), (1._dp, 0._dp), &
                & tout, size(tout, dim=1), psi_Nstates, 1, (0._dp, 0._dp), &
                & psi_Nstates1, 1)
            psi_rhs(i, :) = psi_rhs(i, :) + psi_Nstates1(:)
        end do

        ! Step 4: Multiply by -i to get d(psi)/dt = -i * H * psi
        psi_rhs = -im * psi_rhs

    end subroutine rhs_1d

    !> Perform one full RK4 time step: psi -> psi + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
    !! k1 = rhs(psi, t)
    !! k2 = rhs(psi + k1*dt/2, t + dt/2)
    !! k3 = rhs(psi + k2*dt/2, t + dt/2)
    !! k4 = rhs(psi + k3*dt, t + dt)
    subroutine rk4_step(this, psi, dt, E_field_now, E_field_half, E_field_next, mu_all, adb)
        use global_vars, only: NR, Nstates
        class(rk4_operator_type), intent(inout) :: this
        complex(dp), intent(inout) :: psi(NR, Nstates)
        real(dp), intent(in)      :: dt
        real(dp), intent(in)      :: E_field_now    ! E(t)
        real(dp), intent(in)      :: E_field_half   ! E(t + dt/2)
        real(dp), intent(in)      :: E_field_next   ! E(t + dt)
        real(dp), intent(in)      :: mu_all(Nstates, Nstates, NR)
        real(dp), intent(in)      :: adb(NR, Nstates)

        complex(dp), allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
        complex(dp), allocatable :: psi_tmp(:,:)

        allocate(k1(NR, Nstates), k2(NR, Nstates))
        allocate(k3(NR, Nstates), k4(NR, Nstates))
        allocate(psi_tmp(NR, Nstates))

        ! k1 = rhs(psi, t)        using E(t)
        call this%rhs_1d(psi, k1, E_field_now, mu_all, adb)

        ! k2 = rhs(psi + k1*dt/2, t + dt/2)   using E(t+dt/2)
        psi_tmp = psi + k1 * (0.5_dp * dt)
        call this%rhs_1d(psi_tmp, k2, E_field_half, mu_all, adb)

        ! k3 = rhs(psi + k2*dt/2, t + dt/2)   using E(t+dt/2)
        psi_tmp = psi + k2 * (0.5_dp * dt)
        call this%rhs_1d(psi_tmp, k3, E_field_half, mu_all, adb)

        ! k4 = rhs(psi + k3*dt,   t + dt)     using E(t+dt)
        psi_tmp = psi + k3 * dt
        call this%rhs_1d(psi_tmp, k4, E_field_next, mu_all, adb)

        ! Update: psi_new = psi + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
        psi = psi + (k1 + 2._dp * k2 + 2._dp * k3 + k4) * (dt / 6._dp)

        deallocate(k1, k2, k3, k4, psi_tmp)

    end subroutine rk4_step

    !> Build the dipole coupling matrix (non-exponentiated) for direct application
    !! tout(i,j) = -kap * mu(i,j) * E  (off-diagonal), 0 on diagonal
    subroutine pulse_direct(tout, mu, E)
        use global_vars, only: Nstates, kap
        integer :: i, j
        real(dp), intent(in)  :: mu(Nstates, Nstates), E
        complex(dp), intent(out) :: tout(Nstates, Nstates)

        tout = (0._dp, 0._dp)
        do i = 1, Nstates - 1
            do j = i + 1, Nstates
                tout(i, j) = -kap * mu(i, j) * E
                tout(j, i) = -kap * mu(i, j) * E
            end do
        end do

    end subroutine pulse_direct

    !> Clean up FFTW plans, memory, and kinetic energy array
    subroutine finalize(this)
        use FFTW3
        class(rk4_operator_type), intent(inout) :: this

        call fftw_destroy_plan(this%planF)
        call fftw_destroy_plan(this%planB)
        call fftw_free(this%p_in)
        call fftw_free(this%p_out)
        if (allocated(this%kin_energy)) deallocate(this%kin_energy)

    end subroutine finalize

end module rk4_operator_mod
