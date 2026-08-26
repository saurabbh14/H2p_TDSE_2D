module continuum_1d_mod
    use global_vars, only: dp
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: continuum_1d_type

    !> Standalone type for continuum/absorber propagation in 1D
    !! Provides FFTW plans, kinetic propagator, and absorber application.
    !! Completely independent of any time-propagation scheme.
    type :: continuum_1d_type
        ! Kinetic propagator (free particle, full time step)
        complex(dp), allocatable :: kprop_full(:)
        ! Bare kinetic energy pR²/(2*m_red) — used by the RK4 right-hand side
        real(dp), allocatable :: kin_energy(:)
        ! time-evolution scheme: "split_operator" (default) or "rk4"
        character(20) :: evol_mode = "split_operator"
        ! FFTW plan and memory pointers
        type(C_PTR) :: planF, planB, p_in, p_out
        ! FFTW input/output arrays
        complex(C_DOUBLE_COMPLEX), pointer :: psi_in(:), psi_out(:)
    contains
        procedure :: initialize          ! FFTW setup + kprop_full generation
        procedure :: apply               ! Absorber mask + free propagation of absorbed part
        procedure :: rhs => rhs_outR     ! RHS of the TDSE for the absorbed packet
        procedure :: rk4_step => rk4_outR ! One RK4 step of the absorbed packet
        procedure :: forward_fft         ! Forward FFT of a 1D array (for post-prop analysis)
        procedure :: finalize            ! Clean up FFTW plans, memory, and kprop_full
    end type continuum_1d_type

contains

    !> Initialize FFTW plans and generate kinetic propagator
    !! `propagator` selects the evolution scheme for the absorbed wave packet and
    !! is normally the main propagator name ("split_operator" | "rk4").
    subroutine initialize(this, propagator)
        use global_vars, only: NR, dt, m_red, PR, prop_par_FFTW
        use data_au, only: im
        use FFTW3
        class(continuum_1d_type), intent(inout) :: this
        character(*), intent(in), optional :: propagator

        print*
        print*, "Continuum/Absorber: FFTW initialization ..."

        ! Time-evolution scheme: follow the main propagator; anything unknown
        ! falls back to the (here exact) exponential of the diagonal Hamiltonian.
        this%evol_mode = "split_operator"
        if (present(propagator)) then
            if (trim(adjustl(propagator)) == "rk4") this%evol_mode = "rk4"
        end if

        ! Creating aligned memory for FFTW
        this%p_in = fftw_alloc_complex(int(NR, C_SiZE_T))
        call c_f_pointer(this%p_in, this%psi_in, [NR])
        this%p_out = fftw_alloc_complex(int(NR, C_SiZE_T))
        call c_f_pointer(this%p_out, this%psi_out, [NR])

        call fftw_initialize_threads
        print*, "Continuum/Absorber: FFTW plan creation ..."
        call fftw_create_c2c_plans(this%psi_in, this%psi_out, NR, &
            & this%planF, this%planB, prop_par_FFTW)

        ! Kinetic energy and kinetic propagator for a full time step
        allocate(this%kin_energy(NR))
        this%kin_energy(:) = pR(:) * pR(:) / (2._dp * m_red)
        allocate(this%kprop_full(NR))
        this%kprop_full = exp(-im * dt * this%kin_energy(:))

        if (trim(adjustl(this%evol_mode)) == "rk4") then
            print*, "Continuum/Absorber: evolution mode = RK4"
        else
            print*, "Continuum/Absorber: evolution mode = split-operator"
        end if

        print*, "Continuum/Absorber: Done."

    end subroutine initialize

    !> Apply absorber and propagate the absorbed wavefunction in momentum space.
    !! On input:  psi_ges, psi_outR (accumulated absorbed wavefunction)
    !! On output: psi_ges = psi_ges * abs_func (masked, bound part)
    !!             psi_outR = psi_outR * kprop_full + FFT_to_momentum(psi_ges*(1-abs_func))
    subroutine apply(this, psi_ges, psi_outR, psi_outR_inc, &
                       abs_func, i_cpmR)
        use global_vars, only: NR, Nstates, dt, adb
        use data_au, only: im
        use FFTW3
        class(continuum_1d_type), intent(inout) :: this
        complex(dp), intent(inout) :: psi_ges(NR, Nstates)
        complex(dp), intent(inout) :: psi_outR(NR, Nstates)
        real(dp), intent(inout)    :: psi_outR_inc(NR, Nstates)
        complex(dp), intent(in)    :: abs_func(NR)
        integer, intent(in)        :: i_cpmR

        integer :: J
        complex(dp), allocatable :: psi_outR1(:,:)
        real(dp) :: adb_b(Nstates)

        allocate(psi_outR1(NR, Nstates))

        ! Constant boundary potential of each electronic state at the absorber edge
        do J = 1, Nstates
            adb_b(J) = adb(NR - i_cpmR, J)
        end do

        ! Propagate the already-absorbed part in momentum space.
        ! H_J = pR²/(2*m_red) + adb(NR-i_cpmR, J) is diagonal in pR, so the
        ! exponential below is exact; the RK4 branch exists only to integrate the
        ! continuum with the same scheme as the main (RK4) propagation.
        select case (trim(adjustl(this%evol_mode)))
        case ("rk4")
            call this%rk4_step(psi_outR, adb_b)
        case default
            do J = 1, Nstates
                psi_outR(:, J) = psi_outR(:, J) * this%kprop_full(:) &
                    & * exp(-im * dt * adb_b(J))
            end do
        end select

        psi_outR1 = (0._dp, 0._dp)
        do J = 1, Nstates
            ! Extract new absorbed part: psi * (1 - absorber)
            psi_outR1(:, J) = psi_ges(:, J) * (1._dp - abs_func(:))

            ! Apply mask to keep only the bound part
            psi_ges(:, J) = psi_ges(:, J) * abs_func(:)
        end do

        ! FFT the newly absorbed part to momentum space for spectrum analysis
        do J = 1, Nstates
            this%psi_in = (0._dp, 0._dp)
            this%psi_out = (0._dp, 0._dp)
            this%psi_in(:) = psi_outR1(:, J)
            call fftw_execute_dft(this%planF, this%psi_in, this%psi_out)
            this%psi_in = this%psi_out / sqrt(dble(NR))
            psi_outR1(:, J) = this%psi_in(:)
        end do

        ! Accumulate
        psi_outR = psi_outR + psi_outR1
        psi_outR_inc = psi_outR_inc + abs(psi_outR1)**2

        deallocate(psi_outR1)

    end subroutine apply

    !> Evaluate the RHS of the TDSE for the absorbed (momentum-space) packet:
    !!   psi_rhs(:,J) = -i * ( pR²/(2*m_red) + adb_b(J) ) * psi(:,J)
    !! The Hamiltonian is diagonal in pR, so no FFT is needed here.
    subroutine rhs_outR(this, psi, psi_rhs, adb_b)
        use global_vars, only: NR, Nstates
        use data_au, only: im
        class(continuum_1d_type), intent(inout) :: this
        complex(dp), intent(in)  :: psi(NR, Nstates)
        complex(dp), intent(out) :: psi_rhs(NR, Nstates)
        real(dp), intent(in)     :: adb_b(Nstates)

        integer :: J

        do J = 1, Nstates
            psi_rhs(:, J) = -im * (this%kin_energy(:) + adb_b(J)) * psi(:, J)
        end do

    end subroutine rhs_outR

    !> One full RK4 step of the absorbed packet:
    !! psi -> psi + (k1 + 2*k2 + 2*k3 + k4) * dt/6.
    !! H is time independent here, so all four stages use the same adb_b.
    subroutine rk4_outR(this, psi_outR, adb_b)
        use global_vars, only: NR, Nstates, dt
        class(continuum_1d_type), intent(inout) :: this
        complex(dp), intent(inout) :: psi_outR(NR, Nstates)
        real(dp), intent(in)       :: adb_b(Nstates)

        complex(dp), allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
        complex(dp), allocatable :: psi_tmp(:,:)

        allocate(k1(NR, Nstates), k2(NR, Nstates))
        allocate(k3(NR, Nstates), k4(NR, Nstates))
        allocate(psi_tmp(NR, Nstates))

        ! k1 = rhs(psi)
        call this%rhs(psi_outR, k1, adb_b)
        ! k2 = rhs(psi + k1*dt/2)
        psi_tmp = psi_outR + k1 * (0.5_dp * dt)
        call this%rhs(psi_tmp, k2, adb_b)
        ! k3 = rhs(psi + k2*dt/2)
        psi_tmp = psi_outR + k2 * (0.5_dp * dt)
        call this%rhs(psi_tmp, k3, adb_b)
        ! k4 = rhs(psi + k3*dt)
        psi_tmp = psi_outR + k3 * dt
        call this%rhs(psi_tmp, k4, adb_b)

        psi_outR = psi_outR + (k1 + 2._dp * k2 + 2._dp * k3 + k4) * (dt / 6._dp)

        deallocate(k1, k2, k3, k4, psi_tmp)

    end subroutine rk4_outR

    !> Perform a forward FFT on a single-state array (for post-propagation KER analysis).
    !! Used by post_prop_analysis to transform dissociated wavefunction to momentum space.
    subroutine forward_fft(this, psi_in_col, psi_out_col)
        use global_vars, only: NR
        use FFTW3
        class(continuum_1d_type), intent(inout) :: this
        complex(dp), intent(in)  :: psi_in_col(NR)
        complex(dp), intent(out) :: psi_out_col(NR)

        this%psi_in = (0._dp, 0._dp)
        this%psi_out = (0._dp, 0._dp)
        this%psi_in(:) = psi_in_col(:)
        call fftw_execute_dft(this%planF, this%psi_in, this%psi_out)
        psi_out_col(:) = this%psi_out(:) / sqrt(dble(NR))

    end subroutine forward_fft

    !> Clean up FFTW plans, memory, and kinetic propagator
    subroutine finalize(this)
        use FFTW3
        class(continuum_1d_type), intent(inout) :: this

        call fftw_destroy_plan(this%planF)
        call fftw_destroy_plan(this%planB)
        call fftw_free(this%p_in)
        call fftw_free(this%p_out)
        if (allocated(this%kprop_full)) deallocate(this%kprop_full)
        if (allocated(this%kin_energy)) deallocate(this%kin_energy)

    end subroutine finalize

end module continuum_1d_mod
