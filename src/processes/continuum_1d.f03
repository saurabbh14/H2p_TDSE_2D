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
        ! FFTW plan and memory pointers
        type(C_PTR) :: planF, planB, p_in, p_out
        ! FFTW input/output arrays
        complex(C_DOUBLE_COMPLEX), pointer :: psi_in(:), psi_out(:)
    contains
        procedure :: initialize          ! FFTW setup + kprop_full generation
        procedure :: apply               ! Absorber mask + free propagation of absorbed part
        procedure :: forward_fft         ! Forward FFT of a 1D array (for post-prop analysis)
    end type continuum_1d_type

contains

    !> Initialize FFTW plans and generate kinetic propagator
    subroutine initialize(this)
        use global_vars, only: NR, dt, m_red, PR, prop_par_FFTW
        use data_au, only: im
        use FFTW3
        class(continuum_1d_type), intent(inout) :: this

        print*
        print*, "Continuum/Absorber: FFTW initialization ..."

        ! Creating aligned memory for FFTW
        this%p_in = fftw_alloc_complex(int(NR, C_SiZE_T))
        call c_f_pointer(this%p_in, this%psi_in, [NR])
        this%p_out = fftw_alloc_complex(int(NR, C_SiZE_T))
        call c_f_pointer(this%p_out, this%psi_out, [NR])

        call fftw_initialize_threads
        print*, "Continuum/Absorber: FFTW plan creation ..."
        call fftw_create_c2c_plans(this%psi_in, this%psi_out, NR, &
            & this%planF, this%planB, prop_par_FFTW)

        ! Generate kinetic propagator for full time step
        allocate(this%kprop_full(NR))
        this%kprop_full = exp(-im * dt * pR * pR / (2._dp * m_red))

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

        allocate(psi_outR1(NR, Nstates))

        psi_outR1 = (0._dp, 0._dp)
        do J = 1, Nstates
            ! Propagate existing absorbed part in momentum space
            psi_outR(:, J) = psi_outR(:, J) * this%kprop_full(:) &
                & * exp(-im * dt * adb(NR - i_cpmR, J))

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

end module continuum_1d_mod