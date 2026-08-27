!> Full-2D Imaginary Time Propagation (ITP) module.
!>
!> Computes the lowest bound eigenstates of the full 2D (R, x) Hamiltonian
!> (nuclear coordinate R + electronic coordinate x) by propagating a trial
!> wavefunction in imaginary time.  
!> Results are stored in the global arrays psi_itp_2d(NR, Nx, nstates) and
!> E_itp_2d(nstates) and can be used as the initial state of the real-time
!> evolution (input: initial_state.distribution = "2d itp",
!> initial_state.itp_state = <N>).
module itp_2d_mod
    use varprecision, only: dp
    use global_vars, only: NR, Nx, Nstates, R, x, dR, dx, dt, kappa, &
        & pot, ewf, adb, ITP_par_FFTW, itp_2d_nstates, itp_2d_max_iter, &
        & itp_2d_thresh, itp_2d_dt_scale, itp_2d_guess, itp_2d_save, &
        & itp_2d_read, itp_2d_dir, psi_itp_2d, E_itp_2d
    use data_au, only: au2eV
    use split_operator_2d_mod, only: split_operator_2d_type
    implicit none
    private
    public :: itp_2d_type

    !> Class holding the configuration and working data of a 2D ITP run.
    type :: itp_2d_type
        private
        type(split_operator_2d_type) :: propagator ! FFTW plans + propagators (reused)
        logical :: plans_ready = .false.           ! FFTW plans were created
        real(dp) :: dt2                            ! imaginary time step (dt_scale * dt)
        integer  :: nstates                        ! number of eigenstates to compute
        integer  :: max_iter                       ! max iterations per state
        real(dp) :: thresh                         ! energy convergence threshold (a.u.)
        character(20) :: guess                     ! "auto" | "gaussian" | "ewf"
        logical :: save_results                    ! write eigenstates to disk
        logical :: read_results                    ! reuse existing files if present
    contains
        procedure :: itp_2d_calc         ! high-level driver
        procedure :: read_params
        procedure :: guess_state
        procedure :: compute_eigenstates
        procedure :: check_files_exist
        procedure :: read_eigenstates
        procedure :: write_eigenstates
        procedure :: finalize
    end type itp_2d_type

contains

    !> High-level driver: read parameters, reuse or compute the eigenstates,
    !! store them in the global arrays and (optionally) write them to disk.
    subroutine itp_2d_calc(this)
        class(itp_2d_type), intent(inout) :: this
        logical :: data_exists

        print*
        print*, "2D Imaginary Time Propagation (ITP) ..."
        call this%read_params()

        if (allocated(psi_itp_2d)) deallocate(psi_itp_2d)
        if (allocated(E_itp_2d)) deallocate(E_itp_2d)
        allocate(psi_itp_2d(NR, Nx, this%nstates))
        allocate(E_itp_2d(this%nstates))
        psi_itp_2d = 0._dp
        E_itp_2d   = 0._dp

        data_exists = .false.
        if (this%read_results) data_exists = this%check_files_exist()

        if (data_exists) then
            print*
            print*, "2D ITP data files found. Loading from disk..."
            call this%read_eigenstates()
            print*, "2D ITP data loaded successfully."
        else
            print*
            print*, "Running full 2D ITP ..."
            call this%compute_eigenstates()
            if (this%save_results) call this%write_eigenstates()
        end if

        call this%finalize()
        print*, "Leaving 2D ITP calculations."
    end subroutine itp_2d_calc

    !> Copy the 2D ITP settings from the global (input) variables.
    subroutine read_params(this)
        class(itp_2d_type), intent(inout) :: this
        this%nstates      = max(1, itp_2d_nstates)
        this%max_iter     = max(1, itp_2d_max_iter)
        this%thresh       = itp_2d_thresh
        this%dt2          = itp_2d_dt_scale * dt
        this%guess        = adjustl(trim(itp_2d_guess))
        this%save_results = itp_2d_save
        this%read_results = itp_2d_read
        if (this%dt2 <= 0._dp) then
            print*, "WARNING: [itp_2d] dt_scale <= 0; resetting to 0.1."
            this%dt2 = 0.1_dp * dt
        end if
        print*, "Number of 2D eigenstates to compute:", this%nstates
        print*, "ITP time step: dt_itp = ", sngl(this%dt2), " a.u. (dt_scale * dt)"
        print*, "Convergence threshold:", sngl(this%thresh), " a.u."
        print*, "Max iterations per state:", this%max_iter
        print*, "Initial guess strategy:", trim(this%guess)
    end subroutine read_params

    !> Build the initial guess for eigenstate N and normalize it.
    !!
    !! Both strategies alternate the parity of the startup function with the
    !! state index, exactly like the 1D ITP modules do (`nuclear_wv.f08` uses
    !! `exp(kappa(R-Rin)^2) + (-1)^(V-1) exp(-0.5(R+Rin)^2)` and
    !! `adiabatic.f08` uses two x-centred Gaussians with alternating sign).
    !! Without this the Gram-Schmidt projection would have to build the excited
    !! states out of numerical residue only, which converges very slowly.
    !!
    !! The state index N is mapped onto an electronic index (cycling through the
    !! Nstates electronic surfaces) and a nuclear excitation index:
    !!   N_el  = mod(N-1, Nstates) + 1      (electronic parity: gerade/ungerade)
    !!   v_idx = (N-1) / Nstates            (nuclear parity)
    !!
    !!   "gaussian" : two electronic Gaussians centred on the two nuclei
    !!                (+-mn1/mn2 * R, as in adiabatic.f08) times a nuclear
    !!                Gaussian centred at the R of the global potential minimum
    !!   "ewf"      : adiabatic electronic state ewf(:,:,N_el) times the nuclear
    !!                startup function of nuclear_wv.f08 (centred at the
    !!                equilibrium distance of BO surface N_el)
    !!   "auto"     : "ewf" if adiabatic wavefunctions are available, else "gaussian"
    subroutine guess_state(this, N, psi)
        use global_vars, only: mn1, mn2
        class(itp_2d_type), intent(inout) :: this
        integer, intent(in) :: N
        complex(dp), intent(out) :: psi(NR, Nx)
        integer :: j, N_el, v_idx, iR_in
        integer :: imin(2)
        real(dp) :: Rin, R_eq, norm, el_parity, nu_parity
        real(dp), allocatable :: chi_R(:)
        character(20) :: strategy
        logical :: ewf_available

        strategy      = adjustl(trim(this%guess))
        ewf_available = .false.
        if (allocated(ewf)) ewf_available = (maxval(abs(ewf)) > 0._dp)

        if (strategy == "auto") then
            if (ewf_available) then
                strategy = "ewf"
            else
                strategy = "gaussian"
            end if
        else if (strategy == "ewf" .and. .not. ewf_available) then
            print*, "WARNING: [itp_2d] guess = 'ewf' but no adiabatic wavefunctions"
            print*, "         are available; falling back to a Gaussian guess."
            strategy = "gaussian"
        end if

        ! Electronic / nuclear excitation indices and the associated parities
        N_el      = mod(N - 1, Nstates) + 1
        v_idx     = (N - 1) / Nstates
        el_parity = (-1._dp)**(N_el - 1)
        nu_parity = (-1._dp)**v_idx

        allocate(chi_R(NR))
        psi = (0._dp, 0._dp)
        select case(strategy)
        case("ewf")
            ! Adiabatic electronic state on surface N_el times the nuclear
            ! startup function centred at the equilibrium distance of that surface.
            iR_in = minloc(adb(:, N_el), 1)
            Rin   = R(iR_in)
            chi_R(:) = exp(kappa * (R(:) - Rin)**2) &
                &    + nu_parity * exp(-0.5_dp * (R(:) + Rin)**2)
            do j = 1, Nx
                psi(:, j) = ewf(:, j, N_el) * chi_R(:)
            end do

        case default   ! plain Gaussian guess centred at the potential minimum
            imin = minloc(pot)
            R_eq = R(imin(1))
            chi_R(:) = exp(kappa * (R(:) - R_eq)**2) &
                &    + nu_parity * exp(-0.5_dp * (R(:) + R_eq)**2)
            do j = 1, Nx
                psi(:, j) = ( exp(kappa * (x(j) - mn1 * R(:))**2) &
                    &      + el_parity * exp(kappa * (x(j) + mn2 * R(:))**2) ) &
                    &      * chi_R(:)
            end do
        end select
        deallocate(chi_R)

        call integ_2d(psi, psi, norm)
        if (norm <= 0._dp) then
            print*, "ERROR: 2D ITP initial guess for state", N, "has zero norm."
            stop
        end if
        psi = psi / sqrt(norm)
    end subroutine guess_state

    !> Successively converge the lowest `nstates` eigenstates of the full 2D
    !! Hamiltonian by imaginary-time split-operator propagation.
    subroutine compute_eigenstates(this)
        class(itp_2d_type), intent(inout) :: this
        integer :: N, G, K
        real(dp) :: E, E1, norm_old, norm_new, c
        complex(dp), allocatable :: psi(:,:), psi_old(:,:)
        complex(dp), allocatable :: ref(:,:,:)   ! converged states (R, x, state)
        character(20) :: in_xR

        allocate(psi(NR, Nx), psi_old(NR, Nx))
        allocate(ref(NR, Nx, this%nstates))
        ref = (0._dp, 0._dp)

        ! FFTW plans and imaginary-time propagators (reuse split_operator_2d)
        call this%propagator%fft_initialize(ITP_par_FFTW)
        call this%propagator%itp_initialize(this%dt2)
        this%plans_ready = .true.

        write(in_xR, '(a)') "inner-xR"

        state_loop: do N = 1, this%nstates
            call this%guess_state(N, psi)
            E = 0._dp
            print*
            print*, "2D ITP state", N, ": converging ..."

            iter_loop: do K = 1, this%max_iter
                psi_old = psi
                E1 = E

                ! Project out previously converged states (Gram-Schmidt)
                if (N > 1) then
                    do G = 1, N - 1
                        call integ_2d(ref(:,:,G), psi, c)
                        psi = psi - c * ref(:,:,G)
                    end do
                end if

                ! One imaginary-time split-operator step
                call this%propagator%split_operator_step(psi, in_xR)

                ! Eigenvalue estimate from the norm change over the step
                call integ_2d(psi_old, psi_old, norm_old)
                call integ_2d(psi, psi, norm_new)
                E = -0.5_dp / this%dt2 * log(norm_new / norm_old)

                ! Renormalize
                psi = psi / sqrt(norm_new)

                if (abs(E - E1) <= this%thresh) then
                    print*, "  State", N, "converged after", K, "iterations."
                    print*, "  E =", sngl(E), " a.u. =", sngl(E * au2eV), " eV"
                    ref(:,:,N)  = psi
                    E_itp_2d(N) = E
                    cycle state_loop
                else if (K == this%max_iter) then
                    print*, "ERROR: 2D ITP state", N, "did not converge after", &
                        & K, "iterations."
                    print*, "       |E - E_old| =", sngl(abs(E - E1)), &
                        & " a.u. (thresh =", sngl(this%thresh), " a.u.)"
                    print*, "       Increase [itp_2d] max_iter or decrease dt_scale."
                    stop
                end if
            end do iter_loop
        end do state_loop

        ! Store the converged states in the global array
        do N = 1, this%nstates
            psi_itp_2d(:,:,N) = ref(:,:,N)
        end do

        deallocate(psi, psi_old, ref)
    end subroutine compute_eigenstates

    !> Check that all 2D ITP result files exist and are consistent with the
    !! current grid dimensions and requested number of states.
    function check_files_exist(this) result(exists)
        class(itp_2d_type), intent(inout) :: this
        logical :: exists
        integer :: fd, ios, N, i_dummy
        integer :: nr_chk, nx_chk
        real(dp) :: dummy
        character(500) :: fn

        exists = .true.

        ! 1. Energies file with at least nstates entries (after the header)
        write(fn, '(a,a)') adjustl(trim(itp_2d_dir)), 'energies.out'
        open(newunit=fd, file=trim(fn), status='old', form='formatted', iostat=ios)
        if (ios /= 0) then
            print*, '  Missing: ', trim(fn)
            exists = .false.
        else
            read(fd, *, iostat=ios)  ! skip header line
            do N = 1, this%nstates
                read(fd, *, iostat=ios) i_dummy, dummy
                if (ios /= 0) exit
            end do
            close(fd)
            if (ios /= 0) then
                print*, '  Incomplete: ', trim(fn), ' (expected at least ', &
                    & this%nstates, ' states)'
                exists = .false.
            end if
        end if

        ! 2. Wavefunction binary for every state
        do N = 1, this%nstates
            write(fn, '(a,a,i0,a)') adjustl(trim(itp_2d_dir)), &
                & 'psi_itp_state_', N, '.bin'
            open(newunit=fd, file=trim(fn), status='old', form='unformatted', &
                & access='stream', iostat=ios)
            if (ios /= 0) then
                print*, '  Missing: ', trim(fn)
                exists = .false.
            else
                read(fd, iostat=ios) nr_chk, nx_chk
                close(fd)
                if (ios /= 0 .or. nr_chk /= NR .or. nx_chk /= Nx) then
                    print*, '  Incompatible grid dimensions in: ', trim(fn)
                    print*, '  Expected NR,Nx:', NR, Nx, ' Found:', nr_chk, nx_chk
                    exists = .false.
                end if
            end if
        end do

        if (.not. exists) then
            print*, 'Some 2D ITP data files are missing or invalid. Will run full ITP.'
        end if
    end function check_files_exist

    !> Load eigenstates and energies from itp_2d_dir.
    subroutine read_eigenstates(this)
        class(itp_2d_type), intent(inout) :: this
        integer :: fd, N, i_dummy
        integer :: nr_chk, nx_chk
        real(dp) :: E
        character(500) :: fn

        write(fn, '(a,a)') adjustl(trim(itp_2d_dir)), 'energies.out'
        open(newunit=fd, file=trim(fn), status='old', form='formatted')
        read(fd, *)  ! skip header
        do N = 1, this%nstates
            read(fd, *) i_dummy, E
            E_itp_2d(N) = E
        end do
        close(fd)

        do N = 1, this%nstates
            write(fn, '(a,a,i0,a)') adjustl(trim(itp_2d_dir)), &
                & 'psi_itp_state_', N, '.bin'
            open(newunit=fd, file=trim(fn), status='old', form='unformatted', &
                & access='stream')
            read(fd) nr_chk, nx_chk
            read(fd) psi_itp_2d(:,:,N)
            close(fd)
        end do

        print*
        print*, 'Loaded', this%nstates, '2D ITP state(s) from ', trim(itp_2d_dir)
        do N = 1, this%nstates
            print*, '  State', N, ': E =', sngl(E_itp_2d(N)), ' a.u. =', &
                & sngl(E_itp_2d(N) * au2eV), ' eV'
        end do
    end subroutine read_eigenstates

    !> Write eigenstates and energies to itp_2d_dir.
    subroutine write_eigenstates(this)
        class(itp_2d_type), intent(inout) :: this
        integer :: fd, N
        character(500) :: fn

        print*
        print*, 'Writing 2D ITP eigenstates to ', trim(itp_2d_dir)

        write(fn, '(a,a)') adjustl(trim(itp_2d_dir)), 'energies.out'
        open(newunit=fd, file=trim(fn), status='replace', form='formatted')
        write(fd, '(a)') '# State   E (a.u.)            E (eV)'
        do N = 1, this%nstates
            write(fd, '(i5, 2ES20.10)') N, E_itp_2d(N), E_itp_2d(N) * au2eV
        end do
        close(fd)

        do N = 1, this%nstates
            write(fn, '(a,a,i0,a)') adjustl(trim(itp_2d_dir)), &
                & 'psi_itp_state_', N, '.bin'
            open(newunit=fd, file=trim(fn), status='replace', form='unformatted', &
                & access='stream')
            write(fd) NR, Nx
            write(fd) psi_itp_2d(:,:,N)
            close(fd)
            print*, '  saved: ', trim(fn)
        end do
        print*, 'Done writing 2D ITP data.'
    end subroutine write_eigenstates

    !> Clean up FFTW resources (only if they were created).
    subroutine finalize(this)
        class(itp_2d_type), intent(inout) :: this
        if (this%plans_ready) call this%propagator%finalize()
        this%plans_ready = .false.
    end subroutine finalize

    !> Real part of the grid inner product <a|b> = sum(a* b) * dR * dx.
    !! For (numerically real) ITP states this is the usual overlap integral.
    subroutine integ_2d(a, b, out)
        complex(dp), intent(in) :: a(NR, Nx), b(NR, Nx)
        real(dp), intent(out) :: out
        out = sum(real(a, dp) * real(b, dp) + aimag(a) * aimag(b)) * dR * dx
    end subroutine integ_2d

end module itp_2d_mod
