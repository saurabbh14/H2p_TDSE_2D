module adiabatic_mod
    use iso_c_binding, only: C_DOUBLE, C_PTR, C_SIZE_T
    use varprecision, only: dp
    use global_vars, only: NR, Nx, Nstates, R, x, Px, dR, dx, dt, &
        & Pot, adb, ewf, mu_all, output_data_dir, adiabatic_dir, ITP_par_FFTW
    use data_au, only: au2eV
    use FFTW3
    use omp_lib
    implicit none
    private
    public :: adiabatic_wavefkt_class, adiabatic_wf_calc

    type :: adiabatic_wavefkt_class
        character(10) :: ITP_par_FFTW
        integer :: max_iter
        real(dp) :: thresh
        real(dp), allocatable :: ref(:,:,:)
        real(dp), allocatable :: d(:,:)
    contains
        procedure :: read_params
        procedure :: initialize_wp_params
        procedure :: imaginary_time_propagation => adiabatic_ITP
        procedure :: post_itp_calculations
        procedure :: compute_transition_dipoles
        procedure :: compute_nonadiabatic_couplings
        procedure :: write_output_files
        procedure :: write_bo_surface
        procedure :: write_transition_dipole_files
        procedure :: write_nonadiabatic_coupling_files
        procedure :: write_ewf_binary
        procedure :: check_adiabatic_data_exists
        procedure :: read_adiabatic_data
        procedure :: adiabatic_wf_calc
    end type adiabatic_wavefkt_class

contains

    subroutine adiabatic_wf_calc(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        logical :: data_exists
        call this%read_params()
        call this%initialize_wp_params()

        data_exists = this%check_adiabatic_data_exists()
        if (data_exists) then
            print*, 'Adiabatic data files found. Loading from disk...'
            call this%read_adiabatic_data()
            print*, 'Adiabatic data loaded successfully.'
        else
            print*, 'Adiabatic data files not found or incomplete. Running full ITP...'
            call this%imaginary_time_propagation()
            call this%post_itp_calculations()
            call this%write_output_files()
            call this%write_ewf_binary()
        end if
    end subroutine adiabatic_wf_calc

    subroutine read_params(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        this%ITP_par_FFTW = ITP_par_FFTW
        this%thresh = 1.d-15
        this%max_iter = 1000000
    end subroutine read_params

    subroutine initialize_wp_params(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        allocate(this%ref(Nx, NR, Nstates))
        allocate(this%d(NR, Nstates))
    end subroutine initialize_wp_params

    subroutine adiabatic_ITP(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        integer :: I, J, K, N, G
        type(C_PTR) :: p_in, p_out
        type(C_PTR) :: planF, planB
        real(dp) :: dt2, E, E1, norm
        real(dp), allocatable :: vpropx(:), kprop(:)
        real(C_DOUBLE), allocatable :: psi(:), psi1(:)
        real(C_DOUBLE), pointer :: psi_in(:), psi_out(:)

        allocate(psi(Nx), psi1(Nx), vpropx(Nx), kprop(Nx))
        call fftw_initialize_threads()
        p_in = fftw_alloc_real(int(Nx, C_SIZE_T))
        call c_f_pointer(p_in, psi_in, [Nx])
        p_out = fftw_alloc_real(int(Nx, C_SIZE_T))
        call c_f_pointer(p_out, psi_out, [Nx])
        call fftw_create_r2r_plans(psi_in, psi_out, Nx, planF, planB, this%ITP_par_FFTW)

        dt2 = dt / 10
        do I = 1, Nx
            kprop(I) = exp(-dt2 * (Px(I)**2) * 0.5_dp)
        end do

        do N = 1, Nstates
            write(*,*) 'Imaginary Time Propagation, State', N
            do I = 1, NR
                if (I == 1) then
                    psi(:) = exp(-10._dp * (x(:) - 0.5_dp * R(I))**2) &
                        & + (-1._dp)**(N-1) * exp(-10._dp * (x(:) + 0.5_dp * R(I))**2)
                else
                    psi(:) = this%ref(:, I-1, N)
                end if
                call integ_real(psi, psi, norm)
                psi = psi / sqrt(norm)
                vpropx(:) = exp(-0.5_dp * dt2 * pot(I, :))
                E = 0._dp
                do K = 1, this%max_iter
                    psi1 = psi
                    E1 = E
                    if (N > 1) then
                        do G = 1, N - 1
                            call integ_real(this%ref(1:Nx, I, G), psi, norm)
                            psi = psi - norm * this%ref(:, I, G)
                        end do
                    end if
                    psi = psi * vpropx
                    psi_in = psi
                    call fftw_execute_r2r(planF, psi_in, psi_out)
                    psi = psi_out
                    psi = psi * kprop
                    psi_in = psi
                    call fftw_execute_r2r(planB, psi_in, psi_out)
                    psi = psi_out
                    psi = psi / dble(Nx)
                    psi = psi * vpropx
                    call eigenvalue_real(psi, psi1, E, dt2)
                    call integ_real(psi, psi, norm)
                    psi = psi / sqrt(norm)
                    if (abs(E - E1) <= this%thresh) then
                        adb(I, N) = E
                        this%ref(:, I, N) = psi
                        exit
                    else if (K == this%max_iter) then
                        write(*,*) 'Iteration not converged for state', N, 'at R-index', I
                        stop
                    end if
                end do
            end do
        end do

        ewf = 0._dp
        do I = 1, NR
            do J = 1, Nx
                do N = 1, Nstates
                    ewf(I, J, N) = this%ref(J, I, N)
                end do
            end do
        end do

        call fftw_destroy_plan(planF)
        call fftw_destroy_plan(planB)
        call fftw_free(p_in)
        call fftw_free(p_out)
        deallocate(psi, psi1, vpropx, kprop)
    end subroutine adiabatic_ITP

    subroutine post_itp_calculations(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        call this%compute_transition_dipoles()
        call this%compute_nonadiabatic_couplings()
    end subroutine post_itp_calculations

    subroutine compute_transition_dipoles(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        integer :: I, L, M
        do L = 1, Nstates
            do M = L + 1, Nstates
                do I = 1, NR
                    mu_all(L, M, I) = sum(this%ref(:, I, L) * x(:) * this%ref(:, I, M)) * dx
                end do
            end do
        end do
    end subroutine compute_transition_dipoles

    subroutine compute_nonadiabatic_couplings(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        integer :: I, L
        do I = 1, NR
            do L = 1, Nstates
                call gradient_R(this%ref(:, I, L), this%d(:, L))
            end do
        end do
    end subroutine compute_nonadiabatic_couplings

    subroutine write_output_files(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        call this%write_bo_surface()
        call this%write_transition_dipole_files()
        call this%write_nonadiabatic_coupling_files()
    end subroutine write_output_files

    subroutine write_bo_surface(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        integer :: I
        integer :: fd
        if (.not. allocated(this%ref)) return
        open(newunit=fd, file=adjustl(trim(output_data_dir)) // '/H2+_BO.dat', status='replace', form='formatted')
        write(fd, '(a)') '# Internuclear Distance (a.u.)    Potential Energy (a.u.)'
        do I = 1, NR
            write(fd,*) R(I), adb(I, :)
        end do
        close(fd)
    end subroutine write_bo_surface

    subroutine write_transition_dipole_files(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        integer :: I, L, M
        integer :: fd
        character(200) :: fn
        if (.not. allocated(this%d)) return
        do L = 1, Nstates
            do M = L + 1, Nstates
                write(fn, '(a,a,i0,a,i0,a)') adjustl(trim(output_data_dir)), '/transition-dipole_', L, '-', M, '.out'
                open(newunit=fd, file=trim(fn), status='replace', form='formatted')
                write(fd,*) '# Internuclear Distance (a.u.)    Transition Dipole moment'
                do I = 1, NR
                    write(fd,*) R(I), mu_all(L, M, I)
                end do
                close(fd)
            end do
        end do
    end subroutine write_transition_dipole_files

    subroutine write_nonadiabatic_coupling_files(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        integer :: I, L, M
        integer :: fd
        character(200) :: fn
        do L = 1, Nstates
            do M = L + 1, Nstates
                write(fn, '(a,a,i0,a,i0,a)') adjustl(trim(output_data_dir)), '/Non-adiabatic_coupling_', L, '-', M, '.out'
                open(newunit=fd, file=trim(fn), status='replace', form='formatted')
                write(fd,*) '# Internuclear Distance (a.u.)    Non-adiabatic coupling'
                do I = 1, NR
                    write(fd,*) R(I), sum(this%ref(:, I, L) * this%d(:, M)) * dx
                end do
                close(fd)
            end do
        end do
    end subroutine write_nonadiabatic_coupling_files

    subroutine write_ewf_binary(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        integer :: fd
        character(500) :: fn
        write(fn, '(a,a)') adjustl(trim(adiabatic_dir)), 'ewf.bin'
        open(newunit=fd, file=trim(fn), status='replace', form='unformatted', access='stream')
        write(fd) NR, Nx, Nstates
        write(fd) ewf
        close(fd)
        print*, 'Electronic wavefunctions saved to ', trim(fn)
    end subroutine write_ewf_binary

    function check_adiabatic_data_exists(this) result(exists)
        class(adiabatic_wavefkt_class), intent(in) :: this
        logical :: exists
        integer :: fd, L, M, ios
        integer :: nr_chk, nx_chk, ns_chk
        character(500) :: fn
        logical :: file_ok

        exists = .true.

        ! 1. Check ewf binary file
        write(fn, '(a,a)') adjustl(trim(adiabatic_dir)), 'ewf.bin'
        open(newunit=fd, file=trim(fn), status='old', form='unformatted', access='stream', iostat=ios)
        if (ios /= 0) then
            print*, '  Missing: ', trim(fn)
            exists = .false.
        else
            read(fd, iostat=ios) nr_chk, nx_chk, ns_chk
            close(fd)
            if (ios /= 0 .or. nr_chk /= NR .or. nx_chk /= Nx .or. ns_chk /= Nstates) then
                print*, '  Incompatible grid dimensions in: ', trim(fn)
                print*, '  Expected NR,Nx,Nstates:', NR, Nx, Nstates, &
                    & ' Found:', nr_chk, nx_chk, ns_chk
                exists = .false.
            end if
        end if

        ! 2. Check H2+_BO.dat
        write(fn, '(a,a)') adjustl(trim(output_data_dir)), '/H2+_BO.dat'
        inquire(file=trim(fn), exist=file_ok)
        if (.not. file_ok) then
            print*, '  Missing: ', trim(fn)
            exists = .false.
        end if

        ! 3. Check all transition dipole files
        do L = 1, Nstates
            do M = L + 1, Nstates
                write(fn, '(a,a,i0,a,i0,a)') adjustl(trim(output_data_dir)), &
                    & '/transition-dipole_', L, '-', M, '.out'
                inquire(file=trim(fn), exist=file_ok)
                if (.not. file_ok) then
                    print*, '  Missing: ', trim(fn)
                    exists = .false.
                end if
            end do
        end do

        if (.not. exists) then
            print*, 'Some adiabatic data files are missing or invalid. Will run full ITP.'
        end if
    end function check_adiabatic_data_exists

    subroutine read_adiabatic_data(this)
        class(adiabatic_wavefkt_class), intent(inout) :: this
        integer :: fd, I, L, M, J, N
        integer :: nr_chk, nx_chk, ns_chk
        character(500) :: fn
        real(dp) :: dummy

        ! 1. Read ewf from binary file
        write(fn, '(a,a)') adjustl(trim(adiabatic_dir)), 'ewf.bin'
        open(newunit=fd, file=trim(fn), status='old', form='unformatted', access='stream')
        read(fd) nr_chk, nx_chk, ns_chk
        read(fd) ewf
        close(fd)

        ! Populate this%ref from ewf (ref is indexed as ref(x, R, state))
        do I = 1, NR
            do J = 1, Nx
                do N = 1, Nstates
                    this%ref(J, I, N) = ewf(I, J, N)
                end do
            end do
        end do

        ! 2. Read adb from H2+_BO.dat
        write(fn, '(a,a)') adjustl(trim(output_data_dir)), '/H2+_BO.dat'
        open(newunit=fd, file=trim(fn), status='old', form='formatted')
        read(fd, *)  ! skip header
        do I = 1, NR
            read(fd, *) dummy, adb(I, :)
        end do
        close(fd)

        ! 3. Read mu_all from transition dipole files
        do L = 1, Nstates
            do M = L + 1, Nstates
                write(fn, '(a,a,i0,a,i0,a)') adjustl(trim(output_data_dir)), &
                    & '/transition-dipole_', L, '-', M, '.out'
                open(newunit=fd, file=trim(fn), status='old', form='formatted')
                read(fd, *)  ! skip header
                do I = 1, NR
                    read(fd, *) dummy, mu_all(L, M, I)
                end do
                close(fd)
            end do
        end do

        ! 4. Recompute nonadiabatic couplings from ref
        call this%compute_nonadiabatic_couplings()
    end subroutine read_adiabatic_data

    subroutine eigenvalue_real(A, B, E, dt2)
        use global_vars, only: Nx
        implicit none
        real(dp), intent(in) :: A(Nx), B(Nx), dt2
        real(dp), intent(out) :: E
        real(dp) :: e1, e2, norm
        call integ_real(B, B, norm)
        e1 = norm
        call integ_real(A, A, norm)
        e2 = norm
        E = (-0.5_dp / dt2) * log(e2 / e1)
    end subroutine eigenvalue_real

    subroutine integ_real(A, B, C)
        use global_vars, only: Nx, dx
        implicit none
        real(dp), intent(in) :: A(Nx), B(Nx)
        real(dp), intent(out) :: C
        C = 0._dp
        C = 0.5_dp * (A(1) * B(1) + A(Nx) * B(Nx))
        C = C + sum(A(2:Nx-1) * B(2:Nx-1))
        C = C * dx
    end subroutine integ_real

    subroutine gradient_R(F, d)
        use global_vars, only: NR, dR
        implicit none
        real(dp), intent(in) :: F(NR)
        real(dp), intent(out) :: d(NR)
        integer :: I
        do I = 2, NR - 1
            d(I) = (F(I + 1) - F(I - 1)) / (2._dp * dR)
        end do
        d(1) = (-3._dp * F(1) + 4._dp * F(2) - F(3)) / (2._dp * dR)
        d(NR) = (3._dp * F(NR) - 4._dp * F(NR-1) + F(NR-2)) / (2._dp * dR)
    end subroutine gradient_R
end module adiabatic_mod
