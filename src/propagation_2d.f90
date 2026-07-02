module propagation2d_mod
    use iso_c_binding, only: C_DOUBLE, C_PTR, C_SIZE_T
    use varprecision, only: dp
    use global_vars, only: Nt, Nx, NR, Nstates, guess_vstates, &
        & R, x, dR, dx, dt, m_eff, m_red, pR, Px, kap, lam, time
    use data_au, only: au2eV
    use split_operator_2d_mod, only: split_operator_2d_type
    use rk4_operator_2d_mod, only: rk4_operator_2d_type
    use setpot_mod, only: build_kh_potential_at_time
    use FFTW3
    use omp_lib
    implicit none
    private
    public :: time_prop_2d, propagation_2D

    !> Type to hold all data and methods for 2D time propagation
    type :: time_prop_2d
        logical :: initialized = .false.

        ! Wavefunctions and Propagators
        complex(dp), allocatable :: psi(:,:)
        complex(dp), allocatable :: psi0(:,:)
        complex(dp), allocatable :: psi_end(:,:)

        ! Index for absorber placement
        integer :: i_cpmR, i_cpmx
        ! Absorber function (complex or mask)
        complex(dp), allocatable :: abs_R(:), abs_x(:)
        
        ! Absorbed/Out wavefunctions (x and R dimensions)
        complex(dp), allocatable :: psi_out_x(:,:), psi_out_R(:,:)

        ! Bound and Dissociated wavefunctions
        complex(dp), allocatable :: psi_bound(:,:), psi_diss(:,:)
        complex(dp), allocatable :: psi_diss_g(:,:), psi_diss_u(:,:)

        ! Densities
        real(dp), allocatable :: idensR(:), idensx(:)
        real(dp), allocatable :: idenspR(:), idenspx(:)
        
        ! FFTW Plans (Using integer*8 if using legacy dfftw_ or type(C_PTR))
        integer*8 :: planF, planB, planFd, planBd, planBx, planFx 
        integer*8 :: planBR, planFR

        ! File Handles
        integer :: psi_2d_tk, dens_R_tk, dens_x_tk, norm_2d_tk
        integer :: avgR_2d_tk, avgx_2d_tk, momt_2d_tk, norm_pn_2d_tk
        integer :: field_2d_tk
        integer :: abs_R_tk, abs_x_tk
        integer :: psi_outR_norm_2d_tk, psi_outR_Pdens_2d_tk
        
    contains
        ! Core API equivalent to 1D
        procedure :: propagation_2D
        
        ! Lifecycle & Setup
        procedure :: initialize
        procedure :: ini_dist_choice
        procedure :: absorber_gen
        
        ! File Management
        procedure :: open_files_to_write
        !procedure :: write_headers_to_files
        
        ! Physics / Execution
        procedure :: time_evolution
        !procedure :: expected_position
        !procedure :: post_prop_analysis
        
        ! Cleanup
        procedure :: deallocate_all
    end type time_prop_2d

contains

    subroutine propagation_2D(this, E2, A, alpha_t, propagator_method)
        class(time_prop_2d), intent(inout) :: this
        real(dp), intent(in) :: E2(:), A(:), alpha_t(:)
        character(*), intent(in) :: propagator_method

        if (.not. this%initialized) call this%initialize()
        
        ! I/O Setup
        call this%open_files_to_write()
        !call this%write_headers_to_files()

        ! Setup phase
        call this%ini_dist_choice()
        call this%absorber_gen()
        
        ! Main Loop
        call this%time_evolution(E2, A, alpha_t, propagator_method)
        
        ! Teardown & Analysis
        !call this%post_prop_analysis()
        call this%deallocate_all()
    end subroutine propagation_2D

    subroutine initialize(this)
        use global_vars, only: Nx, NR
        class(time_prop_2d), intent(inout) :: this
        
        ! Allocate main grids based on global dimensions
        allocate(this%psi(NR, Nx))
        allocate(this%psi0(NR, Nx))
        
        allocate(this%psi_out_x(NR, Nx))
        allocate(this%psi_out_R(NR, Nx))
        
        allocate(this%idensR(NR))
        allocate(this%idensx(Nx))
        
        ! Zero out arrays ...
        this%psi = (0.0_dp, 0.0_dp)
        this%initialized = .true.
    end subroutine initialize

    !> Clean up and array deallocation
    subroutine deallocate_all(this)
        class(time_prop_2d), intent(inout):: this

        print*
        print*, "Cleaning up 2d time propagation variables ..."
        if(allocated(this%psi)) deallocate(this%psi)
        if(allocated(this%psi0)) deallocate(this%psi0)
        if(allocated(this%psi_end)) deallocate(this%psi_end)
        if(allocated(this%abs_R)) deallocate(this%abs_R)
        if(allocated(this%abs_x)) deallocate(this%abs_x)
        if(allocated(this%psi_out_R)) deallocate(this%psi_out_R)
        if(allocated(this%psi_out_x)) deallocate(this%psi_out_x)
        if(allocated(this%idensR)) deallocate(this%idensR)
        if(allocated(this%idensx)) deallocate(this%idensx)
        if(allocated(this%idenspR)) deallocate(this%idenspR)
        if(allocated(this%idenspx)) deallocate(this%idenspx)

        print*, "Done."
    end subroutine deallocate_all

    subroutine ini_dist_choice(this)
        use global_vars, only: NR, v_ini, N_ini, Ri_tdse, kappa_tdse, &
            & initial_distribution, R, x, chi0, ewf
        use data_au, only: au2a, au2eV
        class(time_prop_2d), intent(inout) :: this
        character(len=5):: divider
        integer :: i, j 
        real(dp):: norm

        print*
        print*, "Initial wavefunction:"
        this%psi = (0._dp, 0._dp)
        select case(initial_distribution) 
        case("single vibrational state")
            print*, "initial wavefunction is in..."
            print*, N_ini-1, "electronic state and in", v_ini-1, "vibrational state"
            do j = 1, Nx
                this%psi(:,j) = ewf(:,j,N_ini) * chi0(:,v_ini,N_ini)  !   & *exp(im*dpR*(r-ri))
            end do  
        case("gaussian distribution")
            print*, "initial wavefunction is in..."
            print*, N_ini-1, "electronic state and with a Gaussian distribution centered around",&
                & Ri_tdse, "a.u. \n with deviation of", kappa_tdse, "."
            do j = 1, Nx
                this%psi(:,j) = ewf(:,j,N_ini) * exp(kappa_tdse*(R(:)-Ri_tdse)**2) 
            enddo
          
        case default
            print*, "Default case: initial wavefunction is in..."
            print*, N_ini-1, "electronic state and with a Gaussian distribution centered around",&
                & Ri_tdse, "a.u. \n with deviation of", kappa_tdse, "."
            do j = 1, Nx
                this%psi(:,j) = ewf(:,j,N_ini) * exp(kappa_tdse*(R(:)-Ri_tdse)**2) 
            enddo
        end select

        ! Normalization of the initial wavefunction
        call integ_2d(this%psi, norm)
        print*,'norm1 =', sngl(norm)
        this%psi = this%psi / sqrt(norm)

        ! Writing header for initial wavefunction file
        write(this%psi_2d_tk,'(a,a,a,a,a)') "# R-grid(a.u.) ", divider, &
            & "x-grid(a.u.)", divider, &
            & "Electronic Ground state density "
        ! Writing initial wavefunction to file
        do i = 1, NR
            do j = 1, Nx
                write(this%psi_2d_tk,*) R(i), x(j), abs(this%psi(i,j))**2
            end do
            write(this%psi_2d_tk, *)
        end do
        close(this%psi_2d_tk)

        print*, "Wavefunction initialized."
    end subroutine ini_dist_choice

    !> Sets up absorber function for boundary treatment
    subroutine absorber_gen(this)
        use global_vars, only: NR, R, Nx, x, absorber
        use pot_param, only: cpmR, cpmx
        use varprecision, only: dp
        class(time_prop_2d), intent(inout) :: this
        character(len=5):: divider
        real(dp), allocatable:: cof(:), V_abs(:)
        complex(dp), allocatable:: exp_abs(:)
        integer:: i
        real(dp):: c
  
        allocate(this%abs_R(NR), this%abs_x(Nx))

        print*
        print*, "Absorber on R-grid is placed around the number of grid points from the end of the grid:"
        this%i_cpmR = minloc(abs(R(:) - cpmR), 1) - 50
        print*, "i_cpmR = ", this%i_cpmR, ", NR-i_cpmR", NR - this%i_cpmR
        print*, "R(NR-i_cpmR) = ", R(NR - this%i_cpmR)

        print*
        print*, "Absorber on x-grid is placed around the number of grid points from the both ends of the grid:"
        this%i_cpmx = minloc(abs(x(:) - cpmx), 1) - 50
        print*, "i_cpmx = ", this%i_cpmx, ", Nx-i_cpmx", Nx - this%i_cpmx
        print*, "x(Nx-i_cpmx) = ", x(Nx - this%i_cpmx)

        select case(absorber)
        case ("CAP")
            allocate(V_abs(NR), exp_abs(NR))
            call Complex_absorber_function(R, NR, cpmR, V_abs, exp_abs)
            this%abs_R = exp_abs
            deallocate(V_abs, exp_abs)

            allocate(V_abs(Nx), exp_abs(Nx))
            call Complex_absorber_function(x, Nx, cpmx, V_abs, exp_abs)
            this%abs_x = exp_abs
            deallocate(V_abs, exp_abs)
        case("mask")
            allocate(cof(NR))
            c = 1._dp
            call mask_function_ex(R, NR, cpmR, c, cof, sides=1)
            this%abs_R = cof
            deallocate(cof)
            
            allocate(cof(Nx))
            c = 1._dp
            call mask_function_ex(x, Nx, cpmx, c, cof, sides=2)
            this%abs_x = cof
            deallocate(cof)
        case default
            print*, "No absorber selected. Reflections off the grid boundary may occur!"
        end select

        ! writing hearder to the absorber function file
        write(this%abs_R_tk, '(a, a, a)') "# R-grid(a.u.) ", divider, &
            & "Absorber function(arb. units)" 
        ! Writing absorber function to file
        do i = 1, NR
            write(this%abs_R_tk,*) R(i), abs(this%abs_R(i))
        end do
        close(this%abs_R_tk)

        ! writing hearder to the absorber function file
        write(this%abs_x_tk, '(a, a, a)') "# x-grid(a.u.) ", divider, &
            & "Absorber function(arb. units)" 
        ! Writing absorber function to file
        do i = 1, Nx
            write(this%abs_x_tk,*) x(i), abs(this%abs_x(i))
        end do
        close(this%abs_x_tk)

        print*, "Done setting up the absorber."

    end subroutine absorber_gen

    subroutine open_files_to_write(this)
        use global_vars, only: time_prop_dir_2d
        class(time_prop_2d), intent(inout) :: this
        character(150) :: filepath

        ! inintial wavefunction
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "psi0_2d.out"
        open(newunit=this%psi_2d_tk,file=filepath,status='unknown') 
        ! Absorber function
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "absorber_function_x-grid.out"
        open(newunit=this%abs_x_tk,file=filepath,status='unknown') 
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "absorber_function_R-grid.out"
        open(newunit=this%abs_R_tk,file=filepath,status='unknown') 
        
        ! Propagation outputs
        ! time dependent norm
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "norm_2d.out"
        open(newunit=this%norm_2d_tk,file=filepath,status='unknown')
        
        ! time dependent R density 
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "td-density_R.out"
        open(newunit=this%dens_R_tk,file=filepath,status='unknown')
        
        ! time dependent x density
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "td-density_x.out"
        open(newunit=this%dens_x_tk,file=filepath,status='unknown')
        
        ! time dependent momentum density 
        !write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "momt_density_1d_pm3d.out"
        !open(newunit=this%Pdens_1d_tk,file=filepath,status='unknown')
        
        ! time dependent average position 
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "avgR_2d.out"
        open(newunit=this%avgR_2d_tk,file=filepath,status='unknown') 
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "avgx_2d.out"
        open(newunit=this%avgx_2d_tk,file=filepath,status='unknown') 
        
        ! time dependent norm in localized states
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "norm_pn_2d.out"
        open(newunit=this%norm_pn_2d_tk,file=filepath,status='unknown')
        
        ! time dependent electric field 
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "field_2d.out"
        open(newunit=this%field_2d_tk,file=filepath,status='unknown')
        
        ! time dependent absorbed momentum
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "momentum_2d.out"
        open(newunit=this%momt_2d_tk,file=filepath,status='unknown')
        
        ! time dependent vibrational populations 
        !write(filepath,'(a,a)') adjustl(trim(time_prop_dir_2d)), 'vibpop1D_lambda.out'
        !open(newunit=this%vibpop_1d_tk,file=filepath,status='unknown')
        
        ! time dependent norm of absorbed wavepacket
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "psi_outR_norm_2d.out"
        open(newunit=this%psi_outR_norm_2d_tk,file=filepath,status='unknown')
        
        ! time dependent momentum density of absorbed wavepacket 
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_2d)), "psi_outR_momt_density_2d_pm3d.out"
        open(newunit=this%psi_outR_Pdens_2d_tk,file=filepath,status='unknown') 
    
    end subroutine open_files_to_write

    subroutine write_output_files(this)
        class(time_prop_2d), intent(inout) :: this
        ! Post-propagation output is currently handled inside propagation_2D_impl
    end subroutine write_output_files

    subroutine time_evolution(this, E, A, alpha_t, propagator_method)
        use global_vars, only: NR, time, Nt, dp, dR, R, omp_nthreads, pot, x, CalcMode
        use data_au, only: au2fs
        use omp_lib
        class(time_prop_2d), intent(inout) :: this
        type(split_operator_2d_type) :: split_operator_2d
        type(rk4_operator_2d_type)   :: rk4_operator_2d
        integer :: i, j, k
        integer :: max_num_threads
        real(dp) :: evR, evx, epx, epR
        real(dp) :: norm
        real(dp) :: E(Nt), A(Nt)
        real(dp) :: alpha_t(Nt)
        real(dp) :: E_half, A_half       ! fields at t + dt/2 for RK4
        real(dp) :: alpha_half           ! quiver disp at t + dt/2 for RK4-KH
        real(dp), allocatable :: pot_kh(:,:), pot_kh_half(:,:), pot_kh_next(:,:)
        real(dp) :: E_zero, A_zero       ! zero fields for KH mode
        character(*), intent(in) :: propagator_method
        character(20) :: in_xR, out_x, out_R, out_xR
        complex(dp), allocatable :: psi_out_R_tmp(:,:), psi_out_x_tmp(:,:)
        complex(dp), allocatable :: psi_out_xR_tmp(:,:)
        
        allocate(psi_out_R_tmp(NR,Nx), psi_out_x_tmp(NR,Nx))
        allocate(psi_out_xR_tmp(NR,Nx))

        ! Initialize propagator based on selected method
        select case(trim(adjustl(propagator_method)))
        case("split_operator")
            split_operator_2d%gauge = "length"
            call split_operator_2d%fft_initialize()
            call split_operator_2d%split_operator_initialize()
        case("rk4")
            rk4_operator_2d%gauge = "length"
            call rk4_operator_2d%fft_initialize()
        case default
            split_operator_2d%gauge = "length"
            call split_operator_2d%fft_initialize()
            call split_operator_2d%split_operator_initialize()
        end select

        ! Defining simulation regions
        write(in_xR, '(a)') "inner-xR"
        write(out_x, '(a)') "outer-x"
        write(out_R, '(a)') "outer-R"
        write(out_xR, '(a)') "outer-xR"

        print*
        print*,'2D time evolution...'
        print*

        ! Allocate KH potential arrays if using time-dependent KH mode
        if (trim(CalcMode) == "KH_td") then
            allocate(pot_kh(NR, Nx))
            if (trim(adjustl(propagator_method)) == "rk4") then
                allocate(pot_kh_half(NR, Nx))
                allocate(pot_kh_next(NR, Nx))
            end if
            E_zero = 0._dp
            A_zero = 0._dp
        end if

        timeloop: do k = 1, Nt
            if (mod(k,1000) .eq. 0 .and. time(k)*au2fs .lt. 100._dp) then
                print('(a,i0,a,f5.2,a)'), "Progress: time step #", k, ", time:", time(k)*au2fs, " fs"
            elseif (mod(k,1000) .eq. 0 .and. time(k)*au2fs .ge. 100._dp) then
                print('(a,i0,a,f6.2,a)'), "Progress: time step #", k, ", time:", time(k)*au2fs, " fs"
            endif

            evR = 0.d0
            evx = 0.d0
            epR = 0.d0
            epx = 0.d0

            !====================================================================
            ! KH time-dependent mode: compute instantaneous KH potential
            !====================================================================
            if (trim(CalcMode) == "KH_td") then
                select case(trim(adjustl(propagator_method)))
                case("split_operator")
                    call build_kh_potential_at_time(pot_kh, alpha_t(k))
                    call split_operator_2d%vprop_gen_kh(pot_kh)
                    call split_operator_2d%split_operator_step(this%psi, in_xR)
                case("rk4")
                    call build_kh_potential_at_time(pot_kh, alpha_t(k))
                    if (k < Nt) then
                        alpha_half = 0.5_dp * (alpha_t(k) + alpha_t(k+1))
                    else
                        alpha_half = alpha_t(k)
                    end if
                    call build_kh_potential_at_time(pot_kh_half, alpha_half)
                    call build_kh_potential_at_time(pot_kh_next, alpha_t(k+1))
                    call rk4_operator_2d%rk4_step_kh(this%psi, dt, pot_kh, pot_kh_half, pot_kh_next)
                case default
                    call build_kh_potential_at_time(pot_kh, alpha_t(k))
                    call split_operator_2d%vprop_gen_kh(pot_kh)
                    call split_operator_2d%split_operator_step(this%psi, in_xR)
                end select
            else
                !================================================================
                ! Original Lab-frame propagation
                !================================================================
                select case(trim(adjustl(propagator_method)))
                !=============================================================
                case("split_operator")
                !=============================================================
                    call split_operator_2d%vprop_gen_len(E(k), A(k))
                    call split_operator_2d%split_operator_step(this%psi, in_xR)

                !=============================================================
                case("rk4")
                !=============================================================
                    ! Compute fields at halftime: E(t+dt/2), A(t+dt/2)
                    ! and next step: E(t+dt), A(t+dt) for the k4 evaluation
                    if (k < Nt) then
                        E_half = 0.5_dp * (E(k) + E(k+1))
                        A_half = 0.5_dp * (A(k) + A(k+1))
                        call rk4_operator_2d%rk4_step(this%psi, dt, E(k), E_half, E(k+1), &
                            & A(k), A_half, A(k+1), pot)
                    else
                        E_half = E(k)
                        A_half = A(k)
                        call rk4_operator_2d%rk4_step(this%psi, dt, E(k), E_half, E(k), &
                            & A(k), A_half, A(k), pot)
                    end if

                !=============================================================
                case default
                !=============================================================
                    call split_operator_2d%vprop_gen_len(E(k), A(k))
                    call split_operator_2d%split_operator_step(this%psi, in_xR)

                end select
            end if

            call density(this%psi, this%idensR, this%idensx)
            call integ_2D(this%psi, norm)
 
            evR = sum(dble(R(:) * this%idensR(:)))*dR
            evx = sum(dble(x(:) * this%idensx(:)))*dx
            evR = evR / norm
            evx = evx / norm

            ! write time dependent outputs to files
            write(this%avgR_2d_tk,*) time(k) * au2fs, evR !, sngl(epR)
            write(this%avgx_2d_tk,*) time(k) * au2fs, evx !, sngl(epx)
            write(this%norm_2d_tk,*) time(k) * au2fs, norm 
            write(this%field_2d_tk,*) time(k) * au2fs, E(k), A(K)
            
            if(mod(K,50).eq.0) then
                ! R density map  
                do i = 1, NR, 4
                    write(this%dens_R_tk,*) time(k) * au2fs, R(i), this%idensR(i)
                enddo
                write(this%dens_R_tk, *)

                ! x density map  
                do j = Nx/4, 3*Nx/4
                    write(this%dens_x_tk,*) time(k) *au2fs, x(j), this%idensx(j)
                enddo
                write(this%dens_x_tk, *)
            endif

            ! absorbed wavepacket
            do j = 1, Nx
                psi_out_xR_tmp(:,j) = this%psi_out_x(:,j) * (1.0d0 - this%abs_R(:)) ! dissociative ionization
                psi_out_x_tmp(:,j) = this%psi(:,j) * (1.d0 - this%abs_x(j))  ! ionization
                this%psi_out_x(:,j) = this%psi_out_x(:,j) * this%abs_R(:)   ! cutting off the dissociated wavefunction
    
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !%% localization and dissociation 
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                psi_out_R_tmp(:,j) = this%psi(:,j)*(1.0d0 - this%abs_R(j)) ! dissociating localized wavefunction

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !%% Total cut off %%%%%%%%%%%%%%%%%%%%%%
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                this%psi(:,j) = this%psi(:,j) * this%abs_x(j) * this%abs_R(:)   
            enddo

        end do timeloop

        deallocate(psi_out_R_tmp, psi_out_x_tmp, psi_out_xR_tmp)
        if (allocated(pot_kh)) deallocate(pot_kh)
        if (allocated(pot_kh_half)) deallocate(pot_kh_half)
        if (allocated(pot_kh_next)) deallocate(pot_kh_next)
        close(this%avgR_2d_tk)
        close(this%avgx_2d_tk)
        close(this%norm_2d_tk)
        close(this%dens_R_tk)
        close(this%dens_x_tk)
        close(this%field_2d_tk)

        ! Destroy FFTW plans and free memory using encapsulated finalize
        select case(trim(adjustl(propagator_method)))
        case("split_operator")
            call split_operator_2d%finalize()
        case("rk4")
            call rk4_operator_2d%finalize()
        case default
            call split_operator_2d%finalize()
        end select
    end subroutine time_evolution

    !_________________________________________________________

    Function cis(expo)
        use data_au, only: im
        implicit none
 
        complex*16:: cis
        double precision:: expo
 
        cis=cos(expo)+im*sin(expo)
        return
    end function

    ! .................................................................
    subroutine integ_2D(psi, norm)

        use global_vars, only: NR, Nx, dx, dR
        implicit none

        double precision,intent(out):: norm
        complex*16,intent(in):: psi(NR, Nx)

        norm = 0.d0

        norm = sum(abs(psi(:,:))**2)*dx*dR
 
        return
    end subroutine


    ! ..................................................................

    subroutine density(psi,idensR,idensx)

        use global_vars, only: NR, Nx, dx, dR

        implicit none
        integer:: I, J 
        double precision,intent(out)::idensx(Nx), idensR(NR) 
        complex*16,intent(in):: psi(NR,Nx)
 
        idensR = 0.d0
        idensx = 0.d0

        do I = 1, NR
            idensR(I) =sum(abs(psi(I,:))**2)*dx
        end do


        do J = 1, Nx
            idensx(J) =sum(abs(psi(:,J))**2)*dR
        end do

        return
    end subroutine

    !...................................................

    subroutine overlap_2d(psi1, psi2, C)

        use global_vars, only: NR, Nx, dx, dR
        implicit none
        integer:: I
        complex*16:: C, F(NR)
        complex*16,intent(in):: psi1(NR,Nx), psi2(NR,Nx)

        C = (0.d0, 0.d0); F = (0.d0, 0.d0)

        do I = 1, NR
            f(i) = sum(conjg(psi1(i,:)) * psi2(i,:))
        end do

        f = f * dx

        C = sum(f(:))*dR

        return
    end subroutine

    !........................................................................

    subroutine pop_analysis(psi, time, ewf, trans_exp_pop, b)

        use global_vars, only: NR, Nx, Nstates, dx, dR
        use data_au, only: au2fs

        implicit none

        integer:: I, N

        double precision,intent(in):: time, ewf(Nx, NR, Nstates)
        double precision:: B(Nstates)
        complex*16,intent(in):: psi(NR,Nx)
        complex*16:: a(NR,Nstates)
        complex*16::trans_exp_pop(NR,Nx)

        a = (0.d0,0.d0)
        b = 0.d0

        !$OMP parallel do
        do N = 1, Nstates
            !   !$OMP parallel do
            do I = 1, NR
                !    !$OMP parallel do
                !     do J = 1, Nx
                a(I,n) = sum(trans_exp_pop(I,:)*ewf(:,I,n) * psi(I,:))
                !     end do
                !    !$OMP end parallel do
            end do
            !   !$OMP end parallel do
        end do
        !$OMP end parallel do

        a = a * dx

        do N = 1, Nstates
            b(n) = sum(abs(a(:,n))**2)
        end do

        b = b * dR

        write(500,*) time *au2fs, b(1)
        write(501,*) time *au2fs, b(2) 
   
        return
    end subroutine

! .......................................................................

subroutine osc_dipole(psi, d_t1, d_t2,grad)

use global_vars, only: NR, Nx, dx, dR

 implicit none

 integer:: I, J
 double precision,intent(out):: d_t1, d_t2(NR)
 double precision,intent(in)::grad(nr,nx)
 complex*16,intent(in):: psi(NR,Nx)
 
 d_t1 = 0.d0
 d_t2 = 0.d0

 do i = 1, Nr
  do j = 1, Nx
   d_t2(i) = d_t2(i) + abs(psi(i,j))**2 * grad(i,j)    
  end do
   d_t2(i) = -d_t2(i) * dx
 end do 
 
 
 do i = 1, Nr
  d_t1 = d_t1 + d_t2(i)
 end do
 
  d_t1 = d_t1 * dR
  

return 
end subroutine

!........................................................................

    subroutine localization(psi, time, ewf, K, trans_exp_pop)
    
    use global_vars, only: NR, Nx, Nstates, dx, dR
    use data_au, only: au2fs
     implicit none
    
     integer:: I, K, N
    
     double precision,intent(in):: time
     double precision,intent(in):: ewf(Nx,NR,Nstates)
     double precision:: B(Nstates), pl_loc(Nx,nR), neg_loc(nx,Nr)
     complex(kind=kind(0.d0)),intent(in):: psi(NR,Nx)
     complex(kind=kind(0.d0)):: psi2(NR,Nx), a(nr,Nstates) 
     complex*16::trans_exp_pop(NR,Nx)
    
     psi2 = psi
    
     pl_loc = 0.d0
     neg_loc = 0.d0
    
     do i = 1, Nr
    !  do j = 1, nx
       pl_loc(:,i) = 1./sqrt(2.d0)* (ewf(:,i,1) + ewf(:,i,2))
       neg_loc(:,i) = 1./sqrt(2.d0)* (ewf(:,i,1) - ewf(:,i,2))
    !  end do
     end do
   
   
      a = (0.d0,0.d0)
      b = 0.d0
   
   

       do I = 1, NR
     !   do J = 1, Nx
         a(I,1) = sum(pl_loc(:,I) * psi2(I,:) * trans_exp_pop(I,:))
         a(I,2) = sum(neg_loc(:,I) * psi2(I,:) * trans_exp_pop(I,:))
     !   end do
       end do
     
     
      a = a * dx
     

     
      do N = 1, Nstates
    !   do I = 1, NR
         b(n) = sum(abs(a(:,n))**2)
    !   end do
      end do
    
      b = b * dR


     if(mod(K,100).eq.0) then  ! writing out the probability amplitudes
       do I = 1, NR
      
        write(504,*) time*au2fs, R(I), abs(a(i,1))**2
        write(505,*) time*au2fs, R(I), abs(a(i,2))**2
       end do
     
       write(504,*)
       write(505,*)
     end if
   
      write(502,*) sngl(time *au2fs), b(1)
      write(503,*) sngl(time *au2fs), b(2)

    return
    end subroutine


!............... Cut off Functions ................

    !> Generates complex absorber function for boundary
    subroutine complex_absorber_function(R, NR, cpmR, v_abs, f)
        use global_vars, only: dp, dt
        
        integer i, NR
        real(dp):: a, eps, V_abs(NR), n, R0, p
        real(dp):: R(NR), cpmR
        complex(dp):: f(NR)
    
        eps = epsilon(a) 
        print*, "Lower limit of the precision:", eps
        n = 4 ! power of absorber function
        R0 = R(NR)- cpmR ! start of the absorber
        p = 20._dp ! optimal absorption momentum
        a = -log(eps) *(n+1) *p / (2*(R(NR)-R0)**(n+1))
        print*, "Absorber prefactor a:", a

        do i = 1, NR
            if (R(i) .gt. abs(R0)) then
                V_abs(i) = a*(R(i)-R0)**n
            else
                V_abs(i) = 0._dp
            endif
        enddo
        f(:) = exp(-dt *V_abs(:))
    end subroutine

    !> Generates mask absorber function (exponential)
    subroutine mask_function_ex(R, NR, cpmR, c, cof, sides)
        use global_vars, only: dp
    
        integer :: i, NR, sides
        real(dp):: cof(NR),c, R(NR), cpmR
  
        select case (sides)
        case(1) ! left boundary
            do i = 1, NR
                cof(i)=1.0d0/(1.0d0+exp(c*(R(i)-R(NR)+cpmR)))
            end do
        case(2) ! Both boundaries
            do i = 1, NR
                cof(i)=1.0d0/(1.0d0+exp(c*(R(i)-R(NR)+cpmR)))
            end do
            do i = 1, NR/2
                cof(i) = cof(NR-i+1)
            end do
        end select

    end subroutine mask_function_ex

end module propagation2d_mod

