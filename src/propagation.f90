module propagation1d_mod
    use varprecision, only: dp
    use split_operator_mod, only: split_operator_type
    use rk4_operator_mod, only: rk4_operator_type
    use continuum_1d_mod, only: continuum_1d_type
    implicit none
    private
    public :: time_prop, propagation_1D

    !> Type to hold all data and methods for time propagation
    type :: time_prop
            ! Vibrational wavefunctions for all electronic states
            real(dp), allocatable :: chi0(:,:,:) 
            ! Total propagated wavefunction (all electronic states)
            complex(dp), allocatable :: psi_ges(:,:)
            ! Vibrational energies for all states
            real(dp), allocatable :: vib_en(:,:)
            ! Vibrational population amplitudes (projection)
            complex(dp), allocatable :: psi_chi(:)
            ! Index for absorber placement
            integer :: i_cpmR
            ! Absorber function (complex or mask)
            complex(dp), allocatable :: abs_func(:)

            ! Absorbed/dissociated part of wavefunction
            complex(dp), allocatable :: psi_outR(:,:)
            ! Incremental absorbed density (for analysis)
            real(dp), allocatable :: psi_outR_inc(:,:)
        
            ! File handles for output files
            integer:: psi_1d_tk, cof_1d_tk
            integer:: norm_1d_tk, dens_1d_tk, ex_dens_1d_tk, Pdens_1d_tk
            integer:: avgR_1d_tk, avgpR_1d_tk, norm_pn_1d_tk
            integer:: field_1d_tk, accumulation_1d_tk, momt_1d_tk
            integer:: vibpop_1d_tk, psi_outR_norm_1d_tk, psi_outR_Pdens_1d_tk
    contains
    procedure :: initialize
    procedure :: read_pot_files      ! Reads vibrational state and energy files
    procedure :: ini_dist_choice     ! Initializes the wavefunction distribution
    procedure :: time_evolution      ! Main time propagation loop
    procedure :: localized_states_norm, expected_position
    procedure :: vib_pop_analysis    ! Analyzes vibrational populations
    procedure :: propagation_1D      ! High-level wrapper for propagation
    procedure :: absorber_gen        ! Sets up absorber function
    procedure :: open_files_to_write, write_headers_to_files
    procedure :: continuum_prop      ! Handles continuum part of wavefunction
    procedure :: post_propagation_analysis => post_prop_analysis
    procedure :: deallocate_all      ! Cleans up and deallocates arrays
    end type time_prop

contains

    !> High-level wrapper for 1D time propagation
    subroutine propagation_1D(this, E, A, propagator_method)
        use global_vars, only: Nt
        class(time_prop), intent(inout) :: this
        real(dp), intent(in) :: E(Nt), A(Nt)
        character(*), intent(in) :: propagator_method
        print*, "Preparing time propagation ..."
        call this%initialize()
        call this%read_pot_files()
        call this%open_files_to_write()
        call this%write_headers_to_files()
        call this%ini_dist_choice()
        call this%absorber_gen()
        call this%time_evolution(E, A, propagator_method)
        call this%post_propagation_analysis()
        call this%deallocate_all()
    end subroutine propagation_1D

    !> Allocates arrays and initializes variables for propagation
    subroutine initialize(this)
        use global_vars, only: NR, Nstates, guess_vstates
        class(time_prop), intent(inout) :: this
        allocate(this%chi0(NR,guess_vstates,Nstates))
        allocate(this%psi_ges(NR,Nstates))
        allocate(this%psi_outR(NR,Nstates), this%psi_outR_inc(NR,Nstates))
        allocate(this%psi_chi(guess_vstates))
        allocate(this%vib_en(guess_vstates,Nstates))
    end subroutine initialize

    !> Reads vibrational state and energy files from disk
    subroutine read_pot_files(this)
        use global_vars, only: NR, Nstates, Vstates, output_data_dir, chi0
        use varprecision, only: dp, sp
        class(time_prop), intent(inout) :: this
        character(30):: nucl_wp_path
        character(150):: filepath
        integer:: chi0_tk, vstates_tk, vib_en_tk
        integer:: i, N, V, i_dummy
        real(sp):: dummy

        nucl_wp_path = "nuclear_wavepacket_data/"

        ! implementation for reading potential files
        write(filepath,'(a,a,a)') adjustl(trim(output_data_dir)), adjustl(trim(nucl_wp_path)), "Bound-vibstates_in_Nthstates.out"
        call file_status_check(filepath)
        open(newunit=vstates_tk,file=filepath,status='unknown')


        this%chi0 = 0._dp
        this%vib_en = 0._dp
        do N = 1, Nstates 
            print*
            print*, "Reading vibrational states in the Electronic state ", N-1
            print*

            read(vstates_tk,*) i_dummy, Vstates(N)

            write(filepath,'(a,a,a,i0,a)') adjustl(trim(output_data_dir)), adjustl(trim(nucl_wp_path)), &
                & "BO_Electronic-state-g", int(N-1), "_Evib.out"
            call file_status_check(filepath)
            open(newunit=vib_en_tk,file=filepath,status='unknown')
            do V = 1, Vstates(N)
                read(vib_en_tk,*) i_dummy, this%vib_en(V,N)
            enddo
            close(vib_en_tk)

            write(filepath,'(a,a,a,i0,a)') adjustl(trim(output_data_dir)), adjustl(trim(nucl_wp_path)), &
                & "BO_Electronic-state-g", int(N-1), "_vibstates.out"
            call file_status_check(filepath)
            open(newunit=chi0_tk,file=filepath,status='unknown')
            print*, "NR:", NR, "Vstates(N)", Vstates(N)

            do i = 1, NR
                read(chi0_tk,*) dummy, this%chi0(i,1:Vstates(N),N)
                !print*, i, dummy, this%chi0(i,Vstates(N),N)
            enddo 
            close(chi0_tk)
        enddo
        chi0 = this%chi0 ! updating globally
        print*, "Done reading the vibrational states files"
        close(vstates_tk)
    end subroutine read_pot_files

    !> Initializes the wavefunction distribution based on user choice
    subroutine ini_dist_choice(this)
        use global_vars, only: NR, Nstates, Vstates, v_ini, N_ini, Ri_tdse, kappa_tdse, &
            initial_distribution, R
        use data_au, only: au2a, au2eV
        class(time_prop), intent(inout) :: this
        character(len=5):: divider
        integer :: i, N, v
        real(dp):: norm(Nstates)
        real(dp), allocatable, dimension(:) :: vib_dist
        print*
        print*, "Initial wavefunction:"
        this%psi_ges = (0._dp, 0._dp)
        select case(initial_distribution) 
            case("single vibrational state")
                print*, "initial wavefunction is in..."
                print*, N_ini-1, "electronic state and in", v_ini-1, "vibrational state"
                do i = 1, NR
                    this%psi_ges(i,N_ini) = this%chi0(i,v_ini,N_ini) 
                enddo  

            case("gaussian distribution")
                print*, "initial wavefunction is in..."
                print*, N_ini-1, "electronic state and with a Gaussian distribution centered around",&
                    & Ri_tdse/au2a, "a.u. \n with deviation of", kappa_tdse, "."
                do i = 1, NR
                    this%psi_ges(i,1)=exp(kappa_tdse*(R(i)-Ri_tdse/au2a)**2) 
                enddo

            ! case ("input dist")
                !  allocate(v_dist_ini(N_ini))
                !  v_ini_check = 0
                !  do N = 1, N_ini
                !    write(filepath,'(a,a,i0,a)') adjustl(trim(output_data_dir)), "vib_dist_", int(N-1), ".out"
                !    open(newunit=vib_dist_tk,file=filepath,status='unknown')
                !    do 
                !       read(vib_dist_tk,*,iostat=io)
                !       if (io /= 0) exit
                !       v_ini_check = v_ini_check + 1
                !    enddo
                !    if (v_ini /= v_ini_check) then
                !       write(*,'(a,a,i0,a)') "Number of vibstates in 'vib_dist_", N, ".out' not equal to input" 
                !    endif
                !  enddo
                !  allocate(vib_ini(guess_vstates,N_ini), vib_dist(guess_vstates, N_ini))
                !  do N = 1, N_ini
                !     do v = 1, v_dist_ini(N)
                !         read(vib_dist_tk,*) vib_ini(v,N), vib_dist(v,N) 
                !     enddo
                !  enddo
                !  do i = 1, NR
                !    do v = 1, v_ini
                !       do N = 1, N_ini
                !         psi_ges(i,N) = psi_ges(i,N) + vib_dist(vib_ini(v,N),N) * chi0(i,vib_ini(v,N),N)
                !       enddo
                !    enddo
                !  enddo

            case ("Boltzmann dist")
                do N = 1, N_ini
                    call Boltzmann_distribution(N, this%vib_en(:,N)/au2eV, vib_dist)
                    do v = 1, Vstates(N)
                        do i = 1, NR
                            this%psi_ges(i,N) = this%psi_ges(i,N) + vib_dist(v) * this%chi0(i,v,N)
                        enddo
                    enddo
                enddo
                deallocate(vib_dist)
              
            case default
                do i = 1, NR
                    this%psi_ges(i,1) = this%chi0(i,1,1)
                enddo
              
        end select

        ! Normalization of the initial wavefunction
        call integ(this%psi_ges, norm)
        this%psi_ges(:,1) = this%psi_ges(:,1) / sqrt(norm(1))
        !this%psi_ges(:,2) = this%psi_ges(:,2) / sqrt(norm(2))

        ! Writing header for initial wavefunction file
        write(this%psi_1d_tk,'(a,a,a,a,a)') "R-grid(a.u.) ", divider, &
            & "Electronic Ground state density ", divider,&
            & "Electronic First excited state density"
        ! Writing initial wavefunction to file
        do i = 1, NR
            write(this%psi_1d_tk,*) R(i), abs(this%psi_ges(i,1))**2, abs(this%psi_ges(i,2))**2
        end do
        close(this%psi_1d_tk)

        print*, "Wavefunction initialized."
    end subroutine ini_dist_choice

    !> Sets up absorber function for boundary treatment
    subroutine absorber_gen(this)
        use global_vars, only: NR, R, absorber
        use pot_param, only: cpmR
        use varprecision, only: dp
        class(time_prop), intent(inout) :: this
        character(len=5):: divider
        real(dp), allocatable:: cof(:), V_abs(:)
        complex(dp), allocatable:: exp_abs(:)
        integer:: i
  
        allocate(this%abs_func(NR))

        print*
        print*, "Absorber placed around the number of grid points from the end of the grid:"
        this%i_cpmR = minloc(abs(R(:) - cpmR), 1) - 50
        print*, "i_cpmR =", this%i_cpmR, ", NR-i_cpmR", NR - this%i_cpmR
        print*, "R(NR-i_cpmR) =", R(NR - this%i_cpmR)

        select case(absorber)
            case ("CAP")
                allocate(V_abs(NR), exp_abs(NR))
                call Complex_absorber_function(V_abs, exp_abs)
                this%abs_func = exp_abs
            case("mask")
                allocate(cof(NR))
                call mask_function_ex(cof)
                this%abs_func = cof
            case default
                print*, "No absorber selected. Reflections off the grid boundary may occur!"
        end select

        ! writing hearder to the absorber function file
        write(this%cof_1d_tk, '(a, a, a)') "R-grid(a.u.) ", divider, &
            & "Absorber function(arb. units)" 
        ! Writing absorber function to file
        do i = 1, NR
            write(this%cof_1d_tk,*) R(i), abs(this%abs_func(i))
        end do
        close(this%cof_1d_tk)

        print*, "Done setting up the absorber."

    end subroutine absorber_gen

    !> Opens all output files for writing propagation results
    subroutine open_files_to_write(this)
        use global_vars, only: time_prop_dir_1d
        class(time_prop), intent(inout) :: this
        character(150):: filepath

        ! inintial wavefunction
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "psi0_1d.out"
        open(newunit=this%psi_1d_tk,file=filepath,status='unknown') 
        ! Absorber function
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "absorber_function.out"
        open(newunit=this%cof_1d_tk,file=filepath,status='unknown') 
        
        ! Propagation outputs
        ! time dependent norm
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "norm_1d.out"
        open(newunit=this%norm_1d_tk,file=filepath,status='unknown')
        
        ! time dependent ground state density 
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "density_1d_pm3d.out"
        open(newunit=this%dens_1d_tk,file=filepath,status='unknown')
        
        ! time dependent excited state density
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "ex_density_1d_pm3d.out"
        open(newunit=this%ex_dens_1d_tk,file=filepath,status='unknown')
        
        ! time dependent momentum density 
        !write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "momt_density_1d_pm3d.out"
        !open(newunit=this%Pdens_1d_tk,file=filepath,status='unknown')
        
        ! time dependent average position 
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "avgR_1d.out"
        open(newunit=this%avgR_1d_tk,file=filepath,status='unknown') 
        
        ! time dependent norm in localized states
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "norm_pn_1d.out"
        open(newunit=this%norm_pn_1d_tk,file=filepath,status='unknown')
        
        ! time dependent electric field 
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "field_1d.out"
        open(newunit=this%field_1d_tk,file=filepath,status='unknown')
        
        ! time dependent absorbed momentum
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "momentum_1d.out"
        open(newunit=this%momt_1d_tk,file=filepath,status='unknown')
        
        ! time dependent vibrational populations 
        write(filepath,'(a,a)') adjustl(trim(time_prop_dir_1d)), 'vibpop1D_lambda.out'
        open(newunit=this%vibpop_1d_tk,file=filepath,status='unknown')
        
        ! time dependent norm of absorbed wavepacket
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "psi_outR_norm_1d.out"
        open(newunit=this%psi_outR_norm_1d_tk,file=filepath,status='unknown')
        
        ! time dependent momentum density of absorbed wavepacket 
        write(filepath, '(a,a)') adjustl(trim(time_prop_dir_1d)), "psi_outR_momt_density_1d_pm3d.out"
        open(newunit=this%psi_outR_Pdens_1d_tk,file=filepath,status='unknown') 

    end subroutine open_files_to_write

    !> Writes headers to output files for easier post-processing
    subroutine write_headers_to_files(this)
        class(time_prop), intent(inout) :: this
        character(len=5):: divider

        ! time dependent norm
        write(this%norm_1d_tk,'(a,a,a)') "# Time(fs)", divider, &
            & "Norm(∫|psi|^2)->[Nstates Columns]"

        ! time dependent ground state density 
        write(this%dens_1d_tk,'(a,a,a,a,a)') "#(gnuplot pm3d format) Time(fs)", divider, "R(a.u.)", &
            & divider, "Density |psi|^2"
        
        ! time dependent excited state density
        write(this%ex_dens_1d_tk,'(a,a,a,a,a)') "#(gnuplot pm3d format) Time(fs)", divider, "R(a.u.)", &
            & divider, "Density |psi|^2 -> 1:Nstates Columns" 

        ! time dependent average position 
        write(this%avgR_1d_tk,'(a,a,a)') "# Time(fs)", divider, &
            & "Expectation value of R(a.u.)->[Nstates Columns]"
        
        ! time dependent norm in localized states
        write(this%norm_pn_1d_tk,'(a,a,a,a,a)') "# Time(fs)", divider, &
            & "Norm_p(∫|psi_+|^2) ", divider, "Norm_n(∫|psi_-|^2)"
        
        ! time dependent electric field 
        write(this%field_1d_tk,'(a,a,a,a,a)') "# Time(fs)", divider, &
            & "Electric Field(a.u.)", divider, "Vector field(a.u.)"
        
        ! time dependent vibrational populations (only on electonic ground state)
        write(this%vibpop_1d_tk,'(a,a,a,a)') "# Time(fs)", divider, &
            & "Vibrational Poulation density on the Electronic ground state", &
            & "->[Vstates Columns]"
        
        ! time dependent norm of absorbed wavepacket
        write(this%psi_outR_norm_1d_tk,'(a,a,a)') "# Time(fs)", divider, &
            & "Absorbed Norm(∫|psi_absorbed|^2)"

        ! time dependent momentum density of absorbed wavepacket 
        write(this%psi_outR_Pdens_1d_tk,'(a,a,a,a,a)') "# Time(fs)", divider, & 
            & "Momentum(a.u.)", divider, "Density |psi|^2 -> [Nstates Columns]"

    end subroutine write_headers_to_files

    !> Main time propagation loop: evolves wavefunction and writes outputs
    !! Supports both "split_operator" and "rk4" propagator methods
    subroutine time_evolution(this, E, A, propagator_method)
        use global_vars, only: NR, Nstates, time, mu_all, Nt, &
            & dp, dR, guess_vstates, dt, Vstates, R, pR, omp_nthreads, adb
        use data_au, only: au2fs
        use blas_interfaces_module, only: zgemv
        use FFTW3
        use omp_lib
        class(time_prop), intent(inout) :: this
        type(split_operator_type) :: split_operator
        type(rk4_operator_type)   :: rk4_operator
        type(continuum_1d_type)   :: continuum_1d
        integer :: i, j, k, v
        integer :: max_num_threads
        real(dp), allocatable, dimension(:) :: evr, momt
        real(dp), allocatable, dimension(:) :: norm, normPn, norm_outR
        real(dp), allocatable, dimension(:) :: norm_SE_outR
        real(dp) :: E(Nt), A(Nt)
        real(dp) :: E_half       ! E(t + dt/2) for RK4
        character(*), intent(in) :: propagator_method
        complex(dp), allocatable :: tout(:,:)
        complex(dp), allocatable :: psi_Nstates(:), psi_Nstates1(:)
        real(dp), allocatable :: psi_Nstates_real(:), psi_Nstates_imag(:)
        real(dp), allocatable :: vib_pop(:,:)

        allocate(psi_Nstates(Nstates), psi_Nstates1(Nstates))
        allocate(psi_Nstates_real(Nstates), psi_Nstates_imag(Nstates))
        allocate(evr(Nstates), momt(Nstates))
        allocate(norm(Nstates), normPn(Nstates), norm_outR(Nstates))
        allocate(norm_SE_outR(Nstates))
        allocate(tout(Nstates,Nstates))
        allocate(vib_pop(guess_vstates,Nstates))

        ! Initialize continuum/absorber (always needed regardless of propagator)
        call continuum_1d%initialize()

        ! Initialize propagator based on selected method
        select case(trim(adjustl(propagator_method)))
        case("split_operator")
            print*, "Using split-operator propagator ..."
            call split_operator%fft_initialize()
            call split_operator%kprop_gen()
            call split_operator%vprop_gen()
        case("rk4")
            print*, "Using RK4 propagator ..."
            call rk4_operator%fft_initialize()
        case default
            print*, "Warning: Unknown propagator_method '", trim(adjustl(propagator_method)), &
                & "'. Defaulting to split_operator."
            call split_operator%fft_initialize()
            call split_operator%kprop_gen()
            call split_operator%vprop_gen()
        end select

        max_num_threads = omp_get_max_threads()
        print*
        print*, "Setting openmp threads for matrix operations ..."
        if (omp_nthreads > 0 .and. omp_nthreads < max_num_threads) then
            continue
        else
            omp_nthreads = max_num_threads
        endif       
        print*, "Number of openmp threads (if not speciefied then maximum available):", omp_nthreads

        call OMP_SET_NUM_THREADS(omp_nthreads) 
        print*, "Done  setting up Openmp threads."

        print*
        print*,'1D time evolution...'
        print*
    
        this%psi_outR = (0._dp,0._dp)
        this%psi_outR_inc = 0._dp
        
        timeloop: do k = 1, Nt
            
            if (mod(k,1000) .eq. 0 .and. time(k)*au2fs .lt. 100._dp) then
                print('(a,i0,a,f5.2,a)'), "Progress: time step #", k, ", time:", time(k)*au2fs, " fs"
            elseif (mod(k,1000) .eq. 0 .and. time(k)*au2fs .ge. 100._dp) then
                print('(a,i0,a,f6.2,a)'), "Progress: time step #", k, ", time:", time(k)*au2fs, " fs"
            endif
            
            evr = 0._dp
            momt=0._dp
            norm = 0._dp
            normPn = 0._dp
            norm_outR = 0._dp

            select case(trim(adjustl(propagator_method)))
            !=============================================================
            case("split_operator")
            !=============================================================
                ! Kinetic propagation half step
                call split_operator%split_operator_step(this%psi_ges)
                ! Potential propagation full step
                ! Diagonal potential matrix half step
                do j = 1, Nstates         
                    this%psi_ges(:,j) = this%psi_ges(:,j) * split_operator%vprop(:,j)
                end do
                
                ! Off-diagonal potential matrix full step i.e interstates coupling
                !$OMP PARALLEL DO DEFAULT(NONE) FiRSTPRiVATE(tout, psi_Nstates, psi_Nstates1) &
                !$OMP FiRSTPRiVATE(psi_Nstates_real, psi_Nstates_imag) &
                !$OMP SHARED(mu_all, E, this, k, Nstates, NR, dt)
                do i = 1, NR
                    call pulse2(tout, mu_all(:,:,i), E(k))   
                    psi_Nstates(:) = this%psi_ges(i,:)
                    call zgemv('N', int(Nstates), int(Nstates), (1._dp, 0._dp),  &
                        & tout, size(tout,dim=1), psi_Nstates, 1, (0._dp,0._dp), &
                        & psi_Nstates1, 1)
                    this%psi_ges(i,:) = psi_Nstates1(:)
                end do
                !$OMP END PARALLEL DO

                ! Diagonal potential matrix half step
                do j = 1, Nstates         
                    this%psi_ges(:,j) = this%psi_ges(:,j) * split_operator%vprop(:,j)
                end do

                ! Kinetic propagation half step
                call split_operator%split_operator_step(this%psi_ges)

            !=============================================================
            case("rk4")
            !=============================================================
                ! Compute E_field at halftime: E(t + dt/2)
                ! For simplicity, linearly interpolate between E(k) and E(k+1)
                if (k < Nt) then
                    E_half = 0.5_dp * (E(k) + E(k+1))
                else
                    E_half = E(k)  ! last step, use E(k)
                end if
                call rk4_operator%rk4_step(this%psi_ges, dt, E(k), E_half, mu_all, adb)

            end select

            ! On the fly analysis and writing to files
            ! external field
            write(this%field_1d_tk,*) time(k) * au2fs, E(k), A(k)
            ! norm
            call integ(this%psi_ges, norm)
            write(this%norm_1d_tk,*) time(k) * au2fs, norm
            ! localized states norm
            call this%localized_states_norm(normPn)
            write(this%norm_pn_1d_tk,*) time(k) * au2fs, normPn
            ! expected position
            call this%expected_position(norm, evr)
            write(this%avgR_1d_tk,*) time(k) * au2fs, evr(1:Nstates)
            ! vibrational populations (writing only for the ground state)
            call this%vib_pop_analysis(1)
            do v = 1, Vstates(1)
                vib_pop(v,1) = abs(this%psi_chi(v))**2
            enddo
            write(this%vibpop_1d_tk,*) time(k) *au2fs, vib_pop(1:Vstates(1),1)

            ! density
            ! Coordinate space density 
            if(mod(k,100).eq.0) then
                do i = 1, NR, 4   
                    write(this%dens_1d_tk,*) time(k) *au2fs, R(i), abs(this%psi_ges(i,1))**2
                    write(this%ex_dens_1d_tk,*) time(k) *au2fs, R(i), abs(this%psi_ges(i,2:Nstates))**2
                end do 
                write(this%dens_1d_tk,*)
                write(this%ex_dens_1d_tk,*)
            end if

            !-----------------------------------
            ! Cotinuum part treatment
            call this%continuum_prop(continuum_1d)
            
            ! Absorbed density norm
            call integ(this%psi_outR, norm_outR)
            do J = 1, Nstates
                norm_SE_outR(J) = sum(this%psi_outR_inc(:,J))*dR
            enddo
            write(this%psi_outR_norm_1d_tk,*) time(k), norm_outR

            ! Absorbed momentum density
            ! Momentum space density
            if (mod(K,100) .eq. 0) then
                do i = NR/2 +1, NR, 4
                    write(this%psi_outR_Pdens_1d_tk,*) time(k) *au2fs, & 
                        & pR(i), abs(this%psi_outR(i,:))**2
                end do
                do i = 1, NR/2, 4
                    write(this%psi_outR_Pdens_1d_tk,*) time(k) *au2fs, &
                        & pR(i), abs(this%psi_outR(i,:))**2
                end do 
                write(this%psi_outR_Pdens_1d_tk,*)
            endif
        
        end do timeloop

        deallocate(psi_Nstates, psi_Nstates1)
        deallocate(psi_Nstates_real, psi_Nstates_imag)
        deallocate(evr, momt)
        deallocate(norm, normPn, norm_outR)
        deallocate(norm_SE_outR)
        deallocate(tout)
        deallocate(vib_pop)
        ! Destroy FFTW plans and free memory
        ! Continuum/absorber (always)
        call fftw_destroy_plan(continuum_1d%planF)
        call fftw_destroy_plan(continuum_1d%planB)
        call fftw_free(continuum_1d%p_in)
        call fftw_free(continuum_1d%p_out)

        ! Propagator-specific cleanup
        select case(trim(adjustl(propagator_method)))
        case("split_operator")
            call fftw_destroy_plan(split_operator%planF)
            call fftw_destroy_plan(split_operator%planB)
            call fftw_free(split_operator%p_in)
            call fftw_free(split_operator%p_out)
        case("rk4")
            call fftw_destroy_plan(rk4_operator%planF)
            call fftw_destroy_plan(rk4_operator%planB)
            call fftw_free(rk4_operator%p_in)
            call fftw_free(rk4_operator%p_out)
            deallocate(rk4_operator%kin_energy)
        case default
            call fftw_destroy_plan(split_operator%planF)
            call fftw_destroy_plan(split_operator%planB)
            call fftw_free(split_operator%p_in)
            call fftw_free(split_operator%p_out)
        end select
        ! Close all files
        close(this%norm_1d_tk)
        close(this%dens_1d_tk)
        close(this%ex_dens_1d_tk)
        !close(this%Pdens_1d_tk)
        close(this%avgR_1d_tk)
        close(this%norm_pn_1d_tk)
        close(this%field_1d_tk)
        close(this%momt_1d_tk)
        close(this%vibpop_1d_tk)
        close(this%psi_outR_norm_1d_tk)
        close(this%psi_outR_Pdens_1d_tk)
    end subroutine time_evolution

    !> Analyzes norm in localized states
    subroutine localized_states_norm(this, normPn)
        use global_vars, only: NR, Nstates
        class(time_prop), intent(inout) :: this
        real(dp):: normPn(Nstates)
        complex(dp) :: psi_loc(NR,Nstates)
        normPn = 0._dp

        psi_loc(:,1) = 1._dp/sqrt(2._dp)*(this%psi_ges(:,1) + this%psi_ges(:,2))
        psi_loc(:,2) = 1._dp/sqrt(2._dp)*(this%psi_ges(:,1) - this%psi_ges(:,2))

        call integ(psi_loc, normPn)
        
    end subroutine localized_states_norm

    !> Calculates expected position for each state
    subroutine expected_position(this, norm, evr)
        use global_vars, only: Nstates, dR, R
        class(time_prop), intent(inout) :: this
        real(dp):: evR(Nstates), norm(Nstates)
        integer:: N

        evr = 0._dp
        do N = 1, Nstates
            evR(N) = sum(abs(this%psi_ges(:,N))**2 * R(:))
        enddo
        evR = evR * dR

        do N = 1,Nstates
            if (norm(N).ge.1.e-8_dp) then
                evR(N) = evR(N)/norm(N)
            endif
        enddo

    end subroutine expected_position

    !> Analyzes vibrational populations for a given electronic state
    subroutine vib_pop_analysis(this, num_state)
        use global_vars, only: Vstates, dR
        class(time_prop), intent(inout) :: this
        integer:: L
        integer, intent(in) :: num_state ! electronic state for which vibpop is calculated

        this%psi_chi = 0._dp
        do L=1,vstates(num_state)
            this%psi_chi(L) = sum(this%chi0(:,L,num_state) * (this%psi_ges(:,num_state)))
            this%psi_chi(L) = this%psi_chi(L)*dR
        enddo

    end subroutine vib_pop_analysis

    !> Continuum propagation of wavefunction (absorbed/dissociated)
    subroutine continuum_prop(this, continuum_1d)
        class(time_prop), intent(inout) :: this
        class(continuum_1d_type), intent(inout) :: continuum_1d

        call continuum_1d%apply(this%psi_ges, this%psi_outR, this%psi_outR_inc, &
            & this%abs_func, this%i_cpmR)

    end subroutine continuum_prop

    !> Post-propagation analysis: calculates spectra and populations
    subroutine post_prop_analysis(this)
        use global_vars, only: NR, Nstates, dR, vstates, pR, m_red, time_prop_dir_1d
        use FFTW3
        class(time_prop), intent(inout) :: this
        type(continuum_1d_type) :: continuum_1d
        integer:: N, i, j
        character(150):: filepath
        integer:: KER_spectra_tk, KER_spectra_un_tk
        integer:: momt_spectra_tk, momt_spectra_un_tk

        real(dp), allocatable :: norm_bound(:), norm_diss(:), norm_outR(:)
        complex(dp), allocatable :: psi_bound(:,:), psi_diss(:,:)


        allocate(norm_bound(Nstates))
        allocate(norm_diss(Nstates))
        allocate(norm_outR(Nstates))
        allocate(psi_bound(NR,Nstates))
        allocate(psi_diss(NR,Nstates))
        print*
        print*, "Post-propagation analysis..."

        call continuum_1d%initialize()
        ! implementation of post-propagation analysis
        !Final vibrational population in ground state
        do N =1, Nstates
            print*
            print*, "Final vibrational population in the elec. state", int(N-1)
            print*, "Number of vibrational states:", vstates(N)
            psi_bound = 0._dp
            call this%vib_pop_analysis(N)
            
            do j = 1, vstates(N)
                do i = 1, NR
                    psi_bound(i,N) = psi_bound(i,N) + this%psi_chi(J) * this%chi0(i,J,N)
                enddo
                print*, 'Vibpop (',J-1,') =', abs(this%psi_chi(J))**2
            enddo
            norm_bound(N) = sum(abs(psi_bound(:,N))**2) * dR
            print*, 'Total population in state', int(N-1), ":", sum(abs(this%psi_ges(:,N))**2) * dR
            print*, 'Bound population in state', int(N-1), ":", sum(abs(psi_bound(:,N))**2) * dR
            psi_diss(:,N) = this%psi_ges(:,N) - psi_bound(:,N)
            print*, 'Unbound population in state', int(N-1), ":", sum(abs(psi_diss(:,N))**2) * dR
                !& sum(abs(psi_ges(:,N)-psi_bound(:,N))**2)*dR
            print*, 'Unbound population reached to the absorber in state', int(N-1), ":", &
                & sum(abs(this%psi_outR(:,N))**2) * dR
            print*, 'Calculating KER spectra in state ', int(N-1)

            call continuum_1d%forward_fft(psi_diss(:,N), psi_diss(:,N))

            norm_diss(N) = sum(abs(psi_diss(:,N))**2)*dR
            norm_outR(N) = sum(abs(this%psi_outR(:,N))**2)*dR

            print*, "Writing KER spectra in state ", int(N-1)
            ! file for unnormalized KER spectra  
            write(filepath,'(a,a,i0,a)') adjustl(trim(time_prop_dir_1d)), "KER_spectra_from_state_g", &
                &  int(N-1), "_unnormalized.out"
            open(newunit=KER_spectra_un_tk,file=filepath,status='unknown')
            ! file for unnormalized momentum spectra
            write(filepath,'(a,a,i0,a)') adjustl(trim(time_prop_dir_1d)), "momt_spectra_from_state_g",&
                &  int(N-1), "_unnormalized.out"
            open(newunit=momt_spectra_un_tk,file=filepath,status='unknown')
            ! file for normalized KER spectra
            write(filepath,'(a,a,i0,a)') adjustl(trim(time_prop_dir_1d)), "KER_spectra_from_state_g", &
                &  int(N-1), ".out"
            open(newunit=KER_spectra_tk,file=filepath,status='unknown')
            ! file for normalized momentum spectra
            write(filepath,'(a,a,i0,a)') adjustl(trim(time_prop_dir_1d)), "momt_spectra_from_state_g",&
                &  int(N-1), ".out"
            open(newunit=momt_spectra_tk,file=filepath,status='unknown')

            do i=NR/2 +1,NR
                ! Writing unnormalized momentum spectra
                write(momt_spectra_un_tk,*) pR(i), abs(psi_diss(i,N))**2, abs(this%psi_outR(i,N))**2
                ! Writing normalized momentum spectra
                write(momt_spectra_tk,*) pR(i), abs(psi_diss(i,N)/sqrt(norm_diss(N)))**2, &
                    & abs(this%psi_outR(i,N)/sqrt(norm_diss(N)))**2
            enddo
            do i=1,NR/2
                ! Writing unnormalized momentum spectra
                write(momt_spectra_un_tk,*) pR(i), abs(psi_diss(i,N))**2, abs(this%psi_outR(i,N))**2
                ! Writing normalized momentum spectra
                write(momt_spectra_tk,*) pR(i), abs(psi_diss(i,N)/sqrt(norm_diss(N)))**2, &
                    & abs(this%psi_outR(i,N)/sqrt(norm_outR(N)))**2
                ! Writing unnormalized KER spectra
                write(KER_spectra_un_tk,*) pR(i)**2/(2*m_red), m_red*abs(psi_diss(i,N))**2 / pR(i),&
                    & m_red*abs(this%psi_outR(i,N))**2 / pR(i), m_red*abs(this%psi_outR_inc(i,N))/pR(i) 
                ! Writing normalized KER spectra
                write(KER_spectra_tk,*) pR(i)**2/(2*m_red), m_red*abs(psi_diss(i,N)/sqrt(norm_diss(N)))**2 / pR(i), &
                    & m_red*abs(this%psi_outR(i,N)/sqrt(norm_outR(N)))**2 / pR(i) 
            enddo
            close(momt_spectra_un_tk)
            close(KER_spectra_un_tk)
            close(momt_spectra_tk)
            close(KER_spectra_tk)
        enddo

        write(filepath,'(a,a)') adjustl(trim(time_prop_dir_1d)), "Total_KER_spectra.out"
        open(newunit=KER_spectra_un_tk,file=filepath,status='unknown')
        
        write(filepath,'(a,a)') adjustl(trim(time_prop_dir_1d)), "Total_momt_spectra.out"
        open(newunit=momt_spectra_un_tk,file=filepath,status='unknown')
        
        write(filepath,'(a,a)') adjustl(trim(time_prop_dir_1d)), "Total_KER_spectra_normalized.out"
        open(newunit=KER_spectra_tk,file=filepath,status='unknown')
        
        write(filepath,'(a,a)') adjustl(trim(time_prop_dir_1d)), "Total_momt_spectra_normalized.out"
        open(newunit=momt_spectra_tk,file=filepath,status='unknown')

        do i=NR/2 +1,NR
            write(momt_spectra_un_tk,*) pR(i), sum(abs(psi_diss(i,:))**2), sum(abs(this%psi_outR(i,:))**2), &
                & sum(abs(this%psi_outR_inc(i,:)))
            write(momt_spectra_tk,*) pR(i), sum(abs(psi_diss(i,:)/sqrt(norm_diss(:)))**2), &
                & sum(abs(this%psi_outR(i,:)/sqrt(norm_outR(:)))**2)
        enddo
        do i=1,NR/2
            write(momt_spectra_un_tk,*) pR(i), sum(abs(psi_diss(i,:))**2), sum(abs(this%psi_outR(i,:))**2), &
                & sum(abs(this%psi_outR_inc(i,:)))
            write(momt_spectra_tk,*) pR(i), sum(abs(psi_diss(i,:)/sqrt(norm_diss(:)))**2), &
                & sum(abs(this%psi_outR(i,:)/sqrt(norm_outR(:)))**2)
            write(KER_spectra_un_tk,*) pR(i)**2/(2*m_red), m_red*sum(abs(psi_diss(i,:))**2) / pR(i), &
                & m_red*sum(abs(this%psi_outR(i,:))**2) / pR(i), m_red*sum(abs(this%psi_outR_inc(i,:)))/pR(i)
            write(KER_spectra_tk,*) pR(i)**2/(2*m_red), m_red*sum(abs(psi_diss(i,:)/sqrt(norm_diss(:)))**2) /pR(i), &
                & m_red*sum(abs(this%psi_outR(i,:)/sqrt(norm_outR(:)))**2) / pR(i)
        enddo
        close(momt_spectra_un_tk)
        close(KER_spectra_un_tk)
        close(momt_spectra_tk)
        close(KER_spectra_tk)
        ! Destroy FFTW plans and free memory
        call fftw_destroy_plan(continuum_1d%planF)
        call fftw_destroy_plan(continuum_1d%planB)
        call fftw_free(continuum_1d%p_in)
        call fftw_free(continuum_1d%p_out)
        
    end subroutine post_prop_analysis

    !> Clean up and array deallocation
    subroutine deallocate_all(this)
        class(time_prop), intent(inout):: this

        print*
        print*, "Cleaning up time propagation variables ..."
        if(allocated(this%chi0)) deallocate(this%chi0)
        if(allocated(this%psi_ges)) deallocate(this%psi_ges)
        if(allocated(this%vib_en)) deallocate(this%vib_en)
        if(allocated(this%psi_chi)) deallocate(this%psi_chi)
        if(allocated(this%abs_func)) deallocate(this%abs_func)
        if(allocated(this%psi_outR)) deallocate(this%psi_outR)
        if(allocated(this%psi_outR_inc)) deallocate(this%psi_outR_inc)
        print*, "Done."
    end subroutine deallocate_all

    !------------------------------------------------------------------
    !-------- Helper functions ----------------------------------------
    !------------------------------------------------------------------
    !> Checks file status for safe reading
    subroutine file_status_check(filepath)
        integer:: size1, size2
        logical:: EXST
        character(len=*)::filepath

        do 
            inquire(file=adjustl(trim(filepath)), exist=EXST, size=size1)
            if (EXST .and. size1 > 0) then
                call sleep(1)
                inquire(file=adjustl(trim(filepath)), size=size2)
                if (size1 == size2) exit 
            else
                call sleep(1)
            endif
        enddo
        print*, "Status of ", adjustl(trim(filepath)), " checked"
    end subroutine file_status_check

    !> Integrates norm of wavefunction over grid
    subroutine integ(psi, norm)
         
        use global_vars, only:NR, Nstates, dR, dp
        implicit none
        integer J
         
        real(dp) norm(Nstates)
        complex(dp) psi(NR,Nstates)
         
        norm = 0.d0
    
        do J = 1, Nstates
            norm(J)= sum(abs(psi(:,J))**2)   
        end do
    
        norm = norm * dR
        return 
    end subroutine

    !> Calculates Boltzmann distribution for vibrational populations
    subroutine Boltzmann_distribution(N, E, Boltzmann_populations)
        use global_vars, only: temperature, dp, Vstates, guess_vstates
        use data_au, only: kB

        integer:: N, V
        real(dp), intent(in):: E(guess_vstates)
        real(dp), allocatable, intent(out):: Boltzmann_populations(:)
        real(dp):: total_pop
           
        total_pop = sum(exp(-E(:)/(kB*temperature)))
        print*, "total populations =", total_pop
   
        allocate(Boltzmann_populations(Vstates(N)))
        do V = 1, Vstates(N)
            if (V .gt. 1) then
                Boltzmann_populations(V) = exp(-(E(V)-E(1))/(kB*temperature))!/total_pop
            else
                Boltzmann_populations(V)=1._dp
            endif
            print*, "state", V-1, "population ratio =", Boltzmann_populations(V)
        enddo
   
        total_pop = sum(Boltzmann_populations(:))
        print*, "sum of ratios =", total_pop
        Boltzmann_populations = Boltzmann_populations/total_pop
        total_pop = sum(Boltzmann_populations(:))
        print*, "total population =", total_pop
   
        do V =1, Vstates(N)
            print*, "state", V-1, "population probability =", Boltzmann_populations(V)
        enddo  

    end subroutine Boltzmann_distribution

    !> Generates complex absorber function for boundary
    subroutine complex_absorber_function(v_abs, f)
        use global_vars, only: NR, dp, dt, R
        use pot_param, only: cpmR
        
        integer i
        real(dp):: a, eps, V_abs(NR), n, R0, p
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

    !> Generates mask absorber function (cosine)
    subroutine mask_function_cos(cof)
        use global_vars, only:NR, R, dR, dp
        use data_au
        use pot_param
   
        integer :: J
        real(dp):: cof(NR)
   
        do j = 1, NR
            R(J) = R0 + (j - 1) * dR
            if(R(J).lt.(Rend - cpmR)) then
                cof(j) = 1.d0
            else
                cof(j) = cos(((R(J) - Rend + cpmR) / -cpmR) * (0.5d0 * pi))
                cof(j) = cof(j)**2
            end if
        end do
        return
    end subroutine
    !------------------------------------------------
    !> Generates mask absorber function (exponential)
    subroutine mask_function_ex(cof)
        use global_vars, only:NR, R, dR, dp
        use data_au
        use pot_param
    
        integer :: J
        real(dp):: cof(NR),c
  
        c=1.00d0
        do j = 1, NR
            R(J) = R0 + (j - 1) * dR
            cof(j)=1.0d0/(1.0d0+exp(c*(R(J)-Rend+cpmR)))
        end do
  
        return
    end subroutine

    !> Generates pulse matrix for interstate coupling
    subroutine pulse2(tout,mu,E)
        use global_vars, only:dt, Nstates,kap, dp
        use data_au, only:im
        use blas_interfaces_module, only : zgemv, dgemv
   
        integer:: i, J
        real(dp):: u, uv, d, mu(Nstates,Nstates), E
        integer info, Lwork
   
        complex(dp):: tout, b, z, p, pT
        dimension:: u(Nstates,Nstates), d(Nstates), b(Nstates,Nstates), &
            & z(Nstates,Nstates), tout(Nstates,Nstates), &
            & p(Nstates,Nstates), uv(Nstates,Nstates), &
            & pT(Nstates,Nstates)
        real(dp) work(1000) 
        character(len=1):: JOBZ
   
       
        !u(1,1) = 0.d0
        !u(1,2) = -kap*mu(1,2) * E
        !u(2,1) = -kap*mu(2,1) * E
        !u(2,2) = 0.d0
   
        !Diapole matrix
        u=0._dp
        do i=1, Nstates-1
            do J=i+1, Nstates
                u(i,J)= -kap*mu(i,J) * E
                u(J,i)= -kap*mu(i,J) * E
            enddo
        enddo
        !print*, "u"
        !write( * , * ) ((u(i,j),j=1,Nstates), i=1,Nstates )
    
        uv = 0._dp
        uv = u
   
        Lwork=-1
        JOBZ='V'
        call dsyev(JOBZ,'U', int(Nstates), uv, int(Nstates), &
            & d, work, Lwork,info)
        ! call jacobi(u,Nstates,d)
        Lwork = min(1000, int(work(1)))
        
        !---- Solve eigenproblem.
        call dsyev('V', 'U', int(Nstates), uv, int(Nstates), &
            & d, work, Lwork, info)
   
        if( info.gt.0 ) then
            write(*,*)'The algorithm failed to compute eigenvalues.'
            stop
        endif
    
        p = (0._dp,0._dp)
        p = uv ! transfering eigenvectors to a complex array
        !print*,"eigen vector p"
        !write( * , * ) ( (p(i,j),j=1,Nstates), i=1,Nstates )
   
        b= (0._dp,0._dp)
        do J = 1,Nstates  
            b(J,J) = exp(-im * dt * d(J)) ! exponentiate diagonal matrix i.e. eigenvalues
        end do
        !print*, "b"
        !write( * , * ) ((b(i,j),j=1,Nstates), i=1,Nstates )
   
        ! Calculating e^u = P (e^d) P^(-1)
        !z = matmul(p,b)
        call zgemm('N', 'N', int(Nstates), int(Nstates), &
            & int(Nstates),(1._dp,0._dp), p, size(p,dim=1), &
            & b, size(b,dim=1), (0._dp,0._dp), z, size(z,dim=1))
        !print*, "z"
        !write( * , * ) ((z(i,j),j=1,Nstates), i=1,Nstates )
   
        pT = transpose(p)
        !tout = matmul(z,transpose(p))
        call zgemm('N', 'N', int(Nstates), int(Nstates), &
            & int(Nstates),(1._dp,0._dp), z, size(z,dim=1), &
            & pT, size(pT,dim=1), (0._dp,0._dp), tout, size(tout,dim=1))
        !print*, "tout"
        !write( * , * ) ((tout(i,j),j=1,Nstates), i=1,Nstates )
   
    return
    end subroutine pulse2

end module propagation1d_mod


