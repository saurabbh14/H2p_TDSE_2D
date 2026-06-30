!> This module calculates the vibrational energies and states for given electronic states
module nuclear_wavefkt

    use global_vars
    use pot_param
    use data_au
    use FFTW3
    use omp_lib
    implicit none
    private
    public :: nuclear_wavefkt_class, nuclear_wf_calc

    !> Class to hold all data and methods for nuclear wavefunction calculations
    type :: nuclear_wavefkt_class
        integer :: guess_Vstates                ! Number of maximum guesse vibrational states
        real(dp) :: RI, kappa                   ! Initial Gaussian center, kappa parameter
        character(10) :: ITP_par_FFTW           ! FFTW parameter string for parallelization
        logical, allocatable :: Files_exist(:)  ! Flags for existence of output files per state
        integer :: vstates_tk                   ! File handle for vibrational states
        integer, allocatable, dimension(:) :: vib_en_tk, chi0_vib_en_tk, chi0_tk ! File handles per state
        integer :: not_converged_tk             ! File handle for unconverged states
    contains
        procedure :: read_guess_wp_params       ! Read initial parameters from global variables
        procedure :: initialize_wp_params       ! Initialize/allocate arrays and convert units
        procedure :: open_files, open_not_converged_files ! Open output files for results
        procedure :: close_files                ! Close all opened files
        procedure :: imaginary_time_propagation => ITP ! Main vibrational state calculation
        procedure :: nuclear_wf_calc            ! High-level wrapper for calculation
        procedure :: deallocate_all             ! Deallocate arrays and clean up
    end type nuclear_wavefkt_class
  
contains

    !> High-level routine to run the vibrational state calculation
    subroutine nuclear_wf_calc(this)
        class(nuclear_wavefkt_class), intent(inout) :: this
        call this%read_guess_wp_params()
        call this%initialize_wp_params()
        call this%imaginary_time_propagation()
        call this%deallocate_all()
        print*, "Leaving vibrational state calculations ..."
    end subroutine nuclear_wf_calc  

    !> Read initial guess parameters from global variables
    subroutine read_guess_wp_params(this)
        class(nuclear_wavefkt_class), intent(inout) :: this
        this%guess_Vstates = guess_vstates
        this%RI = RI
        this%kappa = kappa
        this%ITP_par_FFTW = ITP_par_FFTW
    end subroutine read_guess_wp_params

    !> Initialize/allocate arrays and convert units for calculation
    subroutine initialize_wp_params(this)
        use data_au, only: au2a
        class(nuclear_wavefkt_class), intent(inout) :: this
        ! Convert initial Gaussian center from Angstrom (or given units) to atomic units
        this%RI = this%RI / au2a
        allocate(this%Files_exist(Nstates))
        this%Files_exist = .false.
        allocate(this%vib_en_tk(Nstates), this%chi0_vib_en_tk(Nstates), this%chi0_tk(Nstates))
    end subroutine initialize_wp_params

    !> Open output files for vibrational state results
    subroutine open_files(this)
        use global_vars, only: Nstates, nucl_wf_dir
        class(nuclear_wavefkt_class), intent(inout) :: this
        logical :: Ext
        integer :: j
        character(150) :: filename

        ! Open file for total vibrational states
        write(filename,'(a,a,i0,a)') adjustl(trim(nucl_wf_dir)), &
                & "Bound-vibstates_in_Nthstates.out"
        open(newunit=this%vstates_tk,file=filename,status='unknown')

        Nloop: do j = 1, Nstates ! varying the different adiabatic states
            ! Check if output file for this state exists
            write(filename,'(a,a,i0,a)') adjustl(trim(nucl_wf_dir)), &
                & "BO_Electronic-state-g", int(j-1), "_vibstates.out"
            inquire(file=filename, Exist=Ext)
            print*, trim(filename), " ", Ext
            this%Files_exist(j) = Ext
            
            if (Ext) cycle ! Skip if file exists

            ! Open files for vibrational energies and wavefunctions
            write(filename,'(a,a,i0,a)') adjustl(trim(nucl_wf_dir)), &
                & "BO_Electronic-state-g", int(j-1), "_Evib.out"
            open(newunit=this%vib_en_tk(j),file=filename,status='unknown')
            write(filename,'(a,a,i0,a)') adjustl(trim(nucl_wf_dir)), &
                & "BO_Electronic-state-g", int(j-1), "_chi0-Evib.out"
            open(newunit=this%chi0_vib_en_tk(j),file=filename,status='unknown')

            write(filename,'(a,a,i0,a)') adjustl(trim(nucl_wf_dir)), &
              & "BO_Electronic-state-g", int(J-1), "_vibstates.out"
            open(newunit=this%chi0_tk(j),file=filename,status='unknown')

        end do Nloop
    end subroutine open_files

    !> Close all open files
    subroutine close_files(this)
        class(nuclear_wavefkt_class), intent(inout) :: this
        integer :: j

        close(this%vstates_tk)

        do j = 1, Nstates
            if (.not. this%Files_exist(j)) then
                close(this%vib_en_tk(j))
                close(this%chi0_vib_en_tk(j))
                close(this%chi0_tk(j))
            end if
        end do
    end subroutine close_files

    !> Open file for unconverged vibrational state wavefunction
    subroutine open_not_converged_files(this, v, j)
        class(nuclear_wavefkt_class), intent(inout) :: this
        integer :: j, v
        character(1000) :: filename

        write(filename,'(a,a,i0,a,i0,a)') adjustl(trim(nucl_wf_dir)), &
            & "Electronic-state-g", int(j-1), "_vib-state", int(v-1), &
            & "_uncoverged-final-wf.out"
        open(newunit=this%not_converged_tk,file=filename,status='unknown')
    end subroutine open_not_converged_files

    !> Main routine for vibrational state calculation using Imaginary Time Propagation (ITP)
    subroutine ITP(this)
        use global_vars, only: NR, Nstates, dp, R, adb, m_red, &
            & Vstates, chi0
        use data_au, only: au2eV, eV2nm
        class(nuclear_wavefkt_class), intent(inout) :: this
        integer :: i, j, V, G, K
        integer:: istep  
        type(C_PTR):: planF, planB, p_in, p_out
        real(dp):: dt2
        real(dp):: E1, norm
        real(dp), allocatable:: E(:,:)
        real(dp):: thresh
        real(dp):: Rin

      !  double precision, allocatable:: chi0(:,:)
      !  double precision:: total_pop, Boltzmann_populations(guess_Vstates) 
      !  double precision:: trans_dipole(guess_Vstates,guess_Vstates) 
 
        ! Working arrays for wavefunction and reference states
        real(dp), allocatable, dimension(:):: vprop
        real(dp), allocatable, dimension(:):: psi1
        real(C_DOUBLE), allocatable, dimension(:):: psi
        real(C_DOUBLE), pointer:: psi_in(:), psi_out(:)
        real(dp), allocatable, dimension(:,:):: ref

        ! Allocate arrays for calculation
        allocate(psi(NR),psi1(NR))
        allocate(vprop(NR), E(guess_vstates,Nstates))
        allocate(ref(NR,guess_vstates))
        
        print*
        print*, "FFTW intialization ..."
        print*
        ! Creating aligned memory for FFTW
        p_in = fftw_alloc_real(int(NR, C_SIZE_T)) 
        call c_f_pointer(p_in,psi_in,[NR])
        p_out = fftw_alloc_real(int(NR, C_SIZE_T)) 
        call c_f_pointer(p_out,psi_out,[NR])

        call fftw_initialize_threads
        print*, "FFTW plan creation..."
        call fftw_create_r2r_plans(psi_in, psi_out, NR, planF, planB, ITP_par_FFTW)
        print*
        print*, "Done setting up FFTW."

        dt2 = dt!*10        ! time step for ITP
        thresh = 1d-18      ! convergence threshold
        istep = 1e8         ! max ITP steps
        print*
        print*, "reduced mass:", m_red
        print*
        print*, 'Start of energy calculation...'
        print*
        E = 0._dp
        call this%open_files()

        Nloop: do J = 1, Nstates ! Loop over electronic states
            if (this%Files_exist(J)) then
                print*, "Vibrational states already computed for state ", J-1, ". Skipping ITP."
                cycle
            endif

            ! Prepare potential for propagation
            do i = 1, NR
                vprop(i) = exp(-0.5_dp * dt2 * adb(i,J))
            end do
            Rin = R(minloc(adb(:,J),1))

            print*, 'Dp =', adb(NR,J)*au2eV, 'eV'
            
            Vloop: do V = 1, guess_Vstates ! loop over the vibrational state

                ! Initialize wavefunction for this vibrational state
                do i = 1, NR  ! symmetry for the startup function
                    psi(i) = exp(this%kappa * (R(I) -Rin )**2) + &
                    &   (-1._dp)**(V - 1) * exp(-0.5_dp * (R(I) + Rin)**2)
                end do

                call integ_r(psi, psi, norm)
                psi = psi / sqrt(norm)
                !.......... Imaginary Time Propagation ........................
                timeloop: do K = 1, istep

                    psi1 = psi        ! storing wave function of iteration step (N - 1)
                    E1 = E(V,J)            ! storing eigenvalue of iteration step (N - 1)
                    if (V.gt.1) then   !projecting out the vibrational ground state...
                        do G = 1, (V - 1)

                        call integ_r(ref(1:NR,G), psi, norm)

                            do i = 1, NR
                                psi(i) = psi(i) - norm * ref(i, G)
                            end do

                        end do
                    end if

                    psi = psi * vprop
                    psi_in = psi
                    call fftw_execute_r2r(planF, psi_in, psi_out)
                    psi = psi_out   
                    psi = psi * exp((-dt2 * pR**2) / (2._dp*m_red)) 
                    psi_in = psi  
                    call fftw_execute_r2r(planB, psi_in, psi_out)
                    psi = psi_out   
                    psi = psi / dble(Nr)         
                    psi = psi * vprop

                    call eigenvalue_r(psi, psi1, E(V,J), dt2)
                    call integ_r(psi, psi, norm)

                    psi = psi / sqrt(norm)

                    ! Check for convergence
                    if(abs(E(V,J) - E1).le.thresh) then
                        print*, 'Surface', J-1, 'Vibrational state', V-1, E(V,J) * au2eV, 'eV'
                        print*, 'Resonance freq for dissociation', (adb(NR,J)-E(V,J))*au2eV, 'eV'
                        print*, 'Resonance freq for dissociation', eV2nm/((adb(NR,J)-E(V,J))*au2eV), 'nm'
                        if (V.gt.1) then
                            do G=1,(V-1)
                                print*, "freq for the transition from", G-1, "to", V-1, ":", &
                                    & eV2nm/(abs(E(G,J)-E(V,J))*au2eV), "nm E=", &
                                    & abs(E(G,J)-E(V,J))*au2eV, "eV"
                            enddo
                        endif
                        write(this%vib_en_tk(j),*) v-1, e(V,J)*au2ev
                        do I = 1, NR
                            ref(I,V) = psi(I)             ! storing as reference for the next loop
                            write(this%chi0_vib_en_tk(j),*) R(I), ref(I,V)+E(V,J)*au2eV     
                        end do
                        write(this%chi0_vib_en_tk(j),*)
                        exit
                    elseif(K .eq. istep) then
                        ! If not converged, write out unconverged wavefunction and stop
                        print*,'Iteration not converged!'
                        print*,'Program stopped!'
                        print*
                        print*,'Surface:', J-1
                        print*,'Vibstate:', V-1
                        print*,'E =', E(V,J) *au2eV
                        print*,'E1 =', E1 *au2eV
                        print*,'thresh =', thresh *au2eV
                        print*,'step =', K
                        call this%open_not_converged_files(v, j)
                        do i = 1, nr
                            write(this%not_converged_tk,*) R(I), psi(i), abs(psi(I))**2
                        end do
                        close(this%not_converged_tk)
                        stop
                    else
                        cycle
                    end if
                end do timeloop
                if (E(V,J) .gt. adb(NR,J)) then
                    Vstates(J) = v-1
                    exit
                endif
            end do Vloop               ! end of vibrational states loop

            print*, "Number of bound vibrational states in the ", int(j-1), &
                & " electronic state:", Vstates(j)
            write(this%vstates_tk,*) j-1, Vstates(j)
            do v = 1, vstates(j)
                do i = 1, NR
                    chi0(i, v, j) = ref(i, v)
                enddo
            enddo
            
            do i = 1, NR
                write(this%chi0_tk(j),*) R(i), chi0(i,1:vstates(j), j)
            enddo
        end do Nloop            ! end of surface loop
        call this%close_files()

        ! Clean up FFTW plans and memory
        call fftw_destroy_plan(planF)
        call fftw_destroy_plan(planB)
        call fftw_free(p_in)
        call fftw_free(p_out)
        deallocate(psi, psi1, vprop, ref)
   
    end subroutine ITP

    !> A subroutine for deallocating all arrays
    subroutine deallocate_all(this)
        class(nuclear_wavefkt_class), intent(inout) :: this

        print*
        print*, "Cleaning up nuclear_wp variables ..."
        if (allocated(this%Files_exist)) deallocate(this%Files_exist)
        if (allocated(this%vib_en_tk)) deallocate(this%vib_en_tk)
        if (allocated(this%chi0_vib_en_tk)) deallocate(this%chi0_vib_en_tk)
        if (allocated(this%chi0_tk)) deallocate(this%chi0_tk)
        print*,"Done"
    end subroutine deallocate_all

    !_________________ Helper Subroutines______________________________________
    
    !> Calculate eigenvalue from two wavefunctions using log ratio
    subroutine eigenvalue_R(A, B, E, dt2)      
      
        use global_vars, only: NR, dp
        real(dp):: E, e1, e2, norm
        real(dp), intent(in):: dt2, A(nr), B(nr)
           
        call integ_r(B, B, norm)  
        e1 = norm
  
        call integ_r(A, A, norm)  
        e2 = norm
   
        E = (-0.5d0/dt2) * log(e2/e1)
 
        return
    end subroutine eigenvalue_R    
  
    !> Integrate product of two arrays over grid                                                          
    subroutine integ_r(A, B, C)

        use global_vars, only:NR, dR, dp   
        integer i 
        real(dp),intent(in):: A(Nr), B(Nr)
        real(dp):: C
  
        C = 0.d0
  
        do i = 1, Nr  
            C = C + A(i) * B(i)   
        end do
  
        C = C * dr
  
        return  
      end subroutine integ_r

end module nuclear_wavefkt

!!! COMMENTED Boltzmann Distribution part
! total_pop = sum(exp(-E(:)/(kB*temperature)))
! print*, "total populations =", total_pop

! do V = 1, Vstates
!   if (V .gt. 1) then
!   Boltzmann_populations(V) = exp(-(E(V)-E(1))/(kB*temperature))!/total_pop
!   else
!   Boltzmann_populations(V)=1.d0
!   endif
!   print*, "state", V-1, "population ratio =", Boltzmann_populations(V)
! enddo
!
! total_pop = sum(Boltzmann_populations(:))
! print*, "sum of ratios =", total_pop
! Boltzmann_populations = Boltzmann_populations/total_pop
! total_pop = sum(Boltzmann_populations(:))
! print*, "total population =", total_pop
!
! do V =1, Vstates
!   print*, "state", V-1, "population probability =", Boltzmann_populations(V)
! enddo
! 
! print*
! do V=1, Vstates
!   write(filename,'(a,i0,a)') "trans_dipole_v_",V,".out"
!   open(unit=trans_dipole_tk, file=filename, status="unknown")
!   !if (V .gt. 1) then
!    do G = V,Vstates
!    trans_dipole(G,V) = abs(sum(ref(:,G)*(kap*(mu_all(:,1,1)+(mr-0.5d0)*R(:))+lam*R(:))*ref(:,V))*dR) !**2
!  !   trans_dipole(G,V) = abs(sum(ref(:,G)*(kap*mu_all(:,1,1)+lam*R(:))*ref(:,V))*dR) !**2
!     print*, 'transition dipole', V, "to", G, ":", trans_dipole(G,V)
!     write(trans_dipole_tk,*) G, trans_dipole(G,V) 
!    enddo
!   !endif
!   print*
! enddo

