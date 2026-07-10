!> Initializes directories, grids, masses and reads potentials/dipoles.
!> Allocations and unit conversions happen here.

module initializer
    implicit none
    private
    
    public :: initializer_setup

contains

    !> Prepare working directories, grids and arrays, read potentials/dipoles.
    subroutine initializer_setup
        use output_dir_mod, only: output_dir_check, setup_output_dir
        use setpot_mod, only: build_2d_potential
        call output_dir_check()
        call setup_output_dir()
        call allocate_arrays()
        
        ! Read potentials and transition dipoles (or generate synthetic ones)
        call sc_read()
        !call bo_pot_read()
        !call trans_dipole_read()

        ! Prepare grids and related parameters
        !call into_atomic_units()
        call r_grid()
        call p_grid()
        call x_grid()
        call px_grid()

        call time_grid()
        call mass_setup()

        ! Build the 2D softcore potential for the electronic coordinate x and internuclear coordinate R
        call build_2d_potential()
    end subroutine initializer_setup

    subroutine allocate_arrays
        use global_vars, only: adb, mu_all, R, pR, x, px, pot, ewf, &
             & NR, Nx, Nstates, zeff, alpha2
        ! Allocate arrays for potentials, dipoles and grids
        allocate(R(NR), x(Nx))
        allocate(pR(NR), px(Nx))
        allocate(zeff(NR), alpha2(NR))
        allocate(pot(NR,Nx), ewf(NR,Nx,Nstates))
        allocate(adb(NR,Nstates),mu_all(Nstates,Nstates,NR))
    end subroutine allocate_arrays

    ! Position grid (internuclear distance)
    subroutine r_grid
        use global_vars, only: R, NR, dR, dpR, dp
        use data_au, only: pi
        use pot_param, only: R0, Rend
        ! Set derived grid parameters and unit conversions
        R0 = R(1)                  ! leftmost grid point (coordinate space)
        Rend = R(NR)               ! rightmost grid point
        dR = R(2) - R(1)           ! grid spacing (assumed uniform)
        dpr = (2._dp * pi) / (dR * NR)  ! momentum-grid spacing via FFT conventions
    end subroutine r_grid

    ! Position grid (internuclear distance)
    subroutine x_grid
        use global_vars, only: x, Nx, dx, dpx, dp, xmin, xmax
        use data_au, only: pi
        use pot_param, only: x0, xend
        integer i
        ! Set derived grid parameters and unit conversions
        x0 = xmin                  ! leftmost grid point (coordinate space)
        xend = xmax               ! rightmost grid point
        dx = (xend - x0) / (Nx-1)           ! grid spacing (assumed uniform)
        dpx = (2._dp * pi) / (dx * Nx)  ! momentum-grid spacing via FFT conventions
        do i = 1, Nx
            x(i) = x0 + (i-1) * dx
        enddo

    end subroutine x_grid
    
    ! Calculating relavent masses in a diatomic system (or a system divided into two mass blocks)
    subroutine mass_setup
        use global_vars, only: m1, m2, mn, m_red, m_eff, mn1, mn2, lam, kap, dp
        use data_au, only: mass
        ! Masses: convert to internal mass units and compute reduced/effective masses
        m1 = m1 * mass
        m2 = m2 * mass
        mn = m1 + m2
        m_red = m1*m2/(m1+m2)
        ! Effective mass used in some reduced-dimension expressions (keeps consistency)
        m_eff = (m1 + m2) / (m1 + m2 + 1.0_dp)
        mn1 = m1 / mn
        mn2 = m2 / mn
        
        ! diapole parameters
        kap = (mn + 2._dp) / (mn + 1._dp)
        lam = (m2 - m1) / mn

    end subroutine mass_setup

    subroutine into_atomic_units
        use global_vars, only: dt
        use data_au, only: au2fs

        ! Convert time step from femtoseconds to atomic units for propagation
        dt = dt / au2fs
                
    end subroutine into_atomic_units
 
    !...................... Impulsgrid......................
    !> Construct the momentum grid consistent with the FFT layout used in split-operator.
    !> Uses the standard ordering: 0, dp, 2dp, ..., (NR/2-1)dp, -NR/2 dp, ..., -dp 

    subroutine p_grid
        use global_vars, only:pR, dpR, NR, R
        use pot_param
        integer:: I 
  
        do I = 1, NR  
            if (I.le.(NR / 2)) then    
                PR(I) = (I - 1) * dpR    
            else    
                PR(I) = - (NR + 1 - I) * dpR    
            end if
        end do
        print'(a,f4.2,a,f6.2)', 'R0=', R(1), ', Rend=', R(NR)
        print'(a,f7.2,a,f7.2)', 'PR0=', PR((NR/2)+1), ', PRend=', PR(NR/2)
        return
    end subroutine  
    subroutine px_grid
        use global_vars, only:px, dpx, Nx, x
        use pot_param
        integer:: I 
  
        do I = 1, Nx  
            if (I.le.(Nx / 2)) then    
                Px(I) = (I - 1) * dpx    
            else    
                Px(I) = - (Nx + 1 - I) * dpx    
            end if
        end do
        print'(a,f7.2,a,f6.2)', 'x0=', x(1), ', xend=', x(Nx)
        print'(a,f6.2,a,f5.2)', 'Px0=', Px((Nx/2)+1), ', Pxend=', Px(Nx/2)
        return
    end subroutine  

    ! Time-grid in a.u.
    subroutine time_grid
        use global_vars, only: dt, Nt, time
        integer:: k
        allocate(time(Nt))
        do k = 1, Nt
            time(k) = k * dt
        end do
    end subroutine time_grid

    !------------------------------------------------------------------------------
    !%%%%%% File IO: reading potentials and dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !------------------------------------------------------------------------------

    !> Read effective charge and softcore parameters from input file or
    !> use static zeff and alpha2 value when sc_kind="static".
    subroutine sc_read
        use global_vars, only:R, NR, sc_params, &
                & input_data_dir, output_data_dir, &
                & sc_kind, dp, Rmin, Rmax, dR, zeff, &
                & alpha2
        use data_au
        use pot_param, only: morse_potential

        character(4000):: filepath
        integer:: I, sc_tk, sc_out_tk

        select case(sc_kind)
            case ("on_grid")
                ! Compose full path and read potential data file with NR lines.
                write(filepath,'(a,a)') adjustl(trim(input_data_dir)), adjustl(trim(sc_params))  
                print*, "sc parameters are in path:", trim(filepath)
                open(newunit=sc_tk,file=adjustl(trim(filepath)),status='unknown')
                read(sc_tk,*)
                do I = 1, NR
                    read(sc_tk,*) R(I), zeff(I), alpha2(I) 
                end do
                close(sc_tk)

                ! Also write a copy into output directory for verification
                write(filepath,'(a,a,a)') adjustl(trim(output_data_dir)), adjustl(trim(sc_params)),&
                        & "_read.out"  
                open(newunit=sc_out_tk,file=adjustl(trim(filepath)),status='unknown')
                do I = 1, NR
                    write(sc_out_tk,*) R(I), zeff(I), alpha2(I) 
                end do
                close(sc_out_tk)

            case ("old_sc")
                ! Compose full path and read potential data file with NR lines.
                write(filepath,'(a,a)') adjustl(trim(input_data_dir)), adjustl(trim(sc_params))  
                print*, "sc parameters are in path:", trim(filepath)
                open(newunit=sc_tk,file=adjustl(trim(filepath)),status='unknown')
                !read(sc_tk,*)
                zeff = 1._dp
                do I = 1, NR
                    read(sc_tk,*) R(I), alpha2(I) 
                    read(sc_tk,*)
                end do
                close(sc_tk)
                R = R / au2a

                ! Also write a copy into output directory for verification
                write(filepath,'(a,a,a)') adjustl(trim(output_data_dir)), adjustl(trim(sc_params)),&
                        & "_read.out"  
                open(newunit=sc_out_tk,file=adjustl(trim(filepath)),status='unknown')
                do I = 1, NR
                    write(sc_out_tk,*) R(I), zeff(I), alpha2(I) 
                end do
                close(sc_out_tk)
 
            case("static")
                dR = (51.2_dp - 0.1_dp) / (NR-1)
                do I = 1, NR
                    R(I) = 0.1_dp + (I-1) * dR
                enddo

                alpha2 = 1._dp
                zeff = 1._dp
                ! Write generated Morse surface to output 
                write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "sc_params_read.out"  
                open(newunit=sc_out_tk,file=adjustl(trim(filepath)),status='unknown')
                do I = 1, NR
                    write(sc_out_tk,*) R(I), zeff(I), alpha2(I) !, sngl(adb(I,2)*au2eV), &
                    ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
                end do
                close(sc_out_tk)
        end select
    end subroutine sc_read

    !> Read electronic Born–Oppenheimer potential surfaces from input file or
    !> construct a Morse potential when bo_pot_kind="Morse".
    subroutine bo_pot_read
        use global_vars, only:R, NR, adb, adb_pot, &
                & input_data_dir, output_data_dir, &
                & bo_pot_kind, dp, Rmin, Rmax, dR
        use data_au
        use pot_param, only: morse_potential

        character(2000):: filepath
        real(dp):: dummy
        integer:: I, pot_tk, pot_out_tk

        select case(bo_pot_kind)
            case ("on_nuclr_grid")
                ! Compose full path and read potential data file with NR lines.
                write(filepath,'(a,a)') adjustl(trim(input_data_dir)), adjustl(trim(adb_pot))  
                print*, "BO Potential surfaces in path:", trim(filepath)
                open(newunit=pot_tk,file=adjustl(trim(filepath)),status='unknown')
                do I = 1, NR
                    read(pot_tk,*) dummy, adb(I,:) !, sngl(adb(I,2)*au2eV), &
                        ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
                end do
                close(pot_tk)

                ! Also write a copy into output directory for verification
                write(filepath,'(a,a,a)') adjustl(trim(output_data_dir)), adjustl(trim(adb_pot)),&
                        & "_read.out"  
                open(newunit=pot_out_tk,file=adjustl(trim(filepath)),status='unknown')
                do I = 1, NR
                    write(pot_out_tk,*) R(I), adb(I,:) !, sngl(adb(I,2)*au2eV), &
                        ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
                end do
                close(pot_out_tk)
 
            case("Morse")
                ! Make R-grid spacing
                dR = (Rmax-Rmin)/(NR-1)
                ! Fill only the ground state with a Morse potential at specific grid points
                do i = 1, NR
                    R(i) = Rmin + (i-1)*dR
                    adb(I,1) = morse_potential(0.17_dp,1.85_dp,0.743_dp/au2a,R(I))
                enddo
   
                ! Write generated Morse surface to output 
                write(filepath,'(a,a)') adjustl(trim(output_data_dir)), "Morse_pot_read.out"  
                open(newunit=pot_out_tk,file=adjustl(trim(filepath)),status='unknown')
                do I = 1, NR
                    write(pot_out_tk,*) R(I), adb(I,:) !, sngl(adb(I,2)*au2eV), &
                    ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
                end do
                close(pot_out_tk)
        end select
    end subroutine bo_pot_read
    
    !> Read transition dipole moments for each electronic-state pair.
    !> Supports optional prefixing and allows switching off specific transitions.
    subroutine trans_dipole_read
        use global_vars, only:R, NR, mu_all, trans_dip_prefix, Nstates, &
            & input_data_dir, output_data_dir, dp, total_trans_off, &
            & trans_off

        integer:: I, L, M, N1, N2
        character(2000):: fn
        character(2):: trans_off_parse(total_trans_off)
        character(2):: tr
        integer:: input_tk, output_tk
        real(dp):: dummy
  
        ! Parse list of transitions to disable (e.g. "12 23")
        trans_off = trim(adjustl(trans_off))
        read(trans_off,*) trans_off_parse

        print*, "Transitions to be switched off"
        do L =1, total_trans_off
            trans_off_parse(L) = trim(adjustl(trans_off_parse(L)))
            print*, "Transition", L,": ", trans_off_parse(L) 
        enddo
   
        print*, "Transition dipoles with file prefix \", trim(trans_dip_prefix), " \."
        mu_all = 0._dp

        ! Two supported modes:
        !  - no prefix: files are named input_data_dir + "<L><M>.dat"
        !  - with prefix: files are input_data_dir + trans_dip_prefix + "<L><M>.dat"
        if (trim(trans_dip_prefix) .eq.'') then
            do L = 1, Nstates
                do M = L+1, Nstates
                    write(fn,fmt='(a,i0,i0,a)') adjustl(trim(input_data_dir)),L,M,'.dat'
                    print*, trim(fn)
                    open(newunit=input_tk, file=adjustl(trim(fn)), form='formatted')
                    do I = 1, NR
                        read(input_tk,*) dummy, mu_all(L,M,I)
                    enddo
                    close(input_tk)
                enddo
            enddo
    
            ! Zero-out any transitions explicitly switched off in input (trans_off)
            do L = 1, total_trans_off
                tr = trans_off_parse(L)
                read(tr(1:1),*) N1
                read(tr(2:2),*) N2 
                mu_all(N1,N2,:) = 0._dp 
            enddo
            mu_all = abs(mu_all) ! ensure positive magnitudes
  
            ! Write read dipoles to output for inspection
            do L = 1, Nstates
                do M = L+1, Nstates
                    write(fn,fmt='(a,i0,i0,a)') adjustl(trim(output_data_dir)),L,M,'_read.out'
                    print*, trim(fn)
                    open(newunit=output_tk, file=adjustl(trim(fn)), form='formatted')
                    do I=1,NR
                        write(output_tk,*) R(I), mu_all(L,M,I)
                    enddo
                    close(output_tk)
                enddo
            enddo

        else
            ! Prefix-mode: read files with a common prefix + state indices
            do L = 1, Nstates
                do M = L+1, Nstates
                    write(fn,fmt='(a,a,i0,i0,a)') adjustl(trim(input_data_dir)), &
                        & adjustl(trim(trans_dip_prefix)),L,M,'.dat'
                    print*, trim(fn)
                    open(newunit=input_tk, file=adjustl(trim(fn)), form='formatted')
                    do I = 1, NR
                        read(input_tk,*) dummy, mu_all(L,M,I)
                    enddo
                    close(input_tk)
                enddo
            enddo
 
            ! Write read dipoles to output directory with prefix
            do L = 1, Nstates
                do M = L+1, Nstates
                    write(fn,fmt='(a,a,i0,i0,a)') adjustl(trim(output_data_dir)), &
                        & adjustl(trim(trans_dip_prefix)),L,M,'_read.out'
                    print*, trim(fn)
                    open(newunit=output_tk, file=adjustl(trim(fn)), form='formatted')
                    do I=1,NR
                        write(output_tk,*) R(I), mu_all(L,M,I)
                    enddo
                    close(output_tk)
                enddo
            enddo
        endif

    end subroutine trans_dipole_read



end module initializer
