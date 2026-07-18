module setpot_mod
    use varprecision, only: dp
    use data_au, only: au2a, au2eV, pi
    use global_vars, only: R, x, NR, Nx, mn1, mn2, output_data_dir, &
        & Pot, zeff, alpha2, alpha0, CalcMode, dt
    implicit none
    private
    public :: build_2d_potential, build_kh_potential_at_time
contains
    subroutine build_2d_potential()
        integer :: I, J, pot_unit, ionic_unit
        integer :: t, T_steps
        real(dp) :: time
        real(dp), allocatable :: v1_t(:), v2_t(:), dx1(:), dx2(:)
        real(dp), allocatable :: v12(:), v1e(:), v2e(:)
        real(dp) :: GM
        character(2000) :: pot_file, ionic_file

        allocate(v1_t(NR), v2_t(NR), dx1(NR), dx2(NR))
        allocate(v12(NR), v1e(NR), v2e(NR))

        select case (trim(CalcMode))
        case ("Lab")
            do J = 1, Nx      
                v12(:) = 1._dp/abs(R(:))
                v1e(:) = -zeff(:) / sqrt((x(J) - mn1*R(:))**2 + alpha2(:))  !proton core
                v2e(:) = -zeff(:) / sqrt((x(J) + mn2*R(:))**2 + alpha2(:))  !proton core  
                pot(:,J) = v12(:) + v1e(:) + v2e(:) 
            end do
 
        case ("KH")
            alpha0 = 2.5_dp
            T_steps = int(2*pi / 1._dp / dt) ! omega = 1.0
            v12(:) = 1._dp/abs(R(:))
            do J = 1, Nx
                v1_t = 0._dp
                v2_t = 0._dp
                do t = 1, T_steps
                    time = t * dt  
                    dx1(:) = x(J) - mn1*R(:) + alpha0 * cos(time)
                    dx2(:) = x(J) + mn2*R(:) + alpha0 * cos(time)
                    v1_t(:) = v1_t(:) + 1._dp/sqrt(dx1(:)*dx1(:) + alpha2(:))
                    v2_t(:) = v2_t(:) + 1._dp/sqrt(dx2(:)*dx2(:) + alpha2(:))
                end do
                v1e(:) = -zeff(:) * v1_t(:) * dt / (2 * pi )
                v2e(:) = -zeff(:) * v2_t(:) * dt / (2 * pi )
                pot(:,J) = v12(:) + v1e(:) + v2e(:)
            end do
 
        case default 
            do J = 1, Nx      
                v12(:) = 1._dp/abs(R(:))
                v1e(:) = -zeff(:) / sqrt((x(J) - mn1*R(:))**2 + alpha2(:))  !proton core
                v2e(:) = -zeff(:) / sqrt((x(J) + mn2*R(:))**2 + alpha2(:))  !proton core  
                pot(:,J) = v12(:) + v1e(:) + v2e(:) 
            end do
        end select

        GM = minval(Pot)
        print*, "2D potential global minimum (eV):", GM * au2eV

        pot_file = adjustl(trim(output_data_dir)) // 'gesamtpotential.out'
        ionic_file = adjustl(trim(output_data_dir)) // 'ionic_pot.out'

        open(newunit=pot_unit, file=pot_file, status='replace', form='formatted')
        do I = 1, NR
            do J = 1, Nx
                if (mod(I,4) == 0 .and. mod(J,8) == 0) then
                    write(pot_unit,*) R(I) , x(J), Pot(I,J) * au2eV
                end if
            end do
            if (mod(I,4) == 0) write(pot_unit,*)
        end do
        close(pot_unit)

        open(newunit=ionic_unit, file=ionic_file, status='replace', form='formatted')
        do I = 1, NR
            write(ionic_unit,*) R(I), (1._dp / abs(R(I))) * au2eV
        end do
        close(ionic_unit)

        deallocate(v12, v1e, v2e)
    end subroutine build_2d_potential

    subroutine build_kh_potential_at_time(pot_KH, alpha_val)
        real(dp), intent(out) :: pot_KH(NR, Nx)
        real(dp), intent(in)  :: alpha_val
        integer :: I, J
        real(dp), allocatable :: v12(:), v1e(:), v2e(:)
        real(dp), allocatable :: dx1(:), dx2(:)

        allocate(v12(NR), v1e(NR), v2e(NR))
        allocate(dx1(NR), dx2(NR))

        ! Nuclear repulsion (unchanged in KH frame)
        v12(:) = 1._dp / abs(R(:))

        do J = 1, Nx
            ! Electronic coordinate shifted by quiver displacement alpha_val
            dx1(:) = (x(J) + alpha_val) - mn1 * R(:)
            dx2(:) = (x(J) + alpha_val) + mn2 * R(:)
            v1e(:) = -zeff(:) / sqrt(dx1(:) * dx1(:) + alpha2(:))
            v2e(:) = -zeff(:) / sqrt(dx2(:) * dx2(:) + alpha2(:))
            pot_KH(:, J) = v12(:) + v1e(:) + v2e(:)
        end do

        deallocate(v12, v1e, v2e, dx1, dx2)
    end subroutine build_kh_potential_at_time
end module setpot_mod
        
            

