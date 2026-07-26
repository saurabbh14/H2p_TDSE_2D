module differentiation
    use global_vars, only: dp
    implicit none

contains

    subroutine central_diff_on_grid(F, Nz, dz, d, points)
        integer K, n_pts
        integer s_loop, mid_pt 
        integer, intent(in) :: Nz
        real(dp), intent(in) :: F(Nz), dz
        real(dp), intent(out) :: d(Nz)
        integer, intent(in), optional :: points
        double precision, allocatable :: central_stencil(:), forward_stencil(:)

        ! Determine stencil size; default to 3-point
        if (present(points)) then
            n_pts = points
        else
            n_pts = 3
        end if

        ! Validate: only odd orders 3,5,7,9 are supported
        if (n_pts /= 3 .and. n_pts /= 5 .and. n_pts /= 7 .and. n_pts /= 9) then
            error stop "differentiation: points must be 3, 5, 7, or 9"
        end if

        allocate(central_stencil(n_pts), forward_stencil(n_pts))

        select case (n_pts)
        case (3)
            central_stencil = (/ -1._dp/2,      0._dp,        1._dp/2     /)
            forward_stencil = (/ -3._dp/2,      2._dp,       -1._dp/2     /)
        case (5)
            central_stencil = (/  1._dp/12,    -2._dp/3,      0._dp,       2._dp/3,    -1._dp/12  /)
            forward_stencil = (/-25._dp/12,     4._dp,       -3._dp,       4._dp/3,    -1._dp/4   /)
        case (7)
            central_stencil = (/ -1._dp/60,     3._dp/20,    -3._dp/4,     0._dp,       3._dp/4,   &
                                -3._dp/20,      1._dp/60  /)
            forward_stencil = (/-49._dp/20,     6._dp,      -15._dp/2,   20._dp/3,   -15._dp/4,  &
                                 6._dp/5,      -1._dp/6   /)
        case (9)
            central_stencil = (/  1._dp/280,   -4._dp/105,   1._dp/5,    -4._dp/5,     0._dp,     &
                                 4._dp/5,      -1._dp/5,     4._dp/105,  -1._dp/280 /)
            forward_stencil = (/-761._dp/280,   8._dp,     -14._dp,      56._dp/3,  -35._dp/2,   &
                                56._dp/5,     -14._dp/3,    8._dp/7,    -1._dp/8   /)
        end select

        d = 0.d0
        mid_pt = n_pts / 2 + 1

        ! Central difference (interior points)
        do K = n_pts/2 + 1, Nz - n_pts/2
            do s_loop = 1, n_pts
                d(K) = d(K) + central_stencil(s_loop) * F(K - mid_pt + s_loop)
            end do
        end do

        ! Forward difference (left boundary)
        do K = 1, n_pts/2
            do s_loop = 1, n_pts
                d(K) = d(K) + forward_stencil(s_loop) * F(K + s_loop - 1)
            end do
        end do

        ! Backward difference (right boundary)
        do K = Nz - n_pts/2 + 1, Nz
            do s_loop = 1, n_pts
                d(K) = d(K) - forward_stencil(s_loop) * F(K - s_loop + 1)
            end do
        end do

        d(:) = d(:) / dz

    end subroutine

end module differentiation