module differentiation
    use global_vars, only: dp
    implicit none

contains

    subroutine central_diff_on_grid(F, Nz, dz, d)
        integer K, points, points_forward
        integer s_loop, mid_pt 
        integer, intent(in) :: Nz
        real(dp), intent(in) :: F(Nz), dz
        real(dp), intent(out) :: d(Nz)
        double precision, allocatable:: central_stencil(:), forward_stencil(:)

        points = 3
        points_forward = 4
        allocate(central_stencil(points), forward_stencil(points_forward))
        central_stencil = (/ -1._dp/2, 0._dp, 1._dp/2 /)
        forward_stencil = (/ -11._dp/6, 3._dp, -3._dp/2, 1._dp/3/)

        d = 0.d0
        mid_pt = points/2 + 1

        ! Central difference
        do K = points/2 +1, Nz - points/2
            !d(K) = (F(K+1) + 0.0 * F(K) - F(K-1)) / (2 * dz) ! 3-point central difference
            do s_loop = 1, points
                d(K) = d(K) + central_stencil(s_loop) * F(K - mid_pt + s_loop)
            enddo
        enddo
        ! Forward difference
        do K = 1, points/2 
            do s_loop = 1, points_forward
                d(K) = d(K) + forward_stencil(s_loop) * F(K + s_loop - 1)
            enddo
            !d(K) = sum(forward_stencil(:) * F(K:K+points))
        enddo   
        ! Backward difference
        do K = Nz - points/2 + 1, Nz 
            do s_loop = 1, points_forward
                d(K) = d(K) - forward_stencil(s_loop) * F(K - s_loop + 1)
            enddo
        enddo   

        d(:) = d(:) / dz

    end subroutine
    
end module differentiation