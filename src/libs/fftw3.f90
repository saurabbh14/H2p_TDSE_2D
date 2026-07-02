module FFTW3
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
  !  include '/usr/include/fftw3.f03'                                        ! Desktop packet
  !  include '/home/me23jok/ProjectX/FFTW3/include/fftw3.f03' ! ARA cluster
  !  include '/usr/local/include/fftw3.f03'
  !  include '/scratch/Saurabh/FFTW3/install/include/fftw3.f03'  ! nias
  !  include '/home/me23jok/fftw-3.3.10/include/fftw3.f03' ! draco



contains

    subroutine fftw_initialize_threads
        integer(C_INT) :: void

        void = fftw_init_threads( )
        if (void == 0) then
            write(*,*) 'Error in fftw_init_threads, quitting'
            stop
        else
            print*, 'fftw threads initiate successfully'
        endif
    end subroutine fftw_initialize_threads

    subroutine fftw_create_r2r_plans(psi_in, psi_out, NR, planF, planB, parallel)
        use omp_lib
        real(C_DOUBLE), intent(inout) :: psi_in(:), psi_out(:)
        integer, intent(in) :: NR
        type(C_PTR), intent(inout) :: planF, planB
        character(len=10), intent(in) :: parallel
        select case(parallel)
        case("parallel")
            call fftw_plan_with_nthreads(omp_get_max_threads())
            planF = fftw_plan_r2r_1d(NR, psi_in, psi_out, FFTW_R2HC, FFTW_MEASURE)
            planB = fftw_plan_r2r_1d(NR, psi_in, psi_out, FFTW_HC2R, FFTW_MEASURE)
            print*, "Created fftw plans with ", omp_get_max_threads(), ' threads.'
        case default
            planF = fftw_plan_r2r_1d(NR, psi_in, psi_out, FFTW_R2HC, FFTW_MEASURE)
            planB = fftw_plan_r2r_1d(NR, psi_in, psi_out, FFTW_HC2R, FFTW_MEASURE)
        end select
    end subroutine fftw_create_r2r_plans

    subroutine fftw_create_c2c_plans(psi_in, psi_out, NR, planF, planB, parallel)
        use omp_lib
        complex(C_DOUBLE_COMPLEX), intent(inout) :: psi_in(:), psi_out(:)
        integer, intent(in) :: NR
        type(C_PTR), intent(inout) :: planF, planB
        character(len=10), intent(in) :: parallel
        select case(parallel)
        case("parallel")
            call fftw_plan_with_nthreads(omp_get_max_threads())
            planF = fftw_plan_dft_1d(NR, psi_in, psi_out, FFTW_FORWARD, FFTW_MEASURE)
            planB = fftw_plan_dft_1d(NR, psi_in, psi_out, FFTW_BACKWARD, FFTW_MEASURE)
            print*, "Created fftw plans with ", omp_get_max_threads(), ' threads.'
        case default
            planF = fftw_plan_dft_1d(NR, psi_in, psi_out, FFTW_FORWARD, FFTW_MEASURE)
            planB = fftw_plan_dft_1d(NR, psi_in, psi_out, FFTW_BACKWARD, FFTW_MEASURE)
        end select
    end subroutine fftw_create_c2c_plans

    subroutine fftw_create_c2c_2d_plans(psi_in, psi_out, NR, Nx, planF, planB, parallel)
        use omp_lib
        complex(C_DOUBLE_COMPLEX), intent(inout) :: psi_in(:,:), psi_out(:,:)
        integer, intent(in) :: NR, Nx
        type(C_PTR), intent(inout) :: planF, planB
        character(len=10), intent(in) :: parallel
        select case(parallel)
        case("parallel")
            call fftw_plan_with_nthreads(omp_get_max_threads()/2)
            planF = fftw_plan_dft_2d(Nx, NR, psi_in, psi_out, FFTW_FORWARD, FFTW_MEASURE)
            planB = fftw_plan_dft_2d(Nx, NR, psi_in, psi_out, FFTW_BACKWARD, FFTW_MEASURE)
            print*, "Created fftw plans with ", omp_get_max_threads()/2, ' threads.'
        case default
            planF = fftw_plan_dft_2d(Nx, NR, psi_in, psi_out, FFTW_FORWARD, FFTW_MEASURE)
            planB = fftw_plan_dft_2d(Nx, NR, psi_in, psi_out, FFTW_BACKWARD, FFTW_MEASURE)
        end select
    end subroutine fftw_create_c2c_2d_plans
end module