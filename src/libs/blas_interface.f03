module blas_interfaces_module
    use VarPrecision, only: wp=>dp, idp
    implicit none 

    interface
        subroutine zgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
            import wp
            character*1, intent(in):: transa
            character*1, intent(in):: transb
            integer, intent(in):: m
            integer, intent(in):: n
            integer, intent(in):: k
            complex(wp), intent(in) :: alpha
            complex(wp), dimension(*), intent(in):: a
            integer, intent(in):: lda
            complex(wp), dimension(*), intent(in):: b
            integer, intent(in):: ldb
            complex(wp), intent(in):: beta
            complex(wp), dimension(*), intent(out):: c
            integer, intent(in):: ldc
        end subroutine zgemm
        
        subroutine dgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
            import wp
            character*1, intent(in):: transa
            character*1, intent(in):: transb
            integer, intent(in):: m
            integer, intent(in):: n
            integer, intent(in):: k
            real(wp), intent(in) :: alpha
            real(wp), dimension(*), intent(in):: a
            integer, intent(in):: lda
            real(wp), dimension(*), intent(in):: b
            integer, intent(in):: ldb
            real(wp), intent(in):: beta
            real(wp), dimension(*), intent(out):: c
            integer, intent(in):: ldc
        end subroutine dgemm
  
        subroutine zgemv (trans, m, n, alpha,a,lda,x,incx,beta,y,incy)
            import wp
            character*1, intent(in):: trans
            integer, intent(in):: m
            integer, intent(in):: n
            complex(wp), intent(in)::  alpha
            complex(wp), dimension(*), intent(in):: a
            integer, intent(in):: lda
            complex(wp), dimension(*), intent(in):: x
            integer, intent(in):: incx
            complex(wp), intent(in):: beta
            complex(wp), dimension(*), intent(out):: y
            integer, intent(in):: incy
        end subroutine zgemv
  
        subroutine dgemv (trans, m, n, alpha,a,lda,x,incx,beta,y,incy)
            import wp
            character*1, intent(in):: trans
            integer, intent(in):: m
            integer, intent(in):: n
            real(wp), intent(in)::  alpha
            real(wp), dimension(*), intent(in):: a
            integer, intent(in):: lda
            real(wp), dimension(*), intent(in):: x
            integer, intent(in):: incx
            real(wp), intent(in):: beta
            real(wp), dimension(*), intent(out):: y
            integer, intent(in):: incy
        end subroutine dgemv
  
        subroutine dsyev ( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO) 
            import wp
            character*1:: JOBZ
            character*1:: UPLO
            integer:: N
            real(wp), dimension(*):: A
            integer:: LDA
            real(wp), dimension(*):: W
            real(wp), dimension(*):: WORK
            integer:: LWORK
            integer:: INFO
        end subroutine dsyev
    end interface
 
contains

    subroutine blas_check
        use VarPrecision, only: wp=>dp
        real(wp):: A(2,3), x(3), y(2)
 
        A = reshape((/1._wp,0._wp,-1._wp,-3._wp,2._wp,1._wp/),(/2,3/))
        call write_matrix(A)
 
        x = (/2._wp,1._wp,0._wp/)
        write(*,*)
        y = matmul(A,x)
        write(*,*) y
        write(*,*)
        y =0._wp
        call dgemv('N', 2, 3, 1._wp, A, size(A,dim=1), x, 1, 0._wp, y, 1)
        write(*,*) y
    end subroutine blas_check

    subroutine write_matrix(a)
        use VarPrecision, only: wp=>dp
        integer :: i, j
        real(wp), dimension(:,:) :: a
        write(*,*)

        do i = lbound(a,1), ubound(a,1)
            write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
        end do
    end subroutine write_matrix

end module blas_interfaces_module