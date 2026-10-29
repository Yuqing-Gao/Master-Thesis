module scalar_1d
    use iso_fortran_env, only: real64
    implicit none

contains

    subroutine setup_scalar_1d_case(A, agg, f, P0, N0)
        real(real64), allocatable, intent(out) :: A(:,:), P0(:,:)
        logical, allocatable, intent(out) :: N0(:,:)
        integer, allocatable, intent(out) :: agg(:), f(:)

        real(real64) :: w(13)
        integer :: n, n_agg, i

        w = [1.0_real64, 1.0_real64, 0.25_real64, 1.0_real64, &
            2.0_real64, 0.5_real64, 0.8_real64, 1.0_real64,  &
            3.0_real64, 1.0_real64, 0.4_real64, 1.0_real64,  &
            1.0_real64]

        n = size(w) - 1     ! number of the interior nodes
        n_agg = n / 2   ! every 2 adjacent nodes form 1 aggregate
        
        allocate(A(n,n), P0(n, n_agg), agg(n), f(n))

        call build_diffusion_1d(w, A)   ! build A

        do i = 1, n_agg; agg(2*i-1) = i; agg(2*i) = i; enddo    ! build aggregate

        f = 1; f(1) = 0; f(n) = 0;  ! build f

        P0 = 0.0_real64
        do i = 1, n; P0(i, agg(i)) = 1; enddo   ! build unsmoothed interpolation P0
        
        N0 = (abs(P0) > 0.0_real64)     ! build N0

    end subroutine setup_scalar_1d_case

    subroutine build_diffusion_1d(w, A)
        real(real64), intent(in) :: w(:)    ! w is the weight of the edges
        real(real64), intent(out) :: A(:,:)

        integer :: n, i

        n = size(w) - 1

        A = 0.0_real64

        do i = 1, n; A(i,i) = w(i) + w(i+1); enddo
        do i = 1, n-1; A(i,i+1) = -w(i+1); A(i+1,i) = -w(i+1); enddo

    end subroutine build_diffusion_1d

end module scalar_1d