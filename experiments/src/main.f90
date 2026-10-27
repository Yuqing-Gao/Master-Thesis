program main
    use iso_fortran_env, only: real64
    use matrix_ops, only: compute_energy, apply_Dinv
    use patterns, only: build_patterns
    use constraints, only: project_direction
    use minimization, only: algorithm3_step, algorithm3_iterate
    implicit none
    
    integer :: n, n_agg  ! n: the size of A; n_agg: the number of aggregations
    integer :: i, n_steps
    real(real64) :: energy, omega
    integer, allocatable :: agg(:), f(:)
    logical, allocatable :: N0(:,:), N1(:,:)
    real(real64), allocatable :: A(:,:), P_hat(:,:), H(:,:), diag_D(:), V(:,:), P_new(:,:)

    n = 4
    n_agg = 2
    n_steps = 10
    omega = 0.5

    allocate(A(n,n), agg(n), P_hat(n,n_agg), H(n,n_agg), V(n,n_agg), P_new(n,n_agg))
    allocate(N0(n,n_agg), N1(n,n_agg))
    allocate(f(n), diag_D(n))

    A = 0.0_real64
    H = 0.0_real64
    V = 0.0_real64
    P_hat = 0.0_real64
    agg = [1, 1, 2, 2]
    f = [1, 1, 1, 1]
    diag_D = 0.0_real64

    ! constrct matrix A
    do i = 1, n
        A(i, i) = 2.0_real64
        if (i /= 1) then
            A(i, i-1) = -1.0_real64
        end if
        if (i /= n) then
            A(i, i+1) = -1.0_real64
        end if
    end do


    ! construct the nonsmoothed interpolation P_hat
    do i = 1, n
        P_hat(i, agg(i)) = 1
    end do

    ! compute the energy of P by J(P) = \sum\sum Pij(AP)ij
    energy = compute_energy(A, P_hat)
    print *, "Energy =", energy

    ! construct diag_D
    do i = 1, n
        diag_D(i) = A(i,i)
    end do

    call build_patterns(A, agg, N0, N1)

    call algorithm3_iterate(A, P_hat, N1, f, diag_D, omega, n_steps, P_new)


end program main