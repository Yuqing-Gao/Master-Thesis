program main
    use iso_fortran_env, only: real64
    use matrix_ops, only: compute_energy, apply_Dinv
    use patterns, only: build_patterns
    use constraints, only: project_direction
    use minimization, only: algorithm3_step, algorithm3_iterate
    use candidate_ranking, only: rank_candidate
    use scalar_1d, only: setup_scalar_1d_case, build_diffusion_1d
    implicit none
    

    integer, allocatable :: agg(:), f(:)
    real(real64), allocatable :: A(:,:), P0(:,:), P_new(:,:), diag_D(:)
    logical, allocatable :: N0(:,:), N1(:,:)
    real(real64) :: omega, tol
    integer :: i, j, k, max_steps
    logical :: converged

    real(real64), allocatable :: score(:), G0(:,:)
    integer, allocatable :: rows(:), cols(:), order(:)
    logical, allocatable :: C(:,:), N_trial(:,:)
    
    call setup_scalar_1d_case(A, agg, f, P0, N0)

    print *, compute_energy(A, P0)

    call build_patterns(A, agg, N0, N1)

    allocate(diag_D(size(A,1))); do i = 1, size(A,1); diag_D(i) = A(i,i); enddo
    
    G0 = matmul (A, P0)

    allocate(C(size(G0,1), size(G0,2)))

    C = N1 .and. (.not. N0)

    allocate(rows(count(C)), cols(count(C)), score(count(C)), order(count(C)))
    score = 0.0_real64
    k = 0
    
    ! TODO: use pack instead of loop?
    do i = 1, size(G0,1)
        do j = 1, size(G0,2)
            if (C(i,j)) then
                k = k + 1
                score(k) = abs(G0(i,j)) / A(i,i)
                rows(k) = i
                cols(k) = j
                print *, "position:", rows(k), cols(k), "score = ", score(k)
            endif
        enddo
    enddo


    call rank_candidate(rows, cols, score, order)
    print *, order

    N_trial = N0
    N_trial(rows(order(1)), cols(order(1))) = .true.
    print *, count(N_trial)


    omega = 0.5_real64
    tol = 1e-8
    max_steps = 1e3

    call algorithm3_iterate(A, P0, N0, f, diag_D, omega, tol, max_steps, P_new, converged)

    call algorithm3_iterate(A, P0, N_trial, f, diag_D, omega, tol, max_steps, P_new, converged)


end program main