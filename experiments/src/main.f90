program main
    use iso_fortran_env, only: real64
    use matrix_ops, only: compute_energy
    use patterns, only: build_patterns, build_ranked_pattern
    use minimization, only: algorithm3_iterate
    use connection_score, only: compute_connection_score
    use energy_score, only: compute_energy_score
    use candidate_ranking, only: rank_candidate
    use scalar_1d, only: setup_scalar_1d_case
    use pattern_experiments, only: exhaustive_result, compute_exhaustive_reference
    implicit none
    

    integer, allocatable :: agg(:), f(:)
    real(real64), allocatable :: A(:,:), P0(:,:), P_new(:,:), diag_D(:)
    logical, allocatable :: N0(:,:), N1(:,:)
    real(real64) :: omega, tol
    integer :: i, max_steps
    logical :: converged
    type(exhaustive_result) :: ref

    real(real64), allocatable :: score(:)
    integer :: budget, m
    real(real64) :: energy
    real(real64), allocatable :: ranked_J(:)
    integer, allocatable :: rows(:), cols(:), order(:)
    logical, allocatable :: C(:,:), N_trial(:,:)
    
    call setup_scalar_1d_case(A, agg, f, P0, N0)

    print *, "Energy of unsmoothed interpolation P0 = ", compute_energy(A, P0)

    call build_patterns(A, agg, N0, N1) ! build unsmoothed pattern N0 and smoothed pattern N1

    allocate(diag_D(size(A,1))); do i = 1, size(A,1); diag_D(i) = A(i,i); enddo
    allocate(C(size(A,1), size(P0,2)))

    C = N1 .and. (.not. N0) ! build candidates C
    m = count(C)
    allocate(rows(m), cols(m), score(m), order(m))
    ! call compute_connection_score(A, P0, C, agg, rows, cols, score)  ! compute the connection score for each candidates and record their positions
    call compute_energy_score(A, P0, C, agg, rows, cols, score)
    call rank_candidate(rows, cols, score, order)   ! rank the scores of all candidates and give them an order

    omega = 0.5_real64
    tol = 1.0e-8_real64
    max_steps = 1000

    allocate(ranked_J(0:m))

    N_trial = N0
    do budget = 0, m
        call build_ranked_pattern(rows, cols, order, budget, N_trial)
        ! print *, "budget =", budget, "allowed positions =", count(N_trial)
        call algorithm3_iterate(A, P0, N_trial, f, diag_D, omega, tol, max_steps, P_new, converged, energy)
        if (.not. converged) then
            print *, "Ranked solve failed at budget =", budget
            error stop
        end if
        ranked_J(budget) = energy
    end do

    print *, "--------------Exhaustive reference--------------"
    call compute_exhaustive_reference(A, P0, N0, f, rows , cols, ref,&
                                    omega=omega, tol=tol, max_steps=max_steps, diag_D=diag_D)
    print *, "failed solves =", ref%failed_count

    print "(A)", "budget   ranked_J         best_J            gap"

    do budget = 0, m
        print "(I6,3(1X,F16.12))", budget, ranked_J(budget), ref%best_J(budget), &
            ranked_J(budget) - ref%best_J(budget)
    end do
    
    ! connection score 
    ! budget   ranked_J         best_J            gap
    !  0   6.950000000000   6.950000000000   0.000000000000
    !  1   5.350000000000   4.700000000000   0.650000000000
    !  2   3.100000000000   3.100000000000   0.000000000000
    !  3   2.778571428571   2.305395683453   0.473175745118
    !  4   2.664285714286   1.866620853949   0.797664860337
    !  5   1.951064265248   1.616620853949   0.334443411299
    !  6   1.804706513958   1.376866398733   0.427840115225
    !  7   1.538709000441   1.242075163399   0.296633837043
    !  8   1.326834976137   1.135304908869   0.191530067268
    !  9   1.076834976137   1.076834976137   0.000000000000
    ! 10   1.024609515445   1.024609515445   0.000000000000

    ! energy score
    ! budget   ranked_J         best_J            gap
    !  0   6.950000000000   6.950000000000   0.000000000000
    !  1   4.700000000000   4.700000000000   0.000000000000
    !  2   4.378571428571   3.100000000000   1.278571428571 ! only when budget = 1 & 2, the 2 scores lead to different gap
    !  3   2.778571428571   2.305395683453   0.473175745118
    !  4   2.664285714286   1.866620853949   0.797664860337
    !  5   1.951064265248   1.616620853949   0.334443411299
    !  6   1.804706513958   1.376866398733   0.427840115225
    !  7   1.538709000441   1.242075163399   0.296633837043
    !  8   1.326834976137   1.135304908869   0.191530067268
    !  9   1.076834976137   1.076834976137   0.000000000000
    ! 10   1.024609515445   1.024609515445   0.000000000000

end program main