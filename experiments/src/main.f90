program main
    use iso_fortran_env, only: real64
    use patterns, only: build_patterns, build_ranked_pattern
    use minimization, only: algorithm3_iterate
    use connection_score, only: compute_connection_score
    use energy_score, only: compute_energy_score
    use coupling_strength, only: compute_coupling_strength
    use candidate_ranking, only: rank_candidate
    use scalar_1d, only: setup_scalar_1d_case
    use scalar_2d, only: setup_scalar_2d_case
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

    integer :: result_unit, score_unit, rank, k
    integer, parameter :: exhaustive_limit = 12
    character(len=*), parameter :: run_name = 'rail371_coupling'
    
    ! call setup_scalar_1d_case(A, agg, f, P0, N0)
    call setup_scalar_2d_case(A, agg, f, P0, N0)

    call build_patterns(A, agg, N0, N1) ! build unsmoothed pattern N0 and smoothed pattern N1

    allocate(diag_D(size(A,1))); do i = 1, size(A,1); diag_D(i) = A(i,i); enddo
    allocate(C(size(A,1), size(P0,2)))

    C = N1 .and. (.not. N0) ! build candidates C
    m = count(C)
    allocate(rows(m), cols(m), score(m), order(m))
    ! call compute_connection_score(A, P0, C, agg, rows, cols, score)  ! compute the connection score for each candidates and record their positions
    ! call compute_energy_score(A, P0, C, agg, rows, cols, score)
    call compute_coupling_strength(A, C, agg, rows, cols, score)
    call rank_candidate(rows, cols, score, order)   ! rank the scores of all candidates and give them an order

    ! write out the candidates and their scores
    open(newunit=score_unit, file='experiments/result/' // run_name // '_scores.txt', status='replace', action='write')
    write(score_unit, '(A)') 'rank candidate_id row col score'

    do rank = 1, m
        k = order(rank)
        write(score_unit, '(I0,3(1X,I0),1X,ES24.16E3)') rank, k, rows(k), cols(k), score(k)
    end do

    close(score_unit)

    ! compute ranked energy under different budgets
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

    if (m <= exhaustive_limit) then
        print *, "--------------Exhaustive reference--------------"
        call compute_exhaustive_reference(A, P0, N0, f, rows , cols, ref,&
                                        omega=omega, tol=tol, max_steps=max_steps, diag_D=diag_D)
        print *, "failed solves =", ref%failed_count

        print "(A)", "budget   ranked_J         best_J            gap"

        do budget = 0, m
            print "(I6,3(1X,F16.12))", budget, ranked_J(budget), ref%best_J(budget), &
                ranked_J(budget) - ref%best_J(budget)
        end do
    else
        open(newunit=result_unit, file='experiments/result/' // run_name // '_energy.txt', status='replace', action='write')

        write(result_unit, '(A,I0)') '# exhaustive skipped, candidates = ', m
        write(result_unit, '(A)') 'budget ranked_J'

        do budget = 0, m
            write(result_unit, '(I6,1X,ES24.16E3)') budget, ranked_J(budget)
        end do
        close(result_unit)
    end if

end program main