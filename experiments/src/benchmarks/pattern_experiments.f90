module pattern_experiments
    use iso_fortran_env, only: real64
    use minimization, only: algorithm3_iterate
    implicit none

    type :: exhaustive_result
        real(real64), allocatable :: best_J(:)  ! best energy under a certain budget
        integer, allocatable :: best_mask(:)    ! a mask that can reach teh best energy
        integer :: failed_count = 0 ! count of cases that not converged
    end type exhaustive_result

contains

    ! exhaustively calculate the best reference for the current case
    subroutine compute_exhaustive_reference(A, P0, N0, f, rows , cols, ref,&
                                            omega, tol, max_steps, diag_D)
        type(exhaustive_result) :: ref
        real(real64), intent(in) :: A(:,:), P0(:,:), diag_D(:)
        integer, intent(in) :: rows(:), cols(:)
        logical, intent(in) :: N0(:,:)
        integer, intent(in) :: f(:), max_steps
        real(real64), intent(in) :: omega, tol
        
        
        real(real64), allocatable :: P_new(:,:)
        real(real64) :: energy
        logical, allocatable :: N_trial(:,:)
        logical :: converged
        integer :: m, k, mask, budget

        m = size(rows)

        allocate(ref%best_J(0:m), ref%best_mask(0:m))

        ref%best_J = huge(1.0_real64)

        do mask = 0 ,2**m-1
            N_trial = N0

            do k = 1, m
                if (btest(mask,k-1)) then
                    N_trial(rows(k), cols(k)) = .true.
                end if
            end do

            call algorithm3_iterate(A, P0, N_trial, f, diag_D, omega, tol, max_steps, P_new, converged, energy)

            budget = popcnt(mask)

            if (converged) then
                ref%best_J(budget) = min(ref%best_J(budget), energy)
            else
                ref%failed_count = ref%failed_count + 1
            end if
        end do

    end subroutine compute_exhaustive_reference

end module pattern_experiments