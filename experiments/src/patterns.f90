module patterns
    use iso_fortran_env, only: real64
    implicit none

contains

    subroutine build_patterns(A, agg, N0, N1)
        real(real64), intent(in) :: A(:,:)
        integer, intent(in) :: agg(:)
        logical, intent(in) :: N0(:,:)
        logical, allocatable, intent(out) :: N1(:,:)
        integer :: i, k

        N1 = N0
        
        ! Now I have only considered ending in full one neighbor pattern, i.e. after 1 time smoothing
        do i = 1, size(agg)
            do k = 1, size(agg)
                ! allow the ith row to connect to their neighbors' aggregates
                if (abs(A(i,k)) > 0.0_real64) then
                    N1(i,agg(k)) = .true.
                end if
            end do
        end do

    end subroutine build_patterns

    subroutine build_ranked_pattern(rows, cols, order, budget, N_trial)
        integer, intent(in) :: rows(:), cols(:), order(:), budget

        logical, intent(out) :: N_trial(:,:)

        ! integer :: i
        
        if (budget > 0) then
            N_trial(rows(order(budget)), cols(order(budget))) = .true.
            ! if (budget == 4) then
            !     do i = 1, budget
            !         print *, rows(order(i)), cols(order(i))
            !     end do
            ! end if
        end if

    end subroutine build_ranked_pattern

end module patterns