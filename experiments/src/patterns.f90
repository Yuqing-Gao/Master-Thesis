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
        
        do i = 1, size(agg)
            do k = 1, size(agg)
                ! allow the ith row to connect to their neighbors' aggregates
                if (abs(A(i,k)) > 0.0_real64) then
                    N1(i,agg(k)) = .true.
                end if
            end do
        end do

    end subroutine build_patterns


end module patterns