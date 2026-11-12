module coupling_strength
    use iso_fortran_env, only: real64
    implicit none

contains

    subroutine compute_coupling_strength(A, C, agg, rows, cols, score)
        real(real64), intent(in) :: A(:,:)
        logical, intent(in) :: C(:,:)
        integer, intent(in) :: agg(:)
        
        integer, intent(out) :: rows(:), cols(:)
        real(real64), intent(out) :: score(:)

        integer :: neighbor, i, j, k
        real(real64) :: row_max, agg_max

        k = 0
        score = 0.0_real64

        do i = 1, size(C,1)
            ! compute max_{k\neq i}(-a_ik)
            row_max = 0.0_real64
            do neighbor = 1, size(A,2)
                if (neighbor == i .or. A(i,neighbor) >= 0.0_real64) cycle
                if (- A(i,neighbor) >= row_max) then
                    row_max = - A(i,neighbor)
                end if
            end do

            do j = 1, size(C,2)
                if (.not. C(i,j)) cycle
                
                k = k+1
                rows(k) = i
                cols(k) = j

                ! compute max_{k\in A_j}c_ik
                agg_max = 0.0_real64
                do neighbor = 1, size(A,2)
                    if (neighbor == i .or. A(i,neighbor) >= 0.0_real64) cycle
                    if (agg(neighbor) /= j) cycle
                    if (- A(i,neighbor) >= agg_max) then
                        agg_max = - A(i,neighbor)
                    end if
                end do
                
                if (row_max > 0.0_real64) then
                    score(k) = agg_max / row_max
                end if
                ! print *, "position:", rows(k), cols(k), "score = ", score(k)
            end do
        end do
                
    end subroutine compute_coupling_strength

end module coupling_strength