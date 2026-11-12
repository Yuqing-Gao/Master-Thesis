module energy_score
    use iso_fortran_env, only: real64
    implicit none

contains

    subroutine compute_energy_score(A, P0, C, agg, rows, cols, score)
        real(real64), intent(in) :: A(:,:), P0(:,:)
        logical, intent(in) :: C(:,:)
        integer, intent(in) :: agg(:)
        
        integer, intent(out) :: rows(:), cols(:)
        real(real64), intent(out) :: score(:)

        real(real64), allocatable :: G0(:,:)
        integer :: i, j, k, coarse_col

        G0 = matmul (A, P0)

        score = 0.0_real64
        k = 0
        
        ! TODO: use pack instead of loop?
        do i = 1, size(G0,1)
            do j = 1, size(G0,2)
                if (C(i,j)) then
                    k = k + 1
                    coarse_col = agg(i)
                    score(k) = (G0(i,coarse_col) - G0(i,j)) ** 2 / (4.0_real64 * A(i,i))
                    rows(k) = i
                    cols(k) = j
                    ! print *, "position:", rows(k), cols(k), "score = ", score(k)
                endif
            enddo
        enddo
    
    end subroutine

end module energy_score