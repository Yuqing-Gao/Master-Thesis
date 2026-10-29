module candidate_ranking
    use iso_fortran_env, only: real64
    implicit none

contains
    
    subroutine rank_candidate(rows, cols, score, order)
        integer, intent(in) :: rows(:), cols(:)
        real(real64), intent(in) :: score(:)

        integer, allocatable, intent(out) :: order(:)

        integer :: m, i, j, best

        m = size(score)

        order = [(i, i = 1, m)]

        ! selection sort
        do i = 1, m-1
            best = i
            do j = i+1, m
                if (precedes(order(j), order(best))) best = j
            enddo
            if (best /= i) order([i, best]) = order ([best, i])
        enddo

        contains
        
        logical function precedes(a, b)     ! if a would precede b
            integer, intent(in) :: a, b

            if (score(a) > score(b)) then
                precedes = .true.
            else if (score(a) < score(b)) then
                precedes = .false.
            else if (rows(a) /= rows(b)) then
                precedes = rows(a) < rows(b)
            else
                precedes = cols(a) < cols(b)
            end if
        end function

    end subroutine rank_candidate

    

end module candidate_ranking