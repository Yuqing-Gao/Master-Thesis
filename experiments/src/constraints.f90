module constraints
    use iso_fortran_env, only: real64
    implicit none

contains

    subroutine project_direction(G, N, f, H)
        real(real64), intent(in) :: G(:,:) ! G: the gradient AP
        logical, intent(in) :: N(:,:) ! N: the allowed pattern
        integer, intent(in) :: f(:) ! f: the indicator of Dirichlet boundary
        
        real(real64), intent(out) :: H(:,:) ! H: the output project direction

        integer :: i, j

        do i = 1, size(G, 1)
            do j = 1, size(G, 2)
                if (N(i,j) .and. f(i)==0) then
                    H(i,j) = G(i,j)
                else if (N(i,j) .and. f(i)==1) then
                    H(i,j) = G(i,j) - sum(G(i,:), mask=N(i,:)) / count(N(i,:))
                else
                    H(i,j) = 0
                end if
            end do
        end do


    end subroutine project_direction
    

end module constraints