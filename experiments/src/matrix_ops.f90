module matrix_ops
    use iso_fortran_env, only: real64
    implicit none

contains

    function compute_energy(A, P) result(value)
        real(real64), intent(in) :: A(:,:), P(:,:)
        real(real64) :: value
        
        value = sum(P * matmul(A, P))
    end function compute_energy

    subroutine apply_Dinv(diag_D, G, V)
        real(real64), intent(in) :: diag_D(:)
        real(real64), intent(in) :: G(:,:) ! G = AP
        real(real64), intent(out) :: V(:,:) ! V = D^{-1}G

        integer :: i, j

        do i = 1, size(G, 1)
            do j = 1, size(G, 2)
                V(i,j) = G(i,j) / diag_D(i)
            end do
        end do
    
    end subroutine apply_Dinv


end module matrix_ops