module minimization
    use iso_fortran_env, only: real64
    use matrix_ops, only: compute_energy, apply_Dinv
    use constraints, only: project_direction
    implicit none

contains

    subroutine algorithm3_step(A, P, N, f, diag_D, omega, P_next, projected_gradient)
        real(real64), intent(in) :: A(:,:), P(:,:), diag_D(:), omega
        integer, intent(in) :: f(:)
        logical, intent(in) :: N(:,:)

        real(real64), intent(out) :: P_next(:,:), projected_gradient

        real(real64) :: G(size(A,1),size(P,2)), V(size(A,1),size(P,2)), H(size(A,1),size(P,2))

        G = matmul(A, P)
        where (.not. N) G = 0.0_real64  ! G = N \odot AP
        call apply_Dinv(diag_D, G, V)   ! V = D^-1 G
        call project_direction(V, N, f, H)  ! H = Z(V)
        projected_gradient = sqrt(sum(H**2))    ! ||H||_F of the current P
        P_next = P - omega * H

    end subroutine algorithm3_step

    subroutine algorithm3_iterate(A, P0, N, f, diag_D, omega, tol, max_steps, P_final, converged)
        real(real64), intent(in) :: A(:,:), P0(:,:), diag_D(:), omega, tol
        integer, intent(in) :: f(:), max_steps
        logical, intent(in) :: N(:,:)

        real(real64), allocatable, intent(out) :: P_final(:,:)
        logical, intent(out) :: converged

        real(real64), allocatable :: P_current(:,:), P_next(:,:)
        real(real64) :: energy, current_projected_gradient
        integer :: iter

        P_current = P0
        P_next = P0
        converged = .false.

        do iter = 0, max_steps
            ! energy = compute_energy(A, P_current)
            call algorithm3_step(A, P_current, N, f, diag_D, omega, P_next, current_projected_gradient)
            ! print *, current_projected_gradient

            if (current_projected_gradient <= tol) then
                converged = .true.
                exit
            end if
            if (iter == max_steps) exit

            P_current = P_next
        end do

        P_final = P_current
        energy = compute_energy(A, P_final)

        if (converged) then
            print *, "converged at iter", iter, "with energy", energy, "and projected gradient", current_projected_gradient
        else
            print *, "max_steps reached without convergence; energy", energy, "projected gradient", current_projected_gradient
        end if

    end subroutine algorithm3_iterate

end module minimization