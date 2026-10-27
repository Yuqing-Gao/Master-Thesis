module minimization
    use iso_fortran_env, only: real64
    use matrix_ops, only: compute_energy, apply_Dinv
    use constraints, only: project_direction
    implicit none

contains

    subroutine algorithm3_step(A, P, N, f, diag_D, omega, P_next)
        real(real64), intent(in) :: A(:,:), P(:,:), diag_D(:), omega
        integer, intent(in) :: f(:)
        logical, intent(in) :: N(:,:)

        real(real64), intent(out) :: P_next(:,:)

        real(real64) :: G(size(A,1),size(P,2)), V(size(A,1),size(P,2)), H(size(A,1),size(P,2))

        G = matmul(A, P)
        where (.not. N) G = 0.0_real64  ! G = N \odot AP
        call apply_Dinv(diag_D, G, V)   ! V = D^-1 G
        call project_direction(V, N, f, H)  ! H = Z(V)
        P_next = P - omega * H

    end subroutine algorithm3_step

    subroutine algorithm3_iterate(A, P0, N, f, diag_D, omega, n_steps, P_final)
        real(real64), intent(in) :: A(:,:), P0(:,:), diag_D(:), omega
        integer, intent(in) :: f(:), n_steps
        logical, intent(in) :: N(:,:)

        real(real64), intent(out) :: P_final(:,:)

        real(real64), allocatable :: P_current(:,:)
        real(real64) :: energy
        integer :: iter

        P_current = P0

        do iter = 1, n_steps
            call algorithm3_step(A, P_current, N, f, diag_D, omega, P_final)
            P_current = P_final
            energy = compute_energy(A, P_current)
            print *, "step = ", iter, "energy =", energy
        end do

    end subroutine algorithm3_iterate

end module minimization