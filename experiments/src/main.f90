program main
    use iso_fortran_env, only: real64
    use matrix_ops, only: compute_energy, apply_Dinv
    use patterns, only: build_patterns
    use constraints, only: project_direction
    use minimization, only: algorithm3_step, algorithm3_iterate
    use scalar_1d, only: setup_scalar_1d_case, build_diffusion_1d
    implicit none
    

    integer, allocatable :: agg(:), f(:)
    real(real64), allocatable :: A(:,:), P0(:,:), P_new(:,:), diag_D(:)
    logical, allocatable :: N0(:,:), N1(:,:)
    real(real64) :: omega, tol
    integer :: i, max_steps
    logical :: converged

    call setup_scalar_1d_case(A, agg, f, P0, N0)

    print *, compute_energy(A, P0)

    call build_patterns(A, agg, N0, N1)

    allocate(diag_D(size(A,1)))
    do i = 1, size(A,1); diag_D(i) = A(i,i); enddo

    omega = 0.5_real64
    tol = 1e-8
    max_steps = 1e3

    call algorithm3_iterate(A, P0, N1, f, diag_D, omega, tol, max_steps, P_new, converged)

end program main