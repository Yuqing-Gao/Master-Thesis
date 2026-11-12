module scalar_2d
    use iso_fortran_env, only: real64
    implicit none

contains

    subroutine read_mtx_dense(filename, matrix)
        character(len=*), intent(in) :: filename
        real(real64), allocatable, intent(out) :: matrix(:,:)

        character(len=256) :: line
        integer :: nrows, ncols, nentries, i, j, k, unit, ios
        real(real64) :: value

        open(newunit=unit, file=filename, status='old', action='read')
        
        ! read dimension and number of non-zero entries
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) error stop 'could not read matrix dimension'

            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '%') cycle
            exit
        end do

        ! line = 371 371 2343
        read(line, *) nrows, ncols, nentries
        allocate(matrix(nrows, ncols))
        matrix = 0.0_real64

        do k = 1, nentries
            read(unit, *) i, j, value
            matrix(i,j) = matrix(i,j) + value
        end do

        close(unit)

    end subroutine read_mtx_dense

    subroutine setup_scalar_2d_case(A, agg, f, P0, N0)
        real(real64), allocatable, intent(out) :: A(:,:), P0(:,:)
        logical, allocatable, intent(out) :: N0(:,:)
        integer, allocatable, intent(out) :: agg(:), f(:)

        real(real64), allocatable :: S(:,:), M_gamma(:,:)
        integer, parameter :: max_agg_size = 4
        integer :: n, n_agg, seed, member, neighbor, best_neighbor, i
        real(real64) :: weight, best_weight
        
        call read_mtx_dense('experiments/data/fenics_rail_371/ODE_371_S.mtx', S)
        call read_mtx_dense('experiments/data/fenics_rail_371/ODE_371_M_GAMMA.mtx', M_gamma)
        A = S + M_gamma

        n = size(A, 1)
        allocate(agg(n), f(n))
        agg = 0
        f = 1
        n_agg = 0

        do seed = 1, n
            if (agg(seed) /= 0) cycle

            n_agg = n_agg + 1
            agg(seed) = n_agg

            do member = 2, max_agg_size
                best_neighbor = 0
                best_weight = 0.0_real64

                do neighbor = 1, n
                    if (agg(neighbor) /= 0) cycle
                    weight = abs(A(seed, neighbor))
                    if (weight > best_weight) then
                        best_weight = weight
                        best_neighbor = neighbor
                    end if
                end do

                if (best_neighbor == 0) exit
                agg(best_neighbor) = n_agg
            end do
        end do

        allocate(P0(n, n_agg))
        P0 = 0.0_real64
        do i = 1, n
            P0(i, agg(i)) = 1.0_real64
        end do

        N0 = (P0 > 0.0_real64)


    end subroutine setup_scalar_2d_case

end module scalar_2d