! Explicit Euler with Laplacian - Fortran OpenMP Implementation
! Usage: ./laplacian_fortran <num_threads> <grid_size>
program laplacian_benchmark
    use omp_lib
    implicit none

    integer :: nx, ny, nz
    integer, parameter :: nsteps = 100
    real(8), parameter :: dx = 0.1d0, dt = 0.0001d0, alpha = 0.1d0
    real(8), allocatable :: u(:,:,:), u_new(:,:,:)
    real(8) :: coeff, start_time, end_time, elapsed_time
    real(8) :: result_sum, result_min, result_max, result_avg
    integer :: i, j, k, step, num_threads
    character(len=32) :: arg

    ! Get number of threads from command line
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) num_threads
    else
        num_threads = 1
    end if

    ! Get grid size from command line
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *) nx
    else
        nx = 100
    end if
    ny = nx
    nz = nx

    call omp_set_num_threads(num_threads)

    ! Allocate arrays
    allocate(u(nx, ny, nz))
    allocate(u_new(nx, ny, nz))

    ! Initialize with a test pattern
    coeff = alpha * dt / (dx * dx)

    !$omp parallel do collapse(3) private(i, j, k)
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                u(i, j, k) = sin(real(i, 8) * 0.1d0) * sin(real(j, 8) * 0.1d0) * sin(real(k, 8) * 0.1d0)
            end do
        end do
    end do
    !$omp end parallel do

    u_new = u

    ! Start timing
    start_time = omp_get_wtime()

    ! Time stepping loop
    do step = 1, nsteps
        !$omp parallel do collapse(3) private(i, j, k)
        do k = 2, nz-1
            do j = 2, ny-1
                do i = 2, nx-1
                    u_new(i, j, k) = u(i, j, k) + coeff * ( &
                        u(i+1, j, k) + u(i-1, j, k) + &
                        u(i, j+1, k) + u(i, j-1, k) + &
                        u(i, j, k+1) + u(i, j, k-1) - &
                        6.0d0 * u(i, j, k))
                end do
            end do
        end do
        !$omp end parallel do

        ! Swap arrays
        !$omp parallel do collapse(3) private(i, j, k)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    u(i, j, k) = u_new(i, j, k)
                end do
            end do
        end do
        !$omp end parallel do
    end do

    ! End timing
    end_time = omp_get_wtime()
    elapsed_time = end_time - start_time

    ! Compute verification metrics
    result_sum = 0.0d0
    result_min = u(1, 1, 1)
    result_max = u(1, 1, 1)

    !$omp parallel do collapse(3) private(i, j, k) reduction(+:result_sum) reduction(min:result_min) reduction(max:result_max)
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                result_sum = result_sum + u(i, j, k)
                result_min = min(result_min, u(i, j, k))
                result_max = max(result_max, u(i, j, k))
            end do
        end do
    end do
    !$omp end parallel do

    result_avg = result_sum / real(nx * ny * nz, 8)

    ! Output results in CSV format for easy parsing
    write(*, '(I4,",",E22.14,",",E22.14,",",E22.14,",",E22.14,",",F12.6)') &
        num_threads, result_sum, result_min, result_max, result_avg, elapsed_time

    deallocate(u, u_new)

end program laplacian_benchmark
