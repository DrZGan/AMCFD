program tdma_benchmark
    use omp_lib
    implicit none

    ! Grid dimensions - 100x100x100 for benchmark
    integer, parameter :: ni = 100, nj = 100, nk = 100
    integer, parameter :: nim1 = ni - 1, njm1 = nj - 1, nkm1 = nk - 1
    integer, parameter :: n_iterations = 10

    ! Arrays
    real(8), allocatable :: enthalpy(:,:,:)
    real(8), allocatable :: ap(:,:,:), ae(:,:,:), aw(:,:,:)
    real(8), allocatable :: an(:,:,:), as(:,:,:), at(:,:,:), ab(:,:,:)
    real(8), allocatable :: su(:,:,:)

    ! Timing
    real(8) :: start_time, end_time, total_time, avg_time
    real(8) :: checksum, max_val, min_val, avg_val
    integer :: iter, num_threads
    character(len=32) :: arg

    ! Get number of threads from command line
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) num_threads
    else
        num_threads = 1
    end if

    call omp_set_num_threads(num_threads)

    ! Allocate arrays
    allocate(enthalpy(ni, nj, nk))
    allocate(ap(ni, nj, nk), ae(ni, nj, nk), aw(ni, nj, nk))
    allocate(an(ni, nj, nk), as(ni, nj, nk), at(ni, nj, nk), ab(ni, nj, nk))
    allocate(su(ni, nj, nk))

    ! Initialize arrays
    call initialize_arrays()

    ! Warmup run
    call tdma_solve()

    ! Re-initialize for fair comparison
    call initialize_arrays()

    ! Benchmark
    total_time = 0.0d0
    do iter = 1, n_iterations
        start_time = omp_get_wtime()
        call tdma_solve()
        end_time = omp_get_wtime()
        total_time = total_time + (end_time - start_time)
    end do

    avg_time = (total_time / n_iterations) * 1000.0d0  ! Convert to ms
    checksum = sum(enthalpy)
    max_val = maxval(enthalpy)
    min_val = minval(enthalpy)
    avg_val = checksum / (ni * nj * nk)

    ! Output in CSV-friendly format: threads, avg_time_ms, checksum, max, min, avg
    write(*,'(I2,A,F12.6,A,E20.12,A,E20.12,A,E20.12,A,E20.12)') &
        num_threads, ',', avg_time, ',', checksum, ',', max_val, ',', min_val, ',', avg_val

    ! Deallocate
    deallocate(enthalpy, ap, ae, aw, an, as, at, ab, su)

contains

    subroutine initialize_arrays()
        integer :: i, j, k

        do k = 1, nk
            do j = 1, nj
                do i = 1, ni
                    enthalpy(i,j,k) = 1.0d0
                    ap(i,j,k) = 6.0d0
                    ae(i,j,k) = 1.0d0
                    aw(i,j,k) = 1.0d0
                    an(i,j,k) = 1.0d0
                    as(i,j,k) = 1.0d0
                    at(i,j,k) = 1.0d0
                    ab(i,j,k) = 1.0d0
                    su(i,j,k) = 0.1d0
                end do
            end do
        end do
    end subroutine initialize_arrays

    subroutine tdma_solve()
        integer :: i, j, k, ksweep, jsweep
        real(8) :: pr(ni), qr(ni)
        real(8) :: d, denom

        do ksweep = 1, 2
            do k = nkm1, 2, -1
                do jsweep = 1, 2
                    !$OMP PARALLEL DO PRIVATE(pr, qr, d, denom, i) SCHEDULE(STATIC)
                    do j = 2, njm1
                        pr(1) = 0.0d0
                        qr(1) = enthalpy(1,j,k)

                        do i = 2, nim1
                            d = at(i,j,k)*enthalpy(i,j,k+1) + ab(i,j,k)*enthalpy(i,j,k-1) + &
                                an(i,j,k)*enthalpy(i,j+1,k) + as(i,j,k)*enthalpy(i,j-1,k) + su(i,j,k)
                            denom = ap(i,j,k) - aw(i,j,k)*pr(i-1)

                            if (denom <= 1e-12 .and. denom >= 0.0d0) denom = denom + 1e-13
                            if (denom >= -1e-12 .and. denom < 0.0d0) denom = denom - 1e-13

                            pr(i) = ae(i,j,k) / denom
                            qr(i) = (d + aw(i,j,k)*qr(i-1)) / denom
                        end do

                        do i = nim1, 2, -1
                            enthalpy(i,j,k) = pr(i)*enthalpy(i+1,j,k) + qr(i)
                        end do
                    end do
                    !$OMP END PARALLEL DO
                end do
            end do
        end do
    end subroutine tdma_solve

end program tdma_benchmark
