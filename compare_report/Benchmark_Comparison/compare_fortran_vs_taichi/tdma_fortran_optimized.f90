program tdma_benchmark_final
    use omp_lib
    implicit none

    ! Grid dimensions
    integer, parameter :: ni = 100, nj = 100, nk = 100
    integer, parameter :: nim1 = ni - 1, njm1 = nj - 1, nkm1 = nk - 1
    integer, parameter :: n_iterations = 10
    integer, parameter :: MAX_THREADS = 16

    ! Arrays
    real(8), allocatable :: enthalpy(:,:,:)
    real(8), allocatable :: ap(:,:,:), ae(:,:,:), aw(:,:,:)
    real(8), allocatable :: an(:,:,:), as(:,:,:), at(:,:,:), ab(:,:,:)
    real(8), allocatable :: su(:,:,:)

    ! Pre-allocated thread-local arrays (avoid repeated allocation)
    real(8) :: pr(ni, 0:MAX_THREADS-1)
    real(8) :: qr(ni, 0:MAX_THREADS-1)

    real(8) :: start_time, end_time, total_time, avg_time, checksum
    integer :: iter, num_threads
    character(len=32) :: arg

    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) num_threads
    else
        num_threads = 1
    end if

    call omp_set_num_threads(num_threads)

    allocate(enthalpy(ni, nj, nk))
    allocate(ap(ni, nj, nk), ae(ni, nj, nk), aw(ni, nj, nk))
    allocate(an(ni, nj, nk), as(ni, nj, nk), at(ni, nj, nk), ab(ni, nj, nk))
    allocate(su(ni, nj, nk))

    call init()
    call solve()  ! Warmup
    call init()

    total_time = 0.0d0
    do iter = 1, n_iterations
        start_time = omp_get_wtime()
        call solve()
        end_time = omp_get_wtime()
        total_time = total_time + (end_time - start_time)
    end do

    avg_time = (total_time / n_iterations) * 1000.0d0
    checksum = sum(enthalpy)

    write(*,'(I2,A,F12.6,A,E20.12)') num_threads, ',', avg_time, ',', checksum

    deallocate(enthalpy, ap, ae, aw, an, as, at, ab, su)

contains

    subroutine init()
        enthalpy = 1.0d0
        ap = 6.0d0
        ae = 1.0d0
        aw = 1.0d0
        an = 1.0d0
        as = 1.0d0
        at = 1.0d0
        ab = 1.0d0
        su = 0.1d0
    end subroutine init

    subroutine solve()
        integer :: i, j, k, ksweep, jsweep, tid
        real(8) :: d, denom

        do ksweep = 1, 2
            do k = nkm1, 2, -1
                do jsweep = 1, 2
                    !$OMP PARALLEL DO PRIVATE(tid, d, denom, i) SCHEDULE(GUIDED, 4)
                    do j = 2, njm1
                        tid = omp_get_thread_num()

                        pr(1, tid) = 0.0d0
                        qr(1, tid) = enthalpy(1,j,k)

                        do i = 2, nim1
                            d = at(i,j,k)*enthalpy(i,j,k+1) + ab(i,j,k)*enthalpy(i,j,k-1) + &
                                an(i,j,k)*enthalpy(i,j+1,k) + as(i,j,k)*enthalpy(i,j-1,k) + su(i,j,k)
                            denom = ap(i,j,k) - aw(i,j,k)*pr(i-1, tid)

                            if (denom <= 1.0d-12 .and. denom >= 0.0d0) denom = denom + 1.0d-13
                            if (denom >= -1.0d-12 .and. denom < 0.0d0) denom = denom - 1.0d-13

                            pr(i, tid) = ae(i,j,k) / denom
                            qr(i, tid) = (d + aw(i,j,k)*qr(i-1, tid)) / denom
                        end do

                        do i = nim1, 2, -1
                            enthalpy(i,j,k) = pr(i, tid)*enthalpy(i+1,j,k) + qr(i, tid)
                        end do
                    end do
                    !$OMP END PARALLEL DO
                end do
            end do
        end do
    end subroutine solve

end program tdma_benchmark_final
