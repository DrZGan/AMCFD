program tdma_benchmark
    use omp_lib
    implicit none

    ! Grid dimensions
    integer, parameter :: ni = 128, nj = 128, nk = 128
    integer, parameter :: nim1 = ni - 1, njm1 = nj - 1, nkm1 = nk - 1
    integer, parameter :: n_iterations = 10

    ! Arrays
    real(8), allocatable :: enthalpy(:,:,:)
    real(8), allocatable :: ap(:,:,:), ae(:,:,:), aw(:,:,:)
    real(8), allocatable :: an(:,:,:), as(:,:,:), at(:,:,:), ab(:,:,:)
    real(8), allocatable :: su(:,:,:)

    ! Timing
    real(8) :: start_time, end_time, total_time
    integer :: iter

    ! Allocate arrays
    allocate(enthalpy(ni, nj, nk))
    allocate(ap(ni, nj, nk), ae(ni, nj, nk), aw(ni, nj, nk))
    allocate(an(ni, nj, nk), as(ni, nj, nk), at(ni, nj, nk), ab(ni, nj, nk))
    allocate(su(ni, nj, nk))

    ! Initialize arrays with some values
    call initialize_arrays()

    ! Warmup run
    call tdma_solve()

    ! Benchmark
    total_time = 0.0d0
    do iter = 1, n_iterations
        start_time = omp_get_wtime()
        call tdma_solve()
        end_time = omp_get_wtime()
        total_time = total_time + (end_time - start_time)
        write(*,'(A,I3,A,F10.6,A)') 'Iteration ', iter, ': ', (end_time - start_time)*1000.0d0, ' ms'
    end do

    write(*,'(A)') '----------------------------------------'
    write(*,'(A,I6,A,I6,A,I6)') 'Grid size: ', ni, ' x ', nj, ' x ', nk
    write(*,'(A,I3)') 'Number of iterations: ', n_iterations
    write(*,'(A,F10.6,A)') 'Average time per iteration: ', (total_time/n_iterations)*1000.0d0, ' ms'
    write(*,'(A,F10.6,A)') 'Total time: ', total_time*1000.0d0, ' ms'
    write(*,'(A,I3)') 'Number of OpenMP threads: ', omp_get_max_threads()

    ! Print checksum for verification
    write(*,'(A,E15.8)') 'Checksum (sum of enthalpy): ', sum(enthalpy)

    ! Deallocate
    deallocate(enthalpy, ap, ae, aw, an, as, at, ab, su)

contains

    subroutine initialize_arrays()
        integer :: i, j, k

        ! Initialize with some non-trivial values
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

        do ksweep = 1, 2                          ! scan k direction twice
            do k = nkm1, 2, -1
                do jsweep = 1, 2                  ! scan j direction twice
                    !$OMP PARALLEL PRIVATE(pr, qr, d, denom)
                    !$OMP DO
                    do j = 2, njm1
                        pr(1) = 0.0d0
                        qr(1) = enthalpy(1,j,k)

                        do i = 2, nim1
                            d = at(i,j,k)*enthalpy(i,j,k+1) + ab(i,j,k)*enthalpy(i,j,k-1) + &
                                an(i,j,k)*enthalpy(i,j+1,k) + as(i,j,k)*enthalpy(i,j-1,k) + su(i,j,k)
                            denom = ap(i,j,k) - aw(i,j,k)*pr(i-1)

                            if (denom <= 1e-12 .and. denom >= 0) denom = denom + 1e-13
                            if (denom >= -1e-12 .and. denom < 0) denom = denom - 1e-13

                            pr(i) = ae(i,j,k) / denom
                            qr(i) = (d + aw(i,j,k)*qr(i-1)) / denom
                        end do

                        ! Back substitution
                        do i = nim1, 2, -1
                            enthalpy(i,j,k) = pr(i)*enthalpy(i+1,j,k) + qr(i)
                        end do
                    end do
                    !$OMP END DO
                    !$OMP END PARALLEL
                end do
            end do
        end do
    end subroutine tdma_solve

end program tdma_benchmark
