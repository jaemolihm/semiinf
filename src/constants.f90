MODULE constants
    USE mpi
    IMPLICIT NONE
    SAVE
    INTEGER, PARAMETER :: DP = selected_real_kind(14, 200)
    REAL(DP), PARAMETER :: PI = 4. * ATAN(1.)
    complex(dp), PARAMETER :: CMPLX_IM = CMPLX(0.0, 1.0, dp)
    complex(dp), PARAMETER :: CMPLX_1 = CMPLX(1.0, 0.0, dp)
    complex(dp), PARAMETER :: CMPLX_0 = CMPLX(0.0, 0.0, dp)

    ! variables for mpi
    INTEGER, PARAMETER :: root_process = 0
    LOGICAL :: is_root
    INTEGER :: my_id, num_procs
CONTAINS
    SUBROUTINE mpi_setup()
        IMPLICIT NONE
        INTEGER :: ierr
        CALL MPI_INIT(ierr)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
        is_root = .false.
        IF (my_id == root_process) is_root = .true.
    END SUBROUTINE
    SUBROUTINE mpi_end()
        IMPLICIT NONE
        INTEGER :: ierr
        CALL MPI_FINALIZE(ierr)
    END SUBROUTINE
    ! SUBROUTINE iden_cmplx(iden, n_size)
    !     complex(dp), allocatable :: iden(:,:)
    !     integer, intent(in) :: n_size
    !     integer :: iw
    !     allocate(iden(n_size, n_size))
    !     iden = CMPLX(0.,0.,DP)
    !     do iw = 1, n_size
    !         iden(iw,iw) = CMPLX(1.,0.,DP)
    !     end do
    ! END SUBROUTINE
    SUBROUTINE io_error(error_msg)
        IMPLICIT NONE
        character(len=*), intent(IN) :: error_msg
        write(*,*) error_msg
        STOP
    END SUBROUTINE

    SUBROUTINE linspace(out_list, start, end, ind_start, ind_end)
        real(DP), intent(IN) :: start, end
        integer, intent(IN) :: ind_start, ind_end
        real(DP), dimension(:), intent(OUT) :: out_list
        real(DP) :: value, step
        integer :: i
        step = (end - start) / (ind_end - ind_start)
        value = start
        DO i = ind_start, ind_end
            out_list(i) = value
            value = value + step
        END DO
    END SUBROUTINE

    SUBROUTINE inv_omega_minus_mat(ndim, mat_in, omega_, mat_out, flag)
        ! mat_out = (omega_ - mat_in)^-1
        ! mat_in and mat_out are square matrices with dimension(ndim, ndim)
        INTEGER, INTENT(IN) :: ndim
        COMPLEX(dp), INTENT(IN) :: omega_
        complex(DP), INTENT(IN) :: mat_in(ndim,ndim)
        complex(DP), INTENT(OUT) :: mat_out(ndim,ndim)
        CHARACTER(len=*), INTENT(IN) :: flag

        complex(DP), allocatable :: temp_mat_inv(:,:)
        INTEGER, allocatable :: ipiv(:)
        INTEGER :: i, info
        ALLOCATE(temp_mat_inv(ndim,ndim))
        ALLOCATE(ipiv(ndim))
        ! temp_mat_inv = omega_ * iden - mat_in
        ! mat_out = identity
        temp_mat_inv = -mat_in
        mat_out = CMPLX_0
        DO i = 1, ndim
            temp_mat_inv(i,i) = omega_ + temp_mat_inv(i,i)
            mat_out(i,i) = CMPLX_1
        END DO
        ! mat_out = (temp_mat_inv)^-1
        CALL ZGESV(ndim,ndim, temp_mat_inv, ndim,ipiv, mat_out, ndim,info)
        IF (info /= 0) CALL io_error('Error in ZGESV: matinv inversion in ' // flag)
        DEALLOCATE(temp_mat_inv)
        DEALLOCATE(ipiv)
    END SUBROUTINE

    subroutine k_operator(nrpts, oper_r, rvec, ndegen, kx_, ky_, oper_k)
        IMPLICIT NONE
        integer, intent(in) :: nrpts
        real(dp), intent(in) :: kx_, ky_
        complex(dp), intent(in) :: oper_r(:,:,:)
        real(dp), intent(in)  :: rvec(:,:), ndegen(:)
        complex(dp), intent(out) :: oper_k(:,:)

        integer :: irpt
        real(dp) :: kx_crystal, ky_crystal
        complex(dp) :: coeff_ik

        kx_crystal = 2.0_dp*PI*kx_
        ky_crystal = 2.0_dp*PI*ky_
        
        oper_k = CMPLX(0.0, 0.0, dp)
        do irpt = 1, nrpts
            coeff_ik = EXP( CMPLX_IM * ( kx_crystal*rvec(1,irpt) + ky_crystal*rvec(2,irpt) ) ) / ndegen(irpt)
            oper_k = oper_k + oper_r(:,:,irpt) * coeff_ik
        end do
    end subroutine
END MODULE
