!------------------------------------------------------------------------
MODULE parallel_include
!------------------------------------------------------------------------
!! Taken from QE/UtilXlib/parallel_include.f90
!! Include files for mpi
!------------------------------------------------------------------------
#if defined (__MPI)
#if defined (__MPI_MODULE)
  USE mpi
#else
  INCLUDE 'mpif.h'
#endif
#else
  ! dummy world and null communicator
  INTEGER, PARAMETER :: MPI_COMM_WORLD =  0
  INTEGER, PARAMETER :: MPI_COMM_NULL  = -1
  INTEGER, PARAMETER :: MPI_COMM_SELF  = -2
#endif
END MODULE parallel_include

!------------------------------------------------------------------------
MODULE comms
!------------------------------------------------------------------------
!! Common variables and subroutines used in semiinf.x
!------------------------------------------------------------------------
  USE parallel_include
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER, PARAMETER :: DP = selected_real_kind(14, 200)
  REAL(DP), PARAMETER :: pi = 3.141592653589793238462643383279502884197169399375105820974944_dp
  COMPLEX(DP), PARAMETER :: ci = (0.d0, 1.d0)
  COMPLEX(DP), PARAMETER :: cone = (1.d0, 0.d0)
  COMPLEX(DP), PARAMETER :: czero = (0.d0, 0.d0)
  !
  ! MPI variables
  INTEGER, PARAMETER :: root_process = 0
  LOGICAL :: is_root = .true.
  !! true if this processor is rooot (my_id == 0)
  INTEGER :: num_procs = 1
  !! Number of total processors
  INTEGER :: my_id = 0
  !! Processor index (from 0 to num_procs - 1)
  INTEGER :: world_comm = 0
  !! communicator
CONTAINS
!
!--------------------------------------------------------------------------
SUBROUTINE mp_setup()
  IMPLICIT NONE
  INTEGER :: ierr
#if defined (__MPI)
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  is_root = .false.
  IF (my_id == root_process) is_root = .true.
  world_comm = MPI_COMM_WORLD
#endif
END SUBROUTINE mp_setup
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
SUBROUTINE mp_end()
  IMPLICIT NONE
  INTEGER :: ierr
#if defined (__MPI)
  CALL MPI_FINALIZE(ierr)
#endif
END SUBROUTINE mp_end
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
SUBROUTINE mp_barrier(gid)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: gid
  INTEGER :: ierr
#if defined (__MPI)
  CALL MPI_BARRIER(gid, ierr)
  IF (ierr /= 0) CALL MPI_ABORT(MPI_COMM_WORLD, 8066, ierr)
#endif
END SUBROUTINE mp_barrier
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
SUBROUTINE io_error(error_msg)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: error_msg
  WRITE(*,*) error_msg
  STOP
END SUBROUTINE io_error
!--------------------------------------------------------------------------
!
! Remains only due to legacy. Not used in code.
! !------------------------------------------------------------------------
! SUBROUTINE linspace(out_list, start, end, ind_start, ind_end)
! !------------------------------------------------------------------------
! !! mimic linspace function of python numpy
! !------------------------------------------------------------------------
!   REAL(DP), INTENT(IN) :: start, end
!   INTEGER, INTENT(IN) :: ind_start, ind_end
!   REAL(DP), INTENT(OUT) :: out_list(:)
!   !
!   REAL(DP) :: value, step
!   INTEGER :: i
!   step = (end - start) / (ind_end - ind_start)
!   value = start
!   DO i = ind_start, ind_end
!       out_list(i) = value
!       value = value + step
!   END DO
! END SUBROUTINE linspace
!
!------------------------------------------------------------------------
SUBROUTINE inv_omega_minus_mat(ndim, mat_in, omega, mat_out, flag)
!------------------------------------------------------------------------
!! Compute mat_out = (omega - mat_in)^-1
!------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: ndim
  COMPLEX(DP), INTENT(IN) :: omega
  COMPLEX(DP), INTENT(IN) :: mat_in(ndim, ndim)
  COMPLEX(DP), INTENT(OUT) :: mat_out(ndim, ndim)
  CHARACTER(LEN=*), INTENT(IN) :: flag
  !! Description to be included in the error message when inversion fails
  !
  COMPLEX(DP), ALLOCATABLE :: mat_temp(:, :)
  INTEGER, ALLOCATABLE :: ipiv(:)
  INTEGER :: i
  INTEGER :: info
  !
  ! Set mat_temp = omega - mat_in
  ALLOCATE(mat_temp(ndim, ndim))
  mat_temp = - mat_in
  mat_out = czero
  DO i = 1, ndim
      mat_temp(i, i) = omega + mat_temp(i, i)
      mat_out(i, i) = cone
  END DO
  !
  ! Invert matrix
  ALLOCATE(ipiv(ndim))
  CALL ZGESV(ndim,ndim, mat_temp, ndim, ipiv, mat_out, ndim, info)
  IF (info /= 0) CALL io_error('Error in ZGESV: matinv inversion in ' // flag)
  !
  DEALLOCATE(mat_temp)
  DEALLOCATE(ipiv)
END SUBROUTINE inv_omega_minus_mat
!
!------------------------------------------------------------------------
SUBROUTINE k_operator(nrpts, oper_r, rvec, ndegen, kx, ky, oper_k)
!------------------------------------------------------------------------
!! Compute fourier transform of an operator
!! oper_k(:,:) = exp(1j*k*rvec) * oper_r(:,:,ir) / ndegen(ir)
!------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nrpts
  REAL(DP), INTENT(IN) :: kx, ky
  COMPLEX(DP), INTENT(IN) :: oper_r(:, :, :)
  REAL(DP), INTENT(IN) :: rvec(:, :)
  REAL(DP), INTENT(IN) :: ndegen(:)
  COMPLEX(DP), INTENT(OUT) :: oper_k(:, :)
  !
  INTEGER :: ir
  COMPLEX(DP) :: coeff_ik
  !
  oper_k = czero
  DO ir = 1, nrpts
      coeff_ik = EXP(ci * 2.d0 * pi * (kx * rvec(1, ir) + ky * rvec(2, ir)))
      coeff_ik = coeff_ik / ndegen(ir)
      oper_k = oper_k + oper_r(:,:,ir) * coeff_ik
  ENDDO
!--------------------------------------------------------------------------
END SUBROUTINE k_operator
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
FUNCTION find_free_unit()
!! Taken from QE/UtilXlib/find_free_unit.f90
  IMPLICIT NONE
  !
  INTEGER :: find_free_unit
  INTEGER :: iunit
  LOGICAL :: opnd
  !
  find_free_unit = -1
  DO iunit = 99, 1, -1
     INQUIRE( UNIT = iunit, OPENED = opnd )
     IF ( .NOT. opnd ) THEN
        find_free_unit = iunit
        RETURN
     END IF
  ENDDO
  CALL io_error('find_free_unit(): free unit not found ?!?')
END FUNCTION find_free_unit
END MODULE comms