!------------------------------------------------------------------------
MODULE iter_bulk
!------------------------------------------------------------------------
!! Driver of iteration for the bulk-only case (no surface modification)
!------------------------------------------------------------------------
  USE comms
  USE parameters
  USE hamiltonian, ONLY : h11, h12, omega
  USE postprocess_green, ONLY : green_s, green_b
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! private temporary variables used in iter_bulk_update and iter_bulk_green
  COMPLEX(DP), ALLOCATABLE :: a_b(:,:), b_b(:,:), e_s(:,:), e_b(:,:), &
    inv_e_b(:,:), a_b_prev(:,:), b_b_prev(:,:), temp_mat(:,:)
  INTEGER, ALLOCATABLE :: ipiv(:)
  ! private variables used in iter_bulk_check_convergence
  REAL(DP) :: conver1, conver2
  PUBLIC :: iter_bulk_main, iter_bulk_allocate, iter_bulk_deallocate
CONTAINS
!------------------------------------------------------------------------
SUBROUTINE iter_bulk_main()
!------------------------------------------------------------------------
!! main driver of the iteration
!------------------------------------------------------------------------
  CALL iter_bulk_initialize()
  ! ... iterate semi-infinite slab
  DO n_iter = 1, MAX_N_ITER
    CALL iter_bulk_update()
    CALL iter_bulk_check_convergence()
    IF (flag_converged) EXIT
  ENDDO
  CALL iter_bulk_green()
!------------------------------------------------------------------------
END SUBROUTINE iter_bulk_main
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_bulk_initialize()
!! initialize iteration paramters
  IMPLICIT NONE
  e_s = h11
  e_b = h11
  a_b = h12
  b_b = TRANSPOSE(CONJG(h12))
!------------------------------------------------------------------------
END SUBROUTINE iter_bulk_initialize
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_bulk_check_convergence()
!------------------------------------------------------------------------
!! check convergences of hopping matrices
!------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: i, j
  conver1 = 0.0_dp
  conver2 = 0.0_dp
  flag_converged = .false.
  DO j = 1, nbulk
    DO i = 1, nbulk
      conver1 = conver1 + SQRT(REAL(a_b(i,j), DP)**2 + AIMAG(a_b(i,j))**2)
      conver2 = conver2 + SQRT(REAL(b_b(i,j), DP)**2 + AIMAG(b_b(i,j))**2)
    END DO
  END DO
  IF (conver1 < hopping_tol .AND. conver2 < hopping_tol) THEN
    flag_converged = .true.
  END IF
!------------------------------------------------------------------------
END SUBROUTINE iter_bulk_check_convergence
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_bulk_update()
!------------------------------------------------------------------------
!! use recurrence relations to update matrices
!------------------------------------------------------------------------
  IMPLICIT NONE
  ! inv_e_b = (omgea - e_b)^-1
  CALL inv_omega_minus_mat(nbulk, e_b, omega, inv_e_b, 'iter_bulk_update for inv_e_b')
  !
  CALL ZCOPY(nbulk*nbulk, a_b, 1, a_b_prev, 1)
  CALL ZCOPY(nbulk*nbulk, b_b, 1, b_b_prev, 1)
  !
  ! temp_mat = a_b_prev * inv_e_b
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,a_b_prev,nbulk,inv_e_b,nbulk,czero,temp_mat,nbulk)
  ! e_s = e_s(prev) + temp_mat * b_b_prev
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,temp_mat,nbulk,b_b_prev,nbulk,cone,e_s,nbulk)
  ! e_b = e_b(prev) + temp_mat * b_b_prev
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,temp_mat,nbulk,b_b_prev,nbulk,cone,e_b,nbulk)
  ! a_b = temp_mat * a_b_prev
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,temp_mat,nbulk,a_b_prev,nbulk,czero,a_b,nbulk)
  !
  ! temp_mat = b_b_prev * inv_e_b
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,b_b_prev,nbulk,inv_e_b,nbulk,czero,temp_mat,nbulk)
  ! e_b = e_b(prev) + temp_mat * a_b_prev
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,temp_mat,nbulk,a_b_prev,nbulk,cone,e_b,nbulk)
  ! b_b = temp_mat * b_b_prev
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,temp_mat,nbulk,b_b_prev,nbulk,czero,b_b,nbulk)
!------------------------------------------------------------------------
END SUBROUTINE iter_bulk_update
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_bulk_green()
!! Calculate Green function after iteration
  IMPLICIT NONE
  CALL inv_omega_minus_mat(nbulk, e_s, omega, green_s, 'iter_bulk_green for green_s')
  CALL inv_omega_minus_mat(nbulk, e_b, omega, green_b, 'iter_bulk_green for green_b')
!------------------------------------------------------------------------
END SUBROUTINE
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_bulk_allocate()
!! Allocate matrices used during the iteration
  IMPLICIT NONE
  INTEGER :: ierr
  ALLOCATE(a_b(nbulk, nbulk), stat=ierr)
  ALLOCATE(b_b(nbulk, nbulk), stat=ierr)
  ALLOCATE(e_s(nbulk, nbulk), stat=ierr)
  ALLOCATE(e_b(nbulk, nbulk), stat=ierr)
  ALLOCATE(inv_e_b(nbulk,nbulk), stat=ierr)
  ALLOCATE(a_b_prev(nbulk,nbulk), stat=ierr)
  ALLOCATE(b_b_prev(nbulk,nbulk), stat=ierr)
  ALLOCATE(temp_mat(nbulk,nbulk), stat=ierr)
  ALLOCATE(ipiv(nbulk), stat=ierr)
  IF (ierr /= 0) CALL io_error('Error during allocation in iter_bulk_allocate')
!------------------------------------------------------------------------
END SUBROUTINE
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_bulk_deallocate()
!! Deallocate matrices when the iteration is done
  IMPLICIT NONE
  INTEGER :: ierr
  DEALLOCATE(a_b, stat=ierr)
  DEALLOCATE(b_b, stat=ierr)
  DEALLOCATE(e_s, stat=ierr)
  DEALLOCATE(e_b, stat=ierr)
  DEALLOCATE(inv_e_b, stat=ierr)
  DEALLOCATE(a_b_prev, stat=ierr)
  DEALLOCATE(b_b_prev, stat=ierr)
  DEALLOCATE(temp_mat, stat=ierr)
  DEALLOCATE(ipiv, stat=ierr)
  IF (ierr /= 0) CALL io_error('Error during deallocation in iter_bulk_deallocate')
!------------------------------------------------------------------------
END SUBROUTINE
!------------------------------------------------------------------------
END MODULE iter_bulk