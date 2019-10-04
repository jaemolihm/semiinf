!------------------------------------------------------------------------
module iter_slab
!------------------------------------------------------------------------
!! Driver of iteration for the bulk-only case (no surface modification)
!------------------------------------------------------------------------
  USE comms, ONLY : DP, io_error, cone, czero
  USE parameters, ONLY : nbulk, nsurf
  USE hamiltonian, ONLY : h00, h01, h11, h12, omega
  USE postprocess_green, ONLY : green_s, green_s1, green_b
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! private temporary variables used in iter_slab_update and iter_slab_green
  COMPLEX(DP), ALLOCATABLE :: a_s(:,:), a_b(:,:), b_s(:,:), b_b(:,:), &
    e_s0(:,:), e_s1(:,:), e_b(:,:), inv_e_b(:,:), a_b_prev(:,:), &
    b_b_prev(:,:), a_s_prev(:,:), b_s_prev(:,:), temp_mat(:,:), &
    temp_mat_sb(:,:)
  INTEGER, ALLOCATABLE :: ipiv_s(:), ipiv_b(:)
  PUBLIC :: iter_slab_main, iter_slab_allocate, iter_slab_deallocate
CONTAINS
!------------------------------------------------------------------------
SUBROUTINE iter_slab_main(flag_converged)
!------------------------------------------------------------------------
!! main driver of the iteration
!------------------------------------------------------------------------
  USE parameters, ONLY : max_n_iter
  IMPLICIT NONE
  !
  LOGICAL, INTENT(INOUT) :: flag_converged
  INTEGER :: n_iter
  flag_converged = .FALSE.
  !
  CALL iter_slab_initialize()
  ! ... iterate semi-infinite slab
  DO n_iter = 1, max_n_iter
    CALL iter_slab_update()
    CALL iter_slab_check_convergence(flag_converged)
    IF (flag_converged) EXIT
  ENDDO
  CALL iter_slab_green()
!------------------------------------------------------------------------
END SUBROUTINE iter_slab_main
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_slab_initialize()
!! initialize iteration paramters
  IMPLICIT NONE
  e_s0 = h00
  e_s1 = h11
  e_b = h11
  a_s = h01
  b_s = TRANSPOSE(CONJG(h01))
  a_b = h12
  b_b = TRANSPOSE(CONJG(h12))
!------------------------------------------------------------------------
END SUBROUTINE iter_slab_initialize
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_slab_check_convergence(flag_converged)
!------------------------------------------------------------------------
!! check convergences of hopping matrices
!------------------------------------------------------------------------
  USE parameters, ONLY : hopping_tol
  IMPLICIT NONE
  !
  LOGICAL, INTENT(OUT) :: flag_converged
  !! Set to True if convergence is reached. False otherwise.
  !
  INTEGER :: i, j
  REAL(DP) :: conver_as, conver_ab, conver_bs, conver_bb
  !
  conver_ab = 0.0_dp
  conver_bb = 0.0_dp
  conver_as = 0.0_dp
  conver_bs = 0.0_dp
  flag_converged = .FALSE.
  DO j = 1, nbulk
    DO i = 1, nbulk
      conver_ab = conver_ab + SQRT(REAL(a_b(i,j), DP)**2+AIMAG(a_b(i,j))**2)
      conver_bb = conver_bb + SQRT(REAL(b_b(i,j), DP)**2+AIMAG(b_b(i,j))**2)
    END DO
  END DO
  DO j = 1, nbulk
    DO i = 1, nsurf
      conver_as = conver_as + SQRT(REAL(a_s(i,j), DP)**2+AIMAG(a_s(i,j))**2)
    END DO
  END DO
  DO j = 1, nsurf
    DO i = 1, nbulk
      conver_bs = conver_bs + SQRT(REAL(b_s(i,j), DP)**2+AIMAG(b_s(i,j))**2)
    END DO
  END DO
  IF (conver_ab < hopping_tol .AND. conver_bb < hopping_tol &
    .AND. conver_as < hopping_tol .AND. conver_bs < hopping_tol) THEN
    flag_converged = .TRUE.
  END IF
!------------------------------------------------------------------------
END SUBROUTINE iter_slab_check_convergence
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_slab_update()
!------------------------------------------------------------------------
!! use recurrence relation to update matrices
!------------------------------------------------------------------------
  USE comms, ONLY : inv_omega_minus_mat
  IMPLICIT NONE
  ! inv_e_b = (omgea - e_b)^-1
  CALL inv_omega_minus_mat(nbulk, e_b, omega, inv_e_b, 'iter_bulk_update for inv_e_b')
  !
  CALL ZCOPY(nsurf*nbulk, a_s, 1, a_s_prev, 1)
  CALL ZCOPY(nbulk*nsurf, b_s, 1, b_s_prev, 1)
  CALL ZCOPY(nbulk*nbulk, a_b, 1, a_b_prev, 1)
  CALL ZCOPY(nbulk*nbulk, b_b, 1, b_b_prev, 1)
  !
  ! temp_mat = inv_e_b * a_b
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,inv_e_b,nbulk,a_b,nbulk,czero,temp_mat,nbulk)
  ! a_s = a_s_prev * temp_mat
  CALL ZGEMM('N','N',nsurf,nbulk,nbulk,cone,a_s_prev,nsurf,temp_mat,nbulk,czero,a_s,nsurf)
  ! a_b = a_b_prev * temp_mat
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,a_b_prev,nbulk,temp_mat,nbulk,czero,a_b,nbulk)
  ! e_b = e_b(prev) + b_b_prev * temp_mat
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,b_b_prev,nbulk,temp_mat,nbulk,cone,e_b,nbulk)
  !
  ! temp_mat = b_b_prev * inv_e_b
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,b_b_prev,nbulk,inv_e_b,nbulk,czero,temp_mat,nbulk)
  ! b_s = temp_mat * b_s_prev
  CALL ZGEMM('N','N',nbulk,nsurf,nbulk,cone,temp_mat,nbulk,b_s_prev,nbulk,czero,b_s,nbulk)
  ! b_b = temp_mat * b_b_prev
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,temp_mat,nbulk,b_b_prev,nbulk,czero,b_b,nbulk)
  !
  ! temp_mat = inv_e_b * b_b_prev
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,inv_e_b,nbulk,b_b_prev,nbulk,czero,temp_mat,nbulk)
  ! e_b = e_b(prev) + a_b_prev * temp_mat
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,a_b_prev,nbulk,temp_mat,nbulk,cone,e_b,nbulk)
  ! e_s1 = e_s1(prev) + a_b_prev * temp_mat
  CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,a_b_prev,nbulk,temp_mat,nbulk,cone,e_s1,nbulk)
  !
  ! temp_mat_sb = a_s_prev * inv_e_b
  CALL ZGEMM('N','N',nsurf,nbulk,nbulk,cone,a_s_prev,nsurf,inv_e_b,nbulk,czero,temp_mat_sb,nsurf)
  ! e_s0 = e_s0(prev) + temp_mat_sb * b_s_prev
  CALL ZGEMM('N','N',nsurf,nsurf,nbulk,cone,temp_mat_sb,nsurf,b_s_prev,nbulk,cone,e_s0,nsurf)
!------------------------------------------------------------------------
END SUBROUTINE iter_slab_update
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_slab_green()
!! Calculate Green function after iteration
  USE comms, ONLY : inv_omega_minus_mat
  IMPLICIT NONE
  CALL inv_omega_minus_mat(nsurf, e_s0, omega, green_s, 'iter_slab_green for green_s')
  CALL inv_omega_minus_mat(nbulk, e_s1, omega, green_s1, 'iter_slab_green for green_s1')
  CALL inv_omega_minus_mat(nbulk, e_b, omega, green_b, 'iter_slab_green for green_b')
END SUBROUTINE iter_slab_green
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_slab_allocate()
!! Allocate matrices used during the iteration
  IMPLICIT NONE
  INTEGER :: ierr
  ALLOCATE(a_s(nsurf, nbulk), STAT=ierr)
  ALLOCATE(b_s(nbulk, nsurf), STAT=ierr)
  ALLOCATE(a_b(nbulk, nbulk), STAT=ierr)
  ALLOCATE(b_b(nbulk, nbulk), STAT=ierr)
  ALLOCATE(e_s0(nsurf, nsurf), STAT=ierr)
  ALLOCATE(e_s1(nbulk, nbulk), STAT=ierr)
  ALLOCATE(e_b(nbulk, nbulk), STAT=ierr)
  ALLOCATE(inv_e_b(nbulk,nbulk), STAT=ierr)
  ALLOCATE(a_s_prev(nsurf,nbulk), STAT=ierr)
  ALLOCATE(b_s_prev(nbulk,nsurf), STAT=ierr)
  ALLOCATE(a_b_prev(nbulk,nbulk), STAT=ierr)
  ALLOCATE(b_b_prev(nbulk,nbulk), STAT=ierr)
  ALLOCATE(temp_mat(nbulk,nbulk), STAT=ierr)
  ALLOCATE(temp_mat_sb(nsurf,nbulk), STAT=ierr)
  ALLOCATE(ipiv_s(nsurf), STAT=ierr)
  ALLOCATE(ipiv_b(nbulk), STAT=ierr)
  if (ierr /= 0) CALL io_error('Error during allocation in iter_slab_allocate')
END SUBROUTINE iter_slab_allocate
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE iter_slab_deallocate()
  !! Deallocate matrices when the iteration is done
  IMPLICIT NONE
  INTEGER :: ierr
  DEALLOCATE(a_s, STAT=ierr)
  DEALLOCATE(b_s, STAT=ierr)
  DEALLOCATE(a_b, STAT=ierr)
  DEALLOCATE(b_b, STAT=ierr)
  DEALLOCATE(e_s0, STAT=ierr)
  DEALLOCATE(e_s1, STAT=ierr)
  DEALLOCATE(e_b, STAT=ierr)
  DEALLOCATE(inv_e_b, STAT=ierr)
  DEALLOCATE(a_s_prev, STAT=ierr)
  DEALLOCATE(b_s_prev, STAT=ierr)
  DEALLOCATE(a_b_prev, STAT=ierr)
  DEALLOCATE(b_b_prev, STAT=ierr)
  DEALLOCATE(temp_mat, STAT=ierr)
  DEALLOCATE(temp_mat_sb, STAT=ierr)
  DEALLOCATE(ipiv_s, STAT=ierr)
  DEALLOCATE(ipiv_b, STAT=ierr)
  if (ierr /= 0) CALL io_error('Error during deallocation in iter_slab_deallocate')
END SUBROUTINE iter_slab_deallocate
!------------------------------------------------------------------------
end module iter_slab
