!------------------------------------------------------------------------
MODULE postprocess_green
!------------------------------------------------------------------------
!! Driver of postprocessing the Green function to calculate dos and spin-dos
!------------------------------------------------------------------------
  USE comms, ONLY : DP, io_error, is_root, cone, czero, pi
  USE parameters, ONLY : seedname, is_ideal_surf, isspin, nbulk, nsurf, green_s, &
    green_s1, green_b, kx, ky, omega
  USE hamiltonian, ONLY : h01, h11, h12
  !
  IMPLICIT NONE
  PRIVATE
  SAVE
  !
  COMPLEX(DP), ALLOCATABLE :: t_mat(:,:), s_mat(:,:), teff_mat(:,:), seff_mat(:,:)
  COMPLEX(DP), ALLOCATABLE :: green_layer(:,:,:)
  !
  PUBLIC :: get_dos_s, get_dos_b, get_dos_nlayer, set_transfer_mat, &
    get_spin_s, pp_green_deallocate_green_layer, pp_green_setup
CONTAINS
!------------------------------------------------------------------------
SUBROUTINE get_dos_s(dos_out)
!! Calculate surface DOS from Green function
  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: dos_out
  INTEGER :: i
  dos_out = 0.0_dp
  DO i = 1, nsurf
    IF ( AIMAG(green_s(i,i)) < 0.d0) THEN
      WRITE(*,'(a,I)') "WARNING: negative DOS in green_b, i = ", i
      ! CALL io_error('Negative DOS in get_dos_s')
    END IF
    dos_out = dos_out + AIMAG(green_s(i,i)) / pi
  END DO
END SUBROUTINE get_dos_s
!------------------------------------------------------------------------
SUBROUTINE get_dos_b(dos_out)
!! Calculate surface DOS from Green function
  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: dos_out
  INTEGER :: i
  dos_out = 0.0_dp
  DO i = 1, nbulk
    IF ( AIMAG(green_b(i,i)) < 0.d0) THEN
      WRITE(*,'(a,I)') "WARNING: negative DOS in green_b, i = ", i
      ! CALL io_error('Negative DOS in get_dos_s')
    END IF
    dos_out = dos_out + AIMAG(green_b(i,i)) / pi
  END DO
END SUBROUTINE get_dos_b
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
! FIXME: to be uncommented
! SUBROUTINE get_dos_up(dos_out)
!     IMPLICIT NONE
!     REAL(DP), INTENT(OUT) :: dos_out
!     INTEGER :: i, j
!     dos_out = 0.0_dp
!     DO j = 1, fin_nup
!         i = fin_ind_up(j)
!         IF ( AIMAG(green_s(i,i)) < 0 ) THEN
!             WRITE(*,*) "negative DOS: green_s basis ", i, ", value: ", AIMAG(green_s(i,i))
!             ! STOP
!         END IF
!         dos_out = dos_out + AIMAG(green_s(i,i)) / pi
!     END DO
! END SUBROUTINE get_dos_up
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
! FIXME: to be uncommented
! SUBROUTINE get_dos_dn(dos_out)
!     IMPLICIT NONE
!     REAL(DP), INTENT(OUT) :: dos_out
!     INTEGER :: i, j
!     dos_out = 0.0_dp
!     DO j = 1, fin_ndn
!         i = fin_ind_dn(j)
!         IF ( AIMAG(green_s(i,i)) < 0 ) THEN
!             WRITE(*,*) "negative DOS: green_s basis ", i, ", value: ", AIMAG(green_s(i,i))
!             ! STOP
!         END IF
!         dos_out = dos_out + AIMAG(green_s(i,i)) / pi
!     END DO
! END SUBROUTINE get_dos_dn
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
! FIXME: to be uncommented
! SUBROUTINE get_spin_up(spin_out)
!     IMPLICIT NONE
!     REAL(DP), INTENT(OUT) :: spin_out(3)
!     COMPLEX(DP), ALLOCATABLE :: spnk(:,:)
!     INTEGER :: i, j, ii, jj, ispin
!     IF (.NOT. isspin) CALL io_error('isspin must be .true. to calculate spinDOS')
!     ALLOCATE(spnk(nsurf,nsurf))
!     spin_out = 0.0_dp
!     DO ispin = 1, 3
!         CALL k_operator(snrpts, spnr(1:nsurf, 1:nsurf,:,ispin), srvec, sndegen, kx, ky, spnk)
!         DO ii = 1, fin_nup
!             i = fin_ind_up(ii)
!             DO jj = 1, fin_nup
!                 j = fin_ind_up(jj)
!                 spin_out(ispin) = spin_out(ispin) + AIMAG(spnk(i,j)*green_s(j,i)) / pi
!             END DO
!         END DO
!     END DO
!     DEALLOCATE(spnk)
! END SUBROUTINE get_spin_up
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
! FIXME: to be uncommented
! SUBROUTINE get_spin_dn(spin_out)
!     IMPLICIT NONE
!     REAL(DP), INTENT(OUT) :: spin_out(3)
!     COMPLEX(DP), ALLOCATABLE :: spnk(:,:)
!     INTEGER :: i, j, ii, jj, ispin
!     IF (.NOT. isspin) CALL io_error('isspin must be .true. to calculate spinDOS')
!     ALLOCATE(spnk(nsurf,nsurf))
!     spin_out = 0.0_dp
!     DO ispin = 1, 3
!         CALL k_operator(snrpts, spnr(1:nsurf, 1:nsurf,:,ispin), srvec, sndegen, kx, ky, spnk)
!         DO ii = 1, fin_ndn
!             i = fin_ind_dn(ii)
!             DO jj = 1, fin_ndn
!                 j = fin_ind_dn(jj)
!                 spin_out(ispin) = spin_out(ispin) + AIMAG(spnk(i,j)*green_s(j,i)) / pi
!             END DO
!         END DO
!     END DO
!     DEALLOCATE(spnk)
! END SUBROUTINE get_spin_dn
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE get_dos_nlayer(nlayer, dos_out)
!! Calculate Green function and DOS for sub-surface layers
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nlayer
  REAL(DP), INTENT(OUT) :: dos_out(nlayer)
  !
  INTEGER :: i, il
  !
  IF (nlayer .lt. 1) RETURN
  !
  ! Use recurrence relation to calculate Green function for sub-surface layers
  CALL set_green_layer(nlayer)
  !
  ! Calculate DOS
  dos_out = 0.d0
  DO il = 1, nlayer
    DO i = 1, nbulk
      IF (AIMAG(green_layer(i, i, il)) < 0.d0) THEN
        WRITE(*,'("WARNING: negative DOS: layer = ",I5," basis ",I5,", value: ",ES9.2)') &
          il, i, AIMAG(green_layer(i,i,il))
        ! STOP
      ENDIF
      dos_out(il) = dos_out(il) + AIMAG(green_layer(i, i, il)) / pi
    ENDDO
  ENDDO
END SUBROUTINE get_dos_nlayer
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE get_spin_s(spin_out)
  USE hamiltonian, ONLY : spn00
  IMPLICIT NONE
  !
  ! TODO: implement inter-layer spin matrix elements
  REAL(DP), INTENT(OUT) :: spin_out(3)
  INTEGER :: i, j, ispin
  IF (.NOT. isspin) CALL io_error('isspin must be .TRUE. to calculate spin-DOS')
  spin_out = 0.0_dp
  DO ispin = 1, 3
    DO i = 1, nsurf
      DO j = 1, nsurf
        spin_out(ispin) = spin_out(ispin) + AIMAG(spn00(i,j,ispin) * green_s(j,i)) / pi
      END DO
    END DO
  END DO
END SUBROUTINE
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE set_green_layer(nl)
!! Calculate Green function for sub-surface layers using recurrence relations
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nl
  !
  INTEGER :: nl_set, il_from, il_to, il
  COMPLEX(DP), ALLOCATABLE :: green_layer_temp(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: temp_mat1(:,:), temp_mat2(:,:), temp_mat3(:,:)
  !
  ! Allocate green_layer appropriate size and set the range of ilayer(il)
  ! to calculate green_layer
  IF (ALLOCATED(green_layer)) THEN
    nl_set = SIZE(green_layer, 3)
    IF (nl <= nl_set) THEN
      ! green_layer already calculated.
      RETURN
    ELSE
      ! green_layer is already calculated, but only for a few layers
      ! extend the size of green_layer from nl_set to nl
      ALLOCATE(green_layer_temp(nbulk, nbulk, nl_set))
      green_layer_temp = green_layer
      DEALLOCATE(green_layer)
      ALLOCATE(green_layer(nbulk, nbulk, nl))
      green_layer(:, :, 1:nl_set) = green_layer_temp
      DEALLOCATE(green_layer_temp)
      il_from = nl_set + 1
      il_to = nl
    ENDIF
  ELSE ! green_layer not allocated
    ALLOCATE(green_layer(nbulk, nbulk, nl))
    il_from = 1
    il_to = nl
  ENDIF
  !
  ! Calculate green_layer using recurrence relations
  DO il = il_from, il_to
    IF (il == 1 .AND. (.NOT. is_ideal_surf)) THEN
      ! green_layer(:,:,1) = teff_mat + ( ( (teff_mat * h01^dagger) * green_s) * h01) * seff_mat
      ALLOCATE(temp_mat1(nbulk,nsurf))
      ALLOCATE(temp_mat2(nbulk,nsurf))
      ALLOCATE(temp_mat3(nbulk,nbulk))
      CALL ZCOPY(nbulk*nbulk, teff_mat, 1, green_layer(:,:,1), 1)
      CALL ZGEMM('N','C',nbulk,nsurf,nbulk,cone, teff_mat, nbulk, h01, nsurf,czero, temp_mat1, nbulk)
      CALL ZGEMM('N','N',nbulk,nsurf,nsurf,cone, temp_mat1, nbulk, green_s, nsurf,czero, temp_mat2, nbulk)
      CALL ZGEMM('N','N',nbulk,nbulk,nsurf,cone, temp_mat2, nbulk, h01, nsurf,czero, temp_mat3, nbulk)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, temp_mat3, nbulk, seff_mat, nbulk,cone, green_layer(:,:,1), nbulk)
      DEALLOCATE(temp_mat1)
      DEALLOCATE(temp_mat2)
      DEALLOCATE(temp_mat3)
    ELSE IF(il == 1 .AND. is_ideal_surf) THEN
      ! temp_mat1 = t_mat * green_s
      ! green_layer(:,:,1) = green_s + temp_mat1 * s_mat
      ALLOCATE(temp_mat1(nbulk,nbulk))
      CALL ZCOPY(nbulk*nbulk, green_s, 1, green_layer(:,:,1), 1)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, t_mat, nbulk, green_s, nbulk,czero, temp_mat1, nbulk)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, temp_mat1, nbulk, s_mat, nbulk,cone, green_layer(:,:,1), nbulk)
      DEALLOCATE(temp_mat1)
    ELSE IF (il > 1 .AND. (.NOT. is_ideal_surf)) THEN
      ! green_layer(:,:,il) = teff_mat + ( (teff_mat * h12^dagger) * green_layer(:,:,il-1)) * s_mat
      ALLOCATE(temp_mat1(nbulk,nbulk))
      ALLOCATE(temp_mat2(nbulk,nbulk))
      CALL ZCOPY(nbulk*nbulk, teff_mat, 1, green_layer(:,:,il), 1)
      CALL ZGEMM('N','C',nbulk,nbulk,nbulk,cone, teff_mat, nbulk, h12, nbulk,czero, temp_mat1, nbulk)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, temp_mat1, nbulk, green_layer(:,:,il-1), nbulk,czero, temp_mat2, nbulk)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, temp_mat2, nbulk, s_mat, nbulk,cone, green_layer(:,:,il), nbulk)
      DEALLOCATE(temp_mat1)
      DEALLOCATE(temp_mat2)
    ELSE ! il > 1 and is_ideal_surf
      ! green_layer(:,:,il) = green_s + ( (t_mat * green_layer(:,:,il-1)) * s_mat )
      ALLOCATE(temp_mat1(nbulk,nbulk))
      CALL ZCOPY(nbulk*nbulk, green_s, 1, green_layer(:,:,il), 1)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, t_mat, nbulk, green_layer(:,:,il-1), nbulk,czero, temp_mat1, nbulk)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, temp_mat1, nbulk, s_mat, nbulk,cone, green_layer(:,:,il), nbulk)
      DEALLOCATE(temp_mat1)
    ENDIF
  ENDDO
END SUBROUTINE set_green_layer
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE pp_green_deallocate_green_layer()
  IMPLICIT NONE
  IF(ALLOCATED(green_layer)) DEALLOCATE(green_layer)
END SUBROUTINE pp_green_deallocate_green_layer
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE set_transfer_mat()
  USE comms, ONLY : inv_omega_minus_mat
  IMPLICIT NONE
  !
  COMPLEX(DP), ALLOCATABLE :: temp_mat(:,:)
  IF (.NOT. is_ideal_surf) THEN
    ! t_mat = green_s1 * h12^dagger
    CALL ZGEMM('N','C',nbulk,nbulk,nbulk,cone,green_s1,nbulk,h12,nbulk,czero,t_mat,nbulk)
    ! s_mat = h12 * green_s1
    CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,h12,nbulk,green_s1,nbulk,czero,s_mat,nbulk)

    ALLOCATE(temp_mat(nbulk,nbulk))
    ! temp_mat = h11 + h12*t_mat
    ! teff_mat = (omega - temp_mat)^-1
    CALL ZCOPY(nbulk*nbulk, h11, 1, temp_mat, 1)
    CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, h12, nbulk, t_mat, nbulk,cone, temp_mat, nbulk)
    CALL inv_omega_minus_mat(nbulk, temp_mat, omega, teff_mat, 'set_transfer_mat for teff_mat')
    ! temp_mat = h11 + s_mat*h12^dagger
    ! seff_mat = (omega - temp_mat)^-1
    CALL ZCOPY(nbulk*nbulk, h11, 1, temp_mat, 1)
    CALL ZGEMM('N','C',nbulk,nbulk,nbulk,cone, s_mat, nbulk, h12, nbulk,cone, temp_mat, nbulk)
    CALL inv_omega_minus_mat(nbulk, temp_mat, omega, seff_mat, 'set_transfer_mat for teff_mat')
    DEALLOCATE(temp_mat)
  ELSE ! is_ideal_surf
    ! t_mat = green_s * h12^dagger
    CALL ZGEMM('N','C',nbulk,nbulk,nbulk,cone,green_s,nbulk,h12,nbulk,czero,t_mat,nbulk)
    ! s_mat = h12 * green_s
    CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,h12,nbulk,green_s,nbulk,czero,s_mat,nbulk)
  ENDIF
END SUBROUTINE set_transfer_mat
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE pp_green_deallocate()
  IMPLICIT NONE
  INTEGER :: ierr
  DEALLOCATE(green_s, STAT=ierr)
  DEALLOCATE(green_b, STAT=ierr)
  IF(ALLOCATED(t_mat)) DEALLOCATE(t_mat,STAT=ierr)
  IF(ALLOCATED(s_mat)) DEALLOCATE(s_mat,STAT=ierr)
  IF(ALLOCATED(teff_mat)) DEALLOCATE(teff_mat,STAT=ierr)
  IF(ALLOCATED(seff_mat)) DEALLOCATE(seff_mat,STAT=ierr)
  if (ierr /= 0) CALL io_error('Error in deallocation in pp_green_allocate')
END SUBROUTINE pp_green_deallocate
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE pp_green_setup()
!! Allocate variables used in each iteration. Read seedname_spnr.dat for
!! the required nr3=0, and save them in spnr, srvec, sndegen, snrpts.
  IMPLICIT NONE
  INTEGER :: ierr
  ALLOCATE(green_s(nsurf,nsurf), STAT=ierr)
  ALLOCATE(green_b(nbulk,nbulk), STAT=ierr)
  ALLOCATE(t_mat(nbulk,nbulk), STAT=ierr)
  ALLOCATE(s_mat(nbulk,nbulk), STAT=ierr)
  IF (.NOT. is_ideal_surf) ALLOCATE(green_s1(nbulk,nbulk), STAT=ierr)
  IF (.NOT. is_ideal_surf) ALLOCATE(teff_mat(nbulk,nbulk))
  IF (.NOT. is_ideal_surf) ALLOCATE(seff_mat(nbulk,nbulk))
  IF (ierr /= 0) CALL io_error('Error during allocation in pp_green_allocate')
END SUBROUTINE pp_green_setup
!------------------------------------------------------------------------
END MODULE postprocess_green