!------------------------------------------------------------------------
MODULE postprocess_green
!------------------------------------------------------------------------
!! Driver of postprocessing the Green function to calculate dos and spin-dos
!------------------------------------------------------------------------
  USE comms
  USE parameters
  USE hamiltonian
  !
  IMPLICIT NONE
  SAVE
  !
  COMPLEX(DP), ALLOCATABLE :: green_s(:,:), green_s1(:,:), green_b(:,:), &
    t_mat(:,:), s_mat(:,:), teff_mat(:,:), seff_mat(:,:)
  COMPLEX(DP), ALLOCATABLE :: green_layer(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: spnr(:,:,:,:)
  REAL(DP), ALLOCATABLE :: srvec(:,:), sndegen(:)
  INTEGER :: snrpts
  !
  ! PUBLIC :: get_dos_s, get_dos_b, pp_green_allocate, pp_green_deallocate
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
    dos_out = dos_out + AIMAG(green_b(i,i)) / PI
  END DO
END SUBROUTINE get_dos_b
!------------------------------------------------------------------------

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
!         dos_out = dos_out + AIMAG(green_s(i,i)) / PI
!     END DO
! END SUBROUTINE
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
!         dos_out = dos_out + AIMAG(green_s(i,i)) / PI
!     END DO
! END SUBROUTINE
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
!                 spin_out(ispin) = spin_out(ispin) + AIMAG(spnk(i,j)*green_s(j,i)) / PI
!             END DO
!         END DO
!     END DO
!     DEALLOCATE(spnk)
! END SUBROUTINE
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
!                 spin_out(ispin) = spin_out(ispin) + AIMAG(spnk(i,j)*green_s(j,i)) / PI
!             END DO
!         END DO
!     END DO
!     DEALLOCATE(spnk)
! END SUBROUTINE



SUBROUTINE get_dos_nlayer(nlayer, dos_out)
  INTEGER, INTENT(IN) :: nlayer
  COMPLEX(DP), INTENT(OUT) :: dos_out(64)
  INTEGER :: i, il, j, cnt
  IF (nlayer .lt. 1) RETURN
  ! recurrence to calculate green_ilayer
  CALL set_green_layer(nlayer)
  dos_out = czero
  ! DO il = 1, nlayer
  !     ! calculate trace and DOS
  !     DO i = 1, nbulk
  !         IF ( AIMAG(green_layer(i,i,il)) < 0 ) THEN
  !             WRITE(*,'("negative DOS: green_",I," basis ",I3,", value: ",ES8.1)')  il, i, AIMAG(green_layer(i,i,il))
  !             ! STOP
  !         END IF
  !         dos_out(i, il) = AIMAG(green_layer(i,i,il)) / PI
  !     END DO
  ! END DO
  cnt = 1
  do i = 41, 48
    do j = 41, 48
      dos_out(cnt) = green_layer(i,j,1)
      cnt = cnt + 1
    end do
  end do
END SUBROUTINE

SUBROUTINE get_spin_s(spin_out)
  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: spin_out(3)
  COMPLEX(DP), ALLOCATABLE :: spnk(:,:)
  INTEGER :: i, j, ispin
  IF (.NOT. isspin) CALL io_error('isspin must be .true. to calculate spinDOS')
  ALLOCATE(spnk(nsurf,nsurf))
  spin_out = 0.0_dp
  DO ispin = 1, 3
    CALL k_operator(snrpts, spnr(1:nsurf, 1:nsurf,:,ispin), srvec, sndegen, kx, ky, spnk)
    DO i = 1, nsurf
      DO j = 1, nsurf
        spin_out(ispin) = spin_out(ispin) + AIMAG(spnk(i,j)*green_s(j,i)) / PI
      END DO
    END DO
  END DO
  DEALLOCATE(spnk)
END SUBROUTINE

SUBROUTINE set_green_layer(nl)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nl
  INTEGER :: nl_set, il_from, il_to, il
  COMPLEX(DP), ALLOCATABLE :: green_layer_temp(:,:,:)
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: temp_mat1, temp_mat2, temp_mat3
  ! allocate green_layer appropriate size
  ! and set the range of ilayer(il) to calculate green_layer
  IF(ALLOCATED(green_layer)) THEN
    nl_set = SIZE(green_layer,3)
    IF (nl <= nl_set) RETURN
    ! ELSE nl > nl_set, extend SIZE(green_layer,3) from nl_set to nl
    ALLOCATE(green_layer_temp(nbulk,nbulk,nl_set))
    green_layer_temp = green_layer
    DEALLOCATE(green_layer)
    ALLOCATE(green_layer(nbulk,nbulk,nl))
    green_layer(:,:,1:nl_set) = green_layer_temp
    DEALLOCATE(green_layer_temp)
    il_from = nl_set+1
    il_to = nl
  ELSE ! green_layer not allocated
    ALLOCATE(green_layer(nbulk,nbulk,nl))
    il_from = 1
    il_to = nl
  ENDIF

  ! calculate green_layer using recurrence relations
  DO il = il_from, il_to
    IF ((il==1) .AND. isslab) THEN
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
    ELSE IF((il==1) .AND. .NOT.(isslab)) THEN
      ! temp_mat1 = t_mat * green_s
      ! green_layer(:,:,1) = green_s + temp_mat1 * s_mat
      ALLOCATE(temp_mat1(nbulk,nbulk))
      CALL ZCOPY(nbulk*nbulk, green_s, 1, green_layer(:,:,1), 1)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, t_mat, nbulk, green_s, nbulk,czero, temp_mat1, nbulk)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, temp_mat1, nbulk, s_mat, nbulk,cone, green_layer(:,:,1), nbulk)
      DEALLOCATE(temp_mat1)
    ELSE IF ((il>1) .AND. isslab) THEN
      ! green_layer(:,:,il) = teff_mat + ( (teff_mat * h12^dagger) * green_layer(:,:,il-1)) * s_mat
      ALLOCATE(temp_mat1(nbulk,nbulk))
      ALLOCATE(temp_mat2(nbulk,nbulk))
      CALL ZCOPY(nbulk*nbulk, teff_mat, 1, green_layer(:,:,il), 1)
      CALL ZGEMM('N','C',nbulk,nbulk,nbulk,cone, teff_mat, nbulk, h12, nbulk,czero, temp_mat1, nbulk)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, temp_mat1, nbulk, green_layer(:,:,il-1), nbulk,czero, temp_mat2, nbulk)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, temp_mat2, nbulk, s_mat, nbulk,cone, green_layer(:,:,il), nbulk)
      DEALLOCATE(temp_mat1)
      DEALLOCATE(temp_mat2)
    ELSE ! il > 1 and .NOT. isslab
      ! green_layer(:,:,il) = green_s + ( (t_mat * green_layer(:,:,il-1)) * s_mat )
      ALLOCATE(temp_mat1(nbulk,nbulk))
      CALL ZCOPY(nbulk*nbulk, green_s, 1, green_layer(:,:,il), 1)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, t_mat, nbulk, green_layer(:,:,il-1), nbulk,czero, temp_mat1, nbulk)
      CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone, temp_mat1, nbulk, s_mat, nbulk,cone, green_layer(:,:,il), nbulk)
      DEALLOCATE(temp_mat1)
    END IF
  END DO
END SUBROUTINE

SUBROUTINE pp_green_deallocate_green_layer()
  IMPLICIT NONE
  IF(ALLOCATED(green_layer)) DEALLOCATE(green_layer)
END SUBROUTINE

SUBROUTINE set_transfer_mat()
  COMPLEX(DP), ALLOCATABLE :: temp_mat(:,:)
  IF (isslab) THEN
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
  ELSE ! .NOT. isslab
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
  IF (isslab) ALLOCATE(green_s1(nbulk,nbulk), STAT=ierr)
  IF (isslab) ALLOCATE(teff_mat(nbulk,nbulk))
  IF (isslab) ALLOCATE(seff_mat(nbulk,nbulk))
  IF (ierr /= 0) CALL io_error('Error during allocation in pp_green_allocate')
  IF (isspin) CALL read_spin(0, spnr, srvec, sndegen, snrpts)
END SUBROUTINE pp_green_setup
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE read_spin(nr3, spnr_nr3, rvec_nr3, ndegen_nr3, irpts_nr3)
!! Read seedname_spnr.dat. Save spin elements <m0|S|nR> with R_z = nr3
!! to variables spnr_nr3.
  USE comms, ONLY : find_free_unit
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nr3
  COMPLEX(DP), ALLOCATABLE, INTENT(OUT) :: spnr_nr3(:,:,:,:)
  REAL(DP), ALLOCATABLE, INTENT(OUT)  :: rvec_nr3(:,:), ndegen_nr3(:)
  INTEGER, INTENT(OUT) :: irpts_nr3
  !
  INTEGER :: iun
  INTEGER :: irpts,iw,ni,nj,ierr,nrpts, ir1,ir2,ir3, num_spnr_wann
  REAL(DP) :: sp1r, sp1i, sp2r, sp2i, sp3r, sp3i
  COMPLEX(DP), ALLOCATABLE :: spnr_nr3_tmp(:,:,:,:)
  REAL(DP), ALLOCATABLE :: ndegen(:), rvec_nr3_tmp(:,:), ndegen_nr3_tmp(:)
  CHARACTER(LEN=80) :: buffer
  LOGICAL :: ir3_checked
  CHARACTER(LEN=256) :: filename
  !
  iun = find_free_unit()
  filename = TRIM(seedname) // '_spnr.dat'
  OPEN(UNIT=iun, FILE=TRIM(filename), ACTION='read', IOSTAT=ierr)
  IF (ierr /= 0) CALL io_error('Error opening ' // TRIM(filename))
  !
  READ(iun,'(80a)', IOSTAT=ierr) buffer
  READ(iun,'(80a)', IOSTAT=ierr) buffer
  READ(buffer, *) num_spnr_wann
  READ(iun,'(80a)', IOSTAT=ierr) buffer
  READ(buffer, *) nrpts
  IF(is_root) WRITE(*,*) 'begin reading ' // TRIM(filename)
  IF(is_root) WRITE(*,*) 'num_spnr_wann = ', num_spnr_wann
  IF(is_root) WRITE(*,*) 'nrpts = ', nrpts
  !
  ALLOCATE(ndegen(nrpts))
  ALLOCATE(spnr_nr3_tmp(num_spnr_wann, num_spnr_wann, nrpts, 3))
  ALLOCATE(rvec_nr3_tmp(2, nrpts))
  ALLOCATE(ndegen_nr3_tmp(nrpts))

  DO irpts = 1, (nrpts-1)/15
    READ(iun, *, IOSTAT=ierr) ndegen((irpts-1)*15+1:(irpts-1)*15+15)
  ENDDO
  READ(iun, *, IOSTAT=ierr) ndegen(((nrpts-1)/15)*15+1:nrpts)
  !
  irpts_nr3 = 0
  DO irpts = 1, nrpts
    ir3_checked = .false.
    DO iw = 1,num_spnr_wann**2
      READ(iun, *, IOSTAT=ierr) ir1, ir2, ir3, ni, nj, &
                                sp1r, sp1i, sp2r, sp2i, sp3r, sp3i
      IF (ierr /= 0) CALL io_error('Error reading file '// TRIM(filename))
      ! Choose lattice vectors with R_3 = nr3
      IF (ir3 == nr3) THEN
        IF (.NOT. ir3_checked) THEN
          irpts_nr3 = irpts_nr3 + 1
          ir3_checked = .true.
        ENDIF
        spnr_nr3_tmp(ni, nj, irpts_nr3, 1) = CMPLX(sp1r, sp1i, DP)
        spnr_nr3_tmp(ni, nj, irpts_nr3, 2) = CMPLX(sp2r, sp2i, DP)
        spnr_nr3_tmp(ni, nj, irpts_nr3, 3) = CMPLX(sp3r, sp3i, DP)
        rvec_nr3_tmp(1, irpts_nr3) = ir1
        rvec_nr3_tmp(2, irpts_nr3) = ir2
        ndegen_nr3_tmp(irpts_nr3) = ndegen(irpts)
      ENDIF
    ENDDO
  ENDDO
  CLOSE(iun)
  !
  ! Change size of array to the read size
  ! This is required because the number of R vectors with R3=nr3 is not known
  ! before reading.
  ALLOCATE(spnr_nr3(num_spnr_wann, num_spnr_wann, irpts_nr3, 3))
  ALLOCATE(rvec_nr3(2, irpts_nr3))
  ALLOCATE(ndegen_nr3(irpts_nr3))
  spnr_nr3 = spnr_nr3_tmp(:,:,1:irpts_nr3,:)
  rvec_nr3 = rvec_nr3_tmp(:,1:irpts_nr3)
  ndegen_nr3 = ndegen_nr3_tmp(1:irpts_nr3)
  DEALLOCATE(spnr_nr3_tmp)
  DEALLOCATE(rvec_nr3_tmp)
  DEALLOCATE(ndegen_nr3_tmp)
  DEALLOCATE(ndegen)
  !
  IF(is_root) write(*,'("Spin_",I1," shape : ",4I5)') nr3, shape(spnr_nr3)
  IF(is_root) write(*,*) 'end reading ' // TRIM(filename)
  !
END SUBROUTINE read_spin
!------------------------------------------------------------------------

END MODULE postprocess_green