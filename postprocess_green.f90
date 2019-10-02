MODULE postprocess_green
    USE constants
    USE parameters
    USE hamiltonian
    IMPLICIT NONE
    COMPLEX(dp), allocatable, dimension(:,:), SAVE, PUBLIC :: green_s, green_s1, green_b
    COMPLEX(dp), allocatable, dimension(:,:), SAVE, PUBLIC :: t_mat, s_mat, teff_mat, seff_mat
    COMPLEX(dp), allocatable, SAVE, PUBLIC :: green_layer(:,:,:)

    COMPLEX(dp), allocatable, SAVE :: spnr(:,:,:,:)
    REAL(dp), allocatable, SAVE :: srvec(:,:), sndegen(:)
    INTEGER :: snrpts

    ! PUBLIC :: get_dos_s, get_dos_b, pp_green_allocate, pp_green_deallocate    
CONTAINS
SUBROUTINE get_dos_s(dos_out)
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: dos_out
    INTEGER :: i
    dos_out = 0.0_dp
    DO i = 1, nsurf
        IF ( AIMAG(green_s(i,i)) < 0) THEN
            WRITE(*,*) "negative DOS: green_s"
            ! STOP
        END IF
        dos_out = dos_out + AIMAG(green_s(i,i)) / PI
    END DO
    RETURN
END SUBROUTINE

! SUBROUTINE get_dos_up(dos_out)
!     IMPLICIT NONE
!     REAL(dp), INTENT(OUT) :: dos_out
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
!     REAL(dp), INTENT(OUT) :: dos_out
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
!     REAL(dp), INTENT(OUT) :: spin_out(3)
!     COMPLEX(dp), ALLOCATABLE :: spnk(:,:)
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
!     REAL(dp), INTENT(OUT) :: spin_out(3)
!     COMPLEX(dp), ALLOCATABLE :: spnk(:,:)
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

SUBROUTINE get_dos_b(dos_out)
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: dos_out
    INTEGER :: i
    dos_out = 0.0_dp
    DO i = 1, nbulk
        IF ( AIMAG(green_b(i,i)) < 0) THEN
            WRITE(*,*) "negative DOS: green_b"
            ! STOP
        END IF
        dos_out = dos_out + AIMAG(green_b(i,i)) / PI
    END DO
END SUBROUTINE

SUBROUTINE get_dos_nlayer(nlayer, dos_out)
    INTEGER, INTENT(IN) :: nlayer
    complex(dp), INTENT(OUT) :: dos_out(64)
    INTEGER :: i, il, j, cnt
    IF (nlayer .lt. 1) RETURN
    ! recurrence to calculate green_ilayer
    CALL set_green_layer(nlayer)
    dos_out = cmplx_0
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
    REAL(dp), INTENT(OUT) :: spin_out(3)
    COMPLEX(dp), ALLOCATABLE :: spnk(:,:)
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
    COMPLEX(dp), ALLOCATABLE :: green_layer_temp(:,:,:)
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:) :: temp_mat1, temp_mat2, temp_mat3
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
            CALL ZGEMM('N','C',nbulk,nsurf,nbulk,CMPLX_1, teff_mat, nbulk, h01, nsurf,CMPLX_0, temp_mat1, nbulk)
            CALL ZGEMM('N','N',nbulk,nsurf,nsurf,CMPLX_1, temp_mat1, nbulk, green_s, nsurf,CMPLX_0, temp_mat2, nbulk)
            CALL ZGEMM('N','N',nbulk,nbulk,nsurf,CMPLX_1, temp_mat2, nbulk, h01, nsurf,CMPLX_0, temp_mat3, nbulk)
            CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1, temp_mat3, nbulk, seff_mat, nbulk,CMPLX_1, green_layer(:,:,1), nbulk)
            DEALLOCATE(temp_mat1)
            DEALLOCATE(temp_mat2)
            DEALLOCATE(temp_mat3)
        ELSE IF((il==1) .AND. .NOT.(isslab)) THEN
            ! temp_mat1 = t_mat * green_s
            ! green_layer(:,:,1) = green_s + temp_mat1 * s_mat
            ALLOCATE(temp_mat1(nbulk,nbulk))
            CALL ZCOPY(nbulk*nbulk, green_s, 1, green_layer(:,:,1), 1)
            CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1, t_mat, nbulk, green_s, nbulk,CMPLX_0, temp_mat1, nbulk)
            CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1, temp_mat1, nbulk, s_mat, nbulk,CMPLX_1, green_layer(:,:,1), nbulk)
            DEALLOCATE(temp_mat1)
        ELSE IF ((il>1) .AND. isslab) THEN
            ! green_layer(:,:,il) = teff_mat + ( (teff_mat * h12^dagger) * green_layer(:,:,il-1)) * s_mat
            ALLOCATE(temp_mat1(nbulk,nbulk))
            ALLOCATE(temp_mat2(nbulk,nbulk))
            CALL ZCOPY(nbulk*nbulk, teff_mat, 1, green_layer(:,:,il), 1)
            CALL ZGEMM('N','C',nbulk,nbulk,nbulk,CMPLX_1, teff_mat, nbulk, h12, nbulk,CMPLX_0, temp_mat1, nbulk)
            CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1, temp_mat1, nbulk, green_layer(:,:,il-1), nbulk,CMPLX_0, temp_mat2, nbulk)
            CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1, temp_mat2, nbulk, s_mat, nbulk,CMPLX_1, green_layer(:,:,il), nbulk)
            DEALLOCATE(temp_mat1)
            DEALLOCATE(temp_mat2)
        ELSE ! il > 1 and .NOT. isslab
            ! green_layer(:,:,il) = green_s + ( (t_mat * green_layer(:,:,il-1)) * s_mat )
            ALLOCATE(temp_mat1(nbulk,nbulk))
            CALL ZCOPY(nbulk*nbulk, green_s, 1, green_layer(:,:,il), 1)
            CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1, t_mat, nbulk, green_layer(:,:,il-1), nbulk,CMPLX_0, temp_mat1, nbulk)
            CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1, temp_mat1, nbulk, s_mat, nbulk,CMPLX_1, green_layer(:,:,il), nbulk)
            DEALLOCATE(temp_mat1)
        END IF
    END DO
END SUBROUTINE

SUBROUTINE pp_green_deallocate_green_layer()
    IMPLICIT NONE
    IF(ALLOCATED(green_layer)) DEALLOCATE(green_layer)
END SUBROUTINE

SUBROUTINE set_transfer_mat()
    COMPLEX(dp), ALLOCATABLE :: temp_mat(:,:)
    IF (isslab) THEN
        ! t_mat = green_s1 * h12^dagger
        CALL ZGEMM('N','C',nbulk,nbulk,nbulk,CMPLX_1,green_s1,nbulk,h12,nbulk,CMPLX_0,t_mat,nbulk)
        ! s_mat = h12 * green_s1
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1,h12,nbulk,green_s1,nbulk,CMPLX_0,s_mat,nbulk)
        
        ALLOCATE(temp_mat(nbulk,nbulk))
        ! temp_mat = h11 + h12*t_mat
        ! teff_mat = (omega - temp_mat)^-1
        CALL ZCOPY(nbulk*nbulk, h11, 1, temp_mat, 1)
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1, h12, nbulk, t_mat, nbulk,CMPLX_1, temp_mat, nbulk)
        CALL inv_omega_minus_mat(nbulk, temp_mat, omega, teff_mat, 'set_transfer_mat for teff_mat')
        ! temp_mat = h11 + s_mat*h12^dagger
        ! seff_mat = (omega - temp_mat)^-1
        CALL ZCOPY(nbulk*nbulk, h11, 1, temp_mat, 1)
        CALL ZGEMM('N','C',nbulk,nbulk,nbulk,CMPLX_1, s_mat, nbulk, h12, nbulk,CMPLX_1, temp_mat, nbulk)
        CALL inv_omega_minus_mat(nbulk, temp_mat, omega, seff_mat, 'set_transfer_mat for teff_mat')
        DEALLOCATE(temp_mat)
    ELSE ! .NOT. isslab
        ! t_mat = green_s * h12^dagger
        CALL ZGEMM('N','C',nbulk,nbulk,nbulk,CMPLX_1,green_s,nbulk,h12,nbulk,CMPLX_0,t_mat,nbulk)
        ! s_mat = h12 * green_s
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1,h12,nbulk,green_s,nbulk,CMPLX_0,s_mat,nbulk)
    ENDIF
END SUBROUTINE

SUBROUTINE pp_green_deallocate()
    IMPLICIT NONE
    INTEGER :: ierr
    DEALLOCATE(green_s, stat=ierr)
    DEALLOCATE(green_b, stat=ierr)
    IF(ALLOCATED(t_mat)) DEALLOCATE(t_mat,stat=ierr)
    IF(ALLOCATED(s_mat)) DEALLOCATE(s_mat,stat=ierr)
    IF(ALLOCATED(teff_mat)) DEALLOCATE(teff_mat,stat=ierr)
    IF(ALLOCATED(seff_mat)) DEALLOCATE(seff_mat,stat=ierr)
    if (ierr /= 0) CALL io_error('Error in deallocation in pp_green_allocate')
END SUBROUTINE


SUBROUTINE pp_green_setup()
    ! allocate variables used in each iteration
    ! read seedname_spnr.dat for the required nr3=0, and save
    ! them in spnr, srvec, sndegen, snrpts
    IMPLICIT NONE
    INTEGER :: ierr
    ! allocation
    ALLOCATE(green_s(nsurf,nsurf), stat=ierr)
    ALLOCATE(green_b(nbulk,nbulk), stat=ierr)
    ALLOCATE(t_mat(nbulk,nbulk), stat=ierr)
    ALLOCATE(s_mat(nbulk,nbulk), stat=ierr)
    IF (isslab) ALLOCATE(green_s1(nbulk,nbulk), stat=ierr)
    IF (isslab) ALLOCATE(teff_mat(nbulk,nbulk))
    IF (isslab) ALLOCATE(seff_mat(nbulk,nbulk))
    if (ierr /= 0) CALL io_error('Error in allocation in pp_green_allocate')
    IF (isspin) CALL read_spin(0, spnr, srvec, sndegen, snrpts)
END SUBROUTINE

subroutine read_spin(nr3,spnr_nr3,rvec_nr3,ndegen_nr3,irpts_nr3)
    ! from seedname_spnr.dat, parse spin elements such that
    ! lattice vector along dir=3 is nr3
    implicit none
    integer, intent(IN) :: nr3
    complex(dp), allocatable, intent(out) :: spnr_nr3(:,:,:,:)
    real(dp), allocatable, intent(out)  :: rvec_nr3(:,:), ndegen_nr3(:)
    integer, intent(out) :: irpts_nr3

    integer :: irpts,iw,ni,nj,ioerror,nrpts, ir1,ir2,ir3, num_spnr_wann
    real(dp) :: sp1r, sp1i, sp2r, sp2i, sp3r, sp3i
    complex(dp), allocatable :: spnr_nr3_tmp(:,:,:,:)
    real(dp), allocatable :: ndegen(:), rvec_nr3_tmp(:,:), ndegen_nr3_tmp(:)
    character(80) :: buffer
    logical :: ir3_checked

    ioerror=0
    irpts_nr3 = 0

    open(unit=10,file=trim(seedname)//'_spnr.dat',action='read',iostat=ioerror)
    if (ioerror/=0) CALL io_error('Error opening '//trim(seedname)//'_spnr.dat : ')

    read(10,'(80a)',iostat=ioerror)buffer
    read(10,'(80a)',iostat=ioerror)buffer
    read(buffer, *) num_spnr_wann
    read(10,'(80a)',iostat=ioerror)buffer
    read(buffer, *) nrpts
    IF(is_root) write(*,*) 'begin reading '//trim(seedname)//'_spnr.dat'
    IF(is_root) write(*,*) 'num_spnr_wann = ', num_spnr_wann
    IF(is_root) write(*,*) 'nrpts = ', nrpts

    allocate(ndegen(nrpts))
    allocate(spnr_nr3_tmp(num_spnr_wann,num_spnr_wann,nrpts,3))
    allocate(rvec_nr3_tmp(2,nrpts))
    allocate(ndegen_nr3_tmp(nrpts))
    
    do irpts = 1, (nrpts-1)/15
        read(10, *) ndegen((irpts-1)*15+1:(irpts-1)*15+15)
    end do
    read(10, *) ndegen(((nrpts-1)/15)*15+1:nrpts)
    do irpts = 1, nrpts
        ir3_checked = .false.
        do iw=1,num_spnr_wann**2
            read (10,*,iostat=ioerror) ir1,ir2,ir3, ni,nj,&
                sp1r, sp1i, sp2r, sp2i, sp3r, sp3i
            if (ioerror/=0) CALL io_error('Error reading file: '//trim(seedname)//'_spnr.dat : ')
            if (ir3 == nr3) then
                if (.not. ir3_checked) then
                    irpts_nr3 = irpts_nr3 + 1
                    ir3_checked = .true.
                end if
                spnr_nr3_tmp(ni,nj,irpts_nr3,1) = CMPLX(sp1r, sp1i, dp)
                spnr_nr3_tmp(ni,nj,irpts_nr3,2) = CMPLX(sp2r, sp2i, dp)
                spnr_nr3_tmp(ni,nj,irpts_nr3,3) = CMPLX(sp3r, sp3i, dp)
                rvec_nr3_tmp(:,irpts_nr3) = (/ ir1, ir2 /)
                ndegen_nr3_tmp(irpts_nr3) = ndegen(irpts)
            end if
        enddo
    enddo

    ! change size of array to known(read) size
    allocate(spnr_nr3(num_spnr_wann,num_spnr_wann,irpts_nr3,3))
    allocate(rvec_nr3(2,irpts_nr3))
    allocate(ndegen_nr3(irpts_nr3))
    spnr_nr3 = spnr_nr3_tmp(:,:,1:irpts_nr3,:)
    rvec_nr3 = rvec_nr3_tmp(:,1:irpts_nr3)
    ndegen_nr3 = ndegen_nr3_tmp(1:irpts_nr3)
    deallocate(spnr_nr3_tmp)
    deallocate(rvec_nr3_tmp)
    deallocate(ndegen_nr3_tmp)
    deallocate(ndegen)

    IF(is_root) write(*,'("Spin_",I1," shape : ",4I5)') nr3, shape(spnr_nr3)
    IF(is_root) write(*,*) 'end reading '//trim(seedname)//'_spnr.dat'
    CLOSE(10)
end subroutine

END MODULE
