!------------------------------------------------------------------------
MODULE hamiltonian
!------------------------------------------------------------------------
!! Module for real space hamintonian
!------------------------------------------------------------------------
  USE comms, ONLY : DP, is_root, io_error
  USE parameters, ONLY : is_ideal_surf, hr_stitching, nsurf, nbulk
  !
  IMPLICIT NONE
  !
  SAVE
  PRIVATE
  !
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: h00(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: h01(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: h11(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: h12(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: spn00(:,:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: spn01(:,:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: spn11(:,:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: spn12(:,:,:)
  ! INTEGER, ALLOCATABLE, PUBLIC :: ind_surf(:)
  ! INTEGER, ALLOCATABLE, PUBLIC :: ind_bulk(:)
  !
  ! private variables for saving values parsed from seedname_hr.dat
  COMPLEX(DP), ALLOCATABLE :: hr0s(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: hr0(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: hr1(:,:,:)
  INTEGER, ALLOCATABLE :: rvec0s(:,:)
  INTEGER, ALLOCATABLE :: rvec0(:,:)
  INTEGER, ALLOCATABLE :: rvec1(:,:)
  INTEGER, ALLOCATABLE :: ndegen0s(:)
  INTEGER, ALLOCATABLE :: ndegen0(:)
  INTEGER, ALLOCATABLE :: ndegen1(:)
  INTEGER :: nrpts0s, nrpts0, nrpts1
  ! private variables for saving values parsed from seedname_spnr.dat
  COMPLEX(DP), ALLOCATABLE :: spnr0s(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: spnr0(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: spnr1(:,:,:,:)
  INTEGER, ALLOCATABLE :: srvec0s(:,:)
  INTEGER, ALLOCATABLE :: srvec0(:,:)
  INTEGER, ALLOCATABLE :: srvec1(:,:)
  INTEGER, ALLOCATABLE :: sndegen0s(:)
  INTEGER, ALLOCATABLE :: sndegen0(:)
  INTEGER, ALLOCATABLE :: sndegen1(:)
  INTEGER :: snrpts0s, snrpts0, snrpts1
  !
  COMPLEX(DP) :: bulk_shift
  !! shift of bulk and slab energy (add to bulk onsite)
  !
  INTEGER :: num_hr_wann ! dimension of hamiltonian in seedname_hr.dat
  !
  PUBLIC :: hamiltonian_setup, hamiltonian_tb_to_k
CONTAINS
!------------------------------------------------------------------------
SUBROUTINE hamiltonian_setup()
!------------------------------------------------------------------------
!! ALLOCATE required variables, and read seedname_hr.dat for the required
!! nr3, and save the information in hr*, rvec* ndegen*, nrpts*
!------------------------------------------------------------------------
  USE parameters, ONLY : isspin, seedname, bulk_rz, ind_1
  IMPLICIT NONE
  INTEGER :: i, ir, ierr
  !
  ! allocation
  IF (.NOT. is_ideal_surf) THEN
    ALLOCATE(h00(nsurf, nsurf), STAT=ierr)
    ALLOCATE(h01(nsurf, nbulk), STAT=ierr)
    ALLOCATE(h11(nbulk, nbulk), STAT=ierr)
    ALLOCATE(h12(nbulk, nbulk), STAT=ierr)
    IF (isspin) THEN
      ALLOCATE(spn00(nsurf, nsurf, 3), STAT=ierr)
      ALLOCATE(spn01(nsurf, nbulk, 3), STAT=ierr)
      ALLOCATE(spn11(nbulk, nbulk, 3), STAT=ierr)
      ALLOCATE(spn12(nbulk, nbulk, 3), STAT=ierr)
    ENDIF
  ELSE
    ALLOCATE(h11(nbulk, nbulk), STAT=ierr)
    ALLOCATE(h12(nbulk, nbulk), STAT=ierr)
    IF (isspin) THEN
      ALLOCATE(spn11(nbulk, nbulk, 3), STAT=ierr)
      ALLOCATE(spn12(nbulk, nbulk, 3), STAT=ierr)
    ENDIF
  ENDIF
  IF (ierr /= 0) CALL io_error('Error during allocation in hamiltonian_setup')
  !
  ! Dummy allocation to avoid passing unallocated arrays as arguments
  ALLOCATE(hr0(1,1,1))
  ALLOCATE(rvec0(1,1))
  ALLOCATE(ndegen0(1))
  ALLOCATE(spnr0(1,1,1,1))
  ALLOCATE(srvec0(1,1))
  ALLOCATE(sndegen0(1))
  IF (is_ideal_surf .OR. (.NOT. is_ideal_surf .AND. hr_stitching)) THEN
    ALLOCATE(hr1(1,1,1))
    ALLOCATE(rvec1(1,1))
    ALLOCATE(ndegen1(1))
    ALLOCATE(spnr1(1,1,1,1))
    ALLOCATE(srvec1(1,1))
    ALLOCATE(sndegen1(1))
  ENDIF
  IF (.NOT. is_ideal_surf .AND. hr_stitching) THEN
    ALLOCATE(hr0s(1,1,1))
    ALLOCATE(rvec0s(1,1))
    ALLOCATE(ndegen0s(1))
    ALLOCATE(spnr0s(1,1,1,1))
    ALLOCATE(srvec0s(1,1))
    ALLOCATE(sndegen0s(1))
  ENDIF
  !
  ! reading seedname_hr.dat
  IF (.NOT. is_ideal_surf) THEN
    IF (hr_stitching) THEN
      CALL read_hr_dat(0, hr0, rvec0, ndegen0, nrpts0, TRIM(seedname)//".bulk")
      CALL read_hr_dat(bulk_rz, hr1, rvec1, ndegen1, nrpts1, TRIM(seedname)//".bulk")
      CALL read_hr_dat(0, hr0s, rvec0s, ndegen0s, nrpts0s, TRIM(seedname)//".slab")
      IF (isspin) THEN
        CALL read_spnr_dat(0, spnr0, srvec0, sndegen0, snrpts0, TRIM(seedname)//".bulk")
        CALL read_spnr_dat(bulk_rz, spnr1, rvec1, ndegen1, nrpts1, TRIM(seedname)//".bulk")
        CALL read_spnr_dat(0, spnr0s, srvec0s, sndegen0s, snrpts0s, TRIM(seedname)//".slab")
      ENDIF
      !
      ! Calculate the required onsite energy shift
      ! bulk_shift = average(hr0s_onsite - hr0_onsite)
      bulk_shift = 0.0_dp
      DO ir = 1, nrpts0s
        IF (rvec0s(1,ir) == 0 .AND. rvec0s(2,ir) == 0) THEN
          DO i = 1, nbulk
            ! FIXME: remove ndegen0s(ir) because it is always 1
            IF (ndegen0s(ir) /= 1) CALL io_error('Error: ndegen0s(ir) /= 1 in hamiltonian_setup')
            bulk_shift = bulk_shift + hr0s(ind_1(i), ind_1(i), ir) / ndegen0s(ir)
          ENDDO
        ENDIF
      ENDDO
      DO ir = 1, nrpts0
        IF (rvec0(1,ir) == 0 .AND. rvec0(2,ir) == 0) THEN
          DO i = 1, nbulk
            ! FIXME: remove ndegen0(ir) because it is always 1
            IF (ndegen0(ir) /= 1) CALL io_error('Error: ndegen0(ir) /= 1 in hamiltonian_setup')
            bulk_shift = bulk_shift - hr0(i, i, ir) / ndegen0(ir)
          ENDDO
        ENDIF
      ENDDO
      bulk_shift = CMPLX(REAL(bulk_shift / nbulk, DP), 0.0_dp, DP)
      IF (is_root) WRITE(*,*) "bulk_shift = ", bulk_shift
      IF (is_root) WRITE(*,*) "Note: this may(should) be small if util_match is used to generate slab hr.dat file"
    ELSE ! .NOT. is_ideal_surf and .NOT. hr_stitching
      CALL read_hr_dat(0, hr0, rvec0, ndegen0, nrpts0, TRIM(seedname))
      IF (isspin) THEN
        CALL read_spnr_dat(0, spnr0, srvec0, sndegen0, snrpts0, TRIM(seedname))
      ENDIF
    ENDIF ! hr_stitching
  ELSE ! is_ideal_surf
    CALL read_hr_dat(0, hr0, rvec0, ndegen0, nrpts0, TRIM(seedname))
    CALL read_hr_dat(bulk_rz, hr1, rvec1, ndegen1, nrpts1, TRIM(seedname))
    IF (isspin) THEN
      CALL read_spnr_dat(0, spnr0, srvec0, sndegen0, snrpts0, TRIM(seedname))
      CALL read_spnr_dat(bulk_rz, spnr1, srvec1, sndegen1, snrpts1, TRIM(seedname))
    ENDIF
  ENDIF
!------------------------------------------------------------------------
END SUBROUTINE hamiltonian_setup
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE hamiltonian_tb_to_k(kx, ky)
!------------------------------------------------------------------------
!! sum the in-plane hamiltonian multiplied with the phase factor
!! to create the tight-binding hamiltonian, hij.
!------------------------------------------------------------------------
  USE comms, ONLY : k_operator
  USE parameters, ONLY : isspin, ind_0, ind_1, ind_2
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: kx
  !! k vector along x axis, in crystal coordinate
  REAL(DP), INTENT(IN) :: ky
  !! k vector along y axis, in crystal coordinate
  !
  COMPLEX(DP), ALLOCATABLE :: htemp(:,:)
  INTEGER :: i
  INTEGER :: is
  !! Spin direction index (1=x, 2=y, 3=z)
  !
  IF (.NOT. is_ideal_surf) THEN
    IF (hr_stitching) THEN
      ALLOCATE(htemp(num_hr_wann,num_hr_wann))
      CALL k_operator(nrpts0, hr0, rvec0, ndegen0, kx, ky, h11)
      CALL k_operator(nrpts1, hr1, rvec1, ndegen1, kx, ky, h12)
      CALL k_operator(nrpts0s, hr0s, rvec0s, ndegen0s, kx, ky, htemp)
      h00 = htemp(ind_0, ind_0)
      h01 = htemp(ind_0, ind_1)
      ! shift diagonal elements of bulk hamiltonian using bulk_shift
      ! to match bulk and slab energy reference
      DO i = 1, nsurf
        h00(i, i) = h00(i,i) - bulk_shift
      ENDDO
      !
      IF (isspin) THEN
        DO is = 1, 3
          CALL k_operator(snrpts0, spnr0(:,:,:,is), srvec0, sndegen0, kx, ky, spn11(:,:,is))
          CALL k_operator(snrpts1, spnr1(:,:,:,is), srvec1, sndegen1, kx, ky, spn12(:,:,is))
          CALL k_operator(snrpts0s, spnr0s(:,:,:,is), srvec0s, sndegen0s, kx, ky, htemp)
          spn00(:,:,is) = htemp(ind_0, ind_0)
          spn01(:,:,is) = htemp(ind_0, ind_1)
        ENDDO
      ENDIF ! isspin
      DEALLOCATE(htemp)
    ELSE ! .NOT. is_ideal_surf and hr_stitching
      ALLOCATE(htemp(num_hr_wann,num_hr_wann))
      CALL k_operator(nrpts0, hr0, rvec0, ndegen0, kx, ky, htemp)
      h00 = htemp(ind_0, ind_0)
      h01 = htemp(ind_0, ind_1)
      h11 = htemp(ind_1, ind_1)
      h12 = htemp(ind_1, ind_2)
      IF (isspin) THEN
        DO is = 1, 3
          CALL k_operator(snrpts0, spnr0(:,:,:,is), srvec0, sndegen0, kx, ky, htemp)
          spn00(:,:,is) = htemp(ind_0, ind_0)
          spn01(:,:,is) = htemp(ind_0, ind_1)
          spn11(:,:,is) = htemp(ind_1, ind_1)
          spn12(:,:,is) = htemp(ind_1, ind_2)
        ENDDO
      ENDIF
      DEALLOCATE(htemp)
    ENDIF ! hr_stitching
  ELSE ! is_ideal_surf
    CALL k_operator(nrpts0, hr0, rvec0, ndegen0, kx, ky, h11)
    CALL k_operator(nrpts1, hr1, rvec1, ndegen1, kx, ky, h12)
    IF (isspin) THEN
      DO is = 1, 3
        CALL k_operator(snrpts0, spnr0(:,:,:,is), srvec0, sndegen0, kx, ky, spn11(:,:,is))
        CALL k_operator(snrpts1, spnr1(:,:,:,is), srvec1, sndegen1, kx, ky, spn12(:,:,is))
      ENDDO
    ENDIF
  ENDIF
!------------------------------------------------------------------------
END SUBROUTINE hamiltonian_tb_to_k
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE read_hr_dat(nr3, hr_nr3, rvec_nr3, ndegen_nr3, nrpts_nr3, prefix)
!------------------------------------------------------------------------
!! Read seedname_hr.dat file, parse hamiltonian elements <m0|H|nR>
!! such that R(3) == nr3
!------------------------------------------------------------------------
  USE comms, ONLY : find_free_unit
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nr3
  COMPLEX(DP), ALLOCATABLE, INTENT(INOUT) :: hr_nr3(:,:,:)
  INTEGER, ALLOCATABLE, INTENT(INOUT) :: rvec_nr3(:,:)
  INTEGER, ALLOCATABLE, INTENT(INOUT) :: ndegen_nr3(:)
  INTEGER, INTENT(OUT) :: nrpts_nr3
  CHARACTER(*), INTENT(IN) :: prefix
  !
  INTEGER :: iun
  INTEGER :: irpts, iw, ni, nj, ierr, nrpts, ir1, ir2, ir3
  REAL(DP) :: hr_re, hr_im
  COMPLEX(DP), ALLOCATABLE :: hr_nr3_tmp(:,:,:)
  INTEGER, ALLOCATABLE :: ndegen(:), rvec_nr3_tmp(:,:), ndegen_nr3_tmp(:)
  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=80) :: buffer
  !
  iun = find_free_unit()
  filename = TRIM(prefix)//'_hr.dat'
  OPEN(UNIT=iun, FILE=TRIM(filename), ACTION='read', IOSTAT=ierr)
  IF (ierr /= 0) CALL io_error('Error opening ' // TRIM(filename))
  !
  READ(iun, '(80a)', IOSTAT=ierr) buffer
  READ(iun, '(80a)', IOSTAT=ierr) buffer
  READ(buffer, *) num_hr_wann
  READ(iun, '(80a)', IOSTAT=ierr) buffer
  READ(buffer, *) nrpts
  IF(is_root) WRITE(*,*) 'begin reading ' // TRIM(filename)
  IF(is_root) WRITE(*,*) 'rvec_z("nr3") = ', nr3
  IF(is_root) WRITE(*,*) 'num_hr_wann = ', num_hr_wann
  IF(is_root) WRITE(*,*) 'nrpts = ', nrpts
  !
  ALLOCATE(ndegen(nrpts))
  ALLOCATE(hr_nr3_tmp(num_hr_wann, num_hr_wann, nrpts))
  ALLOCATE(rvec_nr3_tmp(2, nrpts))
  ALLOCATE(ndegen_nr3_tmp(nrpts))
  !
  DO irpts = 1, (nrpts-1) / 15
    READ(iun, *, IOSTAT=ierr) ndegen((irpts-1)*15+1:(irpts-1)*15+15)
  ENDDO
  READ(iun, *, IOSTAT=ierr) ndegen(((nrpts-1)/15)*15+1:nrpts)
  IF (ierr/=0) CALL io_error('Error reading file: ' // filename)
  !
  nrpts_nr3 = 0
  DO irpts = 1, nrpts
    DO iw = 1, num_hr_wann**2
      READ (iun,'(5I5,2F12.6)', IOSTAT=ierr) ir1, ir2, ir3, ni, nj, hr_re, hr_im
      IF (ierr/=0) CALL io_error('Error reading file: ' // filename)
      IF (ir3 == nr3) THEN
        IF (iw == 1) THEN
          nrpts_nr3 = nrpts_nr3 + 1
          rvec_nr3_tmp(1, nrpts_nr3) = ir1
          rvec_nr3_tmp(2, nrpts_nr3) = ir2
          ndegen_nr3_tmp(nrpts_nr3) = ndegen(irpts)
        ENDIF
        hr_nr3_tmp(ni, nj, nrpts_nr3) = CMPLX(hr_re, hr_im, DP)
      ENDIF
    ENDDO
  ENDDO
  CLOSE(iun)
  !
  ! change size of array to known(read) size
  DEALLOCATE(hr_nr3)
  DEALLOCATE(rvec_nr3)
  DEALLOCATE(ndegen_nr3)
  ALLOCATE(hr_nr3(num_hr_wann, num_hr_wann, nrpts_nr3))
  ALLOCATE(rvec_nr3(2, nrpts_nr3))
  ALLOCATE(ndegen_nr3(nrpts_nr3))
  hr_nr3 = hr_nr3_tmp(:, :, 1:nrpts_nr3)
  rvec_nr3 = rvec_nr3_tmp(:, 1:nrpts_nr3)
  ndegen_nr3 = ndegen_nr3_tmp(1:nrpts_nr3)
  DEALLOCATE(hr_nr3_tmp)
  DEALLOCATE(rvec_nr3_tmp)
  DEALLOCATE(ndegen_nr3_tmp)
  DEALLOCATE(ndegen)
  !
  IF (is_root) write(*,'("Hamiltonian_",I2," shape : ",3I5)') nr3, SHAPE(hr_nr3)
  IF (is_root) write(*,*) 'end reading ' // TRIM(filename)
!------------------------------------------------------------------------
END SUBROUTINE read_hr_dat
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE read_spnr_dat(nr3, spnr_nr3, rvec_nr3, ndegen_nr3, nrpts_nr3, prefix)
!! Read seedname_spnr.dat. Save spin elements <m0|S|nR> with R_z = nr3
!! to variables spnr_nr3.
  USE comms, ONLY : find_free_unit
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nr3
  COMPLEX(DP), ALLOCATABLE, INTENT(INOUT) :: spnr_nr3(:,:,:,:)
  INTEGER, ALLOCATABLE, INTENT(INOUT) :: rvec_nr3(:,:)
  INTEGER, ALLOCATABLE, INTENT(INOUT) :: ndegen_nr3(:)
  INTEGER, INTENT(OUT) :: nrpts_nr3
  CHARACTER(*), INTENT(IN) :: prefix
  !
  INTEGER :: iun
  INTEGER :: irpts,iw,ni,nj,ierr,nrpts, ir1,ir2,ir3, num_spnr_wann
  REAL(DP) :: sp1r, sp1i, sp2r, sp2i, sp3r, sp3i
  COMPLEX(DP), ALLOCATABLE :: spnr_nr3_tmp(:,:,:,:)
  INTEGER, ALLOCATABLE :: ndegen(:), rvec_nr3_tmp(:,:), ndegen_nr3_tmp(:)
  CHARACTER(LEN=80) :: buffer
  LOGICAL :: ir3_checked
  CHARACTER(LEN=256) :: filename
  !
  iun = find_free_unit()
  filename = TRIM(prefix) // '_spnr.dat'
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
  nrpts_nr3 = 0
  DO irpts = 1, nrpts
    ir3_checked = .false.
    DO iw = 1,num_spnr_wann**2
      READ(iun, *, IOSTAT=ierr) ir1, ir2, ir3, ni, nj, &
                                sp1r, sp1i, sp2r, sp2i, sp3r, sp3i
      IF (ierr /= 0) CALL io_error('Error reading file '// TRIM(filename))
      ! Choose lattice vectors with R_3 = nr3
      IF (ir3 == nr3) THEN
        IF (.NOT. ir3_checked) THEN
          nrpts_nr3 = nrpts_nr3 + 1
          ir3_checked = .true.
        ENDIF
        spnr_nr3_tmp(ni, nj, nrpts_nr3, 1) = CMPLX(sp1r, sp1i, DP)
        spnr_nr3_tmp(ni, nj, nrpts_nr3, 2) = CMPLX(sp2r, sp2i, DP)
        spnr_nr3_tmp(ni, nj, nrpts_nr3, 3) = CMPLX(sp3r, sp3i, DP)
        rvec_nr3_tmp(1, nrpts_nr3) = ir1
        rvec_nr3_tmp(2, nrpts_nr3) = ir2
        ndegen_nr3_tmp(nrpts_nr3) = ndegen(irpts)
      ENDIF
    ENDDO
  ENDDO
  CLOSE(iun)
  !
  ! Change size of array to the read size
  ! This is required because the number of R vectors with R3=nr3 is not known
  ! before reading.
  DEALLOCATE(spnr_nr3)
  DEALLOCATE(rvec_nr3)
  DEALLOCATE(ndegen_nr3)
  ALLOCATE(spnr_nr3(num_spnr_wann, num_spnr_wann, nrpts_nr3, 3))
  ALLOCATE(rvec_nr3(2, nrpts_nr3))
  ALLOCATE(ndegen_nr3(nrpts_nr3))
  spnr_nr3 = spnr_nr3_tmp(:, :, 1:nrpts_nr3, :)
  rvec_nr3 = rvec_nr3_tmp(:, 1:nrpts_nr3)
  ndegen_nr3 = ndegen_nr3_tmp(1:nrpts_nr3)
  DEALLOCATE(spnr_nr3_tmp)
  DEALLOCATE(rvec_nr3_tmp)
  DEALLOCATE(ndegen_nr3_tmp)
  DEALLOCATE(ndegen)
  !
  IF(is_root) write(*,'("Spin_",I1," shape : ",4I5)') nr3, shape(spnr_nr3)
  IF(is_root) write(*,*) 'end reading ' // TRIM(filename)
  !
END SUBROUTINE read_spnr_dat
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE hamiltonian_deallocate()
!------------------------------------------------------------------------
!! deallocate hamiltonian variables
!------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: ierr
  DEALLOCATE(h11, STAT=ierr)
  DEALLOCATE(h12, STAT=ierr)
  IF(ALLOCATED(h00)) DEALLOCATE(h00,STAT=ierr)
  IF(ALLOCATED(h01)) DEALLOCATE(h01,STAT=ierr)
  IF(ALLOCATED(hr0)) DEALLOCATE(hr0,STAT=ierr)
  IF(ALLOCATED(hr1)) DEALLOCATE(hr1,STAT=ierr)
  IF(ALLOCATED(rvec0)) DEALLOCATE(rvec0,STAT=ierr)
  IF(ALLOCATED(rvec1)) DEALLOCATE(rvec1,STAT=ierr)
  IF(ALLOCATED(ndegen0)) DEALLOCATE(ndegen0,STAT=ierr)
  IF(ALLOCATED(ndegen1)) DEALLOCATE(ndegen1,STAT=ierr)
  IF (ierr /= 0) CALL io_error('Error in allocation in hamiltonian_ALLOCATE_slab')
!------------------------------------------------------------------------
END SUBROUTINE hamiltonian_deallocate
!------------------------------------------------------------------------
END MODULE hamiltonian
