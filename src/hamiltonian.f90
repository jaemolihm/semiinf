!------------------------------------------------------------------------
MODULE hamiltonian
!------------------------------------------------------------------------
!! Module for real space hamintonian
!------------------------------------------------------------------------
  USE comms, ONLY : DP, is_root, io_error
  USE parameters, ONLY : isslab, hr_stitching, nsurf, nbulk
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
  ! INTEGER, ALLOCATABLE, PUBLIC :: ind_surf(:)
  ! INTEGER, ALLOCATABLE, PUBLIC :: ind_bulk(:)
  !
  ! private variables for saving values parsed from seedname_hr.dat
  COMPLEX(DP), ALLOCATABLE :: hr0s(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: hr0(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: hr1(:,:,:)
  REAL(DP), ALLOCATABLE :: rvec0s(:,:)
  REAL(DP), ALLOCATABLE :: rvec0(:,:)
  REAL(DP), ALLOCATABLE :: rvec1(:,:)
  REAL(DP), ALLOCATABLE :: ndegen0s(:)
  REAL(DP), ALLOCATABLE :: ndegen0(:)
  REAL(DP), ALLOCATABLE :: ndegen1(:)
  COMPLEX(DP) :: bulk_shift ! shift of bulk and slab energy (add to bulk onsite)
  INTEGER :: nrpts0s, nrpts0, nrpts1
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
  USE parameters, ONLY : seedname, bulk_rz, ind_1
  IMPLICIT NONE
  INTEGER :: i, ir, ierr
  !
  ! allocation
  IF (isslab) THEN
    ALLOCATE(h00(nsurf,nsurf), STAT=ierr)
    ALLOCATE(h01(nsurf,nbulk), STAT=ierr)
    ALLOCATE(h11(nbulk,nbulk), STAT=ierr)
    ALLOCATE(h12(nbulk,nbulk), STAT=ierr)
  ELSE
    ALLOCATE(h11(nbulk,nbulk), STAT=ierr)
    ALLOCATE(h12(nbulk,nbulk), STAT=ierr)
  ENDIF
  IF (ierr /= 0) CALL io_error('Error during allocation in hamiltonian_setup')
  !
  ! reading seedname_hr.dat
  IF (isslab) THEN
    IF (hr_stitching) THEN
      CALL read_hamiltonian(0, hr0, rvec0, ndegen0, nrpts0, TRIM(seedname)//".bulk")
      CALL read_hamiltonian(bulk_rz, hr1, rvec1, ndegen1, nrpts1, TRIM(seedname)//".bulk")
      CALL read_hamiltonian(0, hr0s, rvec0s, ndegen0s, nrpts0s, TRIM(seedname)//".slab")
      seedname = TRIM(seedname) // ".slab"
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
    ELSE ! isslab and .NOT. hr_stitching)
      CALL read_hamiltonian(0, hr0, rvec0, ndegen0, nrpts0, TRIM(seedname))
    ENDIF ! hr_stitching
  ELSE ! .NOT. isslab
    CALL read_hamiltonian(0, hr0, rvec0, ndegen0, nrpts0, TRIM(seedname))
    CALL read_hamiltonian(bulk_rz, hr1, rvec1, ndegen1, nrpts1, TRIM(seedname))
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
  USE parameters, ONLY : ind_0, ind_1, ind_2
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: kx
  !! k vector along x axis, in crystal coordinate
  REAL(DP), INTENT(IN) :: ky
  !! k vector along y axis, in crystal coordinate
  !
  COMPLEX(DP), ALLOCATABLE :: htemp(:,:)
  INTEGER :: i
  !
  IF (isslab) THEN
    IF (hr_stitching) THEN
      CALL k_operator(nrpts0, hr0, rvec0, ndegen0, kx, ky, h11)
      CALL k_operator(nrpts1, hr1, rvec1, ndegen1, kx, ky, h12)
      ALLOCATE(htemp(num_hr_wann,num_hr_wann))
      CALL k_operator(nrpts0s, hr0s, rvec0s, ndegen0s, kx, ky, htemp)
      h00 = htemp(ind_0, ind_0)
      h01 = htemp(ind_0, ind_1)
      ! shift diagonal elements of bulk hamiltonian using bulk_shift
      ! to match bulk and slab energy reference
      DO i=1,nsurf
        h00(i,i) = h00(i,i) - bulk_shift
      ENDDO
      DEALLOCATE(htemp)
    ELSE ! isslab and hr_stitching
      ALLOCATE(htemp(num_hr_wann,num_hr_wann))
      CALL k_operator(nrpts0, hr0, rvec0, ndegen0, kx, ky, htemp)
      h00 = htemp(ind_0, ind_0)
      h01 = htemp(ind_0, ind_1)
      h11 = htemp(ind_1, ind_1)
      h12 = htemp(ind_1, ind_2)
      DEALLOCATE(htemp)
    ENDIF ! hr_stitching
  ELSE !.NOT. isslab
    CALL k_operator(nrpts0, hr0, rvec0, ndegen0, kx, ky, h11)
    CALL k_operator(nrpts1, hr1, rvec1, ndegen1, kx, ky, h12)
  ENDIF
!------------------------------------------------------------------------
END SUBROUTINE hamiltonian_tb_to_k
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE read_hamiltonian(nr3, hr_nr3, rvec_nr3, ndegen_nr3, irpts_nr3, prefix)
!------------------------------------------------------------------------
!! Read seedname_hr.dat file, parse hamiltonian elements <m0|H|nR>
!! such that R(3) == nr3
!------------------------------------------------------------------------
  USE comms, ONLY : find_free_unit
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nr3
  COMPLEX(DP), ALLOCATABLE, INTENT(INOUT) :: hr_nr3(:,:,:)
  REAL(DP), ALLOCATABLE, INTENT(INOUT) :: rvec_nr3(:,:)
  REAL(DP), ALLOCATABLE, INTENT(INOUT) :: ndegen_nr3(:)
  INTEGER, INTENT(OUT) :: irpts_nr3
  CHARACTER(*), INTENT(IN) :: prefix
  !
  INTEGER :: iun
  INTEGER :: irpts, iw, ni, nj, ierr, nrpts, ir1, ir2, ir3
  REAL(DP) :: hr_re, hr_im
  COMPLEX(DP), ALLOCATABLE :: hr_nr3_tmp(:,:,:)
  REAL(DP), ALLOCATABLE :: ndegen(:), rvec_nr3_tmp(:,:), ndegen_nr3_tmp(:)
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
  irpts_nr3 = 0
  DO irpts = 1, nrpts
    DO iw = 1, num_hr_wann**2
      READ (iun,'(5I5,2F12.6)', IOSTAT=ierr) ir1, ir2, ir3, ni, nj, hr_re, hr_im
      IF (ierr/=0) CALL io_error('Error reading file: ' // filename)
      IF (ir3 == nr3) THEN
        IF (iw == 1) THEN
          irpts_nr3 = irpts_nr3 + 1
          rvec_nr3_tmp(1, irpts_nr3) = ir1
          rvec_nr3_tmp(2, irpts_nr3) = ir2
          ndegen_nr3_tmp(irpts_nr3) = ndegen(irpts)
        ENDIF
        hr_nr3_tmp(ni, nj, irpts_nr3) = CMPLX(hr_re, hr_im, DP)
      ENDIF
    ENDDO
  ENDDO
  CLOSE(iun)
  !
  ! change size of array to known(read) size
  ALLOCATE(hr_nr3(num_hr_wann, num_hr_wann, irpts_nr3))
  ALLOCATE(rvec_nr3(2, irpts_nr3))
  ALLOCATE(ndegen_nr3(irpts_nr3))
  hr_nr3 = hr_nr3_tmp(:, :, 1:irpts_nr3)
  rvec_nr3 = rvec_nr3_tmp(:, 1:irpts_nr3)
  ndegen_nr3 = ndegen_nr3_tmp(1:irpts_nr3)
  DEALLOCATE(hr_nr3_tmp)
  DEALLOCATE(rvec_nr3_tmp)
  DEALLOCATE(ndegen_nr3_tmp)
  DEALLOCATE(ndegen)
  !
  IF (is_root) write(*,'("Hamiltonian_",I2," shape : ",3I5)') nr3, SHAPE(hr_nr3)
  IF (is_root) write(*,*) 'end reading ' // TRIM(filename)
!------------------------------------------------------------------------
END SUBROUTINE read_hamiltonian
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
