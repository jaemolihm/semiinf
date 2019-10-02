!------------------------------------------------------------------------
MODULE hamiltonian
!------------------------------------------------------------------------
!! Module for real space hamintonian
!------------------------------------------------------------------------
  USE comms
  USE parameters
  !
  IMPLICIT NONE
  !
  SAVE
  !
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: h00(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: h01(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: h11(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: h12(:,:)
  COMPLEX(DP), PUBLIC :: omega
  REAL(DP), PUBLIC :: kx
  REAL(DP), PUBLIC :: ky
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
  ! PUBLIC :: hamiltonian_set_hread, hamiltonian_set_hij
CONTAINS
!------------------------------------------------------------------------
SUBROUTINE hamiltonian_setup()
!------------------------------------------------------------------------
!! ALLOCATE required variables, and read seedname_hr.dat for the required
!! nr3, and save the information in hr*, rvec* ndegen*, nrpts*
!------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: i, ir, ierr
  !
  ! allocation
  IF (isslab) THEN
    ALLOCATE(h00(nsurf,nsurf), stat=ierr)
    ALLOCATE(h01(nsurf,nbulk), stat=ierr)
    ALLOCATE(h11(nbulk,nbulk), stat=ierr)
    ALLOCATE(h12(nbulk,nbulk), stat=ierr)
  ELSE
    ALLOCATE(h11(nbulk,nbulk), stat=ierr)
    ALLOCATE(h12(nbulk,nbulk), stat=ierr)
  ENDIF
  IF (ierr /= 0) CALL io_error('Error in allocation in hamiltonian_setup')
  !
  ! reading seedname_hr.dat
  IF (isslab_hamil) THEN
    CALL read_hamiltonian(0, hr0, rvec0, ndegen0, nrpts0)
  ELSE IF (isslab_match) THEN
    CALL read_hamiltonian(0, hr0, rvec0, ndegen0, nrpts0, trim(seedname)//".bulk")
    CALL read_hamiltonian(bulk_rz, hr1, rvec1, ndegen1, nrpts1, trim(seedname)//".bulk")
    CALL read_hamiltonian(0, hr0s, rvec0s, ndegen0s, nrpts0s, trim(seedname)//".slab")
    seedname = trim(seedname)//".slab"

    ! calculate the required onsite energy shift
    ! bulk_shift = average(hr0s_onsite - hr0_onsite)
    bulk_shift = 0.0_dp
    DO ir = 1, nrpts0s
      IF (rvec0s(1,ir)==0 .and. rvec0s(2,ir)==0) THEN
        DO i = 1, nbulk
          bulk_shift = bulk_shift + hr0s(ind_1(i), ind_1(i), ir) / ndegen0s(ir)
        ENDDO
      ENDIF
    ENDDO
    DO ir = 1, nrpts0
      IF (rvec0(1,ir)==0 .and. rvec0(2,ir)==0) THEN
        DO i = 1, nbulk
          bulk_shift = bulk_shift - hr0(i, i, ir) / ndegen0(ir)
        ENDDO
      ENDIF
    ENDDO
    bulk_shift = CMPLX(REAL(bulk_shift / nbulk, dp), 0.0_dp, dp)
    IF (is_root) WRITE(*,*) "bulk_shift = ", bulk_shift
    IF (is_root) WRITE(*,*) "Note: this may(should) be small if util_match is used to generate slab hr.dat file"
  ELSE
    CALL read_hamiltonian(0, hr0, rvec0, ndegen0, nrpts0)
    CALL read_hamiltonian(1, hr1, rvec1, ndegen1, nrpts1)
  ENDIF
!------------------------------------------------------------------------
END SUBROUTINE hamiltonian_setup
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE hamiltonian_tb_to_k()
!------------------------------------------------------------------------
!! sum the in-plane hamiltonian multiplied with the phase factor
!! to create the tight-binding hamiltonian, hij.
!------------------------------------------------------------------------
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE :: htemp(:,:)
  INTEGER :: i
  !
  IF (isslab_hamil) THEN
    ALLOCATE(htemp(num_hr_wann,num_hr_wann))
    CALL k_operator(nrpts0, hr0, rvec0, ndegen0, kx, ky, htemp)
    IF (isslab) h00 = htemp(ind_0, ind_0)
    IF (isslab) h01 = htemp(ind_0, ind_1)
    h11 = htemp(ind_1, ind_1)
    h12 = htemp(ind_1, ind_2)
    DEALLOCATE(htemp)
  ELSE IF (isbulk_add_onsite) THEN
    CALL k_operator(nrpts0, hr0, rvec0, ndegen0, kx, ky, h11)
    CALL k_operator(nrpts1, hr1, rvec1, ndegen1, kx, ky, h12)
    h01 = h12
    h00 = h11
    do i = 1, nbulk
      h00(i,i) = h00(i,i) + bulk_onsite_energy(i)
    !     print*, i, bulk_onsite_energy(i)
    ENDDO
  ELSE IF (isslab_match) THEN
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
  ELSE !.NOT. isslab, .NOT. isslab_hamil, .NOT. isslab_match
    CALL k_operator(nrpts0, hr0, rvec0, ndegen0, kx, ky, h11)
    CALL k_operator(nrpts1, hr1, rvec1, ndegen1, kx, ky, h12)
  ENDIF
!------------------------------------------------------------------------
END SUBROUTINE hamiltonian_tb_to_k
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE read_hamiltonian(nr3,hr_nr3,rvec_nr3,ndegen_nr3,irpts_nr3,seedname_)
!------------------------------------------------------------------------
!! Read seedname_hr.dat file, parse hamiltonian elements such that
!! lattice vector along dir=3 is nr3
!------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nr3
  COMPLEX(DP), ALLOCATABLE, INTENT(INOUT) :: hr_nr3(:,:,:)
  REAL(DP), ALLOCATABLE, INTENT(INOUT) :: rvec_nr3(:,:)
  REAL(DP), ALLOCATABLE, INTENT(INOUT) :: ndegen_nr3(:)
  INTEGER, INTENT(OUT) :: irpts_nr3
  CHARACTER(*), INTENT(IN), optional :: seedname_

  CHARACTER(LEN=80) :: filename

  INTEGER :: irpts, iw, ni, nj, ioerror, nrpts, ir1, ir2, ir3
  REAL(DP) :: hr_re, hr_im
  COMPLEX(DP), ALLOCATABLE :: hr_nr3_tmp(:,:,:)
  REAL(DP), ALLOCATABLE :: ndegen(:), rvec_nr3_tmp(:,:), ndegen_nr3_tmp(:)
  CHARACTER(LEN=80) :: buffer
  !
  ioerror = 0
  irpts_nr3 = 0
  !
  filename = TRIM(seedname)//'_hr.dat'
  IF (PRESENT(seedname_)) filename = TRIM(seedname_) // '_hr.dat'
  OPEN(UNIT=10, FILE=filename, ACTION='read', IOSTAT=ioerror)
  IF (ioerror /= 0) CALL io_error('Error opening '//filename)
  !
  READ(10, '(80a)', iostat=ioerror)buffer
  READ(10, '(80a)', iostat=ioerror)buffer
  READ(buffer, *) num_hr_wann
  READ(10, '(80a)', iostat=ioerror)buffer
  READ(buffer, *) nrpts
  IF(is_root) WRITE(*,*) 'begin reading ' // filename
  IF(is_root) WRITE(*,*) 'rvec_z("nr3") = ', nr3
  IF(is_root) WRITE(*,*) 'num_hr_wann = ', num_hr_wann
  IF(is_root) WRITE(*,*) 'nrpts = ', nrpts

  ALLOCATE(ndegen(nrpts))
  ALLOCATE(hr_nr3_tmp(num_hr_wann, num_hr_wann, nrpts))
  ALLOCATE(rvec_nr3_tmp(2, nrpts))
  ALLOCATE(ndegen_nr3_tmp(nrpts))

  DO irpts = 1, (nrpts-1) / 15
    READ(10, *) ndegen((irpts-1)*15+1:(irpts-1)*15+15)
  ENDDO
  READ(10, *) ndegen(((nrpts-1)/15)*15+1:nrpts)
  DO irpts = 1, nrpts
    DO iw = 1, num_hr_wann**2
      READ (10,'(5I5,2F12.6)',iostat=ioerror) ir1, ir2, ir3, ni, nj, hr_re, hr_im
      IF (ioerror/=0) CALL io_error('Error reading file: ' // filename)
      IF (ir3 == nr3) THEN
        IF (iw == 1) THEN
          irpts_nr3 = irpts_nr3 + 1
          rvec_nr3_tmp(:, irpts_nr3) = (/ ir1, ir2 /)
          ndegen_nr3_tmp(irpts_nr3) = ndegen(irpts)
        ENDIF
        hr_nr3_tmp(ni, nj, irpts_nr3) = CMPLX(hr_re, hr_im, DP)
      ENDIF
    ENDDO
  ENDDO
  CLOSE(10)

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
  IF (is_root) write(*,*) 'end reading ' // filename
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
  DEALLOCATE(h11, stat=ierr)
  DEALLOCATE(h12, stat=ierr)
  IF(ALLOCATED(h00)) DEALLOCATE(h00,stat=ierr)
  IF(ALLOCATED(h01)) DEALLOCATE(h01,stat=ierr)
  IF(ALLOCATED(hr0)) DEALLOCATE(hr0,stat=ierr)
  IF(ALLOCATED(hr1)) DEALLOCATE(hr1,stat=ierr)
  IF(ALLOCATED(rvec0)) DEALLOCATE(rvec0,stat=ierr)
  IF(ALLOCATED(rvec1)) DEALLOCATE(rvec1,stat=ierr)
  IF(ALLOCATED(ndegen0)) DEALLOCATE(ndegen0,stat=ierr)
  IF(ALLOCATED(ndegen1)) DEALLOCATE(ndegen1,stat=ierr)
  IF (ierr /= 0) CALL io_error('Error in allocation in hamiltonian_ALLOCATE_slab')
!------------------------------------------------------------------------
END SUBROUTINE hamiltonian_deallocate
!------------------------------------------------------------------------
END MODULE hamiltonian
