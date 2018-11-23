module hamiltonian
    use constants
    use parameters
    IMPLICIT NONE
    COMPLEX(dp), allocatable, dimension(:,:), public, save :: h00, h01, h11, h12
    COMPLEX(dp), public, save :: omega
    REAL(dp), public, save :: kx, ky
    ! INTEGER, allocatable, public, save :: ind_surf(:), ind_bulk(:)

    ! private variables for saving values parsed from seedname_hr.dat
    COMPLEX(dp), allocatable, dimension(:,:,:), save :: hr0s, hr0, hr1
    REAL(dp), allocatable, dimension(:,:), save :: rvec0s, rvec0, rvec1
    REAL(dp), allocatable, dimension(:), save :: ndegen0s, ndegen0, ndegen1
    complex(dp), save :: bulk_shift ! shift of bulk and slab energy (add to bulk onsite)
    INTEGER, save :: nrpts0s, nrpts0, nrpts1

    INTEGER, save :: num_hr_wann ! dimension of hamiltonian in seedname_hr.dat

    ! PUBLIC :: hamiltonian_set_hread, hamiltonian_set_hij
contains
SUBROUTINE hamiltonian_setup()
    ! allocate required variables, and
    ! read seedname_hr.dat for the required nr3, and save
    ! them in hr*, rvec* ndegen*, nrpts*
    IMPLICIT NONE
    INTEGER :: i, ir, ierr
    ! allocation
    IF (isslab) THEN
        ALLOCATE(h00(nsurf,nsurf), stat=ierr)
        ALLOCATE(h01(nsurf,nbulk), stat=ierr)
        ALLOCATE(h11(nbulk,nbulk), stat=ierr)
        ALLOCATE(h12(nbulk,nbulk), stat=ierr)
    ELSE
        ALLOCATE(h11(nbulk,nbulk), stat=ierr)
        ALLOCATE(h12(nbulk,nbulk), stat=ierr)
    END IF
    if (ierr /= 0) CALL io_error('Error in allocation in hamiltonian_setup')
    ! reading seedname_hr.dat
    IF (isslab_hamil) THEN
        CALL read_hamiltonian(0, hr0, rvec0, ndegen0, nrpts0)
    ELSE IF (isslab_match) THEN
        CALL read_hamiltonian(0, hr0, rvec0, ndegen0, nrpts0, trim(seedname)//".bulk")
        CALL read_hamiltonian(1, hr1, rvec1, ndegen1, nrpts1, trim(seedname)//".bulk")
        CALL read_hamiltonian(0, hr0s, rvec0s, ndegen0s, nrpts0s, trim(seedname)//".slab")
        seedname = trim(seedname)//".slab"

        ! calculate the required onsite energy shift
        ! bulk_shift = average(hr0s_onsite - hr0_onsite)
        bulk_shift = 0.0_dp
        do ir = 1, nrpts0s
            if (rvec0s(1,ir)==0 .and. rvec0s(2,ir)==0) then
                do i = 1, nbulk
                    bulk_shift = bulk_shift + hr0s(ind_1(i), ind_1(i), ir) / ndegen0s(ir)
                end do
            end if
        end do
        do ir = 1, nrpts0
            if (rvec0(1,ir)==0 .and. rvec0(2,ir)==0) then
                do i = 1, nbulk
                    bulk_shift = bulk_shift - hr0(i, i, ir) / ndegen0(ir)
                end do
            end if
        end do
        bulk_shift = cmplx(REAL(bulk_shift / nbulk, dp), 0.0_dp, dp)
        if (is_root) write(*,*) "bulk_shift = ", bulk_shift
        if (is_root) write(*,*) "Note: this may(should) be small if util_match is used to generate slab hr.dat file"

    ELSE
        CALL read_hamiltonian(0, hr0, rvec0, ndegen0, nrpts0)
        CALL read_hamiltonian(1, hr1, rvec1, ndegen1, nrpts1)
    END IF
END SUBROUTINE

SUBROUTINE hamiltonian_set_hij()
    ! sum the in-plane hamiltonian multiplied with the phase factor
    ! to create the tight-binding hamiltonian, hij.
    IMPLICIT NONE
    COMPLEX(dp), allocatable, dimension(:,:) :: htemp
    INTEGER :: i
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
        end do
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
        END DO
        DEALLOCATE(htemp)
    ELSE !.NOT. isslab, .NOT. isslab_hamil, .NOT. isslab_match
        CALL k_operator(nrpts0, hr0, rvec0, ndegen0, kx, ky, h11)
        CALL k_operator(nrpts1, hr1, rvec1, ndegen1, kx, ky, h12)
    END IF
END SUBROUTINE

subroutine read_hamiltonian(nr3,hr_nr3,rvec_nr3,ndegen_nr3,irpts_nr3,seedname_)
    ! from seedname_hr.dat, parse hamiltonian elements such that
    ! lattice vector along dir=3 is nr3
    implicit none
    integer, intent(IN) :: nr3
    complex(dp), allocatable, intent(out) :: hr_nr3(:,:,:)
    real(dp), allocatable, intent(out)  :: rvec_nr3(:,:), ndegen_nr3(:)
    integer, intent(out) :: irpts_nr3
    character(*), intent(IN), optional :: seedname_

    character(len=80) :: filename

    integer :: irpts,iw,ni,nj,ioerror,nrpts, ir1,ir2,ir3
    real(dp) :: hr_re, hr_im
    complex(dp), allocatable :: hr_nr3_tmp(:,:,:)
    real(dp), allocatable :: ndegen(:), rvec_nr3_tmp(:,:), ndegen_nr3_tmp(:)
    character(80) :: buffer
    logical :: ir3_checked

    ioerror=0
    irpts_nr3 = 0

    filename = trim(seedname)//'_hr.dat'
    IF (present(seedname_)) filename = trim(trim(seedname_)//'_hr.dat')
    open(unit=10,file=filename,action='read',iostat=ioerror)
    if (ioerror/=0) CALL io_error('Error opening '//filename)

    read(10,'(80a)',iostat=ioerror)buffer
    read(10,'(80a)',iostat=ioerror)buffer
    read(buffer, *) num_hr_wann
    read(10,'(80a)',iostat=ioerror)buffer
    read(buffer, *) nrpts
    IF(is_root) write(*,*) 'begin reading '//filename
    IF(is_root) write(*,*) 'rvec_z("nr3") = ', nr3
    IF(is_root) write(*,*) 'num_hr_wann = ', num_hr_wann
    IF(is_root) write(*,*) 'nrpts = ', nrpts

    allocate(ndegen(nrpts))
    allocate(hr_nr3_tmp(num_hr_wann,num_hr_wann,nrpts))
    allocate(rvec_nr3_tmp(2,nrpts))
    allocate(ndegen_nr3_tmp(nrpts))
    
    do irpts = 1, (nrpts-1)/15
        read(10, *) ndegen((irpts-1)*15+1:(irpts-1)*15+15)
    end do
    read(10, *) ndegen(((nrpts-1)/15)*15+1:nrpts)
    do irpts = 1, nrpts
        ir3_checked = .false.
        do iw=1,num_hr_wann**2
            read (10,'(5I5,2F12.6)',iostat=ioerror) ir1,ir2,ir3, ni,nj, hr_re, hr_im
            if (ioerror/=0) CALL io_error('Error reading file: '//filename)
            if (ir3 == nr3) then
                if (.not. ir3_checked) then
                    irpts_nr3 = irpts_nr3 + 1
                    ir3_checked = .true.
                end if
                hr_nr3_tmp(ni,nj,irpts_nr3) = cmplx(hr_re, hr_im, dp)
                rvec_nr3_tmp(:,irpts_nr3) = (/ ir1, ir2 /)
                ndegen_nr3_tmp(irpts_nr3) = ndegen(irpts)
            end if
        enddo
    enddo
    CLOSE(10)

    ! change size of array to known(read) size
    allocate(hr_nr3(num_hr_wann,num_hr_wann,irpts_nr3))
    allocate(rvec_nr3(2,irpts_nr3))
    allocate(ndegen_nr3(irpts_nr3))
    hr_nr3 = hr_nr3_tmp(:,:,1:irpts_nr3)
    rvec_nr3 = rvec_nr3_tmp(:,1:irpts_nr3)
    ndegen_nr3 = ndegen_nr3_tmp(1:irpts_nr3)
    deallocate(hr_nr3_tmp)
    deallocate(rvec_nr3_tmp)
    deallocate(ndegen_nr3_tmp)
    deallocate(ndegen)

    IF (is_root) write(*,'("Hamiltonian_",I1," shape : ",3I5)') nr3, shape(hr_nr3)
    IF(is_root) write(*,*) 'end reading '//filename
end subroutine


SUBROUTINE hamiltonian_deallocate()
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
    IF (ierr /= 0) CALL io_error('Error in allocation in hamiltonian_allocate_slab')
END SUBROUTINE
end module
