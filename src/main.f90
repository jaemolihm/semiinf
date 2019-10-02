program main
    USE comms
    USE parameters
    USE hamiltonian, ONLY : omega, kx, ky
    USE iter_bulk
    USE iter_slab
    USE postprocess_green
    IMPLICIT NONE

    character(len=20) :: buffer_temp
    integer :: ik, irec, ik_start, ik_end, il
    real(DP), allocatable :: dos_s0(:), dos_b(:), dos_spn_s(:,:)
    real(DP), allocatable :: dos_up(:), dos_dn(:), dos_spn_s_up(:,:), dos_spn_s_dn(:,:)
    complex(DP), allocatable :: dos_layer(:,:)
    integer :: n_dos_layer ! this should later go to parameters.f90
    n_dos_layer = 0

    CALL mpi_setup()
    
    ! ... read input file
    CALL GETARG(1, input_filename)
    input_filename = TRIM(input_filename)
    CALL param_read()
    IF (is_root) CALL param_write()

    ! ... allocate array in modules
    IF (isslab) CALL iter_slab_allocate()
    ! IF ((.NOT. isslab) .AND. (.NOT. isfinite)) CALL iter_bulk_allocate()
    IF ((.NOT. isslab)) CALL iter_bulk_allocate()
    CALL hamiltonian_setup()
    CALL pp_green_setup()
    
    ! IF (isfinite .and. isspin) THEN
        ! ALLOCATE(dos_up(num_energy))
        ! ALLOCATE(dos_dn(num_energy))
        ! ALLOCATE(dos_s0(num_energy))
        ! IF (isspin) ALLOCATE(dos_spn_s(num_energy,3))
        ! IF (isspin) ALLOCATE(dos_spn_s_up(num_energy,3))
        ! IF (isspin) ALLOCATE(dos_spn_s_dn(num_energy,3))
        ! inquire (iolength = irec) dos_s0
        ! open (unit = 8, file = "dos_surfup.out", form="unformatted", access="direct", recl=irec)
        ! open (unit = 9, file = "dos_surfdn.out", form="unformatted", access="direct", recl=irec)
        ! open (unit = 10, file = "dos_slab.out", form="unformatted", access="direct", recl=irec)
        ! open (unit = 3, file = "dos_surf_sslab.out", form="unformatted", access="direct", recl=irec*3)
        ! open (unit = 4, file = "dos_surf_sup.out", form="unformatted", access="direct", recl=irec*3)
        ! open (unit = 5, file = "dos_surf_sdn.out", form="unformatted", access="direct", recl=irec*3)
    ! ELSE
        ALLOCATE(dos_s0(num_energy))
        ALLOCATE(dos_b(num_energy))
        IF (n_dos_layer>0) ALLOCATE(dos_layer(num_energy, 64))
        IF (isspin) ALLOCATE(dos_spn_s(num_energy,3))
        inquire (iolength = irec) dos_s0
        open (unit = 9, file = "dos_surf00.out", form="unformatted", access="direct", recl=irec)
        inquire (iolength = irec) dos_b
        open (unit = 10, file = "dos_bulk.out", form="unformatted", access="direct", recl=irec)
        IF (isspin) open (unit = 3, file = "dos_surf_sx.out", form="unformatted", access="direct", recl=irec)
        IF (isspin) open (unit = 4, file = "dos_surf_sy.out", form="unformatted", access="direct", recl=irec)
        IF (isspin) open (unit = 5, file = "dos_surf_sz.out", form="unformatted", access="direct", recl=irec)
        IF (n_dos_layer > 0) THEN
            inquire (iolength = irec) dos_layer(:,:)
            write(buffer_temp,'("dos_surf",I0.2,".out")'), 1
            open (unit=11, file = buffer_temp, form="unformatted", access="direct", recl=irec)
        END IF
        ! DO il = 1, n_dos_layer
            ! write(buffer_temp,'("dos_surf",I0.2,".out")'), il
            ! open (unit=unit_il(il), file = buffer_temp, form="unformatted", access="direct", recl=irec)
        ! END DO
    ! END IF

    CALL distribute_kpoint()    
    ! ... loop over k points
    DO ik = ik_start, ik_end
        CALL run_kpoint()
    END DO
    CLOSE(9)
    CLOSE(10)
    IF (isspin) CLOSE(3)
    IF (isspin) CLOSE(4)
    IF (isspin) CLOSE(5)
    CLOSE(11)
    ! DO il = 1, n_dos_layer
        ! CLOSE(unit_il(il))
    ! END DO

    CALL mpi_end()
    
CONTAINS
SUBROUTINE distribute_kpoint()
    integer :: ierr, an_id, nk_per_proc_quo, nk_per_proc_rem
    
    nk_per_proc_quo = num_kpoint / num_procs
    nk_per_proc_rem = num_kpoint - (num_procs*nk_per_proc_quo)
    ! nk = q * (np-r) + (q+1) * r
    write(*,"(1X,I,' or ',I,' k-points per process')") nk_per_proc_quo, nk_per_proc_quo+1
    IF (my_id < (num_procs - nk_per_proc_rem)) THEN
        ik_end = nk_per_proc_quo * (my_id + 1)
        ik_start = ik_end - nk_per_proc_quo + 1
    ELSE
        ik_end = nk_per_proc_quo * (my_id + 1) + (my_id - num_procs + nk_per_proc_rem + 1)
        ik_start = ik_end - nk_per_proc_quo
    ENDIF
    print*, "number of k point in node ", my_id, " is ", ik_end-ik_start+1
END SUBROUTINE

SUBROUTINE run_kpoint()
    IMPLICIT NONE
    INTEGER :: ie, il
    ! ... set k-dependent hamiltonian
    kx = plot_kpoint(1,ik)
    ky = plot_kpoint(2,ik)
    write(*,'("ik = ", I4, " kx = ", F6.3, " ky = ", F6.3)') ik, kx, ky
    CALL hamiltonian_tb_to_k()
    ! ... loop over energy points
    DO ie = 1, num_energy
        omega = energy(ie) * cone - sigma * ci
        ! run iteration
        IF (isslab) THEN
            CALL iter_slab_main()
        ! ELSE IF (isfinite) THEN
            ! CALL iter_slab_green_finite()
            ! flag_converged = .true.
        ELSE
            CALL iter_bulk_main()
        END IF
        IF (.NOT. flag_converged) WRITE (*,'("DOS NOT converged, s=",ES8.1," N=",I2," ik=",I4," ie=",I5)') sigma, n_iter, ik, ie
!        IF (MOD(ie,100)==0) WRITE (*,'("DOS converged, s=",ES8.1," N=",I2," ik=",I4," ie=",I5)') sigma, n_iter, ik, ie

        ! run pp_green
        ! IF (isfinite) THEN
            ! CALL get_dos_s(dos_s0(ie))
            ! CALL get_dos_up(dos_up(ie))
            ! CALL get_dos_dn(dos_dn(ie))
            ! IF (isspin) CALL get_spin_s(dos_spn_s(ie,:))
            ! IF (isspin) CALL get_spin_up(dos_spn_s_up(ie,:))
            ! IF (isspin) CALL get_spin_dn(dos_spn_s_dn(ie,:))
        ! ELSE
            CALL set_transfer_mat()
            CALL get_dos_s(dos_s0(ie))
            CALL get_dos_b(dos_b(ie))
            IF (n_dos_layer > 0) CALL get_dos_nlayer(n_dos_layer, dos_layer(ie,:))
            IF (isspin) CALL get_spin_s(dos_spn_s(ie,:))
        ! END IF
        CALL pp_green_deallocate_green_layer()
    END DO
    ! IF (isfinite) THEN
        ! WRITE(unit=10, rec=ik) dos_s0
        ! WRITE(unit=8, rec=ik) dos_up
        ! WRITE(unit=9, rec=ik) dos_dn
        ! IF (isspin) WRITE(unit=3, rec=ik) dos_spn_s
        ! IF (isspin) WRITE(unit=4, rec=ik) dos_spn_s_up
        ! IF (isspin) WRITE(unit=5, rec=ik) dos_spn_s_dn
    ! ELSE
        WRITE(unit=9, rec=ik) dos_s0(:)
        WRITE(unit=10, rec=ik) dos_b(:)
        IF (isspin) WRITE(unit=3, rec=ik) dos_spn_s(:,1)
        IF (isspin) WRITE(unit=4, rec=ik) dos_spn_s(:,2)
        IF (isspin) WRITE(unit=5, rec=ik) dos_spn_s(:,3)
        DO il = 1, n_dos_layer
            ! WRITE(unit=unit_il(il), rec=ik) dos_layer(:,il)
            WRITE(unit=11, rec=ik) dos_layer(:,:)
        END DO
    ! END IF
    WRITE(*,'("ik=",I6.5," done")') ik
END SUBROUTINE

! SUBROUTINE write_dos(unit, dos_list)
!     ! USE kinds
!     IMPLICIT NONE
!     integer, INTENT(IN) :: unit, ik
!     real(DP), dimension(:), INTENT(IN) :: dos_list
!     WRITE(unit, rec=ik) dos_list
! END SUBROUTINE
end program main

