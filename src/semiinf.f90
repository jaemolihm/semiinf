!------------------------------------------------------------------------
PROGRAM semiinf
!------------------------------------------------------------------------
!! Main driver of the semiinf.x program
!------------------------------------------------------------------------
  USE comms, ONLY : DP, is_root, mp_setup, mp_barrier, mp_end, world_comm, &
    io_error, ci, cone, find_free_unit
  USE parameters, ONLY : nbulk, isspin, is_ideal_surf, sigma, num_energy, energy, &
    input_filename, param_read, param_write
  USE iter_bulk, ONLY : iter_bulk_allocate
  USE iter_slab, ONLY : iter_slab_allocate
  USE hamiltonian, ONLY : hamiltonian_setup
  USE postprocess_green, ONLY : pp_green_setup
  !
  IMPLICIT NONE
  !
  LOGICAL :: flag_converged
  !! Output of iterations for Green function. True if convergence is reached.
  INTEGER :: n_dos_layer = 0 ! FIXME: this should go to parameters.f90 as input parameter
  INTEGER :: ik, irec, ik_start, ik_end
  INTEGER :: iunsurf, iunbulk, iunsurfsx, iunsurfsy, iunsurfsz, iunlayer
  REAL(DP), ALLOCATABLE :: dos_s0(:), dos_b(:), dos_spn_s(:,:)
  ! REAL(DP), ALLOCATABLE :: dos_up(:), dos_dn(:), dos_spn_s_up(:,:), dos_spn_s_dn(:,:)
  REAL(DP), ALLOCATABLE :: dos_layer(:,:)
  !
  CALL mp_setup()
  IF (is_root) THEN
    WRITE(*, '(1x,a)') "======================================================"
    WRITE(*, '(1x,a)') '                 Welcome to semiinf.x !'
  ENDIF
  !
  ! Read input file
  CALL GETARG(1, input_filename)
  input_filename = TRIM(input_filename)
  CALL param_read()
  IF (is_root) CALL param_write()
  !
  ! Allocate array in modules
  IF (is_ideal_surf) THEN
    CALL iter_bulk_allocate()
  ELSE
    CALL iter_slab_allocate()
  ENDIF
  CALL hamiltonian_setup()
  CALL pp_green_setup()
  !
  ALLOCATE(dos_s0(num_energy))
  ALLOCATE(dos_b(num_energy))
  IF (n_dos_layer > 0) ALLOCATE(dos_layer(num_energy, n_dos_layer))
  IF (isspin) ALLOCATE(dos_spn_s(num_energy, 3))
  !
  ! Open units for binary file output
  INQUIRE(IOLENGTH=irec) dos_s0
  iunsurf = find_free_unit()
  OPEN(UNIT=iunsurf, FILE="dos_surf00.out", FORM="unformatted", ACCESS="direct", RECL=irec)
  !
  INQUIRE(IOLENGTH=irec) dos_b
  iunbulk = find_free_unit()
  OPEN(UNIT=iunbulk, FILE="dos_bulk.out", FORM="unformatted", ACCESS="direct", RECL=irec)
  IF (isspin) THEN
    iunsurfsx = find_free_unit()
    OPEN(UNIT=iunsurfsx, FILE="dos_surf_sx.out", FORM="unformatted", ACCESS="direct", RECL=irec)
    iunsurfsy = find_free_unit()
    OPEN(UNIT=iunsurfsy, FILE="dos_surf_sy.out", FORM="unformatted", ACCESS="direct", RECL=irec)
    iunsurfsz = find_free_unit()
    OPEN(UNIT=iunsurfsz, FILE="dos_surf_sz.out", FORM="unformatted", ACCESS="direct", RECL=irec)
  ENDIF
  IF (n_dos_layer > 0) THEN
    INQUIRE(IOLENGTH=irec) dos_layer(:,:)
    iunlayer = find_free_unit()
    OPEN(UNIT=iunlayer, FILE="dos_layer.out", FORM="unformatted", ACCESS="direct", RECL=irec)
  ENDIF
  !
  CALL distribute_kpoint()
  !
  ! Loop over k points
  DO ik = ik_start, ik_end
    CALL run_kpoint(ik)
  END DO
  CLOSE(iunsurf)
  CLOSE(iunbulk)
  IF (isspin) THEN
    CLOSE(iunsurfsx)
    CLOSE(iunsurfsy)
    CLOSE(iunsurfsz)
  ENDIF
  IF (n_dos_layer > 0) THEN
    CLOSE(iunlayer)
  ENDIF
  !
  CALL mp_barrier(world_comm)
  IF (is_root) THEN
    WRITE(*, '(1x,a)') 'Calculation complete: semiinf exiting'
    WRITE(*, '(1x,a)') "======================================================"
  ENDIF
  !
  CALL mp_end()
  !
CONTAINS
!------------------------------------------------------------------------
SUBROUTINE distribute_kpoint()
!------------------------------------------------------------------------
!! distribute k points to processors for k-point parallelization
!------------------------------------------------------------------------
  USE comms, ONLY : num_procs, my_id, is_root
  USE parameters, ONLY : num_kpoint
  IMPLICIT NONE
  !
  INTEGER :: q, r
  ! nk = q * (np-r) + (q+1) * r
  q = num_kpoint / num_procs
  r = num_kpoint - (num_procs * q)
  IF (is_root) THEN
    WRITE(*, "(X,'There are', I8,' or ',I8,' k-points per processor')") q, q+1
  ENDIF
  !
  IF (my_id < (num_procs - r)) THEN
    ik_end = q * (my_id + 1)
    ik_start = ik_end - q + 1
  ELSE
    ik_end = q * (my_id + 1) + (my_id - num_procs + r + 1)
    ik_start = ik_end - q
  ENDIF
  ! For debugging
  ! PRINT*, "Number of k points in node", my_id, " is", ik_end - ik_start + 1
!------------------------------------------------------------------------
END SUBROUTINE distribute_kpoint
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE run_kpoint(ik)
!------------------------------------------------------------------------
!! Main driver of the calculation for each k point.
!! For given kx and ky, calculate Green function and write DOS to file
!------------------------------------------------------------------------
  USE parameters, ONLY : kpoints, omega
  USE hamiltonian, ONLY : hamiltonian_tb_to_k
  USE iter_bulk, ONLY : iter_bulk_main
  USE iter_slab, ONLY : iter_slab_main
  USE postprocess_green, ONLY : get_dos_s, get_dos_b, get_dos_nlayer, &
    set_transfer_mat, get_spin_s, pp_green_deallocate_green_layer
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik
  !! Index of the k point to be calculated
  REAL(DP) :: kx
  !! k vector along x axis, in crystal coordinate
  REAL(DP) :: ky
  !! k vector along y axis, in crystal coordinate
  INTEGER :: ie, il
  !
  ! Set k-dependent hamiltonian
  kx = kpoints(1, ik)
  ky = kpoints(2, ik)
  WRITE(*, '("ik = ", I4, " kx = ", F6.3, " ky = ", F6.3)') ik, kx, ky
  CALL hamiltonian_tb_to_k(kx, ky)
  !
  ! Loop over energy points
  DO ie = 1, num_energy
    omega = energy(ie) * cone - sigma * ci
    ! Compute Green function by iteration
    IF (is_ideal_surf) THEN
      CALL iter_bulk_main(flag_converged)
    ELSE
      CALL iter_slab_main(flag_converged)
    ENDIF
    IF (.NOT. flag_converged) &
      WRITE(*, '("WARNING: DOS NOT converged, ik=",I8," ie=",I8)') ik, ie
    !
    ! Postprocess Green functions to calculate DOS, spin-DOS, etc.
    CALL set_transfer_mat()
    CALL get_dos_s(dos_s0(ie))
    CALL get_dos_b(dos_b(ie))
    IF (n_dos_layer > 0) CALL get_dos_nlayer(n_dos_layer, dos_layer(ie,:))
    IF (isspin) CALL get_spin_s(dos_spn_s(ie,:))
    !
    CALL pp_green_deallocate_green_layer()
  END DO ! ie
  !
  ! Write DOS to file
  WRITE(UNIT=iunsurf, REC=ik) dos_s0(:)
  WRITE(UNIT=iunbulk, REC=ik) dos_b(:)
  IF (isspin) WRITE(UNIT=iunsurfsx, REC=ik) dos_spn_s(:,1)
  IF (isspin) WRITE(UNIT=iunsurfsy, REC=ik) dos_spn_s(:,2)
  IF (isspin) WRITE(UNIT=iunsurfsz, REC=ik) dos_spn_s(:,3)
  DO il = 1, n_dos_layer
    ! WRITE(unit=unit_il(il), REC=ik) dos_layer(:,il)
    WRITE(UNIT=11, REC=ik) dos_layer(:,:)
  END DO
  !
  WRITE(*,'("ik=",I6.5," done")') ik
!------------------------------------------------------------------------
END SUBROUTINE run_kpoint
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
END PROGRAM semiinf
!------------------------------------------------------------------------