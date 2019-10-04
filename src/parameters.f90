!------------------------------------------------------------------------
!! Adapted by Jae-Mo Lihm from wannier90/src/parameters.F90
!------------------------------------------------------------------------
MODULE parameters
!------------------------------------------------------------------------
!! This module contains parameters to control the actions of semiinf.x.
!! Also have routines to read the parameters and write them out again.
!------------------------------------------------------------------------
  USE comms, ONLY : DP, io_error
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER, PARAMETER :: maxlen = 256
  !
  ! --------------------- input variables ---------------------
  CHARACTER(LEN=256) :: seedname
  !! Prefix of the .dat files.
  LOGICAL :: is_ideal_surf
  !! Is the calculation with surface modification. If .TRUE., all principal
  !! layers, including the surface principal layer, are identical.
  LOGICAL :: hr_stitching
  !! Is the calculation with surface (slab) and bulk matching.
  !! If .true., read seedname.bulk_hr.dat and seedname.slab_hr.dat file
  !! If .false., read seedname_hr.dat file only, and set surface and bulk
  !! regions using ind_0, ind_1, ind_2.
  INTEGER :: nbulk
  !! Number of basis functions for the bulk principal layer
  INTEGER :: nsurf
  !! Number of basis functions for the surface principal layer
  !! If is_ideal_surf, automatically set to nbulk, but not used.
  INTEGER, ALLOCATABLE :: ind_0(:)
  !! (nsurf) Indices of basis functions for surface principal layer.
  !! Used only if is_ideal_surf is false.
  INTEGER, ALLOCATABLE :: ind_1(:)
  !! (nbulk) Indices of basis functions for bulk principal layer.
  !! Used only if is_ideal_surf is false.
  INTEGER, ALLOCATABLE :: ind_2(:)
  !! (nbulk) Indices of basis functions for second bulk principal layer.
  !! Used to set the bulk-bulk interlayer hopping.
  !! Used only if is_ideal_surf is false and hr_stitching is false.
  INTEGER :: bulk_rz
  !! 1 or -1. Determines the surface normal direction. Default is 1.
  !! From bulk _hr.dat, read matrix elements <m0|H|nR> with R_3 == bulk_rz
  !! Used only if is_ideal_surf is false or hr_stitching is true.
  !! If is_ideal_surf is false and hr_stitching is false, the surface normal direction
  !! is determined by ind_0,1,2.
  REAL(DP) :: hopping_tol
  !! Convergence criterion of the iterative surface Green function calculation
  INTEGER :: max_n_iter
  !! Maximum number of iterations for surface Green function
  !
  ! Variables related to energy sampling
  REAL(DP) :: dos_e_min
  !! Minimum energy to calculate DOS. (in eV)
  REAL(DP) :: dos_e_max
  !! Maximum energy to calculate DOS. (in eV)
  REAL(DP) :: dos_e_step
  !! Step size of energy sampling. (in eV)
  REAL(DP) :: sigma
  !! Small imaginary number for green function
  LOGICAL :: isspin
  !! If true, read _spnr.dat file and calculate spin-DOS
  !
  ! Variables related to kpoint
  CHARACTER(4) :: kpoint_type
  !! The way k point is given. 'file' or 'grid' only.
  !! file: the kpoints are read from file.
  !! grid: the kpoints are chosen as a grid in the 2d BZ
  CHARACTER(LEN=256) :: kpoint_file
  !! Filename where the k points are written.
  !! Used only if kpoint_type is file.
  !! The file content must be like
  !! num_kpoint
  !! kpoints(1, 1)  kpoints(2, 1)
  !! ...
  !! kpoints(1, num_kpoint)  kpoints(2, num_kpoint)
  INTEGER :: nkgrid_1
  !! Number of k point grids along crystal lattice vector 1.
  !! Used only if kpoint_type is grid
  INTEGER :: nkgrid_2
  !! Number of k point grids along crystal lattice vector 2.
  !! Used only if kpoint_type is grid
  REAL(DP) :: k1_min
  !! Minimum value of k along crystal lattice vector 1.
  !! Used only if kpoint_type is grid.
  !! Default is -0.5
  REAL(DP) :: k1_max
  !! Maximum value of k along crystal lattice vector 1.
  !! Used only if kpoint_type is grid.
  !! Default is +0.5
  REAL(DP) :: k2_min
  !! Minimum value of k along crystal lattice vector 2.
  !! Used only if kpoint_type is grid.
  !! Default is -0.5
  REAL(DP) :: k2_max
  !! Maximum value of k along crystal lattice vector 2.
  !! Used only if kpoint_type is grid.
  !! Default is +0.5
  !
  ! --------------------- derived variables ---------------------
  ! These variables are not directly read, but derived from the input variables
  INTEGER :: num_energy
  !! Number of energy values to calculate DOS
  REAL(DP), ALLOCATABLE :: energy(:)
  !! (num_energy) List of energy values to calculate DOS (in eV)
  INTEGER :: num_kpoint
  !! number of sampled k-points, in crystal coordinates
  REAL(DP), ALLOCATABLE :: kpoints(:,:)
  !! (3, num_kpoint) k points to calculate the spectral function
  !
  ! --------------------- other variables ---------------------
  ! These variables are not set during the input step.
  CHARACTER(LEN=256) :: input_filename
  !! name of the input file. Read as inline argument.
  COMPLEX(DP), PUBLIC :: omega
  !! Frequency to calculate Green function. omega = energy - i * sigma
  COMPLEX(DP), ALLOCATABLE :: green_s(:,:)
  !! Green function for the surface principal layer
  COMPLEX(DP), ALLOCATABLE :: green_s1(:,:)
  !! Green function for the first sub-surface principal layer
  COMPLEX(DP), ALLOCATABLE :: green_b(:,:)
  !! Green function for a bulk principal layer
  !
  ! --------------------- private variables ---------------------
  ! For reading kpoint path
  INTEGER, PRIVATE :: num_paths
  CHARACTER(LEN=1), ALLOCATABLE, PRIVATE :: bands_label(:)
  INTEGER, PRIVATE :: bands_num_points
  INTEGER, PRIVATE :: bands_num_spec_points
  REAL(DP), ALLOCATABLE, PRIVATE :: bands_spec_points(:,:)
  ! Others
  INTEGER, PRIVATE :: num_lines
  INTEGER, PRIVATE :: itemp1, itemp2
  CHARACTER(LEN=256), ALLOCATABLE, PRIVATE :: in_data(:)
  CHARACTER(LEN=256), PRIVATE :: buffer
CONTAINS
!------------------------------------------------------------------------
SUBROUTINE param_write
  WRITE(*,'("is_ideal_surf = ", L)') is_ideal_surf
  WRITE(*,'("hr_stitching = ", L)') hr_stitching
  WRITE(*,'("nsurf = ", I)') nsurf
  WRITE(*,'("nbulk = ", I)') nbulk
  WRITE(*,'("bulk_rz = ", I)') bulk_rz
  WRITE(*,'("num_energy = ",I)') num_energy
  WRITE(*,'("num_kpoint = ",I)') num_kpoint
END SUBROUTINE
!------------------------------------------------------------------------
SUBROUTINE param_read
!------------------------------------------------------------------------
!! Read input parameters and calculate derived values
!------------------------------------------------------------------------
  USE comms, ONLY : pi, find_free_unit
  IMPLICIT NONE
  LOGICAL :: found
  INTEGER :: ierr, i_temp, ik1, ik2, ik, nlen_temp, iunk
  !
  CALL param_in_file
  !
  CALL param_get_keyword('seedname',found,c_value=seedname)
  !
  ! Set default values
  hopping_tol = 1.D-8
  max_n_iter = 40
  isspin = .FALSE.
  bulk_rz = 1
  ! read parameters
  CALL param_get_keyword('hopping_tol', found, r_value=hopping_tol)
  CALL param_get_keyword('max_n_iter', found, i_value=max_n_iter)
  CALL param_get_keyword('isspin', found, l_value=isspin)
  !
  CALL param_get_keyword('is_ideal_surf', found, l_value=is_ideal_surf)
  IF (.NOT. found) CALL io_error('is_ideal_surf must be set')
  !
  CALL param_get_keyword('hr_stitching', found, l_value=hr_stitching)
  if (.NOT. is_ideal_surf .AND. .NOT. found) &
    CALL io_error('If is_ideal_surf, hr_stitching must be set.')
  CALL param_get_keyword('nbulk', found, i_value=nbulk)
  IF (.NOT. found) CALL io_error('nbulk must be set')
  !
  CALL param_get_keyword('nsurf', found, i_value=nsurf)
  IF (.NOT. is_ideal_surf .AND. .NOT. found) &
    CALL io_error('If not is_ideal_surf, nsurf must be set')
  !
  CALL param_get_keyword('bulk_rz', found, i_value=bulk_rz)
  !
  ! Check validity of input parameters
  if (.NOT. is_ideal_surf .AND. (.NOT. hr_stitching) .AND. (nsurf <= nbulk)) &
    CALL io_error('ERROR: if not is_ideal_surf and not hr_stitching, nsurf must be greater than nbulk')
    ! FIXME : why?
  !
  IF (is_ideal_surf) THEN
    nsurf = nbulk
  ENDIF
  !
  IF (.NOT. is_ideal_surf) THEN
    CALL param_get_range_vector('ind_0', found, nlen_temp, .TRUE.)
    IF (nlen_temp /= nsurf) &
      CALL io_error('ERROR: length of ind_0 must be equal to nsurf')
    ALLOCATE(ind_0(nsurf))
    CALL param_get_range_vector('ind_0', found, nsurf,.FALSE., ind_0)
    !
    CALL param_get_range_vector('ind_1', found, nlen_temp, .TRUE.)
    IF (nlen_temp /= nbulk) &
      CALL io_error('ERROR: length of ind_1 must be equal to nbulk')
    ALLOCATE(ind_1(nbulk))
    CALL param_get_range_vector('ind_1', found, nbulk, .FALSE., ind_1)
  ENDIF
  IF (.NOT. is_ideal_surf .AND. .NOT. hr_stitching) THEN
    CALL param_get_range_vector('ind_2', found, nlen_temp, .TRUE.)
    IF (nlen_temp /= nbulk) &
      CALL io_error('ERROR: length of ind_2 must be equal to nbulk')
    ALLOCATE(ind_2(nbulk))
    CALL param_get_range_vector('ind_2', found, nbulk, .FALSE., ind_2)
  ENDIF
  !
  ! energy
  CALL param_get_keyword('sigma', found, r_value=sigma)
  CALL param_get_keyword('dos_e_min', found, r_value=dos_e_min)
  IF (.NOT. found) CALL io_error("dos_e_min must be set")
  CALL param_get_keyword('dos_e_max', found, r_value=dos_e_max)
  IF (.NOT. found) CALL io_error("dos_e_max must be set")
  CALL param_get_keyword('dos_e_step', found, r_value=dos_e_step)
  IF (.NOT. found) CALL io_error("dos_e_step must be set")
  num_energy = NINT((dos_e_max - dos_e_min) / dos_e_step)
  dos_e_step = (dos_e_max - dos_e_min) / (num_energy-1)
  ALLOCATE(energy(num_energy))
  DO i_temp = 1, num_energy
    energy(i_temp) = dos_e_min + dos_e_step * (i_temp-1)
  ENDDO
  !
  ! kpoints
  CALL param_get_keyword('kpoint_type', found, c_value=kpoint_type)
  IF (.NOT. found) CALL io_error('kpoint_type must be set')
  !
  ! set kpoints path
  IF (kpoint_type == 'file') THEN
    CALL param_get_keyword('kpoint_file', found, c_value=kpoint_file)
    IF (.NOT. found) CALL io_error('If kpoint_type is file, kpoint_file must be set')
    !
    iunk = find_free_unit()
    OPEN(UNIT=iunk, FILE=TRIM(kpoint_file), ACTION='read', IOSTAT=ierr)
    IF (ierr /= 0) CALL io_error('Error opening ' // TRIM(kpoint_file))
    !
    READ(iunk, *, IOSTAT=ierr) num_kpoint
    IF (ierr /= 0) CALL io_error('Error reading num_kpoint from ' // TRIM(kpoint_file))
    !
    ALLOCATE(kpoints(2, num_kpoint))
    DO ik = 1, num_kpoint
      READ(iunk, *, IOSTAT=ierr) kpoints(1,ik), kpoints(2,ik)
      IF (ierr /= 0) CALL io_error('Error reading num_kpoint from ' // TRIM(kpoint_file))
    ENDDO
    !
    CLOSE(iunk)
  ELSE IF (kpoint_type == 'grid') THEN
    k1_min = -0.5d0
    k1_max = 0.5d0
    k2_min = -0.5d0
    k2_max = 0.5d0
    CALL param_get_keyword('nkgrid_1', found, i_value=nkgrid_1)
    IF (.NOT. found) CALL io_error('If kpoint_type is grid, nkgrid_1 must be set')
    CALL param_get_keyword('nkgrid_2', found, i_value=nkgrid_2)
    IF (.NOT. found) CALL io_error('If kpoint_type is grid, nkgrid_2 must be set')
    CALL param_get_keyword('k1_min', found, r_value=k1_min)
    CALL param_get_keyword('k1_max', found, r_value=k1_max)
    CALL param_get_keyword('k2_min', found, r_value=k2_min)
    CALL param_get_keyword('k2_max', found, r_value=k2_max)
    !
    num_kpoint = nkgrid_1 * nkgrid_2
    ALLOCATE(kpoints(2, num_kpoint), STAT=ierr)
    ik = 0
    DO ik1 = 1, nkgrid_1
      DO ik2 = 1, nkgrid_2
        ik = ik + 1
        kpoints(1,ik) = &
            REAL(ik1-1, DP) / REAL(nkgrid_1-1, DP) * (k1_max - k1_min) + k1_min
        kpoints(2,ik) = &
            REAL(ik2-1, DP) / REAL(nkgrid_2-1, DP) * (k2_max - k2_min) + k2_min
      ENDDO
    ENDDO
  ELSE
    CALL io_error('kpoint_type not understood. It must be file or grid.')
  ENDIF
!------------------------------------------------------------------------
END SUBROUTINE param_read
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE param_get_keyword(keyword,found,c_value,l_value,i_value,r_value)
!------------------------------------------------------------------------
!! Finds the value of the required keyword.
!------------------------------------------------------------------------
  IMPLICIT NONE

  character(*),      intent(in)  :: keyword
  !! Keyword to examine
  logical          , intent(out) :: found
  !! Is keyword present
  character(*)     ,optional, intent(inout) :: c_value
  !! Keyword value
  logical          ,optional, intent(inout) :: l_value
  !! Keyword value
  INTEGER          ,optional, intent(inout) :: i_value
  !! Keyword value
  real(kind=dp)    ,optional, intent(inout) :: r_value
  !! Keyword value

  INTEGER           :: kl, in,loop,itmp
  character(len=maxlen) :: dummy

  kl=len_trim(keyword)

  found=.false.

  do loop=1,num_lines
     in=index(in_data(loop),trim(keyword))
     if (in==0 .or. in>1 ) cycle
     itmp=in+len(trim(keyword))
     if (in_data(loop)(itmp:itmp)/='=' &
          .and. in_data(loop)(itmp:itmp)/=':' &
          .and. in_data(loop)(itmp:itmp)/=' ') cycle
     if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
     endif
     found=.true.
     dummy=in_data(loop)(kl+1:)
     in_data(loop)(1:maxlen) = ' '
     dummy=adjustl(dummy)
     if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
        dummy=dummy(2:)
        dummy=adjustl(dummy)
     end if
  end do

  if(found) then
     if( present(c_value) ) c_value=dummy
     if( present(l_value) ) then
        if (index(dummy,'t') > 0) then
           l_value=.true.
        elseif (index(dummy,'f') > 0) then
           l_value=.false.
        else
           call io_error('Error: Problem reading logical keyword '//trim(keyword))
        endif
     endif
     if( present(i_value) ) read(dummy,*,err=220,end=220) i_value
     if( present(r_value) ) read(dummy,*,err=220,end=220) r_value
  end if

  RETURN

220 CALL io_error('Error: Problem reading keyword '//trim(keyword))
!------------------------------------------------------------------------
END SUBROUTINE param_get_keyword
!------------------------------------------------------------------------
SUBROUTINE param_in_file
!------------------------------------------------------------------------
!! Load the *.win file into a character array in_file, ignoring comments and
!! blank lines and converting everything to lowercase characters
!------------------------------------------------------------------------
  USE comms, ONLY : find_free_unit
  IMPLICIT NONE

  INTEGER           :: in_unit,tot_num_lines,ierr,line_counter,loop,in1,in2
  character(len=maxlen) :: dummy
  INTEGER           :: pos
  character, parameter :: TABCHAR = char(9)
  !
  in_unit = find_free_unit()
  open (in_unit, file=trim(input_filename),form='formatted',status='old',err=101)

  num_lines=0;tot_num_lines=0
  do
     read(in_unit, '(a)', iostat = ierr, err= 200, end =210 ) dummy
     ! [GP-begin, Apr13, 2012]: I convert all tabulation characters to spaces
     pos = index(dummy,TABCHAR)
     do while (pos .ne. 0)
        dummy(pos:pos) = ' '
        pos = index(dummy,TABCHAR)
     end do
     ! [GP-end]
     dummy=adjustl(dummy)
     tot_num_lines=tot_num_lines+1
     if( .not.dummy(1:1)=='!'  .and. .not. dummy(1:1)=='#' ) then
        if(len(trim(dummy)) > 0 ) num_lines=num_lines+1
     endif

  end do

101 call io_error('Error: Problem opening input file '//trim(input_filename))
200 call io_error('Error: Problem reading input file '//trim(input_filename))
210 continue
  rewind(in_unit)

  allocate(in_data(num_lines),stat=ierr)
  if (ierr/=0) call io_error('Error allocating in_data in param_in_file')

  line_counter=0
  do loop=1,tot_num_lines
     read(in_unit, '(a)', iostat = ierr, err= 200 ) dummy
     ! [GP-begin, Apr13, 2012]: I convert all tabulation characters to spaces
     pos = index(dummy,TABCHAR)
     do while (pos .ne. 0)
        dummy(pos:pos) = ' '
        pos = index(dummy,TABCHAR)
     end do
     ! [GP-end]
     dummy=utility_lowercase(dummy)
     dummy=adjustl(dummy)
     if( dummy(1:1)=='!' .or.  dummy(1:1)=='#' ) cycle
     if(len(trim(dummy)) == 0 ) cycle
     line_counter=line_counter+1
     in1=index(dummy,'!')
     in2=index(dummy,'#')
     if(in1==0 .and. in2==0)  in_data(line_counter)=dummy
     if(in1==0 .and. in2>0 )  in_data(line_counter)=dummy(:in2-1)
     if(in2==0 .and. in1>0 )  in_data(line_counter)=dummy(:in1-1)
     if(in2>0 .and. in1>0 )   in_data(line_counter)=dummy(:min(in1,in2)-1)
  end do

  CLOSE(in_unit)
!------------------------------------------------------------------------
END SUBROUTINE param_in_file
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
SUBROUTINE param_get_range_vector(keyword,found,length,lcount,i_value)
!------------------------------------------------------------------------
!! Read a range vector eg. 1,2,3,4-10  or 1 3 400:100
!! if(lcount) we return the number of states in length
!------------------------------------------------------------------------
  IMPLICIT NONE

  character(*),      intent(in)    :: keyword
  !! Keyword to examine
  logical          , intent(out)   :: found
  !! Is keyword found
  INTEGER,           intent(inout) :: length
  !! Number of states
  logical,           intent(in)    :: lcount
  !! If T only count states
  INTEGER, optional, intent(out)   :: i_value(length)
  !! States specified in range vector

  INTEGER   :: kl, in,loop,num1,num2,i_punc
  INTEGER   :: counter,i_digit,loop_r,range_size
  character(len=maxlen) :: dummy
  character(len=10), parameter :: c_digit="0123456789"
  character(len=2) , parameter :: c_range="-:"
  character(len=3) , parameter :: c_sep=" ,;"
  character(len=5) , parameter :: c_punc=" ,;-:"
  character(len=5)  :: c_num1,c_num2


  if(lcount .and. present(i_value) ) call io_error('param_get_range_vector: incorrect call')

  kl=len_trim(keyword)

  found=.false.

  do loop=1,num_lines
     in=index(in_data(loop),trim(keyword))
     if (in==0 .or. in>1 ) cycle
     if (found) then
        call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
     endif
     found=.true.
     dummy=in_data(loop)(kl+1:)
     dummy=adjustl(dummy)
     if(.not. lcount) in_data(loop)(1:maxlen) = ' '
     if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
        dummy=dummy(2:)
        dummy=adjustl(dummy)
     end if
  end do

  if(.not. found) return

  counter=0
  if (len_trim(dummy)==0) call io_error('Error: keyword '//trim(keyword)//' is blank')
  dummy=adjustl(dummy)
  do
     i_punc=scan(dummy,c_punc)
     if(i_punc==0) call io_error('Error parsing keyword '//trim(keyword))
     c_num1=dummy(1:i_punc-1)
     read(c_num1,*,err=101,end=101) num1
     dummy=adjustl(dummy(i_punc:))
     !look for range
     if(scan(dummy,c_range)==1) then
        i_digit=scan(dummy,c_digit)
        dummy=adjustl(dummy(i_digit:))
        i_punc=scan(dummy,c_punc)
        c_num2=dummy(1:i_punc-1)
        read(c_num2,*,err=101,end=101) num2
        dummy=adjustl(dummy(i_punc:))
        range_size=abs(num2-num1)+1
        do loop_r=1,range_size
           counter=counter+1
           if(.not. lcount) i_value(counter)=min(num1,num2)+loop_r-1
        end do
     else
        counter=counter+1
        if(.not. lcount) i_value(counter)=num1
     end if

     if(scan(dummy,c_sep)==1) dummy=adjustl(dummy(2:))
     if(scan(dummy,c_range)==1) call io_error('Error parsing keyword '//trim(keyword)//' incorrect range')
     if(index(dummy,' ')==1) exit
  end do

  if(lcount) length=counter
  if(.not.lcount) then
     do loop=1,counter-1
        do loop_r=loop+1,counter
           if(i_value(loop)==i_value(loop_r)) &
              call io_error('Error parsing keyword '//trim(keyword)//' duplicate values')
        end do
      end do
  end if
  !
  RETURN
  !
101 CALL io_error('Error parsing keyword '//trim(keyword))
!------------------------------------------------------------------------
END SUBROUTINE param_get_range_vector
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
FUNCTION utility_lowercase(string)!
!------------------------------------------------------------------------  !                                 !
!! Takes a string and converts to lowercase characters
!------------------------------------------------------------------------
  IMPLICIT NONE

  character(len=*), intent(in) :: string
  character(len=maxlen) :: utility_lowercase

  INTEGER :: iA,iZ,idiff,ipos,ilett

  iA = ichar('A')
  iZ = ichar('Z')
  idiff = iZ-ichar('z')

  utility_lowercase = string

  do ipos=1,len(string)
     ilett = ichar(string(ipos:ipos))
     if ((ilett.ge.iA).and.(ilett.le.iZ)) &
          utility_lowercase(ipos:ipos)=char(ilett-idiff)
  enddo

  utility_lowercase = trim(adjustl(utility_lowercase))
  RETURN
!------------------------------------------------------------------------
END FUNCTION utility_lowercase
!------------------------------------------------------------------------
END MODULE parameters
