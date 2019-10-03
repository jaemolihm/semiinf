!------------------------------------------------------------------------
MODULE parameters
!------------------------------------------------------------------------
!! This module contains parameters to control the actions of wannier90.
!! Also routines to read the parameters and write them out again.
!------------------------------------------------------------------------
  USE comms
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: in_unit = 7
  INTEGER, PARAMETER :: maxlen = 256

  REAL(DP), PUBLIC, SAVE :: hopping_tol
  INTEGER, PUBLIC, SAVE :: MAX_N_ITER

  LOGICAL, PUBLIC, SAVE :: isslab
  !! is the calculation with surface reconstruction
  LOGICAL, PUBLIC, SAVE :: isslab_match
  !! is the calculation with surface + bulk matching
  LOGICAL, PUBLIC, SAVE :: isslab_hamil
  !! is the input hamiltonian from slab calculation
  INTEGER, PUBLIC, SAVE :: nbulk
  !! number of wannier basis for bulk principal layer
  INTEGER, PUBLIC, SAVE :: nsurf
  !! number of wannier basis for surface principal layer
  INTEGER, ALLOCATABLE, PUBLIC, SAVE :: ind_0(:)
  !! indices for hij
  INTEGER, ALLOCATABLE, PUBLIC, SAVE :: ind_1(:)
  !! indices for hij
  INTEGER, ALLOCATABLE, PUBLIC, SAVE :: ind_2(:)
  !! indices for hij
  INTEGER, PUBLIC, SAVE :: bulk_rz
  !! 1 or -1: append bulk layer along +z or -z direction
  !
  ! energy sampling
  INTEGER, PUBLIC, SAVE :: num_energy
  REAL(DP), PUBLIC, SAVE :: dos_e_min, dos_e_max, dos_e_step
  REAL(DP), ALLOCATABLE, PUBLIC, SAVE :: energy(:)
  REAL(DP), PUBLIC, SAVE :: sigma
  !! small imaginary number for green function
  LOGICAL, PUBLIC, SAVE :: isspin
  !! if true, read _spnr.dat file and calculate spin
  !
  ! kpoint common
  LOGICAL, PUBLIC, SAVE :: kpoint_is_path
  !! true if path, false if grid
  INTEGER, PUBLIC, SAVE :: num_kpoint
  !! number of sampled k-points
  REAL(DP), ALLOCATABLE, PUBLIC, SAVE :: plot_kpoint(:,:)
  !! (3, num_kpoint) k points to calculate the spectral function
  ! kpoint path
  INTEGER, PUBLIC :: num_paths
  INTEGER, PUBLIC :: bands_num_points
  INTEGER, PUBLIC :: bands_num_spec_points
  CHARACTER(LEN=1), ALLOCATABLE ::bands_label(:)
  REAL(DP), ALLOCATABLE :: bands_spec_points(:,:)
  !! kpoint grid
  INTEGER :: kpoint_grid_num(2)

  CHARACTER(LEN=maxlen), SAVE :: seedname
  CHARACTER(LEN=maxlen), SAVE :: input_filename

!  ! artificial, surface local e field
!  REAL(DP), PUBLIC :: efield_amp
!  REAL(DP), PUBLIC :: zcenter

  ! values NOT determined from input
  LOGICAL, PUBLIC, SAVE :: flag_converged
  INTEGER, PUBLIC, SAVE :: n_iter

  ! private
  INTEGER :: num_lines
  INTEGER :: itemp1, itemp2
  CHARACTER(LEN=maxlen), ALLOCATABLE :: in_data(:)
  CHARACTER(LEN=maxlen) :: buffer


CONTAINS
  !------------------------------------------------------------------------
  SUBROUTINE param_write
    WRITE(*,'("isslab= ", L)') isslab
    WRITE(*,'("isslab_match= ", L)') isslab_match
    WRITE(*,'("isslab_hamil= ", L)') isslab_hamil
    WRITE(*,'("nsurf= ", I)') nsurf
    WRITE(*,'("nbulk= ", I)') nbulk
    WRITE(*,'("bulk_rz= ", I)') bulk_rz
    WRITE(*,'("num_energy= ",I)') num_energy
    WRITE(*,'("num_kpoint= ",I)') num_kpoint
  END SUBROUTINE
  !------------------------------------------------------------------------
  SUBROUTINE param_read
  !------------------------------------------------------------------------
  !! Read input parameters and calculate derived values
  !------------------------------------------------------------------------
    IMPLICIT NONE
    LOGICAL :: found
    INTEGER :: ierr, i_temp, loop_spts, loop_i, loop_j, counter
    REAL(DP) :: vec(3)
    REAL(DP), ALLOCATABLE :: xval(:), kpath_len(:)
    INTEGER, ALLOCATABLE :: kpath_pts(:)

    CALL param_in_file

    CALL param_get_keyword('seedname',found,c_value=seedname)

!    call param_get_keyword('efield_amp',found,r_value=efield_amp)
!    call param_get_keyword('zcenter',found,r_value=zcenter)

    ! default values
    hopping_tol = 1.D-8
    MAX_N_ITER = 40
    isspin = .false.
    bulk_rz = 1
    ! read parameters
    CALL param_get_keyword('hopping_tol', found, r_value=hopping_tol)
    CALL param_get_keyword('max_n_iter',  found,i_value=MAX_N_ITER)
    CALL param_get_keyword('isspin',      found,l_value=isspin)
    CALL param_get_keyword('isslab',      found,l_value=isslab)
    CALL param_get_keyword('isslab_match',found,l_value=isslab_match)
    CALL param_get_keyword('isslab_hamil',found,l_value=isslab_hamil)
    CALL param_get_keyword('nsurf',       found,i_value=nsurf)
    CALL param_get_keyword('nbulk',       found,i_value=nbulk)
    CALL param_get_keyword('bulk_rz',     found,i_value=bulk_rz)
    !
    ! check validity of the input parameters
    if (isslab .AND. .NOT.(isslab_hamil .OR. isslab_match)) &
      CALL io_error('if isslab, isslab_hamil or isslab_match must be .true.')
    if (isslab .AND. (.NOT. isslab_match) .AND. (nsurf <= nbulk)) &
      CALL io_error('if isslab, nsurf must be GREATER THAN nbulk')
    ! if ((.NOT.isslab_hamil) .AND. (nsurf/=nbulk)) &
    !   CALL io_error('if NOT isslab_hamil, nsurf must be equal to nbulk')
    IF (isslab) THEN
      ALLOCATE(ind_0(nsurf))
      CALL param_get_range_vector('ind_0',found,nsurf,.false.,ind_0)
    ELSE
      nsurf = nbulk
    END IF
    ALLOCATE(ind_1(nbulk))
    CALL param_get_range_vector('ind_1',found,nbulk,.false.,ind_1)
    ALLOCATE(ind_2(nbulk))
    CALL param_get_range_vector('ind_2',found,nbulk,.false.,ind_2)

    ! energy
    CALL param_get_keyword('sigma', found, r_value=sigma)
    CALL param_get_keyword('e_min', found, r_value=dos_e_min)
    CALL param_get_keyword('e_max', found, r_value=dos_e_max)
    CALL param_get_keyword('e_step', found, r_value=dos_e_step)
    num_energy = NINT((dos_e_max - dos_e_min) / dos_e_step)
    dos_e_step = (dos_e_max - dos_e_min) / (num_energy-1)
    ALLOCATE(energy(num_energy))
    DO i_temp = 1, num_energy
      energy(i_temp) = dos_e_min + dos_e_step * (i_temp-1)
    ENDDO

    ! kpoints
    call param_get_keyword('kpoint_type',found,c_value=buffer)
    IF (buffer == 'path') THEN
      kpoint_is_path = .TRUE.
    ELSE IF (buffer == 'grid') THEN
      kpoint_is_path = .FALSE.
    ELSE
      CALL io_error('kpoint_type is path or grid')
    ENDIF

    ! set kpoints path
    IF (kpoint_is_path) THEN
      bands_num_spec_points=0
      CALL param_get_block_length('kpoint_path',found,i_temp)
      IF (found) THEN
        bands_num_spec_points = i_temp * 2
        ALLOCATE(bands_label(bands_num_spec_points),stat=ierr)
        IF (ierr/=0) CALL io_error('Error allocating bands_label in param_read')
        ALLOCATE(bands_spec_points(3,bands_num_spec_points),stat=ierr)
        IF (ierr/=0) CALL io_error('Error allocating bands_spec_points in param_read')
        CALL param_get_keyword_kpath
      ENDIF
      CALL param_get_keyword('bands_num_points',found,i_value=bands_num_points)

      num_paths = bands_num_spec_points / 2
      ALLOCATE(kpath_len(num_paths))
      ALLOCATE(kpath_pts(num_paths))
      ! num_spts=num_paths+1
      DO loop_spts = 1, num_paths
        vec = bands_spec_points(:, 2*loop_spts) - bands_spec_points(:, 2*loop_spts-1)
        ! kpath_len(loop_spts)=sqrt(dot_product(vec,(matmul(recip_metric,vec))))
        kpath_len(loop_spts) = 1.0_dp
        IF (loop_spts == 1) THEN
          kpath_pts(loop_spts)=bands_num_points
        ELSE
          kpath_pts(loop_spts)=nint(real(bands_num_points,dp)*kpath_len(loop_spts)/kpath_len(1))
        ENDIF
      ENDDO
      num_kpoint = SUM(kpath_pts)+1

      ALLOCATE(plot_kpoint(3,num_kpoint), stat=ierr)
      ALLOCATE(xval(num_kpoint), stat=ierr)

      counter=0
      DO loop_spts=1,num_paths
        DO loop_i=1,kpath_pts(loop_spts)
          counter = counter+1
          IF (counter == 1) THEN
            xval(counter)=0.0_dp
          ELSE
            xval(counter)=xval(counter-1)+kpath_len(loop_spts)/real(kpath_pts(loop_spts),dp)
          ENDIF
          plot_kpoint(:,counter)=bands_spec_points(:,2*loop_spts-1)+ &
               (bands_spec_points(:,2*loop_spts)-bands_spec_points(:,2*loop_spts-1))* &
               (real(loop_i-1,dp)/real(kpath_pts(loop_spts),dp))
        ENDDO
      ENDDO
      xval(num_kpoint)=sum(kpath_len)
      plot_kpoint(:,num_kpoint)=bands_spec_points(:,bands_num_spec_points)
    ELSE
      ! .NOT. kpoint_is_path, kpoint_type = 'grid'
      CALL param_get_keyword('grid_1',found,i_value=kpoint_grid_num(1))
      CALL param_get_keyword('grid_2',found,i_value=kpoint_grid_num(2))
      num_kpoint = kpoint_grid_num(1) * kpoint_grid_num(2)
      DO loop_j = 1, kpoint_grid_num(2)
        DO loop_i = 1, kpoint_grid_num(1)
          plot_kpoint(1,(loop_j-1)*kpoint_grid_num(1)+loop_i) = &
              MOD(REAL(loop_i,dp) / REAL(kpoint_grid_num(1),dp) + 0.5_dp, 1.0_dp) - 0.5_dp
          plot_kpoint(2,(loop_j-1)*kpoint_grid_num(1)+loop_i) = &
              MOD(REAL(loop_j,dp) / REAL(kpoint_grid_num(2),dp) + 0.5_dp, 1.0_dp) - 0.5_dp
        END DO
      END DO
      plot_kpoint(3,:) = 0.0_dp
    END IF
  !------------------------------------------------------------------------
  END SUBROUTINE param_read
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE param_get_keyword_kpath
  !------------------------------------------------------------------------
  !!  Fills the kpath data block
  !------------------------------------------------------------------------
    IMPLICIT NONE
    LOGICAL :: found_e,found_s
    CHARACTER(LEN=20) :: keyword
    CHARACTER(LEN=maxlen) :: dummy,end_st,start_st
    INTEGER :: in,ins,ine,loop,i,line_e,line_s,counter

    keyword = "kpoint_path"

    found_s = .false.
    found_e = .false.

    start_st='begin '//trim(keyword)
    end_st='end '//trim(keyword)

    do loop=1,num_lines
      ins=index(in_data(loop),trim(keyword))
      if (ins==0 ) cycle
      in=index(in_data(loop),'begin')
      if (in==0 .or. in>1) cycle
      line_s=loop
      if (found_s) then
        call io_error('Error: Found '//trim(start_st)//' more than once in input file')
      endif
      found_s=.true.
    end do

    do loop=1,num_lines
      ine=index(in_data(loop),trim(keyword))
      if (ine==0 ) cycle
      in=index(in_data(loop),'end')
      if (in==0 .or. in>1) cycle
      line_e=loop
      if (found_e) then
        call io_error('Error: Found '//trim(end_st)//' more than once in input file')
      endif
      found_e=.true.
    end do

    if(.not. found_e) then
      call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if(line_e<=line_s) then
      call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    counter=0
    do loop=line_s+1,line_e-1
      counter=counter+2
      dummy=in_data(loop)
      read(dummy,*,err=240,end=240) bands_label(counter-1),(bands_spec_points(i,counter-1),i=1,3)&
           ,bands_label(counter),(bands_spec_points(i,counter),i=1,3)
    end do

    in_data(line_s:line_e)(1:maxlen) = ' '

    RETURN

240 CALL io_error('param_get_keyword_kpath: Problem reading kpath '//trim(dummy))
  !------------------------------------------------------------------------
  END SUBROUTINE param_get_keyword_kpath
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE param_get_block_length(keyword,found,rows,lunits)
  !------------------------------------------------------------------------
  !! Finds the length of the data block
  !------------------------------------------------------------------------
    IMPLICIT NONE

    character(*),      intent(in)  :: keyword
    !! Keyword to examine
    logical,           intent(out) :: found
    !! Is keyword present
    INTEGER,           intent(out) :: rows
    !! Number of rows
    logical, optional, intent(out) :: lunits
    !! Have we found a unit specification

    INTEGER           :: i,in,ins,ine,loop,line_e,line_s
    logical           :: found_e,found_s
    character(len=maxlen) :: end_st,start_st,dummy
    character(len=2)  :: atsym
    real(kind=dp)     :: atpos(3)

    rows=0
    found_s=.false.
    found_e=.false.

    start_st='begin '//trim(keyword)
    end_st='end '//trim(keyword)

    do loop=1,num_lines
       ins=index(in_data(loop),trim(keyword))
       if (ins==0 ) cycle
       in=index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s=loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s=.true.
    end do


    if(.not. found_s) then
       found=.false.
       return
    end if


    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0 ) cycle
       in=index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e=loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e=.true.
    end do


    if(.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if(line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    rows=line_e-line_s-1

    found=.true.

    if (present(lunits)) then
       dummy=in_data(line_s+1)
       !       write(stdout,*) dummy
       !       write(stdout,*) trim(dummy)
       read(dummy,*,end=555) atsym, (atpos(i),i=1,3)
       lunits=.false.
    endif

    if(rows<=0) then !cope with empty blocks
       found=.false.
       in_data(line_s:line_e)(1:maxlen) = ' '
    end if

    RETURN

555 lunits=.true.

    if(rows<=1) then !cope with empty blocks
       found=.false.
       in_data(line_s:line_e)(1:maxlen) = ' '
    end if
    !
    RETURN
  !------------------------------------------------------------------------
  END SUBROUTINE param_get_block_length
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
    IMPLICIT NONE

    INTEGER           :: in_unit,tot_num_lines,ierr,line_counter,loop,in1,in2
    character(len=maxlen) :: dummy
    INTEGER           :: pos
    character, parameter :: TABCHAR = char(9)

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
