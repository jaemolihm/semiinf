module parameters
  !! This module contains parameters to control the actions of wannier90.
  !! Also routines to read the parameters and write them out again.

  USE comms
  integer, parameter :: in_unit = 7
  integer, parameter :: maxlen = 100

  real(DP), public, save :: hopping_tol
  integer, public, save :: MAX_N_ITER

  LOGICAL, public, save :: isslab ! is the calculation with surface reconstruction
  LOGICAL, public, save :: isslab_match ! is the calculation with surface + bulk matching
  LOGICAL, public, save :: isslab_hamil ! is the input hamiltonian from slab calculation
  LOGICAL, public, save :: isbulk_add_onsite ! is the calculation with bulk hamiltonian with onsite terms on the surface
  integer, public, save :: nbulk ! number of wannier basis for bulk principal layer
  integer, public, save :: nsurf ! number of wannier basis for surface principal layer
  integer, allocatable, public, save :: ind_0(:), ind_1(:), ind_2(:) ! indices for hij
  integer, public, save :: bulk_rz ! 1 or -1: append bulk layer along +z or -z direction

  ! add onsite energy
  real(DP), allocatable, public, save :: bulk_onsite_energy(:)

  ! energy sampling
  integer, public, save :: num_energy
  real(dp), public, save :: dos_e_min, dos_e_max, dos_e_step
  real(DP), allocatable, public, save :: energy(:)
  real(dp), public, save :: sigma ! small imaginary number for green function
  LOGICAL, public, save :: isspin ! if true, read _spnr.dat file and calculate spin

  ! kpoint common
  logical, public, save :: kpoint_is_path ! true if path, false if grid
  integer, public, save :: num_kpoint ! number of sampled k-points
  real(dp), allocatable, public, save :: plot_kpoint(:,:) ! size(3,num_kpoint)
  ! kpoint path
  integer, public :: num_paths
  integer, public :: bands_num_points
  integer, public :: bands_num_spec_points
  character(len=1), allocatable ::bands_label(:)
  real(dp), allocatable :: bands_spec_points(:,:)
  ! kpoint grid
  integer :: kpoint_grid_num(1:2)

  character(len=maxlen), save :: seedname
  character(len=maxlen), save :: input_filename

!  ! artificial, surface local e field
!  real(dp), public :: efield_amp
!  real(dp), public :: zcenter
  
  ! values NOT determined from input
  LOGICAL, public, save :: flag_converged
  INTEGER, public, save :: n_iter

  ! private
  integer :: num_lines
  integer :: itemp1, itemp2
  character(len=maxlen), allocatable :: in_data(:)
  character(len=maxlen) :: buffer


CONTAINS
  SUBROUTINE param_write
    WRITE(*,'("isslab= ", L)') isslab
    WRITE(*,'("isslab_match= ", L)') isslab_match
    WRITE(*,'("isslab_hamil= ", L)') isslab_hamil
    WRITE(*,'("isbulk_add_onsite= ", L)') isbulk_add_onsite
    WRITE(*,'("nsurf= ", I)') nsurf
    WRITE(*,'("nbulk= ", I)') nbulk
    WRITE(*,'("bulk_rz= ", I)') bulk_rz
    WRITE(*,'("num_energy= ",I)') num_energy
    WRITE(*,'("num_kpoint= ",I)') num_kpoint
  END SUBROUTINE
  !==================================================================!
  SUBROUTINE param_read
  !==================================================================!
  !                                                                  !
  !! Read parameters and calculate derived values                    
  !                                                                  !
  !===================================================================  
    ! use w90_constants, only : bohr, eps6
    ! use w90_utility,   only : utility_recip_lattice,utility_metric
    ! use w90_io,        only : io_error,io_file_unit,input_filename,post_proc_flag
    implicit none

    !local variables
    ! real(kind=dp)  :: real_lattice_tmp(3,3)
    ! integer :: nkp,i,j,n,k,itmp,i_temp2,eig_unit,loop,ierr,iv_temp(3), rows
    ! logical :: found2,lunits,chk_found
    ! character(len=6) :: spin_str
    ! real(kind=dp) :: cosa(3),rv_temp(3)

    logical :: found
    integer :: ierr, i_temp, loop_spts, loop_i, loop_j, counter
    real(dP) :: vec(3)
    real(kind=dp), allocatable :: xval(:), kpath_len(:)
    integer, allocatable :: kpath_pts(:)
    ! integer, allocatable, dimension(:,:) :: nnkpts_block
    ! integer, allocatable, dimension(:) :: nnkpts_idx

    call param_in_file

    call param_get_keyword('seedname',found,c_value=seedname)

!    call param_get_keyword('efield_amp',found,r_value=efield_amp)
!    call param_get_keyword('zcenter',found,r_value=zcenter)

    ! iteration parameter
    hopping_tol = 1.D-8
    call param_get_keyword('hopping_tol',found,r_value=hopping_tol)
    MAX_N_ITER = 40
    call param_get_keyword('max_n_iter',found,i_value=MAX_N_ITER)
    
    isspin = .false.
    call param_get_keyword('isspin',found,l_value=isspin)

    ! setting hamiltonian
    call param_get_keyword('isslab',found,l_value=isslab)
    call param_get_keyword('isslab_match',found,l_value=isslab_match)
    call param_get_keyword('isslab_hamil',found,l_value=isslab_hamil)
    call param_get_keyword('isbulk_add_onsite',found,l_value=isbulk_add_onsite)
    call param_get_keyword('nsurf',found,i_value=nsurf)
    call param_get_keyword('nbulk',found,i_value=nbulk)
    bulk_rz = 1
    call param_get_keyword('bulk_rz',found,i_value=bulk_rz)
    if (isslab .AND. .NOT.(isslab_hamil .or. isbulk_add_onsite .or. isslab_match)) &
      call io_error('if isslab, isslab_hamil or isslab_match or isbulk_add_onsite must be .true.')
    if (isbulk_add_onsite .AND. .NOT.(isslab)) &
      call io_error('if isbulk_add_onsite, isslab must be .true.')
    if (isslab .AND. (.NOT. isslab_match) .AND. (nsurf<=nbulk) .AND. .NOT.(isbulk_add_onsite)) &
      call io_error('if isslab, nsurf must be GREATER THAN nbulk')
    ! if ((.NOT.isslab_hamil) .AND. (nsurf/=nbulk)) &
    !   call io_error('if NOT isslab_hamil, nsurf must be equal to nbulk')
    IF (isslab) THEN
      ALLOCATE(ind_0(nsurf))
      call param_get_range_vector('ind_0',found,nsurf,.false.,ind_0)
    ELSE
      nsurf = nbulk
    END IF
    ALLOCATE(ind_1(nbulk))
    call param_get_range_vector('ind_1',found,nbulk,.false.,ind_1)
    ALLOCATE(ind_2(nbulk))
    call param_get_range_vector('ind_2',found,nbulk,.false.,ind_2)

    ! energy
    call param_get_keyword('sigma',found,r_value=sigma)
    call param_get_keyword('e_min',found,r_value=dos_e_min)
    call param_get_keyword('e_max',found,r_value=dos_e_max)
    call param_get_keyword('e_step',found,r_value=dos_e_step)
    num_energy = NINT((dos_e_max - dos_e_min) / dos_e_step) 
    dos_e_step = (dos_e_max - dos_e_min) / (num_energy-1)
    ALLOCATE(energy(num_energy))
    DO i_temp = 1, num_energy
      energy(i_temp) = dos_e_min + dos_e_step * (i_temp-1)
    END DO

    ! kpoints
    call param_get_keyword('kpoint_type',found,c_value=buffer)
    IF (buffer .eq. 'path') THEN
        kpoint_is_path = .TRUE.
    ELSE IF (buffer .eq. 'grid') THEN
        kpoint_is_path = .FALSE.
    ELSE
        call io_error('kpoint_type is path or grid')
    END IF

    ! kpoints path
    IF (kpoint_is_path) THEN
        bands_num_spec_points=0
        call param_get_block_length('kpoint_path',found,i_temp)
        if (found) then
           bands_num_spec_points=i_temp*2
           allocate(bands_label(bands_num_spec_points),stat=ierr)
           if (ierr/=0) call io_error('Error allocating bands_label in param_read')
           allocate(bands_spec_points(3,bands_num_spec_points),stat=ierr)
           if (ierr/=0) call io_error('Error allocating bands_spec_points in param_read')
           call param_get_keyword_kpath
        end if
        call param_get_keyword('bands_num_points',found,i_value=bands_num_points)

        num_paths=bands_num_spec_points/2
        allocate(kpath_len(num_paths))
        allocate(kpath_pts(num_paths))
        ! num_spts=num_paths+1
        do loop_spts=1,num_paths
          vec=bands_spec_points(:,2*loop_spts)-bands_spec_points(:,2*loop_spts-1)
          ! kpath_len(loop_spts)=sqrt(dot_product(vec,(matmul(recip_metric,vec))))
          kpath_len(loop_spts) = 1.0_dp
          if(loop_spts==1) then
            kpath_pts(loop_spts)=bands_num_points
          else
            kpath_pts(loop_spts)=nint(real(bands_num_points,dp)*kpath_len(loop_spts)/kpath_len(1))
          end if
        end do
        num_kpoint=sum(kpath_pts)+1

        allocate(plot_kpoint(3,num_kpoint),stat=ierr)
        allocate(xval(num_kpoint),stat=ierr)

        counter=0
        do loop_spts=1,num_paths
          do loop_i=1,kpath_pts(loop_spts)
            counter=counter+1
            if(counter==1) then
              xval(counter)=0.0_dp
            else
              xval(counter)=xval(counter-1)+kpath_len(loop_spts)/real(kpath_pts(loop_spts),dp)
            endif
            plot_kpoint(:,counter)=bands_spec_points(:,2*loop_spts-1)+ &
                 (bands_spec_points(:,2*loop_spts)-bands_spec_points(:,2*loop_spts-1))* &
                 (real(loop_i-1,dp)/real(kpath_pts(loop_spts),dp))
          end do
        end do
        xval(num_kpoint)=sum(kpath_len)
        plot_kpoint(:,num_kpoint)=bands_spec_points(:,bands_num_spec_points)
    ELSE ! .NOT. kpoint_is_path, kpoint_type = 'grid'
        call param_get_keyword('grid_1',found,i_value=kpoint_grid_num(1))
        call param_get_keyword('grid_2',found,i_value=kpoint_grid_num(2))
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

    IF (isbulk_add_onsite) THEN
      ! a ad-hoc patch for GeTe...
      if (trim(seedname) .ne. 'gete.bulk') &
        call io_error('Error: isbulk_add_onsite NOT IMPLEMENTED except for GeTe.bulk')
      allocate(bulk_onsite_energy(nbulk))
      bulk_onsite_energy(1:8) = 0.157589417
      bulk_onsite_energy(9:16) = 0.152999
      bulk_onsite_energy(17:24) = 0.204782292
      bulk_onsite_energy(25:32) = 0.2090715
      bulk_onsite_energy(33:40) = 0.180956042
      bulk_onsite_energy(41:48) = 0.316390875
    END IF

  END SUBROUTINE

  !===================================!
  SUBROUTINE param_get_keyword_kpath
  !===================================!
  !                                   !
  !!  Fills the kpath data block     
  !                                   !
  !===================================!

    implicit none

    character(len=20) :: keyword
    integer           :: in,ins,ine,loop,i,line_e,line_s,counter
    logical           :: found_e,found_s
    character(len=maxlen) :: dummy,end_st,start_st

    keyword="kpoint_path"

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

    return


240 call io_error('param_get_keyword_kpath: Problem reading kpath '//trim(dummy))

  end SUBROUTINE param_get_keyword_kpath


  !=====================================================!
  SUBROUTINE param_get_block_length(keyword,found,rows,lunits)
  !=====================================================!
  !                                                     !
  !! Finds the length of the data block       
  !                                                     !
  !=====================================================!


    implicit none

    character(*),      intent(in)  :: keyword
    !! Keyword to examine
    logical,           intent(out) :: found
    !! Is keyword present
    integer,           intent(out) :: rows
    !! Number of rows
    logical, optional, intent(out) :: lunits
    !! Have we found a unit specification

    integer           :: i,in,ins,ine,loop,line_e,line_s
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


    return

555 lunits=.true.

    if(rows<=1) then !cope with empty blocks
       found=.false.
       in_data(line_s:line_e)(1:maxlen) = ' '
    end if


    return

  end SUBROUTINE param_get_block_length






  !===========================================================================!
  SUBROUTINE param_get_keyword(keyword,found,c_value,l_value,i_value,r_value)
    !===========================================================================!
    !                                                                           !
    !! Finds the value of the required keyword.
    !                                                                           !
    !===========================================================================!
    implicit none

    character(*),      intent(in)  :: keyword
    !! Keyword to examine
    logical          , intent(out) :: found
    !! Is keyword present
    character(*)     ,optional, intent(inout) :: c_value
    !! Keyword value
    logical          ,optional, intent(inout) :: l_value
    !! Keyword value
    integer          ,optional, intent(inout) :: i_value
    !! Keyword value
    real(kind=dp)    ,optional, intent(inout) :: r_value
    !! Keyword value

    integer           :: kl, in,loop,itmp
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

    return

220 call io_error('Error: Problem reading keyword '//trim(keyword))

  end SUBROUTINE param_get_keyword


  !=======================================!
  SUBROUTINE param_in_file
  !=======================================!
  !! Load the *.win file into a character  
  !! array in_file, ignoring comments and  
  !! blank lines and converting everything 
  !! to lowercase characters               
  !=======================================!
    IMPLICIT NONE

    integer           :: in_unit,tot_num_lines,ierr,line_counter,loop,in1,in2
    character(len=maxlen) :: dummy
    integer           :: pos
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

    close(in_unit)

  end SUBROUTINE param_in_file

    !====================================================================!
    subroutine param_get_range_vector(keyword,found,length,lcount,i_value)
    !====================================================================!
    !!   Read a range vector eg. 1,2,3,4-10  or 1 3 400:100           
    !!   if(lcount) we return the number of states in length            
    !====================================================================!
    implicit none

    character(*),      intent(in)    :: keyword
    !! Keyword to examine
    logical          , intent(out)   :: found
    !! Is keyword found
    integer,           intent(inout) :: length
    !! Number of states
    logical,           intent(in)    :: lcount
    !! If T only count states
    integer, optional, intent(out)   :: i_value(length)
    !! States specified in range vector

    integer   :: kl, in,loop,num1,num2,i_punc
    integer   :: counter,i_digit,loop_r,range_size
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

    return

101 call io_error('Error parsing keyword '//trim(keyword))


   end  subroutine param_get_range_vector


  !=================================!
  function utility_lowercase(string)!
  !=================================!
  !                                 !
  !! Takes a string and converts to
  !!  lowercase characters
  !                                 !
  !=================================!
    implicit none

    character(len=*), intent(in) :: string
    character(len=maxlen) :: utility_lowercase

    integer :: iA,iZ,idiff,ipos,ilett

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

    return
  end function utility_lowercase
end module parameters
