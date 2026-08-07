module turbotmp_helperF

  use MOM_string_functions, only : uppercase
  use MOM_error_infra, only : MOM_err, FATAL
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use iso_c_binding, only : c_int
  use array_mod, only : RealArray_t, LogicalArray_t
  use box_mod, only : Box_t

  implicit none

  logical, parameter :: use_AMREX = .TRUE.
  integer, parameter :: TIMH_runAMREX   = 101, &
                        TIMH_capture    = 102, &
                        TIMH_runFORTRAN = 103

   public :: TIMH_runAMREX, TIMH_capture, TIMH_runFORTRAN
   public :: getenv_mode

   public :: already_recorded, mark_recorded

   integer, parameter :: max_kernels = 1000
   character(len=128), save :: recorded(max_kernels)
   integer, save :: n_recorded = 0
   integer, parameter :: type_read = 1, type_write = 2

   !> Metadata describing a variable in a binary I/O stream
   type :: io_entry
     character(len=:), allocatable :: name         !< The name of a variable on which to perform IO
     character(len=:), allocatable :: type_name    !< The type of the variable
     integer(kind=int64) :: offset                 !< Byte offset in the binary output file
   end type io_entry

   !> An I/O capture capability for individual variables
   type :: io_recorder
     integer :: unit_bin                       !< File unit for  binary file
     integer :: unit_meta                      !< File unit for  metadata file
     type(io_entry), allocatable :: entries(:) !< A description of each of the variables in the file
     integer :: n = 0                          !< Counter for number of variables in a file
     integer :: type                           !< Indicates if an open operation is
                                               !< for reading or writing
     contains
       procedure :: open_write       !< open metadata and binary files for writing
       procedure :: open_read        !< open metadata and binary files for reading
       procedure :: add_entry        !< Add variable descriptor to the metadata file
       procedure :: close            !< close the binary and metadata file

       ! Write side
       procedure :: add_realarray    !< Write a RealArray_t variable to the capture file
       procedure :: add_logicalarray !< Write a LogicalArray_t variable to the capture file
       procedure :: add_box          !< Write a Box_t variable to the capture file
       procedure :: add_real         !< Write a real scalar variable to the capture file
       procedure :: add_integer      !< Write an integer scalar variable to the capture file
       procedure :: add_logical      !< Write a logical scalar variable to the capture file

       ! Read side
       procedure :: get_realarray    !< Read a RealArray_t variable from the capture file
       procedure :: get_logicalarray !< Read a LogicalArray_t variable from the capture file
       procedure :: get_box          !< Read a Box_t variable from the capture file
       procedure :: get_real         !< Read a real scalar variable from the capture file
       procedure :: get_integer      !< Read an integer scalar variable from the capture file
       procedure :: get_logical      !< Read a logical scalar variable from the capture file
       procedure :: load_metadata    !< Read in the metadata file
       procedure :: find_entry       !< Query the locaiton of the variable in the binary
                                     !! capture file
       !> Generic interface to add a variable
       generic   :: add => add_realarray, add_logicalarray, add_box, add_real, add_integer, &
                            add_logical
       !> Generic interface to get a variable
       generic   :: get => get_realarray, get_logicalarray, get_box, get_real, get_integer, &
                            get_logical
   end type io_recorder

contains

  !< read environment variables that control control flow for the Shim layer
  function getenv_mode(name, default) result(mode)
    character(len=*), intent(in) :: name         !< The name of the environment variable
    integer, intent(in), optional :: default     !< The default value if environment variable not set
    integer :: mode

    character(len=:), allocatable :: str
    integer :: length, status
    character(len=256) :: mesg

    ! Get length first
    call get_environment_variable(name, length=length, status=status)

    if (status /= 0 .or. length == 0) then
      if (present(default)) then
        mode = default
      else
        mode = TIMH_runFORTRAN  ! sensible fallback
      end if
      return
    end if

    allocate(character(len=length) :: str)
    call get_environment_variable(name, str)

    ! Normalize (optional but recommended)
    str= uppercase(str)

    select case (trim(str))
    case ("AMREX")
      mode = TIMH_runAMREX
    case ("FORTRAN")
      mode = TIMH_runFORTRAN
    case ("CAPTURE")
      mode = TIMH_capture
    case default
      write(mesg,'("tim_helperF::getenv_mod called with a ",A, &
               & " allowed values: AMREX, FORTRAN, CAPTURE")') TRIM(str)
      call MOM_err(FATAL,mesg)
    end select

  end function getenv_mode


  !< function that determins if a variable name has already been recorded
  logical function already_recorded(name)
    character(len=*), intent(in) :: name        !< Name of the variable to check status
    integer :: i

    already_recorded = .false.
    do i = 1, n_recorded
      if (trim(recorded(i)) == trim(name)) then
        already_recorded = .true.
        return
      end if
    end do
  end function already_recorded

  !< Mark that a particular variable has already been recorded
  subroutine mark_recorded(name)
    character(len=*), intent(in) :: name      !< Name of the variable to mark as areadly captured

    if (.not. already_recorded(name)) then
      if (n_recorded < max_kernels) then
        n_recorded = n_recorded + 1
        recorded(n_recorded) = name
      else
        print *, "Recorder registry full!"
      end if
    end if
  end subroutine mark_recorded

!< Record and write a variable of type real64
subroutine add_real(this, name, val)
  class(io_recorder), intent(inout) :: this    !< The state recorder class
  character(*), intent(in) :: name             !< The name of the variable to write
  real(kind=real64), intent(in) :: val         !< The variable to write

  integer(kind=int64) :: pos

  inquire(unit=this%unit_bin, pos=pos)

  call this%add_entry(name, 'real64', pos)

  write(this%unit_bin) val
end subroutine add_real

!< Read a variable of type real64
subroutine get_real(this, name, val)
  class(io_recorder), intent(inout) :: this !< The state recorder class
  character(*), intent(in) :: name          !< The name of the variable to read
  real(kind=real64), intent(out) :: val     !< The variable to read

  ! local variables
  integer :: idx
  integer(kind=int64) :: pos
  character(len=256) :: mesg

  idx = this%find_entry(name)
  if (idx < 0) then
     write(mesg,'("tim_helperF::get_real variable ",A," not found ")') TRIM(name)
     call MOM_err(FATAL,mesg)
  endif

  pos = this%entries(idx)%offset

  read(this%unit_bin, pos=pos) val
end subroutine get_real

!< Record and write a variable of type integer
subroutine add_integer(this, name, val)
  class(io_recorder), intent(inout) :: this  !< The state recorder class
  character(*), intent(in) :: name           !< The name of the variable to write
  integer, intent(in) :: val                 !< The variable to write

  integer(kind=int64) :: pos

  inquire(unit=this%unit_bin, pos=pos)

  call this%add_entry(name, 'integer', pos)

  write(this%unit_bin) val
end subroutine add_integer

!< Read a variable of type integer
subroutine get_integer(this, name, val)
  class(io_recorder), intent(inout) :: this   !< The state recorder class
  character(*), intent(in) :: name            !< The name of the variable to read
  integer, intent(out) :: val                 !< The variable to read

  ! local variables
  integer :: idx
  integer(kind=int64) :: pos
  character(len=256) :: mesg

  idx = this%find_entry(name)
  if (idx < 0) then
     write(mesg,'("tim_helperF::get_integer variable ",A," not found ")') TRIM(name)
     call MOM_err(FATAL,mesg)
  endif

  pos = this%entries(idx)%offset

  read(this%unit_bin, pos=pos) val
end subroutine get_integer

!< Record and write a variable of type logical
subroutine add_logical(this, name, val)
  class(io_recorder), intent(inout) :: this !< The state recorder class
  character(*), intent(in) :: name          !< The name of the variable to write
  logical, intent(in) :: val                !< The variable to write

  integer(kind=int64) :: pos

  inquire(unit=this%unit_bin, pos=pos)

  call this%add_entry(name, 'logical', pos)

  write(this%unit_bin) val
end subroutine add_logical

!< Read a variable of type logical
subroutine get_logical(this, name, val)
  class(io_recorder), intent(inout) :: this !< The state recorder class
  character(*), intent(in) :: name          !< The name of the variable to read
  logical, intent(out) :: val               !< The variable to read

  ! local variables
  integer :: idx
  integer(kind=int64) :: pos
  character(len=256) :: mesg

  idx = this%find_entry(name)
  if (idx < 0) then
     write(mesg,'("tim_helperF::get_logical variable ",A," not found ")') TRIM(name)
     call MOM_err(FATAL,mesg)
  endif

  pos = this%entries(idx)%offset

  read(this%unit_bin, pos=pos) val

end subroutine get_logical

!< Open a capture file for writing
subroutine open_write(this, binfile, metafile)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: binfile    !< The name of the binary file to open for writing
  character(*), intent(in) :: metafile   !< The name of the metadata file to open for writing

  open(newunit=this%unit_bin, file=binfile, &
       access='stream', form='unformatted', status='replace')

  open(newunit=this%unit_meta, file=metafile, &
       form='formatted', status='replace')

  this%type = type_write
  this%n = 0
end subroutine open_write

!< Open a capture file for reading
subroutine open_read(this, binfile, metafile)
  class(io_recorder), intent(inout) :: this
  character(*), intent(in) :: binfile    !< The name of the binary file to open for reading
  character(*), intent(in) :: metafile   !< The name of the metadata file to open for reading

  open(newunit=this%unit_bin, file=binfile, &
       access='stream', form='unformatted', status='old')

  open(newunit=this%unit_meta, file=metafile, &
       form='formatted', status='old')

  call this%load_metadata()

  this%type = type_read

end subroutine open_read

!< Load a metadata file
subroutine load_metadata(this)
  class(io_recorder), intent(inout) :: this   !< The state recorder class

  ! local variables
  character(len=128) :: name, type_name
  integer(kind=int64) :: offset

  this%n = 0
  if (allocated(this%entries)) deallocate(this%entries)

  do
    read(this%unit_meta, *, end=100) name, type_name, offset

    this%n = this%n + 1
    if (.not. allocated(this%entries)) then
      allocate(this%entries(1))
    else
      this%entries = [this%entries, io_entry(name, type_name, offset)]
      cycle
    endif

    this%entries(1)%name = name
    this%entries(1)%type_name = type_name
    this%entries(1)%offset = offset
  enddo

  100 continue
end subroutine load_metadata

!< Find an entry in the list of variables
function find_entry(this, name) result(idx)
  class(io_recorder), intent(in) :: this       !< The state recorder class
  character(*), intent(in) :: name             !< The name of the variable

  ! local variables
  integer :: idx
  integer :: i

  idx = -1
  do i = 1, this%n
    if (trim(this%entries(i)%name) == trim(name)) then
      idx = i
      return
    endif
  enddo
end function find_entry

!< Add an entry into a list of variables
subroutine add_entry(this, name, type_name, offset)
  class(io_recorder), intent(inout) :: this     !< The state recorder class
  character(*), intent(in) :: name, type_name   !< The name of the variable
  integer(kind=int64), intent(in) :: offset     !< The offset into the binary file

  this%n = this%n + 1

  if (.not. allocated(this%entries)) then
    allocate(this%entries(1))
  else
    this%entries = [this%entries, io_entry(name, type_name, offset)]
    return
  endif

  this%entries(1)%name = name
  this%entries(1)%type_name = type_name
  this%entries(1)%offset = offset
end subroutine add_entry

!< Record and write a variable of type RealArray_t
subroutine add_realarray(this, name, val)
  class(io_recorder), intent(inout) :: this    !< The state recorder class
  character(*), intent(in) :: name             !< The name of the variable
  type(RealArray_t), intent(in) :: val         !< The RealArray_t array to write

  ! local variables
  integer(kind=int64) :: pos

  ! --- Get current file position ---
  inquire(unit=this%unit_bin, pos=pos)

  ! --- Register metadata ---
  call this%add_entry(name, 'RealArray_t', pos)

  ! --- Write binary payload ---
  call val%write_binary(this%unit_bin)

end subroutine add_realarray

!< Read a variable of type RealArray_t
subroutine get_realarray(this, name, val)
  class(io_recorder), intent(inout) :: this  !< The state recorder class
  character(*), intent(in) :: name           !< The name of the variable
  type(RealArray_t), intent(inout) :: val    !< The RealArray_t array to read

  ! local variables
  integer :: idx
  integer(kind=int64) :: pos
  character(len=256) :: mesg

  ! --- Find metadata entry ---
  idx = this%find_entry(name)
  if (idx < 0) then
    write(mesg,'("tim_helperF::get_realarray variable ",A," not found ")') TRIM(name)
    call MOM_err(FATAL,mesg)
  endif

  ! --- Optional type check ---
  if (trim(this%entries(idx)%type_name) /= 'RealArray_t') then
    write(mesg,'("tim_helperF::get_realarray variable ",A," type mismatch ")') TRIM(name)
    call MOM_err(FATAL,mesg)
  endif

  pos = this%entries(idx)%offset

  ! --- Seek to position ---
  read(this%unit_bin, pos=pos)

  ! --- Delegate to type ---
  call val%read_binary(this%unit_bin)

end subroutine get_realarray

!< Record and write a variable of type LogicalArray_t
subroutine add_logicalarray(this, name, val)
  class(io_recorder), intent(inout) :: this    !< The state recorder class
  character(*), intent(in) :: name             !< The name of the variable
  type(LogicalArray_t), intent(in) :: val      !< The LogicalArray_t array to write

  ! local variables
  integer(kind=int64) :: pos

  ! --- Get current file position ---
  inquire(unit=this%unit_bin, pos=pos)

  ! --- Register metadata ---
  call this%add_entry(name, 'LogicalArray_t', pos)

  ! --- Write binary payload ---
  call val%write_binary(this%unit_bin)

end subroutine add_logicalarray

!< Read a variable of type LogicalArray_t
subroutine get_logicalarray(this, name, val)
  class(io_recorder), intent(inout) :: this    !< The state recorder class
  character(*), intent(in) :: name             !< The name of the variable
  type(LogicalArray_t), intent(inout) :: val   !< The LogicalArray_t array to read

  ! local variables
  integer :: idx
  integer(kind=int64) :: pos
  character(len=256) :: mesg

  ! --- Find metadata entry ---
  idx = this%find_entry(name)
  if (idx < 0) then
    write(mesg,'("tim_helperF::get_logicalarray variable ",A," not found ")') TRIM(name)
    call MOM_err(FATAL,mesg)
  endif

  ! --- Optional type check ---
  if (trim(this%entries(idx)%type_name) /= 'LogicalArray_t') then
    write(mesg,'("tim_helperF::get_logicalarray variable ",A," type mismatch ")') TRIM(name)
    call MOM_err(FATAL,mesg)
  endif

  pos = this%entries(idx)%offset

  ! --- Seek to position ---
  read(this%unit_bin, pos=pos)

  ! --- Delegate to type ---
  call val%read_binary(this%unit_bin)

end subroutine get_logicalarray

!< Record and write a variable of type box_t
subroutine add_box(this, name, val)
  class(io_recorder), intent(inout) :: this     !< The state recorder class
  character(*), intent(in) :: name              !< The name of the variable
  type(Box_T), intent(in) :: val                !< The box_T to write

  integer(kind=int64) :: pos

  inquire(unit=this%unit_bin, pos=pos)

  call this%add_entry(name, 'Box_t', pos)
  call val%write_binary(this%unit_bin)
end subroutine add_box

!< Read a variable of type box_t
subroutine get_box(this, name, val)
  class(io_recorder), intent(inout) :: this   !< The state recorder class
  character(*), intent(in) :: name            !< The name of the variable
  type(Box_t), intent(inout) :: val           !< The box_T to read

  integer :: idx
  integer(kind=int64) :: pos
  character(len=256) :: mesg

  idx = this%find_entry(name)
  if (idx < 0) then
    write(mesg,'("tim_helperF::get_box variable ",A," not found ")') TRIM(name)
    call MOM_err(FATAL,mesg)
  endif

  if (trim(this%entries(idx)%type_name) /= 'Box_t') then
    write(mesg,'("tim_helperF::get_box variable ",A," type mismatch ")') TRIM(name)
    call MOM_err(FATAL,mesg)
  endif

  pos = this%entries(idx)%offset

  read(this%unit_bin, pos=pos)

  call val%read_binary(this%unit_bin)

end subroutine get_box

!< Close the capture file
subroutine close(this)
  class(io_recorder), intent(inout) :: this   !< The state recorder class

  ! Local array
  integer :: i

  ! Write out the metadata
  if(this%type.eq.type_write) then
    do i = 1, this%n
      write(this%unit_meta,'(A,1X,A,1X,I0)') &
        this%entries(i)%name, &
        this%entries(i)%type_name, &
        this%entries(i)%offset
    enddo
  endif

  ! Close the files
  close(this%unit_bin)
  close(this%unit_meta)
end subroutine close

end module turbotmp_helperF
