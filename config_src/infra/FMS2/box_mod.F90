module box_mod
  use iso_c_binding, only : c_ptr, c_loc, c_int, c_null_ptr
  use MOM_error_infra, only : MOM_err, FATAL
  implicit none
  private

  public :: Box_t, Box_c

  !< Box struct for C bridge
  type, bind(C) :: Box_C
     type(c_ptr) :: idxS   !< start index of a box
     type(c_ptr) :: idxE   !< end index of a box
  end type Box_C

  !< Class box_t defines an iteration index range
  type :: Box_T
     integer(c_int), allocatable :: idxS(:)  !< Start index of a box
     integer(c_int), allocatable :: idxE(:)  !< End index of a box
  contains
     procedure   :: safe_alloc   !< allocate an index box
     procedure   :: set          !< Sets the index range for the bocx
     procedure   :: free         !< deallocates index box
     procedure   :: to_c         !< Converts Box to C
     procedure   :: grow         !< Increase the bounds of a box in one dimension
                                 !! both extents of box are increased by a fixed amount
     procedure   :: growLo       !< Increase the start extent of the box in one dimension
     procedure   :: growHi       !< Increase the end extent of the box in one dimension
     procedure   :: shrink       !< Decrease the bounds of a box in one dimension
                                 !! both extents of box are decreased by a fixed amount
     procedure   :: write_binary !< Write a box_t to a binary file
     procedure   :: read_binary  !< Read a box_t from a binary file
  end type Box_T

contains

!< Read a box_t from a binary file
subroutine read_binary(this, unit)
  class(Box_t), intent(inout) :: this  !< The box_t variable to read from a binary file
  integer,      intent(in)    :: unit  !< The file unit

  integer :: rank

  ! --- Read rank ---
  read(unit) rank

  ! --- Null case ---
  if (rank == -1) then
    if (allocated(this%idxS)) deallocate(this%idxS)
    if (allocated(this%idxE)) deallocate(this%idxE)
    return
  endif

  ! --- Allocate ---
  if (allocated(this%idxS)) deallocate(this%idxS)
  if (allocated(this%idxE)) deallocate(this%idxE)

  allocate(this%idxS(rank))
  allocate(this%idxE(rank))

  ! --- Read bounds ---
  read(unit) this%idxS
  read(unit) this%idxE

end subroutine read_binary

!< Write a box_t to a binary file
subroutine write_binary(this, unit)
  class(Box_t), intent(in) :: this   !< The box_t variable to write to a binary file
  integer,      intent(in) :: unit   !< The file unit

  integer :: rank

  ! --- Handle unallocated pointers ---
  if (.not. allocated(this%idxS) .or. .not. allocated(this%idxE)) then
    rank = -1
    write(unit) rank
    return
  endif

  ! --- Determine rank ---
  rank = size(this%idxS)

  ! --- Consistency check ---
  if (size(this%idxE) /= rank) then
    call MOM_err(FATAL,"Box_t%write_binary: idxS/idxE size mismatch")
  endif

  ! --- Write rank ---
  write(unit) rank

  ! --- Write bounds arrays ---
  write(unit) this%idxS
  write(unit) this%idxE

end subroutine write_binary

!< Allocates an iteration box
subroutine safe_alloc(this,ndims)
  class(Box_t), intent(inout) :: this   !< The box to be allocated
  integer, intent(in) :: ndims          !< The number of dimension in the box

  ! If already allocated deallocate
  if(allocated(this%idxS)) deallocate(this%idxS)
  if(allocated(this%idxE)) deallocate(this%idxE)

  allocate(this%idxS(ndims), source=0)
  allocate(this%idxE(ndims), source=0)
end subroutine safe_alloc

!< Allocates an iteration box
subroutine free(this)
  class(Box_t), intent(inout) :: this   !< The box to be deallocated

  if(allocated(this%idxS)) deallocate(this%idxS)
  if(allocated(this%idxE)) deallocate(this%idxE)

end subroutine free

!< Set the extents of the iteration box
subroutine set(this,idxS,idxE)
  class(Box_t), intent(inout) :: this        !< The box to set
  integer, dimension(:), intent(in) :: idxS  !< The starting indices
  integer, dimension(:), intent(in) :: idxE  !< The ending indices

  if(allocated(this%idxS) .and. allocated(this%idxS)) then
    this%idxS(:)=idxS(:)
    this%idxE(:)=idxE(:)
  else
    call MOM_err(FATAL, "class(box_t)%set index box must first be allocated")
  endif

end subroutine set

!< Return a new box with expanded iteration extents
function grow(this,dim,n) result(new)
  class(Box_t), intent(in) :: this !< The iteration box to modify
  integer, intent(in)      :: dim  !< The dimension to grow
  integer, intent(in)      :: n    !< The length to grow
  type(Box_t) :: new

  ! Local variables
  integer ::rank

  if(allocated(new%idxS)) deallocate(new%idxS)
  if(allocated(new%idxE)) deallocate(new%idxE)

  rank = SIZE(this%idxS)
  allocate(new%idxS(rank),new%idxE(rank))
  new%idxS(:) = this%idxS(:)
  new%idxE(:) = this%idxE(:)
  new%idxS(dim) = new%idxS(dim)-n
  new%idxE(dim) = new%idxE(dim)+n

end function grow

!< Return a new box with expanded start iteration extents
function growLo(this,dim,n) result(new)
  class(Box_t), intent(in) :: this !< The iteration box to modify
  integer, intent(in)      :: dim  !< The dimension to grow
  integer, intent(in)      :: n    !< The length to grow
  type(Box_t) :: new

  ! Local variables
  integer ::rank

  if(allocated(new%idxS)) deallocate(new%idxS)
  if(allocated(new%idxE)) deallocate(new%idxE)

  rank = SIZE(this%idxS)
  allocate(new%idxS(rank),new%idxE(rank))
  new%idxS(:) = this%idxS(:)
  new%idxE(:) = this%idxE(:)
  new%idxS(dim) = new%idxS(dim)-n
  new%idxE(dim) = new%idxE(dim)

end function growLo

!< Return a new box with expanded end iteration extents
function growHi(this,dim,n) result(new)
  class(Box_t), intent(in) :: this !< The iteration box to modify
  integer, intent(in)      :: dim  !< The dimension to grow
  integer, intent(in)      :: n    !< The length to grow
  type(Box_t) :: new

  ! Local variables
  integer ::rank

  if(allocated(new%idxS)) deallocate(new%idxS)
  if(allocated(new%idxE)) deallocate(new%idxE)

  rank = SIZE(this%idxS)
  allocate(new%idxS(rank),new%idxE(rank))
  new%idxS(:) = this%idxS(:)
  new%idxE(:) = this%idxE(:)
  new%idxS(dim) = new%idxS(dim)
  new%idxE(dim) = new%idxE(dim)+n

end function growHi

!< Return a new box with contracted iteration extents
function shrink(this,dim,n) result(new)
  class(Box_t), intent(in) :: this   !< The iteration box to modify
  integer, intent(in)      :: dim    !< The dimension to shrink
  integer, intent(in)      :: n      !< The length to shrink
  type(Box_t) :: new

  ! Local variables
  integer ::rank

  if(allocated(new%idxS)) deallocate(new%idxS)
  if(allocated(new%idxE)) deallocate(new%idxE)

  rank = SIZE(this%idxS)
  allocate(new%idxS(rank),new%idxE(rank))
  new%idxS(:) = this%idxS(:)
  new%idxE(:) = this%idxE(:)
  new%idxS(dim) = new%idxS(dim)+n
  new%idxE(dim) = new%idxE(dim)-n

end function shrink

!< Convert Fortran box to C
function to_c(this) result(cdesc)
  class(Box_t), target, intent(in) :: this  !< The box to convert
  type(Box_C) :: cdesc              !< C compatible pointers

  integer(c_int), pointer :: p(:)

  cdesc%idxS = c_null_ptr
  cdesc%idxE = c_null_ptr

  if(allocated(this%idxS)) then
     p => this%idxS
     cdesc%idxS = c_loc(p(1))
  endif
  if(allocated(this%idxE)) then
     p => this%idxE
     cdesc%idxE = c_loc(p(1))
   endif
end function to_c

end module box_mod
