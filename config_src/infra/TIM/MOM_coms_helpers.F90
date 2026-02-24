!> Thin interfaces to non-domain-oriented mpp communication subroutines
module MOM_coms_helpers

! This file is part of MOM6. See LICENSE.md for the license.

use iso_fortran_env, only : int32, int64

use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes, mpp_set_root_pe
use mpp_mod, only : mpp_set_current_pelist, mpp_get_current_pelist
use mpp_mod, only : mpp_sync, mpp_sync_self

implicit none ; private

public :: PE_here, root_PE, num_PEs, is_root_PE, set_rootPE, Set_PElist, Get_PElist, sync_PEs

contains

!> Return the ID of the PE for the current process.
function PE_here() result(pe)
  integer :: pe   !< PE ID of the current process
  pe = mpp_pe()
end function PE_here

!> Return the ID of the root PE for the PE list of the current procss.
function root_PE() result(pe)
  integer :: pe   !< root PE ID
  pe = mpp_root_pe()
end function root_PE

!> Return the number of PEs for the current PE list.
function num_PEs() result(npes)
  integer :: npes   !< Number of PEs
  npes = mpp_npes()
end function num_PEs

!> Designate a PE as the root PE
subroutine set_rootPE(pe)
  integer, intent(in) :: pe   !< ID of the PE to be assigned as root
  call mpp_set_root_pe(pe)
end subroutine

!> is_root_pe returns .true. if the current PE is the root PE.
logical function is_root_pe()
  is_root_pe = .false.
  if (PE_here() == root_PE()) is_root_pe = .true.
end function is_root_pe

!> Set the current PE list.  If no list is provided, then the current PE list
!! is set to the list of all available PEs on the communicator.  Setting the
!! list will trigger a rank synchronization unless the `no_sync` flag is set.
subroutine Set_PEList(pelist, no_sync)
  integer, optional, intent(in) :: pelist(:)  !< List of PEs to set for communication
  logical, optional, intent(in) :: no_sync    !< Do not sync after list update.
  call mpp_set_current_pelist(pelist, no_sync)
end subroutine Set_PEList

!> Retrieve the current PE list and any metadata if requested.
subroutine Get_PEList(pelist, name, commID)
  integer,                    intent(out) :: pelist(:) !< List of PE IDs of the current PE list
  character(len=*), optional, intent(out) :: name   !< Name of PE list
  integer,          optional, intent(out) :: commID !< Communicator ID of PE list

  call mpp_get_current_pelist(pelist, name, commiD)
end subroutine Get_PEList

!> Sync the PEs at a defined point in the model
subroutine sync_PEs(pelist)
  integer, optional, intent(in) :: pelist(:)  !< The list of PEs to be synced

  call mpp_sync(pelist)
end subroutine sync_PEs

end module MOM_coms_helpers
