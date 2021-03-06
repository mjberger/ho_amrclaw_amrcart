! ============================================================================
!  Program:     AMRClaw
!  File:        resize_storage.f90
!  Created:     2009-01-21
!  Author:      Kyle Mandli and Marsha Berger
! ============================================================================
!  Description:  Resize the alloc array for AMR storage
! ============================================================================


! NOTE:  Older f90 compilers (e.g. gfortran prior to 4.2?)
! may not implement move_alloc.  If this fails, you may need to use
! resize_storage_static.f90 instead of this routine and set the
! allocation large enough in init_alloc.f90 to avoid running out of space.


subroutine resize_storage(new_size,status)
    
    use amr_module
    implicit none
    
    integer, intent(out) :: status
    integer, intent(in) :: new_size
    
    real(kind=8), allocatable, target, dimension(:) :: new_storage
    

    print *, "Expanding storage from ", memsize," to ", new_size
    allocate(new_storage(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_storage(1:memsize) = alloc     ! new way, use allocatable, not pointer       

    call move_alloc(new_storage,alloc)

    memsize = new_size
    
    return
    
end subroutine resize_storage
