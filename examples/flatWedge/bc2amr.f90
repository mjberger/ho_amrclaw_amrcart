! :::::::::: bc2amr ::::::::::::::::::::::::::::::::::::::::::::::;
!> \callgraph
!! \callergraph
!!  Take a grid patch with mesh widths **hx**,**hy**, of dimensions **nrow** by
!!  **ncol**,  and set the values of any piece of
!!  of the patch which extends outside the physical domain 
!!  using the boundary conditions. 
!!
!!  ### Standard boundary condition choices for amr2ez in clawpack
!!
!!  At each boundary  k = 1 (left),  2 (right),  3 (bottom), 4 (top):
!!
!!  mthbc(k) =  
!!  * 0  for user-supplied BC's (must be inserted!)
!!  * 1  for zero-order extrapolation
!!  * 2  for periodic boundary conditions
!!  * 3  for solid walls, assuming this can be implemented
!!                   by reflecting the data about the boundary and then
!!                   negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
!!                   component of q.
!!  * 4  for sphere bcs (left half maps to right half of same side, and vice versa), as if domain folded in half
!!
!!  The corners of the grid patch are at 
!!     (xlo_patch,ylo_patch)  --  lower left corner
!!     (xhi_patch,yhi_patch) --  upper right corner
!!
!!  The physical domain itself is a rectangle bounded by
!!     (xlower,ylower)  -- lower left corner
!!     (xupper,yupper)  -- upper right corner
!!  
!   This figure below does not work with doxygen
!   the picture is the following: 
!  ____________________________________________________
! 
!                _____________________ (xupper,yupper)
!               |                     |  
!           ____|____ (xhi_patch,yhi_patch)   
!           |   |    |                |
!           |   |    |                |
!           |   |    |                |
!           |___|____|                |
!  (xlo_patch,ylo_patch) |            |
!               |                     |
!               |_____________________|   
!    (xlower,ylower)
!  ____________________________________________________
!!
!!
!>  Any cells that lie outside the physical domain are ghost cells whose
!!  values should be set in this routine.  This is tested for by comparing
!!  xlo_patch with xlower to see if values need to be set at the left
!   as in the figure above, 
!
!>  and similarly at the other boundaries.
!!  Patches are guaranteed to have at least 1 row of cells filled
!!  with interior values so it is possible to extrapolate. 
!!  Fix [trimbd()](@ref trimbd) if you want more than 1 row pre-set.
!!
!!  Make sure the order the boundaries are specified is correct
!!  so that diagonal corner cells are also properly taken care of.
!!
!!  Periodic boundaries are set before calling this routine, so if the
!!  domain is periodic in one direction only you
!!  can safely extrapolate in the other direction. 
!!
!!  Don't overwrite ghost cells in periodic directions!
!!
!! \param val data array for solution \f$q \f$ (cover the whole grid **msrc**)
!! \param aux data array for auxiliary variables 
!! \param nrow number of cells in *i* direction on this grid
!! \param ncol number of cells in *j* direction on this grid
!! \param meqn number of equations for the system
!! \param naux number of auxiliary variables
!! \param hx spacing (mesh size) in *i* direction
!! \param hy spacing (mesh size) in *j* direction
!! \param level AMR level of this grid
!! \param time setting ghost cell values at time **time**
!! \param xlo_patch left bound of the input grid
!! \param xhi_patch right bound of the input grid 
!! \param ylo_patch lower bound of the input grid 
!! \param yhi_patch upper bound of the input grid 
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

subroutine bc2amr(val,aux,nrow,ncol,meqn,naux, hx, hy, level, time,   &
                  xlo_patch, xhi_patch, ylo_patch, yhi_patch) 

    use amr_module, only: mthbc, xlower, ylower, xupper, yupper
    use amr_module, only: xperdom,yperdom,spheredom

    implicit none

    ! Input/Output
    integer, intent(in) :: nrow, ncol, meqn, naux, level
    real(kind=8), intent(in) :: hx, hy, time
    real(kind=8), intent(in) :: xlo_patch, xhi_patch
    real(kind=8), intent(in) :: ylo_patch, yhi_patch
    real(kind=8), intent(in out) :: val(meqn, nrow, ncol)
    real(kind=8), intent(in out) :: aux(naux, nrow, ncol)
    
    ! Local storage
    integer :: i, j, ibeg, jbeg, nxl, nxr, nyb, nyt
    real(kind=8) :: hxmarg, hymarg, pi
    real(kind=8) :: dxhalf, dyhalf,wl,rho,u,v,p
    real(kind=8) :: pl,ul,vl,rhol,pr,ur,vr,rhor
    real(kind=8) :: xcen,ycen,sloc,sloc_top,rhot,ut,vt,pt

    common /moredata/ sloc

    hxmarg = hx * .01d0
    hymarg = hy * .01d0

    pi = 3.14159265358979d0
    pl = 116.5d0
    ul = 8.25d0*cos(30.d0*pi/180.d0)
    vl = -8.25d0*sin(30.d0*pi/180.d0)
    rhol = 8.d0
    rhor = 1.4d0
    ur = 0.0d0
    vr = 0.0d0
    pr = 1.d0

    ! Use periodic boundary condition specialized code only, if only one 
    ! boundary is periodic we still proceed below
    if (xperdom .and. (yperdom .or. spheredom)) then
        return
    end if

    ! Each check has an initial check to ensure that the boundary is a real
    ! boundary condition and otherwise skips the code.  Otherwise 
    !-------------------------------------------------------
    ! Left boundary:
    !-------------------------------------------------------
    if (xlo_patch < xlower-hxmarg) then
        ! number of grid cells from this patch lying outside physical domain:
        nxl = int((xlower + hxmarg - xlo_patch) / hx)

        select case(mthbc(1))
            case(0) ! User defined boundary condition
                do j = 1, ncol
                    do i=1, nxl
                        val(1, i, j) = rhol
                        val(2, i, j) = rhol*ul
                        val(3, i, j) = rhol*vl
                        val(4, i, j) = pl/.4d0 + 0.5d0*rhol*(ul**2 + vl**2)

                    end do
                end do
               
            case(1) ! Zero-order extrapolation
                do j = 1, ncol
                    do i=1, nxl
                        val(:, i, j) = val(:, nxl + 1, j)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do j = 1, ncol
                    do i=1, nxl
                        val(:, i, j) = val(:, 2 * nxl + 1 - i, j)
                    end do
                end do
                ! negate the normal velocity:
                do j = 1, ncol
                    do i=1, nxl
                        val(2, i, j) = -val(2, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop
        end select
    end if

    !-------------------------------------------------------
    ! Right boundary:
    !-------------------------------------------------------
    if (xhi_patch > xupper+hxmarg) then

        ! number of grid cells lying outside physical domain:
        nxr = int((xhi_patch - xupper + hxmarg) / hx)
        ibeg = max(nrow - nxr + 1, 1)

        select case(mthbc(2))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."
            case(1) ! Zero-order extrapolation
                do i = ibeg, nrow
                    do j = 1, ncol
                        val(:, i, j) = val(:, ibeg - 1, j)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do i=ibeg, nrow
                    do j = 1, ncol
                        val(:, i, j) = val(:, 2 * ibeg - 1 - i, j)
                    end do
                end do
                ! negate the normal velocity:
                do i = ibeg, nrow
                    do j = 1, ncol
                        val(2, i, j) = -val(2, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

    !-------------------------------------------------------
    ! Bottom boundary:
    !-------------------------------------------------------
    if (ylo_patch < ylower - hymarg) then

        ! number of grid cells lying outside physical domain:
        nyb = int((ylower + hymarg - ylo_patch) / hy)

        select case(mthbc(3))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."
            
            case(1) ! Zero-order extrapolation
                do j = 1, nyb
                    do i = 1, nrow
                        val(:,i,j) = val(:, i, nyb + 1)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                ! actually only wall after initial shock location
                ! otherwise flow through (extrap)
                do j = 1, nyb
                    do i = 1, nrow
                        val(:,i,j) = val(:, i, 2 * nyb + 1 - j)
                        xcen = xlo_patch  + (dfloat(i)-.5d0)* hx
                        if (xcen .ge. sloc) val(3,i,j) = -val(3,i,j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

    !-------------------------------------------------------
    ! Top boundary:
    !-------------------------------------------------------
    if (yhi_patch > yupper + hymarg) then

        ! number of grid cells lying outside physical domain:
        nyt = int((yhi_patch - yupper + hymarg) / hy)
        jbeg = max(ncol - nyt + 1, 1)

        select case(mthbc(4))
            case(0) ! User defined boundary condition
              do j = jbeg,ncol
              do i = 1, nrow
                ycen = ylo_patch  + (dfloat(j)-.5d0)* hy
                xcen = xlo_patch  + (dfloat(i)-.5d0)* hx
                sloc_top = sloc + (ycen+20.d0*time)/sqrt(3.d0)
!               pointwise approach
                if (xcen < sloc_top) then
                    rhot = 8.d0
                    ut = 8.25d0*cos(pi/6.d0)
                    vt = -8.25d0*sin(pi/6.d0)
                    pt =  116.5d0
                else
                    rhot = 1.4d0
                    ut = 0.d0
                    vt = 0.d0
                    pt = 1.d0
                endif
                val(1,i,j) = rhot
                val(2,i,j) = rhot*ut
                val(3,i,j) = rhot*vt
                val(4,i,j) = pt/.4d0+0.5d0*rhot*(ut**2+vt**2)

!               cell average over  exact shock
!               dxhalf = 0.5d0*hx
!               dyhalf = 0.5d0*hy
!               call cellave(xcen-dxhalf,ycen-dyhalf,hx,hy,wl,time)
!               rho = (1.d0-wl)*rhor + wl*rhol
!               u   = (1.d0-wl)*ur + wl*ul
!               v   = (1.d0-wl)*vr + wl*vl
!               p   = (1.d0-wl)*pr + wl*pl

!               val(1,i,j) = rho
!               val(2,i,j) = rho*u
!               val(3,i,j) = rho*v
!               val(4,i,j) = p/.4d0+0.5d0*rho*(u**2+v**2)
              end do
              end do


            case(1) ! Zero-order extrapolation
                do j = jbeg, ncol
                    do i = 1, nrow
                        val(:, i, j) = val(:, i, jbeg - 1)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do j = jbeg, ncol 
                    do i = 1, nrow
                        val(:, i, j) = val(:, i, 2 * jbeg - 1 - j)
                    end do
                end do
                ! negate the normal velocity:
                do j = jbeg, ncol
                    do i = 1, nrow
                        val(3, i, j) = -val(3, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

end subroutine bc2amr
