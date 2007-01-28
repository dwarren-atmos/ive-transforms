      subroutine move_parcel(retpt,curpt,uwind,vwind,wwind,
     &                       delt,repeat,error)

      implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This is a single forward time step from curpt to retpt
!  Calculation is carried out in transformed space (wwind should be omega)
!
!  real retpt 	physical point after time step
!  real curpt	initial physical point
!  real uwind	x-component of the wind
!  real vwind	x-component of the wind
!  real wwind	x-component of the wind
!  real delt	time step
!  real repeat	if this is the first time this calculation has been done
!  logical error	
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "default.icl"

      real retpt(4), curpt(4), uwind, vwind, wwind, delt, repeat
      logical error

!  local variables
      integer flag(4)
      real realij(4), sigma, sigma_new, zs

!  external functions
      real extrap ! extrap is not part of the transform package
      data flag/1,1,1,1/
!
!  get topography height at curpt
      call phys_2_index(curpt,realij,flag,4,1)
      zs = extrap(%val(zbot_p),nx+1,ny+1,1,1,
     &      realij(1),realij(2),1.,1.,0.)

!  transform physical height to transformed height
      sigma = ztop/(ztop - zs) * ( curpt(3) - zs )

!  step forward in transformed space
      retpt(4) = curpt(4) + delt
      sigma_new = sigma + wwind*delt
      retpt(2)  = curpt(2) + vwind*delt
      retpt(1)  = curpt(1) + uwind*delt

!  Get topography height at new point
      call phys_2_index(retpt,realij,flag,4,1)
      zs = extrap(%val(zbot_p),nx+1,ny+1,1,1,
     &              realij(1),realij(2),1.,1.,0.)

! transfrom back to physical space
      retpt(3) = zs + (ztop - zs)/ztop * sigma_new

      return
      end
