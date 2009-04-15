      subroutine move_parcel(retpt,curpt,uwind,vwind,wwind,
     &                       delt,repeat,error)

      use default_module

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


      real retpt(4), curpt(4), uwind, vwind, wwind, delt, repeat
      logical error

!  local variables
      integer flag(4)
      real realij(4), sigma, sigma_new, zs

!  external functions
      real extrap ! extrap is not part of the transform package
      data flag/1,1,1,1/
!
!
!     assumed use of u,u,omega and plotting in km,km,pressure
!     so wwind in Pa/s and Z is Pa However delt is minutes not seconds
!  step forward in transformed space
      retpt(4) = curpt(4) + delt*60.
      retpt(3) =  curpt(3) + wwind*delt*60.
      retpt(2)  = curpt(2) + vwind*delt*60.
      retpt(1)  = curpt(1) + uwind*delt*60.
c      write(6,*)curpt(1),uwind,delt,curpt(3),wwind
!  Get topography height at new point

      return
      end
