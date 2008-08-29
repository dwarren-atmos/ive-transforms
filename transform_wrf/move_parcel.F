      subroutine move_parcel(retpt,curpt,uwind,vwind,wwind,
     &                       delt,repeat,error)

      use DEFAULT_MODULE

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
      integer flag(4), i
      real realij(4), sigma, sigma_new, zs, mupt, hgtval, hgtval2, 
     &     zpart, phys(4), ind(4), ri, spval

!  external functions
      real interp ! interp is not part of the transform package
      data flag/1,1,1,1/
!
      call getrvar('spval', spval, error)
      if(spval .eq. 0.) spval = -999999.999
!  get topography height at curpt
      call phys_2_index(curpt,realij,flag,4,1)
      ind=realij
!     interp is not part of the transform package
      hgtval=interp(%val(phtot_ptr),nx,ny,nz,ntime,
     &     realij(1),realij(2),realij(3),realij(4),spval)
      write(6,*)'hgtval',hgtval
!  step forward in transformed space
      retpt(1)  = curpt(1) + uwind*delt
      retpt(2)  = curpt(2) + vwind*delt
      retpt(4) = curpt(4) + delt
      if(hgtval .eq. spval) then
         retpt(3)=-1
         write(6,*)'left the area'
         return
      endif
      write(6,*)'wwind*delt',wwind*delt,wwind, delt, hgtval + wwind*delt
      sigma_new = hgtval + wwind*delt
!     find sigma_new in the grid
      retpt(3)=curpt(3)
      if(wwind*delt .eq. 0.) then
         return
      endif
!     get the of the index points
      retpt(3)=sigma_new
      write(6,*)'sigma_new',sigma_new
      call hgt_2_index(retpt,realij,4,spval)
!     convert back to pressure
      do i=1,4
         if(abs(ind(i)-realij(i)) .gt. 1.) then
            write(6,*)'diff at ',i,'is ',abs(ind(i)-realij(i))
         endif
      enddo
      call index_2_phys(retpt,realij,flag,4,1)
!      write(6,*)'(',curpt(1),curpt(2),curpt(3),curpt(4),')',
!     &     '->(',retpt(1),retpt(2),retpt(3),retpt(4),')'
!      write(6,*)curpt
!      write(6,*)retpt
!      call sleep(1)
      return
      end
