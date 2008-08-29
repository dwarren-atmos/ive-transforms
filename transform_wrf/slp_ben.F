
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   slp_ben - function that calculates the sea-level pressure using
!             the Benjamin technique of finding the temperature 2 km
!             above ground and interpolating to sea-level.
!
!   slp = slp_ben(tmpk, pres, qvapor, hght)
!
!      tmpk - temperature array
!      pres - pressure of lowest model level
!    qvapor - water mixing ratio array
!      hght - geopotential height array
!       slp - output sea-level pressure array
!
!     created April 2005 Ryan Torn, U. Washington
!     adapted for ive AUG 2006 AReinecke
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(ive_ptr) function slp_ben(
     >                  ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names)
      use default_module
c
      implicit none
c
      integer ndims,dims(4), time1,time2
      real stag(4),min(4),max(4),miss
      character*(*) fld1,dunits,disunits,dim_names(4)
c      
      integer strbeg,strend
      integer ix,iy,iz,it
      real ignore
      logical new
      integer(ive_ptr) pres, tmpk, qvap
c
c Get wind fields
c
      pres = 0
      pres = getvar('PPP',ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names,new)
      ix=dims(1);iy=dims(2);iz=dims(3);it=dims(4)

      tmpk = 0
      tmpk = getvar('THE',ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names,new)

      if (tmpk.eq.0.or.pres.eq.0) then
        write(6,*) 'Error:  cannot access variables'
        return
      end if
      
c      qvap = 0
c      qvap = getvar('qvap',ndims,dims,stag,min,max,miss,dunits,
c     >                           disunits,dim_names,new)

      slp_ben = 0
      slp_ben = getmem(ix*iy*it)
      if (time_diff.eq.0) then
        write(6,*) 'Error:  cannot allocate space for slp_ben.'
        return
      end if

      call calc_slp_ben(%val(slp_ben),%val(pres),%val(tmpk),%val(qvap)
     &                  ,%val(topo),%val(sigm),ztop
     &                  ,ix,iy,iz,it,miss)

      dims(1)=ix;dims(2)=iy;dims(3)=1;dims(4)=it
      dim_names(1)='nx';dim_names(2)='ny';
      dim_names(3)='one';dim_names(4)='nt';
      min(3)=0 ; max(3)=0
      dunits='hPa'; disunits='hPa';

      return
      end
c
      subroutine calc_slp_ben(slp_ben,pres,tmpk,qvap,topo,sigm
     &               ,ztop,ix,iy,iz,it,t1,t2,miss)

      integer, intent(in) ix,iy,iz,it
      real, intent(in) :: miss, ztop, sigm(iz)
     &    , pers(ix,iy,iz,it), tmpk(ix,iy,iz,it), qvap(ix,iy,iz,it)

      real, intent(inout) :: slp_ben(ix,iy,it)

      integer :: i,j,k

      do i=1,ix ; do j=1,iy

      enddo ; enddo

      return
      end

