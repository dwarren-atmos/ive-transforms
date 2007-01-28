#include "pointer.icl"
c
c converts theta to temperature
c
      subroutine sound_out(ipt,jpt,z1,z2,tpt,
     >                  ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names)
      implicit none
#include "default.icl"
c
      real, intent(in) :: z1, z2, ipt, jpt, tpt
      
      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)
c      
      integer :: ix,iy,iz,it,i
      real :: r_val
      logical :: new
      ive_ptr :: ppp, the, temp, u, v, uin, vin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 1) get theta 
! 2) get pressure
! 3) convert theta to temperature
! 4) get u velocity
! 5) get v velocity
! 6) interpolate u and v to mass points
! 7) pass u,v,temp, and pressure into subroutine that outputs a sounding
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 3) Get U field
      uin = 0
      uin = getvar('U',ndims,dims,stag,min,max,miss,dunits,
     &                           disunits,dim_names,new)

! 4) Get V field
      vin = 0
      vin = getvar('V',ndims,dims,stag,min,max,miss,dunits,
     &                           disunits,dim_names,new)

! 1) Get Theta field
      the = 0
      the = getvar('THE',ndims,dims,stag,min,max,miss,dunits,
     &                           disunits,dim_names,new)

! 2) Get Pressure field
      ppp = 0
      ppp = getvar('ppp',ndims,dims,stag,min,max,miss,dunits,
     &                           disunits,dim_names,new)

      if (uin.eq.0 .or. vin.eq.0 .or. the.eq.0 .or. ppp.eq.0) then 
        write(6,*) 'Cannot access one or more variables'
        return
      end if

      ix=dims(1);iy=dims(2);iz=dims(3);it=dims(4)

!   get space for new variables
c      temp = 0; temp = getmem(ix*iy*iz*it)
      u = 0; u = getmem(ix*iy*iz*it)
      v = 0; v = getmem(ix*iy*iz*it)
!    check to make sure we got all the variables
      if (u.eq.0 .or. v.eq.0 ) then !.or. temp.eq.0) then
        write(6,*) 'Cannot get space for new variables'
        return
      end if

! 5) Get Temperature field
c      call calc_the2temp(%val(temp),%val(the),%val(ppp),
c     &                     ix,iy,iz,it,miss)

! 6) interpolate u and v points to mass thermo points

      do i=1,ix*iy*iz*it
        call s_val(%val(u),i,1.0)
        call s_val(%val(v),i,1.0)
      end do
      call var_at_t(%val(u),%val(uin),1,0,0,ix-1,iy-1,iz,it)
      call var_at_t(%val(v),%val(vin),0,1,0,ix-1,iy-1,iz,it)

      call write_sound_out(%val(u),%val(v),%val(the),%val(ppp),
     &                     %val(zbot_p),%val(sgz),ztop,
     &                     ipt,jpt,z1,z2,tpt,ix,iy,iz,it)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  write_sound_out
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_sound_out(uin,vin,tin,pin,
     &                           topo,sigma,ztop,
     &                           ipt,jpt,z1,z2,tpt,ix,iy,iz,it)

      integer, intent(in) :: ix,iy,iz,it,tpt
      real, intent(in) :: 
     &     uin(ix,iy,iz,it),vin(ix,iy,iz,it),tin(ix,iy,iz,it),
     &     pin(ix,iy,iz,it),
     &     topo(ix,iy),sigma(iz),ztop,
     &     ipt,jpt,z1,z2

      integer, parameter :: npts=1, ndims=1, flag=1
      integer :: i, j, k, n, ii, jj, kk
      real :: uout(iz), vout(iz), pout(iz), tout(iz), zhgt, zbot, zfrac
      real :: zout(iz), horiz_interp_field, xx, yy, lon, lat
      character(len=80) fout
      logical :: fexists
      fout='sound_out.sound'

      ii = int(floor(ipt)); jj = int(floor(jpt));

      call index_2_phys(xx, ipt, flag, ndims, npts)
      call index_2_phys(yy, jpt, flag, ndims, npts)
      call phys_2_lonlat(xx, yy, lon, lat, npts)

      zbot = horiz_interp_field(topo,ipt,jpt,ix,iy)
      zfrac = (ztop - zbot)/ztop

      kk=0
      do k=1,iz
        zhgt = zbot + zfrac*sigma(k)
        if(zhgt.ge.z1.and.zhgt.le.z2) then
          kk=kk+1
          zout(kk) = zhgt
          uout(kk) = horiz_interp_field(uin(:,:,k,tpt),ipt,jpt,ix,iy)
          vout(kk) = horiz_interp_field(vin(:,:,k,tpt),ipt,jpt,ix,iy)
          tout(kk) = horiz_interp_field(tin(:,:,k,tpt),ipt,jpt,ix,iy)
          pout(kk) = horiz_interp_field(pin(:,:,k,tpt),ipt,jpt,ix,iy)
        end if
      end do

      inquire(file=trim(fout),exist=fexists)

      if(fexists) then
        print *,fexists,'T'
        open(unit=10,file=trim(fout),status='old')
      else
        print *,fexists,'F'
        open(unit=10,file=trim(fout),status='new')
      end if
      write(10,'(I2)') kk
      write(10,'(7A10)') 'HEIGHT','U','V','THETA','PPP','LAT','LON'
      do kk=1,kk
       write(10,'(F10.2,4F10.2,2F10.3)') zout(kk),uout(kk),vout(kk),
     &                           tout(kk),pout(kk),lat,lon
      end do
      close(10)
     
      return
      end
