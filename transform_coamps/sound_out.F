#include "pointer.icl"
c
c converts theta to temperature
c
      subroutine sound_out(ipt,jpt,z1,z2,rotang,timelev,otype,
     >                  ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names)
      implicit none
#include "default.icl"
c
      integer, intent(in) :: otype, timelev
      real, intent(in) :: z1, z2, ipt, jpt, rotang
      
      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)
c      
      integer, parameter :: enkf_bogus=1, coamps_ideal=2
      integer :: ix,iy,iz,it,i,memdim,memdimU,memdimV
      real :: r_val, vartmp
      logical :: new
      ive_ptr :: ppp, the, temp, u, v, uin, vin, thein, pppin,utmp,vtmp

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
      thein = 0
      thein = getvar('THE',ndims,dims,stag,min,max,miss,dunits,
     &                           disunits,dim_names,new)

! 2) Get Pressure field
      pppin = 0
      pppin = getvar('PPP',ndims,dims,stag,min,max,miss,dunits,
     &                           disunits,dim_names,new)

      if (uin.eq.0 .or. vin.eq.0 .or. thein.eq.0 .or. pppin.eq.0) then 
        write(6,*) 'Cannot access one or more variables'
        return
      end if

      ix=dims(1);iy=dims(2);iz=dims(3);it=dims(4)
      memdim = ix*iy*iz
      memdimU = (ix-1)*iy*iz
      memdimV = ix*(iy-1)*iz

!   get space for new variables
      utmp = 0; utmp = getmem(memdimU)
      vtmp = 0; vtmp = getmem(memdimV)

      u = 0; u = getmem(memdim)
      v = 0; v = getmem(memdim)
      the = 0; the = getmem(memdim)
      ppp = 0; ppp = getmem(memdim)
!    check to make sure we got all the variables
      if (u.eq.0 .or. v.eq.0 .or. ppp.eq.0 .or. the.eq.0 .or. utmp.eq.0
     &      .or. vtmp.eq.0) then 
        write(6,*) 'Cannot get space for new variables'
        return
      end if

! 5) get desired time slice of data
      do i=1,memdim
        vartmp = r_val(%val(thein),(timelev-1)*memdim+i)
        call s_val(%val(the),i,vartmp)
        vartmp = r_val(%val(pppin),(timelev-1)*memdim+i)
        call s_val(%val(ppp),i,vartmp)
      end do
     
      do i=1,memdimU
        vartmp = r_val(%val(uin),(timelev-1)*memdimU+i)
        call s_val(%val(utmp),i,vartmp)
      end do

      do i=1,memdimV
        vartmp = r_val(%val(vin),(timelev-1)*memdimV+i)
        call s_val(%val(vtmp),i,vartmp)
      end do

! 6) interpolate u and v points to mass thermo points

      do i=1,memdim
        call s_val(%val(u),i,1.0)
        call s_val(%val(v),i,1.0)
      end do
      call var_at_t(%val(u),%val(utmp),1,0,0,ix-1,iy-1,iz,1)
      call var_at_t(%val(v),%val(vtmp),0,1,0,ix-1,iy-1,iz,1)

      select case (otype)
      case (enkf_bogus)
        call wsnd_enkf(%val(u),%val(v),%val(the),%val(ppp),
     &                       %val(zbot_p),%val(sgz),ztop,
     &                       ipt,jpt,z1,z2,timelev,ix,iy,iz)
      case (coamps_ideal)
        call wsnd_coamps_ideal(
     &               %val(u),%val(v),%val(the),%val(ppp),
     &               %val(zbot_p),%val(sgz),ztop,rotang,
     &               ipt,jpt,z1,ztop,timelev,ix,iy,iz)
      end select

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  wsnd_enkf -- writes sounding that can be ingested by coamps's enkf 
c               data assimilation scheme for bogus obs.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wsnd_enkf(uin,vin,tin,pin,topo,sigma,ztop,
     &                     ipt,jpt,z1,z2,timelev,ix,iy,iz)

      integer, intent(in) :: ix,iy,iz,timelev
      real, intent(in) :: 
     &     uin(ix,iy,iz,1),vin(ix,iy,iz,1),tin(ix,iy,iz,1),
     &     pin(ix,iy,iz,1),topo(ix,iy),sigma(iz),ztop,ipt,jpt,z1,z2

      integer, parameter :: npts=1, ndims=1, flag=1
      integer :: i, j, k, n, ii, jj, kk
      real :: uout(iz), vout(iz), pout(iz), tout(iz), zhgt, zbot, zfrac
      real :: zout(iz), horiz_interp_field, xx, yy, lon, lat
      logical :: fexists
      character(len=24), parameter :: fout='./observations_bogus.dat'

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
         uout(kk) = horiz_interp_field(uin(:,:,k,1),ipt,jpt,ix,iy)
         vout(kk) = horiz_interp_field(vin(:,:,k,1),ipt,jpt,ix,iy)
         tout(kk) = horiz_interp_field(tin(:,:,k,1),ipt,jpt,ix,iy)
         pout(kk) = horiz_interp_field(pin(:,:,k,1),ipt,jpt,ix,iy)
       end if
      end do

      inquire(file=trim(fout),exist=fexists)

      if(fexists) then
        open(unit=10,file=trim(fout),status='old')
      else
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  wsnd_coamps_ideal -- writes a sounding that coamps can read in when
c                       doing ideal simulations.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wsnd_coamps_ideal(uin,vin,tin,pin,topo,sigma,ztop,
     &                     rotang,ipt,jpt,z1,z2,timelev,ix,iy,iz)

      integer, intent(in) :: ix,iy,iz,timelev
      real, intent(in) :: 
     &     uin(ix,iy,iz,1),vin(ix,iy,iz,1),tin(ix,iy,iz,1),
     &     pin(ix,iy,iz,1),topo(ix,iy),sigma(iz),
     &     ztop,ipt,jpt,z1,z2,rotang

      integer, parameter :: npts=1, ndims=1, flag=1
      integer :: i, j, k, n, ii, jj, kk
      real :: uout(iz), vout(iz), pout(iz), tout(iz), qout(iz),
     &        zhgt, zbot, zfrac, zll, pi, utmp, vtmp, Arot, Brot
      real :: zout(iz), horiz_interp_field, xx, yy, lon, lat
      character(len=80) fout
      logical :: fexists
      character(len=2) timeout

      if(timelev.lt.10) then
        write(timeout,'(A1,I1)') '0',timelev
      else
        write(timeout,'(I2)') timelev
      end if
      write(fout,'(A14,A2,A3)') './ideal_sound.',timeout,'.in'
      print *,'out file = ',trim(fout)

      pi = 4.0*atan2(1.0,1.0)
      ! initialize output data
      uout(:) = 0.0 ; vout(:) = 0.0 ; zout(:) = 0.0 
      pout(:) = 0.0 ; tout(:) = 0.0 ; qout(:) = 0.0

      ii = int(floor(ipt)); jj = int(floor(jpt));

      call index_2_phys(xx, ipt, flag, ndims, npts)
      call index_2_phys(yy, jpt, flag, ndims, npts)
      call phys_2_lonlat(xx, yy, lon, lat, npts)

      zbot = horiz_interp_field(topo,ipt,jpt,ix,iy)
      zfrac = (ztop - zbot)/ztop

      kk=0
      Arot = cos(rotang*pi/180.) ; Brot = sin(rotang*pi/180.)
      do k=1,iz
       zhgt = zbot + zfrac*sigma(k)
       if(zhgt.ge.z1.and.zhgt.le.z2) then
         kk=kk+1
         zout(kk) = zhgt
         utmp     = horiz_interp_field(uin(:,:,k,1),ipt,jpt,ix,iy)
         vtmp     = horiz_interp_field(vin(:,:,k,1),ipt,jpt,ix,iy)
         tout(kk) = horiz_interp_field(tin(:,:,k,1),ipt,jpt,ix,iy)
         pout(kk) = horiz_interp_field(pin(:,:,k,1),ipt,jpt,ix,iy)
         !qout(kk) = horiz_interp_field(qin(:,:,k,1),ipt,jpt,ix,iy)
         uout(kk) = Arot*utmp + Brot*vtmp
         vout(kk) =-Brot*utmp + Arot*vtmp
       end if
      end do

      ! lowest z level, we will remove this so sounding starts at z=0
      zll = zout(1)
      ! for 2d simulations V=0
      vout(:) = 0.0

      inquire(file=trim(fout),exist=fexists)

      if(fexists) then
        open(unit=10,file=trim(fout),status='old')
      else
        open(unit=10,file=trim(fout),status='new')
      end if

      write(10,'(6A10)') 'HEIGHT','PRES','THETA','U','V','QMIX'
      write(6,'(6A10)') 'HEIGHT','PRES','THETA','U','V','QMIX'

      do k=1,kk
       write(10,'(6F10.2)') 
     &       zout(k)-zll,pout(k),tout(k),uout(k),vout(k),qout(k)
       write(6,'(6F10.2)') 
     &       zout(k)-zll,pout(k),tout(k),uout(k),vout(k),qout(k)
      end do

      close(10)
     
      return
      end
