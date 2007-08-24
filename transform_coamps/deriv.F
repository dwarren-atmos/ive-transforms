      integer(ive_ptr) function deriv(varnm,dir,ifac,ndims,dims,stag,
     &                        min,max,miss,dunits,disunits,dim_names)

      USE DEFAULT_MODULE

      integer ndims,dims(4),ifac
      real stag(4),min(4),max(4),miss
      character(len=*) varnm,dir,dunits,disunits,dim_names(4)
      
      integer idir,iu,iv,iw,it,strbeg,strend,k,nt
      real r_val
      character*(80) tmp
      character*(2) dimu
      integer(ive_ptr) :: var
      logical new,err,requ

      deriv = 0
      nt = ntime

      select case (trim(dir))
        case('X')
          idir = 1
        case('Y')
          idir = 2
        case('Z')
          idir = 3
        case('T')
          idir = 4
        case default
          print *,'Error:  do not know how to compute derivative ',
     &          'for direction ',trim(dir)
      end select

      if ((idir.eq.2).and.(twod.eq.1)) then
        write(6,*) 'Error:  cannot take y derivative for 2d run.'
        return
      end if

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     &                          disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot access variable ',trim(varnm)
        return
      end if

      if (dims(idir).eq.1) then
        print *,'Cannot take derivative of trivial dimension ',trim(dir)
        return
      end if

      deriv = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (deriv.eq.0) then
        print *, 'Error:  cannot allocate space for deriv.'
        return
      end if

      call stagering_calc(iu,iv,iw,it,dims,plmin,min,dx,dy,twod)

      if (ifac.eq.1) then

        if (it.eq.1) then
          call diff1_t(%val(var),%val(deriv),%val(time),nt)
        else if (iw.eq.1) then
    !      call diff1_wpts(%val(var),%val(deriv),nx,ny,nz,nt,
    ! &                    idir,iu,iv,iw)
          call diff1_wpts(%val(var),%val(deriv),%val(dqdx_xy),
     &              %val(dqdx_zt),%val(dqdy_xy),%val(dqdy_zt),
     &              %val(dqdz),dx,dy,%val(time),nx,ny,nz,nt,idir,
     &              iu,iv,iw)
        else
    !      call diff1(%val(var),%val(deriv),nx,ny,nz,nt,idir,iu,iv,iw)
          call diff1(%val(var),%val(deriv),%val(dqdx_xy),
     &             %val(dqdx_zw),%val(dqdy_xy),%val(dqdy_zw),
     &             %val(dqdz),dx,dy,%val(time),nx,ny,nz,nt,idir,
     &             iu,iv,iw)
        end if

        select case (idir)
        case (4)
          dimu = '/s'
        case default
          dimu = '/m'
        end select

      else

        write(6,*) 'Error:  second derivative not yet implemented.'  
        call freemem(deriv)
        return

      end if

      tmp = dunits(strbeg(dunits):strend(dunits))//dimu
      call units_smpfy(tmp,dunits,err)

      tmp = disunits(strbeg(disunits):strend(disunits))//dimu
      call units_smpfy(tmp,disunits,err)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !subroutine diff1_wpts(var,deriv,ix,iy,iz,nt,idir,iu,iv,iw)
      subroutine diff1_wpts(var,deriv,dqdx_xy,dqdx_zt,dqdy_xy,
     &                      dqdy_zt,dqdz,dx,dy,time,ix,iy,iz,nt,idir,
     &                      iu,iv,iw)

      use misc_module
      !use default_module

      integer ix,iy,iz,nt,idir,iu,iv,iw
      real var(ix+1-iu,iy+1-iv,iz+1,nt)
	  real deriv(ix+1-iu,iy+1-iv,iz+1,nt)
      real dqdx_xy(ix+1,iy+1,3),dqdx_zt(iz+1)
      real dqdy_xy(ix+1,iy+1,3),dqdy_zt(iz+1)
      real dqdz(ix+1,iy+1,3),dx,dy,time(nt)

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real rdx,rdy,da,db,ta,tb,dta,dtb      

      rdx = 1./dx
      rdy = 1./dy
      if (idir.eq.1) then ! x-derivative

        do n=1,nt
*
* west and east boundaries
*
        do m=1,4
          if (m.eq.1) i = 1
          if (m.eq.2) i = 2
          if (m.eq.3) i = ix-iu
          if (m.eq.4) i = ix+1-iu
          ip1 = min0(i+1,ix+1-iu)
          im1 = max0(i-1,1)
        do j=1,iy+1-iv

          k = 1 
          deriv(i,j,k,n) = 0.5*rdx
     &       *( xrbnd(i)*(var(ip1,j,k,n)-var(i,j,k,n))
     &        + xlbnd(i)*(var(i,j,k,n)-var(im1,j,k,n)) )
     &                   + dqdx_xy(i,j,1)
     &      * rdzc(k)*dqdx_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))

          k = iz+iw 
          deriv(i,j,k,n) = 0.5*rdx
     &       *( xrbnd(i)*(var(ip1,j,k,n)-var(i,j,k,n))
     &        + xlbnd(i)*(var(i,j,k,n)-var(im1,j,k,n)) )
     &                   + dqdx_xy(i,j,1)
     &      * rdzc(k-1)*dqdx_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n))

          do k=2,iz-1+iw
            deriv(i,j,k,n) = 0.5*rdx
     &         *( xrbnd(i)*(var(ip1,j,k,n)-var(i,j,k,n))
     &          + xlbnd(i)*(var(i,j,k,n)-var(im1,j,k,n)) )
     &                     + 0.5*dqdx_xy(i,j,1)
     &        *( rdzc(k)*dqdx_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     &       + rdzc(k-1)*dqdx_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n)) )
          enddo

        enddo
        enddo
*
* interior
*
        do j=1,iy+1-iv
        do i=3,ix-1-iu

          k = 1 
          deriv(i,j,k,n) = 4./6.*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     &                 - 1./12.*rdx*(var(i+2,j,k,n)-var(i-2,j,k,n))
     &                   + dqdx_xy(i,j,1)
     &      * rdzc(k)*dqdx_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))

          k = iz+iw 
          deriv(i,j,k,n) = 4./6.*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     &                 - 1./12.*rdx*(var(i+2,j,k,n)-var(i-2,j,k,n))
     &                   + dqdx_xy(i,j,1)
     &      * rdzc(k-1)*dqdx_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n))

          do k=2,iz-1+iw
            deriv(i,j,k,n) = 4./6.*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     &                   - 1./12.*rdx*(var(i+2,j,k,n)-var(i-2,j,k,n))
     &                     + 0.5*dqdx_xy(i,j,1)
     &        *( rdzc(k)*dqdx_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     &       + rdzc(k-1)*dqdx_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n)) )
          enddo

        enddo ; enddo

        enddo ! end time loop

      else if (idir.eq.2) then ! y-derivative

        do n=1,nt
*
* north and south boundaries
*
        do m=1,4
          if (m.eq.1) j = 1
          if (m.eq.2) j = 2
          if (m.eq.3) j = iy-iv
          if (m.eq.4) j = iy+1-iv
          jp1 = min0(j+1,iy+1-iv)
          jm1 = max0(j-1,1)
        do i=1,ix+1-iu
         
          k = 1 
          deriv(i,j,k,n) = 0.5*rdy
     &       *( yrbnd(j)*(var(i,jp1,k,n)-var(i,j,k,n))
     &        + ylbnd(j)*(var(i,j,k,n)-var(i,jm1,k,n)) )
     &                   + dqdy_xy(i,j,1)
     &      * rdzc(k)*dqdy_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))

          k = iz+iw 
          deriv(i,j,k,n) = 0.5*rdy
     &       *( yrbnd(j)*(var(i,jp1,k,n)-var(i,j,k,n))
     &        + ylbnd(j)*(var(i,j,k,n)-var(i,jm1,k,n)) )
     &                   + dqdy_xy(i,j,1)
     &      * rdzc(k-1)*dqdy_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n))

          do k=2,iz-1+iw
            deriv(i,j,k,n) = 0.5*rdy
     &       *( yrbnd(j)*(var(i,jp1,k,n)-var(i,j,k,n))
     &        + ylbnd(j)*(var(i,j,k,n)-var(i,jm1,k,n)) )
     &                     + 0.5*dqdy_xy(i,j,1)
     &        *( rdzc(k)*dqdy_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     &       + rdzc(k-1)*dqdy_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n)) )
          enddo
     
        enddo
        enddo
*
* interior
*
        do j=3,iy-1-iv
        do i=1,ix+1-iu

          k = 1 
          deriv(i,j,k,n) = 4./6.*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     &                 - 1./12.*rdy*(var(i,j+2,k,n)-var(i,j-2,k,n))
     &                   + dqdy_xy(i,j,1)
     &      * rdzc(k)*dqdy_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))

          k = iz+iw 
          deriv(i,j,k,n) = 4./6.*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     &                 - 1./12.*rdy*(var(i,j+2,k,n)-var(i,j-2,k,n))
     &                   + dqdy_xy(i,j,1)
     &      * rdzc(k-1)*dqdy_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n))

          do k=2,iz-i+iw
            deriv(i,j,k,n) = 4./6.*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     &                   - 1./12.*rdy*(var(i,j+2,k,n)-var(i,j-2,k,n))
     &                     + 0.5*dqdy_xy(i,j,1)
     &        *( rdzc(k)*dqdy_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     &       + rdzc(k-1)*dqdy_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n)) )

          enddo 
        enddo ; enddo

        enddo ! end time loop

      else if (idir.eq.3) then ! z-derivative
        
        do n=1,nt
        do j=1,iy+1-iv
        do i=1,ix+1-iu

          k = 1
          deriv(i,j,k,n) = dqdz(i,j,1)
     &               * rdzc(k)*(var(i,j,k+1,n)-var(i,j,k,n))

          k = iz+iw
          deriv(i,j,k,n) = dqdz(i,j,1)
     &             * rdzc(k-1)*(var(i,j,k,n)-var(i,j,k-1,n))

          do k=2,iz-1+iw
            deriv(i,j,k,n) = 0.5*dqdz(i,j,1)
     &         *( rdzc(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     &        + rdzc(k-1)*(var(i,j,k,n)-var(i,j,k-1,n)) )      
          enddo

        enddo ; enddo ; enddo

      else if (idir.eq.4) then ! t-derivative
        do k=1,iz+iw
        do j=1,iy+1-iv
        do i=1,ix+1-iu 

          deriv(i,j,k,1) = (var(i,j,k,2)-var(i,j,k,1))
     &                        / (time(2)-time(1))
          deriv(i,j,k,nt) = (var(i,j,k,nt)-var(i,j,k,nt-1))
     &                        / (time(nt)-time(nt-1))

          do n=2,nt-1
            dta = time(n+1)-time(n)
            dtb = time(n)-time(n-1)
            deriv(i,j,k,n) = dtb/(dta+dtb)
     &               * (var(i,j,k,n+1)-var(i,j,k,n))/dta
     &                     + dta/(dta+dtb)
     &               * (var(i,j,k,n)-var(i,j,k,n-1))/dtb
          enddo

        enddo ; enddo ; enddo
      end if

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine diff1_t(var,deriv,time,nt)

      integer nt
      real var(nt),deriv(nt),time(nt)
  
      integer n 
      real dta,dtb     

      deriv(1) = (var(2)-var(1))/(time(2)-time(1))
      deriv(nt) = (var(nt)-var(nt-1))/(time(nt)-time(nt-1))

      do n=2,nt-1
        dta = time(n+1)-time(n)
        dtb = time(n)-time(n-1)
        deriv(n) = dtb/(dta+dtb)*(var(n+1)-var(n))/dta
     &                   + dta/(dta+dtb)*(var(n)-var(n-1))/dtb
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !subroutine diff1(var,deriv,ix,iy,iz,nt,idir,iu,iv,iw)
      subroutine diff1(var,deriv,dqdx_xy,dqdx_zw,dqdy_xy,
     &                dqdy_zw,dqdz,dx,dy,time,ix,iy,iz,nt,idir,iu,iv,iw)

      use misc_module
      !use default_module

      integer ix,iy,iz,nt,idir,iu,iv,iw
      real var(ix+1-iu,iy+1-iv,iz+iw,nt)
	  real deriv(ix+1-iu,iy+1-iv,iz+iw,nt)
      real dqdx_xy(ix+1,iy+1,3),dqdx_zw(iz+1)
      real dqdy_xy(ix+1,iy+1,3),dqdy_zw(iz+1)
      real dqdz(ix+1,iy+1,3),dx,dy,time(nt)

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1,idqdx,ii,jj
      real rdx,rdy,da,db,ta,tb,dta,dtb      

      rdx = 1./dx
      rdy = 1./dy
c
c  Determine if we are at u, v, or mass points.
c
      if (iu.eq.0.and.iv.eq.0) idqdx=1
      if (iu.eq.1.and.iv.eq.0) idqdx=2
      if (iu.eq.0.and.iv.eq.1) idqdx=3
c
c x derivative
c
      if (idir.eq.1) then

        do n=1,nt
c
c east and west boundaries
c
        do k=1,iz+iw
          kp1 = min0(k+1,iz+iw)
          km1 = max0(k-1,1)
          do m=1,4
            if (m.eq.1) then
              i = 1 ; ii = 1
            elseif (m.eq.2) then
              i = 2 ; ii = 2
            elseif (m.eq.3) then 
              i = ix-iu ; ii = ix
            elseif (m.eq.4) then
              i = ix+1-iu ; ii = ix+1
            end if
            ip1 = min0(i+1,ix+1-iu) ; im1 = max0(i-1,1)
            do j=1,iy+1-iv
          
               deriv(i,j,k,n) = 0.5*rdx
     &         *( xrbnd(ii)*(var(ip1,j,k,n)-var(i,j,k,n))
     &         + xlbnd(ii)*(var(i,j,k,n)-var(im1,j,k,n)) )
     &                   + 0.5*dqdx_xy(i,j,idqdx)
     &         *( rdze(kp1)*dqdx_zw(kp1)*zrbnd(k)
     &                       * (var(i,j,kp1,n)-var(i,j,k,n))
     &         + rdze(k)*dqdx_zw(k)*zlbnd(k)
     &                       * (var(i,j,k,n)-var(i,j,km1,n)))
        
        enddo ; enddo ; enddo
c
c interior
c
        do k=1,iz+iw      
          kp1 = min0(k+1,iz+iw)
          km1 = max0(k-1,1)
        do j=1,iy+1-iv
        do i=3,ix-1-iu

          deriv(i,j,k,n) = 4./6.*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     &                 - 1./12.*rdx*(var(i+2,j,k,n)-var(i-2,j,k,n))
     &                   + 0.5*dqdx_xy(i,j,idqdx)
     &          *( rdze(kp1)*dqdx_zw(kp1)*zrbnd(k)
     &                       * (var(i,j,kp1,n)-var(i,j,k,n))
     &             + rdze(k)*dqdx_zw(k)*zlbnd(k)
     &                       * (var(i,j,k,n)-var(i,j,km1,n)) )

        enddo ; enddo ; enddo
        enddo ! end time loop

      else if (idir.eq.2) then ! y-derivative

        do n=1,nt
*
* north and south boundaries
*
        do k=1,iz+iw
          kp1 = min0(k+1,iz+iw)
          km1 = max0(k-1,1)
        do m=1,4
          if (m.eq.1) then
            j = 1 ; jj = 1
          elseif (m.eq.2) then
            j = 2 ; jj = 2
          elseif (m.eq.3) then 
            j = iy-iv ; jj = iy
          elseif (m.eq.4) then
            j = iy+1-iv ;  jj = iy+1
          end if
          jp1 = min0(j+1,iy+1-iv)
          jm1 = max0(j-1,1)
        do i=1,ix+1-iu
          
          deriv(i,j,k,n) = 0.5*rdy
     &       *( yrbnd(jj)*(var(i,jp1,k,n)-var(i,j,k,n))
     &        + ylbnd(jj)*(var(i,j,k,n)-var(i,jm1,k,n)) )
     &                   + 0.5*dqdy_xy(i,j,idqdx)
     &          *( rdze(kp1)*dqdy_zw(kp1)*zrbnd(k)
     &                       * (var(i,j,kp1,n)-var(i,j,k,n))
     &             + rdze(k)*dqdy_zw(k)*zlbnd(k)
     &                       * (var(i,j,k,n)-var(i,j,km1,n)) )
        
        enddo ; enddo ; enddo
*
* interior
*
        do k=1,iz+iw
          kp1 = min0(k+1,iz+iw)
          km1 = max0(k-1,1)
        do j=3,iy-1-iv
        do i=1,ix+1-iu

          deriv(i,j,k,n) = 4./6.*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     &                 - 1./12.*rdy*(var(i,j+2,k,n)-var(i,j-2,k,n))
     &                   + 0.5*dqdy_xy(i,j,idqdx)
     &          *( rdze(kp1)*dqdy_zw(kp1)*zrbnd(k)
     &                       * (var(i,j,kp1,n)-var(i,j,k,n))
     &             + rdze(k)*dqdy_zw(k)*zlbnd(k)
     &                       * (var(i,j,k,n)-var(i,j,km1,n)) )

        enddo ; enddo ; enddo
        enddo ! end time loop

      else if (idir.eq.3) then ! z-derivative
        
        do n=1,nt
        do k=1,iz+iw
          kp1 = min0(k+1,iz+iw)
          km1 = max0(k-1,1)
        do j=1,iy+1-iv
        do i=1,ix+1-iu

          deriv(i,j,k,n) = 0.5*dqdz(i,j,idqdx)
     &      *( rdze(kp1)*zrbnd(k)*(var(i,j,kp1,n)-var(i,j,k,n))
     &         + rdze(k)*zlbnd(k)*(var(i,j,k,n)-var(i,j,km1,n)) )

        enddo ; enddo ; enddo ; enddo

      else if (idir.eq.4) then ! t-derivative

        do k=1,iz+iw
        do j=1,iy+1-iv
        do i=1,ix+1-iu 

          deriv(i,j,k,1) = (var(i,j,k,2)-var(i,j,k,1))
     &                        / (time(2)-time(1))
          deriv(i,j,k,nt) = (var(i,j,k,nt)-var(i,j,k,nt-1))
     &                        / (time(nt)-time(nt-1))

          do n=2,nt-1
            dta = time(n+1)-time(n)
            dtb = time(n)-time(n-1)
            deriv(i,j,k,n) = dtb/(dta+dtb)
     &               * (var(i,j,k,n+1)-var(i,j,k,n))/dta
     &                     + dta/(dta+dtb)
     &               * (var(i,j,k,n)-var(i,j,k,n-1))/dtb
          enddo

        enddo ; enddo ; enddo

      end if

      return
      end

