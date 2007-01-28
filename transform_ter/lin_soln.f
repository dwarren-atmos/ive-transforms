
      program lin_soln

      include '/usr/local/include/netcdf.inc'

      integer ispc
      parameter (ispc=10000000)

      real a(ispc)

      integer idcdf,isclr,isgz,iwgz,izs_p,izs_u,izs_v,itran,iwrk
      integer ivar,ivar2,ik,il,iwrk_x,iwrk_y,nx,ny,nz
      real dx,dy,dz,dzbot,x0,y0
      character*(80) outfl

      real R,g,rhos
      parameter (R=287.04, g=9.806)

      integer nvar,nsclr
      parameter (nvar=12,nsclr=16)

      character*(80) var(nvar),units(nvar),sclr(nsclr)
      data var / 'u','v','w','p','rho','zeta','xi','eta','u0','v0',
     >           'rho0','p0' /
      data units / 'm/s','m/s','m/s','N','kg/m/m/m','m','m','m',
     >             'm/s','m/s','kg/m/m/m','N' /
      data sclr / 'f_cor','ztop','type','h','a1','a2','b1','b2',
     >            'y_len','xc','yc','N','U0','V0','ps','Ts' /

 886  FORMAT(A80)

      isclr = 1

      iunit = 32

      open (unit=iunit,file='lin_soln.in',status='old')

      read(iunit,*) x0,y0
      read(iunit,*) nx,ny,nz
      read(iunit,*) dx,dy
      read(iunit,*) dz,dzbot
      read(iunit,*) a(isclr+3-1)
      read(iunit,*) a(isclr+5-1),a(isclr+7-1),a(isclr+4-1)
      read(iunit,*) a(isclr+6-1),a(isclr+8-1),a(isclr+9-1)
      read(iunit,*) a(isclr+10-1),a(isclr+11-1)
      read(iunit,*) a(isclr)
      read(iunit,*) a(isclr+12-1)
      read(iunit,*) a(isclr+13-1)
      read(iunit,*) a(isclr+14-1)
      read(iunit,*) a(isclr+15-1)
      read(iunit,*) a(isclr+16-1)
      read(iunit,886) outfl

      close (unit=iunit)

      a(isclr+2-1) = dz*(nz-1)
      rhos = a(isclr+15-1)/R/a(isclr+16-1)

      if ((mod(nx,2).eq.0).or.(mod(ny,2).eq.0)) then
        write(6,*)
        write(6,*) 'Error:  nx-1 and ny-1 must be even.'
        write(6,*)
        stop
      end if  

      isgz = isclr + nsclr
      iwgz = isgz + nz
      izs_p = iwgz + nz
      izs_u = izs_p + 2*nx*ny 
      izs_v = izs_u + 2*nx*ny
      itran = izs_v + 2*nx*ny
      iwrk = itran + 2*nx*ny
      itmp = iwrk + 2*nx*ny
      ik = itmp + nx*ny
      il = ik + nx
      iwrk_x = il + ny
      iwrk_y = iwrk_x + 4*nx + 15
      ivar = iwrk_y + 4*ny + 15
      ivar2 = ivar + nx*ny*nz

      if ((ivar2+nx*ny*nz).gt.ispc) then
        write(6,*)
        write(6,*) 'Not enough space!  Increase parameter ispc.'
        write(6,*)
        stop
      end if

      call get_ter(a(isclr),a(izs_p),a(izs_u),a(izs_v),x0,y0,dx,dy,
     >                       nx,ny,nsclr)

      call stretch_it(a(isgz),a(iwgz),dz,dzbot,nz)

      call set_outfl(idcdf,nx,ny,nz,dx,dy,dz,x0,y0,a(isgz),
     >                a(iwgz),a(isclr),sclr,var,units,outfl,
     >                nsclr,nvar)

      call write_var(idcdf,'zbot_p',nx,ny,1,a(izs_p))
      call write_var(idcdf,'zbot_u',nx,ny,1,a(izs_u))
      call write_var(idcdf,'zbot_v',nx,ny,1,a(izs_v))

      call calc_wavnm(a(ik),dx,nx-1) 
      call calc_wavnm(a(il),dy,ny-1)

      do ij=1,nx*ny
        a(itmp+ij-1) = a(izs_p+ij-1)
      enddo

      call trnsfrm(a(izs_p),a(itran),a(iwrk_x),a(iwrk_y),nx,ny,0,0)

      do k=1,nvar
        if ((var(k)(1:1).ne.'u').and.(var(k)(1:1).ne.'v')) then
          if (var(k)(1:1).eq.'w') then
            igz = iwgz
            iw = 1
          else
            igz = isgz
            iw = 0
          end if
          call calc_fld(a(ivar),a(ivar2),a(izs_p),a(itran),a(iwrk),
     >         a(itmp),a(igz),a(iwrk_x),a(iwrk_y),a(ik),a(il),
     >         a(isclr+1-1),a(isclr+2-1),a(isclr+12-1),a(isclr+13-1),
     >         a(isclr+14-1),a(isclr+15-1),rhos,g,dx,dy,var(k),
     >         nx,ny,nz,0,0,iw)
          call write_var(idcdf,var(k),nx,ny,nz,a(ivar))
          lgth = index(var(k),' ')-1
          if (var(k)(1:lgth).eq.'p') then
            call drag_calc(a(ivar),a(izs_p),a(itran),a(iwrk),a(izs_u),
     >               a(izs_v),a(iwrk_x),a(iwrk_y),a(ik),a(il),
     >               a(isclr+1-1),a(isclr+12-1),a(isclr+13-1),
     >               a(isclr+14-1),rhos,dx,dy,nx,ny,nz)
          end if
        end if
      enddo

      do ij=1,nx*ny
        a(itmp+ij-1) = a(izs_u+ij-1)
      enddo

      call trnsfrm(a(izs_u),a(itran),a(iwrk_x),a(iwrk_y),nx,ny,1,0)

      do k=1,nvar
        if (var(k)(1:1).eq.'u') then
          call calc_fld(a(ivar),a(ivar2),a(izs_u),a(itran),a(iwrk),
     >         a(itmp),a(isgz),a(iwrk_x),a(iwrk_y),a(ik),a(il),
     >         a(isclr+1-1),a(isclr+2-1),a(isclr+12-1),a(isclr+13-1),
     >         a(isclr+14-1),a(isclr+15-1),rhos,g,dx,dy,var(k),
     >         nx,ny,nz,1,0,0)
          call write_var(idcdf,var(k),nx,ny,nz,a(ivar))
        end if
      enddo

      do ij=1,nx*ny
        a(itmp+ij-1) = a(izs_v+ij-1)
      enddo

      call trnsfrm(a(izs_v),a(itran),a(iwrk_x),a(iwrk_y),nx,ny,0,1)

      do k=1,nvar
        if (var(k)(1:1).eq.'v') then
          call calc_fld(a(ivar),a(ivar2),a(izs_v),a(itran),a(iwrk),
     >         a(itmp),a(isgz),a(iwrk_x),a(iwrk_y),a(ik),a(il),
     >         a(isclr+1-1),a(isclr+2-1),a(isclr+12-1),a(isclr+13-1),
     >         a(isclr+14-1),a(isclr+15-1),rhos,g,dx,dy,var(k),
     >         nx,ny,nz,0,1,0)
          call write_var(idcdf,var(k),nx,ny,nz,a(ivar))
        end if
      enddo

      call ncclos(idcdf,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem closing netcdf file.'
        stop
      end if

      end

*----------------------------------------------------------------------

      subroutine calc_wavnm(ak,ds,N)

      integer N
      real ak(N),ds

      integer q
      real pi,L

      pi = 2.*asin(1.)
      L = N*ds
      q = N/2

      do j=1,N
        if (j.le.q+1) then
          ak(j) = 2.*pi*(j-1)/L
        else
          ak(j) = 2.*pi*(j-1-N)/L
        end if
      enddo

      return
      end

*-----------------------------------------------------------------------

      subroutine trnsfrm(var,tran,cwrkx,cwrky,nx,ny,iu,iv)

      integer nx,ny,iu,iv
      real cwrkx(4*nx+15),cwrky(4*ny+15)
      complex var(nx-1+iu,ny-1+iv),tran(ny-1,nx-1)

      call cffti(nx-1,cwrkx)
      call cffti(ny-1,cwrky)

      call r_to_c(var,nx,ny,iu,iv)

      do j=1,ny-1
        call cfftf(nx-1,var(1,j),cwrkx)
        do i=1,nx-1
          tran(j,i) = var(i,j)
        enddo
      enddo

      do i=1,nx-1
        call cfftf(ny-1,tran(1,i),cwrky)
        do j=1,ny-1
          var(i,j) = tran(j,i)
        enddo
      enddo

      return
      end

*-----------------------------------------------------------------------

      subroutine r_to_c(var,nx,ny,iu,iv)

      integer nx,ny,iu,iv
      real var(nx-1+iu,ny-1+iv,2)

      integer ptr,q,numx,numy

      numx = nx-1+iu
      numy = ny-1+iv

      do j=1,numy
      do i=1,numx
        var(i,j,2) = var(i,j,1)
      enddo
      enddo

      ptr = numx*numy + 1
 
      do k=1,2
      do j=1,numy
      do i=1,numx

        q = (k-1)*numx*numy + (j-1)*numx + i
        if (mod(q,2).ne.0) then
          var(i,j,k) = var(ptr,1,1)
          ptr = ptr + 1
        else
          var(i,j,k) = 0.
        end if

      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine calc_fld(fld,fld2,hft,tran,wrk,zs,gz,cwrkx,cwrky,
     >                 ak,al,fcor,ztop,N,U0,V0,ps,rhos,g,dx,dy,var,
     >                 nx,ny,nz,iu,iv,iw)

      integer nx,ny,nz,iu,iv,iw
      real fld(nx-1+iu,ny-1+iv,nz-1+iw)
      real fld2(nx-1+iu,ny-1+iv,nz-1+iw)
      real zs(nx-1+iu,ny-1+iv),ak(nx-1),al(ny-1),gz(nz-1+iw)
      real cwrkx(4*nx+15),cwrky(4*ny+15),fcor,ztop,N,U0,V0,rhos
      real ps,g,dx,dy
      complex hft(nx-1+iu,ny-1+iv),tran(ny-1,nx-1),wrk(nx-1,ny-1)
      character*(*) var

      real x,y,z,m,m2,gam,gam2,omga,fac1,fac2
      logical requ
      complex pfac,exfac

      lgth = index(var,' ') - 1

      if (var(lgth:lgth).eq.'0') then
 
        if (var(1:lgth).eq.'u0') then
          do k=1,nz-1  
          do j=1,ny-1
          do i=1,nx
            fld(i,j,k) = U0
          enddo
          enddo
          enddo
        else if (var(1:lgth).eq.'v0') then
          do k=1,nz-1
          do j=1,ny
          do i=1,nx-1
            fld(i,j,k) = V0
          enddo
          enddo
          enddo
        else if (var(1:lgth).eq.'rho0') then
          do k=1,nz-1
          do j=1,ny-1
          do i=1,nx-1
            z = zs(i,j) + (ztop-zs(i,j))/ztop*gz(k)
            fld(i,j,k) = rhos*(1-N**2/g*z)        
          enddo
          enddo
          enddo
        else if (var(1:lgth).eq.'p0') then
          do k=1,nz-1
          do j=1,ny-1
          do i=1,nx-1
            z = zs(i,j) + (ztop-zs(i,j))/ztop*gz(k)
            x = ((i-1)+0.5)*dx
            y = ((j-1)+0.5)*dy
            fld(i,j,k) = ps + rhos*fcor*V0*x - rhos*fcor*U0*y
     >                   - rhos*g*z + 0.5*rhos*(N**2)*(z**2)
          enddo
          enddo
          enddo
        else 
          write(6,*) 'Error:  field not recognized.'
          stop
        end if

        return

      end if

      do k=1,nz-1+iw

        z = gz(k)

        do i=1,nx-1

          do j=1,ny-1

            omga = ak(i)*U0+al(j)*V0 
            if (requ(abs(omga),abs(fcor)).or.(omga.eq.0.)) then
              m2 = 0
              gam2 = 0
            else
              m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >                / (omga**2-fcor**2)
              gam2 = -m2
            end if
  
            if (requ(abs(omga),abs(fcor)).or.(omga.eq.0.)) then

              pfac = 0.

            else if (m2.gt.0.) then

              m = sqrt(m2)*omga/abs(omga)
              exfac = cmplx(cos(m*z),sin(m*z))/(nx-1)/(ny-1)
              if (var(1:lgth).eq.'zeta') then
                pfac = exfac
              else if (var(1:lgth).eq.'w') then
                pfac = cmplx(0.,omga)*exfac
              else if (var(1:lgth).eq.'u') then
                fac1 = - m*fcor*al(j)/(ak(i)**2+al(j)**2)  
                fac2 = - m*ak(i)*omga/(ak(i)**2+al(j)**2)
                pfac = cmplx(fac1,fac2)*exfac
              else if (var(1:lgth).eq.'v') then
                fac1 = m*fcor*ak(i)/(ak(i)**2+al(j)**2)
                fac2 = - m*al(j)*omga/(ak(i)**2+al(j)**2)
                pfac = cmplx(fac1,fac2)*exfac
              else if (var(1:lgth).eq.'p') then
                fac2 = rhos*(N**2-omga**2)/m
                pfac = cmplx(0.,fac2)*exfac
              else if (var(1:lgth).eq.'rho') then
                fac1 = rhos*N**2/g
                pfac = cmplx(fac1,0.)*exfac
              else if (var(1:lgth).eq.'xi') then
                fac1 = - m*ak(i)/(ak(i)**2+al(j)**2)
                fac2 = m*fcor*al(j)/omga/(ak(i)**2+al(j)**2)
                pfac = cmplx(fac1,fac2)*exfac
              else if (var(1:lgth).eq.'eta') then
                fac1 = - m*al(j)/(ak(i)**2+al(j)**2)
                fac2 = - m*fcor*ak(i)/omga/(ak(i)**2+al(j)**2)
                pfac = cmplx(fac1,fac2)*exfac
              else 
                write(6,*) 'Error:  field not recognized.'
                stop
              end if 

            else

              gam = sqrt(gam2)
              exfac = cmplx(exp(-gam*z),0.)/(nx-1)/(ny-1)
              if (var(1:lgth).eq.'zeta') then
                pfac = exfac
              else if (var(1:lgth).eq.'w') then
                pfac = cmplx(0.,omga)*exfac
              else if (var(1:lgth).eq.'u') then
                fac1 = gam*ak(i)*omga/(ak(i)**2+al(j)**2)
                fac2 = - gam*fcor*al(j)/(ak(i)**2+al(j)**2)
                pfac = cmplx(fac1,fac2)*exfac
              else if (var(1:lgth).eq.'v') then
                fac1 = gam*al(j)*omga/(ak(i)**2+al(j)**2)
                fac2 = gam*fcor*ak(i)/(ak(i)**2+al(j)**2)
                pfac = cmplx(fac1,fac2)*exfac
              else if (var(1:lgth).eq.'p') then
                fac1 = rhos*(N**2-omga**2)/gam
                pfac = cmplx(fac1,0.)*exfac
              else if (var(1:lgth).eq.'rho') then
                fac1 = rhos*N**2/g
                pfac = cmplx(fac1,0.)*exfac
              else if (var(1:lgth).eq.'xi') then
                fac1 = - gam*fcor*al(j)/omga/(ak(i)**2+al(j)**2)
                fac2 = - gam*ak(i)/(ak(i)**2+al(j)**2)
                pfac = cmplx(fac1,fac2)*exfac
              else if (var(1:lgth).eq.'eta') then
                fac1 = gam*fcor*ak(i)/omga/(ak(i)**2+al(j)**2)
                fac2 = - gam*al(j)/(ak(i)**2+al(j)**2)
                pfac = cmplx(fac1,fac2)*exfac
              else 
                write(6,*) 'Error:  field not recognized.'
                stop
              end if 

            end if

            tran(j,i) = pfac*hft(i,j)  

          enddo

          call cfftb(ny-1,tran(1,i),cwrky) 
          do j=1,ny-1
            wrk(i,j) = tran(j,i)
          enddo

        enddo

        do j=1,ny-1
          call cfftb(nx-1,wrk(1,j),cwrkx)
          do i=1,nx-1
            fld2(i,j,k) = real(wrk(i,j))
          enddo
        enddo

        if (iu.eq.1) then
          do j=1,ny-1
            fld2(nx,j,k) = fld2(1,j,k)
          enddo
        else if (iv.eq.1) then
          do i=1,nx-1
            fld2(i,ny,k) = fld2(i,1,k)
          enddo
        end if

      enddo

      do j=1,ny-1+iv
      do i=1,nx-1+iu

        z0 = zs(i,j)
        tfac = (ztop-z0)/ztop
        kptr = 1
        
        do k=1,nz-1+iw

          z = z0 + tfac*gz(k)

          kplus = 0
 589      if ((kplus.eq.0).and.(kptr.le.nz-1+iw)) then
            if (gz(kptr).gt.z) then
              kplus = kptr
            else
              kptr = kptr + 1
            end if
            goto 589
          end if

          if (kplus.eq.1) then
            write(6,*)
            write(6,*) 'Error:  problem with terrain interpolation.'
            write(6,*)
            stop
          else if (kplus.eq.0) then
            kplus = nz-1+iw
          end if

          fld(i,j,k) = fld2(i,j,kplus) + 
     >                 (fld2(i,j,kplus)-fld2(i,j,kplus-1))
     >                       /(gz(kplus)-gz(kplus-1))
     >                            *(z-gz(kplus))

        enddo

      enddo
      enddo

      return
      end

*---------------------------------------------------------------------

      subroutine drag_calc(p,hft,tran,wrk,zs_u,zs_v,cwrkx,cwrky,ak,al,
     >                       fcor,N,U0,V0,rhos,dx,dy,nx,ny,nz)

      integer nx,ny,nz
      real p(nx-1,ny-1),zs_u(nx,ny-1),zs_v(nx-1,ny),ak(nx-1),al(ny-1)
      real cwrkx(4*nx+15),cwrky(4*ny+15),fcor,N,U0,V0,rhos,dy
      complex hft(nx-1,ny-1),tran(ny-1,nx-1),wrk(nx-1,ny-1)

      real m,gam,m2,gam2,omga,fac,sumx,sumy
      logical requ
      complex pfac

      do i=1,nx-1

        do j=1,ny-1 

          omga = ak(i)*U0+al(j)*V0
          if (requ(abs(omga),abs(fcor)).or.(omga.eq.0.)) then
            m2 = 0
            gam2 = 0
          else
            m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >              / (omga**2-fcor**2)
            gam2 = -m2
          end if 

          if (requ(abs(omga),abs(fcor)).or.(omga.eq.0.)) then
            pfac = 0.
          else if (m2.gt.0.) then
            m = sqrt(m2)*omga/abs(omga)
            fac = rhos*(N**2-omga**2)/m/(nx-1)/(ny-1)
            pfac = cmplx(0.,fac)
          else
            gam = sqrt(gam2) 
            fac = rhos*(N**2-omga**2)/gam/(nx-1)/(ny-1)
            pfac = cmplx(fac,0.)
          end if

          tran(j,i) = pfac*hft(i,j)

        enddo

        call cfftb(ny-1,tran(1,i),cwrky)
        do j=1,ny-1
          wrk(i,j) = tran(j,i)
        enddo

      enddo

      do j=1,ny-1
        call cfftb(nx-1,wrk(1,j),cwrkx)
        do i=1,nx-1
          p(i,j) = real(wrk(i,j))
        enddo
      enddo

      sumx = 0.
      sumy = 0.

      do i=1,nx-1
      do j=1,ny-1

        sumx = sumx + p(i,j)*dy*(zs_u(i+1,j)-zs_u(i,j))
        sumy = sumy + p(i,j)*dx*(zs_v(i,j+1)-zs_v(i,j))

      enddo
      enddo

      write(6,*) 'Net x Drag = ',sumx
      write(6,*) 'Net y Drag = ',sumy 

      return
      end

*----------------------------------------------------------------------

      subroutine set_outfl(idcdf,nx,ny,nz,dx,dy,dz,x0,y0,sgz,wgz,
     >                    asclr,sclr,var,units,outfl,nsclr,nvar)

      include '/usr/local/include/netcdf.inc'

      integer idcdf,nx,ny,nz,nsclr,nvar
      real dx,dy,dz,x0,y0,sgz(nz),wgz(nz),asclr(nsclr)
      character*(*) sclr(nsclr),var(nvar),units(nvar),outfl

      real xmin,ymin,zmin,tmp
      character*(80) cdlfl,command,command2

      lgth = index(outfl,' ') - 1

      cdlfl(1:len(cdlfl)) = ' '
      cdlfl = outfl(1:lgth)//'.cdl'

      iunit = 48

      open(unit=iunit,file=cdlfl(1:lgth+4),status='new')

      write(iunit,*) 'netcdf lin_soln{'
      write(iunit,*)

      write(iunit,*) 'dimensions:'
      write(iunit,*) 'nx=',nx-1,';'
      write(iunit,*) 'ny=',ny-1,';'
      write(iunit,*) 'nz=',nz-1,';'
      write(iunit,*) 'nxp1=',nx,';'
      write(iunit,*) 'nyp1=',ny,';'
      write(iunit,*) 'nzp1=',nz,';'
      write(iunit,*) 'one=',1,';'
      write(iunit,*) 'time=UNLIMITED;'
      write(iunit,*)

      write(iunit,*) 'variables:'
      write(iunit,*)

      do k=1,nsclr
        lgth = index(sclr(k),' ') - 1
        write(iunit,*) 'float ',sclr(k)(1:lgth),'(one);'
        write(iunit,*) sclr(k)(1:lgth),':no_button=1;'
        write(iunit,*)
      enddo

      write(iunit,*) 'float time(time);'
      write(iunit,*) 'time:units="s";'
      write(iunit,*)
      write(iunit,*) 'float sgz(nz);'
      write(iunit,*) 'sgz:units="m";'
      write(iunit,*) 'sgz:def="staggered grid pts in z";'
      write(iunit,*) 'sgz:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float wgz(nzp1);'
      write(iunit,*) 'wgz:units="m";'
      write(iunit,*) 'wgz:def="unstaggered grid pts in z";'
      write(iunit,*) 'wgz:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_p(ny,nx);'
      write(iunit,*) 'zbot_p:units="m";'
      write(iunit,*) 'zbot_p:def="terrain at p points";'
      write(iunit,*) 'zbot_p:x_min=',x0+0.5*dx,';'
      write(iunit,*) 'zbot_p:y_min=',y0+0.5*dy,';'
      write(iunit,*) 'zbot_p:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_u(ny,nxp1);'
      write(iunit,*) 'zbot_u:units="m";'
      write(iunit,*) 'zbot_u:def="terrain at u points";'
      write(iunit,*) 'zbot_u:x_min=',x0,';'
      write(iunit,*) 'zbot_u:y_min=',y0+0.5*dy,';'
      write(iunit,*) 'zbot_u:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_v(nyp1,nx);'
      write(iunit,*) 'zbot_v:units="m";'
      write(iunit,*) 'zbot_v:def="terrain at v points";'
      write(iunit,*) 'zbot_v:x_min=',x0+0.5*dx,';'
      write(iunit,*) 'zbot_v:y_min=',y0,';'
      write(iunit,*) 'zbot_v:no_button=1;'
      write(iunit,*)

      do k=1,nvar

        xmin = x0 + 0.5*dx
        ymin = y0 + 0.5*dy
        zmin = sgz(1)
 
        lgth = index(var(k),' ') - 1

        if (var(k)(1:1).eq.'u') then
          write(iunit,*) 'float ',var(k)(1:lgth),
     >                     '(time,nz,ny,nxp1);'
          xmin = x0
        else if (var(k)(1:1).eq.'v') then
          write(iunit,*) 'float ',var(k)(1:lgth),
     >                     '(time,nz,nyp1,nx);'
          ymin = y0 
        else if (var(k)(1:1).eq.'w') then
          write(iunit,*) 'float ',var(k)(1:lgth),
     >                     '(time,nzp1,ny,nx);'
          zmin = wgz(1)
        else
          write(iunit,*) 'float ',var(k)(1:lgth),
     >                     '(time,nz,ny,nx);'
        end if

        lgth2 = index(units(k),' ') - 1
        if (lgth2.eq.0) lgth2 = 1
        write(iunit,*) var(k)(1:lgth),':units="',
     >              units(k)(1:lgth2),'";'     

        write(iunit,*) var(k)(1:lgth),':x_min=',xmin,';'
        write(iunit,*) var(k)(1:lgth),':y_min=',ymin,';'
        write(iunit,*) var(k)(1:lgth),':z_min=',zmin,';'

        if (var(k)(lgth:lgth).eq.'0') then
          write(iunit,*) var(k)(1:lgth),':no_button=1;'
        end if

        write(iunit,*) 

      enddo

      write(iunit,*) '//global attributes:'
      write(iunit,*)
      write(iunit,*) ':x_min=',x0,';'
      write(iunit,*) ':x_max=',x0 + (nx-1)*dx,';'
      write(iunit,*) ':x_delta=',dx,';'
      write(iunit,*) ':x_units="m";'
      write(iunit,*) ':x_label="x";'
      write(iunit,*) ':x_display_units="km";'
      write(iunit,*)
      write(iunit,*) ':y_min=',y0,';'
      write(iunit,*) ':y_max=',y0 + (ny-1)*dy,';'
      write(iunit,*) ':y_delta=',dy,';'
      write(iunit,*) ':y_units="m";'
      write(iunit,*) ':y_label="y";'
      write(iunit,*) ':y_display_units="km";'
      write(iunit,*)
      write(iunit,*) ':z_min=',0.,';'
      write(iunit,*) ':z_max=',0. + (nz-1)*dz,';'
      write(iunit,*) ':z_delta=',dz,';'
      write(iunit,*) ':z_units="m";'
      write(iunit,*) ':z_label="z";'
      write(iunit,*) ':z_display_units="km";'
      write(iunit,*)

      lgth = index(outfl,' ') - 1
      write(iunit,*) ':runname="steady linear solution ',
     >                     outfl(1:lgth),'";'
      write(iunit,*)
      write(iunit,*) '}'

      close(iunit)

      outfl = outfl(1:lgth)//'.cdf'

      command(1:len(command)) = ' '
      write(command,*) 'ncgen -o ',outfl(1:lgth+4),' ',
     >                         cdlfl(1:lgth+4)
      call system(command)
      write(command2,*) 'rm -f ',cdlfl(1:lgth+4)
      call system(command2)

      idcdf = ncopn(outfl(1:lgth+4),NCWRITE,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem opening netcdf file.'
        stop
      end if

      idvar = ncvid(idcdf,'time',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting variable id for time.'
        stop
      end if
      call ncvpt1(idcdf,idvar,1,0.,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem writing time.'
        stop
      end if

      do k=1,nsclr
        lgth = index(sclr(k),' ') - 1
        idvar = ncvid(idcdf,sclr(k)(1:lgth),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem getting variable id for ',
     >                                 'scalar.'
          stop
        end if
        call ncvpt1(idcdf,idvar,1,asclr(k),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem writing scalar.'
          stop
        end if
      enddo

      idvar = ncvid(idcdf,'sgz',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for sgz.'
        stop
      end if

      call ncvpt(idcdf,idvar,1,nz-1,sgz,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem writing sgz.'
        stop
      end if

      idvar = ncvid(idcdf,'wgz',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for wgz.'
        stop
      end if

      call ncvpt(idcdf,idvar,1,nz,wgz,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem writing wgz.'
        stop
      end if

      return
      end

*------------------------------------------------------------------------

      subroutine get_ter(asclr,zs_p,zs_u,zs_v,x0,y0,dx,dy,
     >                               nx,ny,nsclr)

      integer nsclr,nx,ny
      real asclr(nsclr),zs_p(nx-1,ny-1),zs_u(nx,ny-1),zs_v(nx-1,ny)
      real x0,y0,dx,dy

      real x,xu,xv,y,yu,yv,pi,ys,yn,ytmp,rad,a,b 

      real type,h,a1,a2,b1,b2,bv,y_len,xc,yc

      type = asclr(3)
      h = asclr(4)
      a1 = asclr(5)
      a2 = asclr(6)
      b1 = asclr(7)
      b2 = asclr(8)
      y_len = asclr(9)
      xc = asclr(10)
      yc = asclr(11)

      pi = 2.*asin(1.)

      if (type.eq.0.) then

        do j = 1,ny-1
        do i = 1,nx-1
         zs_p(i,j) = 0.
         zs_u(i,j) = 0.
         zs_v(i,j) = 0.
        enddo
        enddo
  
        do j=1,ny-1
          zs_u(nx,j) = 0.
        enddo
        do i=1,nx-1
          zs_v(i,ny) = 0.
        enddo
 
      else if (type.eq.1.) then

        do i=1,nx-1
          x  = x0 + dx*(i-1) + 0.5*dx - xc
          xu = x0 + dx*(i-1) - xc
          xv = x
        do j=1,ny-1
          zs_p(i,j) = h*a1**2 / (x**2+a1**2)
          zs_u(i,j) = h*a1**2 / (xu**2+a1**2)
          zs_v(i,j) = zs_p(i,j)
        enddo
        enddo

        xu = x0 + dx*(nx-1) - xc
        do j=1,ny-1
          zs_u(nx,j) = h*a1**2 / (xu**2+a1**2)
        enddo

        do i=1,nx-1
          xv = x0 + dx*(i-1) + 0.5*dx - xc
          zs_v(i,ny) = h*a1**2 / (xv**2+a1**2)
        enddo

      else if (type.eq.2.) then

        do j=1,ny-1
          y  = y0 + dy*(j-1) + 0.5*dy - yc
          yv = y0 + dy*(j-1) - yc
          yu = y
        do i=1,nx-1
          zs_p(i,j) = h*b1**2 / (y**2+b1**2)
          zs_u(i,j) = zs_p(i,j)
          zs_v(i,j) = h*b1**2 / (yv**2+b1**2)
        enddo
        enddo

        do j=1,ny-1
          yu = y0 + dy*(j-1) + 0.5*dy - yc
          zs_u(nx,j) = h*b1**2 / (yu**2+b1**2)
        enddo

        yv = y0 + dy*(ny-1) - yc
        do i=1,nx-1 
          zs_v(i,ny) = h*b1**2 / (yv**2+b1**2)
        enddo

      else if (type.eq.3.) then

        do j=1,ny-1
          y  = y0 + dy*(j-1) + 0.5*dy - yc
          yv = y0 + dy*(j-1) - yc
          yu = y
        do i=1,nx-1
          x  = x0 + dx*(i-1) + 0.5*dx - xc 
          xu = x0 + dx*(i-1) - xc
          xv = x
          zs_p(i,j) = h/(1+(x/a1)**2+(y/b1)**2)**(1.5)
          zs_u(i,j) = h/(1+(xu/a1)**2+(yu/b1)**2)**(1.5)
          zs_v(i,j) = h/(1+(xv/a1)**2+(yv/b1)**2)**(1.5)
        enddo
        enddo

        xu = x0 + dx*(nx-1) - xc
        do j=1,ny-1
          yu = y0 + dy*(j-1) + 0.5*dy - yc
          zs_u(nx,j) = h/(1+(xu/a1)**2+(yu/b1)**2)**(1.5)
        enddo

        yv = y0 + dy*(ny-1) - yc
        do i=1,nx-1
          xv = x0 + (i-1)*dx + 0.5*dx - xc
          zs_v(i,ny) = h/(1+(xv/a1)**2+(yv/b1)**2)**(1.5)
        enddo
          
      else if (type.eq.4.) then

        ys = yc - y_len
        yn = yc + y_len

        do j=1,ny

          b = b1
          bv = b1

          ytmp = y0 + dy*(j-1) + 0.5*dy
          if (ytmp.lt.ys) then
            y = ys - ytmp
            b = b2
          else if (ytmp.gt.yn) then
            y = ytmp - yn
          else
            y = 0.
          end if

          ytmp = ytmp - 0.5*dy
          if (ytmp.lt.ys) then
            yv = ys - ytmp
            bv = b2
          else if (ytmp.gt.yn) then
            yv = ytmp - yn
          else
            yv = 0.
          end if

          yu = y

          do i=1,nx

            x = x0 + dx*(i-1) + 0.5*dx - xc
            xu = x - 0.5*dx
            xv = x

            if (x.le.0.) then 
              a = a1
            else
              a = a2
            end if

            rad = sqrt((x/4./a)**2 + (y/4./b)**2)

            if ((i.le.nx-1).and.(j.le.ny-1)) then
              if (rad.lt.1.) then
                zs_p(i,j) = h/16.*(1+cos(pi*rad))**4
              else
                zs_p(i,j) = 0.
              end if
            end if

            rad = sqrt((xv/4./a)**2 + (yv/4./bv)**2)

            if (i.le.nx-1) then
              if (rad.lt.1.) then
                zs_v(i,j) = h/16.*(1+cos(pi*rad))**4
              else
                zs_v(i,j) = 0.
              end if
            end if

            if (xu.le.0.) then
              a = a1
            else
              a = a2
            end if

            rad = sqrt((xu/4./a)**2 + (yu/4./b)**2)

            if (j.le.ny-1) then
              if (rad.lt.1.) then
                zs_u(i,j) = h/16.*(1+cos(pi*rad))**4
              else
                zs_u(i,j) = 0.
              end if
            end if

          enddo

        enddo

      else if (type.eq.5) then

        do i=1,nx-1
          x  = x0 + dx*(i-1) + 0.5*dx - xc
          xu = x0 + dx*(i-1) - xc
          xv = x
        do j=1,ny-1
          if (x.le.0.) then
            a = a1
          else
            a = a2
          end if
          if (abs(x/4./a).lt.1.) then
            zs_p(i,j) = h/16.*(1.+cos(pi*x/4./a))**4
            zs_v(i,j) = zs_p(i,j)
          else
            zs_p(i,j) = 0.
            zs_v(i,j) = 0.
          end if
          if (xu.le.0.) then
            a = a1
          else
            a = a2
          end if
          if (abs(xu/4./a).lt.1.) then
            zs_u(i,j) = h/16.*(1.+cos(pi*xu/4./a))**4 
          else
            zs_u(i,j) = 0.
          end if 
        enddo
        enddo

        xu = x0 + dx*(nx-1) - xc
        if (xu.le.0.) then
          a = a1
        else
          a = a2
        end if
        do j=1,ny-1
          if (abs(xu/4./a).lt.1) then
            zs_u(nx,j) = h/16.*(1.+cos(pi*xu/4./a))**4 
          else
            zs_u(nx,j) = 0.
          end if
        enddo

        do i=1,nx-1
          xv = x0 + dx*(i-1) + 0.5*dx - xc
          if (xv.le.0.) then
            a = a1
          else
            a = a2
          end if
          if (abs(xv/4./a).lt.1.) then
            zs_v(i,ny) = h/16.*(1.+cos(pi*xv/4./a))**4
          else
            zs_v(i,ny) = 0.
          end if
        enddo

      else
        
 554    FORMAT('Error:  terrain type ',F5.1,' not yet implemented.')
        write(*,554) type
        stop
 
      end if
 
      return
      end

*------------------------------------------------------------------------

      subroutine write_var(idcdf,vname,nx,ny,nz,avar)

      integer idcdf,nx,ny,nz
      real avar(nx*ny*nz)
      character*(*) vname

      integer istrt(4),ilen(4),ist2(2),iln2(2)

      lgth = index(vname,' ') - 1
      if (lgth.eq.-1) lgth = len(vname)

      idvar = ncvid(idcdf,vname(1:lgth),ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting variable id.'
        stop
      end if

      istrt(1) = 1
      istrt(2) = 1
      istrt(3) = 1
      istrt(4) = 1
      ist2(1) = 1
      ist2(2) = 1

      ilen(1) = nx-1
      ilen(2) = ny-1
      ilen(3) = nz-1
      ilen(4) = 1
      iln2(1) = nx-1
      iln2(2) = ny-1

      if (vname(1:4).eq.'zbot') then
        if (vname(lgth:lgth).eq.'u') then
          iln2(1) = nx
        else if (vname(lgth:lgth).eq.'v') then
          iln2(2) = ny
        end if
        call ncvpt(idcdf,idvar,ist2,iln2,avar,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem writing terrain variable.'
          stop
        end if
      else
        if (vname(1:1).eq.'u') then
          ilen(1) = nx
        else if (vname(1:1).eq.'v') then
          ilen(2) = ny
        else if (vname(1:1).eq.'w') then
          ilen(3) = nz
        end if
        call ncvpt(idcdf,idvar,istrt,ilen,avar,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem writing variable.'
          stop
        end if
      end if

      return
      end

*-------------------------------------------------------------------------

      subroutine stretch_it(sgz,wgz,dz,dzbot,nz)

      integer nz
      real sgz(nz),wgz(nz),dz,dzbot

      real ztop,step,stretch,dzmax,mid,diff
      parameter (dzmax=700.)

      ztop = dz*(nz-1)
      step = 1.0
      stretch = 1.0

      if (dz.gt.dzbot) then
        k = 1
 32     if (abs(step).gt.(1.0e-10)) then
          step = step*0.5
          mid = stretch + step
          diff = zhgt(dzbot,mid,nz,dzmax) - ztop
          if (diff.le.0.0) stretch = mid
          k = k+1
          if (k.gt.50) then
            write(6,*) 'Problem calculating stretching...'
            stop
          end if
          goto 32
        end if
      end if

      wgz(1) = 0.0
      do k=2,nz
        wgz(k) = zhgt(dzbot,stretch,k,dzmax)
        sgz(k-1) = (wgz(k)+wgz(k-1))/2.
      enddo
      sgz(nz) = wgz(nz)
 
      return
      end
 
*----------------------------------------------------------------------------

      real function zhgt(dzbot,stretch,kz,dzmax)

      integer kz
      real dzbot,stretch,dzmax

      real sum,dznew

      sum = 0.
      do k=1,kz-1
        dznew = amin1(dzbot*(stretch**(k-1)),dzmax)
        sum = sum + dznew
      enddo

      zhgt = sum

      return
      end 

*---------------------------------------------------------------------------

      logical function requ(val1,val2)

      real val1,val2,eps,diff,big	

      parameter (eps=1.e-4)

      big = amax1(abs(val1),abs(val2))

      if (big.eq.0.) then
        requ = .true.
        return
      end if

      diff = abs((val1-val2)/big)

      if (diff.le.eps) then
        requ = .true.
      else
        requ = .false.
      end if

      return
      end

 
     
