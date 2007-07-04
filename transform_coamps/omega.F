#include "pointer.icl"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ive_ptr function omega(uvarnm,vvarnm,wvarnm,ndims,
     &                       dims,stag,min,max,miss,
     &                       dunits,disunits,dim_names)
      use default_module

      implicit none
c
c
      integer ndims,dims(4)
	  real stag(4),min(4),max(4),miss
      character(len=*) :: uvarnm,vvarnm,wvarnm
      character(len=*) ::dunits,disunits,dim_names(4)
c
      integer ntotal,utotal,vtotal,i
	  logical new
      ive_ptr uvartmp,vvartmp,uvarin,vvarin,wvarin
      ive_ptr uvar,vvar,wvar
      character(len=10) :: freevar
c
c     new = .true.
      ntotal = (nx+1)*(ny+1)*(nz+1)*ntime
      utotal = (nx)*(ny+1)*(nz+1)*ntime
      vtotal = (nx+1)*(ny)*(nz+1)*ntime

c
c  Get U variable
c
c
      uvarin = 0
      uvarin = getvar(uvarnm,ndims,dims,stag,min,max,miss,dunits,
     &          disunits,dim_names,new)
c
      if(uvarin.eq.0) then
        write(6,*)  'Error:  cannot access ',uvarnm
        return
      endif

      uvartmp = 0
      uvartmp = getmem(utotal)
      if(uvartmp.eq.0) then
        write(6,*) 'Error:  Cannot allocate memory in for uvartmp'
        return
      end if
c
      call mult_var(%val(uvartmp),%val(uvarin)
     &             ,%val(dqdx_xy),1,0,0,nx,ny,nz,ntime)
c
      uvar = 0
      uvar = getmem(ntotal)
      if(uvar.eq.0) then
        write(6,*) 'Error:  Cannot allocate memory for ',uvarnm
        return
      end if
c
      do i=1,ntotal
        call s_val(%val(uvar),i,1.)
      end do
c
      call var_at_w(%val(uvar),%val(uvartmp),%val(sgz),%val(wgz),
     &              1,0,0,nx,ny,nz,ntime)
c
c      write(freevar,'(A,A)') 'free = ',uvarnm
c      call driver(freevar,0)
c
c  Get V variable
c
      vvarin = 0
      vvarin = getvar(vvarnm,ndims,dims,stag,min,max,miss,dunits,
     &                  disunits,dim_names,new)
c
      if(vvarin.eq.0) then
        write(6,*)  'Error:  cannot access ',vvarnm
        return
      endif
c
      vvartmp = 0
      vvartmp = getmem(vtotal)
      if(vvartmp.eq.0) then
        write(6,*) 'Error:  Cannot allocate memory in for vvartmp'
        return
      end if
c
      call mult_var(%val(vvartmp),%val(vvarin)
     &             ,%val(dqdy_xy),0,1,0,nx,ny,nz,ntime)
c
      vvar = 0
      vvar = getmem(ntotal)
      if(vvar.eq.0) then
        write(6,*) 'Error:  Cannot allocate memory for ', vvarnm
        return
      end if
c
      do i=1,ntotal
        call s_val(%val(vvar),i,1.)
      end do
c
      call var_at_w(%val(vvar),%val(vvartmp),%val(sgz),%val(wgz),
     &              0,1,0,nx,ny,nz,ntime)
c
c      write(freevar,'(A,A)') 'free = ',vvarnm
c      call driver(freevar,0)

c
c  Get W variable
c
      wvarin = 0
      wvarin = getvar(wvarnm,ndims,dims,stag,min,max,miss,dunits,
     &              disunits,dim_names,new)
c
      if(wvarin.eq.0) then
        write(6,*)  'Error:  cannot access ',wvarnm
        return
      endif
c
      wvar = 0
      wvar = getmem(ntotal)
      if(wvar.eq.0) then
        write(6,*) 'Error:  Cannot allocate memory for wvar'
        return
      end if
c
      call mult_var(%val(wvar),%val(wvarin)
     &             ,%val(dqdz),0,0,1,nx,ny,nz,ntime)
c

      omega = 0
      omega = getmem(ntotal)
      if(omega.eq.0) then
        write(6,*) 'Error:  Cannot allocate memory for omega'
        return
      end if
c
      call omega_calc(%val(omega),%val(uvar),%val(vvar),
     &                %val(wvar),%val(dqdx_xy),%val(dqdy_xy),
     &                %val(dqdx_zw),%val(dqdy_zw),%val(dqdz),
     &                nx,ny,nz,ntime)
c
c      write(freevar,'(A,A)') 'free = ',wvarnm
c      call driver(freevar,0)

      call freemem(uvar) !; call freemem(uvartmp)
      call freemem(vvar) !; call freemem(vvartmp)
c
      dims(1) = nx+1
      dims(2) = ny+1
      dims(3) = nz+1
      dunits = 'm/s'
      disunits = 'm/s'
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine omega_calc(omega,u,v,w,dqdx_xy,dqdy_xy,dqdx_zw,
     &                     dqdy_zw,dqdz_xy,nx,ny,nz,nt)
      implicit none
c
      integer nx,ny,nz,nt
      real omega(nx+1,ny+1,nz+1,nt),u(nx+1,ny+1,nz+1,nt),
     &     v(nx+1,ny+1,nz+1,nt),w(nx+1,ny+1,nz+1,nt),
     &     dqdx_xy(nx+1,ny+1,3),dqdy_xy(nx+1,ny+1,3),
     &     dqdx_zw(nz+1),dqdy_zw(nz+1),dqdz_xy(nx+1,ny+1,3)
c
      integer i,j,k,n
      real dqdx,dqdy,dqdz
c
      do n=1,nt ; do k=1,nz+1 ; do j=1,ny+1 ; do i=1,nx+1
c         dqdx = dqdx_xy(i,j,1)*dqdx_zw(k)
c         dqdy = dqdy_xy(i,j,1)*dqdy_zw(k)
c         dqdz = dqdz_xy(i,j,1)
c
        dqdx = dqdx_zw(k)
        dqdy = dqdy_zw(k)
        dqdz = 1.0
c
c   omega is 0 at the bottom boundry
c
        if(k.eq.1) then
          omega(i,j,k,n) = 0.
        else
          omega(i,j,k,n) = 
     &                  dqdx * u(i,j,k,n) +
     &                  dqdy * v(i,j,k,n) +
     &                  dqdz * w(i,j,k,n)
        end if
c
      end do ; end do ; end do ; end do
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mult_var(var_out,var1,var2,iu,iv,iw,nx,ny,nz,nt)
      implicit none
c
      integer nx,ny,nz,nt,iu,iv,iw
      real var_out(nx+1-iu,ny+1-iv,nz+iw,nt)
      real var1(nx+1-iu,ny+1-iv,nz+iw,nt),var2(nx+1,ny+1,3) 
c
      integer i,j,k,n,ipnt
c
      ipnt = 1
      if(iu.eq.1.and.iv.eq.0) ipnt = 2
      if(iu.eq.0.and.iv.eq.1) ipnt = 3
c
      do n=1,nt ; do k=1,nz+iw ; do j=1,ny+1-iv ; do i=1,nx+1-iu
        var_out(i,j,k,n) = var1(i,j,k,n)*var2(i,j,ipnt)
      end do ; end do ; end do ; end do
c
	return
	end

