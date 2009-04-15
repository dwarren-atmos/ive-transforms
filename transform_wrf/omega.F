!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calculate the vertical velocity on the untransformed grid, omega
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(ive_ptr) function omega(presnm,tempnm,wnm,ndims,
     &                       dims,stag,min,max,miss,
     &                       dunits,disunits,dim_names)
      use default_module

      implicit none

      integer ndims,dims(4)
	  real stag(4),min(4),max(4),miss
      character(len=*) :: presnm,tempnm,wnm
      character(len=80) ::dunits,disunits,dim_names(4)
      character *80 string,tstr,wstr
      integer ntotal,utotal,vtotal,i
	  logical newp,newt,neww
      integer(ive_ptr) pvar,tvar,wvar




c
c
c  Get Temp@P variable
c
      tvar = 0
      tstr=' '
      write(tstr,*)trim(tempnm),'@',trim(presnm)
      tvar = getvar(tstr,ndims,dims,stag,min,max,miss,dunits,
     &                  disunits,dim_names,newt)
c
      if(tvar.eq.0) then
        write(6,*)  'Error:  cannot access ',tempnm
        return
      endif
c
c
c
c  Get W@P variable
c
      wvar = 0
      wstr=' '
      write(wstr,*)trim(wnm),'@',trim(presnm)
      wvar = getvar(wstr,ndims,dims,stag,min,max,miss,dunits,
     &              disunits,dim_names,neww)
c
      if(wvar.eq.0) then
        write(6,*)  'Error:  cannot access ',wnm
        return
      endif

c
c  Get Pressure variable
c
      pvar = 0
      pvar = getvar(presnm,ndims,dims,stag,min,max,miss,dunits,
     &          disunits,dim_names,newp)
      if(pvar.eq.0) then
        write(6,*)  'Error:  cannot access ',presnm
        return
      endif

c
c
      omega = 0 ; omega = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if(omega.le.0) then
        write(6,*) 'Error:  Cannot allocate memory for omega'
        return
      end if

      call omega_calc(%val(omega),%val(pvar),%val(tvar),
     &                %val(wvar),dims(1),dims(2),dims(3),dims(4))
c

c
c     Already correct from pressure
!      dims(1) = nx+1
!      dims(2) = ny+1
!      dims(3) = nz+1
      dunits = 'Pa/s'
      disunits = 'Pa/s'
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine omega_calc(omega,p,t,w,nx,ny,nz,nt)
      implicit none
c
      integer nx,ny,nz,nt
      real omega(nx,ny,nz,nt),p(nx,ny,nz,nt),
     &     t(nx,ny,nz,nt),w(nx,ny,nz,nt)
c
      integer i,j,k,n
c
      do n=1,nt ; do k=1,nz ; do j=1,ny ; do i=1,nx

c
c   omega is 0 at the bottom boundry
c   omega = -Dens*G*w
c   Dens = p/RT where R=8.31432*10 for mb pressure
        if(k.eq.1) then
          omega(i,j,k,n) = 0.
        else
          omega(i,j,k,n) =  -1.*(p(i,j,k,n)/(83.1432 * t(i,j,k,n))) *
     &                   9.8*w(i,j,k,n)
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

