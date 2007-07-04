#include "pointer.icl"
c
      ive_ptr function massflx_lat(x1,x2,y1,y2,k1,k2,
     &                 ndims,dims,stag,min,max,miss,
     &                 dunits,disunits,dim_names)
      use default_module

      implicit none
c
      integer ndims,dims(4)
      real x1,x2,y1,y2,k1,k2,stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
c
      integer nmax,numpts
      parameter(nmax = 1000, numpts=2)
c
      integer tcnt,ucnt,vcnt,indx(numpts,2),sgnni,sgnnj,iu,iv,iw
      integer tindex(nmax,2),uindex(nmax,2),vindex(nmax,2),i,j,k
c
	  logical new
c
      ive_ptr uvar,vvar,rho,rho_calc
c
      massflx_lat = 0
c
      massflx_lat = getmem(ntime)
      if (massflx_lat.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field massflx.'
        return
      end if
c
c  Convert phycical points to index (could use phys_2_index)
c
      indx(1,1) = nint(x1/(dx/1000.))+1
      indx(2,1) = nint(x2/(dx/1000.))+1
      indx(1,2) = nint(y1/(dy/1000.))+1
      indx(2,2) = nint(y2/(dy/1000.))+1
c
      write(6,*) 'Calculating the lateral mass flux between'
      write(6,*) 'x1 = ',x1,'y1 = ',y1
      write(6,*) 'x2 = ',x2,'y2 = ',y2
c
      call tline(indx(1,:),indx(2,:),tindex,tcnt,nmax,sgnni,sgnnj)
      call uvline(vindex,uindex,tindex,tcnt,ucnt,vcnt,nmax,sgnni,sgnnj)
c
      uvar = getvar('UTR',ndims,dims,stag,min,max,miss,dunits,
     &              disunits,dim_names,new)
c
      if(uvar.eq.0) then
        write(6,*) 'Error:  Cannot access field UTR'
        return
      endif
c
      vvar = getvar('VTR',ndims,dims,stag,min,max,miss,dunits,
     &              disunits,dim_names,new)
c
      if(vvar.eq.0) then
        write(6,*) 'Error:  Cannot access field VTR'
        return
      endif
c
      rho = 0
      rho = rho_calc(ndims,dims,stag,min,max,miss,dunits,
     &               disunits,dim_names)
c
      call var_at_u(%val(uvar),%val(rho),0,0,0,nx,ny,nz,ntime)
      call var_at_v(%val(vvar),%val(rho),0,0,0,nx,ny,nz,ntime)
c
c     call freemem(rho)
c
      call massflx_lat_calc(uindex,vindex,k1,k2,
     &     %val(massflx_lat),%val(uvar),%val(vvar),
     &     %val(wgz),%val(dzdq),dx,dy,ucnt,vcnt,
     &     sgnni,sgnnj,nx,ny,nz,ntime,nmax)
c
      dunits = 'kg/s' ; disunits = 'kg/s'
c
      do k = 1,3
        dims(k) = 1
        stag(k) = 0.
        min(k) = 0. ;  max(k) = 0.
        dim_names(k) = 'one'
      enddo
c
      return
      end subroutine
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	subroutine massflx_lat_calc(uindex,vindex,k1,k2,
     &	     massflx,u,v,wgz,dzdq,dx,dy,ucnt,vcnt,
     &           sgnni,sgnnj,nx,ny,nz,nt,nmax)
	implicit none
c
    	integer nmax,nx,ny,nz,nt
	integer uindex(nmax,2),vindex(nmax,2),k1,k2,ucnt,vcnt
	integer sgnni,sgnnj
c
	real massflx(nt),u(nx,ny+1,nz,nt),v(nx+1,ny,nz,nt),wgz(nz+1),
     &     dzdq(nx+1,ny+1,3),dx,dy	
c
	integer i,k,n,iloc,jloc
	real dz
c
	do n=1,nt
	  massflx(n) = 0.
	  do k=k1,k2
	    do i=1,ucnt
	      iloc = uindex(i,1)
	      jloc = uindex(i,2)
	 	dz = ( wgz(k+1) - wgz(k) ) * dzdq(iloc,jloc,2)
	      massflx(n) = massflx(n) - 
     &                 dzdq(iloc,jloc,2)*sgnnj*dz*dy*u(iloc,jloc,k,n)
	    enddo
c	    
	    do i=1,vcnt
	      iloc = vindex(i,1)
	      jloc = vindex(i,2)
	 	dz = ( wgz(k+1) - wgz(k) ) * dzdq(iloc,jloc,3)
	      massflx(n) = massflx(n) + 
     &                 dzdq(iloc,jloc,3)*sgnni*dz*dx*v(iloc,jloc,k,n)
	    enddo
	  enddo
	enddo
c
	return
	end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	subroutine uvline(vindex,uindex,tindex,tcnt,ucnt,vcnt,nmax,
     &			sgnni,sgnnj)
	implicit none
c	
	integer tcnt, nmax,sgnni,sgnnj,
     &	  tindex(nmax,2),vindex(nmax,2),uindex(nmax,2)
c
	integer i,di,dj,ucnt,vcnt,ni,nj,sgnm
c
	logical lsthor
c
	integer heaviside
c
	ni = tindex(tcnt,1)-tindex(1,1)
	nj = tindex(tcnt,2)-tindex(1,2)
c
	sgnni = 1
	sgnnj = 1
	if(ni.ne.0) sgnni = abs(ni)/ni
	if(nj.ne.0) sgnnj = abs(nj)/nj
	if(ni.eq.0) sgnni = -sgnnj
	if(nj.eq.0) sgnnj = sgnni
	sgnm = sgnnj/sgnni
c
	ucnt = 0
	vcnt = 0
c
	do i =1,tcnt-1
	  di = tindex(i,1)-tindex(i+1,1)
	  dj = tindex(i,2)-tindex(i+1,2)
c	  
	  if (di.eq.0.and.dj.eq.0) then
	  elseif (dj.eq.0.and.di.ne.0) then
	    vcnt = vcnt + 1
	    vindex(vcnt,1) = tindex(i,1) + sgnni*heaviside(sgnm)
	    vindex(vcnt,2) = tindex(i,2) - heaviside(sgnni)
	  elseif (di.eq.0.and.dj.ne.0) then
	    ucnt = ucnt + 1
	    uindex(ucnt,1) = tindex(i,1) - heaviside(-sgnnj) 
	    uindex(ucnt,2) = tindex(i,2) + sgnnj*heaviside(-sgnm)
	  endif
	enddo
c
	lsthor = .false.
	if(abs(nj).gt.abs(ni)) lsthor = .true.
c	
	if (lsthor) then
	  vcnt = vcnt + 1
	  vindex(vcnt,1) = tindex(tcnt,1)
	  vindex(vcnt,2) = tindex(tcnt,2) - heaviside(-sgnnj)
	else
	  ucnt = ucnt + 1
	  uindex(ucnt,1) = tindex(tcnt,1) - heaviside(-sgnni)
	  uindex(ucnt,2) = tindex(tcnt,2) 
	endif
c
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine tline(indx1,indx2,tindx,count,nmax,sgnni,sgnnj)
	implicit none 
c
c     gets the end points of a line and the grid spacing between endpoints
c
	integer count,ni,nj,nmax,tindx(nmax,2),indx1(2),indx2(2)
c
	integer gcd,dj,di,sgnni,sgnnj,n,i,j,p
c
      integer gcd_calc
c	 
	ni = indx2(1) - indx1(1)
	nj = indx2(2) - indx1(2)
c
      gcd = gcd_calc(nj,ni)
      dj = abs(nj/gcd)
	di = abs(ni/gcd)
c
	sgnni = 1
	sgnnj = 1
	if (ni.ne.0) sgnni = abs(ni)/ni
	if (nj.ne.0) sgnnj = abs(nj)/nj
	if(ni.eq.0) sgnni = -sgnnj
	if(nj.eq.0) sgnnj = sgnni
c
	ni = sgnni*ni
	nj = sgnnj*nj
c
	n = ni
	if((ni.gt.nj.and.nj.ne.0).or.ni.eq.0) n = nj 	
c	 
	count = 1
	tindx(count,1) = indx1(1)
	tindx(count,2) = indx1(2)
c
	do p=1,n
	  do i = 1,abs(di)
	    count = count + 1
	    tindx(count,1) = indx1(1) + sgnni*(di*(p-1) + i)
	    tindx(count,2) = indx1(2) + sgnnj*(p-1)*dj
	  enddo
	  do j = 1,abs(dj)
	    count = count+1
	    tindx(count,1) = indx1(1) + sgnni*p*di
          tindx(count,2) = indx1(2) + sgnnj*(dj*(p-1) + j)
        enddo
	enddo
c	  
	do i=1,2
	 indx2(i) = tindx(count,i)
	enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	integer function gcd_calc(n1,n2)
c
c 	Finds the gcd of num and denom using euclids algorithim
c
	integer n1,n2
c
	integer p,q,r
c
	if (abs(n1).gt.abs(n2)) then
	  p = abs(n1)
	  q = abs(n2)
	else  
	  q = abs(n1)
	  p = abs(n2)
	end if
c
c  Take care of any zeros
c
      if (p.eq.0.and.q.eq.0) then
	  gcd_calc = 1
	  return
	elseif(p.eq.0) then
	  gcd_calc = q
	  return
	elseif(q.eq.0) then
	  gcd_calc = p
	  return
	endif
c
	do 
	  r = mod(p,q)
	  if(r.eq.1) then
	    gcd_calc = 1
	    return
	  elseif(r.eq.0) then
	    gcd_calc = q
	    return
	  else
	    p = q
	    q = r
	  end if
	end do
c	
	return
	end
c	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	integer function heaviside(i)

	integer i

	if (i.gt.0) heaviside = i
	if (i.le.0) heaviside = 0
c
	return
	end
