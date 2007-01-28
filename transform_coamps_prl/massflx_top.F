#include "pointer.icl"
c
	ive_ptr function massflx_top(x1,x2,x3,x4,k2,top,
     >		     ndims,dims,stag,min,max,miss,
     >		     dunits,disunits,dim_names)
	implicit none
c
#include "default.icl"
c
	integer ndims,dims(4),k2
	real x1(2),x2(2),x3(2),x4(2),stag(4),min(4),max(4),miss
	character*(*) dunits,disunits,dim_names(4)
	logical top
c
	integer nmax
	parameter(nmax = 1000)
	integer indx(4,2),windex(nmax,2),wcnt,ktop,j
	logical new
	ive_ptr wvar,rho,rho_calc,omega
c
      massflx_top = 0
      massflx_top = getmem(ntime)
      if (massflx_top.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field massflx.'
	  return
	  end if
c
c  Convert phycical points to index (could use phys_2_index)
c
	indx(1,1) = nint(x1(1)/(dx/1000.)) + 1
	indx(2,1) = nint(x2(1)/(dx/1000.)) + 1
	indx(3,1) = nint(x3(1)/(dx/1000.)) + 1
	indx(4,1) = nint(x4(1)/(dx/1000.)) + 1
	indx(1,2) = nint(x1(2)/(dy/1000.)) + 1
	indx(2,2) = nint(x2(2)/(dy/1000.)) + 1
	indx(3,2) = nint(x3(2)/(dy/1000.)) + 1
	indx(4,2) = nint(x4(2)/(dy/1000.)) + 1
c
	write(6,*) 'Calculating the vertical mass flux between'
	write(6,*) 'x1 = ',x1(1),'y1 = ',x1(2)
	write(6,*) 'x2 = ',x2(1),'y2 = ',x2(2)
	write(6,*) 'x3 = ',x3(1),'y3 = ',x3(2)
	write(6,*) 'x4 = ',x4(1),'y4 = ',x4(2)
c
	call wpts(indx,windex,wcnt,nx,ny,nmax)
c
	wvar = omega(ndims,dims,stag,min,max,miss,dunits,
     >              disunits,dim_names)
c
      rho = 0
	rho = rho_calc(ndims,dims,stag,min,max,miss,dunits,
     >               disunits,dim_names)
c
	call var_at_w(%val(wvar),%val(rho),%val(sgz),%val(wgz),
     >			  0,0,0,nx,ny,nz,ntime)
c
	call freemem(rho)
c
	ktop = k2
	if(top) then
	  ktop = ktop+1	
	endif
c
	call massflx_top_calc(windex,ktop,%val(massflx_top),
     >     %val(dzdq),%val(wvar),dx,dy,wcnt,nx,ny,nz,ntime,nmax,top)	   
c
	call freemem(wvar)
c
	dunits = 'kg/s'
	disunits = 'kg/s'
c
	do j = 1,3
	  dims(j) = 1
	  stag(j) = 0.
	  min(j) = 0.
	  max(j) = 0.
	  dim_names(j) = 'one'
	enddo
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c
	subroutine massflx_top_calc(windex,ktop,massflx,
     >           dzdq,w,dx,dy,wcnt,nx,ny,nz,nt,nmax,top)	   
	implicit none
c	
      integer nx,ny,nz,nt,nmax
	integer windex(nmax,2),ktop,wcnt
	real massflx(nt),w(nx+1,ny+1,nz+1,nt),dx,dy
	real dzdq(nx+1,ny+1,3)
	logical top
c
	integer n,i,iloc,jloc,sgnflx
c
	sgnflx = 1
	if(top) sgnflx = -1
c
	do n = 1,nt
	  massflx(n) = 0.0
	  do i=1,wcnt
  	    iloc = windex(i,1)
	    jloc = windex(i,2)
	    massflx(n) = massflx(n) + 
     >                dzdq(iloc,jloc,1)*sgnflx*dx*dy*w(iloc,jloc,ktop,n)
	  enddo
	enddo
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c
	subroutine wpts(indx,windex,wcnt,nx,ny,nmax)
	implicit none
c
	integer nmax
	integer indx(4,2),windex(nmax,2),wcnt,nx,ny
c
	integer tcnt(4),i,j,k,tindex1(nmax,2),tindex2(nmax,2)
	integer tindex3(nmax,2),tindex4(nmax,2)
	integer indx1(2),indx2(2),indx3(2),indx4(2),sgnni,sgnnj
	
	logical isin(4)
c
	do i=1,2
	  indx1(i) = indx(1,i)
	  indx2(i) = indx(2,i)
	  indx3(i) = indx(3,i)
	  indx4(i) = indx(4,i)
	enddo
	
c
	call tline(indx1,indx2,tindex1,tcnt(1),nmax,sgnni,sgnnj)
	call tline(indx2,indx3,tindex2,tcnt(2),nmax,sgnni,sgnnj)
	call tline(indx3,indx4,tindex3,tcnt(3),nmax,sgnni,sgnnj)
	call tline(indx4,indx1,tindex4,tcnt(4),nmax,sgnni,sgnnj)
c
	wcnt = 0
	do i = 1,nx
	  do j = 1,ny
c
	    do k=1,4
	      isin(k) = .false.
	    enddo
c
	    do k = 1,tcnt(1)
	      if(i.ge.tindex1(k,1).and.j.ge.tindex1(k,2)) 
     >	  isin(1) = .true.
	    enddo	
c
	    do k = 1,tcnt(3)
	      if(i.le.tindex3(k,1).and.j.le.tindex3(k,2)) 
     >	  isin(3) = .true.
	    enddo	
c
	    do k = 1,tcnt(2)
	      if(j.ge.tindex2(k,2).and.i.le.tindex2(k,1)) 
     >  	  isin(2) = .true.
	    enddo	
c
	    do k = 1,tcnt(4)
	      if(j.le.tindex4(k,2).and.i.ge.tindex4(k,1)) 
     >	  isin(4) = .true.
	    enddo	
c
	    do k = 1,4
	      if(.not.isin(k)) GO TO 40
	    enddo
c
	    wcnt = wcnt + 1
	    windex(wcnt,1) = i
	    windex(wcnt,2) = j
c	    print *,windex(wcnt,1),windex(wcnt,2),wcnt
c
40	    CONTINUE
c
	  enddo
	enddo
c
	return
	end
