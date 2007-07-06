      integer(ive_ptr) function massflx_max(
     >		     ndims,dims,stag,min,max,miss,
     >		     dunits,disunits,dim_names)
      use default_module
	implicit none
c
	integer ndims,dims(4),indx,k1,k2
c	real x1(2),x2(2),x3(2),x4(2)
	real stag(4),min(4),max(4),miss
	character*(*) dunits,disunits,dim_names(4)
c
	integer nt,j
	logical new
	integer(ive_ptr) f1,f2,f3,f4,ftop
	integer(ive_ptr) massflx_lat,massflx_top
c
      massflx_max = 0
c
	f1 = 0
	f1 = getvar('F1',ndims,dims,stag,min,max,miss,dunits,
     >		disunits,dim_names,new)
c
	if(f1.eq.0) then
	  write(*,*) 'Error:  Cannont Get F1'
	  return
	endif
c
	f2 = 0
	f2 = getvar('F2',ndims,dims,stag,min,max,miss,dunits,
     >		disunits,dim_names,new)
c
	if(f2.eq.0) then
	  write(*,*) 'Error:  Cannont Get F2'
	  return
	endif
c
	f3 = 0
	f3 = getvar('F3',ndims,dims,stag,min,max,miss,dunits,
     >		disunits,dim_names,new)
c
	if(f3.eq.0) then
	  write(*,*) 'Error:  Cannont Get F3'
	  return
	endif
c
	f4 = 0
	f4 = getvar('F4',ndims,dims,stag,min,max,miss,dunits,
     >		disunits,dim_names,new)
c
	if(f4.eq.0) then
	  write(*,*) 'Error:  Cannont Get F4'
	  return
	endif
c
	ftop = 0
	ftop = getvar('FTOP',ndims,dims,stag,min,max,miss,dunits,
     >		  disunits,dim_names,new)
c
	if(ftop.eq.0) then
	  write(*,*) 'Error:  Cannont Get FTOP'
	  return
	endif
c
	nt = dims(4)
c
      massflx_max = getmem(nt)
      if (massflx_max.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field massflx.'
	  return
	end if
c
	call massflx_max_calc(%val(massflx_max),%val(f1),
     >			     %val(f2),%val(f3),%val(f4),%val(ftop),nt)
c	
	disunits = ' '
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
	subroutine massflx_max_calc(massflx_max,
     >			           f1,f2,f3,f4,ftop,nt)
	implicit none
c	
      integer nt
	real massflx_max(nt),f1(nt),f2(nt),f3(nt),f4(nt),ftop(nt)
c
	integer n
	real fmax
c
	do n = 1,nt
	   fmax = abs(f1(n)) 
	   if (abs(f2(n)).gt.fmax) fmax = abs(f2(n))
	   if (abs(f3(n)).gt.fmax) fmax = abs(f3(n))
	   if (abs(f4(n)).gt.fmax) fmax = abs(f4(n))
	   if (abs(ftop(n)).gt.fmax) fmax = abs(ftop(n))
	   massflx_max(n) = fmax	   
	enddo
c
	return
	end
