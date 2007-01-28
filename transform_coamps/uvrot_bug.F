
#include "pointer.icl"

      ive_ptr function urot(theta,ndims,dims,stag,min,max,miss,
     >                        dunits,disunits,dim_names)
c
c  Rotates the u and v velocity by counter-clockwise by angle theta
c
      integer ndims,dims(4)
      real theta,stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
c
      integer iu,iv,ijk,nx,ny,nz,nt,i
	real r_val,s_val
	logical new
      ive_ptr uvar,vvar,vvartmp,getvar,getmem
c
      urot = 0
	vvar = 0
c
	iu = 1
	iv = 0
c
c  Get U and V variables
c
	new = .true.
      vvartmp = getvar('vtr',ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new)
      if (vvartmp.eq.0) then
        write(*,*) 'Error:  cannot get field v'
        return
      end if
c
	new = .true.
      uvar = getvar('utr',ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new)
      if (uvar.eq.0) then
        write(*,*) 'Error:  cannot get field u'
        return
      end if
c
	nx = dims(1) + 1
	ny = dims(2)
	nz = dims(3)
	nt = dims(4)

	 print *,nx,ny,nz,nt
c
	vvar = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (vvar.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field vvar.' 
        return
      end if
c
c  Initialize vvar
c
	do i = 1,dims(1)*dims(2)*dims(3)*dims(4)
	  call s_val(%val(vvar),i,1)
	enddo
c
c  V at U points
c
	call varatu(%val(vvar),%val(vvartmp),iu,iv,
     >		  nx,ny,nz,nt)
c
	call freemem(vvartmp)
c
	urot = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (urot.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field urot.' 
        return
      end if
c
c  Do Rotation Calculation
c
	call rot_calc('U',%val(urot),%val(uvar),%val(vvar),theta,
     >		   nx,ny,nz,nt,iu,iv)
c
	call freemem(vvar)
c
      return
      end
c
c ************************************************************
c
      ive_ptr function vrot(theta,ndims,dims,stag,min,max,miss,
     >                        dunits,disunits,dim_names)
c
c  Rotates the u and v velocity by counter-clockwise by angle theta
c
      integer ndims,dims(4)
      real theta,stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
c
      integer iu,iv,ijk,nx,ny,nz,nt,i
	real r_val
	logical new
      ive_ptr uvar,vvar,vvartmp,getvar,getmem
c
      vrot = 0
c
	iu = 1
	iv = 0
c
      vvartmp = getvar('v',ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new)
      if (vvartmp.eq.0) then
        write(*,*) 'Error:  cannot get field v'
        return
      end if
c
      uvar = getvar('u',ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new)
      if (uvar.eq.0) then
        write(*,*) 'Error:  cannot get field u'
        return
      end if
c
	nx = dims(1) + 1;
	ny = dims(2)
	nz = dims(3)
	nt = dims(4)
	 print *,nx,ny,nz,nt

c
	vvar = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (vvar.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field vvar.' 
        return
      end if
c
c  Initialize vvar
c
	do i = 1,dims(1)*dims(2)*dims(3)*dims(4)
	  call s_val(%val(vvar),i,1)
	enddo
c
c  V at U points
c
	call varatu(%val(vvar),%val(vvartmp),iu,iv,
     >		  nx,ny,nz,nt)
c
	vrot = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (vrot.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field vrot.' 
        return
      end if
c
c  Do Rotation Calculation
c
	call rot_calc('V',%val(vrot),%val(uvar),%val(vvar),theta,
     >		   nx,ny,nz,nt,iu,iv)
c
	call freemem(vvar)
c
      return
      end
c
c ******************************************************
c
	subroutine rot_calc(vari,urot,uvar,vvar,theta,
     >			  nx,ny,nz,nt,iu,iv)
c
c	Does rotation caluculation 
c
	integer nx,ny,nz,nt,iu,iv
	real theta
	real uvar(nx-iu,ny-iv,nz,nt),vvar(nx-iu,ny-iv,nz,nt)
	real urot(nx-iu,ny-iv,nz,nt)
	character(1) vari
c	
	integer i,j,k,n  
	real A,B,pi
c
	pi = 4.*atan(1.)
c
c  Convert theta to radians
c
	theta = theta*pi/180
c
	if (vari.eq.'U') then
	  A = cos(theta)
	  B =	sin(theta)
	else if (vari.eq.'V') then
	  A = -sin(theta)
	  B = cos(theta)
	end if
c
	do n=1,nt
	 do k=1,nz
	  do j=1,ny
	   do i=1,nx
	    urot(i,j,k,n) = uvar(i,j,k,n)*A +
     >    	              vvar(i,j,k,n)*B
	   end do
        end do
	 end do
	end do

	return
	end
c
c ***************************************************** 
c
	subroutine varatu(var_u,var,iu,iv,nx,ny,nz,nt)
c
c  This subroutine calculates the variable var at the u velocity points
c  (at present var_u and var should be on the same vertical levels)
c  
	integer iu,iv,nx,ny,nz,nt
	real var(nx-1+iu,ny-1+iv,nz,nt),var_u(nx-1,ny,nz,nt)

	integer n,k,j,i
  	
	do n=1,nt
	 do k=1,nz
	  do i=1,nx
	   do j=2-iv,ny-1+iv
	     
	    var_u(i,j,k,n) = 
     >		 	.25*(var(i,j-1+iv,k,n) + var(i+1-iu,j,k,n)
     >		+     var(i,j,k,n) + var(i+1-iu,j-1+iv,k,n) )	
	   end do
	   if (iv.eq.0) then
	    var_u(i,1,k,n) = 
     >			   3/2*(.5*(var(i,1,k,n)  
     >		      +		var(i+1-iu,1,k,n))) 
     >                  -  1/2*(.5*(var(i,2,k,n)  
     >			+	      var(i+1-iu,2,k,n)))  
	    var_u(i,ny,k,n) = 
     >			   3/2*(.5*(var(i,ny-1,k,n)  
     >			+           var(i+1-iu,ny-1,k,n))) 
     >                  -  1/2*(.5*(var(i,ny-2,k,n)  
     >			+           var(i+1-iu,ny-2,k,n)))  
	   end if
	  end do
       end do
	end do

	return
	end
c
c ***************************************************** 
c
	subroutine varatm(var_m,var,iu,iv,iw,nx,ny,nz,nt)
c
c  This subroutine calculates the variable var at the mass points
c  
	integer iu,iv,iw,nx,ny,nz,nt
	real var(nx-1+iu,ny-1+iv,nz+1-iw,nt),var_m(nx,ny,nz,nt)

	integer n,k,j,i
  	
	do n=1,nt
	 do k=1,nz
	  do i=2,nx-1
	   do j=2,ny-1
	     
	    var_m(i,j,k,n) = 
     >		 	.5*(var(i-1+iu,j-1+iv,k,n)+var(i,j,k+1-iw,n))
	   end do

	   if (iv.eq.0) then
	    var_m(i,1,k,n) = 
     >			   3/2*var(i,1,k,n)  
     >                  -  1/2*var(i,2,k,n)  
	    var_m(i,ny,k,n) = 
     >			   3/2*var(i,ny-1,k,n)  
     >                  -  1/2*var(i,ny-2,k,n)  
	   end if
c
	  end do
c
	  if (iu.eq.0) then
	    do j=2,ny-1
	      var_m(1,j,k,n) = 
     >			   3/2*var(1,j,k,n)  
     >                  -  1/2*var(2,j,k,n)  
	      var_m(1,j,k,n) = 
     >			   3/2*var(nx-1,j,k,n)  
     >                  -  1/2*var(nx-2,j,k,n)  
           end do
	  end if
	    
       end do
	end do

	return
	end
