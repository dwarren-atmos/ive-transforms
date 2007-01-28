
#include "pointer.icl"

      ive_ptr function urot(theta,ndims,dims,stag,min,max,miss,
     >                        dunits,disunits,dim_names)
      implicit none
c
#include "default.icl"
c
c  Rotates the u and v velocity by counter-clockwise by angle theta
c
      integer ndims,dims(4)
      real theta,stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
c
      integer iu,iv,ijk,i,nt
	real r_val
	logical new
      ive_ptr uvar,vvar,vvartmp
c
      urot = 0
	vvar = 0
c
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
	nt = dims(4)
c
	vvar = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (vvar.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field vvar.' 
        return
      end if
c
c  Initialize vvar
c
c	do i = 1,dims(1)*dims(2)*dims(3)*dims(4)
	do i = 1,dims(1)*dims(2)*dims(3)*dims(4)
	  call s_val(%val(vvar),i,1)
	enddo
c
c  V at U points
c
 	iu = 0;
	iv = 1;
	call varatu(%val(vvar),%val(vvartmp),
     >		  nx,ny,nz,nt,iu,iv)
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
      iu = 1;
	iv = 0;
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
      implicit none
c
#include "default.icl"
c
c  Rotates the u and v velocity by counter-clockwise by angle theta
c
      integer ndims,dims(4)
      real theta,stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
c
      integer iu,iv,ijk,i,nt
	real r_val
	logical new
      ive_ptr uvar,vvar,vvartmp
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
	nt = dims(4)

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
      iu = 0;
	iv = 1;
	call varatu(%val(vvar),%val(vvartmp),
     >		  nx,ny,nz,nt,iu,iv)
c
	vrot = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (vrot.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field vrot.' 
        return
      end if
c
c  Do Rotation Calculation
c
      iu = 1;
	iv = 0;
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
	real uvar(nx+1-iu,ny+1-iv,nz,nt),vvar(nx+1-iu,ny+1-iv,nz,nt)
	real urot(nx+1-iu,ny+1-iv,nz,nt)
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
	  do j=1,ny+1-iv
	   do i=1,nx+1-iu
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
	subroutine varatu(var_u,var,nx,ny,nz,nt,iu,iv)
c
c  This subroutine calculates the variable var at the u velocity points
c  (at present var_u and var should be on the same vertical levels)
c  
	integer iu,iv,nx,ny,nz,nt
	real var(nx+1-iu,ny+1-iv,nz,nt),var_u(nx,ny+1,nz,nt)

	integer n,k,j,i
  	
	do n=1,nt
	 do k=1,nz
	  do i=1,nx+1
	   do j=1+iv,ny+1-iv
	     
	    var_u(i,j,k,n) = 
     >		 	.25*(var(i,j-iv,k,n) + var(i,j,k,n)
     >		+     var(i+1-iu,j-iv,k,n) + var(i+1-iu,j,k,n) )	
	   end do
	   if (iv.eq.1) then
	    var_u(i,1,k,n) = 
     >			   3/2*(.5*(var(i,1,k,n)  
     >		      +		var(i+1-iu,1,k,n))) 
     >                  -  1/2*(.5*(var(i,2,k,n)  
     >			+	      var(i+1-iu,2,k,n)))  
	    var_u(i,ny,k,n) = 
     >			   3/2*(.5*(var(i,ny,k,n)  
     >			+           var(i+1-iu,ny,k,n))) 
     >                  -  1/2*(.5*(var(i,ny-1,k,n)  
     >			+           var(i+1-iu,ny-1,k,n)))  
	   end if
	  end do
       end do
	end do

	return
	end
c
c ***************************************************** 
c
	subroutine varatm(var_m,var,nx,ny,nz,nt,iu,iv,iw)
c
c  This subroutine calculates the variable var at the mass points
c  
	integer iu,iv,iw,nx,ny,nz,nt
	real var(nx+1-iu,ny+1-iv,nz+iw,nt),var_m(nx+1,ny+1,nz,nt)

	integer n,k,j,i

  	
	do n=1,nt
	 do k=1,nz
	  do i=1+iu,nx+1-iu
	   do j=1+iv,ny+1-iv
	     
	    var_m(i,j,k,n) = 
     >		 	.5*(var(i-iu,j-iv,k+iw,n)+var(i,j,k,n))
	   end do
        end do
c
	    if (iu.eq.1) then
	      do j=1+iv,ny+1-iv
	        var_m(1,j,k,n) = 
     >			   3/2*var(1,j,k,n)  
     >                  -  1/2*var(2,j,k,n)  
	        var_m(nx+1,j,k,n) = 
     >			   3/2*var(nx  ,j,k,n)  
     >                  -  1/2*var(nx-1,j,k,n)  
            end do
	    end if

	   if (iv.eq.1) then
	     do i=1+iu,nx+1-iu
	      var_m(i,1,k,n) = 
     >			   3/2*var(i,1,k,n)  
     >                  -  1/2*var(i,2,k,n)  
	      var_m(i,ny+1,k,n) = 
     >			   3/2*var(i,ny  ,k,n)  
     >                  -  1/2*var(i,ny-1,k,n)  
           end do
	   end if
       end do
	end do

	return
	end
