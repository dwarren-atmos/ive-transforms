      subroutine var_at_t(var_t,var,iu,iv,iw,nx,ny,nz,nt)
c
c  This subroutine calculates the variable var at the thermo points
c  (at present var_t and var should be on the same vertical levels)
c  
      integer iu,iv,iw,nx,ny,nz,nt
	  real var(nx+1-iu,ny+1-iv,nz+iw,nt)
      real var_t(nx+1,ny+1,nz,nt),tmpvar
c
      integer n,k,j,i
c
      if (iv.eq.1.and.iu.eq.1) then
        write(6,*) 'Error:  You are not on the mesh'
        return
      endif
c
      do n=1,nt; do k=1,nz; 
        do i=1+iu,nx+1-iu; do j=1+iv,ny+1-iv;
          tmpvar = .5*(var(i-iu,j-iv,k,n) + var(i,j,k,n))
          var_t(i,j,k,n) = var_t(i,j,k,n) * tmpvar
      end do; end do;
c
      if(iv.eq.1) then
        do i = 1,nx+1
            var_t(i,1,k,n) = 1.5*var_t(i,2,k,n) - 
     &           0.5*var_t(i,3,k,n)
            var_t(i,ny+1,k,n) = 1.5*var_t(i,ny,k,n) - 
     &           0.5*var_t(i,ny-1,k,n)
        enddo
c
      else if(iu.eq.1) then
        do j = 1,ny+1
            var_t(1,j,k,n) = 1.5*var_t(2,j,k,n) - 
     &          0.5*var_t(3,j,k,n)
                 var_t(nx+1,j,k,n) = 1.5*var_t(nx,j,k,n) - 
     &          0.5*var_t(nx-1,j,k,n)
        enddo
      endif
      end do; end do;
c
      return
      end
c
      subroutine var_at_u(var_u,var,iu,iv,iw,nx,ny,nz,nt)
c
c  This subroutine calculates the variable var at the u velocity points
c  (at present var_u and var should be on the same vertical levels)
c  
      integer iu,iv,iw,nx,ny,nz,nt
      real var(nx+1-iu,ny+1-iv,nz+iw,nt)
      real var_u(nx,ny+1,nz,nt),tmpvar
c
      integer n,k,j,i
c 
      do n=1,nt ; do k=1,nz ; do i=1,nx ; do j=1+iv,ny+1-iv
        tmpvar = 
     &          .25*(var(i,j,k,n) + var(i,j-iv,k,n)
     &          + var(i+1-iu,j,k,n) + var(i+1-iu,j-iv,k,n))
        var_u(i,j,k,n) = var_u(i,j,k,n) * tmpvar
      end do
	   if (iv.eq.1) then
	    var_u(i,1,k,n) = var_u(i,1,k,n) *
     >			   (3/2*(.5*(var(i,1,k,n)  
     >		      +		var(i+1,1,k,n))) 
     >                  -  1/2*(.5*(var(i,2,k,n)  
     >			+	      var(i+1,2,k,n))))  
	    var_u(i,ny+1,k,n) = var_u(i,ny+1,k,n) *
     >			   (3/2*(.5*(var(i,ny,k,n)  
     >			+           var(i+1,ny,k,n))) 
     >                  -  1/2*(.5*(var(i,ny-1,k,n)  
     >			+           var(i+1,ny-1,k,n))))  
	   end if
	  end do
       end do
	end do
c
	return
	end
c
	subroutine var_at_v(var_v,var,iu,iv,iw,nx,ny,nz,nt)
c
c  This subroutine calculates the variable var at the v velocity points
c  (at present var_v and var should be on the same vertical levels)
c  
	integer iu,iv,iw,nx,ny,nz,nt
	real var(nx+1-iu,ny+1-iv,nz+iw,nt)
      real var_v(nx+1,ny,nz,nt),tmpvar

	integer n,k,j,i
  	
	do n=1,nt
	 do k=1,nz
	  do j=1,ny
	   do i=1+iu,nx+1-iu
	    tmpvar = 
     >	     .25*(var(i,j,k,n) + var(i,j+1-iv,k,n)
     >	     +    var(i+iu,j,k,n) + var(i+iu,j+1-iv,k,n))	
	    var_v(i,j,k,n) = var_v(i,j,k,n) * tmpvar
	   end do
c
	   if (iu.eq.1) then
	    var_v(1,j,k,n) = var_v(1,j,k,n)*
     >			   (3/2*(.5*(var(1,j,k,n)  
     >		      +		var(1,j+1,k,n))) 
     >                  -  1/2*(.5*(var(2,j,k,n)  
     >			+	      var(2,j+1,k,n))))  
	    var_v(nx+1,j,k,n) = var_v(nx+1,j,k,n) *
     >			   (3/2*(.5*(var(nx,j,k,n)  
     >			+           var(nx,j+1,k,n))) 
     >                  -  1/2*(.5*(var(nx-1,j,k,n)  
     >			+           var(nx-1,j+1,k,n))))  
	   end if
	  end do
       end do
	end do

	return
	end
c
	subroutine var_at_w(var_w,var,sgz,wgz,iu,iv,iw,nx,ny,nz,nt)
c
c  This subroutine calculates the variable var at the w velocity points
c  
	integer iu,iv,iw,nx,ny,nz,nt
	real var(nx+1-iu,ny+1-iv,nz+iw,nt)
      real var_w(nx+1,ny+1,nz+1,nt)
      real var_t(nx+1,ny+1,nz,nt)
	real sgz(nz),wgz(nz+1)

	integer n,k,j,i
	real zfrac,tmpvar
c
	if (iw.eq.0) then	
c
	do n=1,nt
	 do i=1,nx+1
	  do j=1,ny+1
	   do k=1,nz
	     var_t(i,j,k,n) = 1.0
         enddo
        enddo
       enddo
	enddo
c
	call var_at_t(var_t,var,iu,iv,iw,nx,ny,nz,nt)
c
	do n=1,nt
	 do i=1,nx+1
	  do j=1,ny+1
	   do k=2,nz
	     zfrac = (wgz(k)-sgz(k-1))/(sgz(k)-sgz(k-1))
	     tmpvar = var_t(i,j,k-1,n) + 
     >		 (var_t(i,j,k,n) - var_t(i,j,k-1,n))*zfrac
	     var_w(i,j,k,n) = var_w(i,j,k,n)*tmpvar
	   end do
c
c  Take care of top and bottom boundries
c
	    var_w(i,j,1,n) = var_w(i,j,1,n)*
     >	               (var_w(i,j,2,n) +
     >  	               (var_w(i,j,3,n)-var_w(i,j,2,n)) *
     >			   (wgz(1)-sgz(1))/(sgz(2)-sgz(1)))
c
	    var_w(i,j,nz+1,n) = var_w(i,j,nz+1,n) *
     >	               (var_w(i,j,nz-1,n)  +
     >  	               (var_w(i,j,nz,n)-var_w(i,j,nz-1,n)) *
     >			   (wgz(nz+1)-sgz(nz))/(sgz(nz)-sgz(nz-1)))
	  end do
	 end do
	end do
c
	 else
	  do n = 1,nt
	   do k = 1,nz+1
	    do j = 1,ny+1
	     do i = 1,nx+1
	      var_w(i,j,k,n) = var_w(i,j,k,n)*var(i,j,k,n)
	     enddo
    	    enddo
  	   enddo
       enddo		  
      end if
c
	return
	end
