
#include "pointer.icl"

      ive_ptr function area(varnm,i1,i2,j1,j2,k1,thres,
     >                       ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names)

#include "default.icl"

      integer i1,i2,j1,j2,k1,ndims,dims(4)
      real thres,stag(4),min(4),max(4),miss
      character*(*) varnm,dim_names(4),dunits,disunits

      integer strbeg,strend,k,nt
      logical new
      ive_ptr var,getmem

      area = 0

 79   FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 90   FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)

      write(6,*)
      write(6,*) 'Calculating the area of ',
     >                  varnm(strbeg(varnm):strend(varnm)),
     >	     ' less than ',thres,' at the lowest grid point'
      write(6,*) 'over the domain defined by:  '
      write(6,79) plmin(1)+(i1-1)*dx, plmin(1)+(i2-1)*dx
      write(6,90) plmin(2)+(j1-1)*dy, plmin(2)+(j2-1)*dy
      write(6,*)

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot compute field ',
     >              varnm( strbeg(varnm):strend(varnm) )
        return
      end if

	area = getmem(dims(4));
	if (area.eq.0) then
	  write(*,*) 'Error:  cannot allocate space for field area.'
	  return
	end if

      if ( (j1.lt.1).or.(j2.gt.dims(2)).or.  
     >     (i1.lt.1).or.(i2.gt.dims(1))) then
        write(6,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries.'
        return
      end if

      call area_calc(%val(area),%val(var),thres,i1,i2,j1,j2,k1,
     >	         dims(1),dims(2),dims(3),dims(4),dx,dy)

	do k=1,3
	  dims(k) = 1
	  stag(k) = 0.
	  min(k) = 0.
	  max(k) = 0.
	  dim_names(k) = 'one'
	enddo

      dunits = 'm*m'
      disunits = 'm*m'


      return
      end


	subroutine area_calc(area,var,thres,i1,i2,j1,j2,k1,
     >			   nx,ny,nz,nt,dx,dy)
c
c  Calculate area on lowest grid point that var < thres.
c  returns area.
c
	integer i1,i2,j1,j2,k1,nx,ny,nz,nt
	real area(nt),var(nx,ny,nz,nt),thres,dy,dx

	integer i,j,n
	real count

	do n=1,nt
	  count = 0.
	  do i=i1,i2
  	    do j=j1,j2
	      if(var(i,j,k1,n).le.thres) count = count + 1.
          end do
        end do

	  area(n) = count*dx*dy
 	end do

	return 
	end


