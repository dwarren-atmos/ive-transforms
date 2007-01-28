#include "pointer.icl"

      ive_ptr function y_grid(whch,ndims,dims,stag,min,max,
     >			  miss,dunits,disunits,dim_names)
*
*  grid returns the index space vertical coordinate for each node
*    on the appropriate grid (w or thermodynamic).  As such, a plot of
*    grid with contour interval of unity shows the physical space
*    surfaces corresponding to the nodal surfaces in computational space.
*
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),whch

      integer i,j,k,iptr
      logical new
      ive_ptr ptr,getvar,getmem

      y_grid = 0

      ptr = getvar('ppp',ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new) 
      if (ptr.eq.0) then
        write(*,*) 'Error:  could not find pressure field'
        return
      end if
 
      dims(4) = 1
	dims(3) = 1
      stag(4) = 0.
      min(4) = 0.
      max(4) = 0.
      dunits = ' '
      disunits = ' '
      dim_names(4) = 'one'

      y_grid = getmem(dims(1)*dims(2))
      if (y_grid.eq.0) then
        write(*,*) 'Error:  could not allocate space for field grid'
        return
      end if

      do i = 1,dims(1)
      do j = 1,dims(2)

        iptr = (i-1)*dims(2) + j
        if (whch.eq.'mass') then
          call s_val(%val(y_grid),iptr,real(i))
        else 
          call s_val(%val(y_grid),iptr,real(i)-0.5)
        end if

      enddo
      enddo

      return
      end

