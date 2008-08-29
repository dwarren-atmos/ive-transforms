      integer(ive_ptr) function vert_grid(whch,ndims,dims,stag,min,max,
     >                                   miss,dunits,disunits,dim_names)
      use default_module
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
      integer(ive_ptr) ptr

      vert_grid = 0

      ptr = getvar('w',ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new) 
      if (ptr.eq.0) then
        write(*,*) 'Error:  could not find vertical velocity field'
        return
      end if
 
      dims(4) = 1
      stag(4) = 0.
      min(4) = 0.
      max(4) = 0.
      dunits = ' '
      disunits = ' '
      dim_names(4) = 'one'

      vert_grid = getmem(dims(1)*dims(2)*dims(3))
      if (vert_grid.eq.0) then
        write(*,*) 'Error:  could not allocate space for field grid'
        return
      end if

      do k = 1,dims(3)
      do j = 1,dims(2)
      do i = 1,dims(1)

        iptr = (k-1)*dims(1)*dims(2) + (j-1)*dims(1) + i
        if (whch.eq.'W') then
          call s_val(%val(vert_grid),iptr,real(k))
        else 
          call s_val(%val(vert_grid),iptr,real(k)-0.5)
        end if

      enddo
      enddo
      enddo

      return
      end

