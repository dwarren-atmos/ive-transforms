      integer(ive_ptr) function x_grid(whch,ndims,dims,stag,min,max,
     &                                 miss,dunits,disunits,dim_names)

      use default_module
c
c  grid returns the index space vertical coordinate for each node
c    on the appropriate grid (w or thermodynamic).  As such, a plot of
c    grid with contour interval of unity shows the physical space
c    surfaces corresponding to the nodal surfaces in computational space.
c
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),whch

      integer i,j,k,iptr
      logical new
      integer(ive_ptr) ptr

      x_grid = 0

      ptr = getvar('ppp',ndims,dims,stag,min,max,miss,dunits,
     &                    disunits,dim_names,new) 
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

      x_grid = getmem(dims(1)*dims(2))
      if (x_grid.eq.0) then
        write(*,*) 'Error:  could not allocate space for field grid'
        return
      end if

      do i = 1,dims(1)
      do j = 1,dims(2)

        iptr = (i-1)*dims(1) + j
        if (whch.eq.'mass') then
          call s_val(%val(x_grid),iptr,real(j))
        else 
          call s_val(%val(x_grid),iptr,real(j)-0.5)
        end if

      enddo
      enddo

      return
      end



      integer(ive_ptr) function y_grid(whch,ndims,dims,stag,min,max,
     &                                 miss,dunits,disunits,dim_names)
      use default_module
c
c  grid returns the index space vertical coordinate for each node
c    on the appropriate grid (w or thermodynamic).  As such, a plot of
c    grid with contour interval of unity shows the physical space
c    surfaces corresponding to the nodal surfaces in computational space.
c
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),whch

      integer i,j,k,iptr
      logical new
      integer(ive_ptr) ptr

      y_grid = 0

      ptr = getvar('ppp',ndims,dims,stag,min,max,miss,dunits,
     &                    disunits,dim_names,new) 
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

