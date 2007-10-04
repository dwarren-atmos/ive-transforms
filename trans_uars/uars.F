c     $Id: uars.F,v 1.7 1997-12-18 20:35:10 harry Exp $
c     $Log: uars.F,v $
c     Revision 1.7  1997-12-18 20:35:10  harry
c     Add field_delta.
c
c     Revision 1.6  1997/02/12 17:20:19  harry
c     Change pointer to ive_ptr for F90.
c
c     Revision 1.5  1996/02/22 02:25:44  harry
c     Forgot new parameter coord_dep in new_file.
c
c     Revision 1.4  1995/12/15 23:54:12  harry
c     Add call to derivative in calc_field.  Fix Makefile for Solaris.
c
c Revision 1.3  1995/10/05  22:16:34  harry
c Change uars to reflect changes in new_file.
c
c Revision 1.2  1995/09/05  20:36:28  harry
c Modifications for new slicer.
c
c
c     This is a Fortran version of the default transforms for IVE.
c
c     -----------------------------------------------------------------
c
      subroutine new_file(ncid, exact_coord, n, coord_dep)
c
c     new_file : This routine is called every time a file is read in.
c     The routine can be used to read other information from the open
c     netCDF file.  It also can provide an array of times that limit the
c     values the time slider can be set to.
c     
c     Arguments:
c     ncid	  integer   input  Id of currently opened netCDF file.
c     exact_coord ive_ptr   output Pointer to array of discrete coord values
c				   that the window sliders should use.
c     n		  integer   output Number of values in exact_coord.  If this
c				   is 0 => ignore values in exact_coord
c				   => no restrictions on window slider.
c
c     coord_dep   integer   output States the functional dependencies
c     				   between the coordinates, 1 => dependent,
c     				   0 => independent e.g. if z is a function of
c				   x, y and z, then coord_dep(1,3),
c				   coord_dep(2,3) and coord_dep(3,3)
c				   will be one.  A coordinate is always 
c				   assumed to be dependent on itself.
#include "default.icl"
      integer ncid, n(4), coord_dep(4, 4)
      ive_ptr exact_coord(4)
c
c     common block cdf_dimvar is set up
c     by IVE to contain the values of a variable with the same name as the
c     dimension of the largest variable in the netCDF file.  If the
c     size of the domain is set by attributes, or if a dimension
c     variable does not exist, size will be zero, and data will be NULL.
c
      ive_ptr data(4)
      integer size(4)
      common/cdf_dimvar/data, size


      integer i, j

      do j = 1, 4
         do i = 1, 4
            if (i .eq. j) then
               coord_dep(i,j) = 1
            else
               coord_dep(i,j) = 0
            endif
         enddo
      enddo
      do i = 1, 4
         exact_coord(i) = data(i)
         n(i) = size(i)
      enddo
      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine new_field(name, field, inmax1, inmax2, inmax3, inmax4)
c     
c     new_field : This routine is called every time a new field is set.
c     The routine can be used to store necessary information to perform
c     the other transformations.  It also can use used to reorder the
c     data array.
c     
c     Note: before IVE calls new_field it has already called
c     getvar to read in the data for the new field.  As a
c     result of this call, IVE will have also read in
c     ndims, dims, stagger, phmin, phmax, missing, data_units
c     and display_units---provided this information appears
c     in the netCDF.  The transform routines can access this
c     data by calling getrarr, getlvar or getaarr.
c
c     Arguments:
c     name	character input  The name of the current field.
c     field	real	  input  Values of "field".
c     inmax1	integer	  input  Size of dimensions in "Fortran"
c     - inmax4			 order (i.e., A(inmax1,...inmax4)
c				 will hold the data with the first
c				 index varying most rapidly)
c

#include "default.icl"

      character*(*) name
      integer       inmax1, inmax2, inmax3, inmax4
      real          field(inmax1, inmax2, inmax3, inmax4)

      character*128 dnames(4)
      character*81  data_units, display_units, dim_names(4)
      real          missing, stagger(4), start, r_val
      integer       i, j, min(4), max(4), dims(4), nd, error, ndims
      logical       flag

      real logP(100)
      common /TRANS_COM/ logP

c
c     Save the physical and index limits in the common block in
c     default.icl for use in later calculations by other transform
c     routines.
c
      call getrarr("phmin", phmin, 4, error)
      call getrarr("phmax", phmax, 4, error)
      inmax(1) = inmax1
      inmax(2) = inmax2
      inmax(3) = inmax3
      inmax(4) = inmax4
c
c     Get the names of the dimensions of this field.
c
      call getaarr("data_dim_names", dnames, 4, error)
      call getivar("num_dims", ndims, error)
c
c     Set up P for logarithmic interpolation.
c
      if (inmax(3) .ne. 1) then
         phmin(3) = log10(phmin(3))
         phmax(3) = log10(phmax(3))
      endif
c
c     Get values (if they exist) in the 1D array associated with each
c     dimension to see which dimensions are evenly spaced.
c     delta is zero if they are not.  "delta" is stored in the common
c     block in "default.icl"
c
      do i = 1, ndims
c
c     Check that the data array depth exceeds one along this coordinate
c
         if (inmax(i).ne.1) then
c
c     delta is the distance between grid point in this dimension if
c     the points are evenly spaced.
c
            delta(i) = (phmax(i)-phmin(i))/(inmax(i)-1)
c
c     Get the memory address of the values of the variable with
c     the same name as the dimension.
c     If such a variable does not exist, coord(i) == 0.
c
            flag = .false.
            coord(i) = getvar(dnames(i), nd, dims, stagger, min, max,
     &              missing, data_units, display_units, dim_names, flag)
            if (coord(i) .ne. 0) then
               if (i .eq. 3) then
                  do j = 1, inmax(3)
                     logP(j) = log10(r_val(%val(coord(3)),j))
                  enddo
                  coord(3) = %loc(logP)
               endif
c
c     Check to see if the dimension is evenly spaced by getting the first
c     value and then comparing the other values to the value calculated
c     using the delta computed above.  If the dimension is not evenly
c     spaced, set delta to 0.  The function "r_val" is used to 
c     get the actual data  values from the memory pointer since Fortran
c     cannot directly deal with pointers.
c
               start = r_val(%val(coord(i)), 1)
               do j = 2, inmax(i)
                  if (start+(j-1)*delta(i) .ne.
     &                 r_val(%val(coord(i)), j)) then
                     delta(i) = 0.
                     goto 10
                  endif
               enddo
 10            continue
            endif
         endif
      enddo
      call setrarr("field_delta", delta, 4, error)
      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine index_2_phys(phys, index, flag, ndims, npts)

c     
c     index_2_phys : This routine translates array indices into physical
c     coordinates.
c     
c     Arguments:
c     phys	real	output	Array of physical coordinates.
c     index	real	input	Array of indicies.
c     flag	integer	input	Flags that indicate which coordinates
c				to compute.  If iflag(i) = 1 index(i,j)
c				will be converted to phys(i,j),
c				otherwise the value of phys(i,j)
c				is not changed
c     ndims	integer	input	Number of dimensions in each point.
c     npts	integer	input	Number of points to convert.
c
#include "default.icl"

      integer npts, ndims, flag(ndims)
      real    phys(ndims, npts), index(ndims, npts)

      integer i, j, k
      real    r_val	! function

      do i = 1, npts
         do j = 1, ndims
            if (flag(j) .eq. 1) then
c
c     Trivial case where dimension does not exist (depth of array is 
c     unity along this dimension)
c
               if (inmax(j) .eq. 1) then
                  phys(j,i) = phmin(j)
c
c     Case where interpolation must be performed between the grid
c     point locations given in a 1D "dimension variable" array.
c
               else
                  if (delta(j) .eq. 0.) then
                     k = index(j,i)
                     if (k .lt. inmax(j)) then
                        phys(j,i) = r_val(%val(coord(j)),k) +
     &                       (index(j,i)-k)*(r_val(%val(coord(j)),k+1)
     &                       -r_val(%val(coord(j)),k))
                     else
                        phys(j,i) = phmax(j)
                     endif
c
c     Case of evenly spaced physical coord.
c
                  else
                     phys(j,i) = phmin(j) + (index(j,i)-1)*delta(j)
                  endif
c
c     Convert pressure to lognormal.
c
                  if (j .eq. 3) then
                     phys(j,i) = 10**(phys(j,i))
                  endif
               endif
            endif
         enddo
      enddo
      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine phys_2_index(phys, index, flag, ndims, npts)
c     
c     phys_2_index : This routine translates physical coordinates into 
c     array indices.
c     
c     Arguments:
c     phys	real	input	Array of physical coordinates.
c     index	real	output	Array of indicies.  A negative index =>
c				this point gets the missing data value
c				(e.g. is under the terrain).
c     flag	integer	input	Flags that indicate which coordinates
c				to compute.  If iflag(i) = 1 phys(i,j)
c				will be converted to index(i,j),
c				otherwise the value of index(i,j)
c				is not changed
c     ndims	integer	input	Number of dimensions in each point.
c     npts	integer	input	Number of points to convert.
c
#include "default.icl"

      integer npts, ndims, flag(ndims)
      real    p, phys(ndims, npts), index(ndims, npts)

      integer i, j, k
      real    r_val	! function
      

      do i = 1, npts
         do j = 1, ndims
            if (flag(j) .eq. 1) then
c
c     Interpolate log P.
c
               if (j .eq. 3) then
                  p = log10(phys(j,i))
               else
                  p = phys(j,i)
               endif
c
c     Trivial case where dimension does not exist (depth of array is 
c     unity along this dimension)
c
               if(inmax(j) .eq. 1) then
                  index(j,i) = 1.
c
c     Case where interpolation must be performed between the grid
c     point locations given in a 1D "dimension variable" array.
c
               else if (delta(j) .eq. 0.) then
c
c     Physical coord. increasing with index.
c
                  if (phmin(j) .lt. phmax(j)) then
                     if (p .le. phmin(j)) then
                        index(j,i) = 1.
                     else if (p .ge. phmax(j)) then
                        index(j,i) = inmax(j)
                     else
                        do k = 2, inmax(j)
                           if (p.le.r_val(%val(coord(j)),k))
     &                          goto 10
                        enddo
 10                     index(j,i) = float(k) - (p-
     &                       r_val(%val(coord(j)),k))/
     &                       (r_val(%val(coord(j)),k-1)-
     &                       r_val(%val(coord(j)),k))
                     endif
c
c     Physical coord. decreasing with index.
c
                  else
                     if (p .ge. phmin(j)) then
                        index(j,i) = 1.
                     else if (p .le. phmax(j)) then
                        index(j,i) = inmax(j)
                     else
                        do k = 2, inmax(j)
                           if (p.ge.r_val(%val(coord(j)),k))
     &                          goto 20
                        enddo
 20                     index(j,i) = float(k) - (p-
     &                       r_val(%val(coord(j)),k))/
     &                       (r_val(%val(coord(j)),k-1)-
     &                       r_val(%val(coord(j)),k))
                     endif
                  endif
c
c     Case of evenly spaced physical coord.
c
               else
                  index(j,i) = (p-phmin(j))/delta(j) + 1
               endif
            endif
         enddo
      enddo
      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine horiz_ter(topo, nxw, nyw, stagi, stagj, zero, error)
c
c     horiz_ter : This routine fills an array with terrain heights in the
c     windowed domain for horizontal cross sections.
c     
c     Arguments:
c     topo	real		Array containing the terrain heights (output).
c     nxw	integer		Number of x points.
c     nyw	integer		Number of y points.
c     stagi	real		Grid staggering along x dimension.
c     stagj	real		Grid staggering along y dimension.
c     zero	real		Minimum terrain height.
c     error	integer		Error flag 0 => no errors (output).
c
      integer nxw, nyw, error, i, j
      real    topo(nxw, nyw), stagi(nxw), stagj(nyw), zero

      error = 0
      do j = 1, nyw
         do i = 1, nxw
            topo(i, j) = zero
         enddo
      enddo
      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine vert_ter(xter, yter, topo_ht, nter, pt1, pt2, zero)
c     
c     vert_ter : This routine fills an array with terrain heights in the
c     windowed domain for vertical cross sections.
c
c     Arguments:
c	xter	ive_ptr	output	Pointer to array of x-locations for topo_ht
c	yter	ive_ptr	output	Pointer to array of y-locations for topo_ht
c	topo_ht ive_ptr	output	Pointer to array containing 
c				the terrain heights.
c	nter	integer	output	Number of points in topo_ht (output)
c	pt1	real	input	Minimum physical coordinates for slice.
c	pt2	real	input	Maximum physical coordinates for slice.
c	zero	real	input	Minimum terrain height.
c
      integer nter
      real    pt1(4), pt2(4), zero
      ive_ptr xter, yter, topo_ht, getmem

      nter = 2
      xter = getmem(2)
      yter = getmem(2)
      topo_ht = getmem(2)
      call s_val(%val(xter), 1, pt1(1))
      call s_val(%val(xter), 2, pt2(1))
      call s_val(%val(yter), 1, pt1(2))
      call s_val(%val(yter), 2, pt2(2))
      call s_val(%val(topo_ht), 1, zero)
      call s_val(%val(topo_ht), 2, zero)
      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine phys_2_lonlat(x, y, lon, lat, npts)
c     
c     phys_2_lonlat : This routine converts physical coordinates to longitude-
c     latitude coordinates.
c     
c     Arguments:
c     x		real		Physical x coordinate array.
c     y		real		Physical y coordinate array.
c     lon	real		Longitude array (output).
c     lat	real		Latitude array (output).
c     npts	integer		Number of points to convert.
c
      integer npts
      real    x(npts), y(npts), lon(npts), lat(npts)

      integer i

      do i = 1, npts
         lon(i) = x(i)
         lat(i) = y(i)
      enddo
      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine lonlat_2_phys(x, y, lon, lat, npts)
c     
c     lonlat_2_phys : This routine converts longitude-latitude coordinates
c     to physical coordinates.
c     
c     Arguments:
c     x		real		Physical x coordinate array (output).
c     y		real		Physical y coordinate array (output).
c     lon	real		Longitude array (output)
c     lat	real		Latitude array (output)
c     npts	integer		Number of points to convert
c
      integer npts
      real    x(npts), y(npts), lon(npts), lat(npts)

      integer i

      do i = 1, npts
         x(i) = lon(i)
         y(i) = lat(i)
      enddo
      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine default_map(proj, plon, plat, rota, limit,
     &     plm1, plm2, plm3, plm4, exact_fit)
c
c     default_map : This routine returns the map settings when the
c     default map is specified.
c
c     Arguments:
c     proj	character output EZMAP map projection.
c     plon	real	  output Projection longitude.
c     plat	real	  output Projection latitude.
c     rota	real	  output Projection rotation.
c     limit	character output limit and plm1-4 are the arguments for
c     plm1	real	  output   the MAPSET call which specifies the
c     plm2	real	  output   rectagular portion of the u/v plane
c     plm3	real	  output   to be drawn.
c     plm4	real	  output
c     exact_fit	integer	  output If 1 => the contour plot will fit exactly
c     				   in this map with no transform.  Thus IVE
c     				   will just call conpack within the map
c     				   window without calling phys_2_lonlat.
c
      character*(*) proj, limit
      integer       exact_fit
      real          plon, plat, rota, plm1(2), plm2(2), plm3(2), plm4(2)
c
c     The default is the current window in a Cylindrical Equidistant
c     Projection.
c
      proj = 'CE'
      plon = 0.
      plat = 0.
      rota = 0.
c
c     Setting limt to CL tells IVE to clip to the current window using
c     phys_2_lonlat.
c
      limit = 'CL'
      exact_fit = 0
      return
      end
c
c     -----------------------------------------------------------------
c
      ive_ptr function data_slicer_1d(endpt, ni, da)
c
c     data_slicer_1d_ : This routine returns the computational points at
c     which a 1d slice should be taken.
c
c     Arguments:
c     	endpt	real	input	2 endpoints of the line in physical space
c     	ni	integer	output	Number of points
c	da	integer	input	axis coordinate of average
c     	nda	integer	output	Number of points in the da-direction.
c
      real endpt(4,2)
      integer ni, da, nda
c
c     Let's do the obvious, just divide the slice evenly, using the
c     number of computational points represented by the slice.
c
      integer iflag(4), i, j, k, error
      real cendpt(4,2), dv(4), dda, cpmin(4), cmin, cmax, ceil, floor
      ive_ptr pslice, slice, getmem

      do i = 1, 4
         iflag(i) = 1
      enddo
      call phys_2_index(endpt, cendpt, iflag, 4, 2)
      ni = 0
      do i = 1, 4
         if (i .ne. da) ni = ni + abs(cendpt(i, 2) - cendpt(i, 1))
      enddo
      do i = 1, 4
         dv(i) = (endpt(i, 2) - endpt(i, 1))/(ni-1)
      enddo
c
c     No average case.
c
      if (da .eq. 0) then
         pslice = getmem(4*ni)
         slice = getmem(4*ni)
         do i = 1, ni
            do j = 1, 4
               call s_val(%val(pslice), (i-1)*4+j,
     &              endpt(j,1) + (i-1)*dv(j))
            enddo
         enddo
         j = ni
c
c     Average case
c
      else
         dv(da) = 0
         call getiarr("cpmin", cpmin, 4, error)
         if (cendpt(da, 1) .lt. 0) then
            cmin = cpmin(da)
         else
            cmin = cendpt(da, 1)
         endif
         if (cendpt(da, 2) .lt. 0) then
            cmax = cpmin(da)
         else
            cmax = cendpt(da, 1)
         endif
         if (cmax .eq. cmin) then
            nda = 1
         else
            nda = ceil(cmax) - floor(cmin) + 1
            dda = (endpt(da,2) - endpt(da,1))/(nda-1)
         endif
         pslice = getmem(ni*nda)
         slice = getmem(ni*nda)
         do i = 1, ni
            do k = 1, nda
               do j = 1, 4
                  if (j .ne. da .or. nda .eq. 1) then
                     call s_val(%val(pslice), 4*(nda*(i-1) + (k-1)) + j,
     &                    endpt(j,1) + (i-1)*dv(j))
                  else
                     call s_val(%val(pslice), 4*(nda*(i-1) + (k-1)) + j,
     &                    endpt(j,1) + (k-1)*dda)
                  endif
               enddo
            enddo
         enddo
         j = ni*nda
      endif
      call phys_2_index(pslice, slice, iflag, 4, j)
      call freemem(pslice)
      data_slicer_1d = slice
      return
      end
c
c     -----------------------------------------------------------------
c
      ive_ptr function data_slicer_2d(corner, ii, jj, ni, nj, da, nda)
c
c     data_slicer_2d_ : This routine returns the computational points at
c     which a 2d slice should be taken.
c
c     Arguments:
c     	corner	real	input	4 corner points of slice in physical space.
c     	ii	integer	input	axis coordinate of ni.
c     	jj	integer	input	axis coordinate of nj.
c     	ni	integer	output	Number of points in the x-direction.
c     	nj	integer	output	Number of points in the y-direction.
c     	da	integer	input	axis coordinate of average
c     	nda	integer	output	Number of points in the da-direction.
c
c     The corners will be:     2    3
c                           jj
c                              0    1
c                                ii
c
      real corner(4,4)
      integer ii, jj, ni, nj, da, nda
c
c     Call the internal routine that does 2d physical slicing for
c     vectors.
c
      ive_ptr pslicer2d

      data_slicer_2d = pslicer2d(corner, ii, jj, ni, nj, da, nda)
      return
      end
c
c     -----------------------------------------------------------------
c
      ive_ptr function calc_field(name, ndims, dims, stag, min, max,
     &     missing, data_units, data_display_units, dim_names)
c     
c     calc_field : This routine is used to calculate user-derived field.
c     The return value of the routine is a pointer to the field values.
c     NULL return => cannot calculate.
c     
c     Arguments:
c     name	character	The name of the field to derive.
c     ndims	integer		Number of dimensions in field (output).
c     dims	integer		Number of points in the field in Fortran
c				order (x, y, z, t) (output).
c     stag	real		Grid staggering per dimension (output).
c     min	real		Physical space minimum per dimension (output).
c     max	real		Physical space maximum per dimension (output).
c     missing	real		Missing data value, zero => none (output).
c     data_units
c		character	Units for field (output).
c     data_display_units
c		character	Units to use to display field (output).
c     dim_names character	Names of the dimensions in the field (output).
c
      integer ndims, dims(4)
      real    stag(4), min(4), max(4), missing
      character*(*) name, data_units, data_display_units, dim_names(4)

      ive_ptr derivative

      calc_field = derivative(name, ndims, dims, stag, min, max, 
     &     missing, data_units, data_display_units, dim_names)
      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine heading(where, which, line)
c
c     heading_ : This routine allows the user to change headings on plots
c
c     Arguments:
c     	where   integer		Wherefrom called
c     	which   integer		Which heading (first (1) or second (2) line)
c     	line    character	The line to be changed.
c
      integer where, which
      character*(256) line

      return
      end
c
c     -----------------------------------------------------------------
c
c
c     Local routines
c
      real function r_val(array, j)
c
c     This routine is used to return the value stored in an 1d real array.
c     Since Fortran cannot be used to directly access memory given
c     a pointer, this function can be used to dereference the pointer
c     and extract a value from the corresponding array as follows:
c
c     If p is the pointer to some data array X, you get the j th data 
c     value in X (i.e. X(j)) with the following call:
c
c     r_val(%val(p), j)
c
c     Arguments:
c     array	real	input	The array to extract the value from.
c     j		integer	input	The index of array to use.
c
      integer j
      real array(j)

      r_val = array(j)
      return
      end

      subroutine s_val(array, j, val)
c
c     This routine is used to store a value in a 1d real array.
c     Since Fortran cannot be used to directly access memory given
c     a pointer, this function can be used to dereference the pointer
c     and extract a value from the corresponding array as follows:
c
c     If p is the pointer to some data array X, you set the j th data
c     value in X (i.e. X(j)) to val with the following call:
c
c     call s_val(%val(p), j, val)
c
c     Arguments:
c     array	real	output	The array to store the value in.
c     j		integer	input	The index of array to use.
c     val	real	input	The value to set
c
      integer j
      real array(j), val

      array(j) = val
      return
      end

      real function ceil(x)
c
c     ceil returns the least integral value greater than or equal to x.
c
      real x, floor

      if (x .lt. 0) then
         ceil = -floor(-x)
      else
         ceil = int(x)
         if (ceil .ne. x) then
            ceil = ceil+1
         endif
      endif
      return
      end

      real function floor(x)
c
c     floor returns the greatest integral value less than or equal to x.
c
      real x, ceil

      if (x .lt. 0) then
         floor = -ceil(-x)
      else
         floor = int(x)
      endif
      end
