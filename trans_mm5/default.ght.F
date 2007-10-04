c     $Id: default.F,v 1.11 1997/07/25 21:00:24 harry Exp $
c     $Log: default.F,v $
c     Revision 1.11  1997/07/25 21:00:24  harry
c     Add new transforms: data_slicer_3d and file_coordinate.
c
c     Revision 1.10  1997/02/12 16:56:54  harry
c     Change "pointer" to "ive_ptr" to avoid F90 problems.
c
c     Revision 1.9  1996/06/04 22:52:33  harry
c     Add more comments for 4.0.0
c
c     Revision 1.8  1996/06/04 22:43:14  harry
c     Update for IVE 4.0.0 including more comments.
c
c     Revision 1.7  1996/02/22 02:25:04  harry
c     Forgot new parameter coord_dep in new_file.
c
c     Revision 1.6  1995/12/15 23:55:03  harry
c     Fix Fortran version for missing comma.
c
c     Revision 1.5  1995/12/15 23:49:01  harry
c     Add the call to derivative in calc_field.
c
c Revision 1.4  1995/10/05  22:04:39  harry
c Change transforms to reflect change in new_file transform.
c
c Revision 1.3  1995/09/05  20:35:06  harry
c Modifications for new slicer.
c
c Revision 1.2  1994/08/30  17:26:49  harry
c Add default.icl to new_file.
c
c Revision 1.1  1994/08/30  17:17:57  harry
c Moves to trans_default directory
c
c Revision 1.5  1994/02/24  19:23:04  harry
c Change dimension on dnames in calc_field.
c
c Revision 1.4  1994/02/18  00:31:30  harry
c Change exact_times from integer to pointer.
c
c Revision 1.3  1994/02/17  22:20:18  harry
c Eliminate unused variables, correct include files.
c
c Revision 1.2  1994/02/17  19:16:57  harry
c Change type of calc_field from integer to pointer.
c
c Revision 1.1  1994/02/17  19:09:03  warren
c moved getvar declaration to default.icl and declared it properly (pointer).
c This means default.f -> default.F
c
c Revision 1.6  1994/02/17  19:04:21  harry
c Change read_var calls to getvar.
c
c Revision 1.5  1994/01/07  18:53:04  harry
c Add default_map and make minor corrections.
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
      ive_ptr calc_field
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

c ***
      ive_ptr getmem,thsigma,tfsigma,ttime,txter,tptop, addvar
      ive_ptr tdcor,ttgrid,tps0,ttlp,tts0,txght,tdght,txwspd
      ive_ptr txlat,txlon,tpstx,tfcor,tdter,txtheta,txptot
      logical error,flag

      real r_val,stagger(4),min(4),max(4),missing
      real xstagger(4),xmin(4),xmax(4),xmissing
      real dstagger(4),dmin(4),dmax(4),dmissing
      real*8 slope,intercept
      integer ndims,dims(4),k,nzxght,nhsig
      integer xdims(4),ddims(4)
      character*80 data_units,display_units,dim_names(4),dom_units(4)
      character*80 ddim_names(4),xdim_names(4)

      data hsigma,fsigma,xter,dcor,ps0,ts0,tlp,ptop,xght,dght,tgrid,
     &    xlat,xlon,pstx,xtheta,xptot,xwspd 
     &    /0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./

c ***

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
c ***
c get some domain info from IVE
c
      call getaarr("domain_units",dom_units,4,error)
      if (error) then
        write(6,*) 'Error:  could not get domain_units in new_file.'
        stop
      end if
      call getrarr("grid_delta",delta,4,error)
      if (error) then
        write(6,*) 'Error:  could not get grid spacing in new_file.'
        stop
      end if

c
c set up grid point heigts at cross pts.
c
      if (xght.ne.0.) call freemem(xght)

      txght = getvar('xght',ndims,dims,stagger,min,max,missing,
     >           data_units,display_units,dim_names,flag)

      xdims = dims
      xstagger = stagger
      xmin = min
      xmax = max
      xmissing = missing
      xdim_names = dim_names

      if (txght.eq.0) then
        xght = 0.
        write(6,*) 'Error:  could not find grid point height at ',
     >                'cross pts.'
        stop
      else
        nzxght = dims(1)
        xght = getmem(dims(1)*dims(2)*dims(3))
        xght_min = 1e8
        xght_max = -1e8        
        call convert_units(data_units,dom_units(3),slope,intercept)
        print *,"dim ",dims(1),dims(2),dims(3),delta(1)
        do i=1,dims(1)*dims(2)*dims(3)
          call s_val(%val(xght),i,
     >         real(r_val(%val(txght),i)*slope+intercept))
          if ((r_val(%val(txght),i)*slope+intercept).lt.xght_min)
     >       xght_min = r_val(%val(txght),i)*slope+intercept
          if ((r_val(%val(txght),i)*slope+intercept).gt.xght_max)
     >       xght_max = r_val(%val(txght),i)*slope+intercept
        enddo
      endif
      xght_dims(1) = dims(1)
      xght_dims(2) = dims(2)
      xght_dims(3) = dims(3)
      xght_delta = delta(1)
      xght_mind(1) = min(1)
      xght_mind(2) = min(2)
c
c set up grid point heigts at dot pts.
c
      if (dght.ne.0.) call freemem(dght)

      tdght = getvar('dght',ndims,dims,stagger,min,max,missing,
     >           data_units,display_units,dim_names,flag)

      ddims = dims
      dstagger = stagger
      dmin = min
      dmax = max
      dmissing = missing
      ddim_names = dim_names

      if (tdght.eq.0) then
        dght = 0.
        write(6,*) 'Error:  could not find grid point height at ',
     >                'dot pts.'
        stop
      else
        nzdght = dims(1)
        dght = getmem(dims(1)*dims(2)*dims(3))
        dght_min = 1e8
        dght_max = -1e8        
        call convert_units(data_units,dom_units(3),slope,intercept)
        print *,"dim ",dims(1),dims(2),dims(3),delta(1)
        do i=1,dims(1)*dims(2)*dims(3)
          call s_val(%val(dght),i,
     >         real(r_val(%val(tdght),i)*slope+intercept))
          if ((r_val(%val(tdght),i)*slope+intercept).lt.dght_min)
     >       dght_min = r_val(%val(tdght),i)*slope+intercept
          if ((r_val(%val(tdght),i)*slope+intercept).gt.dght_max)
     >       dght_max = r_val(%val(tdght),i)*slope+intercept
        enddo
      endif
      dght_dims(1) = dims(1)
      dght_dims(2) = dims(2)
      dght_dims(3) = dims(3)
      dght_delta = delta(1)
      dght_mind(1) = min(1)
      dght_mind(2) = min(2)
c
c set up coordinate vector for cross pts
c
      if (hsigma.ne.0) call freemem(hsigma)

      thsigma = getvar('hsigma',ndims,dims,stagger,min,max,missing,
     >                 data_units,display_units,dim_names,flag)
      if (thsigma.eq.0) then
        hsigma = 0
        write(6,*) 'Error:  could not find vertical coordinate ',
     >                     'vector hsigma.'
        stop
      else
        nhsig = dims(1)
        hsigma = getmem(nhsig)
        call convert_units(data_units,dom_units(3),slope,intercept)
        do i = 1,nhsig
          call s_val(%val(hsigma),i,
     >          real(r_val(%val(thsigma),i)*slope+intercept))
          if ((r_val(%val(thsigma),i)*slope+intercept).lt.hsig_min)
     >       hsig_min = r_val(%val(thsigma),i)*slope+intercept
          if ((r_val(%val(thsigma),i)*slope+intercept).gt.hsig_max)
     >       hsig_max = r_val(%val(thsigma),i)*slope+intercept
        enddo
c        hsig_min = r_val(%val(hsigma),1)
c        hsig_max = r_val(%val(hsigma),nhsig)
      end if
      print *,hsig_min,hsig_max,nhsig
c
c set up terrain at cross pts; note that we use only the terrain at 
c   the pressure points when transforming coordinates
c
      if (xter.ne.0.) call freemem(xter)

      txter = getvar('xter',ndims,dims,stagger,min,max,missing,
     >           data_units,display_units,dim_names,flag)
      if (txter.eq.0) then
        xter = 0.
        write(6,*) 'Error:  could not find terrain at ',
     >                'thermodynamic pts.'
        stop

      else
        xter_dims(1) = dims(1)
        xter_dims(2) = dims(2)
        xter_min(1) = min(1)
        xter_min(2) = min(2)
        xter_max(1) = max(1)
        xter_max(2) = max(2)
        xter_delta(1) = delta(1)
        xter_delta(2) = delta(2)
        xter = getmem(dims(1)*dims(2))
        call convert_units(data_units,dom_units(3),slope,intercept)
        do i=1,dims(1)*dims(2)
          call s_val(%val(xter),i,
     >         real(r_val(%val(txter),i)*slope+intercept))
        enddo
        coord_dep(1,3) = 1
        coord_dep(2,3) = 1
      endif
      print *,dims(1),dims(2)
      print *,min(1),min(2)
      print *,max(1),max(2)
      print *,delta(1),delta(2)
c
c set up terrain at dot pts
c
      if (dter.ne.0.) call freemem(dter)

      tdter = getvar('dter',ndims,dims,stagger,min,max,missing,
     >           data_units,display_units,dim_names,flag)
      if (tdter.eq.0) then
        dter = 0.
        write(6,*) 'Error:  could not find terrain at ',
     >                'dot.'
        stop

      else
        dter_dims(1) = dims(1)
        dter_dims(2) = dims(2)
        dter_min(1) = min(1)
        dter_min(2) = min(2)
        dter_max(1) = max(1)
        dter_max(2) = max(2)
        dter_delta(1) = delta(1)
        dter_delta(2) = delta(2)
        dter = getmem(dims(1)*dims(2))
        call convert_units(data_units,dom_units(3),slope,intercept)
        do i=1,dims(1)*dims(2)
          call s_val(%val(dter),i,
     >         real(r_val(%val(tdter),i)*slope+intercept))
        enddo
        coord_dep(1,3) = 1
        coord_dep(2,3) = 1
      endif
      print *,dims(1),dims(2)
      print *,min(1),min(2)
      print *,max(1),max(2)
      print *,delta(1),delta(2)
c
cc set up lat/lon at cross pts
c
      if (xlat.ne.0.) call freemem(xlat)

      txlat = getvar('xlat',ndims,dims,stagger,min,max,missing,
     >           data_units,display_units,dim_names,flag)
      if (txlat.eq.0) then
        xlat = 0.
        write(6,*) 'Error:  could not find latitude at ',
     >                'thermodynamic pts.'
        stop

      else
        xlat = getmem(dims(1)*dims(2))
        do i=1,dims(1)*dims(2)
          call s_val(%val(xlat),i,
     >         real(r_val(%val(txlat),i)*slope+intercept))
        enddo
        coord_dep(1,3) = 1
        coord_dep(2,3) = 1
      endif
c
c set up longitude at cross points
c
      if (xlon.ne.0.) call freemem(xlon)

      txlon = getvar('xlon',ndims,dims,stagger,min,max,missing,
     >           data_units,display_units,dim_names,flag)
      if (txlon.eq.0) then
        xlon = 0.
        write(6,*) 'Error:  could not find longitude at ',
     >                'cross points.'
        stop

      else
        xlon = getmem(dims(1)*dims(2))
        do i=1,dims(1)*dims(2)
          call s_val(%val(xlon),i,
     >         real(r_val(%val(txlon),i)*slope+intercept))
        enddo
        coord_dep(1,3) = 1
        coord_dep(2,3) = 1
      endif
c
c get flag for grid transformation
c
      ttgrid = getvar('tgrid',ndims,dims,stagger,min,max,missing,
     >            data_units,display_units,dim_names,flag)
      if (ttgrid.eq.0) then
        tgrid = 1.
      else
        tgrid = r_val(%val(ttgrid),1)
      end if
      print *,"tg= ",tgrid

c
c set up domain info common block
c
      call getrarr("plmax",plmax,4,error)
      if (error) then
        write(6,*) 'Error:  could not get physical domain maximums',
     >                                ' in new_file.'
        stop
      end if
      call getrarr("plmin",plmin,4,error)
      if (error) then
        write(6,*) 'Error:  could not get physical domain minimums',
     >                                ' in new_file.'
        stop
      end if

      if (plmax(2).eq.plmin(2)) then
        twod = 1
      else
        twod = 0
      end if
      print *,twod,plmax(2),plmin(2)

      dx = delta(1)
      nx = nint((plmax(1)-plmin(1))/dx)+1
      if (twod.ne.1) then
        dy = delta(2)
        ny = nint((plmax(2)-plmin(2))/dy)+1
      else
        dy = 1.
        ny = 1
      end if
      dz = delta(3)
      nz = nhsig
      print *,dx,nx,dy,ny,dz,nz

c
c create additional fields
c Create theta
c
c      if (xtheta.ne.0.) call freemem(xtheta)

c      ndims = 4
cc      print *,xmin(1),xmin(2),xmin(3),xmin(4)
cc      print *,xmax(1),xmax(2),xmax(3),xmax(4)
c      xdim_names(4) = 'time'
c      data_units = 'K'
c      display_units = 'K'

c      txtheta = calc_field('THETA',ndims, xdims, xstagger, xmin, 
c     &    xmax,xmissing, data_units, display_units, xdim_names)

c      if (txtheta.eq.0) then
c        xtheta = 0.
c        write(6,*) 'Error:  could not calculate theta at ',
c     >                'thermodynamic pts.'
c        stop

c      else
c        xtheta = getmem(xdims(1)*xdims(2)*xdims(3)*xdims(4))
c        do i=1,xdims(1)*xdims(2)*xdims(3)*xdims(4)
c          call s_val(%val(xtheta),i,
c     >         real(r_val(%val(txtheta),i)))
c        enddo
c        coord_dep(1,3) = 1
c        coord_dep(2,3) = 1
c      endif

c      txtheta = addvar('THETA','THETA',ndims,xdims,xstagger,xmin,
c     >        xmax,xmissing,data_units,display_units,xdim_names,flag)

c
c Total accumulated precip
c
c      if (xptot.ne.0.) call freemem(xptot)

c      ndims = 4
cc      print *,xmin(1),xmin(2),xmin(3),xmin(4)
cc      print *,xmax(1),xmax(2),xmax(3),xmax(4)
c      xdim_names(4) = 'one'
c      data_units = 'mm'
c      display_units = 'mm'      

c      txptot =  calc_field('PTOT', ndims, xdims, xstagger, xmin, 
c     &    xmax,xmissing, data_units, display_units, xdim_names)

c      if (txptot.eq.0) then
c        xptot = 0.
c        write(6,*) 'Error:  could not calculate precip at ',
c     >                'thermodynamic pts.'
c        stop

c      else
c        xptot = getmem(xdims(1)*xdims(2)*xdims(3)*xdims(4))
c        do i=1,xdims(1)*xdims(2)*xdims(3)*xdims(4)
c          call s_val(%val(xptot),i,
c     >         real(r_val(%val(txptot),i)))
c        enddo
c        coord_dep(1,3) = 1
c        coord_dep(2,3) = 1
c      endif

c      txptot = addvar('PTOT','PTOT',ndims,xdims,xstagger,xmin,
c     >        xmax,xmissing,data_units,display_units,xdim_names,flag)

c
c wind speed
c
c     if (xwspd.ne.0.) call freemem(xwspd)

c      ndims = 4
cc      print *,dmin(1),dmin(2),dmin(3),dmin(4)
cc      print *,dmax(1),dmax(2),dmax(3),dmax(4)
cc      pause

c      ddim_names(4) = 'time'
c      data_units = 'm/s'
c      display_units = 'm/s'      

c      txwspd =  calc_field('WSPD', ndims, ddims, dstagger, dmin, 
c     &    dmax,dmissing, data_units, display_units, ddim_names)

c      if (txwspd.eq.0) then
c        xwspd = 0.
c        write(6,*) 'Error:  could not calculate wind speed at ',
c     >                'dot pts.'
c        stop

c      else
c        xwspd = getmem(ddims(1)*ddims(2)*ddims(3)*ddims(4))
c        do i=1,ddims(1)*ddims(2)*ddims(3)*ddims(4)
c          call s_val(%val(xwspd),i,
c     >         real(r_val(%val(txwspd),i)))
c        enddo
c        coord_dep(1,3) = 1
c        coord_dep(2,3) = 1
c      endif

c      txwspd = addvar('WSPD','WSPD',ndims,ddims,dstagger,dmin,
c     >        dmax,dmissing,data_units,display_units,ddim_names,flag)

c
c ***
c
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

      ive_ptr read_var

      character*(*) name
      integer       inmax1, inmax2, inmax3, inmax4
      real          field(inmax1, inmax2, inmax3, inmax4)

      character*128 dnames(4)
      character*81  data_units, display_units, dim_names(4)
      real          missing, stagger(4), start, r_val
      integer       i, j, min(4), max(4), dims(4), nd, error, ndims
      logical       flag

      print *,'In get_field'

c
c     Save the physical and index limits in the common block in
c     default.icl for use in later calculations by other transform
c     routines.
c
      call getrarr("phmin", phmin, 4, error)
      call getrarr("phmax", phmax, 4, error)
      call getrarr("grid_delta", delta, 4, error)
      inmax(1) = inmax1
      inmax(2) = inmax2
      inmax(3) = inmax3
      inmax(4) = inmax4
cc      print *,inmax1,inmax2,inmax3,inmax4
cc      print *,phmin(1),phmin(2),phmin(3),phmin(4)
cc      print *,phmax(1),phmax(2),phmax(3),phmax(4)
cc      print *,delta(1),delta(2),delta(3),delta(4)
c
c     Get the names of the dimensions of this field.
c
      call getaarr("data_dim_names", dnames, 4, error)
      call getivar("num_dims", ndims, error)
cc      print *,dnames(1),dnames(2),dnames(3),dnames(4)
cc      print *,ndims
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
            if (delta(i) .eq. 0.0) then
               delta(i) = (phmax(i)-phmin(i))/(inmax(i)-1)
cc               print *,phmax(i),phmin(i),inmax(i),i
cc               print *,"j",delta(i)
            endif
c
c     Get the memory address of the values of the variable with
c     the same name as the dimension.
c     If such a variable does not exist, coord(i) == 0.
c
            flag = .false.
            coord(i) = read_var(dnames(i), nd, dims, stagger, min, max,
     &              missing, data_units, display_units, dim_names, flag)

c set coord 3 to sigma half levels

            if (i.eq.3) coord(i) = xght
            print *,"c",coord(i),i
            if (coord(i) .ne. 0) then
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
                     print *,"delta= ",delta(i),i
                      goto 10
                  endif
               enddo
 10            continue
            endif
         endif
      enddo

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
      real  zter
      
      do i = 1, npts
         do j = 1, ndims
            if (flag(j) .eq. 1) then
c
c     Trivial case where dimension does not exist (depth of array is 
c     unity along this dimension)
c
               if (inmax(j) .eq. 1 .or. index(j,i) .lt. 1) then
                  phys(j,i) = phmin(j)
c
c     Case where interpolation must be performed between the grid
c     point locations given in a 1D "dimension variable" array.
c
               else if (delta(j) .eq. 0.) then
                 k = index(j,i)
                 if (npts.eq.int(dght_dims(1)).or.
     &               npts.eq.int(dght_dims(2))) then
c                  print *,"dot"
                  phys(j,i) = get_ght(%val(dght),phys(1, i),
     &              phys(2, i),k,int(dght_dims(1)),int(dght_dims(2)),
     &              int(dght_dims(3)),dght_mind, dght_delta,0)
                 else 
c                  print *,"cross"
                  phys(j,i) = get_ght(%val(xght),phys(1, i),
     &              phys(2, i),k,int(xght_dims(1)),int(xght_dims(2)),
     &              int(xght_dims(3)),xght_mind, xght_delta,1)
                 endif

c                  inxk = index_3d(phys(1, i),phys(2, i),k, 
c     &                  xght_dims ,xght_mind, xght_delta)

c                  print *,k,phys(j,i),j,i,npts
c                  pause
c
c     Case of evenly spaced physical coord.
c
               else
                  phys(j,i) = phmin(j) + (index(j,i)-1)*delta(j)
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

      integer npts, ndims, flag(ndims),type
      real    phys(ndims, npts), index(ndims, npts), testh, testhm1
      real    phmtest

      integer i, j, k
      real    r_val	! function
      

      do i = 1, npts
         do j = 1, ndims
            if (flag(j) .eq. 1) then
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


                 if (npts.eq.int(dght_dims(1)).or.
     &               npts.eq.int(dght_dims(2))) then
                    type = 0
                 else
                    type = 1
                 endif

c
c     Physical coord. increasing with index.
c
                  if ((phmin(j) .lt. phmax(j)).and.(j.eq.3)) then
                  phmtest =  get_ter_pt (%val(xter), phys(1,i), 
     &              phys(2,i), xter_dims(1), xter_dims(2), xter_min,
     &              dter_delta)
                     phmtest = amax1(phmtest,phmin(j))
                     phys(j,i) = amax1(phys(j,i),phmin(j))
c                     print *,phmin(j),phmax(j),phys(j,i),npts,j
c                     print *,phmtest
c                     pause
                     if (phys(j,i) .lt. phmtest) then
                        index(j,i) = -1.
                        goto 5
                     else if (phys(j,i) .ge. phmax(j)) then
                        index(j,i) = inmax(j)
                     else
                        do k = 2, inmax(j)
                 if (type.eq.0) then          
c                   print *,"dot"
                   testh = get_ght(%val(dght),phys(1, i),
     &             phys(2, i),k,int(dght_dims(1)),int(dght_dims(2)),
     &             int(dght_dims(3)),dght_mind, dght_delta,0)
                 else
c                   print *,"cross"
                   testh= get_ght(%val(xght),phys(1, i),
     &             phys(2, i),k,int(xght_dims(1)),int(xght_dims(2)),
     &             int(xght_dims(3)),xght_mind, xght_delta,1)
                 endif
                           if (phys(j,i).le.testh)
     &                          goto 10
                        enddo
c 10                     print *,phys(j,i),testh,k
c                        pause
 10              if (type.eq.0) then
c                   print *,"dot"
                   testhm1 = get_ght(%val(dght),phys(1, i),
     &             phys(2, i),k,int(dght_dims(1)),int(dght_dims(2)),
     &             int(dght_dims(3)),dght_mind, dght_delta,0)
                 else
c                   print *,"cross"
                   testhm1 = get_ght(%val(xght),phys(1, i),
     &             phys(2, i),k,int(xght_dims(1)),int(xght_dims(2)),
     &             int(xght_dims(3)),xght_mind, xght_delta,1)
                 endif

                        index(j,i) = float(k) - (phys(j,i)-
     &                       testh)/(testhm1-testh)
c                        print *,index(j,i)
                     endif

c for other grids

                  elseif ((phmin(j).lt.phmax(j)).and.(j.ne.3)) then
                     if (phys(j,i) .le. phmin(j)) then
                        index(j,i) = 1.
                     else if (phys(j,i) .ge. phmax(j)) then
                        index(j,i) = inmax(j)
                     else
                       k = 2
 15                    if ((phys(j,i).gt.r_val(%val(coord(j)),k)).and.
     >                           (k.lt.inmax(j))) then
                         k = k+1
                         goto 15
                       end if
                       index(j,i) = real(k) - (phys(j,i)-
     &                       r_val(%val(coord(j)),k))/
     &                       (r_val(%val(coord(j)),k-1)-
     &                       r_val(%val(coord(j)),k))
                     endif
c
c     Physical coord. decreasing with index.
c
                  else
                     if (phys(j,i) .ge. phmin(j)) then
                        index(j,i) = 1.
                     else if (phys(j,i) .le. phmax(j)) then
                        index(j,i) = inmax(j)
                     else
                        do k = 2, inmax(j)
                           if (phys(j,i).ge.r_val(%val(coord(j)),k))
     &                          goto 20
                        enddo
 20                     index(j,i) = float(k) - (phys(j,i)-
     &                       r_val(%val(coord(j)),k))/
     &                       (r_val(%val(coord(j)),k-1)-
     &                       r_val(%val(coord(j)),k))
                     endif
                  endif
c
c     Case of evenly spaced physical coord.
c
               else
                  index(j,i) = (phys(j,i)-phmin(j))/delta(j) + 1
               endif
            endif
 5          continue
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

#include "default.icl"

      integer nxw, nyw, error, i, j
      real    topo(nxw, nyw), stagi(nxw), stagj(nyw), zero
      real    r_val	! function
      real    phypt(4)

      error = 0

      do j = 1, nyw
         do i = 1, nxw
            topo(i, j) = zero
c
c     Find the location of (i, j) in physical space. Note that
c     cpmpxy is an IVE routine and not part of the transform
c     package
c
            call cpmpxy (2, real(i), real(j), phypt(1), phypt(2))
c
c     Convert to a point in the data grid.
c
            topo(i, j) =  get_ter_pt (%val(xter), phypt(1), 
     &           phypt(2), xter_dims(1), xter_dims(2), xter_min,
     &           dter_delta)
            topo(i, j) = amax1 (topo(i, j), zero)
         enddo
      enddo


      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine vert_ter(xxter, yter, topo_ht, nter, pt1, pt2, zero)
c     
c     vert_ter : This routine fills an array with terrain heights in the
c     windowed domain for vertical cross sections.
c
c     Arguments:
c	xxter	ive_ptr	output	Pointer to array of x-locations for topo_ht
c	yter	ive_ptr	output	Pointer to array of y-locations for topo_ht
c	topo_ht ive_ptr	output	Pointer to array containing 
c				the terrain heights.
c	nter	integer	output	Number of points in topo_ht (output)
c	pt1	real	input	Minimum physical coordinates for slice.
c	pt2	real	input	Maximum physical coordinates for slice.
c	zero	real	input	Minimum terrain height.
c

#include "default.icl"

      integer nter
      real    pt1(4), pt2(4), zero
      ive_ptr xxter, yter, topo_ht, getmem
c
c
c     Local variable declarations.
c
      integer MAXDIM
      parameter (MAXDIM = 4)
      integer   i, iflag
      real      phypt, min, max
      dimension iflag(MAXDIM),phypt(MAXDIM),min(MAXDIM), max(MAXDIM)
c
      integer stepby
      real    value
c
      print *,"vert_ter"
      if ( xter .ne. 0 ) then
         iflag(1) = 1
         iflag(2) = 1
         iflag(3) = 0
         iflag(4) = 0
         print *,pt1(1),pt1(2),pt1(3),pt1(4)
         print *,pt2(1),pt2(2),pt2(3),pt2(4)
         call phys_2_index(pt1, min, iflag, MAXDIM, 1)
         call phys_2_index(pt2, max, iflag, MAXDIM, 1)
         print *,min(1),max(1),min(2),max(2)
         if (max(1) .eq. min(1)) then
            stepby = 2
         else if (max(2) .eq. min(2)) then
            stepby = 1
         else
            stepby = 0
         endif
c
c     Handle case where slice is along coordinate axis.
c     Step by computational coordinate.
c
         if (stepby .ne. 0) then
            if (max(stepby) .eq. int(max(stepby))) then
               i = int(max(stepby))
            else
               i = int(max(stepby) + 1)
            endif
            nter = i - int(min(stepby)) + 1
            xxter = getmem(nter)
            yter = getmem(nter)
            topo_ht = getmem(nter)
            value =  get_ter_pt (%val(xter), pt1(1), pt1(2),
     &              xter_dims(1), xter_dims(2), xter_min,
     &              xter_delta)
            value = amax1(value, zero)
            print *,value,pt1(1),pt1(2)," 1"
            call s_val(%val(xxter), 1, pt1(1))
            call s_val(%val(yter), 1, pt1(2))
            call s_val(%val(topo_ht), 1, value)
            do i = 2, nter-1
               min(stepby) = int(min(stepby)+1)
               call index_2_phys(phypt, min, iflag, MAXDIM, 1)
               value =  get_ter_pt (%val(xter), phypt(1), phypt(2),
     &              xter_dims(1), xter_dims(2), xter_min,
     &              xter_delta)
               value = amax1(value, zero)
               print *,value,phypt(1),phypt(2),i
               call s_val(%val(xxter), i, phypt(1))
               call s_val(%val(yter), i, phypt(2))
               call s_val(%val(topo_ht), i, value)
            enddo
            value =  get_ter_pt (%val(xter), pt2(1), pt2(2),
     &           xter_dims(1), xter_dims(2), xter_min,
     &           xter_delta)
            value = amax1(value, zero)
            print *,value,pt2(1),pt2(2)," last"
            call s_val(%val(xxter), nter, pt2(1))
            call s_val(%val(yter), nter, pt2(2))
            call s_val(%val(topo_ht), nter, value)
c
c     Handle case of arbitrary slice - step by physical coordinate
c
         else
            print *,'any slice'
            nter = 2 * int(max(2) - min(2) + max(1) - min(1))
            print *,max(2),min(2),max(1),min(1)
            print *,nter
            xxter = getmem(nter+1)
            yter = getmem(nter+1)
            topo_ht = getmem(nter+1)
            value =  get_ter_pt (%val(xter), pt1(1), pt1(2),
     &           xter_dims(1), xter_dims(2), xter_min,
     &           xter_delta)
            value = amax1(value, zero)
            print *,value,pt1(1),pt1(2)," 1"
            call s_val(%val(xxter), 1, pt1(1))
            call s_val(%val(yter), 1, pt1(2))
            call s_val(%val(topo_ht), 1, value)
            do i = 2, nter
               phypt(1) = pt1(1) + i * (pt2(1) - pt1(1)) / real(nter)
               phypt(2) = pt1(2) + i * (pt2(2) - pt1(2)) / real(nter)
               value =  get_ter_pt (%val(xter), phypt(1), phypt(2),
     &              xter_dims(1), xter_dims(2), xter_min,
     &              xter_delta)
               value = amax1(value, zero)
               print *,value,phypt(1),phypt(2),i
               call s_val(%val(xxter), i, phypt(1))
               call s_val(%val(yter), i, phypt(2))
               call s_val(%val(topo_ht), 1, value)
            enddo
         endif
c
c     Handle no terrain case.
c
      else
         nter = 2
         xter = getmem(2)
         yter = getmem(2)
         topo_ht = getmem(2)
         call s_val(%val(xxter), 1, pt1(1))
         call s_val(%val(xxter), 2, pt2(1))
         call s_val(%val(yter), 1, pt1(2))
         call s_val(%val(yter), 2, pt2(2))
         call s_val(%val(topo_ht), 1, zero)
         call s_val(%val(topo_ht), 2, zero)
      endif
 
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
#include "default.icl"

      integer npts
      real    x(npts), y(npts), lon(npts), lat(npts)
      real r_val

      integer i

      print *,'in phys_2_lonlat'
c      pause

      do i = 1, npts
         lon(i) = r_val(%val(xlon),i)
         lat(i) = r_val(%val(xlat),i)
c         print *,i,npts,lon(i),lat(i)
         pause
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

      print *,'in lonlat 2 phys'
      pause

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

      print *,'in slice_1d'
c      pause   

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
      ive_ptr function data_slicer_2d(corner, ii, jj, ni, nj, da, nda,
     &     ri, rj)
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
c       ri	real	input	Multiply normal # x-points by this
c       rj	real	input	Multiply normal # y-points by this
c
c     The corners will be:     2    3
c                           jj
c                              0    1
c                                ii
c
      real corner(4,4), ri, rj
      integer ii, jj, ni, nj, da, nda
c
c     Call the internal routine that does 2d physical slicing for
c     vectors.
c
      ive_ptr pslicer2d

      print *,'in slice_2d'
c      pause

      data_slicer_2d = pslicer2d(corner, ii, jj, ni, nj, da, nda,
     &     ri, rj)
      return
      end
c
c     -----------------------------------------------------------------
c
      ive_ptr function data_slicer_3d(corner, ii, jj, kk, ni, nj, nk,
     &     da, nda)
c
c     data_slicer_3d_ : This routine returns the computational points at
c     which a 3d slice should be taken.
c
c     Arguments:
c     	corner	real	input	2 corner points of slice in physical space.
c     	ii	integer	input	axis coordinate of ni.
c     	jj	integer	input	axis coordinate of nj.
c     	ni	integer	output	Number of points in the x-direction.
c     	nj	integer	output	Number of points in the y-direction.
c     	da	integer	input	axis coordinate of average
c     	nda	integer	output	Number of points in the da-direction.
c       ri	real	input	Multiply normal # x-points by this
c       rj	real	input	Multiply normal # y-points by this
c
c
      real corner(4,2)
      integer ii, jj, kk, ni, nj, nk, da, nda
c
c     Call the internal routine that does 2d physical slicing for
c     vectors.
c
      ive_ptr pslicer3d

      print *,'in slice_3d'
c      pause

      data_slicer_3d = pslicer3d(corner, ii, jj, kk, ni, nj, nk,
     &     da, nda)
      return
      end
c
c     -----------------------------------------------------------------
c
      subroutine heading(where, which, line)
c
c     heading : This routine allows the user to change headings on plots
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
      real function file_coordinate(value, coord, ncid)
c
c     file_coordinate : This routine allows the user to change the current
c     value of the file coordinate from what IVE thinks it should be.
c
c     Arguments:
c     value	float	input	The current value of the file coordinate
c     coord	integer	input	Index of the file coordinate (0-3).
c     ncid	integer	input	Id of currently opened netCDF file.
c
      real value
      integer coord, ncid

      file_coordinate = value
      return
      end
c
c     -----------------------------------------------------------------
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

C ****************************************************************

      real function get_ter_pt(ter,x,y,nx,ny,ter_min,del)

      integer nx,ny
      real ter(nx,ny),ter_min(2),del(2),x,y

      integer i,j,ip1,jp1
      real ri,rj,deli,delj

      ri = (x - ter_min(1))/del(1) + 1.
      rj = (y - ter_min(2))/del(2) + 1.

      i = max(min(int(ri),nx-1),1)
      j = max(min(int(rj),ny-1),1)

      deli = ri - real(i)
      delj = rj - real(j)

      ip1 = min(i+1,nx)
      jp1 = min(j+1,ny) 

c      print *,"ter ",nx,ny,ter(i,j),i,j,ri,rj

      get_ter_pt = (1.-deli)*(1.-delj)*ter(i,j)
     >     +  deli*(1.-delj)*ter(ip1,j)
     >     +  delj*(1.-deli)*ter(i,jp1)
     >     +  deli*delj*ter(ip1,jp1)

      return
      end


C ****************************************************************

      real function get_ght(ght,x,y,k,nx,ny,nz,mind,del,flg)

      integer nx,ny,nz,k,flg
      real ght(nx,ny,nz),mind(2),del,x,y

      integer i,j,ip1,jp1
      real ri,rj,deli,delj

      ri = (x - mind(1))/del + 1.
      rj = (y - mind(2))/del + 1.

      i = max(min(int(ri),nx-1),1)
      j = max(min(int(rj),ny-1),1)

      deli = ri - real(i)
      delj = rj - real(j)

      ip1 = min(i+1,nx)
      jp1 = min(j+1,ny) 

c      print *,"ght ",nx,ny,nz,ght(i,j,k),i,j,ri,rj

      get_ght = (1.-deli)*(1.-delj)*ght(i,j,k)
     >     +  deli*(1.-delj)*ght(ip1,j,k)
     >     +  delj*(1.-deli)*ght(i,jp1,k)
     >     +  deli*delj*ght(ip1,jp1,k)

c      get_ght = ght(i,j,k) 

c      print *,get_ght
c      pause

      return
      end

C ****************************************************************

      integer function index_3d(x,y,k,dims,mind,del)

      real mind(2), del,x,y, dims(3)

      integer i,j,k,nx,ny
      real ri,rj

      nx = int(dims(1))
      ny = int(dims(2))

      ri = (x - mind(1))/del + 1.
      rj = (y - mind(2))/del + 1.


      i = max(min(int(ri),nx-1),1)
      j = max(min(int(rj),ny-1),1)


      index_3d = i+((j-1)*ny)+((k-1)*nx*ny)

      return
      end
