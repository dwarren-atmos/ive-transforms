c     $Id: default.F,v 1.8 1996/11/01 19:41:53 cepi Exp $
c     $Log: default.F,v $
c     Revision 1.8  1996/11/01 19:41:53  cepi
c     Added checks to make sure we don't exceed array bounds when interpolating
c     terrain.
c
c     Revision 1.7  1996/11/01 19:08:15  cepi
c     i)   Got rid of the vertical coordinate array for derivative pts.  From
c          now on we'll assume that all fields are at either the w pts or the
c          thermo pts (in the vertical, at least).  This means we'll have to
c          fill out derived fields with missing values.
c     ii)  Added units conversion for terrain heights.
c     iii) Fixed terrain interpolations in get_ter_pt to use correct
c          set of pts.
c
c     Revision 1.6  1996/10/18 22:02:39  harry
c     Change float calls to real (which is the standard conversion function).
c
c     Revision 1.5  1996/10/18 21:35:05  harry
c     Fix default.F to include the changes from cepi in version 1.2 that I
c     accidentally lost in version 1.3.
c
c     Revision 1.4  1996/10/18 20:22:06  harry
c     slope and intercept should be real*8.
c
c     Revision 1.3
c     Fix bug in get_ter_pt
c
c     Revision 1.2  1996/07/22 22:48:37  cepi
c     New version.  Modified way in which stretching is handled.  Minor
c     correction to terrain interpolations.  Moved calc_field to a
c     separate file.  Some other minor stuff.
c     NOTE:  this version won't work with data files created before
c            july 12 update of usrout.f
c
c     Revision 1.1  1996/07/08 20:48:59  harry
c     Initial check in.
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
c     -----------------------------------------------------------------
c
      subroutine new_file(ncid,exact_coord,n,coord_dep)
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
#include "misc.icl"
      integer ncid,n(4),coord_dep(4,4)
      ive_ptr exact_coord(4)
c
c common block cdf_dimvar is set up by IVE to contain the values of a 
c   variable with the same name as the dimension of the largest variable 
c   in the netCDF file.  If the size of the domain is set by attributes, 
c   or if a dimension variable does not exist, size will be zero, and 
c   data will be NULL.
c
      ive_ptr dimvar_data(4)
      integer dimvar_size(4)
c      common/cdf_dimvar/data,size
c
c common block t_cdfid holds the netCDF id of the current input file
c 
      common/t_cdfid/idcdf
      integer idcdf
c
c numfls keeps track of the number of output files produced by the
c   the transforms for a given data set; here we just initialize 
c   it to zero
c     
      common/flcount/numfls
      integer numfls

      ive_ptr getmem,tsgz,twgz,ttime,tzbot_p,tzbot_u,tzbot_v,tztop
      ive_ptr tfcor,ttgrid
      logical error,flag
      real r_val,tgrid,stagger(4),min(4),max(4),missing,mydelta(4)
      real*8 slope,intercept
      integer i,j,ndims,dims(4),k
      character*80 data_units,display_units,dim_names(4),dom_units(4)

      data wgz,sgz,time,zbot_p,zbot_u,zbot_v/0,0,0,0,0,0/

      nobottom = .false.
      call getiarr('dimvar_size', dimvar_size, 4, error)
      call getparr('dimvar_data', dimvar_data, 4, error)

      idcdf = ncid
      numfls = 0

      do j = 1,4
      do i = 1,4
        if (i .eq. j) then
          coord_dep(i,j) = 1
        else
          coord_dep(i,j) = 0
        endif
      enddo
      enddo

      do i = 1,4
        exact_coord(i) = dimvar_data(i)
        n(i) = dimvar_size(i)
      enddo
c
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
c set up coordinate vector for w pts
c
      if (wgz.ne.0) call freemem(wgz)

      print *,'BEFORE WGZ GETVAR'
      twgz = getvar('wgz',ndims,dims,stagger,min,max,missing,
     >            data_units,display_units,dim_names,flag)
      print *,'AFTER WGZ GETVAR'

      if (twgz.eq.0) then
        wgz = 0
        write(6,*) 'Error: could not find vertical coordinate ',
     >                     'vector wgz.'
        stop

      else

        nwgz = dims(1)
        wgz = getmem(nwgz)
        call convert_units(data_units,dom_units(3),slope,intercept)
        do i = 1,nwgz
          call s_val(%val(wgz),i,
     >          real(r_val(%val(twgz),i)*slope+intercept))
        enddo
        wgz_min = r_val(%val(wgz),1)
        wgz_max = r_val(%val(wgz),nwgz)

      endif
      print *,'WGZ_MIN,WGZ_MAX',wgz_min,wgz_max
c
c set up coordinate vector for thermo pts
c
      if (sgz.ne.0) call freemem(sgz)

      tsgz = getvar('sgz',ndims,dims,stagger,min,max,missing,
     >                 data_units,display_units,dim_names,flag)
      if (tsgz.eq.0) then
        sgz = 0
        write(6,*) 'Error:  could not find vertical coordinate ',
     >                     'vector sgz.'
        stop
      else
        nsgz = dims(1)
        sgz = getmem(nsgz)
        call convert_units(data_units,dom_units(3),slope,intercept)
        do i = 1,nsgz
          call s_val(%val(sgz),i,
     >          real(r_val(%val(tsgz),i)*slope+intercept))
        enddo
        sgz_min = r_val(%val(sgz),1)
        sgz_max = r_val(%val(sgz),nsgz)
      end if
      print *,'SGZ_MIN,SGZ_MAX',sgz_min,sgz_max
c
c set up coordinate vector for times
c
      if (time.ne.0) call freemem(time)

      ttime = dimvar_data(4)
      if (ttime.eq.0) then
        time = 0
        write(6,*) 'Error:  could not find time coordinate vector.'
        stop
      else
        ntime = dimvar_size(4)
        time = getmem(ntime)
        do i = 1,ntime
          call s_val(%val(time),i,r_val(%val(ttime),i))
        enddo
      end if
c
c get total domain depth (not necessarily the same as wgz_max)
c
      mydelta=delta
      write(6,*)
      tztop = getvar('ztop',ndims,dims,stagger,min,max,missing,
     >            data_units,display_units,dim_names,flag)
      delta=mydelta
      if (tztop.eq.0) then
        write(6,*) 'Error:  could not find total domain depth ztop.'
        write(6,*) 'Could be old version of model.  Assuming ztop ',
     >                        'is wgz_max.'
        ztop = wgz_max
      else
        call convert_units(data_units,dom_units(3),slope,intercept)
        ztop = real(r_val(%val(tztop),1)*slope+intercept)
      end if
c
c set up terrain at pressure pts; note that we use only the terrain at 
c   the pressure points when transforming coordinates
c
      if (zbot_p.ne.0) call freemem(zbot_p)

      tzbot_p = getvar('zbot_p',ndims,dims,stagger,min,max,missing,
     >           data_units,display_units,dim_names,flag)
      if (tzbot_p.eq.0) then
        zbot_p = 0
        write(6,*) 'Could not find terrain at pressure pts',
     >                ' Hope your run had no bottom.'
        nobottom = .true.
c        stop
      else
        zbot_dims(1) = dims(1)
        zbot_dims(2) = dims(2)
        zbot_min(1) = min(1)
        zbot_min(2) = min(2)
        zbot_max(1) = max(1)
        zbot_max(2) = max(2)
        zbot_delta(1) = delta(1)
        zbot_delta(2) = delta(2)
        zbot_p = getmem(dims(1)*dims(2))
        call convert_units(data_units,dom_units(3),slope,intercept)
        do i=1,dims(1)*dims(2)
          call s_val(%val(zbot_p),i,
     >         real(r_val(%val(tzbot_p),i)*slope+intercept))
        enddo
        coord_dep(1,3) = 1
        coord_dep(2,3) = 1
      endif
c
c set up terrain at x velocity pts
c
      if (zbot_u.ne.0) call freemem(zbot_u)

      tzbot_u = getvar('zbot_u',ndims,dims,stagger,min,max,missing,
     >           data_units,display_units,dim_names,flag)
      if (tzbot_u.eq.0) then
        zbot_u = 0
        write(6,*) 'Could not find terrain at ',
     >                   'x velocity pts. Hope your run had no bottom'
c        stop
        nobottom = .true.
      else
        zbot_u = getmem(dims(1)*dims(2))
        call convert_units(data_units,dom_units(3),slope,intercept)
        do i=1,dims(1)*dims(2)
          call s_val(%val(zbot_u),i,
     >         real(r_val(%val(tzbot_u),i)*slope+intercept))
        enddo
      endif
c
c set up terrain at y velocity pts
c
      if (zbot_v.ne.0) call freemem(zbot_v)

      tzbot_v = getvar('zbot_v',ndims,dims,stagger,min,max,missing,
     >           data_units,display_units,dim_names,flag)
      if (tzbot_v.eq.0) then
        zbot_v = 0
        write(6,*) 'Could not find terrain at ',
     >                   'y velocity pts. Hope your run had no bottom'
c        stop
        nobottom = .true.
      else
        zbot_v = getmem(dims(1)*dims(2))
        call convert_units(data_units,dom_units(3),slope,intercept)
        do i=1,dims(1)*dims(2)
          call s_val(%val(zbot_v),i,
     >         real(r_val(%val(tzbot_v),i)*slope+intercept))
        enddo
      endif
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

      dx = delta(1)
      nx = nint((plmax(1)-plmin(1))/dx) 
      if (twod.ne.1) then
        dy = delta(2)
        ny = nint((plmax(2)-plmin(2))/dy)
      else
        dy = 1.
        ny = 1
      end if
      dz = wgz_max/nsgz
      nz = nsgz
c
c get coriolis parameter
c
      tfcor = getvar('f_cor',ndims,dims,stagger,min,max,missing,
     >           data_units,display_units,dim_names,flag)
      if (tfcor.eq.0) then
        write(6,*) 'Error:  could not find coriolis parameter.'
c       stop
      else
        fcor = r_val(%val(tfcor),1)
      end if
c
c compute metric arrays (subroutine in metrics.F)
c
      if(.not. nobottom) then
         dqdx_xy = getmem(3*(nx+1)*(ny+1))
         dqdx_zt = getmem(nz+1)
         dqdx_zw = getmem(nz+1)
         dqdy_xy = getmem(3*(nx+1)*(ny+1))
         dqdy_zt = getmem(nz+1)
         dqdy_zw = getmem(nz+1)
         dqdz = getmem(3*(nx+1)*(ny+1))
         dzdq = getmem(3*(nx+1)*(ny+1))
         call metric_calc(%val(zbot_p),%val(zbot_u),%val(zbot_v),
     >        %val(sgz),%val(wgz),%val(dqdx_xy),%val(dqdx_zt),
     >        %val(dqdx_zw),%val(dqdy_xy),%val(dqdy_zt),
     >        %val(dqdy_zw),%val(dqdz),%val(dzdq),ztop,dx,dy,
     >        nx,ny,nz)
      endif
c
c set up boundary arrays (subroutine in metrics.F; arrays in misc.icl)
c
 112  format(5x,A5,I6,5x,A8,I6)

      if ((nx.gt.nxmax).or.(ny.gt.nymax).or.(nz.gt.nzmax)) then
        write(6,*) 'Error:  not enough scratch space.'
        write(6,112) 'nx = ',nx,'nxmax = ',nxmax
        write(6,112) 'ny = ',ny,'nymax = ',nymax
        write(6,112) 'nz = ',nz,'nzmax = ',nzmax
        write(6,*) 'Increase parameters in misc.icl'
      end if

      call bnd_arrays(%val(sgz),%val(wgz),nx,ny,nz)
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
      if (tgrid.ne.0.) then
        trans_on = 1
      else
        trans_on = 0
      end if
c
c set domain offsets to zero
c
      offset(1) = 0.
      offset(2) = 0.

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
      integer inmax1,inmax2,inmax3,inmax4
      real field(inmax1,inmax2,inmax3,inmax4)

      character*128 dnames(4)
      character*81 data_units,display_units,dim_names(4)
      real missing,stagger(4),start,r_val,d,dom_min(4),dom_max(4)
      integer i,j,min(4),max(4),dims(4),nd,error,ndims
      integer k,iflg
      logical flag,requ
      ive_ptr getmem
c
c save the physical and index limits in the common block in
c   default.icl for use in later calculations by other transform
c   routines.
c
      if((len(name).ge. 3).and.
     &     (name(1:3) .eq. 'sgz' .or. name(1:3) .eq. 'wgz')) then
         return
      endif
      call getrarr("phmin",phmin,4,error)
      if (error.ne.0) then
        write(6,*) 'Error:  problem getting physical domain ',
     >                               'for field.'
        stop
      end if
      call getrarr("phmax",phmax,4,error)
      if (error.ne.0) then
        write(6,*) 'Error:  problem getting physical domain ',
     >                               'for field.'
        stop
      end if
      inmax(1) = inmax1
      inmax(2) = inmax2
      inmax(3) = inmax3
      inmax(4) = inmax4
c
c translate variable domain limits if necessary (we assume a domain
c   translation greater than (dx,dy))
c
      if ((abs(phmin(1)-plmin(1)).gt.dx).and.(inmax1.ne.1)) then
        phmin(1) = phmin(1) + offset(1)
        phmax(1) = phmax(1) + offset(1)
      end if
      if ((abs(phmin(2)-plmin(2)).gt.dy).and.(inmax2.ne.1)) then
        phmin(2) = phmin(2) + offset(2)
        phmax(2) = phmax(2) + offset(2)
      end if
c
c check to see whether vertical coordinate should be transformed
c
      if ((inmax(3).ne.1).and.(inmax(1).ne.1).and.((inmax(2).ne.1)
     >             .or.(twod.eq.1)).and.(zbot_p.ne.0).and.
     >             (phmin(3).ge.0.).and.(trans_on.eq.1)) then
        trans_z = 1
      else
        trans_z = 0
      end if   
      
      if (phmin(3).lt.0.) phmin(3) = sgz_min
c
c reset the plot domain if necessary
c
c     if ((phmin(1).lt.plmin(1)).or.(phmin(2).lt.plmin(2)).or.
c    >       (phmax(1).gt.plmax(1)).or.(phmax(2).gt.plmax(2))) then
c
c       call setrarr("plmin",phmin,4,error)
c       if (error.ne.0) then
c         write(6,*) 'Error:  problem setting plot domain.'
c         stop
c       end if
c       call setrarr("plwmin",phmin,4,error)
c       if (error.ne.0) then
c         write(6,*) 'Error:  problem setting plot domain.'
c         stop
c       end if
c       call setrarr("plmax",phmax,4,error)
c       if (error.ne.0) then
c         write(6,*) 'Error:  problem setting plot domain.'
c         stop
c       end if
c       call setrarr("plwmax",phmax,4,error)
c       if (error.ne.0) then
c         write(6,*) 'Error:  problem setting plot domain.'
c         stop
c       end if
c       WRITE(*,*) 'Calling info_2d.'
c       do k=1,4
c         WRITE(*,*) 'phmin(',k,') = ', phmin(k)
c       enddo
c       call info_2d(phmin,phmax,4)
c
c     else
c
c       call getrarr("plmin",dom_min,4,error)
c       if (error.ne.0) then
c         write(6,*) 'Error:  problem setting plot domain.'
c         stop
c       end if
c       call getrarr("plmax",dom_max,4,error)
c       if (error.ne.0) then
c         write(6,*) 'Error:  problem setting plot domain.'
c         stop
c       end if
c
c       k = 1
c       iflg = 0
c567    if ((iflg.eq.0).and.(k.le.4)) then
c         if ((dom_min(k).ne.plmin(k))
c    >                   .or.(dom_max(k).ne.plmax(k))) then
c           iflg = 1
c         else
c           k = k+1
c         end if
c         goto 567
c       end if
c
c       if (iflg.ne.0) then
c
c         call setrarr("plmin",plmin,4,error)
c         if (error.ne.0) then
c           write(6,*) 'Error:  problem setting plot domain.'
c           stop
c         end if
c         call setrarr("plwmin",plmin,4,error)
c         if (error.ne.0) then
c           write(6,*) 'Error:  problem setting plot domain.'
c           stop
c         end if
c         call setrarr("plmax",plmax,4,error)
c         if (error.ne.0) then
c           write(6,*) 'Error:  problem setting plot domain.'
c           stop
c         end if
c         call setrarr("plwmax",plmax,4,error)
c         if (error.ne.0) then
c           write(6,*) 'Error:  problem setting plot domain.'
c           stop
c         end if
c         WRITE(*,*) 'Calling info_2d.'
c         do k=1,4
c           WRITE(*,*) 'plmin(',k,') = ', plmin(k)
c         enddo
c         call info_2d(plmin,plmax,4)
c
c       end if
c
c     end if
c
c set grid intervals assuming evenly spaced coordinates 
c
      do k = 1,4
        if (inmax(k).ne.1) then
          delta(k) = (phmax(k)-phmin(k))/(inmax(k)-1)
        else
          delta(k) = 1.
        end if
      enddo
c
c find the appropriate vertical coordinate vector
c
      if (requ(phmin(3),wgz_min)) then 
        coord(3) = wgz
      else if (requ(phmin(3),sgz_min)) then
        coord(3) = sgz
      else 
        write(*,*) 'Warning:  could not find vertical coordinate ',
     >               'vector.' 
        stop
      endif
c
c check to see if vertical coordinate vector is evenly spaced.  If
c    not, set delta = 0 to indicate uneven spacing. 
c
      if (coord(3).eq.0) then
         write(*,*) ' Warning:  could not find vertical coordinate ',
     >        'vector.  Assuming even spacing in vertical.'
      else
         k = 1
         start = r_val(%val(coord(3)),1)
 538     if (k.le.inmax(3)) then
            if (requ(start+(k-1)*delta(3),r_val(%val(coord(3)),k))) then
               k = k+1
            else
               delta(3) = 0.   
               k = inmax(3)+1
            end if
            goto 538
         end if
      endif
c
c check to see if time coordinate vector is evenly spaced
c
      if (inmax(4).ne.1) then
        coord(4) = time
        k = 1
        start = r_val(%val(coord(4)),1)
 978    if (k.le.inmax(4)) then
          if (requ(start+(k-1)*delta(4),r_val(%val(coord(4)),k))) then
            k = k+1
          else
            delta(4) = 0.
            k = inmax(4)+1
          end if 
          goto 978
        end if
      end if

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
      real    zter
      real    r_val	! function

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
                  if (k .lt. inmax(j)) then
                     phys(j,i) = r_val(%val(coord(j)),k) +
     &                    (index(j,i)-k)*(r_val(%val(coord(j)),k+1)
     &                    -r_val(%val(coord(j)),k))
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
c     Handle terrain
c
               if ((j.eq.3).and.(trans_z.eq.1)) then
                  zter = get_ter_pt(%val(zbot_p), phys(1, i),
     &                 phys(2, i), zbot_dims(1),
     &                 zbot_dims(2), zbot_min, zbot_delta)
                  phys(j, i) = phys(j, i) * (ztop-zter)/ztop + zter
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
      real    phys(ndims, npts), index(ndims, npts)

      integer i, j, k
      real    r_val	! function
      real    zter, p

      do i = 1, npts
         do j = 1, ndims
            if (flag(j) .eq. 1) then
               if ((j.eq.3).and.(trans_z.eq.1)) then
                  zter = get_ter_pt(%val(zbot_p), phys(1, i),
     &                 phys(2, i), zbot_dims(1),
     &                 zbot_dims(2), zbot_min, zbot_delta)
                  p = (phys(3, i) - zter) * ztop / (ztop - zter)
                  if (p .lt. 0) then
                     index(j, i) = -1.0
                     goto 5
                  endif
               else
                  p = phys(j, i)
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
                       k = 2
 10                    if ((p.gt.r_val(%val(coord(j)),k)).and.
     >                           (k.lt.inmax(j))) then
                         k = k+1
                         goto 10
                       end if
                       index(j,i) = real(k) - (p-
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
                       k = 2
 20                    if ((p.lt.r_val(%val(coord(j)),k)).and.
     >                           (k.lt.inmax(j))) then
                         k = k+1
                         goto 20
                       end if
                       index(j,i) = real(k) - (p-
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
                  if (index(j,i) .lt. 1) index(j,i) = 1
                  if (index(j,i) .gt. inmax(j)) index(j,i) = inmax(j)
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

      integer nxw, nyw, error, i, j, iflag(4)
      real    topo(nxw, nyw), stagi(nxw), stagj(nyw), zero
      real    compt(4), phypt(4)
      data iflag/1, 1, 0, 0/

      error = 0
      if (zbot_p .eq. 0) then
         do j = 1, nyw
            do i = 1, nxw
               topo(i, j) = zero
            enddo
         enddo
      else
         do j = 1, nyw
            do i = 1, nxw
c
c     Find the location of (i, j) in physical space. Note that
c     cpmpxy is an IVE routine and not part of the transform
c     package
c
               call cpmpxy (2, real(i), real(j), phypt(1), phypt(2))
c
c     Convert to a point in the data grid.
c
               call phys_2_index (phypt, compt, iflag, 4, 1)
c
c     Interpolate terrain to the required point.
c
               topo(i, j) =  get_ter_pt (%val(zbot_p), phypt(1), 
     &              phypt(2), zbot_dims(1), zbot_dims(2), zbot_min,
     &              zbot_delta)
               topo(i, j) = amax1 (topo(i, j), zero)
            enddo
         enddo
      endif

      end
c
c     -----------------------------------------------------------------
c
      subroutine vert_ter (xter, yter, topoht, nter, pt1, pt2,
     &                  zero)
c-----------------------------------------------------------------------
c     Purpose:
c
c        This routine finds the terrain height and location under each
c        point in a vertical slab.
c
c     Arguments:
c	xter	pointer	output	Pointer to array of x-locations for topo_ht
c	yter	pointer	output	Pointer to array of y-locations for topo_ht
c	topoht	pointer	output	Pointer to array containing 
c				the terrain heights.
c	nter	int	output	Number of points in topo_ht (output)
c	pt1	real	output	Minimum physical coordinates for slice.
c	pt2	real	output	Maximum physical coordinates for slice.
c	zero	real	input	the location of the bottom of the plot
c                               window. This value is used as the minimum
c                               terrain height, rather than 0.0.
c-----------------------------------------------------------------------
c
#include "default.icl"
c
c     Argument declarations.
c
      integer   nter
      real      pt1(4), pt2(4), zero
      ive_ptr   xter, yter, topoht
c
c
c     Local variable declarations.
c
      integer MAXDIM
      parameter (MAXDIM = 4)
      integer   i, iflag
      real      phypt, min, max
      dimension iflag(MAXDIM), phypt(MAXDIM),
     &		min(MAXDIM), max(MAXDIM)
      ive_ptr getmem
c
      integer stepby
      real    value
c
      if ( zbot_p .ne. 0 ) then
         iflag(1) = 1
         iflag(2) = 1
         iflag(3) = 0
         iflag(4) = 0
         call phys_2_index(pt1, min, iflag, MAXDIM, 1)
         call phys_2_index(pt2, max, iflag, MAXDIM, 1)
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
            xter = getmem(nter)
            yter = getmem(nter)
            topoht = getmem(nter)
            value =  get_ter_pt (%val(zbot_p), pt1(1), pt2(2),
     &              zbot_dims(1), zbot_dims(2), zbot_min,
     &              zbot_delta)
            value = amax1(value, zero)
            call put_ter(%val(xter), %val(yter), %val(topoht), 1, pt1,
     &           value)
            do i = 2, nter-1
               min(stepby) = int(min(stepby)+1)
               call index_2_phys(phypt, min, iflag, MAXDIM, 1)
               value =  get_ter_pt (%val(zbot_p), phypt(1), phypt(2),
     &              zbot_dims(1), zbot_dims(2), zbot_min,
     &              zbot_delta)
               value = amax1(value, zero)
               call put_ter(%val(xter), %val(yter), %val(topoht), i,
     &              phypt, value)
            enddo
            value =  get_ter_pt (%val(zbot_p), pt2(1), pt2(2),
     &           zbot_dims(1), zbot_dims(2), zbot_min,
     &           zbot_delta)
            value = amax1(value, zero)
            call put_ter(%val(xter), %val(yter), %val(topoht), nter,
     &           pt2, value)
c
c     Handle case of arbitrary slice - step by physical coordinate
c
         else
            nter = 2 * (max(2) - min(2) + max(1) - min(1))
            xter = getmem(nter+1)
            yter = getmem(nter+1)
            topoht = getmem(nter+1)
            value =  get_ter_pt (%val(zbot_p), pt1(1), pt1(2),
     &           zbot_dims(1), zbot_dims(2), zbot_min,
     &           zbot_delta)
            value = amax1(value, zero)
            call put_ter(%val(xter), %val(yter), %val(topoht), 1, pt1,
     &           value)
            do i = 2, nter
               phypt(1) = pt1(1) + i * (pt2(1) - pt1(1)) / real(nter)
               phypt(2) = pt1(2) + i * (pt2(2) - pt1(2)) / real(nter)
               value =  get_ter_pt (%val(zbot_p), phypt(1), phypt(2),
     &              zbot_dims(1), zbot_dims(2), zbot_min,
     &              zbot_delta)
               value = amax1(value, zero)
               call put_ter(%val(xter), %val(yter), %val(topoht), i,
     &              phypt, value)
            enddo
         endif
c
c     Handle no terrain case.
c
      else
         nter = 2
         xter = getmem(2)
         yter = getmem(2)
         topoht = getmem(2)
         call put_ter(%val(xter), %val(yter), %val(topoht), 1, pt1,
     &        zero)
         call put_ter(%val(xter), %val(yter), %val(topoht), 2, pt2,
     &        zero)
      endif
c
      return
      end

      subroutine put_ter(xter, yter, topoht, i, pt, value)
c
c     Load the terrain arrays
c
      integer i
      real xter(i), yter(i), topoht(i), pt(4), value

      xter(i) = pt(1)
      yter(i) = pt(2)
      topoht(i) = value

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
c     If p is the ive_ptr to some data array X, you get the j th data 
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
c     If p is the ive_ptr to some data array X, you set the j th data
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

      get_ter_pt = (1.-deli)*(1.-delj)*ter(i,j)
     >     +  deli*(1.-delj)*ter(ip1,j)
     >     +  delj*(1.-deli)*ter(i,jp1)
     >     +  deli*delj*ter(ip1,jp1)

      return
      end

      logical function requ(val1,val2)

      real val1,val2,eps,diff,big	

      parameter (eps=1.e-4)

      big = amax1(abs(val1),abs(val2))

      if (big.eq.0.) then
        requ = .true.
        return
      end if

      diff = abs((val1-val2)/big)

      if (diff.le.eps) then
        requ = .true.
      else
        requ = .false.
      end if

      return
      end

      subroutine setstr(outstr,instr)

      character*(*) outstr,instr
      integer ibeg,iend,strbeg,strend

      outstr(1:len(outstr)) = ' '
      ibeg = strbeg(instr)
      iend = strend(instr)
      outstr(1:iend-ibeg+2) = instr(ibeg:iend)//char(0)

      return
      end


      


