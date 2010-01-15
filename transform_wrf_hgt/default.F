c     This is a Fortran version of the default transforms for IVE.
c     -----------------------------------------------------------------
c
      subroutine new_file(ncid,exact_coord,n,coord_dep)

      use DEFAULT_MODULE
c
c     new_file : This routine is called every time a file is read in.
c     The routine can be used to read other information from the open
c     netCDF file.  It also can provide an array of times that limit the
c     values the time slider can be set to.
c     
c     Arguments:
c     input:
c     ncid integer   
c         Id of currently opened netCDF file.
c
c     output: 
c     exact_coord integer(ive_ptr)   
c           Pointer to array of discrete coord values that the window 
c           sliders should use.
c
c     n	integer   
c           Number of values in exact_coord.  If this
c				   is 0 => ignore values in exact_coord
c				   => no restrictions on window slider.
c
c     coord_dep   integer    
c           States the functional dependencies
c     		between the coordinates, 1 => dependent,
c     		0 => independent e.g. if z is a function of
c			x, y and z, then coord_dep(1,3),
c			coord_dep(2,3) and coord_dep(3,3)
c			will be one.  A coordinate is always 
c			assumed to be dependent on itself.

      use DEFAULT_MODULE
      use MISC_MODULE
      USE, INTRINSIC :: ISO_C_BINDING
      include 'netcdf.inc'
      integer :: ncid,n(4),coord_dep(4,4)
      integer(ive_ptr) :: exact_coord(4)
c
c common block cdf_dimvar is set up by IVE to contain the values of a 
c   variable with the same name as the dimension of the largest variable 
c   in the netCDF file.  If the size of the domain is set by attributes, 
c   or if a dimension variable does not exist, size will be zero, and 
c   data will be NULL.
c
      integer(ive_ptr) :: dimvar_data(4)
      integer :: dimvar_size(4)
c      common/cdf_dimvar/data,size
c      Use getavar routines instead
c common block t_cdfid holds the netCDF id of the current input file
c 
      integer strbeg, strend
      integer :: idcdf
      common/t_cdfid/idcdf
c
c numfls keeps track of the number of output files produced by the
c   the transforms for a given data set; here we just initialize 
c   it to zero
c     
      integer :: numfls
      common/flcount/numfls

      integer(ive_ptr) :: tlgz,tsgz,twgz,ttime,thgt_ptr,tzbot_u,tzbot_v
      integer(ive_ptr) :: tztop,ttgrid,theight
      logical :: error,error1,error2,flag
      real :: r_val,tgrid,stagger(4),min(4),max(4),missing,pls(4),junk
      double precision :: extrap
      real*8 :: slope,intercept
      integer :: i,j,ndims,dims(4),k,nxp1,nyp1,ierr,iproj
      character(len=80) :: data_units,display_units,
     &                     dim_names(4),dom_units(4)

      integer :: npts, id
      real :: lat, lon, tfcor
      real, parameter :: omega = 0.00007292

      character(len=80) :: dummy
      double precision djunk(2)
      integer(ive_ptr) :: mu_ptr,mub_ptr, ph_ptr,phb_ptr

      integer :: set_button_name

      call getiarr('dimvar_size', dimvar_size, 4, error)
      call getparr('dimvar_data', dimvar_data, 4, error)
      
      !call set_buttons

      znu_ptr = 0 ; znw_ptr=0 ; ptop_ptr = 0 ; hgt_ptr = 0
      psfc_ptr = 0 ; mutot_ptr = 0 ; phtot_ptr = 0
      xlat = 0 ; xlong = 0 ; ptot_ptr = 0

      idcdf = ncid
      numfls = 0
c
c get some domain info from IVE
c
c     Need Xlat and Xlong for lltpij and ijtoll
c
      flag=.true.
      xlat=getvar('XLAT',ndims,dims,stagger,min,max,
     >     missing,data_units,display_units,dim_names,flag)
      lat_ll=getmem(dims(4))
      lat_ur=getmem(dims(4))
      lon_ll=getmem(dims(4))
      lon_ur=getmem(dims(4))
      call transcopylats(%val(xlat),dims(1),dims(2),dims(3),dims(4),
     &     %val(lat_ll), %val(lat_ur))
c      call driver('free=XLAT',0)
      flag=.true.
      xlong=getvar('XLONG',ndims,dims,stagger,min,max,
     >     missing,data_units,display_units,dim_names,flag)
      call transcopylons(%val(xlong),dims(1),dims(2),dims(3),dims(4),
     &     %val(lon_ll), %val(lon_ur))
c      call driver('free=XLONG',0)
      ierr = set_button_name ('XLAT', '')
      ierr = set_button_name ('XLONG', '')
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

      call getrarr("plmin",plmin,4,error)
      if (error) then
        write(6,*) 'Error:  could not get plmin in new_file.'
        stop
      end if

      call getrarr("plmax",plmax,4,error)
      if (error) then
        write(6,*) 'Error:  could not get plmax in new_file.'
        stop
      end if

      ! set up coordinate dependencies
      coord_dep(:,:) = 0
      do j = 1,4 ; do i = 1,4 
        if (i .eq. j) then ! things are dependent on themselves
          coord_dep(i,j) = 1
        end if
        if  (j .eq. 3) then ! z is dependant on all
           coord_dep(i,j) = 1
        endif
      enddo ; enddo

      ! the vertical coord is a function of everything
      coord_dep(:,3) = 1 


c
c set up coordinate vector for times
c
      if (time.ne.0) call freemem(time)
      time = dimvar_data(4) ; ntime = dimvar_size(4)
      flag=.true.
      znw_ptr = getvar('ZNW',ndims,dims,stagger,min,max,missing,
     >                 data_units,display_units,dim_names,flag)
      if(znw_ptr.ne.0) then
        znw_min = r_val(%val(znw_ptr),1)
        znw_max = r_val(%val(znw_ptr),dims(1))
      end if
      nzp1 = dims(1) ; nz = dims(1)-1

      flag=.true.
      znu_ptr = getvar('ZNU',ndims,dims,stagger,min,max,missing,
     >                 data_units,display_units,dim_names,flag)
      if(znu_ptr.ne.0) then
        nznu=dims(1)
        znu_min = r_val(%val(znu_ptr),1)
        znu_max = r_val(%val(znu_ptr),dims(1))
      end if


      flag=.true.
      ptop_ptr = getvar('P_TOP',ndims,dims,stagger,min,max,missing,
     >                 data_units,display_units,dim_names,flag)
      ptop = r_val(%val(ptop_ptr),1)

      flag=.true.
      psfc_ptr = getvar('PSFC',ndims,dims,stagger,min,max,missing,
     >                 data_units,display_units,dim_names,flag)

      ierr = set_button_name ('PSFC'//char(0), '') ! Terrain
      flag=.true.
      hgt_ptr = getvar('HGT',ndims,dims,stagger,min,max,missing,
     >                 data_units,display_units,dim_names,flag)
      height_t=dims(4)
      ierr = set_button_name ('HGT'//char(0), '') ! Terrain
      zbot_dims(1) = dims(1)
      zbot_dims(2) = dims(2)
      zbot_min = min
      ! get flag for grid tranform (if it exists)
      tgrid=1. ; trans_on = 0
      ttgrid = getvar('tgrid',ndims,dims,stagger,min,max,missing,
     >            data_units,display_units,dim_names,flag)
      if (ttgrid.ne.0) tgrid = r_val(%val(ttgrid),1)
      if (tgrid.ne.0.) trans_on = 1
! try allocating memory to reduce likelihood of failure later
      mutot_ptr=getmem(zbot_dims(1)*zbot_dims(2)*35*5*3)
      if( mutot_ptr .ne. 0) then
         call freemem(mutot_ptr)
         mutot_ptr = 0
      endif
      flag=.true.
      mutot_ptr=getvar('MU_TOT',ndims,dims,stagger,min,max,
     >     missing,data_units,display_units,dim_names,flag)
      if(mutot_ptr .eq. 0) then
         mutot_ptr=addvar('MU+MUB','MU_TOT',ndims,dims,stagger,min,max,
     >        missing,data_units,display_units,dim_names,flag)
         if(mutot_ptr .eq. 0) then
            call sleep(1)
            mutot_ptr=addvar('MU+MUB','MU_TOT',ndims,dims,stagger,min,
     >           max,missing,data_units,display_units,dim_names,flag)
            if(mutot_ptr .eq. 0) then
               write(6,*)'MU_TOT failed',flag
               stop 1
            endif
         endif
      endif
      mu_delta=delta(:) ; mu_dims = dims(:)  
      phys_mumin = min(:) ; phys_mumax = max(:) ;
      

      flag=.true.
      phtot_ptr = getvar('PH_TOT',ndims,dims,stagger,min,max,
     >     missing,data_units,display_units,dim_names,flag)
      if(phtot_ptr .eq. 0) then
         phtot_ptr=addvar('PH+PHB','PH_TOT',ndims,dims,stagger,min,max,
     >        missing,data_units,display_units,dim_names,flag)
         if(phtot_ptr .eq. 0) then
            call sleep(1)
            phtot_ptr = addvar('PH+PHB','PH_TOT',ndims,dims,stagger,min,
     >           max,missing,data_units,display_units,dim_names,flag)
            if(phtot_ptr .eq. 0) then
               write(6,*)'PHTOT failed',flag
               stop 1
            endif
         endif
      endif
      nz = dims(3)
      phdims=dims
      flag=.true.
      ptot_ptr = getvar('Pressure',ndims,dims,stagger,min,max,
     >     missing,data_units,display_units,dim_names,flag)
      if(ptot_ptr .eq. 0) then
         ptot_ptr = addvar('P+PB','Pressure',ndims,dims,stagger,min,max,
     >        missing,data_units,display_units,dim_names,flag)
         if(ptot_ptr .eq. 0) then
            call sleep(1)
            ptot_ptr = addvar('P+PB','Pressure',ndims,dims,stagger,min,
     >           max,missing,data_units,display_units,dim_names,flag)
            if(ptot_ptr .eq. 0) then
               write(6,*)'Pressure failed',flag
               stop 1
            endif
         endif
      endif
      ierr = set_button_name ('MU', '')  
      ierr = set_button_name ('MUB', '')
      ierr = set_button_name ('PH', '')
      ierr = set_button_name ('PHB', '')
      
      mapM_ptr = getvar('MAPFAC_M',ndims,dims,stagger,min,max,missing,
     >                  data_units,display_units,dim_names,flag)
      if(mapM_ptr .eq. 0) then
         call sleep(1)
         mapM_ptr = getvar('MAPFAC_M',ndims,dims,stagger,min,max,
     >        missing,data_units,display_units,dim_names,flag)
         if(mapM_ptr .eq. 0) then
            write(6,*)'MAPFAC_M failed',flag
            stop 1
         endif
      endif
         
      mapU_ptr = getvar('MAPFAC_U',ndims,dims,stagger,min,max,missing,
     >                  data_units,display_units,dim_names,flag)
      if(mapU_ptr .eq. 0) then
         call sleep(1)
         mapU_ptr = getvar('MAPFAC_U',ndims,dims,stagger,min,max,
     >        missing,data_units,display_units,dim_names,flag)
         if(mapU_ptr .eq. 0) then
            write(6,*)'MAPFAC_U failed',flag
            stop 1
         endif
      endif
     
      mapV_ptr = getvar('MAPFAC_V',ndims,dims,stagger,min,max,missing,
     >                  data_units,display_units,dim_names,flag)
      if(mapV_ptr .eq. 0) then
         call sleep(1)
         mapV_ptr = getvar('MAPFAC_V',ndims,dims,stagger,min,max,
     >        missing,data_units,display_units,dim_names,flag)
         if(mapV_ptr .eq. 0) then
            write(6,*)'MAPFAC_V failed',flag
            stop 1
         endif
      endif

      twod = 0
      if (plmax(2).eq.plmin(2)) twod = 1
      
      dx = abs(delta(1)) ; nx = nint((plmax(1)-plmin(1))/dx)
      if (twod.ne.1) then
         dy = abs(delta(2)) ; ny = nint((plmax(2)-plmin(2))/dy)
         if(ny .le. 5) then
            dy=1.; ny=1; twod=1
         endif
      else
         dy = 1. ; ny = 1
      end if
      
      nxp1 = nx + 1 ; nyp1 = ny + 1
      if(coord(1) .ne. 0)call freemem(coord(1))
      if(coord(2) .ne. 0)call freemem(coord(2))
      coord(1) = getmem(nx)
      coord(2) = getmem(ny)
      do i=1,nx
         call s_val(%val(coord(1)),i,plmin(1)+(dx*(i-1)))
      enddo
      do i=1,ny
         call s_val(%val(coord(2)),i,plmin(2)+(dy*(i-1)))
      enddo
      coord(4) = dimvar_data(4)
      exact_coord(4) = coord(4)
        n(1)=0
        n(2)=0
        n(3)=0
        n(4) = dimvar_size(4)
      mumax = 0.0
      mumin = 999999.0
      htop = 0.0
      hbot = 100000.0
      do i=1,nx*ny*ntime
         mumax=amax1(r_val(%val(mutot_ptr),i),mumax)
         mumin=amin1(r_val(%val(mutot_ptr),i),mumin)
      end do
      do i=1,phdims(1)*phdims(2)*phdims(3)*phdims(4)
         junk=r_val(%val(phtot_ptr),i)
         htop = amax1(junk,htop)
         hbot = amin1(junk,hbot)
      enddo
      htop = int(htop/9.80665)
      hbot = int(hbot/9.80665)
      if(hbot .le. .01)hbot=0.0
!       plmax(3) = ptop ; plmin(3) = (mumax+ptop)
!       phmax(3) = ptop ; phmin(3) = (mumax+ptop)
!     Z2 - Z1 = RT/G0ln(p1/p2)
!      htop = (29.3*216)*log(1000/ptop)
       plmax(3) = htop ; plmin(3) = hbot
       phmax(3) = htop ; phmin(3) = hbot
       max(3) = htop
       min(3) = hbot
      dz = (max(3)-min(3))/nz

      call setrarr('phmax',phmax,4,error1)
      call setrarr('phmin',phmin,4,error2)
      if (error1 .or. error2) then
        write(6,*) 'Error:  could not set phmax or phmax in new_file.'
        stop
      end if

      call getaarr('domain_units',domain_units,4,error1)
      call setrarr('plmax',plmax,4,error1)
      call setrarr('plmin',plmin,4,error2)
      if (error1 .or. error2) then
        write(6,*) 'Error:  could not set plmax or plmax in new_file.'
        stop
      end if
      if (error1 .or. error2) then
         write(6,*) 
     &       'Error:  could not set scaled plmax or plmax in new_file.'
        stop
      end if

      call getaarr('domain_units',domain_units,4,error1)
      domain_units(3) = domain_units(2)
      call setaarr('domain_units',domain_units,4,error2)
      if (error1 .or. error2) then
        write(6,*) 'Error:  could not change the domain_units.'
        stop
      end if

      call getaarr('domain_display_units',domain_display_units,4,error1)
      domain_display_units(3) =  domain_display_units(2)
      call setaarr('domain_display_units',domain_display_units,4,error2)
      if (error1 .or. error2) then
        write(6,*) 'Error:  could not change the domain_display_units.'
        stop
      end if

      call getaarr('dim_display_names',dim_display_names,4,error1)
      dim_display_names(3) = 'Z'
      call setaarr('dim_display_names',dim_display_names,4,error1)
      if (error1 .or. error2) then
        write(6,*) 'Error:  could not change the dim_display_names.'
        stop
      end if


      ierr=nf_get_att_double(ncid,0,'TRUELAT1',djunk(1))
      if(ierr.eq.NF_NOERR) then
         trulat1=djunk(1)
      else
         trulat1 =30.0
      endif
      ierr=nf_get_att_double(ncid,0,'TRUELAT2',djunk(1))
      if(ierr.eq.NF_NOERR) then
         trulat2=djunk(1)
      else
         trulat2 =30.0
      endif
      ierr=nf_get_att_double(ncid,0,'STAND_LON',djunk(1))
      if(ierr.eq.NF_NOERR) then
      ierr=nf_get_att_double(ncid,0,'CEN_LON',djunk(1))
      endif
      if(ierr.eq.NF_NOERR) then
         cenlon=djunk(1)
      else
         cenlon =60.0
      endif
      ierr=nf_get_att_double(ncid,0,'CEN_LAT',djunk(1))
      if(ierr.eq.NF_NOERR) then
         cenlat=djunk(1)
      else
         cenlon =60.0
      endif
      ierr=nf_get_att_double(ncid,0,'MAP_PROJ',djunk(1))
      if(ierr.eq.NF_NOERR) then
         iproj=djunk(1)
      else
         iproj = 4
      endif

      if(iproj .eq. 0) then
         projection='CE'
      else if(iproj .eq. 1) then
         projection='LC'
      else if(iproj .eq. 2) then
         projection='ST'
      else if(iproj .eq. 3) then
         projection='ME'
      endif

c      end if
c
c compute metric arrays (subroutine in metrics.F)
c
!      dqdx_xy = getmem(3*(nx+1)*(ny+1))
!      dqdx_zt = getmem(nz+1)
!      dqdx_zw = getmem(nz+1)
!      dqdy_xy = getmem(3*(nx+1)*(ny+1))
!      dqdy_zt = getmem(nz+1)
!      dqdy_zw = getmem(nz+1)
!      dqdz = getmem(3*(nx+1)*(ny+1))
!      dzdq = getmem(3*(nx+1)*(ny+1))
!      call metric_calc(%val(hgt_ptr),%val(hgt_ptr),%val(hgt_ptr),
!     >     %val(znu_ptr),%val(znw_ptr),%val(dqdx_xy),%val(dqdx_zt),
!     >     %val(dqdx_zw),%val(dqdy_xy),%val(dqdy_zt),
!     >     %val(dqdy_zw),%val(dqdz),%val(dzdq),ztop,dx,dy,
!     >     nx,ny,nz)

c set domain offsets to zero
      offset(1) = 0. ;  offset(2) = 0. ;  offset(4) = 0.

      return
      end subroutine

c
c     -----------------------------------------------------------------
c
      subroutine new_field(name, field, inmax1, inmax2, inmax3, inmax4)

      use DEFAULT_MODULE
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
      character*(*) name
      integer inmax1,inmax2,inmax3,inmax4
      real field(inmax1,inmax2,inmax3,inmax4)

      character*128 dnames(4)
      character*81 data_units,display_units,dim_names(4)
      real missing,stagger(4),start,r_val,d,dom_min(4),dom_max(4)
      integer i,j,min(4),max(4),dims(4),nd,error,ndims
      integer k,iflg, strbeg,strend
      logical flag,requ
      real :: plmin_scaled(4), plmax_scaled(4)
      integer(ive_ptr) :: dimvar_data(4)
      integer :: dimvar_size(4)
c
c save the physical and index limits in the common block in
c   default.icl for use in later calculations by other transform
c   routines.
c
!      if(len(name) .lt.81) then
!         write(6,*)'len(name)',len(name)
!         write(6,*)'NEW_FIELD name = ',trim(name(1:len(name)))
!      else
!         write(6,*)'len(name',len(name)
!         write(6,*)'NEW_FIELD name = ',name(1:80)
!      endif
      call setrvar('spval',-10.0**27,error)
      call setlvar('reset_spval', .true ., error)
      if(trim(name(1:len(name))) .eq. 'ZNW' .or. 
     &     trim(name(1:len(name))) .eq. 'ZNU' .or.
     &     trim(name(1:len(name))) .eq. 'XLAT' .or. 
     &     trim(name(1:len(name))) .eq. 'XLONG' .or.
     &     trim(name(1:len(name))) .eq. 'MU' .or.
     &     trim(name(1:len(name))) .eq. 'MUB' .or.
     &     trim(name(1:len(name))) .eq. 'PH' .or.
     &     trim(name(1:len(name))) .eq. 'PHB' 
     &     ) return

      call getrarr("plmax_scaled",plmax_scaled,4,error)
      if (error.ne.0) then
        write(6,*) 'Error:  problem getting physical domain ',
     >                               'for field.'
        stop
      end if

      call getrarr("plmin_scaled",plmin_scaled,4,error)
      if (error.ne.0) then
        write(6,*) 'Error:  problem getting physical domain ',
     >                               'for field.'
        stop
      end if

c      print *,'PLMIN/MAX_SCLAED: ', plmin_scaled, plmax_scaled

      call getrarr("phmin",phmin,4,error)
      if (error.ne.0) then
        write(6,*) 'Error:  problem getting physical domain ',
     >                               'for field.'
        stop
      end if

      call getrarr("phmax",phmax,4,error)
      if (error.ne.0) then
        write(6,*) 'Error:  problem getting physical domain ',
     >                               'maximum for field.'
        stop
      end if


      inmax(1) = inmax1 ; inmax(2) = inmax2 
      inmax(3) = inmax3 ; inmax(4) = inmax4
      if(trim(name(1:len(name))) .eq. 'HOR_MASSFLX' ) inmax(3)=1
      if(trim(name(1:len(name))) .eq. 'VER_MASSFLX' ) inmax(1)=1

c      print *,'INMAX: ',inmax
!      print *,'PLMIN: ',PLMIN
!      print *,'PLMAX: ',PLMAX

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

      if (inmax4.ne.1.and.abs(phmin(4)-plmin(4)).gt.0.0) then
        phmin(4) = phmin(4) + offset(4)
        phmax(4) = phmax(4) + offset(4)
      end if

        trans_eta = 1
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


      coord(3) = 0
c
c check to see if time coordinate vector is evenly spaced
c
      if (inmax(4).ne.1) then
         call getparr('dimvar_data', dimvar_data, 4, error)
          coord(4) = dimvar_data(4) !should = time
        k = 1
        start = r_val(%val(coord(4)),1)
        do while(k.le.inmax(4))
          if (requ(start+(k-1)*delta(4),r_val(%val(coord(4)),k))) then
            k = k+1
          else
            delta(4) = 0.
            k = inmax(4)+1
          end if 
        end do
      end if
       
c      print *,' '

      return
      end subroutine
c
c     -----------------------------------------------------------------
c
      subroutine index_2_phys(phys, index, flag, ndims, npts)

      use DEFAULT_MODULE

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

      integer npts, ndims, flag(ndims)
      real    phys(ndims, npts), index(ndims, npts), myind(ndims)

      integer i, j, k, iflag(4)
      logical error
      real    mupt, coord_k, coord_kp1, ind2(ndims), mh, spval
      real junk, geo_to_geom
      real    r_val ! function
      double precision extrap

      data iflag/0, 0, 0, 1/
      
      call getrvar('spval', spval, error)
      if(error)spval = 0.0
      do i = 1, npts
         do j = 1, ndims
            if (flag(j) .eq. 1) then
c     Handle vertical transformation
               if ( j.eq.3 ) then
                  myind=index(:,i)
!                  myind(3)=myind(3)+1
                  phys(j, i) = extrap(%val(phtot_ptr),phdims(1),
     &                 phdims(2),phdims(3),phdims(4),myind(1),
     &                 myind(2),myind(3),myind(4),spval)
                  phys(j, i)=phys(j, i)/9.80665
                  write(6,*)'phys(3)',phys(j, i), myind(3)
c     Case where we are out of the data array
               else if ( index(j,i) .lt. 1 ) then
                  phys(j,i) = spval !missing
c     Trivial case where dimension does not exist (depth of array is 
c     unity along this dimension)
               else if (inmax(j) .eq. 1 .or. index(j,i) .lt. 1) then
                  phys(j,i) = phmin(j)

c     Case where interpolation must be performed between the grid
c     point locations given in a 1D "dimension variable" array.
               else if (delta(j) .eq. 0.) then

                  k = index(j,i)
                  if (k .lt. inmax(j)) then
                     coord_k   = r_val(%val(coord(j)),k)
                     coord_kp1 = r_val(%val(coord(j)),k+1)
                     phys(j,i) = coord_k  +
     &                          (index(j,i) - k)*(coord_kp1 - coord_k)
                  else
                     phys(j,i) = phmax(j)
                  endif

c     Case of evenly spaced physical coord.
               else 
                  phys(j,i) = phmin(j) + (index(j,i)-1)*delta(j)
               endif

            endif
         enddo
      enddo
      
      return
      end subroutine
c
c     -----------------------------------------------------------------
c
      subroutine phys_2_index(phys, index, flag, ndims, npts)

      use DEFAULT_MODULE
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

      integer npts, ndims, flag(ndims)
      real    phys(ndims, npts), index(ndims, npts)

      integer i, j, k, l, myflag(4)
      logical error
      real    r_val, spval
      real    mupt, p, mh, myin(4), myph(4)
      double precision extrap
      data myflag /1,1,0,1/
      call getrvar('spval', spval, error)
      if(error)spval = 0.0
!     can't do without all others
      if(flag(3) .eq. 1) then
         call phys_2_index(phys, index, myflag, ndims, npts)
         if(inmax(3) .eq. 1) then
            index(3,:) = 1
         else
            do i = 1, npts
               mh = extrap(%val(hgt_ptr),nx,ny,1,height_t,
     &              index(1,i),index(2,i),1.,index(4,i),spval)
               if(phys(3, i) .lt. mh .or. phys(3, i) .lt. hbot) then
                  index(3,i)=-1 !under ground
               else if(phys(3, i) .ge. htop) then
                  index(3,i)=nz+1 !in space
               else
                  myin=index(:,i)
                  call hgt_2_index(phys(3,i),myin,4,spval)
                  index(:,i)=myin
               endif
            enddo
         endif
!!!!  Don't care about height
      else
         do i = 1, npts
            do j = ndims, 1, -1
               if (flag(j) .eq. 1) then
                  p = phys(j, i)
                  if(p .le. phmin(j)) then
                     index(j,i) =1
                  elseif(p .ge. phmax(j)) then
                     index(j,i) = inmax(j)
                     
c     Trivial case where dimension does not exist (depth of array is 
c     unity along this dimension)
                  elseif(inmax(j) .eq. 1) then
                     index(j,i) = 1.
                     
c     Case where interpolation must be performed between the grid
c     point locations given in a 1D "dimension variable" array.
                  else if (delta(j) .eq. 0.) then
c     
c     Physical coord. increasing with index.
                     if (phmin(j) .lt. phmax(j)) then
                        k = 2
 30                     if ((p.gt.r_val(%val(coord(j)),k)).and.
     >                       (k.lt.inmax(j))) then
                           k = k+1
                           goto 30
                        end if
                        
                        index(j,i) = real(k) - (p-
     &                       r_val(%val(coord(j)),k))/
     &                       (r_val(%val(coord(j)),k-1)-
     &                       r_val(%val(coord(j)),k))
                        
                     else       ! Physical coord. decreasing with index.
                        k = 2
 40                     if ((p.lt.r_val(%val(coord(j)),k)).and.
     >                       (k.lt.inmax(j))) then
                           k = k+1
                           goto 40
                        end if
                        
                        index(j,i) = real(k) - (p-
     &                       r_val(%val(coord(j)),k))/
     &                       (r_val(%val(coord(j)),k-1)-
     &                       r_val(%val(coord(j)),k))
                     endif
                     
                  else          ! Case of evenly spaced physical coord.
                     index(j,i) = (p-phmin(j))/delta(j) + 1
                  endif
               endif
            enddo
         enddo
      endif
      return
      end subroutine
c
c     -----------------------------------------------------------------
c
      subroutine horiz_ter(topo, nxw, nyw, stagi, stagj, zero, error)

      use DEFAULT_MODULE
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

      integer nxw, nyw, error, i, j, iflag(4)
      real    topo(nxw, nyw), stagi(nxw), stagj(nyw), zero
      real    compt(4), phypt(4)
      data iflag/1, 1, 0, 0/

      error = 0
      if (hgt_ptr .eq. 0) then
         topo(:,:) = zero
      else
         do j = 1, nyw
            do i = 1, nxw
c
c     Find the location of (i, j) in physical space. Note that
c     cpmpxy is an IVE routine and not part of the transform
c     package
c
               call cpmpxy(2, real(i), real(j), phypt(1), phypt(2))
c
c     Convert to a point in the data grid.
c
               call phys_2_index(phypt, compt, iflag, 4, 1)
c
c     Interpolate terrain to the required point.
c
               topo(i, j) =  get_ter_pt (%val(hgt_ptr), compt(1), 
     &              compt(2), zbot_dims(1), zbot_dims(2), zbot_min,
     &              zbot_delta)
               topo(i, j) = amax1 (topo(i, j), zero)
            enddo
         enddo
      endif

      end subroutine
c
c     -----------------------------------------------------------------
c
      subroutine vert_ter (xter, yter, topoht, nter, pt1, pt2,
     &                  zero)

      use DEFAULT_MODULE
c-----------------------------------------------------------------------
c     Purpose:
c
c        This routine finds the terrain height and location under each
c        point in a vertical slab.
c
c     Arguments:
c	xter	integer(ive_ptr)	output	Pointer to array of x-locations for topo_ht
c	yter	integer(ive_ptr)	output	Pointer to array of y-locations for topo_ht
c	topoht	integer(ive_ptr)	output	Pointer to array containing 
c				the terrain heights.
c	nter	int	output	Number of points in topo_ht (output)
c	pt1	real	output	Minimum physical coordinates for slice.
c	pt2	real	output	Maximum physical coordinates for slice.
c	zero	real	input	the location of the bottom of the plot
c                               window. This value is used as the minimum
c                               terrain height, rather than 0.0.
c-----------------------------------------------------------------------
c
c
c     Argument declarations.
c
      integer   nter
      real      pt1(4), pt2(4), zero
      integer(ive_ptr)   xter, yter, topoht
c
c
c     Local variable declarations.
c
      integer MAXDIM
      parameter (MAXDIM = 4)
      integer   i, iflag
      logical error
      real      phypt, min, max, spval, mytime
      dimension iflag(MAXDIM), phypt(MAXDIM),
     &          min(MAXDIM), max(MAXDIM)
c
      integer stepby
      real    value
      double precision extrap

c
      call getrvar('spval', spval, error)
      if(error) call write_message('cant get spval in vert_ter')
      if ( hgt_ptr .ne. 0 ) then
         iflag(1) = 1
         iflag(2) = 1
         iflag(3) = 0
         iflag(4) = 1
         call phys_2_index(pt1, min, iflag, MAXDIM, 1)
         call phys_2_index(pt2, max, iflag, MAXDIM, 1)
         if (max(1) .eq. min(1)) then
            stepby = 2
         else if (max(2) .eq. min(2)) then
            stepby = 1
         else
            stepby = 0
         endif
         if(min(4) .lt. height_t) then
            mytime = min(4)
         else
            mytime= height_t
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

            value = extrap(%val(hgt_ptr),nx,ny,1,height_t,
     &           min(1),min(2),1.,mytime,spval)
c            value =  get_ter_pt (%val(hgt_ptr), min(1), min(2),
c     &              zbot_dims(1), zbot_dims(2), zbot_min,
c     &              zbot_delta)
            if(value .eq. spval)value=0
c            value = amax1(value, zero)
            call put_ter(%val(xter), %val(yter), %val(topoht), 1, pt1,
     &           value)
            do i = 2, nter-1
               min(stepby) = int(min(stepby)+1)
               call index_2_phys(phypt, min, iflag, MAXDIM, 1)
               value = extrap(%val(hgt_ptr),nx,ny,1,height_t,
     &              min(1),min(2),1.,mytime,spval)
!               value =  get_ter_pt (%val(hgt_ptr), min(1), min(2),
!     &              zbot_dims(1), zbot_dims(2), zbot_min,
!     &              zbot_delta)
               if(value .eq. spval)value=0
c               value = amax1(value, zero)
               call put_ter(%val(xter), %val(yter), %val(topoht), i,
     &              phypt, value)
            enddo
            value = extrap(%val(hgt_ptr),nx,ny,1,height_t,
     &           max(1),max(2),1.,mytime,spval)
!            value =  get_ter_pt (%val(hgt_ptr), max(1), max(2),
!     &           zbot_dims(1), zbot_dims(2), zbot_min,
!           &           zbot_delta)
            if(value .eq. spval)value=0
c            value = amax1(value, zero)
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
            value = extrap(%val(hgt_ptr),nx,ny,1,height_t,
     &           min(1),min(2),1.,mytime,spval)
!            value =  get_ter_pt (%val(hgt_ptr), min(1), min(2),
!     &           zbot_dims(1), zbot_dims(2), zbot_min,
!           &           zbot_delta)
            if(value .eq. spval)value=101325.
c            value = amax1(value, zero)
            call put_ter(%val(xter), %val(yter), %val(topoht), 1, pt1,
     &           value)
            do i = 2, nter
               phypt(1) = pt1(1) + i * (pt2(1) - pt1(1)) / real(nter)
               phypt(2) = pt1(2) + i * (pt2(2) - pt1(2)) / real(nter)
               call phys_2_index(phypt, min, iflag, MAXDIM, 1)
            value = extrap(%val(hgt_ptr),nx,ny,1,height_t,
     &           min(1),min(2),1.,mytime,spval)
!               value =  get_ter_pt (%val(hgt_ptr), min(1), min(2),
!     &              zbot_dims(1), zbot_dims(2), zbot_min,
!     &              zbot_delta)
            if(value .eq. spval)value=101325.
c            value = amax1(value, zero)
            call put_ter(%val(xter), %val(yter), %val(topoht), i,
     &           phypt, value)
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
      end subroutine
c
      subroutine put_ter(xter, yter, topoht, i, pt, value)
c
c     Load the terrain arrays
c
      integer i
      real xter(i), yter(i), topoht(i), pt(4), value
c
      xter(i) = pt(1)
      yter(i) = pt(2)
      topoht(i) = value
!      write(6,*) xter(i), yter(i), topoht(i)
c
      return
      end subroutine
c
c     -----------------------------------------------------------------
c
      subroutine phys_2_lonlat(x, y, lon, lat, npts)

      use DEFAULT_MODULE
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
c
c
      integer i
      logical error
      real grdi(npts),grdj(npts), mytime
      double precision interp
c
c
c      write(6,*)'in phys_2_latlon'
      call getrvar('tloc',mytime,error)
      if(error)mytime=1.0

      do i=1,npts
         grdi(i) =  x(i)/dx
         grdj(i) =  y(i)/dx
         lat(i) = interp(%val(xlat), nx, ny, 1, ntime,
     &        grdi(i),grdj(i),1,mytime)
      lon(i) = interp(%val(xlong), nx, ny, 1, ntime,
     &        grdi(i),grdj(i),1,mytime)
      enddo
c	 
      return
      end subroutine
c
c     -----------------------------------------------------------------
c
c
      subroutine lonlat_2_phys(x, y, lon, lat, npts)
      use DEFAULT_MODULE

c     
c     lonlat_2_phys : This routine converts longitude-latitude coordinates
c     to physical coordinates.
c     
c     Arguments:
c     x		real		Physical x coordinate array (output).
c     y		real		Physical y coordinate array (output).
c     lon	real		Longitude array (input)
c     lat	real		Latitude array (input)
c     npts	integer		Number of points to convert
c
      integer npts
      real    x(npts), y(npts), lon(npts), lat(npts)
c
      integer k
      integer iflag(4)
      logical error
      real mytime
      real phys(4,npts),grid(4,npts)
c
!      data iflag/1,1,0,0/

c
c
c     Only thing we know is that we increase in x and y
c
      call getrvar('tloc',mytime,error)
      grid(4,:)=mytime
      grid(3,:)=1.
      if(error)mytime=1.0
      do k = 1, npts
         call ll2ij(grid(1,x),grid(2,k),lat(k),lon(k),%val(xlat),
     &        %val(xlong), nx, ny)
      enddo
c
      call index_2_phys(phys,grid,iflag,4,npts)
c
      do k=1,npts
        x(k) = phys(1,k)
        y(k) = phys(2,k)
      end do
c
      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine physuv_2_lonlatuv_trans(udat, vdat, nu, nv, imap)

      integer, intent(in) :: nu, nv, imap
      real, intent(inout) :: udat(nu,nv), vdat(nu,nv)

      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      use DEFAULT_MODULE
      character*(*) proj, limit
      integer       exact_fit
      real          plon, plat, rota, plm1(2), plm2(2), plm3(2), plm4(2)

      real mytime, r_val
      integer imytime
      logical error
c
c     The default is the current window in a Lambert Conic
c     Projection.
c
      proj(1:len(proj))=' '
      limit(1:len(limit))=' '
      proj(1:2) = projection(1:2)
      if(projection .eq. 'LC') then
         rota = trulat2
         plon = cenlon
         plat = trulat1
      else if(projection .eq. 'ST') then
!         plat=sign(90.,cenlon)
         plat=trulat1
         plon=cenlon
!         rota=-1*rrota(ipl)
         rota=0.
      else
         plat=0.
         rota=0.
         plon = cenlon
      endif
c
c     Setting limt to CL tells IVE to clip to the current window using
c     phys_2_lonlat.
c
      call getrvar('time', mytime, error)
      if(.not. error) then
         do imytime=1,ntime
            if(mytime .eq. r_val(%val(time),imytime)) then
               plm1=r_val(%val(lat_ll),imytime)
               plm2=r_val(%val(lon_ll),imytime)
               plm3=r_val(%val(lat_ur),imytime)
               plm4=r_val(%val(lon_ur),imytime)
            endif
         enddo
      else
         plm1=0.
         plm2=0.
         plm3=0.
         plm4=0.
      endif
      limit(1:2) = 'CO'
      exact_fit = 1
      !exact_fit = 0
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
c     a integer(ive_ptr), this function can be used to dereference the pointer
c     and extract a value from the corresponding array as follows:
c
c     If p is the integer(ive_ptr) to some data array X, you get the j th data 
c     value in X (i.e. X(j)) with the following call:
c
c     r_val(%val(p), j)
c
c     Arguments:
c     array	real	input	The array to extract the value from.
c     j		integer	input	The index of array to use.
c
      real array
      integer j
      dimension array(j)

      r_val = array(j)
      return
      end function
c
      real function r4_val(array, i,j,k,l,ni,nj,nk,nl)
c
c     This routine is used to return the value stored in an 4d real array.
c     Since Fortran cannot be used to directly access memory given
c     a integer(ive_ptr), this function can be used to dereference the pointer
c     and extract a value from the corresponding array as follows:
c
c     If p is the integer(ive_ptr) to some data array X, you get the j th data 
c     value in X (i.e. X(j)) with the following call:
c
c     r4_val(%val(p), i,j,k,l)
c
c     Arguments:
c     array	real	input	The array to extract the value from.
c     j		integer	input	The index of array to use.
c
      real array
      integer i,j,k,l,ni,nj,nk,nl
      dimension array(ni,nj,nk,nl)

      r4_val = array(i,j,k,l)
      return
      end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine s_val(array, j, val)
!
!  This routine is used to store a value in a 1d real array.
!  Since Fortran cannot be used to directly access memory given
!  a integer(ive_ptr), this function can be used to dereference the
!  pointer and extract a value from the corresponding array as follows:
!
!  If p is the integer(ive_ptr) to some data array X, you set the j the
!  data value in X (i.e. X(j)) to val with the following call:
!
!  call s_val(%val(p), j, val)
!
!  Arguments:
!     array	real	output	The array to store the value in.
!     j		integer	input	The index of array to use.
!     val	real	input	The value to set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine s_val(array, j, val)

      integer j
      real array(j), val

      array(j) = val
      return
      end
c
c  ***********************************************************
c
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
c
c   ***************************************************
c
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
c
c  *****************************************************
c
      real function get_2Dpt_time(dat,xx,yy,tt,nx,ny,nt,dat_min,delta)

      integer, intent(in) :: nx,ny,nt
      real, intent(in) :: dat(nx,ny,nt), dat_min(2), delta(2)
     &                   ,xx, yy, tt

      integer i,j,n,ip1,jp1,np1,itmp,jtmp
      real ri,rj,deli,delj,deln
      real get_2Dpt_time_np1

c      print *,'xx,yy,tt: ',xx,yy,tt
c      print *,'nx,ny,nt: ',nx,ny,nt
c      print *,'dat_min: ',dat_min
c      print *,'delta: ',delta


      ri = (xx - dat_min(1))/delta(1) + 1.
      rj = (yy - dat_min(2))/delta(2) + 1.

      i = max(min(int(ri),nx-1),1)  
      j = max(min(int(rj),ny-1),1)
!      n = floor(tt)
      n = int(tt)

c      print *,'i,j,n: ',i,j,n
c      print *,'ri,rj: ',ri,rj

      deli = ri - real(i) 
      delj = rj - real(j)
      deln = tt - real(n)

c      print *,deli,delj,deln
c      print *,dat(i,j,n)

      ip1 = min(i+1,nx)  
      jp1 = min(j+1,ny) 
      np1 = min(n+1,nt)

      get_2Dpt_time = (1.-deli)*(1.-delj)*dat(i,j,n)
     >       +  deli*(1.-delj)*dat(ip1,j,n)
     >       +  delj*(1.-deli)*dat(i,jp1,n)
     >       +  deli*delj*dat(ip1,jp1,n)
   
      if(deln.ne.0.0) then

        get_2Dpt_time_np1 = (1.-deli)*(1.-delj)*dat(i,j,n)
     >         +  deli*(1.-delj)*dat(ip1,j,n)
     >         +  delj*(1.-deli)*dat(i,jp1,n)
     >         +  deli*delj*dat(ip1,jp1,n)
        
        get_2Dpt_time = (1.0 - deln)*get_2Dpt_time 
     &                + deln*get_2Dpt_time_np1
      end if

      return
      end
c
c  *****************************************************
c
      real function get_ter_pt(ter,x,y,nx,ny,ter_min,del)
c
c     ter(nx,ny)   input terrain field
c     x            input real x val (comp point)
c     y            input real y val (comp point)
c     nx           input number in x dim
c     ny           input number in x dim
c     ter_min(2)   input currently unused
c     del input(2) input currently unused
c
      integer nx,ny
      real ter(nx,ny),ter_min(2),del(2),x,y
      double precision interp

      integer i,j,ip1,jp1
      real ri,rj,deli,delj

      i=max(int(x),1)
      j=max(int(y),1)
      i=min(int(x),nx)
      j=min(int(y),ny)
      deli = x - real(i)
      delj = y - real(j)
      ip1=min(i+1,nx)
      jp1=min(j+1,ny)



      get_ter_pt = (1.-deli)*(1.-delj)*ter(i,j)
     >     +  deli*(1.-delj)*ter(ip1,j)
     >     +  delj*(1.-deli)*ter(i,jp1)
     >     +  deli*delj*ter(ip1,jp1)
     
c      write(6,*)ter(int(x),int(y)),x,y,nx,ny
      return
      end
c
c **********************************************
c
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
c
c **********************************************
c
      subroutine setstr(outstr,instr)
c
      character*(*) outstr,instr
      integer ibeg,iend,strbeg,strend
c
      outstr(1:len(outstr)) = ' '
      ibeg = strbeg(instr)
      iend = strend(instr)
      outstr(1:iend-ibeg+2) = instr(ibeg:iend)//char(0)
c
      return
      end subroutine
c
      subroutine ref_latlon(iref,jref,dx,dy)
c
c  Index of reference lat/lon for lambert conic projection
c	
c grid 3 resolution is 3km, grid2 is 9 km and grid1 27 km
c
      integer ii21,jj21,ii32,jj32,irefi,jrefi,iref,jref
      real dx,dy

      data irefi/46/
      data jrefi/46/
      data ii21/28/
      data jj21/22/
      data ii32/41/
      data jj32/41/
c
      if (dx.eq.27000) then
        iref = irefi
        jref = jrefi
      elseif (dx.eq.9000) then
        iref=(irefi-ii21)*3+1
        jref=(jrefi-jj21)*3+1
      elseif (dx.eq.3000) then
        iref=(irefi-ii21)*9-(ii32-1)*3+1
        jref=(jrefi-jj21)*9-(jj32-1)*3+1
      end if
c
      return
      end subroutine


      subroutine hgt_2_index(hgt, index, ndims, spval)

      use DEFAULT_MODULE
c     
c     phys_2_index : This routine translates physical coordinates into 
c     array indices.
c     
c     Arguments:
c     phys	real	input	Array of physical coordinates.
c     index	real	output	Array of indicies.  A negative index =>
c				this point gets the missing data value
c				(e.g. is under the terrain).
c     ndims	integer	input	Number of dimensions in each point.
c     npts	integer	input	Number of points to convert.
c

      integer  ndims, flag(ndims)
      real    hgt, index(ndims)

      integer i, j, k
      real    r_val, spval
      double precision extrap
      
      integer iflag(4),count
      data iflag/1,1,0,1/
      real bottom, top, part, retpres

      k=1
      top=extrap(%val(phtot_ptr),phdims(1),phdims(2),phdims(3),
     &     phdims(4),index(1),index(2), real(k), index(4),spval)
      if(top .gt. hgt) then
          index(3)=1
          return
       endif
      do while(top .lt. hgt)
         k=k+1
         top=extrap(%val(phtot_ptr),phdims(1),phdims(2),phdims(3),
     &        phdims(4),index(1),index(2), real(k), index(4),spval)
      enddo
      if(top .eq. hgt) then
         index(3)=real(k)
         return
      elseif(k .gt. phdims(3)) then
         bottom=extrap(%val(phtot_ptr),phdims(1),phdims(2),phdims(3),
     &        phdims(4),index(1),index(2), phdims(3), index(4),spval)
         part=(hgt-bottom)/(top-bottom)
         index(3)=real(k-1)+part
         return
      else
         bottom=extrap(%val(phtot_ptr),phdims(1),phdims(2),phdims(3),
     &        phdims(4),index(1),index(2), real(k-1), index(4),spval)
         part=(hgt-bottom)/(top-bottom)
         index(3)=real(k-1)+part
         return
      endif
      end


      subroutine pres_2_index(pres, index, ndims, spval)

      use DEFAULT_MODULE
c     
c     phys_2_index : This routine translates physical coordinates into 
c     array indices.
c     
c     Arguments:
c     phys	real	input	Array of physical coordinates.
c     index	real	output	Array of indicies.  A negative index =>
c				this point gets the missing data value
c				(e.g. is under the terrain).
c     ndims	integer	input	Number of dimensions in each point.
c     npts	integer	input	Number of points to convert.
c

      integer  ndims
      real    pres, index(ndims)

      integer i, j, k
      real    r_val, spval
      double precision interp
      
      real bottom, top, part
         
      if(nz .eq. 1) then
         index(3) = 1
         return
      endif
      k=nz
      bottom=interp(%val(ptot_ptr),nx,ny,nz,ntime,
     &     index(1),index(2), real(k), index(4),spval)
      do while(bottom .lt. pres .and. k .gt. 1)
         k=k-1
         bottom=interp(%val(ptot_ptr),nx,ny,nz,ntime,
     &        index(1),index(2), real(k), index(4),spval)
      enddo
!      write(6,*)pres,bottom,k
      if((bottom .gt. pres .and. k .eq. nz) .or.
     &     ((k.eq.1) .and. bottom .lt. pres)) then
         index(3)=-1
         return
      endif

      if(bottom .eq. pres) then
         index(3)=real(k)
         return
      else
         top=interp(%val(ptot_ptr),nx,ny,nz,ntime,
     &        index(1),index(2), real(k+1), index(4),spval)
         part=(pres-bottom)/(top-bottom)
         index(3)=real(k)+part
      endif
      return
      end

         
      subroutine transcopylats(lats, x,y,z,t, lat_ll, lat_ur)
      real lats, lat_ll, lat_ur
      integer x,y,z,t
      dimension lats(x,y,z,t), lat_ll(t), lat_ur(t)

      lat_ll(:) = lats(1,1,1,:)
      lat_ur(:) = lats(x,y,z,:)
      end

      subroutine transcopylons(lons, x,y,z,t, lon_ll, lon_ur)
      real lons, lon_ll, lon_ur
      integer x,y,z,t
      dimension lons(x,y,z,t), lon_ll(t), lon_ur(t)

      lon_ll(:) = lons(1,1,1,:)
      lon_ur(:) = lons(x,y,z,:)
      end

      real function geo_to_geom(geohgt, x, y, t)
      use default_module
      real geohgt, x, y, t, somigliana,eccentricity
      real termg, termr, termrg, slat, zges, sin2, spval
      double precision interp
      
!     Convert geopotential height at layer midpoints to geometric 
!     height using equations (17, 20, 23) in MJ Mahoney's note 
!     "A discussion of various measures of altitude" (2001).  
!     Available on the web at
!     http://mtp.jpl.nasa.gov/notes/altitude/altitude.html
!     From: http://www.emc.ncep.noaa.gov/gmb/treadon/gsi/gsi_browser/html_code/gsi/setupw.f90.html
!     termg  = equation 17
!     termr  = equation 21
!     termrg = first term in the denominator of equation 23
!     zges  = equation 23
      spval=0.0
      somigliana = (6378.1370e3/6356.7523142e3)*
     &     (9.8321849378/9.7803253359) -1.0
      eccentricity=sqrt((6378.1370e3)**2 -(6356.7523142e3)**2)/
     &     6378.1370e3
      slat = interp(%val(xlat), nx,ny,1,ntime,x,y,1,t,spval)
      slat = slat * (3.14159/180.0)
      sin2  = sin(slat)*sin(slat)
      termg = 9.7803253359 *
     &     ((1.0+somigliana*sin2)/
     &     sqrt(1.0-eccentricity*eccentricity*sin2))
      termr =  6378.1370e3 /(1.0 + 
     &     (6378.1370e3-6356.7523142e3)/6378.1370e3 + 
     &     ((7.292115e-5**2 * 6378.1370e3**2 * 6356.7523142e3)/ 
     &     3.986004418e14) - 
     &     (2.0*(6378.1370e3-6356.7523142e3)/6378.1370e3)*sin2)

      termrg = (termg/9.80665)*termr
      geo_to_geom = (termr*geohgt)/(termrg-geohgt) ! eq (23)
      return
      end
      
      
