#include "pointer.icl"
    
      ive_ptr function calc_field(name,ndims,dims,stag,min,max,
     &             missing,data_units,data_display_units,dim_names)
c     
c     calc_field : This routine is used to calculate user-derived field.
c     The return value of the routine is a ive_ptr to the field values.
c     NULL return => cannot calculate.
c     
c     Arguments:
c     name	character	The name of the field to derive.
c     ndims	integer		Number of dimensions in field (output).
c     dims	integer		Number of nodes along each dimension
c				in Fortran order (nx,ny,nz,nt) (output).
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
      integer :: ndims,dims(4),i
      real    :: stag(4),min(4),max(4),missing
      character(len=*) :: name,data_units,data_display_units
     &                   ,dim_names(4)
c
      integer ibeg,iend,idvar,ircode,strbeg,strend,ierr,ncvid
      integer iarg(10),k,rcnt,icnt,ccnt
      real rarg(10),r_val
      character(len=80) :: fld,arg1,arg2,arg3,carg(10)
     &                    ,str1,tstfld,mkfield,fld_cut
      character(len=1) :: comp
      logical isfunc,err,new,top
c
      common/t_cdfid/ idcdf
      integer idcdf
c
      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2),
     &           t_dims(1,10,2),numout(10)
      character (len=80) :: flnam(10)

      real x,y,lat,lon
      integer npts
c
      ive_ptr vert_grid,deriv,urot,vrot,buoy_freq,massflx_lat
      ive_ptr massflx_top,massflx_max,rho_calc,omega,test
      ive_ptr interp_pres,geopot,x_grid,y_grid,deriv_nsew
      ive_ptr geowind,urot2,windspeed,reldiff,time_diff
      ive_ptr variance, maxmin_reg, time_mean, the2temp
      ive_ptr interp_height, ens_cov, kgain
c
      call setstr(fld,name)
c
      data_units(1:len(data_units)) = ' '
      data_display_units(1:len(data_display_units)) = ' '
      iend = index(fld,'[',.true.) - 1
      if(iend.lt.0) then
        fld_cut = trim(fld)
      else
c        fld_cut(1:iend) = fld(1:iend)
        fld_cut = fld(1:iend)
      end if

      select case (trim(fld_cut))
c
c  List of available fields
c
      case ('LIST')
        call list
        calc_field = 1
c
c  List of available fields
c
      case ('SNDOUT_ENKF')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err.or.(rcnt.ne.4.and.icnt.ne.1)) then
           write(*,'(2A)') 'Usage: ',
     &     'SNDOUT_ENKF[REAL IPT:REAL JPT: REAL z1:REAL Z2:INT TIME]'
           return
         else
           call sound_out(rarg(1),rarg(2),rarg(3),rarg(4),0.0,
     &                iarg(1),1,
     &                ndims,dims,stag,min,max,missing,
     &                data_units,data_display_units,dim_names)
           return
         end if

      case ('SNDOUT_COAMPS_IDEAL')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err.or.(rcnt.ne.4.and.icnt.ne.1)) then
           write(*,'(2A)') 'Usage: ',
     &     'SNDOUT_COAMPS_IDEAL[REAL IPT:REAL JPT: ' //
     &     'REAL z1:INT TIME:REAL ROTANG]'
           return
         else
           call sound_out(rarg(1),rarg(2),rarg(3),0.0,rarg(4),
     &                iarg(1),2,
     &                ndims,dims,stag,min,max,missing,
     &                data_units,data_display_units,dim_names)
           return
         end if
c
c
c  Grid heights
c
      case ('GRID')
       call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
       calc_field = vert_grid(carg(1),ndims,dims,stag,min,max,missing,
     &                 data_units,data_display_units,dim_names)
c      else if (isfunc(fld,'GRID',arg1,arg2,arg3)) then
c        calc_field = vert_grid(arg1,ndims,dims,stag,min,max,missing,
c     &                 data_units,data_display_units,dim_names)
c
c      else if (isfunc(fld,'XGRID',arg1,arg2,arg3)) then
c        calc_field = x_grid(arg1,ndims,dims,stag,min,max,missing,
c     &                 data_units,data_display_units,dim_names)
c
c      else if (isfunc(fld,'YGRID',arg1,arg2,arg3)) then
c        calc_field = y_grid(arg1,ndims,dims,stag,min,max,missing,
c     &                 data_units,data_display_units,dim_names)
c
c  Translate domain
c
      case ('TRANSLATE')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (rcnt.ne.2) then
          write(6,*) 'Error:  expected two real arguments.'
          return
        else
          call translate_domain(rarg(1),rarg(2))
          calc_field = 1
          return
        end if
c
c  Turn vertical coordinate transformation on or off
c
      case ('TRANS_ON')
        call trans_switch(1)
        calc_field = 1
        return

      case ('TRANS_OFF')
        call trans_switch(0)
        calc_field = 1
        return
c
c  Test field
c
      case ('TEST')
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c          return
c        else
c          tstfld = 'test'
c          WRITE(*,*) mkfield(tstfld,carg,iarg,rarg,ccnt,icnt,rcnt,4)
        x = 122850
        y = 193200
        npts = 1
c
          return
c        end if
c
c  Derivative (fourth order diff)
c
       case ('D')
         call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err) then
           return
         else if (ccnt.ne.2) then
           write(*,*) 'Error:  expected two string arguments.'
           return
         else
           calc_field = deriv(carg(1),carg(2),1,ndims,dims,stag,min,
     &                   max,missing,data_units,data_display_units,
     &                   dim_names)
         return
         end if
c
       case ('DLL')
         call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err) then
           return
         else if (ccnt.ne.2) then
           write(*,*) 'Error:  expected two string arguments.'
           return
         else
           calc_field = deriv_nsew(carg(1),carg(2),ndims,dims,stag,min,
     &                   max,missing,data_units,data_display_units,
     &                   dim_names)
           return
         end if
c
c  Rotate u wind through angle theta
c
      case ('UROT')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (rcnt.ne.1) then
          write(6,*) 'Error:  expected one real arguments.' 
          write(6,*) 'The angle of rotation.'
          return
        else if (ccnt.ne.2) then
          write(6,*) 'Error:  expected two string arguments.'
          write(6,*) 'U and V wind components.'
          return
        else
          calc_field = urot2(rarg(1),carg(1),carg(2),
     &                ndims,dims,stag,min,max,missing,
     &                data_units,data_display_units,dim_names)
          return
        end if
c
c  Calculate Windspeed
c
      case ('WSPD')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.2) then
          write(6,*) 'Error:  expected two character args.'
          return
        else
          calc_field = windspeed(carg(1),carg(2),ndims,dims,stag,
     &                 min,max,missing,data_units,
     &                 data_display_units,dim_names)
          return
        end if
c
c  Calculate relative diff (fld1 - fld2) / fld2
c
      case ('RELDIFF')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.2.and.rcnt.ne.1) then
          write(6,*) 'Error:  expected two character args and 1 real.'
          return
        else
          calc_field = reldiff(carg(1),carg(2),rarg(1),ndims,dims,stag,
     &                 min,max,missing,data_units,
     &                 data_display_units,dim_names)
          return
        end if
!
! Difference between time level t1 and t2
! timediff = fld(t1) - fld(t2)
!
      case ('TIMEDIFF')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.1.and.icnt.ne.2) then
          write(6,*) 'Error:  expected one character args and 2 ints.'
          return
        else
          calc_field=time_diff(carg(1),iarg(1),iarg(2),ndims,dims,stag,
     &                 min,max,missing,data_units,
     &                 data_display_units,dim_names)
          return
        end if
!
!  Calculate the time mean (assuming equally spaced time intervals)
!
      case ('TIMEMEAN')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.1) then
          write(6,*) 'USAGE: TIMEMEAN[CHAR FLD]'
          return
        else
          calc_field=time_mean(trim(carg(1)),
     &                 ndims,dims,stag,min,max,missing,data_units,
     &                 data_display_units,dim_names)
          return
        end if
c
c
c  RHO
c
      case ('RHO')
        calc_field = rho_calc(ndims,dims,stag,min,max,missing,
     &                      data_units,data_display_units,dim_names)
        return
c
c  OMEGA
c
      case ('OMEGA')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if(ccnt.ne.3) then
          print *,'ERROR: USAGE omega[U,V,W]'
          return
        else  
          calc_field = omega(carg(1),carg(2),carg(3),
     &                      ndims,dims,stag,min,max,missing,
     &                      data_units,data_display_units,dim_names)
          return
        endif
c    
c  Buoyancy frequency
c    
      case ('NSPEC')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.1) then
          write(6,*) 'Error:  expected two character args.'
          return
        else
          calc_field = buoy_freq(carg(1),
     &                   ndims,dims,stag,min,max,missing,data_units,
     &                   data_display_units,dim_names,0)
          return
        end if

      case ('N')
        calc_field = buoy_freq('the',ndims,dims,stag,min,max,missing,
     &                 data_units,data_display_units,dim_names,0)
        return

      case ('N2')
        calc_field = buoy_freq('the',ndims,dims,stag,min,max,missing,
     &                 data_units,data_display_units,dim_names,1)
        return
c
c  MASS FLUX
c
      case ('MASSFLX_LAT')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (rcnt.ne.4) then
          write(6,*) 'Error:  expected four real arguments.'
          return
        else if (icnt.ne.2) then
          write(6,*) 'Error:  expected two integer arguments.'
          return
        else
          calc_field = massflx_lat(rarg(1),rarg(3),rarg(2),rarg(4),
     &                 iarg(1),iarg(2),ndims,dims,stag,min,max,
     &                 missing,data_units,data_display_units,
     &                 dim_names)
          return
        end if
c
      case ('MASSFLX_TOP')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          write(6,*) 'Error:  Unknown'
          return
        else if (rcnt.ne.8) then
          write(6,*) 'Error:  expected eight real arguments.'
          return
        else if (icnt.ne.1) then
          write(6,*) 'Error:  expected one integer arguments.'
          return
        else
          top = .true.
          calc_field = massflx_top(rarg(1:2),rarg(3:4),rarg(5:6),
     &                 rarg(7:8),iarg(1),top,
     &                 ndims,dims,stag,min,max,missing,data_units,
     &                 data_display_units,dim_names)
        return
      end if
c
      case ('MASSFLX_BOT')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
        write(6,*) 'Error:  Unknown'
          return
        else if (rcnt.ne.8) then
          write(6,*) 'Error:  expected eight real arguments.'
          return
        else if (icnt.ne.1) then
          write(6,*) 'Error:  expected one integer arguments.'
          return
        else
          top = .false.
          calc_field = massflx_top(rarg(1:2),rarg(3:4),rarg(5:6),
     &                 rarg(7:8),iarg(1),top,
     &                 ndims,dims,stag,min,max,missing,data_units,
     &                 data_display_units,dim_names)
       return
       end if
c
      case ('MASSFLX_MAX')
          calc_field = massflx_max(
     &                 ndims,dims,stag,min,max,missing,data_units,
     &                 data_display_units,dim_names)
      return

      case ('MAX_SPEC')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
        write(6,*) 'Error:  Unknown'
          return
        else if (rcnt.ne.1) then
          write(6,*) 'Error:  expected 1 real argument.'
          return
        else if (ccnt.ne.1) then
          write(6,*) 'Error:  expected one character argument.'
          return
        else
          calc_field = maxmin_reg(trim(carg(1)),rarg(1),
     &                 ndims,dims,stag,min,max,missing,data_units,
     &                 data_display_units,dim_names)
       return
       end if
c
c  Calculates the temperature from theta and pressure
c
      case ('THE2TEMP')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err.or.ccnt.ne.2) then
           write(*,*) 'Usage: THE2TEMP[THETA FIELD:PRESSURE FIELD]'
           return
         else
           calc_field = the2temp(carg(1),carg(2),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
         return
         end if
c
c  Calculates the variance of the time dimension of a field
c
      case ('VARIANCE')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err.or.ccnt.ne.1) then
           write(*,*) 'Usage: variance[fldname]'
           return
         else
           calc_field = variance(carg(1),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
         return
         end if
c
c  Interpolate a field to a height above ground level
c
      case ('INTERPZAGL')
         call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err.or.ccnt.ne.1.or.rcnt.ne.1) then
           write(*,*) 'Usage: interp[string field: real z (agl)]'
           return
         else
           calc_field = interp_height(carg(1),rarg(1),.true.,
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
           do i=1,4
             print *,min(i),max(i)
           end do
         return
         end if
c
c
c  Interpolate a field to a presure level
c
      case ('INTERPP')
         call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err.or.ccnt.ne.1.or.rcnt.ne.1) then
           write(*,*) 'Usage: interp[string field: real plevel]'
           return
         else
           calc_field = interp_pres(carg(1),rarg(1),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
         return
         end if
c
      case ('GEOPOT')
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err.or.rcnt.ne.1) then
           write(*,*) 'Usage: geopot[real plevel]'
           return
         else
           calc_field = geopot(rarg(1),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
         return
         end if
c
c  U and V components of the geostrophic wind
c
      case ('UG')
         call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         calc_field = geowind('U',
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
         return
c
      case ('VG')
         call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         calc_field = geopot('V',
     &                ndims,dims,stag,min,max,missing,data_units,
     &                data_display_units,dim_names)
         return
c
c  Calculates the covariance between a point and a field
c
      case ('ENS_COV')
         call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err.or.ccnt.lt.2.or.ccnt.gt.3.or.rcnt.ne.3) then
           write(*,*) 'Usage: ',
     &      'ens_cov[string var1: string var2: real x0: real y0: ',
     &              'real z0 : {lmapon T|F F} ]'
           return
         else
           if(ccnt.eq.2) then
             calc_field = ens_cov(carg(1),carg(2),
     &                    rarg(1),rarg(2),rarg(3),.false.,
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
           else
             calc_field = ens_cov(carg(1),carg(2),
     &                    rarg(1),rarg(2),rarg(3),.true.,
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
           endif
           return
         endif
c
c  Calculates the covariance between a point and a field
c
      case ('KGAIN')
         call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
         if (err.or.ccnt.lt.2.or.ccnt.gt.3.or.rcnt.ne.4) then
           write(*,*) 'Usage: ',
     &      'ens_cov[string var1: string var2: real R^2: ',
     &              'real x0: real y0: real z0: ',
     &              'lmapon T|F F} ]'
           return
         else
           if(ccnt.eq.2) then
             calc_field = kgain(carg(1),carg(2),
     &                    rarg(1),rarg(2),rarg(3),rarg(4),.false.,
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
           else
             calc_field = kgain(carg(1),carg(2),
     &                    rarg(1),rarg(2),rarg(3),rarg(4),.true.,
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
           endif
           return
         endif
c
c
c  Variable not found
c
      case default

        calc_field = 0
        return

      end select
c
c  Write out any scalars
c
      if (calc_field.ne.0) then
        if ((dims(1).eq.1).and.(dims(2).eq.1).and.(dims(3).eq.1)
     &                   .and.(dims(4).eq.1)) then
          write(6,*)
          write(6,*) 'Field ',fld(1:iend),' = ',
     &                    r_val(%val(calc_field),1)
          write(6,*)
        end if
      end if

    
      return
      end function
