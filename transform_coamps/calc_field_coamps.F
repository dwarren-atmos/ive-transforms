      integer(ive_ptr) function calc_field(name,ndims,dims,stag,min,max,
     &             missing,data_units,data_display_units,dim_names)

      use default_module
      use misc_module
c     
c     calc_field : This routine is used to calculate user-derived field.
c     The return value of the routine is a integer(ive_ptr) to the field values.
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
c      type(meta_ive) :: meta
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
      logical isfunc,err,new,top,lmap
c
      common/t_cdfid/ idcdf
      integer idcdf
c
      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2),
     &           t_dims(1,10,2),numout(10)
      character (len=80) :: flnam(10)
      character (len=1) :: ctmp

      real x,y,lat,lon
      integer npts
c
      integer(ive_ptr) vert_grid,deriv,urot,vrot,buoy_freq,massflx_lat,
     &                 massflx_top,massflx_max,rho_calc,omega,test,
     &                 interp_pres,findz,x_grid,y_grid,
     &                 geowind,urot2,windspeed,reldiff,time_diff,
     &                 variance, maxmin_reg, time_mean, the2temp,
     &                 interp_height, kgain, localize, vort_z, vort_x,
     &                 vort_y,pvort,sigdiff,rmserr
c
      call setstr(fld,name)
c
      write(6,*)'du:',len(data_units),' ddu: ',len(data_display_units)
      data_units(1:len(data_units)) = ' '
      data_display_units(1:len(data_display_units)) = ' '
!      data_units = ' '
!      data_display_units = ' '
      iend = index(fld,'[',.true.) - 1

      if(iend.lt.0) then
        fld_cut = trim(fld)
      else
        fld_cut = fld(1:iend)
      end if

      call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
      if (err) return

      calc_field = 0
      select case (trim(fld_cut))
c
c  List of available fields
c
      case ('LIST')
        call list
        return
c
c  List of available fields
c
      case ('SNDOUT_ENKF')
         call calc_field_message(fld_cut,rcnt,icnt,ccnt,4,1,0,err)
         if (err) then
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
         call calc_field_message(fld_cut,rcnt,icnt,ccnt,4,1,0,err)
         if (err) then
           write(*,'(3A)') 'Usage: ',
     &     'SNDOUT_COAMPS_IDEAL[REAL IPT:REAL JPT: ',
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
       call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,1,err)
       if(err) return

       calc_field = vert_grid(carg(1),ndims,dims,stag,min,max,missing,
     &                 data_units,data_display_units,dim_names)
c
c  Translate domain
c
      case ('TRANSLATE')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,2,0,0,err)
        if(err) return

        call translate_domain(rarg(1),rarg(2))
        !calc_field = 1
        return
c
c  Translate time
c
      case ('TRANS_TIME')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,0,0,err)
        if(err) return
        call translate_time(rarg(1))

        !calc_field = 1
        return
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
c  Derivative (fourth order diff)
c
       case ('D')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
        if(err) return

        calc_field = deriv(carg(1),carg(2),1,ndims,dims,stag,min,
     &                   max,missing,data_units,data_display_units,
     &                   dim_names)
        return
c
c  Rotate u wind through angle theta
c
      case ('UROT')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,0,2,err)
        if(err) return

        calc_field = urot2(rarg(1),carg(1),carg(2),
     &                ndims,dims,stag,min,max,missing,
     &                data_units,data_display_units,dim_names)
        return
c
c  Calculate Windspeed
c
      case ('WSPD')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
        if (err) return

        calc_field = windspeed(carg(1),carg(2),ndims,dims,stag,
     &                 min,max,missing,data_units,
     &                 data_display_units,dim_names)
        return
c
c  Calculate relative diff (fld1 - fld2) / fld2
c
      case ('RELDIFF')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,0,2,err)
        if (err) return

        calc_field = reldiff(carg(1),carg(2),rarg(1),ndims,dims,stag,
     &                 min,max,missing,data_units,
     &                 data_display_units,dim_names)
        return
!
! Difference between time level t1 and t2
! timediff = fld(t1) - fld(t2)
!
      case ('TIMEDIFF')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,2,1,err)
        if (err) return
        calc_field=time_diff(carg(1),iarg(1),iarg(2),ndims,dims,stag,
     &                 min,max,missing,data_units,
     &                 data_display_units,dim_names)
        return
!
!  Calculate the time mean (assuming equally spaced time intervals)
!
      case ('TIMEMEAN')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,1,err)
        if (err) return

        calc_field=time_mean(trim(carg(1)),
     &                 ndims,dims,stag,min,max,missing,data_units,
     &                 data_display_units,dim_names)
        return
c
c  RHO
c
      case ('RHO')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
        if (err) return

        calc_field = rho_calc(carg(1),carg(2),
     &                        ndims,dims,stag,min,max,missing,
     &                        data_units,data_display_units,dim_names)
        return
c
c  OMEGA
c
      case ('OMEGA')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,3,err)
        if (err) return

        calc_field = omega(carg(1),carg(2),carg(3),
     &                      ndims,dims,stag,min,max,missing,
     &                      data_units,data_display_units,dim_names)
        return
c    
c  Buoyancy frequency
c    
      case ('NSPEC')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,1,err)
        if (err) return

        calc_field = buoy_freq(carg(1),
     &                   ndims,dims,stag,min,max,missing,data_units,
     &                   data_display_units,dim_names,0)
        return

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
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,4,2,0,err)
        if (err) return

        calc_field = massflx_lat(rarg(1),rarg(3),rarg(2),rarg(4),
     &                 iarg(1),iarg(2),ndims,dims,stag,min,max,
     &                 missing,data_units,data_display_units,
     &                 dim_names)
        return
c
      case ('MASSFLX_TOP')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,8,1,0,err)
        if (err) return

        top = .true.
        calc_field = massflx_top(rarg(1:2),rarg(3:4),rarg(5:6),
     &                 rarg(7:8),iarg(1),top,
     &                 ndims,dims,stag,min,max,missing,data_units,
     &                 data_display_units,dim_names)
        return
c
      case ('MASSFLX_BOT')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,8,1,0,err)
        if (err) return

        top = .false.
        calc_field = massflx_top(rarg(1:2),rarg(3:4),rarg(5:6),
     &                 rarg(7:8),iarg(1),top,
     &                 ndims,dims,stag,min,max,missing,data_units,
     &                 data_display_units,dim_names)
       return
c
      case ('MASSFLX_MAX')
        calc_field = massflx_max(
     &                 ndims,dims,stag,min,max,missing,data_units,
     &                 data_display_units,dim_names)
        return

      case ('MAX_SPEC')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,0,1,err)
        if(err) return

        calc_field = maxmin_reg(trim(carg(1)),rarg(1),
     &                 ndims,dims,stag,min,max,missing,data_units,
     &                 data_display_units,dim_names)
        return
c
c  Calculates the temperature from theta and pressure
c
      case ('THE2TEMP')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
        if(err) return

        calc_field = the2temp(carg(1),carg(2),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
        return
c
c  Calculates the variance of the time dimension of a field
c
      case ('VARIANCE')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,1,err)
        if(err) return

        calc_field = variance(carg(1),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
         return
c
c  Interpolate a field to a height (z0<=0 --> abs(z0) agl )
c
      case ('INTERPZ')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,0,1,err)
        if(err) return

        calc_field = interp_height(carg(1),rarg(1),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
        return
c
c
c  Interpolate a field to a presure level
c
      case ('INTERPP')
        if(ccnt.eq.1) then
          ccnt=2 ; carg(2) = 'PPP'
        end if
        if(icnt.eq.0) then
          icnt=1 ; iarg(1) = 0
        end if

        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,1,2,err)
        if(err) return

        calc_field = interp_pres(carg(1),carg(2),rarg(1),iarg(1),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
        return
c
      case ('GEOPOT')
        if(ccnt.eq.0) then
          ccnt=1 ; carg(1) = 'PPP'
        end if

        if(icnt.eq.0) then
          icnt=1 ; iarg(1) = 0
        end if

        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,1,1,err)
        if(err) return

        calc_field = findz(carg(1),rarg(1),.false.,.true.,iarg(1),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
        return

c
c  U and V components of the geostrophic wind
c
      case ('UG')
        calc_field = geowind('U',
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
        return
c
      case ('VG')
        calc_field = geowind('V',
     &                ndims,dims,stag,min,max,missing,data_units,
     &                data_display_units,dim_names)
        return
c
c  Calculates the kalman gain between a point and a field
c
      case ('KGAIN')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,4,0,2,err)
        if(err) return
        
        ctmp = 'F' ; lmap = .false.
        if(ccnt.eq.3) ctmp = carg(3)
        if(ctmp.eq.'T'.or.ctmp.eq.'t') lmap = .true.
        calc_field = kgain(carg(1),carg(2),
     &                    rarg(1),rarg(2),rarg(3),rarg(4),lmap,fld_cut,
     &                    0,0,
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
        return
      case ('KGAIN_SPEC')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,2,0,3,err)
        if(err) return
        
        ctmp = 'F' ; lmap = .false.
        if(ccnt.eq.4) ctmp = carg(4)
        if(ctmp.eq.'T'.or.ctmp.eq.'t') lmap = .true.

        call get_stnloc(rarg(3),rarg(4),carg(3))
        print '(A,F9.4,A1,F8.5,A)',
     &        'draw_point_map=',rarg(3),',',rarg(4),',dot,2' 

        calc_field = kgain(carg(1),carg(2),
     &                    rarg(1),rarg(3),rarg(4),rarg(2),lmap,fld_cut,
     &                    0,0,
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
        return
c
c  Calculates the covariance between a point and a field
c
      case ('ENS_COV')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,3,0,2,err)
        if(err) return

        ctmp = 'F' ; lmap = .false.
        if(ccnt.eq.3) ctmp = carg(3)
        if(ctmp.eq.'T'.or.ctmp.eq.'t') lmap = .true.

        calc_field = kgain(carg(1),carg(2),0.0,
     &                    rarg(1),rarg(2),rarg(3),lmap,fld_cut,0,0,
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
        return
c
c  Calculates the covariance between a point and a field
c
      case ('ENS_COR')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,3,0,2,err)
        if(err) return

        ctmp = 'F' ; lmap = .false.
        if(ccnt.eq.3) ctmp = carg(3)
        if(ctmp.eq.'T'.or.ctmp.eq.'t') lmap = .true.
        calc_field = kgain(carg(1),carg(2),0.0,
     &                    rarg(1),rarg(2),rarg(3),lmap,fld_cut,0,0,
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
        return

c
c  Calculates gaspari-coehn localization radius
c
      case ('GCLOC')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,3,0,1,err)
        if(err) return

        calc_field = localize(rarg(2),rarg(1),rarg(3),carg(1),
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)

        return

      case ('GCLOC_SPEC')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,0,2,err)
        if(err) return

        call get_stnloc(rarg(3),rarg(2),carg(1))

        calc_field = localize(rarg(2),rarg(3),rarg(1),carg(2),
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)

        return
c
c  Calculates the hierarchical localization factor
c
      case ('HIERLOC')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,4,2,2,err)
        if(err) return

        ctmp = 'F' ; lmap = .false.
        if(ccnt.eq.3) ctmp = carg(3)
        if(ctmp.eq.'T'.or.ctmp.eq.'t') lmap = .true.
        calc_field = kgain(carg(1),carg(2),
     &                    rarg(1),rarg(2),rarg(3),rarg(4),lmap,fld_cut,
     &                    iarg(1),iarg(2),
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)

        return

      case ('HIERLOC_SPEC')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,2,2,3,err)
        if(err) return

        ctmp = 'F' ; lmap = .false.
        if(ccnt.eq.4) ctmp = carg(4)
        if(ctmp.eq.'T'.or.ctmp.eq.'t') lmap = .true.

        call get_stnloc(rarg(3),rarg(4),carg(3))
        print '(A,F9.4,A1,F8.5,A)',
     &        'draw_point_map=',rarg(3),',',rarg(4),',dot,2' 

        print *, rarg
        print *, iarg
        print *, carg

        calc_field = kgain(carg(1),carg(2),
     &                    rarg(1),rarg(3),rarg(4),rarg(2),lmap,fld_cut,
     &                    iarg(1),iarg(2),
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)

        return

      case ('VORT_X') ! x component of vorticity
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
        if(err) return

        calc_field = vort_x(carg(1),carg(2),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
        return

      case ('VORT_Y') ! y component of vorticity
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
        if(err) return

        calc_field = vort_y(carg(1),carg(2),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
        return

      case ('VORT_Z') ! z component of vorticity
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
        if(err) return

        calc_field = vort_z(carg(1),carg(2),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
        return

      case ('PVORT') ! potential vorticity
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,5,err)
        if(err) then
          call make_help_widget('Expected: pvort[U:V:W:THE:PRE]')
          return
        end if

        calc_field = pvort(carg(1),carg(2),carg(3),carg(4),carg(5),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
        return

      case ('SIGDIFF') ! Significant difference between two means
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,1,4,err)
        if(err) then
          call make_help_widget('Expected:'
     &         //     'SIGDIFF[VAR1:VAR2:STD1:STD2,CONFLEV,N]')
          return
        end if

        calc_field = sigdiff(carg(1),carg(2),carg(3),carg(4),
     &                       rarg(1),iarg(1),iarg(1),
     &                       ndims,dims,stag,min,max,missing,data_units,
     &                       data_display_units,dim_names)
        return

      case ('RMSERR') ! RMS Difference between two means

        if(icnt.lt.1) then
         icnt = 1 ; iarg(1) = 0
        end if

        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,1,2,err)
        if(err) then
         call make_help_widget('Expected: RMSERR[VAR1:VAR2:[FULL_DOM]]')
         return
        end if

        calc_field = rmserr(carg(1),carg(2),iarg(1),.false.,
     &                       ndims,dims,stag,min,max,missing,data_units,
     &                       data_display_units,dim_names)
        return

      case ('ABSERR') ! ABSOLUTE Difference between two means

        if(icnt.eq.0) then
         icnt = 1 ; iarg(1) = 0
        end if

        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,1,2,err)
        if(err) then
         call make_help_widget('Expected: ABSERR[VAR1:VAR2:[FULL_DOM]]')
         return
        end if

        calc_field = rmserr(carg(1),carg(2),iarg(1),.true.,
     &                       ndims,dims,stag,min,max,missing,data_units,
     &                       data_display_units,dim_names)
        return
c
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
          write(6,*) 'Field ',trim(fld_cut),' = ',
     &                    r_val(%val(calc_field),1)
          write(6,*)
        end if
      end if

    
      return
      end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This Subroutine checks for the correct usage of the ive function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine calc_field_message(fld_cut,rcnt,icnt,ccnt
     &                             ,rcnt_exp,icnt_exp,ccnt_exp,err)

      integer, intent(in) :: rcnt,icnt,ccnt,rcnt_exp,icnt_exp,ccnt_exp
      character(len=*), intent(in) :: fld_cut
      logical, intent(inout) :: err
      
      err = .false.

      if(rcnt.ne.rcnt_exp.or.icnt.ne.icnt_exp.or.ccnt.lt.ccnt_exp) then
        print '(A,1x,A,A1)','-Incorrect usage of',trim(fld_cut),'.'
        err=.true.
      end if

      if(rcnt.ne.rcnt_exp) then
        print '(A,1x,I1,1x,A,1x,I1,A1)',
     &    '---Expecting',rcnt_exp,'real argument received',rcnt,'.'
      end if

      if(icnt.ne.icnt_exp) then
        print '(A,1x,I1,1x,A,1x,I1,A1)',
     &    '---Expecting',icnt_exp,'integer argument received',icnt,'.'
      end if

      if(ccnt.lt.ccnt_exp) then
        print '(A,1x,I1,1x,A,1x,I1,A1)',
     &    '---Expecting',ccnt_exp,'character argument received',ccnt,'.'
      end if

      end subroutine
