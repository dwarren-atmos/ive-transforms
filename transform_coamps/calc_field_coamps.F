      integer(ive_ptr) function calc_field(name,ndims,dims,stag,min,max,
     &             missing,data_units,data_display_units,dim_names)

      use default_module
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
     &                 interp_pres,geopot,x_grid,y_grid,
     &                 geowind,urot2,windspeed,reldiff,time_diff,
     &                 variance, maxmin_reg, time_mean, the2temp,
     &                 interp_height, kgain, localize
c
      call setstr(fld,name)
c
      data_units(1:len(data_units)) = ' '
      data_display_units(1:len(data_display_units)) = ' '
      iend = index(fld,'[',.true.) - 1


      if(iend.lt.0) then
        fld_cut = trim(fld)
      else
        fld_cut = fld(1:iend)
      end if

      if(iend.ge.0) then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) return
      end if

      calc_field = 1
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
        calc_field = 1
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
         if (ccnt.ne.2) then
           write(*,*) 'Error:  expected two string arguments.'
           return
         else
           calc_field = deriv(carg(1),carg(2),1,ndims,dims,stag,min,
     &                   max,missing,data_units,data_display_units,
     &                   dim_names)
         return
         end if
c
c  Rotate u wind through angle theta
c
      case ('UROT')
        if (rcnt.ne.1) then
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
        if (ccnt.ne.2) then
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
        if (ccnt.ne.2.and.rcnt.ne.1) then
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
        if (ccnt.ne.1.and.icnt.ne.2) then
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
        if (ccnt.ne.1) then
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
        if (ccnt.ne.1) then
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
        if (rcnt.ne.4) then
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
        if (rcnt.ne.8) then
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
        if (rcnt.ne.8) then
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
        if (rcnt.ne.1) then
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
         if (ccnt.ne.2) then
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
         if (ccnt.ne.1) then
           write(*,*) 'Usage: variance[fldname]'
           return
         else
           calc_field = variance(carg(1),
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
         return
         end if
c
c  Interpolate a field to a height (z0<=0 --> abs(z0) agl )
c
      case ('INTERPZ')
         if (ccnt.ne.1.or.rcnt.ne.1) then
           write(*,*) 'Usage: interp[string field: real z (agl)]'
           return
         else
           calc_field = interp_height(carg(1),rarg(1),
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
         if (ccnt.ne.1.or.rcnt.ne.1) then
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
         if (rcnt.ne.1) then
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
         calc_field = geowind('U',
     &                  ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
         return
c
      case ('VG')
         calc_field = geopot('V',
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
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
        print *,ndims,missing
        do i = 1,4
          print *,i,dims(i),stag(i),min(i),max(i),dim_names(i)
        end do
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
     &                    rarg(1),rarg(2),rarg(3),lmap,fld_cut,
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
     &                    rarg(1),rarg(2),rarg(3),lmap,fld_cut,
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
        return

c
c  Calculates gaspari-coehn localization radius
c
      case ('GCLOC')
        call calc_field_message(fld_cut,rcnt,icnt,ccnt,3,0,1,err)
        if(err) return

        calc_field = localize(rarg(1),rarg(2),rarg(3),trim(carg(1)),
     &                    ndims,dims,stag,min,max,missing,data_units,
     &                    data_display_units,dim_names)
        do i = 1,4
          print *,i,dims(i),stag(i),min(i),max(i),dim_names(i)
        end do
        return

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
