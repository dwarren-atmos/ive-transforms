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
      character(len=*) :: name
      character(len=80):: data_units,data_display_units
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
     &                 vort_y,pvort,sigdiff
c
      call setstr(fld,name)
c
      data_units(1:len(data_units)) = ' '
      data_display_units(1:len(data_display_units)) = ' '
!      data_units = ' '
!
      data_display_units = ' '
      iend = index(fld,'[',BACK=.false.) - 1
c .true = rindex .false = index

      fld_cut=' '
      if(iend.lt.0) then
        fld_cut = trim(fld(1:len(trim(fld))-1))
      else
        fld_cut = fld(1:iend)
      end if
      carg=' '
      if(iend.ge.0) then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) return
      end if

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
c
c  Grid heights
c
!      case ('GRID')
!       call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,1,err)
!       if(err) return

!       calc_field = vert_grid(carg(1),ndims,dims,stag,min,max,missing,
!     &                 data_units,data_display_units,dim_names)
c
c  Translate domain
c
!      case ('TRANSLATE')
!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,2,0,0,err)
!        if(err) return

!        call translate_domain(rarg(1),rarg(2))
!        !calc_field = 1
!        return
c
c  Translate time
c
!      case ('TRANS_TIME')
!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,0,0,err)
!        if(err) return
!        call translate_time(rarg(1))

        !calc_field = 1
!        return
c
c  Turn vertical coordinate transformation on or off
c
!      case ('TRANS_ON')
!        call trans_switch(1)
!        calc_field = 1
!        return

!      case ('TRANS_OFF')
!        call trans_switch(0)
!        calc_field = 1
!        return
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
!      case ('UROT')
!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,0,2,err)
!        if(err) return

!        calc_field = urot2(rarg(1),carg(1),carg(2),
!     &                ndims,dims,stag,min,max,missing,
!     &                data_units,data_display_units,dim_names)
!        return
c
c  OMEGA
c
      case ('OMEGA')
!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,3,err)
!        if (err) return

!        calc_field = omega(carg(1),carg(2),carg(3),
!     &                      ndims,dims,stag,min,max,missing,
!     &                      data_units,data_display_units,dim_names)
        calc_field = omega('Pressure','T','W',
     &                      ndims,dims,stag,min,max,missing,
     &                      data_units,data_display_units,dim_names)
        return
c    
c  Buoyancy frequency
c    
!      case ('NSPEC')
!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,1,err)
!        if (err) return

!        calc_field = buoy_freq(carg(1),
!     &                   ndims,dims,stag,min,max,missing,data_units,
!     &                   data_display_units,dim_names,0)
!        return

!      case ('N')
!        calc_field = buoy_freq('the',ndims,dims,stag,min,max,missing,
!     &                 data_units,data_display_units,dim_names,0)
!        return

!      case ('N2')
!        calc_field = buoy_freq('the',ndims,dims,stag,min,max,missing,
!     &                 data_units,data_display_units,dim_names,1)
!        return
c
!      case ('THE2TEMP')
!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
!        if(err) return

!        calc_field = the2temp(carg(1),carg(2),
!     &                  ndims,dims,stag,min,max,missing,data_units,
!     &                  data_display_units,dim_names)
!        return

!      case ('VORT_X') ! x component of vorticity
!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
!        if(err) return

!        calc_field = vort_x(carg(1),carg(2),
!     &                  ndims,dims,stag,min,max,missing,data_units,
!     &                  data_display_units,dim_names)
!        return

!      case ('VORT_Y') ! y component of vorticity
!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
!        if(err) return

!        calc_field = vort_y(carg(1),carg(2),
!     &                  ndims,dims,stag,min,max,missing,data_units,
!     &                  data_display_units,dim_names)
!        return

!      case ('VORT_Z') ! z component of vorticity
!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,0,0,2,err)
!        if(err) return
!
!        calc_field = vort_z(carg(1),carg(2),
!     &                  ndims,dims,stag,min,max,missing,data_units,
!     &                  data_display_units,dim_names)
!        return

      case ('PVORT') ! potential vorticity
!                    In WRF we know the names:U,V,T,F etc
!                    So they are hard coded unlike coamps etc.
        calc_field = pvort(ndims,dims,stag,min,max,missing,data_units,
     &                  data_display_units,dim_names)
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
