#include "pointer.icl"
    
      ive_ptr function calc_field(name,ndims,dims,stag,min,max,
     >             missing,data_units,data_display_units,dim_names)
c     
c     calc_field : This routine is used to calculate user-derived field.
c     The return value of the routine is a pointer to the field values.
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

      integer ndims,dims(4) 
      real    stag(4),min(4),max(4),missing
      character*(*) name,data_units,data_display_units,dim_names(4)

      integer ibeg,iend,idvar,ircode,strbeg,strend,ierr,ncvid
      integer iarg(10),k,rcnt,icnt,ccnt
      real rarg(10),r_val
      character*80 fld,arg1,arg2,arg3,carg(10),str1
      character*1 comp
      logical isfunc,err,new

      common/t_cdfid/ idcdf
      integer idcdf

      ive_ptr deriv,tptr,getvar,def_cmp,mix_wpts,mix_tpts,rho,
     >        press_drag,drag_int,per_vol,pert,vol_int,run_int,
     >        area_int_z,area_int_x,area_int_y,bar,grid

      calc_field = 0

      call setstr(fld,name)

      data_units(1:len(data_units)) = ' '
      data_display_units(1:len(data_display_units)) = ' '
      iend = strend(fld)
*
*  List of available fields
*
      if (fld(1:iend).eq.'LIST') then
        calc_field = 0
        call list
      end if
*
*  Derivative
*
      if (isfunc(fld,'D',arg1,arg2,arg3)) then
        calc_field = deriv(arg1,arg2,ndims,dims,stag,min,max,
     >                 missing,data_units,data_display_units,
     >                 dim_names)
        return
      end if
*
*  Grid heights
*
      if (isfunc(fld,'GRID',arg1,arg2,arg3)) then
        calc_field = grid(arg1,ndims,dims,stag,min,max,missing,
     >                 data_units,data_display_units,dim_names)
        return
      end if
*
*  Buoyancy frequency
*
      if (fld(1:iend).eq.'N') then
        tptr = getvar('thet',ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names,new)
        calc_field = deriv('thet','z',ndims,dims,stag,min,max,
     >           missing,data_units,data_display_units,dim_names)
        call buoy_freq(%val(calc_field),%val(tptr),dims(1),dims(2),
     >           dims(3),dims(4),data_units,data_display_units,
     >           missing)
        if (new) then
          call freemem(tptr)
        end if
        return
      end if
*
*  Component of the deformation tensor
*
      if (isfunc(fld,'DEFCMP',arg1,arg2,arg3)) then
        calc_field = def_cmp(arg1,arg2,ndims,dims,stag,min,max,
     >           missing,data_units,data_display_units,dim_names)
        return
      end if
*
*  Square of deformation tensor
*
      if (fld(1:iend).eq.'DEF2_WPTS') then
        calc_field = mix_wpts('def2',ndims,dims,stag,min,max,
     >              missing,data_units,data_display_units,dim_names)
        return
      else if ((fld(1:iend).eq.'DEF2_TPTS').or.
     >             (fld(1:iend).eq.'DEF2')) then
        calc_field = mix_tpts('def2',ndims,dims,stag,min,max,
     >              missing,data_units,data_display_units,dim_names)
        return
      end if
*
*  Richardson number
*
      if (fld(1:iend).eq.'RI_INV_WPTS') then
        calc_field = mix_wpts('richinv',ndims,dims,stag,min,max,
     >            missing,data_units,data_display_units,dim_names)
        return
      else if ((fld(1:iend).eq.'RI_INV_TPTS').or.
     >             (fld(1:iend).eq.'RI_INV')) then
        calc_field = mix_tpts('richinv',ndims,dims,stag,min,max,
     >            missing,data_units,data_display_units,dim_names)
        return
      end if
*
*  Mixing coefficient
*
      if (fld(1:iend).eq.'KM_WPTS') then
        calc_field = mix_wpts('km',ndims,dims,stag,min,max,
     >             missing,data_units,data_display_units,dim_names)
        return
      else if ((fld(1:iend).eq.'KM_TPTS').or.
     >              (fld(1:iend).eq.'KM_IVE')) then
        calc_field = mix_tpts('km',ndims,dims,stag,min,max,
     >             missing,data_units,data_display_units,dim_names)
        return
      end if
*
*  Total density
*
      if (fld(1:iend).eq.'RHO') then
        calc_field = rho(.false.,ndims,dims,stag,min,max,missing,
     >                    data_units,data_display_units,dim_names)
        return
      else if (fld(1:iend).eq.'RHO0') then
        calc_field = rho(.true.,ndims,dims,stag,min,max,missing,
     >                    data_units,data_display_units,dim_names)
        return
      end if
*
*  Pressure drag, x and y components
*
      if (fld(1:4).eq.'DRAG') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.4) then
          write(*,*) 'Error:  expected four integer arguments, ',
     >                           'got ',icnt
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one character argument, ',
     >                           'got ',ccnt
          return
        end if
        calc_field = press_drag(carg(1),iarg(1),iarg(2),iarg(3),
     >                  iarg(4),ndims,dims,stag,min,max,missing,
     >                  data_units,data_display_units,dim_names)
        return
      end if 
*
*  Conversion of per unit mass quantity to per unit volume
*
      if (isfunc(fld,'PER_VOL',arg1,arg2,arg3)) then
        calc_field = per_vol(arg1,.false.,ndims,dims,stag,min,max,
     >           missing,data_units,data_display_units,dim_names)
        return
      else if (isfunc(fld,'PER_VOL_BASE',arg1,arg2,arg3)) then
        calc_field = per_vol(arg1,.true.,ndims,dims,stag,min,max,
     >           missing,data_units,data_display_units,dim_names)
        return
      end if
*
*  Perturbation from basic state
*
      if (fld(1:4).eq.'PERT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.2) then
          write(*,*) 'Error:  expected two string arguments, ',
     >                           'got ',ccnt
          return
        else
          calc_field = pert(carg(1),carg(2),ndims,dims,stag,min,max,
     >               missing,data_units,data_display_units,dim_names)
          return
        end if
      end if
*
*  Volume integral
*
      if (fld(1:7).eq.'VOL_INT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.6) then
          write(*,*) 'Error:  expected 6 integer arguments; ',
     >                   'got ',icnt
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected 1 string argument; ',
     >                    'got ',ccnt
          return
        else
          calc_field = vol_int(carg(1),iarg(1),iarg(2),iarg(3),
     >               iarg(4),iarg(5),iarg(6),ndims,dims,stag,
     >               min,max,missing,data_units,
     >               data_display_units,dim_names)
          return
        end if 
      end if
*
*  Running time integral
*
      if (fld(1:7).eq.'RUN_INT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((rcnt.eq.0).and.(icnt.eq.0)) then
          rarg(1) = 0.
          rarg(2) = 0.
        else if ((rcnt.eq.0).and.(icnt.eq.2)) then
          rarg(1) = real(iarg(1))
          rarg(2) = real(iarg(2))
        else if (rcnt.ne.2) then
          write(*,*) 'Error:  expected two real arguments.',
     >          ' (Please include decimal pts to indicate reals).'
          return
        end if
        calc_field = run_int(carg(1),rarg(1),rarg(2),ndims,dims,
     >                      stag,min,max,missing,data_units,
     >                      data_display_units,dim_names)
        return
      end if
*
*  Area integral at constant z
*
      if (fld(1:10).eq.'AREA_INT_Z') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected five integer arguments,',
     >                     ' got ',icnt
          return
        else
          calc_field = area_int_z(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >                     iarg(4),iarg(5),ndims,dims,stag,min,max,
     >                     missing,data_units,data_display_units,
     >                                     dim_names)
          return
        end if    
      end if
*
*  Area integral at constant x
*
      if (fld(1:10).eq.'AREA_INT_X') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected five integer arguments,',
     >                     ' got ',icnt
          return
        else
          calc_field = area_int_x(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >                     iarg(4),iarg(5),ndims,dims,stag,min,max,
     >                     missing,data_units,data_display_units,
     >                                     dim_names)
          return
        end if    
      end if
*
*  Area integral at constant y
*
      if (fld(1:10).eq.'AREA_INT_Y') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected five integer arguments,',
     >                     ' got ',icnt
          return
        else
          calc_field = area_int_y(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >                     iarg(4),iarg(5),ndims,dims,stag,min,max,
     >                     missing,data_units,data_display_units,
     >                                     dim_names)
          return
        end if
      end if
*
*  Area integral at constant z as function of height and time
*
      if (fld(1:3).eq.'BAR') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.4) then
          write(*,*) 'Error:  expected four integer arguments, ',
     >                         'got ',icnt
          return
        else
          calc_field = bar(carg,ccnt,iarg(1),iarg(2),iarg(3),iarg(4),
     >                   ndims,dims,stag,min,max,missing,data_units,
     >                   data_display_units,dim_names)
          return
        end if
      end if
        
      return
      end


c     WRITE(*,*) 'ndims = ',ndims
c     do k=1,ndims
c       WRITE(*,*) 'dims(',k,') = ',dims(k)
c       WRITE(*,*) 'stag(',k,') = ',stag(k)
c       WRITE(*,*) 'min(',k,') = ',min(k)
c       WRITE(*,*) 'max(',k,') = ',max(k)
c       WRITE(*,*) 'dim_names(',k,') = ',
c    >    dim_names(k)(strbeg(dim_names(k)):strend(dim_names(k)))
c     enddo
c     WRITE(*,*) 'missing = ',missing
c     WRITE(*,*) 'data_units = ',
c    >        data_units(strbeg(data_units):strend(data_units))
c     WRITE(*,*) 'data_display_units = ',
c    >   data_display_units(strbeg(data_display_units):
c    >                  strend(data_display_units)) 



