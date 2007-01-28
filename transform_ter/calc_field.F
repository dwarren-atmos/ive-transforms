#include "pointer.icl"
    
      ive_ptr function calc_field(name,ndims,dims,stag,min,max,
     >             missing,data_units,data_display_units,dim_names)
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
      integer ndims,dims(4) 
      real    stag(4),min(4),max(4),missing
      character*(*) name,data_units,data_display_units,dim_names(4)

      integer ibeg,iend,idvar,ircode,strbeg,strend,ierr,ncvid
      integer iarg(10),k,rcnt,icnt,ccnt
      real rarg(10),r_val
      character*80 fld,arg1,arg2,arg3,carg(10),str1,tstfld,mkfield
      character*1 comp
      logical isfunc,err,new

      common/t_cdfid/ idcdf
      integer idcdf

      common/t_out/id_out,x_dims,y_dims,z_dims,t_dims,numout,flnam
      integer id_out(10),x_dims(3,10,2),y_dims(3,10,2),z_dims(3,10,2),
     >           t_dims(1,10,2),numout(10)
      character*(80) flnam(10)

      ive_ptr deriv,tptr,getvar,def2,rho,press_drag,press_drag_spec,
     >        drag_int,per_vol,pert,vol_int,run_int,area_int_z,
     >        area_int_x,area_int_y,bar,grid,ref,abs_ive 
      ive_ptr diff_area_int,z_ave,ru_pr,baseline,div,
     >        vol_ave,area_ave_x,area_ave_y,area_ave_z,zf_top,
     >        drg_bot,h_ave,volume,uabs,uflx,quot,prod
      ive_ptr max_min,ufcv,symm_swtch,area_int_z_vp,
     >        diff_int_vp,vol_int_vp,vol_ave_vp,area_int_x_vp
      ive_ptr area_int_y_vp,drag_vp,drag_ave_vp,scalar,area_ave_z_vp,
     >        z_coord,dist,net_disp,pvs,tangent,vortx,vorty,vortz,
     >        t_ave,buoy_freq,spctrm,wndspeed,area_int_xz0_y,nrmlz
      ive_ptr area_int_yz0_x,vol_int_xyz0,maketime,ave,drag_density,
     >        area_int_xz0,area_int_yz0,flux_density_yz0_x,pbase,
     >        r0b,rb,bernoulli,q_ive,r_ive,s_ive,surf,upert,v2z
      ive_ptr massflx,massflx_vp,massflx_mat,upstrm_hgt,area,
     >        thres_vol
      calc_field = 0

      call setstr(fld,name)

      data_units(1:len(data_units)) = ' '
      data_display_units(1:len(data_display_units)) = ' '
      iend = strend(fld)

c     WRITE(*,*) 'In calc_field, fld is ',fld(1:iend)

  
*
*  List of available fields
*
      if (fld(1:iend).eq.'LIST') then
        calc_field = 0
        call list
*
*  Derivative (fourth order diff)
*
      else if (fld(1:2).eq.'D[') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.2) then
          write(*,*) 'Error:  expected two string arguments.'
          return
        else
          calc_field = deriv(carg(1),carg(2),1,ndims,dims,stag,min,
     >                   max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Second derivative (curvature)
*
      else if (fld(1:2).eq.'D2') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.2) then
          write(*,*) 'Error:  expected two string arguments.'
          return
        else
          calc_field = deriv(carg(1),carg(2),2,ndims,dims,stag,min,
     >                   max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Derivative of product at thermodynamic pts
*
c      else if (fld(1:3).eq.'D''[') then
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c          return
c        else if (ccnt.lt.2) then
c          write(*,*) 'Error:  expected at least two string arguments.'
c          return
c        else
c          calc_field = deriv_pr(carg(1),carg(ccnt),ccnt-1,ndims,dims,
c     >                  stag,min,max,missing,data_units,
c     >                  data_display_units,dim_names)
c        end if
*
*  Use a given field as vertical coordinate
*
      else if (fld(1:7).eq.'Z_COORD') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.lt.2) then
          write(*,*) 'Error:  expected at least two string arguments.'
          return
        else
          if (icnt.eq.0) then
            iarg(1) = 0
          end if
          calc_field = z_coord(carg(1),carg(2),ccnt-1,iarg(1),
     >                  ndims,dims,stag,min,max,missing,data_units,
     >                  data_display_units,dim_names)
        end if
*
*  Tangent of vector to level surface of given field
*
      else if (fld(1:7).eq.'TANGENT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.5) then
          write(*,*) 'Error:  expected five string arguments.'
          return
        else 
          calc_field = tangent(carg(1),carg(2),carg(3),carg(4),
     >               carg(5),ndims,dims,stag,min,max,missing,
     >               data_units,data_display_units,dim_names)
        end if
*
*  Power spectrum
*
      else if (fld(1:6).eq.'SPCTRM') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
        else
          calc_field = spctrm(carg(1),ndims,dims,stag,min,max,
     >                   missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Grid heights
*
      else if (isfunc(fld,'GRID',arg1,arg2,arg3)) then
        calc_field = grid(arg1,ndims,dims,stag,min,max,missing,
     >                 data_units,data_display_units,dim_names)
*
*  Buoyancy frequency
*
      else if (fld(1:iend).eq.'N') then
        calc_field = buoy_freq(ndims,dims,stag,min,max,missing,
     >                 data_units,data_display_units,dim_names,0)
      else if (fld(1:iend).eq.'N2') then
        calc_field = buoy_freq(ndims,dims,stag,min,max,missing,
     >                 data_units,data_display_units,dim_names,1)
*
*  Deformation (squared)
*
      else if (fld(1:iend).eq.'DEF2') then
        calc_field = def2(ndims,dims,stag,min,max,missing,
     >                data_units,data_display_units,dim_names)
*
*  Wind speed
*
      else if (fld(1:iend).eq.'SPEED') then
        calc_field = wndspeed(0,ndims,dims,stag,min,max,missing,
     >                data_units,data_display_units,dim_names)
      else if (fld(1:iend).eq.'SPEED_PERT') then
        calc_field = wndspeed(1,ndims,dims,stag,min,max,missing,
     >                data_units,data_display_units,dim_names)
*
*  Total density
*
      else if (fld(1:iend).eq.'RHO') then
        calc_field = rho(.false.,ndims,dims,stag,min,max,missing,
     >                 data_units,data_display_units,dim_names,.false.)
      else if (fld(1:iend).eq.'RHO0') then
        calc_field = rho(.true.,ndims,dims,stag,min,max,missing,
     >                 data_units,data_display_units,dim_names,.false.)
      else if (fld(1:iend).eq.'RHOB') then
        calc_field = rho(.false.,ndims,dims,stag,min,max,missing,
     >                 data_units,data_display_units,dim_names,.true.)
      else if (fld(1:iend).eq.'RHOB0') then
        calc_field = rho(.true.,ndims,dims,stag,min,max,missing,
     >                 data_units,data_display_units,dim_names,.true.)

*
*  Pressure drag, vertical profile
*
      else if (fld(1:7).eq.'DRAG_VP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected five integer arguments, ',
     >                           'got ',icnt
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one character argument, ',
     >                           'got ',ccnt
          return
        end if
        calc_field = drag_vp(carg(1),iarg(1),iarg(2),iarg(3),
     >        iarg(4),iarg(5),ndims,dims,stag,min,max,missing,
     >        data_units,data_display_units,dim_names)
*
*  Average pressure drag, vertical profile
*
      else if (fld(1:11).eq.'DRAG_AVE_VP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected five integer arguments, ',
     >                           'got ',icnt
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one character argument, ',
     >                           'got ',ccnt
          return
        end if
        calc_field = drag_ave_vp(carg(1),iarg(1),iarg(2),iarg(3),
     >        iarg(4),iarg(5),ndims,dims,stag,min,max,missing,
     >        data_units,data_display_units,dim_names)
*
*  Drag/Length as function of y and t
*
      else if (fld(1:12).eq.'DRAG_DENSITY') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.2) then
          write(*,*) 'Error:  expected two integer arguments.'
          return
        end if
        calc_field = drag_density(iarg(1),iarg(2),ndims,dims,stag,
     >                min,max,missing,data_units,data_display_units,
     >                dim_names)
*
*  Pressure drag (special)
*
      else if (fld(1:9).eq.'DRAG_SPEC') then
        calc_field = press_drag_spec(ndims,dims,stag,min,max,
     >                   missing,data_units,data_display_units,
     >                   dim_names)
*
*  Pressure drag 
*
      else if (fld(1:4).eq.'DRAG') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.6) then
          write(*,*) 'Error:  expected six integer arguments, ',
     >                           'got ',icnt
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one character argument, ',
     >                           'got ',ccnt
          return
        end if
        calc_field = press_drag(carg(1),iarg(1),iarg(2),iarg(3),
     >        iarg(4),iarg(5),iarg(6),ndims,dims,stag,min,max,missing,
     >        data_units,data_display_units,dim_names)
*
*  Conversion of per unit mass quantity to per unit volume
*
      else if (isfunc(fld,'PER_VOL',arg1,arg2,arg3)) then
        calc_field = per_vol(arg1,0,'total',ndims,dims,stag,min,max,
     >           missing,data_units,data_display_units,dim_names)
      else if (isfunc(fld,'PER_VOL_BASE',arg1,arg2,arg3)) then
        calc_field = per_vol(arg1,0,'base',ndims,dims,stag,min,max,
     >           missing,data_units,data_display_units,dim_names)
      else if (fld(1:11).eq.'PER_VOL_REF') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((ccnt.ne.1).or.(icnt.ne.1)) then
          write(*,*) 'Error:  expected form ',
     >                  'PER_VOL_REF[string:integer].'
          return
        else
          calc_field = per_vol(carg(1),iarg(1),'ref',ndims,dims,
     >                  stag,min,max,missing,data_units,
     >                  data_display_units,dim_names)
        end if
*
*  Perturbation from basic state
*
      else if (fld(1:4).eq.'PERT') then
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

        end if
*
*  Perturbation u (including startup)
*
      else if (fld(1:5).eq.'UPERT') then
        calc_field = upert(ndims,dims,stag,min,max,missing,
     >                      data_units,data_display_units,
     >                      dim_names)
*
*  Field at a fixed time
*
      else if (fld(1:3).eq.'REF') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((ccnt.ne.1).or.(icnt.ne.1)) then
          write(*,*) 'Error:  expected the form REF[string:integer].'
          return
        else
          calc_field = ref(carg(1),iarg(1),ndims,dims,stag,min,max,
     >              missing,data_units,data_display_units,dim_names)
        end if
*
*  Field at untransformed heights
*
      else if (fld(1:3).eq.'V2Z') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.1) then
          write(6,*) 'Error:  expected the form V2Z[string].'
          return
        else
          calc_field = v2z(carg(1),ndims,dims,stag,min,max,missing,
     >                  data_units,data_display_units,dim_names)
        end if
*
*  Absolute value
*
      else if (isfunc(fld,'ABS',arg1,arg2,arg3)) then
        calc_field = abs_ive(arg1,ndims,dims,stag,min,max,missing,
     >            data_units,data_display_units,dim_names)
*
*  Minimum value of field
*
      else if (fld(1:7).eq.'MIN_VAL') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
        else
          calc_field = max_min(carg(1),0,0,0,1,1,1,1,ndims,dims,stag,
     >                   min,max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Maximum value of field
*
      else if (fld(1:7).eq.'MAX_VAL') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
        else
          calc_field = max_min(carg(1),1,0,0,1,1,1,1,ndims,dims,stag,
     >                   min,max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Minimum value of field at surface
*
      else if (fld(1:8).eq.'MIN_SURF') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
        else
          calc_field = max_min(carg(1),0,2,0,1,1,1,1,ndims,dims,stag,
     >                   min,max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Maximum value of field at surface
*
      else if (fld(1:8).eq.'MAX_SURF') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
        else
          calc_field = max_min(carg(1),1,2,0,1,1,1,1,ndims,dims,stag,
     >                   min,max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Minimum value of field in a subdomain
*
      else if (fld(1:7).eq.'MIN_SUB') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
      else if (icnt.ne.4) then
      write(*,*) 'Error:  expected four integer arguments.',
     >                'got ', icnt
      return
        else
          calc_field = max_min(carg(1),0,0,1,
     >                  iarg(1),iarg(2),iarg(3),iarg(4),ndims,dims,stag,
     >                  min,max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Maximum value of field in a subdomain
*
      else if (fld(1:7).eq.'MAX_SUB') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
      else if (icnt.ne.4) then
      write(*,*) 'Error:  expected four integer arguments.',
     >                'got ', icnt
      return
        else
          calc_field = max_min(carg(1),1,0,1,
     >      iarg(1),iarg(2),iarg(3),iarg(4),ndims,dims,stag,
     >                   min,max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Minimum value of field at surface in a subdomain
*
      else if (fld(1:11).eq.'MINSURF_SUB') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
	else if (icnt.ne.4) then
	  write(*,*) 'Error:  expected four integer arguments.',
     >                'got ', icnt
	  return
        else
          calc_field = max_min(carg(1),0,1,2,
     >	               iarg(1),iarg(2),iarg(3),iarg(4),ndims,dims,stag,
     >                   min,max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Maximum value of field at surface in a subdomain
*
      else if (fld(1:11).eq.'MAXSURF_SUB') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
	else if (icnt.ne.4) then
	  write(*,*) 'Error:  expected four integer arguments.',
     >                'got ', icnt
	  return
        else
          calc_field = max_min(carg(1),1,1,1,
     >	               iarg(1),iarg(2),iarg(3),iarg(4),ndims,dims,stag,
     >                   min,max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Minimum value of field at lowest grid point in a subdomain
*
      else if (fld(1:11).eq.'MINLGRD_SUB') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
	else if (icnt.ne.4) then
	  write(*,*) 'Error:  expected four integer arguments.',
     >                'got ', icnt
	  return
        else
          calc_field = max_min(carg(1),0,2,1,
     >	               iarg(1),iarg(2),iarg(3),iarg(4),ndims,dims,stag,
     >                   min,max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Maximum value of field at lowest grid point in a subdomain
*
      else if (fld(1:11).eq.'MAXLGRD_SUB') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then 
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
	else if (icnt.ne.4) then
	  write(*,*) 'Error:  expected four integer arguments.',
     >                'got ', icnt
	  return
        else
          calc_field = max_min(carg(1),1,2,1,
     >	               iarg(1),iarg(2),iarg(3),iarg(4),ndims,dims,stag,
     >                   min,max,missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Normalization (by maximum value)
*
      else if (fld(1:5).eq.'NRMLZ') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((ccnt.ne.1).or.(rcnt.ne.1)) then
          write(6,*) 'Error:  expected one string and one real ',
     >                           'argument.'
          return
        else
          calc_field = nrmlz(carg(1),rarg(1),ndims,dims,stag,min,max,
     >                   missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Volume integral, vertical profile
*
      else if (fld(1:10).eq.'VOL_INT_VP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected five integer arguments; ',
     >                   'got ',icnt
          return
        else if (ccnt.lt.1) then
          write(*,*) 'Error:  must specify at least one field.'
          return
        else
          calc_field = vol_int_vp(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >               iarg(4),iarg(5),ndims,dims,stag,min,max,missing,
     >               data_units,data_display_units,dim_names)
        end if 
*
*  Volume average, vertical profile
* 
      else if (fld(1:10).eq.'VOL_AVE_VP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected five integer arguments; ',
     >                   'got ',icnt
          return
        else if (ccnt.lt.1) then
          write(*,*) 'Error:  must specify at least one field.'
          return
        else
          calc_field = vol_ave_vp(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >               iarg(4),iarg(5),ndims,dims,stag,min,max,missing,
     >               data_units,data_display_units,dim_names)
        end if 
*
*  Volume integral over x,y and z0 (i.e., between material surfaces)
*
      else if (fld(1:12).eq.'VOL_INT_XYZ0') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((icnt.ne.4).or.(rcnt.ne.2)) then
          write(6,*) 'Error:  expected 4 integer arguments and ',
     >                     '2 real arguments.'
          return
        else if (ccnt.lt.1) then
          write(6,*) 'Error:  must specify at least one field.'
          return
        else
          calc_field = vol_int_xyz0(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >               iarg(4),rarg(1),rarg(2),ndims,dims,stag,min,max,
     >               missing,data_units,data_display_units,dim_names)
        end if 
*
*  Volume integral
*
      else if (fld(1:7).eq.'VOL_INT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.6) then
          write(*,*) 'Error:  expected 6 integer arguments; ',
     >                   'got ',icnt
          return
        else if (ccnt.lt.1) then
          write(*,*) 'Error:  must specify at least one field.'
          return
        else
          calc_field = vol_int(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >               iarg(4),iarg(5),iarg(6),ndims,dims,stag,
     >               min,max,missing,data_units,
     >               data_display_units,dim_names)
        end if 
*
*  Volume average
*
      else if (fld(1:7).eq.'VOL_AVE') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.6) then
          write(*,*) 'Error:  expected 6 integer arguments; ',
     >                   'got ',icnt
          return
        else if (ccnt.lt.1) then
          write(*,*) 'Error:  must specify at least one field.'
          return
        else
          calc_field = vol_ave(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >               iarg(4),iarg(5),iarg(6),ndims,dims,stag,
     >               min,max,missing,data_units,
     >               data_display_units,dim_names)
        end if 
*
*  Volume of field in subdomain the is within some threshold
*
	else if (fld(1:9).eq.'THRES_VOL') then
	 call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
	 if (err) then
	   return
	 else if (icnt.ne.4) then
         write(*,*) 'Error:  expected 4 integer arguments; ',
     >                         'got ',icnt
         return
	 else if (rcnt.ne.2) then
         write(*,*) 'Error:  expected 2 real arguments; ',
     >                         'got ',rcnt
         return
	 else if (ccnt.ne.1) then
         write(*,*) 'Error:  expected 1 string arguments; ',
     >                         'got ',ccnt
         return
	 else
	   calc_field=thres_vol(carg(1),rarg(1),rarg(2),
     >			 iarg(1),iarg(2),iarg(3),iarg(4),
     >                   ndims,dims,stag,min,max,missing,data_units,
     >                   data_display_units,dim_names)
	 end if
*
*  Area at lowest grid point of given field that passes threshold
*
	else if (fld(1:9).eq.'GRND_AREA') then
	 call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
	 if (err) then
	   return
	 else if (icnt.ne.4) then
         write(*,*) 'Error:  expected 4 integer arguments; ',
     >                         'got ',icnt
         return
	 else if (rcnt.ne.1) then
         write(*,*) 'Error:  expected 4 real arguments; ',
     >                         'got ',rcnt
         return
	 else if (ccnt.ne.1) then
         write(*,*) 'Error:  expected 1 string arguments; ',
     >                         'got ',ccnt
         return
	 else
	   calc_field=area(carg,iarg(1),iarg(2),iarg(3),iarg(4),1,rarg,
     >                     ndims,dims,stag,min,max,missing,data_units,
     >                     data_display_units,dim_names)
	 end if
*
*  Area average at constant x
*
      else if (fld(1:10).eq.'AREA_AVE_X') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected 5 integer arguments; ',
     >                          'got ',icnt
          return
        else if (ccnt.lt.1) then
          write(*,*) 'Error:  must specify at least one field.'
          return
        else
          calc_field = area_ave_x(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >            iarg(4),iarg(5),ndims,dims,stag,min,max,missing,
     >            data_units,data_display_units,dim_names)
        end if
*
*  Area average at constant y
*
      else if (fld(1:10).eq.'AREA_AVE_Y') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected 5 integer arguments; ',
     >                          'got ',icnt
          return
        else if (ccnt.lt.1) then
          write(*,*) 'Error:  must specify at least one field.'
          return
        else
          calc_field = area_ave_y(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >            iarg(4),iarg(5),ndims,dims,stag,min,max,missing,
     >            data_units,data_display_units,dim_names)
        end if
*
*  Area average at constant z, vertical profile
*
      else if (fld(1:13).eq.'AREA_AVE_Z_VP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((icnt.ne.4).and.(icnt.ne.5)) then
          write(*,*) 'Error:  expected four integer arguments,',
     >                          ' got ',icnt
          return
        else
          if (icnt.eq.4) iarg(5) = 0
          calc_field = area_ave_z_vp(carg,ccnt,iarg(1),iarg(2),
     >                  iarg(3),iarg(4),iarg(5),ndims,dims,stag,min,
     >                  max,missing,data_units,data_display_units,
     >                                  dim_names)
        end if
*
*  Area average at constant z
*
      else if (fld(1:10).eq.'AREA_AVE_Z') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected 5 integer arguments; ',
     >                          'got ',icnt
          return
        else if (ccnt.lt.1) then
          write(*,*) 'Error:  must specify at least one field.'
          return
        else
          calc_field = area_ave_z(carg,ccnt,iarg(1),iarg(2),iarg(3),
     >            iarg(4),iarg(5),ndims,dims,stag,min,max,missing,
     >            data_units,data_display_units,dim_names)
        end if
*
*  Running time integral
*
      else if (fld(1:7).eq.'RUN_INT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          write(*,*) 'Error:  could not get arguments.'
          return
        else if ((icnt.eq.0).and.(ccnt.eq.1)) then
          iarg(1) = 0
        else if ((icnt.ne.1).or.(ccnt.ne.1)) then
          write(*,*) 'Error:  expected one string and one ',
     >                    'integer argument.'
          return
        end if
        calc_field = run_int(carg(1),iarg(1),ndims,dims,
     >                      stag,min,max,missing,data_units,
     >                      data_display_units,dim_names)
*
*  Area integral at constant z, vertical profile
*
      else if (fld(1:13).eq.'AREA_INT_Z_VP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((icnt.ne.4).and.(icnt.ne.5)) then
          write(*,*) 'Error:  expected four integer arguments,',
     >                          ' got ',icnt
          return
        else
          if (icnt.eq.4) iarg(5) = 0
          calc_field = area_int_z_vp(carg,ccnt,iarg(1),iarg(2),
     >                  iarg(3),iarg(4),iarg(5),ndims,dims,stag,min,
     >                  max,missing,data_units,data_display_units,
     >                                  dim_names)
        end if
*
*  Area integral at constant z
*
      else if (fld(1:10).eq.'AREA_INT_Z') then
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
        end if    
*
*  Ratio of Mass out fixed height top of control volume to 
*  mass in upstream side
*
      else if (fld(1:11).eq.'MASSFLX_TOP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.6) then
          write(*,*) 'Error:  expected six integer arguments,',
     >                     ' got ',icnt
          return
        else
          calc_field = massflx(iarg(1),iarg(2),iarg(3),
     >                     iarg(4),iarg(5),iarg(6),.false.,ndims,
     >                     dims,stag,min,max,missing,data_units,
     >                     data_display_units,dim_names)
        end if    
*
*  Ratio of Mass out lateral sides of control volume to 
*  mass in upstream side between ridged surfaces
*
      else if (fld(1:11).eq.'MASSFLX_LAT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.6) then
          write(*,*) 'Error:  expected six integer arguments,',
     >                     ' got ',icnt
          return
        else
          calc_field = massflx(iarg(1),iarg(2),iarg(3),
     >                     iarg(4),iarg(5),iarg(6),.true.,ndims,
     >                     dims,stag,min,max,missing,data_units,
     >                     data_display_units,dim_names)
        end if    
*
*  Ratio of Mass out lateral sides of control volume to 
*  mass in upstream side between to material surfaces
*
      else if (fld(1:11).eq.'MASSFLX_MAT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.4) then
          write(*,*) 'Error:  expected four integer arguments,',
     >                     ' got ',icnt
          return
        else if (rcnt.ne.2) then
          write(*,*) 'Error:  expected two real arguments,',
     >                     ' got ',rcnt
          return
        else
          calc_field = massflx_mat(iarg(1),iarg(2),iarg(3),
     >                     iarg(4),rarg(1),rarg(2),ndims,
     >                     dims,stag,min,max,missing,data_units,
     >                     data_display_units,dim_names)
        end if    
*
*  VP of Ratio of Mass out lateral sides of control volume to 
*  mass in upstream side
*
      else if (fld(1:10).eq.'MASSFLX_VP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.7) then
          write(*,*) 'Error:  expected seven integer arguments,',
     >                     ' got ',icnt
          return
        else
          calc_field = massflx_vp(iarg(1),iarg(2),iarg(3),
     >                     iarg(4),iarg(5),iarg(6),iarg(7),ndims,
     >                     dims,stag,min,max,missing,data_units,
     >                     data_display_units,dim_names)
        end if    
*
*  Area integral at constant x, vertical profile
*
      else if (fld(1:13).eq.'AREA_INT_X_VP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.4) then
          write(*,*) 'Error:  expected four integer arguments,',
     >                     ' got ',icnt
          return
        else if (ccnt.lt.1) then
          write(*,*) 'Error:  need to specify at least one field.'
          return
        else
          calc_field = area_int_x_vp(carg,ccnt,iarg(1),iarg(2),
     >            iarg(3),iarg(4),ndims,dims,stag,min,max,missing,
     >            data_units,data_display_units,dim_names)
        end if    
       
*
*  Area integral over y and z0 (between material surfaces) at
*    constant x
*
      else if (fld(1:14).eq.'AREA_INT_YZ0_X') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((rcnt.ne.2).or.(icnt.ne.3)) then
          write(6,*) 'Error:  expected three integer arguments ',
     >                         'and two real arguments.'
          return
        else if (ccnt.lt.1) then
          write(6,*) 'Error:  must specify field(s) to integrate.'
          return
        else
          calc_field = area_int_yz0_x(carg,ccnt,iarg(1),iarg(2),
     >                   rarg(1),rarg(2),iarg(3),ndims,dims,stag,
     >                   min,max,missing,data_units,
     >                   data_display_units,dim_names)
        end if    
*
*  Area integral over x and z0 (between material surfaces) at
*    constant y
*
      else if (fld(1:14).eq.'AREA_INT_XZ0_Y') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((rcnt.ne.2).or.(icnt.ne.3)) then
          write(6,*) 'Error:  expected three integer arguments ',
     >                         'and two real arguments.'
          return
        else if (ccnt.lt.1) then
          write(6,*) 'Error:  must specify field(s) to integrate.'
          return
        else
          calc_field = area_int_xz0_y(carg,ccnt,iarg(1),iarg(2),
     >                   rarg(1),rarg(2),iarg(3),ndims,dims,stag,
     >                   min,max,missing,data_units,
     >                   data_display_units,dim_names)
        end if    
*
*  Area integral over x and z0 (between material surfaces) as
*    function of y
*
      else if (fld(1:12).eq.'AREA_INT_XZ0') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((rcnt.ne.2).or.(icnt.ne.2)) then
          write(6,*) 'Error:  expected two integer arguments ',
     >                         'and two real arguments.'
          return
        else if (ccnt.lt.1) then
          write(6,*) 'Error:  must specify field(s) to integrate.'
          return
        else
          calc_field = area_int_xz0(carg,ccnt,iarg(1),iarg(2),
     >                   rarg(1),rarg(2),ndims,dims,stag,min,max,
     >                   missing,data_units,data_display_units,
     >                   dim_names)
        end if    
*
*  Area integral over y and z0 (between material surfaces) as
*    function of x
*
      else if (fld(1:12).eq.'AREA_INT_YZ0') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((rcnt.ne.2).or.(icnt.ne.2)) then
          write(6,*) 'Error:  expected two integer arguments ',
     >                         'and two real arguments.'
          return
        else if (ccnt.lt.1) then
          write(6,*) 'Error:  must specify field(s) to integrate.'
          return
        else
          calc_field = area_int_yz0(carg,ccnt,iarg(1),iarg(2),
     >                   rarg(1),rarg(2),ndims,dims,stag,min,max,
     >                   missing,data_units,data_display_units,
     >                   dim_names)
        end if    
*
*  Flux/Length over y and z0 (between material surfaces) as function
*    of y at constant x
*
      else if (fld(1:18).eq.'FLUX_DENSITY_YZ0_X') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((rcnt.ne.2).or.(icnt.ne.1)) then
          write(6,*) 'Error:  expected one integer arguments ',
     >                         'and two real arguments.'
          return
        else if (ccnt.lt.1) then
          write(6,*) 'Error:  must specify field(s) to integrate.'
          return
        else
          calc_field = flux_density_yz0_x(carg,ccnt,rarg(1),rarg(2),
     >                   iarg(1),ndims,dims,stag,min,max,missing,
     >                   data_units,data_display_units,dim_names)
        end if    
*
*  Area integral at constant x
*
      else if (fld(1:10).eq.'AREA_INT_X') then
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
        end if    
*
*  Area integral at constant y, vertical profile
*
      else if (fld(1:13).eq.'AREA_INT_Y_VP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.4) then
          write(*,*) 'Error:  expected four integer arguments,',
     >                     ' got ',icnt
          return
        else if (ccnt.lt.1) then
          write(*,*) 'Error:  need to specify at least one field.'
          return
        else
          calc_field = area_int_y_vp(carg,ccnt,iarg(1),iarg(2),
     >            iarg(3),iarg(4),ndims,dims,stag,min,max,missing,
     >            data_units,data_display_units,dim_names)
        end if
*
*  Area integral at constant y
*
      else if (fld(1:10).eq.'AREA_INT_Y') then
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
        end if
*
*  Horizontal average over fixed domain
*
c      else if (fld(1:6).eq.'H_MEAN') then
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c          return
c        else if (icnt.ne.5) then
c          write(*,*) 'Error:  expected five integer arguments.'
c          return
c        else if (ccnt.lt.1) then
c          write(*,*) 'Error:  expected at least one string argument.'
c          return
c        else
c          calc_field = h_mean('N',carg,ccnt,iarg(1),iarg(2),iarg(3),
c     >           iarg(4),iarg(5),ndims,dims,stag,min,max,missing,
c     >           data_units,data_display_units,dim_names)
c        end if
*
*  Vertical derivative of horizontal average over fixed domain
*
c      else if (fld(1:10).eq.'DH_MEAN_DZ') then
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c          return
c        else if (icnt.ne.5) then
c          write(*,*) 'Error:  expected five integer arguments.'
c          return
c        else if (ccnt.lt.1) then
c          write(*,*) 'Error:  expected at least one string argument.'
c          return
c        else 
c          calc_field = h_mean('Z',carg,ccnt,iarg(1),iarg(2),iarg(3),
c     >           iarg(4),iarg(5),ndims,dims,stag,min,max,missing,
c     >           data_units,data_display_units,dim_names)
c        end if
*
*  Line integral at constant x and z
*
c      else if (fld(1:10).eq.'LINE_INT_X') then
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c          return
c        else if (icnt.ne.4) then
c          write(*,*) 'Error:  expected four integer arguments.'
c          return
c        else if (ccnt.lt.1) then
c          write(*,*) 'Error:  expected at least one string argument.'
c          return
c        else
c          calc_field = line_int_x(carg,ccnt,iarg(1),iarg(2),iarg(3),
c     >            iarg(4),ndims,dims,stag,min,max,missing,data_units,
c     >            data_display_units,dim_names)
c        end if
*
*  Line integral at constant y and z
*
c      else if (fld(1:10).eq.'LINE_INT_Y') then
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c          return
c        else if (icnt.ne.4) then
c          write(*,*) 'Error:  expected four integer arguments.'
c          return
c        else if (ccnt.lt.1) then
c          write(*,*) 'Error:  expected at least one string argument.'
c          return
c        else
c          calc_field = line_int_y(carg,ccnt,iarg(1),iarg(2),iarg(3),
c     >            iarg(4),ndims,dims,stag,min,max,missing,data_units,
c     >            data_display_units,dim_names)
c        end if
*
*  Difference between line integrals (for area-average budgets)
*
c      else if (fld(1:13).eq.'DIFF_LINE_INT') then
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c          return
c        else if (icnt.ne.5) then 
c          write(*,*) 'Error:  expected five integer arguments.'
c          return
c        else if (ccnt.lt.2) then
c          write(*,*) 'Error:  expected at least two string arguments.'
c          return
c        else
c          calc_field = diff_line_int(carg(1),carg(2),ccnt-1,iarg(1),
c     >           iarg(2),iarg(3),iarg(4),iarg(5),ndims,dims,stag,
c     >           min,max,missing,data_units,data_display_units,
c     >           dim_names)
c        end if
*
*  Line integral of pressure over y along terrain
*
c      else if (fld(1:13).eq.'LINE_INT_DRGX') then
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c          return
c        else if (icnt.ne.5) then
c          write(*,*) 'Error:  expected four integer arguments.'
c          return
c        else 
c          calc_field = line_int_drgx(iarg(1),iarg(2),iarg(3),iarg(4),
c     >           iarg(5),ndims,dims,stag,min,max,missing,data_units,
c     >           data_display_units,dim_names)
c        end if
*
*  Flux convergence into volume
*
      else if (fld(1:13).eq.'DIFF_AREA_INT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.6) then
          write(*,*) 'Error:  expected six integer arguments; ',
     >                      'got ',icnt
          return
        else if (ccnt.le.1) then
          write(*,*) 'Error:  expected more than one string ',
     >                     'argument.'
          return
        else
          calc_field = diff_area_int(carg(1),carg(2),ccnt-1,
     >           iarg(1),iarg(2),iarg(3),iarg(4),iarg(5),iarg(6),
     >           ndims,dims,stag,min,max,missing,data_units,
     >           data_display_units,dim_names)
        end if
*
*  Flux convergence into volume, vertically profiled
*
      else if (fld(1:11).eq.'DIFF_INT_VP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected five integer arguments; ',
     >                      'got ',icnt
          return
        else if (ccnt.le.1) then
          write(*,*) 'Error:  expected more than one string ',
     >                     'argument.'
          return
        else
          calc_field = diff_int_vp(carg(1),carg(2),ccnt-1,
     >           iarg(1),iarg(2),iarg(3),iarg(4),iarg(5),
     >           ndims,dims,stag,min,max,missing,data_units,
     >           data_display_units,dim_names)
        end if
*
*  Flux of momentum through top of volume (for budgets)
*
      else if (fld(1:6).eq.'ZF_TOP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected five integer arguments; ',
     >                            'got ',icnt
          return
        else
          calc_field = zf_top(iarg(1),iarg(2),iarg(3),iarg(4),
     >          iarg(5),ndims,dims,stag,min,max,missing,
     >          data_units,data_display_units,dim_names)
        end if
*
* Drag at bottom of volume (for budgets)
*
      else if (fld(1:7).eq.'DRG_BOT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.5) then
          write(*,*) 'Error:  expected five integer arguments; ',
     >                            'got ',icnt
          return
        else
          calc_field = drg_bot(iarg(1),iarg(2),iarg(3),iarg(4),
     >          iarg(5),ndims,dims,stag,min,max,missing,
     >          data_units,data_display_units,dim_names)
        end if
*
*  Average momentum loss per unit volume due to vertical fluxes/drag
*
c      else if (fld(1:6).eq.'RU_DRG') then
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c           return
c        else if (icnt.ne.7) then
c          write(*,*) 'Error:  expected seven integer arguments; ',
c     >                  'got ',icnt
c          return
c        else
c          calc_field = ru_drg(iarg(1),iarg(2),iarg(3),iarg(4),iarg(5),
c     >                      iarg(6),iarg(7),ndims,dims,stag,min,max,
c     >                      missing,data_units,data_display_units,
c     >                      dim_names)
c        end if
*
*  Volume average over z as function of x, y, and t
*
      else if (fld(1:5).eq.'Z_AVE') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.2) then
          write(*,*) 'Error:  expected two integer arguments; ',
     >                          'got ',icnt
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument; ',
     >                         'got ',ccnt
          return
        else
          calc_field = z_ave(carg(1),iarg(1),iarg(2),ndims,dims,stag,
     >            min,max,missing,data_units,data_display_units,
     >                 dim_names)
        end if
*
*  Volume average over x as function of y, z, and t
*
      else if (fld(1:5).eq.'X_AVE') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.2) then
          write(*,*) 'Error:  expected two integer arguments; ',
     >                           'got ',icnt
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument; ',
     >                          'got ',ccnt
          return
        else
          calc_field = h_ave('X',carg(1),iarg(1),iarg(2),ndims,dims,
     >                   stag,min,max,missing,data_units,
     >                   data_display_units,dim_names)
        end if
*
*  Volume average over y as function of x, z, and t
*
      else if (fld(1:5).eq.'Y_AVE') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.2) then
          write(*,*) 'Error:  expected two integer arguments; ',
     >                           'got ',icnt
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument; ',
     >                          'got ',ccnt
          return
        else
          calc_field = h_ave('Y',carg(1),iarg(1),iarg(2),ndims,dims,
     >                   stag,min,max,missing,data_units,
     >                   data_display_units,dim_names)
        end if
*
*  Time average
*
      else if (fld(1:5).eq.'T_AVE') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.2) then
          write(*,*) 'Error:  expected two integer arguments; ',
     >                           'got ',icnt
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument; ',
     >                          'got ',ccnt
          return
        else 
          calc_field = t_ave(carg(1),iarg(1),iarg(2),ndims,dims,
     >                    stag,min,max,missing,data_units,
     >                    data_display_units,dim_names)
        end if
*
*  General pointwise average
*
      else if (fld(1:3).eq.'AVE') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((ccnt.ne.2).or.(icnt.ne.2)) then
          write(*,*) 'Error:  expected two string arguments ',
     >                   'and two integer arguments.'
          return
        else
          calc_field = ave(carg(1),carg(2),iarg(1),iarg(2),ndims,
     >                  dims,stag,min,max,missing,data_units,
     >                  data_display_units,dim_names)
        end if
*
*  Perturbation momentum density
*
      else if (fld(1:5).eq.'RU_PR') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.1) then
          write(*,*) 'Error:  expected one integer argument; ',
     >                          'got ',icnt
          return
        else
          calc_field = ru_pr(iarg(1),ndims,dims,stag,min,max,
     >                   missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Zonal momentum flux
*
      else if (fld(1:4).eq.'UFLX') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  expected one string argument.'
          return
        else
          calc_field = uflx(carg(1),ndims,dims,stag,min,max,missing,
     >                    data_units,data_display_units,dim_names)
        end if
*
*  Absolute momentum, zonal component
*
      else if (fld(1:4).eq.'UABS') then
        calc_field = uabs(ndims,dims,stag,min,max,missing,
     >          data_units,data_display_units,dim_names)
*
*  Quotient of two fields at thermo pts
*
      else if (fld(1:4).eq.'QUOT') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.2) then
          write(*,*) 'Error:  expected two string arguments.'
          return
        else
          calc_field = quot(carg(1),carg(2),ndims,dims,stag,min,max,
     >                   missing,data_units,data_display_units,
     >                   dim_names)
        end if
*
*  Product of fields at thermo pts
*
      else if (fld(1:4).eq.'PROD') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.lt.1) then
          write(*,*) 'Error:  expected at least one string argument.'
          return
        else
          calc_field = prod(carg,ccnt,ndims,dims,stag,min,max,missing,
     >                  data_units,data_display_units,dim_names)
        end if
*
*  Sum of fields at thermo pts
*
c      else if (fld(1:3).eq.'SUM') then
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c          return
c        else if (ccnt.lt.1) then
c          write(*,*) 'Error:  expected at least one string argument.'
c          return
c        else
c          calc_field = sum(carg,ccnt,ndims,dims,stag,min,max,missing,
c     >                      data_units,data_display_units,dim_names)
c        end if
*
*  Line  
*
      else if (fld(1:4).eq.'LINE') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if ((ccnt.ne.1).or.(rcnt.ne.1)) then
          write(*,*) 'Error:  expected one string argument and ',
     >                          'one real argument.'
          return
        else
          calc_field = baseline(carg(1),rarg(1),ndims,dims,stag,
     >                   min,max,missing,data_units,
     >                   data_display_units,dim_names)
        end if
*
*  Symmetry switch
*
      else if (fld(1:10).eq.'SYMM_SWTCH') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.2) then
          write(*,*) 'Error:  expected two string arguments.'
          return
        else
          calc_field = symm_swtch(carg(2),carg(1),ndims,dims,stag,
     >                   min,max,missing,data_units,
     >                   data_display_units,dim_names)
        end if
*
*  Fluid volume (taking out terrain volume)
*
      else if (fld(1:6).eq.'VOLUME') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (icnt.ne.6) then
          write(*,*) 'Error:  expected six integer arguments.'
          return
        else
          calc_field = volume(iarg(1),iarg(2),iarg(3),iarg(4),
     >          iarg(5),iarg(6),ndims,dims,stag,min,max,missing,
     >          data_units,data_display_units,dim_names)
        end if
*
*  Divergence of velocity field (see div.F for options)
*
      else if (fld(1:3).eq.'DIV') then
        if (fld(4:4).eq.'[') then
          call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
          if (err) then
            return
          else if (icnt.gt.1) then
            write(*,*) 'Error:  expected one integer argument.'
            return
          else
            calc_field = div(iarg(1),ndims,dims,stag,min,max,missing,
     >               data_units,data_display_units,dim_names)
          end if
        else
          calc_field = div(0,ndims,dims,stag,min,max,missing,
     >                  data_units,data_display_units,dim_names)
        end if
*
*  Zonal momentum flux convergence
*
      else if (fld(1:4).eq.'UFCV') then
        calc_field = ufcv(ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names)
*
*  Volume density of pv substance
*
      else if (fld(1:3).eq.'PVS') then
        calc_field = pvs(ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names)
*
*  x component of vorticity
*
      else if (fld(1:5).eq.'VORTX') then
        calc_field = vortx(ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names)
*
*  x component of vorticity (4th order, therm pts)
*
      else if (fld(1:5).eq.'Q_IVE') then
        calc_field = q_ive(ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names)
*
*  y component of vorticity
*
      else if (fld(1:5).eq.'VORTY') then
        calc_field = vorty(ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names)
*
*  y component of vorticity (4th order, therm pts)
*
      else if (fld(1:5).eq.'R_IVE') then
        calc_field = r_ive(ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names)
*
*  z component of vorticity
*
      else if (fld(1:5).eq.'VORTZ') then
        calc_field = vortz(ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names)
*
*  z component of vorticity (4th order, therm pts)
*
      else if (fld(1:5).eq.'S_IVE') then
        calc_field = s_ive(ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names)
*
*  Bernoulli function (Boussinesq)
*
      else if (fld(1:9).eq.'BERNOULLI') then
        calc_field = bernoulli(ndims,dims,stag,min,max,missing,
     >             data_units,data_display_units,dim_names)
*
*  Basic state pressure (Boussinesq calculation)
*
      else if (fld(1:5).eq.'PBASE') then
        calc_field = pbase(ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names)
*
*  Basic state density (Boussinesq calculation)
*  
      else if (fld(1:3).eq.'R0B') then
        calc_field = r0b(ndims,dims,stag,min,max,missing,
     >        data_units,data_display_units,dim_names)
*
*  Perturbation density (Boussinesq)
*
      else if (fld(1:2).eq.'RB') then
        calc_field = rb(ndims,dims,stag,min,max,missing,
     >       data_units,data_display_units,dim_names)
*
*  Net Lagrangian displacement
*
      else if (fld(1:8).eq.'NET_DISP') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  need one string argument.'
          return
        else
          calc_field = net_disp(carg(1),ndims,dims,stag,min,max,
     >            missing,data_units,data_display_units,dim_names)
        end if
*
*  Uptream Height
*
      else if (fld(1:10).eq.'UPSTRM_HGT') then
          calc_field = upstrm_hgt(ndims,dims,stag,min,max,
     >            missing,data_units,data_display_units,dim_names)
*
*  Geometric distance from origin (x,y,z)
*
      else if (fld(1:5).eq.'DIST[') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.1) then
          write(*,*) 'Error:  need one string argument.'
          return
        else
          calc_field = dist(carg(1),ndims,dims,stag,min,max,missing,
     >           data_units,data_display_units,dim_names)
        end if
*
*  Extrapolation to surface
*
      else if (fld(1:4).eq.'SURF') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.1) then
          write(6,*) 'Error:  expected one string argument.'
          return
        else
          calc_field = surf(carg(1),ndims,dims,stag,min,max,
     >                  missing,data_units,data_display_units,
     >                  dim_names)
        end if
*
*  Copy time-independent field into time-dependent field
*
      else if (fld(1:8).eq.'MAKETIME') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (ccnt.ne.1) then
          write(6,*) 'Error:  expected one string argument.'
          return
        else
          calc_field = maketime(carg(1),ndims,dims,stag,min,max,
     >                     missing,data_units,data_display_units,
     >                     dim_names)
        end if
*
*  Make scalar
*
      else if (fld(1:6).eq.'MKSCLR') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (rcnt.ne.1) then
          write(6,*) 'Error:  need one real argument.'
          return
        else if (ccnt.ne.1) then
          write(6,*) 'Error:  need one string argument.'
          return
        else
          calc_field = scalar(rarg(1),carg(1),ndims,dims,stag,
     >             min,max,missing,data_units,data_display_units,
     >             dim_names)
        end if
*
*  Write mat file
*
c      else if (fld(1:8).eq.'MAKE_MAT') then
c        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
c        if (err) then
c          return
c        else if (ccnt.lt.2) then
c          write(*,*) 'Error:  expected at least two string ',
c     >                           'arguments.'
c          return
c        else if (rcnt.ne.1) then
c          write(*,*) 'Error:  expected one real argument.'
c          return
c        else if ((icnt.eq.0).or.(icnt.eq.2)) then
c          write(*,*) 'Error:  expected either one or three ',
c     >                  'integer arguments.'
c          return 
c        else
c          if (icnt.eq.1) then
c            iarg(2) = 1
c            iarg(3) = 1
c          end if
c          call make_mat(carg(1),carg(2),ccnt-1,rarg(1),
c     >                    iarg(1),iarg(2),iarg(3))
c          calc_field = 0
c          return
c        end if
*
*  Translate domain
*
      else if (fld(1:9).eq.'TRANSLATE') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else if (rcnt.ne.2) then
          write(6,*) 'Error:  expected two real arguments.'
          return
        else
          call translate_domain(rarg(1),rarg(2))
          calc_field = 0
          return
        end if
*
*  Turn vertical coordinate transformation on or off
*
      else if (fld(1:8).eq.'TRANS_ON') then
        call trans_switch(1)
        calc_field = 0
        return
      else if (fld(1:9).eq.'TRANS_OFF') then
        call trans_switch(0)
        calc_field = 0
        return
*
*  Output sounding
*
      else if (fld(1:9).eq.'SOUND_OUT') then
          call sound_out
          calc_field = 0
          return
c        end if
*
*  Test field
*
      else if (fld(1:4).eq.'TEST') then
        call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
        if (err) then
          return
        else
          tstfld = 'test'
          WRITE(*,*) mkfield(tstfld,carg,iarg,rarg,ccnt,icnt,rcnt,4)
          return
        end if
*
*  Variable not found
*
      else
        calc_field = 0
        return
      end if
*
*  Write out any scalars
*
      if (calc_field.ne.0) then
        if ((dims(1).eq.1).and.(dims(2).eq.1).and.(dims(3).eq.1)
     >                   .and.(dims(4).eq.1)) then
          write(6,*)
          write(6,*) 'Field ',fld(1:iend),' = ',
     >                    r_val(%val(calc_field),1)
          write(6,*)
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

