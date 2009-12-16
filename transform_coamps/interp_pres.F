      integer(ive_ptr) function interp_pres(varnm,pvarnam,p0,nflt,
     &        ndims,dims,stag,min,max,miss,dunits,
     &        disunits,dim_names)
c
      use default_module
      implicit none
c
      integer ndims,dims(4),nflt
      real stag(4),min(4),max(4),miss
      character*(*) varnm,pvarnam,dunits,disunits,dim_names(4)
c      
      integer iu,iv,iw,it,strbeg,strend,nxp1,nyp1,nt,i
      real p0,r_val
      character*(80) outdunits,outdisunits
      logical new
      integer(ive_ptr) vartmp,pvar,var
      integer, parameter :: ifin = 3
c
      nxp1 = nx+1 ; nyp1 = ny+1
c
      vartmp = 0
      vartmp = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new)
      print *,r_val(%val(vartmp),202500),miss
c
c Get variable to interp_presolate
c
      if (vartmp.eq.0) then
        write(6,*) 'Error:  cannot access variable ',
     &               varnm(strbeg(varnm):strend(varnm))
        return
      end if
c
      nt = dims(4)
      outdunits = dunits ; outdisunits = disunits
c
      call stagering_calc(iu,iv,iw,it,dims,plmin,min,dx,dy,twod)
c
c  Interpolate to the mass points
c
      var = 0 ; var = getmem(nxp1*nyp1*nz*nt)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot allocate space for var.'
        return
      end if
c
c  set var to 1.
c
      do i=1,nxp1*nyp1*nz*nt
         call s_val(%val(var),i,1.)
      end do
c
      call var_at_t(%val(var),%val(vartmp),iu,iv,iw,nx,ny,nz,nt)
c
c Get pressure field
c
      pvar = 0
      pvar = getvar(pvarnam,ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new)
c
      if (pvar.eq.0) then
        write(6,*) 'Error:  cannot access variable PPP'
        return
      end if
c
      interp_pres = 0 ; interp_pres = getmem(nxp1*nyp1*nt)
      if (interp_pres.eq.0) then
        write(6,*) 'Error:  cannot allocate space for interp_pres.'
        return
      end if
c
      call vert_interp_field(%val(var),nxp1,nyp1,nz,nt,
     &   %val(zbot_p),ifin,%val(interp_pres),
     &   1,p0,2,miss,%val(pvar))

      do i=1,nflt
        call filt9_iso(%val(interp_pres),nxp1,nyp1,ntime)
      end do
c
      call freemem(var)
c
      dims(1) = nxp1; dims(2) = nyp1; dims(3) = 1; dims(4) = nt
      dim_names(1) = 'nx'; dim_names(2) = 'ny'
      dim_names(3) = 'one'; dim_names(4) = 'time'
      min(3) = 0; max(3) = 0
      dunits = outdunits ; disunits = outdisunits
c
      return
      end
