      integer(ive_ptr) function interp_height(varnm,z0,
     &         ndims,dims,stag,min,max,miss,dunits,
     &         disunits,dim_names)

      use default_module

      implicit none

      integer :: ndims,dims(4)
      real :: stag(4),min(4),max(4),miss
      character (len=*) :: varnm,dunits,disunits,dim_names(4)
      logical :: lagl

      integer :: iu,iv,iw,it,strbeg,strend,nxp1,nyp1,nt,i,k,mn
      real :: z0,r_val, sigm,zsfc
      character (len=80) :: outdunits,outdisunits
      logical :: new
      integer(ive_ptr) :: vartmp,var
      real, allocatable :: zhgt(:,:,:)

      integer, parameter :: ifin = 2

      nxp1 = nx+1; nyp1 = ny+1;
      mn = nxp1*nyp1
      vartmp = 0
      vartmp = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     $                          disunits,dim_names,new)
      if (vartmp.eq.0) then
        write(6,*) 'Error:  cannot access variable ',trim(varnm)
        interp_height=0
        return
      end if

      lagl=.false.
      if(z0.le.0) lagl=.true.
      z0 = abs(z0) 

      nt = dims(4)
      outdunits = dunits ; outdisunits = disunits

      call stagering_calc(iu,iv,iw,it,dims,plmin,min,dx,dy,twod)
c
c  Interpolate to the mass points
c
      var = 0 ; var = getmem(mn*nz*nt)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot allocate space for var.'
        interp_height=0
        return
      end if
c
c  set var to 1.
c
      do i=1,nxp1*nyp1*nz*nt
        call s_val(%val(var),i,1.)
      end do

      call var_at_t(%val(var),%val(vartmp),iu,iv,iw,nx,ny,nz,nt)

      interp_height = 0 ; interp_height = getmem(mn*nt)
      if (interp_height.eq.0) then
        write(6,*) 'Error:  cannot allocate space for interp_height.'
        return
      end if

       
      allocate(zhgt(mn,nz,nt))
      do k=1,nz
        sigm = r_val(%val(sgz),k)
        do i=1,mn 
          zsfc = r_val(%val(zbot_p),i);
          if(lagl) then
            zhgt(i,k,:) = (ztop - zsfc)/ztop * sigm
          else
            zhgt(i,k,:) = (ztop - zsfc)/ztop * sigm + zsfc
          endif
        end do 
      end do

      call vert_interp_field(%val(var),nxp1,nyp1,nz,nt,
     &          %val(zbot_p),ifin,%val(interp_height),
     &          1,z0,2,miss,zhgt)
      deallocate(zhgt)

      call freemem(var)
      dims(1) = nxp1 ; dims(2) = nyp1 ; dims(3) = 1 ; dims(4) = nt
      dim_names(1) = 'nx' ; dim_names(2) = 'ny'
      dim_names(3) = 'one' ; dim_names(4) = 'time'
      min(3) = z0 ;  max(3) = z0

      return
      end
