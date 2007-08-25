!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   function vort_z(uvar_name, vvar_name, ...)
!
!   Calculates the vertical component of vorticity 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(ive_ptr) function vort_z(uvar_name,vvar_name,
     &                                 ndims,dims,stag,min,max,miss,
     &                                 dunits,disunits,dim_names)

      USE DEFAULT_MODULE

      character(len=*), intent(in) :: uvar_name, vvar_name

      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout):: dunits,disunits,dim_names(4)
      
      integer :: i, mnk
      integer(ive_ptr) :: d_uvar, d_vvar, d_uvar_t, d_vvar_t
      real :: r_val
      logical :: new,error

      integer(ive_ptr) :: deriv

      vort_z=0 
      mnk = ntime*nz*(ny+1)*(nx+1)

      d_uvar = deriv(uvar_name,'Y',1,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (d_uvar.le.0) return

      d_vvar = deriv(vvar_name,'X',1,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (d_vvar.le.0) return

      ! get memory for variables
      d_uvar_t = getmem(mnk) 
      if(d_uvar_t.le.0) then
        print *,'Error: cannont allocate memory for d_uvar_t in vort_z'
        vort_z=0
        return
      end if
      
      d_vvar_t = getmem(mnk)
      if(d_vvar_t.le.0) then
        print *,'Error: cannont allocate memory for d_vvar_t in vort_z'
        vort_z=0
        return
      end if


      vort_z = getmem(mnk)
      if(d_uvar_t.le.0 .or. d_uvar_t.le.0 .or. vort_z.le.0) then
        print *,'Error: cannont allocate memory for vort_z in vort_z'
        vort_z=0
        return
      end if

      call var_at_t(%val(d_uvar_t),%val(d_uvar),1,0,0,nx,ny,nz,ntime)
      call var_at_t(%val(d_vvar_t),%val(d_vvar),0,1,0,nx,ny,nz,ntime)

      do i=1,mnk
        call s_val(%val(vort_z),i,
     &             r_val(%val(d_vvar_t),i)-r_val(%val(d_uvar_t),i))
      end do

      dunits='1/s' ; disunits='1/s'
      dims(1) = nx+1 ; dims(2) = ny+1 ; dims(3) = nz
      dim_names(1)='nx' ; dim_names(2)='ny' ; dim_names(3)='nz' 
      stag(1) = 0.0 ; stag(2) = 0.0 ; stag(3) = 0.0
      min(1) = plmin(1) ; min(2) = plmin(2)
      max(1) = plmax(1) ; max(2) = plmax(2)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   function vort_x(vvar_name, wvar_name, ...)
!
!   Calculates the x component of vorticity 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer(ive_ptr) function vort_x(vvar_name,wvar_name,
     &                                 ndims,dims,stag,min,max,miss,
     &                                 dunits,disunits,dim_names)

      USE DEFAULT_MODULE

      character(len=*), intent(in) :: wvar_name, vvar_name

      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout):: dunits,disunits,dim_names(4)
      
      integer :: i, mnk
      integer(ive_ptr) :: d_wvar, d_vvar, d_wvar_t, d_vvar_t
      real :: r_val
      logical :: new,error

      integer(ive_ptr) :: deriv

      vort_x=0 
      mnk = ntime*nz*(ny+1)*(nx+1)

      d_wvar = deriv(wvar_name,'Y',1,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (d_wvar.le.0) return

      d_vvar = deriv(vvar_name,'Z',1,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (d_vvar.le.0) return

      ! get memory for variables
      d_wvar_t = getmem(mnk) ; d_vvar_t = getmem(mnk)
      vort_x = getmem(mnk)
      if(d_vvar_t.le.0 .or. d_wvar_t.le.0 .or. vort_x.le.0) then
        print *,'Error: cannont allocate nesecary memory in vort_x'
        vort_x=0
        return
      end if

      call var_at_t(%val(d_vvar_t),%val(d_vvar),0,1,0,nx,ny,nz,ntime)
      call var_at_t(%val(d_wvar_t),%val(d_wvar),0,0,1,nx,ny,nz,ntime)

      do i=1,mnk
        call s_val(%val(vort_x),i,
     &             r_val(%val(d_wvar_t),i)-r_val(%val(d_vvar_t),i))
      end do

      dunits='1/s' ; disunits='1/s'
      dims(1) = nx+1 ; dims(2) = ny+1 ; dims(3) = nz
      dim_names(1)='nx' ; dim_names(2)='ny' ; dim_names(3)='nz' 
      stag(1) = 0.0 ; stag(2) = 0.0 ; stag(3) = 0.0
      min(1) = plmin(1) ; min(2) = plmin(2)
      max(1) = plmax(1) ; max(2) = plmax(2)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   function vort_y(uvar_name, wvar_name, ...)
!
!   Calculates the y component of vorticity 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer(ive_ptr) function vort_y(uvar_name,wvar_name,
     &                                 ndims,dims,stag,min,max,miss,
     &                                 dunits,disunits,dim_names)

      USE DEFAULT_MODULE

      character(len=*), intent(in) :: wvar_name, uvar_name

      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout):: dunits,disunits,dim_names(4)
      
      integer :: i, mnk
      integer(ive_ptr) :: d_wvar, d_uvar, d_wvar_t, d_uvar_t
      real :: r_val
      logical :: new,error

      integer(ive_ptr) :: deriv

      vort_y=0 
      mnk = ntime*nz*(ny+1)*(nx+1)

      d_wvar = deriv(wvar_name,'X',1,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (d_wvar.le.0) return

      d_uvar = deriv(uvar_name,'Z',1,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (d_uvar.le.0) return

      ! get memory for variables
      d_wvar_t = getmem(mnk) ; d_uvar_t = getmem(mnk)
      vort_y = getmem(mnk)
      if(d_uvar_t.le.0 .or. d_wvar_t.le.0 .or. vort_y.le.0) then
        print *,'Error: cannont allocate nesecary memory in vort_y'
        vort_y=0
        return
      end if

      call var_at_t(%val(d_uvar_t),%val(d_uvar),1,0,0,nx,ny,nz,ntime)
      call var_at_t(%val(d_wvar_t),%val(d_wvar),0,0,1,nx,ny,nz,ntime)

      do i=1,mnk
        call s_val(%val(vort_y),i,
     &             r_val(%val(d_uvar_t),i)-r_val(%val(d_wvar_t),i))
      end do

      dunits='1/s' ; disunits='1/s'
      dims(1) = nx+1 ; dims(2) = ny+1 ; dims(3) = nz
      dim_names(1)='nx' ; dim_names(2)='ny' ; dim_names(3)='nz' 
      stag(1) = 0.0 ; stag(2) = 0.0 ; stag(3) = 0.0
      min(1) = plmin(1) ; min(2) = plmin(2)
      max(1) = plmax(1) ; max(2) = plmax(2)

      return
      end

