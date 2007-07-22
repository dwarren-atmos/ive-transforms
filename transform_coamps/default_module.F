      MODULE DEFAULT_MODULE
c
c          coord  an array containing the memory address of each array
c                 holding the "dimension variables" The use
c                 of dimension variables is optional, but
c                 allows one to specify irregular grid intervals
c                 along each spatial cooridinate
c          inmax  an array containing the maximum array indices of the 
c                 data 
c          phmin  an array of the minimum physical-space locations for
c                 each dimension
c          phmax  an array of the maximum physical-space locations for
c                 each dimension
c          delta  the uniform physical-space grid interval (if one exists) 
c                 along dimension
c
     
! Was previously in pointer.icl
#  ifdef SIXTY_FOUR_bptr
      integer, parameter :: ive_ptr=8
#  else
      integer, parameter :: ive_ptr=4
#  endif

! Was previously in default.icl
      integer(ive_ptr) :: getvar, getmem, coord(4), time, wgz, sgz, lgz
      integer(ive_ptr) :: zbot_p, zbot_u, zbot_v
      integer(ive_ptr) :: fcor
      integer(ive_ptr) :: dqdx_xy, dqdx_zt, dqdx_zw
      integer(ive_ptr) :: dqdy_xy, dqdy_zt, dqdy_zw
      integer(ive_ptr) :: dqdz, dzdq
      integer(ive_ptr) :: hmap, grdrot

      integer :: trans_z,trans_on
      integer :: nx,ny,nz,twod
      integer :: ntime, nwgz, nsgz, nlgz
      integer :: zbot_dims(2)
      integer :: igrid, jgrid, iref, jref

      real :: inmax(4),phmin(4),phmax(4),delta(4)
      real :: offset(2)
      real :: dx,dy,dz,plmin(4),plmax(4)
      real :: wgz_max,wgz_min,sgz_max,sgz_min,lgz_min,lgz_max
      real :: zbot_delta(2),zbot_min(4),zbot_max(4),ztop
      real :: reflat,reflon,stdlt1,stdlt2,stdlon
      real :: get_ter_pt

      external :: getvar, getmem, get_ter_pt

      integer, parameter :: MAX_LEN = 80, dim_len = 4
      TYPE meta_ive 
        integer :: ndims, dims(dim_len) 
        real :: stag(dim_len), max_val(dim_len), min_val(dim_len), miss
        character(len=MAX_LEN) :: dunits, disunits, dim_names(dim_len)
      END TYPE

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Sets the meta data associated with an ive variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_var_meta(stagpts,meta)
c     &                 ,ndims,dims,stag,min,max
c     &                 ,miss,dim_names)

      !USE DEFAULT_MODULE

      character(len=*), intent(in) :: stagpts

      type(meta_ive), intent(out) :: meta

      !integer, intent(out) :: ndims, dims(4)
      !real, intent(out) :: stag(4), min(4), max(4), miss
      !character(len=*), intent(out) :: dim_names(4)

      real :: r_val

      ! SET THE DEFAULTS OF ALL OUTPUT VARIABLES
      meta%dim_names(1) = 'nx' ; meta%dim_names(2) = 'ny'
      meta%dim_names(3) = 'nz' ; meta%dim_names(4) = 'time'
      meta%dims(1) = nx+1 ; meta%dims(2) = ny+1 
      meta%dims(3) = nz ; meta%dims(4) = ntime

      meta%max(1) = real(nx)*dx ; max(2) = real(ny)*dy
      meta%min(3) = sgz_min ; meta%max(3) = (ztop)/real(nz+1)*real(nz) 
      !meta%max(3) = sgz_max 
      max(4) = r_val(%val(time),ntime)
      min(:) = 0.0 ! set default min value to 0
      metat%stag(:) = 0.0 ; meta%ndims = 4 ; meta%miss = -9999.0 

      select case (trim(stagpts))
      case ('U')
        meta%dims(1) = nx ; meta%dim_names(1) = 'nxm1' 
        meta%min(1) = 0.5*dx ; meta%max(1) = (real(nx)-0.5)*dx
        meta%stag(1) = 0.0

      case ('V')
        meta%dims(2) = ny ; meta%dim_names(2) = 'nym1'
        meta%min(2) = 0.5*dy ; meta%max(2) = (real(ny)-0.5)*dy
        meta%stag(2) = 0.0

      case ('W')
        meta%dims(3) = nz + 1 ; meta%dim_names(3) = 'nzp1' 
        meta%min(3) = 0.0 ; meta%max(4) = ztop
        meta%stag(3) = 0.0

      end select

      return
      end subroutine
      END MODULE
