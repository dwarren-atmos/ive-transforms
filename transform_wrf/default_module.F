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
      save
#  ifdef SIXTY_FOUR_bptr
      integer, parameter :: ive_ptr=8
#  else
      integer, parameter :: ive_ptr=4
#  endif

      integer(ive_ptr) :: getvar, addvar, getmem, coord(4), time
      integer(ive_ptr) :: znu_ptr, znw_ptr, ptop_ptr, psfc_ptr
      integer(ive_ptr) :: mutot_ptr, phtot_ptr, hgt_ptr
      integer(ive_ptr) :: mapM_ptr, mapU_ptr, mapV_ptr
      integer(ive_ptr) :: lat_ll, lat_ur, lon_ll, lon_ur
      integer(ive_ptr) :: fcor
      integer(ive_ptr) :: dqdx_xy, dqdx_zt, dqdx_zw
      integer(ive_ptr) :: dqdy_xy, dqdy_zt, dqdy_zw
      integer(ive_ptr) :: dqdz, dzdq
      integer(ive_ptr) :: grdrot

      integer :: getvid
      integer :: trans_eta,trans_on
      integer :: nx,ny,nz,twod, nzp1
      integer :: ntime, nznw, nznu
      integer :: zbot_dims(2)
      integer :: igrid, jgrid, iref, jref

      real :: znu_min,znu_max,znw_min,znw_max
      real :: inmax(4),phmin(4),phmax(4),delta(4)
      real :: offset(4)
      real :: dx,dy,dz,plmin(4),plmax(4)
      real :: wgz_max,wgz_min,sgz_max,sgz_min,lgz_min,lgz_max
      real :: zbot_delta(2),zbot_min(4),zbot_max(4),ztop
      real :: phys_mumin(4), phys_mumax(4), mu_dims(4), mu_delta(4)
      real :: ptop
      real :: reflat,reflon,stdlt1,stdlt2,stdlon
      real :: get_ter_pt, get_2Dpt_time
      real :: mumax, mumin
      real :: cenlat, cenlon, trulat1, trulat2, standlon

      character(len=80) :: domain_units(4), domain_display_units(4)
      character(len=80) :: dim_display_names(4)
      character(len=2)  :: projection

      external :: getvar, addvar,getmem,get_ter_pt,getvid,get_2Dpt_time

      integer, parameter :: MAX_LEN = 80, dim_len = 4
      TYPE meta_ive 
        integer :: ndims, dims(dim_len) 
        real :: stag(dim_len), max(dim_len), min(dim_len), miss
        character(len=MAX_LEN) :: dunits, disunits, dim_names(dim_len)
      END TYPE

      real, parameter :: pref=10**5, Rd=287.04, cp=1004., cv=717.,
     &                   kappa = Rd/cp

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Sets the button names and button display
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_buttons

      integer ierr, set_button_name

      ! raw netcdf variables to display
      ierr = set_button_name ('U', 'U')
      ierr = set_button_name ('V', 'V')
      ierr = set_button_name ('W', 'W')
      ierr = set_button_name ('QVAPOR', 'QVAPOR')

      ierr = set_button_name ('HGT', '') ! Terrain

      ierr = set_button_name ('PH', '')  ! Pert geopotential
      ierr = set_button_name ('PHB', '') ! Base geopotential

      ierr = set_button_name ('MU', '')  ! Pert dry-air mass
      ierr = set_button_name ('MUB', '') ! Base dry-air mass

      ierr = set_button_name ('P', '')   ! Pert pressure
      ierr = set_button_name ('PB', '')  ! Base pressure
      ierr = set_button_name ('P_TOP', '') ! Top of model pressure

      ierr = set_button_name ('ZNU', '') ! Mass level eta-values
      ierr = set_button_name ('ZNW', '') ! W level eta-valuse

      ierr = set_button_name ('SMSTAV', '')
      ierr = set_button_name ('SFROFF', '')
      ierr = set_button_name ('UDROFF', '')
      ierr = set_button_name ('CANWAT', '')
      ierr = set_button_name ('ACSNOM', '')
      ierr = set_button_name ('SNOW', '')
      ierr = set_button_name ('WEASD', '')
      ierr = set_button_name ('THZ0', '')
      ierr = set_button_name ('QZ0', '')
      ierr = set_button_name ('UZ0', '')
      ierr = set_button_name ('VZ0', '')
      ierr = set_button_name ('AKHS', '')
      ierr = set_button_name ('AKMS', '')
      ierr = set_button_name ('FNM', '')
      ierr = set_button_name ('FNP', '')
      ierr = set_button_name ('RDNW', '')
      ierr = set_button_name ('RDN', '')
      ierr = set_button_name ('DNW', '')
      ierr = set_button_name ('DN', '')

      ierr = set_button_name ('RU_TEND', '')
      ierr = set_button_name ('RU_TENDF', '')
      ierr = set_button_name ('RV_TEND', '')
      ierr = set_button_name ('RV_TENDF', '')
      ierr = set_button_name ('RW_TEND', '')
      ierr = set_button_name ('RW_TENDF', '')
      ierr = set_button_name ('T_TEND', '')
      ierr = set_button_name ('T_TENDF', '')

      ierr = set_button_name ('NEST_POS', '')
      ierr = set_button_name ('SR', '')
      ierr = set_button_name ('POTEVP', '')
      ierr = set_button_name ('SNOPCX', '')
      ierr = set_button_name ('SOILTB', '')
      ierr = set_button_name ('CFN', '')
      ierr = set_button_name ('CFN1', '')
      ierr = set_button_name ('RDX', '')
      ierr = set_button_name ('RDY', '')
      ierr = set_button_name ('RESM', '')
      ierr = set_button_name ('CF1', '')
      ierr = set_button_name ('CF2', '')
      ierr = set_button_name ('CF3', '')
      ierr = set_button_name ('XTIME', '')
      ierr = set_button_name ('LANDMASK', '')
      ierr = set_button_name ('XICE', '')
      ierr = set_button_name ('MAPFAC_U', '')
      ierr = set_button_name ('MAPFAC_V', '')
      ierr = set_button_name ('SINALPHA', '')
      ierr = set_button_name ('COSALPHA', '')
      ierr = set_button_name ('TSK', '')

      ierr = set_button_name ('LAT_LL_T', '')
      ierr = set_button_name ('LAT_UL_T', '')
      ierr = set_button_name ('LAT_UR_T', '')
      ierr = set_button_name ('LAT_LR_T', '')

      ierr = set_button_name ('LON_LL_T', '')
      ierr = set_button_name ('LON_UL_T', '')
      ierr = set_button_name ('LON_UR_T', '')
      ierr = set_button_name ('LON_LR_T', '')

      ierr = set_button_name ('LAT_LL_U', '')
      ierr = set_button_name ('LAT_UL_U', '')
      ierr = set_button_name ('LAT_UR_U', '')
      ierr = set_button_name ('LAT_LR_U', '')

      ierr = set_button_name ('LON_LL_U', '')
      ierr = set_button_name ('LON_UL_U', '')
      ierr = set_button_name ('LON_UR_U', '')
      ierr = set_button_name ('LON_LR_U', '')

      ierr = set_button_name ('LAT_LL_V', '')
      ierr = set_button_name ('LAT_UL_V', '')
      ierr = set_button_name ('LAT_UR_V', '')
      ierr = set_button_name ('LAT_LR_V', '')

      ierr = set_button_name ('LON_LL_V', '')
      ierr = set_button_name ('LON_UL_V', '')
      ierr = set_button_name ('LON_UR_V', '')
      ierr = set_button_name ('LON_LR_V', '')

      ierr = set_button_name ('LAT_LL_D', '')
      ierr = set_button_name ('LAT_UL_D', '')
      ierr = set_button_name ('LAT_UR_D', '')
      ierr = set_button_name ('LAT_LR_D', '')

      ierr = set_button_name ('LON_LL_D', '')
      ierr = set_button_name ('LON_UL_D', '')
      ierr = set_button_name ('LON_UR_D', '')
      ierr = set_button_name ('LON_LR_D', '')

      ierr = set_button_name ('RAINC', '')
      ierr = set_button_name ('RAINNC', '')
      ierr = set_button_name ('RANCV', '')
      ierr = set_button_name ('RAINBL', '')
      ierr = set_button_name ('ZETATOP', '')

      ierr = set_button_name ('XLONG', '') ! Longitude points
      ierr = set_button_name ('XLAT', '')  ! Latitude points

      ierr = set_button_name ('XLONG_U', '') ! U-longitude points
      ierr = set_button_name ('XLAT_U', '')  ! U-latitude points

      ierr = set_button_name ('XLONG_V', '') ! V-longitude points
      ierr = set_button_name ('XLAT_V', '')  ! V-latitude points

      ierr = set_button_name ('MAPFAC_M', '')

      ierr = set_button_name ('TIMES', '')
      ierr = set_button_name ('SST', '')
      ierr = set_button_name ('E', '')
      ierr = set_button_name ('F', '')
      ierr = set_button_name ('ITIMESTEP', '')
      ierr = set_button_name ('IVGTYP', '')
      ierr = set_button_name ('ISLTYP', '')
      ierr = set_button_name ('VEGFRA', '')
      ierr = set_button_name ('LU_INDEX', '')
      ierr = set_button_name ('XLAND', '')
      ierr = set_button_name ('T', '')

      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Sets the meta data associated with an ive variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_var_meta(stagpts
     &                 ,ndims,dims,stag,min,max
     &                 ,miss,dim_names)

      !USE DEFAULT_MODULE

      character(len=*), intent(in) :: stagpts

      !type(meta_ive), intent(out) :: meta

      integer, intent(out) :: ndims, dims(4)
      real, intent(out) :: stag(4), min(4), max(4), miss
      character(len=*), intent(out) :: dim_names(4)

      real :: r_val

      ! SET THE DEFAULTS OF ALL OUTPUT VARIABLES
      dim_names(1) = 'nx' ; dim_names(2) = 'ny'
      dim_names(3) = 'nz' ; dim_names(4) = 'time'
      dims(1) = nx+1 ; dims(2) = ny+1 
      dims(3) = nz ; dims(4) = ntime

      min(:) = 0.0 ! set default min value to 0
      max(1) = real(nx)*dx ; max(2) = real(ny)*dy
      min(3) = sgz_min ; max(3) = (ptop)/real(nz+1)*real(nz) 
      !max(3) = sgz_max 
      max(4) = r_val(%val(time),ntime)
      stag(:) = 0.0 ; ndims = 4 ; miss = -9999.0 

      select case (trim(stagpts))
      case ('U')
        dims(1) = nx ; dim_names(1) = 'nxm1' 
        min(1) = 0.5*dx ; max(1) = (real(nx)-0.5)*dx
        stag(1) = 0.0

      case ('V')
        dims(2) = ny ; dim_names(2) = 'nym1'
        min(2) = 0.5*dy ; max(2) = (real(ny)-0.5)*dy
        stag(2) = 0.0

      case ('W')
        dims(3) = nz + 1 ; dim_names(3) = 'nzp1' 
        min(3) = 0.0 ; max(3) = ptop
        stag(3) = 0.0

      end select

      return
      end subroutine

      END MODULE
