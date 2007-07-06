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

      END MODULE
