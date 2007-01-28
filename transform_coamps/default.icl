#include "pointer.icl"
c
c          COMMON BLOCK FOR USE WITH THE DEFAULT TRANSFORMS
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
      ive_ptr getvar
      external getvar
	  ive_ptr getmem
	  external getmem
      real get_ter_pt
      external get_ter_pt

      ive_ptr coord(4)
      real inmax(4),phmin(4),phmax(4),delta(4)	
      common/t_default/coord,inmax,phmin,phmax,delta

      integer trans_z,trans_on
      real offset(2)
      common/t_flgs/trans_z,trans_on,offset

      integer nx,ny,nz,twod
      real dx,dy,dz,plmin(4),plmax(4)
      common/t_domain/nx,ny,nz,twod,dx,dy,dz,plmin,plmax
c      real dx,dy,dz,plmin(4),plmax(4),fcor
c      common/t_domain/nx,ny,nz,twod,dx,dy,dz,plmin,plmax,fcor
      
      ive_ptr time
      integer ntime
      common/timeblk/time,ntime

      ive_ptr wgz,sgz,lgz
      integer nwgz,nsgz,nlgz
      real wgz_max,wgz_min,sgz_max,sgz_min,lgz_min,lgz_max
      common/t_gz/wgz,sgz,wgz_max,wgz_min,sgz_max,sgz_min,
     >            nwgz,nsgz,lgz,nlgz,lgz_min,lgz_max

      ive_ptr zbot_p,zbot_u,zbot_v
      integer zbot_dims(2)
      real zbot_delta(2),zbot_min(4),zbot_max(4),ztop
      common/t_ter/zbot_p,zbot_u,zbot_v,zbot_delta,zbot_min,
     >             zbot_max,ztop,zbot_dims

      ive_ptr fcor
	common/f_cor/ fcor

      ive_ptr dqdx_xy,dqdx_zt,dqdx_zw
      ive_ptr dqdy_xy,dqdy_zt,dqdy_zw
      ive_ptr dqdz,dzdq
      common/t_metric/dqdx_xy,dqdx_zt,dqdx_zw,dqdy_xy,dqdy_zt,
     >                   dqdy_zw,dqdz,dzdq

      ive_ptr hmap,grdrot
      integer igrid,jgrid,iref,jref
	real reflat,reflon,stdlt1,stdlt2,stdlon
	common/map_proj/igrid,reflat,reflon,stdlt1,stdlt2,stdlon,
     >			hmap,grdrot,iref,jref
