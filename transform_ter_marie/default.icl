#include "pointer.icl"
c
c          COMMON BLOCK FOR USE WITH THE DEFAULT TRANSFORMS
c
c          coord  an array containing the memory address of each array
c                 holding the "dimension variables" The use
c                 of dimension variables is optional, but
c                 allows one to specify irregular grid intervals
c                 along each spatial cooridinate
c
c          inmax  an array containing the maximum array indices of the 
c                 data 
c
c          phmin  an array of the minimum physical-space locations for
c                 each dimension
c
c          phmax  an array of the maximum physical-space locations for
c                 each dimension
c
c          delta  the uniform physical-space grid interval (if one exists) 
c                 along dimension
c       
c
	ive_ptr coord(4)
	real    inmax(4), phmin(4), phmax(4), delta(4)	

	common/t_default/coord, inmax, phmin, phmax, delta

	ive_ptr getvar
	external getvar
	real get_ter_pt
	external get_ter_pt

	ive_ptr wgz,sgz,zbot_p,time
	real	wgz_min,wgz_max,sgz_min,sgz_max
        real    zbot_delta(2),ztop,zbot_min(4),zbot_max(4)
	integer	nwgz,nsgz,ntime,zbot_dims(2)

        common/timeblk/time,ntime

        common/t_gz/wgz,sgz,wgz_max,wgz_min,sgz_max,sgz_min,
     >                 nwgz,nsgz

        common/t_ter/zbot_p,zbot_delta,zbot_min,zbot_max,ztop,
     >               zbot_dims


