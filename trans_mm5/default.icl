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

        integer nx,ny,nz,twod
        real dx,dy,dz,plmin(4),plmax(4)
        common/t_domain/nx,ny,nz,twod,dx,dy,dz,plmin,plmax

        ive_ptr xter, dter, ter_p
        real get_ter_pt, get_ght
        integer xter_dims(2), index_3d, dter_dims(2)
        real xter_delta(2),xter_min(4),xter_max(4),ztop
        real dter_delta(2),dter_min(4),dter_max(4),trans_z
        common/n_ter/xter,xter_delta,xter_min,ter_p,
     >           xter_max,ztop,xter_dims,trans_z,
     >      dter,dter_delta,dter_min,dter_max,dter_dims

        ive_ptr xght,dght,hsigma,fsigma,dcor,ps0,ts0,tlp,ptop
        ive_ptr xlat,xlon,pstx,tgrid,xtheta,xptot,xwspd,zght
        integer nzhght, nzdght
        real xght_min,xght_max,hsig_min,hsig_max,offset(2)
        real xght_delta,xght_dims(3),xght_mind(2)
        real dght_delta,dght_dims(3),dght_mind(2)
        real dght_min,dght_max
        common/t_xght/xght,hsigma,fsigma,dcor,ps0,ts0,tlp,ptop,
     >      dght,xlat,xlon,pstx,zght,offset,tgrid,xght_delta,xght_dims,
     >      xght_mind,dght_delta,dght_dims,dght_mind,
     >      dght_min,dght_max

