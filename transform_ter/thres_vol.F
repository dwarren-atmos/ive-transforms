c
#include "pointer.icl"
c
	ive_ptr function thres_vol(varnm,thres,zstop,i1,i2,j1,j2,
     >                 ndims,dims,stag,min,max,miss,
     >                 dunits,disunits,dim_names)
c
#include "default.icl"
c
	integer ndims,dims(4),i1,i2,j1,j2
	real stag(4),min(4),max(4),miss,thres,zstop
	character*(*) dunits,disunits,dim_names(4),varnm
	logical new,requ
c
	integer iu,iv,iw,nt,k,strbeg,strend	
	ive_ptr zs,var,getmem
c
      thres_vol = 0
	var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >             disunits,dim_names,new)
	if (var.eq.0) then
	   write(6,*) 'Error:  cannot get ',
     >              varnm( strbeg(varnm):strend(varnm) )
	   return
	end if
c
	thres_vol = getmem(dims(4))
	if (thres_vol.eq.0) then
	   write(6,*) 'Error:  cannot allocate memory for thres_vol '
	   return
	end if
c
	zs = zbot_p
	iu = 0
	iv = 0
	iw = 0
	if (requ(min(1),plmin(1))) then
	  zs = zbot_u
	  iu = 1
	else if ((requ(min(2),plmin(2))).and.(twod.ne.1)) then
	  zs = zbot_v
	  iv = 0
	else if (requ(min(3),plmin(3))) then
	  iw = 0
      end if
c
      if ( (j1.lt.1).or.(j2.gt.dims(2)).or.
     >     (i1.lt.1).or.(i2.gt.dims(1))) then
		write(6,*) 'Error:  specified limits', 
     >		     'exceed domain boundaries.'
	  return
	end if
c
	nt = dims(4)
	call calc_volume(%val(thres_vol),%val(var),thres,%val(zs),
     >     ztop,i1,i2,j1,j2,zstop,dx,dy,%val(sgz),%val(wgz),
     >     iu,iv,iw,nx,ny,nz,nt)
c
	do k=1,3
	  dims(k) = 1
	  stag(k) = 0.
	  min(k) = 0.
	  max(k) = 0.
	  dim_names(k) = 'one'
	enddo
c
	dunits = 'm*m*m'
	disunits = 'm*m*m'
c
	return
	end
c
c This subroutine calculates the volume within some larger defined 
c volume where var is within some threshold
c
	subroutine calc_volume(volume,var,thres,zs,zt,
     >           istart,istop,jstart,jstop,zstop,dx,dy,
     >           sgz,wgz,iu,iv,iw,nx,ny,nz,nt)
c
	integer nx,ny,nz,nt,istart,istop,jstart,jstop,iu,iv,iw
	real volume,var(nx+iu,ny+iv,nz+iw,nt),thres,
     >     zs(nx+iu,ny+iv,nz+iw),zt,
     >     zstop,sgz(nz),wgz(nz+1),dx,dy	
c
      integer i,j,n
      real zs1d,var1d(nz+iw),dv
c
	volume = 0
c	
	do n=1,nt
	 do i = istart,istop
	  do j = jstart,jstop
          call v_to_1d(zs,zs1d,i,j,n,iu,iv,iw,nx,ny,1,nt)
          call v_to_1d(var,var1d,i,j,n,iu,iv,iw,nx,ny,nz,nt)
	    call calc_dv(dv,var1d,thres,zt,zs1d,zstop,
     >		     wgz,sgz,dx,dy,iw,nz)
	    volume = volume + dv
        enddo
	 enddo
	enddo
c
	return
	end
c
c This subroutine calculates the volume in a column where var is 
c within some threshold.
c 
c BEWARE: DOES NOT SUPPORT VARIABLES AT W POINTS YET!
c
	subroutine calc_dv(dv,var,thres,zt,zs,z2,wgz,sgz,dx,dy,iw,nz)

	
	integer iw,nz
	real dv,var(nz),thres,zt,zs,z2,wgz(nz+1),sgz(nz),dx,dy

	integer k
	real trans,dz
	
	trans = (zt - zs)/zt
	dv = 0.
c
	do k = 1,nz
	  if ( z2.lt.(zs + trans*sgz(k)) ) go to 10
c
	  if (var(k).lt.thres) then
	    dz = trans*(wgz(k+1) - wgz(k))
	    dv = dv + dx*dy*dz
	  endif
	enddo
10    continue	
c
	return
	end
c
c This subroutine copies var(i,j,k,n) at point i,j,n to var1d(k) 
c  
	subroutine v_to_1d(var,var1d,i0,j0,n0,iu,iv,iw,nx,ny,nz,nt)

	integer i0,j0,n0,nx,ny,nz,nt,iu,iv,iw
	real var(nx+iu,ny+iv,nz+iw,nt),var1d(nz+iw)

	integer k,limit

	limit = nz+iw
	if (nz.eq.1) limit = 1

	do k=1,limit
	  var1d(k) = var(i0,j0,k,n0)
	enddo

	return
	end
