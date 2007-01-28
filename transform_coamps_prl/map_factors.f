	subroutine map_factors(hmap,grdrot,fm,nx,ny,dx,dy,iref,jref,
     >			  igrid,reflat,reflon,stdlt1,stdlt2,stdlon)
	implicit none
c
c  nx,ny are number of thermo points here
c
	integer nx,ny,igrid,iref,jref
	real dx,dy,reflat,reflon,stdlt1,stdlt2,stdlon
	real hmap(nx,ny,3)
c
	integer istr1,istr2,strgrd,i,j
	real hxm(nx,ny),hym(nx,ny)
      real hxv(nx,ny),hyv(nx,ny),hxu(nx,ny),hyu(nx,ny)
c
	real grdi(nx,ny),grdj(nx,ny),grdlat(nx,ny),grdlon(nx,ny),
     >     fm(nx,ny),xpos(nx,ny),ypos(nx,ny),distx(nx,ny),
     >     disty(nx,ny),dxav1(nx,ny),dxav2(nx,ny),dxm(nx,ny),
     >     dxu(nx,ny),dxv(nx,ny),dyav1(nx,ny),dyav2(nx,ny),
     >     dym(nx,ny),dyu(nx,ny),dyv(nx,ny),grdrot(nx,ny),
     >     flat,gcon
c
c  Map factors
c    1: therm pts
c    2: u pts
c    3: v pts
c
	flat = 0.
	istr1 = 1
	istr2 = 1
	strgrd = 1
	call grdij(nx,ny,grdi,grdj)
	call grid(igrid,reflat,reflon,iref,jref,stdlt1,stdlt2,
     >	    stdlon,dx,dx,dy,dy,grdi,grdj,nx,ny,
     >	    grdlat,grdlon,gcon,fm,flat,hxm,hym,xpos,ypos,
     >	    distx,disty,dxav1,dxav2,dxm,dxu,dxv,dyav1,
     >	    dyav2,dym,dyu,dyv,grdrot,istr1,istr2,strgrd)
c
	call hm2uv(hxm,hym,hxu,hyu,hxv,hyv,nx,ny)
c
	do i=1,nx
	  do j = 1,ny
	    hmap(i,j,1) = hxm(i,j)
	    if(i.lt.nx) hmap(i,j,2) = hxu(i,j)
	    if(j.lt.ny) hmap(i,j,3) = hxv(i,j)
	  enddo
	enddo

	    
	return
	end
