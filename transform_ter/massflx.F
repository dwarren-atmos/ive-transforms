#include "pointer.icl"

      ive_ptr function massflx(i1,i2,j1,j2,k1,k2,
     >                lateral,ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,k1,k2,nt
      real stag(4),min(4),max(4),miss
      character*(*) dim_names(4),dunits,disunits     
      logical lateral,oned

      real r_val
      ive_ptr getmem,fu,fwt,fwb,fvn,fvs
      ive_ptr area_int_x,area_int_z,area_int_y

      massflx = 0

      oned = .false.
      if (j1.eq.j2) then
       oned = .true.
       j1 = j1-1
      end if

      fu = area_int_x('u',1,j1,j2,k1,k2,0,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)

      nt = ntime 
      massflx = getmem(nt)
      
      if( lateral ) then
        if (oned) then
          fvn = area_int_y('dvdy',1,i1,i2,k1,k2,j2,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)

          call massflx_lat_calc_oned(%val(massflx),%val(fu),
     >                     %val(fvn),0,1,nt,dx)

        else
          fvn = area_int_y('v',1,i1,i2,k1,k2,j2,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)

          fvs = area_int_y('v',1,i1,i2,k1,k2,j1,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)

          call massflx_lat_calc(%val(massflx),%val(fu),
     >                     %val(fvn),%val(fvs),0,1,nt)
        endif
      else
        fwt = area_int_z('w',1,i1,i2,j1,j2,k2,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)

        fwb = area_int_z('w',1,i1,i2,j1,j2,k1,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)

        call massflx_top_calc(%val(massflx),%val(fu),
     >                     %val(fwt),%val(fwb),nt)

      end if

      return
      end

ccc
ccc Calculates mass going around mountain between two material surfaces
ccc
      ive_ptr function massflx_mat(i1,i2,j1,j2,z0,z1,
     >                ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2
      real stag(4),min(4),max(4),miss,z0,z1
      character*(*) dim_names(4),dunits,disunits     

      real r_val,fuinf
      integer nt,iw
	logical new,oned
      ive_ptr getmem,freemem,fu2,u0
      ive_ptr area_int_yz0_x
      ive_ptr area_int_xz0_y
c
      nt = dims(4)
      massflx_mat = 0
      massflx_mat = getmem(nt)
      if (massflx_mat.eq.0) then
        write(6,*) 'Error:  cannot allocate space for massflx_mat.'
	return
      end if
c
	u0 = 0
	u0 = getmem(dims(1)*dims(2)*dims(3))
      u0 = getvar('u0',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
    	if (u0.eq.0) then
	  write(*,*) 'Error:  cannot compute field u0'
        return
   	end if
c
c Assumses U0 is constant
c
	if (j1.ne.j2) then
	  fuinf = (j2-j1)*delta(2)*(z1-z0)*r_val(%val(u0),1)
c
	  fu2 = area_int_yz0_x('u',1,j1,j2,z0,z1,i2,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)
	  oned = .false.
	else
	  fuinf = (z1-z0)*r_val(%val(u0),1)
c
	  fu2 = area_int_yz0_x('u',1,j1-1,j2,z0,z1,i2,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)
	  oned = .true.
	end if
c
c	call freemem(u0)
      call massflx_lat_calc_sp(%val(massflx_mat),%val(fu2),
     >                     fuinf,0,1,nt,dy,oned)
c	 
      return
      end

ccc
ccc Calculates vertical profile of mass flux out lateral sides
ccc
      ive_ptr function massflx_vp(i1,i2,j1,j2,k1,k2,
     >                dpth,ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,k1,k2,dpth
      real stag(4),min(4),max(4),miss
      character*(*) dim_names(4),dunits,disunits     

      real r_val
      integer iw,nt
      ive_ptr getmem,fu,fvn,fvs
      ive_ptr area_int_x,area_int_y_vp

      if (mod(dpth,2).eq.0) then
        iw = 1
      else
        iw = 0
      end if

      nt = dims(4)
      massflx_vp = 0
      massflx_vp = getmem((nz+iw)*nt)
      if (massflx_vp.eq.0) then
        write(6,*) 'Error:  cannot allocate space for massflx_vp.'
	return
      end if

      fu = area_int_x('u',1,j1,j2,k1,k2,1,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)
      
      fvn = area_int_y_vp('v',1,i1,i2,dpth,j2,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)

      fvs = area_int_y_vp('v',1,i1,i2,dpth,j1,ndims,dims,
     >                         stag,min,max,miss,dunits,
     >                         disunits,dim_names)

      call massflx_lat_calc(%val(massflx_vp),%val(fu),
     >                     %val(fvn),%val(fvs),iw,nz,nt)

      return
      end
c
c Calculates ratio of mass entering a control volume to that going
c out the fixed lateral sides.
c
      subroutine massflx_lat_calc(massflx,fu,fvn,fvs,
     >             iw,nz,nt)    
      
      integer iw,nz,nt
      real massflx(nt),fu(nt),fvn(nt),fvs(nt)

      integer n,j
	
      do n=1,nt
      	
      	massflx(n) = ( fvn(n) - fvs(n) )/fu(n) 
	      print *,n,massflx(n)

      end do	

      end
c
c Calculates ratio of mass entering a control volume to that going
c out the fixed lateral sides for degenerate case of beta = 1 mountain.
c
      subroutine massflx_lat_calc_oned(massflx,fu,fvn,
     >             iw,nz,nt,dx)    
      
      integer iw,nz,nt
      real massflx(nz+iw,nt),fu(nt),fvn(nz+iw,nt),fvs(nz+iw,nt)
	real dx

      integer n,j
	
      do n=1,nt
       do j=1,nz+iw
      	
      	massflx(j,n) = ( fvn(j,n) )/( fu(n)/dx ) 

       end do
      end do	

      end
c
c Calculates ratio of mass entering a control volume to that going
c out the fixed lateral sides.
c
      subroutine massflx_lat_calc_sp(massflx,fu2,fusp,
     >             iw,nz,nt,dy,oned)    
      
      integer iw,nz,nt
      real massflx(nz+iw,nt),fu2(nt),fusp,dy
	logical oned

      integer n,j
	
	if (oned) then
        do n=1,nt
         do j=1,nz+iw
      	massflx(j,n) = 1.0 -  (fu2(n)/dy)/fusp 
         end do
        end do	
	else
        do n=1,nt
         do j=1,nz+iw
      	massflx(j,n) = 1.0 -  fu2(n)/fusp 
         end do
        end do	
	end if

      end
c
c Calculates ratio of mass entering a control volume to that going
c out the fixed height top.
c
      subroutine massflx_top_calc(massflx,fu,fwt,fwb,
     >             nt)    
      
      integer nt
      real massflx(nt),fu(nt),fwt(nt),fwb(nt)

      integer n
	
      do n=1,nt
      	
      	massflx(n) = (fwt(n)-fwb(n))/fu(n) 

      end do	

      end
