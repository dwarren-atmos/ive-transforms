
#include "pointer.icl"

      ive_ptr function press_drag(whch,i1,i2,j1,j2,ndims,dims,stag,
     >              min,max,miss,dunits,disunits,dim_names)
*
*  press_drag computes the horizontal stress exerted by the topography 
*    on the fluid.  The stress is calculated over the domain defined by
*    i1*dx to i2*dx along the x-axis, and j1*dy to j2*dy along the 
*    y-axis.
*
      integer ndims,dims(4),i1,i2,j1,j2
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),whch

      integer k,iu,iv,istrt,istop,jstrt,jstop,strbeg,strend,twod
      real dx,dy,area,del,dom_min(4),dom_max(4)
      character*1 comp
      ive_ptr p,p0,zs,getvar,getmem
      logical err,new

      press_drag = 0

      call getrarr('plmin',dom_min,4,err)
      if (err) then
        write(*,*) 'Error:  could not find physical space domain ',
     >                            'minimums.'
        return
      end if
      call getrarr('plmax',dom_max,4,err)
      if (err) then
        write(*,*) 'Error:  could not find physical space domain ',
     >                              'maximums.'
        return
      end if

      if (dom_min(2).eq.dom_max(2)) then
        twod = 1
      else 
        twod = 0
      end if

      iu = 0
      iv = 0

      comp = whch(strbeg(whch):strbeg(whch)) 

      if ((comp.eq.'Y').and.(twod.eq.1)) then
        write(*,*) 'Error:  cannot calculate meridional pressure ',
     >                         'drag for 2D run.'
        return
      end if

      if (comp.eq.'X') then
        zs = getvar('zbot_u',ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,err)
        iu = 1
      else if (comp.eq.'Y') then
        zs = getvar('zbot_v',ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,err)
        iv = 1
      else
        write(*,*) 'Error:  cannot calculate pressure drag for ',
     >                         'component ',comp(1:1)
        return
      end if

      p0 = getvar('p0',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      p = getvar('p',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)

      dx = (max(1)-min(1))/(dims(1)-1)
      if (twod.eq.1) then
        j1 = 0
        j2 = 1
        dy = 1.
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if
 
      area = dims(1)*dx*dims(2)*dy

      istrt = i1+1
      istop = i2
      jstrt = j1+1
      jstop = j2

      press_drag = getmem(dims(4))

      if (comp.eq.'X') then
        del = dy
      else 
        del = dx
      end if 

      call drag_calc(%val(press_drag),%val(p),%val(p0),%val(zs),del,
     >                area,istrt,istop,jstrt,jstop,dims(1),dims(2),
     >                dims(3),dims(4),iu,iv) 

      do k=1,3
        dims(k) = 1
        stag(k) = 0.
        min(k) = 0.
        max(k) = 0.
      enddo

      dunits = 'Pa'
      disunits = 'Pa'

      return
      end

************************************************************************

      subroutine drag_calc(drag,p,p0,zs,del,area,istrt,istop,jstrt,
     >                               jstop,nx,ny,nz,nt,iu,iv)

      integer nx,ny,nz,nt,istrt,istop,jstrt,jstop,iu,iv
      real p(nx,ny,nz,nt),p0(nx,ny,nz),zs(nx+iu,ny+iv),drag(nt),
     >                           del,area

      integer i,j,n

      do n=1,nt

        drag(n) = 0

        do i=istrt,istop
        do j=jstrt,jstop
          drag(n) = drag(n) - 100.*(p(i,j,1,n)+p0(i,j,1))
     >                      *(zs(i+iu,j+iv)-zs(i,j))*del
        enddo
        enddo

        drag(n) = drag(n)/area

      enddo

      return
      end
  

