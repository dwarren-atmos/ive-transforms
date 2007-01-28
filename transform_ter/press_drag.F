
#include "pointer.icl"

      ive_ptr function press_drag(whch,i1,i2,j1,j2,k1,k2,ndims,dims,
     >         stag,min,max,miss,dunits,disunits,dim_names)

#include "default.icl"
*
*  press_drag computes the horizontal stress exerted by the topography 
*    on the fluid.  The stress is calculated over the domain defined by
*    i1*dx to i2*dx along the x-axis, and j1*dy to j2*dy along the 
*    y-axis.
*
      integer ndims,dims(4),i1,i2,j1,j2,k1,k2
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),whch

      integer k,iu,iv,istrt,istop,jstrt,jstop,strbeg,strend
      real del,r_val,z1,z2
      character*1 comp
      character*(80) tmp
      ive_ptr p,zs,getmem
      logical err,new

      press_drag = 0

      iu = 0
      iv = 0

      comp = whch(strbeg(whch):strbeg(whch)) 

      if ((comp.eq.'Y').and.(twod.eq.1)) then
        write(*,*) 'Error:  cannot calculate meridional pressure ',
     >                         'drag for 2D run.'
        return
      end if

      if (comp.eq.'X') then
        zs = zbot_u
        iu = 1
      else if (comp.eq.'Y') then
        zs = zbot_v
        iv = 1
      else
        write(*,*) 'Error:  cannot calculate pressure drag for ',
     >                         'component ',comp(1:1)
        return
      end if

      p = getvar('p',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                         dim_names,new)
      if (p.eq.0) then
        write(*,*) 'Error:  unable to access pressure.'
        return
      end if

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if
 
      istrt = i1 + 1
      istop = i2 
      jstrt = j1 + 1
      jstop = j2 
      z1 = k1*dz
      z2 = k2*dz

      if ((z1.ge.z2).or.(istrt.gt.istop).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                        'integration.'
        return
      else if ((istrt.lt.1).or.(istop.gt.dims(1)).or.(jstrt.lt.1)
     >              .or.(jstop.gt.dims(2)).or.(z1.lt.0)) then
        write(*,*) 'Error:  specified integration domain exceeds ',
     >                   'model domain.'
        return
      end if

 79   FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 80   FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)
 81   FORMAT(10x,'z1 = ',F10.2,2x,'z2 = ',F10.2)

      write(*,*)
      write(*,*) 'Calculating pressure drag for dimension ',comp(1:1)
      write(*,*) 'over the domain defined by: '
      write(*,79) plmin(1)+i1*dx, plmin(1)+i2*dx
      write(*,80) plmin(2)+j1*dy, plmin(2)+j2*dy
      write(*,81) z1,z2
      write(*,*)

      press_drag = getmem(dims(4))

      if (comp.eq.'X') then
        del = dy
      else 
        del = dx
      end if 

      call drag_calc(%val(press_drag),%val(p),%val(zs),z1,z2,del,
     >                istrt,istop,jstrt,jstop,dims(1),dims(2),
     >                dims(3),dims(4),iu,iv) 

      do k=1,3
        dims(k) = 1
        stag(k) = 0.
        min(k) = 0.
        max(k) = 0.
        dim_names(k) = 'one'
      enddo

      tmp = dunits(strbeg(dunits):strend(dunits))//'*m*m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '

      tmp = disunits(strbeg(disunits):strend(disunits))//'*m*m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '

      return
      end

*------------------------------------------------------------------------

      subroutine drag_calc(drag,p,zs,z1,z2,del,istrt,istop,jstrt,
     >                               jstop,nx,ny,nz,nt,iu,iv)

      integer nx,ny,nz,nt,istrt,istop,jstrt,jstop,iu,iv
      real p(nx,ny,nz,nt),zs(nx+iu,ny+iv),drag(nt),del,z1,z2 

      integer i,j,n
      real delz,zt,zb

      do n=1,nt

        drag(n) = 0

        do i=istrt,istop
        do j=jstrt,jstop

          zt = max(zs(i,j),zs(i+iu,j+iv))
          zb = min(zs(i,j),zs(i+iu,j+iv))
          
          if ((zt.gt.z1).and.(zb.lt.z2)) then

            delz = min(zt,z2)-max(zb,z1)
            if (zs(i+iu,j+iv).lt.zs(i,j)) delz = -delz
            drag(n) = drag(n) - 0.5*(3.*p(i,j,1,n)-p(i,j,2,n))
     >                    *delz*del
   
          end if   

        enddo
        enddo

      enddo

      return
      end

*---------------------------------------------------------------------------

      ive_ptr function drag_density(k1,k2,ndims,dims,stag,min,max,
     >                      miss,dunits,disunits,dim_names)

#include "default.icl"

      integer k1,k2,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer j,strbeg,strend
      character*(1) comp
      character*(80) tmp
      logical err
      ive_ptr press_drag,drag,getmem

      comp = 'X'

      drag_density = getmem((ny+1)*ntime)
      if (drag_density.eq.0) then
        write(6,*) 'Error:  unable to allocate space for ',
     >                          'drag_density.'
        return
      end if

      drag = getmem(ntime)
      if (drag.eq.0) then
        write(6,*) 'Error:  unable to allocate space for drag.'
        call freemem(drag_density)
        return
      end if

      do j=2,ny
      
        drag = press_drag(comp,1,nx,j-2,j,k1,k2,ndims,dims,stag,min,
     >                         max,miss,dunits,disunits,dim_names)
        if (drag.eq.0) then
          call freemem(drag_density)
          call freemem(drag)
          return
        end if

        call drag_density_calc(%val(drag_density),%val(drag),dy,miss,
     >                                     j,ny,ntime)

      enddo

      miss = -777

      call drag_density_calc(%val(drag_density),%val(drag),dy,miss,
     >                                      1,ny,ntime)
      call drag_density_calc(%val(drag_density),%val(drag),dy,miss,
     >                                     ny+1,ny,ntime)

      dims(2) = ny+1
      stag(2) = 0.
      min(2) = plmin(2)
      max(2) = plmax(2)
      dim_names(2) = 'nyp1'
      
      tmp = dunits//'/m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '

      tmp = disunits//'/m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '

      call freemem(drag)

      return
      end

      subroutine drag_density_calc(drag_density,drag,dy,miss, 
     >                                  j,ny,ntime)

      integer j,ny,ntime
      real drag_density(ny+1,ntime),drag(ntime),dy,miss

      integer n

      if ((j.eq.1).or.(j.eq.ny+1)) then
        do n=1,ntime
          drag_density(j,n) = miss
        enddo
      else
        do n=1,ntime
          drag_density(j,n) = drag(n)/(2.*dy)
        enddo
      end if

      return
      end

*------------------------------------------------------------------------
      
      ive_ptr function press_drag_spec(ndims,dims,stag,min,max,miss,
     >                              dunits,disunits,dim_names)

#include "default.icl"
 
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer strbeg,strend,nt,k
      character*(80) tmp
      ive_ptr p,getmem
      logical err,new

      nt = ntime

      press_drag_spec = 0

      p = getvar('p',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                         dim_names,new)
      if (p.eq.0) then
        write(*,*) 'Error:  unable to access pressure.'
        return
      end if

      write(*,*)
      write(*,*) 'Calculating pressure drag.'
      write(*,*)

      press_drag_spec = getmem(dims(4))

      call drag_spec_calc(%val(press_drag_spec),%val(p),%val(zbot_u),
     >              %val(zbot_p),%val(sgz),ztop,dx,dy,nx,ny,nz,nt)

      do k=1,3
        dims(k) = 1
        stag(k) = 0.
        min(k) = 0.
        max(k) = 0.
        dim_names(k) = 'one'
      enddo

      tmp = dunits(strbeg(dunits):strend(dunits))//'*m*m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '

      tmp = disunits(strbeg(disunits):strend(disunits))//'*m*m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '

      return
      end

*------------------------------------------------------------------------

      subroutine drag_spec_calc(drag,p,zs_u,zs_p,sgz,ztop,dx,dy,
     >                               nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real drag(nt),p(nx,ny,nz,nt),zs_u(nx+1,ny),zs_p(nx,ny)
      real sgz(nz),ztop,dx,dy

      integer i,j,n
      real z0,z1,z2,z3,z4,ps,dhdx,rdx

      rdx = 1./dx

      do n=1,nt

        drag(n) = 0

        do i=2,nx-1
        do j=1,ny

            z3 = zs_p(i,j) + (ztop-zs_p(i,j))/ztop*sgz(3)
            z2 = zs_p(i,j) + (ztop-zs_p(i,j))/ztop*sgz(2)
            z1 = zs_p(i,j) + (ztop-zs_p(i,j))/ztop*sgz(1)
            z0 = zs_p(i,j)

            ps = (z0-z1)*(z0-z2)
     >                    / (z3-z1)/(z3-z2)*p(i,j,3,n)
     >         + (z0-z1)*(z0-z3)
     >                    / (z2-z1)/(z2-z3)*p(i,j,2,n)
     >         + (z0-z2)*(z0-z3)
     >                    / (z1-z2)/(z1-z3)*p(i,j,1,n)

            dhdx = 4./3.*rdx*(zs_u(i+1,j)-zs_u(i,j))
     >                  - 1./6.*rdx*(zs_p(i+1,j)-zs_p(i-1,j))

            drag(n) = drag(n) - ps*dhdx*dx*dy
   
        enddo
        enddo

      enddo

      return
      end

      
   
      
      
  

