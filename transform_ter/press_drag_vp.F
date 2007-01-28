
#include "pointer.icl"

      ive_ptr function drag_vp(whch,i1,i2,j1,j2,dpth,ndims,dims,
     >         stag,min,max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,dpth
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),whch

      integer k,iu,iv,iw,istrt,istop,jstrt,jstop,strbeg,strend
      real del
      character*1 comp
      character*(80) tmp
      ive_ptr p,zs,getmem,wrksp
      logical err,new

      drag_vp = 0

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

      if ((istrt.gt.istop).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                        'integration.'
        return
      else if ((istrt.lt.1).or.(istop.gt.dims(1)).or.(jstrt.lt.1)
     >              .or.(jstop.gt.dims(2)).or.(dpth.gt.nz)) then
        write(*,*) 'Error:  specified integration domain exceeds ',
     >                   'model domain.'
        return
      end if

 79   FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 80   FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)
 81   FORMAT(10x,'depth = ',I3)
 91   FORMAT(10x,'z1 = ',F10.2,2x,'z2 = ',A3)

      write(*,*)
      write(*,*) 'Calculating pressure drag for dimension ',comp(1:1)
      write(*,*) 'over the domain defined by: '
      write(*,79) plmin(1)+i1*dx, plmin(1)+i2*dx
      write(*,80) plmin(2)+j1*dy, plmin(2)+j2*dy
      if (dpth.eq.0) then
        write(*,91) 0.,'z'
      else
        write(*,81) dpth
      end if
      write(*,*)

      if (mod(dpth,2).eq.0) then
        iw = 1
      else 
        iw = 0
      end if

      drag_vp = getmem((dims(3)+iw)*dims(4))
      if (drag_vp.eq.0) then
        write(6,*) 'Error:  could not allocate space for drag_vp.'
        return
      end if

      wrksp = getmem(dims(4))
      if (wrksp.eq.0) then
        write(6,*) 'Error:  could not allocate space for wrksp.'
        return
      end if

      if (comp.eq.'X') then
        del = dy
      else 
        del = dx
      end if 

      call drag_vp_calc(%val(drag_vp),%val(p),%val(zs),%val(wgz),
     >                %val(wrksp),del,istrt,istop,jstrt,jstop,dims(1),
     >                dims(2),dims(3),dims(4),iu,iv,iw,dpth) 

      do k=1,2
        dims(k) = 1
        stag(k) = 0.
        min(k) = 0.
        max(k) = 0.
        dim_names(k) = 'one'
      enddo

      if (iw.eq.1) then
        dims(3) = dims(3) + 1
        stag(3) = 0.
        min(3) = wgz_min
        max(3) = wgz_max
        dim_names(3) = 'nzp1'
      else
        stag(3) = 0.5
        min(3) = sgz_min
        max(3) = sgz_max
        dim_names(3) = 'nz'
      end if

      tmp = dunits(strbeg(dunits):strend(dunits))//'*m*m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '

      tmp = disunits(strbeg(disunits):strend(disunits))//'*m*m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '

      call freemem(wrksp)

      return
      end

*------------------------------------------------------------------------

      subroutine drag_vp_calc(drag,p,zs,wgz,wrk,del,istrt,istop,jstrt,
     >                      jstop,nx,ny,nz,nt,iu,iv,iw,dpth)

      integer nx,ny,nz,nt,istrt,istop,jstrt,jstop,iu,iv,iw,dpth
      real drag(nz+iw,nt),p(nx,ny,nz,nt),zs(nx+iu,ny+iv),wgz(nz+1)
      real wrk(nt),del

      integer n,k,kbot,ktop,idel
      real z1,z2

      idel = dpth/2

      do k=1,nz+iw

        if (dpth.eq.0) then
          kbot = 1
          ktop = k
        else if (iw.eq.1) then
          kbot = max(1,k-idel)
          ktop = min(nz+1,k+idel)
        else
          kbot = max(1,k-idel)
          ktop = min(nz+1,k+1+idel)
        end if

        z1 = wgz(kbot)
        z2 = wgz(ktop)
   
        call drag_calc(wrk,p,zs,z1,z2,del,istrt,istop,jstrt,jstop,
     >                                 nx,ny,nz,nt,iu,iv)
   
        do n=1,nt
          drag(k,n) = wrk(n)
        enddo

      enddo

      return
      end

*-------------------------------------------------------------------------

      ive_ptr function drag_ave_vp(whch,i1,i2,j1,j2,dpth,ndims,dims,
     >         stag,min,max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,dpth
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),whch

      integer strbeg,strend,iw,iarg(5),nt
      real area,rarg(1)
      character*(80) intstr,mkfield,carg(1),tmp,ustr
      logical err,new
      ive_ptr getmem,drag_int

      drag_ave_vp = 0

      nt = ntime

      carg(1)(1:len(carg(1))) = ' '
      carg(1) = whch(strbeg(whch):strend(whch))

      iarg(1) = i1
      iarg(2) = i2
      iarg(3) = j1
      iarg(4) = j2
      iarg(5) = dpth

      intstr = mkfield('DRAG_VP',carg,iarg,rarg,1,5,0,0)

      drag_int = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                              disunits,dim_names,new)
      if (drag_int.eq.0) then
        write(6,*) 'Error:  could not compute integrals.'
        return
      end if

      area = real(i2-i1)*dx*real(j2-j1)*dy

      if (mod(dpth,2).eq.0) then
        iw = 1
      else
        iw = 0
      end if

      drag_ave_vp = getmem((nz+iw)*nt)
      if (drag_ave_vp.eq.0) then
        write(6,*) 'Error:  unable to allocate space for drag_ave_vp.'
        return
      end if

      call drag_ave_vp_calc(%val(drag_ave_vp),%val(drag_int),%val(wgz),
     >                               area,nz,nt,iw,dpth)

      if (dpth.eq.0) then
        ustr = '/m/m'
      else
        ustr = '/m/m/m'
      end if

      tmp = dunits(strbeg(dunits):strend(dunits))//
     >                 ustr(strbeg(ustr):strend(ustr))
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '

      tmp = disunits(strbeg(disunits):strend(disunits))//
     >                 ustr(strbeg(ustr):strend(ustr))
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '

      return
      end 

*------------------------------------------------------------------------

      subroutine drag_ave_vp_calc(drag_ave,drag,wgz,area,nz,nt,iw,dpth)

      integer nz,nt,iw,dpth
      real drag_ave(nz+iw,nt),drag(nz+iw,nt),wgz(nz+1),area

      integer n,k,kbot,ktop,del
      real vol,delz

      del = dpth/2

      do k=1,nz+iw

        if (dpth.eq.0) then
          kbot = 1
          ktop = k
        else if (iw.eq.1) then
          kbot = max(1,k-del)
          ktop = min(nz+1,k+del)
        else
          kbot = max(1,k-del)
          ktop = min(nz+1,k+1+del)
        end if

        if (dpth.eq.0) then
          delz = 1.
        else
          delz = wgz(ktop)-wgz(kbot)
        end if
        vol = area*delz

        do n=1,nt
          drag_ave(k,n) = drag(k,n)/vol
        enddo

      enddo

      return
      end
 



