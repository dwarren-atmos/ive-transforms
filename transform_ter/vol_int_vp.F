
#include "pointer.icl"

      ive_ptr function vol_int_vp(varnm,numvars,i1,i2,j1,j2,dpth,
     >         ndims,dims,stag,min,max,miss,dunits,disunits,
     >         dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,dpth,numvars
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dunits,disunits,dim_names(4)

      integer istrt,istop,jstrt,jstop,iu,iv,iw,k,strbeg,strend
      integer ijk,nt
      real zstrt,zstop
      character*(80) tmp,tmp2
      logical err,new,requ
      ive_ptr getmem,var,var_t,wrksp

      vol_int_vp = 0

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if

      istrt = i1 + 1
      istop = i2 
      jstrt = j1 + 1
      jstop = j2 

 101  FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 104  FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)
 105  FORMAT(10x,'depth = ',I3)
 108  FORMAT(10x,'z1 = ',F10.2,2x,'z2 = ',A3)

      tmp = varnm(1)
      do k=2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
      enddo
    
      write(*,*)
      write(*,*) 'Calculating vertical profile of volume integral of ',
     >              'field ',tmp(strbeg(tmp):strend(tmp))
      write(*,*) 'over the domain defined by:  '
      write(*,101) plmin(1)+i1*dx, plmin(1)+i2*dx
      write(*,104) plmin(2)+j1*dy, plmin(2)+j2*dy
      if (dpth.eq.0) then
        write(*,108) 0.,'z'
      else
        write(*,105) dpth
      end if
      write(*,*)

      if ((istrt.lt.1).or.(istop.gt.nx).or.(jstrt.lt.1).or.
     >               (jstop.gt.ny).or.(dpth.gt.nz)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries.' 
        return
      else if ((istrt.gt.istop).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                        'integration.'
        return
      end if 

      var = getvar(varnm(1),ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot access field ',
     >              varnm(1)(strbeg(varnm(1)):strend(varnm(1)))
        return
      end if

      iu = 0
      iv = 0
      iw = 0
      if (requ(plmin(1),min(1))) then
        iu = 1
      else if (requ(plmin(2),min(2)).and.(twod.ne.1)) then
        iv = 1
      else if (requ(plmin(3),min(3))) then
        iw = 1
      end if

      dims(1) = dims(1) - iu
      dims(2) = dims(2) - iv
      dims(3) = dims(3) - iw

      var_t = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (var_t.eq.0) then
        write(*,*) 'Error:  could not allocate space for var_t.'
        return
      end if

      do ijk=1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(var_t),ijk,1.)
      enddo

      call var_2_tpts(%val(var_t),%val(var),dims(1),dims(2),dims(3),
     >                        dims(4),iu,iv,iw)

      nt = dims(4)

      k = 2

 873  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  cannot access field ',
     >            varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          call freemem(var_t)
          return
        else if (nt.ne.dims(4)) then
          write(*,*) 'Error:  fields ',
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >                        ' and ',
     >        varnm(1)(strbeg(varnm(1)):strend(varnm(1))),
     >           ' have different time dimensions.'
          call freemem(var_t)
          return
        end if
      
        iu = 0
        iv = 0
        iw = 0
        if (requ(plmin(1),min(1))) then
          iu = 1
        else if (requ(plmin(2),min(2)).and.(twod.ne.1)) then
          iv = 1
        else if (requ(plmin(3),min(3))) then
          iw = 1
        end if

        dims(1) = dims(1) - iu
        dims(2) = dims(2) - iv
        dims(3) = dims(3) - iw

        call var_2_tpts(%val(var_t),%val(var),dims(1),dims(2),dims(3),
     >                        dims(4),iu,iv,iw)

        k = k+1

        goto 873
   
      end if

      if (mod(dpth,2).eq.0) then
        iw = 1
      else 
        iw = 0
      end if

      vol_int_vp = getmem((dims(3)+iw)*dims(4))
      if (vol_int_vp.eq.0) then
        write(*,*) 'Error:  could not allocate space for vol_int.'
        call freemem(var_t)
        return
      end if

      wrksp = getmem(nt)
      if (wrksp.eq.0) then
        write(*,*) 'Error:  could not allocate space for wrksp.'
        call freemem(var_t)
        call freemem(vol_int_vp)
        return
      end if

      call vol_int_vp_calc(%val(vol_int_vp),%val(var_t),%val(zbot_p),
     >         %val(sgz),%val(wgz),%val(wrksp),ztop,dx,dy,istrt,
     >         istop,jstrt,jstop,dims(1),dims(2),dims(3),dims(4),
     >         dpth,iw,err)

      do k=1,2
        dims(k) = 1
        min(k) = 0.
        max(k) = 0.
        stag(k) = 0.
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
 
      tmp = dunits(strbeg(dunits):strend(dunits))//'*m*m*m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '

      tmp = disunits(strbeg(disunits):strend(disunits))//'*m*m*m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '

      call freemem(var_t)
      call freemem(wrksp)

      return
      end

*-------------------------------------------------------------------------

      subroutine vol_int_vp_calc(vol_int,v_t,zs_p,sgz,wgz,wrk,ztop,
     >                  dx,dy,istrt,istop,jstrt,jstop,nx,ny,nz,nt,
     >                  dpth,iw,err)

      integer nx,ny,nz,nt,istrt,istop,jstrt,jstop,dpth,iw
      real vol_int(nz+iw,nt),v_t(nx,ny,nz,nt),zs_p(nx,ny)
      real sgz(nz),wgz(nz+1),wrk(nt),ztop,dx,dy
      logical err

      integer n,k,kbot,ktop,del
      real z1,z2

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

        z1 = wgz(kbot)
        z2 = wgz(ktop)

        call vol_int_calc(wrk,v_t,zs_p,sgz,wgz,ztop,dx,dy,z1,z2,
     >              istrt,istop,jstrt,jstop,nx,ny,nz,nt,err)
        if (err) then
          write(6,*) 'Error:  problem calculating volume integral.'
          write(6,*) 'Stopping run.'
          stop
        end if
        
        do n=1,nt
          vol_int(k,n) = wrk(n)
        enddo

      enddo

      return
      end

*-----------------------------------------------------------------------

      ive_ptr function vol_ave_vp(varnm,numvars,i1,i2,j1,j2,dpth,
     >         ndims,dims,stag,min,max,miss,dunits,disunits,
     >         dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,dpth,numvars
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dunits,disunits,dim_names(4)

      integer strbeg,strend,n,iw,iarg(5),nt
      real area,rarg(1)
      character*(80) intstr,mkfield,carg(10),tmp,ustr
      logical new,err
      ive_ptr vol_int,getmem

      vol_ave_vp = 0
  
      nt = ntime

      do n=1,numvars
        carg(n)(1:len(carg(n))) = ' '
        carg(n) = varnm(n)(strbeg(varnm(n)):strend(varnm(n)))
      enddo

      iarg(1) = i1
      iarg(2) = i2
      iarg(3) = j1
      iarg(4) = j2
      iarg(5) = dpth

      intstr = mkfield('VOL_INT_VP',carg,iarg,rarg,numvars,5,0,0)

      vol_int = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                                disunits,dim_names,new)
      if (vol_int.eq.0) then
        write(6,*) 'Error:  could not compute integrals.'
        return
      end if

      area = real(i2-i1)*dx*real(j2-j1)*dy

      if (mod(dpth,2).eq.0) then
        iw = 1
      else
        iw = 0
      end if

      vol_ave_vp = getmem((nz+iw)*nt)
      if (vol_ave_vp.eq.0) then
        write(6,*) 'Error:  unable to allocate space for vol_ave_vp.'
        return
      end if

      call vol_ave_vp_calc(%val(vol_ave_vp),%val(vol_int),%val(wgz),
     >                         area,nz,nt,iw,dpth)

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

*---------------------------------------------------------------------------

      subroutine vol_ave_vp_calc(vol_ave,vol_int,wgz,area,nz,nt,
     >                                      iw,dpth)

      integer nz,nt,iw,dpth
      real vol_ave(nz+iw,nt),vol_int(nz+iw,nt),wgz(nz+1),area

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
          vol_ave(k,n) = vol_int(k,n)/vol
        enddo

      enddo

      return
      end 
      

     
      
