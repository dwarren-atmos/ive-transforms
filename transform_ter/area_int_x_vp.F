
#include "pointer.icl"

      ive_ptr function area_int_x_vp(varnm,numvars,j1,j2,dpth,i,ndims,
     >                   dims,stag,min,max,miss,dunits,disunits,
     >                   dim_names)

#include "default.icl"

      integer numvars,j1,j2,dpth,i,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dim_names(4),dunits,disunits

      integer jstrt,jstop,iu,iv,iw,k,strbeg,strend,nt,ifix
      character*(80) tmp,tmp2,vunits,vdisunits
      logical new,err,requ
      ive_ptr getmem,var_x,var,wrksp

      area_int_x_vp = 0
      vunits(1:len(vunits)) = ' '
      vdisunits(1:len(vdisunits)) = ' '
      err = .false.

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if

      jstrt = j1 + 1
      jstop = j2 
      ifix = i + 1

 79   FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)
 80   FORMAT(10x,'at x = ',F10.2)
 90   FORMAT(10x,'depth = ',I3)
 91   FORMAT(10x,'z1 = ',F10.2,2x,'z2 = ',A3)

      tmp = varnm(1)
      do k=2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
      enddo

      write(*,*)
      write(*,*) 'Calculating the vertical profile of the area ',
     >            'integral of field ',tmp(strbeg(tmp):strend(tmp))
      write(*,*) 'over the domain defined by:  '
      write(*,79) plmin(2)+j1*dy, plmin(2)+j2*dy
      if (dpth.eq.0) write(*,91) 0.,'z'
      write(*,80) plmin(1)+(ifix-1)*dx 
      if (dpth.ne.0) write(*,90) dpth
      write(*,*)

      if ((ifix.lt.2).or.(ifix.gt.nx).or.(jstrt.lt.1).or.
     >          (jstop.gt.ny).or.(dpth.gt.nz)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries.'
        return
      else if (jstrt.gt.jstop) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                     'integration.'
        return
      end if

      var = getvar(varnm(1),ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot compute field ',
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

      dims(1) = dims(1)-iu
      dims(2) = dims(2)-iv
      dims(3) = dims(3)-iw
      nt = dims(4)

      var_x = getmem(dims(2)*dims(3)*nt)
      if (var_x.eq.0) then
        write(*,*) 'Error:  could not allocate space for var_x.'
        return
      end if
      
      do k=1,dims(2)*dims(3)*nt
        call s_val(%val(var_x),k,1.)
      enddo

      call v_at_x(%val(var_x),%val(var),ifix,dims(1),dims(2),dims(3),
     >                            dims(4),iu,iv,iw)

      vunits = dunits
      vdisunits = disunits

      k = 2

 348  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  cannot compute field ',
     >             varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          call freemem(var_x)
          return
        else if (dims(4).ne.nt) then
          write(*,*) 'Error:  fields ',
     >        varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                        ' and ',
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >         'have different time dimensions.'
          call freemem(var_x)
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

        dims(1) = dims(1)-iu
        dims(2) = dims(2)-iv
        dims(3) = dims(3)-iw

        call v_at_x(%val(var_x),%val(var),ifix,dims(1),dims(2),
     >                   dims(3),dims(4),iu,iv,iw)

        tmp = vunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >               dunits(strbeg(dunits):strend(dunits))
        tmp = vdisunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vdisunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >        disunits(strbeg(disunits):strend(disunits))

        k = k+1

        goto 348

      end if

      if (mod(dpth,2).eq.0) then
        iw = 1
      else
        iw = 0
      end if    

      area_int_x_vp = getmem((dims(3)+iw)*nt)
      if (area_int_x_vp.eq.0) then
        write(6,*) 'Error:  cannot allocate space for area_int_x_vp.'
        call freemem(var_x)
        return
      end if

      wrksp = getmem(nt)
      if (wrksp.eq.0) then
        write(6,*) 'Error:  cannot allocate space for wrksp.'
        call freemem(var_x)
        call freemem(area_int_x_vp)
        return
      end if

      call area_int_x_vp_calc(%val(area_int_x_vp),%val(var_x),
     >              %val(zbot_u),%val(sgz),%val(wgz),%val(wrksp),
     >              ztop,dx,dy,jstrt,jstop,ifix,dims(1),dims(2),
     >              dims(3),dims(4),dpth,iw,err)

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

      vunits = vunits(strbeg(vunits):strend(vunits))//'*m*m'
      vdisunits = vdisunits(strbeg(vdisunits):strend(vdisunits))
     >                  //'*m*m'

      call units_smpfy(vunits,dunits,err)
      if (err) dunits = ' '
      call units_smpfy(vdisunits,disunits,err)
      if (err) disunits = ' '

      call freemem(var_x)
      call freemem(wrksp)

      return
      end

*-------------------------------------------------------------------------

      subroutine area_int_x_vp_calc(area_int,v_x,zs_u,sgz,wgz,wrk,
     >              ztop,dx,dy,jstrt,jstop,ifix,nx,ny,nz,nt,dpth,
     >                              iw,err)

      integer nx,ny,nz,nt,jstrt,jstop,ifix,dpth,iw
      real area_int(nz+iw,nt),v_x(ny,nz,nt),zs_u(nx+1,ny)
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

        call area_int_x_calc(wrk,v_x,zs_u,sgz,wgz,ztop,dx,dy,z1,z2,
     >                   jstrt,jstop,ifix,nx,ny,nz,nt,err)
        if (err) then
          write(6,*) 'Error:  problem calculating area integral.'
          write(6,*) 'Stopping run.'
          stop
        end if

        do n=1,nt
          area_int(k,n) = wrk(n)
        enddo

      enddo

      return
      end

      

