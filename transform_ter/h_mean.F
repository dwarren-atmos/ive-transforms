
#include "pointer.icl"

      ive_ptr function h_mean(whch,varnm,numvars,i1,i2,j1,j2,dpth,
     >      ndims,dims,stag,min,max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer numvars,i1,i2,j1,j2,dpth,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) whch,varnm(numvars),dunits,disunits,dim_names(4)

      integer strbeg,strend,iu,iv,iw,k,nt,n
      integer istrt,istop,jstrt,jstop
      real r_val
      character*(80) vunits,vdisunits,tmp,tmp2
      logical new,err,requ
      ive_ptr getmem,v_gz,var,sgz_tmp,wgz_tmp,gz

      h_mean = 0

      if ((whch(1:1).ne.'Z').and.(whch(1:1).ne.'N')) then
        write(*,*) 'Error:  not sure what to do in h_mean.  whch is ',
     >                        whch(1:1)
        return
      end if

      wgz_tmp = getmem(nz+1)
      if (wgz_tmp.eq.0) then
        write(*,*) 'Error:  unable to allocate space for wgz_tmp.'
        return
      end if
      do k=1,nz+1
        call s_val(%val(wgz_tmp),k,r_val(%val(wgz),k))
      enddo

      sgz_tmp = getmem(nz)
      if (sgz_tmp.eq.0) then
        write(*,*) 'Error:  unable to allocate space for sgz_tmp.'
        call freemem(wgz_tmp)
        return
      end if
      do k=1,nz
        call s_val(%val(sgz_tmp),k,r_val(%val(sgz),k))
      enddo

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if

      istrt = i1+1
      istop = i2
      jstrt = j1+1
      jstop = j2

 545  FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 548  FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)

      tmp = varnm(1)
      do n=2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >     varnm(n)(strbeg(varnm(n)):strend(varnm(n)))
      enddo

      write(*,*)
      if (whch.eq.'N') then
        write(*,*) 'Calculating the horizontal mean of field ',
     >             tmp(strbeg(tmp):strend(tmp))
      else
        write(*,*) 'Calculating vertical derivative of horizontal ',
     >                'mean of field ',
     >             tmp(strbeg(tmp):strend(tmp))
      end if
      write(*,*) 'over the domain defined by:  '
      write(*,545) plmin(1)+i1*dx, plmin(1)+i2*dx
      write(*,548) plmin(2)+j1*dy, plmin(2)+j2*dy
      write(*,*) 

      if ((istrt.gt.istop).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain ',
     >                'of integration.'
        call freemem(wgz_tmp)
        call freemem(sgz_tmp)
        return
      else if ((istrt.lt.1).or.(istop.gt.nx).or.(jstrt.lt.1)
     >           .or.(jstop.gt.ny)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >               'boundaries.'
        call freemem(wgz_tmp)
        call freemem(sgz_tmp)
        return
      end if

      var = getvar(varnm(1),ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot access variable ',
     >          varnm(1)(strbeg(varnm(1)):strend(varnm(1)))
        call freemem(sgz_tmp)
        call freemem(wgz_tmp)
        return
      end if
      vunits = dunits
      vdisunits = disunits
      nt = dims(4)

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

      if (whch.eq.'N') then
        v_gz = getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (v_gz.eq.0) then
          write(*,*) 'Error:  unable to allocate space for v_gz.'
          call freemem(sgz_tmp)
          call freemem(wgz_tmp)
          return
        end if
        do k=1,dims(1)*dims(2)*dims(3)*dims(4)
          call s_val(%val(v_gz),k,1.)
        enddo
      else
        v_gz = getmem(dims(1)*dims(2)*(dims(3)+1)*dims(4))
        if (v_gz.eq.0) then
          write(*,*) 'Error:  unable to allocate space for v_gz.'
          call freemem(sgz_tmp)
          call freemem(wgz_tmp)
          return
        end if
        do k=1,dims(1)*dims(2)*(dims(3)+1)*dims(4)
          call s_val(%val(v_gz),k,1.)
        enddo
      end if

      if (iw.eq.1) then
        gz = wgz
      else
        gz = sgz
      end if

      miss = 8675309

      if (whch.eq.'N') then
        call v_at_sgz(%val(v_gz),%val(var),%val(zbot_p),%val(gz),
     >        %val(sgz_tmp),ztop,miss,dims(1),dims(2),dims(3),
     >        dims(4),iu,iv,iw)
      else
        call v_2_wgz(%val(v_gz),%val(var),%val(zbot_p),%val(gz),
     >        %val(wgz_tmp),ztop,miss,dims(1),dims(2),dims(3),
     >        dims(4),iu,iv,iw)
      end if

      k = 2

 547  if (k.le.numvars) then
   
        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  cannot access variable ',
     >            varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          call freemem(v_gz)
          call freemem(wgz_tmp)
          call freemem(sgz_tmp)
          return
        else if (dims(4).ne.nt) then
          write(*,*) 'Error:  variables ',
     >      varnm(1)(strbeg(varnm(1)):strend(varnm(1))),' and ',
     >      varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          write(*,*) 'have different time dimensions.'
          call freemem(v_gz)
          call freemem(wgz_tmp)
          call freemem(sgz_tmp)
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

        if (iw.eq.1) then
          gz = wgz
        else
          gz = sgz
        end if

        miss = 8675309

        if (whch.eq.'N') then
          call v_at_sgz(%val(v_gz),%val(var),%val(zbot_p),%val(gz),
     >        %val(sgz_tmp),ztop,miss,dims(1),dims(2),dims(3),
     >        dims(4),iu,iv,iw)
        else
          call v_2_wgz(%val(v_gz),%val(var),%val(zbot_p),%val(gz),
     >        %val(wgz_tmp),ztop,miss,dims(1),dims(2),dims(3),
     >        dims(4),iu,iv,iw)
        end if

        tmp = vunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >               dunits(strbeg(dunits):strend(dunits))
        tmp = vdisunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vdisunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >            disunits(strbeg(disunits):strend(disunits))

        k = k+1

        goto 547

      end if

      h_mean = getmem(dims(3)*dims(4))
      if (h_mean.eq.0) then
        write(*,*) 'Error:  unable to allocate space for h_mean.'
        call freemem(v_gz)
        call freemem(wgz_tmp)
        call freemem(sgz_tmp)
        return
      end if

      if (whch.eq.'N') then
        call h_mean_calc(%val(h_mean),%val(v_gz),%val(zbot_p),
     >          %val(zbot_u),%val(zbot_v),%val(sgz),%val(wgz),miss,
     >          dx,dy,istrt,istop,jstrt,jstop,dpth,dims(1),dims(2),
     >          dims(3),dims(4),twod)
      else
        call dh_mean_dz_calc(%val(h_mean),%val(v_gz),%val(zbot_p),
     >       %val(zbot_u),%val(zbot_v),%val(wgz),miss,dx,dy,istrt,
     >       istop,jstrt,jstop,dpth,dims(1),dims(2),dims(3),dims(4),
     >       twod)
        vunits = vunits(strbeg(vunits):strend(vunits))//'/m'
        vdisunits = vdisunits(strbeg(vdisunits):strend(vdisunits))
     >                           //'/m'
      end if
 
      dims(1) = 1
      dims(2) = 1 
      stag(1) = 0.
      stag(2) = 0.
      stag(3) = 0.5
      min(1) = 0.
      min(2) = 0.
      min(3) = r_val(%val(sgz),1)
      max(1) = 0.
      max(2) = 0.
      max(3) = r_val(%val(sgz),dims(3))
      dim_names(1) = 'one'
      dim_names(2) = 'one'
      dim_names(3) = 'nz'

      call units_smpfy(vunits,dunits,err)
      if (err) dunits = ' '
      call units_smpfy(vdisunits,disunits,err)
      if (err) disunits =  ' '

      call freemem(v_gz)
      call freemem(wgz_tmp)
      call freemem(sgz_tmp)

      return
      end

*--------------------------------------------------------------------------

      subroutine h_mean_calc(h_mean,v_sgz,zs_p,zs_u,zs_v,sgz,wgz,miss,
     >                dx,dy,istrt,istop,jstrt,jstop,dpth,nx,ny,nz,nt,
     >                twod)

      integer istrt,istop,jstrt,jstop,dpth,nx,ny,nz,nt,twod
      real h_mean(nz,nt),v_sgz(nx,ny,nz,nt),zs_p(nx,ny)
      real zs_u(nx+1,ny),zs_v(nx,ny+1-twod),sgz(nz),wgz(nz+1)
      real miss,dx,dy

      integer k,kk,n,NZMAX,NTMAX,k1,k2,delk,num
      parameter (NZMAX=151,NTMAX=51)
      real area,area_dum,sum(NZMAX,NTMAX),delz

      area = (istop-istrt+1)*dx*(jstop-jstrt+1)*dy

      do n=1,nt
      do k=1,nz

        if (twod.ne.1) then
          call area_int_z_calc2(sum(k,n),area_dum,v_sgz(1,1,k,n),zs_p,
     >         zs_u,zs_v,sgz(k),dx,dy,miss,istrt,istop,jstrt,jstop,
     >         nx,ny)
        else
          call area_int_z_calc1(sum(k,n),area_dum,v_sgz(1,1,k,n),zs_p,
     >         zs_u,sgz(k),dx,miss,istrt,istop,nx)
        end if

        if (area_dum.ne.0) then
          sum(k,n) = sum(k,n)/area
        else
          sum(k,n) = miss
        end if

      enddo
      enddo

      delk = (dpth-1)/2

      do n=1,nt
      do k=1,nz

        h_mean(k,n) = 0. 
        k1 = max(k-delk,1)
        k2 = min(k+delk,nz)
        delz = 0.
      
        do kk=k1,k2
          if (sum(kk,n).ne.miss) then
            h_mean(k,n) = h_mean(k,n) + sum(kk,n)*(wgz(kk+1)-wgz(kk))
            delz = delz + wgz(kk+1)-wgz(kk)
          end if
        enddo

        if (delz.ne.0) then
          h_mean(k,n) = h_mean(k,n)/delz
        else
          h_mean(k,n) = miss
        end if

      enddo 
      enddo

      return
      end

*----------------------------------------------------------------------------

      subroutine dh_mean_dz_calc(h_mean,v_wgz,zs_p,zs_u,zs_v,wgz,
     >         miss,dx,dy,istrt,istop,jstrt,jstop,dpth,nx,ny,nz,nt,
     >         twod)

      integer istrt,istop,jstrt,jstop,dpth,nx,ny,nz,nt,twod
      real h_mean(nz,nt),v_wgz(nx,ny,nz+1,nt),zs_p(nx,ny)
      real zs_u(nx+1,ny),zs_v(nx,ny+1-twod),wgz(nz+1),miss,dx,dy

      integer k,n,NZMAX,NTMAX,delk,k1,k2
      parameter (NZMAX=151,NTMAX=51)
      real sum(NZMAX,NTMAX),area,area_dum

      area = (istop-istrt+1)*dx*(jstop-jstrt+1)*dy

      do n=1,nt
      do k=1,nz+1

        if (twod.ne.1) then
          call area_int_z_calc2(sum(k,n),area_dum,v_wgz(1,1,k,n),zs_p,
     >         zs_u,zs_v,wgz(k),dx,dy,miss,istrt,istop,jstrt,jstop,
     >         nx,ny)
        else
          call area_int_z_calc1(sum(k,n),area_dum,v_wgz(1,1,k,n),zs_p,
     >         zs_u,wgz(k),dx,miss,istrt,istop,nx)
        end if

        if (area_dum.ne.0) then
          sum(k,n) = sum(k,n)/area
        else
          sum(k,n) = miss
        end if
     
      enddo
      enddo

      delk = (dpth-1)/2

      do n=1,nt
      do k=1,nz

        k1 = max(k-delk,1)
        k2 = min(k+delk+1,nz+1)
     
        if ((sum(k2,n).ne.miss).and.(sum(k1,n).ne.miss)) then
          h_mean(k,n) = (sum(k2,n)-sum(k1,n))/(wgz(k2)-wgz(k1))
        else
          h_mean(k,n) = miss
        end if

      enddo
      enddo

      return
      end
