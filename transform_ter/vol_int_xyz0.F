
#include "pointer.icl"

      ive_ptr function vol_int_xyz0(varnm,numvars,i1,i2,j1,j2,z01,z02,
     >               ndims,dims,stag,min,max,miss,dunits,disunits,
     >               dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,numvars
      real stag(4),min(4),max(4),miss,z01,z02
      character*(*) varnm(numvars),dunits,disunits,dim_names(4)

      integer istrt,istop,jstrt,jstop,iu,iv,iw,k,strbeg,strend,ijk,nt
      real r_val
      character*(80) tmp,tmp2,vunits,vdisunits
      logical err,new,units_err,requ
      ive_ptr getmem,var,var_t,zstrt,zstop

      vol_int_xyz0 = 0

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
 103  FORMAT(10x,'z01 = ',F10.2,2x,'z02 = ',F10.2)

      tmp = varnm(1)
      do k=2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
      enddo
    
      write(6,*)
      write(6,*) 'Calculating volume integral of field ',
     >               tmp(strbeg(tmp):strend(tmp))
      write(6,*) 'over the domain defined by:  '
      write(6,101) plmin(1)+i1*dx, plmin(1)+i2*dx
      write(6,104) plmin(2)+j1*dy, plmin(2)+j2*dy
      write(6,103) z01,z02
      write(6,*)

      if ((z01.lt.0.).or.(z02.gt.wgz_max).or.(istrt.lt.1).or.
     >         (istop.gt.nx).or.(jstrt.lt.1).or.
     >               (jstop.gt.ny)) then
        write(6,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries' 
        return
      else if ((z01.ge.z02).or.(istrt.gt.istop)
     >                   .or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                        'integration.'
        return
      end if 

      var = getvar('zeta',ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot access field zeta.'
        return
      end if

      zstrt = getmem(dims(1)*dims(2)*dims(4))
      if (zstrt.eq.0) then
        write(6,*) 'Error:  cannot allocate space for zstrt.'
        return
      end if

      zstop = getmem(dims(1)*dims(2)*dims(4))
      if (zstop.eq.0) then
        write(6,*) 'Error:  cannot allocate space for zstop.'
        call freemem(zstrt)
        return
      end if

      call zrange_xy(%val(var),%val(zstrt),%val(zstop),%val(zbot_p),
     >                 %val(sgz),z01,z02,ztop,dims(1),dims(2),dims(3),
     >                 dims(4),err)
      if (err) then
        call freemem(zstrt)
        call freemem(zstop)
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
      nt = dims(4)

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
 
      vunits = dunits
      vdisunits = disunits

      k = 2

 873  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
        if (var.eq.0) then
          write(6,*) 'Error:  cannot access field ',
     >            varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          call freemem(zstrt)
          call freemem(zstop)
          call freemem(var_t)
          return
        else if (dims(4).ne.nt) then
          write(6,*) 'Error:  fields ',
     >        varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                       ' and ',
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >         ' have different time dimensions.'
          call freemem(zstrt)
          call freemem(zstop)
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

        tmp = vunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >                 dunits(strbeg(dunits):strend(dunits))
        tmp = vdisunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vdisunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >        disunits(strbeg(disunits):strend(disunits)) 

        k = k+1

        goto 873
   
      end if

      vol_int_xyz0 = getmem(nt)
      if (vol_int_xyz0.eq.0) then
        write(*,*) 'Error:  could not allocate space for ',
     >                             'vol_int_xyz0.'
        call freemem(zstrt)
        call freemem(zstop)
        call freemem(var_t)
        return
      end if

      call vol_int_xyz0_calc(%val(vol_int_xyz0),%val(var_t),
     >         %val(zstrt),%val(zstop),%val(zbot_p),%val(sgz),
     >         %val(wgz),ztop,dx,dy,istrt,istop,jstrt,jstop,
     >         dims(1),dims(2),dims(3),dims(4),err)

      if (err) then
        call freemem(var_t)
        call freemem(vol_int_xyz0)
        vol_int_xyz0 = 0
        return
      end if

      do k=1,3
        dims(k) = 1
        min(k) = 0.
        max(k) = 0.
        stag(k) = 0.
        dim_names(k) = 'one'
      enddo
 
      tmp = dunits(strbeg(dunits):strend(dunits))//'*m*m*m'
      call units_smpfy(tmp,dunits,units_err)
      if (units_err) dunits = tmp

      tmp = disunits(strbeg(disunits):strend(disunits))//'*m*m*m'
      call units_smpfy(tmp,disunits,units_err)
      if (units_err) disunits = tmp

      call freemem(var_t)

      return
      end

*------------------------------------------------------------------------

      subroutine zrange_xy(zeta,zstrt,zstop,zs_p,sgz,z01,z02,ztop,
     >                              nx,ny,nz,nt,err)

      integer nx,ny,nz,nt
      real zeta(nx,ny,nz,nt),zstrt(nx,ny,nt),zstop(nx,ny,nt)
      real zs_p(nx,ny),sgz(nz),z01,z02,ztop
      logical err

      integer i,j,k,n,kptr,kplus
      real za,zb,z0a,z0b,tfac,z0_node,zs

      err = .false.

      do n=1,nt
      do j=1,ny
      do i=1,nx

        zs = zs_p(i,j)
        tfac = (ztop-zs)/ztop

        kptr = 1
        kplus = 0

        if (z01.eq.0.) then
     
          zstrt(i,j,n) = 0.

        else

 532      if ((kplus.eq.0).and.(kptr.le.nz)) then
            z0_node = zs + tfac*sgz(kptr) - zeta(i,j,kptr,n)
            if (z0_node.ge.z01) then
              kplus = kptr
            else
              kptr = kptr+1
            end if
            goto 532
          end if

          if (kplus.eq.0) then
            write(6,*) 'Error:  z01 surface not found.'
            err = .true.
            return
          end if

          if (kplus.eq.1) then
            za = zs + tfac*sgz(1)
            zb = zs + tfac*sgz(2)
            z0a = za - zeta(i,j,1,n)
            z0b = zb - zeta(i,j,2,n)
            zstrt(i,j,n) = za + (zb-za)/(z0b-z0a)*(z01-z0a)
          else
            za = zs + tfac*sgz(kplus-1)
            zb = zs + tfac*sgz(kplus)
            z0a = za - zeta(i,j,kplus-1,n)
            z0b = zb - zeta(i,j,kplus,n)
            zstrt(i,j,n) = zb + (zb-za)/(z0b-z0a)*(z01-z0b)
          end if

          if (zstrt(i,j,n).lt.zs) zstrt(i,j,n) = 0.

        end if

        kplus = 0

 535    if ((kplus.eq.0).and.(kptr.le.nz)) then
          z0_node = zs + tfac*sgz(kptr) - zeta(i,j,kptr,n)
          if (z0_node.ge.z02) then
            kplus = kptr
          else
            kptr = kptr+1
          end if
          goto 535
        end if

        if (kplus.eq.0) then
          write(6,*) 'Error:  z02 surface not found.'
          err = .true.
          return
        end if

        if (kplus.eq.1) then
          za = zs + tfac*sgz(1)
          zb = zs + tfac*sgz(2)
          z0a = za - zeta(i,j,1,n)
          z0b = zb - zeta(i,j,2,n)
          zstop(i,j,n) = za + (zb-za)/(z0b-z0a)*(z02-z0a)
        else
          za = zs + tfac*sgz(kplus-1)
          zb = zs + tfac*sgz(kplus)
          z0a = za - zeta(i,j,kplus-1,n)
          z0b = zb - zeta(i,j,kplus,n)
          zstop(i,j,n) = zb + (zb-za)/(z0b-z0a)*(z02-z0b)
        end if

        if (zstop(i,j,n).le.zs) zstop(i,j,n) = 0.

      enddo
      enddo
      enddo

      return
      end

*-------------------------------------------------------------------------

      subroutine vol_int_xyz0_calc(vol_int_xyz0,var,zstrt,zstop,
     >             zs_p,sgz,wgz,ztop,dx,dy,istrt,istop,jstrt,jstop,
     >             nx,ny,nz,nt,err)

      integer istrt,istop,jstrt,jstop,nx,ny,nz,nt
      real vol_int_xyz0(nt),var(nx,ny,nz,nt),zs_p(nx,ny)
      real zstrt(nx,ny,nt),zstop(nx,ny,nt),sgz(nz),wgz(nz+1)
      real ztop,dx,dy
      logical err

      integer i,j,k,n,ktop,kbot
      real dztop,dzbot,z1,z2,z3,z4,z5,z6,tfac,vbot,vtop,zs,ztmp

      do n=1,nt

      vol_int_xyz0(n) = 0.

      do j=jstrt,jstop
      do i=istrt,istop

        zs = zs_p(i,j)
        tfac = (ztop-zs)/ztop

        if (zstop(i,j,n).gt.zs) then

          if (zstrt(i,j,n).lt.zs) then

            call bot_top(zs,wgz,ztop,zstrt(i,j,n),zstop(i,j,n),
     >                             kbot,ktop,nz)

            if (ktop.eq.1) then
 
              z1 = zs + tfac*sgz(1)
              z2 = zs + tfac*sgz(2)
              dztop = zstop(i,j,n) - zs
              vtop = var(i,j,1,n) + (var(i,j,2,n)-var(i,j,1,n))
     >                 / (z2-z1)*(zs+0.5*dztop-z1)
              vbot = 0.
              dzbot = 0.

            else

              vbot = var(i,j,1,n)
              dzbot = tfac*wgz(2)
              z4 = zs + tfac*sgz(ktop-1)
              z5 = zs + tfac*wgz(ktop)
              z6 = zs + tfac*sgz(ktop)
              dztop = zstop(i,j,n) - z5
              if (zstop(i,j,n).eq.ztop) then
                vtop = var(i,j,ktop,n)
              else
                vtop = var(i,j,ktop,n) 
     >               + (var(i,j,ktop,n)-var(i,j,ktop-1,n))/(z6-z4)
     >                      * (zstop(i,j,n)-0.5*dztop-z6)
              end if

            end if

          else

            call bot_top(zs,wgz,ztop,zstrt(i,j,n),zstop(i,j,n),
     >                                kbot,ktop,nz)

            if ((kbot.eq.1).and.(ktop.eq.1)) then

              z1 = zs + tfac*sgz(1)
              z2 = zs + tfac*sgz(2)
              dztop = zstop(i,j,n) - zstrt(i,j,n)
              vtop = var(i,j,1,n) + (var(i,j,2,n)-var(i,j,1,n))
     >                / (z2-z1)*(zstop(i,j,n)-0.5*dztop-z1)
              vbot = 0.
              dzbot = 0.

            else if ((kbot.eq.nz).and.(ktop.eq.nz)) then

              z1 = zs + tfac*sgz(nz-1)
              z2 = zs + tfac*sgz(nz)
              dztop = zstop(i,j,n) - zstrt(i,j,n)
              vtop = var(i,j,nz,n) + (var(i,j,nz,n)-var(i,j,nz-1,n))
     >                   / (z2-z1)*(zstop(i,j,n)-0.5*dztop-z2)
              vbot = 0.
              dzbot = 0.

            else if (kbot.eq.ktop) then

              dztop = zstop(i,j,n) - zstrt(i,j,n)
              ztmp = zstop(i,j,n) - dztop/2.
              z6 = zs + tfac*sgz(ktop)
              if (ztmp.gt.z6) then
                z1 = z6
                z2 = zs + tfac*sgz(ktop+1)
                vtop = var(i,j,ktop,n)
     >              + (var(i,j,ktop+1,n)-var(i,j,ktop,n))/(z2-z1)
     >                          * (ztmp-z1)
              else
                z1 = zs + tfac*sgz(ktop-1)
                z2 = z6
                vtop = var(i,j,ktop,n)
     >              + (var(i,j,ktop,n)-var(i,j,ktop-1,n))/(z2-z1)
     >                          * (ztmp-z2)
              end if
              vbot = 0.
              dzbot = 0.

            else 

              z1 = zs + tfac*sgz(kbot)
              z2 = zs + tfac*wgz(kbot+1)
              z3 = zs + tfac*sgz(kbot+1)
              z4 = zs + tfac*sgz(ktop-1)
              z5 = zs + tfac*wgz(ktop)
              z6 = zs + tfac*sgz(ktop)
              dzbot = z2 - zstrt(i,j,n)
              vbot = var(i,j,kbot,n)   
     >            + (var(i,j,kbot+1,n)-var(i,j,kbot,n))/(z3-z1)
     >                        * (zstrt(i,j,n)+0.5*dzbot-z1)
              dztop = zstop(i,j,n) - z5
              if (zstop(i,j,n).eq.ztop) then
                vtop = var(i,j,ktop,n)
              else
                vtop = var(i,j,ktop,n)  
     >              + (var(i,j,ktop,n)-var(i,j,ktop-1,n))/(z6-z4)
     >                     * (zstop(i,j,n)-0.5*dztop-z6)
              end if

            end if

          end if

          vol_int_xyz0(n) = vol_int_xyz0(n) + vtop*dztop*dx*dy
     >                        + vbot*dzbot*dx*dy

          do k=kbot+1,ktop-1
            vol_int_xyz0(n) = vol_int_xyz0(n) + var(i,j,k,n)*dx*dy
     >                       *tfac*(wgz(k+1)-wgz(k))
          enddo

        end if

      enddo
      enddo
 
      enddo

      return
      end


