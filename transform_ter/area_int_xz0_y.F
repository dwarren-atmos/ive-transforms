
#include "pointer.icl"

      ive_ptr function area_int_xz0_y(varnm,numvars,i1,i2,z01,z02,j,
     >                       ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names)

#include "default.icl"

      integer numvars,i1,i2,j,ndims,dims(4)
      real stag(4),min(4),max(4),miss,z01,z02
      character*(*) varnm(numvars),dim_names(4),dunits,disunits

      integer istrt,istop,iu,iv,iw,k,strbeg,strend,nt,jfix
      real r_val,y
      character*(80) tmp,tmp2,vunits,vdisunits
      logical new,err,units_err,requ
      ive_ptr getmem,var_y,var,zstrt,zstop

      area_int_xz0_y = 0
      vunits(1:len(vunits)) = ' '
      vdisunits(1:len(vdisunits)) = ' '
      err = .false.

      if (twod.eq.1) j = 0

      istrt = i1 + 1
      istop = i2 
      jfix = j + 1

 79   FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 80   FORMAT(10x,'z01 = ',F10.2,2x,'z02 = ',F10.2)
 90   FORMAT(10x,'at y = ',F10.2)

      tmp = varnm(1)
      do k=2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
      enddo

      write(6,*)
      write(6,*) 'Calculating the area integral of field ',
     >                  tmp(strbeg(tmp):strend(tmp))
      write(6,*) 'over the domain defined by:  '
      write(6,79) plmin(1)+i1*dx, plmin(1)+i2*dx
      write(6,80) z01,z02
      write(6,90) plmin(2)+(jfix-1)*dy
      write(6,*)

      if (((jfix.lt.2).and.(twod.ne.1))
     >    .or.(jfix.gt.ny).or.(istrt.lt.1).or.
     >        (istop.gt.nx).or.(z01.lt.0.).or.(z02.gt.wgz_max)) then
        write(6,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries.'
        return
      else if ((z01.ge.z02).or.(istrt.gt.istop)) then
        write(6,*) 'Error:  improper specification of domain of ',
     >                     'integration.'
        return
      end if

      var = getvar('zeta',ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot access field zeta.'
        return
      end if 
     
      zstrt = getmem(dims(1)*dims(4))
      if (zstrt.eq.0) then
        write(6,*) 'Error:  cannot allocate space for zstrt.'
        return
      end if

      zstop = getmem(dims(1)*dims(4))
      if (zstop.eq.0) then
        write(6,*) 'Error:  cannot allocate space for zstop.'
        call freemem(zstrt)
        return
      end if 

      call zrange_y(%val(var),%val(zstrt),%val(zstop),%val(zbot_v),
     >                %val(sgz),z01,z02,ztop,jfix,dims(1),dims(2),
     >                dims(3),dims(4),twod,err)
      if (err) then
        call freemem(zstrt)
        call freemem(zstop)
        return
      end if

      var = getvar(varnm(1),ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot compute field ',
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

      var_y = getmem(dims(1)*dims(3)*nt)
      if (var_y.eq.0) then
        write(6,*) 'Error:  could not allocate space for var_y.'
        return
      end if
      
      do k=1,dims(1)*dims(3)*nt
        call s_val(%val(var_y),k,1.)
      enddo

      call v_at_y(%val(var_y),%val(var),jfix,dims(1),dims(2),dims(3),
     >                            dims(4),iu,iv,iw,twod)

      vunits = dunits
      vdisunits = disunits

      k = 2

 348  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
        if (var.eq.0) then
          write(6,*) 'Error:  cannot compute field ',
     >             varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          call freemem(zstrt)
          call freemem(zstop)
          call freemem(var_y)
          return
        else if (dims(4).ne.nt) then
          write(6,*) 'Error:  fields ',
     >        varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                        ' and ',
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >         'have different time dimensions.'
          call freemem(zstrt)
          call freemem(zstop)
          call freemem(var_y)
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

        call v_at_y(%val(var_y),%val(var),jfix,dims(1),dims(2),
     >                   dims(3),dims(4),iu,iv,iw,twod)

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

      area_int_xz0_y = getmem(nt)
      if (area_int_xz0_y.eq.0) then
        write(6,*) 'Error:  cannot allocate space for area_int_xz0_y.'
        call freemem(zstrt)
        call freemem(zstop)
        call freemem(var_y)
      end if

      call area_int_xz0_y_calc(%val(area_int_xz0_y),%val(var_y),
     >              %val(zstrt),%val(zstop),%val(zbot_v),%val(sgz),
     >              %val(wgz),ztop,dx,dy,istrt,istop,jfix,dims(1),
     >              dims(2),dims(3),dims(4),err,twod)

      if (err) then
        call freemem(zstrt)
        call freemem(zstop)
        call freemem(var_y)
        call freemem(area_int_xz0_y)
        area_int_xz0_y = 0
        return
      end if

      do k=1,3
        dims(k) = 1
        stag(k) = 0.
        min(k) = 0.
        max(k) = 0. 
        dim_names(k) = 'one'
      enddo

      vunits = vunits(strbeg(vunits):strend(vunits))//'*m*m'
      vdisunits = vdisunits(strbeg(vdisunits):strend(vdisunits))
     >                  //'*m*m' 

      call units_smpfy(vunits,dunits,units_err)
      if (units_err) dunits = ' '
      call units_smpfy(vdisunits,disunits,units_err)
      if (units_err) disunits = ' '

      call freemem(zstrt)
      call freemem(zstop)
      call freemem(var_y)

      return
      end

*-------------------------------------------------------------------------

      subroutine zrange_y(zeta,zstrt,zstop,zs_v,sgz,z01,z02,ztop,jfix,
     >                           nx,ny,nz,nt,twod,err)

      integer jfix,nx,ny,nz,nt,twod
      real zeta(nx,ny,nz,nt),zstrt(nx,nt),zstop(nx,nt)
      real zs_v(nx,ny+1-twod),sgz(nz),z01,z02,ztop
      logical err

      integer i,k,n,kptr,kplus
      real za,zb,z0a,z0b,tfac,z0_node,zs

      err = .false.

      do n=1,nt
      do i=1,nx

        zs = zs_v(i,jfix)
        tfac = (ztop-zs)/ztop

        kptr = 1
        kplus = 0

        if (z01.eq.0.) then

          zstrt(i,n) = 0. 

        else

 332      if ((kplus.eq.0).and.(kptr.le.nz)) then
            z0_node = zs + tfac*sgz(kptr) 
     >         - 0.5*(zeta(i,jfix,kptr,n)+zeta(i,jfix-1+twod,kptr,n))
            if (z0_node.ge.z01) then
              kplus = kptr
            else
              kptr = kptr+1
            end if
            goto 332
          end if

          if (kplus.eq.0) then
            write(6,*) 'Error:  z01 surface not found.'
            err = .true.
            return
          end if

          if (kplus.eq.1) then 
            za = zs + tfac*sgz(1)
            zb = zs + tfac*sgz(2)
            z0a = za - 0.5*(zeta(i,jfix,1,n)
     >                            + zeta(i,jfix-1+twod,1,n))
            z0b = zb - 0.5*(zeta(i,jfix,2,n)
     >                            + zeta(i,jfix-1+twod,2,n))
            zstrt(i,n) = za + (zb-za)/(z0b-z0a)*(z01-z0a)
          else
            za = zs + tfac*sgz(kplus-1)
            zb = zs + tfac*sgz(kplus)
            z0a = za - 0.5*(zeta(i,jfix,kplus-1,n)
     >                            + zeta(i,jfix-1+twod,kplus-1,n))
            z0b = zb - 0.5*(zeta(i,jfix,kplus,n)
     >                            + zeta(i,jfix-1+twod,kplus,n))
            zstrt(i,n) = zb + (zb-za)/(z0b-z0a)*(z01-z0b)
          end if

          if (zstrt(i,n).lt.zs) zstrt(i,n) = 0.

        end if

        kplus = 0

 335    if ((kplus.eq.0).and.(kptr.le.nz)) then
          z0_node = zs + tfac*sgz(kptr)   
     >       - 0.5*(zeta(i,jfix,kptr,n)+zeta(i,jfix-1+twod,kptr,n))
          if (z0_node.ge.z02) then
            kplus = kptr
          else
            kptr = kptr+1
          end if
          goto 335
        end if

        if (kplus.eq.0) then
          write(6,*) 'Error:  z02 surface not found.'
          err = .true.
          return
        end if

        if (kplus.eq.1) then 
          za = zs + tfac*sgz(1)
          zb = zs + tfac*sgz(2)
          z0a = za - 0.5*(zeta(i,jfix,1,n)
     >                          + zeta(i,jfix-1+twod,1,n))
          z0b = zb - 0.5*(zeta(i,jfix,2,n)
     >                          + zeta(i,jfix-1+twod,2,n))
          zstop(i,n) = za + (zb-za)/(z0b-z0a)*(z02-z0a)
        else
          za = zs + tfac*sgz(kplus-1)
          zb = zs + tfac*sgz(kplus)
          z0a = za - 0.5*(zeta(i,jfix,kplus-1,n)
     >                          + zeta(i,jfix-1+twod,kplus-1,n))
          z0b = zb - 0.5*(zeta(i,jfix,kplus,n)
     >                          + zeta(i,jfix-1+twod,kplus,n)) 
          zstop(i,n) = zb + (zb-za)/(z0b-z0a)*(z02-z0b)
        end if

        if (zstop(i,n).le.zs) zstop(i,n) = 0.

      enddo 
      enddo

      return
      end

*-------------------------------------------------------------------------

      subroutine area_int_xz0_y_calc(area_int_xz0_y,var_y,zstrt,
     >            zstop,zs_v,sgz,wgz,ztop,dx,dy,istrt,istop,jfix,
     >            nx,ny,nz,nt,err,twod)

      integer istrt,istop,jfix,nx,ny,nz,nt,twod
      real area_int_xz0_y(nt),var_y(nx,nz,nt),zstrt(nx,nt)
      real zstop(nx,nt),zs_v(nx,ny+1-twod),sgz(nz),wgz(nz+1)
      real dx,dy,ztop
      logical err

      integer i,k,n,ktop,kbot
      real dztop,dzbot,z1,z2,z3,z4,z5,z6,tfac,vbot,vtop,zs,ztmp

      do n=1,nt

      area_int_xz0_y(n) = 0.

      do i=istrt,istop

        zs = zs_v(i,jfix)
        tfac = (ztop-zs)/ztop

        if (zstop(i,n).gt.zs) then

          if (zstrt(i,n).lt.zs) then

            call bot_top(zs,wgz,ztop,zstrt(i,n),zstop(i,n),
     >                             kbot,ktop,nz)

            if (ktop.eq.1) then
 
              z1 = zs + tfac*sgz(1)
              z2 = zs + tfac*sgz(2)
              dztop = zstop(i,n) - zs
              vtop = var_y(i,1,n) + (var_y(i,2,n)-var_y(i,1,n))
     >                      / (z2-z1)*(zs+0.5*dztop-z1)
              vbot = 0.
              dzbot = 0.

            else

              vbot = var_y(i,1,n)
              dzbot = tfac*wgz(2)
              z4 = zs + tfac*sgz(ktop-1)
              z5 = zs + tfac*wgz(ktop)
              z6 = zs + tfac*sgz(ktop)
              dztop = zstop(i,n) - z5
              if (zstop(i,n).eq.ztop) then
                vtop = var_y(i,ktop,n)
              else
                vtop = var_y(i,ktop,n) 
     >               + (var_y(i,ktop,n)-var_y(i,ktop-1,n))/(z6-z4)
     >                      * (zstop(i,n)-0.5*dztop-z6)
              end if

            end if

          else

            call bot_top(zs,wgz,ztop,zstrt(i,n),zstop(i,n),
     >                            kbot,ktop,nz)

            if ((kbot.eq.1).and.(ktop.eq.1)) then

              z1 = zs + tfac*sgz(1)
              z2 = zs + tfac*sgz(2)
              dztop = zstop(i,n) - zstrt(i,n)
              vtop = var_y(i,1,n) + (var_y(i,2,n)-var_y(i,1,n))
     >                   / (z2-z1)*(zstop(i,n)-0.5*dztop-z1)
              vbot = 0.
              dzbot = 0.

            else if ((kbot.eq.nz).and.(ktop.eq.nz)) then

              z1 = zs + tfac*sgz(nz-1)
              z2 = zs + tfac*sgz(nz)
              dztop = zstop(i,n) - zstrt(i,n)
              vtop = var_y(i,nz,n) + (var_y(i,nz,n)-var_y(i,nz-1,n))
     >                   / (z2-z1)*(zstop(i,n)-dztop/2.-z2)
              vbot = 0.
              dzbot = 0.

            else if (kbot.eq.ktop) then

              dztop = zstop(i,n) - zstrt(i,n)
              ztmp = zstop(i,n) - 0.5*dztop
              z6 = zs + tfac*sgz(ktop)
              if (ztmp.gt.z6) then
                z1 = z6
                z2 = zs + tfac*sgz(ktop+1)
                vtop = var_y(i,ktop,n)
     >              + (var_y(i,ktop+1,n)-var_y(i,ktop,n))/(z2-z1)
     >                          * (ztmp-z1)
              else
                z1 = zs + tfac*sgz(ktop-1)
                z2 = z6
                vtop = var_y(i,ktop,n)
     >              + (var_y(i,ktop,n)-var_y(i,ktop-1,n))/(z2-z1)
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
              dzbot = z2 - zstrt(i,n)
              vbot = var_y(i,kbot,n)   
     >            + (var_y(i,kbot+1,n)-var_y(i,kbot,n))/(z3-z1)
     >                        * (zstrt(i,n)+0.5*dzbot-z1)
              dztop = zstop(i,n) - z5
              if (zstop(i,n).eq.ztop) then
                vtop = var_y(i,ktop,n)
              else
                vtop = var_y(i,ktop,n)  
     >              + (var_y(i,ktop,n)-var_y(i,ktop-1,n))/(z6-z4)
     >                        * (zstop(i,n)-0.5*dztop-z6)
              end if

            end if

          end if

          area_int_xz0_y(n) = area_int_xz0_y(n) + vtop*dztop*dx
     >                        + vbot*dzbot*dx

          do k=kbot+1,ktop-1
            area_int_xz0_y(n) = area_int_xz0_y(n) + var_y(i,k,n)*dx
     >                       *tfac*(wgz(k+1)-wgz(k))
          enddo

        end if

      enddo

      enddo

      return
      end

