
#include "pointer.icl"

      ive_ptr function area_int_yz0_x(varnm,numvars,j1,j2,z01,z02,i,
     >                       ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names)

#include "default.icl"
 
      integer numvars,j1,j2,i,ndims,dims(4)
      real stag(4),min(4),max(4),miss,z01,z02
      character*(*) varnm(numvars),dim_names(4),dunits,disunits

      integer jstrt,jstop,iu,iv,iw,k,strbeg,strend,nt,ifix,j
      real r_val,x
      character*(80) tmp,tmp2,vunits,vdisunits
      logical new,err,units_err,requ
      ive_ptr getmem,var_x,var,zstrt,zstop

      area_int_yz0_x = 0
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
 80   FORMAT(10x,'z01 = ',F10.2,2x,'z02 = ',F10.2)
 90   FORMAT(10x,'at x = ',F10.2)

      tmp = varnm(1)
      do k=2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
      enddo

      write(*,*)
      write(*,*) 'Calculating the area integral of field ',
     >                  tmp(strbeg(tmp):strend(tmp))
      write(*,*) 'over the domain defined by:  '
      write(*,79) plmin(2)+j1*dy, plmin(2)+j2*dy
      write(*,80) z01,z02
      write(*,90) plmin(1)+(ifix-1)*dx 
      write(*,*)

      if ((ifix.lt.2).or.(ifix.gt.nx).or.(jstrt.lt.1).or.
     >          (jstop.gt.ny).or.(z01.lt.0).or.(z02.gt.wgz_max)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries.'
        return
      else if ((z01.ge.z02).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                     'integration.'
        return
      end if

      var = getvar('zeta',ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot access field zeta.'
        return
      end if

      zstrt = getmem(dims(2)*dims(4))
      if (zstrt.eq.0) then
        write(6,*) 'Error:  cannot allocate space for zstrt.'
        return
      end if

      zstop = getmem(dims(2)*dims(4))
      if (zstop.eq.0) then
        write(6,*) 'Error:  cannot allocate space for zstop.'
        call freemem(zstrt)
        return
      end if

      call zrange_x(%val(var),%val(zstrt),%val(zstop),%val(zbot_u),
     >                %val(sgz),z01,z02,ztop,ifix,dims(1),dims(2),
     >                dims(3),dims(4),err)
      if (err) then
        call freemem(zstrt)
        call freemem(zstop)
        return
      end if

      var = getvar(varnm(1),ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot compute field ',
     >              varnm(1)(strbeg(varnm(1)):strend(varnm(1)))
        return
      end if

      IF (abs(min(1)-plmin(1)).gt.dx) min(1) = min(1) + offset(1) 
      IF (abs(min(2)-plmin(2)).gt.dy) min(2) = min(2) + offset(2)

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
          write(6,*) 'Error:  cannot compute field ',
     >             varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          call freemem(zstrt)
          call freemem(zstop)
          call freemem(var_x)
          return
        else if (dims(4).ne.nt) then
          write(6,*) 'Error:  fields ',
     >        varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                        ' and ',
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >         'have different time dimensions.'
          call freemem(zstrt)
          call freemem(zstop)
          call freemem(var_x)
          return
        end if

        IF (abs(min(1)-plmin(1)).gt.dx) min(1) = min(1) + offset(1) 
        IF (abs(min(2)-plmin(2)).gt.dy) min(2) = min(2) + offset(2)
      
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

c     WRITE(6,*)
c     WRITE(6,*) 'var = ',varnm(k)
c     WRITE(6,*) 'iu = ',iu
c     WRITE(6,*) 'iv = ',iv
c     WRITE(6,*) 'iw = ',iw
c     WRITE(6,*) 'min(1) = ',min(1)
c     WRITE(6,*) 'plmin(1) = ',plmin(1)
c     WRITE(6,*)
   

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

      area_int_yz0_x = getmem(nt)
      if (area_int_yz0_x.eq.0) then
        write(6,*) 'Error:  cannot allocate space for area_int_yz0_x.'
        call freemem(zstrt)
        call freemem(zstop)
        call freemem(var_x)
        return
      end if

      call area_int_yz0_x_calc(%val(area_int_yz0_x),%val(var_x),
     >              %val(zstrt),%val(zstop),%val(zbot_u),%val(sgz),
     >              %val(wgz),ztop,dx,dy,jstrt,jstop,ifix,dims(1),
     >              dims(2),dims(3),dims(4),err)

      if (err) then
        call freemem(zstrt)
        call freemem(zstop)
        call freemem(var_x)
        call freemem(area_int_yz0_x)
        area_int_yz0_x = 0
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
      call freemem(var_x)

      return
      end

*-------------------------------------------------------------------------

      subroutine zrange_x(zeta,zstrt,zstop,zs_u,sgz,z01,z02,ztop,ifix,
     >                           nx,ny,nz,nt,err)

      integer ifix,nx,ny,nz,nt
      real zeta(nx,ny,nz,nt),zstrt(ny,nt),zstop(ny,nt)
      real zs_u(nx+1,ny),sgz(nz),z01,z02,ztop
      logical err

      integer j,k,n,kptr,kplus
      real za,zb,z0a,z0b,tfac,z0_node,zs

      err = .false.
 
      do n=1,nt
      do j=1,ny

        zs = zs_u(ifix,j)
        tfac = (ztop-zs)/ztop

        kptr = 1
        kplus = 0
c	  
c  Finds z01 surface
c
        if (z01.eq.0.) then
 
          zstrt(j,n) = 0.

        else

 442      if ((kplus.eq.0).and.(kptr.le.nz)) then
            z0_node = zs + tfac*sgz(kptr)
     >         - 0.5*(zeta(ifix,j,kptr,n)+zeta(ifix-1,j,kptr,n))
            if (z0_node.ge.z01) then
              kplus = kptr
            else
              kptr = kptr+1
            end if
            goto 442
          end if

          if (kplus.eq.0) then
            write(6,*) 'Error:  z01 surface not found.'
            err = .true.
            return
          end if

          if (kplus.eq.1) then
            za = zs + tfac*sgz(1)
            zb = zs + tfac*sgz(2)
            z0a = za - 0.5*(zeta(ifix,j,1,n)+zeta(ifix-1,j,1,n))
            z0b = zb - 0.5*(zeta(ifix,j,2,n)+zeta(ifix-1,j,2,n))
            zstrt(j,n) = za + (zb-za)/(z0b-z0a)*(z01-z0a)
          else
            za = zs + tfac*sgz(kplus-1)
            zb = zs + tfac*sgz(kplus)
            z0a = za - 0.5*(zeta(ifix,j,kplus-1,n)
     >                               + zeta(ifix-1,j,kplus-1,n))
            z0b = zb - 0.5*(zeta(ifix,j,kplus,n)
     >                               + zeta(ifix-1,j,kplus,n))
            zstrt(j,n) = zb + (zb-za)/(z0b-z0a)*(z01-z0b)
          end if

          if (zstrt(j,n).lt.zs) zstrt(j,n) = 0.

        end if
c
        kplus = 0
c
c  Finds z02 surface
c
 435    if ((kplus.eq.0).and.(kptr.le.nz)) then
          z0_node = zs + tfac*sgz(kptr)
     >       - 0.5*(zeta(ifix,j,kptr,n)+zeta(ifix-1,j,kptr,n))
          if (z0_node.ge.z02) then
            kplus = kptr
          else
            kptr = kptr+1
          end if
          goto 435
        end if
c
        if (kplus.eq.0) then
          write(6,*) 'Error:  z02 surface not found.'
          err = .true.
          return
        end if
c
c  Linear Interpolation to z02
c
        if (kplus.eq.1) then
          za = zs + tfac*sgz(1)
          zb = zs + tfac*sgz(2)
          z0a = za - 0.5*(zeta(ifix,j,1,n)+zeta(ifix-1,j,1,n))
          z0b = zb - 0.5*(zeta(ifix,j,2,n)+zeta(ifix-1,j,2,n))
          zstop(j,n) = za + (zb-za)/(z0b-z0a)*(z02-z0a)
        else
          za = zs + tfac*sgz(kplus-1)
          zb = zs + tfac*sgz(kplus)
          z0a = za - 0.5*(zeta(ifix,j,kplus-1,n)
     >                             + zeta(ifix-1,j,kplus-1,n))
          z0b = zb - 0.5*(zeta(ifix,j,kplus,n)
     >                             + zeta(ifix-1,j,kplus,n))
          zstop(j,n) = zb - (zb-za)/(z0b-z0a)*(z0b-z02)
        end if
c
c  If surface is below ground set to zero
c
        if (zstop(j,n).lt.zs) zstop(j,n) = 0.
c
c	  if(j.eq.1.and.n.eq.1) then
c	    write(*,501) 'j','zstop','zs','zb','za','z0b','z0a',
c     >			'slp','diff'
c	  end if
c
c	  if(j.ge.51.and.j.le.70) then
c	     write(*,500) kplus,zstop(j,n),zs,zb,za,z0b,z0a,
c     >                  (zb-za)/(z0b-z0a),z0b-z02
c	  end if
c
      enddo
      enddo
c	
500	FORMAT(I5,F9.2,F8.1,F9.2,F9.2,F9.2,F9.2,F9.4,F8.1)	
501	FORMAT(A5,A9,A8,A9,A9,A9,A9,A9,A8)
      return
      end

*-------------------------------------------------------------------------

      subroutine area_int_yz0_x_calc(area_int_yz0_x,var_x,zstrt,
     >            zstop,zs_u,sgz,wgz,ztop,dx,dy,jstrt,jstop,ifix,
     >            nx,ny,nz,nt,err)

      integer jstrt,jstop,ifix,nx,ny,nz,nt
      real area_int_yz0_x(nt),var_x(ny,nz,nt),zstrt(ny,nt)
      real zstop(ny,nt),zs_u(nx+1,ny),sgz(nz),wgz(nz+1)
      real dx,dy,ztop
      logical err

      integer j,k,n,ktop,kbot
      real dztop,dzbot,z1,z2,z3,z4,z5,z6,tfac,vbot,vtop,zs,ztmp

      do n=1,nt

      area_int_yz0_x(n) = 0.

      do j=jstrt,jstop

        zs = zs_u(ifix,j)
        tfac = (ztop-zs)/ztop

        if (zstop(j,n).gt.zs) then

          if (zstrt(j,n).lt.zs) then

            call bot_top(zs,wgz,ztop,zstrt(j,n),zstop(j,n),
     >                             kbot,ktop,nz)

            if (ktop.eq.1) then
 
              z1 = zs + tfac*sgz(1)
              z2 = zs + tfac*sgz(2)
              dztop = zstop(j,n) - zs
              vtop = var_x(j,1,n) + (var_x(j,2,n)-var_x(j,1,n))
     >                 / (z2-z1)*(zs+0.5*dztop-z1)
              vbot = 0.
              dzbot = 0.

            else

              vbot = var_x(j,1,n)
              dzbot = tfac*wgz(2)
              z4 = zs + tfac*sgz(ktop-1)
              z5 = zs + tfac*wgz(ktop)
              z6 = zs + tfac*sgz(ktop)
              dztop = zstop(j,n) - z5
              if (zstop(j,n).eq.ztop) then
                vtop = var_x(j,ktop,n)
              else
                vtop = var_x(j,ktop,n) 
     >               + (var_x(j,ktop,n)-var_x(j,ktop-1,n))/(z6-z4)
     >                      * (zstop(j,n)-0.5*dztop-z6)
              end if

            end if

          else

            call bot_top(zs,wgz,ztop,zstrt(j,n),zstop(j,n),
     >                            kbot,ktop,nz)

            if ((kbot.eq.1).and.(ktop.eq.1)) then

              z1 = zs + tfac*sgz(1)
              z2 = zs + tfac*sgz(2)
              dztop = zstop(j,n) - zstrt(j,n)
              vtop = var_x(j,1,n) + (var_x(j,2,n)-var_x(j,1,n))
     >                   / (z2-z1)*(zstop(j,n)-0.5*dztop-z1)
              vbot = 0.
              dzbot = 0.

            else if ((kbot.eq.nz).and.(ktop.eq.nz)) then

              z1 = zs + tfac*sgz(nz-1)
              z2 = zs + tfac*sgz(nz)
              dztop = zstop(j,n) - zstrt(j,n)
              vtop = var_x(j,nz,n) + (var_x(j,nz,n)-var_x(j,nz-1,n))
     >                   / (z2-z1)*(zstop(j,n)-0.5*dztop-z2)
              vbot = 0.
              dzbot = 0.

            else if (kbot.eq.ktop) then

              dztop = zstop(j,n) - zstrt(j,n)
              ztmp = zstop(j,n) - 0.5*dztop
              z6 = zs + tfac*sgz(ktop)
              if (ztmp.gt.z6) then
                z1 = z6
                z2 = zs + tfac*sgz(ktop+1)
                vtop = var_x(j,ktop,n)
     >              + (var_x(j,ktop+1,n)-var_x(j,ktop,n))/(z2-z1)
     >                          * (ztmp-z1)
              else
                z1 = zs + tfac*sgz(ktop-1)
                z2 = z6
                vtop = var_x(j,ktop,n)
     >              + (var_x(j,ktop,n)-var_x(j,ktop-1,n))/(z2-z1)
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
              dzbot = z2 - zstrt(j,n)
              vbot = var_x(j,kbot,n)   
     >            + (var_x(j,kbot+1,n)-var_x(j,kbot,n))/(z3-z1)
     >                        * (zstrt(j,n)+0.5*dzbot-z1)
              dztop = zstop(j,n) - z5
              if (zstop(j,n).eq.ztop) then
                vtop = var_x(j,ktop,n)
              else
                vtop = var_x(j,ktop,n)  
     >              + (var_x(j,ktop,n)-var_x(j,ktop-1,n))/(z6-z4)
     >                         * (zstop(j,n)-0.5*dztop-z6)
              end if

            end if

          end if

          area_int_yz0_x(n) = area_int_yz0_x(n) + vtop*dztop*dy
     >                        + vbot*dzbot*dy

          do k=kbot+1,ktop-1
            area_int_yz0_x(n) = area_int_yz0_x(n) + var_x(j,k,n)*dy
     >                       *tfac*(wgz(k+1)-wgz(k))
          enddo

        end if

      enddo

      enddo

      return
      end


        


