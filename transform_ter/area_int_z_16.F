
#include "pointer.icl"

      ive_ptr function area_int_z(varnm,numvars,i1,i2,j1,j2,k1,
     >                ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,k1,numvars
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dim_names(4),dunits,disunits     

      integer istrt,istop,jstrt,jstop,iu,iv,iw,k
      integer strbeg,strend,nt,n,j
      real z,type,h,a1,a2,b1,b2,y_len,xc,yc
      character*(80) tmp,tmp2,vunits,vdisunits
      logical err,new,errchck,units_err,requ
      ive_ptr getmem,gz,var_z,var_zs,var,area,var_c,zs_q

      area_int_z = 0
      vunits = ' '
      vdisunits = ' '
      err = .false.

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if

      istrt = i1 + 1
      istop = i2 
      jstrt = j1 + 1
      jstop = j2 
      z = k1*dz

 86   FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 87   FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)
 88   FORMAT(10x,'at z = ',F10.2)

      tmp = varnm(1)
      do n = 2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >          varnm(n)(strbeg(varnm(n)):strend(varnm(n)))
      enddo

      write(*,*) 
      write(*,*) 'Calculating the area integral of field ',
     >               tmp(strbeg(tmp):strend(tmp))
      write(*,*) 'over the domain defined by:  '
      write(*,86) plmin(1)+i1*dx, plmin(1)+i2*dx  
      write(*,87) plmin(2)+j1*dy, plmin(2)+j2*dy
      write(*,88) z 
      write(*,*)

      if ((istrt.lt.1).or.(istop.gt.nx).or.(jstrt.lt.1).or.
     >                (jstop.gt.ny).or.(z.lt.0).or.
     >                   (z.gt.wgz_max)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                          'boundaries'
        return
      else if ((istrt.gt.istop).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                          'integration.'
        return
      end if

      err = errchck(%val(zbot_p),%val(sgz),ztop,z,nx,ny,nsgz)
      if (err) then
        write(*,*)  'Error:  specified value for z not valid. '
        write(*,*)  '        z must not exceed the highest ',
     >              'thermodynamic height for all x,y.'
        return
      end if

      var = getvar(varnm(1),ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot compute field ',
     >         varnm(1)(strbeg(varnm(1)):strend(varnm(1)))
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

      if (iw.eq.1) then
        gz = wgz
      else
        gz = sgz
      end if

      nt = dims(4)
      var_z = getmem(dims(1)*dims(2)*nt)
      if (var_z.eq.0) then
        write(*,*) 'Error:  could not allocate space for var_z.'
        return
      end if

      var_zs = getmem(dims(1)*dims(2)*nt)
      if (var_zs.eq.0) then
        write(*,*) 'Error:  could not allocate space for var_zs.'
        return
      end if

      do j=1,dims(1)*dims(2)*nt
        call s_val(%val(var_z),j,1.)
        call s_val(%val(var_zs),j,1.)
      enddo

      miss = -777.

      call v_2_z(%val(var_z),%val(var),%val(zbot_p),%val(gz),z,ztop,
     >               miss,dims(1),dims(2),dims(3)+iw,nt,iu,iv,iw)
 
      call v_2_zs(%val(var_zs),%val(var),dims(1),dims(2),dims(3),
     >                             nt,iu,iv,iw)

      vunits = dunits
      vdisunits = disunits

      k = 2

 658  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                              disunits,dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  cannot compute field ',
     >           varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          call freemem(var_z)
          call freemem(var_zs)
          return
        else if (nt.ne.dims(4)) then
          write(*,*) 'Error:  fields ',
     >         varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                     ' and ',
     >         varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >              ' have different time domains.'
          call freemem(var_z)
          call freemem(var_zs)
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
  
        if (iw.eq.1) then
          gz = wgz
        else
          gz = sgz
        end if

        miss = -777.

        call v_2_z(%val(var_z),%val(var),%val(zbot_p),%val(gz),z,ztop,
     >               miss,dims(1),dims(2),dims(3)+iw,nt,iu,iv,iw)

        call v_2_zs(%val(var_zs),%val(var),dims(1),dims(2),dims(3),
     >                                 nt,iu,iv,iw)

        tmp = vunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >                  dunits(strbeg(dunits):strend(dunits))
        tmp = vdisunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vdisunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >              disunits(strbeg(disunits):strend(disunits))

        k = k+1

        goto 658

      end if

      area = getmem(dims(1)*dims(2))
      if (area.eq.0) then
        write(*,*) 'Error:  could not allocate space for area.'
        call freemem(var_z)
        call freemem(var_zs)
        return
      end if

      var_c = getmem(dims(1)*dims(2)*nt)
      if (var_c.eq.0) then
        write(*,*) 'Error:  could not allocate space for var_c.'
        call freemem(var_z)
        call freemem(var_zs)
        call freemem(area)
        return
      end if 

      zs_q = getmem((4*dims(1)+1)*(4*dims(2)+1))
      if (zs_q.eq.0) then
        write(*,*) 'Error:  could not allocate space for zs_q.'
        call freemem(var_z)
        call freemem(var_zs)
        call freemem(area)
        call freemem(var_c)
      end if

      call ter_params(type,h,a1,a2,b1,b2,y_len,xc,yc)

      call zs_qud(%val(zs_q),dx,dy,plmin(1),plmin(2),type,h,
     >               a1,a2,b1,b2,y_len,xc,yc,dims(1),dims(2))

      call int_z_setup(%val(area),%val(var_c),%val(zs_q),
     >            %val(var_zs),z,dx,dy,dims(1),dims(2),nt)

      call freemem(var_zs)
      call freemem(zs_q)

      area_int_z = getmem(nt)
      if (area_int_z.eq.0) then
        write(*,*) 'Error:  could not allocate space for area_int_z.'
        call freemem(var_z)
        call freemem(area)
        call freemem(var_c)
        return
      end if

      call area_int_z_calc(%val(area_int_z),%val(var_z),%val(area),
     >              %val(var_c),%val(zbot_p),z,miss,istrt,istop,
     >              jstrt,jstop,dims(1),dims(2),nt,twod)

      do k=1,3
        dims(k) = 1
        stag(k) = 0.
        min(k) = 0.
        max(k) = 0.
        dim_names(k) = 'one'
      enddo

      vunits = vunits(strbeg(vunits):strend(vunits))//'*m*m'
      vdisunits = vdisunits(strbeg(vdisunits):strend(vdisunits))
     >              //'*m*m'

      call units_smpfy(vunits,dunits,units_err)
      if (units_err) dunits = ' '
      call units_smpfy(vdisunits,disunits,units_err)
      if (units_err) disunits = ' '

      call freemem(var_z)
      call freemem(area)
      call freemem(var_c)

      return
      end 

*---------------------------------------------------------------------

      subroutine v_2_z(var_z,var,zs,gz,z,ztop,miss,nx,ny,nz,nt,
     >                        iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real var_z(nx,ny,nt),var(nx+iu,ny+iv,nz,nt),zs(nx,ny)
      real gz(nz),z,ztop,miss
      
      integer k,kplus,i,j,n
      real z1,z2,var1,var2,tfac,znode

      do n=1,nt
      do j=1,ny
      do i=1,nx

        if (z.ge.zs(i,j)) then

          tfac = (ztop-zs(i,j))/ztop

          kplus = 0
          k = 1

 876      if ((kplus.eq.0).and.(k.le.nz)) then
            znode = zs(i,j) + tfac*gz(k)
            if (znode.gt.z) then
              kplus = k
            else 
              k = k+1
            end if
            goto 876
          end if

          if (kplus.eq.0) then
            write(*,*) 'Error:  kplus zero in v_2_z.'
            stop
          else if (kplus.eq.1) then
            z1 = zs(i,j) + tfac*gz(1)
            z2 = zs(i,j) + tfac*gz(2)
            var1 = (var(i,j,1,n)+var(i+iu,j+iv,1,n))/2.
            var2 = (var(i,j,2,n)+var(i+iu,j+iv,2,n))/2.
          else
            z1 = zs(i,j) + tfac*gz(kplus-1)
            z2 = zs(i,j) + tfac*gz(kplus)
            var1 = (var(i,j,kplus-1,n)+var(i+iu,j+iv,kplus-1,n))/2.
            var2 = (var(i,j,kplus,n)+var(i+iu,j+iv,kplus,n))/2.
          end if

          var_z(i,j,n) = var_z(i,j,n)
     >             *(var1 + (var2-var1)*(z-z1)/(z2-z1))

        else

          var_z(i,j,n) = miss

        end if

      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine v_2_zs(v_zs,var,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real v_zs(nx,ny,nt),var(nx+iu,ny+iv,nz+iw,nt)

      integer i,j,n

      if (iw.eq.1) then

        do n=1,nt
        do j=1,ny
        do i=1,nx
          v_zs(i,j,n) = v_zs(i,j,n)*var(i,j,1,n)
        enddo
        enddo
        enddo

      else

        do n=1,nt
        do j=1,ny
        do i=1,nx
          v_zs(i,j,n) = v_zs(i,j,n)
     >        * 0.25*(3.*(var(i,j,1,n)+var(i+iu,j+iv,1,n)) 
     >                  - var(i,j,2,n) - var(i+iu,j+iv,2,n))
        enddo
        enddo
        enddo

      end if

      return
      end

*-------------------------------------------------------------------------

      subroutine zs_qud(zs_q,dx,dy,x0,y0,type,h,a1,a2,b1,b2,y_len,
     >                                  xc,yc,nx,ny)

      integer nx,ny
      real zs_q(4*nx+1,4*ny+1),dx,dy,x0,y0,type,h,a1,a2,b1,b2
      real y_len,xc,yc

      integer i,j
      real x,y,ytmp,a,b,rad,pi,ly
      logical requ

      ly = y_len
      pi = 2.*asin(1.)

      do j=1,4*ny+1

        ytmp = y0 + 0.25*dy*real(j-1) - yc

        b = b1

        if (ytmp.lt.-ly) then
          y = -ly - ytmp
          b = b2
        else if (ytmp.gt.ly) then
          y = ytmp - ly
        else
          y = 0.
        end if

        do i=1,4*nx+1

          x = x0 + 0.25*dx*real(i-1) - xc

          if (x.le.0.) then
            a = a1
          else
            a = a2
          end if

          rad = sqrt((x/4./a)**2 + (y/4./b)**2)

          if (rad.lt.1.) then
            zs_q(i,j) = h/16.*(1+cos(pi*rad))**4
          else
            zs_q(i,j) = 0.
          end if

        enddo

      enddo

      return
      end

*--------------------------------------------------------------------------

      subroutine int_z_setup(area,v_c,zs,v_zs,z,dx,dy,nx,ny,nt)

      integer nx,ny,nt
      real zs(4*nx+1,4*ny+1),v_zs(nx,ny,nt),area(nx,ny)
      real v_c(nx,ny,nt),z,dx,dy

      integer sa,sb
      real xa,ya,xb,yb,z0(5)

      integer numpts
      real zmid,zave(4),x(8),y(8),len(0:7),hlen

      integer i,j,n,q,k,flga,flgb,cseg
      real zsmax,qdx,qdy,hdx,hdy,part,xc,yc,slen,w1,w2,w3,w4 
      logical requ

      zsmax = 0.
      do j=1,4*ny+1
      do i=1,4*nx+1
        if (zs(i,j).gt.zsmax) zsmax = zs(i,j)
      enddo
      enddo

      if (z.gt.zsmax) then
        do j=1,ny
        do i=1,nx
          area(i,j) = dx*dy
          do n=1,nt
            v_c(i,j,n) = 0.
          enddo
        enddo
        enddo
        return
      end if

      qdx = 0.25*dx
      qdy = 0.25*dy
      hdx = 0.5*dx
      hdy = 0.5*dy

      do j=1,ny
      do i=1,nx

        area(i,j) = 0.

        do q=1,8
          x(q) = -1.
          y(q) = -1.
        enddo
        numpts = 0

        do k = 1,16

          call side_info(zs,z,dx,dy,nx,ny,i,j,k,sa,sb,z0,xa,ya,xb,yb)

          do n=1,4
            zave(n) = 0.5*(z0(n)+z0(n+1))
          enddo
          zmid = 0.25*(z0(1)+z0(2)+z0(3)+z0(4))

          if (sa.eq.0) then

            part = 0.0625*dx*dy
            if (zmid.gt.0.) part = 0.

          else 
  
            if ((sa.eq.1).and.(sb.eq.3)) then

              part = 0.5*qdx*(ya+yb)
              if (zave(4).gt.0.) part = qdx*qdy - part
  
            else if ((sa.eq.2).and.(sb.eq.4)) then
         
              part = 0.5*(xa+xb)*qdy
              if (zave(1).gt.0.) part = qdx*qdy - part
  
            else if ((sa.eq.1).and.(sb.eq.2)) then
     
              part = 0.5*(qdy-ya)*xb
              if (z0(2).gt.0.) part = qdx*qdy - part
  
            else if ((sa.eq.2).and.(sb.eq.3)) then
  
              part = 0.5*(qdx-xa)*(qdy-yb)
              if (z0(3).gt.0.) part = qdx*qdy - part

            else if ((sa.eq.3).and.(sb.eq.4)) then
   
              part = 0.5*(qdx-xb)*ya
              if (z0(4).gt.0.) part = qdx*qdy - part
  
            else if ((sa.eq.1).and.(sb.eq.4)) then
  
              part = 0.5*ya*xb
              if (z0(1).gt.0.) part = qdx*qdy - part
  
            else
  
              write(6,*) 'Error:  case not recognized in int_z_setup. '
              write(6,*) 'Stopping run.'
              stop
  
            end if

            xa = xa + int((k-1)/4)*qdx
            xb = xb + int((k-1)/4)*qdx 

            ya = ya + mod(k-1,4)*qdy
            yb = yb + mod(k-1,4)*qdy

            flga = 0
            flgb = 0
            n = 1
 124        if (n.le.numpts) then
              if (requ(xa,x(n)).and.requ(ya,y(n))) flga = 1
              if (requ(xb,x(n)).and.requ(yb,y(n))) flgb = 1
              n = n+1
              goto 124
            end if

            if (flga.eq.0) then
              numpts = numpts+1 
              if (numpts.gt.8) then
                write(6,*) 'Error:  too many intersection points in ',
     >                                'int_z_setup.'
                write(6,*) 'Stopping run.'
                stop
              end if
              x(numpts) = xa
              y(numpts) = ya
            end if

            if (flgb.eq.0) then
              numpts = numpts+1
              if (numpts.gt.8) then
                write(6,*) 'Error:  too many intersection points in ',
     >                                'int_z_setup.'
                write(6,*) 'Stopping run.'
                stop
              end if
              x(numpts) = xb
              y(numpts) = yb
            end if
              
          end if

          area(i,j) = area(i,j) + part

        enddo

        if ((numpts.eq.0).and.requ(area(i,j),dx*dy)) then

          do n=1,nt
            v_c(i,j,n) = 0.
          enddo

        else if (numpts.eq.0) then
 
          do n=1,nt
            v_c(i,j,n) = v_zs(i,j,n)
          enddo

        else

          do k = 1,numpts-1
          do n = 1,numpts-k
            if (x(n+1).lt.x(n)) then
              xa = x(n)
              ya = y(n)
              x(n) = x(n+1)
              y(n) = y(n+1)
              x(n+1) = xa
              y(n+1) = ya
            end if
          enddo
          enddo

          len(0) = 0. 
          do k = 1,numpts-1
            len(k) = sqrt((x(k)-x(k+1))**2 + (y(k)-y(k+1))**2)
     >                              + len(k-1)
          enddo
          hlen = 0.5*len(numpts-1)

          cseg = 0
          k = 1
 456      if ((cseg.eq.0).and.(k.le.numpts-1)) then
            if (hlen.lt.len(k)) then
              cseg = k
            else
              k = k+1
            end if
            goto 456 
          end if

          if (cseg.eq.0) then
            write(6,*) 'Error:  did not find contour segment ',
     >                            'in int_z_setup.'
            write(6,*) 'Stopping run.'
            stop
          end if 

          slen = sqrt((x(cseg)-x(cseg+1))**2 
     >                     + (y(cseg)-y(cseg+1))**2)
          xc = x(cseg) + (x(cseg+1)-x(cseg))/slen
     >                  *(hlen-len(cseg-1))
          yc = y(cseg) + (y(cseg+1)-y(cseg))/slen
     >                  *(hlen-len(cseg-1))

          if ((xc.lt.0.).or.(xc.gt.dx).or.
     >             (yc.lt.0.).or.(yc.gt.dy)) then
            write(6,*) 'Error:  xc or yc out of bounds.'
            write(6,*) 'Stopping run.'
            stop
          else if ((xc.lt.hdx).and.(yc.lt.hdy)) then
            q = 1
            xc = hdx + xc
            yc = hdy + yc
          else if (xc.lt.hdx) then
            q = 2
            xc = hdx + xc
            yc = yc - hdy
          else if (yc.lt.hdy) then
            q = 4
            xc = xc - hdx
            yc = hdy + yc
          else
            q = 3
            xc = xc - hdx
            yc = yc - hdy
          end if
     
          w1 = (dx-xc)*(dy-yc)/dx/dy
          w2 = (dx-xc)*yc/dx/dy
          w3 = xc*yc/dx/dy
          w4 = xc*(dy-yc)/dx/dy
 
          if (q.eq.1) then
            do n=1,nt
              v_c(i,j,n) = w1*v_zs(i-1,j-1,n) + w2*v_zs(i-1,j,n)
     >                   + w3*v_zs(i,j,n) + w4*v_zs(i,j-1,n)
            enddo
          else if (q.eq.2) then
            do n=1,nt
              v_c(i,j,n) = w1*v_zs(i-1,j,n) + w2*v_zs(i-1,j+1,n)
     >                   + w3*v_zs(i,j+1,n) + w4*v_zs(i,j,n)
            enddo
          else if (q.eq.3) then
            do n=1,nt
              v_c(i,j,n) = w1*v_zs(i,j,n) + w2*v_zs(i,j+1,n)
     >                   + w3*v_zs(i+1,j+1,n) + w4*v_zs(i+1,j,n)
            enddo
          else 
            do n=1,nt
              v_c(i,j,n) = w1*v_zs(i,j-1,n) + w2*v_zs(i,j,n)
     >                   + w3*v_zs(i+1,j,n) + w4*v_zs(i+1,j-1,n)
            enddo
          end if 

        end if
           
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine side_info(zs,z,dx,dy,nx,ny,i,j,q,sa,sb,z0,
     >                              xa,ya,xb,yb)

      integer nx,ny,i,j,q,sa,sb
      real zs(4*nx+1,4*ny+1),z,dx,dy,z0(5),xa,ya,xb,yb

      integer seg1,seg2,k,flg,i1,j1
      real x(4),y(4)

      i1 = 4*i - 3 + int((q-1)/4)
      j1 = 4*j - 3 + mod(q-1,4)

      z0(1) = zs(i1,j1) - z
      z0(2) = zs(i1,j1+1) - z
      z0(3) = zs(i1+1,j1+1) - z
      z0(4) = zs(i1+1,j1) - z
      z0(5) = z0(1)

      flg = 0
      k = 1
 222  if ((flg.eq.0).and.(k.le.4)) then
        if ((z0(k).eq.0.).and.(z0(k+1).eq.0.)) then
          flg = 1
        else 
          k = k+1
        end if
        goto 222
      end if

      if (flg.eq.1) then
        sa = 0
        sb = 0
        return
      end if

      x(1) = 0.
      x(2) = 0.
      x(3) = 0.25*dx
      x(4) = 0.25*dx

      y(1) = 0.
      y(2) = 0.25*dy
      y(3) = 0.25*dy
      y(4) = 0.
 
      seg1 = 0
      k = 1
 344  if ((seg1.eq.0).and.(k.le.4)) then
        if ((z0(k)*z0(k+1).lt.0.).or.(z0(k).eq.0.)) then
          seg1 = k
        else
          k = k+1
        end if
        goto 344
      end if

      seg2 = 0
      k = seg1+1
 349  if ((seg2.eq.0).and.(k.le.4)) then
        if ((z0(k)*z0(k+1).lt.0.).or.(z0(k).eq.0.)) then
          seg2 = k
        else  
          k = k+1
        end if
        goto 349
      end if

      if (seg1.eq.0) then
        sa = 0
        sb = 0
        return
      else if (seg2.eq.0) then
        sa = 0
        sb = 0
        return
      else if (seg1.eq.1) then
        sa = 1
        xa = x(seg1)  
        ya = y(seg1) + (0.-z0(seg1))/(z0(seg1+1)-z0(seg1))*0.25*dy
      else if (seg1.eq.2) then
        sa = 2
        xa = x(seg1) + (0.-z0(seg1))/(z0(seg1+1)-z0(seg1))*0.25*dx
        ya = y(seg1)
      else if (seg1.eq.3) then
        sa = 3
        xa = x(seg1)
        ya = y(seg1) - (0.-z0(seg1))/(z0(seg1+1)-z0(seg1))*0.25*dy
      else 
        write(6,*) 'Error:  found only one intersection in ',
     >                  'side_info.  Stopping run.'
        stop
      end if

      if (seg2.eq.2) then
        sb = 2
        xb = x(seg2) + (0.-z0(seg2))/(z0(seg2+1)-z0(seg2))*0.25*dx
        yb = y(seg2)
      else if (seg2.eq.3) then
        sb = 3
        xb = x(seg2)
        yb = y(seg2) - (0.-z0(seg2))/(z0(seg2+1)-z0(seg2))*0.25*dy
      else if (seg2.eq.4) then
        sb = 4
        xb = x(seg2) - (0.-z0(seg2))/(z0(seg2+1)-z0(seg2))*0.25*dx
        yb = y(seg2)
      else
        write(6,*) 'Error:  found only one intersection in ',
     >                  'side_info.  Stopping run.'
        stop
      end if

      return
      end

*-----------------------------------------------------------------------

      subroutine area_int_z_calc(area_int,v_z,area,v_c,zs_p,z,miss,
     >                istrt,istop,jstrt,jstop,nx,ny,nt,twod)

      integer nx,ny,nt,twod,istrt,istop,jstrt,jstop
      real area_int(nt),v_z(nx,ny,nt),area(nx,ny)
      real v_c(nx,ny,nt),zs_p(nx,ny),z,miss

      integer n
      real tot_area,sum

      do n=1,nt

        if (twod.eq.1) then
          write(*,*) 'Error:  int_z_1d not yet implemented.' 
          stop
        else
          call int_z_2d(sum,tot_area,v_z(1,1,n),area,v_c(1,1,n),
     >                     zs_p,z,istrt,istop,jstrt,jstop,nx,ny)
        end if

        if (tot_area.ne.0.) then
          area_int(n) = sum
        else
          area_int(n) = miss
        end if

      enddo
  
      return
      end

*--------------------------------------------------------------------------

      subroutine int_z_2d(sum,tot_area,v_z,area,v_c,zs_p,z,istrt,
     >                             istop,jstrt,jstop,nx,ny)

      integer nx,ny,istrt,istop,jstrt,jstop
      real sum,tot_area,v_z(nx,ny),area(nx,ny),v_c(nx,ny)
      real zs_p(nx,ny),z

      integer i,j

      sum = 0.
      tot_area = 0.

      do j=jstrt,jstop
      do i=istrt,istop

        if (area(i,j).ne.0.) then

          if (z.ge.zs_p(i,j)) then
            sum = sum + v_z(i,j)*area(i,j)
          else
            sum = sum + v_c(i,j)*area(i,j)
          end if

          tot_area = tot_area + area(i,j)

        end if

      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine ter_params(type,h,a1,a2,b1,b2,y_len,xc,yc)

      real type,h,a1,a2,b1,b2,y_len,xc,yc

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss,r_val
      character*(80) dunits,disunits,dim_names(4)
      logical new

      ive_ptr var,getvar

      var = getvar('type',ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not get type.'
        write(6,*) 'Stopping run.'
        stop
      end if

      type = r_val(%val(var),1) 

      var = getvar('h',ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not get h.'
        write(6,*) 'Stopping run.'
        stop
      end if

      h = r_val(%val(var),1) 

      var = getvar('a1',ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not get a1.'
        write(6,*) 'Stopping run.'
        stop
      end if

      a1 = r_val(%val(var),1) 

      var = getvar('a2',ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not get a2.'
        write(6,*) 'Stopping run.'
        stop
      end if

      a2 = r_val(%val(var),1) 

      var = getvar('b1',ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not get b1.'
        write(6,*) 'Stopping run.'
        stop
      end if

      b1 = r_val(%val(var),1) 

      var = getvar('b2',ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not get b2.'
        write(6,*) 'Stopping run.'
        stop
      end if

      b2 = r_val(%val(var),1) 

      var = getvar('y_len',ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not get y_len.'
        write(6,*) 'Stopping run.'
        stop
      end if

      y_len = r_val(%val(var),1) 

      var = getvar('xc',ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not get xc.'
        write(6,*) 'Stopping run.'
        stop
      end if

      xc = r_val(%val(var),1) 

      var = getvar('yc',ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not get yc.'
        write(6,*) 'Stopping run.'
        stop
      end if

      yc = r_val(%val(var),1) 

      return
      end

*------------------------------------------------------------------------

      logical function errchck(zs,gz,ztop,z,nx,ny,nz)

      integer nx,ny,nz
      real zs(nx,ny),gz(nz),ztop,z,zmin,zmax,tfac

      integer i,j

      errchck = .false.

      do j=1,ny 
      do i=1,nx

        tfac = (ztop-zs(i,j))/ztop
        zmin = 0.
        zmax = zs(i,j) + tfac*gz(nz)

        if ((z.lt.zmin).or.(z.gt.zmax)) then
          errchck = .true.
          return
        end if

      enddo
      enddo

      return
      end

c
c           var = 0.
c           num = 0
c           if (z.ge.zs_p(i-1,j-1)) then
c              var = var + v_z(i-1,j-1)
c              num = num + 1
c           end if
c           if (z.ge.zs_p(i-1,j)) then
c              var = var + v_z(i-1,j)
c              num = num + 1
c           end if
c           if (z.ge.zs_p(i-1,j+1)) then
c              var = var + v_z(i-1,j+1)
c              num = num + 1
c           end if
c           if (z.ge.zs_p(i,j+1)) then
c              var = var + v_z(i,j+1)
c              num = num + 1
c           end if
c           if (z.ge.zs_p(i+1,j+1)) then
c              var = var + v_z(i+1,j+1)
c              num = num + 1
c           end if
c           if (z.ge.zs_p(i+1,j)) then
c              var = var + v_z(i+1,j)
c              num = num + 1
c           end if
c           if (z.ge.zs_p(i+1,j-1)) then
c              var = var + v_z(i+1,j-1)
c              num = num + 1
c           end if
c           if (z.ge.zs_p(i,j-1)) then
c              var = var + v_z(i,j-1)
c              num = num + 1
c           end if
c
