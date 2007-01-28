
#include "pointer.icl"

      ive_ptr function area_int_y(varnm,numvars,i1,i2,k1,k2,j,ndims,
     >                   dims,stag,min,max,miss,dunits,disunits,
     >                   dim_names)

#include "default.icl"
*
*  area_int_y calculates the area integral at constant y of the product
*    of the variables in varnm(numvars).  The domain of integration is 
*    the rectangular region at y = j*dy enclosed by i1*dx to i2*dx along
*    the x-axis and k1*dz to k2*dz along the z-axis.  Here dz is the
*    average grid spacing in the vertical (without terrain).
*
      integer numvars,i1,i2,k1,k2,j,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dim_names(4),dunits,disunits

      integer istrt,istop,iu,iv,iw,k,strbeg,strend,nt,jfix
      real zstrt,zstop,r_val,y,dzbot
      character*(80) tmp,tmp2,vunits,vdisunits
      logical new,err,units_err,requ
      ive_ptr getmem,var_y,var

      area_int_y = 0
      vunits(1:len(vunits)) = ' '
      vdisunits(1:len(vdisunits)) = ' '
      err = .false.

      dzbot = r_val(%val(wgz),2) - r_val(%val(wgz),1)

      if (twod.eq.1) j = 0

      istrt = i1 + 1
      istop = i2 
c      zstrt = k1*dz
c      zstop = k2*dz
      zstrt = k1*dzbot
      zstop = k2*dzbot
      jfix = j + 1

 79   FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 80   FORMAT(10x,'z1 = ',F10.2,2x,'z2 = ',F10.2)
 90   FORMAT(10x,'at y = ',F10.2)

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
      write(*,79) plmin(1)+i1*dx, plmin(1)+i2*dx
      write(*,80) zstrt, zstop
      write(*,90) plmin(2)+(jfix-1)*dy
      write(*,*)

      if (((jfix.lt.2).and.(twod.ne.1))
     >     .or.(jfix.gt.ny+1).or.(istrt.lt.1).or.
     >            (istop.gt.nx).or.(zstrt.lt.0).or.
     >              (zstop.gt.wgz_max)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries.'
        return
      else if ((zstop.lt.zstrt).or.(istrt.gt.istop)) then
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

      var_y = getmem(dims(1)*dims(3)*nt)
      if (var_y.eq.0) then
        write(*,*) 'Error:  could not allocate space for var_y.'
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
          write(*,*) 'Error:  cannot compute field ',
     >             varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          call freemem(var_y)
          return
        else if (dims(4).ne.nt) then
          write(*,*) 'Error:  fields ',
     >        varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                        ' and ',
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >         'have different time dimensions.'
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

      area_int_y = getmem(nt)
      if (area_int_y.eq.0) then
        write(6,*) 'Error:  cannot allocate space for area_int_y.'
        call freemem(var_y)
      end if

      call area_int_y_calc(%val(area_int_y),%val(var_y),%val(zbot_v),
     >              %val(sgz),%val(wgz),ztop,dx,dy,zstrt,zstop,istrt,
     >              istop,jfix,dims(1),dims(2),dims(3),dims(4),
     >                                     err,twod)

      if (err) then
        call freemem(var_y)
        call freemem(area_int_y)
        area_int_y = 0
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

      call freemem(var_y)

      return
      end

*-------------------------------------------------------------------------

      subroutine v_at_y(var_y,var,jfix,nx,ny,nz,nt,iu,iv,iw,twod)

      integer iu,iv,iw,nx,ny,nz,nt,jfix,twod
      real var_y(nx,nz,nt),var(nx+iu,ny+iv,nz+iw,nt)

      integer n,k,i

      if (twod.eq.1) then
        do n=1,nt
        do k=1,nz
        do i=1,nx
          var_y(i,k,n) = var_y(i,k,n) * (var(i,jfix,k,n) +
     >           var(i+iu,jfix,k+iw,n))/2.
        enddo
        enddo
        enddo
      else if (jfix.eq.ny+1) then
        do n=1,nt
        do k=1,nz
        do i=1,nx
          var_y(i,k,n) = var_y(i,k,n) * (var(i,jfix+iv-1,k,n) + 
     >             var(i+iu,jfix+iv-1,k+iw,n))/2.
        enddo
        enddo
        enddo
      else
        do n=1,nt
        do k=1,nz
        do i=1,nx
          var_y(i,k,n) = var_y(i,k,n) *
     >        (var(i,jfix,k,n) + var(i+iu,jfix,k+iw,n)
     >       + var(i,jfix+iv-1,k,n) + var(i+iu,jfix+iv-1,k+iw,n))/4.
        enddo
        enddo
        enddo
      end if
 
      return
      end

*-------------------------------------------------------------------------

      subroutine area_int_y_calc(area_int_y,var_y,zs,sgz,wgz,ztop, 
     >            dx,dy,zstrt,zstop,istrt,istop,jfix,nx,ny,nz,nt,err,
     >                                twod)

      integer istrt,istop,jfix,nx,ny,nz,nt,twod
      real area_int_y(nt),var_y(nx,nz,nt),zs(nx,ny+1-twod),sgz(nz),
     >        wgz(nz+1),dx,dy,zstrt,zstop,ztop
      logical err

      integer i,k,n,ktop,kbot
      real dztop,dzbot,z1,z2,z3,z4,z5,z6,tfac,vbot,vtop,z0,ztmp

      do n=1,nt

      area_int_y(n) = 0.

      do i=istrt,istop

        z0 = zs(i,jfix)
        tfac = (ztop - z0)/ztop

        if (zstop.gt.z0) then

          if (zstrt.lt.z0) then

            call bot_top(z0,wgz,ztop,zstrt,zstop,kbot,ktop,nz)

            if (ktop.eq.1) then
 
              z1 = z0 + tfac*sgz(1)
              z2 = z0 + tfac*sgz(2)
              dztop = zstop - z0
              vtop = var_y(i,1,n) + (var_y(i,2,n)-var_y(i,1,n))
     >                 * (z0 + dztop/2. - z1)/(z2-z1)
              vbot = 0.
              dzbot = 0.

            else

              vbot = var_y(i,1,n)
              dzbot = tfac*wgz(2)
              z4 = z0 + tfac*sgz(ktop-1)
              z5 = z0 + tfac*wgz(ktop)
              z6 = z0 + tfac*sgz(ktop)
              dztop = zstop - z5
              if (zstop.eq.ztop) then
                vtop = var_y(i,ktop,n)
              else
                vtop = var_y(i,ktop,n) 
     >               + (var_y(i,ktop,n)-var_y(i,ktop-1,n))/(z6-z4)
     >                      * (zstop - dztop/2. - z6)
              end if

            end if

          else

            call bot_top(z0,wgz,ztop,zstrt,zstop,kbot,ktop,nz)

            if ((kbot.eq.1).and.(ktop.eq.1)) then

              z1 = z0 + tfac*sgz(1)
              z2 = z0 + tfac*sgz(2)
              dztop = zstop - zstrt
              vtop = var_y(i,1,n) + (var_y(i,2,n)-var_y(i,1,n))
     >                   *(zstop - dztop/2. - z1)/(z2-z1)
              vbot = 0.
              dzbot = 0.

            else if ((kbot.eq.nz).and.(ktop.eq.nz)) then

              z1 = z0 + tfac*sgz(nz-1)
              z2 = z0 + tfac*sgz(nz)
              dztop = zstop - zstrt
              vtop = var_y(i,nz,n) + (var_y(i,nz,n)-var_y(i,nz-1,n))
     >                   *(zstop - dztop/2. - z2)/(z2-z1)
              vbot = 0.
              dzbot = 0.

            else if (kbot.eq.ktop) then

              dztop = zstop - zstrt
              ztmp = zstop - dztop/2.
              z6 = z0 + tfac*sgz(ktop)
              if (ztmp.gt.z6) then
                z1 = z6
                z2 = z0 + tfac*sgz(ktop+1)
                vtop = var_y(i,ktop,n)
     >              + (var_y(i,ktop+1,n)-var_y(i,ktop,n))/(z2-z1)
     >                          * (ztmp-z1)
              else
                z1 = z0 + tfac*sgz(ktop-1)
                z2 = z6
                vtop = var_y(i,ktop,n)
     >              + (var_y(i,ktop,n)-var_y(i,ktop-1,n))/(z2-z1)
     >                          * (ztmp-z2)
              end if
              vbot = 0.
              dzbot = 0.

            else 

              z1 = z0 + tfac*sgz(kbot)
              z2 = z0 + tfac*wgz(kbot+1)
              z3 = z0 + tfac*sgz(kbot+1)
              z4 = z0 + tfac*sgz(ktop-1)
              z5 = z0 + tfac*wgz(ktop)
              z6 = z0 + tfac*sgz(ktop)
              dzbot = z2 - zstrt
              vbot = var_y(i,kbot,n)   
     >            + (var_y(i,kbot+1,n)-var_y(i,kbot,n))/(z3-z1)
     >                        * (zstrt + dzbot/2. - z1)
              dztop = zstop - z5
              if (zstop.eq.ztop) then
                vtop = var_y(i,ktop,n)
              else
                vtop = var_y(i,ktop,n)  
     >              + (var_y(i,ktop,n)-var_y(i,ktop-1,n))/(z6-z4)
     >                     * (zstop - dztop/2. - z6)
              end if

            end if

          end if

          area_int_y(n) = area_int_y(n) + vtop*dztop*dx
     >                        + vbot*dzbot*dx

          do k=kbot+1,ktop-1
            area_int_y(n) = area_int_y(n) + var_y(i,k,n)*dx
     >                       *tfac*(wgz(k+1)-wgz(k))
          enddo

        end if

      enddo

      enddo

      return
      end


        


