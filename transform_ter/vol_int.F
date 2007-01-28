
#include "pointer.icl"
      ive_ptr function vol_int(varnm,numvars,i1,i2,j1,j2,k1,k2,ndims,
     >          dims,stag,min,max,missing,dunits,disunits,dim_names)

#include "default.icl"
*
*  vol_int computes the volume integral of variable varnm.  The domain of
*    integration is the rectangular region enclosed by i1*dx to i2*dx
*    along the x-axis, j1*dy to j2*dy along the y-axis, and k1*dz to
*    k2*dz along the z-axis.  Here dz is the average grid spacing in
*    the vertical (without terrain).
*
      integer ndims,dims(4),i1,i2,j1,j2,k1,k2,numvars
      real stag(4),min(4),max(4),missing
      character*(*) varnm(numvars),dunits,disunits,dim_names(4)

      integer istrt,istop,jstrt,jstop,iu,iv,iw,k,strbeg,strend,ijk
      real zstrt,zstop,r_val
      character*(80) tmp,tmp2
      logical err,new,units_err,requ
      ive_ptr getmem,var,var_t

      vol_int = 0

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if

      istrt = i1 + 1
      istop = i2 
      jstrt = j1 + 1
      jstop = j2 

      zstrt = k1*dz
      zstop = k2*dz

 101  FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 104  FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)
 103  FORMAT(10x,'z1 = ',F10.2,2x,'z2 = ',F10.2)

      tmp = varnm(1)
      do k=2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
      enddo
    
      write(*,*)
      write(*,*) 'Calculating volume integral of field ',
     >               tmp(strbeg(tmp):strend(tmp))
      write(*,*) 'over the domain defined by:  '
      write(*,101) plmin(1)+i1*dx, plmin(1)+i2*dx
      write(*,104) plmin(2)+j1*dy, plmin(2)+j2*dy
      write(*,103) zstrt, zstop
      write(*,*)

      if ((zstrt.lt.0.).or.(zstop.gt.wgz_max).or.(istrt.lt.1).or.
     >         (istop.gt.nx).or.(jstrt.lt.1).or.
     >               (jstop.gt.ny)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries' 
        return
      else if ((zstop.lt.zstrt).or.(istrt.gt.istop)
     >                       .or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                        'integration.'
        return
      end if 

      var = getvar(varnm(1),ndims,dims,stag,min,max,missing,dunits,
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

      k = 2

 873  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,missing,dunits,
     >                        disunits,dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  cannot access field ',
     >            varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
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

      vol_int = getmem(dims(4))
      if (vol_int.eq.0) then
        write(*,*) 'Error:  could not allocate space for vol_int.'
        call freemem(var_t)
        return
      end if

      call vol_int_calc(%val(vol_int),%val(var_t),%val(zbot_p),
     >         %val(sgz),%val(wgz),ztop,dx,dy,zstrt,zstop,istrt,
     >         istop,jstrt,jstop,dims(1),dims(2),dims(3),dims(4),err)

      if (err) then
        call freemem(var_t)
        call freemem(vol_int)
        vol_int = 0
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

      subroutine var_2_tpts(var_t,var,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real var_t(nx,ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        var_t(i,j,k,n) = var_t(i,j,k,n)
     >          *(var(i+iu,j+iv,k+iw,n)+var(i,j,k,n))/2.

      enddo
      enddo
      enddo
      enddo
  
      return
      end

*-------------------------------------------------------------------------

      subroutine vol_int_calc(vol_int,var,zs,sgz,wgz,ztop,dx,dy,
     >                         zstrt,zstop,istrt,istop,jstrt,jstop,
     >                         nx,ny,nz,nt,err)

      integer istrt,istop,jstrt,jstop,nx,ny,nz,nt
      real vol_int(nt),var(nx,ny,nz,nt),zs(nx,ny),sgz(nz),wgz(nz+1)
      real ztop,dx,dy,zstrt,zstop

      integer i,j,k,n,ktop,kbot
      real dztop,dzbot,z1,z2,z3,z4,z5,z6,tfac,vbot,vtop,z0,ztmp
      logical err

      do n=1,nt

      vol_int(n) = 0.

      do j=jstrt,jstop
      do i=istrt,istop

        z0 = zs(i,j)
        tfac = (ztop - z0)/ztop

        if (zstop.gt.z0) then

          if (zstrt.lt.z0) then

            call bot_top(z0,wgz,ztop,zstrt,zstop,kbot,ktop,nz)

            if (ktop.eq.1) then
 
              z1 = z0 + tfac*sgz(1)
              z2 = z0 + tfac*sgz(2)
              dztop = zstop - z0
              vtop = var(i,j,1,n) + (var(i,j,2,n)-var(i,j,1,n))
     >                 * (z0 + dztop/2. - z1)/(z2-z1)
              vbot = 0.
              dzbot = 0.

            else

              vbot = var(i,j,1,n)
              dzbot = tfac*wgz(2)
              z4 = z0 + tfac*sgz(ktop-1)
              z5 = z0 + tfac*wgz(ktop)
              z6 = z0 + tfac*sgz(ktop)
              dztop = zstop - z5
              if (zstop.eq.ztop) then
                vtop = var(i,j,ktop,n)
              else
                vtop = var(i,j,ktop,n) 
     >               + (var(i,j,ktop,n)-var(i,j,ktop-1,n))/(z6-z4)
     >                      * (zstop - dztop/2. - z6)
              end if

            end if

          else

            call bot_top(z0,wgz,ztop,zstrt,zstop,kbot,ktop,nz)

            if ((kbot.eq.1).and.(ktop.eq.1)) then

              z1 = z0 + tfac*sgz(1)
              z2 = z0 + tfac*sgz(2)
              dztop = zstop - zstrt
              vtop = var(i,j,1,n) + (var(i,j,2,n)-var(i,j,1,n))
     >                   *(zstop - dztop/2. - z1)/(z2-z1)
              vbot = 0.
              dzbot = 0.

            else if ((kbot.eq.nz).and.(ktop.eq.nz)) then

              z1 = z0 + tfac*sgz(nz-1)
              z2 = z0 + tfac*sgz(nz)
              dztop = zstop - zstrt
              vtop = var(i,j,nz,n) + (var(i,j,nz,n)-var(i,j,nz-1,n))
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
                vtop = var(i,j,ktop,n)
     >              + (var(i,j,ktop+1,n)-var(i,j,ktop,n))/(z2-z1)
     >                          * (ztmp-z1)
              else
                z1 = z0 + tfac*sgz(ktop-1)
                z2 = z6
                vtop = var(i,j,ktop,n)
     >              + (var(i,j,ktop,n)-var(i,j,ktop-1,n))/(z2-z1)
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
              vbot = var(i,j,kbot,n)   
     >            + (var(i,j,kbot+1,n)-var(i,j,kbot,n))/(z3-z1)
     >                        * (zstrt + dzbot/2. - z1)
              dztop = zstop - z5
              if (zstop.eq.ztop) then
                vtop = var(i,j,ktop,n)
              else
                vtop = var(i,j,ktop,n)  
     >              + (var(i,j,ktop,n)-var(i,j,ktop-1,n))/(z6-z4)
     >                     * (zstop - dztop/2. - z6)
              end if

            end if

          end if

          vol_int(n) = vol_int(n) + vtop*dztop*dx*dy
     >                        + vbot*dzbot*dx*dy

          do k=kbot+1,ktop-1
            vol_int(n) = vol_int(n) + var(i,j,k,n)*dx*dy
     >                       *tfac*(wgz(k+1)-wgz(k))
          enddo

        end if

      enddo
      enddo
 
      enddo

      return
      end

*-------------------------------------------------------------------------

      subroutine bot_top(zs,wgz,ztop,zstrt,zstop,kbot,ktop,nz)   
*
*    bot_top finds the grid elements intersected by the horizontal 
*    planes defined by zstrt and zstop
*
      integer kbot,ktop,nz
      real zs,wgz(nz+1),ztop,zstrt,zstop

      integer k
      real zgrd

      kbot = 0
      ktop = 0

      k = 1

 786  if ((kbot.eq.0).and.(k.le.nz)) then
        zgrd = zs + (ztop-zs)*wgz(k+1)/ztop
        if (zgrd.gt.zstrt) then
          kbot = k
        else
          k = k+1
        end if
        goto 786
      end if

 897  if ((ktop.eq.0).and.(k.le.nz)) then
        zgrd = zs + (ztop-zs)*wgz(k+1)/ztop
        if (zgrd.ge.zstop) then
          ktop = k
        else
          k = k+1
        end if
        goto 897
      end if

      if (ktop.eq.0) ktop = nz

      return
      end 
        


          
      
