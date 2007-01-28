
#include "pointer.icl"

      ive_ptr function area_int_z_vp(varnm,numvars,i1,i2,j1,j2,is,
     >                ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,is,numvars
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dim_names(4),dunits,disunits     

      integer istrt,istop,jstrt,jstop,iu,iv,iw,k
      integer strbeg,strend,nt,n,j
      real type,h,a1,a2,b1,b2,y_len,xc,yc
      character*(80) tmp,tmp2,vunits,vdisunits
      logical new,units_err,requ
      ive_ptr getmem,gz,var_gz,var_zs,var,area,var_c,zs_d

      area_int_z_vp = 0
      vunits = ' '
      vdisunits = ' '

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if

      istrt = i1 + 1
      istop = i2 
      jstrt = j1 + 1
      jstop = j2 

 86   FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 87   FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)

      tmp = varnm(1)
      do n = 2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >          varnm(n)(strbeg(varnm(n)):strend(varnm(n)))
      enddo

      write(*,*) 
      write(*,*) 'Calculating vertical profile of the area integral ',
     >               'of field ',tmp(strbeg(tmp):strend(tmp))
      write(*,*) 'over the domain defined by:  '
      write(*,86) plmin(1)+i1*dx, plmin(1)+i2*dx  
      write(*,87) plmin(2)+j1*dy, plmin(2)+j2*dy
      write(*,*)

      if ((istrt.lt.1).or.(istop.gt.nx).or.(jstrt.lt.1).or.
     >                (jstop.gt.ny)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                          'boundaries'
        return
      else if ((istrt.gt.istop).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                          'integration.'
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

      if (is.eq.0) then
        var_gz = getmem(dims(1)*dims(2)*(dims(3)+1)*nt)
        if (var_gz.eq.0) then
          write(*,*) 'Error:  could not allocate space for var_gz.'
          return
        end if
        do j=1,dims(1)*dims(2)*(dims(3)+1)*nt
          call s_val(%val(var_gz),j,1.)
        enddo
      else
        var_gz = getmem(dims(1)*dims(2)*dims(3)*nt)
        if (var_gz.eq.0) then
          write(*,*) 'Error:  could not allocate space for var_gz.'
          return
        end if
        do j=1,dims(1)*dims(2)*dims(3)*nt
          call s_val(%val(var_gz),j,1.)
        enddo
      end if

      var_zs = getmem(dims(1)*dims(2)*nt)
      if (var_zs.eq.0) then
        write(*,*) 'Error:  could not allocate space for var_zs.'
        return
      end if

      do j=1,dims(1)*dims(2)*nt
        call s_val(%val(var_zs),j,1.)
      enddo

      miss = -777.

      if (is.eq.0) then
        call v_2_wgz(%val(var_gz),%val(var),%val(zbot_p),%val(gz),
     >           %val(wgz),ztop,miss,dims(1),dims(2),dims(3),nt,
     >           iu,iv,iw)
      else
        call v_2_sgz(%val(var_gz),%val(var),%val(zbot_p),%val(gz),
     >           %val(sgz),ztop,miss,dims(1),dims(2),dims(3),nt,
     >           iu,iv,iw)
      end if
 
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
          call freemem(var_gz)
          call freemem(var_zs)
          return
        else if (nt.ne.dims(4)) then
          write(*,*) 'Error:  fields ',
     >         varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                     ' and ',
     >         varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >              ' have different time domains.'
          call freemem(var_gz)
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

        if (is.eq.0) then
          call v_2_wgz(%val(var_gz),%val(var),%val(zbot_p),%val(gz),
     >             %val(wgz),ztop,miss,dims(1),dims(2),dims(3),nt,
     >             iu,iv,iw)
        else
          call v_2_sgz(%val(var_gz),%val(var),%val(zbot_p),%val(gz),
     >             %val(sgz),ztop,miss,dims(1),dims(2),dims(3),nt,
     >             iu,iv,iw)
        end if

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
        call freemem(var_gz)
        call freemem(var_zs)
        return
      end if

      var_c = getmem(dims(1)*dims(2)*nt)
      if (var_c.eq.0) then
        write(*,*) 'Error:  could not allocate space for var_c.'
        call freemem(var_gz)
        call freemem(var_zs)
        call freemem(area)
        return
      end if 

      zs_d = getmem((2*dims(1)+1)*(2*dims(2)+1))
      if (zs_d.eq.0) then
        write(*,*) 'Error:  could not allocate space for zs_d.'
        call freemem(var_gz)
        call freemem(var_zs)
        call freemem(area)
        call freemem(var_c)
      end if

      call ter_params(type,h,a1,a2,b1,b2,y_len,xc,yc)

      call zs_dbl(%val(zs_d),dx,dy,plmin(1),plmin(2),type,h,
     >               a1,a2,b1,b2,y_len,xc,yc,dims(1),dims(2))

      if (is.eq.0) then
        area_int_z_vp = getmem((dims(3)+1)*nt)
      else
        area_int_z_vp = getmem(dims(3)*nt)
      end if
      if (area_int_z_vp.eq.0) then
        write(*,*) 'Error:  could not allocate space for ',
     >                   'area_int_z_vp.'
        call freemem(var_gz)
        call freemem(var_zs)
        call freemem(zs_d)
        call freemem(area)
        call freemem(var_c)
        return
      end if
    
      call area_int_z_vp_calc(%val(area_int_z_vp),%val(var_gz),
     >              %val(area),%val(var_c),%val(var_zs),%val(zs_d),
     >              %val(zbot_p),%val(wgz),%val(sgz),miss,dx,dy,
     >              istrt,istop,jstrt,jstop,dims(1),dims(2),dims(3),
     >              nt,twod,1-is)

      do k=1,2
        dims(k) = 1
        stag(k) = 0.
        min(k) = 0.
        max(k) = 0.
        dim_names(k) = 'one'
      enddo
 
      if (is.eq.1) then
        stag(3) = 0.5
        min(3) = sgz_min
        max(3) = sgz_max
        dim_names(3) = 'nz'
      else
        dims(3) = dims(3) + 1
        stag(3) = 0.
        min(3) = wgz_min
        max(3) = wgz_max
        dim_names(3) = 'nzp1'
      end if

      vunits = vunits(strbeg(vunits):strend(vunits))//'*m*m'
      vdisunits = vdisunits(strbeg(vdisunits):strend(vdisunits))
     >              //'*m*m'

      call units_smpfy(vunits,dunits,units_err)
      if (units_err) dunits = ' '
      call units_smpfy(vdisunits,disunits,units_err)
      if (units_err) disunits = ' '

      call freemem(var_gz)
      call freemem(area)
      call freemem(var_c)
      call freemem(zs_d)
      call freemem(var_zs)

      return
      end 

*-----------------------------------------------------------------------

      subroutine v_2_wgz(v_wgz,var,zs,gz,wgz,ztop,miss,nx,ny,nz,nt,
     >                        iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real v_wgz(nx,ny,nz+1,nt),var(nx+iu,ny+iv,nz+iw,nt),zs(nx,ny)
      real gz(nz+iw),wgz(nz+1),ztop,miss

      integer i,j,k,n,kplus,ks
      real z1,z2,var1,var2,tfac,znode,z0

      do n=1,nt
      do j=1,ny
      do i=1,nx

        z0 = zs(i,j)
        tfac = (ztop-z0)/ztop

        ks = 1

        do k=1,nz+1

          if (wgz(k).lt.z0) then
 
            v_wgz(i,j,k,n) = miss

          else

            kplus = 0
         
 338        if ((kplus.eq.0).and.(ks.le.(nz+iw))) then
              znode = z0 + tfac*gz(ks) 
              if (znode.gt.wgz(k)) then
                kplus = ks
              else
                ks = ks+1
              end if
              goto 338
            end if

            if (kplus.eq.0) kplus = nz+iw

            if (kplus.eq.1) then
              z1 = z0 + tfac*gz(1)
              z2 = z0 + tfac*gz(2)
              var1 = 0.5*(var(i,j,1,n)+var(i+iu,j+iv,1,n))
              var2 = 0.5*(var(i,j,2,n)+var(i+iu,j+iv,2,n))
              v_wgz(i,j,k,n) = v_wgz(i,j,k,n) *
     >            (var1 + (var2-var1)/(z2-z1)*(wgz(k)-z1))
            else
              z1 = z0 + tfac*gz(kplus-1)
              z2 = z0 + tfac*gz(kplus)
              var1 = 0.5*(var(i,j,kplus-1,n)+var(i+iu,j+iv,kplus-1,n))
              var2 = 0.5*(var(i,j,kplus,n)+var(i+iu,j+iv,kplus,n))
              v_wgz(i,j,k,n) = v_wgz(i,j,k,n) *
     >              (var2 + (var2-var1)/(z2-z1)*(wgz(k)-z2))
            end if

          end if

        enddo

      enddo
      enddo
      enddo

      return
      end

*-------------------------------------------------------------------------

      subroutine v_2_sgz(v_sgz,var,zs,gz,sgz,ztop,miss,nx,ny,nz,nt,
     >                           iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real v_sgz(nx,ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt),zs(nx,ny)
      real gz(nz+iw),sgz(nz),ztop,miss

      integer i,j,k,n,kptr,kplus
      real z1,z2,var1,var2,tfac,znode,z0

      do n=1,nt

      do j=1,ny
      do i=1,nx

        z0 = zs(i,j)
        tfac = (ztop - z0)/ztop
        kptr = 1

        do k=1,nz

          if (sgz(k).lt.z0) then

            v_sgz(i,j,k,n) = miss

          else
 
            kplus = 0

 454        if ((kplus.eq.0).and.(kptr.le.nz+iw)) then
              znode = z0 + tfac*gz(kptr)
              if (znode.ge.sgz(k)) then
                kplus = kptr
              else 
                kptr = kptr+1
              end if
              goto 454
            end if

            if (kplus.eq.0) kplus = nz

            if (kplus.eq.1) then
              z1 = z0 + tfac*gz(1)
              z2 = z0 + tfac*gz(2)
              var1 = (var(i,j,1,n)+var(i+iu,j+iv,1,n))/2.
              var2 = (var(i,j,2,n)+var(i+iu,j+iv,2,n))/2.
              v_sgz(i,j,k,n) = v_sgz(i,j,k,n) * ( var1 +
     >              (var2-var1)/(z2-z1)*(sgz(k)-z1) )
            else
              z1 = z0 + tfac*gz(kplus-1)
              z2 = z0 + tfac*gz(kplus)
              var1 = (var(i,j,kplus-1,n)+var(i+iu,j+iv,kplus-1,n))/2.
              var2 = (var(i,j,kplus,n)+var(i+iu,j+iv,kplus,n))/2.
              v_sgz(i,j,k,n) = v_sgz(i,j,k,n) * ( var2 +
     >              (var2-var1)/(z2-z1)*(sgz(k)-z2) )
            end if

          end if

        enddo

      enddo
      enddo

      enddo

      return
      end

*---------------------------------------------------------------------

      subroutine area_int_z_vp_calc(area_int,v_gz,area,v_c,v_zs,zs_d,
     >              zs_p,wgz,sgz,miss,dx,dy,istrt,istop,jstrt,jstop,
     >              nx,ny,nz,nt,twod,iw)

      integer nx,ny,nz,nt,istrt,istop,jstrt,jstop,twod,iw
      real area_int(nz+iw,nt),v_gz(nx,ny,nz+iw,nt),area(nx,ny)
      real v_c(nx,ny,nt),v_zs(nx,ny,nt),zs_d(2*nx+1,2*ny+1)
      real zs_p(nx,ny),wgz(nz+1),sgz(nz),miss,dx,dy

      integer k,n
      real sum,tot_area

      do k=1,nz+iw

        if (iw.eq.1) then
          call int_z_setup(area,v_c,zs_d,v_zs,wgz(k),dx,dy,nx,ny,nt)
        else
          call int_z_setup(area,v_c,zs_d,v_zs,sgz(k),dx,dy,nx,ny,nt)
        end if

        do n=1,nt
        
          if (twod.eq.1) then
            write(*,*) 'Error:  int_z_1d not yet implemented.'
          else
            if (iw.eq.1) then
              call int_z_2d(sum,tot_area,v_gz(1,1,k,n),area,
     >                v_c(1,1,n),zs_p,wgz(k),istrt,istop,jstrt,jstop,
     >                nx,ny)
            else
              call int_z_2d(sum,tot_area,v_gz(1,1,k,n),area,
     >                v_c(1,1,n),zs_p,sgz(k),istrt,istop,jstrt,jstop,
     >                nx,ny)
            end if
          end if 

          if (tot_area.ne.0.) then
            area_int(k,n) = sum
          else
            area_int(k,n) = miss
          end if

        enddo

      enddo
 
      return
      end

