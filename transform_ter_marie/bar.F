
#include "pointer.icl"

      ive_ptr function bar(varnm,numvars,i1,i2,j1,j2,ndims,dims,
     >             stag,min,max,miss,dunits,disunits,dim_names)
*
*  bar computes the area integral at constant z of the product of the
*    fields listed in varnm(numvars) over the region defined by i1*dx to 
*    i2*dx along the x-axis and j1*dy to j2*dy along the y-axis.  The 
*    integral is taken at the heights defined by the w-nodes in stretched
*    space (i.e., without the terrain), and the resulting function thus
*    depends on height as well as time.
*
      integer numvars,i1,i2,j1,j2,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),varnm(numvars)

      integer istrt,istop,jstrt,jstop,kmin,kmax,iu,iv,iw,k,n,j,
     >              strbeg,strend,nz,nt,twod
      real ztop,dom_min(4),dom_max(4),dx,dy,r_val,zmin,zmax,
     >              max_z,min_z
      character*(80) tmp,tmp2,vunits,vdisunits
      logical err,new,units_err
      ive_ptr getvar,getmem,zs,gz,wgz,sgz,var_z,var,wgz_tmp

      bar = 0
      vunits(1:len(vunits)) = ' '
      vdisunits(1:len(vdisunits)) = ' '
      err = .false.

      call getrarr('plmin',dom_min,4,err)
      if (err) then
        write(*,*)  'Error:  could not find physical space domain',
     >                               ' minimums.'
        return
      end if
      call getrarr('plmax',dom_max,4,err)
      if (err) then
        write(*,*)  'Error:  could not find physical space domain',
     >                               ' maximums.'
        return
      end if

      if (dom_min(2).eq.dom_max(2)) then
        twod = 1
      else 
        twod = 0
      end if

      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names,new)
      ztop = r_val(%val(wgz),dims(1))

      wgz_tmp = getmem(dims(1))
      do j=1,dims(1)
        call s_val(%val(wgz_tmp),j,r_val(%val(wgz),j))
      enddo

      sgz = getvar('sgz',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
      nz = dims(1)

      zs = getvar('zbot_p',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)

      dx = (max(1)-min(1))/(dims(1)-1)
      if (twod.eq.1) then
        j1 = 0
        j2 = 1
        dy = 1.
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      istrt = i1+1
      istop = i2
      jstrt = j1+1
      jstop = j2

 87   FORMAT(10x,A5,F10.2,2x,A5,F10.2)

      tmp = varnm(1)
      do j=2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >          varnm(j)(strbeg(varnm(j)):strend(varnm(j)))
      enddo

      write(*,*)
      write(*,*)  'Calculating the area integral of field ',
     >               tmp(strbeg(tmp):strend(tmp))
      write(*,*)  'over the domain defined by:  '
      write(*,87) 'x1 = ',i1*dx,'x2 = ',i2*dx
      write(*,87) 'y1 = ',j1*dy,'y2 = ',j2*dy

      if ((istrt.lt.0).or.(istop.gt.dims(1)).or.(jstrt.lt.0).or.
     >                    (jstop.gt.dims(2))) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                    'boundaries.'
        return
      else if ((istrt.gt.istop).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                    'integration.'
        return
      end if

      zmin = max_z(%val(zs),r_val(%val(sgz),1),ztop,dims(1),dims(2))
      zmax = min_z(%val(zs),r_val(%val(sgz),nz),ztop,dims(1),dims(2))

      kmin = 0
      k = 2
  
 345  if ((kmin.eq.0).and.(k.le.nz)) then
        if (r_val(%val(wgz),k).ge.zmin) then
          kmin = k
        else
          k = k+1
        end if
        goto 345
      end if

      if (kmin.eq.0) then
        write(*,*) 'Error:  could not find lowest integration ',
     >                      'surface.'
        return
      end if

      kmax = 0

 97   if ((kmax.eq.0).and.(k.le.nz+1)) then
        if (r_val(%val(wgz),k).gt.zmax) then
          kmax = k-1
        else
          k = k+1
        end if
        goto 97      
      end if

      if (kmax.eq.0) then
        write(*,*) 'Error:  could not find highest integration ',
     >                      'surface.'
        return
      else if (kmax.lt.kmin) then
        write(*,*) 'Error:  could not find any heights for ',
     >                               'integration.'
        return
      end if

      var = getvar(varnm(1),ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names,new)
      if (var.eq.0) then 
        write(*,*) 'Error:  cannot compute field ',
     >           varnm(1)(strbeg(varnm(1)):strend(varnm(1)))
        return
      end if

      iu = 0
      iv = 0
      iw = 0
      if (dom_min(1).ge.min(1)) then
        iu = 1
      else if (dom_min(2).ge.min(2)) then
        iv = 1 - twod
      else if (dom_min(3).ge.min(3)) then
        iw = 1
      end if

      dims(1) = dims(1)-iu
      dims(2) = dims(2)-iv
      dims(3) = dims(3)-iw
      nt = dims(4)

      if (iw.eq.1) then
        gz = wgz
      else
        gz = sgz
      end if

      var_z = getmem(dims(1)*dims(2)*dims(3)*dims(4))

      do j=1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(var_z),j,1.)
      enddo

      call v_z(%val(var_z),%val(var),%val(zs),%val(gz),
     >           %val(wgz_tmp),ztop,kmin,kmax,dims(1),dims(2),
     >           dims(3),dims(4),iu,iv,iw)

      vunits = dunits
      vdisunits = disunits

      k = 2
      
 934  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  could not compute field ',
     >         varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          return
        else if (nt.ne.dims(4)) then
          write(*,*) 'Error:  fields ',
     >           varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                          ' and ',
     >           varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >                 ' have different time domains.'
          return
        end if
        
        iu = 0
        iv = 0
        iw = 0
        if (dom_min(1).ge.min(1)) then
          iu = 1
        else if (dom_min(2).ge.min(2)) then
          iv = 1 - twod
        else if (dom_min(3).ge.min(3)) then
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

        call v_z(%val(var_z),%val(var),%val(zs),%val(gz),
     >           %val(wgz_tmp),ztop,kmin,kmax,dims(1),dims(2),
     >           dims(3),dims(4),iu,iv,iw)
 
        tmp = vunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >               dunits(strbeg(dunits):strend(dunits))
        tmp = vdisunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vdisunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >            disunits(strbeg(disunits):strend(disunits))

        k = k+1

        goto 934

      end if

      bar = getmem((dims(3)+1)*dims(4))

      miss = -55555555

      call bar_calc(%val(bar),%val(var_z),dx,dy,miss,istrt,istop,
     >            jstrt,jstop,kmin,kmax,dims(1),dims(2),dims(3),
     >            dims(4))

      do k=1,2
        dims(k) = 1
        min(k) = 0.
        max(k) = 0.
        stag(k) = 0.
        dim_names(k) = 'one'
      enddo

      dims(3) = dims(3)+1
      min(3) = r_val(%val(wgz),1)
      max(3) = r_val(%val(wgz),nz+1)
      stag(3) = 0.

      call units_smpfy(vunits,dunits,units_err)
      if (units_err) dunits = vunits
      call units_smpfy(vdisunits,disunits,units_err)
      if (units_err) disunits = vdisunits

      call freemem(wgz_tmp)
      call freemem(var_z)

      return
      end

**********************************************************************

      subroutine v_z(var_z,var,zs,gz,wgz,ztop,kmin,kmax,nx,ny,nz,
     >                          nt,iu,iv,iw)
*
*  v_z interpolates variable var to the horizontal planes defined by
*    wgz(kmin),wgz(kmin+1),...,wgz(kmax).  We then take the product of 
*    the interpolated values with whatever values already reside in
*    var_z, so as to obtain the product of multiple functions at the
*    desired heights.
*
      integer kmin,kmax,nx,ny,nz,nt,iu,iv,iw
      real var_z(nx,ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt)
      real zs(nx,ny),gz(nz+iw),wgz(nz+1),ztop

      integer i,j,k,n,kplus,ks
      real z1,z2,var1,var2,tfac,znode

      do n=1,nt
      do i=1,nx
      do j=1,ny

        tfac = (ztop-zs(i,j))/ztop
        ks = 1

        do k=kmin,kmax

          kplus = 0
 348      if ((kplus.eq.0).and.(ks.le.(nz+iw))) then
            znode = zs(i,j) + tfac*gz(ks)
            if (znode.gt.wgz(k)) then
              kplus = ks
            else
              ks = ks+1
            end if
            goto 348
          end if

          if (kplus.eq.0) kplus = nz

          z1 = zs(i,j) + tfac*gz(kplus-1)
          z2 = zs(i,j) + tfac*gz(kplus)
          var1 = (var(i+iu,j+iv,kplus-1,n) + var(i,j,kplus-1,n))/2.
          var2 = (var(i+iu,j+iv,kplus,n) + var(i,j,kplus,n))/2.
   
          var_z(i,j,k,n) = var_z(i,j,k,n) *
     >         (var1 + (var2-var1)/(z2-z1)*(wgz(k)-z1))

        enddo

      enddo
      enddo
      enddo

      return
      end

***************************************************************************

      subroutine bar_calc(bar,var_z,dx,dy,miss,istrt,istop,jstrt,
     >                        jstop,kmin,kmax,nx,ny,nz,nt)

      integer nx,ny,nz,nt,istrt,istop,jstrt,jstop,kmin,kmax
      real bar(nz+1,nt),var_z(nx,ny,nz,nt),dx,dy,miss

      integer i,j,k,n
      real area

      area = nx*dx*ny*dy

      do n=1,nt

        do k=kmin,kmax

          bar(k,n) = 0.

          do j=jstrt,jstop
          do i=istrt,istop
            bar(k,n) = bar(k,n) + var_z(i,j,k,n)
          enddo
          enddo

          bar(k,n) = bar(k,n)*dx*dy/area

        enddo

        do k=1,kmin-1
          bar(k,n) = miss
        enddo

        do k=kmax+1,nz+1
          bar(k,n) = miss
        enddo

      enddo

      return
      end

***************************************************************************

      real function max_z(zs,zeta,ztop,nx,ny)

      integer nx,ny
      real zs(nx,ny),zeta,ztop

      integer i,j
      real big_z,new_z

      big_z = zs(1,1) + (ztop-zs(1,1))*zeta/ztop

      do i=1,nx
      do j=1,ny

        new_z = zs(i,j) + (ztop-zs(i,j))*zeta/ztop
        if (new_z.gt.big_z) then
          big_z = new_z
        end if

      enddo
      enddo

      max_z = big_z
 
      return
      end 

***********************************************************************

      real function min_z(zs,zeta,ztop,nx,ny)

      integer nx,ny
      real zs(nx,ny),zeta,ztop

      integer i,j
      real small_z,new_z

      small_z = zs(1,1) + (ztop-zs(1,1))*zeta/ztop

      do i=1,nx
      do j=1,ny

        new_z = zs(i,j) + (ztop-zs(i,j))*zeta/ztop
        if (new_z.lt.small_z) then
          small_z = new_z
        end if

      enddo
      enddo

      min_z = small_z

      return
      end

********************************************************************** 

