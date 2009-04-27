
#include "pointer.icl"

      ive_ptr function area_int_z(varnm,numvars,i1,i2,j1,j2,k1,
     >                ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names)
*
*  area_int_z calculates the area integral at constant z, normalized
*    by the horizontal area of the total domain, of the product
*    of the variables in array varnm(numvars).  Each of the variables
*    in varnm is interpolated to the constant z surface, and the 
*    product of the variables is then taken (i.e., the variables are
*    interpolated separately, rather than interpolating the product).
*    The domain of integration is the rectangular region at height
*    k*dz enclosed by i1*dx to i2*dx along the x-axis and j1*dy to
*    j2*dy along the y-axis.  Here dz is the average grid spacing
*    in the vertical (without terrain).
*    

      integer ndims,dims(4),i1,i2,j1,j2,k1,numvars
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dim_names(4),dunits,disunits     

      integer istrt,istop,jstrt,jstop,iu,iv,iw,k
      integer strbeg,strend,nz,nt,n,j
      real ztop,dom_min(4),dom_max(4),dx,dy,z,r_val,dz
      character*(80) tmp,tmp2,vunits,vdisunits
      logical err,new,errchck,units_err,twod
      ive_ptr getvar,getmem,zs,gz,wgz,sgz,var_z,var

      area_int_z = 0
      vunits = ' '
      vdisunits = ' '
      err = .false.

      call getrarr('plmin',dom_min,4,err)
      if (err) then
        write(*,*) 'Error:  could not find physical space domain ',
     >                     'minimums.'
        return
      end if
      call getrarr('plmax',dom_max,4,err)
      if (err) then
        write(*,*) 'Error:  could not find physical space domain ',
     >                      'maximums.'
        return
      end if

      if (dom_min(2).eq.dom_max(2)) then
        twod = .true.
      else 
        twod = .false.
      end if 

      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
      ztop = r_val(%val(wgz),dims(1))

      sgz = getvar('sgz',ndims,dims,stag,min,max,miss,dunits,
     >                             disunits,dim_names,new)
      nz = dims(1)
      
      dz = ztop/dims(1)

      zs = getvar('zbot_p',ndims,dims,stag,min,max,miss,dunits,
     >                              disunits,dim_names,new)

      dx = (max(1)-min(1))/(dims(1)-1)
      if (twod) then
        dy = 1. 
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      if (twod) then
        j1 = 0
        j2 = 1
      end if

      istrt = i1+1
      istop = i2
      jstrt = j1+1
      jstop = j2
      z = k1*dz

 87   FORMAT(10x,A5,F10.2,2x,A5,F10.2)
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
      write(*,87) 'x1 = ',i1*dx,'x2 = ',i2*dx
      write(*,87) 'y1 = ',j1*dy,'y2 = ',j2*dy
      write(*,88) z 

      if ((istrt.lt.1).or.(istop.gt.dims(1)).or.(jstrt.lt.1).or.
     >                (jstop.gt.dims(2))) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                          'boundaries'
        return
      else if ((istrt.gt.istop).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                          'integration.'
        return
      end if

      err = errchck(%val(zs),%val(sgz),ztop,z,dims(1),dims(2),nz)
      if (err) then
        write(*,*)  'Error:  specified value for z not valid. ',
     >           'z must exceed the minimum thermodynamic height '
        write(*,*)  'and must not exceed the highest thermodynamic',
     >                    ' height for all x,y.'
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
      if (dom_min(1).ge.min(1)) then
        iu = 1
      else if ((dom_min(2).ge.min(2)).and.(.not.twod)) then
        iv = 1
      else if (dom_min(3).ge.min(3)) then
        iw = 1
      end if

      dims(1) = dims(1) - iu
      dims(2) = dims(2) - iv

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

      do j=1,dims(1)*dims(2)*nt
        call s_val(%val(var_z),j,1.)
      enddo

      call v_at_z(%val(var_z),%val(var),%val(zs),%val(gz),z,ztop,
     >                 dims(1),dims(2),dims(3),nt,iu,iv,iw)

      vunits = dunits
      vdisunits = disunits

      k = 2

 658  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                              disunits,dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  cannot compute field ',
     >           varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          return
        else if (nt.ne.dims(4)) then
          write(*,*) 'Error:  fields ',
     >         varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                     ' and ',
     >         varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >              ' have different time domains.'
          return
        end if

        iu = 0
        iv = 0
        iw = 0
        if (dom_min(1).ge.min(1)) then
          iu = 1
        else if ((dom_min(2).ge.min(2)).and.(.not.twod)) then
          iv = 1
        else if (dom_min(3).ge.min(3)) then
          iw = 1
        end if

        dims(1) = dims(1)-iu
        dims(2) = dims(2)-iv
  
        if (iw.eq.1) then
          gz = wgz
        else
          gz = sgz
        end if

        call v_at_z(%val(var_z),%val(var),%val(zs),%val(gz),z,ztop,
     >                dims(1),dims(2),dims(3),dims(4),iu,iv,iw)

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

      area_int_z = getmem(nt)
      if (area_int_z.eq.0) then
        write(*,*) 'Error:  could not allocate space for area_int_z.'
        call freemem(var_z)
        return
      end if

      call area_int_z_calc(%val(area_int_z),%val(var_z),dx,dy,istrt,
     >                       istop,jstrt,jstop,dims(1),dims(2),nt)

      do k=1,3
        dims(k) = 1
        stag(k) = 0.
        min(k) = 0.
        max(k) = 0.
        dim_names(k) = 'one'
      enddo

      call units_smpfy(vunits,dunits,units_err)
      if (units_err) dunits = vunits
      call units_smpfy(vdisunits,disunits,units_err)
      if (units_err) disunits = vdisunits

      call freemem(var_z)

      return
      end 

********************************************************************

      subroutine v_at_z(var_z,var,zs,gz,z,ztop,nx,ny,nz,nt,
     >                        iu,iv,iw)
*
*  v_at_z interpolates variable var to the points defined by the 
*    intersections of the vertical lines through the thermodynamic
*    points with the horizontal plane defined by z.  We then take
*    the product of the interpolated values with whatever values 
*    already reside in var_z (so as to obtain the product of multiple
*    functions at the desired points). 
*     
      integer nx,ny,nz,nt,iu,iv,iw
      real var_z(nx,ny,nt),var(nx+iu,ny+iv,nz,nt),zs(nx,ny)
      real gz(nz),z,ztop
      
      integer k,kplus,i,j,n
      real z1,z2,var1,var2,tfac,znode

      do n=1,nt
      do j=1,ny
      do i=1,nx

        tfac = (ztop-zs(i,j))/ztop

        kplus = 0
        k = 1

 876    if ((kplus.eq.0).and.(k.le.nz)) then
          znode = zs(i,j) + tfac*gz(k)
          if (znode.gt.z) then
            kplus = k
          else 
            k = k+1
          end if
          goto 876
        end if

        if (kplus.eq.0) then
          kplus = nz 
        else if (kplus.eq.1) then
          write(*,*) 'Danger Will Robinson!  kplus = 1 in v_at_z!'
        end if

        z1 = zs(i,j) + tfac*gz(kplus-1)
        z2 = zs(i,j) + tfac*gz(kplus)
        var1 = (var(i,j,kplus-1,n)+var(i+iu,j+iv,kplus-1,n))/2.
        var2 = (var(i,j,kplus,n)+var(i+iu,j+iv,kplus,n))/2.

        var_z(i,j,n) = var_z(i,j,n)
     >             *(var1 + (var2-var1)*(z-z1)/(z2-z1))

      enddo
      enddo
      enddo

      return
      end

**********************************************************************

      subroutine area_int_z_calc(area_int,var_z,dx,dy,istrt,istop,
     >               jstrt,jstop,nx,ny,nt)

      integer nx,ny,nt,istrt,istop,jstrt,jstop
      real area_int(nt),var_z(nx,ny,nt),dx,dy

      integer i,j,n
      real dom_area

      dom_area = nx*dx*ny*dy

      do n=1,nt

        area_int(n) = 0.

        do j=jstrt,jstop
        do i=istrt,istop

          area_int(n) = area_int(n) + var_z(i,j,n)*dx*dy

        enddo
        enddo

        area_int(n) = area_int(n)/dom_area

      enddo

      return
      end

************************************************************************

      logical function errchck(zs,gz,ztop,z,nx,ny,nz)

      integer nx,ny,nz
      real zs(nx,ny),gz(nz),ztop,z,zmin,zmax,tfac

      integer i,j

      errchck = .false.

      do j=1,ny 
      do i=1,nx

        tfac = (ztop-zs(i,j))/ztop
        zmin = zs(i,j) + tfac*gz(1)
        zmax = zs(i,j) + tfac*gz(nz)

        if ((z.lt.zmin).or.(z.gt.zmax)) then
          errchck = .true.
          return
        end if

      enddo
      enddo

      return
      end

