
#include "pointer.icl"

      ive_ptr function area_int_y(varnm,numvars,i1,i2,k1,k2,j,ndims,
     >                   dims,stag,min,max,miss,dunits,disunits,
     >                   dim_names)
*
*  area_int_y calculates the area integral at constant y, normalized by
*    the horizontal area of the total domain, of the product of the 
*    variables in varnm(numvars).  The domain of integration is the
*    rectangular region at y = j*dy enclosed by i1*dx to i2*dx along
*    the x-axis and k1*dz to k2*dz along the z-axis.  Here dz is the
*    average grid spacing in the vertical (without terrain).
*
      integer numvars,i1,i2,k1,k2,j,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dim_names(4),dunits,disunits

      integer istrt,istop,iu,iv,iw,k,strbeg,strend,nt,jfix,twod
      real zstrt,zstop,ztop,dom_min(4),dom_max(4),dx,dy,dz,r_val,y
      character*(80) tmp,tmp2,vunits,vdisunits
      logical new,err,units_err
      ive_ptr getvar,getmem,zs,sgz,wgz,var_y,var

      area_int_y = 0
      vunits(1:len(vunits)) = ' '
      vdisunits(1:len(vdisunits)) = ' '
      err = .false.

      call getrarr('plmin',dom_min,4,err)
      if (err) then
        write(*,*) 'Error:  could not find physical space domain ',
     >                                    'minimums.'
        return
      end if
      call getrarr('plmax',dom_max,4,err)
      if (err) then
        write(*,*) 'Error:  could not find physical space domain ',
     >                                     'maximums.'
        return
      end if

      if (dom_min(2).eq.dom_max(2)) then
        twod = 1
      else 
        twod = 0
      end if

      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                        dim_names,new)
      ztop = r_val(%val(wgz),dims(1))

      sgz = getvar('sgz',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                        dim_names,new)
      dz = ztop/dims(1)

      zs = getvar('zbot_v',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)

      dx = (max(1)-min(1))/(dims(1)-1)
      if (twod.eq.1) then
        dy = 1.
        j = 0
      else      
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      istrt = i1+1
      istop = i2
      zstrt = k1*dz
      zstop = k2*dz
      y = j*dy
      jfix = j+1

 79   FORMAT(10x,A5,F10.2,2x,A5,F10.2)
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
      write(*,79) 'x1 = ',i1*dx,'x2 = ',i2*dx
      write(*,79) 'z1 = ',zstrt,'z2 = ',zstop
      write(*,90) y

      if (((j.lt.1).and.(.not.(twod.eq.1)))
     >           .or.((j.gt.(dims(2)-2)).and.(.not.(twod.eq.1))).or.
     >     (istrt.lt.1).or.(istop.gt.dims(1)).or.(zstrt.lt.0).or.
     >              (zstop.gt.ztop)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries.'
        return
      else if (((zstop-zstrt).le.(4.*dz)).or.(istrt.gt.istop)) then
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
          return
        else if (dims(4).ne.nt) then
          write(*,*) 'Error:  fields ',
     >        varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                        ' and ',
     >        varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >         'have different time dimensions.'
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

      call area_int_y_calc(%val(area_int_y),%val(var_y),%val(zs),
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

      call units_smpfy(vunits,dunits,units_err)
      if (units_err) dunits = vunits
      call units_smpfy(vdisunits,disunits,units_err)
      if (units_err) disunits = vdisunits

      call freemem(var_y)

      return
      end

**************************************************************************

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

**************************************************************************

      subroutine area_int_y_calc(area_int_y,var_y,zs,sgz,wgz,ztop, 
     >            dx,dy,zstrt,zstop,istrt,istop,jfix,nx,ny,nz,nt,err,
     >                                twod)

      integer istrt,istop,jfix,nx,ny,nz,nt,twod
      real area_int_y(nt),var_y(nx,nz,nt),zs(nx,ny+1-twod),sgz(nz),
     >        wgz(nz+1),dx,dy,zstrt,zstop,ztop
      logical err

      integer i,k,n,ktop,kbot
      real area,dztop,dzbot,z1,z2,z3,z4,z5,z6,tfac,vbot,vtop

      area = nx*dx*ny*dy

      do n=1,nt

        area_int_y(n) = 0.

        do i=istrt,istop 

          call bot_top(zs(i,jfix),wgz,ztop,zstrt,zstop,kbot,ktop,nz)

          tfac = (ztop - zs(i,jfix))/ztop

          z1 = zs(i,jfix) + tfac*sgz(kbot)
          z2 = zs(i,jfix) + tfac*wgz(kbot+1)
          z3 = zs(i,jfix) + tfac*sgz(kbot+1)
          z4 = zs(i,jfix) + tfac*sgz(ktop-1)
          z5 = zs(i,jfix) + tfac*wgz(ktop)
          z6 = zs(i,jfix) + tfac*sgz(ktop)

          dzbot = z2 - amax1(zs(i,jfix),zstrt)
         
          if (zs(i,jfix).ge.zstrt) then
            vbot = var_y(i,kbot,n)
          else
            vbot = var_y(i,kbot,n) +
     >              (var_y(i,kbot+1,n)-var_y(i,kbot,n))
     >            *(z2-dzbot/2.-z1)/(z3-z1)
          end if

          dztop = zstop - z5
 
          if (zstop.ge.ztop) then
            vtop = var_y(i,ktop,n)
          else
            vtop = var_y(i,ktop,n) +
     >              (var_y(i,ktop,n)-var_y(i,ktop-1,n))
     >            *(zstop-dztop/2.-z6)/(z6-z4)
          end if

          area_int_y(n) = area_int_y(n) +
     >            vtop*dztop*dx + vbot*dzbot*dx

          do k=kbot+1,ktop-1
            area_int_y(n) = area_int_y(n) +
     >            var_y(i,k,n)*dx*tfac*(wgz(k+1)-wgz(k))
          enddo

        enddo

        area_int_y(n) = area_int_y(n)/area

      enddo

      return
      end


        


