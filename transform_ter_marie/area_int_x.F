
#include "pointer.icl"

      ive_ptr function area_int_x(varnm,numvars,j1,j2,k1,k2,i,ndims,
     >                   dims,stag,min,max,miss,dunits,disunits,
     >                   dim_names)
*
*  area_int_x calculates the area integral at constant x, normalized by
*    the horizontal area of the total domain, of the product of the 
*    variables in varnm(numvars).  The domain of integration is the
*    rectangular region at x = i*dx enclosed by j1*dy to j2*dy along
*    the y-axis and k1*dz to k2*dz along the z-axis.  Here dz is the
*    average grid spacing in the vertical (without terrain).
*
      integer numvars,j1,j2,k1,k2,i,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dim_names(4),dunits,disunits

      integer jstrt,jstop,iu,iv,iw,k,strbeg,strend,nt,ifix
      real zstrt,zstop,ztop,dom_min(4),dom_max(4),dx,dy,dz,r_val,x
      character*(80) tmp,tmp2,vunits,vdisunits
      logical new,err,units_err,twod
      ive_ptr getvar,getmem,zs,sgz,wgz,var_x,var

      area_int_x = 0
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
        twod = .true.
      else 
        twod = .false.
      end if

      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                        dim_names,new)
      ztop = r_val(%val(wgz),dims(1))

      sgz = getvar('sgz',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                        dim_names,new)
      dz = ztop/dims(1)

      zs = getvar('zbot_u',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)

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

      jstrt = j1+1
      jstop = j2
      zstrt = k1*dz
      zstop = k2*dz
      x = i*dx
      ifix = i+1

 79   FORMAT(10x,A5,F10.2,2x,A5,F10.2)
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
      write(*,79) 'y1 = ',j1*dy,'y2 = ',j2*dy
      write(*,79) 'z1 = ',zstrt,'z2 = ',zstop
      write(*,90) x

      if ((i.lt.1).or.(i.gt.(dims(1)-2)).or.(jstrt.lt.1).or.
     >          (jstop.gt.dims(2)).or.(zstrt.lt.0).or.
     >              (zstop.gt.ztop)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries.'
        return
      else if (((zstop-zstrt).le.(4.*dz)).or.(jstrt.gt.jstop)) then
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
      else if ((dom_min(2).ge.min(2)).and.(.not.twod)) then
        iv = 1
      else if (dom_min(3).ge.min(3)) then
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
        else if ((dom_min(2).ge.min(2)).and.(.not.twod)) then
          iv = 1
        else if (dom_min(3).ge.min(3)) then
          iw = 1
        end if

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

      area_int_x = getmem(nt)

      call area_int_x_calc(%val(area_int_x),%val(var_x),%val(zs),
     >              %val(sgz),%val(wgz),ztop,dx,dy,zstrt,zstop,jstrt,
     >              jstop,ifix,dims(1),dims(2),dims(3),dims(4),err)

      if (err) then
        call freemem(var_x)
        call freemem(area_int_x)
        area_int_x = 0
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

      call freemem(var_x)

      return
      end

**************************************************************************

      subroutine v_at_x(var_x,var,ifix,nx,ny,nz,nt,iu,iv,iw)

      integer iu,iv,iw,nx,ny,nz,nt,ifix
      real var_x(ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt)
      integer n,k,j

      do n=1,nt
      do k=1,nz
      do j=1,ny

        var_x(j,k,n) = var_x(j,k,n) *
     >        (var(ifix,j,k,n) + var(ifix,j+iv,k+iw,n)
     >     + var(ifix+iu-1,j,k,n) + var(ifix+iu-1,j+iv,k+iw,n))/4.
     
      enddo
      enddo
      enddo
 
      return
      end

**************************************************************************

      subroutine area_int_x_calc(area_int_x,var_x,zs,sgz,wgz,ztop, 
     >            dx,dy,zstrt,zstop,jstrt,jstop,ifix,nx,ny,nz,nt,err)

      integer jstrt,jstop,ifix,nx,ny,nz,nt
      real area_int_x(nt),var_x(ny,nz,nt),zs(nx+1,ny),sgz(nz),
     >        wgz(nz+1),dx,dy,zstrt,zstop,ztop
      logical err

      integer j,k,n,ktop,kbot
      real area,dztop,dzbot,z1,z2,z3,z4,z5,z6,tfac,vbot,vtop

      area = nx*dx*ny*dy

      do n=1,nt

        area_int_x(n) = 0.

        do j=jstrt,jstop 

          call bot_top(zs(ifix,j),wgz,ztop,zstrt,zstop,kbot,ktop,nz)

          tfac = (ztop - zs(ifix,j))/ztop

          z1 = zs(ifix,j) + tfac*sgz(kbot)
          z2 = zs(ifix,j) + tfac*wgz(kbot+1)
          z3 = zs(ifix,j) + tfac*sgz(kbot+1)
          z4 = zs(ifix,j) + tfac*sgz(ktop-1)
          z5 = zs(ifix,j) + tfac*wgz(ktop)
          z6 = zs(ifix,j) + tfac*sgz(ktop)

          dzbot = z2 - amax1(zs(ifix,j),zstrt)
         
          if (zs(ifix,j).ge.zstrt) then
            vbot = var_x(j,kbot,n)
          else
            vbot = var_x(j,kbot,n) +
     >              (var_x(j,kbot+1,n)-var_x(j,kbot,n))
     >            *(z2-dzbot/2.-z1)/(z3-z1)
          end if

          dztop = zstop - z5
 
          if (zstop.ge.ztop) then
            vtop = var_x(j,ktop,n)
          else
            vtop = var_x(j,ktop,n) +
     >              (var_x(j,ktop,n)-var_x(j,ktop-1,n))
     >            *(zstop-dztop/2.-z6)/(z6-z4)
          end if

          area_int_x(n) = area_int_x(n) +
     >            vtop*dztop*dy + vbot*dzbot*dy

          do k=kbot+1,ktop-1
            area_int_x(n) = area_int_x(n) +
     >            var_x(j,k,n)*dy*tfac*(wgz(k+1)-wgz(k))
          enddo

        enddo

        area_int_x(n) = area_int_x(n)/area

      enddo

      return
      end


        


