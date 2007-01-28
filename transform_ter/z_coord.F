

#include "pointer.icl"

      ive_ptr function z_coord(basnm,varnm,numvars,iflg,ndims,dims,
     >           stag,min,max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),numvars,iflg
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),basnm,varnm(numvars)

      integer strbeg,strend,iu,iv,iw,iub,ivb,iwb,ijk,nt,k
      character*(80) tmp,tmp2,vunits,vdisunits
      logical new,err,requ
      ive_ptr getmem,base,var,b0

      z_coord = 0

      nt = ntime

      tmp = varnm(1)
      do k=2,numvars
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        tmp = tmp2(strbeg(tmp2):strend(tmp2))//
     >      varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
      enddo

      write(6,*)
      write(6,*) 'Calculating field ', tmp(strbeg(tmp):strend(tmp))
      write(6,*) '   with vertical coordinate ',
     >                     basnm(strbeg(basnm):strend(basnm))
      write(6,*)

      base = getvar(basnm,ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names,new)
      if (base.eq.0) then
        write(6,*) 'Error:  cannot access field ',
     >                basnm(strbeg(basnm):strend(basnm))
        return
      end if

      iub = 0
      ivb = 0
      iwb = 0
      if (requ(plmin(1),min(1))) then
        iub = 1
      else if (requ(plmin(2),min(2)).and.(twod.ne.1)) then
        ivb = 1
      else if (requ(plmin(3),min(3))) then
        iwb = 1
      end if

      dims(1) = dims(1) - iub
      dims(2) = dims(2) - ivb
      dims(3) = dims(3) - iwb

      b0 = getmem(dims(3))
      if (b0.eq.0) then
        write(6,*) 'Error:  could not allocate space for b0.'
        return
      end if

      call get_b0(%val(b0),%val(base),dims(1),dims(2),dims(3),
     >                          dims(4),iub,ivb,iwb,err)
      if (err) then
        call freemem(b0)
        return
      end if

      var = getvar(varnm(1),ndims,dims,stag,min,max,miss,dunits,
     >               disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not access variable ',
     >               varnm(1)(strbeg(varnm(1)):strend(varnm(1)))
        call freemem(b0)
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

      z_coord = getmem(dims(1)*dims(2)*dims(3)*dims(4)) 
      if (z_coord.eq.0) then
        write(6,*) 'Error:  could not allocate space for z_coord.'
        call freemem(b0)
        return
      end if

      do ijk=1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(z_coord),ijk,1.)
      enddo

      miss = -0.007

      call var_2_bpts(%val(z_coord),%val(var),%val(base),%val(zbot_p),
     >                  %val(wgz),%val(sgz),%val(b0),ztop,miss,
     >                  dims(1),dims(2),dims(3),dims(4),iu,iv,iw,
     >                  iub,ivb,iwb,iflg)

      vunits = dunits
      vdisunits = disunits

      k = 2

 788  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new)
        if (var.eq.0) then
          write(6,*) 'Error:  could not access variable ',
     >          varnm(k)(strbeg(varnm(k)):strend(varnm(k)))
          call freemem(z_coord)
          call freemem(b0)
          return
        else if (nt.ne.dims(4)) then
          write(6,*) 'Error:  fields ',
     >          varnm(k-1)(strbeg(varnm(k-1)):strend(varnm(k-1))),
     >                         ' and ',
     >          varnm(k)(strbeg(varnm(k)):strend(varnm(k))),
     >                 ' have different time domains.'
          call freemem(z_coord)
          call freemem(b0)
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

        miss = -0.007

        call var_2_bpts(%val(z_coord),%val(var),%val(base),
     >            %val(zbot_p),%val(wgz),%val(sgz),%val(b0),ztop,
     >            miss,dims(1),dims(2),dims(3),dims(4),iu,iv,iw,
     >            iub,ivb,iwb,iflg)

        tmp = vunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >                dunits(strbeg(dunits):strend(dunits))

        tmp = vdisunits
        tmp2 = tmp(strbeg(tmp):strend(tmp))//'*'
        vdisunits = tmp2(strbeg(tmp2):strend(tmp2))//
     >         disunits(strbeg(disunits):strend(disunits))

        k = k + 1

        goto 788

      end if

      stag(1) = 0.5
      stag(2) = 0.5
      stag(3) = 0.5
      min(1) = plmin(1) + 0.5*dx
      min(2) = plmin(2) + 0.5*dy
      min(3) = -1.
      max(1) = plmax(1) - 0.5*dx
      max(2) = plmax(2) - 0.5*dy
      max(3) = sgz_max
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'
      dim_names(3) = 'nz'
 
      call units_smpfy(vunits,dunits,err)
      if (err) dunits = ' '

      call units_smpfy(vdisunits,disunits,err)
      if (err) disunits = ' '

      return
      end     

*------------------------------------------------------------------------

      subroutine get_b0(b0,base,nx,ny,nz,nt,iub,ivb,iwb,err)

      integer nx,ny,nz,nt,iub,ivb,iwb
      real b0(nz),base(nx+iub,ny+ivb,nz+iwb,nt)
      logical err

      real tst
      integer k

      err = .false.
     
      do k=1,nz
        b0(k) = 0.5*(base(1,1,k,1)+base(1+iub,1+ivb,k+iwb,1))
      enddo

      k = 1
 605  if ((err.eq..false.).and.(k.le.nz-1)) then
        tst = (b0(k+1)-b0(k))*(b0(nz)-b0(1))
        if (tst.le.0.) then
          write(6,*) 'Error:  initial base field is not monotonic.'
          err = .true.
        else
          k = k + 1
        end if
        goto 605
      end if

      return
      end 

*-------------------------------------------------------------------------

      subroutine var_2_bpts(v_b,var,base,zs_p,wgz,sgz,b0,ztop,miss,
     >                        nx,ny,nz,nt,iu,iv,iw,iub,ivb,iwb,iflg)

      integer nx,ny,nz,nt,iu,iv,iw,iub,ivb,iwb,iflg
      real v_b(nx,ny,nz,nt),base(nx+iub,ny+ivb,nz+iwb,nt),b0(nz)
      real var(nx+iu,ny+iv,nz+iw,nt),zs_p(nx,ny),wgz(nz+1),sgz(nz)
      real ztop,miss

      integer i,j,k,n,inc,bptr,bmin,vptr,vmin,kmin
      real zb,z1,z2,b1,b2,v1,v2,b_s,tfac,z0,z_node,v

      integer nzmax
      parameter (nzmax=201)
      real vgz(nzmax),bgz(nzmax)

      if (iwb.eq.0) then
        do k=1,nz
          bgz(k) = sgz(k)
        enddo
      else
        do k=1,nz+1
          bgz(k) = wgz(k)
        enddo
      end if

      if (iw.eq.0) then
        do k=1,nz
          vgz(k) = sgz(k)
        enddo
      else
        do k=1,nz+1
          vgz(k) = wgz(k)
        enddo
      end if

      if (b0(nz).gt.b0(1)) then
        inc = 1
      else
        inc = 0
      end if

      do n=1,nt
      do j=1,ny
      do i=1,nx

        z0 = zs_p(i,j)
        tfac = (ztop-z0)/ztop

        bptr = nz+iwb
        vptr = nz+iw

        do k=nz,1,-1

          if (inc.eq.1) then
     
            kmin = 0
 455        if ((kmin.eq.0).and.(bptr.ge.1)) then
              b_s = 0.5*(base(i,j,bptr,n)+base(i+iub,j+ivb,bptr,n))
              if (b0(k).gt.b_s) then
                kmin = bptr
              else
                bptr = bptr - 1
              end if
              goto 455
            end if

          else

            kmin = 0
 475        if ((kmin.eq.0).and.(bptr.ge.1)) then
              b_s = 0.5*(base(i,j,bptr,n)+base(i+iub,j+ivb,bptr,n))
              if (b0(k).lt.b_s) then
                kmin = bptr
              else
                bptr = bptr - 1
              end if
              goto 475
            end if

          end if

          if (kmin.eq.nz+iwb) then
            b1 = 0.5*(base(i,j,nz+iwb-1,n)
     >                     +base(i+iub,j+ivb,nz+iwb-1,n))
            b2 = 0.5*(base(i,j,nz+iwb,n)+base(i+iub,j+ivb,nz+iwb,n))
            if ((inc.eq.1).and.(b2.le.b1)) then
              zb = z0 - 1.
            else if ((inc.eq.0).and.(b1.le.b2)) then
              zb = z0 - 1.
            else
              z1 = z0 + tfac*bgz(nz+iwb-1)  
              z2 = z0 + tfac*bgz(nz+iwb)
              zb = z2 + (z2-z1)/(b2-b1)*(b0(k)-b2)
            end if
          else if (kmin.eq.0) then
            b1 = 0.5*(base(i,j,1,n)+base(i+iub,j+ivb,1,n))
            b2 = 0.5*(base(i,j,2,n)+base(i+iub,j+ivb,2,n))
            if ((inc.eq.1).and.(b2.le.b1)) then
              zb = z0 - 1.
            else if ((inc.eq.0).and.(b1.le.b2)) then
              zb = z0 - 1.
            else
              z1 = z0 + tfac*bgz(1)
              z2 = z0 + tfac*bgz(2)
              zb = z1 + (z2-z1)/(b2-b1)*(b0(k)-b1)
            end if
            if (iflg.eq.1) then
              if (zb.lt.z0) zb = z0
            end if
          else
            b1 = 0.5*(base(i,j,kmin,n)+base(i+iub,j+ivb,kmin,n))
            b2 = 0.5*(base(i,j,kmin+1,n)+base(i+iub,j+ivb,kmin+1,n))
            z1 = z0 + tfac*bgz(kmin)
            z2 = z0 + tfac*bgz(kmin+1)
            zb = z1 + (z2-z1)/(b2-b1)*(b0(k)-b1)
          end if

          if (zb.lt.z0-1.) then

            if (iflg.eq.2) then
              v_b(i,j,k,n) = 0.
            else 
              v_b(i,j,k,n) = miss
            end if
 
          else

            kmin = 0
 505        if ((kmin.eq.0).and.(vptr.ge.1)) then
              z_node = z0 + tfac*vgz(vptr)
              if (zb.gt.z_node) then
                kmin = vptr
              else
                vptr = vptr - 1
              end if
              goto 505
            end if

            if (kmin.eq.0) kmin = 1
  
            if (kmin.eq.nz+iw) then
              v1 = 0.5*(var(i,j,nz+iw-1,n)+var(i+iu,j+iv,nz+iw-1,n))
              v2 = 0.5*(var(i,j,nz+iw,n)+var(i+iu,j+iv,nz+iw,n))
              z1 = z0 + tfac*vgz(nz+iw-1)
              z2 = z0 + tfac*vgz(nz+iw)
              v_b(i,j,k,n) = v_b(i,j,k,n)*
     >                 (v2 + (v2-v1)/(z2-z1)*(zb-z2))
            else
              v1 = 0.5*(var(i,j,kmin,n)+var(i+iu,j+iv,kmin,n))
              v2 = 0.5*(var(i,j,kmin+1,n)+var(i+iu,j+iv,kmin+1,n))
              z1 = z0 + tfac*vgz(kmin)
              z2 = z0 + tfac*vgz(kmin+1)
              v_b(i,j,k,n) = v_b(i,j,k,n)*
     >                 (v1 + (v2-v1)/(z2-z1)*(zb-z1))
            end if

          end if

        enddo

      enddo
      enddo
      enddo

      return
      end

            
