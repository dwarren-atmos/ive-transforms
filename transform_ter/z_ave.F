
#include "pointer.icl"

      ive_ptr function z_ave(fld,z1,z2,ndims,dims,stag,min,max,
     >                    miss,dunits,disunits,dim_names)
*
*  z_ave calculates the volume average of fld over z as a function of
*    x,y and t.
*

#include "default.icl"

      integer ndims,dims(4)
      real miss,stag(4),min(4),max(4),z1,z2
      character*(*) fld,dunits,disunits,dim_names(4)

      integer strbeg,strend,iu,iv,iw,istrt,istop,jstrt,jstop,ijk
      real r_val
      logical err,new,requ
      ive_ptr getmem,var_t,var

      z_ave = 0

      var = getvar(fld,ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new) 
      if (var.eq.0) then
        write(*,*) 'Error:  could not compute field ',
     >                 fld(strbeg(fld):strend(fld))
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
        write(*,*) 'Error:  unable to allocate space for var_t.'
        return
      end if
      do ijk=1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(var_t),ijk,1.)
      enddo

      call var_2_tpts(%val(var_t),%val(var),dims(1),dims(2),dims(3),
     >                              dims(4),iu,iv,iw)
      
      istrt = 1
      istop = dims(1)
      jstrt = 1
      jstop = dims(2) 

 101  FORMAT(10x,'z1 = ',F10.2,2x,'z2 = ',F10.2)

      write(*,*)
      write(*,*) 'Calculating volume average over z of field ',
     >                     fld(strbeg(fld):strend(fld))
      write(*,*) 'over the domain defined by:  '
      write(*,101) z1,z2
      write(*,*)

      if ((z1.lt.0.).or.(z2.gt.ztop)) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                  'boundaries.'
        call freemem(var_t)
        return
      else if (z2.le.z1) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                       'integration.' 
        call freemem(var_t)
        return
      end if

      z_ave = getmem(dims(1)*dims(2)*dims(4))
      if (z_ave.eq.0) then
        write(*,*) 'Error:  unable to allocate space for z_ave.'
        call freemem(var_t)
      end if

      miss = -777

      call z_ave_calc(%val(z_ave),%val(var_t),%val(zbot_p),%val(sgz),
     >             %val(wgz),ztop,dx,dy,z1,z2,miss,istrt,istop,
     >             jstrt,jstop,dims(1),dims(2),dims(3),dims(4))

      dims(3) = 1
      min(3) = 0.
      max(3) = 0.
      stag(3) = 0.
      dim_names(3) = 'one'

      min(1) = plmin(1) + 0.5*dx
      min(2) = plmin(2) + 0.5*dy
      max(1) = plmax(1) - 0.5*dx
      max(2) = plmax(2) - 0.5*dy
      stag(1) = 0.5
      stag(2) = 0.5
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'

      call freemem(var_t)

      return
      end

*---------------------------------------------------------------------

      subroutine z_ave_calc(z_ave,var,zs,sgz,wgz,ztop,dx,dy,zstrt,
     >      zstop,miss,istrt,istop,jstrt,jstop,nx,ny,nz,nt)

      integer istrt,istop,jstrt,jstop,nx,ny,nz,nt
      real  z_ave(nx,ny,nt),var(nx,ny,nz,nt),zs(nx,ny)
      real  sgz(nz),wgz(nz+1),ztop,dx,dy,zstrt,zstop,miss

      integer i,j,n,kbot,ktop,k
      real tfac,dzbot,dztop,vbot,vtop,z1,z2,z3,z4,z5,z6,volume,
     >                  z0,ztmp
      
      do n=1,nt

      do j=jstrt,jstop
      do i=istrt,istop

        z_ave(i,j,n) = 0.

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

          z_ave(i,j,n) = z_ave(i,j,n) + vtop*dztop*dx*dy
     >                        + vbot*dzbot*dx*dy

          do k=kbot+1,ktop-1
            z_ave(i,j,n) = z_ave(i,j,n) + var(i,j,k,n)*dx*dy
     >                       *tfac*(wgz(k+1)-wgz(k))
          enddo

          volume = (zstop - amax1(zstrt,z0))*dx*dy

          z_ave(i,j,n) = z_ave(i,j,n)/volume

        else

          z_ave(i,j,n) = miss

        end if

      enddo
      enddo
 
      enddo

      return
      end
