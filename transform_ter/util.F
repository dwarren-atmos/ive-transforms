
#include "pointer.icl"

      ive_ptr function v2z(fld,ndims,dims,stag,min,max,miss,dunits,
     >                             disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),fld

      integer strbeg,strend
      logical new,requ
      ive_ptr getmem,var,gz,zs

      v2z = 0

      var = getvar(fld,ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not access variable ',
     >                 fld(strbeg(fld):strend(fld))
        return
      end if

      v2z = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (v2z.eq.0) then
        write(6,*) 'Error:  could not allocate space for v2z.'
        return
      end if

      zs = zbot_p
      gz = sgz
      if (requ(min(1),plmin(1))) then
        zs = zbot_u
      else if (requ(min(2),plmin(2)).and.(twod.ne.1)) then
        zs = zbot_v
      else if (requ(min(3),plmin(3))) then
        gz = wgz
      end if

      miss = 7.77e5
      call v2z_calc(%val(v2z),%val(var),%val(zs),%val(gz),ztop,miss,
     >                    dims(1),dims(2),dims(3),dims(4))

      return
      end

*--------------------------------------------------------------------------

      subroutine v2z_calc(vz,var,zs,gz,ztop,miss,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real vz(nx,ny,nz,nt),var(nx,ny,nz,nt),zs(nx,ny),gz(nz)
      real ztop,miss
 
      integer i,j,k,n,kptr,kplus
      real z,z1,z2,z3,z4,tfac,znode,z0

      do n=1,nt
      do j=1,ny
      do i=1,nx

        z0 = zs(i,j)
        tfac = (ztop-z0)/ztop
        kptr = 1

        do k=1,nz
  
          z = gz(k)

          if (z.lt.z0) then

            vz(i,j,k,n) = miss

          else
 
            kplus = 0
 454        if ((kplus.eq.0).and.(kptr.le.nz)) then
              znode = z0 + tfac*gz(kptr)
              if (znode.ge.z) then
                kplus = kptr
              else 
                kptr = kptr+1
              end if
              goto 454
            end if

            if ((kplus.eq.1).or.(kplus.eq.2)) then
              z3 = z0 + tfac*gz(3)
              z2 = z0 + tfac*gz(2)
              z1 = z0 + tfac*gz(1)
              vz(i,j,k,n) = (z-z3)*(z-z2)/(z1-z3)/(z1-z2)
     >                              * var(i,j,1,n)
     >                  + (z-z3)*(z-z1)/(z2-z3)/(z2-z1)
     >                            * var(i,j,2,n)
     >                + (z-z2)*(z-z1)/(z3-z2)/(z3-z1)
     >                          * var(i,j,3,n)
            else if ((kplus.eq.0).or.(kplus.eq.nz)) then
              z3 = z0 + tfac*gz(nz)
              z2 = z0 + tfac*gz(nz-1)
              z1 = z0 + tfac*gz(nz-2)
              vz(i,j,k,n) = (z-z3)*(z-z2)/(z1-z3)/(z1-z2)
     >                              * var(i,j,nz-2,n)
     >                  + (z-z3)*(z-z1)/(z2-z3)/(z2-z1)
     >                            * var(i,j,nz-1,n)
     >                + (z-z2)*(z-z1)/(z3-z2)/(z3-z1)
     >                          * var(i,j,nz,n)
            else
              z4 = z0 + tfac*gz(kplus+1)
              z3 = z0 + tfac*gz(kplus)
              z2 = z0 + tfac*gz(kplus-1)
              z1 = z0 + tfac*gz(kplus-2)
              vz(i,j,k,n) = (z-z4)*(z-z3)*(z-z2)
     >              /(z1-z4)/(z1-z3)/(z1-z2)*var(i,j,kplus-2,n)
     >                    + (z-z4)*(z-z3)*(z-z1)
     >              /(z2-z4)/(z2-z3)/(z2-z1)*var(i,j,kplus-1,n)
     >                    + (z-z4)*(z-z2)*(z-z1)
     >              /(z3-z4)/(z3-z2)/(z3-z1)*var(i,j,kplus,n)
     >                    + (z-z3)*(z-z2)*(z-z1)
     >              /(z4-z3)/(z4-z2)/(z4-z1)*var(i,j,kplus+1,n)
            end if

          end if

        enddo

      enddo
      enddo
      enddo

      return
      end 

