
#include "pointer.icl"

      ive_ptr function per_vol(varnm,n,whch,ndims,dims,stag,min,
     >             max,miss,dunits,disunits,dim_names)
*
*  per_vol converts a per unit mass quantity to a per unit volume
*    quantity.  Specifically, it returns the product of variable
*    varnm and the mass density at the thermodynamic pts.
*    NOTE:  we assume that varnm is defined at the u,v,w or
*           thermodynamic pts
*
      integer ndims,dims(4),n
      real stag(4),min(4),max(4),miss
      character*(*) varnm,dunits,disunits,dim_names(4),whch

      integer iu,iv,iw,k,strend,strbeg
      real rmin(3),rmax(3),rstag(3)
      character*(80) tmp
      logical new,err
      ive_ptr getvar,getmem,rho,var,ref

      per_vol = 0
      
      if (whch(1:4).eq.'base') then
        rho = getvar('RHO0',ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new)
      else if (whch(1:3).eq.'ref') then
        rho = ref('RHO',n,ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names)
      else
        rho = getvar('RHO',ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new)
      end if
      if (rho.eq.0) then
        write(*,*) 'Error:  cannot calculate appropriate mass', 
     >                       ' density.'
        return
      end if

      do k=1,3
        rmin(k) = min(k)
        rmax(k) = max(k)
        rstag(k) = stag(k)
      enddo

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot get variable ',
     >                 varnm(strbeg(varnm):strend(varnm))
        return
      end if

      iu = 0
      iv = 0
      iw = 0
      if (min(1).lt.rmin(1)) then
        iu = 1
      else if (min(2).lt.rmin(2)) then
        iv = 1
      else if (min(3).lt.rmin(3)) then
        iw = 1
      end if

      dims(1) = dims(1) - iu
      dims(2) = dims(2) - iv
      dims(3) = dims(3) - iw

      do k=1,3
        min(k) = rmin(k)
        max(k) = rmax(k)
        stag(k) = rstag(k)
      enddo

      tmp = dunits(strbeg(dunits):strend(dunits))//'*kg/m/m/m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = tmp(strbeg(tmp):strend(tmp))
      tmp = disunits(strbeg(disunits):strend(disunits))//'*kg/m/m/m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = tmp(strbeg(tmp):strend(tmp))

      per_vol = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (per_vol.eq.0) then
        write(*,*) 'Error:  could not allocate space for per_vol.'
        return
      end if

      if ((whch(1:4).eq.'base').or.(whch(1:3).eq.'ref')) then
        call per_vol_calc_base(%val(per_vol),%val(var),%val(rho),
     >             dims(1),dims(2),dims(3),dims(4),iu,iv,iw)
      else 
        call per_vol_calc(%val(per_vol),%val(var),%val(rho),dims(1),
     >             dims(2),dims(3),dims(4),iu,iv,iw)
      end if

      if (whch(1:3).eq.'ref') then
        call freemem(rho)
      end if

      return 
      end
       
*************************************************************************

      subroutine per_vol_calc(per_vol,var,rho,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real per_vol(nx,ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt),
     >                 rho(nx,ny,nz,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        per_vol(i,j,k,n) = (var(i+iu,j+iv,k+iw,n)+var(i,j,k,n))/2.
     >                        * rho(i,j,k,n)

      enddo
      enddo
      enddo
      enddo

      return
      end

*************************************************************************

      subroutine per_vol_calc_base(per_vol,var,r0,nx,ny,nz,nt,
     >                                  iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real per_vol(nx,ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt),
     >                  r0(nx,ny,nz)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        per_vol(i,j,k,n) = (var(i+iu,j+iv,k+iw,n)+var(i,j,k,n))/2.
     >                      * r0(i,j,k)

      enddo
      enddo
      enddo
      enddo

      return
      end

