#include "pointer.icl"
c 
      ive_ptr function kgain(obs_name,upd_name,R2,x0,y0,z0,mapon,
     &         ndims,dims,stag,min,max,miss,dunits,
     &         disunits,dim_names)

      use default_module

      implicit none

      real, intent(in) :: x0, y0, z0, R2
      character(len=*), intent(in) :: obs_name, upd_name 
      logical :: mapon

      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)

      integer :: ix, iy, iz, Ne, mnk
      integer :: iu, iv, iw, it
      real :: ir, jr
      character (len=80) :: outunits,outdisunits
      logical :: new
      ive_ptr :: xvar,xprior

      real, allocatable :: ye(:)

      ! Find location to interpolate to in the horizontal 
      if(mapon) then
        call lonlat_2_phys(ir,jr,x0,y0,1)
      else
        ir=x0 ; jr=y0
      endif

      ! Get first variable now.
      xvar = 0
      xvar = getvar(obs_name,ndims,dims,stag,min,max,miss,dunits,
     $                          disunits,dim_names,new)
      if (xvar.eq.0) then
        write(6,*) 'Error:  cannot access variable ',trim(obs_name)
        return
      end if
      ix = dims(1) ; iy = dims(2) ; iz = dims(3) ; Ne = dims(4)

      kgain = 0 ; kgain = getmem(ix*iy*iz)
      if(kgain.eq.0) then
        write(*,*) 'Error:  Cannot allocate space for kgain'
        return
      endif

      call stagering_calc(iu,iv,iw,it,dims,plmin,min,dx,dy,twod)

      allocate(ye(Ne))
      if(iu.eq.1) then
        call point_val(ye,%val(xvar),%val(zbot_u),%val(sgz),
     &               ztop,ir,jr,z0,ix,iy,iz,Ne)
      elseif(iv.eq.1) then
        call point_val(ye,%val(xvar),%val(zbot_v),%val(sgz),
     &               ztop,ir,jr,z0,ix,iy,iz,Ne)
      elseif(iw.eq.1) then
        call point_val(ye,%val(xvar),%val(zbot_p),%val(wgz),
     &               ztop,ir,jr,z0,ix,iy,iz,Ne)
      else
        call point_val(ye,%val(xvar),%val(zbot_p),%val(sgz),
     &               ztop,ir,jr,z0,ix,iy,iz,Ne)
      endif

      ! Get second variable now.
      xprior = 0
      xprior = getvar(upd_name,ndims,dims,stag,min,max,miss,dunits,
     $                          disunits,dim_names,new)
      if (xprior.eq.0) then
        write(6,*) 'Error:  cannot access variable ',trim(upd_name)
        return
      end if
      outunits = trim(dunits) ; outdisunits = trim(disunits)
      outunits = trim(dunits) ; outdisunits = trim(disunits)
      ix = dims(1) ; iy = dims(2) ; iz = dims(3) ; Ne = dims(4)

      mnk = ix*iy*iz
      call kgain_calc(%val(kgain),ye,R2,%val(xprior),mnk,Ne)

      dims(1) = ix ; dims(2) = iy ; dims(3) = iz ; dims(4) = 1
      dim_names(4) = 'one'
      dunits = outunits ; disunits = outdisunits

      return
      end

      subroutine kgain_calc(kgain,ye,xprior,R2,mnk,Ne)

      integer, intent(in) :: mnk,Ne
      real, intent(in) :: ye(Ne), xprior(mnk,Ne), R2
      real, intent(out) :: kgain(mnk)

      integer :: i,n
      real :: mean_fac,cov_fac,ye_var,cov,kdenom
      real :: ye_pert(Ne), xprior_pert(Ne)
      real :: ye_mean, xprior_mean

      mean_fac = 1.0/real(Ne)
      cov_fac = 1.0/real(Ne-1)

      kgain(:) = 0.0

      ye_mean = sum(ye(:))*mean_fac
      ye_pert(:) = ye(:) - ye_mean
      ye_var = sum(ye_pert(:)*ye_pert(:))*cov_fac
      kdenom = 1/(ye_var + R2)

      do i=1,mnk
        xprior_mean = sum(xprior(i,:))*mean_fac
        xprior_pert(:) = xprior(i,:) - xprior_mean

        cov = sum(ye_pert(:)*xprior_pert(:))*cov_fac
        kgain(i) = cov*kdenom
      enddo

      return
      end
        
