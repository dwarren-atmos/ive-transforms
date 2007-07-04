#include "pointer.icl"
c 
      ive_ptr function ens_cov(varnm1,varnm2,x0,y0,z0,mapon,
     &         ndims,dims,stag,min,max,miss,dunits,
     &         disunits,dim_names)

      use default_module

      implicit none

      real, intent(in) :: x0, y0, z0
      character(len=*), intent(in) :: varnm1, varnm2
      logical :: mapon

      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)

      integer :: ix, iy, iz, Ne, mnk
      integer :: iu, iv, iw, it
      real :: ir, jr
      character (len=80) :: outunits,outdisunits
      logical :: new
      ive_ptr :: var1,var2

      real, allocatable :: var1pt(:)

      ! Find location to interpolate to in the horizontal 
      if(mapon) then
        call lonlat_2_phys(ir,jr,x0,y0,1)
      else
        ir=x0 ; jr=y0
      endif

      ! Get first variable now.
      var1 = 0
      var1 = getvar(varnm1,ndims,dims,stag,min,max,miss,dunits,
     $                          disunits,dim_names,new)
      if (var1.eq.0) then
        write(6,*) 'Error:  cannot access variable ',trim(varnm1)
        return
      end if
      outunits = trim(dunits) ; outdisunits = trim(disunits)
      ix = dims(1) ; iy = dims(2) ; iz = dims(3) ; Ne = dims(4)

      ens_cov = 0 ; ens_cov = getmem(ix*iy*iz)
      if(ens_cov.eq.0) then
        write(*,*) 'Error:  Cannot allocate space for ens_cov'
        return
      endif

      call stagering_calc(iu,iv,iw,it,dims,plmin,min,dx,dy,twod)

      allocate(var1pt(Ne))
      if(iu.eq.1) then
        call point_val(var1pt,%val(var1),%val(zbot_u),%val(sgz),
     &               ztop,ir,jr,z0,ix,iy,iz,Ne)
      elseif(iv.eq.1) then
        call point_val(var1pt,%val(var1),%val(zbot_v),%val(sgz),
     &               ztop,ir,jr,z0,ix,iy,iz,Ne)
      elseif(iw.eq.1) then
        call point_val(var1pt,%val(var1),%val(zbot_p),%val(wgz),
     &               ztop,ir,jr,z0,ix,iy,iz,Ne)
      else
        call point_val(var1pt,%val(var1),%val(zbot_p),%val(sgz),
     &               ztop,ir,jr,z0,ix,iy,iz,Ne)
      endif

      ! Get second variable now.
      var2 = 0
      var2 = getvar(varnm2,ndims,dims,stag,min,max,miss,dunits,
     $                          disunits,dim_names,new)
      if (var2.eq.0) then
        write(6,*) 'Error:  cannot access variable ',trim(varnm2)
        return
      end if
      !outunits = trim(dunits) ; outdisunits = trim(disunits)
      ix = dims(1) ; iy = dims(2) ; iz = dims(3) ; Ne = dims(4)

      mnk = ix*iy*iz
      call ens_cov_calc(%val(ens_cov),var1pt,%val(var2),mnk,Ne)

      dims(1) = ix ; dims(2) = iy ; dims(3) = iz ; dims(4) = 1
      !dim_names(1) = 'one' ; dim_names(2) = 'one'
      !dim_names(3) = 'one' ; 
      dim_names(4) = 'one'
      !min(1) = 0. ;  max(1) = 0.
      !min(2) = 0. ;  max(2) = 0.  
      !min(3) = 0. ;  max(3) = 0.
      dunits = outunits ; disunits = outdisunits

      return
      end

      subroutine ens_cov_calc(ens_cov,fld1,fld2,mnk,Ne)

      integer, intent(in) :: mnk,Ne
      real, intent(in) :: fld1(Ne), fld2(mnk,Ne)
      real, intent(out) :: ens_cov(mnk)

      integer :: i,n
      real :: Ne_inv,var1,var2

      Ne_inv = 1.0/real(Ne)

      ens_cov(:) = 0.0

      do i=1,mnk
        var1 = 0.0 ; var2 = 0.0
        do n=1,Ne
          ens_cov(i) = ens_cov(i) + fld1(n)*fld2(i,n)
          var1 = var1 + fld1(n)**2
          var2 = var2 + fld2(i,n)**2
        enddo
        ens_cov(i) = ens_cov(i) / sqrt(var1*var2)
      enddo

      return
      end
        
