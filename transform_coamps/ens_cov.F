#include "pointer.icl"
c 
      ive_ptr function ens_cov(varnm1,varnm2,x0,y0,z0,mapon,
     &         ndims,dims,stag,min,max,miss,dunits,
     &         disunits,dim_names)

      implicit none

#include "default.icl"

      real, intent(in) :: x0, y0, z0
      character(len=*), intent(in) :: varnm1, varnm2
      logical :: mapon

      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)

      integer :: ix, iy, iz, Ne
      real :: ir, jr
      character (len=80) :: outunits,outdisunits
      logical :: new
      ive_ptr :: var1,var2

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


      ens_cov = 0 ; ens_cov = getmem(Ne)
      !ens_cov = getmem(ix*iy*iz)
      if(ens_cov.eq.0) then
        write(*,*) 'Error:  Cannot allocate space for ens_cov'
        return
      endif


      call point_val(%val(ens_cov),%val(var1),%val(zbot_u),%val(sgz),
     &               ztop,ir,jr,z0,ix,iy,iz,Ne)



      dims(1) = 1 ; dims(2) = 1 ; dims(3) = 1 ; dims(4) = Ne
      dim_names(1) = 'one' ; dim_names(2) = 'one'
      dim_names(3) = 'one' ; dim_names(4) = 'time'
      min(3) = z0 ;  max(3) = z0
      dunits = outunits ; disunits = outdisunits

      return
      end
