      integer(ive_ptr) function kgain(obs_name,upd_name,R,x0,y0,z0,
     &                                mapon,ndims,dims,stag,min,max,
     &                                miss,dunits,disunits,dim_names)

      use default_module

      implicit none

      real, intent(in) :: x0, y0, z0, R
      character(len=*), intent(in) :: obs_name, upd_name 
      logical :: mapon

      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)

      integer :: ix, iy, iz, Ne, mnk
      integer :: iu, iv, iw, it,kk, iflag(2)
      real :: ir, jr, phys(2), grid(2)
      character (len=80) :: outunits,outdisunits
      logical :: new
      integer(ive_ptr) :: xvar,xprior

      real, allocatable :: ye(:)

      data iflag/1,1/

      ! Find location to interpolate to in the horizontal 
      if(mapon) then
        call lonlat_2_phys(phys(1),phys(2),x0,y0,1)
        call phys_2_index(phys,grid,iflag,2,1)
        ir=grid(1) ; jr=grid(2)
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
      ix = dims(1) ; iy = dims(2) ; iz = dims(3) ; Ne = dims(4)

      mnk = ix*iy*iz
      kgain = 0 ; kgain = getmem(mnk)
      if(kgain.eq.0) then
        write(*,*) 'Error:  Cannot allocate space for kgain'
        return
      endif

      call kgain_calc(%val(kgain),ye,%val(xprior),R,mnk,Ne)

      dims(1) = ix ; dims(2) = iy ; dims(3) = iz ; dims(4) = 1
      dim_names(4) = 'one'
      dunits = outunits ; disunits = outdisunits

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calculates the kalman gain
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine kgain_calc(kgain,ye,xprior,R,mnk,Ne)

      integer, intent(in) :: mnk,Ne
      real, intent(in) :: ye(Ne), xprior(mnk,Ne), R
      real, intent(out) :: kgain(mnk)

      integer :: i,n
      real :: mean_fac,cov_fac,ye_var,cov,kdenom
      real :: ye_pert(Ne), xprior_pert(Ne)
      real :: ye_mean, xprior_mean, R2

      mean_fac = 1.0/real(Ne)
      cov_fac = 1.0/real(Ne-1)

      kgain(:) = 0.0
      R2 = R**2

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
        
