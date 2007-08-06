      integer(ive_ptr) function kgain(obs_name,upd_name,R,x0,y0,z0,
     &                                mapon,calc_type,Nbin,nfilt,
     &                                ndims,dims,stag,min,max,
     &                                miss,dunits,disunits,dim_names)

      use default_module
      !use module_numerics

      implicit none

      real, intent(in) :: x0, y0, z0, R
      integer, intent(in) :: Nbin, nfilt 
      character(len=*), intent(in) :: calc_type
      character(len=*), intent(in) :: obs_name, upd_name
      logical :: mapon

      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)

      integer :: ix, iy, iz, Ne, mnk
      integer :: iu, iv, iw, it,kk, iflag(2)
      real :: ir, jr, phys(2), grid(2)
      character (len=80) :: outunits,outdisunits
      logical :: new, lkgain, lcov, lcor, lhier
      integer(ive_ptr) :: xvar,xprior

      real, allocatable :: ye(:)

      data iflag/1,1/

      select case (trim(calc_type))

      case('kgain','k','gain','KGAIN')
        lkgain = .true. ; lcov = .false. ; lcor = .false.
        lhier = .false.
      case('cov','covariance','ENS_COV')
        lkgain = .false. ; lcov = .true. ; lcor = .false.
        lhier = .false.
      case('cor','correlation','ENS_COR')
        lkgain = .false. ; lcov = .false. ; lcor = .true.
        lhier = .false.
      case('hier','hier_loc','hloc','HIERLOC')
        lkgain = .false. ; lcov = .false. ; lcor = .false. 
        lhier = .true.
      case default
        print *,'UNKOWN FIELD TO CALCULATE, RETURNING'
        kgain = 0
      end select

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

      if(lkgain) then
        call kgain_calc(%val(kgain),ye,%val(xprior),R,mnk,Ne)
      elseif(lcov) then
        call covcor_calc(%val(kgain),ye,%val(xprior),mnk,Ne,lcor)
      elseif(lcor) then
        call covcor_calc(%val(kgain),ye,%val(xprior),mnk,Ne,lcor)
      elseif(lhier) then
        call hierloc_calc(%val(kgain),ye,%val(xprior),R,mnk,Ne,Nbin)
        call filt9(%val(kgain),ix,iy,iz,nfilt)
      end if


      dims(1) = ix ; dims(2) = iy ; dims(3) = iz ; dims(4) = 1
      dim_names(4) = 'one' ; max(4) = 0.0 ; min(4) = 0.0
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

      mean_fac = 1.0/real(Ne) ; cov_fac = 1.0/real(Ne-1)

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
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calculates covariance or correlation between a point and a field
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine covcor_calc(covcor,ye,xe,mnk,Ne,lcor)

      integer, intent(in) :: mnk,Ne
      real, intent(in) :: ye(Ne), xe(mnk,Ne)
      logical, intent(in) :: lcor

      real, intent(out) :: covcor(mnk)

      integer :: i,n
      real :: ye_std,xe_std
      real :: ye_pert(Ne), xe_pert(Ne)
      real :: ye_mean, xe_mean,mean_fac,cov_fac

      mean_fac = 1.0/real(Ne) ; cov_fac = 1.0/real(Ne-1)

      covcor(:) = 0.0

      ye_mean = sum(ye(:))*mean_fac
      ye_pert(:) = ye(:) - ye_mean
      ye_std = sqrt(sum(ye_pert(:)*ye_pert(:)))

      if(lcor) then
        do i=1,mnk
          xe_mean = sum(xe(i,:))*mean_fac
          xe_pert(:) = xe(i,:) - xe_mean
          xe_std = sqrt(sum(xe_pert(:)*xe_pert(:)))
          covcor(i) = sum(ye_pert(:)*xe_pert(:))/(xe_std*ye_std)
        enddo
      else
        do i=1,mnk
          xe_mean = sum(xe(i,:))*mean_fac
          xe_pert(:) = xe(i,:) - xe_mean
          covcor(i) = sum(ye_pert(:)*xe_pert(:))*cov_fac
        enddo
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calculates hierarchial localization factor 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine hierloc_calc(hloc,ye,xe,R,mnk,Ne,Nbin)


      integer, intent(in) :: mnk,Ne,Nbin
      real, intent(in) :: ye(Ne), xe(mnk,Ne), R

      real, intent(out) :: hloc(mnk)

      integer :: i,n,Nmem,bin
      real :: ye_mean, ye_var, xe_mean, mean_fac, cov_fac, kdenom, R2
      real :: bin_fac, beta, beta_sum(mnk), beta_var(mnk), hloc_tmp

      integer, allocatable :: pnt_mem(:)
      real, allocatable :: ye_bin(:), ye_pert(:),
     &                     xe_bin(:,:), xe_pert(:)

      R2 = R**2
      Nmem = int(Ne/Nbin)  
      allocate(pnt_mem(Nmem))
      allocate(ye_bin(Nmem)) ; allocate(ye_pert(Nmem))
      allocate(xe_bin(mnk,Nmem)) ; allocate(xe_pert(Nmem))
      
      mean_fac = 1.0/real(Nmem) ; cov_fac = 1.0/real(Nmem-1)
      bin_fac = 1.0/real(Nbin-1)

      beta_sum(:) = 0.0 ; beta_var(:) = 0.0
      do bin=1,Nbin       

        do i=1,Nmem
          pnt_mem(i) = Nmem*(bin-1)+i
        end do

        do n=1,Nmem
          ye_bin(n) = ye(pnt_mem(n))
          xe_bin(:,n) = xe(:,pnt_mem(n))
        end do

        ye_mean = sum(ye_bin(:))*mean_fac
        ye_pert(:) = ye_bin(:) - ye_mean
        ye_var = sum(ye_pert(:)*ye_pert(:))*cov_fac
        kdenom = 1/(ye_var + R2)

        do i=1,mnk
          xe_mean = sum(xe_bin(i,:))*mean_fac
          xe_pert(:) = xe_bin(i,:) - xe_mean
          beta = sum(ye_pert(:)*xe_pert(:))*cov_fac*kdenom
          beta_sum(i) = beta + beta_sum(i)
          beta_var(i) = beta**2 + beta_var(i)
        end do

      end do

      do i=1,mnk
        hloc_tmp = (beta_sum(i)**2/beta_var(i) - 1)*bin_fac 
        hloc(i) = max(hloc_tmp,0.0)
      end do

      deallocate(pnt_mem,ye_bin,ye_pert,xe_bin,xe_pert)

      return
      end
