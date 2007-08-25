!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   function vort_z(uvar_name, vvar_name, ...)
!
!   Calculates the vertical component of vorticity 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(ive_ptr) function pvort(unam,vnam,wnam,tnam,pnam,
     &                                ndims,dims,stag,min,max,miss,
     &                                dunits,disunits,dim_names)

      USE DEFAULT_MODULE

      character(len=*), intent(in) :: unam, vnam, wnam, tnam, pnam

      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout):: dunits,disunits,dim_names(4)
      
      integer :: i, ii, mnk, nt
      integer(ive_ptr) :: xvort,yvort,zvort,dxthe,dythe,dzthe,rho
      real :: r_val
      logical :: new

      integer(ive_ptr) :: deriv,rho_calc,vort_x,vort_y,vort_z

      nt=ntime ; mnk = nt*nz*(ny+1)*(nx+1)

      pvort=0 ; pvort = getmem(mnk)
      if(pvort.le.0) then
        print *,'Error: cannont allocate nesecary memory in pvort'
        return
      end if

      do i=1,mnk ! initialize pvort to 0
        call s_val(%val(pvort),i,0.0)
      end do

      ! Get dtheta/dx
      dxthe = deriv(tnam,'X',1,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (dxthe.le.0) return

      ! Get the x-vorticity
      xvort = vort_x(vnam,wnam,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (xvort.le.0) return

      ! pv calc component wise and then clear mem (saves memorey)
      call sum_prod(%val(pvort),%val(dxthe),%val(xvort),nx+1,ny+1,nz,nt)
      call freemem(dxthe) ; call freemem(xvort)

      ! Get dtheta/dy
      dythe = deriv(tnam,'Y',1,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (dythe.le.0) return

      ! Get the y-vorticity
      yvort = vort_y(unam,wnam,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (yvort.le.0) return

      ! pv calc component wise and then clear mem (saves memorey)
      call sum_prod(%val(pvort),%val(dythe),%val(yvort),nx+1,ny+1,nz,nt)
      call freemem(dythe) ; call freemem(yvort)

      ! Get dtheta/dz
      dzthe = deriv(tnam,'Z',1,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (dzthe.le.0) return

      ! Get the z-vorticity
      zvort = vort_z(unam,vnam,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (zvort.le.0) return

      ! pv calc component wise and then clear mem (saves memorey)
      call abs_vort(%val(zvort),%val(zvort),%val(fcor),nx+1,ny+1,nz,nt)
      call sum_prod(%val(pvort),%val(dzthe),%val(zvort),nx+1,ny+1,nz,nt)
      call freemem(dzthe) ; call freemem(zvort)


      ! Get density
      rho = rho_calc(pnam,tnam,ndims,dims,stag,
     &      min,max,miss,dunits,disunits,dim_names)
      if (rho.le.0) return

      do i=1,mnk ! replace rho with 1/rho
        call s_val(%val(rho),i,1/r_val(%val(rho),i))
      end do

      call sum_prod(%val(pvort),%val(pvort),%val(rho),nx+1,ny+1,nz,nt)
      call freemem(rho)

c      call pvort_calc(%val(pvort),%val(xvort),%val(yvort),%val(zvort),
c     &                %val(dxthe),%val(dythe),%val(dzthe),%val(rho),
c     &                %val(fcor),(nx+1),(ny+1),nz,ntime)

      dunits='K*m*m/kg/s' ; disunits='10^-6*K*m*m/kg/s'

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Does the potential vorticity calculation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pvort_calc(pv,xvort,yvort,zvort,dxthe,dythe,dzthe,
     &                      rho,f,ix,iy,iz,it)

      integer, intent(in) :: ix,iy,iz,it
      real, intent(in) :: xvort(ix,iy,iz,it), yvort(ix,iy,iz,it),
     &                    zvort(ix,iy,iz,it), dxthe(ix,iy,iz,it),
     &                    dythe(ix,iy,iz,it), dzthe(ix,iy,iz,it),
     &                    rho(ix,iy,iz,it), f(ix,iy)

      real, intent(out) :: pv(ix,iy,iz,it)

      integer i,j,k,n

      do n=1,it ; do k=1,iz ; do j=1,iy ; do i=1,ix

        pv(i,j,k,n) = 1/rho(i,j,k,n) * (
     &                xvort(i,j,k,n)*dxthe(i,j,k,n) +
     &                yvort(i,j,k,n)*dythe(i,j,k,n) +
     &               (zvort(i,j,k,n) + f(i,j))*dzthe(i,j,k,n) )
      end do ; end do ; end do ; end do

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! abs_vort(abs_vort,rel_vort,fcor,ix,iy,iz,it)
!
! calculates the absoulute vorticity given the relative and planetary
! contributions.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine abs_vort(avort,rvort,fcor,ix,iy,iz,it)

      integer, intent(in) :: ix, iy, iz, it
      real, intent(in) :: rvort(ix,iy,iz,it), fcor(ix,iy)
      real, intent(inout) :: avort(ix,iy,iz,it)

      integer i, j, k, n

      do n=1,it ; do k=1,iz ; do j=1,iy ; do i=1,ix
        avort(i,j,k,n) = rvort(i,j,k,n) + fcor(i,j)
      end do ; end do ; end do ; end do

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine sum_prod(var_sum,var1,var2,ix,iy,iz,it)    
!
!  Adds the product of two variables to a third
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sum_prod(varsum,var1,var2,ix,iy,iz,it)

      integer, intent(in) :: ix,iy,iz,it
      real, intent(in) :: var1(ix,iy,iz,it), var2(ix,iy,iz,it)
      real, intent(inout) :: varsum(ix,iy,iz,it)

      integer i, j, k, n

      do n=1,it ; do k=1,iz ; do j=1,iy ; do i=1,ix
        varsum(i,j,k,n) = varsum(i,j,k,n) + var1(i,j,k,n)*var2(i,j,k,n)
      end do ; end do ; end do ; end do

      return
      end

