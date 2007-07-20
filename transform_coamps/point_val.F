!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This routine calculates the value at a spatial point for all times.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine point_val(varout,varin,zbot,sigm,ztop,
     &                        ir,jr,z0,ix,iy,iz,it)

      integer, intent(in) :: ix,iy,iz,it
      real, intent(in) :: ir,jr,z0,ztop
      real, intent(in) :: varin(ix,iy,iz,it), zbot(ix,iy), sigm(iz)
      real, intent(out) :: varout(it)

      real ::  varz(iz), phys_hgt(iz), zbot_pnt, z00
      integer :: k, n
      logical :: lagl

      integer, parameter :: ifin=2
      real, parameter :: miss=-9999.0

      real :: horiz_interp_field

      lagl = .false. 
      if(z0.le.0.0) lagl=.true.
      z00 = abs(z0)

      varout(:) = miss
      zbot_pnt = horiz_interp_field(zbot,ir,jr,ix,iy)

      phys_hgt = (ztop-zbot_pnt)/ztop * sigm(:) + zbot_pnt

      if(lagl) z00 = z00 + zbot_pnt

      !print '(A,1x,3(F8.2,1x))','Interpolating field to:',ir,jr,z00

      do n=1,it

        ! First get the vertical profile at the desired point
        do k=1,iz
          varz(k) = horiz_interp_field(varin(:,:,k,n),ir,jr,ix,iy)
        enddo

        ! Second interpolate to the vertical level desired.
        call vert_interp_field(varz,1,1,iz,1,zbot_pnt,ifin,
     &                         varout(n),1,z00,2,miss,phys_hgt)
      enddo

      return
      end
