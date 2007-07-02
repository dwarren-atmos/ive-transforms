!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This routine calculates the value at a spatial point for all times.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine point_val(varout,varin,zbot,sigma,ztop,
     &                        ir,jr,z0,ix,iy,iz,it)

      integer, intent(in) ix,iy,iz
      real, intent(in) :: ir,jr,z0
      real, intent(in) :: varin(ix,iy,iz), zbot(ix,iy), sigma(iz)
      real, intent(out) :: varout(it)

      real ::  varz(1,1,iz), phys_hgt(iz), zbot_pnt
      integer :: k, n

      integer, parameter :: ifin=2
      real, parameter :: miss=-9999.0

      point_val(:) = miss
      zbot_pnt = horiz_interp_field(zbot,ir,jr,ix,iy)
      phys_hgt = (ztop-zbot_pnt)/ztop * sigm(:) + zbot_pnt

      do n=1,it

        ! First get the vertical profile at the desired point
        do k=1,iz
          varz(k) = horiz_interp_field(varin(:,:,k),ir,jr,ix,iy)
        enddo

        ! Second interpolate to the vertical level desired.
        call vert_interp_field(varz,1,1,iz,1,zbot_pnt,ifin,
     &                         point_val(n),1,z0,2,miss,phys_hgt)
      enddo

      return
      end
