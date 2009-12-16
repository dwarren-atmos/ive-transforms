!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!   interp_hght - function that interpolates a variable in the
!                 vertical based on height.
!
!    varble - input variable to interpolate
!      hght - vector of geopotential heights
!     level - vertical level desired
!        iz - number of vertical grid points
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      function interp_hght(varble, hght, level, iz)

      integer, intent(in) :: iz
      real, intent(in)    :: varble(iz), hght(iz), level

      integer :: k, klev
      real :: m, interp_hght

      if ((hght(1) .le. level) .AND. (hght(iz) .ge. level)) then

        ! search for appropriate level
        do k = 2, iz
          if ( hght(k) > level ) then
            klev = k - 1
            exit
          endif
        enddo

        ! linearly interpolate
        m = (varble(klev+1) - varble(klev)) / (hght(klev+1)-hght(klev))
        interp_hght = m * (level - hght(klev)) + varble(klev)

      else  ! level outside domain

        interp_hght = missing

      endif

      return
      end function

