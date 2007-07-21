      integer(ive_ptr) function gc_local(olat,olon,locrad,stagpts,
     &                                ndims,dims,stag,min,max,
     &                                miss,dunits,disunits,dim_names)

      use default_module
      use module_mapping

      implicit none

      real, intent(in) :: olat, olon, locrad
      character(len=1), intent(in) :: stagpts

      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)

      integer :: ix, iy, iz, Ne, mnk
      integer :: ii, jj, kk, iptr, iflag(2)
      real :: phys(2), grid(2), xlon, xlat, dist, gc_tmp

      real, external :: gc_local_calc

      data iflag/1,1/

      print *,olat,olon,locrad

      dim_names(3) = 'nz' ; 
      max(1) = real(nx-1)*dx ; max(2) = real(nx)*dy ; min(:) = 0.0
      min(1:2) = 0.0 ; min(3) = sgz_min ; min(4) = 0.0
      select case (stagpts)
      case ('U')
        ix = nx ; iy = ny+1 ; iz = nz
        dim_names(1) = 'nxm1' ; dim_names(2) = 'ny'  
        min(1) = 0.5*dx ; max = real(ix-1)*dx
      case ('V')
        ix = nx+1 ; iy = ny ; iz = nz
        dim_names(1) = 'nx' ; dim_names(2) = 'nym1'
        min(2) = 0.5*dy ; max = real(iy-1)*dy
      case ('T')
        ix = nx+1 ; iy = ny+1 ; iz = nz
        dim_names(1) = 'nx' ; dim_names(2) = 'ny'
      case ('W')
        ix = nx+1 ; iy = ny+1 ; iz = nz+1
        dim_names(1) = 'nx' ; dim_names(2) = 'ny'
        dim_names(3) = 'nzp1' ; 
      end select




      mnk = ix*iy
      gc_local = 0 ; gc_local = getmem(mnk)
      if(gc_local.eq.0) then
        write(*,*) 'Error:  Cannot allocate space for kgain'
        return
      endif

      do ii=1,ix ; do jj=1,iy

        grid(2) = real(ii) ; grid(1) = real(jj)
        call index_2_phys(phys,grid,iflag,2,1)
        call phys_2_lonlat(phys(1),phys(2),xlon,xlat,1)
        dist = earth_dist(olon, olat, xlon, xlat)

        gc_tmp = gc_local_calc(dist, locrad)

        !do kk=1,iz
          iptr = (ii-1)*ix + jj
          call s_val(%val(gc_local),iptr,gc_tmp)
        !end do
      end do ; end do 
        print *, iptr,mnk,gc_tmp

      ndims = 4
      dims(1) = ix ; dims(2) = iy ; dims(3) = ix ; dims(4) = 1
      dim_names(4) = 'one'
      max(3:4) = 0.0 ; min(:) = 0.0 ; stag(:) = 0.0
      miss = -9999.0

      dunits = ' ' ; disunits = ' '

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   gc_local - function that calculates the Gaspari-Cohn localization 
!              5th order piecewise function
!
!             dist - horizontal distance between points
!  localize_radius - radius to scale function
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function gc_local_calc(dist, localize_radius)

      real, intent(in) :: dist, localize_radius
      real :: a, b, gcv, gc_local_calc
      
      a = localize_radius / 2.
      gcv = dist / a ! gaspari-cohn variable

      if ( dist .le. a ) then
        gc_local_calc = -0.25*(gcv)**5 + 0.5*(gcv)**4 +
     &             (5./8.)*(gcv)**3 - (5./3.)*(gcv)**2 + 1.
      else
        gc_local_calc = (1./12.)*(gcv)**5 - 0.5*(gcv)**4 +
     &             (5./8.)*(gcv)**3 + (5./3.)*(gcv)**2 -
     &             5.*(gcv) + 4 - (2./3.)*(gcv)**(-1)
      endif

      return
      end function

