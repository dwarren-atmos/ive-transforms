      integer(ive_ptr) function localize(olat,olon,locrad,stagpts,meta)
c     &                                ndims,dims,stag,min,max,
c     &                                miss,dunits,disunits,dim_names)

      use default_module
      use module_mapping

      implicit none

      real, intent(in) :: olat, olon, locrad
      character(len=1), intent(in) :: stagpts

      type(meta_ive), intent(out) :: meta
!      integer, intent(out) :: ndims,dims(4)
!      real, intent(out) :: stag(4),min(4),max(4),miss
!      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)

      integer :: ix, iy, iz, it, Ne, mnk
      integer :: ii, jj, kk, iptr, iflag(2)
      real :: phys(2), grid(2), xlon, xlat, dist, local_tmp

      real, external :: gc_local_calc

      data iflag/1,1/


      call set_var_meta(stagpts,meta)
c     &                 ,ndims,dims,stag,min,max
c     &                 ,miss,dim_names)
      ix = meta%dims(1) ; iy = meta%dims(2) ; iz = meta%dims(3)

      print *,'max/min(3)',meta%max(3),meta%min(3)

      meta%max(4)=0.0 ; meta%dim_names(4)='one' ; meta%dims(4)=1
      meta%disunits = ' ' ; meta%dunits = ' '

      mnk = ix*iy*iz
      localize = 0 ; localize = getmem(mnk)
      if(localize.eq.0) then
        write(*,*) 'Error:  Cannot allocate space for localize'
        return
      endif

      do ii=1,ix ; do jj=1,iy

        grid(1) = real(ii) ; grid(2) = real(jj)
        call index_2_phys(phys,grid,iflag,2,1)
        call phys_2_lonlat(phys(1),phys(2),xlon,xlat,1)
        dist = earth_dist(olon, olat, xlon, xlat)
        local_tmp = gc_local_calc(dist, locrad)

        do kk=1,iz
          call set_ptr(%val(localize),local_tmp,ii,jj,kk,1,ix,iy,iz,1)
        end do

      end do ; end do 

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Sets the i,j,k,n value of a pointer
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_ptr(varout,varin,ii,jj,kk,nn,ix,iy,iz,it)

      integer, intent(in) :: ii,jj,kk,nn,ix,iy,iz,it
      real, intent(in) :: varin
      real, intent(inout) :: varout(ix,iy,iz,it)

      varout(ii,jj,kk,nn) = varin

      return
      end subroutine

