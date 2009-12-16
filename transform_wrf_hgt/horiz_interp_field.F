!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! horiz_interp_field()
! This routine interpolates a field in the horizontal.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function horiz_interp_field(fldvar,ipt,jpt,ix,iy)

      integer, intent(in) :: ix,iy
      real, intent(in) :: ipt,jpt,fldvar(ix,iy)

      integer :: ii,jj
      real :: a, b, c, d
      real :: horiz_interp_field

      ii = int(floor(ipt)) ; jj = int(floor(jpt))

      a = real(ii+1) - ipt ; b = ipt - real(ii)
      c = real(jj+1) - jpt ; d = jpt - real(jj)

      horiz_interp_field = c*(a*fldvar(ii,jj) + b*fldvar(ii+1,jj))
     &                 + d*(a*fldvar(ii,jj+1) + b*fldvar(ii+1,jj+1))

      return
      end




