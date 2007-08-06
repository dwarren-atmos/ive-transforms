      MODULE MODULE_NUMERICS

      CONTAINS

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        9 point smoother/desmoother
c
c        input variables:
c
c          z1     : first word address of array to be filtered
c          nfpass : number of passes through filter
c          isd=1  : 9 point smoother
c          isd=2  : 9 point smoother-desmoother
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine filt9(z1,ix,iy,iz,nfpass)
      implicit none

      integer, parameter :: isd = 1

      integer, intent(in) :: ix,iy,iz,nfpass
      real, intent(inout) :: z1(ix,iy,iz)

      integer :: i,j,k,nf,np,imin,imax,jmin,jmax
      real :: s1,z2(ix,iy,iz)

      imin = 2 ; imax = ix ; jmin = 2 ; jmax = iy

      if (nfpass.le.0) return

      z2(:,:,:) = z1(:,:,:)

      ! start loop for number of passes
      do np=1,nfpass
        do nf=1,isd
          if (nf.eq.1) then
            s1=0.25
          else
            s1=-0.25
          endif

          do k=1,iz
            !filter in the x-direction
            do j = jmin-1, jmax+1 ; do i = imin, imax
                 z2(i,j,k) = z1(i,j,k) + 
     &              s1*( z1(i-1,j,k) + z1(i+1,j,k) - 2.0*z1(i,j,k) )
            end do ; end do
            z1(:,:,:)=z2(:,:,:)

            !filter in the y-direction
            do j = jmin, jmax ; do i = imin-1, imax+1
               z2(i,j,k)=z1(i,j,k)
     c              +(z1(i,j-1,k)+z1(i,j+1,k)-2.0*z1(i,j,k))*s1
            end do ; end do
            z1(:,:,:) = z2(:,:,:)

          end do ! end k loop
        end do ! end nf loop
      end do ! end np loop

      return
      end subroutine

      END MODULE
