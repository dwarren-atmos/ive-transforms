!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Interpolates variable from stagerd mesh to thermo points  
!
!  var_t- output field at thermo points (should be unity on input)
!  var  - input field on stagard mesh
!  iu   - 1 -> u momentum points
!  iv   - 1 -> v momentum points
!  iw   - 1 -> w momentum points
!  nx   - number of u momentum points
!  ny   - number of v momentum points
!  nz   - number of vertical thermo points
!  nt   - number of time levels
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine var_at_t(var_t,var,iu,iv,iw,nx,ny,nz,nt)

      implicit none

      integer, intent(in) :: iu,iv,iw,nx,ny,nz,nt
	  real, intent(in) ::  var(nx+1-iu,ny+1-iv,nz+iw,nt)

      real, intent(inout) ::  var_t(nx+1,ny+1,nz,nt)
  
      real :: tmpvar
      integer n,k,j,i
c
      var_t(:,:,:,:) = 1.0

      if (iv.eq.1.and.iu.eq.1) then
        write(6,*) 'Error:  You are not on the mesh'
        return
      endif
c
      do n=1,nt; do k=1,nz; 
        do i=1+iu,nx+1-iu; do j=1+iv,ny+1-iv;
          tmpvar = .5*(var(i-iu,j-iv,k,n) + var(i,j,k,n))
          var_t(i,j,k,n) = var_t(i,j,k,n) * tmpvar
      end do; end do;

      if(iv.eq.1) then
        do i = 1,nx+1
            var_t(i,1,k,n) = 1.5*var_t(i,2,k,n) - 
     &           0.5*var_t(i,3,k,n)
            var_t(i,ny+1,k,n) = 1.5*var_t(i,ny,k,n) - 
     &           0.5*var_t(i,ny-1,k,n)
        enddo
c
      else if(iu.eq.1) then
        do j = 1,ny+1
            var_t(1,j,k,n) = 1.5*var_t(2,j,k,n) - 
     &          0.5*var_t(3,j,k,n)
                 var_t(nx+1,j,k,n) = 1.5*var_t(nx,j,k,n) - 
     &          0.5*var_t(nx-1,j,k,n)
        enddo
      endif

      end do; end do;
c
      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Interpolates variable from stagerd mesh u momentum points
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine var_at_u(var_u,var,iu,iv,iw,nx,ny,nz,nt)

      implicit none
c  
      integer iu,iv,iw,nx,ny,nz,nt
      real var(nx+1-iu,ny+1-iv,nz+iw,nt)
      real var_u(nx,ny+1,nz,nt),tmpvar
c
      integer n,k,j,i
c 
      var_u(:,:,:,:) = 1.0
      do n=1,nt ; do k=1,nz ; do i=1,nx ; 
        do j=1+iv,ny+1-iv
          tmpvar = 
     &          .25*(var(i,j,k,n) + var(i,j-iv,k,n)
     &          + var(i+1-iu,j,k,n) + var(i+1-iu,j-iv,k,n))
          var_u(i,j,k,n) = var_u(i,j,k,n) * tmpvar
        end do

        if (iv.eq.1) then
          var_u(i,1,k,n) = var_u(i,1,k,n) *
     >             (3/2*(.5*(var(i,1,k,n)  
     >       +      var(i+1,1,k,n))) 
     >                  -  1/2*(.5*(var(i,2,k,n)  
     >       +      var(i+1,2,k,n))))  
          var_u(i,ny+1,k,n) = var_u(i,ny+1,k,n) *
     >           (3/2*(.5*(var(i,ny,k,n)  
     >       +           var(i+1,ny,k,n))) 
     >                  -  1/2*(.5*(var(i,ny-1,k,n)  
     >       +           var(i+1,ny-1,k,n))))  
         end if
      end do ; end do ; end do
c
      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Interpolates variable from stagerd mesh v momentum points
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine var_at_v(var_v,var,iu,iv,iw,nx,ny,nz,nt)

      implicit none

      integer, intent(in) :: iu,iv,iw,nx,ny,nz,nt
	  real, intent(in) :: var(nx+1-iu,ny+1-iv,nz+iw,nt)
      real, intent(inout) :: var_v(nx+1,ny,nz,nt)

      real :: tmpvar
      integer :: n,k,j,i

      var_v(:,:,:,:) = 1.0
      do n=1,nt ; do k=1,nz ; do j=1,ny

        do i=1+iu,nx+1-iu
          tmpvar = 
     &       .25*(var(i,j,k,n) + var(i,j+1-iv,k,n)
     &    +   var(i+iu,j,k,n) + var(i+iu,j+1-iv,k,n))
              var_v(i,j,k,n) = var_v(i,j,k,n) * tmpvar
        enddo
c
        if (iu.eq.1) then

          var_v(1,j,k,n) = var_v(1,j,k,n)*
     &           (3/2*(.5*(var(1,j,k,n)  
     &                   + var(1,j+1,k,n))) 
     &          -  1/2*(.5*(var(2,j,k,n)  
     &                    + var(2,j+1,k,n))))  

          var_v(nx+1,j,k,n) = var_v(nx+1,j,k,n) *
     &              (3/2*(.5*(var(nx,j,k,n)  
     &                      + var(nx,j+1,k,n))) 
     &            -  1/2*(.5*(var(nx-1,j,k,n)  
     &                      + var(nx-1,j+1,k,n))))  

        end if
      end do ; end do ; end do

      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Interpolates variable from stagerd mesh w momentum points
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine var_at_w(var_w,var,sgz,wgz,iu,iv,iw,nx,ny,nz,nt)

      implicit none

      integer iu,iv,iw,nx,ny,nz,nt
	  real var(nx+1-iu,ny+1-iv,nz+iw,nt)
      real var_w(nx+1,ny+1,nz+1,nt)
      real var_t(nx+1,ny+1,nz,nt)
	  real sgz(nz),wgz(nz+1)

      integer n,k,j,i
	  real zfrac,tmpvar
c
      var_w(:,:,:,:) = 1.0
      if (iw.eq.0) then
c
        do n=1,nt ; do i=1,nx+1 ; do j=1,ny+1
          do k=1,nz
            var_t(i,j,k,n) = 1.0
          enddo
        enddo ; enddo ; enddo
c
        call var_at_t(var_t,var,iu,iv,iw,nx,ny,nz,nt)
c
        do n=1,nt ; do i=1,nx+1 ; do j=1,ny+1
          do k=2,nz
            zfrac = (wgz(k)-sgz(k-1))/(sgz(k)-sgz(k-1))
            tmpvar = var_t(i,j,k-1,n) + 
     >              (var_t(i,j,k,n) - var_t(i,j,k-1,n))*zfrac
            var_w(i,j,k,n) = var_w(i,j,k,n)*tmpvar
          end do
c
c  Take care of top and bottom boundries
c

          var_w(i,j,1,n) = var_w(i,j,1,n)*
     >                    (var_w(i,j,2,n) +
     >                    (var_w(i,j,3,n)-var_w(i,j,2,n)) *
     >                    (wgz(1)-sgz(1))/(sgz(2)-sgz(1)))
c
          var_w(i,j,nz+1,n) = var_w(i,j,nz+1,n) *
     >                       (var_w(i,j,nz-1,n)  +
     >                       (var_w(i,j,nz,n)-var_w(i,j,nz-1,n)) *
     >                       (wgz(nz+1)-sgz(nz))/(sgz(nz)-sgz(nz-1)))
        end do ; end do ; end do
c
      else
        do n = 1,nt ; do k = 1,nz+1 ; do j = 1,ny+1 ; do i = 1,nx+1
           var_w(i,j,k,n) = var_w(i,j,k,n)*var(i,j,k,n)
        enddo ; enddo ; enddo ; enddo
      end if
c
      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calculates the stagering for the given boundaries
!  iu = 1 -> u points
!  iv = 1 -> v points
!  iw = 1 -> w points
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine stagering_calc(iu,iv,iw,it,dims,plmin,min,dx,dy,twod)

      implicit none

      integer, intent(in) :: twod
      real, intent(in) :: dims(4),min(4),plmin(4),dx,dy
      integer, intent(out) :: iu,iv,iw,it

      logical :: requ
c
      it = 0 ; iu = 0 ; iv = 0 ; iw = 0
      if ((dims(1).eq.1).and.(dims(2).eq.1).and.(dims(3).eq.1)) then 
        it = 1
      else if (requ(plmin(1),min(1)-dx/2.)) then
        iu = 1
      else if (requ(plmin(2),min(2)-dy/2.).and.(twod.ne.1)) then
        iv = 1
      else if (requ(plmin(3),min(3))) then
        iw = 1
      end if

      return
      end subroutine
