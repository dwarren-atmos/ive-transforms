!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  function regrid(var1_name,var2_name, ...)
!
!  regrids var1 to var2's grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(ive_ptr) function regrid(var1_name,var2_name,
     &                           ndims,dims,stag,min,max,
     &                           miss,dunits,disunits,dim_names,
     &                           error)
      use default_module
      implicit none

      character(len=*), intent(in) :: var1_name, var2_name
c
      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(out) :: dunits,disunits,dim_names(4)
      logical new,error

      CHARACTER*80 driverstuff
      integer :: i, j, k, l,ndims1,dims1(4)
      integer(ive_ptr) :: var1,var2
      real :: stag1(4),min1(4),max1(4),miss1,val
      double precision extrap
      character(len=80) :: dunits1,disunits1,dim_names1(4)
      logical :: new1, new2
      regrid = 0

      var1 = getvar(var1_name,
     &     ndims1,dims1,stag1,min1,max1,miss1,
     &     dunits1,disunits1,dim_names1,new1)

      var2 = getvar(var2_name,
     &     ndims,dims,stag,min,max,miss,
     &     dunits,disunits,dim_names,new)

      if(ndims .ne. ndims1) then
         error = .true.
         return
      endif
      stag=(dims1-dims)/2.!pos if shrinking neg if growing
      regrid=getmem(dims(1)*dims(2)*dims(3)*dims(4))
      call regrid_fill(%val(var1), %val(regrid), stag, dims(1), dims(2),
     &     dims(3), dims(4), dims1(1), dims1(2), dims1(3), dims1(4),
     &     miss)
      return
      end

      subroutine regrid_fill(from, to, stag, nx, ny, nz, nt, 
     &     nx1, ny1, nz1, nt1, miss)
      real:: from, to, stag, miss
      integer :: nx,ny,nz,nt,nx1,ny1,nz1,nt1
      dimension stag(4)
      dimension from( nx1, ny1, nz1, nt1)
      dimension to( nx, ny, nz, nt)

      double precision extrap
      integer:: i, j, k, l
      do l=1, nt
         do k=1, nz
            do j=1, ny
               do i=1, nx
                  to(i,j,k,l)=extrap(from,nx1,ny1,nz1,nt1,
     &                 real(i)-stag(1),real(j)-stag(2),
     &                 real(k)-stag(3),real(l)-stag(4),miss)
               enddo
            enddo
         enddo
      enddo
      return
      end
