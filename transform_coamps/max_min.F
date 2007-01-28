#include "pointer.icl"

      ive_ptr function maxmin_reg(varnm,zlid,
     &                            ndims,dims,stag,min,max,miss,
     &                            dunits,disunits,dim_names)
#include "default.icl"

      real, intent(in) :: zlid
      character(len=*), intent(in) :: varnm

      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*) , intent(inout) :: dunits,disunits,
     &                                    dim_names(4)


c      integer, parameter :: npts = 218
      integer, parameter :: npts = 155
      integer :: ipts(npts), jpts(npts), k
      logical :: new
      ive_ptr :: var

      data ipts/
     &      77,78,79,80,81,78,79,80,81,82,78,79,80,81,82,78,79,80,81,82,
     &      78,79,80,81,82,79,80,81,82,83,79,80,81,82,83,79,80,81,82,83,
     &      79,80,81,82,83,80,81,82,83,84,80,81,82,83,84,80,81,82,83,84,
     &      81,82,83,84,85,81,82,83,84,85,81,82,83,84,85,81,82,83,84,85,
     &      81,82,83,84,85,81,82,83,84,85,81,82,83,84,85,82,83,84,85,86,
     &      82,83,84,85,86,83,84,85,86,87,83,84,85,86,87,84,85,86,87,88,
     &      84,85,86,87,88,85,86,87,88,89,86,87,88,89,90,87,88,89,90,91,
     &      87,88,89,90,91,87,88,89,90,91,87,88,89,90,91/

c     &     78,79,80,81,82,78,79,80,81,82,83,78,79,80,81,82,83,84,79,80,
c     &     81,82,83,84,79,80,81,82,83,84,79,80,81,82,83,84,79,80,81,82,
c     &     83,84,80,81,82,83,84,80,81,82,83,84,80,81,82,83,84,80,81,82,
c     &     83,84,81,82,83,84,85,81,82,83,84,85,86,81,82,83,84,85,86,81,
c     &     82,83,84,85,86,87,81,82,83,84,85,86,87,81,82,83,84,85,86,87,
c     &     88,81,82,83,84,85,86,87,88,82,83,84,85,86,87,88,82,83,84,85,
c     &     86,87,88,89,83,84,85,86,87,88,89,83,84,85,86,87,88,89,90,84,
c     &     85,86,87,88,89,90,84,85,86,87,88,89,90,91,85,86,87,88,89,90,
c     &     91,85,86,87,88,89,90,91,92,86,87,88,89,90,91,92,93,94,87,88,
c     &     89,90,91,92,93,94,95,87,88,89,90,91,92,93,94,95,96,87,88,89,
c     &     90,91,92,93,94,95,96,87,88,89,90,91,92,93,94,95,96,97/
      
      data jpts/
     & 100,100,100,100,100,99,99,99,99,99,98,98,98,98,98,97,97,97,97,97,
     &      96,96,96,96,96,95,95,95,95,95,94,94,94,94,94,93,93,93,93,93,
     &      92,92,92,92,92,91,91,91,91,91,90,90,90,90,90,89,89,89,89,89,
     &      88,88,88,88,88,87,87,87,87,87,86,86,86,86,86,85,85,85,85,85,
     &      84,84,84,84,84,83,83,83,83,83,82,82,82,82,82,81,81,81,81,81,
     &      80,80,80,80,80,79,79,79,79,79,78,78,78,78,78,77,77,77,77,77,
     &      76,76,76,76,76,75,75,75,75,75,74,74,74,74,74,73,73,73,73,73,
     &      72,72,72,72,72,71,71,71,71,71,70,70,70,70,70/

c     &100,100,100,100,100,99,99,99,99,99,99,98,98,98,98,98,98,98,97,97,
c     &     97,97,97,97,96,96,96,96,96,96,95,95,95,95,95,95,94,94,94,94,
c     &     94,94,93,93,93,93,93,92,92,92,92,92,91,91,91,91,91,90,90,90,
c     &     90,90,89,89,89,89,89,88,88,88,88,88,88,87,87,87,87,87,87,86,
c     &     86,86,86,86,86,86,85,85,85,85,85,85,85,84,84,84,84,84,84,84,
c     &     84,83,83,83,83,83,83,83,83,82,82,82,82,82,82,82,81,81,81,81,
c     &     81,81,81,81,80,80,80,80,80,80,80,79,79,79,79,79,79,79,79,78,
c     &     78,78,78,78,78,78,77,77,77,77,77,77,77,77,76,76,76,76,76,76,
c     &     76,75,75,75,75,75,75,75,75,74,74,74,74,74,74,74,74,74,73,73,
c     &     73,73,73,73,73,73,73,72,72,72,72,72,72,72,72,72,72,71,71,71,
c     &     71,71,71,71,71,71,71,70,70,70,70,70,70,70,70,70,70,70/


      var = getvar(trim(varnm),ndims,dims,stag,min,max,miss,dunits,
     &                    disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot get field ', trim(varnm)
        return
      end if
      
      maxmin_reg = 0
      maxmin_reg = getmem(dims(4))

      call mean_calc(%val(maxmin_reg),%val(var),zlid,
     &                  ipts,jpts,npts,%val(sgz),%val(zbot_u),ztop,
     &                  dims(1),dims(2),dims(3),dims(4))

      do k = 1,3
        dims(k)=1 ; stag(k)=0.0
        min(k)=0.0 ;  max(k)=0.0
        dim_names(k) = 'one'
      enddo


      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE TO DO CALCULATION OF MAX OVER A REGION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine max_calc(extrema,fldvar,zlid,ipts,jpts,npts,
     &                       sigm,topo,ztop,ix,iy,iz,it) 


      integer, intent(in) :: ix,iy,iz,it,npts
      integer, intent(in) :: ipts(npts), jpts(npts)
      real, intent(in) :: fldvar(ix,iy,iz,it),zlid,
     &                    sigm(iz),topo(ix,iy),ztop

      real, intent(out) :: extrema(it)

      integer ii,kk,n,i,j
      real :: zfrac, level

      do n=1,it
        extrema(n) = -99999.0
        do ii = 1,npts
          i = ipts(ii) ; j = jpts(ii)
          zfrac = (ztop-topo(i,j))/ztop 
          do kk=1,iz
            level = zfrac*sigm(kk)
              
            if(n.eq.1.and.ii.lt.75) then
             print '(4I3,2x,6F14.4)',n,i,j,kk,
     &       zfrac,sigm(kk),ztop,topo(i,j),level,fldvar(i,j,kk,n)
            end if
            if(level.le.zlid) then
              if(extrema(n).lt.fldvar(i,j,kk,n)) then
                extrema(n)=fldvar(i,j,kk,n)
              end if
            end if
          end do
        enddo
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE TO DO CALCULATION OF MEAN OVER A REGION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mean_calc(meanval,fldvar,zlid,ipts,jpts,npts,
     &                       sigm,topo,ztop,ix,iy,iz,it) 


      integer, intent(in) :: ix,iy,iz,it,npts
      integer, intent(in) :: ipts(npts), jpts(npts)
      real, intent(in) :: fldvar(ix,iy,iz,it),zlid,
     &                    sigm(iz),topo(ix,iy),ztop

      real, intent(out) :: meanval(it)

      integer ii,kk,n,i,j
      real :: zfrac, level, DZ, zvol

      do n=1,it
        meanval(n) = 0.0 ; zvol = 0.0
        do ii = 1,npts
          i = ipts(ii) ; j = jpts(ii)
          zfrac = (ztop-topo(i,j))/ztop 
          do kk=1,iz
            level = zfrac*sigm(kk)
            DZ = zfrac*(sigm(kk+1) - sigm(kk))
            if(level.le.zlid) then
                zvol = zvol + DZ
                meanval(n)=meanval(n)+fldvar(i,j,kk,n)*DZ
            end if
          end do
        enddo
        meanval(n) = meanval(n)/zvol
      enddo

      return
      end
