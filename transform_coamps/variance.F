#include "pointer.icl"
c
c  Calculates the variance over the time dimension
c
      ive_ptr function variance(fldnam,ndims,dims,stag,min,max,miss
     &                         ,dunits,disunits,dim_names)

      use default_module
c
      IMPLICIT NONE
c
      integer :: ndims,dims(4)
      real :: stag(4),min(4),max(4),miss
      character(len=*) :: fldnam, dunits,disunits,dim_names(4)
      
      integer :: nt,mnk
      logical :: new
      ive_ptr fldval
c
c Get pressure field
c
      fldval = 0
      fldval = getvar(fldnam,ndims,dims,stag,min,max,miss,dunits,
     &                          disunits,dim_names,new)
      if (fldval.eq.0) then
        write(6,*) 'Error:  cannot access variable ',fldnam
        return
      end if
c
      mnk = dims(1)*dims(2)*dims(3)
      nt = dims(4)
c
      variance = 0
      variance = getmem(mnk)
      if (variance.eq.0) then
        write(6,*) 'Error:  cannot allocate space for geopot.'
        return
      end if
c
      call calc_variance(%val(variance),%val(fldval),miss,mnk,nt)
c
      dims(4) = 1
      dim_names(4) = 'one'
      min(4) = 0
      max(4) = 0
c
      return
      end
c
c
      subroutine calc_variance(var,fld,missing,mnk,nt)

      integer, intent(in) :: mnk, nt
      real, intent(in) :: fld(mnk,nt), missing
      real, intent(inout) :: var(mnk)
c
      real :: scr(mnk,nt), mean(mnk)
      integer :: n
c
      scr(:,:) = fld(:,:)
      mean(:) = 0.0

      ! calculate the mean and subtract
      do n=1,nt
        mean(:) = mean(:) + scr(:,n)/real(nt)        
      enddo

      do n=1,nt
        scr(:,n) = scr(:,n) - mean(:)
      enddo

      !calculate variance
      var(:) = 0.0
      do n=1,nt
        var(:) = var(:) + scr(:,n)*scr(:,n)/real(nt-1)
      enddo

      return
      end

