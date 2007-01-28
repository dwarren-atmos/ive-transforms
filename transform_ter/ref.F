
#include "pointer.icl"

      ive_ptr function ref(varnm,n,ndims,dims,stag,min,max,miss,
     >                           dunits,disunits,dim_names)
*
*  ref returns the field varnm evaluated at time (n-1)*dt [i.e., the
*    returned field has dimensions (dims(1),dims(2),dims(3),1)].
*
      integer ndims,dims(4),n
      real stag(4),min(4),max(4),miss
      character*(*) varnm,dunits,disunits,dim_names(4)

      integer strbeg,strend,nfix
      real tfix,r_val
      logical new
      ive_ptr var,getvar,getmem

      common/timeblk/time,ntime
      ive_ptr time
      integer ntime

      ref = 0

      if ((n.gt.ntime-1).or.(n.lt.0)) then
        write(*,*) 'Error:  specified time exceeds largest time.'
        return
      end if

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot get field ',
     >                varnm(strbeg(varnm):strend(varnm))
        return
      end if

      nfix = n+1
      tfix = r_val(%val(time),nfix) 

 14   FORMAT(5x,'at time ',F8.2)

      write(*,*)
      write(*,*) 'Calculating reference field from field ',
     >                varnm(strbeg(varnm):strend(varnm))
      write(*,14) tfix
      write(*,*)
     
      ref = getmem(dims(1)*dims(2)*dims(3))
      if (ref.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field ref.'
        return
      end if

      call ref_calc(%val(ref),%val(var),dims(1),dims(2),dims(3),
     >                         dims(4),nfix)

      dims(4) = 1
      stag(4) = 0.
      min(4) = tfix
      max(4) = tfix
      dim_names(4) = 'one'

      return
      end
  
*----------------------------------------------------------------

      subroutine ref_calc(ref,var,nx,ny,nz,nt,nfix)

      integer nx,ny,nz,nt,nfix
      real ref(nx,ny,nz),var(nx,ny,nz,nt)

      integer i,j,k

      do k=1,nz
      do j=1,ny
      do i=1,nx

        ref(i,j,k) = var(i,j,k,nfix)

      enddo
      enddo
      enddo

      return
      end

      
