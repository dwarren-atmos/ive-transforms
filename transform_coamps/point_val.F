#include "pointer.icl"

      real function point_val(varnm,ipt,jpt,zlev,tpt
     &                            ndims,dims,stag,min,max,miss,
     &                            dunits,disunits,dim_names)
#include "default.icl"

      integer, intent(in) :: ipt,jpt,tpt
      real, intent(in) :: zlev
      character(len=*), intent(in) :: varnm

      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*) , intent(inout) :: dunits,disunits,
     &                                    dim_names(4)

      integer :: k
      logical :: new
      ive_ptr :: var


      var = getvar(trim(varnm),ndims,dims,stag,min,max,miss,dunits,
     &                    disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot get field ', trim(varnm)
        return
      end if
      
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
