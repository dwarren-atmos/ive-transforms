#include "pointer.icl"
c
c   Calculates the difference between two time levels in a variable
c
      ive_ptr function time_diff(fld1,time1,time2,
     >                  ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names)
      implicit none
c
#include "default.icl"
c
      integer ndims,dims(4), time1,time2
      real stag(4),min(4),max(4),miss
      character*(*) fld1,dunits,disunits,dim_names(4)
c      
      integer strbeg,strend
      integer ix,iy,iz,it
      real ignore
      logical new
      ive_ptr var1
c
c Get wind fields
c
      var1 = 0
      var1 = getvar(fld1,ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names,new)
      ix=dims(1);iy=dims(2);iz=dims(3);it=dims(4)
      if (var1.eq.0) then
        write(6,*) 'Error:  cannot access variable fld1'
        return
      end if

      time_diff = 0
      time_diff = getmem(ix*iy*iz)
      if (time_diff.eq.0) then
        write(6,*) 'Error:  cannot allocate space for time_diff.'
        return
      end if

      call calc_timediff(%val(time_diff),%val(var1),ix,iy,iz,it
     &                  ,time1,time2,miss)

      dims(1)=ix;dims(2)=iy;dims(3)=iz;dims(4)=1
      dim_names(1)='nx';dim_names(2)='ny';
      dim_names(3)='nz';dim_names(4)='one';
      min(4)=0 ; max(4)=0
      dunits=dunits; disunits=disunits;

      return
      end
c
      subroutine calc_timediff(timediff,var1,ix,iy,iz,it,t1,t2,miss)

      integer t1,t2,ix,iy,iz,it,i,j,k
      real var1(ix,iy,iz,it),timediff(ix,iy,iz),miss

      timediff(:,:,:) = var1(:,:,:,t1) - var1(:,:,:,t2)

      return
      end

