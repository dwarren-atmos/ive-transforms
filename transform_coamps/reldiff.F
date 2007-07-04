#include "pointer.icl"
c
c  Calculates the relative difference between field 1 and 2
c  reldiff = (fld1 - fld2)/fld2
c  if fld2 < ignore label missing
c
      ive_ptr function reldiff(fld1,fld2,ignore,
     >                  ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names)
      use default_module

      implicit none

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) fld1,fld2,dunits,disunits,dim_names(4)
c      
      integer strbeg,strend
      integer nx1,ny1,nz1,nt1,nx2,ny2,nz2,nt2,len
      real ignore
      logical new
      ive_ptr var1,var2
c
c Get wind fields
c
      var1 = 0; var2=0;
      var1 = getvar(fld1,ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names,new)
      nx1 = dims(1);ny1=dims(2);nz1=dims(3);nt1=dims(4)
      var2 = getvar(fld2,ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new)
      nx2 = dims(1);ny2=dims(2);nz2=dims(3);nt2=dims(4)
      if (var1.eq.0.or.var2.eq.0) then
        write(6,*) 'Error:  cannot access variable fld2'
        return
      end if

      if(nx1.ne.nx2.or.ny1.ne.ny2.or.nz1.ne.nz2.or.nt1.ne.nt2) then
        write(6,*) 'Error:  fld1 and 2 must be the same size'
        return
      end if

      reldiff = 0
      reldiff = getmem(nx1*ny1*nz1*nt1)
      if (reldiff.eq.0) then
        write(6,*) 'Error:  cannot allocate space for reldiff.'
        return
      end if

      len = nx1*ny1*nz1*nt1
      call calc_reldiff(%val(reldiff),%val(var1),%val(var2),ignore
     >                     ,len,miss)

      dims(1)=nx1;dims(2)=ny1;dims(3)=nz1;dims(4)=nt1
      dim_names(1)='nx';dim_names(2)='ny';
      dim_names(3)='nz';dim_names(4)='time';
      dunits=' '; disunits=' ';

      return
      end
c
      subroutine calc_reldiff(reldiff,var1,var2,ignore,len,miss)

      integer len,i
      real var1(len),var2(len),reldiff(len),miss,ignore
      real absdiff

c      print *,'ignore = ',ignore
      do i=1,len
       absdiff = abs(var1(i) - var2(i))
c       if(abs(var2(i)).lt.ignore.or.absdiff.lt.ignore) then
c       if(abs(var2(i)).lt.1.0.or.absdiff.lt.ignore) then
       if(abs(var2(i)).lt.5.0) then
c         reldiff(i) = miss
         reldiff(i) = 0.0
c         print *,i,' missing'
       else
c         reldiff(i)=(abs(var1(i))-abs(var2(i)))/var2(i)
         reldiff(i)=abs(var1(i))/var2(i)
       end if
      end do

      return
      end

