
c
c   Calculates the time mean of a field
c   (ASSUMES EVENLY SPACED TIME INTERVALS)
c
      integer(ive_ptr) function time_mean(fld,
     >                  ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names)
      use default_module
c
      implicit none

      character(len=*), intent(in) :: fld
      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)
c      
      integer :: ix,iy,iz,it
      logical :: new
      integer(ive_ptr) :: var
c
c Get wind fields
c
      print '(2A)','FLD ',fld
      var = 0
      var = getvar(fld,ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names,new)
      ix=dims(1);iy=dims(2);iz=dims(3);it=dims(4)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot access variable ',fld
        return
      end if

      time_mean = 0
      time_mean = getmem(ix*iy*iz)
      if (time_mean.eq.0) then
        write(6,*) 'Error:  cannot allocate space for time_diff.'
        return
      end if

      call calc_time_mean(%val(time_mean),%val(var),ix,iy,iz,it,miss)

      dims(1)=ix;dims(2)=iy;dims(3)=iz;dims(4)=1
      dim_names(1)='nx';dim_names(2)='ny';
      dim_names(3)='nz';dim_names(4)='one';
      min(4)=0 ; max(4)=0
      dunits=dunits; disunits=disunits;

      

      return
      end
c
      subroutine calc_time_mean(time_mean,var,ix,iy,iz,it,miss)

      integer, intent(in) :: ix,iy,iz,it
      real, intent(in) :: var(ix,iy,iz,it),miss
      real, intent(out) :: time_mean(ix,iy,iz)

      integer :: i,j,k,n

      time_mean(:,:,:) = 0.0
      do i=1,ix ; do j=1,iy ; do k=1,iz ; 
        do n=1,it
          if(var(i,j,k,n).eq.miss) then
            time_mean(i,j,k) = miss
            exit
          else
            time_mean(i,j,k) = time_mean(i,j,k) + var(i,j,k,n)/real(it)
          end if
        end do
      end do ; end do ; end do

      return
      end

