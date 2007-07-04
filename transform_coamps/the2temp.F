#include "pointer.icl"
c
c converts theta to temperature
c
      ive_ptr function the2temp(fld_the,fld_pre,
     >                  ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names)
      use default_module
c
      implicit none

      character(len=*), intent(in) :: fld_the, fld_pre
      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout) :: dunits,disunits,dim_names(4)
c      
      integer :: ix,iy,iz,it
      logical :: new
      ive_ptr :: ppp, the

c Get Theta field
      the = 0
      the = getvar(fld_the,ndims,dims,stag,min,max,miss,dunits,
     &                           disunits,dim_names,new)
      ix=dims(1);iy=dims(2);iz=dims(3);it=dims(4)
      if (the.eq.0) then
        write(6,*) 'Error:  cannot access variable ',fld_the
        return
      end if

c Get Pressure field
      ppp = 0
      ppp = getvar(fld_pre,ndims,dims,stag,min,max,miss,dunits,
     &                           disunits,dim_names,new)
      if (ppp.eq.0) then
        write(6,*) 'Error:  cannot access variable ',fld_pre
        return
      end if

      the2temp = 0
      the2temp = getmem(ix*iy*iz*it)
      if (the2temp.eq.0) then
        write(6,*) 'Error:  cannot allocate space for the2temp.'
        return
      end if

      call calc_the2temp(%val(the2temp),%val(the),%val(ppp),
     &                     ix,iy,iz,it,miss)

      dims(1)=ix; dims(2)=iy; dims(3)=iz; dims(4)=it
      dim_names(1)='nx';dim_names(2)='ny';
      dim_names(3)='nz';dim_names(4)='time';
c      min(4)= ; max(4)=0
      dunits='Kelvin'; disunits='Kelvin';

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  calc_the2temp
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_the2temp(temp,thet,pres,ix,iy,iz,it,miss)

      integer, intent(in) :: ix,iy,iz,it
      real, intent(in) :: thet(ix,iy,iz,it),pres(ix,iy,iz,it),miss
      real, intent(out) :: temp(ix,iy,iz,it)

      integer :: i,j,k,n
      real, parameter :: Rd=287, Cp=1004.6, p00=1000
      real, parameter :: kappa=Rd/Cp

      temp(:,:,:,:) = 0.0
      temp(:,:,:,:) = thet(:,:,:,:)*( pres(:,:,:,:)/p00 )**kappa

      return
      end
