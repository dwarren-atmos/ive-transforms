ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Rotates the u and v velocity by counter-clockwise by angle theta
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "pointer.icl"

      ive_ptr function urot2(theta,varnm1,varnm2,ndims,dims,stag,
     &                       min,max,miss,dunits,disunits,dim_names)
      IMPLICIT NONE

#include "default.icl"

      integer :: ndims,dims(4)
      real :: theta,stag(4),min(4),max(4),miss
      character(len=*) :: varnm1,varnm2
      character(len=*) :: dunits,disunits,dim_names(4)
c
      integer iu,iv,ijk,i,nt
      real r_val
      logical new
      ive_ptr uvar,vvar,vvartmp
c
      urot2 = 0
      vvar = 0

      vvartmp = getvar(varnm2,ndims,dims,stag,min,max,miss,dunits,
     &                 disunits,dim_names,new)
      if (vvartmp.eq.0) then
        write(*,*) 'Error:  cannot get field v'
        return
      end if
c
      uvar = getvar(varnm1,ndims,dims,stag,min,max,miss,dunits,
     &               disunits,dim_names,new)
      if (uvar.eq.0) then
        write(*,*) 'Error:  cannot get field u'
        return
      end if
c
      nt = dims(4)
c
      vvar = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (vvar.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field vvar.' 
        return
      end if
c
c  Initialize vvar
c
      do i = 1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(vvar),i,1)
      enddo
c
c  V at U points
c
      iu = 0;
      iv = 1;
      call varatu(%val(vvar),%val(vvartmp),
     &            nx,ny,nz,nt,iu,iv)
c
c	call freemem(vvartmp)
c
      urot2 = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (urot2.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field urot2.' 
        return
      end if
c
c  Do Rotation Calculation
c
      iu = 1;
      iv = 0;
      call rot_calc('U',%val(urot2),%val(uvar),%val(vvar),theta,
     &              nx,ny,nz,nt,iu,iv)

      call freemem(vvar)

      return
      end
c
c ************************************************************
c
      ive_ptr function vrot2(theta,ndims,dims,stag,min,max,miss,
     &                        dunits,disunits,dim_names)
      implicit none

#include "default.icl"
c
c  Rotates the u and v velocity by counter-clockwise by angle theta
c
      integer ndims,dims(4)
      real theta,stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
c
      integer iu,iv,ijk,i,nt
      real r_val
      logical new
      ive_ptr uvar,vvar,vvartmp
c
      vrot2 = 0
      iu = 1 ; iv = 0
c
      vvartmp = getvar('v',ndims,dims,stag,min,max,miss,dunits,
     &                    disunits,dim_names,new)
      if (vvartmp.eq.0) then
        write(*,*) 'Error:  cannot get field v'
        return
      end if
c
      uvar = getvar('u',ndims,dims,stag,min,max,miss,dunits,
     &                    disunits,dim_names,new)
      if (uvar.eq.0) then
        write(*,*) 'Error:  cannot get field u'
        return
      end if
c
      nt = dims(4)

      vvar = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (vvar.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field vvar.' 
        return
      end if
c
c  Initialize vvar
c
      do i = 1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(vvar),i,1)
      enddo
c
c  V at U points
c
      iu = 0 ; iv = 1
      call varatu(%val(vvar),%val(vvartmp),
     &            nx,ny,nz,nt,iu,iv)
c
      vrot2 = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (vrot2.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field vrot2.' 
        return
      end if
c
c  Do Rotation Calculation
c
      iu = 1 ; iv = 0
      call rot_calc('V',%val(vrot2),%val(uvar),%val(vvar),theta,
     &              nx,ny,nz,nt,iu,iv)
c
      call freemem(vvar)
c
      return
      end
c
c ******************************************************
