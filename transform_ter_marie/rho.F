#include "pointer.icl"

      ive_ptr function rho(base,ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)
*
*  Calculates total density.  Kind of prosaic.
* 
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      logical new,base
      ive_ptr qv,p,thet,getvar,getmem

      rho = 0

      if (.not.base) then
      
        qv = getvar('qv',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
        if (qv.eq.0) then
          write(*,*) 'Warning:  cannot find field qv.  ',
     >                  'Calculating dry density.'
        end if

        p = getvar('p+p0',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
        if (p.eq.0) then
          write(*,*) 'Error:  cannot calculate total pressure.'
          return
        end if

        thet = getvar('thet',ndims,dims,stag,min,max,miss,
     >                         dunits,disunits,dim_names,new)
        if (thet.eq.0) then
          write(*,*) 'Error:  cannot find field thet.'
          return
        end if

      else 

        qv = getvar('qv0',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
        if (qv.eq.0) then
          write(*,*) 'Warning:  cannot find field qv0.  ',
     >                  'Calculating dry density.'
        end if

        p = getvar('p0',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
        if (p.eq.0) then
          write(*,*) 'Error:  cannot find field p0.'
          return
        end if

        thet = getvar('thet0',ndims,dims,stag,min,max,miss,
     >                         dunits,disunits,dim_names,new)
        if (thet.eq.0) then
          write(*,*) 'Error:  cannot find field thet0.'
          return
        end if

      end if

      rho = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (rho.eq.0) then
        write(*,*) 'Error:  could not allocate space for rho.'
        return
      end if

      dunits = 'kg/m/m/m'
      disunits = 'kg/m/m/m'

      if (qv.eq.0) then
        call rho_dry(%val(rho),%val(p),%val(thet),dims(1),dims(2),
     >                             dims(3),dims(4))
      else 
        call rho_moist(%val(rho),%val(qv),%val(p),%val(thet),dims(1),
     >                        dims(2),dims(3),dims(4))
      end if

      return
      end

**************************************************************************

      subroutine rho_dry(rho,p,thet,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real rho(nx,ny,nz,nt),p(nx,ny,nz,nt),thet(nx,ny,nz,nt)

      integer i,j,k,n
      real kappa,Rd

      PARAMETER (kappa = 287.04/1004., Rd = 287.04)

      do n = 1,nt
      do k = 1,nz
      do j = 1,ny
      do i = 1,nx

        rho(i,j,k,n) =  100.*p(i,j,k,n) / (Rd*thet(i,j,k,n)
     >                  *(p(i,j,k,n)/1000.)**kappa)

      enddo
      enddo
      enddo
      enddo

      return
      end

**************************************************************************

      subroutine rho_moist(rho,qv,p,thet,nx,ny,nz,nt)
*
*  Note:  if this doesn't work, check units on qv
*
      integer nx,ny,nz,nt
      real rho(nx,ny,nz,nt),p(nx,ny,nz,nt),thet(nx,ny,nz,nt)
      real qv(nx,ny,nz,nt)

      integer i,j,k,n
      real kappa,Rd

      PARAMETER (kappa = 287.04/1004., Rd = 287.04)

      do n = 1,nt
      do k = 1,nz
      do j = 1,ny
      do i = 1,nx

        rho(i,j,k,n) = 100.*p(i,j,k,n) / (Rd*thet(i,j,k,n)
     >                  *(p(i,j,k,n)/1000.)**kappa
     >                  *(1 + 0.61*qv(i,j,k,n)))

      enddo
      enddo
      enddo
      enddo

      return
      end

