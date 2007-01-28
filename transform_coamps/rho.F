
#include "pointer.icl"

      ive_ptr function rho_calc(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)
	implicit none
c
c  Calculates total density.  Kind of prosaic.
c 
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
      logical new
      ive_ptr p,thet,getvar,getmem
c
      p = getvar('ppp',ndims,dims,stag,min,max,miss,dunits,
     >          disunits,dim_names,new)
c
      if (p.eq.0) then
         write(*,*) 'Error:  cannot find field ppp.'
         return
      end if
c
      rho_calc = 0
      rho_calc = getmem(dims(1)*dims(2)*dims(3)*dims(4))
c
      if (rho_calc.eq.0) then
        write(*,*) 'Error:  could not allocate space for rho.'
        return
      end if
c
      call rho_I(%val(rho_calc),%val(p),dims(1),dims(2),dims(3),
     >              dims(4))
c
      call driver('free = PPP',0)
c
      thet = getvar('THE',ndims,dims,stag,min,max,miss,
     >                         dunits,disunits,dim_names,new)
c
      if (thet.eq.0) then
        write(*,*) 'Error:  cannot find field thet.'
        return
      end if
c
      call rho_II(%val(rho_calc),%val(thet),dims(1),dims(2),dims(3),
     >                            dims(4))
c
      call driver('free = THE',0)
c
      dunits = 'kg/m/m/m'
      disunits = 'kg/m/m/m'
c
      return
      end

*---------------------------------------------------------------------------

      subroutine rho_I(rho,p,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real rho(nx,ny,nz,nt),p(nx,ny,nz,nt)

      integer i,j,k,n
      real kappa,pref,g,thets,cp,R,rkappa,ppascal

      PARAMETER (pref = 100000., g=9.806, thets = 273., 
     >           cp = 1004., R = 287.04, kappa = R/cp, 
     >           rkappa = cp/R)

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

	  ppascal = 100*p(i,j,k,n)
        rho(i,j,k,n) = ppascal/((ppascal/pref)**kappa)

      enddo
      enddo
      enddo
      enddo

      return
      end

*--------------------------------------------------------------------------

      subroutine rho_II(rho,thet,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real rho(nx,ny,nz,nt),thet(nx,ny,nz,nt)

      integer i,j,k,n
      real Rd

      PARAMETER (Rd = 287.04)

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx
         
        rho(i,j,k,n) = rho(i,j,k,n)/(Rd*thet(i,j,k,n))

      enddo
      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

