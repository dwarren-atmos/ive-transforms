#include "pointer.icl"

      ive_ptr function rho(base,ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names,bouss)
*
*  Calculates total density.  Kind of prosaic.
* 
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
      character therm0*5,therm*5

      logical new,base,bouss
      ive_ptr p,pprm,p0,thet,getvar,getmem

      rho = 0

      if(bouss) then
         therm = 'b'
         therm0 = 'b0'
      else
         therm = 'thet'
         therm0 = 'thet0'
      end if

      if (.not.base) then

        p0 = getvar('p0',ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
        if (p0.eq.0) then
          write(*,*) 'Error:  unable to access basic state pressure.'
          return
        end if

        pprm = getvar('p',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
        if (pprm.eq.0) then
          write(*,*) 'Error:  cannot access perturbation pressure.'
          return
        end if

        p = getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (p.eq.0) then
          write(*,*) 'Error:  unable to allocate space for ',
     >                          'total pressure.'
        end if

        call tot_press(%val(p),%val(pprm),%val(p0),dims(1),dims(2),
     >                               dims(3),dims(4))

      else 

        p = getvar('p0',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
        if (p.eq.0) then
          write(*,*) 'Error:  cannot find field p0.'
          return
        end if

      end if

      rho = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (rho.eq.0) then
        write(*,*) 'Error:  could not allocate space for rho.'
        if (.not.base) call freemem(p)
        return
      end if

      call rho_I(%val(rho),%val(p),dims(1),dims(2),dims(3),
     >              dims(4),bouss)

      if (.not.base) then

        call freemem(p)
        call driver('free = p',0)
        call driver('free = p0',0)

        thet = getvar(therm,ndims,dims,stag,min,max,miss,
     >                         dunits,disunits,dim_names,new)
        if (thet.eq.0) then
          write(*,*) 'Error:  cannot find field thet.'
          return
        end if

      else

        call driver('free = p0',0)

        thet = getvar(therm0,ndims,dims,stag,min,max,miss,
     >                         dunits,disunits,dim_names,new)
        if (thet.eq.0) then
          write(*,*) 'Error:  cannot find field thet0.'
          return
        end if

      end if

      call rho_II(%val(rho),%val(thet),dims(1),dims(2),dims(3),
     >                            dims(4))

      if (.not.base) then

        call driver('free = '//therm,0)

      else

        call driver('free = '//therm0,0)

      end if

  

      dunits = 'kg/m/m/m'
      disunits = 'kg/m/m/m'
      
      


      return
      end

*------------------------------------------------------------------------

      subroutine tot_press(p,pprm,p0,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real p(nx,ny,nz,nt),pprm(nx,ny,nz,nt),p0(nx,ny,nz)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        p(i,j,k,n) = ( pprm(i,j,k,n) + p0(i,j,k) )

      enddo
      enddo
      enddo
      enddo

      return
      end

*---------------------------------------------------------------------------

      subroutine rho_I(rho,p,nx,ny,nz,nt,bouss)

      integer nx,ny,nz,nt
      real rho(nx,ny,nz,nt),p(nx,ny,nz,nt)

      integer i,j,k,n
      real kappa,pref,g,thets,cp,R,rkappa

      logical bouss

      PARAMETER (pref = 100000., g=9.806, thets = 273., 
     >           cp = 1004., R = 287.04, kappa = R/cp, 
     >           rkappa = cp/R)

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        if(bouss) then
        rho(i,j,k,n) = g*pref/thets*(p(i,j,k,n)/(cp*thets))**(rkappa-1)
        else
        rho(i,j,k,n) = p(i,j,k,n)/((p(i,j,k,n)/pref)**kappa)
        end if

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

