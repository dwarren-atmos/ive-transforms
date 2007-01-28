
#include "pointer.icl"

      ive_ptr function spctrm(fld,ndims,dims,stag,min,max,miss,
     >                      dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),fld

      integer strbeg,strend
      logical new,requ,err
      character*(80) tmp
      ive_ptr getmem,var,wrk,cwrkx

      spctrm = 0

      var = getvar(fld,ndims,dims,stag,min,max,miss,dunits,disunits, 
     >                                    dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  could not access field ',
     >                 fld(strbeg(fld):strend(fld))
        return
      end if
    
      if (mod(dims(1),2).ne.0) then
        write(6,*) 'Error:  test version.  Assumes nx even.'
        write(6,*) 'Stopping.'
        stop 
      end if

      wrk = getmem(2*dims(1))
      if (wrk.eq.0) then
        write(6,*) 'Error:  could not allocate space for wrk.'
        return
      end if

      cwrkx = getmem(4*dims(1)+15)
      if (cwrkx.eq.0) then
        write(6,*) 'Error: could not allocate space for cwrkx.'
        call freemem(wrk)
        return
      end if

      spctrm = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (spctrm.eq.0) then
        write(6,*) 'Error:  could not allocate space for spctrm.'
        call freemem(wrk)
        call freemem(cwrkx)
        return
      end if

      call spctrm_calc(%val(spctrm),%val(var),%val(wrk),
     >            %val(cwrkx),dims(1),dims(2),dims(3),dims(4))

      call freemem(wrk)
      call freemem(cwrkx)

      end

*--------------------------------------------------------------------------

      subroutine spctrm_calc(spctrm,var,wrk,cwrkx,nx,ny,nz,nt) 

      integer nx,ny,nz,nt
      real spctrm(nx,ny,nz,nt),var(nx,ny,nz,nt),cwrkx(4*nx+15)
      complex wrk(nx)

      integer i,j,k,n,qx
      real a,b,amp,phase

      qx = nx/2
 
      call cffti(nx,cwrkx)

      do n=1,nt
      do k=1,nz
      do j=1,ny

        do i=1,nx
          wrk(i) = cmplx(var(i,j,k,n),0.)
        enddo
        call cfftf(nx,wrk,cwrkx)
        do i=1,nx
          wrk(i) = wrk(i)/nx
        enddo

        spctrm(1,j,k,n) = real(wrk(1))

        do i=2,qx+1
          a = real(wrk(i))
          b = aimag(wrk(i))
          amp = 2*sqrt(a**2+b**2)
          phase = atan(b/a)
          spctrm(i,j,k,n) = amp
          if (i.lt.qx+1) then
            spctrm(i+qx,j,k,n) = phase
          end if
        enddo

      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

c     logical function requ(val1,val2)

c     real val1,val2,eps,diff,big	

c     parameter (eps=1.e-4)

c     big = amax1(abs(val1),abs(val2))

c     if (big.eq.0.) then
c       requ = .true.
c       return
c     end if

c     diff = abs((val1-val2)/big)

c     if (diff.le.eps) then
c       requ = .true.
c     else
c       requ = .false.
c     end if

c     return
c     end

