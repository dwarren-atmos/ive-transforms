
#include "pointer.icl"

      ive_ptr function run_int(varnm,var0,int0,ndims,dims,stag,
     >                    min,max,miss,dunits,disunits,dim_names)
*
*  run_int computes the running time integral of variable varnm.
*    Here var0 specifies the value of the variable at time zero (zero 
*    by default) and int0 specifies the value of the integral
*    at time zero (also zero by default).  Note that var0 is only
*    used if min(t) is not zero.
*
      integer ndims,dims(4)
      real var0,int0,stag(4),min(4),max(4),miss
      character*(*) varnm,dunits,disunits,dim_names(4)

      integer k,t,strbeg,strend
      real dt
      character*(80) tmp
      logical new,units_err  
      ive_ptr getvar,getmem,var

      run_int = 0

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names,new)

      t = 0
      k = 1

 089  if ((t.eq.0).and.(k.le.4)) then
        if (dims(k).ne.1) then
          t = k
        else
          k = k+1
        end if
        goto 089
      end if

      if (t.eq.0) then
        write(*,*) 'Error:  cannot find non-trivial dimension.'
        return
      else if (t.lt.4) then 
        do k=t+1,4
          if (dims(k).ne.1) then
            write(*,*) 'Error:  variable ',
     >                 varnm(strbeg(varnm):strend(varnm)),
     >             ' has more than one non-trivial dimension.'
            return
          end if
        enddo
      end if

      dt = (max(t)-min(t))/(dims(t)-1)

      run_int = getmem(dims(t))
      if (run_int.eq.0) then
        write(*,*) 'Error:  could not allocate space for run_int.'
        return
      end if

 345  FORMAT(1x,'dimension ',I1,' with dt = ',F8.2)

      write(*,*)
      write(*,*) 'Calculating running integral of ',
     >         varnm(strbeg(varnm):strend(varnm)),' over '
      write(*,345) t,dt

      miss = .060458

      call run_int_calc(%val(run_int),%val(var),var0,int0,min(t),
     >                                dt,miss,dims(t))

      tmp = dunits(strbeg(dunits):strend(dunits))//'*s'
      call units_smpfy(tmp,dunits,units_err)
      if (units_err) dunits = tmp
      tmp = disunits(strbeg(disunits):strend(disunits))//'*s'
      call units_smpfy(tmp,disunits,units_err)
      if (units_err) disunits = tmp

      return
      end

**************************************************************************

      subroutine run_int_calc(run_int,var,var0,int0,t0,dt,miss,nt)

      integer nt
      real run_int(nt),var(nt),var0,int0,t0,dt,miss

      integer n,nstrt

      if (t0.eq.0) then
        run_int(1) = int0
        nstrt = 1
      else if (t0.gt.0) then
        run_int(1) = int0 + (var(1)+var0)*t0/2.
        nstrt = 1
      else 
        nstrt = 1
 777    if ((t0+(nstrt-1)*dt).le.(-dt/2.)) then
          run_int(nstrt) = miss
          nstrt = nstrt + 1
          goto 777
        end if
        run_int(nstrt) = int0
      end if

      do n=nstrt+1,nt
      
        run_int(n) = run_int(n-1) + (var(n)+var(n-1))*dt/2.

      enddo

      return
      end
     
     
      
