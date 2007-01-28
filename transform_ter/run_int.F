
#include "pointer.icl"

      ive_ptr function run_int(varnm,istrt,ndims,dims,stag,
     >                    min,max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),istrt
      real stag(4),min(4),max(4),miss
      character*(*) varnm,dunits,disunits,dim_names(4)

      integer k,strbeg,strend
      real r_val
      character*(80) tmp
      logical new,err
      ive_ptr getmem,var

      run_int = 0

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  unable to access variable ',
     >          varnm(strbeg(varnm):strend(varnm))
        return
      end if

      if (dims(4).ne.ntime) then
        write(*,*) 'Error:  variable ',
     >               varnm(strbeg(varnm):strend(varnm)),
     >                 ' does not have all time levels.'
        return
      else if ((istrt.lt.0).or.(istrt.gt.(ntime-2))) then
        write(*,*) 'Error:  specified start time exceeds ',
     >                          'time domain.'
        return
      end if

      run_int = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (run_int.eq.0) then
        write(*,*) 'Error:  could not allocate space for run_int.'
        return
      end if

 991  FORMAT(5x,'starting at time ',F8.2)

      write(*,*)
      write(*,*) 'Calculating running integral of ',
     >         varnm(strbeg(varnm):strend(varnm))
      write(*,991) r_val(%val(time),istrt+1)
      write(*,*)

      miss = .060458

      call run_int_calc(%val(run_int),%val(var),%val(time),miss,
     >                  istrt,dims(1),dims(2),dims(3),dims(4))

      tmp = dunits(strbeg(dunits):strend(dunits))//'*s'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '

      tmp = disunits(strbeg(disunits):strend(disunits))//'*s'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' ' 

      return
      end

*-------------------------------------------------------------------------

      subroutine run_int_calc(run_int,var,time,miss,istrt,
     >                               nx,ny,nz,nt)

      integer nx,ny,nz,nt,istrt
      real run_int(nx,ny,nz,nt),var(nx,ny,nz,nt),time(nt),miss

      integer i,j,k,n

      do k=1,nz
      do j=1,ny
      do i=1,nx

        do n=1,istrt
          run_int(i,j,k,n) = miss
        enddo

        run_int(i,j,k,istrt+1) = 0.

        do n=istrt+2,nt
          run_int(i,j,k,n) = run_int(i,j,k,n-1) + 
     >         0.5*(var(i,j,k,n)+var(i,j,k,n-1))
     >                 *(time(n)-time(n-1))
        enddo

      enddo
      enddo 
      enddo

      return
      end

*---------------------------------------------------------------------

      ive_ptr function t_ave(varnm,istrt,istop,ndims,dims,stag,
     >              min,max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer istrt,istop,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),varnm

      integer k,strbeg,strend,iflg,isize
      real r_val
      logical new,err
      ive_ptr getmem,var

      t_ave = 0

      if ((istrt.lt.0).or.(istop.gt.ntime-1)) then
        write(*,*) 'Error:  specified times exceed time limits.'
        return
      else if (istrt.ge.istop) then
        write(*,*) 'Error:  stop time must exceed start time.'
        return
      end if

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                 disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  unable to access variable ',
     >                varnm(strbeg(varnm):strend(varnm))
        return
      end if

      if (dims(4).ne.ntime) then
        write(*,*) 'Error:  variable ',
     >                varnm(strbeg(varnm):strend(varnm)),
     >                    ' does not have all time levels.'
        return
      end if

      isize = dims(1)*dims(2)*dims(3)

      t_ave = getmem(isize)
      if (t_ave.eq.0) then
        write(*,*) 'Error:  could not allocate space for t_ave.'
        return
      end if

 591  FORMAT(5x,'start time ',F8.2)
 592  FORMAT(5x,'end time   ',F8.2)

      write(*,*)
      write(*,*) 'Calculating time average of ',
     >      varnm(strbeg(varnm):strend(varnm)),' with '
      write(*,591) r_val(%val(time),istrt+1) 
      write(*,592) r_val(%val(time),istop+1) 
      write(*,*)

      miss = .060458

      call t_ave_calc(%val(t_ave),%val(var),%val(time),istrt,istop,
     >                     dims(1),dims(2),dims(3),dims(4))

      k = 1
      iflg = 0
 599  if ((iflg.eq.0).and.(k.le.3)) then
        if (dims(k).ne.1) then
          iflg = 1
        else
          k = k+1
        end if
        goto 599 
      end if

      if (iflg.eq.0) then
        write(*,*) '   Average = ',r_val(%val(t_ave),1)
        write(*,*)
        call freemem(t_ave)
        t_ave = 0 
      else   
        dims(4) = 1
        stag(4) = 0.
        min(4) = 0.
        max(4) = 0.
        dim_names(4) = 'one'
      end if

      return
      end

*----------------------------------------------------------------------

      subroutine t_ave_calc(t_ave,var,time,istrt,istop,nx,ny,nz,nt)

      integer nx,ny,nz,nt,istrt,istop
      real t_ave(nx,ny,nz),var(nx,ny,nz,nt),time(nt)

      integer i,j,k,n
      real int,delt

      delt = time(istop+1)-time(istrt+1)

      do k=1,nz
      do j=1,ny
      do i=1,nx

        int = 0.
        do n=istrt+2,istop+1
          int = int + 0.5*(var(i,j,k,n)+var(i,j,k,n-1))
     >                  *(time(n)-time(n-1))/delt
        enddo

        t_ave(i,j,k) = int

      enddo
      enddo
      enddo

      return
      end 

      
