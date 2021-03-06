c-----------------------------------------------------------------------
c     This file contains utilities
c     $Id: utilities.f,v 1.1 1994/11/14 22:37:31 warren Exp $
c     history
c       $Log: utilities.f,v $
c Revision 1.1  1994/11/14  22:37:31  warren
c Christoph Schaer's European Model transforms.
c
c Revision 1.1  1994/05/26  11:51:59  schaer
c Initial revision
c
c-----------------------------------------------------------------------


      subroutine exttime(art,ar,nx,ny,nz,nt,it)
c     =========================================
      integer nx,ny,nz,nt,it,i,j,k
      real    ar(nx,ny,nz,nt),art(nx,ny,nz)
      do i=1,nx
      do j=1,ny
      do k=1,nz
        art(i,j,k)=ar(i,j,k,it)
      enddo
      enddo
      enddo
      end


      subroutine extlevel(arl,ar,nx,ny,nz,nt,ik)
c     ==========================================
      integer nx,ny,nz,nt,ik,i,j,l
      real    ar(nx,ny,nz,nt),arl(nx,ny,nt)
      do i=1,nx
      do j=1,ny
      do l=1,nt
        arl(i,j,l)=ar(i,j,ik,l)
      enddo
      enddo
      enddo
      end


      real function str2nmb(chrs,ierr)
c---------------------------------------------------------------
c     This function takes the string chrs and converts it into a
c     real number. If the string is not a number, the routine
c     returns with str2numb=0 and ierr=1
c---------------------------------------------------------------

      character*(*) chrs
      integer       ierr

      integer       i,ibeg,iend,cp,ce,cf,strbeg,strend

c     Check whether string is a number
      ierr=0
      ibeg=strbeg(chrs)
      iend=strend(chrs)
      if (ibeg.gt.iend) return
      cp=0
      ce=0
      cf=1
      do i=ibeg,iend
c       check for + and -
        if ((cf.eq.1).and.
     &    ((chrs(i:i).eq.'-').or.(chrs(i:i).eq.'+'))) then
          cf=0
        else
          cf=0        
          if (chrs(i:i).eq.'.') cp=cp+1
c         check for format with exponents
          if ((chrs(i:i).eq.'E').or.(chrs(i:i).eq.'e')) then
            ce=ce+1
            cf=1
            if (cp.le.1) cp=1
          else
            if (((chrs(i:i).lt.'0').or.(chrs(i:i).gt.'9')).and.
     &        (chrs(i:i).ne.'.')) goto 100
          endif
        endif
      enddo

      if ((cp.gt.1).or.(ce.gt.1)) goto 100

c     convert the string into a number
      read (chrs,*,err=100) str2nmb
      return

c     error-exit
 100  continue
      str2nmb=0.
      ierr=1
      end


      integer function str2int(chrs,ierr)
c---------------------------------------------------------------
c     This function takes the string chrs and converts it into an
c     integer number. If the string is not an integer, the routine
c     returns with str2int=0 and ierr=1
c---------------------------------------------------------------
      real           r,str2nmb
      character*(*)  chrs
      integer        ierr
      str2int=0
      r=str2nmb(chrs,ierr)
      if (ierr.ne.0) return
      str2int=nint(r)
      if (abs(r-real(str2int)).lt.1.e-2) return
      str2int=0
      ierr=1
      end


      subroutine minmax(a,d,umin,umax,nx,ny,nz,nt)
c     ============================================
c     Overwrite upper-left and lower-right corner with the values
c     provided in umin and umax. In addition, the new value is bound into
c     the interval [umin,umax]

c     declaration of parameters
      integer    nx,ny,nz,nt
      real       a(nx,ny,nz,nt),d(nx,ny,nz,nt),umin,umax

c     declaration of variables
      integer    i,j,k,l

      print *,'entering minmax'

c     start computation
      do l=1,nt
        do k=1,nz
          do i=1,nx
            do j=1,ny
              d(i,j,k,l)=amax1(umin,amin1(umax,a(i,j,k,l)))
            enddo
          enddo
          d(1,ny,k,l)=umin
          d(nx,1,k,l)=umax
        enddo
      enddo
      end


      subroutine dostagx(a,s,nx,ny,nz,nt,misdat)
c     ==========================================
c     Computes an x-staggered version of array 'a'.

c     declaration of arguments
      integer    nx,ny,nz,nt
      real       a(nx+1,ny,nz,nt),s(nx,ny,nz,nt),misdat
c     declaration of variables
      integer    i,j,k,l

c     start computation
      do l=1,nt
        do k=1,nz
          do j=1,ny
            do i=1,nx
              if ( ((a(i,j,k,l).ne.misdat).and.(a(i+1,j,k,l).ne.misdat))
     &          .or. (misdat.eq.0.) ) then
                s(i,j,k,l)=0.5*(a(i,j,k,l)+a(i+1,j,k,l))
              else
                s(i,j,k,l)=misdat
              endif
            enddo
          enddo
        enddo
      enddo
      end




      subroutine remmisdat(ar,nx,ny,nz,nt,misdat)
c     ===========================================
c     this subroutine removes the missing data of a field by appropriate
c     vertical inter- and extrapolation

c     arguments
      integer   nx,ny,nz,nt
      real      ar(nx,ny,nz,nt),misdat

c     variables
      integer   i,j,k,l,kl,ku,kh
      logical   missing

      print *,nx,ny,nz,nt,misdat

      do i=1,nx
      do j=1,ny
      do l=1,nt
        kl=-1
        ku=-1
        missing=.false.
        do k=1,nz
          if (ar(i,j,k,l).ne.misdat) then
            if (.not.missing) then
              kl=k
            else
              ku=k
            endif
          else
            missing=.true.
            if (k.eq.1)  kl=0
            if (k.eq.nz) ku=nz+1
          endif
          if ((missing).and.(kl.ne.-1).and.(ku.ne.-1)) then
c           interpolate / extrapolate the values between kl and ku
            if ((kl.eq.0).and.(ku.le.nz)) then
              do kh=kl+1,ku-1
                ar(i,j,kh,l)=ar(i,j,ku,l)
              enddo
            else if ((kl.ge.1).and.(ku.eq.nz+1)) then
              do kh=kl+1,ku-1
                ar(i,j,kh,l)=ar(i,j,kl,l)
              enddo
            else if ((kl.ge.1).and.(ku.le.nz)) then
              do kh=kl+1,ku-1
                ar(i,j,kh,l)=
     &            (ar(i,j,kl,l)*real(ku-kh)+ar(i,j,ku,l)*real(kh-kl))
     &                    /real(ku-kl)
              enddo
            endif
            kl=ku
            ku=-1
            missing=.false.
          endif
        enddo
      enddo
      enddo
      enddo
      end


      subroutine blnkbnd(ar,misdat,nlev,nx,ny,nz,nt)
c     ==============================================
c     overwrites the lowevermost nlev levels with misdat
      integer   nx,ny,nz,nt,nlev,i,j,k,l
      real      ar(nx,ny,nz,nt),misdat
      do l=1,nt
      do k=1,nlev
      do j=1,ny
      do i=1,nx
        ar(i,j,k,l)=misdat
      enddo
      enddo
      enddo
      enddo
      end


      subroutine cregrid (fld1,ar,ie,je,datmin,datmax,misdat)
c     ======================================================
      include 'rotpol.icl'
c     creates grid-array
      integer       ie,je ,i,j, strbeg,strend
      character*(*) fld1
      character*(20) fld
      real          ar(ie,je),datmin(2),datmax(2),dx,dy,misdat
      real          lon,lat,lonr,latr

      dx=(datmax(1)-datmin(1))/real(ie-1)
      dy=(datmax(2)-datmin(2))/real(je-1)

      if ((fld1.eq.'COS').or.(fld1.eq.'SIN')) then
        fld='RPHI'
      else if ((fld1.eq.'COSR').or.(fld1.eq.'TAN')) then
        fld='RPHIM'
      else
        fld=fld1(1:strend(fld1))
      endif

      do i=1,ie
       do j=1,je
        if (fld.eq.'GRIDX') ar(i,j)=real(i)
        if (fld.eq.'GRIDY') ar(i,j)=real(j)
        if (fld.eq.'RLON')  ar(i,j)=datmin(1)+real(i-1)*dx
        if (fld.eq.'RLAT')  ar(i,j)=datmin(2)+real(j-1)*dy
        if (fld.eq.'RLAM')  ar(i,j)=(datmin(1)+real(i-1)*dx)*zpir18
	if ((fld.eq.'RPHI').or.(fld.eq.'RPHIM')) then
	  if (abs(pollam-1.).lt.0.1) then
	    ar(i,j)=0.
	  else
	    ar(i,j)=(datmin(2)+real(j-1)*dy)*zpir18
	  endif
	endif
        if ((fld.eq.'LAT').or.(fld.eq.'LON')
     &  .or.(fld.eq.'PHI').or.(fld.eq.'LAM')) then
          lonr=datmin(1)+real(i-1)*dx
          latr=datmin(2)+real(j-1)*dy
          call t_ph2ll(latr,lonr,lat,lon)
          if (fld.eq.'LON') ar(i,j)=lon
          if (fld.eq.'LAT') ar(i,j)=lat
          if (fld.eq.'LAM') ar(i,j)=lon*zpir18
	  if (fld.eq.'PHI') then
	    if (abs(pollam-1.).lt.0.1) then
	      ar(i,j)=polphi*zpir18
	    else
	      ar(i,j)=lat*zpir18
	    endif
	  endif
        endif
       enddo
c      insert missing data if too close to the pole
       if (fld.eq.'RPHIM') then
         if (abs(datmin(2)+90.).lt.1.e-3) then
           misdat=99.99
           ar(i,1)=misdat
         endif
         if (abs(datmax(2)-90.).lt.1.e-3) then
           misdat=99.99
           ar(i,je)=misdat
         endif
       endif
      enddo

      if ((fld1.eq.'COS').or.(fld1.eq.'SIN').or.
     &    (fld1.eq.'TAN').or.(fld1.eq.'COSR')) then
        do i=1,ie
        do j=1,je
          if (ar(i,j).ne.99.99) then
            if (fld1.eq.'COS') ar(i,j)=cos(ar(i,j))
            if (fld1.eq.'SIN') ar(i,j)=sin(ar(i,j))
            if (fld1.eq.'TAN') ar(i,j)=tan(ar(i,j))
            if (fld1.eq.'COSR') ar(i,j)=1./cos(ar(i,j))
          endif
        enddo 
        enddo
      endif
      end


      subroutine addmiss(ar,n,val,misdat)
c-----------------------------------------------------------------------
c     ersetzt die Werte 'val' in ar(n) mit misdat.
c-----------------------------------------------------------------------
      integer     n,i
      real        ar(n),misd,val,misdat,eps

c     define a missing data value
      if (misdat.eq.0.) then
        misd=999.99
      else
        misd=misdat
      endif

c     find the extremas of the array
      eps=1.e-4
      do i=1,n
        if (abs(ar(i)-val).lt.eps) then
          ar(i)=misd
          misdat=misd
        endif
      enddo
      end


      subroutine levelar(ar,nx,ny,nz,nt)
c     ==================================
      integer  nx,ny,nz,nt,i,j,k,l
      real     ar(nx,ny,nz,nt)
      do i=1,nx
       do j=1,ny
        do k=1,nz
         do l=1,nt
           ar(i,j,k,l)=float(k)
         end do
        end do
       end do
      end do
      end



      subroutine anadata (ar,nn,misdat)
c     =================================
      integer   nn
      real      ar(nn),misdat

      integer   powmin, powmax
      parameter (powmin=-10,powmax=10)
      integer   i,lrgpow,smlpow,power(powmin:powmax)
      integer   missing,undef,neg,pos,zero,oom,pmin,pmax
      real      absmin,absmax,minval,maxval,log10,exp10,x,ln10

c     define statement-function for logarithm and exponent with base 10
      log10(x)=log(x)/ln10
      exp10(x)=exp(ln10*x)

c     initialize constants
      ln10=log(10.)
      absmin=exp10(float(powmin))
      absmax=exp10(float(powmax+1))

c     initialize counters
      missing=0
      undef=0
      zero=0
      neg=0
      pos=0
      lrgpow=0
      smlpow=0
      do i=powmin,powmax
        power(i)=0
      enddo
      maxval=-1.e20
      minval= 1.e20

      do i=1,nn
        if ((ar(i).eq.misdat).and.(misdat.ne.0.)) then
          missing=missing+1
        else if ((ar(i).lt.1.).or.(ar(i).gt.-1.)) then
c         do evaluation of extremum
          minval=amin1(minval,ar(i))
          maxval=amax1(maxval,ar(i))
c         do zero/neg/pos classification
          if (ar(i).eq.0.) zero=zero+1
          if (ar(i).gt.0.) pos=pos+1
          if (ar(i).lt.0.) neg=neg+1
c         do order of magnitude statistics
          if (abs(ar(i)).lt.absmin) then
            smlpow=smlpow+1
          else if (abs(ar(i)).ge.absmax) then
            lrgpow=lrgpow+1
          else
            x=abs(ar(i))
            oom=int(log10(x)+1000.)-1000
            power(oom)=power(oom)+1
          endif
        else
c         it depends on the compiler and the compiling-options, whether an
c         undefined value is accounted for in this loop
          undef=undef+1
        endif
      enddo
           
c     printout
      print 800,'number of data-points:',nn
      print 800,'missing:              ',missing
      print 800,'undefined:            ',undef
      print 800,'negative values:      ',neg
      print 800,'zero     values:      ',zero
      print 800,'positive values:      ',pos
      print 801,'minimum:              ',minval
      print 801,'maximum:              ',maxval
      pmin=powmax
      pmax=powmin
      do i=powmin,powmax
        if (power(i).gt.0) then
          if (i.lt.pmin) pmin=i
          if (i.gt.pmax) pmax=i
        endif
      enddo
      if (pmax.ge.pmin) then
        print *,'order of magnitude:'
        pmin=max0(pmin-1,powmin)
        pmax=min0(pmax+1,powmax)
        if (pmin.eq.powmin) print 811,powmin,smlpow
        do i=pmin,pmax
          print 810,i,i+1,power(i)
        enddo
        if (pmax.eq.powmax) print 812,powmax+1,lrgpow
      endif
      
  800 format(a,i12)
  801 format(a,e12.3)
  810 format('1.E',i3,' <=|x|< 1.E',i3,': ',i12)
  811 format(   '         |x|< 1.E',i3,': ',i12)
  812 format('1.E',i3,' <=|x|:         '   ,i12)

      end




      logical function isfunc4(expr,fun,arg1,arg2,arg3,arg4)
c-----------------------------------------------------------------------
c     Purpose:
c        This function is identical to isfunc, except that 4 arguments
c        can be specified. It checks whether the expression 'expr' 
c        represents a function involving the function-name 'fun'. 
c        The arguments of the function must be delimited by (), [], 
c        or {}. Valid argument-delimiters are ',' and ':'. A maximum 
c        of 4 arguments can be specified. If less then 4
c        arguments are specified, the unspecified arguments are
c        returned as '#'. 
c     Arguments:
c        expr    chrs  input   the expression to be tested
c        fun     chrs  input   the function-name to be screend for
c        arg1    chrs  output  the first  argument
c        arg2    chrs  output  the second argument
c        arg3    chrs  output  the third  argument
c        arg4    chrs  output  the fourth  argument
c        isfunc4 logi  output  .true. if 'expr' is of the form
c                              fun(arg1,arg2,arg3), where either (), []
c                              or {} can be used to delimit the arguments,
c                              and either ',' or ':' to separate the
c                              arguments.
c     Examples:
c        isfunc4('ABS[expression]','ABS',..)   returns isfunc=.true.
c                                                     arg1='expression'
c                                                     arg2='#'
c                                                     arg3='#'
c        isfunc4('ABS(exp1,exp2)' ,'ABS',..)   returns isfunc=.true.
c        isfunc4('ABS()'          ,'ABS',..)   returns isfunc=.true.
c        isfunc4('ABS(exp1)*(c+d)','ABS',..)   returns isfunc=.false.
c        isfunc4('ABS'            ,'ABS',..)   returns isfunc=.false.
c        isfunc4('ABSOLUTELY'     ,'ABS',..)   returns isfunc=.false.
c-----------------------------------------------------------------------

c     parameter declaration
      integer   LENVAR
      parameter (LENVAR=80)

c     argument declaration
      character*(*)  expr,fun,arg1,arg2,arg3,arg4

c     variable declaration
      integer        strbeg,strend
      integer        ibeg,iend,fbeg,fend,flen,abeg,aend,alen1,alen2
      logical        err
      character*(LENVAR) hargs,harg1,harg2,harg3,harg4
      character*(1)  del

c     initially set result to .false.
      isfunc4=.false.

c     check for length of string
      ibeg=strbeg(expr)
      iend=strend(expr)
      fbeg=strbeg(fun)
      fend=strend(fun)
      flen=fend-fbeg+1
      if (iend-ibeg+1.lt.flen+2) return

c     check for presence of function-name
      if (expr(ibeg:ibeg+flen-1).ne.(fun(fbeg:fend))) return

c     check for consistent structure of arguments, alen1 is the length of
c     the arguments including enclosing (),{},[]. alen2 without enclosing
c     brackets.
      alen1=iend-(ibeg+flen)+1
      call bracket(expr(ibeg+flen:iend),abeg,aend,err)
      alen2=aend-abeg+1
      if ((err).or.(alen1.eq.alen2)) return

c     the expression 'expr' is a function involving 'fun'
      isfunc4=.true.

c     parse out arguments
      hargs=expr(ibeg+flen+1:iend-1)
      call fnddel(hargs,':,',harg3,harg4,del)
      if (del.eq.'#') then
c       the function has just 1 argument
        arg1=harg3
        arg2='#'
        arg3='#'
        arg4='#'
      else
c       the function has 2 or more arguments
        hargs=harg3
        call fnddel(hargs,':,',harg2,harg3,del)
        if (del.eq.'#') then
c         the function has 2 arguments
          arg1=harg2
          arg2=harg4
          arg3='#'
          arg4='#'
        else
c         the function has 3 or 4 arguments
          hargs=harg2
          call fnddel(hargs,':,',harg1,harg2,del)
          if (del.eq.'#') then
c           the function has 3 arguments
            arg1=harg1
            arg2=harg3
            arg3=harg4
            arg4='#'
          else
c           the function has 4 arguments
            arg1=harg1
            arg2=harg2
            arg3=harg3
            arg4=harg4
          endif
        endif
      endif

      end


      subroutine copylevel (ar2d,ar3d,ie,je,ke,le,k)
c------------------------------------------------------------------------------
c     arguments and local variables
      integer  ie,je,ke,le,i,j,k,l
      real     ar3d(ie,je,ke,le),ar2d(ie,je,le)
      real     rmin,rmax
      rmin= 1.e20
      rmax=-1.e20
    
      do i=1,ie
      do j=1,je
      do l=1,le
        ar3d(i,j,k,l)=ar2d(i,j,l)
        rmin=amin1(ar2d(i,j,l),rmin)
        rmax=amax1(ar2d(i,j,l),rmax)
      enddo
      enddo
      enddo
      end


      subroutine mask(ar,size,rnumb1,rnumb2,misdat)
c----------------------------------------------------------------------------
c     overwrites all values of ar with misdat, except those which
c     satisfy rnumb1<ar(i)<rnumb2. If misdat=0, then misdat will
c     be redefined
c----------------------------------------------------------------------------
      integer  size,i
      real     ar(size),rnumb1,rnumb2,misdat
      if (misdat.eq.0.) misdat=-999.99
      print *,size,rnumb1,rnumb2,misdat
      do i=1,size
        if ((ar(i).le.rnumb1).or.(ar(i).ge.rnumb2)) ar(i)=misdat
      enddo
      end


      subroutine timediff(date1,date2,diff)
C     =====================================
C     Calculates the time difference in hours for the two dates specified
C     by the two arrays date1 and date2. They are expected to contain the
C     following date information:
C     year      month   day     time    step.
C
C     date1     array specifying the first date
C     date2     array specifying the second date
C     diff      time differenc between date1 and date2 in hours
C
C     Warning:  ihdiff is equal to 0 for date1/2 = 880203_12_00 and
C               880202_12_24 !!!

      integer   date1(5),date2(5)
      integer   idays(12)       ! array containing the days of the monthes
      real      diff
      integer   ixday,ihdiff,iddiff,j

      data idays/31,28,31,30,31,30,31,31,30,31,30,31/

C     Determine if the period between date1 and date2 contains a Feb.29

      ixday=0   ! extra day flag

      if (((mod(date1(1),4).eq.0).and.(date1(1).ne.0)).or.
     >    ((mod(date2(1),4).eq.0).and.(date2(1).ne.0))) then
        if (((date1(2).le.2).and.(date2(2).gt.2)).or.
     >      ((date1(2).gt.2).and.(date2(2).le.2))) then
          ixday=1       ! set extra day flag to 1
          idays(2)=29   ! february has 29 days
        endif
      endif

      ihdiff=0  ! diff. in hours between date1/date2
      iddiff=0  ! diff. in days  between date1/date2

      if (date1(1).gt.date2(1)) then            ! compare years
        do j=date2(1),date1(1)-1
          iddiff=iddiff+365+ixday
        enddo
      else if (date1(1).lt.date2(1)) then
        do j=date1(1),date2(1)-1
          iddiff=iddiff-365-ixday
        enddo
      endif

      if (date1(2).gt.date2(2)) then            ! compare monthes
        do j=date2(2),date1(2)-1
          iddiff=iddiff+idays(j)
        enddo
      else if (date1(2).lt.date2(2)) then
        do j=date1(2),date2(2)-1
          iddiff=iddiff-idays(j)
        enddo
      endif

      iddiff=iddiff+date1(3)-date2(3)
      ihdiff=iddiff*24+date1(4)-date2(4)+date1(5)-date2(5)

      diff=real(ihdiff)

      return
      end




      integer function time2index(time)
c------------------------------------------------------------------------------
c     finds the time-index for a given time. If 'time' is not found in
c     the array '', 0 is returned.
c------------------------------------------------------------------------------
      include 'attributes.icl'

      real     time,eps
      integer  it

      if (ntimes.eq.0) then
        time2index=0.
        return
      endif

      eps=amax1(abs(timeval(ntimes)-timeval(1))/real(ntimes*100),1.e-5)
      time2index=0.
      do it=1,ntimes
        if (abs(timeval(it)-time).lt.eps) then
          time2index=it
          return
        endif
      enddo
      end




      subroutine initatt(ndims,dims, stag, datmin, datmax, misdat,
     &                               data_units, data_display_units)
c-------------------------------------------------------------------
c     Purpose: 
c        Initializes the attributes of a variable as required for
c        a thermodynamic variable. The information is taken from
c        the constants common block.
c-------------------------------------------------------------------
      include       'constants.icl'
      include       'attributes.icl'

c     argument declaration
      integer       MAXDIM
      parameter     (MAXDIM=4)
      integer       ndims, dims(MAXDIM)
      real          stag(MAXDIM), datmin(MAXDIM), datmax(MAXDIM),misdat
      character*80  data_units, data_display_units

      ndims=4
      dims(1)=nx
      dims(2)=ny
      dims(3)=nz
      dims(4)=max0(1,max0(ntimes,nt))
      stag(1)=0.
      stag(2)=0.
      stag(3)=0.
      misdat=0.
      datmin(1)=lonmin
      datmax(1)=lonmax
      datmin(2)=latmin
      datmax(2)=latmax
      datmin(3)=presssrf
      datmax(3)=presstop
      datmin(4)=timeval(1)
      datmax(4)=timeval(ntimes)
      data_units=" "
      data_display_units=" "
      end 
