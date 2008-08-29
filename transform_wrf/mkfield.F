
      character*(*) function mkfield(func,carg,iarg,rarg,numc,
     >                                   numi,numr,ndec)
*
*  mkfield constructs an IVE command string which can be used in a
*    call to getvar
*
      integer numc,numi,numr,iarg(numi),ndec
      real rarg(numr)
      character*(*) func,carg(numc)

      integer strbeg,strend,k
      character*(80) int2str,real2str,rstr,istr

      mkfield(1:len(mkfield)) = ' '
      mkfield = func(strbeg(func):strend(func))//'['

      do k=1,numc
        mkfield = mkfield(strbeg(mkfield):strend(mkfield))//
     >             carg(k)(strbeg(carg(k)):strend(carg(k)))
        mkfield = mkfield(strbeg(mkfield):strend(mkfield))//':'
      enddo

      do k=1,numr
        rstr = real2str(rarg(k),ndec)
        mkfield = mkfield(strbeg(mkfield):strend(mkfield))//
     >                     rstr(strbeg(rstr):strend(rstr))
        mkfield = mkfield(strbeg(mkfield):strend(mkfield))//':'
      enddo

      do k=1,numi
        istr = int2str(iarg(k))
        mkfield = mkfield(strbeg(mkfield):strend(mkfield))//
     >                istr(strbeg(istr):strend(istr))
        mkfield = mkfield(strbeg(mkfield):strend(mkfield))//':'
      enddo

      mkfield(strend(mkfield):strend(mkfield)) = ']'

      return
      end

*------------------------------------------------------------------------

      character*(*) function int2str(iarg)

      integer iarg

      integer ptr,k,fac,numplcs,num,intgr
      character*(1) chint,int2ch

      int2str(1:len(int2str)) = ' '
      ptr = 1

      if (iarg.lt.0) then
        int2str(1:1) = '-'
        ptr = ptr+1
      else if (iarg.eq.0) then
        int2str(1:1) = '0'
        return
      end if

      intgr = abs(iarg)
      numplcs = 1
      fac = 10

 89   if ((intgr/fac).ne.0) then
        numplcs = numplcs+1
        fac = 10*fac
        goto 89  
      end if

      do k=numplcs,2,-1
        fac = fac/10
        num = intgr/fac
        chint = int2ch(num)
        int2str(ptr:ptr) = chint(1:1)
        ptr = ptr+1
        intgr = intgr - num*fac
      enddo

      chint = int2ch(intgr)
      int2str(ptr:ptr) = chint(1:1)

      return
      end

*------------------------------------------------------------------------

      character*(*) function real2str(rarg,ndec)

      real rarg
      integer ndec

      real rl,fac
      integer ptr,k,numplcs,num
      character*(1) chint,int2ch

      real2str(1:len(real2str)) = ' '
      ptr = 1

      if (rarg.lt.0) then
        real2str(1:1) = '-'
        ptr = ptr+1
      else if ((rarg.lt.1.).and.(rarg.gt.(-1.))) then
        real2str(1:2) = '0.'
        ptr = ptr+2
      end if

      rl = abs(rarg)
      numplcs = 0
      fac = 1.

 79   if (int(rl/fac).ne.0) then
        numplcs = numplcs+1
        fac = 10.*fac
        goto 79
      end if

      do k=1,numplcs+ndec
        fac = fac/10.
        num = int(rl/fac)
        chint = int2ch(num)
        real2str(ptr:ptr) = chint(1:1)
        ptr = ptr+1
        if ((fac.gt.0.9).and.(fac.lt.(1.1))) then
          real2str(ptr:ptr) = '.'
          ptr = ptr+1
        end if
        rl = rl - num*fac
      enddo

      return
      end

*------------------------------------------------------------------------

      character*(*) function int2ch(iarg)

      integer iarg

      if (iarg.eq.0) then
        int2ch = '0'
      else if (iarg.eq.1) then
        int2ch = '1'
      else if (iarg.eq.2) then
        int2ch = '2'
      else if (iarg.eq.3) then
        int2ch = '3'
      else if (iarg.eq.4) then
        int2ch = '4'
      else if (iarg.eq.5) then
        int2ch = '5'
      else if (iarg.eq.6) then
        int2ch = '6'
      else if (iarg.eq.7) then
        int2ch = '7'
      else if (iarg.eq.8) then
        int2ch = '8'
      else if (iarg.eq.9) then
        int2ch = '9'
      else
        write(*,*) 'Error:  argument to int2ch not a single digit ',
     >                            'integer.'
        stop
      end if

      return
      end





