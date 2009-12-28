      subroutine s2pint(din,dout,zin,zout,kin,kout,len,nt,
     >                  missing,value)

c rcs keywords: $RCSfile: s2pint.f,v $ 
c               $Revision: 1.1.1.1 $ $Date: 2001/04/10 21:59:36 $
c
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer kin
      integer kout
      integer len
      integer nt
      integer n
      integer missing
c
      real din (len,kin,nt)
      real dout (len,kout,nt)
      real value
      real zin (len,kin,nt)
      real zout (kout)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer ki
      integer ko
c
      real temp
c
c***********************************************************************
c          end of definitions
c***********************************************************************
c***********************************************************************
c          interpolate from sigma levels to p-levels
c          set up interpolation loops
c***********************************************************************
c
      do n=1,nt
      do ko=1,kout
        do i=1,len
c          if (zout(ko).le.zin(i,1,n)) dout(i,ko,n)=din(i,1,n)
c***********************************************************************
c         assign under ground (missing) values
c***********************************************************************
c          if (zout(ko).ge.zin(i,kin,n).or.zout(ko).le.zin(i,1,n))
          if (zout(ko).le.zin(i,kin,n).or.zout(ko).ge.zin(i,1,n))
     >    then
            if(missing.eq.1) then
              dout(i,ko,n)=din(i,kin,n)
            else
              dout(i,ko,n)=value
            endif
          endif

        enddo
      enddo
      do ki=2,kin
        do ko=1,kout
          do i=1,len
            if ((zout(ko).lt.zin(i,ki-1,n)).and.
     1          (zout(ko).ge.zin(i,ki,n)))
     2       then
              temp=(alog(zout(ko))-alog(zin(i,ki-1,n)))
     1            /(alog(zin(i,ki,n))-alog(zin(i,ki-1,n)))
              dout(i,ko,n)=din(i,ki-1,n)+
     1                    temp*(din(i,ki,n)-din(i,ki-1,n))
            endif
          enddo
        enddo
      enddo

      enddo
c
c***********************************************************************
c
      return
      end
