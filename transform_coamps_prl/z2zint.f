      subroutine z2zint(din,dout,zin,zout,kin,kout,len
     1  ,missing,value)

c rcs keywords: $RCSfile: z2zint.f,v $ 
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
c
      real din    (len,kin)
      real dout   (len,kout)
      real zin    (len,kin)
c      real zout   (kout)
      real zout   (100)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer ki
      integer ko
      integer missing
c
      real temp
      real value

c
c***********************************************************************
c          interpolate from sigma-z levels to constant z levels
c***********************************************************************
c
      do i=1,len
        do ko=1,kout
          if (zout(ko).ge.zin(i,1)) dout(i,ko)=din(i,1)
c          if (zout(ko).le.zin(i,kin)) dout(i,ko)=din(i,kin)
c***********************************************************************
c         assign under ground (missing) values
c***********************************************************************
          if (zout(ko).le.zin(i,kin)) then
            if(missing.eq.1) then
              dout(i,ko)=din(i,kin)
            else
              dout(i,ko)=value
            endif
          endif
c
        enddo
      enddo
c
      do ki=2,kin
        do ko=1,kout
          do i=1,len
            if ((zout(ko).lt.zin(i,ki-1)).and.(zout(ko).ge.zin(i,ki)))
     1       then
              temp=(zout(ko)-zin(i,ki-1))/(zin(i,ki)-zin(i,ki-1))
              dout(i,ko)=din(i,ki-1)+temp*(din(i,ki)-din(i,ki-1))
            endif
          enddo
        enddo
      enddo
c
c***********************************************************************
c
      return
      end
