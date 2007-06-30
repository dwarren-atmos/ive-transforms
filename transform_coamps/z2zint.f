      subroutine z2zint(din,dout,zin,zout,kin,kout,len,nt
     1  ,missing,value)

c rcs keywords: $RCSfile: z2zint.f,v $ 
c               $Revision: 1.1.1.1 $ $Date: 2001/04/10 21:59:36 $
c
      implicit none

c***********************************************************************
c           parameters:
c***********************************************************************

      integer kin, kout, len, nt

      real din(len,kin,nt)
      real dout(len,kout,nt)
      real zin(len,kin,nt)
      real zout(kout)

c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************

      integer i,n
      integer ki
      integer ko
      integer missing

      real temp
      real value

c
c***********************************************************************
c          interpolate from sigma-z levels to constant z levels
c***********************************************************************
c
      do n=1,nt ; do i=1,len

        do ko=1,kout

          if(zout(ko).ge.zin(i,1,n) .and. zout(ko).le.zin(i,kin,n)) then
             dout(i,ko,n)=din(i,1,n)
          else
            if(missing.eq.1) then
              dout(i,ko,n)=din(i,1,n)
            else
              dout(i,ko,n)=value
            endif
            goto 65
          end if

          do ki=1,kin-1

            if((zout(ko).lt.zin(i,ki+1,n)).and.
     &        (zout(ko).ge.zin(i,ki,n))) then
              temp=(zout(ko)-zin(i,ki,n))/(zin(i,ki+1,n)-zin(i,ki,n))
              dout(i,ko,n)=din(i,ki,n)+temp*(din(i,ki+1,n)-din(i,ki,n))
              goto 65
            endif

          end do

65        continue
        end do
      end do ; end do
c
c***********************************************************************
c
      return
      end
