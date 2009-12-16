      subroutine grdij(m,n,grdi,grdj)

c rcs keywords: $RCSfile: grdij.f,v $ 
c               $Revision: 1.1.1.1 $ $Date: 2001/04/10 21:59:35 $c
c***********************************************************************
c
      implicit none
c
c***********************************************************************
c
      integer i
      integer j
      integer m
      integer n
c
      real grdi   (m,n)
      real grdj   (m,n)
c
c***********************************************************************
c
      do j=1,n
        do i=1,m
          grdi(i,j)=float(i)
          grdj(i,j)=float(j)
        enddo
      enddo
c
c***********************************************************************
c
      return
      end
