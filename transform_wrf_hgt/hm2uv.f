      subroutine hm2uv(hxm,hym,hxu,hyu,hxv,hyv,m,n)

c rcs keywords: $RCSfile: hm2uv.f,v $ 
c               $Revision: 1.1.1.1 $ $Date: 2001/04/10 21:59:35 $c
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer m
      integer n
c
      real hxm    (m,n)
      real hxu    (m,n)
      real hxv    (m,n)
      real hym    (m,n)
      real hyu    (m,n)
      real hyv    (m,n)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer j
      integer m1
      integer n1
c
c***********************************************************************
c          local constants
c***********************************************************************
c
      m1=m-1
      n1=n-1
c
c***********************************************************************
c          convert the map factor in x and y direction at mass points 
c          to the map factor in x and y direction at u wind points
c***********************************************************************
c
      do j=1,n
        do i=1,m1
          hxu(i,j)=(hxm(i,j)+hxm(i+1,j))*0.5
          hyu(i,j)=(hym(i,j)+hym(i+1,j))*0.5
        enddo
      enddo
      do j=1,n
        hxu(m,j)=hxu(m1,j)
        hyu(m,j)=hyu(m1,j)
      enddo
c
c***********************************************************************
c          convert the map factor in x and y direction at mass points 
c          to the map factor in x and y direction at v wind points
c***********************************************************************
c
      do j=1,n1
        do i=1,m
          hxv(i,j)=(hxm(i,j)+hxm(i,j+1))*0.5
          hyv(i,j)=(hym(i,j)+hym(i,j+1))*0.5
        enddo
      enddo
      do i=1,m
        hxv(i,n)=hxv(i,n1)
        hyv(i,n)=hyv(i,n1)
      enddo
c
c***********************************************************************
c
      return
      end
