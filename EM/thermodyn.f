c-----------------------------------------------------------------------
c     This file contains thermodynamic computations
c     $Id: thermodyn.f,v 1.1 1994/11/14 22:37:28 warren Exp $
c     history
c       $Log: thermodyn.f,v $
c Revision 1.1  1994/11/14  22:37:28  warren
c Christoph Schaer's European Model transforms.
c
c Revision 1.1  1994/05/26  11:51:59  schaer
c Initial revision
c
c-----------------------------------------------------------------------


      subroutine pottemp(pt,t,ie,je,ke,le,stagz)
c     ==========================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     pt(ie,je,ke,le),t(ie,je,ke,le),getps,stagz

c     variable declaration
      integer  i,j,k,l
      real     rdcp,tzero,psrf
      data     rdcp,tzero /0.286,273.15/

c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation of potential temperature
      do l=1,le
      do i=1,ie
      do j=1,je
        psrf=getps(i,j,l)
        do k=1,ke
c         distinction of staggered and unstaggered temperatuer-variable
          if (stagz.eq.-0.5) then
c           distinction of temperature in K and deg C
            if (t(i,j,k,l).lt.100.) then
              pt(i,j,k,l)=(t(i,j,k,l)+tzero)*( (1000./prlay(k))**rdcp )
            else
              pt(i,j,k,l)=t(i,j,k,l)*( (1000./prlay(k))**rdcp )
            endif
          else
c           distinction of temperature in K and deg C
            if (t(i,j,k,l).lt.100.) then
              pt(i,j,k,l)=(t(i,j,k,l)+tzero)*( (1000./prlev(k))**rdcp )
            else
              pt(i,j,k,l)=t(i,j,k,l)*( (1000./prlev(k))**rdcp )
            endif
          endif
        enddo
      enddo
      enddo
      enddo
      end

      subroutine equpot(ap,t,qd,ie,je,ke,le)
c     ======================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     ap(ie,je,ke,le),t(ie,je,ke,le),getps
      real     qd(ie,je,ke,le)
      
c     variable declaration
      integer  i,j,k,l
      real     rdcp,tzero,psrf
      data     rdcp,tzero /0.286,273.15/

c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation of potential temperature
      do l=1,le
        do i=1,ie
          do j=1,je
            psrf=getps(i,j,l)
            do k=1,ke
              ap(i,j,k,l) = (t(i,j,k,l)+tzero)*(1000./prlay(k))
     +           **(0.2854*(1.0-0.28*qd(i,j,k,l)))*exp(
     +           (3.376/(2840.0/(3.5*alog(t(i,j,k,l)+tzero)-alog(
     +           100.*prlay(k)*max(1.0E-10,qd(i,j,k,l))/(0.622+0.378*
     +           qd(i,j,k,l)))-0.1998)+55.0)-0.00254)*1.0E3*
     +           max(1.0E-10,qd(i,j,k,l))*(1.0+0.81*qd(i,j,k,l)))
            enddo
          enddo
        enddo
      enddo
      end

      subroutine zlayer(ap,t,z,ie,je,ke,le)
c     ======================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     ap(ie,je,ke,le),t(ie,je,ke,le),getps
      real     z(ie,je,ke,le)
      
c     variable declaration
      integer  i,j,k,l
      real     rdg,tzero,psrf
      data     rdg,tzero /29.271,273.15/

c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation of height of main levels
      do l=1,le
	do i=1,ie
	  do j=1,je
	    psrf=getps(i,j,l)
	    do k=1,ke-1
	      ap(i,j,k,l) = (z(i,j,k,l) - rdg*(t(i,j,k,l)+tzero)*
     +                       alog(prlay(k)/prlev(k)))/1000.
	    enddo
	    ap(i,j,ke,l) = (z(i,j,ke-1,l) + rdg*(t(i,j,ke,l)+tzero)*
     +                     alog(prlev(ke-1)/prlay(ke)))/1000.
	  enddo
	enddo
      enddo
      end

      subroutine relhum(rh,t,qd,qw,ie,je,ke,le)
c     =========================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     rh(ie,je,ke,le),t(ie,je,ke,le),getps
      real     qd(ie,je,ke,le),qw(ie,je,ke,le)

c     variable declaration
      integer  i,j,k,l
      real     rdcp,tzero,psrf
      real     b1,b2w,b3,b4w,r,rd,gqd,ge
      data     rdcp,tzero /0.286,273.15/
      data     b1,b2w,b3,b4w,r,rd /6.1078, 17.2693882, 273.16, 35.86,
     &                  287.05, 461.51/

     
c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation of potential temperature
      do l=1,le
        do i=1,ie
          do j=1,je
            psrf=getps(i,j,l)
            do k=1,ke
              ge = b1*exp(b2w*(t(i,j,k,l))/(t(i,j,k,l)+b3-b4w))
              gqd= r/rd*ge/(prlay(k)-(1.-r/rd)*ge)
              rh(i,j,k,l)=(qd(i,j,k,l)+qw(i,j,k,l))/gqd
            enddo
          enddo
        enddo
      enddo
      end

      subroutine pstonn(pnn,t,zb,ie,je,ke,le)
c     =======================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     pnn(ie,je,le),t(ie,je,ke,le),getps
      real     zb(ie,je)

c     variable declaration
      integer  i,j,l
      real     rdcp,tzero,psrf,r,g
      real     ztstar,zalpha,zt0
      data     rdcp,tzero,r,g /0.286,273.15,287.05,9.80665/
     
c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation of potential temperature
      do l=1,le
        do i=1,ie
          do j=1,je
            psrf=getps(i,j,l)
            if (zb(i,j).lt.1.) then
              pnn(i,j,l)=getps(i,j,l)
            else
              ztstar = (t(i,j,1,l)+tzero)*(1. + 0.0065*r/g*
     &                   (getps(i,j,l)/prlay(1) - 1.0))
              zalpha = 0.0065*r
              zt0    = ztstar + 0.0065*zb(i,j)
              if (zt0.gt.290.5) then
                if (ztstar.gt.290.5) then
                  zalpha = 0.0
                  ztstar = 0.5*(ztstar+290.5)
                else
                  zalpha = r*(290.5-ztstar)/zb(i,j)
                endif
              else if (ztstar.lt.255.) then
                ztstar = 0.5*(255.0+ztstar)
              endif
              pnn(i,j,l) = getps(i,j,l)* exp(g*zb(i,j)/(r*ztstar)*
     &               (1.0 - 0.5*(zalpha*zb(i,j)/(r*ztstar)) +
     &                    0.333*(zalpha*zb(i,j)/(r*ztstar))**2))
            endif
          enddo
        enddo
      enddo
      end

      subroutine zboden(zb,t,z,ie,je,ke,le)
c     =======================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     zb(ie,je),t(ie,je,ke,le),getps
      real     z(ie,je,ke,le),akbk

c     variable declaration
      integer  i,j
      real     tzero,r,g
      data     tzero,r,g /273.15,287.05,9.80665/
     
c     computation of surface height
      do j=1,je
        do i=1,ie
          akbk = aklev(1)/getps(i,j,1)+bklev(1)
          zb(i,j) = z(i,j,1,1) + alog(akbk)*r/g*
     +                (t(i,j,1,1)+tzero)
        enddo
      enddo
      end

      subroutine pressure(pr,stag3,ie,je,ke,le)
c     =========================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     pr(ie,je,ke,le),getps,stag3

c     variable declaration
      integer  i,j,k,l
      real     psrf

c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation pressure
      do l=1,le
        do i=1,ie
          do j=1,je
            psrf=getps(i,j,l)
            do k=1,ke
              if (stag3.eq.0.) then
                pr(i,j,k,l)=prlev(k)
              else
                pr(i,j,k,l)=prlay(k)
              endif
            enddo
          enddo
        enddo
      enddo
      end
