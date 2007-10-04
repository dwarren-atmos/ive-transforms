c--------------------------------------------------------------------------
c     This file contains routines for the computation of means and extremes
c     $Id: means.f,v 1.1 1994/11/14 22:37:20 warren Exp $
c     history
c       $Log: means.f,v $
c Revision 1.1  1994/11/14  22:37:20  warren
c Christoph Schaer's European Model transforms.
c
c Revision 1.3  1994/09/28  08:06:09  schaer
c Neue Funktionen TSUM und TACCU (Chrigel Frei).
c
c Revision 1.2  1994/05/26  12:17:49  schaer
c Verbesserung Missing-Data Checking in DMEAN.
c
c Revision 1.1  1994/05/26  11:51:59  schaer
c Initial revision
c
c--------------------------------------------------------------------------

      subroutine tmean(a,m,nx,ny,nz,nt,misdat,sflag)
c     ==============================================
c     Computes time mean. Missing values in 'a' will only yield a missing value
c     in 'd' if more than 5% of the values required for the mean are missing.
c
c     if sflag=0, compute mean
c     if sflag=1, compute sum

c     declaration of parameters
      integer    nx,ny,nz,nt,sflag
      real       a(nx,ny,nz,nt),m(nx,ny,nz),misdat

c     declaration of variables
      integer    i,j,k,l,cnt
      real       sum

      print *,'Computing TMEAN for an array dimensioned',nx,ny,nz,nt
c      print *,'misdat=',misdat

c     start computation
      do i=1,nx
        do k=1,nz
          do j=1,ny
            sum=0.
            cnt=0
            do l=1,nt
              if ((a(i,j,k,l).ne.misdat).or.(misdat.eq.0.)) then
                sum=sum+a(i,j,k,l)
                cnt=cnt+1
              endif
            enddo
            if (sflag.eq.1) then
c             calculate the sum
              if (float(cnt)/float(nt).ge.0.95) then
                m(i,j,k)=sum
              else
                m(i,j,k)=misdat
              endif
            elseif (sflag.eq.0) then
c             calculate the mean
              if (float(cnt)/float(nt).ge.0.95) then
                m(i,j,k)=sum/float(cnt)
              else
                m(i,j,k)=misdat
              endif
            else
              print *,'invalid setting for sflag'
              return
            endif
          enddo
        enddo
      enddo
c      print *,'misdat=',misdat
      end


      subroutine taccu(a,tac,nx,ny,nz,nt,misdat)
c     ==========================================
c     Computes time accumulation tac of field a. 
c     Missing values in 'a' will yield a missing value in tac.


c     declaration of parameters
      integer    nx,ny,nz,nt
      real       a(nx,ny,nz,nt),tac(nx,ny,nz,nt),misdat

c     declaration of variables
      integer    i,j,k,l

      print *,'Computing TACCU for an array dimensioned',nx,ny,nz,nt
c      print *,'misdat=',misdat

c     check for time dependence of the field
      if (nt.le.1) then
        print *,'field to accumulate is not time dependent'
        return
      endif

c     initialise first time level
      do i=1,nx
        do k=1,nz
          do j=1,ny
            tac(i,j,k,1)=a(i,j,k,1)
          enddo
        enddo
      enddo

c     start computation
      do l=2,nt
        do i=1,nx
          do k=1,nz
            do j=1,ny
              if ((tac(i,j,k,l-1).ne.misdat).or.(misdat.eq.0.)) then
                if ((a(i,j,k,l).ne.misdat).or.(misdat.eq.0.)) then
                  tac(i,j,k,l) = tac(i,j,k,l-1)+a(i,j,k,l)
                else
                  tac(i,j,k,l) = misdat
                endif
              else
                tac(i,j,k,l) = misdat
              endif
            enddo
          enddo
        enddo
      enddo
c      print *,'misdat=',misdat
      end


      subroutine xmean(a,m,nx,ny,nz,nt,misdat)
c     ========================================
c     Computes zonal mean. Missing values in 'a' will only yield a missing value
c     in 'd' if more than 5% of the values required for the mean are missing.

c     declaration of parameters
      integer    nx,ny,nz,nt
      real       a(nx,ny,nz,nt),m(ny,nz,nt),misdat

c     declaration of variables
      integer    i,j,k,l,cnt
      real       sum

      print *,'Computing XMEAN for an array dimensioned',nx,ny,nz,nt

c     start computation
      do l=1,nt
        do k=1,nz
          do j=1,ny
            sum=0.
            cnt=0.
            do i=1,nx
              if ((a(i,j,k,l).ne.misdat).or.(misdat.eq.0.)) then
                sum=sum+a(i,j,k,l)
                cnt=cnt+1
              endif
            enddo
            if (float(cnt)/float(nx).ge.0.95) then
              m(j,k,l)=sum/float(cnt)
            else
              m(j,k,l)=misdat
            endif
          enddo
        enddo
      enddo
      end


      subroutine dmean (a,af,nx,ny,nz,nt,i1,i2,j1,j2,
     &    datmin,datmax,misdat,
     &    icheck,ac,nxc,nyc,nzc,ntc,misdatc)
c----------------------------------------------------------------------------
c     this subroutine computes the horizontal mean of a field in a previously
c     specified subdomain. The subdomain is specified by i1,i2,j1,j2.
c     If icheck=1, then only the gridpoints which satisfy ac(i,j,k,l)<>misdatc
c     will contribute to the horizontal mean. If icheck=0, ac is ignored.
c----------------------------------------------------------------------------

      include  'rotpol.icl'

c     argument declarations
      integer  nx,ny,nz,nt,i1,i2,j1,j2,icheck,nxc,nyc,nzc,ntc
      real     a(nx,ny,nz,nt),af(1,1,nz,nt),datmin(2),datmax(2),misdat
      real     ac(nxc,nyc,nzc,ntc),misdatc

c     variable declarations
      integer  ii,jj,k,l
      real     cnti,cnth,cnt,sum,lat
      logical  missing

      do l=1,nt
      do k=1,nz
        cnt=0.
        sum=0.
        do jj=j1,j2
          if (ny.gt.1) then
            lat=datmin(2)+real(jj-1)*(datmax(2)-datmin(2))/real(ny-1)
          else
            lat=datmin(2)
          endif
          cnti=cos(lat*zpir18)
          do ii=i1,i2
            cnth=cnti
            if ((ii.eq.i1).or.(ii.eq.i2)) cnth=cnth*0.5
            if ((jj.eq.j1).or.(jj.eq.j2)) cnth=cnth*0.5
            missing=((misdat.ne.0.).and.(a(ii,jj,k,l).eq.misdat))
            if ((icheck.eq.1).and.(misdatc.ne.0.))
     &        missing=missing.or.
     &                (ac(ii,jj,min0(k,nzc),min0(l,ntc)).eq.misdatc)
            if (.not.missing) then
              cnt=cnt+cnth
              sum=sum+cnth*a(ii,jj,k,l)
            endif
          enddo
        enddo
        if (cnt.gt.0.) then
          af(1,1,k,l)=sum/cnt
          print 100,k,l,af(1,1,k,l)
 100      format ('    k=',i3,'   l=',i3,'   domain-average=',f10.4)
        else
          if (misdat.eq.0.) misdat=misdatc
          af(1,1,k,l)=misdat
        endif

      enddo
      enddo
      end
  

      subroutine tmax(a,m,it1,it2,nx,ny,nz,nt,misdat)
c     ===============================================
c     Computes the temporal max of a field'

c     declaration of parameters
      integer    nx,ny,nz,nt,it1,it2
      real       a(nx,ny,nz,nt),m(nx,ny,nz),misdat

c     declaration of variables
      integer    i,j,k,l
      real       defval
      data       defval/-1.e15/

c     initialize field
      do k=1,nz
      do j=1,ny
      do i=1,nx
        m(i,j,k)=defval
      enddo
      enddo
      enddo

c     start computation
      print *,'TMAX mit l=',it1,it2
      print *,nx,ny,nz,nt
      do l=it1,it2
        print *,l
        do k=1,nz
        do j=1,ny
        do i=1,nx
          if ((a(i,j,k,l).ne.misdat).or.(misdat.eq.0.)) then
            m(i,j,k)=amax1(m(i,j,k),a(i,j,k,l))
          endif
        enddo
        enddo
        enddo
      enddo

c     clean up
      if (misdat.ne.0.) then
        do k=1,nz
        do j=1,ny
        do i=1,nx
          if (m(i,j,k).eq.defval) m(i,j,k)=misdat
        enddo
        enddo
        enddo
      endif
      end
