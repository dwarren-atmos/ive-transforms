c--------------------------------------------------------------------------
c     SUBROUTINES FOR COMPUTATION OF DERIVATIVES ON LON/LAT GRIDS
c--------------------------------------------------------------------------
c     der1d:   vertical derivative (with respect to 'P' or 'T')
c     der1dm:  dito, but missing data checking in addition
c     der2d:   horizontal derivative on data-surfaces (either 'X' or 'Y')
c     der2dm:  dito, but missing data checking in addition
c     der3d    horizontal derivative on sigma-surfaces
c
c     hasmiss: checks whether a field has missing data (logical function) 
c--------------------------------------------------------------------------



      subroutine der1dm(a,d,dir,ie,je,ke,le,
     &                                datmin,datmax,stag,misa,misd)
c--------------------------------------------------------------------------
c     Purpose: VERTICAL DERIVATIVE WITH MISSING DATA
c       Compute the vertical derivative with full missing data checking.
c       The derivative is taken from array 'a' in the direction of 'dir',
c       where 'dir' is either 'P','T'. The result is stored in array 'd'. 
c          The missing data values for 'a' and 'd' are 'misa' (input) and
c       'misd' (input/output), respectively. The flag misd is overwritten
c       with misd=0 on output if array a is free of missing data.
c           Depending on the value of 'derptype', 3 point weighted centered 
c       differencing (derptype=1) or simple centered differencing 
c       (dertype=2) is used.
c           The vertical level-structure of the data is of the form 
c       p=ak(k)+bk(k)*ps, where the surface pressure ps is obtained
c       through a call to the subroutine getps.
c     History:
c       dertype=2:  Christoph Schaer
c       dertype=1:  Andrea Rossa
c--------------------------------------------------------------------------

c     declaration of arguments
      integer       MAXDIM
      parameter     (MAXDIM=4)
      integer       ie,je,ke,le
      real          a(ie,je,ke,le),d(ie,je,ke,le)
      real          misa,misd
      real          datmin(MAXDIM),datmax(MAXDIM),stag(MAXDIM)
      character*(*) dir

c     common blocks for grid definition
      include 'constants.icl'

c     variable declaration
      integer       i,j,k,l
      real          dpu,dpl,quot,fac,psrf,getps
      real          ak(MAXK),bk(MAXK),as(MAXK),bs(MAXK)
      logical       hasmiss

c     catch cases without missing data
      if (.not.hasmiss(a,ie*je*ke*le,misa)) then
        call der1d(a,d,dir,ie,je,ke,le,datmin,datmax,stag)
        misd=0
        return
      endif

c     specify scaling factor associated with derivation with respect
c     to pressure
      if (dir.eq.'P') then
        fac=0.01
      else
        fac=1.
      endif

c     define ak and bk appropriate for current field 
      if (stag(3).eq.0.) then
        do k=1,nz
          ak(k)=aklev(k)
          bk(k)=bklev(k)
          as(k)=aslev(k)
          bs(k)=bslev(k)
        enddo
      else
        do k=1,nz
          ak(k)=aklay(k)
          bk(k)=bklay(k)
          as(k)=aslay(k)
          bs(k)=bslay(k)
        enddo
      endif

c     compute vertical 3 point derivative
      if (derptype.eq.1) then 
c       ---------------------------
c       3-point vertical derivative
c       ---------------------------
        do l=1,le
         do j=1,je
          do i=1,ie
c           get surface pressure at current grid-point
            psrf=getps(i,j,l)
c           points at k=1
            if ((a(i,j,1,l).eq.misa).or.(a(i,j,2,l).eq.misa)) then
               d(i,j,1,l)=misd
            else
              dpu=(ak(1)+bk(1)*psrf)-(ak(2)+bk(2)*psrf)
              d(i,j,1,l)=(a(i,j,1,l)-a(i,j,2,l))*fac/dpu  
            endif
c           points at 1<k<ke
            do k=2,ke-1
             if ((a(i,j,k-1,l).eq.misa).or.(a(i,j,k,l).eq.misa).or.
     &           (a(i,j,k+1,l).eq.misa)) then
               d(i,j,k,l)=misd
             else
               dpu       =(ak(k)+bk(k)*psrf)-(ak(k+1)+bk(k+1)*psrf)
               dpl       =(ak(k-1)+bk(k-1)*psrf)-(ak(k)+bk(k)*psrf)
               quot      =dpu/dpl
               d(i,j,k,l)=(quot*(a(i,j,k-1,l)-a(i,j,k,l))
     &           +1./quot*(a(i,j,k,l)-a(i,j,k+1,l)))*fac/(dpu+dpl)
             endif
           enddo
c          points at k=ke
           if ((a(i,j,ke-1,l).eq.misa).or.(a(i,j,ke,l).eq.misa)) then
             d(i,j,ke,l)=misd
           else
             dpl=(ak(ke-1)+bk(ke-1)*psrf)-(ak(ke)+bk(ke)*psrf)
             d(i,j,ke,l)=(a(i,j,ke-1,l)-a(i,j,ke,l))*fac/dpl 
           endif
          enddo
         enddo
        enddo
      else
c       ---------------------
c       centered differencing
c       ---------------------
        do l=1,le
         do  j=1,je
          do i=1,ie
c           get surface pressure at current grid-point
            psrf=getps(i,j,l)
c           points at 1<k<ke
            do k=2,ke-1
              if ((a(i,j,k-1,l).eq.misa).or.(a(i,j,k+1,l).eq.misa)) then
                d(i,j,k,l)=misd
              else
                d(i,j,k,l)=(a(i,j,k-1,l)-a(i,j,k+1,l))/2.
     &            /(as(k)+bs(k)*psrf)*fac
              endif
            enddo
c           points at k=1 and k=ke
            k=1
            if ((a(i,j,k,l).eq.misa).or.(a(i,j,k+1,l).eq.misa)) then
              d(i,j,k,l)=misd
            else
              d(i,j,k,l)=(a(i,j,k,l)-a(i,j,k+1,l))
     &           /(as(k)+bs(k)*psrf)*fac
            endif
            k=ke
            if ((a(i,j,k-1,l).eq.misa).or.(a(i,j,k,l).eq.misa)) then
              d(i,j,k,l)=misd
            else
              d(i,j,k,l)=(a(i,j,k-1,l)-a(i,j,k,l))
     &           /(as(k)+bs(k)*psrf)*fac
            endif
          enddo
         enddo
        enddo

      endif
      end


      subroutine der1d(a,d,dir,ie,je,ke,le,datmin,datmax,stag)
c--------------------------------------------------------------------------
c     Purpose: VERTICAL DERIVATIVE WITHOUT MISSING DATA
c       Compute the vertical derivative without missing data checking.
c       The derivative is taken from array 'a' in the direction of 'dir',
c       where 'dir' is either 'P','T'. The result is stored in array 'd'. 
c           Depending on the value of 'derptype', 3 point weighted centered 
c       differencing (derptype=1) or simple centered differencing 
c       (dertype=2) is used.
c           The vertical level-structure of the data is of the form 
c       p=ak(k)+bk(k)*ps, where the surface pressure ps is obtained
c       through a call to the subroutine getps.
c     History:
c       dertype=2:  Christoph Schaer
c       dertype=1:  Andrea Rossa
c--------------------------------------------------------------------------

c     declaration of arguments
      integer       MAXDIM
      parameter     (MAXDIM=4)
      integer       ie,je,ke,le
      real          a(ie,je,ke,le),d(ie,je,ke,le)
      real          datmin(MAXDIM),datmax(MAXDIM),stag(MAXDIM)
      character*(*) dir

c     common blocks for grid definition
      include 'constants.icl'

c     variable declaration
      integer       i,j,k,l
      real          dpu,dpl,quot,fac,psrf,getps
      real          ak(MAXK),bk(MAXK),as(MAXK),bs(MAXK)

c     specify scaling factor associated with derivation with respect
c     to pressure
      if (dir.eq.'P') then
        fac=0.01
      else
        fac=1.
      endif

c     define ak and bk appropriate for current field 
      if (stag(3).eq.0.) then
        do k=1,nz
          ak(k)=aklev(k)
          bk(k)=bklev(k)
          as(k)=aslev(k)
          bs(k)=bslev(k)
        enddo
      else
        do k=1,nz
          ak(k)=aklay(k)
          bk(k)=bklay(k)
          as(k)=aslay(k)
          bs(k)=bslay(k)
        enddo
      endif

c     compute vertical 3 point derivative
      if (derptype.eq.1) then 
c       ---------------------------
c       3-point vertical derivative
c       ---------------------------
        do l=1,le
         do j=1,je
          do i=1,ie
c           get surface pressure at current grid-point
            psrf=getps(i,j,l)
c           points at k=1
            dpu=(ak(1)+bk(1)*psrf)-(ak(2)+bk(2)*psrf)
            d(i,j,1,l)=(a(i,j,1,l)-a(i,j,2,l))*fac/dpu  
c           points at 1<k<ke
            do k=2,ke-1
              dpu       =(ak(k)+bk(k)*psrf)-(ak(k+1)+bk(k+1)*psrf)
              dpl       =(ak(k-1)+bk(k-1)*psrf)-(ak(k)+bk(k)*psrf)
              quot      =dpu/dpl
              d(i,j,k,l)=(quot*(a(i,j,k-1,l)-a(i,j,k,l))
     &           +1./quot*(a(i,j,k,l)-a(i,j,k+1,l)))*fac/(dpu+dpl)
            enddo
c           points at k=ke
            dpl=(ak(ke-1)+bk(ke-1)*psrf)-(ak(ke)+bk(ke)*psrf)
            d(i,j,ke,l)=(a(i,j,ke-1,l)-a(i,j,ke,l))*fac/dpl 
          enddo
         enddo
        enddo

      else
c       ---------------------
c       centered differencing
c       ---------------------
        do l=1,le
         do  j=1,je
          do i=1,ie
c           get surface pressure at current grid-point
            psrf=getps(i,j,l)
c           points at 1<k<ke
            do k=2,ke-1
              d(i,j,k,l)=(a(i,j,k-1,l)-a(i,j,k+1,l))/2.
     &           /(as(k)+bs(k)*psrf)*fac
            enddo
c           points at k=1 and k=ke
            k=1
            d(i,j,k,l)=(a(i,j,k,l)-a(i,j,k+1,l))
     &           /(as(k)+bs(k)*psrf)*fac
            k=ke
            d(i,j,k,l)=(a(i,j,k-1,l)-a(i,j,k,l))
     &           /(as(k)+bs(k)*psrf)*fac
          enddo
         enddo
        enddo

      endif
      end


      subroutine der2dm(a,d,cl,dir,ie,je,ke,datmin,datmax,misa,misd)
c-----------------------------------------------------------------------
c     Purpose: HORIZONTAL DERIVATIVE ON DATA-SURFACES WITH MISSING DATA
c       Compute the horizontal derivative with full missing data checking.
c       The derivative is taken from array 'a' in the direction of 'dir',
c       where 'dir' is either 'X','Y'. The result is stored in array 'd'. 
c          The missing data values for 'a' and 'd' are 'misa' (input) 
c       and 'misd' (input/output), respectively. The flag misd is
c       overwritten with misd=0 on output if array a is free of missing
c       data.
c          The routine accounts for derivatives at the pole and periodic
c       boundaries in the longitudinal direction (depending on
c       the value of datmin, datmax). If the data-set does not reach to 
c       the pole, a one-sided derivative is taken. Pole-treatment is only 
c       carried out if the data-set covers 360 deg in longitude, and it 
c       requires that ie=4*ii+1, where ii is an integer.
c-----------------------------------------------------------------------

c     declaration of arguments
      integer       MAXDIM
      parameter     (MAXDIM=4)
      integer       ie,je,ke
      real          a(ie,je,ke),d(ie,je,ke),cl(ie,je)
      real          misa,misd
      real          datmin(MAXDIM),datmax(MAXDIM)
      character*(*) dir

c     local variable declaration
      integer       i,j,k,ip1,im1,jp1,jm1,ip,im,j1,j2
      real          dlat,dlon,coslat,dx,dy,dxr,dyr
      integer       northpl,southpl,lonper
      logical       hasmiss

c     rerd and circ are the mean radius and diameter of the earth in meter
      real          rerd,circ,pi
      data          rerd,circ,pi / 6.37e6,4.e7,3.141592654/

c     catch cases without missing data
      if (.not.hasmiss(a,ie*je*ke,misa)) then
        call der2d(a,d,cl,dir,ie,je,ke,datmin,datmax)
        misd=0
        return
      endif

c     compute flags for pole and periodic treatment
      southpl=0
      northpl=0
      lonper =0
      j1=1
      j2=je
      if (abs(datmax(1)-datmin(1)-360.).lt.1.e-3) then
        lonper=1
        if (abs(datmin(2)+90.).lt.1.e-3) then
          southpl=1
          j1=2
        endif
        if (abs(datmax(2)-90.).lt.1.e-3) then
          northpl=1
          j2=je-1
        endif
      endif

      dlon=((datmax(1)-datmin(1))/float(ie-1)) *pi/180.
      dlat=((datmax(2)-datmin(2))/float(je-1)) *pi/180.

c      print *,'Computing derivative ',dir(1:1),
c     &        ' of an array dimensioned ',ie,je,ke

      if (dir(1:1).eq.'X') then
c       -----------------------------
c       derivation in the x-direction
c       -----------------------------
        do k=1,ke

c         do gridpoints at j1<=j<=j2
          do j=j1,j2
           coslat=cl(1,j)

c          do regular gridpoints at 1<i<ie, 1<j<je
           dx =rerd*coslat*dlon
           dxr=1./(2.*dx)
           do i=2,ie-1
            ip1=i+1
            im1=i-1
            if ((a(ip1,j,k).eq.misa).or.(a(im1,j,k).eq.misa)) then
              d(i,j,k)=misd
            else
              d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
            endif
           enddo ! i-loop
c          completed regular gridpoints at 1<i<ie, 1<j<je

c          do gridpoints at i=1, i=ie, 1<j<je
           if (lonper.eq.1) then
c           use periodic boundaries
            i=1
            ip1=2
            im1=ie-1
            if ((a(ip1,j,k).eq.misa).or.(a(im1,j,k).eq.misa)) then
              d(i,j,k)=misd
            else
              d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
            endif
            d(ie,j,k)=d(1,j,k)
           else
c           use one-sided derivatives
            dxr=1./dx
            i=1
            ip1=2
            im1=1
            if ((a(ip1,j,k).eq.misa).or.(a(im1,j,k).eq.misa)) then
              d(i,j,k)=misd
            else
              d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
            endif
            i=ie
            ip1=ie
            im1=ie-1
            if ((a(ip1,j,k).eq.misa).or.(a(im1,j,k).eq.misa)) then
              d(i,j,k)=misd
            else
              d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
            endif

           endif
c          completed gridpoints at i=1, i=ie, j1<=j<=j2

          enddo ! j-loop
c         completed gridpoints at 1<j<je

c         do gridpoints at j=je
          if (northpl.eq.1) then
c          for these gridpoints, the derivative in the x-direction is a
c          derivative in the y-direction at another pole-gridpoint
           dy =rerd*dlat
           dyr=1./(2.*dy)
           j=je
           jp1=je-1
           jm1=je-1
           do i=1,ie
            ip=mod(i-1+  (ie-1)/4,ie)+1
            im=mod(i-1+3*(ie-1)/4,ie)+1
            if ((a(ip,jp1,k).eq.misa).or.(a(im,jm1,k).eq.misa)) then
              d(i,j,k)=misd
            else
              d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
            endif
           enddo    ! i-loop
c          completed gridpoints at j=je
          endif
c         do gridpoints at j=1
          if (southpl.eq.1) then
           dy =rerd*dlat
           dyr=1./(2.*dy)
           j=1
           jp1=2
           jm1=2
           do i=1,ie
            ip=mod(i-1+  (ie-1)/4,ie)+1
            im=mod(i-1+3*(ie-1)/4,ie)+1
            if ((a(ip,jp1,k).eq.misa).or.(a(im,jm1,k).eq.misa)) then
              d(i,j,k)=misd
            else
              d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
            endif
           enddo    ! i-loop
          endif
c         completed gridpoints at j=1

        enddo    ! k-loop

      else if (dir(1:1).eq.'Y') then
c       -----------------------------
c       derivation in the y-direction
c       -----------------------------
        dy =dlat*rerd
        dyr=1./(2.*dy)
        do k=1,ke
         do i=1,ie

c          do regular gridpoints
           do j=2,je-1
            jp1=j+1
            jm1=j-1
            if ((a(i,jp1,k).eq.misa).or.(a(i,jm1,k).eq.misa)) then
              d(i,j,k)=misd
            else
              d(i,j,k)=dyr*(a(i,jp1,k)-a(i,jm1,k))
            endif
           enddo

c          do gridpoints at j=je
           if (northpl.eq.1) then
c           pole-treatment
            j=je
            jm1=j-1
            jp1=j-1
            ip=mod(i-1+(ie-1)/2,ie)+1
            im=i
            if ((a(ip,jp1,k).eq.misa).or.(a(im,jm1,k).eq.misa)) then
              d(i,j,k)=misd
            else
             d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
            endif
           else
c           one-sided derivative
            j=je
            jm1=je-1
            jp1=je
            if ((a(i,jp1,k).eq.misa).or.(a(i,jm1,k).eq.misa)) then
              d(i,j,k)=misd
            else
              d(i,j,k)=2.*dyr*(a(i,jp1,k)-a(i,jm1,k))
            endif
           endif
c          completed gridpoints at j=je

c          do gridpoints at j=1
           if (southpl.eq.1) then
c           pole-treatment
            j=1
            jm1=2
            jp1=2
            ip=i
            im=mod(i-1+(ie-1)/2,ie)+1
            if ((a(ip,jp1,k).eq.misa).or.(a(im,jm1,k).eq.misa)) then
              d(i,j,k)=misd
            else
              d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
            endif
           else
c           one-sided derivative
            j=1
            jm1=1
            jp1=2
            if ((a(i,jp1,k).eq.misa).or.(a(i,jm1,k).eq.misa)) then
              d(i,j,k)=misd
            else
              d(i,j,k)=2.*dyr*(a(i,jp1,k)-a(i,jm1,k))
            endif
           endif
c          completed gridpoints at j=1

         enddo
        enddo


      endif
      end




      subroutine der2d(a,d,cl,dir,ie,je,ke,datmin,datmax)
c-----------------------------------------------------------------------
c     Purpose: HORIZONTAL DERIVATIVE ON DATA-SURFACES WITHOUT MISSING DATA
c       Compute the horizontal derivative without missing data checking.
c       The derivative is taken from array 'a' in the direction of 'dir',
c       where 'dir' is either 'X','Y'. The result is stored in array 'd'. 
c          The routine accounts for derivatives at the pole and periodic
c       boundaries in the longitudinal direction (depending on
c       the value of datmin, datmax). If the data-set does not reach to 
c       the pole, a one-sided derivative is taken. Pole-treatment is only 
c       carried out if the data-set covers 360 deg in longitude, and it 
c       requires that ie=4*ii+1, where ii is an integer.
c-----------------------------------------------------------------------

c     declaration of arguments
      integer       MAXDIM
      parameter     (MAXDIM=4)
      integer       ie,je,ke
      real          a(ie,je,ke),d(ie,je,ke),cl(ie,je)
      real          datmin(MAXDIM),datmax(MAXDIM)
      character*(*) dir

c     local variable declaration
      integer       i,j,k,ip1,im1,jp1,jm1,ip,im,j1,j2
      real          dlat,dlon,coslat,dx,dy,dxr,dyr
      integer       northpl,southpl,lonper

c     rerd and circ are the mean radius and diameter of the earth in meter
      real          rerd,circ,pi
      data          rerd,circ,pi / 6.37e6,4.e7,3.141592654/

c     compute flags for pole and periodic treatment
      southpl=0
      northpl=0
      lonper =0
      j1=1
      j2=je
      if (abs(datmax(1)-datmin(1)-360.).lt.1.e-3) then
        lonper=1
        if (abs(datmin(2)+90.).lt.1.e-3) then
          southpl=1
          j1=2
        endif
        if (abs(datmax(2)-90.).lt.1.e-3) then
          northpl=1
          j2=je-1
        endif
      endif

      dlon=((datmax(1)-datmin(1))/float(ie-1)) *pi/180.
      dlat=((datmax(2)-datmin(2))/float(je-1)) *pi/180.

c      print *,'Computing derivative ',dir(1:1),
c     &        ' of an array dimensioned ',ie,je,ke

      if (dir(1:1).eq.'X') then
c       -----------------------------
c       derivation in the x-direction
c       -----------------------------
        do k=1,ke

c         do gridpoints at j1<=j<=j2
          do j=j1,j2
           coslat=cl(1,j)

c          do regular gridpoints at 1<i<ie, 1<j<je
           dx =rerd*coslat*dlon
           dxr=1./(2.*dx)
           do i=2,ie-1
            ip1=i+1
            im1=i-1
            d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
           enddo ! i-loop
c          completed regular gridpoints at 1<i<ie, 1<j<je

c          do gridpoints at i=1, i=ie, 1<j<je
           if (lonper.eq.1) then
c           use periodic boundaries
            i=1
            ip1=2
            im1=ie-1
            d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
            d(ie,j,k)=d(1,j,k)
           else
c           use one-sided derivatives
            dxr=1./dx
            i=1
            ip1=2
            im1=1
            d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
            i=ie
            ip1=ie
            im1=ie-1
            d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
           endif
c          completed gridpoints at i=1, i=ie, j1<=j<=j2

          enddo ! j-loop
c         completed gridpoints at 1<j<je

c         do gridpoints at j=je
          if (northpl.eq.1) then
c          for these gridpoints, the derivative in the x-direction is a
c          derivative in the y-direction at another pole-gridpoint
           dy =rerd*dlat
           dyr=1./(2.*dy)
           j=je
           jp1=je-1
           jm1=je-1
           do i=1,ie
            ip=mod(i-1+  (ie-1)/4,ie)+1
            im=mod(i-1+3*(ie-1)/4,ie)+1
            d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
           enddo    ! i-loop
c          completed gridpoints at j=je
          endif
c         do gridpoints at j=1
          if (southpl.eq.1) then
           dy =rerd*dlat
           dyr=1./(2.*dy)
           j=1
           jp1=2
           jm1=2
           do i=1,ie
            ip=mod(i-1+  (ie-1)/4,ie)+1
            im=mod(i-1+3*(ie-1)/4,ie)+1
            d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
           enddo    ! i-loop
          endif
c         completed gridpoints at j=1

        enddo    ! k-loop

      else if (dir(1:1).eq.'Y') then
c       -----------------------------
c       derivation in the y-direction
c       -----------------------------
        dy =dlat*rerd
        dyr=1./(2.*dy)
        do k=1,ke
         do i=1,ie

c          do regular gridpoints
           do j=2,je-1
            jp1=j+1
            jm1=j-1
            d(i,j,k)=dyr*(a(i,jp1,k)-a(i,jm1,k))
           enddo

c          do gridpoints at j=je
           if (northpl.eq.1) then
c           pole-treatment
            j=je
            jm1=j-1
            jp1=j-1
            ip=mod(i-1+(ie-1)/2,ie)+1
            im=i
            d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
           else
c           one-sided derivative
            j=je
            jm1=j-1
            jp1=j
            d(i,j,k)=2.*dyr*(a(i,jp1,k)-a(i,jm1,k))
           endif
c          completed gridpoints at j=je

c          do gridpoints at j=1
           if (southpl.eq.1) then
c           pole-treatment
            j=1
            jm1=2
            jp1=2
            ip=i
            im=mod(i-1+(ie-1)/2,ie)+1
            d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
           else
c           one-sided derivative
            j=1
            jm1=1
            jp1=2
            d(i,j,k)=2.*dyr*(a(i,jp1,k)-a(i,jm1,k))
           endif
c          completed gridpoints at j=1

         enddo
        enddo


      endif
      end



      subroutine der3d(a,d,cl,dir,ie,je,ke,le,ps,dps,datmin,datmax,stag)
c-----------------------------------------------------------------------
c     Purpose: HORIZONTAL DERIVATIVE ON PRESSURE-SURFACES WITHOUT MISSING DATA
c       The derivative is taken from array 'a' in the direction of 'dir',
c       where 'dir' is either 'X','Y'. The result is stored in array 'd'. 
c          The routine accounts for derivatives at the pole and periodic
c       boundaries in the longitudinal direction (depending on
c       the value of datmin, datmax). If the data-set does not reach to 
c       the pole, a one-sided derivative is taken. Pole-treatment is only 
c       carried out if the data-set covers 360 deg in longitude, and it 
c       requires that ie=4*ii+1, where ii is an integer.
c     History:
c       Daniel Luethi
c-----------------------------------------------------------------------

c     declaration of arguments
      integer       MAXDIM
      parameter     (MAXDIM=4)
      integer       ie,je,ke,le
      real          a(ie,je,ke,le),d(ie,je,ke,le),cl(ie,je)
      real          ps(ie,je,le),dps(ie,je,le)
      real          datmin(MAXDIM),datmax(MAXDIM),stag(MAXDIM)
      character*(*) dir

c     common blocks for grid definition
      include 'constants.icl'

c     variable declaration
      integer       i,j,k,l
      real          ak(MAXK),bk(MAXK),as(MAXK),bs(MAXK)

c     define ak and bk appropriate for current field 
      if (stag(3).eq.0.) then
        do k=1,nz
          ak(k)=aklev(k)
          bk(k)=bklev(k)
          as(k)=aslev(k)
          bs(k)=bslev(k)
        enddo
      else
        do k=1,nz
          ak(k)=aklay(k)
          bk(k)=bklay(k)
          as(k)=aslay(k)
          bs(k)=bslay(k)
        enddo
      endif

c     compute horizontal derivatives on sigma surfaces
      call der2d(a,d,cl,dir,ie,je,ke*le,datmin,datmax)

c     apply correction for horizontal derivative on p-surfaces
      do l=1,le
       do j=1,je
        do i=1,ie
          do k=2,ke-1
           d(i,j,k,l)=d(i,j,k,l)+bk(k)*dps(i,j,l)/2./(as(k)+
     &                 bs(k)*ps(i,j,l))*(a(i,j,k+1,l)-a(i,j,k-1,l))
          enddo
          d(i,j,1,l)=d(i,j,1,l)+bk(1)*dps(i,j,l)/(as(1)+
     &                 bs(1)*ps(i,j,l))*(a(i,j,2,l)-a(i,j,1,l))
          k=ke
          d(i,j,k,l)=d(i,j,k,l)+bk(k)*dps(i,j,l)/(as(k)+
     &                 bs(k)*ps(i,j,l))*(a(i,j,k,l)-a(i,j,k-1,l))
        enddo
       enddo
      enddo
      end


      logical function hasmiss(a,n,misdat) 
c--------------------------------------------------------------------------
c     Purpose:
c       Checks whether array a(n) has missing data. 
c--------------------------------------------------------------------------
      integer   n,i
      real      a(n),misdat
      hasmiss=.false.
      if (misdat.eq.0.) return
      do i=1,n
        if (a(i).eq.misdat) then
          hasmiss=.true.
          return
        endif
      enddo
      end
