c-----------------------------------------------------------------------
c     This file contains filter-routines
c     $Id: filter.f,v 1.1 1994/11/14 22:37:14 warren Exp $
c     history
c       $Log: filter.f,v $
c Revision 1.1  1994/11/14  22:37:14  warren
c Christoph Schaer's European Model transforms.
c
c Revision 1.1  1994/05/26  11:51:59  schaer
c Initial revision
c
c-----------------------------------------------------------------------


      subroutine filt4d (a,af,nx,ny,nz,nt,datmin,datmax,fil,misdat,ierr)
c     ==================================================================
c     this subroutine applies the diffusion-operator onto the horizontal
c     levels of a 4d-array

c     argument declarations
      integer  nx,ny,nz,nt,ierr
      real     a(nx,ny,nz*nt),af(nx,ny,nz*nt),fil,misdat
      real        datmin(2),datmax(2)

c     variable declarations
      integer  k,ptrf1,ptrf2,i,nfil
      integer  getmem
      real     rfil
      integer  iperx,ipery,inpol,ispol

      nfil=int(fil-0.00001)
      rfil=fil-float(nfil)

c     get work-memory
      ptrf1=getmem((nx+1)*ny)
      ptrf2=getmem(nx*(ny+1))
      if ((ptrf1.eq.0).or.(ptrf2.eq.0)) then
        ierr=1
        write (*,*) ' Cannot get memory in filt4d'
        return
      else
        ierr=0
      endif

c     initialize flags for special boundaries
      iperx=0 
      ipery=0
      ispol=0
      inpol=0
      if (abs(abs(datmax(1)-datmin(1))-360.).lt.0.01) iperx=1
      if (abs(datmin(2)+90.).lt.0.01) ispol=1
      if (abs(datmax(2)-90.).lt.0.01) inpol=1

c     call the 2d filter routine for every level
      do k=1,nz*nt
        print *,'Calling filt2d r for k=',k
        call filt2d(a(1,1,k),af(1,1,k),%val(ptrf1),%val(ptrf2),nx,ny,
     &                 rfil,misdat,iperx,ipery,ispol,inpol)
        do i=1,nfil
          print *,'Calling filt2d i for k=',k
          call filt2d(af(1,1,k),af(1,1,k),%val(ptrf1),%val(ptrf2),nx,ny,
     &                 1.  ,misdat,iperx,ipery,ispol,inpol)
        enddo
      enddo

c     return work-memory
      call freemem(ptrf1)
      call freemem(ptrf2)
      end


      subroutine filt2d (a,af,f1,f2,nx,ny,fil,misdat,
     &                                     iperx,ipery,ispol,inpol)
c     =============================================================
c     Apply a conservative diffusion operator onto the 2d field a,
c     with full missing data checking.
c
c     a     real   inp  array to be filtered, dimensioned (nx,ny)
c     af    real   out  filtered array, dimensioned (nx,ny), can be
c                       equivalenced with array a in the calling routine
c     f1    real        workarray, dimensioned (nx+1,ny)
c     f2    real        workarray, dimensioned (nx,ny+1)
c     fil   real   inp  filter-coeff., 0<afil<=1. Maximum filtering with afil=1
c                       corresponds to one application of the linear filter.
c     misdat real  inp  missing-data value, a(i,j)=misdat indicates that
c                       the corresponding value is not available. The
c                       misdat-checking can be switched off with with misdat=0.
c     iperx int    inp  periodic boundaries in the x-direction (1=yes,0=no)
c     ipery int    inp  periodic boundaries in the y-direction (1=yes,0=no)
c     inpol int    inp  northpole at j=ny  (1=yes,0=no)
c     ispol int    inp  southpole at j=1   (1=yes,0=no)
c
c     Christoph Schaer, 1993

c     argument declaration
      integer     nx,ny
      real        a(nx,ny),af(nx,ny),f1(nx+1,ny),f2(nx,ny+1),fil,misdat
      integer     iperx,ipery,inpol,ispol

c     local variable declaration
      integer     i,j,is
      real        fh

c     compute constant fh
      fh=0.125*fil

c     compute fluxes in x-direction
      if (misdat.eq.0.) then
        do j=1,ny
        do i=2,nx
          f1(i,j)=a(i-1,j)-a(i,j)
        enddo
        enddo
      else
        do j=1,ny
        do i=2,nx
          if ((a(i,j).eq.misdat).or.(a(i-1,j).eq.misdat)) then
            f1(i,j)=0.
          else
            f1(i,j)=a(i-1,j)-a(i,j)
          endif
        enddo
        enddo
      endif
      if (iperx.eq.1) then
c       do periodic boundaries in the x-direction
        do j=1,ny
          f1(1,j)=f1(nx,j)
          f1(nx+1,j)=f1(2,j)
        enddo
      else
c       set boundary-fluxes to zero
        do j=1,ny
          f1(1,j)=0.
          f1(nx+1,j)=0.
        enddo
      endif

c     compute fluxes in y-direction
      if (misdat.eq.0.) then
        do j=2,ny
        do i=1,nx
          f2(i,j)=a(i,j-1)-a(i,j)
        enddo
        enddo
      else
        do j=2,ny
        do i=1,nx
          if ((a(i,j).eq.misdat).or.(a(i,j-1).eq.misdat)) then
            f2(i,j)=0.
          else
            f2(i,j)=a(i,j-1)-a(i,j)
          endif
        enddo
        enddo
      endif
c     set boundary-fluxes to zero
      do i=1,nx
        f2(i,1)=0.
        f2(i,ny+1)=0.
      enddo
      if (ipery.eq.1) then
c       do periodic boundaries in the x-direction
        do i=1,nx
          f2(i,1)=f2(i,ny)
          f2(i,ny+1)=f2(i,2)
        enddo
      endif
      if (iperx.eq.1) then
        if (ispol.eq.1) then
c         do south-pole
          is=(nx-1)/2
          do i=1,nx
            f2(i,1)=-f2(mod(i-1+is,nx)+1,2)
          enddo
        endif
        if (inpol.eq.1) then
c         do north-pole
          is=(nx-1)/2
          do i=1,nx
            f2(i,ny+1)=-f2(mod(i-1+is,nx)+1,ny)
          enddo
        endif
      endif

c     compute flux-convergence -> filter
      if (misdat.eq.0.) then
        do j=1,ny
        do i=1,nx
            af(i,j)=a(i,j)+fh*(f1(i,j)-f1(i+1,j)+f2(i,j)-f2(i,j+1))
        enddo
        enddo
      else
        do j=1,ny
        do i=1,nx
          if (a(i,j).eq.misdat) then
            af(i,j)=misdat
          else
            af(i,j)=a(i,j)+fh*(f1(i,j)-f1(i+1,j)+f2(i,j)-f2(i,j+1))
          endif
        enddo
        enddo
      endif
      end

