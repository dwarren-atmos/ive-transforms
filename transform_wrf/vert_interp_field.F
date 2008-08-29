      subroutine vert_interp_field(fin,m,n,kk,nt,zsfc,ifin,fout,
     &                        kkf,zcoord,missing,value,zin)

c
c***********************************************************************
c  Function: perform vertical interpolation on height, pressure or theta
c            surfaces
c  
c  INPUT -
c  fin:    input 3D sigma array to be interpolated
c  m:      number of grid points in the x direction
c  n:      number of grid points in the y direction
c  kk:     total number of sigma levels for mass variables
c  zsfc:   terrain height
c  ifin:   vertical interpolation method 
c          1 - no interpolation (sigma)
c          2 - interpolated to height levels
c          3 - interpolated to pressure levels
c          4 - interpolated to theta levels
c  kkf:    number of vertical levels for output field fout
c  zcoord: vertical coordinate values for output
c  fout:   3 D array contains interpolated output values
c  missing: defined underground missing values
c          1 - use the lowest sigma values
c          2 - user define missing values
c value:   values for underground points (if missing=2)
c
      implicit none

      integer i
      integer idbms
      integer ifin
      integer istat
      integer j
      integer k
      integer kk
      integer kkf
      integer klvls
      integer l
      integer len
      integer m
      integer mn
      integer n
      integer nt
      integer missing

      real factor
      real fin(m,n,kk,nt)
      real fout(m,n,kkf,nt)
      real zin(m,n,kk,nt)
      real rlev1
      real rlev2
      real value
      real zcoord(kkf)
      real zsfc(m,n)

      character*1   cfluid
      character*4   cparm
      character*6   cparm_new
      character*2   lvltyp

      mn=m*n
      if(ifin.eq.1.or.ifin.eq.5) then
c
c***********************************************************************
c    no interpolation, assign fin to be fout
c***********************************************************************
c
         fout(:,:,:,:)=fin(:,:,:,:)
c
c***********************************************************************
c  interpolate to height surfaces
c***********************************************************************
c
      else if(ifin.eq.2) then
        call z2zint(fin,fout,zin,zcoord,kk,kkf,mn,nt,missing,value)
c
c***********************************************************************
c          interpolate to pressure surfaces
c          pin:  total pressure
c***********************************************************************
      else if(ifin.eq.3) then
           call s2pint(fin,fout,zin,zcoord,kk,kkf,mn,nt,missing,value)
c
c***********************************************************************
c          interpolate to theta surface
c***********************************************************************
c
      else if(ifin.eq.4) then
           call z2zint(fin,fout,zin,zcoord,kk,kkf,mn,nt,missing,value)

      endif


      return
      end
