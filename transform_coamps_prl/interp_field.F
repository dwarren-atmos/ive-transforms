      subroutine interp_field(fin,m,n,kk,nt,zsfc,ifin
     1                 ,fout,kkf,zcoord,missing,value,zin)

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
c  inest:  domain id 
c  ihour:  model forecast hour
c  minute: model forecast minutes
c  isec:   model forecast second
c  dsetff: directory path name where coamps flat files are located
c  cdtg:   date-time-group of coamps forecast
c  kkf:    number of vertical levels for output field fout
c  zcoord: vertical coordinate values for output
c  fout:   3 D array contains interpolated output values
c  ifprint:
c          0 - minimum print out
c          1 - more diagnostic print out 
c missing: defined underground missing values
c          1 - use the lowest sigma values
c          2 - user define missing values
c value:   values for underground points (if missing=2)
c
      implicit none

      integer i
      integer idbms
      integer ifin
      integer ifprint
      integer ihour
      integer inest
      integer isec
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
      integer minute
      integer missing

      real factor
      real fin(m,n,kk,nt)
      real fout(m,n,kkf,nt)
      real zin(m,n,kk,nt)
      real rlev1
      real rlev2
      real value
      real zcoord(100)
      real zsfc(m,n)

      character*10  cdtg
      character*1   cfluid
      character*4   cparm
      character*6   cparm_new
      character*100 dsetff
      character*2   lvltyp


      mn=m*n
      if(ifin.eq.1.or.ifin.eq.5) then
c
c***********************************************************************
c    no interpolation, assign fin to be fout
c***********************************************************************
c
       do l=1,nt
         do k=1,kk
            do i=1,m
               do j=1,n
                  fout(i,j,k,l)=fin(i,j,k,l)
               enddo
            enddo
         enddo
       enddo
c
c***********************************************************************
c  interpolate to height surfaces
c***********************************************************************
c
      else if(ifin.eq.2) then
 
	 
        call z2zint(fin,fout,zin,zcoord,kk,kkf,mn,missing,value)
c
c***********************************************************************
c          interpolate to pressure surfaces
c          pin:  total pressure
c***********************************************************************
      else if(ifin.eq.3) then
           call s2pint(fin,fout,zin,zcoord,kk,kkf,mn,nt,missing,value)
c           call maxmin(fin,m,n,kkf)
c
c***********************************************************************
c          interpolate to theta surface
c***********************************************************************
c
      else if(ifin.eq.4) then
           call z2zint(fin,fout,zin,zcoord,kk,kkf,mn,missing,value)

      endif
	 

      return
      end
