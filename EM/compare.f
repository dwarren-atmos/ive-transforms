c-----------------------------------------------------------------------
c     This file contains routines for the comparision of fields
c     $Id: compare.f,v 1.1 1994/11/14 22:36:57 warren Exp $
c     history
c       $Log: compare.f,v $
c Revision 1.1  1994/11/14  22:36:57  warren
c Christoph Schaer's European Model transforms.
c
c Revision 1.1  1994/05/26  11:51:59  schaer
c Initial revision
c
c-----------------------------------------------------------------------



      subroutine rmsdiff(a1,a2,af,nx,ny,nz,nt,i1,i2,j1,j2,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutines determines the root-mean-square difference
c        between two fields a1 and a2 on a subdomain (i1..i2,j1..j2) 
c        of the total horizontal extent (1..nx,1..ny). Missing data
c        is treated correctly.
c     Arguments:
c        a1  real arr  input   first field
c        a2  real arr  input   second field
c        af  real arr output   array containing rms-differences
c        nx   integer  input   first dimension of a1 and a2
c        ny   integer  input   second dimension of a1 and a2
c        nz   integer  input   number of vertical levels for a1 and a2
c        nt   integer  input   number of timesteps for a1 and a2
c        i1   integer  input   lower bound for first index 
c        i2   integer  input   upper bound for first index
c        j1   integer  input   lower bound for second index
c        j2   integer  input   upper bound for second index
c        misdat  real  input   missing data value
c-----------------------------------------------------------------------
      integer   nx, ny, nz, nt, i1, i2, j1, j2
      real      a1(nx,ny,nz,nt), a2(nx,ny,nz,nt), af(1,1,nz,nt)
      real      misdat

      integer   kk, ll
      real      rmserr

      do ll=1,nt
        do kk=1,nz
          af(1,1,kk,ll)=rmserr(a1(1,1,kk,ll),a2(1,1,kk,ll),
     &                          nx,ny,i1,i2,j1,j2,misdat)
        enddo
      enddo

      return
      end

      subroutine persdiff(a1,af,nx,ny,nz,nt,i1,i2,j1,j2,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutines determines the root-mean-square difference
c        between the starting point and any consecutive timestep
c        from a subdomain (i1..i2,j1..j2) from the field a1
c        with the total horizontal extent (1..nx,1..ny). Missing data
c        is treated correctly.
c     Arguments:
c        a1  real arr  input   first field
c        af  real arr output   array containing rms-differences
c        nx   integer  input   first dimension of a1 and a2
c        ny   integer  input   second dimension of a1 and a2
c        nz   integer  input   number of vertical levels for a1 and a2
c        nt   integer  input   number of timesteps for a1 and a2
c        i1   integer  input   lower bound for first index 
c        i2   integer  input   upper bound for first index
c        j1   integer  input   lower bound for second index
c        j2   integer  input   upper bound for second index
c        misdat  real  input   missing data value
c-----------------------------------------------------------------------
      integer   nx, ny, nz, nt, i1, i2, j1, j2
      real      a1(nx,ny,nz,nt), af(1,1,nz,nt)
      real      misdat

      integer   kk, ll
      real      rmserr

      do kk=1,nz
        af(1,1,kk,1)=0.
      enddo
      do ll=2,nt
        do kk=1,nz
          af(1,1,kk,ll)=rmserr(a1(1,1,kk,1),a1(1,1,kk,ll),
     &                          nx,ny,i1,i2,j1,j2,misdat)
        enddo
      enddo

      return
      end


      subroutine varblty(a1,af,nx,ny,nz,nt,i1,i2,j1,j2,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutines determines the root-mean-square difference
c        between the starting point and any consecutive timestep
c        from a subdomain (i1..i2,j1..j2) from the field a1
c        with the total horizontal extent (1..nx,1..ny). Missing data
c        is treated correctly.
c     Arguments:
c        a1  real arr  input   first field
c        af  real arr output   array containing rms-differences
c        nx   integer  input   first dimension of a1 and a2
c        ny   integer  input   second dimension of a1 and a2
c        nz   integer  input   number of vertical levels for a1 and a2
c        nt   integer  input   number of timesteps for a1 and a2
c        i1   integer  input   lower bound for first index 
c        i2   integer  input   upper bound for first index
c        j1   integer  input   lower bound for second index
c        j2   integer  input   upper bound for second index
c        misdat  real  input   missing data value
c-----------------------------------------------------------------------
      integer   nx, ny, nz, nt, i1, i2, j1, j2
      real      a1(nx,ny,nz,nt), af(1,1,nz,nt)
      real      misdat

      integer   kk, ll
      real      rmserr

      do kk=1,nz
        af(1,1,kk,1)=0.
      enddo
      do ll=2,nt
        do kk=1,nz
          af(1,1,kk,ll)=rmserr(a1(1,1,kk,ll-1),a1(1,1,kk,ll),
     &                          nx,ny,i1,i2,j1,j2,misdat)
        enddo
      enddo

      return
      end


      real function rmserr(a1,a2,nx,ny,i1,i2,j1,j2,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This function determines the root-mean-square difference
c        between two fields a1 and a2 on a subdomain (i1..i2,j1..j2) 
c        of the total horizontal extent (1..nx,1..ny). Missing data
c        is treated correctly.
c     Arguments:
c        a1  real arr  input   first field
c        a2  real arr  input   second field
c        nx   integer  input   first dimension of a1 and a2
c        ny   integer  input   second dimension of a1 and a2
c        i1   integer  input   lower bound for first index 
c        i2   integer  input   upper bound for first index
c        j1   integer  input   lower bound for second index
c        j2   integer  input   upper bound for second index
c        misdat  real  input   missing data value
c-----------------------------------------------------------------------

      integer  nx, ny, i1, i2, j1, j2
      real     a1(nx,ny), a2(nx,ny), misdat

      real     sum, cnt
      integer  ii, jj

      cnt=0.
      sum=0.
      do jj=j1,j2
        do ii=i1,i2
          if ((misdat.eq.0).or.((a1(ii,jj).ne.misdat).and.(
     &                           a2(ii,jj).ne.misdat))) then
            sum=sum+(a1(ii,jj)-a2(ii,jj))**2
            cnt=cnt+1.
          endif
        enddo
      enddo

      if (cnt.gt.0) then
        rmserr=sqrt(sum/cnt)
      else
        rmserr=misdat
      endif
 
      return
      end
