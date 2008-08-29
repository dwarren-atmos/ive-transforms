        subroutine uvg2uv (u, v, kk, mn, rot,utru,vtru)
c
c****************************************************************** 
c arguments:
c
c        u:   u-component of winds in grid coordinate
c
c        v:   v-component of winds in grid coordinate 
c
c        mn:  number of points in the x and y plane
c
c        kk:  number of points in the z (or k) direction
c
c        utru: u-component of winds in world (lat/long) coordinate
c
c        vtru: v-component of winds in world (lat/long) coordinate
c
c  convert grid u and v to real u and v
c  assuming grid u = real u and grid v = real v along
c  the standard longitude and rot is the rotation array
c
c  Author: Pedro T.H. Tsai   30 Nov 1998
c
c********************************************************************
c
        implicit none
        integer mn, i, k, kk
        real u(mn,kk), v(mn,kk), rot(mn), utru(mn,kk), vtru(mn,kk)
        real r2d, d2r, ff, dd, ndd, uu, vv
  
        d2r = atan(1.0)/45.0
        r2d = 1.0/d2r
        do 10 k = 1, kk
          do 10 i = 1, mn
c  calc rotation angle between grid lon and standard lon
            uu   = u(i,k)
            vv   = v(i,k)
            ff   = sqrt((uu*uu) + (vv*vv))
            if (uu .eq. 0.0) uu = 1.0E-6
            dd   = 270.0 - (atan2(vv,uu)*r2d)
            dd   = amod(dd, 360.0)
            ndd  = dd - rot(i) + 360.0
            ndd  = amod(ndd, 360.0)
            utru(i,k) = -sin(ndd*d2r)*ff
            vtru(i,k) = -cos(ndd*d2r)*ff
10    continue
  
      return
      end
