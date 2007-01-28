	subroutine stagering_calc(iu,iv,iw,it,dims,plmin,min,dx,dy,twod)
c
	implicit none
c
c  Calculates what the stagering is
c
      integer iu,iv,iw,it,twod
	real dims(4),min(4),plmin(4),dx,dy
	logical requ
c
      it = 0
      iu = 0
      iv = 0
      iw = 0
      if ((dims(1).eq.1).and.(dims(2).eq.1).and.(dims(3).eq.1)) then 
        it = 1
      else if (requ(plmin(1),min(1)-dx/2.)) then
        iu = 1
      else if (requ(plmin(2),min(2)-dy/2.).and.(twod.ne.1)) then
        iv = 1
      else if (requ(plmin(3),min(3))) then
        iw = 1
      end if

	return
	end
