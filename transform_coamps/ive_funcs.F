
      subroutine translate_domain(delx,dely)

      use default_module

      real delx,dely

      integer k,ierr
      real*8 slope(4),intercept(4)
      real plmin_scaled(4),plmax_scaled(4)
      character*(80) dom_units(4),display_units(4),ctmp

 76   FORMAT(3x,'x offset:  ',F8.1,' km')
 77   FORMAT(3x,'y offset:  ',F8.1,' km')

      write(6,*)
      write(6,*) 'Translating domain:  '
      write(6,76) delx
      write(6,77) dely
      write(6,*)

      offset(1) = 1000.*delx
      offset(2) = 1000.*dely

      plmin(1) = plmin(1) + offset(1)
      plmax(1) = plmax(1) + offset(1)
      plmin(2) = plmin(2) + offset(2)
      plmax(2) = plmax(2) + offset(2)

      call setrarr("plmin",plmin,4,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem resetting domain limits.'
        stop
      end if
      call setrarr("plmax",plmax,4,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem resetting domain limits.'
        stop
      end if

      call getaarr("domain_units",dom_units,4,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting domain units.'
        stop
      end if
      call getaarr("domain_display_units",display_units,4,ierr) 
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting domain display units.'
        stop
      end if

      do k=1,4
        call convert_units(dom_units(k),display_units(k),slope(k),
     >                             intercept(k))
        plmin_scaled(k) = plmin(k)
        plmax_scaled(k) = plmax(k)
        call scale(plmin_scaled(k),1,slope(k),intercept(k),0.0)
        call scale(plmax_scaled(k),1,slope(k),intercept(k),0.0)
      enddo

      call setrarr("plmin_scaled",plmin_scaled,4,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem resetting scaled domain limits.' 
        stop
      end if
      call setrarr("plmax_scaled",plmax_scaled,4,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem resetting scaled domain limits.'
        stop
      end if

      call info_2d(plmin_scaled,plmax_scaled,4)

      do k=1,2
        zbot_min(k) = zbot_min(k) + offset(k)
        zbot_max(k) = zbot_max(k) + offset(k)
      enddo

      return
      end

*---------------------------------------------------------------------------

      subroutine trans_switch(iflg)

      use default_module

      integer iflg

      write(6,*)
      if (iflg.eq.1) then
        write(6,*) 'Vertical coordinate transformation on.'
        trans_on = 1
      else 
        write(6,*) 'Vertical coordinate transformation off.'
        trans_on = 0
      end if
      write(6,*)

      return
      end

