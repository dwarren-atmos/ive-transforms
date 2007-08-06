      MODULE MISC_MODULE

      USE DEFAULT_MODULE

      integer, parameter :: nxmax=501, nymax=501, nzmax=501
      real :: rdzc(nzmax),rdze(nzmax),xlbnd(nxmax),xrbnd(nxmax)
     &       ,ylbnd(nymax),yrbnd(nymax),zlbnd(nzmax),zrbnd(nzmax)

      CONTAINS

      subroutine translate_time(time_shift)

      real, intent(in) :: time_shift

      integer :: ierr

      plmin(4) = plmin(4) + time_shift
      plmax(4) = plmax(4) + time_shift

      print *,plmin(4),plmax(4),time_shift

      call setrarr("plmin",plmin,4,ierr) ; if(ierr) goto 25
      call setrarr("plmax",plmin,4,ierr) ; if(ierr) goto 30
      return

25    print *, 'Problems with serarr(plmin)'
      return

30    print *, 'Problems with serarr(plmax)'
      return

      end subroutine
      END MODULE
