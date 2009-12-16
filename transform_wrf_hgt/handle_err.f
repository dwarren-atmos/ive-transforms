      subroutine handle_err(str)
c
c   subroutine to print netcdf error
c
      implicit none

      integer len
      character*(*) str

      write (6,*) str
      write (6,100)
 100  format(' error, program stop')
      stop
      return
      end
