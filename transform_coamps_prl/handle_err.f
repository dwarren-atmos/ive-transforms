      subroutine handle_err(str)

      implicit none

      character(len=*) str

      write (6,'(A,2x,A)') 'ERROR:',str
      write (6,100)
 100  format('Program stoping')

      stop
      end
