      subroutine list

      integer numlines,k
      
      PARAMETER (numlines = 5 )
      
      character*(65) line(numlines)
      
      data line/
     >     'LIST            generates this list of fields',
     >     'D[fld:dir]      derivative of field ''fld'' in ',
     >     '                direction ''dir'' ',
     >     'OMEGA[U:V:W]    Calculate the vertical velocity on ',
     >     '                the untransformed grid, omega.'/
      
      write(*,*)
      write(*,*) '------------------------------------------------'
      write(*,*) 'List of Available Fields: '
      write(*,*)
      do k = 1,numlines
         write(*,*) line(k)
      enddo
      write(*,*) '------------------------------------------------'
      
      return
      end
