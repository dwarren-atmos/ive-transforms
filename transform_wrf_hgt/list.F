      subroutine list

      integer numlines,k
      
      PARAMETER (numlines = 6 )
      
      character*(65) line(numlines)
      
      data line/
     >     'LIST            generates this list of fields',
     >     'D[fld:dir]      derivative of field ''fld'' in ',
     >     '                direction ''dir'' ',
     >     'VER_MASFLX      Calculate the vertical mass flux',
     >     'HOR_MASSFLX     Calculate the horizontal mass flux',
     >     'RHO             Calculate the total density' /
      
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
