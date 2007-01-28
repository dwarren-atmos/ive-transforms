      subroutine list

      integer numlines,k

      PARAMETER (numlines = 60 )

      character*(65) line(numlines)

      data line/
     >  'LIST:', 
     >  'Generates this list of user defined functions',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'GRID[which]:',   
     >  'Returns index space heights for each node ',
     >  '(e.g., a plot of GRID[p] with contour ',
     >  'interval of unity shows the heights of the',
     >  'thermodynamic nodes)',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'TRANS_OFF(ON):', 
     >  'Turns the vertical transform off (on)',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'D[string Var,string Dir]:',    
     >  'Derivative of Variable Var in direction Dir',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'DLL[string Var,string Dir]:',    
     >  'Derivative of Variable Var the N-S or E-W Dir',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'U(V)ROT[Deg]:',  
     >  'Rotate u-wind by angle Deg clockwise from x>0',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'N,N2:',
     >  'BF-Frequency and BF-Frequency squared',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'RHO:',
     >  'Total Density',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'OMEGA:',
     >  'vertical velocity on untransformed grid',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'INERP[string var:real plevel]',
     >  'Interpolates the field var to the pressure level plevel',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'GEOPOT[real plevel]',
     >  'Finds the geopotential height of the pressure level plevel',
     >  ' ',
     >  '############################################# ',
     >  ' ',
     >  'UG (VG)',
     >  'Calculates the U and V geostrophic winds'/


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
