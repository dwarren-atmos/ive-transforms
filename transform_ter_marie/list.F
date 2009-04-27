

      subroutine list

      integer numlines,k

      PARAMETER (numlines = 51 )

      character*(65) line(numlines)

      data line/
     >  'LIST            generates this list of fields',
     >  'D[fld:dir]      derivative of field ''fld'' in ',
     >  '                direction ''dir'' ',
     >  'GRID[which]     returns index space heights for each node ',
     >  '                (e.g., a plot of GRID[p] with contour ',
     >  '                interval of unity shows the heights of the',
     >  '                thermodynamic nodes)',
     >  'N               bouyancy frequency',
     >  'DEFCMP[i1,i2]   i1,i2 component of the deformation tensor',
     >  'DEF2            square of deformation tensor',
     >  'RI_INV          inverse Richardson number',
     >  'KM_IVE          mixing coefficient',
     >  'RHO, RHO0       total density, basic state total density',
     >  'DRAG[cmp,i1,i2,j1,j2]  computes the ''cmp'' component of the',
     >  '                pressure drag over the domain defined',
     >  '                by i1*dx to i2*dx, j1*dy to j2*dy',
     >  'PER_VOL[fld]    takes the product of field ''fld'' and RHO',
     >  'PER_VOL_BASE[fld] takes the product of field ''fld'' ',
     >  '                and RHO0',
     >  'PERT[fld,base]  takes the perturbation of field ''fld'' from',
     >  '                the basic state field ''base'' ',
     >  'VOL_INT[fld1,fld2,...,i1,i2,j1,j2,k1,k2]  computes the ',
     >  '                volume integral of the product of fields ',
     >  '                ''fld1'',''fld2'',... over the domain ',
     >  '                defined by i1*dx to i2*dx, j1*dy to j2*dy,',
     >  '                k1*dz to k2*dz ',
     >  'RUN_INT[fld,var0,int0]   computes the running time integral',
     >  '                of field ''fld''.  Here ''var0'' is the ',
     >  '                value of the field at t=0 and ''int0'' is ',
     >  '                the value of the integral at t=0.',
     >  'AREA_INT_Z[fld1,fld2,...,i1,i2,j1,j2,k]  computes the area',
     >  '                integral at constant z of the product of ',
     >  '                the fields ''fld1'',''fld2'',... over the',
     >  '                domain defined by i1*dx to i2*dx, j1*dy to',
     >  '                j2*dy, at k*dz',
     >  'AREA_INT_X[fld1,fld2,...,j1,j2,k1,k2,i]  computes the area',
     >  '                integral at constant x of the product of ',
     >  '                the fields ''fld1'',''fld2'',... over the',
     >  '                domain defined by j1*dy to j2*dy, k1*dz to',
     >  '                k2*dz, at i*dx',
     >  'AREA_INT_Y[fld1,fld2,...,i1,i2,k1,k2,i]  computes the area',
     >  '                integral at constant y of the product of ',
     >  '                the fields ''fld1'',''fld2'',... over the',
     >  '                domain defined by i1*dx to i2*dx, k1*dz to',
     >  '                k2*dz, at j*dy',
     >  'BAR[fld1,fld2,...,i1,i2,j1,j2]   computes the area integral',
     >  '                at constant z of the product of the fields ',
     >  '                ''fld1'',''fld2'',... over the domain ',
     >  '                defined by i1*dx to i2*dx, j1*dy to j2*dy.',
     >  '                BAR returns this integral as a function of ',
     >  '                height and time.'/


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

      
      
