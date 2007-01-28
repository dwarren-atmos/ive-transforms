

      program put_sclr

      include '/usr/local/include/netcdf.inc'

      integer maxsclr,numsclr

      parameter (maxsclr = 10)

      integer idcdf,id_one,idvar(maxsclr),lgth,lendef(maxsclr),k
      character*(80) outfile,varnm(maxsclr),units(maxsclr)
      character*(80) def(maxsclr)
      real val(maxsclr)

 23   FORMAT(A80)
 24   FORMAT('Variable ',I3)

      write(6,*) 'Output file:  '
      read(6,23) outfile

      write(6,*) 'Number of scalars:'
      read(6,*) numsclr
    
      if (numsclr.gt.10) then
        write(6,*) 'Error:  too many scalars.'
        stop
      end if 
      
      do k=1,numsclr

        write(6,*)
        write(6,24) k
        write(6,*)
 
        write(6,*) 'Variable name:  '
        read(6,23) varnm(k)

        write(6,*) 'Value: '
        read(6,*) val(k)

        write(6,*) 'Units:  '
        read(6,23) units(k)

        write(6,*) 'Def:  '
        read(6,23) def(k)

        write(6,*) 'Def length: '
        read(6,*) lendef(k)

      enddo 

      lgth = index(outfile,' ') - 1

      idcdf = ncopn(outfile(1:lgth),NCWRITE,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  cannot open file ',outfile(1:lgth)
        stop
      end if

      id_one = ncdid(idcdf,'one',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  failed to find id for dimension one.'
        call ncclos(idcdf,ierr)
        stop
      end if

      call ncredf(idcdf,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  failed to put file into definition mode.'
        call ncclos(idcdf,ierr)
        stop
      end if

      do k=1,numsclr

        lgth = index(varnm(k),' ') - 1

        idvar(k) = ncvdef(idcdf,varnm(k)(1:lgth),NCFLOAT,1,id_one,
     >                                 ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  failed to define variable.'
          call ncclos(idcdf,ierr)
          stop
        end if

        lgth = index(units(k),' ') - 1

        call ncaptc(idcdf,idvar(k),'units',NCCHAR,lgth,
     >                    units(k)(1:lgth),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  failed to put attribute units.'
          call ncclos(idcdf,ierr)
          stop
        end if

        call ncaptc(idcdf,idvar(k),'def',NCCHAR,lendef(k),
     >                   def(k)(1:lendef(k)),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  failed to put attribute def.'
          call ncclos(idcdf,ierr)
          stop
        end if

        call ncapt(idcdf,idvar(k),'no_button',NCLONG,1,1,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  failed to put attribute no_button.'
          call ncclos(idcdf,ierr)
          stop
        end if

      enddo

      call ncendf(idcdf,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  failed to exit definition mode.'
        call ncclos(idcdf,ierr)
        stop
      end if

      do k=1,numsclr

        call ncvpt1(idcdf,idvar(k),1,val(k),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  failed to write data value.'
          call ncclos(idcdf,ierr)
          stop
        end if
     
      enddo

      call ncclos(idcdf,ierr)

      end


      




       
