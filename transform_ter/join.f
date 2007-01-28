
      program join

      include '/usr/local/include/netcdf.inc'

      integer maxsrc
      parameter (maxsrc=5)

      integer nfl,idsrc(maxsrc),idnew,nt
      character*(80) srcfl(maxsrc),newfl

      integer ispc
      parameter (ispc=2000000)

      real a(ispc)

      integer itime,isrcfl,insrc,ivar

      integer maxvar,maxdim,maxgatt,maxvatt
      parameter (maxvar=40,maxdim=10,maxgatt=25,maxvatt=10)

      common/newfl_info/ndims,dim_size,dim_names,nvars,v_names,
     >          ngatts,gatt_names,igatt,rgatt,cgatt,nvdims,vdims,
     >          nvatts,vatt_names,ivatt,rvatt,cvatt,rec_dim,
     >          rec_var
      integer ndims,ngatts,nvars,nvatts(maxvar)
      integer dim_size(maxdim),nvdims(maxvar),vdims(4,maxvar)
      integer igatt(maxgatt),ivatt(maxvatt,maxvar)
      integer rec_dim,rec_var
      real rgatt(maxgatt),rvatt(maxvatt,maxvar)
      character*(80) dim_names(maxdim),v_names(maxvar)
      character*(80) gatt_names(maxgatt,2)
      character*(80) vatt_names(maxvatt,maxvar,2)
      character*(80) cgatt(maxgatt),cvatt(maxvatt,maxvar)

      integer ierr,k,n,q,p,ifnd
      integer itmp,size,idnt,maxt,numt
      character*(80) ctmp

 101  FORMAT(A80)

      open(unit=12,file='join.in',status='old')

      read(12,101) newfl
      read(12,*) nfl
      if (nfl.gt.maxsrc) then
        write(6,*) 'Error:  too many files.  Increase maxsrc.' 
        stop
      end if
      do k=1,nfl
        read(12,101) srcfl(k)
        lgth = index(srcfl(k),' ')-1
        srcfl(k) = srcfl(k)(1:lgth)//'.cdf'
        idsrc(k) =  ncopn(srcfl(k)(1:lgth+4),NCNOWRIT,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem opening source file.'
          stop
        end if
      enddo

      close(unit=12)

      call get_stuff(idsrc,nfl)

      call set_newfl(newfl,idnew,srcfl,nfl)

      maxt = 0
      do k=1,nfl
        call ncdinq(idsrc(k),rec_dim,ctmp,numt,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem inquiring about dimension time.'
          stop
        end if
        maxt = maxt + numt
      enddo

      itime = 1
      isrcfl = itime + maxt
      insrc = isrcfl + maxt
      ivar = insrc + maxt

      call set_times(a(itime),a(isrcfl),a(insrc),idsrc,nt,maxt,nfl)

      size = 1
      itmp = 0
      do k=1,ndims
        lgth = index(dim_names(k),' ')-1
        if ((dim_names(k)(1:lgth).eq.'nxp1').or.
     >          (dim_names(k)(1:lgth).eq.'nyp1').or.
     >               (dim_names(k)(1:lgth).eq.'nzp1')) then
          itmp = itmp + 1
          size = size*dim_size(k)
        end if
      enddo

      if (itmp.lt.3) then
        write(6,*)
        write(6,*) 'Error:  did not find three spatial dimensions.'
        write(6,*) '        Proceeding without memory check.'
        write(6,*)
        size = ispc-ivar+1
      else if (ispc-ivar+1.lt.size) then
        write(6,*) 'Error:  not enough memory for variables.  ',
     >                 'Increase parameter ispc.'
        stop
      end if

      do k=1,nvars

        q = 1
        ifnd = 0
 87     if ((ifnd.eq.0).and.(q.le.nvdims(k))) then
          if (vdims(q,k).eq.rec_dim) then
            ifnd = 1
          else 
            q = q+1
          end if
          goto 87
        end if

        if (ifnd.eq.0) then
          call joinvar(a(ivar),a(isrcfl),a(insrc),idnew,idsrc,size,
     >                                    nfl,nt,k,0)
        end if

      enddo

      do n=1,nt
      
        call ncvpt1(idnew,rec_var,n,a(itime+n-1),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem writing time.'
          stop
        end if

        do k=1,nvars

          if (k.ne.rec_var) then

            q = 1
            ifnd = 0
 897        if ((ifnd.eq.0).and.(q.le.nvdims(k))) then
              if (vdims(q,k).eq.rec_dim) then
                ifnd = 1
              else 
                q = q+1
              end if
              goto 897
            end if

            if (ifnd.eq.1) then
              call joinvar(a(ivar),a(isrcfl),a(insrc),idnew,idsrc,
     >                            size,nfl,nt,k,n)
            end if

          end if

        enddo

      enddo

      call ncclos(idnew,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem closing new file.'
        stop
      end if 

      do k=1,nfl
        call ncclos(idsrc(k),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem closing source file.'
          stop
        end if
      enddo

      end
      
*--------------------------------------------------------------------

      subroutine set_times(times,srcfl,nsrc,idsrc,nt,maxt,nfl)

      integer nfl,maxt,idsrc(nfl),srcfl(maxt),nsrc(maxt)
      real times(maxt)

      integer ierr,k,q,p,itmp,idtime,idnt,numt
      real rtmp
      character*(80) ctmp

      nt = 0

      do k=1,nfl

        idtime = ncvid(idsrc(k),'time',ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem getting id for variable time.'
          stop
        end if
        idnt = ncdid(idsrc(k),'time',ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem getting id for dimension time.'
          stop
        end if
        call ncdinq(idsrc(k),idnt,ctmp,numt,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem inquiring about dimension time.'
          stop
        end if

        do q=1,numt

          call ncvgt1(idsrc(k),idtime,q,rtmp,ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem reading variable time.'
            stop
          end if

          p = 1
          ifnd = 0
 343      if ((ifnd.eq.0).and.(p.le.nt)) then
            if (rtmp.eq.times(p)) then
              ifnd = 1
            else
              p = p+1
            end if
            goto 343
          end if

          if (ifnd.eq.0) then
            nt = nt+1
            times(nt) = rtmp
            srcfl(nt) = k
            nsrc(nt) = q
          end if

        enddo

      enddo
         
      do k=1,nt-1
      do q=nt-1,k,-1

        if (times(q).gt.times(q+1)) then
          rtmp = times(q)
          times(q) = times(q+1)
          times(q+1) = rtmp
          itmp = srcfl(q)
          srcfl(q) = srcfl(q+1)
          srcfl(q+1) = itmp
          itmp = nsrc(q)
          nsrc(q) = nsrc(q+1)
          nsrc(q+1) = itmp
        end if

      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine get_stuff(idsrc,nfl)

      include '/usr/local/include/netcdf.inc'

      integer nfl,idsrc(nfl)

      integer maxvar,maxdim,maxgatt,maxvatt
      parameter (maxvar=40,maxdim=10,maxgatt=25,maxvatt=10)

      common/newfl_info/ndims,dim_size,dim_names,nvars,v_names,
     >          ngatts,gatt_names,igatt,rgatt,cgatt,nvdims,vdims,
     >          nvatts,vatt_names,ivatt,rvatt,cvatt,rec_dim,
     >          rec_var
      integer ndims,ngatts,nvars,nvatts(maxvar)
      integer dim_size(maxdim),nvdims(maxvar),vdims(4,maxvar)
      integer igatt(maxgatt),ivatt(maxvatt,maxvar)
      integer rec_dim,rec_var
      real rgatt(maxgatt),rvatt(maxvatt,maxvar)
      character*(80) dim_names(maxdim),v_names(maxvar)
      character*(80) gatt_names(maxgatt,2)
      character*(80) vatt_names(maxvatt,maxvar,2)
      character*(80) cgatt(maxgatt),cvatt(maxvatt,maxvar)

      integer k,p,q,vtyp,ierr,itmp(6),ifnd,lgth,lgth2
      integer atyp,alen
      real*8 atmp
      character*(80) ctmp(2)

      call ncinq(idsrc(1),ndims,nvars,ngatts,rec_dim,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem with file inquiry.'
        stop
      end if

      do k=2,nfl
        call ncinq(idsrc(k),itmp(1),itmp(2),itmp(3),itmp(4),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem with file inquiry.'
          stop
        else if (ndims.ne.itmp(1)) then
          write(6,*) 'Error:  files have different numbers of ',
     >                              'dimensions.'
          stop
        else if (nvars.ne.itmp(2)) then 
          write(6,*) 'Error:  files have different numbers of ',
     >                              'variables.'
          stop
        else if (ngatts.ne.itmp(3)) then
          write(6,*) 'Error:  files have different numbers of ',
     >                              'global attributes.'
          stop
        end if
      enddo

      if (nvars.gt.maxvar) then
        write(6,*) 'Error:  too many variables.  Increase maxvar.'
        stop
      else if (ndims.gt.maxdim) then
        write(6,*) 'Error:  too many dimensions.  Increase maxdim.'
        stop
      else if (ngatts.gt.maxgatt) then
        write(6,*) 'Error:  too many global attributes.  Increase ',
     >                       'maxgatt.'
        stop
      end if

      do p=1,ndims

        call ncdinq(idsrc(1),p,dim_names(p),dim_size(p),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem with dimension inquiry.'
          stop
        end if
        lgth = index(dim_names(p),' ')-1

        do k=2,nfl
          call ncdinq(idsrc(k),p,ctmp(1),itmp(1),ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem with dimension inquiry.'
            stop
          end if
          lgth2 = index(ctmp(1),' ')-1
          if (ctmp(1)(1:lgth2).ne.dim_names(p)(1:lgth)) then
            write(6,*) 'Error:  files have different dimensions.'
            stop
          else if ((itmp(1).ne.dim_size(p)).and.(p.ne.rec_dim)) then
            write(6,*) 'Error:  files have different dimensions.'
            stop
          end if
        enddo

      enddo 

      do p=1,nvars

        call ncvinq(idsrc(1),p,v_names(p),vtyp,nvdims(p),vdims(1,p),
     >                            nvatts(p),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem with variable inquiry.'
          stop
        else if (nvatts(p).gt.maxvatt) then
          write(6,*) 'Error:  too many variable attributes. ',
     >                 'Increase maxvatt.'
          stop
        end if 
        lgth = index(v_names(p),' ')-1

        do k=2,nfl

          ifnd = 0
          q = 1
 1234     if ((ifnd.eq.0).and.(q.le.nvars)) then
            call ncvinq(idsrc(k),q,ctmp(1),vtyp,itmp(5),itmp(1),
     >                             itmp(6),ierr)
            if (ierr.ne.0) then
              write(6,*) 'Error:  problem with variable inquiry.'
              stop
            end if
            lgth2 = index(ctmp(1),' ')-1
            if (ctmp(1)(1:lgth2).eq.v_names(p)(1:lgth)) then
              ifnd = 1
            else
              q = q+1
            end if
            goto 1234
          end if

          if ((ifnd.eq.0).or.(itmp(5).ne.nvdims(p))) then
            write(6,*) 'Error:  files have different variables.'
            stop
          else
            q = 1
 2233       if (q.le.nvdims(p)) then
              if (itmp(q).ne.vdims(q,p)) then
                write(6,*) 'Error:  files have different variables.'
                stop
              else
                q = q+1
              end if
              goto 2233
            end if
          end if

        enddo   

        do k=1,nvatts(p)

          call ncanam(idsrc(1),p,k,vatt_names(k,p,1),ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem getting attribute name.'
            stop
          end if

          call ncainq(idsrc(1),p,vatt_names(k,p,1),atyp,alen,ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem with attribute inquiry.'
            stop
          end if

          if (atyp.eq.NCCHAR) then
            vatt_names(k,p,2) = 'c'
            call ncagtc(idsrc(1),p,vatt_names(k,p,1),cvatt(k,p),80,
     >                                      ierr)
            if (ierr.ne.0) then
              write(6,*) 'Error:  problem reading string attribute.'
              stop
            end if
          else if (atyp.eq.NCLONG) then
            vatt_names(k,p,2) = 'i'
            call ncagt(idsrc(1),p,vatt_names(k,p,1),ivatt(k,p),ierr)
            if (ierr.ne.0) then
              write(6,*) 'Error:  problem reading integer attribute.'
              stop
            end if
          else if (atyp.eq.NCFLOAT) then
            vatt_names(k,p,2) = 'r'
            call ncagt(idsrc(1),p,vatt_names(k,p,1),rvatt(k,p),ierr)
            if (ierr.ne.0) then
              write(6,*) 'Error:  problem reading real attribute.'
              stop
            end if
          else if (atyp.eq.NCDOUBLE) then
            vatt_names(k,p,2) = 'r'
            call ncagt(idsrc(1),p,vatt_names(k,p,1),atmp,ierr)
            if (ierr.ne.0) then
              write(6,*) 'Error:  problem reading real attribute.'
              stop
            end if
            rvatt(k,p) = real(atmp)
          else
            write(6,*) 'Error:  attribute type not recognized.' 
            stop 
          end if

        enddo

      enddo

      do p=1,ngatts

        call ncanam(idsrc(1),NCGLOBAL,p,gatt_names(p,1),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem getting global attribute name.'
          stop
        end if 

        call ncainq(idsrc(1),NCGLOBAL,gatt_names(p,1),atyp,
     >                                 alen,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem inquiring about global ',
     >                                'attribute.'
          stop
        end if

        if (atyp.eq.NCCHAR) then
          gatt_names(p,2) = 'c'
          call ncagtc(idsrc(1),NCGLOBAL,gatt_names(p,1),
     >                     cgatt(p),80,ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem reading string global ',
     >                              'attribute.'
            stop
          end if
        else if (atyp.eq.NCLONG) then
          gatt_names(p,2) = 'i'
          call ncagt(idsrc(1),NCGLOBAL,gatt_names(p,1),
     >                    igatt(p),ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem reading integer global ',
     >                              'attribute.'
            stop
          end if
        else if (atyp.eq.NCFLOAT) then
          gatt_names(p,2) = 'r'
          call ncagt(idsrc(1),NCGLOBAL,gatt_names(p,1),
     >                    rgatt(p),ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem reading real global ',
     >                              'attribute.'
            stop
          end if
        else if (atyp.eq.NCDOUBLE) then
          gatt_names(p,2) = 'r'
          call ncagt(idsrc(1),NCGLOBAL,gatt_names(p,1),
     >                    atmp,ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem reading real global ',
     >                              'attribute.'
            stop
          end if
          rgatt(p) = real(atmp)
        else
          write(6,*) 'Error:  global attribute type not recognized.'
          stop
        end if

      enddo

      if (rec_dim.ne.-1) then

        lgth = index(dim_names(rec_dim),' ')-1
        k = 1
        ifnd = 0
 2020   if ((ifnd.eq.0).and.(k.le.nvars)) then
          lgth2 = index(v_names(k),' ')-1
          if (v_names(k)(1:lgth2).eq.dim_names(rec_dim)(1:lgth)) then 
            ifnd = 1
            rec_var = k
          else
            k = k+1
          end if
          goto 2020
        end if
    
        if (ifnd.eq.0) then
          write(6,*)
          write(6,*) 'Warning:  no coordinate variable found to ',
     >                  'match record dimension.'
          write(6,*)
          rec_var = -1
        end if

      else
   
        rec_var = -1

      end if
  
      return
      end

*----------------------------------------------------------------------

      subroutine set_newfl(newfl,idnew,srcfl,nfl)

      include '/usr/local/include/netcdf.inc'

      integer idnew,nfl
      character*(*) newfl,srcfl(nfl)

      integer maxvar,maxdim,maxgatt,maxvatt
      parameter (maxvar=40,maxdim=10,maxgatt=25,maxvatt=10)

      common/newfl_info/ndims,dim_size,dim_names,nvars,v_names,
     >          ngatts,gatt_names,igatt,rgatt,cgatt,nvdims,vdims,
     >          nvatts,vatt_names,ivatt,rvatt,cvatt,rec_dim,
     >          rec_var
      integer ndims,ngatts,nvars,nvatts(maxvar)
      integer dim_size(maxdim),nvdims(maxvar),vdims(4,maxvar)
      integer igatt(maxgatt),ivatt(maxvatt,maxvar)
      integer rec_dim,rec_var
      real rgatt(maxgatt),rvatt(maxvatt,maxvar)
      character*(80) dim_names(maxdim),v_names(maxvar)
      character*(80) gatt_names(maxgatt,2)
      character*(80) vatt_names(maxvatt,maxvar,2)
      character*(80) cgatt(maxgatt),cvatt(maxvatt,maxvar)

      integer k,q,p,slen,iunit,lgth,lgth2
      character*(80) outstr,cdl_name,cdf_name,command,command2
      character*(1) num

      cdl_name(1:len(cdl_name)) = ' '
      lgth = index(newfl,' ')-1
      cdl_name = newfl(1:lgth)//'.cdl'

      write(6,*)
      write(6,*) 'JOIN_OUT:  creating new file'
      write(6,*) 'JOIN_OUT:  cdl file is ',cdl_name(1:lgth+4)

      iunit = 777
      open(unit=iunit,file=cdl_name(1:lgth+4),status='new')

      write(iunit,*) 'netcdf JoinOut{'
      write(iunit,*) 

      write(iunit,*) 'dimensions:'
      write(iunit,*)
      do k=1,ndims
        lgth2 = index(dim_names(k),' ')-1
        if (k.eq.rec_dim) then
          write(iunit,*) dim_names(k)(1:lgth2),'=UNLIMITED;'
        else
          write(iunit,*) dim_names(k)(1:lgth2),'=',dim_size(k),';'
        end if
      enddo
      write(iunit,*)

      write(iunit,*) 'variables:'
      write(iunit,*)

      do k=1,nvars

        lgth = index(v_names(k),' ')-1
        outstr(1:len(outstr)) = ' '
        outstr = 'float '//v_names(k)(1:lgth)
        slen = 6 + lgth
        outstr = outstr(1:slen)//'('
        slen = slen + 1
        do q=nvdims(k),1,-1
          lgth2 = index(dim_names(vdims(q,k)),' ')-1
          outstr = outstr(1:slen)//dim_names(vdims(q,k))(1:lgth2)
          slen = slen + lgth2
          outstr = outstr(1:slen)//','
          slen = slen + 1
        enddo
        outstr = outstr(1:slen-1)//');'

        write(iunit,*) outstr(1:slen+1)

        do p=1,nvatts(k)

          outstr(1:len(outstr)) = ' '
          outstr = v_names(k)(1:lgth)//':'
          lgth2 = index(vatt_names(p,k,1),' ')-1
          outstr = outstr(1:lgth+1)//vatt_names(p,k,1)(1:lgth2)
          slen = lgth+1+lgth2
          outstr = outstr(1:slen)//'='
          if (vatt_names(p,k,2)(1:1).eq.'c') then
            lgth2 = index(cvatt(p,k),'  ')-1
            if (lgth2.eq.0) lgth2 = 1
            write(iunit,*) outstr(1:slen+1),'"',cvatt(p,k)(1:lgth2),
     >                             '";'
          else if (vatt_names(p,k,2)(1:1).eq.'i') then 
            write(iunit,*) outstr(1:slen+1),ivatt(p,k),';'
          else if (vatt_names(p,k,2)(1:1).eq.'r') then
            write(iunit,*) outstr(1:slen+1),rvatt(p,k),';'
          end if

        enddo 

        write(iunit,*)

      enddo

      write(iunit,*) '//global attributes:'
      write(iunit,*)

      do p=1,ngatts

        lgth = index(gatt_names(p,1),' ')-1
        if (gatt_names(p,1).ne.'runname') then
          outstr(1:len(outstr)) = ' '
          outstr = ':'//gatt_names(p,1)(1:lgth)
          outstr = outstr(1:lgth+1)//'='
          if (gatt_names(p,2)(1:1).eq.'c') then
            lgth2 = index(cgatt(p),'  ')-1
            write(iunit,*) outstr(1:lgth+2),'"',cgatt(p)(1:lgth2),'";'
          else if (gatt_names(p,2)(1:1).eq.'i') then
            write(iunit,*) outstr(1:lgth+2),igatt(p),';'
          else if (gatt_names(p,2)(1:1).eq.'r') then
            write(iunit,*) outstr(1:lgth+2),rgatt(p),';'
          end if
        end if

      enddo

 123  FORMAT(I1)

      lgth = index(newfl,' ')-1
      write(iunit,*) ':runname="union of runs ',newfl(1:lgth),'";'
      do p=1,nfl
        lgth2 = index(srcfl(p),' ')-5
        write(num,123) p
        write(iunit,*) ':source',num,'="',srcfl(p)(1:lgth2),'";'
      enddo
      write(iunit,*) 
      write(iunit,*) '}'

      close(iunit)

      cdf_name = newfl(1:lgth)//'.cdf'
      command(1:len(command)) = ' '
      write(6,*) 'JOIN_OUT:  netcdf file is ',cdf_name(1:lgth+4)
      write(6,*)
      write(command,*) 'ncgen -o',cdf_name(1:lgth+4),' ',
     >                     cdl_name(1:lgth+4)
      call system(command)
      write(command2,*) 'rm -f ',cdl_name(1:lgth+4)
      call system(command2)     

      idnew = ncopn(cdf_name(1:lgth+4),NCWRITE,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem opening new data file.'
        stop
      end if

      return 
      end 

*----------------------------------------------------------------------

      subroutine joinvar(var,srcfl,nsrc,idnew,idsrc,size,nfl,nt,k,n)

      integer size,nfl,nt,srcfl(nt),nsrc(nt),idnew,idsrc(nfl),k,n
      real var(size)

      integer maxvar,maxdim,maxgatt,maxvatt
      parameter (maxvar=40,maxdim=10,maxgatt=25,maxvatt=10)

      common/newfl_info/ndims,dim_size,dim_names,nvars,v_names,
     >          ngatts,gatt_names,igatt,rgatt,cgatt,nvdims,vdims,
     >          nvatts,vatt_names,ivatt,rvatt,cvatt,rec_dim,
     >          rec_var
      integer ndims,ngatts,nvars,nvatts(maxvar)
      integer dim_size(maxdim),nvdims(maxvar),vdims(4,maxvar)
      integer igatt(maxgatt),ivatt(maxvatt,maxvar)
      integer rec_dim,rec_var
      real rgatt(maxgatt),rvatt(maxvatt,maxvar)
      character*(80) dim_names(maxdim),v_names(maxvar)
      character*(80) gatt_names(maxgatt,2)
      character*(80) vatt_names(maxvatt,maxvar,2)
      character*(80) cgatt(maxgatt),cvatt(maxvatt,maxvar)

      integer idvar,ist(4),iln(4),ierr,q,p,ifnd,itmp
      character*(80) ctmp

      if (n.eq.0) then
        q = 1
      else
        q = srcfl(n)
      end if

      idvar = ncvid(idsrc(q),v_names(k),ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for variable.'
        stop
      end if

      do p=1,nvdims(k)
        if (vdims(p,k).eq.rec_dim) then
          ist(p) = nsrc(n)
          iln(p) = 1
        else
          ist(p) = 1
          iln(p) = dim_size(vdims(p,k))
        end if
      enddo

      call ncvgt(idsrc(q),idvar,ist,iln,var,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem reading variable.'
        stop
      end if
 
      ifnd = 0
      p = 1
 759  if ((ifnd.eq.0).and.(p.le.nvdims(k))) then
        if (vdims(p,k).eq.rec_dim) then
          ifnd = 1
          ist(p) = n
          iln(p) = 1
        else
          p = p+1
        end if
        goto 759
      end if

      call ncvpt(idnew,k,ist,iln,var,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem writing variable.'
        stop
      end if

      return
      end


