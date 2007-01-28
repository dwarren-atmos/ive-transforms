
      program make_an

      include '/usr/local/include/netcdf.inc'

      integer mxsrc,nsrc
      parameter (mxsrc=10)

      integer icnst(0:mxsrc),rcnst(0:mxsrc),isrc(2,mxsrc)
      integer jsrc(2,mxsrc),ian(2,mxsrc),jan(2,mxsrc),idcdf(0:mxsrc)
      integer itime,ntime,isclr(mxsrc),nsclr(mxsrc),sgz,wgz
      integer srctm,srcntm,fptr,sspc,aspc,mspc,ksrc,kvar,i_val
      logical sclr_err

      integer mxsclr
      parameter (mxsclr=15)

      character*(80) srcfl(mxsrc),anfl,sclrnm(mxsclr,mxsrc)

      integer ispc
      parameter (ispc=12000000)

      real a(ispc)

      integer nvar
      parameter (nvar=5)

      integer van(nvar),vsrc

      character*(80) var(nvar),units(nvar),tvar(3)
      data var / 'u','v','w','p','thet' /
      data units / 'm/s','m/s','m/s','N/m/m','K' /

 886  FORMAT(A80)

      call ncpopt(NCVERBOS)

      itime = 1

      iunit = 32

      open (unit=iunit,file='make_an.in',status='old')

      read(iunit,*) nsrc
      do k=1,nsrc
        read(iunit,886) srcfl(k)
      enddo
      read(iunit,*) ntime
      do k=1,ntime
        read(iunit,*) a(k)
      enddo

      icnst(0) = itime+ntime
      rcnst(0) = icnst(0)+3

      call an_param(a(icnst(0)),a(rcnst(0)),iunit)

      read(iunit,886) anfl

      close (unit=iunit)

      do k=1,nsrc
        lgth = index(srcfl(k),' ')-1 
        srcfl(k) = srcfl(k)(1:lgth)//'.cdf'
        idcdf(k) = ncopn(srcfl(k)(1:lgth+4),NCNOWRIT,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem opening file ',
     >                     srcfl(k)(1:lgth+4)
          stop
        end if
        icnst(k) = rcnst(k-1)+9
        rcnst(k) = icnst(k)+3
        call get_atts(idcdf(k),a(icnst(k)),a(rcnst(k)))
      enddo

      fptr = rcnst(nsrc)+9

      call srt_fls(idcdf(1),srcfl(1),a(1),ispc,icnst(0),rcnst(0),nsrc)

      a(icnst(0)+2) = a(icnst(1)+2)
      a(rcnst(0)+2) = a(rcnst(1)+2)
      a(rcnst(0)+5) = a(rcnst(1)+5)
      a(rcnst(0)+8) = a(rcnst(1)+8)

      mspc = 0
      do k=0,nsrc
        sspc = i_val(a(icnst(k)))*i_val(a(icnst(k)+1))
     >                      *i_val(a(icnst(k)+2)) 
        if (k.eq.0) then
          aspc = sspc
        else
          if (sspc.gt.mspc) mspc = sspc
          isrc(1,k) = fptr
          isrc(2,k) = isrc(1,k) + i_val(a(icnst(0)))
          jsrc(1,k) = isrc(2,k) + i_val(a(icnst(0)))
          jsrc(2,k) = jsrc(1,k) + i_val(a(icnst(0)+1))
          ian(1,k) = jsrc(2,k) + i_val(a(icnst(0)+1))
          ian(2,k) = ian(1,k) + i_val(a(icnst(0)))
          jan(1,k) = ian(2,k) + i_val(a(icnst(0)))
          jan(2,k) = jan(1,k) + i_val(a(icnst(0)+1))
          fptr = jan(2,k) + i_val(a(icnst(0)+1))
        end if
      enddo

      vsrc = fptr
      van(1) = vsrc+mspc
      do k=2,nvar
        van(k) = van(k-1)+aspc
      enddo
      fptr = van(nvar) + aspc
       
 449  FORMAT('Error:  grids ',I2,' and ',I2,' have different 
     >                            scalars.') 

      do k=1,nsrc
        isclr(k) = fptr
        call get_sclrs(idcdf(k),a(isclr(k)),sclrnm(1,k),nsclr(k),
     >                                   mxsclr)
        if (k.gt.1) then
          if (sclr_err(a(isclr(1)),a(isclr(k)),
     >        sclrnm(1,1),sclrnm(1,k),nsclr(1),nsclr(k))) then
            write(6,449) k,1
            stop
          end if
        end if
        fptr = fptr+nsclr(k)
      enddo

      sgz = fptr
      wgz = sgz + i_val(a(icnst(1)+2))
      srctm = wgz + i_val(a(icnst(1)+2))

      call get1da(idcdf(1),a(sgz),a(wgz),a(srctm),a(icnst(1)+2),
     >                              srcntm)

      fptr = srctm + srcntm

      call set_anfl(idcdf(0),a(icnst(0)),a(icnst(0)+1),a(icnst(0)+2),
     >                a(rcnst(0)),a(sgz),a(wgz),a(isclr(1)),
     >                sclrnm(1,1),var,units,anfl,srcfl,nsrc,
     >                nsclr(1),nvar)               

      do k=1,nsrc
        call intrp_info(a(rcnst(0)),a(rcnst(k)),a(isrc(1,k)),
     >                    a(isrc(2,k)),a(jsrc(1,k)),a(jsrc(2,k)),
     >                    a(ian(1,k)),a(ian(2,k)),a(jan(1,k)),
     >                    a(jan(2,k)),a(icnst(0)),a(icnst(0)+1),
     >                    a(icnst(k)),a(icnst(k)+1))
      enddo 

      tvar(1) = 'zbot_p'
      tvar(2) = 'zbot_u'
      tvar(3) = 'zbot_v'

      do kvar=1,3
        do k=1,nsrc
          call read_var(idcdf(k),tvar(kvar),a(icnst(k)),a(icnst(k)+1),
     >                           a(icnst(k)+2),a(vsrc),1)
          call interp_var(a(vsrc),a(van(1)),a(isrc(1,k)),
     >                a(isrc(2,k)),a(jsrc(1,k)),a(jsrc(2,k)),
     >                a(ian(1,k)),a(ian(2,k)),a(jan(1,k)),
     >                a(jan(2,k)),a(icnst(0)),a(icnst(0)+1),1,
     >                a(icnst(k)),a(icnst(k)+1),1,tvar(kvar))
        enddo
        call write_var(idcdf(0),tvar(kvar),a(icnst(0)),a(icnst(0)+1),
     >                a(icnst(0)+2),a(van(1)),1)
      enddo

      do k=1,nsrc
      do kvar=1,nvar
        lgth = index(var(kvar),' ')-1
        if (var(kvar)(lgth:lgth).eq.'0') then
          call read_var(idcdf(k),var(kvar),a(icnst(k)),a(icnst(k)+1),
     >                         a(icnst(k)+2),a(vsrc),1)
          call interp_var(a(vsrc),a(van(kvar)),a(isrc(1,k)),
     >                    a(isrc(2,k)),a(jsrc(1,k)),a(jsrc(2,k)),
     >                    a(ian(1,k)),a(ian(2,k)),a(jan(1,k)),
     >                    a(jan(2,k)),a(icnst(0)),a(icnst(0)+1),
     >                    a(icnst(0)+2),a(icnst(k)),a(icnst(k)+1),
     >                    a(icnst(k)+2),var(kvar))
        end if
      enddo
      enddo

      do kvar=1,nvar
        lgth = index(var(kvar),' ')-1
        if (var(kvar)(lgth:lgth).eq.'0') then
          call write_var(idcdf(0),var(kvar),a(icnst(0)),
     >             a(icnst(0)+1),a(icnst(0)+2),a(van(kvar)),1)
        end if
      enddo

      idtime = ncvid(idcdf(0),'time',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for time variable.'
        stop
      end if

      do n=1,ntime

        call ncvpt1(idcdf(0),idtime,n,a(itime+n-1),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem writing time.'
          stop
        end if

        ktime = int_time(a(itime+n-1),a(srctm),srcntm)

        do k=1,nsrc
        do kvar=1,nvar
          lgth = index(var(kvar),' ')-1
          if (var(kvar)(lgth:lgth).ne.'0') then
            call read_var(idcdf(k),var(kvar),a(icnst(k)),a(icnst(k)+1),
     >                           a(icnst(k)+2),a(vsrc),ktime)
            call interp_var(a(vsrc),a(van(kvar)),a(isrc(1,k)),
     >                      a(isrc(2,k)),a(jsrc(1,k)),a(jsrc(2,k)),
     >                      a(ian(1,k)),a(ian(2,k)),a(jan(1,k)),
     >                      a(jan(2,k)),a(icnst(0)),a(icnst(0)+1),
     >                      a(icnst(0)+2),a(icnst(k)),a(icnst(k)+1),
     >                      a(icnst(k)+2),var(kvar))
          end if
        enddo
        enddo

        do kvar=1,nvar
          lgth = index(var(kvar),' ')-1
          if (var(kvar)(lgth:lgth).ne.'0') then
            call write_var(idcdf(0),var(kvar),a(icnst(0)),
     >               a(icnst(0)+1),a(icnst(0)+2),a(van(kvar)),n)
          end if
        enddo

      enddo

      do k=0,nsrc
        call ncclos(idcdf(k),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem closing file.'
          stop
        end if
      enddo

      end

*------------------------------------------------------------------------

      subroutine an_param(icnst,rcnst,iunit)

      integer icnst(3),iunit
      real rcnst(9)

      read(iunit,*) icnst(1)
      read(iunit,*) icnst(2)
      read(iunit,*) rcnst(1)
      read(iunit,*) rcnst(2)
      read(iunit,*) rcnst(4)
      read(iunit,*) rcnst(5)

      rcnst(7) = rcnst(4) + (icnst(1)-1)*rcnst(1)
      rcnst(8) = rcnst(5) + (icnst(2)-1)*rcnst(2)

      return
      end

*-------------------------------------------------------------------------

      subroutine get_atts(idcdf,icnst,rcnst)

      integer idcdf,icnst(3)
      real rcnst(9)

      real*8 aval

      call ncagt(idcdf,NCGLOBAL,'x_delta',aval,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting x_delta.'
        stop
      end if
      rcnst(1) = real(aval)
      call ncagt(idcdf,NCGLOBAL,'y_delta',aval,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting y_delta.'
        stop
      end if
      rcnst(2) = real(aval)
      call ncagt(idcdf,NCGLOBAL,'z_delta',aval,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting z_delta.'
        stop
      end if
      rcnst(3) = real(aval)
      call ncagt(idcdf,NCGLOBAL,'x_min',aval,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting x_min.'
        stop
      end if
      rcnst(4) = real(aval)
      call ncagt(idcdf,NCGLOBAL,'y_min',aval,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting y_min.'
        stop
      end if
      rcnst(5) = real(aval)
      call ncagt(idcdf,NCGLOBAL,'z_min',aval,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting z_min.'
        stop
      end if
      rcnst(6) = real(aval)
      call ncagt(idcdf,NCGLOBAL,'x_max',aval,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting x_max.'
        stop
      end if
      rcnst(7) = real(aval)
      call ncagt(idcdf,NCGLOBAL,'y_max',aval,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting y_max.'
        stop
      end if
      rcnst(8) = real(aval)
      call ncagt(idcdf,NCGLOBAL,'z_max',aval,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting z_max.'
        stop
      end if
      rcnst(9) = real(aval)

      icnst(1) = nint((rcnst(7)-rcnst(4))/rcnst(1)) + 1
      icnst(2) = nint((rcnst(8)-rcnst(5))/rcnst(2)) + 1
      icnst(3) = nint((rcnst(9)-rcnst(6))/rcnst(3)) + 1

      return
      end

*-------------------------------------------------------------------------

      subroutine srt_fls(idcdf,srcfl,a,ispc,icnst,rcnst,nsrc)

      integer ispc,nsrc,idcdf(nsrc),icnst(0:nsrc),rcnst(0:nsrc)
      real a(ispc)
      character*(*) srcfl(nsrc)

      character*(80) ctmp
      logical sam_dpth

      do i=1,nsrc-1
      do j=nsrc-1,i,-1

        if (a(rcnst(j+1)).gt.a(rcnst(j))) then

          ctmp = srcfl(j)
          srcfl(j) = srcfl(j+1)
          srcfl(j+1) = ctmp

          itmp = idcdf(j)
          idcdf(j) = idcdf(j+1)
          idcdf(j+1) = itmp

          itmp = rcnst(j)
          rcnst(j) = rcnst(j+1)
          rcnst(j+1) = itmp

          itmp = icnst(j)
          icnst(j) = icnst(j+1)
          icnst(j+1) = itmp

        end if

      enddo
      enddo

 533  FORMAT('Error:  grids ',I2,' and ',I2,' not properly nested.')
 534  FORMAT('Error:  grids ',I2,' and ',I2,' have different depths.')

      if ((a(rcnst(0)+3).lt.a(rcnst(1)+3))
     >    .or.(a(rcnst(0)+4).lt.a(rcnst(1)+4))
     >       .or.(a(rcnst(0)+6).gt.a(rcnst(1)+6))
     >          .or.(a(rcnst(0)+7).gt.a(rcnst(1)+7))) then
        write(6,533) k,1
        stop
      end if

      do k=2,nsrc

        if ((a(rcnst(k)+3).lt.a(rcnst(1)+3))
     >      .or.(a(rcnst(k)+4).lt.a(rcnst(1)+4))
     >         .or.(a(rcnst(k)+6).gt.a(rcnst(1)+6))
     >            .or.(a(rcnst(k)+7).gt.a(rcnst(1)+7))) then
          write(6,533) k,1
          stop
        end if

        if (.not.sam_dpth(a(icnst(k)+2),a(icnst(1)+2),a(rcnst(k)+2),
     >                     a(rcnst(1)+2))) then
          write(6,534) k,1
          stop
        end if

      enddo

      return
      end

*-----------------------------------------------------------------------

      logical function sam_dpth(nzk,nz1,dzk,dz1)

      integer nzk,nz1
      real dzk,dz1

      logical requ

      sam_dpth = .true.
      if (nzk.ne.nz1) sam_dpth = .false.
      if (.not.requ(dzk,dz1)) sam_dpth = .false.

      return
      end 

*---------------------------------------------------------------------------

      subroutine get_sclrs(idcdf,asclr,sclrnm,nsclr,mxsclr)

      integer idcdf,nsclr,mxsclr
      real asclr(mxsclr)
      character*(80) sclrnm(mxsclr)

      integer ndims,nvars,natts,idrec,did1
      integer vtyp,did(4)
      character*(80) vnm

      call ncinq(idcdf,ndims,nvars,natts,idrec,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem with file inquiry.'
        stop
      end if

      did1 = ncdid(idcdf,'one',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for dimension one.'
        stop
      end if

      nsclr = 0

      do k=1,nvars

        call ncvinq(idcdf,k,vnm,vtyp,ndims,did,natts,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem with variable inquiry.'
          stop
        end if

        iflg = 1
        n = 1
 458    if ((iflg.eq.1).and.(n.le.ndims)) then
          if (did(n).ne.did1) then
            iflg = 0
          else 
            n = n+1
          end if
          goto 458
        end if

        if (iflg.eq.1) then
          nsclr = nsclr+1
          if (nsclr.gt.mxsclr) then
            write(6,*) 'Error:  too many scalars.  Increase mxsclr.'
            stop
          end if
          lgth = index(vnm,' ')-1 
          sclrnm(nsclr) = vnm(1:lgth)
          call ncvgt1(idcdf,k,1,asclr(nsclr),ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem reading scalar.'
            stop
          end if
        end if

      enddo

      return
      end

*-----------------------------------------------------------------------

      logical function sclr_err(asclr1,asclr,sclrnm1,sclrnm,
     >                             nsclr1,nsclr)

      integer nsclr1,nsclr
      real asclr1(nsclr1),asclr(nsclr)
      character*(80) sclrnm1(nsclr1),sclrnm(nsclr)

      logical requ

      sclr_err = .false.

      if (nsclr1.ne.nsclr) then

        sclr_err = .true.

      else

        k = 1
 48     if ((sclr_err.eq..false.).and.(k.le.nsclr1)) then
          
          lgth1 = index(sclrnm1(k),' ')-1

          n = 1
          ifnd = 0
 112      if ((ifnd.eq.0).and.(n.le.nsclr)) then
            lgth = index(sclrnm(n),' ')-1
            if (sclrnm(n)(1:lgth).eq.sclrnm1(k)(1:lgth1)) then
              ifnd = 1
            else
              n = n+1
            end if
            goto 112
          end if

          if (ifnd.eq.0) then
            sclr_err = .true.
          else if (.not.requ(asclr1(k),asclr(n))) then
            sclr_err = .true.
          else
            k = k+1
          end if

          goto 48

        end if

      end if

      return
      end

*-------------------------------------------------------------------------

      subroutine set_anfl(idcdf,nx,ny,nz,rcnst,sgz,wgz,asclr,sclrnm,
     >                     var,units,anfl,srcfl,nsrc,nsclr,nvar)

      include '/usr/local/include/netcdf.inc'

      integer idcdf,nx,ny,nz,nsclr,nvar,nsrc
      real  rcnst(9),sgz(nz-1),wgz(nz),asclr(nsclr)
      character*(*) sclrnm(nsclr),var(nvar),units(nvar),anfl
      character*(*) srcfl(nsrc)

      real xmin,ymin,zmin 
      character*(80) cdlfl,command,command2  
      character*(1) num

      lgth = index(anfl,' ')-1

      cdlfl(1:len(cdlfl)) = ' '
      cdlfl = anfl(1:lgth)//'.cdl'

      iunit = 48

      open(unit=iunit,file=cdlfl(1:lgth+4),status='new')

      write(iunit,*) 'netcdf an_out{'
      write(iunit,*)

      write(iunit,*) 'dimensions:'
      write(iunit,*) 'nx=',nx-1,';'
      write(iunit,*) 'ny=',ny-1,';'
      write(iunit,*) 'nz=',nz-1,';'
      write(iunit,*) 'nxp1=',nx,';'
      write(iunit,*) 'nyp1=',ny,';'
      write(iunit,*) 'nzp1=',nz,';'
      write(iunit,*) 'time=UNLIMITED;'
      write(iunit,*) 'one=',1,';'
      write(iunit,*)

      write(iunit,*) 'variables:'
      write(iunit,*)

      do k=1,nsclr
        lgth = index(sclrnm(k),' ')-1
        write(iunit,*) 'float ',sclrnm(k)(1:lgth),'(one);'
        write(iunit,*) sclrnm(k)(1:lgth),':no_button=1;'
        write(iunit,*)
      enddo

      write(iunit,*) 'float time(time);'
      write(iunit,*) 'time:units="s";'
      write(iunit,*)
      write(iunit,*) 'float sgz(nz);'
      write(iunit,*) 'sgz:units="m";'
      write(iunit,*) 'sgz:def="staggered grid pts in z";'
      write(iunit,*) 'sgz:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float wgz(nzp1);'
      write(iunit,*) 'wgz:units="m";'
      write(iunit,*) 'wgz:def="unstaggered grid pts in z";'
      write(iunit,*) 'wgz:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_p(ny,nx);'
      write(iunit,*) 'zbot_p:units="m";'
      write(iunit,*) 'zbot_p:def="terrain at p points";'
      write(iunit,*) 'zbot_p:x_min=',rcnst(4)+0.5*rcnst(1),';'
      write(iunit,*) 'zbot_p:y_min=',rcnst(5)+0.5*rcnst(2),';'
      write(iunit,*) 'zbot_p:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_u(ny,nxp1);'
      write(iunit,*) 'zbot_u:units="m";'
      write(iunit,*) 'zbot_u:def="terrain at u points";'
      write(iunit,*) 'zbot_u:x_min=',rcnst(4),';'
      write(iunit,*) 'zbot_u:y_min=',rcnst(5)+0.5*rcnst(2),';'
      write(iunit,*) 'zbot_u:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_v(nyp1,nx);'
      write(iunit,*) 'zbot_v:units="m";'
      write(iunit,*) 'zbot_v:def="terrain at v points";'
      write(iunit,*) 'zbot_v:x_min=',rcnst(4)+0.5*rcnst(1),';'
      write(iunit,*) 'zbot_v:y_min=',rcnst(5),';'
      write(iunit,*) 'zbot_v:no_button=1;'
      write(iunit,*)

      do k=1,nvar

        xmin = rcnst(4) + 0.5*rcnst(1)
        ymin = rcnst(5) + 0.5*rcnst(2)      
        zmin = sgz(1)
 
        lgth = index(var(k),' ')-1

        if (var(k)(lgth:lgth).eq.'0') then
          if (var(k)(1:1).eq.'u') then
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(one,nz,ny,nxp1);'
            xmin = rcnst(4)
          else if (var(k)(1:1).eq.'v') then
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(one,nz,nyp1,nx);'
            ymin = rcnst(5)
          else if (var(k)(1:1).eq.'w') then
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(one,nzp1,ny,nx);'
            zmin = wgz(1)
          else
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(one,nz,ny,nx);'
          end if
        else
          if (var(k)(1:1).eq.'u') then
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(time,nz,ny,nxp1);'
            xmin = rcnst(4)
          else if (var(k)(1:1).eq.'v') then
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(time,nz,nyp1,nx);'
            ymin = rcnst(5)
          else if (var(k)(1:1).eq.'w') then
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(time,nzp1,ny,nx);'
            zmin = wgz(1)
          else
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(time,nz,ny,nx);'
          end if
        end if

        lgth2 = index(units(k),' ')-1
        if (lgth2.eq.0) lgth2 = 1
        write(iunit,*) var(k)(1:lgth),':units="',
     >              units(k)(1:lgth2),'";'     

        write(iunit,*) var(k)(1:lgth),':x_min=',xmin,';'
        write(iunit,*) var(k)(1:lgth),':y_min=',ymin,';'
        write(iunit,*) var(k)(1:lgth),':z_min=',zmin,';'

        if (var(k)(lgth:lgth).eq.'0') then
          write(iunit,*) var(k)(1:lgth),':no_button=1;'
        end if

        write(iunit,*) 

      enddo

      write(iunit,*) '//global attributes:'
      write(iunit,*)
      write(iunit,*) ':x_min=',rcnst(4),';'
      write(iunit,*) ':x_max=',rcnst(7),';'
      write(iunit,*) ':x_delta=',rcnst(1),';'
      write(iunit,*) ':x_units="m";'
      write(iunit,*) ':x_label="x";'
      write(iunit,*) ':x_display_units="km";'
      write(iunit,*)
      write(iunit,*) ':y_min=',rcnst(5),';'
      write(iunit,*) ':y_max=',rcnst(8),';'
      write(iunit,*) ':y_delta=',rcnst(2),';'
      write(iunit,*) ':y_units="m";'
      write(iunit,*) ':y_label="y";'
      write(iunit,*) ':y_display_units="km";'
      write(iunit,*)
      write(iunit,*) ':z_min=',rcnst(6),';'
      write(iunit,*) ':z_max=',rcnst(9),';'
      write(iunit,*) ':z_delta=',rcnst(3),';'
      write(iunit,*) ':z_units="m";'
      write(iunit,*) ':z_label="z";'
      write(iunit,*) ':z_display_units="km";'
      write(iunit,*)

 123  FORMAT(I1)

      lgth = index(anfl,' ')-1
      write(iunit,*) ':runname="a posteriori analysis ',
     >                     anfl(1:lgth),'";'
      do k=1,nsrc
        lgths = index(srcfl(k),' ')-5
        write(num,123) k
        write(iunit,*) ':source',num,'="',srcfl(k)(1:lgths),'";'
      enddo
      write(iunit,*)
      write(iunit,*) '}'

      close(iunit)

      anfl = anfl(1:lgth)//'.cdf'

      command(1:len(command)) = ' '
      write(command,*) 'ncgen -o ',anfl(1:lgth+4),' ',
     >                         cdlfl(1:lgth+4)
      call system(command)
      write(command2,*) 'rm -f ',cdlfl(1:lgth+4)
      call system(command2)

      idcdf = ncopn(anfl(1:lgth+4),NCWRITE,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem opening netcdf file.'
        stop
      end if

      do k=1,nsclr
        lgth = index(sclrnm(k),' ')-1
        idvar = ncvid(idcdf,sclrnm(k)(1:lgth),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem getting variable id for ',
     >                                 'scalar.'
          stop
        end if
        call ncvpt1(idcdf,idvar,1,asclr(k),ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem writing scalar.'
          stop
        end if
      enddo

      idvar = ncvid(idcdf,'sgz',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for sgz.'
        stop
      end if

      call ncvpt(idcdf,idvar,1,nz-1,sgz,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem writing sgz.'
        stop
      end if

      idvar = ncvid(idcdf,'wgz',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for wgz.'
        stop
      end if

      call ncvpt(idcdf,idvar,1,nz,wgz,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem writing wgz.'
        stop
      end if

      return
      end

*----------------------------------------------------------------------

      subroutine intrp_info(rcnst0,rcnst,isrc,isrcu,jsrc,jsrcv,
     >                   ian,ianu,jan,janv,nx0,ny0,nx,ny)

      integer nx0,ny0,nx,ny,isrc(nx0),isrcu(nx0),jsrc(ny0),jsrcv(ny0)
      real rcnst0(9),rcnst(9),ian(nx0),ianu(nx0),jan(ny0),janv(ny0)

      real x,xu,y,yv,xsrc,xsrcu,ysrc,ysrcv

      k = 1
      ku = 1

      do i=1,nx0

        x = rcnst0(4) + (i-0.5)*rcnst0(1)
        xu = rcnst0(4) + (i-1)*rcnst0(1)

        if ((x.lt.rcnst(4)+0.5*rcnst(1))
     >                  .or.(x.gt.rcnst(7)-0.5*rcnst(1))) then
          isrc(i) = 0
        else
          ifnd = 0
 870      if ((ifnd.eq.0).and.(k.lt.nx-1)) then
            xsrc = rcnst(4) + (k-0.5)*rcnst(1)
            if ((xsrc.le.x).and.(xsrc+rcnst(1).ge.x)) then
              ifnd = 1
            else
              k = k+1 
            end if
            goto 870
          end if
          if (ifnd.eq.0) then
            write(6,*) 'Error:  no source point found in ',
     >                                 'interp_info.'
            stop 
          else
            isrc(i) = k
            ian(i) = (x-xsrc)/rcnst(1)
          end if
        end if

        if ((xu.lt.rcnst(4)).or.(xu.gt.rcnst(7))) then
          isrcu(i) = 0
        else
          ifnd = 0
 872      if ((ifnd.eq.0).and.(ku.lt.nx)) then
            xsrcu = rcnst(4) + (ku-1)*rcnst(1)
            if ((xsrcu.le.xu).and.(xsrcu+rcnst(1).ge.xu)) then
              ifnd = 1
            else
              ku = ku+1 
            end if
            goto 872
          end if
          if (ifnd.eq.0) then
            write(6,*) 'Error:  no source point found in ',
     >                                    'intrp_info.'
            stop
          else
            isrcu(i) = ku
            ianu(i) = (xu-xsrcu)/rcnst(1)
          end if
        end if

      enddo

      k = 1
      kv = 1

      do j=1,ny0

        y = rcnst0(5) + (j-0.5)*rcnst0(2)
        yv = rcnst0(5) + (j-1)*rcnst0(2)

        if ((y.lt.rcnst(5)+0.5*rcnst(2))
     >              .or.(y.gt.rcnst(8)-0.5*rcnst(2))) then
          jsrc(j) = 0
        else
          ifnd = 0
 874      if ((ifnd.eq.0).and.(k.lt.ny-1)) then
            ysrc = rcnst(5) + (k-0.5)*rcnst(2)
            if ((ysrc.le.y).and.(ysrc+rcnst(2).ge.y)) then
              ifnd = 1
            else
              k = k+1
            end if
            goto 874
          end if
          if (ifnd.eq.0) then
            write(6,*) 'Error:  no source point found in ',
     >                             'interp_info.'
            stop
          else
            jsrc(j) = k
            jan(j) = (y-ysrc)/rcnst(2)
          end if
        end if
        
        if ((yv.lt.rcnst(5)).or.(yv.gt.rcnst(8))) then
          jsrcv(j) = 0
        else
          ifnd = 0
 876      if ((ifnd.eq.0).and.(kv.lt.ny)) then
            ysrcv = rcnst(5) + (kv-1)*rcnst(2)
            if ((ysrcv.le.yv).and.(ysrcv+rcnst(2).ge.yv)) then
              ifnd = 1
            else
              kv = kv+1
            end if
            goto 876
          end if
          if (ifnd.eq.0) then
            write(6,*) 'Error:  no source point found in ',
     >                             'interp_info.'
            stop
          else
            jsrcv(j) = kv
            janv(j) = (yv-ysrcv)/rcnst(2)
          end if
        end if

      enddo

      return
      end

*----------------------------------------------------------------------
 
      subroutine read_var(idcdf,var,nx,ny,nz,vsrc,itime)

      integer idcdf,nx,ny,nz,itime
      real vsrc(nx*ny*nz)
      character*(*) var

      integer istrt(4),ilen(4),ist2(2),iln2(2)

      lgth = index(var,' ')-1
      idvar = ncvid(idcdf,var(1:lgth),ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting variable id.'
        stop
      end if

      istrt(1) = 1
      istrt(2) = 1
      istrt(3) = 1
      istrt(4) = itime
      ist2(1) = 1
      ist2(2) = 1

      ilen(1) = nx-1
      ilen(2) = ny-1
      ilen(3) = nz-1
      ilen(4) = 1
      iln2(1) = nx-1
      iln2(2) = ny-1

      if (var(1:4).eq.'zbot') then
        if (var(lgth:lgth).eq.'u') then
          iln2(1) = nx
        else if (var(lgth:lgth).eq.'v') then
          iln2(2) = ny
        end if
        call ncvgt(idcdf,idvar,ist2,iln2,vsrc,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem reading terrain variable.'
          stop
        end if
      else
        if (var(1:1).eq.'u') then
          ilen(1) = nx
        else if (var(1:1).eq.'v') then
          ilen(2) = ny
        else if (var(1:1).eq.'w') then
          ilen(3) = nz
        end if
        call ncvgt(idcdf,idvar,istrt,ilen,vsrc,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem reading variable.'
          stop
        end if
      end if 

      return
      end

*---------------------------------------------------------------------

      subroutine interp_var(vsrc,van,isrc,isrcu,jsrc,jsrcv,
     >                        ian,ianu,jan,janv,nx0,ny0,nz0,
     >                        nx,ny,nz,var)

      integer nx0,ny0,nz0,nx,ny,nz,isrc(nx0),isrcu(nx0)
      integer jsrc(ny0),jsrcv(ny0)
      real ian(nx0),ianu(nx0),jan(ny0),janv(ny0)
      real vsrc(nx*ny*nz),van(nx0*ny0*nz0)
      character*(*) var

      integer nxs,nys,nzs,nxa,nya,nza

      nxs = nx-1
      nys = ny-1
      nzs = nz-1
      nxa = nx0-1
      nya = ny0-1
      nza = nz0-1

      lgth = index(var,' ')-1

      if (var(1:4).eq.'zbot') then
        nzs = 1
        nza = 1  
        if (var(lgth:lgth).eq.'u') then
          nxs = nx
          nxa = nx0
        else if (var(lgth:lgth).eq.'v') then
          nys = ny
          nya = ny0
        end if
      else
        if (var(1:1).eq.'u') then
          nxs = nx
          nxa = nx0
        else if (var(1:1).eq.'v') then
          nys = ny
          nya = ny0
        else if (var(1:1).eq.'w') then
          nzs = nz
          nza = nz0
        end if
      end if

      if (nxs.ne.nx-1) then
        call interp(vsrc,van,isrcu,jsrc,ianu,jan,nxa,nya,nza,
     >                           nxs,nys,nzs)
      else if (nys.ne.ny-1) then
        call interp(vsrc,van,isrc,jsrcv,ian,janv,nxa,nya,nza,
     >                           nxs,nys,nzs)
      else
        call interp(vsrc,van,isrc,jsrc,ian,jan,nxa,nya,nza,
     >                           nxs,nys,nzs)
      end if

      return
      end

*--------------------------------------------------------------------------

      subroutine interp(vsrc,van,isrc,jsrc,ian,jan,nxa,nya,nza,
     >                               nxs,nys,nzs) 

      integer nxa,nya,nza,nxs,nys,nzs,isrc(nxa),jsrc(nya)
      real ian(nxa),jan(nya),vsrc(nxs,nys,nzs),van(nxa,nya,nza)

      real wt(4)   

      do k=1,nza
      do j=1,nya
      do i=1,nxa

        if ((isrc(i).ne.0).and.(jsrc(j).ne.0)) then
          wt(1) = (1.-ian(i))*(1.-jan(j))
          wt(2) = (1.-ian(i))*jan(j)
          wt(3) = ian(i)*jan(j)
          wt(4) = ian(i)*(1.-jan(j))
          van(i,j,k) = wt(1)*vsrc(isrc(i),jsrc(j),k)
     >            + wt(2)*vsrc(isrc(i),jsrc(j)+1,k)
     >            + wt(3)*vsrc(isrc(i)+1,jsrc(j)+1,k)
     >            + wt(4)*vsrc(isrc(i)+1,jsrc(j),k)
        end if

      enddo
      enddo
      enddo

      return
      end       

*--------------------------------------------------------------------

      subroutine write_var(idcdf,var,nx,ny,nz,van,itime)

      integer idcdf,nx,ny,nz,itime
      real van(nx*ny*nz)
      character*(*) var

      integer istrt(4),ilen(4),ist2(2),iln2(2)

      lgth = index(var,' ')-1

      idvar = ncvid(idcdf,var(1:lgth),ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting variable id.'
        stop
      end if

      istrt(1) = 1
      istrt(2) = 1
      istrt(3) = 1
      istrt(4) = itime
      ist2(1) = 1
      ist2(2) = 1

      ilen(1) = nx-1
      ilen(2) = ny-1
      ilen(3) = nz-1
      ilen(4) = 1
      iln2(1) = nx-1
      iln2(2) = ny-1

      if (var(1:4).eq.'zbot') then
        if (var(lgth:lgth).eq.'u') then
          iln2(1) = nx
        else if (var(lgth:lgth).eq.'v') then
          iln2(2) = ny
        end if
        call ncvpt(idcdf,idvar,ist2,iln2,van,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem writing terrain variable.'
          stop
        end if
      else
        if (var(1:1).eq.'u') then
          ilen(1) = nx
        else if (var(1:1).eq.'v') then
          ilen(2) = ny
        else if (var(1:1).eq.'w') then
          ilen(3) = nz
        end if
        call ncvpt(idcdf,idvar,istrt,ilen,van,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem writing variable.'
          stop
        end if
      end if

      return
      end

*------------------------------------------------------------------------
 
      integer function i_val(inval)

      integer inval

      i_val = inval

      return
      end

*------------------------------------------------------------------------

      subroutine get1da(idcdf,sgz,wgz,time,nz,nt)

      integer nz,nt,idcdf
      real sgz(nz),wgz(nz),time(*)

      character*(80) dum

      idvar = ncvid(idcdf,'sgz',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for sgz.'
        stop
      end if

      call ncvgt(idcdf,idvar,1,nz-1,sgz,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem reading sgz.'
        stop
      end if

      idvar = ncvid(idcdf,'wgz',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for wgz.'
        stop
      end if

      call ncvgt(idcdf,idvar,1,nz,wgz,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem reading wgz.'
        stop
      end if

      idtime = ncdid(idcdf,'time',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for dimension time.'
        stop
      end if

      call ncdinq(idcdf,idtime,dum,nt,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem inquiring about time dimension.'
        stop
      end if

      idvar = ncvid(idcdf,'time',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting id for time variable.'
        stop
      end if

      call ncvgt(idcdf,idvar,1,nt,time,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem reading time variable.'
        stop
      end if

      return
      end

*--------------------------------------------------------------------------

      integer function int_time(outtm,srctm,ntime)

      integer ntime
      real outtm,srctm(ntime)

      logical requ

      int_time = 0

      k = 1
      ifnd = 0
 224  if ((ifnd.eq.0).and.(k.le.ntime)) then
        if (requ(outtm,srctm(k))) then
          ifnd = 1
        else
          k = k+1
        end if
        goto 224 
      end if

      if (ifnd.eq.0) then
        write(6,*) 'Error:  requested time not found.'
        stop
      else
        int_time = k
      end if

      return
      end
      
*-------------------------------------------------------------------------

      logical function requ(val1,val2)

      real val1,val2,eps,diff,big	

      parameter (eps=1.e-4)

      big = amax1(abs(val1),abs(val2))

      if (big.eq.0.) then
        requ = .true.
        return
      end if

      diff = abs((val1-val2)/big)

      if (diff.le.eps) then
        requ = .true.
      else
        requ = .false.
      end if

      return
      end
       
      



          

      

      
