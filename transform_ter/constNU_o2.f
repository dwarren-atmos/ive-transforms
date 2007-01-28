
      program constNU_o2

      include '/usr/local/include/netcdf.inc'

      integer ispc,nspc
      parameter (ispc=275000000,nspc=20)

      real a(ispc)

      integer maxslb
      parameter (maxslb=501)

      integer idcdf1,idcdf2,isclr,isgz,iwgz,izs_p,izs_u,izs_v
      integer itran(3),iwrk(3),iter(3),ik,il,iwrk_x,iwrk_y,mptr(nspc)
      integer fptr,fptr1,fptr2,fptr3,fptr4,iptr,igu,iu1,igv,iv1,igw
      integer iw1,igb,ib1,ib0,id11,id12,id13,id22,id23,id33,ikm,idgu
      integer idgv,idgb,ius,ivs,iws,idgus,idgvs,igws,igbs,igxi,iget
      integer igze,iw2,iwf,idwf,ipf,inf1,inf2,inf3,inf4,inf5,iuf,ivf
      integer iu,iv,iw,i1,i2,izsr,izsf,idz,ifld,ift,isgzt,iwgzt
      integer iGs,igh,ircof,isz(maxslb),iwz(maxslb),nslb,rslb,ntop
      integer igut,igvt,igus,igvs,iuzt,ivzt,iwzt,icalc,is,ifld1,ifld2
      integer nx,ny,nz,nzt,nxout,nyout,nzout,nvout,itopo,ivrt,nvrt
      integer iq1,ir1,is1,iq2,ir2,is2,is3,igs2,iwrk4,itran4
      real dx,dy,x0,y0,zt,ztop,pi,vel,mix,xout,yout,twod
      real time1,time2,rtmp(2)
      character*(80) outfl1,outfl2,ctmp

      integer nvar,nsclr
      parameter (nvar=28,nsclr=22)

      real dz
      common /dzblk/ dz

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      character*(80) var(nvar),units(nvar),sclr(nsclr)
      data var / 'u1','u2','v1','v2','w1','w2','p1','p2','km1','b1',
     >           'b2','zeta1','zeta2','xi1','xi2','eta1','eta2',
     >           'u0','v0','b0','p0','q1','q2','r1','r2','s1',
     >           's2','s3' /
      data units / 'm/s','m/s','m/s','m/s','m/s','m/s','m*m/s/s',
     >             'm*m/s/s','m*m/s','m/s/s','m/s/s','m','m',
     >             'm','m','m','m','m/s','m/s','m/s/s','m*m/s/s',
     >             '1/s','1/s','1/s','1/s','1/s','1/s','1/s' /
      data sclr / 'f_cor','ztop','type','h','a1','a2','b1','b2',
     >            'y_len','xc','yc','N','U0','V0','ps','Ts',
     >            'rhos','zt','mix','Cm','rprandl','twod' /

 886  FORMAT(A80)

      time1 = etime(rtmp)

      if (nspc.lt.20) then
        write(6,*)
        write(6,*) 'Error:  need space for at least 20 variables.'
        write(6,*) '        Increase nspc.'
        write(6,*)
        stop
      end if 

      isclr = 1

      iunit = 32

      open (unit=iunit,file='constNU_o2.in',status='old')

      read(iunit,*) x0,y0
      read(iunit,*) nx,ny,nz
      read(iunit,*) dx,dy,dz
      read(iunit,*) xout,yout
      read(iunit,*) nxout,nyout,nzout
      read(iunit,*) a(isclr+22-1)
      read(iunit,*) a(isclr+3-1)
      read(iunit,*) a(isclr+5-1),a(isclr+7-1),a(isclr+4-1)
      read(iunit,*) a(isclr+6-1),a(isclr+8-1),a(isclr+9-1)
      read(iunit,*) a(isclr+10-1),a(isclr+11-1)
      read(iunit,*) a(isclr+1-1)
      read(iunit,*) a(isclr+12-1)
      read(iunit,*) a(isclr+13-1)
      read(iunit,*) a(isclr+14-1)
      read(iunit,*) a(isclr+15-1)
      read(iunit,*) a(isclr+16-1)
      read(iunit,*) a(isclr+18-1)
      read(iunit,*) a(isclr+19-1)
      read(iunit,*) a(isclr+20-1)
      read(iunit,*) a(isclr+21-1)
      read(iunit,*) ivrt
      read(iunit,*) itopo
      read(iunit,886) outfl2

      close (unit=iunit)

      if (ivrt.eq.1) then
        if (a(isclr+1-1).eq.0.) then
          nvout = nvar
          nvrt = 7
        else
          nvout = nvar-1
          nvrt = 6
        end if
      else
        nvout = nvar-7
        nvrt = 0
      end if

      mix = a(isclr+19-1)
      twod = a(isclr+22-1)
      if (itopo.eq.1) then
        lgth = index(outfl2,' ') - 1
        outfl1 = outfl2(1:lgth)//'_tmp' 
      else
        outfl1 = outfl2
      end if

      if ((mod(nx,2).eq.0).or.(mod(ny,2).eq.0)) then
        write(6,*)
        write(6,*) 'Error:  nx-1 and ny-1 must be even.'
        write(6,*)
        stop
      else if ((mix.ne.0.).and.(nz.lt.7)) then
        write(6,*)
        write(6,*) 'Error:  need a layer depth of at least seven ',
     >             '        when using mixing.'
        write(6,*)
        stop
      else if ((twod.ne.0.).and.(ny.ne.7)) then
        write(6,*) 
        write(6,*) 'Error:  set ny to seven when solving in 2D.'
        write(6,*)
        stop
      end if  

      iout = nint(xout/dx) + 1
      jout = nint(yout/dy) + 1
      if ((iout+nxout-1.gt.nx).or.(jout+nyout-1.gt.ny)) then
        write(6,*) 
        write(6,*) 'Error:  output domain exceeds calculation domain.'
        write(6,*)
        stop
      end if
      xout = (iout-1)*dx
      yout = (jout-1)*dy

      pi = 2.*asin(1.)
      vel = sqrt(a(isclr+13-1)**2.+a(isclr+14-1)**2.)
      a(isclr+2-1) = dz*(nzout-1)
      a(isclr+17-1) = a(isclr+15-1)/287.04/a(isclr+16-1)
      zt = a(isclr+18-1)*2.*pi*vel/a(isclr+12-1)
      nzt = nint(zt/dz)+1
      zt = (nzt-1)*dz

      if (nzt.lt.nzout) then
        write(6,*) 'Error:  upper radiation boundary less ',
     >                      'than output domain height.'
        stop
      end if

      call set_sclrblk(a(isclr),nsclr)

      isgzt = isclr + nsclr
      iwgzt = isgzt + nzt + 2
      isgz = iwgzt + nzt + 2
      iwgz = isgz + nzout
      izs_p = iwgz + nzout
      izs_u = izs_p + 2*nx*ny 
      izs_v = izs_u + 2*nx*ny
      igus = izs_v + 2*nx*ny
      igvs = igus + 2*nx*ny
      igut = igvs + 2*nx*ny
      igvt = igut + 2*nx*ny
      itran(1) = igvt + 2*nx*ny
      itran(2) = itran(1) + 2*nx*ny
      itran(3) = itran(2) + 2*nx*ny
      iwrk(1) = itran(3) + 2*nx*ny
      iwrk(2) = iwrk(1) + 2*nx*ny
      iwrk(3) = iwrk(2) + 2*nx*ny
      iter(1) = iwrk(3) + 2*nx*ny
      iter(2) = iter(1) + nx*ny
      iter(3) = iter(2) + nx*ny
      ik = iter(3) + nx*ny
      il = ik + nx
      iwrk_x = il + ny
      iwrk_y = iwrk_x + 4*nx + 15
      ircof = iwrk_y + 4*ny + 15
      igh = ircof + 2*nx*ny
      iGs = igh + 2*nx*ny
     
      mptr(1) = iGs + 2*nx*ny
      k = 2 
 89   if (k.le.nspc) then
        mptr(k) = mptr(k-1) + nx*ny*nz
        k = k+1
        goto 89
      end if

      fptr = mptr(nspc) + nx*ny*nz
      if (fptr.gt.ispc) then
        write(6,*) 
        write(6,*) 'Error:  not enough space!  Increase ispc.'
        write(6,*) '        needed:  ',fptr
        write(6,*) '        allocated:  ',ispc
        write(6,*)
        stop
      end if

      call get_ter(a(isclr),a(izs_p),a(izs_u),a(izs_v),x0,y0,dx,dy,
     >                       nx,ny,nsclr)

      do ij=1,nx*ny
        a(iter(1)+ij-1) = a(izs_p+ij-1)
        a(iter(2)+ij-1) = a(izs_u+ij-1)
        a(iter(3)+ij-1) = a(izs_v+ij-1)
      enddo

      a(iwgz+1-1) = 0.
      a(isgz+1-1) = 0.5*dz
      do k=2,nzout
        a(isgz+k-1) = a(isgz+k-2) + dz
        a(iwgz+k-1) = a(iwgz+k-2) + dz
      enddo

      call set_outfl(idcdf1,nxout,nyout,nzout,dx,dy,dz,xout,yout,
     >                a(isgz),a(iwgz),a(isclr),0.,sclr,var,units,
     >                outfl1,nsclr,nvout)

      call write_2d(idcdf1,a(izs_p),a(itran(1)),'zbot_p',iout,jout,
     >                  nx,ny,nxout,nyout,0,0)
      call write_2d(idcdf1,a(izs_u),a(itran(1)),'zbot_u',iout,jout,
     >                  nx,ny,nxout,nyout,1,0)
      call write_2d(idcdf1,a(izs_v),a(itran(1)),'zbot_v',iout,jout,
     >                  nx,ny,nxout,nyout,0,1)

      call calc_wavnm(a(ik),dx,nx-1) 
      call calc_wavnm(a(il),dy,ny-1)

      call cffti(nx-1,a(iwrk_x))
      call cffti(ny-1,a(iwrk_y))
      call trnsfrm2d(a(izs_p),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,0,0)
      call trnsfrm2d(a(izs_u),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,1,0)
      call trnsfrm2d(a(izs_v),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,0,1)

      fptr1 = mptr(14)
      fptr2 = mptr(16)
      fptr3 = mptr(18)

      ius = mptr(1)
      ivs = mptr(2)
      iws = mptr(3)

      do k=1,3
        if (k.eq.1) then
          iptr = ius
          ctmp = 'u1'
        else if (k.eq.2) then 
          iptr = ivs
          ctmp = 'v1'
        else 
          iptr = iws
          ctmp = 'w1'
        end if
        call calc_o1(a(iptr),a(fptr1),a(izs_p),a(itran(1)),a(iwrk(1)),
     >            0.,a(iwrk_x),a(iwrk_y),a(ik),a(il),dx,dy,ctmp,
     >            nx,ny,1,0,0)
      enddo

      idgus = mptr(4)
      idgvs = mptr(6)
      igws = mptr(8)
      igbs = mptr(10)

      do ij=1,nx*ny
        a(idgus+ij-1) = 0.
        a(idgvs+ij-1) = 0.
        a(igws+ij-1) = 0.
        a(igbs+ij-1) = 0.
        a(igus+ij-1) = 0.
        a(igvs+ij-1) = 0.
      enddo

      do k=1,6
        if (k.eq.1) then
          iptr = idgus
          ctmp = 'dgu'
          idz = 1
        else if (k.eq.2) then 
          iptr = idgvs
          ctmp = 'dgv'
          idz = 1
        else if (k.eq.3) then 
          iptr = igws
          ctmp = 'gw'
          idz = 0
        else if (k.eq.4) then
          iptr = igbs
          ctmp = 'gb'
          idz = 0
        else if (k.eq.5) then
          iptr = igus
          ctmp = 'gu'
          idz = 0
        else
          iptr = igvs
          ctmp = 'gv'
          idz = 0
        end if
        call calc_advc(ctmp,a(iptr),a(ius),a(ivs),a(iws),a(fptr1),
     >         a(fptr2),a(fptr3),a(izs_p),a(itran(1)),a(itran(2)),
     >         a(itran(3)),a(iwrk(1)),a(iwrk(2)),a(iwrk(3)),
     >         a(iwrk_x),a(iwrk_y),a(ik),a(il),0.,nx,ny,1,0,0,idz)
        call trnsfrm2d(a(iptr),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,
     >                                       0,0)
      enddo

      call calc_Gfac(a(iGs),a(idgus),a(idgvs),a(igws),a(igbs),
     >                a(ik),a(il),nx,ny)
      call calc_gh(a(igh),a(ius),a(ivs),a(iter(1)),a(izs_p),
     >         a(itran(1)),a(itran(2)),a(itran(3)),a(iwrk(1)),
     >         a(iwrk(2)),a(iwrk(3)),a(iwrk_x),a(iwrk_y),a(ik),
     >         a(il),nx,ny)
      call trnsfrm2d(a(igh),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,0,0)

      iuzt = mptr(1)
      ivzt = mptr(2)
      iwzt = mptr(3)

      do k=1,3
        if (k.eq.1) then
          iptr = iuzt
          ctmp = 'u1'
        else if (k.eq.2) then 
          iptr = ivzt
          ctmp = 'v1'
        else 
          iptr = iwzt
          ctmp = 'w1'
        end if
        call calc_o1(a(iptr),a(fptr1),a(izs_p),a(itran(1)),a(iwrk(1)),
     >            zt,a(iwrk_x),a(iwrk_y),a(ik),a(il),dx,dy,ctmp,
     >            nx,ny,1,0,0)
      enddo

      do ij=1,nx*ny
        a(igut+ij-1) = 0.
        a(igvt+ij-1) = 0.
      enddo

      do k=1,2
        if (k.eq.1) then
          iptr = igut
          ctmp = 'gu'
          idz = 0
        else if (k.eq.2) then 
          iptr = igvt
          ctmp = 'gv'
          idz = 0
        end if
        call calc_advc(ctmp,a(iptr),a(iuzt),a(ivzt),a(iwzt),a(fptr1),
     >         a(fptr2),a(fptr3),a(izs_p),a(itran(1)),a(itran(2)),
     >         a(itran(3)),a(iwrk(1)),a(iwrk(2)),a(iwrk(3)),
     >         a(iwrk_x),a(iwrk_y),a(ik),a(il),zt,nx,ny,1,0,0,idz)
      enddo

      k=1
 91   if (k.le.nzt+2) then
        if (k.le.nzout) then
          a(iwgzt+k-1) = a(iwgz+k-1)
          a(isgzt+k-1) = a(isgz+k-1)
        else
          a(iwgzt+k-1) = a(iwgzt+k-2) + dz
          a(isgzt+k-1) = a(isgzt+k-2) + dz 
        end if
        k = k+1
        goto 91 
      end if

      rslb = 0
      if (mix.ne.0.) then
        isz(1) = isgzt + min(nz-6,nzout+2-6)
        iwz(1) = iwgzt + min(nz-6,nzout+2-6)
 97     if (isz(rslb+1)-isgzt+1.lt.nzt-3) then
          rslb = rslb + 1
          if (rslb.gt.maxslb-1) then
            write(6,*) 'Error:  too many slabs for radiation ',
     >                           'calculation!'
            stop
          end if
          isz(rslb+1) = isz(rslb) + nz-6
          iwz(rslb+1) = iwz(rslb) + nz-6
          goto 97
        end if
        ntop = nzt+2 - (isz(rslb)-isgzt+1) + 1
      else
        isz(1) = isgzt + min(nz,nzout-1)
        iwz(1) = iwgzt + min(nz,nzout-1)
 99     if (isz(rslb+1)-isgzt+1.lt.nzt) then
          rslb = rslb + 1
          if (rslb.gt.maxslb-1) then
            write(6,*) 'Error:  too many slabs for radiation ',
     >                          'calculation!'
            stop
          end if
          isz(rslb+1) = isz(rslb) + nz
          iwz(rslb+1) = iwz(rslb) + nz
          goto 99
        end if
        ntop = nzt - (isz(rslb)-isgzt+1)
      end if

      do ij=1,nx*ny
        a(ircof+ij-1) = 0.
      enddo

      idgu = mptr(1)
      iu1 = mptr(2)
      idgv = mptr(3)
      iv1 = mptr(4)
      igw = mptr(5)
      iw1 = mptr(6)
      igb = mptr(7)
      ib1 = mptr(8)
      idgb = mptr(9)
      igu = mptr(11)
      igv = mptr(12)
      id11 = mptr(13)
      id12 = mptr(14)
      id13 = mptr(15)
      id22 = mptr(16)
      id23 = mptr(17)
      id33 = mptr(18)
      ikm = mptr(19)

      if (rslb.gt.0) then
        call rad_coeff(a,ispc,ircof,iu1,iv1,iw1,ib1,igu,igv,igw,igb,
     >                  idgu,idgv,idgb,id11,id12,id13,id22,id23,id33,
     >                  ikm,fptr1,fptr2,fptr3,izs_p,izs_u,izs_v,
     >                  itran,iwrk,iwrk_x,iwrk_y,ik,il,isclr,isz(1),
     >                  iwz(1),rslb,ntop,dx,dy,dz,nx,ny,nz,1)
        if (mix.ne.0.) then
          call mix_zt(a(igut),a(igvt),a(igu),a(igv),nx,ny,ntop)
        end if
      end if

      isz(1) = isgzt
      iwz(1) = iwgzt
      nslb = 0
      if (mix.ne.0.) then
 102    if (isz(nslb+1)-isgzt+1.lt.nzout-3) then
          nslb = nslb+1
          isz(nslb+1) = isz(nslb) + nz-6
          iwz(nslb+1) = iwz(nslb) + nz-6
          goto 102
        end if
        ntop = nzout+2 - (isz(nslb)-isgzt+1) + 1
      else
 104    if (isz(nslb+1)-isgzt+1.lt.nzout) then
          nslb = nslb+1
          isz(nslb+1) = isz(nslb) + nz
          iwz(nslb+1) = iwz(nslb) + nz
          goto 104
        end if
        ntop = nzout - (isz(nslb)-isgzt+1)
      end if

      is = 1
 107  if (is.le.nslb) then
 
        if (is.eq.nslb) then
          nzs = ntop
        else
          nzs = nz
        end if
 
        idgu = mptr(1)
        iu1 = mptr(2)
        idgv = mptr(3)
        iv1 = mptr(4)
        igw = mptr(5)
        iw1 = mptr(6)
        igb = mptr(7)
        ib1 = mptr(8)
        idgb = mptr(9)
        igu = mptr(11)
        igv = mptr(12)
        id11 = mptr(13)
        id12 = mptr(14)
        id13 = mptr(15)
        id22 = mptr(16)
        id23 = mptr(17)
        id33 = mptr(18)
        ikm = mptr(19)

        fptr1 = mptr(14)
        fptr2 = mptr(16)
        fptr3 = mptr(18)

        do ijk=1,nx*ny*nzs
          a(igu+ijk-1) = 0.
          a(igv+ijk-1) = 0.
          a(igw+ijk-1) = 0.
          a(igb+ijk-1) = 0.
          a(idgu+ijk-1) = 0.
          a(idgv+ijk-1) = 0.
          a(idgb+ijk-1) = 0.
          a(ikm+ijk-1) = 0.
        enddo

        do k=1,4
          iu = 0
          iv = 0
          igz = isz(is)
          izsf = izs_p
          if (k.eq.1) then
            ctmp = 'u1'
            iu = 1
            izsf = izs_u
            iptr = iu1
          else if (k.eq.2) then
            ctmp = 'v1'
            iv = 1
            izsf = izs_v
            iptr = iv1
          else if (k.eq.3) then
            ctmp = 'w1'
            igz = iwz(is)
            iptr = iw1
          else
            ctmp = 'b1'
            iptr = ib1
          end if
          call calc_o1(a(iptr),a(fptr1),a(izsf),a(itran(1)),
     >                 a(iwrk(1)),a(igz),a(iwrk_x),a(iwrk_y),
     >                 a(ik),a(il),dx,dy,ctmp,nx,ny,nzs,iu,iv)
        enddo

        if (mix.ne.0.) then
 
          call calc_def(a(id11),a(id12),a(id13),a(id22),a(id23),
     >                 a(id33),a(iu1),a(iv1),a(iw1),a(isz(is)),
     >                 a(iwz(is)),dx,dy,nx,ny,nzs)
          call calc_km(a(id11),a(id12),a(id13),a(id22),a(id23),
     >                 a(id33),a(ib1),a(ikm),a(isz(is)),a(isclr+5-1),
     >                 a(isclr+20-1),a(isclr+21-1),nx,ny,nzs)
          call u_mix(a(id11),a(id12),a(id13),a(igu),a(idgu),
     >                    a(iwz(is)),a(isz(is)),dx,dy,nx,ny,nzs)
          call v_mix(a(id12),a(id22),a(id23),a(igv),a(idgv),
     >                    a(iwz(is)),a(isz(is)),dx,dy,nx,ny,nzs)
          call w_mix(a(id13),a(id23),a(id33),a(igw),a(isz(is)),dx,dy,
     >                             nx,ny,nzs)
          call b_mix(a(ikm),a(ib1),a(igb),a(idgb),a(isz(is)),
     >                  a(iwz(is)),a(isclr+21-1),dx,dy,nx,ny,nzs)
          if (rslb.eq.0) then
            call mix_zt(a(igut),a(igvt),a(igu),a(igv),nx,ny,nzs)
          end if

          if (is.ne.1) then
            nzs = nzs-6
            idgu = idgu + 3*(nx-1)*(ny-1)
            idgv = idgv + 3*(nx-1)*(ny-1)
            idgb = idgb + 3*(nx-1)*(ny-1)
            igw = igw + 3*(nx-1)*(ny-1)
            igb = igb + 3*(nx-1)*(ny-1)
            igu = igu + 3*(nx-1)*(ny-1)
            igv = igv + 3*(nx-1)*(ny-1)
            iu1 = iu1 + 3*nx*(ny-1)
            iv1 = iv1 + 3*(nx-1)*ny
            iw1 = iw1 + 3*(nx-1)*(ny-1)
            ib1 = ib1 + 3*(nx-1)*(ny-1)
            ikm = ikm + 3*(nx-1)*(ny-1)
            isz(is) = isz(is) + 3
            iwz(is) = iwz(is) + 3
          else
            nzs = nzs-3
          end if
          if (twod.ne.0.) then
            call set_2d(a(idgu),a(idgv),a(idgb),a(igw),a(igb),
     >                     a(igu),a(igv),a(ikm),a(igut),a(igvt),
     >                     nx,ny,nzs)
          end if

        end if

        do k=1,nvout-nvrt
          lgth = index(var(k),' ')-1
          if ((var(k)(lgth:lgth).eq.'1')
     >                 .or.(var(k)(lgth:lgth).eq.'0')) then
            iu = 0
            iv = 0
            izsf = izs_p
            izsr = iter(1)
            igz = isz(is)
            iptr = fptr2
            icalc = 1
            if (var(k)(1:1).eq.'u') then
              iu = 1
              izsf = izs_u
              izsr = iter(2)
              if (var(k)(1:lgth).eq.'u1') then
                iptr = iu1
                icalc = 0
              end if 
            else if (var(k)(1:1).eq.'v') then
              iv = 1
              izsf = izs_v
              izsr = iter(3)
              if (var(k)(1:lgth).eq.'v1') then
                iptr = iv1
                icalc = 0
              end if
            else if (var(k)(1:1).eq.'w') then
              igz = iwz(is)
              if (var(k)(1:lgth).eq.'w1') then
                iptr = iw1
                icalc = 0
              end if
            else if (var(k)(1:lgth).eq.'b1') then
              iptr = ib1
              icalc = 0
            else if (var(k)(1:lgth).eq.'km1') then
              iptr = ikm
              icalc = 0
            end if
            i1 = iu
            i2 = iv
            if (icalc.eq.1) then
              call calc_o1(a(iptr),a(fptr1),a(izsf),a(itran(1)),
     >                 a(iwrk(1)),a(igz),a(iwrk_x),a(iwrk_y),a(ik),
     >                 a(il),dx,dy,var(k),nx,ny,nzs,iu,iv)
            end if
            call write_3d(idcdf1,a(iptr),a(fptr1),a(izsr),a(igz),
     >             a(isclr+2-1),mix,var(k),nx,ny,nz,nzs,nxout,
     >             nyout,nzout,iout,jout,is,iu,iv,i1,i2)
          end if
        enddo

        iu1 = mptr(2)
        iv1 = mptr(4)
        iw1 = mptr(6)

        do k=1,3
          if (k.eq.1) then
            iptr = iu1
            ctmp = 'u1'
          else if (k.eq.2) then 
            iptr = iv1
            ctmp = 'v1'
          else 
            iptr = iw1
            ctmp = 'w1'
          end if
          call calc_o1(a(iptr),a(fptr1),a(izs_p),a(itran(1)),
     >            a(iwrk(1)),a(isz(is)),a(iwrk_x),a(iwrk_y),a(ik),
     >            a(il),dx,dy,ctmp,nx,ny,nzs,0,0)
        enddo

        do k=1,6+ivrt
          idz = 0
          if (k.eq.1) then
            iptr = idgu
            ctmp = 'dgu'
            idz = 1
          else if (k.eq.2) then
            iptr = idgv
            ctmp = 'dgv'
            idz = 1
          else if (k.eq.3) then
            iptr = igu
            ctmp = 'gu'
          else if (k.eq.4) then
            iptr = igv
            ctmp = 'gv'
          else if (k.eq.5) then
            iptr = igw
            ctmp = 'gw'
          else if (k.eq.6) then
            iptr = igb
            ctmp = 'gb'
          else if (k.eq.7) then
            iptr = idgb
            ctmp = 'dgb'
            idz = 1
          end if
          call calc_advc(ctmp,a(iptr),a(iu1),a(iv1),a(iw1),a(fptr1),
     >         a(fptr2),a(fptr3),a(izs_p),a(itran(1)),a(itran(2)),
     >         a(itran(3)),a(iwrk(1)),a(iwrk(2)),a(iwrk(3)),
     >         a(iwrk_x),a(iwrk_y),a(ik),a(il),a(isz(is)),nx,ny,nzs,
     >         0,0,idz)
c         call write_3d(idcdf1,a(iptr),a(fptr1),a(iter(1)),
c    >           a(isz(is)),a(isclr+2-1),mix,ctmp,nx,ny,nz,nzs,
c    >           nxout,nyout,nzout,iout,jout,is,0,0,0,0)     
        enddo

        call trnsfrm3d(a(igw),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,
     >                            nzs,0,0)
        call trnsfrm3d(a(igb),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,
     >                            nzs,0,0)
        call trnsfrm3d(a(idgu),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,
     >                            nzs,0,0)
        call trnsfrm3d(a(idgv),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,
     >                            nzs,0,0)

        if (is.eq.1) then
          itmp = nzs
          call rad_coeff(a,ispc,ircof,0,0,0,0,0,0,igw,igb,idgu,idgv,
     >                  0,0,0,0,0,0,0,0,0,0,0,izs_p,izs_u,izs_v,
     >                  itran,iwrk,iwrk_x,iwrk_y,ik,il,isclr,isz(1),
     >                  iwz(1),1,itmp,dx,dy,dz,nx,ny,nzs,0)
          call trnsfrm2d(a(igut),a(itran(1)),a(iwrk_x),a(iwrk_y),
     >                                 nx,ny,0,0)
          call trnsfrm2d(a(igvt),a(itran(1)),a(iwrk_x),a(iwrk_y),
     >                                 nx,ny,0,0)
          call rad_bc(a(ircof),a(igut),a(igvt),a(ik),a(il),zt,nx,ny)
          iwrk4 = igut
          itran4 = igvt
          call drag_calc(a(igh),a(ircof),a(izs_p),a(igus),a(igvs),
     >                  a(fptr1),a(fptr2),a(itran(1)),a(itran(2)), 
     >                  a(itran(3)),a(itran4),a(iwrk(1)),a(iwrk(2)),
     >                  a(iwrk(3)),a(iwrk4),a(iter(1)),a(iwrk_x),
     >                  a(iwrk_y),a(ik),a(il),dx,dy,nx,ny)
          ic1 = igus
          ic2 = igvs
          ic1w = iwrk4
          ic2w = itran4
          do ij=1,2*nx*ny
            a(ic1+ij-1) = 0.
            a(ic2+ij-1) = 0.
            a(ic1w+ij-1) = 0.
            a(ic2w+ij-1) = 0.
          enddo
        end if

        fptr1 = 0
        fptr2 = 0
        fptr3 = 0

        iw2 = mptr(13)
        iwf = mptr(15)
        idwf = mptr(17)
        fptr = mptr(19)

        call calc_w2(a(iw2),a(iwf),a(idwf),a(idgu),a(idgv),a(igb),
     >               a(igw),a(igh),a(ircof),a(iGs),a(ic1),a(ic2),
     >               a(ic1w),a(ic2w),a(itran(1)),a(iwrk(1)),
     >               a(iwrk(2)),a(iwrk_x),a(iwrk_y),a(ik),a(il),
     >               a(isz(is)),a(iwz(is)),dz,nx,ny,nzs,is)

        call write_3d(idcdf1,a(iw2),a(fptr),a(iter(1)),a(iwz(is)),
     >             a(isclr+2-1),mix,'w2',nx,ny,nz,nzs,nxout,nyout,
     >             nzout,iout,jout,is,0,0,0,0)

        if (ivrt.eq.1) then

          iq2 = mptr(13)
          ir2 = mptr(14)

          do k=1,2
            ift = fptr
            if (k.eq.1) then
              ctmp = 'q2'
              ifld = iq2
              inf1 = iwf
              inf2 = igb
              inf3 = igw
              inf4 = idgu
              inf5 = idgv
            else
              ctmp = 'r2'
              ifld = ir2
              inf1 = iwf
              inf2 = igb
              inf3 = igw
              inf4 = idgu
              inf5 = idgv
            end if
            call calc_o2(ctmp,a(ifld),a(ift),a(inf1),a(inf2),a(inf3),
     >                  a(inf4),a(inf5),a(itran(1)),a(iwrk(1)),
     >                  a(iwrk_x),a(iwrk_y),a(ik),a(il),nx,ny,nzs)
            call write_3d(idcdf1,a(ifld),a(fptr),a(iter(1)),
     >               a(isz(is)),a(isclr+2-1),mix,ctmp,nx,ny,nz,nzs,
     >               nxout,nyout,nzout,iout,jout,is,0,0,0,0)
          enddo

          iq1 = mptr(1)
          ir1 = mptr(2)
          is1 = mptr(3)
          is2 = mptr(4)
          igs2 = mptr(5)

          do ijk=1,nx*ny*nzs
            a(igs2+ijk-1) = 0.
          enddo

          do k=1,3
            if (k.eq.1) then
              iptr = iq1
              ctmp = 'q1'
            else if (k.eq.2) then
              iptr = ir1
              ctmp = 'r1'
            else
              iptr = is1
              ctmp = 's1'
            end if
            call calc_o1(a(iptr),a(fptr),a(izs_p),a(itran(1)),
     >            a(iwrk(1)),a(isz(is)),a(iwrk_x),a(iwrk_y),a(ik),
     >            a(il),dx,dy,ctmp,nx,ny,nzs,0,0)
            call write_3d(idcdf1,a(iptr),a(fptr),a(iter(1)),
     >               a(isz(is)),a(isclr+2-1),mix,ctmp,nx,ny,nz,nzs,
     >               nxout,nyout,nzout,iout,jout,is,0,0,0,0)
          enddo

          call calc_gs(a(igs2),a(iq1),a(ir1),a(is1),a(iq2),a(ir2),
     >                    a(is2),a(iwf),a(igb),a(idwf),a(idgb),
     >                    a(fptr),a(izs_p),a(itran(1)),a(iwrk(1)),
     >                    a(iwrk_x),a(iwrk_y),a(ik),a(il),a(isz(is)),
     >                    nx,ny,nzs,2)

          call trnsfrm3d(a(igs2),a(itran(1)),a(iwrk_x),a(iwrk_y),
     >                              nx,ny,nzs,0,0)

          ctmp = 's2'
          idum1 = mptr(7)
          idum2 = mptr(15)

          call calc_o2(ctmp,a(is2),a(fptr),a(igs2),a(idwf),a(idgb),
     >                a(idum1),a(idum2),a(itran(1)),a(iwrk(1)),
     >                a(iwrk_x),a(iwrk_y),a(ik),a(il),nx,ny,nzs)
          call write_3d(idcdf1,a(is2),a(fptr),a(iter(1)),
     >             a(isz(is)),a(isclr+2-1),mix,ctmp,nx,ny,nz,nzs,
     >             nxout,nyout,nzout,iout,jout,is,0,0,0,0)

          if (a(isclr+1-1).eq.0.) then

            is3 = mptr(4)
            igs3 = mptr(5)
          
            call calc_gs(a(igs3),a(iq1),a(ir1),a(is1),a(iq2),a(ir2),
     >                    a(is2),a(iwf),a(igb),a(idwf),a(idgb),
     >                    a(fptr),a(izs_p),a(itran(1)),a(iwrk(1)),
     >                    a(iwrk_x),a(iwrk_y),a(ik),a(il),a(isz(is)),
     >                    nx,ny,nzs,3)

            call trnsfrm3d(a(igs3),a(itran(1)),a(iwrk_x),a(iwrk_y),
     >                               nx,ny,nzs,0,0)

            ctmp = 's3'
            idum = mptr(7)
            fptr1 = mptr(1)
            fptr2 = mptr(9)
            fptr3 = mptr(13)
         
            call calc_o2(ctmp,a(is3),a(fptr),a(igs3),a(fptr1),
     >                a(fptr2),a(fptr3),a(idum),a(itran(1)),
     >                a(iwrk(1)),a(iwrk_x),a(iwrk_y),a(ik),a(il),
     >                nx,ny,nzs)
            call write_3d(idcdf1,a(is3),a(fptr),a(iter(1)),
     >             a(isz(is)),a(isclr+2-1),mix,ctmp,nx,ny,nz,nzs,
     >             nxout,nyout,nzout,iout,jout,is,0,0,0,0)

          end if

        end if

        do ijk=1,nx*ny*nzs
          a(mptr(1)+ijk-1) = a(igu+ijk-1)
          a(mptr(3)+ijk-1) = a(igv+ijk-1) 
        enddo

        igu = mptr(1)
        igv = mptr(3)

        call trnsfrm3d(a(igu),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,
     >                           nzs,0,0)
        call trnsfrm3d(a(igv),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,
     >                           nzs,0,0)

        fptr1 = mptr(9)
        fptr2 = mptr(11)
        fptr3 = mptr(13)

        do k=1,4
          i1 = 0
          i2 = 0 
          izsr = iter(1)
          if (k.eq.1) then
            ctmp = 'b2'
            ifld = mptr(5)
            ift = mptr(19)
            inf1 = iwf
            inf2 = igb
            inf3 = fptr1
            inf4 = fptr2
            inf5 = fptr3
          else if (k.eq.2) then
            ctmp = 'p2'
            ipf = mptr(7)
            ifld = mptr(5)
            ift = ipf
            inf1 = idwf
            inf2 = igu
            inf3 = igv
            inf4 = fptr1
            inf5 = fptr2
          else if (k.eq.3) then
            ctmp = 'u2'
            iuf = mptr(17)
            ifld = mptr(5)
            ift = iuf
            inf1 = ipf
            inf2 = igu
            inf3 = igv
            inf4 = fptr1
            inf5 = fptr2
            i1 = 1
            izsr = iter(2)
          else if (k.eq.4) then
            ctmp = 'v2'
            ivf = mptr(19)
            ifld = mptr(5)
            ift = ivf
            inf1 = ipf
            inf2 = igu
            inf3 = igv
            inf4 = fptr1
            inf5 = fptr2
            i2 = 1
            izsr = iter(3)
          end if
          call calc_o2(ctmp,a(ifld),a(ift),a(inf1),a(inf2),a(inf3),
     >                  a(inf4),a(inf5),a(itran(1)),a(iwrk(1)),
     >                  a(iwrk_x),a(iwrk_y),a(ik),a(il),nx,ny,nzs)
          call write_3d(idcdf1,a(ifld),a(fptr1),a(izsr),a(isz(is)),
     >               a(isclr+2-1),mix,ctmp,nx,ny,nz,nzs,nxout,nyout,
     >               nzout,iout,jout,is,0,0,i1,i2)
        enddo

        iu1 = mptr(2)
        iv1 = mptr(4)
        iw1 = mptr(6)

        do k=1,3
          if (k.eq.1) then
            iptr = iu1
            ctmp = 'u1'
          else if (k.eq.2) then 
            iptr = iv1
            ctmp = 'v1'
          else 
            iptr = iw1
            ctmp = 'w1'
          end if
          call calc_o1(a(iptr),a(fptr1),a(izs_p),a(itran(1)),
     >            a(iwrk(1)),a(isz(is)),a(iwrk_x),a(iwrk_y),a(ik),
     >            a(il),dx,dy,ctmp,nx,ny,nzs,0,0)
        enddo

        igxi = mptr(1)
        iget = mptr(3)
        igze = mptr(5)

        do ijk=1,nx*ny*nzs
          a(igxi+ijk-1) = 0.
          a(iget+ijk-1) = 0.
          a(igze+ijk-1) = 0.
        enddo

        do k=1,3
          if (k.eq.1) then
            iptr = igxi
            ctmp = 'gxi'
          else if (k.eq.2) then
            iptr = iget
            ctmp = 'get'
          else if (k.eq.3) then
            iptr = igze
            ctmp = 'gze'
          end if
          call calc_advc(ctmp,a(iptr),a(iu1),a(iv1),a(iw1),a(fptr1),
     >         a(fptr2),a(fptr3),a(izs_p),a(itran(1)),a(itran(2)),
     >         a(itran(3)),a(iwrk(1)),a(iwrk(2)),a(iwrk(3)),
     >         a(iwrk_x),a(iwrk_y),a(ik),a(il),a(isz(is)),nx,ny,nzs,
     >         0,0,0)
c         call write_3d(idcdf1,a(iptr),a(fptr1),a(iter(1)),
c    >           a(isz(is)),a(isclr+2-1),mix,ctmp,nx,ny,nz,nzs,
c    >           nxout,nyout,nzout,iout,jout,is,0,0,0,0)     
        enddo

        call trnsfrm3d(a(igxi),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,
     >                           nzs,0,0)
        call trnsfrm3d(a(iget),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,
     >                           nzs,0,0)
        call trnsfrm3d(a(igze),a(itran(1)),a(iwrk_x),a(iwrk_y),nx,ny,
     >                           nzs,0,0)

        ifld = mptr(7) 
        ift = mptr(13) 
        fptr1 = mptr(9)
        fptr2 = mptr(11)
        fptr3 = 0

        do k=1,3
          if (k.eq.1) then
            idum = mptr(19)
            ctmp = 'xi2'
            inf1 = iuf
            inf2 = igxi
            inf3 = fptr1
            inf4 = fptr2
            inf5 = idum
          else if (k.eq.2) then
            fptr3 = mptr(17)
            ctmp = 'eta2'
            inf1 = ivf
            inf2 = iget
            inf3 = fptr1
            inf4 = fptr2
            inf5 = fptr3
          else if (k.eq.3) then
            ctmp = 'zeta2'
            inf1 = iwf
            inf2 = igze
            inf3 = fptr1
            inf4 = fptr2
            inf5 = fptr3
          end if
          call calc_o2(ctmp,a(ifld),a(ift),a(inf1),a(inf2),a(inf3),
     >                  a(inf4),a(inf5),a(itran(1)),a(iwrk(1)),
     >                  a(iwrk_x),a(iwrk_y),a(ik),a(il),nx,ny,nzs)
          call write_3d(idcdf1,a(ifld),a(fptr1),a(iter(1)),a(isz(is)),
     >             a(isclr+2-1),mix,ctmp,nx,ny,nz,nzs,nxout,nyout,
     >             nzout,iout,jout,is,0,0,0,0)
        enddo

        is = is+1
        goto 107

      end if

      iptr = mptr(1)
      fptr = mptr(3)
      ctmp = 'w1'
      call calc_o1(a(iptr),a(fptr),a(izs_p),a(itran(1)),a(iwrk(1)),
     >                a(iwgz+nzout-1),a(iwrk_x),a(iwrk_y),a(ik),a(il),
     >                dx,dy,ctmp,nx,ny,1,0,0)
      call write_3d(idcdf1,a(iptr),a(fptr),a(iter(1)),a(iwgz+nzout-1),
     >                a(isclr+2-1),mix,ctmp,nx,ny,nz,1,nxout,nyout,
     >                nzout,iout,jout,is,0,0,0,0)

      iwf = mptr(5)
      idwf = mptr(7)
      idgu = mptr(9)
      idgv = mptr(11)
      igb = mptr(13)
      igw = mptr(15)
      ctmp = 'w2'
      do ij=1,2*nx*ny
        a(idgu+ij-1) = 0.
        a(idgv+ij-1) = 0.
        a(igb+ij-1) = 0.
        a(igw+ij-1) = 0.
      enddo
      call calc_w2(a(iptr),a(iwf),a(idwf),a(idgu),a(idgv),a(igb),
     >                  a(igw),a(igh),a(ircof),a(iGs),a(ic1),a(ic2),
     >                  a(ic1w),a(ic2w),a(itran(1)),a(iwrk(1)),
     >                  a(iwrk(2)),a(iwrk_x),a(iwrk_y),a(ik),a(il),
     >                  a(isgzt+nzout-1),a(iwgzt+nzout-1),dz,
     >                  nx,ny,1,is)
      call write_3d(idcdf1,a(iptr),a(fptr),a(iter(1)),a(iwgz+nzout-1),
     >                a(isclr+2-1),mix,ctmp,nx,ny,nz,1,nxout,nyout,
     >                nzout,iout,jout,is,0,0,0,0)

      if (itopo.eq.1) then

        lgth = index(outfl1,' ') - 1

        call set_outfl(idcdf2,nxout,nyout,nzout,dx,dy,dz,xout,yout,
     >                  a(isgz),a(iwgz),a(isclr),1.,sclr,var,units,
     >                  outfl2,nsclr,nvout)

        call write_2d(idcdf2,a(iter(1)),a(itran(1)),'zbot_p',iout,
     >                   jout,nx,ny,nxout,nyout,0,0)
        call write_2d(idcdf2,a(iter(2)),a(itran(1)),'zbot_u',iout,
     >                   jout,nx,ny,nxout,nyout,1,0)
        call write_2d(idcdf2,a(iter(3)),a(itran(1)),'zbot_v',iout,
     >                   jout,nx,ny,nxout,nyout,0,1)

        ifld1 = mptr(1)
        ifld2 = ifld1 + nxout*nyout*nzout
        fptr = ifld2 + nxout*nyout*nzout

        if (fptr.gt.ispc) then
          write(6,*)
          write(6,*) 'Error:  not enough space for interpolations!'
          write(6,*) '        needed:  ',fptr
          write(6,*) '        allocated:  ',ispc
          write(6,*) 'Data on constant geometric height surfaces is ',
     >                       'in ',outfl1(1:lgth),'.'
          write(6,*)
          stop
        end if

        call interp(idcdf1,idcdf2,a,ispc,ifld1,ifld2,iter,isgz,iwgz,
     >                 isclr,nx,ny,nxout,nyout,nzout,iout,jout)

        call ncclos(idcdf2,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem closing data file.'
          stop
        end if

      end if

      call ncclos(idcdf1,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem closing data file.'
        stop
      end if

      if (itopo.eq.1) then
        ctmp(1:len(ctmp)) = ' '
        ctmp = 'rm -f '//outfl1(1:lgth)
        call system(ctmp)
      end if

 45   FORMAT(A80)

      time2 = etime(rtmp)

      write(6,*)
      write(6,*) 'Elapsed Time:  ',time2-time1
      write(6,*)

      end

*----------------------------------------------------------------------

      subroutine rad_coeff(a,ispc,ircof,iu1,iv1,iw1,ib1,igu,igv,igw,
     >                       igb,idgu,idgv,idgb,id11,id12,id13,id22,
     >                       id23,id33,ikm,fptr1,fptr2,fptr3,izs_p,
     >                       izs_u,izs_v,itran,iwrk,iwrk_x,iwrk_y,
     >                       ik,il,isclr,isz,iwz,nslb,ntop,dx,dy,dz,
     >                       nx,ny,nz,icalc)

      integer ispc
      real a(ispc)

      integer ircof,iu1,iv1,iw1,ib1,igu,igv,igw,igb,idgu,idgv,idgb
      integer id11,id12,id13,id22,id23,id33,ikm,fptr1,fptr2,fptr3
      integer izs_p,izs_u,izs_v,itran(3),iwrk(3),iwrk_x,iwrk_y,ik,il
      integer isclr,nslb,isz(nslb),iwz(nslb),ntop,nx,ny,nz,icalc
      real dx,dy,dz

      integer ijk,is,j,k
      integer nzs,iu,iv,iw,igz,idz,iGfac,isave(4)
      real mix,twod
      character*(80) ctmp

      mix = a(isclr+19-1)
      twod = a(isclr+22-1)

      isave(1) = idgu
      isave(2) = idgv
      isave(3) = igw
      isave(4) = igb

      is = 1
 99   if (is.le.nslb) then

        if ((is.eq.nslb).and.(icalc.eq.1)) then
          nzs = ntop 
        else
          nzs = nz
        end if

        if (icalc.eq.1) then

          do ijk=1,nx*ny*nz
            a(igu+ijk-1) = 0.
            a(igv+ijk-1) = 0.
            a(igw+ijk-1) = 0.
            a(igb+ijk-1) = 0.
            a(idgu+ijk-1) = 0.
            a(idgv+ijk-1) = 0.
            a(idgb+ijk-1) = 0.
          enddo

          if (mix.ne.0.) then
 
            do j=1,4
              iu = 0
              iv = 0
              igz = isz(is)
              izs = izs_p
              if (j.eq.1) then
                ctmp = 'u1'
                iptr = iu1
                iu = 1
                izs = izs_u
              else if (j.eq.2) then
                ctmp = 'v1'
                iptr = iv1
                iv = 1
                izs = izs_v
              else if (j.eq.3) then 
                ctmp = 'w1'
                iptr = iw1
                igz = iwz(is)
              else
                ctmp = 'b1'
                iptr = ib1
              end if
              call calc_o1(a(iptr),a(fptr1),a(izs),a(itran(1)),
     >           a(iwrk(1)),a(igz),a(iwrk_x),a(iwrk_y),a(ik),a(il),
     >           dx,dy,ctmp,nx,ny,nzs,iu,iv)
            enddo 
 
            call calc_def(a(id11),a(id12),a(id13),a(id22),a(id23),
     >                      a(id33),a(iu1),a(iv1),a(iw1),a(isz(is)),        
     >                      a(iwz(is)),dx,dy,nx,ny,nzs)
            call calc_km(a(id11),a(id12),a(id13),a(id22),a(id23),
     >                      a(id33),a(ib1),a(ikm),a(isz(is)),
     >                      a(isclr+5-1),a(isclr+20-1),a(isclr+21-1),
     >                      nx,ny,nzs)
            call u_mix(a(id11),a(id12),a(id13),a(igu),a(idgu),
     >                      a(iwz(is)),a(isz(is)),dx,dy,nx,ny,nzs)
            call v_mix(a(id12),a(id22),a(id23),a(igv),a(idgv),
     >                      a(iwz(is)),a(isz(is)),dx,dy,nx,ny,nzs)
            call w_mix(a(id13),a(id23),a(id33),a(igw),a(isz(is)),
     >                              dx,dy,nx,ny,nzs)
            call b_mix(a(ikm),a(ib1),a(igb),a(idgb),a(isz(is)),
     >                    a(iwz(is)),a(isclr+21-1),dx,dy,nx,ny,nzs)

            nzs = nzs-6
            idgu = idgu + 3*(nx-1)*(ny-1)
            idgv = idgv + 3*(nx-1)*(ny-1)
            igw = igw + 3*(nx-1)*(ny-1)
            igb = igb + 3*(nx-1)*(ny-1)
            isz(is) = isz(is) + 3
            iwz(is) = iwz(is) + 3
            if (twod.ne.0.) then
              call set_2d(a(idgu),a(idgv),a(idgb),a(igw),a(igb),
     >                      a(igu),a(igv),a(ikm),a(iwrk(1)),
     >                      a(iwrk(2)),nx,ny,nzs)
            end if

          end if

          do k=1,3
            if (k.eq.1) then
              iptr = iu1
              ctmp = 'u1'
            else if (k.eq.2) then
              iptr = iv1
              ctmp = 'v1'
            else
              iptr = iw1
              ctmp = 'w1'
            end if
            call calc_o1(a(iptr),a(fptr1),a(izs_p),a(itran(1)),
     >               a(iwrk(1)),a(isz(is)),a(iwrk_x),a(iwrk_y),a(ik),
     >               a(il),dx,dy,ctmp,nx,ny,nzs,0,0)
          enddo

          do k=1,4
            if (k.eq.1) then
              iptr = idgu
              ctmp = 'dgu'
              idz = 1
            else if (k.eq.2) then
              iptr = idgv
              ctmp = 'dgv'
              idz = 1
            else if (k.eq.3) then
              iptr = igw
              ctmp = 'gw'
              idz = 0
            else
              iptr = igb
              ctmp = 'gb'
              idz = 0
            end if
            call calc_advc(ctmp,a(iptr),a(iu1),a(iv1),a(iw1),a(fptr1),
     >                 a(fptr2),a(fptr3),a(izs_p),a(itran(1)),
     >                 a(itran(2)),a(itran(3)),a(iwrk(1)),a(iwrk(2)),
     >                 a(iwrk(3)),a(iwrk_x),a(iwrk_y),a(ik),a(il),
     >                 a(isz(is)),nx,ny,nzs,0,0,idz) 
          enddo

          call trnsfrm3d(a(igw),a(itran(1)),a(iwrk_x),a(iwrk_y),
     >                          nx,ny,nzs,0,0)
          call trnsfrm3d(a(igb),a(itran(1)),a(iwrk_x),a(iwrk_y),
     >                          nx,ny,nzs,0,0)
          call trnsfrm3d(a(idgu),a(itran(1)),a(iwrk_x),a(iwrk_y),
     >                          nx,ny,nzs,0,0)
          call trnsfrm3d(a(idgv),a(itran(1)),a(iwrk_x),a(iwrk_y),
     >                          nx,ny,nzs,0,0)

        end if

        iGfac = iwrk(1) 

        call calc_rad_coeff(a(ircof),a(idgu),a(idgv),a(igw),a(igb),
     >                          a(iGfac),a(ik),a(il),a(isz(is)),
     >                          dz,nx,ny,nzs)

        if (mix.ne.0.) then
          idgu = isave(1)
          idgv = isave(2)
          igw = isave(3)
          igb = isave(4)
        end if

        is = is+1
        goto 99

      end if

      return
      end

*---------------------------------------------------------------------

      subroutine calc_rad_coeff(rcof,dgu,dgv,gw,gb,Gfac,ak,al,
     >                                sgz,dz,nx,ny,nz)

      integer nx,ny,nz
      real ak(nx-1),al(ny-1),sgz(nz),dz
      complex rcof(nx-1,ny-1),dgu(nx-1,ny-1,nz),dgv(nx-1,ny-1,nz)
      complex gw(nx-1,ny-1,nz),gb(nx-1,ny-1,nz),Gfac(nx-1,ny-1)
        
      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      integer i,j,k
      real omga,m2,gam2,m,gam
      complex exfac
      logical mzero2,requ

      do k=1,nz

        tau = sgz(k)

        call calc_Gfac(Gfac,dgu(1,1,k),dgv(1,1,k),gw(1,1,k),gb(1,1,k),
     >                             ak,al,nx,ny)

        do j=1,ny-1
        do i=1,nx-1

          omga = ak(i)*U0+al(j)*V0
          if (requ(abs(omga),abs(fcor))) then
            m2 = 0.
            gam2 = 0.
          else 
            m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >                 / (omga**2-fcor**2)
            gam2 = -m2
          end if

          if (mzero2(omga,N,fcor,m2)) then

            rcof(i,j) = 0.

          else if (m2.gt.0.) then

            m = sqrt(m2)*omga/abs(omga)
            exfac = cmplx(cos(m*tau),sin(m*tau))
            rcof(i,j) = rcof(i,j) + cmplx(0.,0.5/m)*Gfac(i,j)
     >                          * exfac*dz

          else 

            gam = sqrt(gam2)
            exfac = cmplx(exp(-gam*tau),0.)
            rcof(i,j) = rcof(i,j) + cmplx(0.5/gam,0.)*Gfac(i,j)
     >                      * exfac*dz

          end if

        enddo
        enddo

      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine rad_bc(rcof,gu,gv,ak,al,zt,nx,ny)

      integer nx,ny
      real ak(nx-1),al(ny-1),zt
      complex rcof(nx-1,ny-1),gu(nx-1,ny-1),gv(nx-1,ny-1)

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      integer i,j
      real omga,m2,gam2,m,gam,fac(4)
      logical mzero2,requ
      complex exfac

      do j=1,ny-1
      do i=1,nx-1

        omga = ak(i)*U0+al(j)*V0
        if (requ(abs(omga),abs(fcor))) then
          m2 = 0.
          gam2 = 0.
        else
          m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >               / (omga**2-fcor**2)
          gam2 = -m2
        end if

        if (mzero2(omga,N,fcor,m2)) then

          rcof(i,j) = 0.

        else if (m2.gt.0.) then

          m = sqrt(m2)*omga/abs(omga)
          exfac = cmplx(cos(m*zt),sin(m*zt))
          fac(1) = 0.5*fcor*al(j)/m/(omga**2-fcor**2)
          fac(2) = -0.5*omga*ak(i)/m/(omga**2-fcor**2)
          fac(3) = -0.5*fcor*ak(i)/m/(omga**2-fcor**2)
          fac(4) = -0.5*omga*al(j)/m/(omga**2-fcor**2)
          rcof(i,j) = rcof(i,j) + cmplx(fac(1),fac(2))*gu(i,j)*exfac
     >                   + cmplx(fac(3),fac(4))*gv(i,j)*exfac

        else

          gam = sqrt(gam2)
          exfac = cmplx(exp(-gam*zt),0.)
          fac(1) = -0.5*omga*ak(i)/gam/(omga**2-fcor**2)
          fac(2) = -0.5*fcor*al(j)/gam/(omga**2-fcor**2)
          fac(3) = -0.5*omga*al(j)/gam/(omga**2-fcor**2)
          fac(4) = 0.5*fcor*ak(i)/gam/(omga**2-fcor**2)
          rcof(i,j) = rcof(i,j) + cmplx(fac(1),fac(2))*gu(i,j)*exfac
     >                   + cmplx(fac(3),fac(4))*gv(i,j)*exfac

        end if

      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine set_2d(dgu,dgv,dgb,gw,gb,gu,gv,km,gut,gvt,nx,ny,nz)

      integer nx,ny,nz
      real dgu(nx-1,ny-1,nz),dgv(nx-1,ny-1,nz),dgb(nx-1,ny-1,nz)
      real gw(nx-1,ny-1,nz),gb(nx-1,ny-1,nz),gu(nx-1,ny-1,nz)
      real gv(nx-1,ny-1,nz),km(nx-1,ny-1,nz)
      real gut(nx-1,ny-1),gvt(nx-1,ny-1)

      integer i,j,k

      do j=1,ny-1
      do i=1,nx-1

        gut(i,j) = gut(i,3)
        gvt(i,j) = gvt(i,3)
        do k=1,nz
          dgu(i,j,k) = dgu(i,3,k)
          dgv(i,j,k) = dgv(i,3,k)
          dgb(i,j,k) = dgb(i,3,k)
          gw(i,j,k) = gw(i,3,k)
          gb(i,j,k) = gb(i,3,k)
          gu(i,j,k) = gu(i,3,k)
          gv(i,j,k) = gv(i,3,k)
          km(i,j,k) = km(i,3,k)
        enddo

      enddo
      enddo

      return
      end
 
*----------------------------------------------------------------------

      subroutine mix_zt(gut,gvt,gu,gv,nx,ny,ntop)

      integer nx,ny,ntop
      real gut(nx-1,ny-1),gvt(nx-1,ny-1)
      real gu(nx-1,ny-1,ntop),gv(nx-1,ny-1,ntop)

      integer i,j

      do j=1,ny-1
      do i=1,nx-1
        gut(i,j) = gut(i,j) + 0.5*(gu(i,j,ntop-2)+gu(i,j,ntop-3))
        gvt(i,j) = gvt(i,j) + 0.5*(gv(i,j,ntop-2)+gv(i,j,ntop-3))
      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine calc_wavnm(ak,ds,n)

      integer n
      real ak(n),ds

      integer q
      real pi,L

      pi = 2.*asin(1.)
      L = n*ds
      q = n/2

      do j=1,n
        if (j.le.q+1) then
          ak(j) = 2.*pi*(j-1)/L
        else
          ak(j) = 2.*pi*(j-1-n)/L
        end if
      enddo

      return
      end

*-----------------------------------------------------------------------

      subroutine trnsfrm2d(fld,tran,cwrkx,cwrky,nx,ny,iu,iv)

      integer nx,ny,iu,iv
      real cwrkx(4*nx+15),cwrky(4*ny+15)
      complex fld(nx-1+iu,ny-1+iv),tran(ny-1,nx-1)

      call r_to_c(fld,nx,ny,iu,iv)

      do j=1,ny-1
        call cfftf(nx-1,fld(1,j),cwrkx)
        do i=1,nx-1
          tran(j,i) = fld(i,j)
        enddo
      enddo

      do i=1,nx-1
        call cfftf(ny-1,tran(1,i),cwrky)
        do j=1,ny-1
          fld(i,j) = tran(j,i)/(nx-1)/(ny-1)
        enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine trnsfrm3d(fld,tran,cwrkx,cwrky,nx,ny,nz,iu,iv)

      integer nx,ny,nz,iu,iv
      real cwrkx(4*nx+15),cwrky(4*ny+15)
      complex fld(nx-1+iu,ny-1+iv,nz),tran(ny-1,nx-1)

      integer k

      call space_out(fld,nx-1+iu,ny-1+iv,nz)

      do k=1,nz
        call trnsfrm2d(fld(1,1,k),tran,cwrkx,cwrky,nx,ny,iu,iv)
      enddo

      return
      end

*--------------------------------------------------------------------------

      subroutine space_out(fld,n1,n2,n3)

      integer n1,n2,n3
      real fld(n1,n2,2*n3)

      integer i,j,k,q

      do k=n3,2,-1

        q = 2*k-1

        do i=1,n1
        do j=1,n2
          fld(i,j,q) = fld(i,j,k)
        enddo
        enddo

      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine r_to_c(fld,nx,ny,iu,iv)

      integer nx,ny,iu,iv
      real fld(nx-1+iu,ny-1+iv,2)

      integer ptr,q,numx,numy

      numx = nx-1+iu
      numy = ny-1+iv

      do j=1,numy
      do i=1,numx
        fld(i,j,2) = fld(i,j,1)
      enddo
      enddo

      ptr = numx*numy + 1
 
      do k=1,2
      do j=1,numy
      do i=1,numx

        q = (k-1)*numx*numy + (j-1)*numx + i
        if (mod(q,2).ne.0) then
          fld(i,j,k) = fld(ptr,1,1)
          ptr = ptr + 1
        else
          fld(i,j,k) = 0.
        end if

      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine calc_Gfac(Gfac,dgu,dgv,gw,gb,ak,al,nx,ny)      

      integer nx,ny
      real ak(nx-1),al(ny-1)
      complex Gfac(nx-1,ny-1),dgu(nx-1,ny-1),dgv(nx-1,ny-1)
      complex gw(nx-1,ny-1),gb(nx-1,ny-1)

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      integer i,j
      real omga,m2
      logical mzero2,requ
      complex fac(4)

      do j=1,ny-1
      do i=1,nx-1

        omga = ak(i)*U0+al(j)*V0

        if (requ(abs(omga),abs(fcor))) then
          m2 = 0.
        else
          m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >             / (omga**2-fcor**2)
        end if        

        if (mzero2(omga,N,fcor,m2)) then
          Gfac(i,j) = 0.
        else
          fac(1) = cmplx(omga*ak(i),fcor*al(j))
          fac(2) = cmplx(omga*al(j),-fcor*ak(i))
          fac(3) = cmplx(-ak(i)**2-al(j)**2,0.)
          fac(4) = cmplx(0.,-omga*(ak(i)**2+al(j)**2))
          Gfac(i,j) = (fac(1)*dgu(i,j) + fac(2)*dgv(i,j)
     >                + fac(3)*gb(i,j) + fac(4)*gw(i,j))
     >                   /(omga**2-fcor**2)
        end if

      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine calc_o1(fld,ft,hft,tran,wrk,gz,cwrkx,cwrky,ak,al,
     >                       dx,dy,var,nx,ny,nz,iu,iv)

      integer nx,ny,nz,iu,iv
      real fld(nx-1+iu,ny-1+iv,nz)
      real ak(nx-1),al(ny-1),gz(nz),cwrkx(4*nx+15),cwrky(4*ny+15)
      real dx,dy
      complex ft(nx-1,ny-1,nz)
      complex hft(nx-1+iu,ny-1+iv),tran(ny-1,nx-1),wrk(nx-1,ny-1)
      character*(*) var

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      real x,y,z,m,m2,gam,gam2,omga
      logical mzero1,mzero2,requ
      complex pfac,exfac

      lgth = index(var,' ') - 1

      if (var(lgth:lgth).eq.'0') then
 
        if (var(1:lgth).eq.'u0') then
          do k=1,nz  
          do j=1,ny-1
          do i=1,nx
            fld(i,j,k) = U0
          enddo
          enddo
          enddo
        else if (var(1:lgth).eq.'v0') then
          do k=1,nz
          do j=1,ny
          do i=1,nx-1
            fld(i,j,k) = V0
          enddo
          enddo
          enddo
        else if (var(1:lgth).eq.'b0') then
          do k=1,nz
          do j=1,ny-1
          do i=1,nx-1
            z = gz(k)
            fld(i,j,k) = -g*(1.-N**2/g*z)        
          enddo
          enddo
          enddo
        else if (var(1:lgth).eq.'p0') then
          do k=1,nz
          do j=1,ny-1
          do i=1,nx-1
            z = gz(k)
            x = ((i-1)+0.5)*dx
            y = ((j-1)+0.5)*dy
            fld(i,j,k) = ps/rhos + fcor*V0*x - fcor*U0*y
     >                   - g*z + 0.5*(N**2)*(z**2)
          enddo
          enddo
          enddo
        else 
          write(6,*) 'Error:  field ',var(1:lgth),' not recognized.'
          stop
        end if

        return

      end if

      do j=1,ny-1
      do i=1,nx-1

        omga = ak(i)*U0+al(j)*V0 
        if (requ(abs(omga),abs(fcor))) then
          m2 = 0.
          gam2 = 0.
        else
          m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >               / (omga**2-fcor**2)
          gam2 = -m2
        end if
  
        if (mzero1(omga,N,fcor,m2)) then

          do k=1,nz
            ft(i,j,k) = 0.
          enddo

        else if (m2.gt.0.) then

          m = sqrt(m2)*omga/abs(omga)
          call calc_pfac(pfac,ak,al,m2,m,gam2,0.,N,omga,fcor,
     >                            var,i,j,nx,ny,0)
          do k=1,nz
            z = gz(k)
            exfac = cmplx(cos(m*z),sin(m*z))
            ft(i,j,k) = pfac*exfac*hft(i,j)
          enddo

        else

          gam = sqrt(gam2)
          call calc_pfac(pfac,ak,al,m2,0.,gam2,gam,N,omga,fcor,
     >                              var,i,j,nx,ny,0)
          do k=1,nz
            z = gz(k)
            exfac = cmplx(exp(-gam*z),0.)
            ft(i,j,k) = pfac*exfac*hft(i,j)
          enddo

        end if

      enddo
      enddo

      do k=1,nz

        do i=1,nx-1
          do j=1,ny-1
            tran(j,i) = ft(i,j,k)
          enddo
          call cfftb(ny-1,tran(1,i),cwrky) 
          do j=1,ny-1
            wrk(i,j) = tran(j,i)
          enddo
        enddo

        do j=1,ny-1
          call cfftb(nx-1,wrk(1,j),cwrkx)
          do i=1,nx-1
            fld(i,j,k) = real(wrk(i,j))
          enddo
        enddo

        if (iu.eq.1) then
          do j=1,ny-1
            fld(nx,j,k) = fld(1,j,k)
          enddo
        else if (iv.eq.1) then
          do i=1,nx-1
            fld(i,ny,k) = fld(i,1,k)
          enddo
        end if

      enddo

      return
      end

*---------------------------------------------------------------------

      subroutine calc_pfac(pfac,ak,al,m2,m,gam2,gam,N,omga,fcor,var, 
     >                                  i,j,nx,ny,idz)

      integer i,j,nx,ny,idz
      real ak(nx-1),al(ny-1),m2,m,gam2,gam,N,omga,fcor
      character*(*) var
      complex pfac

      integer lgth
      real fac1,fac2
      
      lgth = index(var,' ') - 1

      if (m2.gt.0.) then
      
        if (var(1:1).eq.'w') then
          if (idz.eq.0) then
            pfac = cmplx(0.,omga)
          else
            pfac = cmplx(-m*omga,0.)
          end if
        else if (var(1:1).eq.'u') then
          if (idz.eq.0) then
            fac1 = -m*fcor*al(j)/(ak(i)**2+al(j)**2) 
            fac2 = -m*ak(i)*omga/(ak(i)**2+al(j)**2) 
          else
            fac1 = m2*ak(i)*omga/(ak(i)**2+al(j)**2)
            fac2 = -m2*fcor*al(j)/(ak(i)**2+al(j)**2)
          end if
          pfac = cmplx(fac1,fac2)
        else if (var(1:1).eq.'v') then
          if (idz.eq.0) then
            fac1 = m*fcor*ak(i)/(ak(i)**2+al(j)**2)
            fac2 = -m*al(j)*omga/(ak(i)**2+al(j)**2)
          else
            fac1 = m2*al(j)*omga/(ak(i)**2+al(j)**2)
            fac2 = m2*fcor*ak(i)/(ak(i)**2+al(j)**2)
          end if     
          pfac = cmplx(fac1,fac2)
        else if (var(1:1).eq.'p') then
          if (idz.eq.0) then
            fac1 = 0.
            fac2 = (N**2-omga**2)/m
          else
            fac1 = -(N**2-omga**2)
            fac2 = 0.
          end if
          pfac = cmplx(fac1,fac2)
        else if (var(1:1).eq.'b') then
          if (idz.eq.0) then
            pfac = cmplx(-N**2,0.)
          else
            pfac = cmplx(0.,-m*N**2)
          end if
        else if (var(1:1).eq.'q') then
          fac1 = -al(j)*omga*(ak(i)**2+al(j)**2+m2)
     >                 / (ak(i)**2+al(j)**2)
          fac2 = -fcor*ak(i)*m2/(ak(i)**2+al(j)**2)
          pfac = cmplx(fac1,fac2) 
        else if (var(1:1).eq.'r') then
          fac1 = ak(i)*omga*(ak(i)**2+al(j)**2+m2)
     >                 / (ak(i)**2+al(j)**2)
          fac2 = -fcor*al(j)*m2/(ak(i)**2+al(j)**2)
          pfac = cmplx(fac1,fac2)
        else if (var(1:1).eq.'s') then
          pfac = cmplx(0.,fcor*m)
        else if (var(1:2).eq.'xi') then
          fac1 = -m*ak(i)/(ak(i)**2+al(j)**2)
          fac2 = m*fcor*al(j)/omga/(ak(i)**2+al(j)**2)
          pfac = cmplx(fac1,fac2)
        else if (var(1:2).eq.'et') then
          fac1 = -m*al(j)/(ak(i)**2+al(j)**2)
          fac2 = -m*fcor*ak(i)/omga/(ak(i)**2+al(j)**2)
          pfac = cmplx(fac1,fac2)
        else if (var(1:2).eq.'ze') then
          pfac = cmplx(1.,0.)
        else 
          write(6,*) 'Error:  field ',var(1:lgth),
     >                      ' not recognized.'
          stop
        end if 

      else

        if (var(1:1).eq.'w') then
          if (idz.eq.0) then
            pfac = cmplx(0.,omga)
          else
            pfac = cmplx(0.,-gam*omga)
          end if
        else if (var(1:1).eq.'u') then
          if (idz.eq.0) then
            fac1 = gam*ak(i)*omga/(ak(i)**2+al(j)**2)
            fac2 = -gam*fcor*al(j)/(ak(i)**2+al(j)**2)
          else
            fac1 = -gam2*ak(i)*omga/(ak(i)**2+al(j)**2)
            fac2 = gam2*fcor*al(j)/(ak(i)**2+al(j)**2)
          end if
          pfac = cmplx(fac1,fac2)
        else if (var(1:1).eq.'v') then
          if (idz.eq.0) then
            fac1 = gam*al(j)*omga/(ak(i)**2+al(j)**2)
            fac2 = gam*fcor*ak(i)/(ak(i)**2+al(j)**2)
          else
            fac1 = -gam2*al(j)*omga/(ak(i)**2+al(j)**2)
            fac2 = -gam2*fcor*ak(i)/(ak(i)**2+al(j)**2)
          end if
          pfac = cmplx(fac1,fac2)
        else if (var(1:1).eq.'p') then
          if (idz.eq.0) then
            fac1 = (N**2-omga**2)/gam
            fac2 = 0.
          else
            fac1 = -(N**2-omga**2)
            fac2 = 0.
          end if
          pfac = cmplx(fac1,fac2)
        else if (var(1:1).eq.'b') then
          if (idz.eq.0) then
            pfac = cmplx(-N**2,0.)
          else
            pfac = cmplx(gam*N**2,0.)
          end if
        else if (var(1:1).eq.'q') then
          fac1 = -al(j)*omga*(ak(i)**2+al(j)**2-gam2)
     >                   / (ak(i)**2+al(j)**2)
          fac2 = fcor*ak(i)*gam2/(ak(i)**2+al(j)**2)
          pfac = cmplx(fac1,fac2)
        else if (var(1:1).eq.'r') then
          fac1 = ak(i)*omga*(ak(i)**2+al(j)**2-gam2)
     >                   / (ak(i)**2+al(j)**2)
          fac2 = fcor*al(j)*gam2/(ak(i)**2+al(j)**2)
          pfac = cmplx(fac1,fac2)
        else if (var(1:1).eq.'s') then
          pfac = cmplx(-fcor*gam,0.)
        else if (var(1:2).eq.'xi') then
          fac1 = -gam*fcor*al(j)/omga/(ak(i)**2+al(j)**2)
          fac2 = -gam*ak(i)/(ak(i)**2+al(j)**2)
          pfac = cmplx(fac1,fac2)
        else if (var(1:2).eq.'et') then
          fac1 = gam*fcor*ak(i)/omga/(ak(i)**2+al(j)**2)
          fac2 = -gam*al(j)/(ak(i)**2+al(j)**2)
          pfac = cmplx(fac1,fac2)
        else if (var(1:2).eq.'ze') then
          pfac = cmplx(1.,0.)
        else 
          write(6,*) 'Error:  field ',var(1:lgth),' not recognized.'
          stop
        end if 

      end if

      return
      end

*---------------------------------------------------------------------

      subroutine calc_gh(gh,us,vs,zs_p,hft,tran1,tran2,tran3,wrk1,
     >                     wrk2,wrk3,cwrkx,cwrky,ak,al,nx,ny)

      integer nx,ny
      real gh(nx-1,ny-1),us(nx-1,ny-1),vs(nx-1,ny-1)
      real zs_p(nx-1,ny-1)
      real ak(nx-1),al(ny-1),cwrkx(4*nx+15),cwrky(4*ny+15)
      complex hft(nx-1,ny-1)
      complex tran1(ny-1,nx-1),tran2(ny-1,nx-1),tran3(ny-1,nx-1)
      complex wrk1(nx-1,ny-1),wrk2(nx-1,ny-1),wrk3(nx-1,ny-1)

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      integer i,j
      real dwdz,dhdx,dhdy
      real m,m2,gam,gam2,omga,fac1,fac2
      character*1 var
      logical mzero2,requ
      complex pfac

      var(1:1) = 'w'

      do i=1,nx-1

        do j=1,ny-1

          omga = ak(i)*U0+al(j)*V0
          if (requ(abs(omga),abs(fcor))) then
            m2 = 0.
            gam2 = 0.
          else
            m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >                 / (omga**2-fcor**2)
            gam2 = -m2
          end if

          if (mzero2(omga,N,fcor,m2)) then
            pfac = 0.
          else if (m2.gt.0.) then
            m = sqrt(m2)*omga/abs(omga)
            call calc_pfac(pfac,ak,al,m2,m,gam2,0.,N,omga,fcor,var,
     >                               i,j,nx,ny,1) 
          else 
            gam = sqrt(gam2)
            call calc_pfac(pfac,ak,al,m2,0.,gam2,gam,N,omga,fcor,var,
     >                               i,j,nx,ny,1)
          end if

          tran1(j,i) = cmplx(0.,ak(i))*hft(i,j)
          tran2(j,i) = cmplx(0.,al(j))*hft(i,j)
          tran3(j,i) = pfac*hft(i,j)

        enddo

        call cfftb(ny-1,tran1(1,i),cwrky)
        call cfftb(ny-1,tran2(1,i),cwrky)
        call cfftb(ny-1,tran3(1,i),cwrky)
        do j=1,ny-1
          wrk1(i,j) = tran1(j,i)
          wrk2(i,j) = tran2(j,i)
          wrk3(i,j) = tran3(j,i)
        enddo

      enddo

      do j=1,ny-1
        call cfftb(nx-1,wrk1(1,j),cwrkx)
        call cfftb(nx-1,wrk2(1,j),cwrkx)
        call cfftb(nx-1,wrk3(1,j),cwrkx)
        do i=1,nx-1
          dhdx = real(wrk1(i,j))
          dhdy = real(wrk2(i,j))
          dwdz = real(wrk3(i,j))
          gh(i,j) = us(i,j)*dhdx + vs(i,j)*dhdy - dwdz*zs_p(i,j)
        enddo
      enddo

      return
      end

*---------------------------------------------------------------------

      subroutine calc_gs(gs,q1,r1,s1,q2,r2,s2,wf,gb,dwf,dgb,ft,hft,
     >                     tran,wrk,cwrkx,cwrky,ak,al,gz,nx,ny,nz,
     >                     iflg)

      integer nx,ny,nz,iflg
      real gs(nx-1,ny-1,nz),q1(nx-1,ny-1,nz),r1(nx-1,ny-1,nz)
      real s1(nx-1,ny-1,nz),q2(nx-1,ny-1,nz),r2(nx-1,ny-1,nz)
      real s2(nx-1,ny-1,nz),cwrkx(4*nx+15),cwrky(4*ny+15)
      real ak(nx-1),al(ny-1),gz(nz)
      complex wf(nx-1,ny-1,nz),gb(nx-1,ny-1,nz)
      complex dwf(nx-1,ny-1,nz),dgb(nx-1,ny-1,nz) 
      complex ft(nx-1,ny-1,nz),hft(nx-1,ny-1)
      complex tran(ny-1,nx-1),wrk(nx-1,ny-1)

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g
 
      integer i,j,k,trm,idz
      real m,m2,gam,gam2,omga,fac1,fac2,fac3,fac4,z
      real dbdx,dbdy,dbdz,q,r,s
      character*(8) ctmp
      logical mzero2,requ
      complex pfac,exfac

      if ((iflg.ne.2).and.(iflg.ne.3)) then
        write(6,*) 'Error:  do not know which field to calculate ',
     >                         'in calc_gs.'
        write(6,*) 'Stopping run.'
        stop
      end if

      ctmp = 'b1'
      idz = 0

      do trm=1,3

        if (trm.eq.3) idz = 1 

        do j=1,ny-1       
        do i=1,nx-1

          omga = ak(i)*U0+al(j)*V0
          if (requ(abs(omga),abs(fcor))) then
            m2 = 0.
            gam2 = 0.
          else
            m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >                  / (omga**2-fcor**2)
            gam2 = -m2
          end if

          if (mzero2(omga,N,fcor,m2)) then

            do k=1,nz
              ft(i,j,k) = 0.
            enddo
        
          else if (m2.gt.0.) then
     
            m = sqrt(m2)*omga/abs(omga)
            call calc_pfac(pfac,ak,al,m2,m,gam2,0.,N,omga,fcor,ctmp,
     >                              i,j,nx,ny,idz)
            if (trm.eq.1) then
              pfac = cmplx(0.,ak(i))*pfac
            else if (trm.eq.2) then
              pfac = cmplx(0.,al(j))*pfac
            end if

            do k=1,nz
              z = gz(k)
              exfac = cmplx(cos(m*z),sin(m*z))
              ft(i,j,k) = pfac*exfac*hft(i,j)
            enddo

          else

            gam = sqrt(gam2)
            call calc_pfac(pfac,ak,al,m2,0.,gam2,gam,N,omga,fcor,
     >                             ctmp,i,j,nx,ny,idz)
            if (trm.eq.1) then
              pfac = cmplx(0.,ak(i))*pfac
            else if (trm.eq.2) then
              pfac = cmplx(0.,al(j))*pfac
            end if

            do k=1,nz
              z = gz(k)
              exfac = cmplx(exp(-gam*z),0.)
              ft(i,j,k) = pfac*exfac*hft(i,j)
            enddo

          end if

        enddo
        enddo

        do k=1,nz

          do i=1,nx-1
            do j=1,ny-1
              tran(j,i) = ft(i,j,k)
            enddo
            call cfftb(ny-1,tran(1,i),cwrky)
            do j=1,ny-1
              wrk(i,j) = tran(j,i)
            enddo
          enddo

          if (trm.eq.1) then

            do j=1,ny-1
              call cfftb(nx-1,wrk(1,j),cwrkx)
              do i=1,nx-1
                dbdx = real(wrk(i,j))
                if (iflg.eq.2) then
                  q = q1(i,j,k)
                else 
                  q = q2(i,j,k)
                end if
                gs(i,j,k) = gs(i,j,k) - q*dbdx
              enddo
            enddo

          else if (trm.eq.2) then 

            do j=1,ny-1
              call cfftb(nx-1,wrk(1,j),cwrkx)
              do i=1,nx-1
                dbdy = real(wrk(i,j))
                if (iflg.eq.2) then
                  r = r1(i,j,k)
                else 
                  r = r2(i,j,k)
                end if
                gs(i,j,k) = gs(i,j,k) - r*dbdy
              enddo
            enddo

          else
                  
            do j=1,ny-1
              call cfftb(nx-1,wrk(1,j),cwrkx)
              do i=1,nx-1
                dbdz = real(wrk(i,j))
                if (iflg.eq.2) then
                  s = s1(i,j,k)
                else 
                  s = s2(i,j,k)
                end if
                gs(i,j,k) = gs(i,j,k) - s*dbdz
              enddo
            enddo

          end if

        enddo

      enddo

      if (iflg.eq.3) then

        do trm = 1,3

          do j=1,ny-1
          do i=1,nx-1

            omga = ak(i)*U0+al(j)*V0
            if (requ(abs(omga),abs(fcor))) then
              m2 = 0.
              gam2 = 0.
            else
              m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >                    / (omga**2-fcor**2)
              gam2 = -m2
            end if

            if (mzero2(omga,N,fcor,m2)) then
            
              do k=1,nz
                ft(i,j,k) = 0.
              enddo

            else

              fac1 = 0.
              fac2 = N**2/omga
              fac3 = 0.
              fac4 = -1./omga 

              if (trm.eq.3) then
                do k=1,nz
                  ft(i,j,k) = cmplx(fac1,fac2)*dwf(i,j,k)
     >                            + cmplx(fac3,fac4)*dgb(i,j,k)
                enddo
              else
                do k=1,nz
                  ft(i,j,k) = cmplx(fac1,fac2)*wf(i,j,k)
     >                            + cmplx(fac3,fac4)*gb(i,j,k)
                  if (trm.eq.1) then
                    ft(i,j,k) = cmplx(0.,ak(i))*ft(i,j,k)
                  else
                    ft(i,j,k) = cmplx(0.,al(j))*ft(i,j,k)
                  end if
                enddo
              end if

            end if

          enddo
          enddo

          do k=1,nz

            do i=1,nx-1
              do j=1,ny-1
                tran(j,i) = ft(i,j,k)
              enddo
              call cfftb(ny-1,tran(1,i),cwrky)
              do j=1,ny-1
                wrk(i,j) = tran(j,i)
              enddo
            enddo

            if (trm.eq.1) then

              do j=1,ny-1
                call cfftb(nx-1,wrk(1,j),cwrkx)
                do i=1,nx-1
                  dbdx = real(wrk(i,j))
                  q = q1(i,j,k)
                  gs(i,j,k) = gs(i,j,k) - q*dbdx
                enddo
              enddo

            else if (trm.eq.2) then 

              do j=1,ny-1
                call cfftb(nx-1,wrk(1,j),cwrkx)
                do i=1,nx-1
                  dbdy = real(wrk(i,j))
                  r = r1(i,j,k)
                  gs(i,j,k) = gs(i,j,k) - r*dbdy
                enddo
              enddo

            else

              do j=1,ny-1
                call cfftb(nx-1,wrk(1,j),cwrkx)
                do i=1,nx-1
                  dbdz = real(wrk(i,j))
                  s = s1(i,j,k)
                  gs(i,j,k) = gs(i,j,k) - s*dbdz
                enddo
              enddo

            end if

          enddo

        enddo

      end if

      return
      end

*---------------------------------------------------------------------

      subroutine calc_advc(var,gs,u,v,w,ft1,ft2,ft3,hft,tran1,tran2,
     >                     tran3,wrk1,wrk2,wrk3,cwrkx,cwrky,ak,al,gz,
     >                     nx,ny,nz,iu,iv,idz)

      integer nx,ny,nz,iu,iv,idz
      real gs(nx-1+iu,ny-1+iv,nz),u(nx-1+iu,ny-1+iv,nz)
      real v(nx-1+iu,ny-1+iv,nz),w(nx-1+iu,ny-1+iv,nz)
      real cwrkx(4*nx+15),cwrky(4*ny+15)
      real ak(nx-1),al(ny-1),gz(nz)
      complex ft1(nx-1,ny-1,nz),ft2(nx-1,ny-1,nz)
      complex ft3(nx-1,ny-1,nz),hft(nx-1+iu,ny-1+iv)
      complex tran1(ny-1,nx-1),tran2(ny-1,nx-1),tran3(ny-1,nx-1)
      complex wrk1(nx-1,ny-1),wrk2(nx-1,ny-1),wrk3(nx-1,ny-1)
      character*(*) var

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      integer i,j,k,trm
      real dsdx,dsdy,dsdz,du,ds
      real m,m2,gam,gam2,omga,fac1,fac2,z
      character*(8) vtmp,ctmp
      logical mzero2,requ
      complex pfac,exfac,pfac1,pfac2

      lgth = index(var,' ') - 1
      if (lgth.eq.0) lgth = len(var)
      if (var(1:1).eq.'g') then
        vtmp(1:lgth-1) = var(2:lgth)
      else if (var(1:2).eq.'dg') then
        vtmp(1:lgth-2) = var(3:lgth)
      else
        write(6,*) 'Error:  do not know what to calculate in ',
     >                               'calc_advc.'
        write(6,*) 'Stopping run.'
        stop
      end if
        
      do j=1,ny-1
      do i=1,nx-1
 
        omga = ak(i)*U0+al(j)*V0
        if (requ(abs(omga),abs(fcor))) then
          m2 = 0.
          gam2 = 0.
        else
          m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >               / (omga**2-fcor**2)
          gam2 = -m2
        end if

        if (mzero2(omga,N,fcor,m2)) then

          do k=1,nz
            ft1(i,j,k) = 0.
            ft2(i,j,k) = 0.
            ft3(i,j,k) = 0.
          enddo 
         
        else if (m2.gt.0.) then

          m = sqrt(m2)*omga/abs(omga)
          call calc_pfac(pfac,ak,al,m2,m,gam2,0.,N,omga,fcor,vtmp,
     >                             i,j,nx,ny,idz)
          do k=1,nz
            z = gz(k)
            exfac = cmplx(cos(m*z),sin(m*z))
            ft1(i,j,k) = cmplx(0.,ak(i))*pfac*exfac*hft(i,j)
            ft2(i,j,k) = cmplx(0.,al(j))*pfac*exfac*hft(i,j)
            ft3(i,j,k) = cmplx(0.,m)*pfac*exfac*hft(i,j)
          enddo

        else

          gam = sqrt(gam2)
          call calc_pfac(pfac,ak,al,m2,0.,gam2,gam,N,omga,fcor,vtmp,
     >                               i,j,nx,ny,idz)
          do k=1,nz
            z = gz(k)
            exfac = cmplx(exp(-gam*z),0.)
            ft1(i,j,k) = cmplx(0.,ak(i))*pfac*exfac*hft(i,j)
            ft2(i,j,k) = cmplx(0.,al(j))*pfac*exfac*hft(i,j)
            ft3(i,j,k) = cmplx(-gam,0.)*pfac*exfac*hft(i,j)
          enddo

        end if

      enddo
      enddo

      do k=1,nz

        do i=1,nx-1
          do j=1,ny-1
            tran1(j,i) = ft1(i,j,k)
            tran2(j,i) = ft2(i,j,k)
            tran3(j,i) = ft3(i,j,k)
          enddo
          call cfftb(ny-1,tran1(1,i),cwrky)
          call cfftb(ny-1,tran2(1,i),cwrky)
          call cfftb(ny-1,tran3(1,i),cwrky)
          do j=1,ny-1
            wrk1(i,j) = tran1(j,i)
            wrk2(i,j) = tran2(j,i)
            wrk3(i,j) = tran3(j,i)
          enddo
        enddo

        do j=1,ny-1
          call cfftb(nx-1,wrk1(1,j),cwrkx)
          call cfftb(nx-1,wrk2(1,j),cwrkx)
          call cfftb(nx-1,wrk3(1,j),cwrkx)
          do i=1,nx-1
            dsdx = real(wrk1(i,j))
            dsdy = real(wrk2(i,j))
            dsdz = real(wrk3(i,j))
            gs(i,j,k) = gs(i,j,k) - u(i,j,k)*dsdx - v(i,j,k)*dsdy
     >                                   - w(i,j,k)*dsdz
          enddo
        enddo

        if (iu.eq.1) then
          do j=1,ny-1
            gs(nx,j,k) = gs(1,j,k)
          enddo
        else if (iv.eq.1) then
          do i=1,nx-1
            gs(i,ny,k) = gs(i,1,k)
          enddo
        end if

      enddo

      if (idz.ne.0) then

        do trm=1,3
            
          do j=1,ny-1
          do i=1,nx-1

            omga = ak(i)*U0+al(j)*V0
            if (requ(abs(omga),abs(fcor))) then
              m2 = 0.
              gam2 = 0.
            else
              m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >                   / (omga**2-fcor**2)
              gam2 = -m2
            end if

            if (mzero2(omga,N,fcor,m2)) then

              do k=1,nz
                ft1(i,j,k) = 0.
                ft2(i,j,k) = 0.
              enddo

            else if (m2.gt.0.) then

              m = sqrt(m2)*omga/abs(omga)
              if (trm.eq.1) then
                ctmp(1:1) = 'u'
                call calc_pfac(pfac1,ak,al,m2,m,gam2,0.,N,omga,fcor,
     >                              ctmp,i,j,nx,ny,1)
                call calc_pfac(pfac2,ak,al,m2,m,gam2,0.,N,omga,fcor,
     >                              vtmp,i,j,nx,ny,0)
                pfac2 = cmplx(0.,ak(i))*pfac2
              else if (trm.eq.2) then
                ctmp(1:1) = 'v'
                call calc_pfac(pfac1,ak,al,m2,m,gam2,0.,N,omga,fcor,
     >                              ctmp,i,j,nx,ny,1)
                call calc_pfac(pfac2,ak,al,m2,m,gam2,0.,N,omga,fcor,
     >                              vtmp,i,j,nx,ny,0)
                pfac2 = cmplx(0.,al(j))*pfac2
              else
                ctmp(1:1) = 'w'
                call calc_pfac(pfac1,ak,al,m2,m,gam2,0.,N,omga,fcor,
     >                              ctmp,i,j,nx,ny,1)
                call calc_pfac(pfac2,ak,al,m2,m,gam2,0.,N,omga,fcor,
     >                              vtmp,i,j,nx,ny,1)
              end if

              do k=1,nz
                z = gz(k)
                exfac = cmplx(cos(m*z),sin(m*z))
                ft1(i,j,k) = pfac1*exfac*hft(i,j)
                ft2(i,j,k) = pfac2*exfac*hft(i,j)
              enddo

            else

              gam = sqrt(gam2)
              if (trm.eq.1) then
                ctmp(1:1) = 'u'
                call calc_pfac(pfac1,ak,al,m2,0.,gam2,gam,N,omga,fcor,
     >                                    ctmp,i,j,nx,ny,1)
                call calc_pfac(pfac2,ak,al,m2,0.,gam2,gam,N,omga,fcor,
     >                                    vtmp,i,j,nx,ny,0)
                pfac2 = cmplx(0.,ak(i))*pfac2
              else if (trm.eq.2) then
                ctmp(1:1) = 'v'
                call calc_pfac(pfac1,ak,al,m2,0.,gam2,gam,N,omga,fcor,
     >                                    ctmp,i,j,nx,ny,1)
                call calc_pfac(pfac2,ak,al,m2,0.,gam2,gam,N,omga,fcor,
     >                                    vtmp,i,j,nx,ny,0)
                pfac2 = cmplx(0.,al(j))*pfac2
              else
                ctmp(1:1) = 'w'
                call calc_pfac(pfac1,ak,al,m2,0.,gam2,gam,N,omga,fcor,
     >                                    ctmp,i,j,nx,ny,1)
                call calc_pfac(pfac2,ak,al,m2,0.,gam2,gam,N,omga,fcor,
     >                                    vtmp,i,j,nx,ny,1)
              end if

              do k=1,nz
                z = gz(k)
                exfac = cmplx(exp(-gam*z),0.)
                ft1(i,j,k) = pfac1*exfac*hft(i,j)
                ft2(i,j,k) = pfac2*exfac*hft(i,j)
              enddo

            end if

          enddo
          enddo

          do k=1,nz

            do i=1,nx-1
              do j=1,ny-1
                tran1(j,i) = ft1(i,j,k)
                tran2(j,i) = ft2(i,j,k)
              enddo
              call cfftb(ny-1,tran1(1,i),cwrky)
              call cfftb(ny-1,tran2(1,i),cwrky)
              do j=1,ny-1
                wrk1(i,j) = tran1(j,i)
                wrk2(i,j) = tran2(j,i)
              enddo
            enddo 

            do j=1,ny-1
              call cfftb(nx-1,wrk1(1,j),cwrkx)
              call cfftb(nx-1,wrk2(1,j),cwrkx)
              do i=1,nx-1
                du = real(wrk1(i,j))
                ds = real(wrk2(i,j))
                gs(i,j,k) = gs(i,j,k) - du*ds 
              enddo
            enddo

          enddo

        enddo

      end if

      return
      end

*----------------------------------------------------------------------

      subroutine calc_w2(w2,wf,dwf,dgu,dgv,gb,gw,gh,rcof,Gs,c1,c2,
     >                      c1w,c2w,tran,wrk,Gfac,cwrkx,cwrky,ak,al,
     >                      sgz,wgz,dz,nx,ny,nz,is)

      integer nx,ny,nz,is
      real cwrkx(4*nx+15),cwrky(4*ny+15),ak(nx-1),al(ny-1)
      real w2(nx-1,ny-1,nz),sgz(nz),wgz(nz),dz
      complex wf(nx-1,ny-1,nz),dwf(nx-1,ny-1,nz)
      complex dgu(nx-1,ny-1,nz),dgv(nx-1,ny-1,nz)
      complex gb(nx-1,ny-1,nz),gw(nx-1,ny-1,nz)
      complex gh(nx-1,ny-1),rcof(nx-1,ny-1),Gs(nx-1,ny-1)
      complex tran(ny-1,nx-1),wrk(nx-1,ny-1),Gfac(nz)
      complex c1(nx-1,ny-1),c2(nx-1,ny-1),c1w(nx-1,ny-1)
      complex c2w(nx-1,ny-1)

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      real m,m2,gam,gam2,omga,z,tau,dz1
      complex fac(4),exfac1,exfac2,exfac3,exfac4,Gave
      logical mzero2,requ
      integer i,j

      do j=1,ny-1   
      do i=1,nx-1

        omga = ak(i)*U0+al(j)*V0
        if (requ(abs(omga),abs(fcor))) then
          m2 = 0.
          gam2 = 0.
        else
          m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >               / (omga**2-fcor**2)
          gam2 = -m2
          fac(1) = cmplx(omga*ak(i),fcor*al(j))
          fac(2) = cmplx(omga*al(j),-fcor*ak(i))
          fac(3) = cmplx(-ak(i)**2-al(j)**2,0.)
          fac(4) = cmplx(0.,-omga*(ak(i)**2+al(j)**2))
          do k=1,nz
            Gfac(k) = (fac(1)*dgu(i,j,k) + fac(2)*dgv(i,j,k)
     >                 + fac(3)*gb(i,j,k) + fac(4)*gw(i,j,k))
     >                     /(omga**2-fcor**2)
          enddo
        end if

        if (mzero2(omga,N,fcor,m2)) then

          do k=1,nz
           
            wf(i,j,k) = 0.         

          enddo

        else if (m2.gt.0.) then

          m = sqrt(m2)*omga/abs(omga)

          do k=1,nz

            z = wgz(k)
            exfac1 = cmplx(cos(m*z),sin(m*z))
            exfac2 = cmplx(cos(m*z),-sin(m*z))

            wf(i,j,k) = (gh(i,j)+c1w(i,j)-rcof(i,j))*exfac1
     >                     + (rcof(i,j)-c2w(i,j))*exfac2

            tau = wgz(k) + 0.5*dz
            exfac3 = cmplx(cos(m*tau),sin(m*tau))
            exfac4 = cmplx(cos(m*tau),-sin(m*tau)) 
            c1w(i,j) = c1w(i,j) + cmplx(0.,0.5/m)*Gfac(k)
     >                          * exfac4*dz
            c2w(i,j) = c2w(i,j) + cmplx(0.,0.5/m)*Gfac(k)
     >                          * exfac3*dz

          enddo

        else 

          gam = sqrt(gam2)

          do k=1,nz

            z = wgz(k)
            exfac1 = cmplx(exp(-gam*z),0.)
            exfac2 = cmplx(exp(gam*z),0.)

            wf(i,j,k) = (gh(i,j)+c1w(i,j)-rcof(i,j))*exfac1
     >                     + (rcof(i,j)-c2w(i,j))*exfac2

            tau = wgz(k) + 0.5*dz
            exfac3 = cmplx(exp(-gam*tau),0.)
            exfac4 = cmplx(exp(gam*tau),0.)
            c1w(i,j) = c1w(i,j) + cmplx(0.5/gam,0.)*Gfac(k)
     >                          * exfac4*dz
            c2w(i,j) = c2w(i,j) + cmplx(0.5/gam,0.)*Gfac(k)
     >                          * exfac3*dz

          enddo

        end if

      enddo
      enddo

      do k=1,nz

        do i=1,nx-1
          do j=1,ny-1
            tran(j,i) = wf(i,j,k)
          enddo
          call cfftb(ny-1,tran(1,i),cwrky)
          do j=1,ny-1
            wrk(i,j) = tran(j,i)
          enddo
        enddo

        do j=1,ny-1
          call cfftb(nx-1,wrk(1,j),cwrkx)
          do i=1,nx-1
            w2(i,j,k) = real(wrk(i,j))
          enddo
        enddo

      enddo

      do j=1,ny-1   
      do i=1,nx-1

        omga = ak(i)*U0+al(j)*V0
        if (requ(abs(omga),abs(fcor))) then
          m2 = 0.
          gam2 = 0.
        else
          m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >               / (omga**2-fcor**2)
          gam2 = -m2
          fac(1) = cmplx(omga*ak(i),fcor*al(j))
          fac(2) = cmplx(omga*al(j),-fcor*ak(i))
          fac(3) = cmplx(-ak(i)**2-al(j)**2,0.)
          fac(4) = cmplx(0.,-omga*(ak(i)**2+al(j)**2))
          do k=1,nz
            Gfac(k) = (fac(1)*dgu(i,j,k) + fac(2)*dgv(i,j,k)
     >                  + fac(3)*gb(i,j,k) + fac(4)*gw(i,j,k))
     >                     /(omga**2-fcor**2)
          enddo
        end if

        if (mzero2(omga,N,fcor,m2)) then

          do k=1,nz
           
            wf(i,j,k) = 0.         
            dwf(i,j,k) = 0.

          enddo

        else if (m2.gt.0.) then

          m = sqrt(m2)*omga/abs(omga)

          z = sgz(1)
          exfac1 = cmplx(cos(m*z),sin(m*z))
          exfac2 = cmplx(cos(m*z),-sin(m*z))
          if (is.eq.1) then
            tau = 0.25*dz
            dz1 = 0.5*dz
          else
            tau = sgz(1) - 0.5*dz
            dz1 = dz
          end if
          Gave = 0.5*(Gs(i,j)+Gfac(1))
          exfac3 = cmplx(cos(m*tau),sin(m*tau))
          exfac4 = cmplx(cos(m*tau),-sin(m*tau))
          c1(i,j) = c1(i,j) + cmplx(0.,0.5/m)*Gave*exfac4*dz1
          c2(i,j) = c2(i,j) + cmplx(0.,0.5/m)*Gave*exfac3*dz1
          
          wf(i,j,1) = (gh(i,j)+c1(i,j)-rcof(i,j))*exfac1
     >                    + (rcof(i,j)-c2(i,j))*exfac2
          dwf(i,j,1) = cmplx(0.,m)*(gh(i,j)+c1(i,j)-rcof(i,j))*exfac1
     >                    + cmplx(0.,-m)*(rcof(i,j)-c2(i,j))*exfac2

          do k=2,nz

            z = sgz(k)
            exfac1 = cmplx(cos(m*z),sin(m*z))
            exfac2 = cmplx(cos(m*z),-sin(m*z))
            tau = 0.5*(sgz(k)+sgz(k-1))
            Gave = 0.5*(Gfac(k)+Gfac(k-1))
            exfac3 = cmplx(cos(m*tau),sin(m*tau))
            exfac4 = cmplx(cos(m*tau),-sin(m*tau))
            c1(i,j) = c1(i,j) + cmplx(0.,0.5/m)*Gave*exfac4
     >                       * (sgz(k)-sgz(k-1))
            c2(i,j) = c2(i,j) + cmplx(0.,0.5/m)*Gave*exfac3
     >                       * (sgz(k)-sgz(k-1))

            wf(i,j,k) = (gh(i,j)+c1(i,j)-rcof(i,j))*exfac1
     >                     + (rcof(i,j)-c2(i,j))*exfac2
            dwf(i,j,k) = cmplx(0.,m)*(gh(i,j)+c1(i,j)-rcof(i,j))
     >                                   *exfac1
     >               + cmplx(0.,-m)*(rcof(i,j)-c2(i,j))*exfac2 

          enddo

          Gs(i,j) = Gfac(nz)

        else 

          gam = sqrt(gam2)

          z = sgz(1)
          exfac1 = cmplx(exp(-gam*z),0.)
          exfac2 = cmplx(exp(gam*z),0.)
          if (is.eq.1) then
            tau = 0.25*dz
            dz1 = 0.5*dz
          else
            tau = sgz(1) - 0.5*dz
            dz1 = dz
          end if
          Gave = 0.5*(Gs(i,j)+Gfac(1))
          exfac3 = cmplx(exp(-gam*tau),0.)
          exfac4 = cmplx(exp(gam*tau),0.)
          c1(i,j) = c1(i,j) + cmplx(0.5/gam,0.)*Gave*exfac4*dz1
          c2(i,j) = c2(i,j) + cmplx(0.5/gam,0.)*Gave*exfac3*dz1
          
          wf(i,j,1) = (gh(i,j)+c1(i,j)-rcof(i,j))*exfac1
     >                    + (rcof(i,j)-c2(i,j))*exfac2
          dwf(i,j,1) = cmplx(-gam,0.)*(gh(i,j)+c1(i,j)-rcof(i,j))
     >                                 * exfac1
     >                + cmplx(gam,0.)*(rcof(i,j)-c2(i,j))*exfac2
          
          do k=2,nz

            z = sgz(k)
            exfac1 = cmplx(exp(-gam*z),0.)
            exfac2 = cmplx(exp(gam*z),0.)
            tau = 0.5*(sgz(k)+sgz(k-1))
            Gave = 0.5*(Gfac(k)+Gfac(k-1))
            exfac3 = cmplx(exp(-gam*tau),0.)
            exfac4 = cmplx(exp(gam*tau),0.)
            c1(i,j) = c1(i,j) + cmplx(0.5/gam,0.)*Gave*exfac4
     >                         * (sgz(k)-sgz(k-1))
            c2(i,j) = c2(i,j) + cmplx(0.5/gam,0.)*Gave*exfac3
     >                         * (sgz(k)-sgz(k-1))

            wf(i,j,k) = (gh(i,j)+c1(i,j)-rcof(i,j))*exfac1 
     >                       + (rcof(i,j)-c2(i,j))*exfac2
            dwf(i,j,k) = cmplx(-gam,0.)*(gh(i,j)+c1(i,j)-rcof(i,j))
     >                                  *exfac1
     >                   + cmplx(gam,0.)*(rcof(i,j)-c2(i,j))*exfac2 

          enddo

          Gs(i,j) = Gfac(nz)

        end if

      enddo
      enddo

      return
      end            

*------------------------------------------------------------------------

      subroutine calc_o2(var,fld,ft,inf1,inf2,inf3,inf4,inf5,tran,
     >                         wrk,cwrkx,cwrky,ak,al,nx,ny,nz)

      integer nx,ny,nz
      real fld(nx-1,ny-1,nz),cwrkx(4*nx+15),cwrky(4*ny+15)
      real ak(nx-1),al(ny-1)
      complex ft(nx-1,ny-1,nz),inf1(nx-1,ny-1,nz)
      complex inf2(nx-1,ny-1,nz),inf3(nx-1,ny-1,nz)
      complex inf4(nx-1,ny-1,nz),inf5(nx-1,ny-1,nz)
      complex tran(ny-1,nx-1),wrk(nx-1,ny-1)
      character*(*) var

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      real omga,m2,fac(10)
      integer i,j,k
      logical mzero2,requ
  
      lgth = index(var,' ') - 1

      do j=1,ny-1
      do i=1,nx-1

        omga = ak(i)*U0+al(j)*V0
        if (requ(abs(omga),abs(fcor))) then
          m2 = 0.
        else
          m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >               / (omga**2-fcor**2)
        end if

        if (mzero2(omga,N,fcor,m2)) then
        
          do k=1,nz 
            ft(i,j,k) = 0.
          enddo

        else

          if (var(1:1).eq.'p') then
            fac(1) = 0.
            fac(2) = -(omga**2-fcor**2)/omga/(ak(i)**2+al(j)**2)
            fac(3) = fcor*al(j)/omga/(ak(i)**2+al(j)**2)
            fac(4) = -ak(i)/(ak(i)**2+al(j)**2)
            fac(5) = -fcor*ak(i)/omga/(ak(i)**2+al(j)**2)
            fac(6) = -al(j)/(ak(i)**2+al(j)**2)
            fac(7) = 0.
            fac(8) = 0.
            fac(9) = 0.
            fac(10) = 0.
          else if (var(1:1).eq.'b') then
            fac(1) = 0.
            fac(2) = N**2/omga
            fac(3) = 0.
            fac(4) = -1./omga
            fac(5) = 0.
            fac(6) = 0.
            fac(7) = 0.
            fac(8) = 0.
            fac(9) = 0.
            fac(10) = 0.
          else if (var(1:1).eq.'u') then 
            fac(1) = -ak(i)*omga/(omga**2-fcor**2)
            fac(2) = fcor*al(j)/(omga**2-fcor**2)
            fac(3) = 0.
            fac(4) = -omga/(omga**2-fcor**2)
            fac(5) = -fcor/(omga**2-fcor**2)
            fac(6) = 0.
            fac(7) = 0.
            fac(8) = 0.
            fac(9) = 0.
            fac(10) = 0.
          else if (var(1:1).eq.'v') then
            fac(1) = -al(j)*omga/(omga**2-fcor**2)
            fac(2) = -fcor*ak(i)/(omga**2-fcor**2)
            fac(3) = fcor/(omga**2-fcor**2)
            fac(4) = 0.
            fac(5) = 0.
            fac(6) = -omga/(omga**2-fcor**2)
            fac(7) = 0.
            fac(8) = 0.
            fac(9) = 0.
            fac(10) = 0.
          else if (var(1:1).eq.'q') then
            fac(1) = -fcor*ak(i)*m2/omga/(ak(i)**2+al(j)**2)
            fac(2) = m2*al(j)/(ak(i)**2+al(j)**2) + al(j)
            fac(3) = fcor*ak(i)/omga/(omga**2-fcor**2)
            fac(4) = -al(j)/(omga**2-fcor**2)
            fac(5) = al(j)*omga/(omga**2-fcor**2)
            fac(6) = fcor*ak(i)/(omga**2-fcor**2)
            fac(7) = -fcor/(omga**2-fcor**2)
            fac(8) = 0.
            fac(9) = 0.
            fac(10) = omga/(omga**2-fcor**2) 
          else if (var(1:1).eq.'r') then
            fac(1) = -fcor*al(j)*m2/omga/(ak(i)**2+al(j)**2)
            fac(2) = -m2*ak(i)/(ak(i)**2+al(j)**2) - ak(i)
            fac(3) = fcor*al(j)/omga/(omga**2-fcor**2)
            fac(4) = ak(i)/(omga**2-fcor**2)
            fac(5) = -ak(i)*omga/(omga**2-fcor**2)
            fac(6) = fcor*al(j)/(omga**2-fcor**2)
            fac(7) = 0.
            fac(8) = -omga/(omga**2-fcor**2)
            fac(9) = -fcor/(omga**2-fcor**2)
            fac(10) = 0.
          else if (var(1:2).eq.'s2') then
            fac(1) = 1./N**2
            fac(2) = 0.
            fac(3) = 0.
            fac(4) = -fcor/omga
            fac(5) = 0.
            fac(6) = fcor/omga/N**2
            fac(7) = 0.
            fac(8) = 0.
            fac(9) = 0.
            fac(10) = 0.
          else if (var(1:2).eq.'s3') then
            fac(1) = 1./N**2
            fac(2) = 0.
            fac(3) = 0.
            fac(4) = 0.
            fac(5) = 0.
            fac(6) = 0.
            fac(7) = 0.
            fac(8) = 0.
            fac(9) = 0.
            fac(10) = 0.
          else 
            fac(1) = 0.
            fac(2) = -1./omga
            fac(3) = 0.
            fac(4) = -1./omga
            fac(5) = 0.
            fac(6) = 0.
            fac(7) = 0.
            fac(8) = 0.
            fac(9) = 0.
            fac(10) = 0.
          end if
 
          do k=1,nz
            ft(i,j,k) = cmplx(fac(1),fac(2))*inf1(i,j,k) 
     >              + cmplx(fac(3),fac(4))*inf2(i,j,k)
     >                + cmplx(fac(5),fac(6))*inf3(i,j,k)
     >              + cmplx(fac(7),fac(8))*inf4(i,j,k)
     >                + cmplx(fac(9),fac(10))*inf5(i,j,k)
          enddo

        end if

      enddo
      enddo

      do k=1,nz

        do i=1,nx-1
          do j=1,ny-1
            tran(j,i) = ft(i,j,k)
          enddo
          call cfftb(ny-1,tran(1,i),cwrky)
          do j=1,ny-1
            wrk(i,j) = tran(j,i)
          enddo
        enddo

        do j=1,ny-1
          call cfftb(nx-1,wrk(1,j),cwrkx)
          do i=1,nx-1
            fld(i,j,k) = real(wrk(i,j))
          enddo
        enddo

      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine drag_calc(gh,rcof,hft,gus,gvs,p2,dwdz,tran1,tran2,
     >                      tran3,tran4,wrk1,wrk2,wrk3,wrk4,zs_p,
     >                      cwrkx,cwrky,ak,al,dx,dy,nx,ny) 

      integer nx,ny
      real zs_p(nx-1,ny-1),ak(nx-1),al(ny-1)
      real cwrkx(4*nx+15),cwrky(4*ny+15),dx,dy
      real p2(nx-1,ny-1)
      complex gh(nx-1,ny-1),rcof(nx-1,ny-1),hft(nx-1,ny-1)
      complex gus(nx-1,ny-1),gvs(nx-1,ny-1),dwdz(nx-1,ny-1)
      complex tran1(ny-1,nx-1),tran2(ny-1,nx-1)
      complex tran3(ny-1,nx-1),tran4(ny-1,nx-1)
      complex wrk1(nx-1,ny-1),wrk2(nx-1,ny-1)
      complex wrk3(nx-1,ny-1),wrk4(nx-1,ny-1)
 
      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      real m,gam,m2,gam2,omga,fac,sum1x,sum1y,sum2x,sum2y
      real p1,dp1dz,dhdx,dhdy
      character*(8) ctmp
      logical requ,mzero2
      complex pfac1,pfac2 

      ctmp = 'p2'

      do j=1,ny-1
      do i=1,nx-1

        omga = ak(i)*U0+al(j)*V0
        if (requ(abs(omga),abs(fcor))) then
          m2 = 0.
          gam2 = 0.
        else
          m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >               / (omga**2-fcor**2)
          gam2 = -m2
        end if

        if (mzero2(omga,N,fcor,m2)) then
          dwdz(i,j) = 0.
        else if (m2.gt.0.) then
          m = sqrt(m2)*omga/abs(omga)
          dwdz(i,j) = cmplx(0.,m)*(gh(i,j)-rcof(i,j))
     >                 + cmplx(0.,-m)*rcof(i,j)
        else
          gam = sqrt(gam2)
          dwdz(i,j) = cmplx(-gam,0.)*(gh(i,j)-rcof(i,j))
     >                 + cmplx(gam,0.)*rcof(i,j)
        end if

      enddo
      enddo

      call calc_o2(ctmp,p2,wrk2,dwdz,gus,gvs,wrk3,wrk4,tran1,wrk1,
     >                  cwrkx,cwrky,ak,al,nx,ny,1)

      ctmp = 'p1'

      do i=1,nx-1
 
        do j=1,ny-1 
 
          omga = ak(i)*U0+al(j)*V0
          if (requ(abs(omga),abs(fcor))) then
            m2 = 0.
            gam2 = 0.
          else
            m2 = (N**2-omga**2)*(ak(i)**2+al(j)**2)
     >              / (omga**2-fcor**2)
            gam2 = -m2
          end if 
 
          if (mzero2(omga,N,fcor,m2)) then
            pfac1 = 0.
            pfac2 = 0.
          else if (m2.gt.0.) then
            m = sqrt(m2)*omga/abs(omga)
            call calc_pfac(pfac1,ak,al,m2,m,gam2,0.,N,omga,fcor,ctmp,
     >                                  i,j,nx,ny,0)
            call calc_pfac(pfac2,ak,al,m2,m,gam2,0.,N,omga,fcor,ctmp,
     >                                  i,j,nx,ny,1)
          else
            gam = sqrt(gam2) 
            call calc_pfac(pfac1,ak,al,m2,0.,gam2,gam,N,omga,fcor,
     >                              ctmp,i,j,nx,ny,0)
            call calc_pfac(pfac2,ak,al,m2,0.,gam2,gam,N,omga,fcor,
     >                              ctmp,i,j,nx,ny,1)
          end if
 
          tran1(j,i) = pfac1*hft(i,j)
          tran2(j,i) = pfac2*hft(i,j)
          tran3(j,i) = cmplx(0.,ak(i))*hft(i,j)
          tran4(j,i) = cmplx(0.,al(j))*hft(i,j)
 
        enddo
 
        call cfftb(ny-1,tran1(1,i),cwrky)
        call cfftb(ny-1,tran2(1,i),cwrky)
        call cfftb(ny-1,tran3(1,i),cwrky)
        call cfftb(ny-1,tran4(1,i),cwrky)
        do j=1,ny-1
          wrk1(i,j) = tran1(j,i)
          wrk2(i,j) = tran2(j,i)
          wrk3(i,j) = tran3(j,i)
          wrk4(i,j) = tran4(j,i)
        enddo
 
      enddo

      sum1x = 0.
      sum2x = 0.
      sum1y = 0.
      sum2y = 0.
 
      do j=1,ny-1

        call cfftb(nx-1,wrk1(1,j),cwrkx)
        call cfftb(nx-1,wrk2(1,j),cwrkx)
        call cfftb(nx-1,wrk3(1,j),cwrkx)
        call cfftb(nx-1,wrk4(1,j),cwrkx)

        do i=1,nx-1

          p1 = real(wrk1(i,j))
          dp1dz = real(wrk2(i,j))
          dhdx = real(wrk3(i,j))
          dhdy = real(wrk4(i,j))

          sum1x = sum1x + p1*dhdx*dx*dy
          sum1y = sum1y + p1*dhdy*dx*dy 

          sum2x = sum2x + (p1 + dp1dz*zs_p(i,j) + p2(i,j))
     >                           * dhdx*dx*dy
          sum2y = sum2y + (p1 + dp1dz*zs_p(i,j) + p2(i,j))
     >                           * dhdy*dx*dy

        enddo

      enddo

      write(6,*) 
      write(6,*) 'Net x drag, first order = ',sum1x
      write(6,*) 'Net y drag, first order = ',sum1y 
      write(6,*) 
      write(6,*) 'Net x drag, second order = ',sum2x
      write(6,*) 'Net y drag, second order = ',sum2y
      write(6,*)
 
      return
      end

*----------------------------------------------------------------------

      subroutine calc_def(d11,d12,d13,d22,d23,d33,u,v,w,sgz,wgz,
     >                                dx,dy,nx,ny,nz)

      integer nx,ny,nz
      real d11(nx,ny,nz),d12(nx,ny,nz),d13(nx,ny,nz)
      real d22(nx,ny,nz),d23(nx,ny,nz),d33(nx,ny,nz)
      real u(nx,ny-1,nz),v(nx-1,ny,nz),w(nx-1,ny-1,nz)
      real sgz(nz),wgz(nz),dx,dy

      integer i,j,k,ip1,im1,kp1,km1,jp1,jm1,nzmax
      parameter (nzmax=201)
      real ux,uy,uz,vx,vy,vz,wx,wy,wz
      real rdx,rdy,wdz(nzmax),sdz(nzmax)

      if (nz.gt.nzmax) then
        write(6,*) 'Error:  nz too big!  Increase nzmax in calc_def.'
        stop
      end if

      do k=1,nz-1
        wdz(k) = wgz(k+1)-wgz(k) 
        sdz(k) = sgz(k+1)-sgz(k)
      enddo

      wdz(nz) = 1.
      sdz(nz) = 1.

      rdx = 1./dx
      rdy = 1./dy

      do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
      do j=1,ny-1
        jp1 = min(j+1,ny-1)
        jm1 = max(j-1,1)
      do i=1,nx-1
        ip1 = min(i+1,nx-1)
        im1 = max(i-1,1)

        ux = rdx*(u(i+1,j,k)-u(i,j,k))
        vy = rdy*(v(i,j+1,k)-v(i,j,k))
        wz = (w(i,j,kp1)-w(i,j,k))/wdz(k)
        d11(i,j,k) = (4.*ux - 2.*(vy+wz))/3.
        d22(i,j,k) = (4.*vy - 2.*(ux+wz))/3.
        d33(i,j,k) = (4.*wz - 2.*(ux+vy))/3.

        uy = rdy*(u(i,j,k)-u(i,jm1,k))
        vx = rdx*(v(i,j,k)-v(im1,j,k))
        d12(i,j,k) = uy + vx

        uz = (u(i,j,k)-u(i,j,km1))/sdz(km1)
        wx = rdx*(w(i,j,k)-w(im1,j,k))
        d13(i,j,k) = uz + wx

        vz = (v(i,j,k)-v(i,j,km1))/sdz(km1)
        wy = rdy*(w(i,j,k)-w(i,jm1,k))
        d23(i,j,k) = vz + wy

      enddo
      enddo
      enddo

      do k=1,nz
      do j=1,ny-1
        d13(1,j,k) = d13(2,j,k)
        d13(nx,j,k) = d13(nx-1,j,k)
        d12(1,j,k) = d12(2,j,k)
        d12(nx,j,k) = d12(nx-1,j,k)
      enddo
      enddo

      do k=1,nz
      do i=1,nx
        d12(i,1,k) = d12(i,2,k)
        d12(i,ny,k) = d12(i,ny-1,k)
        if (i.lt.nx) then
          d23(i,1,k) = d23(i,2,k)
          d23(i,ny,k) = d23(i,ny-1,k)
        end if
      enddo
      enddo

      do j=1,ny
      do i=1,nx
        if (j.lt.ny) then
          d13(i,j,1) = 0.
        end if
        if (i.lt.nx) then
          d23(i,j,1) = 0.
        end if
      enddo
      enddo

      return
      end 

*----------------------------------------------------------------------

      subroutine calc_km(d11,d12,d13,d22,d23,d33,b1,km,sgz,a1,Cm,
     >                        rprandl,nx,ny,nz)

      integer nx,ny,nz
      real d11(nx,ny,nz),d12(nx,ny,nz),d13(nx,ny,nz)
      real d22(nx,ny,nz),d23(nx,ny,nz),d33(nx,ny,nz)
      real b1(nx-1,ny-1,nz),km(nx-1,ny-1,nz),sgz(nz)
      real a1,Cm,rprandl

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      integer i,j,k,km1,im1,jm1,ii,jj
      real aved12,aved13,aved23,def2,N2,Cd,zl,zu,pi

      pi = 2.*asin(1.)
      zl = 0.5*pi*U0/N
      zu = 2.5*pi*U0/N

      Cd = (Cm**2)/a1*(U0/N)**3

      do k=1,nz-1
        km1 = max(k-1,1)
      do j=1,ny-1
      do i=1,nx-1

        if ((sgz(k).gt.zl).and.(sgz(k).lt.zu)) then

          aved12 = 0.25*(d12(i,j,k)+d12(i+1,j,k)+d12(i,j+1,k)
     >                         +d12(i+1,j+1,k))
          aved13 = 0.25*(d13(i,j,k)+d13(i+1,j,k)+d13(i,j,k+1)
     >                         +d13(i+1,j,k+1))
          aved23 = 0.25*(d23(i,j,k)+d23(i,j+1,k)+d23(i,j,k+1)
     >                         +d23(i,j+1,k+1)) 

          def2 = (d11(i,j,k)**2 + d22(i,j,k)**2 + d33(i,j,k)**2)/2.
     >          + aved12**2 + aved13**2 + aved23**2

          N2 = N**2 + (b1(i,j,k+1)-b1(i,j,km1))/(sgz(k+1)-sgz(km1))

          km(i,j,k) = Cd*sqrt(amax1(0.,(def2-rprandl*N2)))

        else

          km(i,j,k) = 0.

        end if

      enddo
      enddo
      enddo

      do i=1,nx-1
      do j=1,ny-1
        km(i,j,nz) = km(i,j,nz-1)
      enddo
      enddo

      do k=1,nz
        km1 = max(k-1,1)
      do j=1,ny
        jm1 = max(j-1,1)
        jj = min(j,ny-1)
      do i=1,nx
        im1 = max(i-1,1)
        ii = min(i,nx-1)

        if ((i.lt.nx).and.(j.lt.ny)) then
          d11(i,j,k) = km(i,j,k)*d11(i,j,k)
          d22(i,j,k) = km(i,j,k)*d22(i,j,k)
          d33(i,j,k) = km(i,j,k)*d33(i,j,k)
        end if

        d12(i,j,k) = 0.25*(km(ii,jj,k)+km(im1,jj,k)+km(ii,jm1,k)
     >                       +km(im1,jm1,k))*d12(i,j,k)
        if (j.lt.ny) then
          d13(i,j,k) = 0.25*(km(ii,j,k)+km(im1,j,k)+km(ii,j,km1)
     >                       +km(im1,j,km1))*d13(i,j,k)
        end if
        if (i.lt.nx) then
          d23(i,j,k) = 0.25*(km(i,jj,k)+km(i,jm1,k)+km(i,jj,km1)
     >                       +km(i,jm1,km1))*d23(i,j,k)
        end if

      enddo
      enddo
      enddo

      return
      end 

*----------------------------------------------------------------------

      subroutine u_mix(d11,d12,d13,gu,dgu,wgz,sgz,dx,dy,nx,ny,nz)

      integer nx,ny,nz
      real d11(nx,ny,nz),d12(nx,ny,nz),d13(nx,ny,nz)
      real gu(nx-1,ny-1,nz),dgu(nx-1,ny-1,nz)
      real wgz(nz),sgz(nz),dx,dy

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      integer i,j,k,kp1,km1,nzmax
      parameter (nzmax=201)
      real rdx,rdy,rdz,zrbnd(nzmax),zlbnd(nzmax),dz(nzmax),zl,pi

      if (nz-1.gt.nzmax) then
        write(6,*) 'Error:  nz too big!  Increase nzmax in u_mix.'
        stop
      end if
 
      pi = 2.*asin(1.)
      zl = 0.5*pi*U0/N

      do k=1,nz-1
        zrbnd(k) = 1.
        zlbnd(k) = 1.
        dz(k) = sgz(k+1)-sgz(k)
      enddo

      zrbnd(nz-1) = 0.
      zrbnd(1) = 2.
      zlbnd(nz-1) = 2.
      zlbnd(1) = 0.

      rdx = 1./dx
      rdy = 1./dy

      do k=1,nz-1

        if (sgz(k).gt.zl) then

          rdz = 1./(wgz(k+1)-wgz(k))

          do j=1,ny-1
          do i=2,nx-2
            gu(i,j,k) = gu(i,j,k) + 0.5*rdx*
     >                       (d11(i+1,j,k)-d11(i-1,j,k))
          enddo
          enddo 

          do j=2,ny-2
          do i=1,nx-1
            gu(i,j,k) = gu(i,j,k) + 0.5*rdy*(d12(i,j+1,k)-d12(i,j,k)
     >                            + d12(i+1,j+1,k)-d12(i+1,j,k))
          enddo
          enddo

          do j=1,ny-1
          do i=1,nx-1
            gu(i,j,k) = gu(i,j,k) + 0.5*rdz*(d13(i,j,k+1)-d13(i,j,k)
     >                            + d13(i+1,j,k+1)-d13(i+1,j,k))
          enddo
          enddo

        end if

      enddo

      do k=1,nz-1
        kp1 = min(k+1,nz-1)
        km1 = max(k-1,1)
      do j=1,ny-1
      do i=1,nx-1
        dgu(i,j,k) = dgu(i,j,k) + 0.5*zrbnd(k)*(gu(i,j,kp1)
     >                              - gu(i,j,k))/dz(k)
     >                          + 0.5*zlbnd(k)*(gu(i,j,k)
     >                           - gu(i,j,km1))/dz(km1)
      enddo
      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine v_mix(d12,d22,d23,gv,dgv,wgz,sgz,dx,dy,nx,ny,nz)

      integer nx,ny,nz
      real d12(nx,ny,nz),d22(nx,ny,nz),d23(nx,ny,nz)
      real gv(nx-1,ny-1,nz),dgv(nx-1,ny-1,nz)
      real wgz(nz),sgz(nz),dx,dy

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      integer i,j,k,kp1,km1,nzmax
      parameter (nzmax=201)
      real rdx,rdy,rdz,zrbnd(nzmax),zlbnd(nzmax),dz(nzmax),zl,pi

      if (nz-1.gt.nzmax) then
        write(6,*) 'Error:  nz too big!  Increase nzmax in v_mix.'
        stop
      end if

      pi = 2.*asin(1.)
      zl = 0.5*pi*U0/N

      do k=1,nz-1
        zrbnd(k) = 1.
        zlbnd(k) = 1.
        dz(k) = sgz(k+1)-sgz(k)
      enddo

      zrbnd(nz-1) = 0.
      zrbnd(1) = 2.
      zlbnd(nz-1) = 2.
      zlbnd(1) = 0.

      rdx = 1./dx
      rdy = 1./dy

      do k=1,nz-1

        if (sgz(k).gt.zl) then

          rdz = 1./(wgz(k+1)-wgz(k))

          do j=1,ny-1
          do i=2,nx-2
            gv(i,j,k) = gv(i,j,k) + 0.5*rdx*(d12(i+1,j,k)-d12(i,j,k)
     >                            + d12(i+1,j+1,k)-d12(i,j+1,k))
          enddo
          enddo

          do j=2,ny-2
          do i=1,nx-1
            gv(i,j,k) = gv(i,j,k) + 0.5*rdy*
     >                        (d22(i,j+1,k)-d22(i,j-1,k))
          enddo
          enddo

          do j=1,ny-1
          do i=1,nx-1
            gv(i,j,k) = gv(i,j,k) + 0.5*rdz*(d23(i,j,k+1)-d23(i,j,k)
     >                            + d23(i,j+1,k+1)-d23(i,j+1,k))
          enddo
          enddo

        end if

      enddo
  
      do k=1,nz-1
        kp1 = min(k+1,nz-1)
        km1 = max(k-1,1)
      do j=1,ny-1
      do i=1,nx-1
        dgv(i,j,k) = dgv(i,j,k) + 0.5*zrbnd(k)*(gv(i,j,kp1)
     >                              - gv(i,j,k))/dz(k)
     >                          + 0.5*zlbnd(k)*(gv(i,j,k)
     >                           - gv(i,j,km1))/dz(km1)
      enddo
      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine w_mix(d13,d23,d33,gw,sgz,dx,dy,nx,ny,nz)

      integer nx,ny,nz
      real d13(nx,ny,nz),d23(nx,ny,nz),d33(nx,ny,nz)
      real gw(nx-1,ny-1,nz),sgz(nz),dx,dy

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      integer i,j,k,kp1,km1,nzmax
      parameter (nzmax=201)
      real rdx,rdy,zrbnd(nzmax),zlbnd(nzmax),dz(nzmax),zl,pi

      if (nz-1.gt.nzmax) then
        write(6,*) 'Error:  nz too big!  Increase nzmax in w_mix.'
        stop
      end if

      pi = 2.*asin(1.)
      zl = 0.5*pi*U0/N

      do k=1,nz-1
        zrbnd(k) = 1.
        zlbnd(k) = 1.
        dz(k) = sgz(k+1)-sgz(k)
      enddo

      zrbnd(nz-1) = 0.
      zrbnd(1) = 2.
      zlbnd(nz-1) = 2.
      zlbnd(1) = 0. 

      rdx = 1./dx
      rdy = 1./dy

      do k=1,nz-1
        kp1 = min(k+1,nz-1)
        km1 = max(k-1,1)

        if (sgz(k).gt.zl) then

          do j=1,ny-1
          do i=2,nx-2
            gw(i,j,k) = gw(i,j,k) + 0.5*rdx*(d13(i+1,j,k)-d13(i,j,k)
     >                            + d13(i+1,j,k+1)-d13(i,j,k+1))
          enddo
          enddo

          do j=2,ny-2
          do i=1,nx-1
            gw(i,j,k) = gw(i,j,k) + 0.5*rdy*(d23(i,j+1,k)-d23(i,j,k)
     >                            + d23(i,j+1,k+1)-d23(i,j,k+1))
          enddo
          enddo

          do j=1,ny-1
          do i=1,nx-1
            gw(i,j,k) = gw(i,j,k) + 0.5*zrbnd(k)*(d33(i,j,kp1)
     >                             - d33(i,j,k))/dz(k)
     >                        + 0.5*zlbnd(k)*(d33(i,j,k)
     >                           - d33(i,j,km1))/dz(km1)
          enddo
          enddo

        end if

      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine b_mix(km,b1,gb,dgb,sgz,wgz,rprandl,dx,dy,nx,ny,nz)

      integer nx,ny,nz
      real km(nx-1,ny-1,nz),b1(nx-1,ny-1,nz),gb(nx-1,ny-1,nz)
      real dgb(nx-1,ny-1,nz),sgz(nz),wgz(nz),rprandl,dx,dy

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      integer i,j,k,nzmax,kp1,km1
      parameter (nzmax=201)
      real rdx,rdy,dz(nzmax),zrbnd(nzmax),zlbnd(nzmax),zl,pi

      if (nz-1.gt.nzmax) then
        write(6,*) 'Error:  nz too big!  Increase nzmax in b_mix.'
        stop
      end if

      pi = 2.*asin(1.)
      zl = 0.5*pi*U0/N

      do k=1,nz-1
        zrbnd(k) = 1.
        zlbnd(k) = 1.
        dz(k) = sgz(k+1)-sgz(k)
      enddo

      zrbnd(nz-1) = 0.
      zrbnd(1) = 2.
      zlbnd(nz-1) = 2.
      zlbnd(1) = 0.

      rdx = 1./dx
      rdy = 1./dy

      do k=1,nz-1

        if (sgz(k).gt.zl) then

          do j=1,ny-1
          do i=2,nx-2
            gb(i,j,k) = gb(i,j,k) + rprandl*rdx*(
     >                         0.5*(km(i+1,j,k)+km(i,j,k))
     >                    *rdx*(b1(i+1,j,k)-b1(i,j,k))
     >                        - 0.5*(km(i,j,k)+km(i-1,j,k))
     >                    *rdx*(b1(i,j,k)-b1(i-1,j,k)))
          enddo
          enddo 

          do j=2,ny-2
          do i=1,nx-1
            gb(i,j,k) = gb(i,j,k) + rprandl*rdy*(
     >                         0.5*(km(i,j+1,k)+km(i,j,k))
     >                    *rdy*(b1(i,j+1,k)-b1(i,j,k))
     >                        - 0.5*(km(i,j,k)+km(i,j-1,k))
     >                    *rdy*(b1(i,j,k)-b1(i,j-1,k)))
          enddo
          enddo

          if (k.eq.1) then
            do j=1,ny-1
            do i=1,nx-1
              gb(i,j,k) = gb(i,j,k) + rprandl*0.5
     >                        *(km(i,j,k+1)+km(i,j,k))
     >                     *(b1(i,j,k+1)-b1(i,j,k))/dz(1)/wgz(2)
            enddo 
            enddo
          else 
            do j=1,ny-1
            do i=1,nx-1
              gb(i,j,k) = gb(i,j,k) + rprandl*(
     >                        0.5*(km(i,j,k+1)+km(i,j,k))
     >                   *(b1(i,j,k+1)-b1(i,j,k))/dz(k)
     >                      - 0.5*(km(i,j,k)+km(i,j,k-1))
     >                   *(b1(i,j,k)-b1(i,j,k-1))/dz(k-1))
     >                      /(wgz(k+1)-wgz(k))
            enddo
            enddo
          end if

        end if

      enddo

      do k=1,nz-1
        kp1 = min(k+1,nz-1)
        km1 = max(k-1,1)
      do j=1,ny-1
      do i=1,nx-1
        dgb(i,j,k) = dgb(i,j,k) + 0.5*zrbnd(k)*(gb(i,j,kp1)
     >                              - gb(i,j,k))/dz(k) 
     >                          + 0.5*zlbnd(k)*(gb(i,j,k)
     >                           - gb(i,j,km1))/dz(km1)
      enddo
      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine interp(idcdf1,idcdf2,a,ispc,ifld1,ifld2,iter,
     >                 isgz,iwgz,isclr,nx,ny,nxout,nyout,nzout,
     >                 iout,jout)

      include '/usr/local/include/netcdf.inc'

      integer ispc
      real a(ispc)

      integer idcdf1,idcdf2,ifld1,ifld2,iter(3),isgz,iwgz,isclr
      integer nx,ny,nxout,nyout,nzout,iout,jout

      integer idvar1,idvar2,oneid,igz,izs,iu,iv,iw
      integer ndims,nvars,ngatts,rec_dim,ierr,k,q
      integer vtyp,nvdims,vdims(4),nvatts,ist(4),iln(4)
      real ztop
      character*(80) vname

      ztop = a(isclr+2-1)

      oneid = ncdid(idcdf1,'one',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  failed to find id for dimension one.' 
        stop 
      end if

      call ncinq(idcdf1,ndims,nvars,ngatts,rec_dim,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem with file inquiry.'
        stop
      end if

      k = 1
 567  if (k.le.nvars) then

        vname(1:len(vname)) = ' '
        do q=1,4
          vdims(q) = 0
        enddo
        call ncvinq(idcdf1,k,vname,vtyp,nvdims,vdims,nvatts,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Error:  problem with variable inquiry.'
          stop
        end if

        if ((nvdims.eq.4).and.
     >          ((vdims(1).ne.oneid).and.(vdims(2).ne.oneid)
     >                  .and.(vdims(3).ne.oneid))) then

          idvar1 = k
          idvar2 = ncvid(idcdf2,vname,ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem getting variable id.'
            stop
          end if

          izs = iter(1)
          igz = isgz
          iu = 0
          iv = 0
          iw = 0
          if (vname(1:1).eq.'u') then
            iu = 1
            izs = iter(2)
          else if (vname(1:1).eq.'v') then
            iv = 1
            izs = iter(3)
          else if (vname(1:1).eq.'w') then
            iw = 1
            igz = iwgz
          end if

          ist(1) = 1
          ist(2) = 1
          ist(3) = 1
          ist(4) = 1

          iln(1) = nxout-1+iu
          iln(2) = nyout-1+iv
          iln(3) = nzout-1+iw
          iln(4) = 1

          call ncvgt(idcdf1,k,ist,iln,a(ifld1),ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem reading variable.'
            stop
          end if

          call interp_calc(a(ifld1),a(ifld2),a(izs),a(igz),ztop,
     >                          nx,ny,nxout,nyout,nzout,
     >                          iu,iv,iw,iout,jout)

          call ncvpt(idcdf2,idvar2,ist,iln,a(ifld2),ierr)
          if (ierr.ne.0) then
            write(6,*) 'Error:  problem writing variable.'
            stop
          end if

        end if

        k = k+1
        goto 567

      end if

      return 
      end

*---------------------------------------------------------------------

      subroutine interp_calc(fld1,fld2,zs,gz,ztop,nx,ny,nxout,nyout,
     >                              nzout,iu,iv,iw,iout,jout)

      integer nx,ny,nxout,nyout,nzout,iu,iv,iw,iout,jout
      real fld1(nxout-1+iu,nyout-1+iv,nzout-1+iw)
      real fld2(nxout-1+iu,nyout-1+iv,nzout-1+iw)
      real zs(nx-1+iu,ny-1+iv),gz(nzout-1+iw),ztop

      integer i,j,k,kptr,kplus,ktop
      real z0,tfac,z1,z2,z3,z4

      do j=1,nyout-1+iv
      do i=1,nxout-1+iu

        z0 = zs(iout+i-1,jout+j-1)
        tfac = (ztop-z0)/ztop
        kptr = 1

        do k=1,nzout-1+iw
 
          z = z0 + tfac*gz(k)

          kplus = 0
 789      if ((kplus.eq.0).and.(kptr.le.nzout-1+iw)) then
            if (gz(kptr).gt.z) then
              kplus = kptr
            else
              kptr = kptr + 1
            end if
            goto 789
          end if

          if ((kplus.eq.1).or.(kplus.eq.2)) then
            z3 = gz(3)
            z2 = gz(2)
            z1 = gz(1)
            fld2(i,j,k) = (z-z3)*(z-z2)/(z1-z3)/(z1-z2)*fld1(i,j,1)
     >                + (z-z3)*(z-z1)/(z2-z3)/(z2-z1)*fld1(i,j,2)
     >              + (z-z2)*(z-z1)/(z3-z2)/(z3-z1)*fld1(i,j,3)
          else if ((kplus.eq.0).or.(kplus.eq.nzout-1+iw)) then
            ktop = nzout-1+iw
            z3 = gz(ktop)
            z2 = gz(ktop-1)
            z1 = gz(ktop-2)
            fld2(i,j,k) = (z-z3)*(z-z2)/(z1-z3)/(z1-z2)
     >                           * fld1(i,j,ktop-2) 
     >                + (z-z3)*(z-z1)/(z2-z3)/(z2-z1)
     >                         * fld1(i,j,ktop-1)
     >              + (z-z2)*(z-z1)/(z3-z2)/(z3-z1)
     >                       * fld1(i,j,ktop)
          else 
            z4 = gz(kplus+1)
            z3 = gz(kplus)
            z2 = gz(kplus-1)
            z1 = gz(kplus-2) 
            fld2(i,j,k) = (z-z4)*(z-z3)*(z-z2)
     >              /(z1-z4)/(z1-z3)/(z1-z2)*fld1(i,j,kplus-2)
     >                  + (z-z4)*(z-z3)*(z-z1)
     >              /(z2-z4)/(z2-z3)/(z2-z1)*fld1(i,j,kplus-1)
     >                  + (z-z4)*(z-z2)*(z-z1)
     >              /(z3-z4)/(z3-z2)/(z3-z1)*fld1(i,j,kplus)
     >                  + (z-z3)*(z-z2)*(z-z1)
     >              /(z4-z3)/(z4-z2)/(z4-z1)*fld1(i,j,kplus+1)
          end if

        enddo

      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine set_outfl(idcdf,nx,ny,nz,dx,dy,dz,x0,y0,sgz,wgz,
     >                asclr,tgrid,sclr,var,units,outfl,nsclr,nvout)

      include '/usr/local/include/netcdf.inc'

      integer idcdf,nx,ny,nz,nsclr,nvout
      real dx,dy,dz,x0,y0,sgz(nz),wgz(nz),asclr(nsclr),tgrid
      character*(*) sclr(nsclr),var(nvout),units(nvout),outfl

      real xmin,ymin,zmin,tmp
      character*(80) cdlfl,command,command2

      lgth = index(outfl,' ') - 1

      cdlfl(1:len(cdlfl)) = ' '
      cdlfl = outfl(1:lgth)//'.cdl'

      iunit = 48

      open(unit=iunit,file=cdlfl(1:lgth+4),status='new')

      write(iunit,*) 'netcdf lin_soln{'
      write(iunit,*)

      write(iunit,*) 'dimensions:'
      write(iunit,*)
      write(iunit,*) 'nx=',nx-1,';'
      write(iunit,*) 'ny=',ny-1,';'
      write(iunit,*) 'nz=',nz-1,';'
      write(iunit,*) 'nxp1=',nx,';'
      write(iunit,*) 'nyp1=',ny,';'
      write(iunit,*) 'nzp1=',nz,';'
      write(iunit,*) 'one=',1,';'
      write(iunit,*) 'time=UNLIMITED;'
      write(iunit,*)

      write(iunit,*) 'variables:'
      write(iunit,*)
      write(iunit,*) 'float tgrid(one);'
      write(iunit,*) 'tgrid:no_button=1;'
      write(iunit,*)

      do k=1,2    
        lgth = index(sclr(k),' ') - 1
        write(iunit,*) 'float ',sclr(k)(1:lgth),'(one);'
        write(iunit,*) sclr(k)(1:lgth),':no_button=1;'
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
      write(iunit,*) 'zbot_p:x_min=',x0+0.5*dx,';'
      write(iunit,*) 'zbot_p:y_min=',y0+0.5*dy,';'
      write(iunit,*) 'zbot_p:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_u(ny,nxp1);'
      write(iunit,*) 'zbot_u:units="m";'
      write(iunit,*) 'zbot_u:def="terrain at u points";'
      write(iunit,*) 'zbot_u:x_min=',x0,';'
      write(iunit,*) 'zbot_u:y_min=',y0+0.5*dy,';'
      write(iunit,*) 'zbot_u:no_button=1;'
      write(iunit,*)
      write(iunit,*) 'float zbot_v(nyp1,nx);'
      write(iunit,*) 'zbot_v:units="m";'
      write(iunit,*) 'zbot_v:def="terrain at v points";'
      write(iunit,*) 'zbot_v:x_min=',x0+0.5*dx,';'
      write(iunit,*) 'zbot_v:y_min=',y0,';'
      write(iunit,*) 'zbot_v:no_button=1;'
      write(iunit,*)

      do k=1,nvout

        xmin = x0 + 0.5*dx
        ymin = y0 + 0.5*dy
        zmin = sgz(1)
 
        lgth = index(var(k),' ') - 1

        if (var(k)(lgth:lgth).eq.'0') then

          if (var(k)(1:lgth).eq.'u0') then
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(one,nz,ny,nxp1);'
            xmin = x0
          else if (var(k)(1:lgth).eq.'v0') then
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(one,nz,nyp1,nx);'
            ymin = y0 
          else if (var(k)(1:lgth).eq.'w0') then
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
            xmin = x0
          else if (var(k)(1:1).eq.'v') then
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(time,nz,nyp1,nx);'
            ymin = y0 
          else if (var(k)(1:1).eq.'w') then
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(time,nzp1,ny,nx);'
            zmin = wgz(1)
          else
            write(iunit,*) 'float ',var(k)(1:lgth),
     >                       '(time,nz,ny,nx);'
          end if

        end if

        lgth2 = index(units(k),' ') - 1
        if (lgth2.eq.0) lgth2 = 1
        write(iunit,*) var(k)(1:lgth),':units="',
     >              units(k)(1:lgth2),'";'     

        write(iunit,*) var(k)(1:lgth),':x_min=',xmin,';'
        write(iunit,*) var(k)(1:lgth),':y_min=',ymin,';'
        write(iunit,*) var(k)(1:lgth),':z_min=',zmin,';'

        if ((var(k)(lgth:lgth).ne.'1')
     >           .and.(var(k)(lgth:lgth).ne.'2')
     >                .and.(var(k)(lgth:lgth).ne.'3')) then
          write(iunit,*) var(k)(1:lgth),':no_button=1;'
        end if

        write(iunit,*) 

      enddo

      write(iunit,*) '//global attributes:'
      write(iunit,*)
      write(iunit,*) ':x_min=',x0,';'
      write(iunit,*) ':x_max=',x0 + (nx-1)*dx,';'
      write(iunit,*) ':x_delta=',dx,';'
      write(iunit,*) ':x_units="m";'
      write(iunit,*) ':x_label="x";'
      write(iunit,*) ':x_display_units="km";'
      write(iunit,*)
      write(iunit,*) ':y_min=',y0,';'
      write(iunit,*) ':y_max=',y0 + (ny-1)*dy,';'
      write(iunit,*) ':y_delta=',dy,';'
      write(iunit,*) ':y_units="m";'
      write(iunit,*) ':y_label="y";'
      write(iunit,*) ':y_display_units="km";'
      write(iunit,*)
      write(iunit,*) ':z_min=',0.,';'
      write(iunit,*) ':z_max=',0. + (nz-1)*dz,';'
      write(iunit,*) ':z_delta=',dz,';'
      write(iunit,*) ':z_units="m";'
      write(iunit,*) ':z_label="z";'
      write(iunit,*) ':z_display_units="km";'
      write(iunit,*)

      lgth = index(outfl,' ') - 1
      write(iunit,*) ':runname="steady 2nd-order solution ',
     >                     outfl(1:lgth),'";'
      write(iunit,*)
      write(iunit,*) '}'

      close(iunit)

      outfl = outfl(1:lgth)//'.cdf'

      command(1:len(command)) = ' '
      write(command,*) 'ncgen -o ',outfl(1:lgth+4),' ',
     >                         cdlfl(1:lgth+4)
      call system(command)
      write(command2,*) 'rm -f ',cdlfl(1:lgth+4)
      call system(command2)

      idcdf = ncopn(outfl(1:lgth+4),NCWRITE,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem opening netcdf file.'
        stop
      end if

      idvar = ncvid(idcdf,'tgrid',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting variable id for tgrid.'
        stop
      end if
      call ncvpt1(idcdf,idvar,1,tgrid,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem writing tgrid.'
        stop
      end if

      idvar = ncvid(idcdf,'time',ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting variable id for time.'
        stop
      end if
      call ncvpt1(idcdf,idvar,1,0.,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem writing time.'
        stop
      end if

      do k=1,2    
        lgth = index(sclr(k),' ') - 1
        idvar = ncvid(idcdf,sclr(k)(1:lgth),ierr)
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

*------------------------------------------------------------------------

      subroutine get_ter(asclr,zs_p,zs_u,zs_v,x0,y0,dx,dy,
     >                               nx,ny,nsclr)

      integer nsclr,nx,ny
      real asclr(nsclr),zs_p(nx-1,ny-1),zs_u(nx,ny-1),zs_v(nx-1,ny)
      real x0,y0,dx,dy

      real x,xu,xv,y,yu,yv,pi,ys,yn,ytmp,rad,a,b 

      real type,h,a1,a2,b1,b2,bv,y_len,xc,yc

      type = asclr(3)
      h = asclr(4)
      a1 = asclr(5)
      a2 = asclr(6)
      b1 = asclr(7)
      b2 = asclr(8)
      y_len = asclr(9)
      xc = asclr(10)
      yc = asclr(11)

      pi = 2.*asin(1.)

      if (type.eq.0.) then

        do j = 1,ny-1
        do i = 1,nx-1
         zs_p(i,j) = 0.
         zs_u(i,j) = 0.
         zs_v(i,j) = 0.
        enddo
        enddo
  
        do j=1,ny-1
          zs_u(nx,j) = 0.
        enddo
        do i=1,nx-1
          zs_v(i,ny) = 0.
        enddo
 
      else if (type.eq.1.) then

        do i=1,nx-1
          x  = x0 + dx*(i-1) + 0.5*dx - xc
          xu = x0 + dx*(i-1) - xc
          xv = x
        do j=1,ny-1
          zs_p(i,j) = h*a1**2 / (x**2+a1**2)
          zs_u(i,j) = h*a1**2 / (xu**2+a1**2)
          zs_v(i,j) = zs_p(i,j)
        enddo
        enddo

        xu = x0 + dx*(nx-1) - xc
        do j=1,ny-1
          zs_u(nx,j) = h*a1**2 / (xu**2+a1**2)
        enddo

        do i=1,nx-1
          xv = x0 + dx*(i-1) + 0.5*dx - xc
          zs_v(i,ny) = h*a1**2 / (xv**2+a1**2)
        enddo

      else if (type.eq.2.) then

        do j=1,ny-1
          y  = y0 + dy*(j-1) + 0.5*dy - yc
          yv = y0 + dy*(j-1) - yc
          yu = y
        do i=1,nx-1
          zs_p(i,j) = h*b1**2 / (y**2+b1**2)
          zs_u(i,j) = zs_p(i,j)
          zs_v(i,j) = h*b1**2 / (yv**2+b1**2)
        enddo
        enddo

        do j=1,ny-1
          yu = y0 + dy*(j-1) + 0.5*dy - yc
          zs_u(nx,j) = h*b1**2 / (yu**2+b1**2)
        enddo

        yv = y0 + dy*(ny-1) - yc
        do i=1,nx-1 
          zs_v(i,ny) = h*b1**2 / (yv**2+b1**2)
        enddo

      else if (type.eq.3.) then

        do j=1,ny-1
          y  = y0 + dy*(j-1) + 0.5*dy - yc
          yv = y0 + dy*(j-1) - yc
          yu = y
        do i=1,nx-1
          x  = x0 + dx*(i-1) + 0.5*dx - xc 
          xu = x0 + dx*(i-1) - xc
          xv = x
          zs_p(i,j) = h/(1+(x/a1)**2+(y/b1)**2)**(1.5)
          zs_u(i,j) = h/(1+(xu/a1)**2+(yu/b1)**2)**(1.5)
          zs_v(i,j) = h/(1+(xv/a1)**2+(yv/b1)**2)**(1.5)
        enddo
        enddo

        xu = x0 + dx*(nx-1) - xc
        do j=1,ny-1
          yu = y0 + dy*(j-1) + 0.5*dy - yc
          zs_u(nx,j) = h/(1+(xu/a1)**2+(yu/b1)**2)**(1.5)
        enddo

        yv = y0 + dy*(ny-1) - yc
        do i=1,nx-1
          xv = x0 + (i-1)*dx + 0.5*dx - xc
          zs_v(i,ny) = h/(1+(xv/a1)**2+(yv/b1)**2)**(1.5)
        enddo
          
      else if (type.eq.4.) then

        ys = yc - y_len
        yn = yc + y_len

        do j=1,ny

          b = b1
          bv = b1

          ytmp = y0 + dy*(j-1) + 0.5*dy
          if (ytmp.lt.ys) then
            y = ys - ytmp
            b = b2
          else if (ytmp.gt.yn) then
            y = ytmp - yn
          else
            y = 0.
          end if

          ytmp = ytmp - 0.5*dy
          if (ytmp.lt.ys) then
            yv = ys - ytmp
            bv = b2
          else if (ytmp.gt.yn) then
            yv = ytmp - yn
          else
            yv = 0.
          end if

          yu = y

          do i=1,nx

            x = x0 + dx*(i-1) + 0.5*dx - xc
            xu = x - 0.5*dx
            xv = x

            if (x.le.0.) then 
              a = a1
            else
              a = a2
            end if

            rad = sqrt((x/4./a)**2 + (y/4./b)**2)

            if ((i.le.nx-1).and.(j.le.ny-1)) then
              if (rad.lt.1.) then
                zs_p(i,j) = h/16.*(1+cos(pi*rad))**4
              else
                zs_p(i,j) = 0.
              end if
            end if

            rad = sqrt((xv/4./a)**2 + (yv/4./bv)**2)

            if (i.le.nx-1) then
              if (rad.lt.1.) then
                zs_v(i,j) = h/16.*(1+cos(pi*rad))**4
              else
                zs_v(i,j) = 0.
              end if
            end if

            if (xu.le.0.) then
              a = a1
            else
              a = a2
            end if

            rad = sqrt((xu/4./a)**2 + (yu/4./b)**2)

            if (j.le.ny-1) then
              if (rad.lt.1.) then
                zs_u(i,j) = h/16.*(1+cos(pi*rad))**4
              else
                zs_u(i,j) = 0.
              end if
            end if

          enddo

        enddo

      else if (type.eq.5) then

        do i=1,nx-1
          x  = x0 + dx*(i-1) + 0.5*dx - xc
          xu = x0 + dx*(i-1) - xc
          xv = x
        do j=1,ny-1
          if (x.le.0.) then
            a = a1
          else
            a = a2
          end if
          if (abs(x/4./a).lt.1.) then
            zs_p(i,j) = h/16.*(1.+cos(pi*x/4./a))**4
            zs_v(i,j) = zs_p(i,j)
          else
            zs_p(i,j) = 0.
            zs_v(i,j) = 0.
          end if
          if (xu.le.0.) then
            a = a1
          else
            a = a2
          end if
          if (abs(xu/4./a).lt.1.) then
            zs_u(i,j) = h/16.*(1.+cos(pi*xu/4./a))**4 
          else
            zs_u(i,j) = 0.
          end if 
        enddo
        enddo

        xu = x0 + dx*(nx-1) - xc
        if (xu.le.0.) then
          a = a1
        else
          a = a2
        end if
        do j=1,ny-1
          if (abs(xu/4./a).lt.1) then
            zs_u(nx,j) = h/16.*(1.+cos(pi*xu/4./a))**4 
          else
            zs_u(nx,j) = 0.
          end if
        enddo

        do i=1,nx-1
          xv = x0 + dx*(i-1) + 0.5*dx - xc
          if (xv.le.0.) then
            a = a1
          else
            a = a2
          end if
          if (abs(xv/4./a).lt.1.) then
            zs_v(i,ny) = h/16.*(1.+cos(pi*xv/4./a))**4
          else
            zs_v(i,ny) = 0.
          end if
        enddo

      else if (type.eq.6.) then

        do j=1,ny-1
          y  = y0 + dy*(j-1) + 0.5*dy - yc
          yv = y0 + dy*(j-1) - yc
          yu = y
        do i=1,nx-1
          x  = x0 + dx*(i-1) + 0.5*dx - xc 
          xu = x0 + dx*(i-1) - xc
          xv = x
          zs_p(i,j) = h*(a1**2+b1**2)/(x**2+y**2+a1**2+b1**2)
          zs_u(i,j) = h*(a1**2+b1**2)/(xu**2+yu**2+a1**2+b1**2)
          zs_v(i,j) = h*(a1**2+b1**2)/(xv**2+yv**2+a1**2+b1**2)
        enddo
        enddo

        xu = x0 + dx*(nx-1) - xc
        do j=1,ny-1
          yu = y0 + dy*(j-1) + 0.5*dy - yc
          zs_u(nx,j) = h*(a1**2+b1**2)/(xu**2+yu**2+a1**2+b1**2)
        enddo

        yv = y0 + dy*(ny-1) - yc
        do i=1,nx-1
          xv = x0 + (i-1)*dx + 0.5*dx - xc
          zs_v(i,ny) = h*(a1**2+b1**2)/(xv**2+yv**2+a1**2+b1**2)
        enddo

      else
        
 554    FORMAT('Error:  terrain type ',F5.1,' not yet implemented.')
        write(*,554) type
        stop
 
      end if
 
      return
      end

*------------------------------------------------------------------------

      subroutine write_2d(idcdf,fld,fld2,vname,iout,jout,nx,ny,
     >                            nxout,nyout,iu,iv)

      integer idcdf,iout,jout,nx,ny,nxout,nyout,iu,iv
      real fld(nx-1+iu,ny-1+iv),fld2(nxout-1+iu,nyout-1+iv)
      character*(*) vname

      integer istrt(2),ilen(2),idvar

      lgth = index(vname,' ') - 1
      if (lgth.eq.-1) lgth = len(vname)

      idvar = ncvid(idcdf,vname(1:lgth),ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting variable id.'
        stop
      end if

      do j=1,nyout-1+iv
      do i=1,nxout-1+iu
        fld2(i,j) = fld(iout+i-1,jout+j-1)
      enddo
      enddo

      istrt(1) = 1
      istrt(2) = 1
      
      ilen(1) = nxout-1+iu
      ilen(2) = nyout-1+iv

      call ncvpt(idcdf,idvar,istrt,ilen,fld2,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem writing 2d variable.'
        stop
      end if

      return
      end 

*------------------------------------------------------------------------

      subroutine write_3d(idcdf,fld,fld2,zs,gz,ztop,mix,vname,nx,ny,
     >                           nz,nzs,nxout,nyout,nzout,iout,jout,
     >                           is,iu,iv,iu2,iv2)

      integer idcdf,nx,ny,nz,nzs,nxout,nyout,nzout,iout,jout,is
      integer iu,iv,iu2,iv2
      real fld(nx-1+iu,ny-1+iv,nzs)
      real fld2(nxout-1+iu2,nyout-1+iv2,nzs)
      real zs(nx-1+iu2,ny-1+iv2),gz(nzs),ztop
      character*(*) vname

      integer istrt(4),ilen(4),idvar
      integer i,j,k,kptr,kplus,ii,im1,jj,jm1
      real z0,tfac,fp,fm

c     do j=1,ny-1+iv2
c       jj = min(j,ny-1)
c       jm1 = max(j-1,1)
c     do i=1,nx-1+iu2
c       ii = min(i,nx-1)
c       im1 = max(i-1,1)
c
c       z0 = zs(i,j)
c       tfac = (ztop-z0)/ztop
c       kptr = 1
c        
c       do k=1,nz-1+iw
c
c         z = z0 + tfac*gz(k)
c
c         kplus = 0
c589      if ((kplus.eq.0).and.(kptr.le.nz-1+iw)) then
c           if (gz(kptr).gt.z) then
c             kplus = kptr
c           else
c             kptr = kptr + 1
c           end if
c           goto 589
c         end if
c
c         if (kplus.eq.1) then
c           write(6,*)
c           write(6,*) 'Error:  problem with terrain interpolation.'
c           write(6,*)
c           stop
c         else if (kplus.eq.0) then
c           kplus = nz-1+iw
c         end if
c
c         if (iu2-iu.eq.1) then
c           fp = 0.5*(fld(ii,j,kplus)+fld(im1,j,kplus))
c           fm = 0.5*(fld(ii,j,kplus-1)+fld(im1,j,kplus-1))
c         else if (iv2-iv.eq.1) then
c           fp = 0.5*(fld(i,jj,kplus)+fld(i,jm1,kplus))
c           fm = 0.5*(fld(i,jj,kplus-1)+fld(i,jm1,kplus-1))
c         else 
c           fp = fld(i,j,kplus)
c           fm = fld(i,j,kplus-1)
c         end if
c
c         fld2(i,j,k) = fp + (fp-fm)/(gz(kplus)-gz(kplus-1))
c    >                            *(z-gz(kplus))
c
c       enddo
c
c     enddo
c     enddo

      do k=1,nzs
      do j=1,nyout-1+iv2
      do i=1,nxout-1+iu2

        ii = min(iout+i-1,nx-1)
        im1 = max(iout+i-2,1)
        jj = min(jout+j-1,ny-1)
        jm1 = max(jout+j-2,1)

        if (iu2-iu.eq.1) then
          fld2(i,j,k) = 0.5*(fld(ii,jout+j-1,k)+fld(im1,jout+j-1,k))
        else if (iv2-iv.eq.1) then
          fld2(i,j,k) = 0.5*(fld(iout+i-1,jj,k)+fld(iout+i-1,jm1,k))
        else
          fld2(i,j,k) = fld(iout+i-1,jout+j-1,k)
        end if

      enddo
      enddo
      enddo

      lgth = index(vname,' ') - 1
      if (lgth.eq.-1) lgth = len(vname)
      idvar = ncvid(idcdf,vname(1:lgth),ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem getting variable id.'
        stop
      end if

      istrt(1) = 1
      istrt(2) = 1
      if (mix.ne.0.) then
        if (is.eq.1) then
          istrt(3) = 1
        else
          istrt(3) = 3 + (is-1)*(nz-6) + 1
        end if
      else
        istrt(3) = (is-1)*nz + 1
      end if
      if (istrt(3).gt.nzout-1) istrt(3) = nzout
      istrt(4) = 1

      ilen(1) = nxout-1+iu2
      ilen(2) = nyout-1+iv2
      ilen(3) = nzs
      ilen(4) = 1

      call ncvpt(idcdf,idvar,istrt,ilen,fld2,ierr)
      if (ierr.ne.0) then
        write(6,*) 'Error:  problem writing variable.'
        stop
      end if

      return
      end

*---------------------------------------------------------------------------

      logical function mzero1(omga,N,fcor,m2)

      real omga,N,fcor,m2

      real dz
      common /dzblk/ dz

      real kz,pi
      logical requ

      pi = 2.*asin(1.)
      kz = sqrt(abs(m2))

      mzero1 = .false.
      
      if (requ(abs(omga),abs(fcor)).or.(omga.eq.0.)
     >      .or.((m2.gt.0.).and.(kz.gt.(2.*pi/3./dz)))
     >         .or.((m2.lt.0.).and.(kz.gt.(1./3./dz)))) then
        mzero1 = .true.
      end if

      return
      end

*---------------------------------------------------------------------------

      logical function mzero2(omga,N,fcor,m2)

      real omga,N,fcor,m2

      real dz
      common /dzblk/ dz

      real kz,pi
      logical requ

      pi = 2.*asin(1.)
      kz = sqrt(abs(m2))

      mzero2 = .false.
      
      if (requ(abs(omga),abs(fcor)).or.(omga.eq.0.)
     >      .or.((m2.gt.0.).and.(kz.gt.(2.*pi/5./dz)))
     >         .or.((m2.lt.0.).and.(kz.gt.(1./5./dz)))) then
        mzero2 = .true.
      end if

      return
      end

*---------------------------------------------------------------------------

      logical function requ(val1,val2)

      real val1,val2,eps,diff,big	

      parameter (eps=1.e-3)

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

*--------------------------------------------------------------------------

      subroutine set_sclrblk(asclr,nsclr)

      integer nsclr
      real asclr(nsclr)

      real fcor,N,U0,V0,ps,rhos,g
      common /sclrblk/ fcor,N,U0,V0,ps,rhos,g

      fcor = asclr(1)
      N = asclr(12)
      U0 = asclr(13)
      V0 = asclr(14)
      ps = asclr(15)
      rhos = asclr(17)
      g = 9.806

      return
      end

*----------------------------------------------------------------------
      
      subroutine invert(fld,ft,tran,wrk,cwrkx,cwrky,nx,ny,nz)

      integer nx,ny,nz
      real fld(nx-1,ny-1,nz-1),cwrkx(4*nx+15),cwrky(4*ny+15)
      complex ft(nx-1,ny-1,nz-1),tran(ny-1,nx-1),wrk(nx-1,ny-1)

      integer i,j,k

      do k=1,nz-1

        do i=1,nx-1
          do j=1,ny-1
            tran(j,i) = ft(i,j,k)
          enddo
          call cfftb(ny-1,tran(1,i),cwrky)
          do j=1,ny-1
            wrk(i,j) = tran(j,i)
          enddo
        enddo

        do j=1,ny-1
          call cfftb(nx-1,wrk(1,j),cwrkx)
          do i=1,nx-1
            fld(i,j,k) = real(wrk(i,j))
          enddo
        enddo

      enddo

      return
      end

*---------------------------------------------------------------------

      subroutine stretch_it(sgz,wgz,dz,dzbot,nz)

      integer nz
      real sgz(nz),wgz(nz),dz,dzbot

      real ztop,step,stretch,dzmax,mid,diff
      parameter (dzmax=700.)

      ztop = dz*(nz-1)
      step = 1.0
      stretch = 1.0

      if (dz.gt.dzbot) then
        k = 1
 32     if (abs(step).gt.(1.0e-10)) then
          step = step*0.5
          mid = stretch + step
          diff = zhgt(dzbot,mid,nz,dzmax) - ztop
          if (diff.le.0.0) stretch = mid
          k = k+1
          if (k.gt.50) then
            write(6,*) 'Problem calculating stretching...'
            stop
          end if
          goto 32
        end if
      end if

      wgz(1) = 0.0
      do k=2,nz
        wgz(k) = zhgt(dzbot,stretch,k,dzmax)
        sgz(k-1) = (wgz(k)+wgz(k-1))/2.
      enddo
      sgz(nz) = wgz(nz)
 
      return
      end
 
*----------------------------------------------------------------------------

      real function zhgt(dzbot,stretch,kz,dzmax)

      integer kz
      real dzbot,stretch,dzmax

      real sum,dznew

      sum = 0.
      do k=1,kz-1
        dznew = amin1(dzbot*(stretch**(k-1)),dzmax)
        sum = sum + dznew
      enddo

      zhgt = sum

      return
      end 

*--------------------------------------------------------------------------

      subroutine write_it(ctmp,var,nx,ny)
 
      integer nx,ny
      complex var(nx-1,ny-1)
      character*(*) ctmp

      lgth = index(ctmp,' ') - 1

      WRITE(*,*)
      do i=1,nx-1
        WRITE(*,*) 'i = ',i,'  ',ctmp(1:lgth),' = ',var(i,1)
      enddo
      WRITE(*,*)
 
      return
      end
      

