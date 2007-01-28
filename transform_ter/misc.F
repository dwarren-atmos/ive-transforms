
#include "pointer.icl"

      ive_ptr function ru_pr(itime,ndims,dims,stag,min,max,miss,
     >                dunits,disunits,dim_names)

      integer ndims,dims(4),itime
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      character*(80) com1,com2,carg(5),mkfield
      integer iarg(5),strbeg,strend
      real rarg(5)
      logical new
      ive_ptr ru_0,ru,getvar

      ru_pr = 0

      carg(1) = 'U'
      iarg(1) = itime

      com1 = mkfield('REF',carg,iarg,rarg,1,1,0,0)

      carg(1) = com1(strbeg(com1):strend(com1))
      iarg(1) = itime

      com2 = mkfield('PER_VOL_REF',carg,iarg,rarg,1,1,0,0)

      carg(1) = 'PER_VOL[U]'
      carg(2) = com2(strbeg(com2):strend(com2))

      com1 = mkfield('PERT',carg,iarg,rarg,2,0,0,0)

      ru_pr = getvar(com1,ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
      if (ru_pr.eq.0) then
        write(*,*) 'Error:  unable to calculate ru_pr.'
        return
      end if
  
      return
      end 

*-----------------------------------------------------------------------

      ive_ptr function baseline(whch,val,ndims,dims,stag,min,max,miss,
     >                dunits,disunits,dim_names)

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss,val
      character*(*) dunits,disunits,dim_names(4),whch

      integer strbeg,strend,k,iwhch
      real del,r_val
      logical err,new
      ive_ptr getmem,dummy,getvar

      common/timeblk/time,ntime
      integer ntime
      ive_ptr time

      baseline = 0
 
      if (whch(1:1).eq.'X') then
        iwhch = 1
        dummy = getvar('zbot_u',ndims,dims,stag,min,max,miss,dunits,
     >                             disunits,dim_names,new)
      else if (whch(1:1).eq.'Y') then
        iwhch = 2
        dummy = getvar('zbot_v',ndims,dims,stag,min,max,miss,dunits,
     >                             disunits,dim_names,new)
      else if (whch(1:1).eq.'Z') then
        iwhch = 3
        dummy = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,
     >                             disunits,dim_names,new)
        dims(3) = dims(1)
        min(3) = r_val(%val(dummy),1)
        max(3) = r_val(%val(dummy),dims(1))
        stag(3) = 0.
        miss = 0.007
        dim_names(3) = 'nzp1'
      else 
        iwhch = 4
        dims(4) = ntime
        min(4) = r_val(%val(time),1)
        max(4) = r_val(%val(time),ntime)
        stag(4) = 0.
        miss = 0.007 
        dim_names(4) = 'time'
      end if

      ndims = 4

      if ((dummy.eq.0).and.(iwhch.ne.4)) then
        write(*,*) 'Error:  unable to find appropriate dummy field.'
        return
      end if

      dunits(1:len(dunits)) = ' '
      disunits(1:len(disunits)) = ' ' 
      miss = 0.007
 
      do k=1,4
        if (k.ne.iwhch) then
          dims(k) = 1
          stag(k) = 0.
          min(k) = 0.
          max(k) = 0.
          dim_names(k) = 'one'
        end if
      enddo

      baseline = getmem(dims(iwhch))
      if (baseline.eq.0) then
        write(*,*) 'Error:  unable to allocate space for line.'
        return
      end if

      do k=1,dims(iwhch)
        call s_val(%val(baseline),k,val)
      enddo

      return
      end 

*-------------------------------------------------------------------------

      ive_ptr function volume(i1,i2,j1,j2,k1,k2,ndims,dims,stag,min,
     >                    max,miss,dunits,disunits,dim_names)

      integer i1,i2,j1,j2,k1,k2,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer twod,k
      real dx,dy,dz,ter_vol,dom_min(4),dom_max(4),rvol,r_val,z1,z2
      logical new,err
      ive_ptr getvar,zs,getmem,wgz

      volume = 0
     
      call getrarr('plmin',dom_min,4,err)
      if (err) then
        write(*,*) 'Error:  unable to access physical domain ',
     >                          'minimums.'
        return
      end if
      call getrarr('plmax',dom_max,4,err)
      if (err) then
        write(*,*) 'Error:  unable to access physical domain ',
     >                          'maximums.'
        return
      end if

      if (dom_min(2).eq.dom_max(2)) then
        twod = 1
      else
        twod = 0
      end if

      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,
     >                 disunits,dim_names,new)
      if (wgz.eq.0) then
        write(*,*) 'Error:  unable to access wgz.'
        return
      end if

      dz = r_val(%val(wgz),dims(1))/(dims(1)-1) 
      z1 = real(k1)*dz
      z2 = real(k2)*dz

      zs = getvar('zbot_p',ndims,dims,stag,min,max,miss,dunits,
     >                 disunits,dim_names,new)
      if (zs.eq.0) then
        write(*,*) 'Error:  cannot access terrain array zbot_p.'
        return
      end if

      dx = (max(1)-min(1))/(dims(1)-1)
      if (twod.eq.1) then
        dy = 1.
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      volume = getmem(1)

      rvol = real(i2-i1)*dx*real(j2-j1)*dy*(z2-z1)
     >   - ter_vol(%val(zs),dx,dy,i1,i2,j1,j2,z1,z2,dims(1),dims(2)) 

      call s_val(%val(volume),1,rvol)

      ndims = 4
      do k=1,4
        dims(k) = 1
        min(k) = 0.
        max(k) = 0.
        stag(k) = 0.
        dim_names(k) = 'one'
      enddo
      dunits = 'm*m*m'
      disunits = 'm*m*m'

      return
      end

*------------------------------------------------------------------------- 

      real function ter_vol(zbot_p,dx,dy,i1,i2,j1,j2,z1,z2,nx,ny)

      integer i1,i2,j1,j2,nx,ny
      real zbot_p(nx,ny),dx,dy,z1,z2

      integer istrt,istop,jstrt,jstop,i,j
      real sum

      istrt = i1+1
      istop = i2
      jstrt = j1+1
      jstop = j2

      if ((istrt.gt.istop).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error: improper domain specification in ',
     >                        'ter_vol.'
        write(*,*) 'Stopping run.'
        stop 
      end if

      sum = 0.

      do i=istrt,istop
      do j=jstrt,jstop
        if (zbot_p(i,j).gt.z1) then
          sum = sum + amin1(zbot_p(i,j),z2) - z1
        end if
      enddo
      enddo

      ter_vol = sum*dx*dy

      return
      end

*---------------------------------------------------------------------

      real function ter_area_x(zbot_u,dy,j1,j2,i1,z1,z2,nx,ny)

      integer nx,ny,j1,j2,i1
      real zbot_u(nx+1,ny),dy,z1,z2

      integer ifix,jstrt,jstop,j
      real sum

      ifix = i1+1
      jstrt = j1+1
      jstop = j2

      sum = 0.

      do j=jstrt,jstop
        if (zbot_u(ifix,j).gt.z1) then
          sum = sum + amin1(zbot_u(ifix,j),z2) - z1
        end if
      enddo

      ter_area_x = sum*dy

      return
      end

*----------------------------------------------------------------------

      real function ter_area_y(zbot_v,dx,i1,i2,j1,z1,z2,nx,ny)

      integer nx,ny,i1,i2,j1
      real zbot_v(nx,ny+1),dx,z1,z2

      integer jfix,istrt,istop,i
      real sum

      jfix = j1+1
      istrt = i1+1
      istop = i2

      sum = 0.

      do i=istrt,istop
        if (zbot_v(i,jfix).gt.z1) then
          sum = sum + amin1(zbot_v(i,jfix),z2) - z1
        end if
      enddo

      ter_area_y = sum*dx

      return
      end 

*-------------------------------------------------------------------------

      real function ter_area_z(zs,dx,dy,i1,i2,j1,j2,z,nx,ny)

      integer i1,i2,j1,j2,nx,ny
      real zs(nx,ny),dx,dy,z

      integer istrt,istop,jstrt,jstop,i,j
      real sum

      istrt = i1+1
      istop = i2
      jstrt = j1+1
      jstop = j2

      sum = 0.

      do j=jstrt,jstop
      do i=istrt,istop
        if (zs(i,j).ge.z) then
          sum = sum + dx*dy
        end if
      enddo
      enddo

      ter_area_z = sum

      return
      end 

*----------------------------------------------------------------------------

      ive_ptr function uabs(ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names)  

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      logical new
      ive_ptr getmem,eta,u

      uabs = 0

      eta = getvar('eta',ndims,dims,stag,min,max,miss,dunits,
     >               disunits,dim_names,new)
      if (eta.eq.0) then
        write(*,*) 'Error:  unable to access field eta.'
        return
      end if

      uabs = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (uabs.eq.0) then
        write(*,*) 'Error:  unable to allocate space for uabs.'
        return
      end if

      call uabs_I(%val(uabs),%val(eta),fcor,dims(1),dims(2),dims(3),
     >                                  dims(4))

      call driver('free = eta',0)

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                   dim_names,new)
      if (u.eq.0) then
        write(*,*) 'Error:  unable to access zonal velocity.'
        return
      end if

      dims(1) = dims(1) - 1

      call uabs_II(%val(uabs),%val(u),dims(1),dims(2),dims(3),dims(4))

      stag(1) = 0.5
      min(1) = plmin(1) + 0.5*dx
      max(1) = plmax(1) - 0.5*dx
      dim_names(1) = 'nx'

      return
      end

*-------------------------------------------------------------------------

      subroutine uabs_I(uabs,eta,fcor,nx,ny,nz,nt)

      integer  nx,ny,nz,nt
      real uabs(nx,ny,nz,nt),eta(nx,ny,nz,nt),fcor

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        uabs(i,j,k,n) = - fcor*eta(i,j,k,n)

      enddo
      enddo
      enddo
      enddo

      return
      end

*-------------------------------------------------------------------------

      subroutine uabs_II(uabs,u,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real uabs(nx,ny,nz,nt),u(nx+1,ny,nz,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        uabs(i,j,k,n) = 0.5*(u(i,j,k,n)+u(i+1,j,k,n)) + uabs(i,j,k,n)

      enddo
      enddo
      enddo
      enddo

      return
      end

*--------------------------------------------------------------------------

      ive_ptr function quot(fld1,fld2,ndims,dims,stag,min,max,miss,
     >            dunits,disunits,dim_names)

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),fld1,fld2

      character*(80) dunits1,disunits1,dunits2,disunits2,tmp

      integer strbeg,strend,twod,iu1,iv1,iw1,iu2,iv2,iw2
      real dom_min(4),dom_max(4),dx,dy,r_val
      logical new,err
      ive_ptr getvar,getmem,var1,var2,zs,sgz

      quot = 0

      call getrarr('plmin',dom_min,4,err)
      if (err) then
        write(*,*) 'Error:  unable to access physical domain ',
     >                          'minimums.'
        return
      end if
      call getrarr('plmax',dom_max,4,err)
      if (err) then
        write(*,*) 'Error:  unable to access physical domain ',
     >                          'maximums.'
        return
      end if

      if (dom_min(2).eq.dom_max(2)) then
        twod = 1
      else
        twod = 0
      end if

      zs = getvar('zbot_p',ndims,dims,stag,min,max,miss,dunits,
     >                 disunits,dim_names,new)
      if (zs.eq.0) then
        write(*,*) 'Error:  cannot access terrain array zbot_p.'
        return
      end if

      dx = (max(1)-min(1))/(dims(1)-1)
      if (twod.eq.1) then
        dy = 1.
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      sgz = getvar('sgz',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      if (sgz.eq.0) then
        write(*,*) 'Error:  unable to access sgz.'
        return
      end if

      var1 = getvar(fld1,ndims,dims,stag,min,max,miss,dunits1,
     >                           disunits1,dim_names,new)
      if (var1.eq.0) then
        write(*,*) 'Error:  unable to access field ',
     >                fld1(strbeg(fld1):strend(fld1))
        return
      end if

      iu1 = 0
      iv1 = 0
      iw1 = 0
      if (dom_min(1).ge.min(1)) then
        iu1 = 1
      else if ((dom_min(2).ge.min(2)).and.(twod.ne.1)) then
        iv1 = 1
      else if (dom_min(3).ge.min(3)) then
        iw1 = 1
      end if

      var2 = getvar(fld2,ndims,dims,stag,min,max,miss,dunits2,
     >                        disunits2,dim_names,new)
      if (var2.eq.0) then
        write(*,*) 'Error:  unable to access field ',
     >                        fld2(strbeg(fld2):strend(fld2))
        return
      end if

      iu2 = 0
      iv2 = 0
      iw2 = 0
      if (dom_min(1).ge.min(1)) then
        iu2 = 1
      else if ((dom_min(2).ge.min(2)).and.(twod.ne.1)) then
        iv2 = 1
      else if (dom_min(3).ge.min(3)) then
        iw2 = 1
      end if

      dims(1) = dims(1) - iu2
      dims(2) = dims(2) - iv2
      dims(3) = dims(3) - iw2

      quot = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (quot.eq.0) then
        write(*,*) 'Error:  cannot allocate memory for quot.'
        return
      end if

      miss = 8675309
      call quot_calc(%val(quot),%val(var1),%val(var2),miss,dims(1),
     >              dims(2),dims(3),dims(4),iu1,iv1,iw1,iu2,iv2,iw2)

      min(1) = min(1) + real(iu2)*dx/2.
      max(1) = max(1) - real(iu2)*dx/2.
      min(2) = min(2) + real(iv2)*dy/2.
      max(2) = max(2) - real(iv2)*dy/2.
      min(3) = r_val(%val(sgz),1)
      max(3) = r_val(%val(sgz),dims(3))
      stag(1) = 0.5
      stag(2) = 0.5
      stag(3) = 0.5
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'
      dim_names(3) = 'nz'

      tmp = dunits1(strbeg(dunits1):strend(dunits1))//'/('
      tmp = tmp(strbeg(tmp):strend(tmp))//
     >             dunits2(strbeg(dunits2):strend(dunits2))
      tmp = tmp(strbeg(tmp):strend(tmp))//')' 
      call units_smpfy(tmp,dunits,err)
      if (err) dunits(1:len(dunits)) = ' '

      tmp = disunits1(strbeg(disunits1):strend(disunits1))//'/('
      tmp = tmp(strbeg(tmp):strend(tmp))//
     >          disunits2(strbeg(disunits2):strend(disunits2))
      tmp = tmp(strbeg(tmp):strend(tmp))//')'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits(1:len(disunits)) = ' '      

      return
      end 
       
*------------------------------------------------------------------------

      subroutine quot_calc(quot,var1,var2,miss,nx,ny,nz,nt,
     >                  iu1,iv1,iw1,iu2,iv2,iw2)

      integer nx,ny,nz,nt,iu1,iv1,iw1,iu2,iv2,iw2
      real quot(nx,ny,nz,nt),var1(nx+iu1,ny+iv1,nz+iw1,nt)
      real var2(nx+iu2,ny+iv2,nz+iw2,nt),miss

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        if (0.5*(var2(i,j,k,n)+var2(i+iu2,j+iv2,k+iw2,n)).ne.0.) then
          quot(i,j,k,n) = 
     >          0.5*(var1(i,j,k,n)+var1(i+iu1,j+iv1,k+iw1,n))
     >        / (0.5*(var2(i,j,k,n)+var2(i+iu2,j+iv2,k+iw2,n)))
        else
          quot(i,j,k,n) = miss
        end if

      enddo
      enddo
      enddo
      enddo

      return
      end 

*---------------------------------------------------------------------

      ive_ptr function prod(varnm,numvars,ndims,dims,stag,min,max,
     >                   miss,dunits,disunits,dim_names)

      integer ndims,dims(4),numvars
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),varnm(numvars)

      integer strbeg,strend,twod,iu,iv,iw,ijk,k
      real dom_min(4),dom_max(4),dx,dy,r_val
      logical new,err
      character*(80) tmp,vunits,vdisunits
      ive_ptr getvar,getmem,var,zs,sgz

      prod = 0

      call getrarr('plmin',dom_min,4,err)
      if (err) then
        write(*,*) 'Error:  unable to access physical domain ',
     >                          'minimums.'
        return
      end if
      call getrarr('plmax',dom_max,4,err)
      if (err) then
        write(*,*) 'Error:  unable to access physical domain ',
     >                          'maximums.'
        return
      end if

      if (dom_min(2).eq.dom_max(2)) then
        twod = 1
      else
        twod = 0
      end if

      zs = getvar('zbot_p',ndims,dims,stag,min,max,miss,dunits,
     >                 disunits,dim_names,new)
      if (zs.eq.0) then
        write(*,*) 'Error:  cannot access terrain array zbot_p.'
        return
      end if

      dx = (max(1)-min(1))/(dims(1)-1)
      if (twod.eq.1) then
        dy = 1.
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      sgz = getvar('sgz',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      if (sgz.eq.0) then
        write(*,*) 'Error:  unable to access sgz.'
        return
      end if

      var = getvar(varnm(1),ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  unable to access field ',
     >               varnm(1)(strbeg(varnm(1)):strend(varnm(1)))
        return
      end if

      vunits = dunits
      vdisunits = disunits

      iu = 0
      iv = 0
      iw = 0
      if (dom_min(1).ge.min(1)) then
        iu = 1
      else if ((dom_min(2).ge.min(2)).and.(twod.ne.1)) then
        iv = 1
      else if (dom_min(3).ge.min(3)) then
        iw = 1
      end if

      dims(1) = dims(1) - iu
      dims(2) = dims(2) - iv
      dims(3) = dims(3) - iw

      prod = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (prod.eq.0) then
        write(*,*) 'Error:  unable to allocate memory for prod.'
        return
      end if

      do ijk=1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(prod),ijk,1.)
      enddo

      call var_2_tpts2(%val(prod),%val(var),dims(1),dims(2),dims(3),
     >             dims(4),iu,iv,iw)

      k = 2

 976  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                 disunits,dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  unable to access field ',
     >             varnm(2)(strbeg(varnm(2)):strend(varnm(2)))
          call freemem(prod)
          return
        end if

        tmp = vunits(strbeg(vunits):strend(vunits))//'*'
        vunits = tmp(strbeg(tmp):strend(tmp))//
     >                  dunits(strbeg(dunits):strend(dunits))
        tmp = vdisunits(strbeg(vdisunits):strend(vdisunits))//'*'
        vdisunits = tmp(strbeg(tmp):strend(tmp))//
     >                disunits(strbeg(disunits):strend(disunits))

        iu = 0
        iv = 0
        iw = 0
        if (dom_min(1).ge.min(1)) then
          iu = 1
        else if ((dom_min(2).ge.min(2)).and.(twod.ne.1)) then
          iv = 1
        else if (dom_min(3).ge.min(3)) then
          iw = 1
        end if

        dims(1) = dims(1) - iu
        dims(2) = dims(2) - iv
        dims(3) = dims(3) - iw

        call var_2_tpts2(%val(prod),%val(var),dims(1),dims(2),dims(3),
     >                     dims(4),iu,iv,iw)

        k = k+1

        goto 976

      end if

      min(1) = min(1) + real(iu)*dx/2.
      max(1) = max(1) - real(iu)*dx/2.
      min(2) = min(2) + real(iv)*dy/2.
      max(2) = max(2) - real(iv)*dy/2.
      min(3) = r_val(%val(sgz),1)
      max(3) = r_val(%val(sgz),dims(3))
      stag(3) = 0.5
      stag(2) = 0.5
      stag(1) = 0.5
      dim_names(1) = 'nx'
      dim_names(3) = 'nz'
      dim_names(2) = 'ny'

      call units_smpfy(vunits,dunits,err)
      if (err) dunits(1:len(dunits)) = ' '

      call units_smpfy(vdisunits,disunits,err)
      if (err) disunits(1:len(disunits)) = ' '

      return
      end

*-----------------------------------------------------------------------

      subroutine var_2_tpts2(var_t,var,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real var_t(nx,ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        var_t(i,j,k,n) = var_t(i,j,k,n) *
     >          0.5*(var(i,j,k,n)+var(i+iu,j+iv,k+iw,n))

      enddo
      enddo
      enddo
      enddo

      return
      end

*-------------------------------------------------------------------------

      subroutine var_2_tpts3(var_t,var,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real var_t(nx,ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        var_t(i,j,k,n) = var_t(i,j,k,n) + 
     >          0.5*(var(i,j,k,n)+var(i+iu,j+iv,k+iw,n))

      enddo
      enddo
      enddo
      enddo

      return
      end 

*-------------------------------------------------------------------------

      subroutine zs_2_upts(zs_u,zs_p,nx,ny)

      integer nx,ny
      real zs_u(nx+1,ny),zs_p(nx,ny)

      integer i,j

      do j=1,ny
        zs_u(1,j) = zs_p(1,j)
        do i=2,nx
          zs_u(i,j) = 0.5*(zs_p(i,j)+zs_p(i-1,j))
        enddo
        zs_u(nx+1,j) = zs_p(nx,j)
      enddo

      return
      end

*--------------------------------------------------------------------------- 

      subroutine zs_2_vpts(zs_v,zs_p,nx,ny)

      integer nx,ny
      real zs_v(nx,ny+1),zs_p(nx,ny)

      integer i,j

      do i=1,nx
        zs_v(i,1) = zs_p(i,1)
        do j=2,ny
          zs_v(i,j) = 0.5*(zs_p(i,j)+zs_p(i,j-1))
        enddo
        zs_v(i,ny+1) = zs_p(i,ny)
      enddo

      return
      end

*---------------------------------------------------------------------

      ive_ptr function sum(varnm,numvars,ndims,dims,stag,min,max,
     >                   miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),numvars
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),varnm(numvars)

      integer strbeg,strend,iu,iv,iw,ijk,k
      real r_val
      logical new,err,requ
      character*(80) tmp,vunits,vdisunits
      ive_ptr getmem,var

      sum = 0

      var = getvar(varnm(1),ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  unable to access field ',
     >               varnm(1)(strbeg(varnm(1)):strend(varnm(1)))
        return
      end if

      vunits = dunits
      vdisunits = disunits

      if (abs(min(1)-plmin(1)).gt.dx) min(1) = min(1) + offset(1)
      if (abs(min(2)-plmin(2)).gt.dy) min(2) = min(2) + offset(2)

      iu = 0
      iv = 0
      iw = 0
      if (requ(min(1),plmin(1))) then
        iu = 1
      else if (requ(min(2),plmin(2)).and.(twod.ne.1)) then
        iv = 1
      else if (requ(min(3),plmin(3))) then
        iw = 1
      end if

      dims(1) = dims(1) - iu
      dims(2) = dims(2) - iv
      dims(3) = dims(3) - iw

      sum = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (sum.eq.0) then
        write(*,*) 'Error:  unable to allocate memory for sum.'
        return
      end if

      do ijk=1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(sum),ijk,0.)
      enddo

      call var_2_tpts3(%val(sum),%val(var),dims(1),dims(2),dims(3),
     >             dims(4),iu,iv,iw)

      k = 2

 976  if (k.le.numvars) then

        var = getvar(varnm(k),ndims,dims,stag,min,max,miss,dunits,
     >                 disunits,dim_names,new)
        if (var.eq.0) then
          write(*,*) 'Error:  unable to access field ',
     >             varnm(2)(strbeg(varnm(2)):strend(varnm(2)))
          call freemem(sum)
          return
        end if

        if (dunits(strbeg(dunits):strend(dunits)).ne.
     >            vunits(strbeg(vunits):strend(vunits))) then
          vunits = ' '
        end if
        if (disunits(strbeg(disunits):strend(disunits)).ne.
     >       vdisunits(strbeg(vdisunits):strend(vdisunits))) then
          vdisunits = ' '
        end if

        if (abs(min(1)-plmin(1)).gt.dx) min(1) = min(1) + offset(1)
        if (abs(min(2)-plmin(2)).gt.dy) min(2) = min(2) + offset(2)

        iu = 0
        iv = 0
        iw = 0
        if (requ(min(1),plmin(1))) then
          iu = 1
        else if (requ(min(2),plmin(2)).and.(twod.ne.1)) then
          iv = 1
        else if (requ(min(3),plmin(3))) then
          iw = 1
        end if

        dims(1) = dims(1) - iu
        dims(2) = dims(2) - iv
        dims(3) = dims(3) - iw

        call var_2_tpts3(%val(sum),%val(var),dims(1),dims(2),dims(3),
     >                     dims(4),iu,iv,iw)

        k = k+1

        goto 976

      end if

      min(1) = plmin(1) + 0.5*dx
      max(1) = plmax(1) - 0.5*dx
      if (twod.ne.1) then
        min(2) = plmin(2) + 0.5*dy
        max(2) = plmax(2) - 0.5*dy
      else
        min(2) = 0.
        max(2) = 0.
      end if
      min(3) = sgz_min
      max(3) = sgz_max
      stag(1) = 0.5
      stag(2) = 0.5
      stag(3) = 0.5
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'
      dim_names(3) = 'nz'

      call units_smpfy(vunits,dunits,err)
      if (err) dunits(1:len(dunits)) = ' '

      call units_smpfy(vdisunits,disunits,err)
      if (err) disunits(1:len(disunits)) = ' '

      return
      end

*----------------------------------------------------------------------

      ive_ptr function max_min(varnm,imax,isurf,isubdom,i1,i2,j1,j2,
     >            ndims,dims,stag,min,max,miss,dunits,disunits,
     >            dim_names)

#include "default.icl"

      integer imax,isurf,isubdom,i1,i2,j1,j2,k1,k2,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) varnm,dunits,disunits,dim_names(4)

      integer strbeg,strend,k,iw,iloc,jloc,kloc
      real r_val,dzbot,zstrt,zstop
      character*(3) ctmp
      logical new,requ
      ive_ptr getmem,var

      k1 = 0;
      k2 = 64;

      max_min = 0

      if (imax.eq.1) then
        ctmp = 'max'
      else
        ctmp = 'min'
      end if

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >               disunits,dim_names,new)

      dzbot = r_val(%val(wgz),2) - r_val(%val(wgz),1)
      zstrt = (k1)*dzbot;
      zstop = (k2)*dzbot

 86   FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)
 87   FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)
 88   FORMAT(10x,'z1 = ',F10.2,2x,'z2 = ',F10.2)

      if (var.eq.0) then
        write(*,*) 'Error:  unable to access variable ',
     >       varnm(strbeg(varnm):strend(varnm))
        return
      end if

      iw = 0
      if (requ(plmin(3),min(3))) iw = 1

      write(6,*) 
      write(6,*) 'Calculating ',ctmp(1:3),' value for field ',
     >             varnm(strbeg(varnm):strend(varnm))
      if (isurf.eq.1) then
        write(6,*) '     at the surface.  '
	else if (isurf.eq.2) then
	  write(6,*) '     at the lowest grid point.'
      end if
      if (isubdom.eq.0) then
	  write(6,*) '     over whole domain.  '
      else
	  write(6,*) '     over domain,  '
        write(*,86) plmin(1)+i1*dx, plmin(1)+i2*dx  
        write(*,87) plmin(2)+j1*dy, plmin(2)+j2*dy
c        write(*,88) zstrt, zstop
      end if
        write(6,*)

      max_min = getmem(dims(4))
      if (max_min.eq.0) then
        write(*,*) 'Error:  unable to allocate space for max_min.'
        return
      end if

c      call max_min_calc(%val(max_min),%val(var),imax,isurf,isubdom,iw,
c     >                 i1,i2,j1,j2,dims(1),dims(2),dims(3),dims(4))

      call max_min_calc(%val(max_min),%val(var),imax,isurf,isubdom,iw,
     >                 i1,i2,j1,j2,iloc,jloc,kloc,
     >                 dims(1),dims(2),dims(3),dims(4))

      do k=1,3
        dims(k) = 1
        stag(k) = 0.
        min(k) = 0.
        max(k) = 0.
        dim_names(k) = 'one'
      enddo

 90   FORMAT(10x,'x extreme = ',F10.2)
 91   FORMAT(10x,'y extreme = ',F10.2)
 92   FORMAT(10x,'z extreme = ',I10)

      if (dims(4).eq.1) then
        if (imax.eq.1) then
          write(6,*) '   Max Val = ',r_val(%val(max_min),1)
        else
          write(6,*) '   Min Val = ',r_val(%val(max_min),1)
        end if
        write(6,*)
        call freemem(max_min)
        max_min = 0

        write(*,90) plmin(1)+iloc*dx  
        write(*,91) plmin(2)+jloc*dy
        write(*,92)  kloc
	write(*,*)
	
      end if
        
      return
      end

*-------------------------------------------------------------------------

      subroutine max_min_calc(max_min,var,imax,isurf,isubdom,iw,
     >                         i1,i2,j1,j2,iloc,jloc,kloc,
     >                         nx,ny,nz,nt) 

      integer imax,isurf,isubdom,iw,i1,j1,i2,j2,nx,ny,nz,nt
      integer iloc,jloc,kloc
      integer numi1,numj1,numi2,numj2
      real max_min(nt),var(nx,ny,nz,nt)

      integer i,j,k,n,numz
      real extreme,val

      if (isurf.eq.1.or.isurf.eq.2) then
        numz = 1
      else
        numz = nz
      end if

      if (isubdom.eq.0) then
        numi1 = 1
      	numi2 = nx
	numj1 = 1
        numj2 = ny
      else
        numi1 = i1
      	numi2 = i2
	numj1 = j1
        numj2 = j2
      end if

      if (imax.eq.1) then

        do n=1,nt
          extreme = var(numi1,numj1,1,n)
          do k=1,numz
          do j=numj1,numj2
          do i=numi1,numi2
            if ((isurf.eq.1).and.(iw.eq.0)) then
              val = 0.5*(3.*var(i,j,1,n)-var(i,j,2,n))
            else
              val = var(i,j,k,n)
            end if 
            if (val.gt.extreme) then
	      extreme = val
	      iloc = i
	      jloc = j
	      kloc = k
            end if
          enddo
          enddo
          enddo
          max_min(n) = extreme 
        enddo

      else

        do n=1,nt
          extreme = var(numi1,numj1,1,n)
          do k=1,numz
          do j=numj1,numj2
          do i=numi1,numi2
            if ((isurf.eq.1).and.(iw.eq.0)) then
              val = 0.5*(3.*var(i,j,1,n)-var(i,j,2,n))
            else
              val = var(i,j,k,n)
            end if
            if (val.lt.extreme) then
	      extreme = val
	      iloc = i
	      jloc = j
	      kloc = k
	    end if
          enddo
          enddo
          enddo
          max_min(n) = extreme
        enddo

      end if

      return
      end

*------------------------------------------------------------------------

      ive_ptr function symm_swtch(fld,whch,ndims,dims,stag,min,max,
     >                   miss,dunits,disunits,dim_names)

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),fld,whch

      integer strbeg,strend
      logical new
      ive_ptr var,getmem,getvar

      symm_swtch = 0

      if ((whch(1:1).ne.'X').and.(whch(1:1).ne.'Y')) then
        write(*,*) 'Error:  do not know how to do symmetry switch ',
     >                   'for direction ',whch(1:1)
        return
      end if

      var = getvar(fld,ndims,dims,stag,min,max,miss,dunits,disunits,
     >                             dim_names,new) 
      if (var.eq.0) then
        write(*,*) 'Error:  unable to access field ',
     >                  fld(strbeg(fld):strend(fld))
        return
      end if

      symm_swtch = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (symm_swtch.eq.0) then
        write(*,*) 'Error:  unable to allocate space for ',
     >                        'symm_swtch.'
        return
      end if 

      call symm_swtch_calc(%val(symm_swtch),%val(var),whch,dims(1),
     >               dims(2),dims(3),dims(4))

      return
      end
 
*------------------------------------------------------------------------

      subroutine symm_swtch_calc(symm_swtch,var,whch,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real symm_swtch(nx,ny,nz,nt),var(nx,ny,nz,nt)
      character*(*) whch

      integer i,j,k,n,nx2,ny2

      if (whch(1:1).eq.'X') then

        nx2 = nx/2

        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx2
          symm_swtch(i,j,k,n) = var(nx-i+1,j,k,n)
          symm_swtch(nx-i+1,j,k,n) = var(i,j,k,n)
        enddo
        enddo
        enddo
        enddo

      else 

        ny2 = ny/2

        do n=1,nt
        do k=1,nz
        do j=1,ny2
        do i=1,nx
          symm_swtch(i,j,k,n) = var(i,ny-j+1,k,n)
          symm_swtch(i,ny-j+1,k,n) = var(i,j,k,n)
        enddo
        enddo
        enddo
        enddo

      end if

      return
      end

*----------------------------------------------------------------------

      ive_ptr function scalar(val,vunits,ndims,dims,stag,min,max,miss,
     >               dunits,disunits,dim_names)

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss,val
      character*(*) dunits,disunits,dim_names(4),vunits

      integer strbeg,strend,q
      ive_ptr getmem

      scalar = 0

      scalar = getmem(1)
      if (scalar.eq.0) then
        write(6,*) 'Error:  could not allocate space for scalar.'
        return
      end if

      call s_val(%val(scalar),1,val)

      ndims = 4
      miss = -123455.
      dunits(1:len(dunits)) = ' '
      disunits(1:len(disunits)) = ' '
      dunits = vunits(strbeg(vunits):strend(vunits))
      disunits = vunits(strbeg(vunits):strend(vunits))

      do q=1,4
        dims(q) = 1
        stag(q) = 0.
        min(q) = 0.
        max(q) = 0.
        dim_names(q) = 'one'
      enddo

      return
      end

*---------------------------------------------------------------------------

      ive_ptr function dist(dir,ndims,dims,stag,min,max,miss,dunits,
     >                             disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
      character*(1) dir

      integer nt
      real r_val
      ive_ptr getmem

      if ((dir.ne.'X').and.(dir.ne.'Y').and.(dir.ne.'Z')) then
        write(6,*) 'Error:  do not know how to compute dist for ',
     >                        'direction ',dir
        dist = 0
        return
      end if

      nt = ntime

      dist = getmem(nx*ny*nz*nt)
      if (dist.eq.0) then
        write(6,*) 'Error:  unable to allocate space for dist.'
        return
      end if
 
      call dist_calc(%val(dist),%val(zbot_p),%val(sgz),ztop,plmin(1),
     >                   plmin(2),dx,dy,dir,nx,ny,nz,nt)

      ndims = 4
      dims(1) = nx
      dims(2) = ny
      dims(3) = nz
      dims(4) = nt
      stag(1) = 0.5
      stag(2) = 0.5
      stag(3) = 0.5
      stag(4) = 0.0
      min(1) = plmin(1) + 0.5*dx
      min(2) = plmin(2) + 0.5*dy
      min(3) = sgz_min
      min(4) = r_val(%val(time),1)
      max(1) = plmax(1) - 0.5*dx
      max(2) = plmax(2) - 0.5*dy
      max(3) = sgz_max
      max(4) = r_val(%val(time),nt)
      miss = -0.007
      dunits = 'm'
      disunits = 'm'
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'
      dim_names(3) = 'nz'
      dim_names(4) = 'time'

      return
      end

      subroutine dist_calc(dist,zs_p,sgz,ztop,x0,y0,dx,dy,dir,
     >                              nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real dist(nx,ny,nz,nt),zs_p(nx,ny),sgz(nz),ztop,x0,y0,dx,dy
      character*(1) dir

      integer i,j,k,n
      real z0,tfac

      if (dir.eq.'X') then

        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx
          dist(i,j,k,n) = x0 + (i-1)*dx
        enddo
        enddo
        enddo
        enddo

      else if (dir.eq.'Y') then

        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx
          dist(i,j,k,n) = y0 + (j-1)*dy
        enddo
        enddo
        enddo
        enddo

      else

        do n=1,nt
        do j=1,ny 
        do i=1,nx
          z0 = zs_p(i,j)
          tfac = (ztop-z0)/ztop
          do k=1,nz
            dist(i,j,k,n) = z0 + tfac*sgz(k)
          enddo
        enddo
        enddo
        enddo

      end if

      return
      end 

*---------------------------------------------------------------------------

      ive_ptr function upstrm_hgt(ndims,dims,stag,min,max,miss,
     >                        dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      ive_ptr getmem,zeta
      logical new

      upstrm_hgt = 0

      zeta = getvar('zeta',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new)
      if (zeta.eq.0) then
        write(6,*) 'Error:  could not access zeta.'
        return
      end if
 
      upstrm_hgt = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (zeta.eq.0) then
        write(6,*) 'Error:  unable to allocate space for zeta.'
        return
      end if

      call upstrm_hgt_calc(%val(upstrm_hgt),%val(zeta),%val(sgz),
     >           %val(zbot_p),ztop,dims(1),dims(2),dims(3),dims(4))

      return
      end
*--------------------------------------------------------------------

      subroutine upstrm_hgt_calc(upstrm_hgt,zeta,sgz,zs_p,ztop,
     1		              nx,ny,nz,nt)


      integer nx,ny,nz,nt
      real upstrm_hgt(nx,ny,nz,nt),zeta(nx,ny,nz,nt),sgz(nz)
	real zs_p(nx,ny),zs,tfac,ztop

      integer i,j,k,n

	print *,nx,ny

      do n=1,nt
        do k=1,nz
          do j=1,ny
            do i=1,nx
               zs = zs_p(i,j)
               tfac = (ztop-zs)/ztop
               upstrm_hgt(i,j,k,n) = zs + tfac*sgz(k) - zeta(i,j,k,n)
            enddo
          enddo
        enddo
      enddo

      return 
      end 

*-----------------------------------------------------------------------

      ive_ptr function net_disp(dir,ndims,dims,stag,min,max,miss,
     >                        dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
      character*(1) dir

      integer nt,iu,iv
      logical new
      ive_ptr getmem,disp_pr,v0

      net_disp = 0

      nt = ntime
      iu = 0
      iv = 0

      if (dir.eq.'X') then

        v0 = getvar('u0',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                           dim_names,new)
        if (v0.eq.0) then
          write(6,*) 'Error:  could not access u0.'
          return
        end if

        disp_pr = getvar('xi',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
        if (disp_pr.eq.0) then
          write(6,*) 'Error:  could not access xi.'
          return
        end if

        iu = 1

      else if (dir.eq.'Y') then

        v0 = getvar('v0',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                           dim_names,new)
        if (v0.eq.0) then
          write(6,*) 'Error:  could not access v0.'
          return
        end if

        disp_pr = getvar('eta',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
        if (disp_pr.eq.0) then
          write(6,*) 'Error:  could not access eta.'
          return
        end if

        iv = 1

      else
      
        write(6,*) 'Error:  do not know how to compute net_disp for ',
     >                        'direction ',dir
        return

      end if
 
      net_disp = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (net_disp.eq.0) then
        write(6,*) 'Error:  unable to allocate space for net_disp.'
        return
      end if

      call net_disp_calc(%val(net_disp),%val(disp_pr),%val(v0),
     >           %val(time),dims(1),dims(2),dims(3),dims(4),iu,iv)

      return
      end

*--------------------------------------------------------------------

      subroutine net_disp_calc(disp,disp_pr,v0,time,nx,ny,nz,nt,iu,iv)

      integer nx,ny,nz,nt,iu,iv
      real disp(nx,ny,nz,nt),disp_pr(nx,ny,nz,nt),v0(nx+iu,ny+iv,nz)
      real time(nt)

      integer i,j,k,n
      real disp0,vel,pi,delt,stime

      vel = v0(1,1,1)
      pi = 2.*asin(1.)

      do n=1,nt

        if (time(n).lt.0.) then
          delt = time(n) - time(1)
          stime = -time(1)
          disp0 = 0.5*vel*delt - 0.5*vel*sin(delt*pi/stime)*stime/pi
        else
          delt = time(n)
          stime = -time(1)
          disp0 = 0.5*vel*stime + vel*delt
        end if

        do k=1,nz
        do j=1,ny
        do i=1,nx

          disp(i,j,k,n) = disp0 + disp_pr(i,j,k,n)

        enddo
        enddo
        enddo

      enddo

      return 
      end 

*-------------------------------------------------------------------------

      ive_ptr function wndspeed(ipert,ndims,dims,stag,min,max,miss,
     >                                dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),ipert
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,iu,iv,iw,size
      logical new
      ive_ptr getmem,u,v,w,u0,v0

      wndspeed = getmem(nx*ny*nz*ntime)
      if (wndspeed.eq.0) then
        write(6,*) 'Error:  could not allocate space for wndspeed.'
        return
      end if

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                               dim_names,new)
      if (u.eq.0) then
        write(6,*) 'Error:  could not access u.'
        return
      end if
      u0 = getvar('u0',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                               dim_names,new)
      if (u0.eq.0) then
        write(6,*) 'Error:  could not access u0.'
        return
      end if
      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                               dim_names,new)
      if (v.eq.0) then
        write(6,*) 'Error:  could not access v.'
        return
      end if
      v0 = getvar('v0',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                               dim_names,new)
      if (v0.eq.0) then
        write(6,*) 'Error:  could not access v0.'
        return
      end if
      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                               dim_names,new)
      if (v.eq.0) then
        write(6,*) 'Error:  could not access w.'
        return
      end if

      call wndspeed_calc(%val(wndspeed),%val(u),%val(v),%val(w),
     >             %val(u0),%val(v0),nx,ny,nz,ntime,ipert,twod)

      dims(3) = dims(3)-1
      stag(3) = 0.5
      min(3) = sgz_min
      max(3) = sgz_max
      dim_names(3) = 'nz'

      return
      end

*--------------------------------------------------------------------------

      subroutine wndspeed_calc(speed,u,v,w,u0,v0,nx,ny,nz,nt,
     >                                   ipert,twod)

      integer nx,ny,nz,nt,ipert,twod
      real speed(nx,ny,nz,nt),u(nx+1,ny,nz,nt),v(nx,ny+1-twod,nz,nt)
      real w(nx,ny,nz+1,nt),u0(nx+1,ny,nz),v0(nx,ny+1-twod,nz)

      real tmp1,tmp2,tmp3
      integer i,j,k,n

      if (ipert.ne.0) then

        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx

          tmp1 = 0.5*(u(i,j,k,n)+u(i+1,j,k,n)
     >                       - u0(i,j,k)-u0(i+1,j,k))
          tmp2 = 0.5*(v(i,j,k,n)+v(i,j+1-twod,k,n)
     >                       - v0(i,j,k)-v0(i,j+1-twod,k))
          tmp3 = 0.5*(w(i,j,k,n)+w(i,j,k+1,n))

          speed(i,j,k,n) = sqrt(tmp1**2 + tmp2**2 + tmp3**2)

        enddo
        enddo
        enddo
        enddo

      else

        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx

          tmp1 = 0.5*(u(i,j,k,n)+u(i+1,j,k,n))
          tmp2 = 0.5*(v(i,j,k,n)+v(i,j+1-twod,k,n))
          tmp3 = 0.5*(w(i,j,k,n)+w(i,j,k+1,n))

          speed(i,j,k,n) = sqrt(tmp1**2 + tmp2**2 + tmp3**2)

        enddo
        enddo
        enddo
        enddo

      end if

      return
      end

*---------------------------------------------------------------------

      ive_ptr function maketime(varnm,ndims,dims,stag,min,max,miss,
     >                   dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) varnm,dunits,disunits,dim_names(4)

      integer strbeg,strend
      logical new
      ive_ptr getmem,var

      maketime = 0

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  unable to access field ',
     >               varnm(strbeg(varnm):strend(varnm))
        return
      else if (dims(4).ne.1) then
        write(6,*) 'Error:  field ',
     >                 varnm(strbeg(varnm):strend(varnm)),
     >                     ' already time-dependent.'
        return
      end if

      maketime = getmem(dims(1)*dims(2)*dims(3)*ntime)
      if (maketime.eq.0) then
        write(6,*) 'Error:  unable to allocate space for maketime.'
        return
      end if

      call maketime_calc(%val(maketime),%val(var),dims(1),dims(2),
     >                            dims(3),ntime)

      dims(4) = ntime

      return
      end

      subroutine maketime_calc(maketime,var,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real maketime(nx,ny,nz,nt),var(nx,ny,nz)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        maketime(i,j,k,n) = var(i,j,k)

      enddo
      enddo
      enddo
      enddo

      return
      end

*---------------------------------------------------------------------

      ive_ptr function ave(varnm,whch,n1,n2,ndims,dims,stag,min,max,
     >                           miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),n1,n2
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),varnm,whch

      integer strbeg,strend
      character*(80) new
      ive_ptr getmem,var

      ave = 0

      if ((whch(1:1).ne.'X').and.(whch(1:1).ne.'Y')
     >        .and.(whch(1:1).ne.'Z').and.(whch(1:1).ne.'T')) then
        write(6,*) 'Error:  do not know direction over which to ',
     >                            'compute average.'
        return
      end if

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  unable to access field ',
     >                      varnm(strbeg(varnm):strend(varnm))
        return
      end if

      if (whch(1:1).eq.'X') then

        if ((n1.lt.1).or.(n2.gt.dims(1))) then
          write(6,*) 'Error:  specified averaging domain exceeds ',
     >                           'model domain.'
          return
        else if (n1.ge.n2) then
          write(6,*) 'Error:  improper specification of averaging ',
     >                            'domain.'
          return
        end if

        ave = getmem(dims(2)*dims(3)*dims(4))
        if (ave.eq.0) then
          write(6,*) 'Error unable to allocate space for ave.'
          return
        end if

        call xave_calc(%val(ave),%val(var),miss,dims(1),dims(2),
     >                    dims(3),dims(4),n1,n2)

        dims(1) = 1
        stag(1) = 0. 
        min(1) = 0.
        max(1) = 0.
        dim_names(1) = 'one'

      else if (whch(1:1).eq.'Y') then

        if ((n1.lt.1).or.(n2.gt.dims(2))) then
          write(6,*) 'Error:  specified averaging domain exceeds ',
     >                           'model domain.'
          return
        else if (n1.ge.n2) then
          write(6,*) 'Error:  improper specification of averaging ',
     >                            'domain.'
          return
        end if

        ave = getmem(dims(1)*dims(3)*dims(4))
        if (ave.eq.0) then
          write(6,*) 'Error unable to allocate space for ave.'
          return
        end if

        call yave_calc(%val(ave),%val(var),miss,dims(1),dims(2),
     >                    dims(3),dims(4),n1,n2)

        dims(2) = 1
        stag(2) = 0. 
        min(2) = 0.
        max(2) = 0.
        dim_names(2) = 'one'

      else if (whch(1:1).eq.'Z') then

        if ((n1.lt.1).or.(n2.gt.dims(3))) then
          write(6,*) 'Error:  specified averaging domain exceeds ',
     >                           'model domain.'
          return
        else if (n1.ge.n2) then
          write(6,*) 'Error:  improper specification of averaging ',
     >                            'domain.'
          return
        end if

        ave = getmem(dims(1)*dims(2)*dims(4))
        if (ave.eq.0) then
          write(6,*) 'Error unable to allocate space for ave.'
          return
        end if

        call zave_calc(%val(ave),%val(var),miss,dims(1),dims(2),
     >                    dims(3),dims(4),n1,n2)

        dims(3) = 1
        stag(3) = 0. 
        min(3) = 0.
        max(3) = 0.
        dim_names(3) = 'one'

      else

        if ((n1.lt.1).or.(n2.gt.dims(4))) then
          write(6,*) 'Error:  specified averaging domain exceeds ',
     >                           'model domain.'
          return
        else if (n1.ge.n2) then
          write(6,*) 'Error:  improper specification of averaging ',
     >                            'domain.'
          return
        end if

        ave = getmem(dims(1)*dims(2)*dims(3))
        if (ave.eq.0) then
          write(6,*) 'Error unable to allocate space for ave.'
          return
        end if

        call tave_calc(%val(ave),%val(var),miss,dims(1),dims(2),
     >                    dims(3),dims(4),n1,n2)

        dims(4) = 1
        stag(4) = 0. 
        min(4) = 0.
        max(4) = 0.
        dim_names(4) = 'one'

      end if

      return
      end

*------------------------------------------------------------------------

      subroutine xave_calc(ave,var,miss,nx,ny,nz,nt,n1,n2)

      integer nx,ny,nz,nt,n1,n2
      real ave(ny,nz,nt),var(nx,ny,nz,nt),miss

      integer i,j,k,n,cnt
      real sum

      do n=1,nt
      do k=1,nz
      do j=1,ny

         sum = 0.
         cnt = 0

         do i=n1,n2
           if (var(i,j,k,n).ne.miss) then
             sum = sum + var(i,j,k,n)
             cnt = cnt + 1 
           end if
         enddo

         if (cnt.ne.0) ave(j,k,n) = sum/cnt

      enddo
      enddo
      enddo
 
      return
      end
      
      subroutine yave_calc(ave,var,miss,nx,ny,nz,nt,n1,n2)

      integer nx,ny,nz,nt,n1,n2
      real ave(nx,nz,nt),var(nx,ny,nz,nt),miss

      integer i,j,k,n,cnt
      real sum

      do n=1,nt
      do k=1,nz
      do i=1,nx

         sum = 0.
         cnt = 0

         do j=n1,n2
           if (var(i,j,k,n).ne.miss) then
             sum = sum + var(i,j,k,n)
             cnt = cnt + 1 
           end if
         enddo

         if (cnt.ne.0) ave(i,k,n) = sum/cnt

      enddo
      enddo
      enddo
 
      return
      end
      
      subroutine zave_calc(ave,var,miss,nx,ny,nz,nt,n1,n2)

      integer nx,ny,nz,nt,n1,n2
      real ave(nx,ny,nt),var(nx,ny,nz,nt),miss

      integer i,j,k,n,cnt
      real sum

      do n=1,nt
      do j=1,ny
      do i=1,nx

         sum = 0.
         cnt = 0

         do k=n1,n2
           if (var(i,j,k,n).ne.miss) then
             sum = sum + var(i,j,k,n)
             cnt = cnt + 1 
           end if
         enddo

         if (cnt.ne.0) ave(i,j,n) = sum/cnt

      enddo
      enddo
      enddo
 
      return
      end
      
      subroutine tave_calc(ave,var,miss,nx,ny,nz,nt,n1,n2)

      integer nx,ny,nz,nt,n1,n2
      real ave(nx,ny,nz),var(nx,ny,nz,nt),miss

      integer i,j,k,n,cnt
      real sum

      do k=1,nz
      do j=1,ny 
      do i=1,nx

         sum = 0.
         cnt = 0

         do n=n1,n2
           if (var(i,j,k,n).ne.miss) then
             sum = sum + var(i,j,k,n)
             cnt = cnt + 1 
           end if
         enddo

         if (cnt.ne.0) ave(i,j,k) = sum/cnt

      enddo
      enddo
      enddo
 
      return
      end

*------------------------------------------------------------------------

      ive_ptr function pbase(ndims,dims,stag,min,max,miss,dunits,
     >                              disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      ive_ptr getmem

      pbase = getmem(nx*ny*nz)
      if (pbase.eq.0) then
        write(6,*) 'Error:  cannot allocate memory for pbase.'
        return
      end if

      call pbase_calc(%val(pbase),%val(zbot_p),%val(sgz),ztop,
     >                               nx,ny,nz)

      ndims = 4
      dims(1) = nx
      dims(2) = ny
      dims(3) = nz
      dims(4) = 1
      stag(1) = 0.5
      stag(2) = 0.5
      stag(3) = 0.5
      stag(4) = 0.0
      min(1) = plmin(1) + 0.5*dx
      min(2) = plmin(2) + 0.5*dy
      min(3) = sgz_min
      min(4) = 0.
      max(1) = plmax(1) - 0.5*dx
      max(2) = plmax(2) - 0.5*dy
      max(3) = sgz_max
      max(4) = 0.
      miss = -777.
      dunits = 'N/m/m'
      disunits = 'N/m/m'
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'
      dim_names(3) = 'nz'
      dim_names(4) = 'one'

      return
      end

      subroutine pbase_calc(pbase,zs_p,sgz,ztop,nx,ny,nz)

      integer nx,ny,nz
      real pbase(nx,ny,nz),zs_p(nx,ny),sgz(nz),ztop

      integer i,j,k
      real z0,tfac,z

      real ps,rhos,g,N2
      parameter (ps = 100000., rhos = 1.2754, g = 9.806, N2 = 0.000144)

      do j=1,ny 
      do i=1,nx

        z0 = zs_p(i,j)
        tfac = (ztop-z0)/ztop

        do k=1,nz
          z = z0 + tfac*sgz(k)
          pbase(i,j,k) = ps - rhos*g*z + 0.5*rhos*N2*(z**2)
        enddo

      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      ive_ptr function r0b(ndims,dims,stag,min,max,miss,dunits,
     >                              disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      ive_ptr getmem

      r0b = getmem(nx*ny*nz)
      if (r0b.eq.0) then
        write(6,*) 'Error:  cannot allocate memory for r0b.'
        return
      end if

      call r0b_calc(%val(r0b),%val(zbot_p),%val(sgz),ztop,nx,ny,nz)

      ndims = 4
      dims(1) = nx
      dims(2) = ny
      dims(3) = nz
      dims(4) = 1
      stag(1) = 0.5
      stag(2) = 0.5
      stag(3) = 0.5
      stag(4) = 0.0
      min(1) = plmin(1) + 0.5*dx
      min(2) = plmin(2) + 0.5*dy
      min(3) = sgz_min
      min(4) = 0.
      max(1) = plmax(1) - 0.5*dx
      max(2) = plmax(2) - 0.5*dy
      max(3) = sgz_max
      max(4) = 0.
      miss = -777.
      dunits = 'kg/m/m/m'
      disunits = 'kg/m/m/m'
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'
      dim_names(3) = 'nz'
      dim_names(4) = 'one'

      return
      end

      subroutine r0b_calc(r0b,zs_p,sgz,ztop,nx,ny,nz)

      integer nx,ny,nz
      real r0b(nx,ny,nz),zs_p(nx,ny),sgz(nz),ztop

      integer i,j,k
      real z0,tfac,z

      real rhos,g,N2
      parameter (rhos = 1.27613, g = 9.806, N2 = 0.000144)

      do j=1,ny 
      do i=1,nx

        z0 = zs_p(i,j)
        tfac = (ztop-z0)/ztop

        do k=1,nz
          z = z0 + tfac*sgz(k)
          r0b(i,j,k) = rhos*(1.-N2/g*z) 
        enddo

      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      ive_ptr function rb(ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      character*(3) ctmp
      logical new
      ive_ptr getmem,b,b0,r0b

      rb = 0

      ctmp = 'R0B'
      r0b = getvar(ctmp,ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new)
      if (r0b.eq.0) then
        write(6,*) 'Error:  cannot access variable r0b.'
        return
      end if

      ctmp = 'b0'
      b0 = getvar(ctmp,ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new)
      if (b0.eq.0) then
        write(6,*) 'Error:  cannot access variable b0.'
        return
      end if

      ctmp = 'b'
      b = getvar(ctmp,ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new)
      if (b.eq.0) then
        write(6,*) 'Error:  cannot access variable b.'
        return
      end if

      rb = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (rb.eq.0) then
        write(6,*) 'Error:  cannot allocate space for rb.'
        return
      end if

      call rb_calc(%val(rb),%val(b),%val(b0),%val(r0b),dims(1),
     >                      dims(2),dims(3),dims(4))

      dunits = 'kg/m/m/m'
      disunits = 'kg/m/m/m'

      return
      end

      subroutine rb_calc(rb,b,b0,r0b,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real rb(nx,ny,nz,nt),b(nx,ny,nz,nt),b0(nx,ny,nz)
      real r0b(nx,ny,nz)

      integer i,j,k,n
      real rhos,g
      parameter (rhos = 1.27613, g = 9.806) 

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        rb(i,j,k,n) = - rhos/g*(b(i,j,k,n)-b0(i,j,k)) + r0b(i,j,k)

      enddo
      enddo
      enddo
      enddo

      return
      end

*--------------------------------------------------------------------------

      ive_ptr function surf(fld,ndims,dims,stag,min,max,miss,dunits,
     >                             disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) fld,dunits,disunits,dim_names(4)

      integer nt,strbeg,strend
      logical new,requ
      ive_ptr var,getmem,zs,gz

      surf = 0

      nt = ntime

      var = getvar(fld,ndims,dims,stag,min,max,miss,dunits,disunits,
     >                        dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  unable to access field ',
     >                fld(strbeg(fld):strend(fld))
        return
      end if

      surf = getmem(dims(1)*dims(2)*nt)

      zs = zbot_p
      gz = sgz
      if (requ(min(1),plmin(1))) then
        zs = zbot_u
      else if ((requ(min(2),plmin(2))).and.(twod.ne.1)) then
        zs = zbot_v
      else if (requ(min(3),plmin(3))) then 
        gz = wgz
      end if

      call surf_calc(%val(surf),%val(var),%val(zs),%val(gz),ztop,
     >                             dims(1),dims(2),dims(3),nt)

      dims(3) = 1
      stag(3) = 0.0
      min(3) = wgz_min
      max(3) = wgz_min
      dim_names(3) = 'one'

      return
      end

*---------------------------------------------------------------------------

      subroutine surf_calc(surf,var,zs,gz,ztop,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real surf(nx,ny,nt),var(nx,ny,nz,nt),zs(nx,ny),gz(nz),ztop

      integer i,j,n
      real z0,z1,z2,z3,vs

      do j=1,ny
      do i=1,nx

        z3 = zs(i,j) + (ztop-zs(i,j))/ztop*gz(3)
        z2 = zs(i,j) + (ztop-zs(i,j))/ztop*gz(2)
        z1 = zs(i,j) + (ztop-zs(i,j))/ztop*gz(1)
        z0 = zs(i,j)

        do n=1,nt

          surf(i,j,n) = (z0-z1)*(z0-z2)/(z3-z1)/(z3-z2)*var(i,j,3,n)
     >             + (z0-z1)*(z0-z3)/(z2-z1)/(z2-z3)*var(i,j,2,n)
     >          + (z0-z2)*(z0-z3)/(z1-z2)/(z1-z3)*var(i,j,1,n)

        enddo

      enddo
      enddo

      return
      end

*--------------------------------------------------------------------------

      ive_ptr function upert(ndims,dims,stag,min,max,miss,dunits,
     >                             disunits,dim_names)
   
#include "default.icl" 

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer nt
      logical new
      ive_ptr u,u0,getmem

      upert = 0

      nt = ntime

      u0 = getvar('u0',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                             dim_names,new)
      if (u0.eq.0) then
        write(6,*) 'Error:  unable to access u0.'
        return
      end if

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                              dim_names,new)
      if (u.eq.0) then
        write(6,*) 'Error:  unable to access u.'
        return
      end if

      upert = getmem((nx+1)*ny*nz*nt)
      if (upert.eq.0) then
        write(6,*) 'Error:  unable to allocate space for upert.'
        return
      end if

      call upert_calc(%val(upert),%val(u),%val(u0),%val(time),
     >                          nx,ny,nz,nt)

      return
      end

*--------------------------------------------------------------------------

      subroutine upert_calc(upert,u,u0,time,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real upert(nx+1,ny,nz,nt),u(nx+1,ny,nz,nt),u0(nx+1,ny,nz)
      real time(nt)

      integer i,j,k,n
      real stup_time,pi,ubase,delt

      pi = 2.*asin(1.)
      stup_time = abs(time(1))

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx+1

        delt = time(n) - time(1)
        if (delt.lt.stup_time) then
          ubase = 0.5*u0(i,j,k)*(1.-cos(delt/stup_time*pi))
        else
          ubase = u0(i,j,k)
        end if

        upert(i,j,k,n) = u(i,j,k,n) - ubase

      enddo
      enddo
      enddo
      enddo

      return
      end 

*--------------------------------------------------------------------------

      ive_ptr function nrmlz(fld,t0,ndims,dims,stag,min,max,miss,
     >                          dunits,disunits,dim_names)

#include "default.icl"

      real t0
      character*(*) fld 

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
    
      integer strbeg,strend,it0,ifnd
      real r_val
      logical new 
      ive_ptr getmem,var

      nrmlz = 0

      var = getvar(fld,ndims,dims,stag,min,max,miss,dunits,disunits,
     >                        dim_names,new) 
      if (var.eq.0) then
        write(6,*) 'Error:  unable to access field ',
     >                fld(strbeg(fld):strend(fld))       
        return
      end if

      it0 = 1
      ifnd = 0
 444  if ((it0.le.ntime).and.(ifnd.eq.0)) then
        if (r_val(%val(time),it0).ge.t0) then
          ifnd = 1
        else 
          it0 = it0 + 1
        end if
        goto 444
      end if

      if (ifnd.eq.0) then
        write(6,*) 'Error:  cannot find specified time.'
        return
      else if (it0.gt.dims(4)) then
        write(6,*) 'Error:  specified time exceeds limits for ',
     >                    fld(strbeg(fld):strend(fld))
        return
      end if

      nrmlz = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (nrmlz.eq.0) then
        write(6,*) 'Error:  cannot allocate space for nrmlz.'
        return
      end if

      call nrmlz_calc(fld,%val(var),%val(nrmlz),%val(time),it0,
     >                 dims(1),dims(2),dims(3),dims(4),ntime)

      return
      end

*---------------------------------------------------------------------------

      subroutine nrmlz_calc(fld,var,nrmlz,time,it0,nx,ny,nz,nt,ntime)

      integer it0,nx,ny,nz,nt,ntime
      real var(nx,ny,nz,nt),nrmlz(nx,ny,nz,nt),time(ntime)
      character*(*) fld

      integer i,j,k,n,strbeg,strend
      real maxval

      maxval = 0.

      do n=it0,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx
        if (abs(var(i,j,k,n)).gt.maxval) then
          maxval = abs(var(i,j,k,n))  
        end if
      enddo
      enddo
      enddo
      enddo

      write(6,*)
      write(6,*) 'Max val for ',fld(strbeg(fld):strend(fld)),' = ',
     >                maxval
      write(6,*)

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx
        nrmlz(i,j,k,n) = var(i,j,k,n)/maxval
      enddo
      enddo
      enddo
      enddo

      return
      end

