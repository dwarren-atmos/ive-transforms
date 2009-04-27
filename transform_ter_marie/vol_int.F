
#include "pointer.icl"

      ive_ptr function vol_int(varnm,i1,i2,j1,j2,k1,k2,ndims,dims,
     >               stag,min,max,missing,dunits,disunits,dim_names)
*
*  vol_int computes the volume integral, normalized by the horizontal
*    area of the total domain, of variable varnm.  The domain of
*    integration is the rectangular region enclosed by i1*dx to i2*dx
*    along the x-axis, j1*dy to j2*dy along the y-axis, and k1*dz to
*    k2*dz along the z-axis.  Here dz is the average grid spacing in
*    the vertical (without terrain).
*
      integer ndims,dims(4),i1,i2,j1,j2,k1,k2
      real stag(4),min(4),max(4),missing
      character*(*) varnm,dunits,disunits,dim_names(4)

      integer istrt,istop,jstrt,jstop,iu,iv,iw,k,strbeg,strend
      real ztop,dom_min(4),dom_max(4),dx,dy,zstrt,zstop,r_val
      character*(80) tmp,tmp2
      logical err,new,units_err,twod
      ive_ptr getmem,getvar,zs,sgz,wgz,var,var_t

      vol_int = 0

      call getrarr('plmin',dom_min,4,err)
      if (err) then
        write(*,*) 'Error:  could not find physical space domain ',
     >                   'minimums.'
        return
      end if
      call getrarr('plmax',dom_max,4,err)
      if (err) then
        write(*,*) 'Error:  could not find physical space domain ',
     >                    'maximums.'
        return
      end if

      if (dom_min(2).eq.dom_max(2)) then
        twod = .true.
      else
        twod = .false.
      end if

      wgz = getvar('wgz',ndims,dims,stag,min,max,missing,dunits,
     >                      disunits,dim_names,new)
      ztop = r_val(%val(wgz),dims(1))

      sgz = getvar('sgz',ndims,dims,stag,min,max,missing,dunits,
     >                      disunits,dim_names,new)

      zs = getvar('zbot_p',ndims,dims,stag,min,max,missing,dunits,
     >                      disunits,dim_names,new)

      dx = (max(1)-min(1))/(dims(1)-1)
      if (twod) then
        dy = 1.
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      var = getvar(varnm,ndims,dims,stag,min,max,missing,dunits,
     >                      disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot calculate field ',
     >                  varnm(strbeg(varnm):strend(varnm))
        return
      end if

      iu = 0
      iv = 0
      iw = 0
      if (dom_min(1).ge.min(1)) then
        iu = 1
      else if ((dom_min(2).ge.min(2)).and.(.not.twod)) then
        iv = 1
      else if (dom_min(3).ge.min(3)) then
        iw = 1
      end if

      dims(1) = dims(1) - iu
      dims(2) = dims(2) - iv
      dims(3) = dims(3) - iw
*
*  Interpolate variable varnm onto thermodynamic grid
*
      var_t = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (var_t.eq.0) then
        write(*,*) 'Error:  could not allocate space for var_t.'
        return
      end if

      call var_2_tpts(%val(var_t),%val(var),dims(1),dims(2),dims(3),
     >                        dims(4),iu,iv,iw)
*
*  Determine grid cells over which to integrate
*
      if (twod) then
        j1 = 0
        j2 = 1
      end if

      istrt = i1+1
      istop = i2
      jstrt = j1+1
      jstop = j2
*
*  Set upper and lower vertical limits in physical space
*
      zstrt = k1*ztop/dims(3)
      zstop = k2*ztop/dims(3)

 101  FORMAT(10x,A5,F10.2,2x,A5,F10.2)
    
      write(*,*)
      write(*,*) 'Calculating volume integral of field ',
     >                  varnm(strbeg(varnm):strend(varnm))
      write(*,*) 'over the domain defined by:  '
      write(*,101) 'x1 = ',i1*dx,'x2 = ',i2*dx
      write(*,101) 'y1 = ',j1*dy,'y2 = ',j2*dy
      write(*,101) 'z1 = ',zstrt,'z2 = ',zstop

      if ((zstrt.lt.0.).or.(zstop.gt.ztop).or.(istrt.lt.1).or.
     >         (istop.gt.dims(1)).or.(jstrt.lt.1).or.
     >         (jstop.gt.dims(2))) then
        write(*,*) 'Error:  specified limits exceed domain ',
     >                       'boundaries' 
        call freemem(var_t)
        return
      else if (((zstop-zstrt).le.(4.*ztop/dims(3))).or.
     >               (istrt.gt.istop).or.(jstrt.gt.jstop)) then
        write(*,*) 'Error:  improper specification of domain of ',
     >                        'integration.'
        call freemem(var_t)
        return
      end if 
        
      vol_int = getmem(dims(4))
      if (vol_int.eq.0) then
        write(*,*) 'Error:  could not allocate space for vol_int.'
        call freemem(var_t)
        return
      end if

      call vol_int_calc(%val(vol_int),%val(var_t),%val(zs),%val(sgz),
     >            %val(wgz),ztop,dx,dy,zstrt,zstop,istrt,istop,jstrt,
     >            jstop,dims(1),dims(2),dims(3),dims(4),err)

      if (err) then
        call freemem(var_t)
        call freemem(vol_int)
        vol_int = 0
        return
      end if

      do k=1,3
        dims(k) = 1
        min(k) = 0.
        max(k) = 0.
        stag(k) = 0.
        dim_names(k) = 'one'
      enddo
 
      tmp = dunits(strbeg(dunits):strend(dunits))//'*m'
      call units_smpfy(tmp,dunits,units_err)
      if (units_err) dunits = tmp

      tmp = disunits(strbeg(disunits):strend(disunits))//'*m'
      call units_smpfy(tmp,disunits,units_err)
      if (units_err) disunits = tmp

      call freemem(var_t)

      return
      end

************************************************************************

      subroutine var_2_tpts(var_t,var,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real var_t(nx,ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        var_t(i,j,k,n) = (var(i+iu,j+iv,k+iw,n) + var(i,j,k,n))/2.

      enddo
      enddo
      enddo
      enddo
  
      return
      end

*************************************************************************

      subroutine vol_int_calc(vol_int,var,zs,sgz,wgz,ztop,dx,dy,
     >                         zstrt,zstop,istrt,istop,jstrt,jstop,
     >                         nx,ny,nz,nt,err)

      integer istrt,istop,jstrt,jstop,nx,ny,nz,nt
      real vol_int(nt),var(nx,ny,nz,nt),zs(nx,ny),sgz(nz),wgz(nz+1)
      real ztop,dx,dy,zstrt,zstop

      integer i,j,k,n,ktop,kbot
      real area,dztop,dzbot,z1,z2,z3,z4,z5,z6,tfac,vbot,vtop
      logical err

      area = nx*dx*ny*dy

      do n=1,nt

        vol_int(n) = 0.

        do j = jstrt,jstop
        do i = istrt,istop

          call bot_top(zs(i,j),wgz,ztop,zstrt,zstop,kbot,ktop,nz)
          if (kbot.eq.0) then
            write(*,*) 'Error:  could not find lowest grid cell ',
     >                        'for volume integral.'
            err = .true.
            return  
          else if (kbot.eq.ktop) then
            write(*,*) 'Error:  volume of integration does not span ',
     >                         'one grid cell in depth.'
            err = .true.
            return  
          end if

          tfac = (ztop - zs(i,j))/ztop

          z1 = zs(i,j) + tfac*sgz(kbot)
          z2 = zs(i,j) + tfac*wgz(kbot+1)
          z3 = zs(i,j) + tfac*sgz(kbot+1)
          z4 = zs(i,j) + tfac*sgz(ktop-1)
          z5 = zs(i,j) + tfac*wgz(ktop)
          z6 = zs(i,j) + tfac*sgz(ktop) 
*
*  Find depths of top and bottom integration volume elements.  
*    Interpolate variable to the center of these elements.
*
          dzbot = z2 - amax1(zs(i,j),zstrt)

          if (zs(i,j).ge.zstrt) then
            vbot = var(i,j,kbot,n)
          else
            vbot = var(i,j,kbot,n) + (var(i,j,kbot+1,n)
     >                  - var(i,j,kbot,n)) 
     >             * ((z2-dzbot/2.) - z1)/(z3-z1)
          end if

          dztop = zstop - z5

          if (zstop.ge.ztop) then
            vtop = var(i,j,ktop,n)
          else
            vtop = var(i,j,ktop,n) + (var(i,j,ktop,n) 
     >                   - var(i,j,ktop-1,n))
     >             * ((zstop-dztop/2.) - z6)/(z6-z4)
          end if

          vol_int(n) = vol_int(n) + vbot*dx*dy*dzbot 
     >                    + vtop*dx*dy*dztop
*
*  Add on contributions of intervening volume elements
*
          do k = kbot+1,ktop-1
            vol_int(n) = vol_int(n) + var(i,j,k,n)*dx*dy
     >                      * tfac*(wgz(k+1)-wgz(k))
          enddo

        enddo
        enddo

        vol_int(n) = vol_int(n)/area

      enddo

      return
      end

**************************************************************************

      subroutine bot_top(zs,wgz,ztop,zstrt,zstop,kbot,ktop,nz)   
*
*  bot_top finds the grid elements intersected by the horizontal 
*    planes defined by zstrt and zstop
*
      integer kbot,ktop,nz
      real zs,wgz(nz+1),ztop,zstrt,zstop

      integer k
      real zgrd

      kbot = 0
      ktop = 0

      k = 1

 786  if ((kbot.eq.0).and.(k.lt.nz)) then
        zgrd = zs + (ztop-zs)*wgz(k+1)/ztop
        if (zgrd.gt.zstrt) then
          kbot = k
        else
          k = k+1
        end if
        goto 786
      end if

 897  if ((ktop.eq.0).and.(k.le.nz)) then
        zgrd = zs + (ztop-zs)*wgz(k+1)/ztop
        if (zgrd.ge.zstop) then
          ktop = k
        else
          k = k+1
        end if
        goto 897
      end if

      if (ktop.eq.0) ktop = nz

      return
      end 
        


          
      
