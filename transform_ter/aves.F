
#include "pointer.icl"

      ive_ptr function vol_ave(varnm,numvars,i1,i2,j1,j2,k1,k2,ndims,
     >          dims,stag,min,max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,k1,k2,numvars
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dunits,disunits,dim_names

      integer strbeg,strend,iarg(6),k
      real volume,r_val,ter_vol,rarg(1)
      character*(80) tmp,carg(10),intstr,mkfield
      logical err,new
      ive_ptr vol_int

      vol_ave = 0

      if (twod.eq.1) then
        j2 = 1
        j1 = 0
      end if

      volume = real(i2-i1)*dx*real(j2-j1)*dy*real(k2-k1)*dz
     >    - ter_vol(%val(zbot_p),dx,dy,i1,i2,j1,j2,k1*dz,
     >                              k2*dz,nx,ny)

      do k=1,numvars
        carg(k) = varnm(k)
      enddo
      iarg(1) = i1
      iarg(2) = i2
      iarg(3) = j1
      iarg(4) = j2
      iarg(5) = k1
      iarg(6) = k2

      intstr = mkfield('VOL_INT',carg,iarg,rarg,numvars,6,0,0)

      vol_ave = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names,new)
      if (vol_ave.eq.0) then
        write(*,*) 'Error:  cannot compute volume integral.'
        return
      end if

      call vol_ave_calc(%val(vol_ave),volume,dims(4))

      tmp = dunits(strbeg(dunits):strend(dunits))//'/m/m/m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = tmp

      tmp = disunits(strbeg(disunits):strend(disunits))//'/m/m/m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = tmp

      return
      end

*----------------------------------------------------------------------

      subroutine vol_ave_calc(vol_ave,volume,nt)

      integer nt,k
      real vol_ave(nt),volume

      do k=1,nt
        vol_ave(k) = vol_ave(k)/volume
      enddo
 
      return
      end
        
*---------------------------------------------------------------------

      ive_ptr function area_ave_x(varnm,numvars,j1,j2,k1,k2,i,ndims,
     >             dims,stag,min,max,miss,dunits,disunits,dim_names)

      integer ndims,dims(4),j1,j2,k1,k2,i,numvars
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dunits,disunits,dim_names(4)

      integer strbeg,strend,nx,ny,iarg(5),k
      real area,dx,dy,dz,r_val,ztop,rarg(1),ter_area_x 
      character*(80) tmp,carg(10),mkfield,intstr
      logical err,new
      ive_ptr zbot_u,wgz,getvar

      area_ave_x = 0

      zbot_u = getvar('zbot_u',ndims,dims,stag,min,max,miss,dunits,
     >             disunits,dim_names,new)
      if (zbot_u.eq.0) then
        write(*,*) 'Error:  cannot find field zbot_u.'
        return
      end if

      nx = dims(1)-1
      ny = dims(2)

      dx = (max(1)-min(1))/(dims(1)-1)
      if (max(2).eq.min(2)) then
        dy = 1.
        j2 = 1
        j1 = 0
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,
     >            disunits,dim_names,new)
      if (wgz.eq.0) then
        write(*,*) 'Error:  cannot find field wgz.'
        return
      end if

      ztop = r_val(%val(wgz),dims(1))

      dz = ztop/(dims(1)-1)

      area = real(j2-j1)*dy*real(k2-k1)*dz 
     >  - ter_area_x(%val(zbot_u),dy,j1,j2,i,k1*dz,k2*dz,nx,ny)

      do k=1,numvars
        carg(k) = varnm(k)
      enddo
      iarg(1) = j1
      iarg(2) = j2
      iarg(3) = k1
      iarg(4) = k2
      iarg(5) = i

      intstr = mkfield('AREA_INT_X',carg,iarg,rarg,numvars,5,0,0)
     
      area_ave_x = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new)
      if (area_ave_x.eq.0) then
        write(*,*) 'Error:  cannot compute area integral.'
        return
      end if

      call area_ave_x_calc(%val(area_ave_x),area,dims(4))
   
      tmp = dunits(strbeg(dunits):strend(dunits))//'/'
      tmp = tmp(strbeg(tmp):strend(tmp))//'m*m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = tmp

      tmp = disunits(strbeg(disunits):strend(disunits))//'/'
      tmp = tmp(strbeg(tmp):strend(tmp))//'m*m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = tmp

      return
      end

      subroutine area_ave_x_calc(area_ave,area,nt)

      integer nt,k
      real area_ave(nt),area
 
      do k=1,nt
        area_ave(k) = area_ave(k)/area
      enddo

      return
      end

*------------------------------------------------------------------------
      
      ive_ptr function area_ave_y(varnm,numvars,i1,i2,k1,k2,j,ndims,
     >             dims,stag,min,max,miss,dunits,disunits,dim_names)

      integer ndims,dims(4),i1,i2,k1,k2,j,numvars
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dunits,disunits,dim_names(4)

      integer strbeg,strend,nx,ny,iarg(5),k
      real area,dx,dy,dz,r_val,ztop,rarg(1),ter_area_y
      character*(80) tmp,carg(10),mkfield,intstr
      logical err,new
      ive_ptr zbot_v,wgz,getvar

      area_ave_y = 0

      zbot_v = getvar('zbot_v',ndims,dims,stag,min,max,miss,dunits,
     >             disunits,dim_names,new)
      if (zbot_v.eq.0) then
        write(*,*) 'Error:  cannot find field zbot_v.'
        return
      end if

      nx = dims(1)
      ny = dims(2)-1

      dx = (max(1)-min(1))/(dims(1)-1)
      if (max(2).eq.min(2)) then
        dy = 1.
        j = 1
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,
     >            disunits,dim_names,new)
      if (wgz.eq.0) then
        write(*,*) 'Error:  cannot find field wgz.'
        return
      end if

      ztop = r_val(%val(wgz),dims(1))

      dz = ztop/(dims(1)-1)

      area = real(i2-i1)*dx*real(k2-k1)*dz 
     >  - ter_area_y(%val(zbot_v),dx,i1,i2,j,k1*dz,k2*dz,nx,ny)

      do k=1,numvars
        carg(k) = varnm(k)
      enddo
      iarg(1) = i1
      iarg(2) = i2
      iarg(3) = k1
      iarg(4) = k2
      iarg(5) = j

      intstr = mkfield('AREA_INT_Y',carg,iarg,rarg,numvars,5,0,0)
     
      area_ave_y = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new)
      if (area_ave_y.eq.0) then
        write(*,*) 'Error:  cannot compute area integral.'
        return
      end if

      call area_ave_y_calc(%val(area_ave_y),area,dims(4))
   
      tmp = dunits(strbeg(dunits):strend(dunits))//'/'
      tmp = tmp(strbeg(tmp):strend(tmp))//'m*m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = tmp

      tmp = disunits(strbeg(disunits):strend(disunits))//'/'
      tmp = tmp(strbeg(tmp):strend(tmp))//'m*m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = tmp

      return
      end

      subroutine area_ave_y_calc(area_ave,area,nt)

      integer nt,k
      real area_ave(nt),area
 
      do k=1,nt
        area_ave(k) = area_ave(k)/area
      enddo

      return
      end

*------------------------------------------------------------------

      ive_ptr function area_ave_z(varnm,numvars,i1,i2,j1,j2,k,ndims,
     >              dims,stag,min,max,miss,dunits,disunits,dim_names)

      integer ndims,dims(4),i1,i2,j1,j2,k,numvars
      real stag(4),min(4),max(4),miss
      character*(*) varnm(numvars),dunits,disunits,dim_names(4)

      integer strbeg,strend,iarg(5),n,nx,ny
      real area,dx,dy,rarg(1),ztop,dz,r_val,ter_area_z
      character*(80) tmp,carg(10),mkfield,intstr
      logical err,new
      ive_ptr zbot_p,getvar,wgz

      area_ave_z = 0

      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,
     >                 disunits,dim_names,new)
      if (wgz.eq.0) then
        write(*,*) 'Error:  unable to access wgz.'
        return
      end if

      ztop = r_val(%val(wgz),dims(1))
      dz = ztop/(dims(1)-1)

      zbot_p = getvar('zbot_p',ndims,dims,stag,min,max,miss,dunits,
     >             disunits,dim_names,new)
      if (zbot_p.eq.0) then
        write(*,*) 'Error:  unable to access zbot_p.'
        return
      end if

      nx = dims(1)
      ny = dims(2)

      dx = (max(1)-min(1))/(dims(1)-1)
      if (max(2).eq.min(2)) then
        dy = 1.
        j2 = 1
        j1 = 0
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      area = real(i2-i1)*dx*real(j2-j1)*dy 
     >        - ter_area_z(%val(zbot_p),dx,dy,i1,i2,j1,j2,k*dz,nx,ny)

      do n=1,numvars
        carg(n) = varnm(n)
      enddo
      iarg(1) = i1
      iarg(2) = i2
      iarg(3) = j1
      iarg(4) = j2
      iarg(5) = k

      intstr = mkfield('AREA_INT_Z',carg,iarg,rarg,numvars,5,0,0)

      area_ave_z = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names,new)
      if (area_ave_z.eq.0) then
        write(*,*) 'Error:  cannot compute area integral.'
        return
      end if

      call area_ave_z_calc(%val(area_ave_z),area,dims(4))
 
      tmp = dunits(strbeg(dunits):strend(dunits))//'/'
      tmp = tmp(strbeg(tmp):strend(tmp))//'m*m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = tmp

      tmp = disunits(strbeg(disunits):strend(disunits))//'/'
      tmp = tmp(strbeg(tmp):strend(tmp))//'m*m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = tmp

      return
      end

      subroutine area_ave_z_calc(area_ave,area,nt)

      integer nt,k
      real area_ave(nt),area

      do k=1,nt
        area_ave(k) = area_ave(k)/area
      enddo
  
      return
      end
     
*-----------------------------------------------------------------------

      ive_ptr function area_ave_z_vp(varnm,numvars,i1,i2,j1,j2,dpth,
     >                     ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),numvars,i1,i2,j1,j2,dpth
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),varnm(numvars)

      integer strbeg,strend,iarg(5),n,nt
      real area,rarg(1)
      logical new,err
      character*(80) tmp,mkfield,intstr,carg(10)
      ive_ptr int,getmem

      area_ave_z_vp = 0

      nt = ntime 

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if

      if ((i1.ge.i2).or.(j1.ge.j2)) then
        write(6,*) 'Error:  improper specification of integration ',
     >                                 'domain.'
        return
      end if

      area = real(i2-i1)*dx*real(j2-j1)*dy

      do n=1,numvars
        carg(n)(1:len(carg(n))) = ' '
        carg(n) = varnm(n)(strbeg(varnm(n)):strend(varnm(n)))
      enddo

      iarg(1) = i1
      iarg(2) = i2
      iarg(3) = j1
      iarg(4) = j2
      iarg(5) = 0

      intstr = mkfield('AREA_INT_Z_VP',carg,iarg,rarg,numvars,5,0,0)

      int = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new)
      if (int.eq.0) then
        write(6,*) 'Error:  unable to calculate integrals.'
        return
      end if

      area_ave_z_vp = getmem((nz+1)*dims(4))
      if (area_ave_z_vp.eq.0) then
        write(6,*) 'Error:  unable to allocate space for ',
     >                           'area_ave_z_vp.'
        return
      end if

      call area_ave_z_vp_calc(%val(area_ave_z_vp),%val(int),area,
     >                                    nz,nt)

      tmp = dunits(strbeg(dunits):strend(dunits))//'/m/m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '

      tmp = disunits(strbeg(disunits):strend(disunits))//'/m/m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '

      return
      end

      subroutine area_ave_z_vp_calc(area_ave,area_int,area,nz,nt)

      integer nz,nt
      real area_ave(nz+1,nt),area_int(nz+1,nt),area

      integer k,n

      do n=1,nt
      do k=1,nz+1
        area_ave(k,n) = area_int(k,n)/area
      enddo
      enddo

      return
      end

