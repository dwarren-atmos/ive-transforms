
#include "pointer.icl"

      ive_ptr function diff_int_vp(dir,varnm,numvars,i1,i2,j1,j2,dpth,
     >                  ndims,dims,stag,min,max,miss,dunits,disunits,
     >                  dim_names)

#include "default.icl"

      integer ndims,dims(4),numvars,i1,i2,j1,j2,dpth
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),dir,varnm(numvars)

      integer strbeg,strend,iarg(5),n,nt,iw
      real area,rarg(1)
      logical new,err
      character*(80) tmp,mkfield,intstr,carg(10),ustr
      ive_ptr int1,int2,getmem

      nt = ntime

      if (mod(dpth,2).eq.0) then
        diff_int_vp = getmem((nz+1)*nt) 
        iw = 1
      else
        diff_int_vp = getmem(nz*nt)
        iw = 0
      end if

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if

      if ((i1.ge.i2).or.(j1.ge.j2)) then
        write(*,*) 'Error:  improper specification of integration ',
     >                                 'domain.'
        call freemem(diff_int_vp)
        return
      end if

      area = real(i2-i1)*dx*real(j2-j1)*dy

      do n=1,numvars
        carg(n)(1:len(carg(n))) = ' '
        carg(n) = varnm(n)(strbeg(varnm(n)):strend(varnm(n)))
      enddo

      if (dir(1:1).eq.'X') then

        iarg(1) = j1
        iarg(2) = j2
        iarg(3) = dpth
        iarg(4) = i1

        intstr = mkfield('AREA_INT_X_VP',carg,iarg,rarg,numvars,4,0,0)

        int1 = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
        if (int1.eq.0) then
          write(6,*) 'Error:  unable to calculate integrals.'
          call freemem(diff_int_vp)
          return
        end if

        iarg(4) = i2
      
        intstr = mkfield('AREA_INT_X_VP',carg,iarg,rarg,numvars,4,0,0)

        int2 = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
        if (int2.eq.0) then
          write(6,*) 'Error:  unable to calculate integrals.'
          call freemem(diff_int_vp)
          return
        end if

        call diff_int_vp_xy(%val(diff_int_vp),%val(int1),%val(int2),
     >                           %val(wgz),area,nz,nt,dpth,iw)

      else if (dir(1:1).eq.'Y') then
 
        iarg(1) = i1
        iarg(2) = i2
        iarg(3) = dpth
        iarg(4) = j1

        intstr = mkfield('AREA_INT_Y_VP',carg,iarg,rarg,numvars,4,0,0)

        int1 = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
        if (int1.eq.0) then
          write(6,*) 'Error:  unable to calculate integrals.'
          call freemem(diff_int_vp)
          return
        end if

        iarg(4) = j2
      
        intstr = mkfield('AREA_INT_Y_VP',carg,iarg,rarg,numvars,4,0,0)

        int2 = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
        if (int2.eq.0) then
          write(6,*) 'Error:  unable to calculate integrals.'
          call freemem(diff_int_vp)
          return
        end if

        call diff_int_vp_xy(%val(diff_int_vp),%val(int1),%val(int2),
     >                           %val(wgz),area,nz,nt,dpth,iw)
     
      else if (dir(1:1).eq.'Z') then

        iarg(1) = i1
        iarg(2) = i2
        iarg(3) = j1
        iarg(4) = j2
        iarg(5) = 0

        intstr = mkfield('AREA_INT_Z_VP',carg,iarg,rarg,numvars,5,0,0)

        int1 = getvar(intstr,ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
        if (int1.eq.0) then
          write(6,*) 'Error:  unable to calculate integrals.'
          call freemem(diff_int_vp)
          return
        end if

        call diff_int_vp_z(%val(diff_int_vp),%val(int1),%val(wgz),area,
     >                             nz,nt,dpth,iw)

      else

        write(6,*) 'Error:  unable to calculate difference for ',
     >               'direction ',dir(strbeg(dir):strend(dir))
        call freemem(diff_int_vp)
        return

      end if

      if (iw.eq.1) then
        dims(3) = nz+1
        stag(3) = 0.
        min(3) = wgz_min
        max(3) = wgz_max
        dim_names(3) = 'nzp1'
      else
        dims(3) = nz
        stag(3) = 0.5
        min(3) = sgz_min
        max(3) = sgz_max
        dim_names(3) = 'nz'
      end if

      if (dpth.eq.0) then
        ustr = '/m/m'
      else
        ustr = '/m/m/m'
      end if

      tmp = dunits(strbeg(dunits):strend(dunits))//
     >                 ustr(strbeg(ustr):strend(ustr))
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '

      tmp = disunits(strbeg(disunits):strend(disunits))//
     >                 ustr(strbeg(ustr):strend(ustr))
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '
 
      return
      end

*-------------------------------------------------------------------------

      subroutine diff_int_vp_xy(diff_int,int1,int2,wgz,area,nz,nt,
     >                                       dpth,iw)

      integer nz,nt,dpth,iw
      real diff_int(nz+iw,nt),int1(nz+iw,nt),int2(nz+iw,nt),wgz(nz+1)
      real area

      integer n,k,ktop,kbot,del
      real vol,delz

      del = dpth/2

      do k=1,nz+iw

        if (dpth.eq.0) then
          kbot = 1
          ktop = k
        else if (iw.eq.1) then
          kbot = max(1,k-del)
          ktop = min(nz+1,k+del)
        else
          kbot = max(1,k-del)
          ktop = min(nz+1,k+1+del)
        end if

        if (dpth.eq.0) then
          delz = 1.
        else
          delz = wgz(ktop)-wgz(kbot)
        end if
        vol = area*delz

        do n=1,nt
          diff_int(k,n) = int1(k,n)/vol - int2(k,n)/vol
        enddo

      enddo

      return
      end

*--------------------------------------------------------------------------

      subroutine diff_int_vp_z(diff_int,area_int,wgz,area,nz,nt,
     >                                      dpth,iw)        

      integer nz,nt,dpth,iw
      real diff_int(nz+iw,nt),area_int(nz+1,nt),wgz(nz+1),area

      integer n,k,kbot,ktop,del
      real vol,delz

      del = dpth/2

      do k=1,nz+iw

        if (dpth.eq.0) then
          kbot = 1
          ktop = k
        else if (iw.eq.1) then
          kbot = max(1,k-del)
          ktop = min(nz+1,k+del)
        else
          kbot = max(1,k-del)
          ktop = min(nz+1,k+1+del)
        end if

        if (dpth.eq.0) then
          delz = 1.
        else
          delz = wgz(ktop)-wgz(kbot)
        end if 
        vol = area*delz

        do n=1,nt
          diff_int(k,n) = area_int(kbot,n)/vol - area_int(ktop,n)/vol
        enddo

      enddo

      return
      end


