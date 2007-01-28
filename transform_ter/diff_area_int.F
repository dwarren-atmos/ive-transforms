
#include "pointer.icl"

      ive_ptr function diff_area_int(dir,varnm,numvars,i1,i2,j1,j2,
     >             k1,k2,ndims,dims,stag,min,max,miss,dunits,
     >             disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer numvars,i1,i2,j1,j2,k1,k2    
      character*(*) dir,varnm(numvars)

      real volume,r_val,ter_vol
      integer strbeg,strend
      logical new,err
      character*(80) tmp
      ive_ptr int1,int2,getmem,area_int_x,area_int_y,area_int_z

      diff_area_int = 0

      if ((i1.ge.i2).or.(j1.ge.j2).or.(k1.ge.k2)) then
        write(*,*) 'Error:  improper specification of surfaces ',
     >                         'of integration.'
        return
      end if

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if

      volume = real(i2-i1)*dx*real(j2-j1)*dy*real(k2-k1)*dz
     >    - ter_vol(%val(zbot_p),dx,dy,i1,i2,j1,j2,k1*dz,k2*dz,nx,ny)

      if (dir(1:1).eq.'X') then
        int1 = area_int_x(varnm,numvars,j1,j2,k1,k2,i1,ndims,dims,
     >               stag,min,max,miss,dunits,disunits,dim_names)
        if (int1.eq.0) then
          write(*,*) 'Error:  could not compute area integral'
          return
        end if
        int2 = area_int_x(varnm,numvars,j1,j2,k1,k2,i2,ndims,dims,
     >               stag,min,max,miss,dunits,disunits,dim_names)
        if (int2.eq.0) then
          write(*,*) 'Error:  could not compute area integral'
          call freemem(int1)
          return
        end if
      else if (dir(1:1).eq.'Y') then
        int1 = area_int_y(varnm,numvars,i1,i2,k1,k2,j1,ndims,dims,
     >               stag,min,max,miss,dunits,disunits,dim_names)
        if (int1.eq.0) then
          write(*,*) 'Error:  could not compute area integral'
          return
        end if
        int2 = area_int_y(varnm,numvars,i1,i2,k1,k2,j2,ndims,dims,
     >               stag,min,max,miss,dunits,disunits,dim_names)
        if (int2.eq.0) then
          write(*,*) 'Error:  could not compute area integral'
          call freemem(int1)
          return
        end if
      else if (dir(1:1).eq.'Z') then
        int1 = area_int_z(varnm,numvars,i1,i2,j1,j2,k1,ndims,dims,
     >               stag,min,max,miss,dunits,disunits,dim_names)
        if (int1.eq.0) then
          write(*,*) 'Error:  could not compute area integral'
          return
        end if
        int2 = area_int_z(varnm,numvars,i1,i2,j1,j2,k2,ndims,dims,
     >               stag,min,max,miss,dunits,disunits,dim_names)
        if (int2.eq.0) then
          write(*,*) 'Error:  could not compute area integral'
          call freemem(int1)
          return
        end if
      else
        write(*,*) 'Error:  don''t know how to compute area ',
     >         'integrals for direction ',
     >             dir(strbeg(dir):strend(dir))
        return
      end if

      diff_area_int = getmem(dims(4))
      if (diff_area_int.eq.0) then
        write(*,*) 'Error:  cannot allocate space for diff_area_int.'
        call freemem(int1)
        call freemem(int2)
      end if

      call diff_area_int_calc(%val(diff_area_int),%val(int1),
     >                 %val(int2),volume,dims(4))

      tmp = dunits(strbeg(dunits):strend(dunits))//'/m/m/m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '
     
      tmp = disunits(strbeg(disunits):strend(disunits))//'/m/m/m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '

      call freemem(int1)
      call freemem(int2)

      return
      end

*-----------------------------------------------------------------

      subroutine diff_area_int_calc(diff_area_int,int1,int2,
     >                     volume,nt)

      integer nt
      real diff_area_int(nt),int1(nt),int2(nt),volume

      integer n

      do n=1,nt
        diff_area_int(n) = int1(n)/volume - int2(n)/volume
      enddo

      return
      end

*------------------------------------------------------------------

      ive_ptr function zf_top(i1,i2,j1,j2,k,ndims,dims,stag,min,
     >              max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer i1,i2,j1,j2,k
        
      integer iarg(5),strbeg,strend
      real r_val,volume,rarg(1),ter_vol
      logical new,err
      character*(80) carg(5),mkfield,zintstr,tmp
      
      zf_top = 0 

      if (twod.eq.1) then
        j2 = 1
        j1 = 0
      end if

      volume = real(i2-i1)*dx*real(j2-j1)*dy*real(k)*dz
     >    - ter_vol(%val(zbot_p),dx,dy,i1,i2,j1,j2,0.,k*dz,nx,ny)

      carg(1) = 'RHO'
      if (fcor.ne.0.) then
        carg(2) = 'UABS'
      else
        carg(2) = 'U'
      end if
      carg(3) = 'W'
      iarg(1) = i1
      iarg(2) = i2
      iarg(3) = j1
      iarg(4) = j2
      iarg(5) = k

      zintstr = mkfield('AREA_INT_Z',carg,iarg,rarg,3,5,0,0)

      zf_top = getvar(zintstr,ndims,dims,stag,min,max,miss,dunits,
     >            disunits,dim_names,new)
      if (zf_top.eq.0) then 
        write(*,*) 'Error:  unable to calculate area integral.'
        return
      end if

      call zf_top_calc(%val(zf_top),volume,dims(4))

      tmp = dunits(strbeg(dunits):strend(dunits))//'/m/m/m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '

      tmp = disunits(strbeg(disunits):strend(disunits))//'/m/m/m' 
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '

      return
      end 

*----------------------------------------------------------------------

      subroutine zf_top_calc(zf,volume,nt)

      integer nt
      real zf(nt),volume

      integer n

      do n=1,nt
        zf(n) = - zf(n)/volume
      enddo

      return
      end

*---------------------------------------------------------------------

      ive_ptr function drg_bot(i1,i2,j1,j2,k,ndims,dims,stag,min,max,
     >           miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),i1,i2,j1,j2,k
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer iarg(6),strbeg,strend
      real r_val,volume,rarg(1),ter_vol
      logical new,err
      character*(80) carg(1),mkfield,drgstr,tmp

      drg_bot = 0
      
      if (twod.eq.1) then
        j2 = 1
        j1 = 0
      end if

      volume = real(i2-i1)*dx*real(j2-j1)*dy*real(k)*dz
     >    - ter_vol(%val(zbot_p),dx,dy,i1,i2,j1,j2,0.,k*dz,nx,ny)

      carg(1) = 'X'
      iarg(1) = i1
      iarg(2) = i2
      iarg(3) = j1
      iarg(4) = j2
      iarg(5) = 0
      iarg(6) = k

      drgstr = mkfield('DRAG',carg,iarg,rarg,1,6,0,0)

      drg_bot = getvar(drgstr,ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (drg_bot.eq.0) then
        write(*,*) 'Error:  cannot calculate drag.'
        return
      end if

      call drg_bot_calc(%val(drg_bot),volume,dims(4))

      tmp = dunits(strbeg(dunits):strend(dunits))//'/m/m/m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '
 
      tmp = disunits(strbeg(disunits):strend(disunits))//'/m/m/m' 
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '

      return
      end

*--------------------------------------------------------------------

      subroutine drg_bot_calc(drg,volume,nt)

      integer nt
      real drg(nt),volume
     
      integer k

      do k=1,nt
        drg(k) = drg(k)/volume
      enddo

      return
      end 

