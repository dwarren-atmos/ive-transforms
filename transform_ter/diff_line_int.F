
#include "pointer.icl"

      ive_ptr function diff_line_int(dir,varnm,numvars,i1,i2,j1,j2,
     >             dpth,ndims,dims,stag,min,max,miss,dunits,disunits,
     >             dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer numvars,i1,i2,j1,j2,dpth
      character*(*) dir,varnm(numvars)

      real area
      integer strbeg,strend
      character*(80) tmp
      logical new,err
      ive_ptr int1,int2,getmem,line_int_x,line_int_y

      diff_line_int = 0

      if ((i1.ge.i2).or.(j1.ge.j2)) then
        write(*,*) 'Error:  improper specification of surfaces ',
     >                         'of integration.'
        return
      end if

      if (twod.eq.1) then
        j1 = 0
        j2 = 1
      end if

      area = real(i2-i1)*dx*real(j2-j1)*dy

      if (dir(1:1).eq.'X') then
        int1 = line_int_x(varnm,numvars,j1,j2,i1,dpth,ndims,dims,
     >               stag,min,max,miss,dunits,disunits,dim_names)
        if (int1.eq.0) then
          write(*,*) 'Error:  could not compute line integral'
          return
        end if
        int2 = line_int_x(varnm,numvars,j1,j2,i2,dpth,ndims,dims,
     >               stag,min,max,miss,dunits,disunits,dim_names)
        if (int2.eq.0) then
          write(*,*) 'Error:  could not compute line integral'
          call freemem(int1)
          return
        end if
      else if (dir(1:1).eq.'Y') then
        int1 = line_int_y(varnm,numvars,i1,i2,j1,dpth,ndims,dims,
     >               stag,min,max,miss,dunits,disunits,dim_names)
        if (int1.eq.0) then
          write(*,*) 'Error:  could not compute line integral'
          return
        end if
        int2 = line_int_y(varnm,numvars,i1,i2,j2,dpth,ndims,dims,
     >               stag,min,max,miss,dunits,disunits,dim_names)
        if (int2.eq.0) then
          write(*,*) 'Error:  could not compute line integral'
          call freemem(int1)
          return
        end if
      else
        write(*,*) 'Error:  don''t know how to compute line ',
     >         'integrals for direction ',
     >             dir(strbeg(dir):strend(dir))
        return
      end if

      diff_line_int = getmem(dims(3)*dims(4))

      call diff_line_int_calc(%val(diff_line_int),%val(int1),
     >                 %val(int2),area,dims(3),dims(4))

      tmp = dunits(strbeg(dunits):strend(dunits))//'/m/m'
      call units_smpfy(tmp,dunits,err)
      if (err) dunits = ' '
 
      tmp = disunits(strbeg(disunits):strend(disunits))//'/m/m'
      call units_smpfy(tmp,disunits,err)
      if (err) disunits = ' '      

      call freemem(int1)
      call freemem(int2)

      return
      end

*-----------------------------------------------------------------

      subroutine diff_line_int_calc(diff_line_int,int1,int2,
     >                     area,nz,nt)

      integer nz,nt
      real diff_line_int(nz,nt),int1(nz,nt),int2(nz,nt),area

      integer n,k

      do n=1,nt
      do k=1,nz
        diff_line_int(k,n) = int1(k,n)/area - int2(k,n)/area
      enddo
      enddo

      return
      end

