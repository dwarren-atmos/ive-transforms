
#include "pointer.icl"

      ive_ptr function h_ave(whch,fld,n1,n2,ndims,dims,stag,min,
     >                     max,miss,dunits,disunits,dim_names)
*
* h_ave calculates the volume average of fld i) over x as a function
*   of y,z and t; or ii) over y as a function of x,z and t
*

#include "default.icl"

      integer n1,n2,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) whch,fld,dunits,disunits,dim_names(4)

      integer strbeg,strend,iu,iv,iw,istrt,istop,jstrt,jstop,k
      real r_val
      logical new,err,requ
      ive_ptr getmem,v_sgz,var,gz,sgz_tmp

      h_ave = 0

      if ((whch.ne.'X').and.(whch.ne.'Y')) then
        write(*,*) 'Error:  do not know over which direction ',
     >                          'to average.'
        return
      else if ((whch.eq.'Y').and.(twod.eq.1)) then
        write(*,*) 'Error:  cannot average over y for 2d run.'
        return
      end if

      sgz_tmp = getmem(nz)
      if (sgz_tmp.eq.0) then
        write(*,*) 'Error:  unable to allocate space for sgz_tmp.'
        return
      end if
      do k=1,nz
        call s_val(%val(sgz_tmp),k,r_val(%val(sgz),k))
      enddo

      var = getvar(fld,ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot find field ' ,
     >                         fld(strbeg(fld):strend(fld))
        call freemem(sgz_tmp)
        return
      end if

      iu = 0
      iv = 0
      iw = 0
      if (requ(plmin(1),min(1))) then
        iu = 1
      else if (requ(plmin(2),min(2))) then
        iv = 1
      else if (requ(plmin(3),min(3))) then
        iw = 1
      end if

      dims(1) = dims(1) - iu
      dims(2) = dims(2) - iv
      dims(3) = dims(3) - iw

      v_sgz = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (v_sgz.eq.0) then
        write(*,*) 'Error:  unable to allocate space for v_sgz.'
        call freemem(sgz_tmp)
        return
      end if
      do k = 1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(v_sgz),k,1.)
      enddo

      if (iw.eq.1) then
        gz = wgz
      else
        gz = sgz
      end if

      miss = -777
      call v_2_sgz(%val(v_sgz),%val(var),%val(zbot_p),%val(gz),
     >           %val(sgz_tmp),ztop,miss,dims(1),dims(2),dims(3),
     >           dims(4),iu,iv,iw)

      if (whch.eq.'X') then
      
        istrt = n1+1
        istop = n2 
        
 123    FORMAT(10x,'x1 = ',F10.2,2x,'x2 = ',F10.2)

        write(*,*) 
        write(*,*) 'Calculating volume average over x of field ',
     >                   fld(strbeg(fld):strend(fld))
        write(*,*) 'over the domain defined by:  '
        write(*,123) plmin(1)+n1*dx,plmin(1)+n2*dx
        write(*,*)

        if ((n1.lt.0).or.(n2.gt.dims(1))) then
          write(*,*) 'Error:  specified limits exceed domain ',
     >                            'boundaries.'
          call freemem(v_sgz)
          call freemem(sgz_tmp)
          return
        else if (n1.ge.n2) then
          write(*,*) 'Error:  improper specification of domain ',
     >                          'of integration.' 
          call freemem(v_sgz)
          call freemem(sgz_tmp)
          return
        end if

        h_ave = getmem(dims(2)*dims(3)*dims(4))
        if (h_ave.eq.0) then
          write(*,*) 'Error:  unable to allocate space for h_ave.'
          call freemem(v_sgz)
          call freemem(sgz_tmp)
          return
        end if

        call x_ave_calc(%val(h_ave),%val(v_sgz),%val(zbot_p),
     >              %val(sgz),miss,istrt,istop,dims(1),dims(2),
     >              dims(3),dims(4))

        dims(1) = 1
        min(1) = 0.
        max(1) = 0.
        stag(1) = 0.
        dim_names(1) = 'one'

        min(2) = plmin(2) + 0.5*dy
        max(2) = plmax(2) - 0.5*dy
        min(3) = sgz_min
        max(3) = sgz_max
        stag(2) = 0.5
        stag(3) = 0.5
        dim_names(2) = 'ny'
        dim_names(3) = 'nz'
   
      else

        jstrt = n1+1
        jstop = n2
      
124     FORMAT(10x,'y1 = ',F10.2,2x,'y2 = ',F10.2)

        write(*,*) 
        write(*,*) 'Calculating volume average over y of field ',
     >                          fld(strbeg(fld):strend(fld))
        write(*,*) 'over the domain defined by: '
        write(*,124) plmin(2)+n1*dy,plmin(2)+n2*dy
        write(*,*)

        if ((n1.lt.0).or.(n2.gt.dims(2))) then
          write(*,*) 'Error:  specified limits exceed domain ',
     >                            'boundaries.'
          call freemem(v_sgz)
          call freemem(sgz_tmp)
          return
        else if (n2.le.n1) then
          write(*,*) 'Error:  improper specification of domain ',
     >                         'of integration.'
          call freemem(v_sgz)
          call freemem(sgz_tmp)
          return
        end if

        h_ave = getmem(dims(1)*dims(3)*dims(4))
        if (h_ave.eq.0) then
          write(*,*) 'Error:  unable to allocate space for h_ave.'
          call freemem(v_sgz)
          call freemem(sgz_tmp)
          return
        end if

        call y_ave_calc(%val(h_ave),%val(v_sgz),%val(zbot_p),
     >              %val(sgz),miss,jstrt,jstop,dims(1),dims(2),
     >              dims(3),dims(4))

        dims(2) = 1
        min(2) = 0.
        max(2) = 0.
        stag(2) = 0.
        dim_names(2) = 'one'

        min(1) = plmin(1) + 0.5*dx
        max(1) = plmax(1) - 0.5*dx
        min(3) = sgz_min 
        max(3) = sgz_max
        stag(1) = 0.5
        stag(3) = 0.5
        dim_names(1) = 'nx'
        dim_names(3) = 'nz'

      end if

      call freemem(v_sgz)
      call freemem(sgz_tmp)

      return
      end

*------------------------------------------------------------------------

      subroutine x_ave_calc(x_ave,v_sgz,zs,sgz,miss,istrt,istop,
     >                           nx,ny,nz,nt)   

      integer istrt,istop,nx,ny,nz,nt
      real v_sgz(nx,ny,nz,nt),x_ave(ny,nz,nt),zs(nx,ny),sgz(nz),
     >            miss

      integer n,j,i,k
      real volume

      do n=1,nt
      do k=1,nz
      do j=1,ny

        x_ave(j,k,n) = 0.
        volume = 0.

        do i = istrt,istop
          if (sgz(k).ge.zs(i,j)) then
            x_ave(j,k,n) = x_ave(j,k,n) + v_sgz(i,j,k,n)
            volume = volume + 1.
          end if
        enddo

        if (volume.eq.0.) then
          x_ave(j,k,n) = miss
        else
          x_ave(j,k,n) = x_ave(j,k,n)/volume
        end if

      enddo
      enddo
      enddo

      return
      end          

*-----------------------------------------------------------------------

      subroutine y_ave_calc(y_ave,v_sgz,zs,sgz,miss,jstrt,jstop,
     >                           nx,ny,nz,nt)   

      integer jstrt,jstop,nx,ny,nz,nt
      real v_sgz(nx,ny,nz,nt),y_ave(nx,nz,nt),zs(nx,ny),sgz(nz),
     >                 miss

      integer n,i,j,k
      real volume

      do n=1,nt
      do k=1,nz
      do i=1,nx

        y_ave(i,k,n) = 0.
        volume = 0.

        do j = jstrt,jstop
          if (sgz(k).ge.zs(i,j)) then
            y_ave(i,k,n) = y_ave(i,k,n) + v_sgz(i,j,k,n)
            volume = volume + 1.
          end if
        enddo

        if (volume.eq.0.) then
          y_ave(i,k,n) = miss
        else
          y_ave(i,k,n) = y_ave(i,k,n)/volume
        end if

      enddo
      enddo
      enddo

      return
      end          

