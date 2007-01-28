
#include "pointer.icl"

      ive_ptr function filter(fld,ndims,dims,stag,min,max,miss,
     >                           dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      character*(*) fld

      integer iu,iv,iw,ijk,nt,strbeg,strend
      real fdamp,dt,r_val
      logical new
      ive_ptr getmem,var,b0,fdamp_ptr,dt_ptr

      filter = 0

      nt = ntime

      fdamp_ptr = getvar('fdamp',ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names,new)
      if (fdamp_ptr.eq.0) then
        write(6,*) 'Variable fdamp not found.  Specify damping ',
     >                    'coefficient:'
        read(6,*) fdamp
      else
        fdamp = r_val(%val(fdamp_ptr),1)
      end if

      dt_ptr = getvar('dt',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      if (dt_ptr.eq.0) then
        write(6,*) 'Variable dt not found.  Specify time step:'
        read(6,*) dt
      else
        dt = r_val(%val(dt_ptr),1)
      end if 

      if (fld(strbeg(fld):strend(fld)).eq.'B') then
        b0 = getvar('b0',ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new)
        if (b0.eq.0) then
          write(6,*) 'Error:  cannot access field b0.'
          return
        end if
      end if

      var = getvar(fld,ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot access field ',
     >               fld(strbeg(fld):strend(fld))
        return
      end if

      if (fld(strbeg(fld):strend(fld)).eq.'B') 
     >              call sub_b0(%val(var),%val(b0),nx,ny,nz,nt)

      iu = 0
      iv = 0
      iw = 0
      if (dims(1).eq.nx+1) then
        iu = 1
      else if (dims(2).eq.ny+1) then
        iv = 1
      else if (dims(3).eq.nz+1) then
        iw = 1
      end if

      filter = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (filter.eq.0) then
        write(6,*) 'Error:  cannot allocate space for filter.'
        return
      end if

      do ijk=1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(filter),ijk,0.)
      enddo

      write(6,*)
      write(6,*) 'Computing fourth-derivative filter on field ',
     >              fld(strbeg(fld):strend(fld)),' with'
      write(6,*) '     fdamp = ',fdamp
      write(6,*) '     dt = ',dt
      write(6,*)

      call filter_calc(%val(filter),%val(var),fdamp,dt,rdx,rdy,
     >                      nx,ny,nz,nt,iu,iv,iw)

      if (fld(strbeg(fld):strend(fld)).eq.'B') 
     >              call add_b0(%val(var),%val(b0),nx,ny,nz,nt)

      dunits = dunits(strbeg(dunits):strend(dunits))//'/s'
      disunits = disunits(strbeg(disunits):strend(disunits))//'/s'

      return
      end

*-----------------------------------------------------------------------

      subroutine filter_calc(fs,var,fdamp,dt,rdx,rdy,nx,ny,nz,nt,
     >                                    iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real fs(nx+iu,ny+iv,nz+iw,nt),var(nx+iu,ny+iv,nz+iw,nt)
      real fdamp,dt,rdx,rdy

      integer i,j,k,n
      real f4damp,f2damp

      f4damp = 0.0625*fdamp/2./dt
      f2damp = 0.25*fdamp/2./dt

      do n=1,nt
      do k=1+iw,nz
*
* x pass
*
      do j=1+iv,ny
      do i=3,nx+iu-2
        fs(i,j,k,n) = fs(i,j,k,n) 
     >    + f4damp*(-var(i+2,j,k,n)+4.*var(i+1,j,k,n)
     >       - 6.*var(i,j,k,n)+4.*var(i-1,j,k,n)-var(i-2,j,k,n))
      enddo
      enddo
*
* boundaries
*
      do j=1+iv,ny
        i = 2
        fs(i,j,k,n) = fs(i,j,k,n)
     >     + f2damp*(var(i+1,j,k,n)-2.*var(i,j,k,n)+var(i-1,j,k,n))
        i = nx+iu-1
        fs(i,j,k,n) = fs(i,j,k,n)
     >     + f2damp*(var(i+1,j,k,n)-2.*var(i,j,k,n)+var(i-1,j,k,n))
      enddo
*
* y pass
*
      do j=3,ny+iv-2
      do i=1+iu,nx
        fs(i,j,k,n) = fs(i,j,k,n)
     >    + f4damp*(-var(i,j+2,k,n)+4.*var(i,j+1,k,n)
     >       - 6.*var(i,j,k,n)+4.*var(i,j-1,k,n)-var(i,j-2,k,n))
      enddo
      enddo
*
* boundaries
*
      do i=1+iu,nx
        j = 2
        fs(i,j,k,n) = fs(i,j,k,n)
     >     + f2damp*(var(i,j+1,k,n)-2.*var(i,j,k,n)+var(i,j-1,k,n))
        j = ny+iv-1
        fs(i,j,k,n) = fs(i,j,k,n)
     >     + f2damp*(var(i,j+1,k,n)-2.*var(i,j,k,n)+var(i,j-1,k,n))
      enddo
*
* end outer loops
*
      enddo
      enddo

      return
      end
