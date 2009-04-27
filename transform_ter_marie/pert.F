
#include "pointer.icl"

      ive_ptr function pert(varnm,basnm,ndims,dims,stag,min,max,
     >                  missing,dunits,disunits,dim_names)
*
*  pert calculates the deviation of variable varnm(nx,ny,nz,nt)
*    from the basic state variable basnm(nx,ny,nz,1).  The two
*    variables must be defined on the same grid.
*
      integer ndims,dims(4)
      real stag(4),min(4),max(4),missing
      character*(*) varnm,basnm,dunits,disunits,dim_names(4)

      integer strbeg,strend,k
      real bmin(3)
      logical new
      ive_ptr var,base,getvar,getmem

      pert = 0

      base = getvar(basnm,ndims,dims,stag,min,max,missing,dunits,
     >                 disunits,dim_names,new)
      if (base.eq.0) then
        write(*,*) 'Error:  cannot compute field ',
     >                basnm(strbeg(basnm):strend(basnm))
        return
      else if (dims(4).ne.1) then
        write(*,*) 'Error:  field ',
     >            basnm(strbeg(basnm):strend(basnm)),
     >              ' is not a basic state field.'
        return
      end if

      do k=1,3
        bmin(k) = min(k)
      enddo

      var = getvar(varnm,ndims,dims,stag,min,max,missing,dunits,
     >                   disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot compute field ',
     >                varnm(strbeg(varnm):strend(varnm))
        return
      end if

      k = 1
 45   if (k.le.3) then
        if (bmin(k).ne.min(k)) then
          write(*,*) 'Error:  fields ',
     >            basnm(strbeg(basnm):strend(basnm)),' and ',
     >              varnm(strbeg(varnm):strend(varnm)),
     >           ' are defined on different grids.'
          return
        end if
        k = k+1
        goto 45
      end if

      pert = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (pert.eq.0) then
        write(*,*) 'Error:  could not allocate space for pert.'
        return
      end if
 
      call pert_calc(%val(pert),%val(var),%val(base),dims(1),
     >                     dims(2),dims(3),dims(4))

      return
      end

************************************************************************

      subroutine pert_calc(pert,var,base,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real pert(nx,ny,nz,nt),var(nx,ny,nz,nt),base(nx,ny,nz)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        pert(i,j,k,n) = var(i,j,k,n) - base(i,j,k)

      enddo
      enddo
      enddo
      enddo
  
      return
      end 
 
    

      
