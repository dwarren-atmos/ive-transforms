      integer(ive_ptr) function abs_ive(varnm,ndims,dims,stag,min,max,
     &                          miss,dunits,disunits,dim_names)
      use default_module
*
*  abs_ive returns the absolute value of field varnm
*
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) varnm,dunits,disunits,dim_names(4)

      integer strbeg,strend
      logical new
      integer(ive_ptr) var,getvar,getmem

      abs_ive = 0

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new)
      if (var.eq.0) then
        write(*,*) 'Error:  cannot get field ',
     >                   varnm(strbeg(varnm):strend(varnm))
        return
      end if

      abs_ive = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (abs_ive.eq.0) then
        write(*,*) 'Error:  cannot allocate space for field abs.' 
        return
      end if

      call abs_calc(%val(abs_ive),%val(var),dims(1),dims(2),dims(3),
     >                            dims(4))

      return
      end

*-------------------------------------------------------------------

      subroutine abs_calc(abs_ive,var,nx,ny,nz,nt)
  
      integer nx,ny,nz,nt
      real abs_ive(nx,ny,nz,nt),var(nx,ny,nz,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        abs_ive(i,j,k,n) = abs(var(i,j,k,n))

      enddo
      enddo
      enddo
      enddo

      return
      end

   
