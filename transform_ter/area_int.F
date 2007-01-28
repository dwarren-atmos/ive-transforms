
#include "pointer.icl"

      ive_ptr function area_int_xz0(varnm,numvars,i1,i2,z01,z02,
     >                        ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names)

#include "default.icl"

      integer numvars,i1,i2,ndims,dims(4)
      real stag(4),min(4),max(4),miss,z01,z02
      character*(*) varnm(numvars),dim_names(4),dunits,disunits

      integer j,strbeg,strend
      character*(80) tmp
      logical err
      ive_ptr area_int_xz0_y,area_int,getmem

      area_int_xz0 = getmem((ny+1)*ntime)
      if (area_int_xz0.eq.0) then
        write(6,*) 'Error:  unable to allocate space for ',
     >                      'area_int_xz0.'
        return
      end if

      area_int = getmem(ntime)
      if (area_int.eq.0) then
        write(6,*) 'Error:  unable to allocate space for area_int.'
        call freemem(area_int_xz0)
        return
      end if

      do j=2,ny

        area_int = area_int_xz0_y(varnm,numvars,i1,i2,z01,z02,j-1,
     >                   ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names)
        if (area_int.eq.0) then
          call freemem(area_int_xz0)
          call freemem(area_int)
          return
        end if

        call area_int_xz0_calc(%val(area_int_xz0),%val(area_int),
     >                              miss,j,ny,ntime)

      enddo

      miss = -777

      call area_int_xz0_calc(%val(area_int_xz0),%val(area_int),miss,
     >                                   1,ny,ntime)
      call area_int_xz0_calc(%val(area_int_xz0),%val(area_int),miss,
     >                                 ny+1,ny,ntime)

      dims(2) = ny+1
      stag(2) = 0.
      min(2) = plmin(2)
      max(2) = plmax(2)
      dim_names(2) = 'nyp1'

      call freemem(area_int)

      return
      end

      subroutine area_int_xz0_calc(area_int_xz0,area_int,miss,
     >                            j,ny,ntime)
 
      integer ny,ntime,j
      real area_int_xz0(ny+1,ntime),area_int(ntime),miss

      integer n

      if ((j.eq.1).or.(j.eq.ny+1)) then
        do n=1,ntime
          area_int_xz0(j,n) = miss
        enddo
      else
        do n=1,ntime
          area_int_xz0(j,n) = area_int(n)
        enddo
      end if

      return
      end

*-------------------------------------------------------------------------

      ive_ptr function area_int_yz0(varnm,numvars,j1,j2,z01,z02,
     >                        ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names)

#include "default.icl"

      integer numvars,j1,j2,ndims,dims(4)
      real stag(4),min(4),max(4),miss,z01,z02
      character*(*) varnm(numvars),dim_names(4),dunits,disunits

      integer i,strbeg,strend
      character*(80) tmp
      logical err
      ive_ptr area_int_yz0_x,area_int,getmem

      area_int_yz0 = getmem((nx+1)*ntime)
      if (area_int_yz0.eq.0) then
        write(6,*) 'Error:  unable to allocate space for ',
     >                      'area_int_yz0.'
        return
      end if

      area_int = getmem(ntime)
      if (area_int.eq.0) then
        write(6,*) 'Error:  unable to allocate space for area_int.'
        call freemem(area_int_yz0)
        return
      end if

      do i=2,nx

        area_int = area_int_yz0_x(varnm,numvars,j1,j2,z01,z02,i-1,
     >                   ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names)
        if (area_int.eq.0) then
          call freemem(area_int_yz0)
          call freemem(area_int)
          return
        end if

        call area_int_yz0_calc(%val(area_int_yz0),%val(area_int),
     >                              miss,i,nx,ntime)

      enddo

      miss = -777

      call area_int_yz0_calc(%val(area_int_yz0),%val(area_int),miss,
     >                                   1,nx,ntime)
      call area_int_yz0_calc(%val(area_int_yz0),%val(area_int),miss,
     >                                 nx+1,nx,ntime)

      dims(1) = nx+1
      stag(1) = 0.
      min(1) = plmin(1)
      max(1) = plmax(1)
      dim_names(1) = 'nxp1'

      call freemem(area_int)

      return
      end

      subroutine area_int_yz0_calc(area_int_yz0,area_int,miss,
     >                            i,nx,ntime)
 
      integer nx,ntime,i
      real area_int_yz0(nx+1,ntime),area_int(ntime),miss

      integer n

      if ((i.eq.1).or.(i.eq.nx+1)) then
        do n=1,ntime
          area_int_yz0(i,n) = miss
        enddo
      else
        do n=1,ntime
          area_int_yz0(i,n) = area_int(n)
        enddo
      end if

      return
      end

*------------------------------------------------------------------------

      ive_ptr function flux_density_yz0_x(varnm,numvars,z01,z02,i,
     >                  ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names)

#include "default.icl"

      integer numvars,i,ndims,dims(4)
      real stag(4),min(4),max(4),miss,z01,z02
      character*(*) varnm(numvars),dim_names(4),dunits,disunits

      integer j,strbeg,strend
      character*(80) tmp
      logical err
      ive_ptr area_int_yz0_x,area_int,getmem

      flux_density_yz0_x = getmem((ny+1)*ntime)
      if (flux_density_yz0_x.eq.0) then
        write(6,*) 'Error:  unable to allocate space for ',
     >                         'flux_density_yz0_x.'
        return
      end if

      area_int = getmem(ntime)
      if (area_int.eq.0) then
        write(6,*) 'Error:  unable to allocate space for area_int.'
        call freemem(flux_density_yz0_x)
        return
      end if

      do j=2,ny

        area_int = area_int_yz0_x(varnm,numvars,j-2,j,z01,z02,i,
     >                       ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names)
        if (area_int.eq.0) then
          call freemem(flux_density_yz0_x)
          call freemem(area_int)
          return
        end if

        call flux_density_yz0_x_calc(%val(flux_density_yz0_x),
     >             %val(area_int),dy,miss,j,ny,ntime)

      enddo

      miss = -777

      call flux_density_yz0_x_calc(%val(flux_density_yz0_x),
     >               %val(area_int),dy,miss,1,ny,ntime)
      call flux_density_yz0_x_calc(%val(flux_density_yz0_x),
     >               %val(area_int),dy,miss,ny+1,ny,ntime)

      dims(2) = ny+1
      stag(2) = 0.
      min(2) = plmin(2)
      max(2) = plmax(2)
      dim_names(2) = 'nyp1'
    
      call freemem(area_int)
 
      return
      end

      subroutine flux_density_yz0_x_calc(flux_density_yz0_x,area_int,
     >                           dy,miss,j,ny,ntime)

      integer j,ny,ntime
      real flux_density_yz0_x(ny+1,ntime),area_int(ntime),dy,miss
           
      integer n

      if ((j.eq.1).or.(j.eq.ny+1)) then
        do n=1,ntime
          flux_density_yz0_x(j,n) = miss
        enddo
      else
        do n=1,ntime
          flux_density_yz0_x(j,n) = area_int(n)/(2.*dy) 
        enddo
      end if

      return
      end
