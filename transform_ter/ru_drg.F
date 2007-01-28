
#include "pointer.icl"

      ive_ptr function ru_drg(i1,i2,j1,j2,k1,k2,n0,ndims,dims,stag,
     >                  min,max,miss,dunits,disunits,dim_names)
*
*  ru_drg calculates the average perturbation momentum per unit volume
*    that would result if all the momentum loss predicted by the convergence
*    of the vertical momentum flux (or drag) on the specified domain
*    (i1,i2,j1,j2,k1,k2) were in fact lossed on the specified domain 
*    (i.e., as if the horizontal pressure gradients and momentum fluxes 
*    were zero.)
*
      integer i1,i2,j1,j2,k1,k2,n0,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer n,strbeg,strend,iarg(10),nx,ny
      real dom_min(4),dom_max(4),area,volume,dz,dx,dy,r_val,fcv,
     >        rarg(10),fac
      logical new,err,rotating
      character*(80) varnm(10),tmp,zfstr,carg(10),mkfield
      ive_ptr zf1,zf2,rtmp,ter,run_int,getvar,getmem,wgz,fcor

      ru_drg = 0

      fcor = getvar('f_cor',ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names,new)
      if (fcor.eq.0) then
        write(*,*) 'Error:  cannot find coriolis parameter.'
        return
      end if

      if (r_val(%val(fcor),1).eq.0.) then
        rotating = .false.
      else
        rotating = .true.
      end if

      if (k1.eq.0) then

        iarg(1) = i1
        iarg(2) = i2
        iarg(3) = j1
        iarg(4) = j2
        iarg(5) = k2

        zfstr = mkfield('DRG_BOT',carg,iarg,rarg,0,5,0,0)

        zf1 = run_int(zfstr,n0,ndims,dims,stag,min,max,miss,dunits,
     >                            disunits,dim_names)
        if (zf1.eq.0) then
          write(*,*) 'Error:  cannot calculate running integral ',
     >                              'of pressure drag.'
          return
        end if

        zfstr = mkfield('ZF_TOP',carg,iarg,rarg,0,5,0,0)

        zf2 = run_int(zfstr,n0,ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        if (zf2.eq.0) then
          write(*,*) 'Error:  cannot calculate running integral ',
     >                          'of vertical momentum flux.'
          call freemem(zf1)
          return
        end if

        miss = 0.060458

        ru_drg = getmem(dims(4))

        do n=1,n0
          call s_val(%val(ru_drg),n,miss)
        enddo

        do n=n0+1,dims(4)
          fcv = r_val(%val(zf1),n) + r_val(%val(zf2),n)
          call s_val(%val(ru_drg),n,fcv)
        enddo

        call freemem(zf1)
        call freemem(zf2)

      else

        carg(1) = 'Z'
        carg(2) = 'RHO'
        if (rotating) then
          carg(3) = 'UABS'
        else
          carg(3) = 'U'
        end if
        carg(4) = 'W'
        iarg(1) = i1
        iarg(2) = i2
        iarg(3) = j1
        iarg(4) = j2
        iarg(5) = k1
        iarg(6) = k2

        zfstr = mkfield('DIFF_AREA_INT',carg,iarg,rarg,4,6,0,0)

        ru_drg = run_int(zfstr,n0,ndims,dims,stag,min,max,miss,
     >                   dunits,disunits,dim_names)
        if (ru_drg.eq.0) then 
          write(*,*) 'Error:  unable to calculate running integral ',
     >                      'of vertical momentum fluxes.'
          return
        end if

      end if
 
      return
      end

