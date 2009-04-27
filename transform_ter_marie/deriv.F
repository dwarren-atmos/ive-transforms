
#include "pointer.icl"

      ive_ptr function deriv(varnm,dir,ndims,dims,stag,min,max,
     >                        miss,dunits,disunits,dim_names)
*
*  deriv takes a variable name and a direction as input and returns
*    a pointer to the gradient of the variable in the specified
*    direction.  In this routine we assume that the input variable is
*    defined on one of the standard grids (i.e., u,v,w, or thermo).
*
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) varnm,dir,dunits,disunits,dim_names(4)
      
      integer ndims2,dims2(4)
      real stag2(4),min2(4),max2(4),miss2
      character*(80) dunits2,disunits2,dim_names2(4)
      
      integer idir,iu,iv,iw,k,strbeg,strend,ibeg,iend
      integer strcasecmp
      real dom_min(4),del(4),ztop,r_val
      real*8 slope,intercept
      character*80 dim_dis(4),dom_units(4),str1,str2,ternm,gznm
      ive_ptr vptr,tptr,ttptr,zptr,tzptr,getvar,getmem
      logical new,err

      deriv = 0
      call getaarr('dim_display_names',dim_dis,4,err)
      if (err.eq..true.) then
        write(*,*) 'Error:  could not access display names for ',
     >                         'coordinates'
        return
      end if
      call getaarr('domain_units',dom_units,4,err)
      if (err.eq..true.) then
        write(*,*) 'Error:  could not access units of ',
     >                         'domain specification.'
        return
      end if
      call getrarr('plmin',dom_min,4,err)
      if (err.eq..true.) then
        write(*,*) 'Error:  could not access physical space ',
     >                       'coordinate minimums.'
        return
      end if
 
      call setstr(str1,dir)
      call setstr(str2,dim_dis(1))

      k = 1
 79   if ((strcasecmp(str1,str2).ne.0).and.(k.lt.4)) then
        k = k+1
        call setstr(str2,dim_dis(k))
        goto 79
      end if

      if (k.eq.4) then
        if (strcasecmp(str1,str2).ne.0) then
          write(*,*) 'Error:  cannot find dimension ', 
     >                    dir(strbeg(dir):strend(dir))
          return
        end if
      end if

      idir = k

      tzptr = getvar('wgz',ndims2,dims2,stag2,min2,max2,miss2,
     >                     dunits2,disunits2,dim_names2,new)
      if (tzptr.eq.0) then
        write(*,*) 'Error:  cannot find variable wgz.'
        return
      end if
      call convert_units(dunits2,dom_units(3),slope,intercept)
      ztop = real(r_val(%val(tzptr),dims2(1))*slope + intercept)

      vptr = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new)
      if (vptr.eq.0) then
        write(*,*) 'Error:  cannot find variable ',
     >               varnm(strbeg(varnm):strend(varnm))
        return
      end if

      if (dims(idir).eq.1) then
        write(*,*) 'Error:  cannot take derivative of trivial ',
     >                 'dimension ',dir(strbeg(dir):strend(dir))
        return
      end if

      del(1) = (max(1)-min(1))/(dims(1)-1)
      del(2) = (max(2)-min(2))/(dims(2)-1) 
      del(4) = (max(4)-min(4))/(dims(4)-1)

      iu = 0
      iv = 0
      iw = 0
      if (dom_min(1).ge.min(1)) then
        iu = 1
        ternm(1:6) = 'zbot_u'
        gznm(1:3) = 'sgz'
      else if (dom_min(2).ge.min(2)) then
        iv = 1
        ternm(1:6) = 'zbot_v'
        gznm(1:3) = 'sgz'
      else if (dom_min(3).ge.min(3)) then
        iw = 1
        ternm(1:6) = 'zbot_p'
        gznm(1:3) = 'wgz'
      else 
        ternm(1:6) = 'zbot_p'
        gznm(1:3) = 'sgz'
      end if

      ttptr = getvar(ternm(1:6),ndims2,dims2,stag2,min2,max2,miss2,
     >                     dunits2,disunits2,dim_names2,new)
      if (ttptr.eq.0) then
        write(*,*) 'Error:  could not find terrain array ',ternm(1:6)
        return
      end if
      tptr = getmem(dims2(1)*dims2(2))
      call convert_units(dunits2,dom_units(3),slope,intercept)
      do k=1,dims2(1)*dims2(2)
        call s_val(%val(tptr),k,
     >               real(r_val(%val(ttptr),k)*slope + intercept))
      enddo

      tzptr = getvar(gznm(1:3),ndims2,dims2,stag2,min2,max2,miss2,
     >                   dunits2,disunits2,dim_names2,new)
      if (tzptr.eq.0) then
        write(*,*) 'Error:  could not find vertical coordinate ',
     >                       'array ',gznm(1:3)
        return
      end if
      zptr = getmem(dims2(1))
      call convert_units(dunits2,dom_units(3),slope,intercept)
      do k=1,dims2(1)
        call s_val(%val(zptr),k,
     >               real(r_val(%val(tzptr),k)*slope + intercept))
      enddo

      ibeg = strbeg(dom_units(idir))
      iend = strend(dom_units(idir))
      dunits2 = dunits(strbeg(dunits):strend(dunits))//'/'
      dunits2 = dunits2(strbeg(dunits2):strend(dunits2))//
     >                        dom_units(idir)(ibeg:iend)
      call units_smpfy(dunits2,dunits,err)
      if (err) dunits = dunits2

      call getaarr('domain_display_units',dom_units,4,err)
      ibeg = strbeg(dom_units(idir))
      iend = strend(dom_units(idir))
      disunits2 = disunits(strbeg(disunits):strend(disunits))//'/'
      disunits2 = disunits2(strbeg(disunits2):strend(disunits2))//
     >                        dom_units(idir)(ibeg:iend)
      call units_smpfy(disunits2,disunits,err)
      if (err) disunits = disunits2

      deriv = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (deriv.eq.0) then
        write(*,*) 'Error:  cannot allocate space for deriv.'
        return
      end if

      call diff(%val(deriv),%val(vptr),%val(tptr),%val(zptr),ztop,
     >           del(1),del(2),del(4),dims(1),dims(2),dims(3),
     >           dims(4),idir)

      call freemem(tptr)
      call freemem(zptr)

      return
      end




