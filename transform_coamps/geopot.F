!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Finds the first physical height of the given field.
! The search starts from either the top or the bottom of the domain.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(ive_ptr) function findz(fldnam,p0,ltop,lpre,
     &  ndims,dims,stag,min,max,miss,dunits,
     &  disunits,dim_names)

      use default_module
      implicit none

      character(len=*), intent(in) :: fldnam
      real, intent(in) :: p0
      logical, intent(in) :: ltop,lpre
c
      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(out) :: dunits,disunits,dim_names(4)
c      
      integer :: i,strbeg,strend,ix,iy,iz,it,mn
      real :: r_val, p0_wrk
      logical :: new
      integer(ive_ptr) var, zbot_wrk
c
      var = 0
      var = getvar(fldnam,ndims,dims,stag,min,max,miss,
     &             dunits,disunits,dim_names,new)
      if (var.le.0) then
        write(6,*) 'Error:  cannot access variable ',trim(fldnam)
        return
      end if
c
      ix=dims(1) ; iy=dims(2) ; iz=dims(3) ; it=ntime ; mn = ix*iy
c
      findz = 0 ; findz = getmem(mn*it)
      if (findz.le.0) then
        write(6,*) 'Error:  cannot allocate space for findz.'
        return
      end if
c
      p0_wrk = p0
c      if(lpre) then
c        p0_wrk = alog(p0_wrk)
c        do i=1,mn*it
c          call s_val(%val(var),i,alog(r_val(%val(var),i)))
c        end do
c      end if          

      zbot_wrk = zbot_p
      !zbot_wrk = zbot_u
      !zbot_wrk = zbot_v

      if(.not.ltop) then
      call findz_bot(%val(findz),p0_wrk,%val(var),%val(sgz),
     &                 %val(zbot_wrk),ztop,miss,mn,iz,it)
      else
        print *,'Search from top is not yet implemented.'
        call freemem(findz) ; findz=0
        return
      end if

c      if(lpre) then
c        do i=1,mn*iz*it
c          call s_val(%val(var),i,exp(r_val(%val(var),i)))
c        end do
c      end if          

c
      dims(1) = ix ; dims(2) = iy ; dims(3) = 1 ; dims(4) = ntime
      dim_names(3) = 'one' ; min(3) = 0 ; max(3) = 0
      dunits = 'm' ; disunits = 'decameters'
c
      return
      end

      subroutine findz_bot(fout,f0,fin,sgz,zsurf,ztop,missing,
     &                       mn,nz,nt)

      integer :: mn, nz, nt
	  real :: f0, fin(mn,nz,nt), fout(mn,nt), sgz(nz), zsurf(mn),
     &     ztop, missing
c
      integer i,k,n,k1,k2
	  real temp, zfrac
c
      do n=1,nt ! outer time loop
c
c   Assign missing date 
c
        do i=1,mn
          if(f0.gt.fin(i,1,n).or.f0.lt.fin(i,nz,n)) then
            fout(i,n) = missing
          end if
        end do
c
        do k=2,nz
          do i=1,mn
            if ( f0.lt.fin(i,k-1,n).and.f0.ge.fin(i,k,n)) then
              zfrac = (ztop - zsurf(i))/ztop
              k1 = k-1 ; k2 = k
              temp =  ( alog(f0) - alog(fin(i,k1,n)) )/
     &                ( alog(fin(i,k2,n)) - alog(fin(i,k1,n)) )
              temp = sgz(k1) + temp*( sgz(k2) - sgz(k1) ) 
              fout(i,n) = temp*zfrac + zsurf(i)
            end if
          end do
        end do
      end do

      return
      end subroutine

