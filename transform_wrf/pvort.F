!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   function vort_z(uvar_name, vvar_name, ...)
!
!   Calculates the vertical component of vorticity 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(ive_ptr) function pvort(ndims,dims,stag,min,max,miss,
     &                                dunits,disunits,dim_names)

      USE DEFAULT_MODULE
c     have nx,ny,nz,mapM,mapU,mapV,ptot in default 

      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=80), intent(inout):: dunits,disunits,dim_names(4)
      
      integer :: i, ii, nt
      real :: r_val
      logical :: new, newu, newv, newt, newp, newc

      integer(ive_ptr) :: cor, myu, myv, myt, myp

      nt=ntime 
!
      newu = .false.; newv = .false.; newt = .false.; newp = .false.
      newc=.false.
      myu = 0; myv = 0; myt = 0; myp = 0
      pvort=0
!     We know total pressure exists from new_file_trans
      cor = getvar('F',ndims,dims,stag,min,max,
     >     miss,dunits,disunits,dim_names,newc)
      if(cor.eq.0) then
        print *,'Error: cannot allocate nesecary memory in pvort'
        return
      end if
      myu = getvar('U@Pressure',ndims,dims,stag,min,max,
     >     miss,dunits,disunits,dim_names,newu)
      if(myu.eq.0) then
        print *,'Error: cannot allocate nesecary memory in pvort'
        if(newc)call freemem(cor)
        return
      end if
      myv = getvar('V@Pressure',ndims,dims,stag,min,max,
     >     miss,dunits,disunits,dim_names,newv)
      if(myv.eq.0) then
        print *,'Error: cannont allocate nesecary memory in pvort'
        if(newu)call driver('free=U@Pressure',.false.)
        return
      end if
      myt = getvar('T@Pressure',ndims,dims,stag,min,max,
     >     miss,dunits,disunits,dim_names,newt)
      if(myt.eq.0) then
        print *,'Error: cannont allocate nesecary memory in pvort'
        if(newc)call freemem(cor)
        if(newu)call freemem(myu)
        if(newv)call freemem(myv)
        return
      end if
      myp = getvar('Pressure',ndims,dims,stag,min,max,
     >     miss,dunits,disunits,dim_names,newp)
      if(myp.eq.0) then
        print *,'Error: cannont allocate nesecary memory in pvort'
        if(newc)call freemem(cor)
        if(newu)call freemem(myu)
        if(newv)call freemem(myv)
        if(newt)call freemem(myt)
        return
      end if


      pvort = getmem(nx*ny*nz*nt)
      write(6,*)'alloc ',nx*ny*nz*nt,'for pvort'
      if(pvort.eq.0) then
        print *,'Error: cannont allocate nesecary memory in pvort'
        return
        if(newc)call freemem(cor)
        if(newu)call freemem(myu)
        if(newv)call freemem(myv)
        if(newt)call freemem(myt)
      end if
!      do i=1,mnk ! initialize pvort to 0
!        call s_val(%val(pvort),i,0.0)
!      end do

      call calc_pvort(%val(pvort), 
     &     %val(myu), %val(myv), %val(myp), %val(myt),
     &     %val(mapU_ptr), %val(mapV_ptr), %val(mapM_ptr),
     &     %val(cor),
     &     nx, ny, nz, nt, dx, dy)

      write(6,*) 'set dunits etc'
      dunits='K*m*m/kg/s' ; disunits='10^-6*K*m*m/kg/s'
      

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Does the potential vorticity calculation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_pvort(pv,u,v,p,t,mu,mv,mf,
     &                      cor,ix,iy,iz,it,dx,dy)

      integer, intent(in) :: ix,iy,iz,it
      real, intent(in) :: u(ix,iy,iz,it), v(ix,iy,iz,it),
     &                    p(ix,iy,iz,it), t(ix,iy,iz,it),
     &                    mu(ix,iy,it), mv(ix,iy,it),
     &                    mf(ix,iy,it), cor(ix,iy,it),
     &                    dx,dy

      real, intent(out) :: pv(ix,iy,iz,it)

      integer i,j,k,n,ip1,im1,jp1,jm1,kp1,km1
      real mm, dsx, dsy, dudy, dvdx, avort, dp, dudp, dvdp,
     &     dtdp, dtdx, dtdy

!Init it first
!      pv=0.0
      do n=1,it 
         do k=1,iz
            kp1 = min(k+1, iz)
            km1 = max(k-1, 1)
            do j=1,iy 
               jp1 = min(j+1, iy)
               jm1 = max(j-1, 1)
               dsy =(jp1-jm1)*dy
               do i=1,ix
                  ip1 = min(i+1, ix)
                  im1 = max(i-1, 1)
                  dsx =(ip1-im1)*dx
                  mm = mf(i,j,n)**2
                  dudy = 0.50*(u(i,jp1,k,n)/mu(i,jp1,n) +
     &                 u(ip1,jp1,k,n)/mu(ip1,jp1,n) -
     &                 u(i,jm1,k,n)/mu(i,jm1,n) -
     &                 u(ip1,jm1,k,n)/mu(ip1,jm1,n))/dsy*mm
                  dvdx = 0.50*(v(ip1,j,k,n)/mv(ip1,j,n) +
     &                 v(ip1,jp1,k,n)/mv(ip1,jp1,n) -
     &                 v(im1,j,k,n)/mv(im1,j,n) -
     &                 v(im1,jp1,k,n)/mv(im1,jp1,n))/dsx*mm
                  avort = dvdx - dudy + cor(i,j,n)

                  dp = p(i,j,kp1,n) - p(i,j,km1,n)
                  dudp = 0.50*(u(i,j,kp1,n)+u(ip1,j,kp1,n) -
     &                 u(i,j,km1,n) - u(ip1,j,km1,n))/dp
                  dvdp = 0.50*(v(i,j,kp1,n)+v(i,jp1,kp1,n) -
     &                 v(i,j,km1,n) - v(i,jp1,km1,n))/dp
                  dtdp = (t(i,j,kp1,n) - t(i,j,km1,n))/dp
                  dtdx = (t(ip1,j,k,n) - t(im1,j,k,n))/
     &                 dsx*mf(i,j,n)
                  dtdy = (t(i,jp1,k,n) - t(i,jm1,k,n))/
     &                 dsy*mf(i,j,n)

                  pv(i,j,k,n) = -9.810*(dtdp*avort - dvdp*dtdx +
     &                 dudp*dtdy)*1000000.0

      end do ; end do ; end do ; end do

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine sum_prod(var_sum,var1,var2,ix,iy,iz,it)    
!
!  Adds the product of two variables to a third
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sum_prod(varsum,var1,var2,ix,iy,iz,it)

      integer, intent(in) :: ix,iy,iz,it
      real, intent(in) :: var1(ix,iy,iz,it), var2(ix,iy,iz,it)
      real, intent(inout) :: varsum(ix,iy,iz,it)

      integer i, j, k, n

      do n=1,it ; do k=1,iz ; do j=1,iy ; do i=1,ix
        varsum(i,j,k,n) = varsum(i,j,k,n) + var1(i,j,k,n)*var2(i,j,k,n)
      end do ; end do ; end do ; end do

      return
      end

