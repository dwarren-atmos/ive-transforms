*
*     Get stuff to make mass flux
*

      integer(ive_ptr) function get_vmassflux( RVAR, WVAR,
     &     ndims,dims,stag,min,max,miss,dunits,disunits,dim_names,new,
     &     err)

      use default_module
      implicit none

      character(len=*), intent(in) :: RVAR, WVAR
c
      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(out) :: dunits,disunits,dim_names(4)
      logical new, err

      integer :: i, mnk
      integer(ive_ptr) :: r,w,mf

      get_vmassflux = 0
      err = .false.
      r = getvar(RVAR, ndims,dims,stag,min,max,miss,dunits,
     >     disunits,dim_names,new)
      if(err) return
      w = getvar(WVAR, ndims,dims,stag,min,max,miss,dunits,
     >     disunits,dim_names,new)
      if(err) return
      
      mf=getmem(dims(2)*dims(3)*dims(4))
      call ver_mass_flux(dims(1),dims(2),dims(3),dims(4),
     &     %val(w),%val(r),%val(mf))
      get_vmassflux=mf
      dims(1)=1
      return
      end


      integer(ive_ptr) function get_hmassflux( RVAR, UVAR,
     &     ndims,dims,stag,min,max,miss,dunits,disunits,dim_names,new,
     &     err)

      use default_module
      implicit none

      character(len=*), intent(in) :: RVAR, UVAR
c
      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(out) :: dunits,disunits,dim_names(4)
      logical new, err

      integer :: i, mnk
      integer(ive_ptr) :: r,u,mf

      get_hmassflux = 0
      err = .false.
      r = getvar(RVAR, ndims,dims,stag,min,max,miss,dunits,
     >     disunits,dim_names,new)
      if(err) return
      u = getvar(UVAR, ndims,dims,stag,min,max,miss,dunits,
     >     disunits,dim_names,new)
      if(err) return
      
      mf=getmem(dims(1)*dims(2)*dims(4))
      call hor_mass_flux(dims(1),dims(2),dims(3),dims(4),
     &     %val(u),%val(r),%val(mf))
      get_hmassflux=mf
      dims(3)=1
      return
      end


*
*   Calculate vertical mass flux and horizental mass flux
*
       subroutine ver_mass_flux(mnx,mny,mnz,mnt,watrho,rho,ver_massflx)

      use default_module
       integer mnx,mny,mnz,mnt
       real watrho(mnx,mny,mnz,mnt),rho(mnx,mny,mnz,mnt),
     &         ver_massflx(1,mny,mnz,mnt)
       integer i,j,k,n

       ver_massflx = 0
       do n=1,mnt
          do j=1,mny
             do i=1,mnx
                do k=1,mnz
                   ver_massflx(1,j,k,n)=ver_massflx(1,j,k,n)+
     &                  rho(i,j,k,n)*(watrho(i,j,k,n))*dx
                enddo
             enddo
          enddo
       enddo
       return
       end
*
*
*
      subroutine hor_mass_flux(mnx,mny,mnz,mnt,uatrho,rho,hor_massflx)


      use default_module
       integer mnx,mny,mnz,mnt
       real uatrho(mnx,mny,mnz,mnt),rho(mnx,mny,mnz,mnt),
     &         hor_massflx(mnx,mny,1,mnt)
       integer i,j,k,n

       hor_massflx = 0
       do n=1,mnt
          do k=1,mnz
             do j=1,mny
                do i=1,mnx
                   hor_massflx(i,j,1,n)=hor_massflx(i,j,1,n)+
     &                  rho(i,j,k,n)*(uatrho(i,j,k,n))*dz
                enddo
             enddo
          enddo
          write(6,*)hor_massflx(5,1,1,n)
       enddo
       return
       end
