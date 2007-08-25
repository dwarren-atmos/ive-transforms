!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  function rho_calc(pvar_name,tvar_name, ...)
!
!  Calculates Density
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(ive_ptr) function rho_calc(pvar_name,tvar_name,
     &                           ndims,dims,stag,min,max,
     &                           miss,dunits,disunits,dim_names)
      use default_module
      implicit none

      character(len=*), intent(in) :: pvar_name, tvar_name 
c 
      integer, intent(out) :: ndims,dims(4)
      real, intent(out) :: stag(4),min(4),max(4),miss
      character(len=*), intent(out) :: dunits,disunits,dim_names(4)
      logical new

      integer :: i, mnk
      integer(ive_ptr) :: pvar,tvar
      real :: the,pre,rho,r_val,slope,intercept

      rho_calc = 0
      mnk=(nx+1)*(ny+1)*nz*ntime

      pvar = getvar(pvar_name,
     &              ndims,dims,stag,min,max,miss,
     &              dunits,disunits,dim_names,new)
      !call convert_units(dunits,'Pa',slope,intercept)

      tvar = getvar(tvar_name,
     &              ndims,dims,stag,min,max,miss,
     &              dunits,disunits,dim_names,new)
  
      if (tvar.le.0 .or. pvar.le.0) then
         write(6,40) trim(pvar_name),trim(tvar_name)
         return
      end if
40    FORMAT('Error getting either ',A,' or ',A,'.')

      rho_calc = getmem(mnk) ; if(rho_calc.le.0) return

      do i=1,mnk
        the = r_val(%val(tvar),i) 
        pre = r_val(%val(pvar),i)
        pre = pre*100.0 ! convert from hPa to Pa
        !pre = pre*slope + interpcept
        rho = pref/(Rd*the)*(pre/pref)**(cv/cp)
        call s_val(%val(rho_calc),i,rho)
      end do

      call freemem(pvar) ; call freemem(tvar)

      dunits = 'kg/m/m/m' ; disunits = 'kg/m/m/m'

      return
      end

