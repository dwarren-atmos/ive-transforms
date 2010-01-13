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

      CHARACTER*80 driverstuff
      integer :: i, mnk
      integer(ive_ptr) :: pvar,tvar
      real :: the,pre,rho,r_val,slope,intercept, mod
      logical :: pnew, tnew
      rho_calc = 0

      pvar = getvar(pvar_name,
     &              ndims,dims,stag,min,max,miss,
     &              dunits,disunits,dim_names,pnew)
      !call convert_units(dunits,'Pa',slope,intercept)

      tvar = getvar(tvar_name,
     &              ndims,dims,stag,min,max,miss,
     &              dunits,disunits,dim_names,tnew)
  
      mnk=dims(1)*dims(2)*dims(3)*dims(4)

      if (tvar.le.0 .or. pvar.le.0) then
         write(6,40) trim(pvar_name),trim(tvar_name)
         return
      end if
40    FORMAT('Error getting either ',A,' or ',A,'.')

      rho_calc = getmem(mnk) ; if(rho_calc.le.0) return

      do i=1,mnk
        pre = r_val(%val(pvar),i)
!        pre = pre*100.0         ! convert from hPa to Pa
!     already Pa in wrf no hPa
!     convert to actual temperatureper wrf manual
        the = (r_val(%val(tvar),i)+300)
        !pre = pre*slope + interpcept
        rho = ((pre**(1-(Rd/cp)))*((1000**(Rd/cp))))/(Rd*the)
!        if(mod(i,100))write(6,*)'rho=',rho
        call s_val(%val(rho_calc),i,rho)
      end do
      new = .true.

      dunits = 'kg/m/m/m' ; disunits = 'kg/m/m/m'

      return
      end

