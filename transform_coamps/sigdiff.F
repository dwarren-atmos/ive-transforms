!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   function sigdiff(varnam1,varnam2,stdnam1,stdnam2,conflev ...)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(ive_ptr) function sigdiff(varnam1,varnam2,stdnam1,
     &                                  stdnam2,conflev,N1,N2,
     &                                  ndims,dims,stag,min,max,miss,
     &                                  dunits,disunits,dim_names)

      USE DEFAULT_MODULE

      character(len=*), intent(in) :: varnam1,varnam2,stdnam1,stdnam2
      real, intent(in) :: conflev
      integer, intent(in) :: N1,N2

      integer, intent(inout) :: ndims,dims(4)
      real, intent(inout) :: stag(4),min(4),max(4),miss
      character(len=*), intent(inout):: dunits,disunits,dim_names(4)
      
      integer :: dims1(4), dims2(4)
      character(len=80) :: dunits1,disunits1,dim_names1(4)
     &                   , dunits2,disunits2,dim_names2(4)

      integer :: i, ii, mnk, ix, iy, iz, it
      integer(ive_ptr) :: var1,var2,std1,std2
      real :: r_val, tcrit, x1, x2, s1, s2, tval, sig
      logical :: new

      ! external function declaration
      real :: tinv
       
      sigdiff=0

      ! get the 2 variables and check for equal dims and units.
      print *,'GETTING ',trim(varnam1)
      var1 = getvar(varnam1,ndims,dims1,stag,min,max,miss,dunits1,
     &              disunits1,dim_names1,new)
      if(var1.eq.0) then
        print *,'Cannot get ',trim(varnam1)
        return
      end if

      print *,'GETTING ',trim(varnam2)
      var2 = getvar(varnam2,ndims,dims2,stag,min,max,miss,dunits2,
     &              disunits2,dim_names2,new)
      if(var2.eq.0) then
        print *,'Cannot get ',trim(varnam2)
        return
      end if

      do i=1,3
        if(dims1(i).ne.dims2(i)) then
          write(6,100) trim(dim_names1(i)),'dimension',trim(varnam1)
     &                ,trim(dim_names2(i)),'dimension',trim(varnam2)
          return
        end if
      end do

      if(dunits1.ne.dunits2) then
        write(6,100) trim(dunits1),'units',trim(varnam1)
     &              ,trim(dunits2),'units',trim(varnam2)
        return
      end if

      print *,'GETTING ',trim(stdnam1)
      std1 = getvar(stdnam1,ndims,dims1,stag,min,max,miss,dunits1,
     &             disunits1,dim_names1,new)
      if(std1.eq.0) then
        print *,'Cannot get ',trim(stdnam1)
        return
      end if

      print *,'GETTING ',trim(stdnam2)
      std2 = getvar(stdnam2,ndims,dims2,stag,min,max,miss,dunits2,
     &              disunits2,dim_names2,new)

      if(std2.eq.0) then
        print *,'Cannot get ',trim(stdnam2)
        return
      end if

      do i=1,3
        if(dims1(i).ne.dims2(i)) then
          write(6,100) trim(dim_names1(i)),'dimension',trim(stdnam1)
     &                ,trim(dim_names2(i)),'dimension',trim(stdnam2)
          return
        end if
      end do

      if(dunits1.ne.dunits2) then
        write(6,100) trim(dunits1),'units',trim(stdnam1)
     &              ,trim(dunits2),'units',trim(stdnam2)
        return
      end if

100   format(A,1x,A,' for ',A,' is not the same as ',
     &       A,1x,A,' for ',A,'.')

      dims(:)=dims1(:) ; dim_names(:)=dim_names2(:)
      ix=dims(1) ; iy=dims(2) ; iz=dims(3) ; it=ntime
      mnk = it*iz*iy*ix

      sigdiff = getmem(mnk)
      if(sigdiff.le.0) then
        print *,'Error: cannont allocate nesecary memory in sigdiff'
        return
      end if

      !tcrit = tinv(conflev,N1+N2-2)
       tcrit=2.1009 ! alpha=0.05, nu=18

      do i=1,mnk 

        x1 = r_val(%val(var1),i) ; x2 = r_val(%val(var2),i)
        s1 = r_val(%val(std1),i) ; s2 = r_val(%val(std2),i)

        sig = sqrt((real(N1)*s1*s1 + real(N2)*s2*s2)/real(N1+N2-2))
        tval = abs(x1-x2)/(sig*sqrt(1/real(N1) + 1/real(N2)))

        if(tval < tcrit) then
          call s_val(%val(sigdiff),i,-1.0)
        else
          call s_val(%val(sigdiff),i,1.0)
        end if
      end do


      dunits='' ; disunits=''

      return
      end
