      subroutine run_sub(name)

      use default_module
      use misc_module

      character(len=*), intent(in) :: name

      integer :: ibeg, iend, iarg(10), rcnt, icnt, ccnt
      real :: rarg(10), r_val
      character(len=80) :: fld, carg(10), fld_cut
      logical :: err

      call setstr(fld,name)
      iend = index(fld,'[',.true.) - 1

      if(iend.lt.0) then
        fld_cut = trim(fld)
        !print '(3A,I2)','lt0: ',trim(fld_cut),':',iend
      else
        fld_cut = fld(1:iend)
        !print '(3A,I2)','ge0: ',trim(fld_cut),':',iend
      end if

      call get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
      if (err) return
      !print '(A,1x,3(I1,1x))',trim(fld_cut),rcnt,icnt,ccnt

      select case (trim(fld_cut))

c      case ('SNDOUT_ENKF')
c          call calc_field_message(fld_cut,rcnt,icnt,ccnt,4,1,0,err)
c         if (err) return

c         call sound_out(rarg(1),rarg(2),rarg(3),rarg(4),0.0,
c     &                iarg(1),1)
c         return

c      case ('SNDOUT_COAMPS_IDEAL')
c        call calc_field_message(fld_cut,rcnt,icnt,ccnt,4,1,0,err)
c        if (err) return

c        call sound_out(rarg(1),rarg(2),rarg(3),0.0,rarg(4),
c     &                iarg(1),2)
c        return

!      case ('TRANSLATE') ! translate the horizontal spatial domain
!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,2,0,0,err)
!        if(err) return

!        call translate_domain(rarg(1),rarg(2))
!        return

!      case ('TRANS_TIME','trans_time') ! translate the time domain
!!        call calc_field_message(fld_cut,rcnt,icnt,ccnt,1,0,0,err)
!!        if(err) return
!        call translate_time(rarg(1))

!        return

!      case ('TRANS_ON') ! Turn on vertical transform
!        call trans_switch(1)
!        return

!      case ('TRANS_OFF') ! Turn off vertical transform
!        call trans_switch(0)
!        return

      case ('TEST','test')
        print *,'THIS IS A TEST'
        return

      case default

        print *,'Unknwon routine ',trim(fld_cut),'.'
        return

      end select

      return
      end 
