
      subroutine units_smpfy(str1,str2,err)
*
*  units_smpfy takes a unit string as input and attempts to simplify
*    the string by cancelling like factors.  This routine does not
*    do any unit conversions before simplifying.
*
      character*(*) str1,str2

      integer top,bot,mult,div,list,lgth(2),swap(10),pcnt,ptr,iend,
     >        itmp,tptr,tcnt,ifac,char2int,k,expfac,exlist(10),fac,
     >        strbeg,strend,ptr_t,ptr_b,same
      character*(80) instr,outstr,u_top,u_bot,utmp 
      character*(10) ulist(2,40)
      logical op_sym,err

      err = .false.

      instr(1:len(instr)) = ' '
      outstr(1:len(outstr)) = ' '
      call setstr(instr,str1)
      iend = strend(instr)

      top = 1
      bot = 2
      lgth(top) = 0
      lgth(bot) = 0

      do k=1,40
        ulist(top,k)(1:len(ulist(top,k))) = ' '
        ulist(bot,k)(1:len(ulist(bot,k))) = ' '
      enddo

      mult = top
      div = bot
      list = mult

      pcnt = 0
      expfac = 1
      
      do k=1,10
        swap(k) = 0 
        exlist(k) = 1
      enddo

      ptr = 1

 111  if (ptr.le.iend) then

        if (instr(ptr:ptr).eq.'(') then

          pcnt = pcnt+1

          tptr = ptr+1
          tcnt = 0

 222      if (((instr(tptr:tptr).ne.')').and.(tptr.le.iend)).or.
     >                  ((tcnt.ne.0).and.(tptr.le.iend))) then
            if (instr(tptr:tptr).eq.'(') then
              tcnt = tcnt+1
            else if (instr(tptr:tptr).eq.')') then
              tcnt = tcnt-1
            end if
            if (tcnt.lt.0) then
              write(*,*) 'Warning:  error in simplifying units; ',
     >                            'tcnt less than zero.'
              err = .true.
              return
            end if
            tptr = tptr+1
            goto 222
          end if

          if (tptr.gt.iend) then
            write(*,*) 'Warning:  error in simplifying units; ',
     >                          'could not find end parenthesis.'
            err = .true. 
            return
          end if

          if ((tptr.lt.iend).and.(instr(tptr+1:tptr+1).eq.'^')) then
            exlist(pcnt) = char2int(instr(tptr+2:tptr+2))
            if (exlist(pcnt).le.0) then
              write(*,*) 'Warning:  error in simplifying units; ',
     >                       'could not find exponent after ^.'
              err = .true.
              return
            else
              expfac = expfac*exlist(pcnt)
            end if
          else
            exlist(pcnt) = 1
          end if

          if (list.eq.div) then
            itmp = mult
            mult = div
            div = itmp
            swap(pcnt) = 1
          else
            swap(pcnt) = 0
          end if

          list = mult
          ptr = ptr+1

        else if (instr(ptr:ptr).eq.')') then

 333      if (instr(ptr:ptr).eq.')') then
            if (pcnt.eq.0) then
              write(*,*) 'Warning:  error in simplifying units; ',
     >                       'unbalanced end parenthesis.'
              err = .true.
              return
            end if
            if (swap(pcnt).eq.1) then
              itmp = mult
              mult = div
              div = itmp
            end if
            expfac = expfac/exlist(pcnt)
            if ((ptr.lt.iend).and.(instr(ptr+1:ptr+1).eq.'^')) then
              ptr = ptr+2
            end if
            pcnt = pcnt-1
            ptr = ptr+1
            goto 333
          end if

          if (instr(ptr:ptr).eq.'*') then
            list = mult
            ptr = ptr+1
          else if (instr(ptr:ptr).eq.'/') then
            list = div
            ptr = ptr+1
          else if (ptr.le.iend) then
            write(*,*) 'Warning:  error in simplifying units; ',
     >                     'improper syntax in input unit string.'
            err = .true.
            return
          end if

        else

          tptr = ptr

 444      if ((.not.op_sym(instr(tptr:tptr))).and.(tptr.le.iend)) then
            tptr = tptr+1
            goto 444
          end if
          if (tptr.eq.ptr) then
            write(*,*) 'Warning:  error in simplifying units; ',
     >                       'null unit.'
            err = .true.
            return
          end if

          if (instr(tptr:tptr).eq.'^') then
            ifac = char2int(instr(tptr+1:tptr+1))
            if (ifac.le.0) then
              write(*,*) 'Warning:  error in simplifying units; ',
     >                        'could not find exponent after ^.'
              err = .true.
              return
            else
              fac = expfac*ifac
            end if
            do k=1,fac
              lgth(list) = lgth(list)+1
              ulist(list,lgth(list)) = instr(ptr:tptr-1)
            enddo 
            ptr = tptr+2
            if (instr(ptr:ptr).eq.'*') then
              list = mult
              ptr = ptr+1
            else if (instr(ptr:ptr).eq.'/') then
              list = div
              ptr = ptr+1
            else if ((ptr.le.iend).and.(instr(ptr:ptr).ne.')')) then
              write(*,*) 'Warning:  error in simplifying units; ',
     >                         'improper syntax in input string.'
              err = .true.
              return
            end if
          else
            do k=1,expfac
              lgth(list) = lgth(list)+1
              ulist(list,lgth(list)) = instr(ptr:tptr-1)
            enddo
            ptr = tptr
            if (instr(ptr:ptr).eq.'*') then
              list = mult
              ptr = ptr+1
            else if (instr(ptr:ptr).eq.'/') then
              list = div
              ptr = ptr+1
            else if ((ptr.le.iend).and.(instr(ptr:ptr).ne.')')) then
              write(*,*) 'Warning:  error in simplifying units; ',
     >                         'improper syntax in input string.'
              err = .true.
              return
            end if
          end if

        end if

        goto 111

      end if

      ptr_t = 1

 555  if (ptr_t.le.lgth(top)) then
 
        u_top = ulist(top,ptr_t)
        ptr_b = 1
        same = 0

 666    if ((same.eq.0).and.(ptr_b.le.lgth(bot))) then
          u_bot = ulist(bot,ptr_b)
          if (u_top(strbeg(u_top):strend(u_top)).eq.
     >                   u_bot(strbeg(u_bot):strend(u_bot))) then
            same = 1
            lgth(top) = lgth(top)-1 
            do k=ptr_t,lgth(top)
              ulist(top,k) =  ulist(top,k+1)
            enddo
            lgth(bot) = lgth(bot)-1
            do k=ptr_b,lgth(bot)
              ulist(bot,k) = ulist(bot,k+1)
            enddo
          else
            ptr_b = ptr_b+1
          end if
          goto 666
        end if

        if (same.eq.0) ptr_t = ptr_t+1

        goto 555

      end if

      if (lgth(top).gt.0) then
        outstr = ulist(top,1)
        do k=2,lgth(top)
          utmp = outstr(strbeg(outstr):strend(outstr))//'*'
          u_top = ulist(top,k)
          outstr = utmp(strbeg(utmp):strend(utmp))//
     >                  u_top(strbeg(u_top):strend(u_top))
        enddo
      else
        outstr(1:1) = '1'
      end if

      if (lgth(bot).gt.0) then
        do k=1,lgth(bot)
          utmp = outstr(strbeg(outstr):strend(outstr))//'/'
          u_bot = ulist(bot,k)
          outstr = utmp(strbeg(utmp):strend(utmp))//
     >                 u_bot(strbeg(u_bot):strend(u_bot))
        enddo
      end if

      if (outstr(1:len(outstr)).eq.'1') outstr(1:len(outstr)) = ' '

      call setstr(str2,outstr)

      return
      end

******************************************************************

      logical function op_sym(chr)

      character*(1) chr 

      if ((chr.eq.'*').or.(chr.eq.'/').or.(chr.eq.'^')
     >              .or.(chr.eq.')').or.(chr.eq.' ')) then
        op_sym = .true.
      else
        op_sym = .false.
      end if

      return
      end 
 
*****************************************************************
