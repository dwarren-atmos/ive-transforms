      subroutine get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
*
*  get_args searches the command string fld for the arguments of a
*    function.  This routine returns an array of arguments for each
*    type (character, integer, real).  For the purposes of this routine,
*    the presence of a decimal distinguishes a real argument from
*    an integer argument.
*
      integer, parameter :: argmax=10
      integer :: iarg(argmax),rcnt,icnt,ccnt
      real :: rarg(argmax)
      character(len=*) :: fld,carg(argmax)
      logical :: err

      integer i1,iend,k,bcnt,str2int
      real str2real
      character*80 cstr
      logical delimit,charstr,realstr,lendbrk,lbegbrk

      err = .false.
      carg(:) = ' ' ; iarg(:) = -1 ; rarg(:) = -999.
      rcnt = 0 ; icnt = 0 ; ccnt = 0

      call setstr(cstr,fld)
      iend = index(cstr,char(0)) - 1

      if (cstr(iend:iend).ne.']') then
        lendbrk=.false.
      end if
*
*  Find open bracket.
*
      i1 = 0
      k = 1
 543  if ((cstr(k:k).ne.'[').and.(k.lt.iend)) then
        k = k+1
        goto 543
      end if

      if (k.eq.iend) then
        lbegbrk=.false.
      end if

      if(lbegbrk.and..not.lendbrk) then
        write(*,*) 'Error:  no ending bracket.'
        err = .true.
        return
      else if(.not.lbegbrk.and.lendbrk) then
        write(*,*) 'Error:  no opening bracket.'
        err = .true.
        return
      end if
*
*  Find arguments.  Check to determine type and assign to the 
*    appropriate array.
*
 643  if (k.lt.iend) then
     
        bcnt = 0
        i1 = k
        k = k+1
        
 754    if (((.not.delimit(cstr(k:k))).and.(k.lt.iend)).or.
     >               ((bcnt.ne.0).and.(k.lt.iend))) then

          if (cstr(k:k).eq.'[') then
            bcnt = bcnt+1
          else if (cstr(k:k).eq.']') then
            bcnt = bcnt-1
          end if
     
          if (bcnt.lt.0) then
            write(*,*)  'Error:  unbalanced brackets.'
            err = .true.
            return
          end if

          k = k+1        
 
          goto 754

        end if

        if (k.eq.i1+1) then
          write(*,*) 'Error:  null argument.'
          err = .true.
          return
        else if (bcnt.ne.0) then
          write(*,*) 'Error:  unbalanced brackets.'
          err = .true.
          return
        end if


        if (charstr(cstr(i1+1:k-1))) then

          ccnt = ccnt+1
          if (ccnt.gt.argmax) then
            write(*,*) 'Ooops!  Too many string arguments!'
            err = .true.
            return
          end if
          call setstr(carg(ccnt),cstr(i1+1:k-1))

        else if (realstr(cstr(i1+1:k-1))) then

          rcnt = rcnt+1
          if (rcnt.gt.argmax) then
            write(*,*) 'Ooops!  Too many real arguments!'
            err = .true.
            return
          end if
          rarg(rcnt) = str2real(cstr(i1+1:k-1),err)

          if (err) return

        else 

          icnt = icnt+1
          if (icnt.gt.argmax) then
            write(*,*) 'Ooops!  Too many integer arguments!'
            err = .true.
            return
          end if
          iarg(icnt) = str2int(cstr(i1+1:k-1),err)
          if (err) return

        end if

        goto 643

      end if

      return 
      end

**********************************************************************

      integer function str2int(str,err)
*
*  str2int attempts to convert the string str to the corresponding
*    integer.
*
      character*(*) str
      logical err

      integer strbeg,strend,ibeg,iend,lgth,k,ichar,istr
      integer char2int
      logical minus

      if(trim(str).eq.'.FALSE.') then
        str2int = 0 ; return
      elseif(trim(str).eq.'.TRUE.') then 
        str2int = 1 ; return
      endif

      ibeg = strbeg(str) ; iend = strend(str)

      if (str(strbeg(str):strbeg(str)).eq.'-') then
        minus = .true.
        ibeg = ibeg+1
      else
        minus = .false.
      end if

      lgth = iend-ibeg+1
      istr = 0

      do k=ibeg,iend

        ichar = char2int(str(k:k))

        if (ichar.lt.0) then
          write(*,*) 'Error:  ',str(ibeg:iend),' is not an integer.'
          err = .true.
          str2int = -1
          return
        end if

        istr = istr + ichar*(10**(lgth-(k-ibeg+1)))

      enddo 

      if (minus) then
        str2int = -istr
      else
        str2int = istr
      end if

      return
      end

*******************************************************************

      integer function char2int(chr)
*
*  char2int attempts to convert the character char to the corresponding
*    integer
*
      character*1 chr
      
      select case(trim(chr))
      case('0')
        char2int = 0
      case('1')
        char2int = 1
      case('2')
        char2int = 2
      case('3')
        char2int = 3
      case('4')
        char2int = 4
      case('5')
        char2int = 5
      case('6')
        char2int = 6
      case('7')
        char2int = 7
      case('8')
        char2int = 8
      case('9')
        char2int = 9
      case default
        char2int = -1
      end select

      return
      end

***********************************************************************

      real function str2real(str,err)
*
*  str2real attempts to convert the string str to the corresponding
*    real.
*
      character*(*) str
      logical err

      integer strbeg,strend,ibeg,iend,idec,k,ichar,char2int
      real rstr
      logical minus
  
      ibeg = strbeg(str)
      iend = strend(str)

      rstr = 0.

      if (str(ibeg:ibeg).eq.'-') then
        minus = .true.
        ibeg = ibeg+1
      else
        minus = .false.
      end if

      idec = 0
      k = ibeg
 646  if ((idec.eq.0).and.(k.lt.iend)) then
        if (str(k:k).eq.'.') then
          idec = k-ibeg+1
        else
          k = k+1
        end if  
        goto 646
      end if

      if (k.eq.iend) idec = iend-ibeg+1

      do k = ibeg,iend

        if ((k-ibeg+1).ne.idec) then

          ichar = char2int(str(k:k))

          if (ichar.lt.0) then
            write(*,*) 'Error:  ',str(ibeg:iend),' is not a ',
     >                  'valid argument.' 
            err = .true.
            str2real = -999.
            return
          end if

          if ((k-ibeg+1).lt.idec) then
 
            rstr = rstr + real(ichar)*(10.**(idec-(k-ibeg+2)))

          else if ((k-ibeg+1).gt.idec) then

            rstr = rstr + real(ichar)*(10.**(idec-(k-ibeg+1))) 
    
          end if

        end if

      enddo

      if (minus) then
        str2real = -rstr
      else
        str2real = rstr
      end if

      return
      end

***********************************************************************

      logical function delimit(chr)

      character*1 chr

      if (chr.eq.':') then
        delimit = .true.
      else if (chr.eq.';') then
        delimit = .true.
      else
        delimit = .false.
      end if

      return
      end

************************************************************************

      logical function charstr(str)

      character*(*) str

      integer strbeg,strend,k
      logical notnumber

      charstr = .false.
      k = strbeg(str)

      if(trim(str).eq.'.FALSE.' .or. trim(str).eq.'.TRUE.') return

 9867 if (k.le.strend(str)) then
     
        if (notnumber(str(k:k))) then
          charstr = .true.
          return
        else
          k = k+1
        end if

        goto 9867 
 
      end if

      return
      end

**********************************************************************

      logical function notnumber(chr)

      character*1 chr

      if ((chr.ne.' ').and.(chr.ne.'1').and.(chr.ne.'2').and.
     >     (chr.ne.'3').and.(chr.ne.'4').and.(chr.ne.'5').and.
     >     (chr.ne.'6').and.(chr.ne.'7').and.(chr.ne.'8').and.
     >     (chr.ne.'9').and.(chr.ne.'0').and.(chr.ne.'.').and.
     >     (chr.ne.'-')) then
        
        notnumber = .true.

      else
 
        notnumber = .false.

      end if

      return
      end

*************************************************************************88

      logical function realstr(str)      

      character*(*) str

      integer k,strbeg,strend

      realstr = .false.
      k = strbeg(str)

      if(trim(str).eq.'.FALSE.' .or. trim(str).eq.'.TRUE.') return

 98   if (k.le.strend(str)) then

        if (str(k:k).eq.'.') then
          realstr = .true.
          return
        else
          k = k+1
        end if

        goto 98

      end if

      return
      end


      



        
        

