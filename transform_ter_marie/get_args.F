
      subroutine get_args(fld,rarg,iarg,carg,rcnt,icnt,ccnt,err)
*
*  get_args searches the command string fld for the arguments of a
*    function.  This routine returns an array of arguments for each
*    type (character, integer, real).  For the purposes of this routine,
*    the presence of a decimal distinguishes a real argument from
*    an integer argument.
*
      integer iarg(10),rcnt,icnt,ccnt
      real rarg(10)
      character*(*) fld,carg(10)
      logical err

      integer i1,iend,k,bcnt,str2int
      real str2real
      character*80 cstr
      logical delimit,charstr,realstr

      err = .false.
      do k=1,10
        carg(k) = ' '
        iarg(k) = -1
        rarg(k) = -999.
      enddo
      rcnt = 0
      icnt = 0
      ccnt = 0

      call setstr(cstr,fld)
      iend = index(cstr,char(0)) - 1
      if (cstr(iend:iend).ne.']') then
        write(*,*) 'Error:  no ending bracket.'
        err = .true.
        return
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
          if (ccnt.gt.10) then
            write(*,*) 'Ooops!  Too many string arguments!'
            err = .true.
            return
          end if
          call setstr(carg(ccnt),cstr(i1+1:k-1))

        else if (realstr(cstr(i1+1:k-1))) then
 
          rcnt = rcnt+1
          if (rcnt.gt.10) then
            write(*,*) 'Ooops!  Too many real arguments!'
            err = .true.
            return
          end if
          rarg(rcnt) = str2real(cstr(i1+1:k-1),err)
          if (err) return

        else 

          icnt = icnt+1
          if (icnt.gt.10) then
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

      ibeg = strbeg(str)
      iend = strend(str)
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

      str2int = istr

      return
      end

*******************************************************************

      integer function char2int(chr)
*
*  char2int attempts to convert the character char to the corresponding
*    integer
*
      character*1 chr

      if (chr.eq.'1') then
        char2int = 1
      else if (chr.eq.'2') then
        char2int = 2
      else if (chr.eq.'3') then
        char2int = 3
      else if (chr.eq.'4') then
        char2int = 4
      else if (chr.eq.'5') then
        char2int = 5
      else if (chr.eq.'6') then
        char2int = 6
      else if (chr.eq.'7') then
        char2int = 7
      else if (chr.eq.'8') then
        char2int = 8
      else if (chr.eq.'9') then
        char2int = 9
      else if (chr.eq.'0') then
        char2int = 0
      else
        char2int = -1
      end if

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
  
      ibeg = strbeg(str)
      iend = strend(str)

      rstr = 0.

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

      if (k.eq.iend) idec = iend

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

      str2real = rstr

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
     >     (chr.ne.'9').and.(chr.ne.'0').and.(chr.ne.'.')) then
        
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


      



        
        

