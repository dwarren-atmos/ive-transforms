C datelib.f
C it contains routines to manage the date:
c
c id=id_after_00(y,m,d) calculates the days after 1900/12/31
c date_from_days(y,m,d,id) returns the date when id is given.
c day=day_name(id) returns the name of the day (1=monday...)
c datestr(y,m,d,date) given y,m,d it returns YYMMDD
c hourstr(h,hour) given h, it returns HH
c hmstr(h,min,hour) given h,min it returns HH:MM
c
c------------------------------------------------------------------------------
C Testtreiber fuer Datum-Routinen:
c to run it, remove all left c down to the line end of remove.
c      program datetest
c      integer iy,im,id,n,i,iday    
c      external id_after_00
c      external day_name
c      integer day_name
c      print *,'year:'
c      read *,iy
c      print *,'month:'
c      read *,im
c      print *,'day:'
c      read *,id
c      print *,'date: ',iy,im,id
c      n=id_after_00(iy,im,id)
c      print *,'days after 1900: ',n
c      call date_from_days(iy,im,id,n)
c      iday=day_name(n)
c      print *,'date: ',iy,im,id,iday
c      print *,' '
c      do i=1,10
c      call date_from_days(iy,im,id,n+i)
c      iday=day_name(n+i)
c      print *,n+i,'-> date: ',iy,im,id,iday
c      enddo
c      end
c end of remove


      subroutine date_from_days(iy,im,id,nnn)
c     ======================================
c     returns the date iy,im,id for a given day number nn after 1900/12/31

      integer iy,im,id,nnn,nn,nm,nmm,idofy
      external idofy

      nn=nnn

c     calculate the year
      iy=1900
   10 continue
        nm=idofy(iy,12,31)
        if (nn.gt.nm) then
          nn=nn-nm
          iy=iy+1
        else
          goto 15
        endif
      goto 10

c     calculate the month
   15 continue
      im=1
      nmm=0
   20 continue
        nm=idofy(iy,im,99)
        if (nn.gt.nm) then
          im=im+1
          nmm=nm
        else
          goto 25
        endif
      goto 20

c     calculate the day
   25 continue
      id=nn-nmm

      return
      end

c--------------------------------------------------------------------------

      integer function id_after_00(iy,im,id)
c     ======================================
c     determines the number of days after 1900/12/31

      integer i,idofy,nn,iy,im,id
      external idofy

      nn=0

      if (iy.lt.1900) then
         if (iy.lt.100.) then
c david: add 1900, interpret as year since 1900:
            iy=iy+1900
         else
            print *, 'function value of id_after_00 not determined'
            stop
         endif
      endif
 
      do i=1900,iy
        if (i.lt.iy) then
          nn=nn+idofy(i,12,31)
        else
          nn=nn+idofy(iy,im,id)
        endif
      enddo

      id_after_00 = nn

      return
      end

c--------------------------------------------------------------------------

      integer function idofy(iyr,imon,iday)
c     =====================================
c     returns the day number of the date iyr/imon/iday
c
      integer  iyr,imon,iday,nr,im,id

      im = imon
      id = iday

c     in case the last of the month is to be returned
      if (id.gt.31) then
        im=im+1
        id=0
      endif

c     in case the last of the year is to be returned
      if (im.gt.12) then
        im=12
        id=31
      endif

      if (im.eq.1) then
        nr=id
      elseif (im.eq.2) then
        nr=id+31
      elseif (im.eq.3) then
        nr=id+59
      elseif (im.eq.4) then
        nr=id+90
      elseif (im.eq.5) then
        nr=id+120
      elseif (im.eq.6) then
        nr=id+151
      elseif (im.eq.7) then
        nr=id+181
      elseif (im.eq.8) then
        nr=id+212
      elseif (im.eq.9) then
        nr=id+243
      elseif (im.eq.10) then
        nr=id+273
      elseif (im.eq.11) then
        nr=id+304
      elseif (im.eq.12) then
        nr=id+334
      endif

c     schalttag:
      if ((mod(iyr,4).eq.0).and.(im.ge.3)) nr=nr+1

      idofy = nr

      return
      end

c--------------------------------------------------------------------------

      integer function day_name(id)
c     =====================================
c     returns the day name number of the day id after 1900/12/31
c     dayname=1 Monday
c     dayname=2 Tuesday etc. =7 Sunday
c     not working because of problems with mod.
c
      integer  id

      day_name = mod(id-2,7)+1
      
      return
      end


      subroutine datestr(y,m,d,date)
c-------------------------------------------------------------------------
c given y,m,d it returns YYMMDD
c david 941108
      integer   y,m,d
      character *(80) date
      character *2 dummy1,dummy2,dummy3
      character    dummy
c
c     function forward-def:
c      integer   strbeg,strend
c
      if (y.lt.10) then 
         write(dummy,'(i1)') y
         dummy1='0' // dummy
      else 
         write(dummy1,'(i2)') y 
      endif
      if (m.lt.10) then
         write(dummy,'(i1)') m
         dummy2='0' // dummy
      else
         write(dummy2,'(i2)') m
      endif
      if (d.lt.10) then
         write(dummy,'(i1)') d
         dummy3='0' // dummy
      else
         write(dummy3,'(i2)') d
      endif
      date=dummy1 // dummy2 // dummy3

      return
      end !datestr


      subroutine hmstr(h,m,hour)
c------------------------------------------------------------------------
c given h,min it returns HH:MM
c david 941108
      integer h,m
      character*(80) hour,min
      integer   strbeg,strend

      call int2str2(h,hour)
      call int2str2(m,min)
      hour=hour(strbeg(hour):strend(hour))//':'
     &     //min(strbeg(min):strend(min))

      return
      end !of hmstr


      subroutine shortyear(y)
c------------------------------------------------------------------------
c remove the leading 19 from 19XX
c david 941108
      integer y
      y=y-1900

      return
      end !end of shortyear


      subroutine int2str2(i,str)
c-------------------------------------------------------------------------
c given i, it returns a string with leading 0 if necessary:
c example: 5 -> '05' , 67 -> '67'
c david 941109
      integer  i
      character*(80) str
      character *2 dummy1
      character    dummy

      if (i.lt.10) then 
         write(dummy,'(i1)') i
         dummy1='0' // dummy
      else 
         write(dummy1,'(i2)') i
      endif
      str=dummy1

      return
      end !of int2str2
