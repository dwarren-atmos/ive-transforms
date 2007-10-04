
c -----------------------------------------------------------------------------
c     HEADING-ROUTINE fuer EM, ECMWF Daten auf Druck und Modellflaechen
c -----------------------------------------------------------------------------
c     $Id: heading.f,v 1.1 1994/11/14 22:37:15 warren Exp $
c     $Log: heading.f,v $
c Revision 1.1  1994/11/14  22:37:15  warren
c Christoph Schaer's European Model transforms.
c
c Revision 1.8  1994/11/11  09:02:12  schaer
c Aenderung fuer IVE-3-2-beta
c Verbessertes Heading
c Startzeit neu auch auf Datenfile
c
c Revision 1.7  1994/03/07  08:50:53  schaer
c Aenderung bei der Bestimmung des runnames, cleanup.
c
c Revision 1.6  1994/02/18  07:43:32  schaer
c Anpassung an z-Koordinaten im Plot-Heading.
c
c Revision 1.5  1993/12/02  13:15:32  schaer
c Neue Version von IVE mit neuen Transformationen.
c
c Revision 1.4  1993/10/08  12:08:04  schaer
c Allow Simulation-Names with len(name)<>3.
c
c Revision 1.3  1993/10/07  15:30:31  schaer
c Genauere Angaben zum gewaehlten Level.
c
c Revision 1.2  1993/10/07  13:36:18  schaer
c Added new heading-routines from David Bresch.
c
c Revision 1.2  1993/09/15  14:27:06  dbr
c fheading.f created. Labels defined by user
c 3 poss:Forecast,Verification,Analysis
c
c Revision 1.1  1993/09/08  11:09:43  schaer
c Initial revision
c
c
c -----------------------------------------------------------------------------


      subroutine HEADING(where,which,line)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine allows the user to change the contents of the
c        heading title of the plot.
c        see the line: set timemode (DONE BY THE USER) 
c        about simple changes.
c     Arguments:
c        where  integer input 1:called from lab2dc (line contour)
c                             2:called from lab2ds (filled contour)
c                             3:called from drlineplt (line plot)
c        which  integer input 1:first line to be edited
c                             2:second line to be edited
c        line   string  i/o   content of first/second line.
c     History:
c        written by David N. Bresch 8.9.93 sma.ch
c        two smaller corrections , CS ETHZ, 7.10.93
c-----------------------------------------------------------------------
c
      include 'attributes.icl'
      include 'constants.icl'
c
c     Argument declarations.
c
      integer         where,which
      character *(80) line
c
c     Internal var decalrations.
c
      integer         ibeg,iend,i,j,timemode,type,lock(4)
      real            scaled_loc(4),plwmin(4),t,z
      character *(80) dtstr,rnname,zstr
      logical         error,prntz

c     Function type definition:
      integer         strbeg, strend
c     float to char and int to char
      character *(80) ftoa,itoa

c     set timemode (DONE BY USER via field=heading or automatically):
c     =0:Original IVE-headings
c     =1:Forecast: APPEND timestep to starttime
c     =2:Verification: ADD timestep to starttime and print VT
c     =3:Analyse: ADD timestep to starttime
c     =4:Forecast with day
c     =5:Analyse with day

c     type means:
c     =0:append :94/10/23/12+12 -> 94/10/93/12+12
c     =1:add timestep to starttime: 94/10/23/12+12 -> 94/10/24/00
c     =2:don't add or append

c     Neu:Entscheide aufrgrund des Wertes von headtyp:
      if (headtyp.eq.0) then
c     IVE-default:
        timemode=0
      elseif (headtyp.eq.1) then
c     Automatische Bestimmung:
c     Aenderung CS, 7.10.93: Automatische Bestimmung von timemode
        if ((dattyp.eq.52).or.(dattyp.eq.2)) then
c     Forecast
          timemode=1
          type=0
        else if ((dattyp.eq.58).or.(dattyp.eq.1)) then
c     Analyse
          timemode=3
          type=1
        else
          timemode=0
        endif
c     Ende Aenderung CS
      else if (headtyp.eq.2) then
c     Analyse:
          timemode=3
          type=1
      else if (headtyp.eq.3) then
c     Forecast:
          timemode=1
          type=0
      else if (headtyp.eq.4) then
c     Analyse day:
          timemode=4
          type=3
      else if (headtyp.eq.5) then
c     Forecast day:
          timemode=5
          type=2
      endif

      if (timemode.eq.0) return

c     get the timestep:
      call getrarr('scaled_loc',scaled_loc,4,error)
      call getrarr('plwmin',plwmin,4,error)
      call getiarr('lock',lock,4,error)
c     get the right t location:
      if (ntimes.eq.1) then
c       there is only one time on the file, get its value (Aenderung CS 7.10.)
        t=timeval(1)
      else if (lock(4).eq.1) then
c       t is locked, get loc value:            
        t=scaled_loc(4)
      else
c       get lower limit of the t range:
        t=plwmin(4)
      endif
c
c     Ergaenzung CS, 7.10.93
c     get the right z-location      
      if (trdims(3).eq.1) then
c       this is a 1d-field, don't print z
        prntz=.false.
      else if (lock(3).eq.1) then
c       z is locked, get loc value:            
        prntz=.true.
        z=scaled_loc(3)
      else
c       z is unlocked
        prntz=.false.
      endif
c     Ende Ergaenzung CS 

c     Aenderung CS: get runname from data-filenam rather than constants-filenam
      rnname=runname(strbeg(runname):strend(runname))
      ibeg = 1
      iend = strend(rnname)
      call upcase (rnname(1:iend), iend)
c     cut the _HST extension
      do i=2,iend-3
        if (rnname(i:i+3).eq.'_HST') then
          rnname=rnname(1:i-1)
          goto21
        endif
      enddo
 21   continue

c     Aenderung CS: get time from file-name for analysis files. This is a bug
c     fix for inconsistently written analysis-cdf-files. These files
c     are recognized through a leading digit in the simulation name
      if ((line(1:1).ge.'0').and.(line(1:1).le.'9')) then
c       find group associated with time
        ibeg=0
        iend=0
        do i=1,strend(line)-1
          if ((ibeg.eq.0).and.(line(i:i).eq.'_')
     &      .and.(line(i+1:i+1).ge.'0').and.(line(i+1:i+1).le.'9')) then
            ibeg=i+1
            iend=strend(line)
            do j=ibeg,strend(line)
              if ((line(j:j).lt.'0').or.(line(j:j).gt.'9')) then
                iend=j-1
                goto 25
              endif
            enddo
            goto 25
          endif
        enddo
 25     if (iend.gt.ibeg) then
          t=0.
          do i=ibeg,iend
            t=t*10.+real(ichar(line(i:i))-ichar('0'))
          enddo
          starth=0
        endif
      endif
c
c     calculate the new date when timestep is ADDED to starttime 
c     (type=1) or just APPEND the timestep (type=0):
c     (the addition of 0.0 is necessary for type match)
c      if (type.lt.2) then
        call udtime(type,starty,startm,startd,starth+0.0,t,dtstr)
c      endif
c
c     decide, if first or second line:
      if (which.eq.1) then
c     first line to be edited:
        if ((timemode.eq.1).or.(timemode.eq.3)) then
c     Forecast or Analyse (Unterschied in type):
          line=rnname(strbeg(rnname):strend(rnname))//' '
     &      //dtstr(strbeg(dtstr):strend(dtstr))
        else if (timemode.eq.2) then
c     Verification time, write VT
          line=rnname(strbeg(rnname):strend(rnname))//' VT: '
     &      //dtstr(strbeg(dtstr):strend(dtstr))
        else if (timemode.eq.4) then
c     Analyse mit day:
          line=rnname(strbeg(rnname):strend(rnname))//' '
     &      //dtstr(strbeg(dtstr):strend(dtstr))
        else if (timemode.eq.5) then
c     Forecast mit day:
          line=rnname(strbeg(rnname):strend(rnname))//' '
     &      //dtstr(strbeg(dtstr):strend(dtstr))
        endif          
      else
c     second line to be edited

c     Aenderung CS (7.10.83)
c     allgemeine Behandlung der Level-Angaben
c     reduziere gegenwaertigen String um "at time ..."
        ibeg=1
        iend=strend(line)
        do i=1,strend(line)-3
          if ((line(i:i).eq.' ').and.(line(i+1:i+1).eq.'a')
     &       .and.(line(i+2:i+2).eq.'t')
     &       .and.(line(i+3:i+3).eq.' ')) iend=i-1
        enddo
        iend=max(ibeg,iend)
        line=line(ibeg:iend)
              
c     verfasse Level-Angabe
        if (prntz) then   
          if (abs(z-real(nint(z))).lt.0.01) then
            zstr=itoa(nint(z))
          else
            zstr=ftoa(z)
        endif
                if (disptype.eq.0) then
                  line=line(ibeg:iend) // ' at K=' // 
     &                         zstr(strbeg(zstr):strend(zstr))
                else if (disptype.eq.1) then
                  line=line(ibeg:iend) // ' at ' // 
     &                         zstr(strbeg(zstr):strend(zstr)) // 'hPa'
                else if (disptype.eq.2) then
                  line=line(ibeg:iend) // ' at ' // 
     &                         zstr(strbeg(zstr):strend(zstr)) // 'K'
                else if (disptype.eq.3) then
                  line=line(ibeg:iend) // ' at ' // 
     &                         zstr(strbeg(zstr):strend(zstr)) // 'km'
                else 
                  line=line(ibeg:iend) // ' at ' // 
     &                         zstr(strbeg(zstr):strend(zstr)) // '?'
                endif
              endif
c             ende Aenderung CS

            endif
      end


      subroutine udtime(type,y,m,d,h,t,dtstr)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine updates date and time, depending on type and
c        returns a string containing it.
c        type   integer input 0:append timestep time separated by '+'
c                             1:add timestep
c                             2:fc with day
c                             3:ana with day
c        y      integer input year
c        m      integer input month
c        d      integer input day
c        h      real    input hour
c        t      real    input actual hour offset
c        dtstr  character output the string
c     History:
c        written by David N. Bresch 8.9.93 sma.ch
c-----------------------------------------------------------------------
c
c     Variables:
c
      integer         type,y,m,d
      real            t,h
      character *(80) dtstr,daystr
      character *(80) time
c
c     Function type definition:
      integer         strbeg, strend
      character *(80) ftoa
c
c     Internal vars:
c
      integer         ny,nm,nd,ml,ltot
      real            nh
c
c     convert time increment to string:
      if (t.eq.0.0) then
         time='0.0'
      else
         time=ftoa(t)
      endif
c     check, if time has 2 digits before comma, add one, if not:
      if (t.lt.10.0) then
         time='0'//time(strbeg(time):strend(time))
      endif
c
c     decide, if adding or appending:
      if (type.eq.0) then
c     append time:
         daystr=' '
         call datestring(y,m,d,h,daystr,dtstr)
c
c        alte version DB
c               dtstr=dtstr(strbeg(dtstr):strend(dtstr))//'Z'//
c     &               '+'//time(strbeg(time):strend(time))//'h'
               
c        neue version CS: Eliminiert 'h' und 'z'. Im weiteren werden
c        Kommastellen unterdrueckt falls nicht erforderlich
         if (abs(t-real(nint(t))).lt.0.05) then
            dtstr=dtstr(strbeg(dtstr):strend(dtstr))//
     &            '+'//time(strbeg(time):strend(time)-2)
         else
            dtstr=dtstr(strbeg(dtstr):strend(dtstr))//
     &            '+'//time(strbeg(time):strend(time))
         endif
c
       else
c           add time:
               ny=y
               nm=m
               nd=d
               nh=h+t
c              process hours:
 10            if (nh.ge.24.0) then
                        nh=nh-24.0
                        nd=nd+1
                        go to 10
                     endif
c              get the lenght of the month:
               call getmlength(nm,ml,ltot)
 20            if (nd.gt.ml) then
                        nd=nd-ml
                        nm=nm+1
c                       get new month length:
                        call getmlength(nm,ml,ltot)
                        go to 20
                     endif
c              process month:    
 30            if (nm.ge.12) then
                        nm=nm-12
                        ny=ny+1
                        go to 30
                     endif 
c              decide, which format:
               if (type.eq.1) then
c           analyse:
                  daystr=' '
                  call datestring(ny,nm,nd,nh,daystr,dtstr)
               else if (type.eq.2) then
c           forecast with day:
                  call new_daystring(ny,nm,nd,nh,daystr)
                  !daystr='hallo fc'
                  call datestring(ny,nm,nd,nh,daystr,dtstr)
                  if (abs(t-real(nint(t))).lt.0.05) then
                     dtstr=dtstr(strbeg(dtstr):strend(dtstr))//
     &             ' (+'//time(strbeg(time):strend(time)-2)//')'
                  else
                     dtstr=dtstr(strbeg(dtstr):strend(dtstr))//
     &             ' (+'//time(strbeg(time):strend(time))//')'
                  endif  
               else if (type.eq.3) then
c           analyse with day:
                  call new_daystring(ny,nm,nd,nh,daystr)
                  !daystr='hallu ana'
                  call datestring(ny,nm,nd,nh,daystr,dtstr)                  
                  if (abs(t-real(nint(t))).lt.0.05) then
                     dtstr=dtstr(strbeg(dtstr):strend(dtstr))
                  else
                     dtstr=dtstr(strbeg(dtstr):strend(dtstr))
                  endif  
               else
c           wie type=1 at the moment.
                  daystr=' '
                  call datestring(ny,nm,nd,nh,daystr,dtstr)
               endif
c              dtstr=dtstr(strbeg(dtstr):strend(dtstr))//'Z'
               dtstr=dtstr(strbeg(dtstr):strend(dtstr))
            endif

      return
      end !udtime


      subroutine datestring(y,m,d,h,dstr,str)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine converts date and time to yy/mm/dd tt
c        or to yy/mm/dd day tt if daystr is given.
c     Arguments:
c        y      integer input year
c        m      integer input month
c        d      integer input day
c        h      real    input hour
c        dtstr  character output the string
c     History:
c        written by David N. Bresch 8.9.93 sma.ch
c-----------------------------------------------------------------------
c
c     Variables:
c
      integer         y,m,d
      real            h
      character *(80) str,dstr
c
c     Function type definition:
      integer         strbeg, strend
      character *(80) ftoa,itoa
c
c     Internal vars:
      character *(80) year
      character *(80) month
      character *(80) day
      character *(80) hour
c
c     convert integer->string
      year=itoa(y)
      if (strend(year).lt.2) then
                year='0'//year(strbeg(year):strend(year))
            endif
      month=itoa(m)
      if (strend(month).lt.2) then
                month='0'//month(strbeg(month):strend(month))
            endif

      day=itoa(d)
      if (strend(day).lt.2) then
                day='0'//day(strbeg(day):strend(day))
            endif
c     If there is no time, set it to zero:
      if (h.eq.0.0) then
               hour='0'
            else
               hour=ftoa(h)
            endif
c     check, if time has 2 digits before comma, add one, if not:
      if (h.lt.10.0) then
                hour='0'//hour(strbeg(hour):strend(hour))
            endif

c     remove unnecessary digits
      if (abs(h-real(nint(h))).lt.0.01) then
        hour=hour(1:2)
      endif
c
c here, the string format yy/mm/dd hh is defined:
c format yy/mm/dd hh
      if (dstr(strbeg(dstr):strend(dstr)).ne.' ') then
         str=year(strbeg(year):strend(year))//
     &      '/'//month(strbeg(month):strend(month))//
     &      '/'//day(strbeg(day):strend(day))//
     &      ' '//dstr(strbeg(dstr):strend(dstr))//
     &      ' '//hour(strbeg(hour):strend(hour))
      else
c without daystring:
         str=year(strbeg(year):strend(year))//
     &      '/'//month(strbeg(month):strend(month))//
     &      '/'//day(strbeg(day):strend(day))//
     &      ' '//hour(strbeg(hour):strend(hour))
      endif
c
      return
      end


      subroutine new_daystring(y,m,d,h,daystr)
c---------------------------------------------------------------------------
c     purpose:
c        this routine returns the day to a given date, using
c        an eternal calendar.
c     arguments:
c        y      integer input year
c        m      integer input month
c        d      integer input day
c        h      real    input hour
c        daytstr  character output the string as 'Mon' 
c
c     history:
c        written by david n. bresch 941107
c        it replaces the older daystring routine
c-----------------------------------------------------------------------
c
c     Variables:
c
      integer         y,m,d
      real            h
      character *(80) daystr
c
      integer         id,day
      integer         id_after_00
      integer         day_name
      CHARACTER*3 CDAY(7)
c
      CDAY(1)='Mon'
      CDAY(2)='Tue'
      CDAY(3)='Wed'
      CDAY(4)='Thu'
      CDAY(5)='Fri'
      CDAY(6)='Sat'
      CDAY(7)='Sun'

      id=id_after_00(y,m,d)
      day=day_name(id)
      
      daystr=CDAY(day)

      return
      end !of new_daystring


      subroutine daystring(y,m,d,h,daystr)
c---------------------------------------------------------------------------
c     purpose:
c        this routine returns the day to a given date
c     arguments:
c        y      integer input year
c        m      integer input month
c        d      integer input day
c        h      real    input hour
c        daytstr  character output the string as ' Mon ' 
c
c     history:
c        written by david n. bresch 941026
c-----------------------------------------------------------------------
c
c     Variables:
c
      integer         y,m,d
      real            h
      character *(80) daystr
c
c     Internal vars: (c for current or actual)
      integer         cy,cm,cd,cwd,inc,day,l,ltot,ltot2
      integer         dy,dd
      CHARACTER*3 CDAY(7)
c
      CDAY(1)='Mon'
      CDAY(2)='Tue'
      CDAY(3)='Wed'
      CDAY(4)='Thu'
      CDAY(5)='Fri'
      CDAY(6)='Sat'
      CDAY(7)='Sun'
c
c     get the current time and the day from the system:
c     the current  time is referred as ctime
c     the time passed to this routine as time.
c     UNIX UNIX UNIX UNIX UNIX
      call system ('/usr/bin/rm -f fort.98')
c     get the current weekday:
      call system ("date '+%w' > fort.98")
      REWIND(98)
      READ(98,'(I1)') cwd
      CLOSE(98)
      call system ('/usr/bin/rm -f fort.98')
c     get the current year:
      call system ("date '+%y' > fort.98")
      REWIND(98)
      READ(98,'(I2)') cy
      CLOSE(98)
      call system ('/usr/bin/rm -f fort.98')
c     get the current month:
      call system ("date '+%m' > fort.98")
      REWIND(98)
      READ(98,'(I2)') cm
      CLOSE(98)
      call system ('/usr/bin/rm -f fort.98')
c     get the current day:
      call system ("date '+%d' > fort.98")
      REWIND(98)
      READ(98,'(I2)') cd
      CLOSE(98)
      call system ('/usr/bin/rm -f fort.98')
c     UNIX UNIX UNIX UNIX 
c     Calculate the difference of days between system time and
c     the time passed to this routine:
      dd=d-cd
      dy=y-cy
c     get the lenght of every month between ctime and time:
c     difference of full years:
      call getmlength(m,l,ltot)      
      call getmlength(cm,l,ltot2)
      inc=ltot-ltot2
      inc=365*dy+inc+dd+cwd
      day=mod(inc,7)
  20   if (day.lt.1) then !avoid negative days
         day=day+7
         go to 20
      endif
 30   if (day.gt.7) then !avoid to large a day
         day=day-7
         go to 30
      endif
      if ((day.gt.0).and.(day.lt.8)) then
         daystr=CDAY(day)
      else
         daystr='error'
      endif

      return
      end !of daystring


      subroutine getmlength(m,l,ltot)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine returns the length of a month in days.
c        (doesn't care about Feb 29)
c     Arguments:
c        m      integer input month
c        l      integer output number of days
c        ltot   integer output number of days to begin of month since jan.
c     History:
c        written by David N. Bresch 8.9.93 sma.ch
c-----------------------------------------------------------------------
c
c     Variables:
c
      integer         m,l,ltot
c
      if (m.eq.1) then
c           Jan
               l=31
               ltot=0
            elseif (m.eq.2) then
c           Feb
               l=28
               ltot=31
            elseif (m.eq.3) then
c           Mar
               l=31
               ltot=59
            elseif (m.eq.4) then
c           Apr
               l=30
               ltot=90
            elseif (m.eq.5) then
c           Mai
               l=31
               ltot=120
            elseif (m.eq.6) then
c           Jun
               l=30
               ltot=151
            elseif (m.eq.7) then
c           Jul
               l=31
               ltot=181
            elseif (m.eq.8) then
c           Aug
               l=31
               ltot=212
            elseif (m.eq.9) then
c           Sep
               l=30
               ltot=243
            elseif (m.eq.10) then
c           Oct
               l=31
               ltot=273
            elseif (m.eq.11) then 
c           Nov
               l=30
               ltot=304
            elseif (m.eq.12) then
c           Dec
               l=31
               ltot=334
            endif
c
      return
      end
