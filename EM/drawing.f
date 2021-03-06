c-----------------------------------------------------------------------
c     This file contains auxiliary plotting modules
c     $Id: drawing.f,v 1.1 1994/11/14 22:37:04 warren Exp $
c     history
c       $Log: drawing.f,v $
c Revision 1.1  1994/11/14  22:37:04  warren
c Christoph Schaer's European Model transforms.
c
c Revision 1.1  1994/05/26  11:51:59  schaer
c Initial revision
c
c-----------------------------------------------------------------------


      subroutine t_drawfile(arg1,arg2,arg3)
c---------------------------------------------------------------------------
c     purpose:
c        this routine is used to put user-text from a file on plots.
c        if there is a map set, the user passes lat and lon, if
c        not, values between 0 and 1 (left to right, bottom to top).
c     arguments:
c        arg1    char  input   the filename of the file which contains the
C			       plot information.
c        arg2    char  input   which parameter from the file in arg1
c        arg3    char  input   not used.
c     history:
c        written by Christoph Frei & David Bresch 050194
c-----------------------------------------------------------------------
c     argument declarations.
      character *(80)  arg1,arg2,arg3
      
c     internal vars:
      integer   ibeg,iend,argno
      character*16 text(4)
      character*8  char8y,char8x
      character*80 ueber,chary,charx
c
c     function forward declarations:
      integer   strbeg,strend
c
c     convert arg2 into an integer, specifying the column number:
      read(arg2(1:strend(arg2)),'(i5)') argno

c     convert arg1 in a filename:
c     IVE passes only UPCASE strings, convert it:
      ibeg=strbeg(arg1)
      iend=strend(arg1)
      call locase(arg1(ibeg:iend),iend-ibeg+1)
      print *,'Drawn from file: ',arg1(strbeg(arg1):strend(arg1))
c
c     open the file with filename arg1:
      open(31,file=arg1,status='OLD')

c     loop over every line in file:
c
c     read the position (chary=latitude, charx=longitude) and text:
 888  read (31,99,end=999) char8y,char8x,text(1),text(2),text(3),text(4)
  99  format (2A8,4A16)
c     convert from char*8 to char*80:
      charx=char8x
      chary=char8y
c     draw the caracters in text on the plot at location x,y:
c     chary:latitude, charx:longitude.
      ueber=text(argno)
      call t_drawtext(charx,chary,ueber,0)
c     next:
      goto 888

c     end loop
 999  continue
      close(31)
      return
      end !of drawfile

      subroutine t_set_rangecheck(arg1,arg2,arg3)
c---------------------------------------------------------------------------
c     purpose:
c        this routine sets the rangecheck flag:
c        1:the parameter range when drawing user-defined on the plot
c        is checked, 0:not. 2: draw and symb coords are interpreted as
c        relative screen coords, ranging from 0 to 1 in plot region.
c     arguments:
c        arg1   char  input   containig rangecheck (int) as string
c                             rangecheck=1: the range is checked.
c        arg2   char  input   not used
c        arg3   char  input   not used
c     history:
c        written by david n. bresch 940120
c     arg2 and arg3 are free to be used!
c-----------------------------------------------------------------------
      include 'text.icl'
c subroutine arguments:
      character*(80)  arg1,arg2,arg3
c internal vars:
c
c forward definitions:
      integer         strend
c
c process arg1: convert it to int:
      if (arg1.ne.'#') then
         read(arg1(1:strend(arg1)),'(i3)') rangecheck
      endif
c  
      return
      end ! of set_rangecheck

      subroutine t_set_drawparam(arg1,arg2,arg3)
c---------------------------------------------------------------------------
c     purpose:
c        this routine sets the parameters for drawing etc. on the plot
c     arguments:
c        arg1   char  input   containig the size (real) as string
c        arg2   char  input   containig the color (int) as string
c        arg3   char  input   containig center, orient (real) as string
c                             separated by '/', like:
c                             -1.0/0.0
c                             (the code below can easily be extended to
c                             up to six parameters separated by ',')
c     history:
c        written by david n. bresch 940118
c     to simplify the adaption of these routine to specific user
c     needs, the locations where a new variable has to be inserted are
c     indicated (the new variable is called my_var), see also the 
c     file that contains the var-declarations: text.icl.
c-----------------------------------------------------------------------
      include 'text.icl'
c subroutine arguments:
      character*(80)  arg1,arg2,arg3
c internal vars:
      integer         maxarr
      parameter(maxarr=6) !change this if you need more than 6 arguments
                          !in arg3, separated by '/'.
      integer         ierr,p,l,begin
      real            rarg3(maxarr) !maximal maxarr arguments in last parameter
c
c forward definitions:
      integer         strbeg,strend
      real            str2nmb
c
c set default parameters for the arguments (add my_var here):
      data defcolor /1/
      data center /-1.0/
      data orient /0.0/
      data tsize /0.75/ !the same size as the minor text on plots
      data rangecheck /1/
c
c fill rarg3(*) with dummy values to be able to check them later
c if they are filled with good arguments:
      do p=1,maxarr
         rarg3(p)=-99.9
      enddo
c
c process arg1: convert it to real:
      if (arg1.ne.'#') then
         tsize=str2nmb(arg1(1:strend(arg1)),ierr)
      endif
c
c process arg2: convert it to int:
      if (arg2.ne.'#') then
         read(arg2(1:strend(arg2)),'(i3)') defcolor
      endif
c
c process arg3: search for '/' which separates numbers and fill each number
c               into rarg2(i)
      if (arg3.eq.'#') then
         return !return immediately
      endif
         p=1
         begin=1
         do l = strbeg(arg3), strend(arg3)
            if (arg3(l:l).eq.'/') then
               rarg3(p)=str2nmb(arg3(begin:l-1),ierr)
               p=p+1
               begin=l+1
               if (p.gt.maxarr) then 
                  go to 99
               endif
            endif
         enddo
         rarg3(p)=str2nmb(arg3(begin:strend(arg3)),ierr)
 99      continue
c set the values of all parameters passed in arg3 (add here, if you are
c adding a new parameter to arg3):
         if (rarg3(1).ne.-99.9) then
            center=rarg3(1)
         endif
         if (rarg3(2).ne.-99.9) then
            orient=rarg3(2)
         endif
c change the following lines to add your new var:
c         if (rarg3(3).ne.-99.9) then
c           my_var=rarg3(3)
c         endif
c
      return
      end ! of t_set_drawparam


      subroutine t_drawtext(arg1,arg2,text,norm)
c---------------------------------------------------------------------------
c     purpose:
c        this routine is used to put user-text on plots.
c        the coord-system is the one shown when drawing axes. If there
c        is a map, global geographical coords are used.
c        other purposes: see the parameter norm.
c     arguments:
c        arg1    char  input   the latitude of the text position. (as char)
c        arg2    char  input   the longitude of the text position. (as char)
c        text    char  input   the text itself.
c        norm    int   input   =0(default):use data coords resp. map
c                              =1:use normalized coords, interprete
c                              arg and arg2 as coords ranging from 0 to 1
c                              (in this case, no rangecheck is done)
c     history:
c        written by david n. bresch 931102
c        this routine plots directly!
c-----------------------------------------------------------------------
c
      include 'text.icl'
c
c     argument declarations.
      real             userx,usery
      character *(80)  text,arg1,arg2
      integer           norm
c
      logical          mapflg
c
c     local variable declarations.
c
      integer          ierror, line_index, linlog,  
     &                 qual, text_index
      real             vpl, vpr, vpb, vpt, wdl, wdr, wdb, wdt
      real             szsf,xpos, ypos
c
      real             umin,umax,vmin,vmax,u,v
      integer          found
c
c     external function declarations.
      integer         strbeg, strend
      real            str2nmb
c
c     internal vars:
      integer         ierr
      logical         in_region,lerror
c
      real            plwmin(4),plwmax(4),x1,y1
      real            xmin,xmax,ymin,ymax,dummy
      integer         lock(4),firstaxis,error
c
c     convert characters to real
c
      userx=str2nmb(arg1(1:strend(arg1)),ierr)
      usery=str2nmb(arg2(1:strend(arg2)),ierr)
c
c     if rangecheck is =2, set norm to 1 to use normalized coords:
      if (rangecheck.eq.2) then
         norm=1
      endif
c
c     set text and line colors to the background color.
c
      call gqtxci (ierror, text_index)
      call gqplci (ierror, line_index)
      call gstxci (defcolor) !was 1
      call gsplci (defcolor) !was 1
      call gstxfp (1, 1)
c
c     get current set parameters.
c
      call gsclip (0)
      call getset (vpl, vpr, vpb, vpt, wdl, wdr, wdb, wdt, linlog)
      szsf = tsize * (vpr - vpl) * 0.019
c
c     call set so we can use fractional coordinates.
c
      call set (vpl, vpr, vpb, vpt, 0.0, 1.0, 0.0, 1.0, 1)
c
      call pcgeti ('QU - quality flag', qual)
c     
c     use medium quality characters.
c
      call pcseti ('QU - quality flag', 1)
c
c     get the mapflag:
      call getlvar('mapflg',mapflg,ierror)
c     if z or t are free, don't use map-conversion:
      call getiarr('lock',lock,4,error)
      if ((lock(3).eq.0).or.(lock(4).eq.0)) then
         mapflg=.false.
      endif
c     test, if a map has been underlied:mapflg or normal coords are used:
      if (mapflg.or.(norm.eq.1)) then !begins MAP clause
c         print *,'mapflg is TRUE or norm =1'
         if (norm.eq.0) then
c        convert coordinates:
c        userx and usery contain in this case longitude and latitude resp.
c        the foll. 2 lines are used if rotated coord. are given by the user:
c           call ph2ll(usery,userx,lat,lon)
c           call maptrn(lat,lon,u,v)
            call maptrn(usery,userx,u,v)
c        get the umin,umax ... values:
            call map_clip(umin,umax,vmin,vmax,found)
c
c     test, if the points are in the plotted region:
            if (rangecheck.eq.1) then
               call get_region(1,lerror)
               if (in_region(1,u,v).eq.0)  then
c                 print *,'point not in region'
                  norm=0
                  return
               endif
            endif
c
c        scale the position to screen-coordinates:
            xpos=(u-umin)/(umax-umin)
            ypos=(v-vmin)/(vmax-vmin)
         else
            xpos=userx
            ypos=usery
         endif !of map convert clause
c
c     call set so we can use fractional coordinates.
c
         call set (vpl, vpr, vpb, vpt, 0.0, 1.0, 0.0, 1.0, 1)
c
         call plchlq (xpos, ypos, text(strbeg(text):strend(text)), szsf,
     &            orient,center)
c
      else !end of MAP clause, begin of NO MAP clause
c
c     copy userx to x1 because of code in drawline
c
      x1=userx
      y1=usery
c
c     get the boundary values (for range checking and window settting):
c
         call getrarr('plwmin',plwmin,4,error)
         call getrarr('plwmax',plwmax,4,error)
c         
c     decide, which axes are active (only two values in lock are ne 0):
         firstaxis=1
         if (lock(1).eq.0) then
c     x-axis active:
            xmin=plwmin(1)
            xmax=plwmax(1)
            firstaxis=0
         endif
         if (lock(2).eq.0) then
c     y-axis active:
            if (firstaxis.eq.1) then
               xmin=plwmin(2)
               xmax=plwmax(2)
               firstaxis=0
            else
               ymin=plwmin(2)
               ymax=plwmax(2)
            endif
         endif
         if (lock(3).eq.0) then
c     z-axis active:
            if (firstaxis.eq.1) then
               xmin=plwmin(3)
               xmax=plwmax(3)
               firstaxis=0
            else
               ymin=plwmin(3)
               ymax=plwmax(3)
            endif
         endif
         if (lock(4).eq.0) then
c     t-axis active:
            if (firstaxis.eq.1) then
               xmin=plwmin(4)
               xmax=plwmax(4)
               firstaxis=0
            else
               ymin=plwmin(4)
               ymax=plwmax(4)
            endif
         endif !of axis check
c     call set to adjust coordinate frame:
      call set (vpl, vpr, vpb, vpt, xmin, xmax, ymin, ymax, linlog)
c
c     if the coordinates are in reverse order, change:
         if (xmax.lt.xmin) then
            dummy=xmin
            xmin =xmax
            xmax =dummy
         endif
         if (ymax.lt.ymin) then
            dummy=ymin
            ymin =ymax
            ymax =dummy
         endif
c
c     test user coordinates: if one is out of range, don't plot anything:
c     new:test only, if the rangecheck-flag is set:
         if (rangecheck.eq.1) then
            if ((x1.lt.xmin).or.(x1.gt.xmax)) then
               norm=0
               return
            endif
            if ((y1.lt.ymin).or.(y1.gt.ymax)) then
               norm=0
               return
            endif
         endif !of rangecheck
c
c     call set so we can use fractional coordinates.
c
c         call set (vpl, vpr, vpb, vpt, 0.0, 1.0, 0.0, 1.0, 1)
c     convert x1 and y1 into positions fractional:
c
c         xpos=x1/(xmax-xmin)
c         ypos=y1/(ymax-ymin)
         xpos=x1
         ypos=y1
c
c     place label just above top plot border.(folowing 2 lines)
c
         call plchlq (xpos, ypos, text(strbeg(text):strend(text)), szsf,
     &            orient,center)
c
      endif !of NO MAP clause
c
c     reset quality flag, and call set with original values.
c
      call pcseti ('QU - quality flag', qual) 
      call set (vpl, vpr, vpb, vpt, wdl, wdr, wdb, wdt, linlog)
c
c     reset text and line color.
c
      call gstxci (text_index)
      call gsplci (line_index)
c
c     release all output.
c
      call plotit (0, 0, 0)
c
c     set norm to zero:
      norm=0
      return
      end ! of t_drawtext
c
      subroutine t_drawsymb(arg1,arg2,arg3)
c---------------------------------------------------------------------------
c     purpose:
c        this routine is used to draw symbols on plots.
c        normally, the slider coordinates (the same that are appearing
c        when drawing the axes) are used.
c        if a map is set, the global geographic coordinates are used.
c     arguments:
c        arg1   char  input   horizontal axis or longitude
c        arg2   char  input   vertical axis or latitude
c        arg3   char  input   the symbol type:
c                             0(default):cross,1:triangle
c     history:
c        written by david n. bresch 931124
c-----------------------------------------------------------------------
c
      include 'text.icl'
c
c     argument declarations.
      character *(80)  arg1,arg2,arg3
c
c     local variable declarations.
c
      real            lat,lon,locxsize,locysize
      integer         type,lock(4),firstaxis,error
      real            plwmin(4),plwmax(4),xmin,xmax,ymin,ymax

c     internal vars:
      real            x1,x2,y1,y2,ierr

c     external function declarations.
c
      integer         strend,nint
      real            str2nmb
c
c     convert characters to real
c
      lat=str2nmb(arg2(1:strend(arg2)),ierr)
      lon=str2nmb(arg1(1:strend(arg1)),ierr)
      type=nint(str2nmb(arg3(1:strend(arg3)),ierr))
c
c     set the size:
c     according to the active coordinate-range and map (if any):
      call getrarr('plwmin',plwmin,4,error)
      call getrarr('plwmax',plwmax,4,error)
      call getiarr('lock',lock,4,error)
c
c     decide, which axes are active (only two values in lock are ne 0):
      firstaxis=1
      if (lock(1).eq.0) then
c     x-axis active:
         xmin=plwmin(1)
         xmax=plwmax(1)
         firstaxis=0
      endif
      if (lock(2).eq.0) then
c     x-axis active:
         if (firstaxis.eq.1) then
            xmin=plwmin(2)
            xmax=plwmax(2)
            firstaxis=0
         else
            ymin=plwmin(2)
            ymax=plwmax(2)
         endif
      endif
      if (lock(3).eq.0) then
c     z-axis active:
         if (firstaxis.eq.1) then
            xmin=plwmin(3)
            xmax=plwmax(3)
            firstaxis=0
         else
            ymin=plwmin(3)
            ymax=plwmax(3)
         endif
      endif
      if (lock(4).eq.0) then
c     t-axis active:
         if (firstaxis.eq.1) then
            xmin=plwmin(4)
            xmax=plwmax(4)
            firstaxis=0
         else
            ymin=plwmin(4)
            ymax=plwmax(4)
         endif
      endif
c
c     set the scaled size:
c
      locxsize=tsize*(xmax-xmin)/100.
      locysize=tsize*(ymax-ymin)/100.
c
      if (type.eq.0) then
c     print a cross:
         x1=lon-locxsize
         y1=lat-locysize
         x2=lon+locxsize
         y2=lat+locysize
c     call the line-drawing routine:
c     hint:(lon,lat),(lon,lat),color
         call t_drawline(x1,y1,x2,y2,defcolor)
         x1=lon-locxsize
         y1=lat+locysize
         x2=lon+locxsize
         y2=lat-locysize
         call t_drawline(x1,y1,x2,y2,defcolor)
      elseif (type.eq.1) then
c     print a triangle:
         x1=lon-locxsize
         y1=lat-locysize
         x2=lon
         y2=lat+locysize
         call t_drawline(x1,y1,x2,y2,defcolor)
         x1=lon
         y1=lat+locysize
         x2=lon+locxsize
         y2=lat-locysize
         call t_drawline(x1,y1,x2,y2,defcolor)
         x1=lon+locxsize
         y1=lat-locysize
         x2=lon-locxsize
         y2=lat-locysize
         call t_drawline(x1,y1,x2,y2,defcolor)
cc    start here to define your own symbol:
cc    the center of your symbol should be x=lon,y=lat
cc    the unitlenght in x-direction is locxsize, in y-dir locysize
c     elseif (type.eq.2) then
c     print a user defined symbol
c     first line (x1,y1)-(x2,y2)
c        x1=
c        y1=
c        x2=
c        y2=
c        call t_drawline(x1,y1,x2,y2,defcolor)
c     copy the last five lines for the next line to be drawn... 
c     end of definition of new symbols.     
      endif
      return
      end ! of t_drawsymb

      subroutine t_drawline(x1,y1,x2,y2,linecol)
c---------------------------------------------------------------------------
c     purpose:
c        this routine is used to draw a line on a plot.
c        from (x1,y1) to (x2,y2)
c        (lo-level routine, don't change it see the calling routines).
c        if there is no map, the routine uses coordinates
c        ranging from 0 to 1 in both directions.
c     arguments:
c        x1     real  input   horizontal axis or longitude start
c        y1     real  input   vertical axis or latitude start
c        x2                   end
c        y2
c        linecol integer inp  the color of the line:
c                             1:black,2:white,3:red,4:orange,5:yellow
c                             6:green,7:blue,8:violet,9:black,10:grey
c     history:
c        written by david n. bresch 940120
c        this routine plots directly!
c-----------------------------------------------------------------------
c
      include 'text.icl'
c
c     argument declarations.
      real             x1,y1,x2,y2
      integer          linecol
c
      logical          mapflg,in_region,error
c
c     local variable declarations.
c
      integer          ierror, line_index, linlog,  
     &                 qual, text_index
      real             vpl, vpr, vpb, vpt, wdl, wdr, wdb, wdt
      integer          lineclr1,linepat1,found
      real             xpts(2),ypts(2)
      real             xmin, xmax, ymin, ymax,dummy
      real             plwmin(4),plwmax(4)
      integer          lock(4),firstaxis
      real             umin,umax,vmin,vmax,u1,u2,v1,v2
c
c     set text and line colors to the background color.
c
      call gqtxci (ierror, text_index)
      call gqplci (ierror, line_index)
      call gstxci (1)
      call gsplci (1)
      call gstxfp (1, 1)
c
c     get current set parameters.
c
      call gsclip (0)
      call getset (vpl, vpr, vpb, vpt, wdl, wdr, wdb, wdt, linlog)
c
      call pcgeti ('QU - quality flag', qual)
c     
c     use medium quality characters.
c
      call pcseti ('QU - quality flag', 1)
c
c     set the linecolor:
      call getivar ('hicolor', lineclr1, error)
c     set user-linecolor, if defined
      if (linecol.ne.0) then
         lineclr1=linecol
      endif
      call gsplci (lineclr1)
      call getivar ('hclpat', linepat1, error)
      call dashdb (linepat1)
c
c     get the mapflag:
      call getlvar('mapflg',mapflg,ierror)
c     if z or t are free, don't use map-conversion:
      call getiarr('lock',lock,4,error)
      if ((lock(3).eq.0).or.(lock(4).eq.0)) then
         mapflg=.false.
      endif
c     test, if a map has been underlied:mapflg
      if ((mapflg).or.(rangecheck.eq.2)) then !begin of MAP clause
c     the follwing lines up to the label NO MAP are used to plot with map:
         if (rangecheck.ne.2) then !=2:use normalized coords 0-1
c        get the umin,umax ... values:
            call map_clip(umin,umax,vmin,vmax,found)
c        convert coordinates:
c        userx and usery contain in this case longitude and latitude resp.
c        the foll. 2 lines are used if rotated coord. are given by the user:
            call maptrn(y1,x1,u1,v1)
            call maptrn(y2,x2,u2,v2)
            call get_region(1,error)
c
c        test, if the points are in the plotted region:
            if (rangecheck.eq.1) then
               if    ( (in_region(1,u1,v1).eq.0) .or. 
     &            (in_region(1,u2,v2).eq.0) ) then
                  return
               endif
            endif
c
c     scale the position to screen-coordinates:
            xpts(1)=(u1-umin)/(umax-umin)
            ypts(1)=(v1-vmin)/(vmax-vmin)
            xpts(2)=(u2-umin)/(umax-umin)
            ypts(2)=(v2-vmin)/(vmax-vmin)
         else
c     use normalized coords:
            xpts(1)=x1
            xpts(2)=x2
            ypts(1)=y1
            ypts(2)=y2
         endif !of rangecheck.ne.2 clause.
c
c     draw the line:
c     call set so we can use fractional coordinates.
c
         call set (vpl, vpr, vpb, vpt, 0.0, 1.0, 0.0, 1.0, 1)
c
         call curved (xpts, ypts, 2)
c
      else !end of MAP clause, begin of NO MAP clause
c     
c     get the boundary values (for range checking and window settting):
c
         call getrarr('plwmin',plwmin,4,error)
         call getrarr('plwmax',plwmax,4,error)
         call getiarr('lock',lock,4,error)
c         
c     decide, which axes are active (only two values in lock are ne 0):
         firstaxis=1
         if (lock(1).eq.0) then
c     x-axis active:
            xmin=plwmin(1)
            xmax=plwmax(1)
            firstaxis=0
         endif
         if (lock(2).eq.0) then
c     y-axis active:
            if (firstaxis.eq.1) then
               xmin=plwmin(2)
               xmax=plwmax(2)
               firstaxis=0
            else
               ymin=plwmin(2)
               ymax=plwmax(2)
            endif
         endif
         if (lock(3).eq.0) then
c     z-axis active:
            if (firstaxis.eq.1) then
               xmin=plwmin(3)
               xmax=plwmax(3)
               firstaxis=0
            else
               ymin=plwmin(3)
               ymax=plwmax(3)
            endif
         endif
         if (lock(4).eq.0) then
c     t-axis active:
            if (firstaxis.eq.1) then
               xmin=plwmin(4)
               xmax=plwmax(4)
               firstaxis=0
            else
               ymin=plwmin(4)
               ymax=plwmax(4)
            endif
         endif !of axis check

c     call set to adjust coordinate frame:
      call set (vpl, vpr, vpb, vpt, xmin, xmax, ymin, ymax, linlog)
c
c     if the coordinates are in reverse order, change:
         if (xmax.lt.xmin) then
            dummy=xmin
            xmin =xmax
            xmax =dummy
         endif
         if (ymax.lt.ymin) then
            dummy=ymin
            ymin =ymax
            ymax =dummy
         endif
c
c     test user coordinates: if one is out of range, don't plot anything:
c     new:test only, if the rangecheck-flag is set:
         if (rangecheck.eq.1) then
            if ((x1.lt.xmin).or.(x1.gt.xmax)) then
               return
            endif
            if ((x2.lt.xmin).or.(x2.gt.xmax)) then
               return
            endif
            if ((y1.lt.ymin).or.(y1.gt.ymax)) then
               return
            endif
            if ((y2.lt.ymin).or.(y2.gt.ymax)) then
               return
            endif
         endif !of rangecheck
c
c     plot the line:
         call line ( x1, y1, x2, y2 )
c
      endif !of NO MAP clause
c
c     reset quality flag, and call set with original values.
c
      call pcseti ('QU - quality flag', qual) 
      call set (vpl, vpr, vpb, vpt, wdl, wdr, wdb, wdt, linlog)
c
c     reset text and line color.
c
      call gstxci (text_index)
      call gsplci (line_index)
c
c     release all output.
c
      call plotit (0, 0, 0)
c
      return
      end ! of t_drawline

      subroutine t_proj_wind(a,b,c,d,m,nx,ny,nz,nt,misdat)
c----------------------------------------------------------------------------
c
c     computes the projection of the horizontal wind field to the
c     inclined cross-section plane.
c     a contains u
c     b contains v
c     c contains the inclination of the plane (counterclockwise, radian)
c     d integer: 0, if projection on the plane, 1 if projection to the
c                perpendicular plane.
c
c     written by david bresch 11/93
c---------------------------------------------------------------------------
c
c     declaration of parameters
      integer    nx,ny,nz,nt,d
      real       a(nx,ny,nz,nt),b(nx,ny,nz,nt),m(nx,ny,nz,nt),
     &           misdat,c

c     declaration of variables
      integer    i,j,k,l
      real       ang
      real       pt1(4),pt2(4)
      integer    err

c     definition of function types:
      real       atan,cos

      print *,'computing vector projection'
      print *,'  for an array dimensioned',nx,ny,nz,nt

c     test, if the user provided an angle (c not zero):
      if (c.eq.0.0) then
c        get the cross-sections-angle:
         call getrarr('point_1',pt1,4,err)
         call getrarr('point_2',pt2,4,err)
c         print *,'pt1: ',pt1(1),pt1(2),pt1(3),pt1(4)
c         print *,'pt2: ',pt2(1),pt2(2),pt2(3),pt2(4)
         if ((pt2(1)-pt1(1)).ne.0.0) then
c        atan is only from -pi to pi, test x:
            if ( (pt2(1)-pt1(1)) .lt.0.0) then
               c=atan((pt2(2)-pt1(2))/(pt2(1)-pt1(1)))+3.1415926
            else
               c=atan((pt2(2)-pt1(2))/(pt2(1)-pt1(1)))
            endif
         else
            c=3.1415926/2.0
         endif
         print *,'  for an angle (computed) of',c
      else
         print *,'  for an angle (user defined) of',c
      endif
c
c     start computation
      do l=1,nt
        do k=1,nz
          do j=1,ny
            do i=1,nx
               if ((a(i,j,k,l).ne.0.0).and.(a(i,j,k,l).ne.misdat).and.
     &            (b(i,j,k,l).ne.misdat)) then
c                 atan is only from -pi to pi, test x:
                  if (a(i,j,k,l).lt.0.0) then
                     ang=atan(b(i,j,k,l)/a(i,j,k,l))+3.1415926
                  else
                     ang=atan(b(i,j,k,l)/a(i,j,k,l))
                  endif
c                 check, if projection on the plane or perpendicular:
                  if (d.eq.0) then
                     m(i,j,k,l)=(a(i,j,k,l)*a(i,j,k,l)+
     &                          b(i,j,k,l)*b(i,j,k,l))
     &                          **(0.5)*cos(ang-c)
                  else
                     m(i,j,k,l)=-(a(i,j,k,l)*a(i,j,k,l)+
     &                          b(i,j,k,l)*b(i,j,k,l))
     &                          **(0.5)*sin(ang-c)
                  endif
              else
                  m(i,j,k,l)=misdat
              endif
            enddo
          enddo
        enddo
      enddo
      end !of t_proj_wind
