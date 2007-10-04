      integer function privat (varnam, ndims, dims, stag, datmin, 
     &                           datmax, misdat, 
     &                           data_units, data_display_units) 
c-----------------------------------------------------------------------
c     $id: privat.f,v 1.2 1993/09/02 07:25:49 schaer exp $
c     purpose:
c        this is a user-written function that is used to calculate 
c        user-derived fields. 
c     arguments:
c        varnam  char  input   the name of the variable to be derived.
c        ndims   int   output  the number of dimensions of the variable.
c                              assumed to be 4, which is okay for any
c                              variable with not more than 4 dimensions.
c        dims    int   output  the number of data points along each 
c                              dimension. 
c                              this array must be supplied by the user.
c                              for example, if dims(1) = nx, 
c                              dims(2) = ny, dims(3), this is a 3d data
c                              set that is nx x ny x nz.
c        stag    real  output  the grid staggering along each dimension.
c                              this array must be supplied by the user.
c        datmin  real  output  the location of the origin of the data
c                              in physical space.
c                              this array must be supplied by the user.
c        datmax  real  output  the extent of the data in physical space
c                              along each dimension.
c                              this array must be supplied by the user.
c        misdat  real  output  the missing data value. any point whose
c                              value is misdat will be ignored by the
c                              plotting routines.
c                              this value must be supplied by the user.
c     note : on output the value of the function is equal to either 0 or
c            some integer that points to memory where the derived field
c            is stored. if the function value is 0, it indicates to
c            the calling routine that there was some problem in 
c            calculating the field.
c     history:
c	$log: privat.f,v $
c revision 1.2  1993/09/02  07:25:49  schaer
c kleinere aenderungen damit es schoener aussieht.
c
c revision 1.1  1993/08/27  07:01:24  schaer
c initial revision
c
c---------------------------------------------------------------------- 
c
c     do compiler dependent definition of recursive variables
c#ifdef ultrix
c
c#else
      implicit automatic (a-z)
      static diagfld
c#endif

c     include the common-block with the level information
      include 'constants.icl'

c     argument declarations. maxdim is the maximum data dimension allowed,
c     lenvar is the maximum length of the variable-name.
      integer   maxdim, lenvar
      parameter (maxdim=4, lenvar=80)

      character*(80)   data_units, data_display_units
      character*(lenvar)  varnam
      integer          ndims, dims(maxdim)
      real             stag(maxdim), misdat,
     &                 datmin(maxdim), datmax(maxdim)

c     local variable declarations
      character*(lenvar) fld1,fld2,arg1,arg2,arg3,dvar
      real             misdatu,misdatv,misdvel,misdatd,str2nmb
      integer          ptrt,ptrvel,ptrd,ptrrr,size,ptrdps,ierr
      logical          errflg,localu,localv

c     local variable declarations.
      character*(lenvar) fld
      integer          i, ibeg, iend, ibeg1, iend1, ptr, dims1(maxdim)
      logical          err, local
      real             delta1,delta2

c     variables for wind:
      integer        ptr1,ptr2,ptrd,size,d
      real           c

c     external function declarations.
      integer        getmem, addvar, getvar, strbeg, strend, gcdfvar
      logical        isfunc

c     initially, set pointer 'privat' to unused, indicating a
c     failure in computing the requested field.
      privat=0

c     define field-name to be computed (removes leading and trailing blanks)
      ibeg = strbeg (varnam)
      iend = strend (varnam)
      fld=varnam(ibeg:iend)
      iend = strend (fld)

c     start computation of fields

      if (fld.eq.'PVALL') then

c          addvar und getvar sind in ~schaer/uwgap7/uars/getvar.f

           privat=getvar('-(vort+f)*9.80665*d[theta:p]',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if (privat.eq.0) return
c          assign units (pvu's)
           data_units='k m**2 / (kg s)'
           data_display_units='1.e-6  k m**2 / (kg s)'

      endif

      if (isfunc(fld,'WIND',arg1,arg2,arg3)) then
c     calculates the projection of the wind vector to the current
c     cross-section plane (if the passed argument is positive,
c     it is just used to number the field, if negative, it is
c     converted to positive and interpreted as the desired angle)
c     the second argument is not requested, but when set to 1, the
c     vectors are interpolated to a plane perpendicular to the
c     cross-section.
c     the formula is (resp. -sin for second argument =1):
c     wind=((u*u+v*v)^(.5)*cos(atan(v/u)-arg1))
      if (arg1.eq.'#') then
        print *,'enter the cross-section angle in ',fld(1:iend)
        return
      endif
c     addvar und getvar sind in ~schaer/uwgap7/uars/getvar.f
c
c     arg1 contains the angle of the plane (as string, in degree),
c     convert it to real:
           c=str2nmb(arg1(1:strend(arg1)),ierr)
c           read(arg1(1:strend(arg1)),'(f5.1)') c
           c=c*0.0174533   ! convert degree to rad
c     test, if c is positive (in this case, the angle is 
c     calculated in proj_wind) or negative (then it is the angle):
           if (c.ge.0.0) then
               c=0.0
           else
c     c is the angle, convert it to positive:
               c=-c
           endif
c     convert the second argument:
      if (arg2.ne.'#') then
           d=str2nmb(arg2(1:strend(arg2)),ierr)
c           read(arg2(1:strend(arg2)),'(i3)') d
      else
           d=0
      endif
c
c     get the u-component:
           ptr1=getvar('u',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if (ptr1.eq.0) return
c
c     get the v-component:
           ptr2=getvar('v',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if (ptr2.eq.0) return
c
c       define the dimension of the field and get the memory for it
        size=1
        do i = 1, ndims
          size = size * dims(i)
        enddo
        ptrd  = getmem (size)
        if (ptrd.eq.0) return

c       call the subroutine to compute the projected wind:
           call proj_wind(%val(ptr1),%val(ptr2),c,d,%val(ptrd),
     &             dims(1),dims(2),dims(3),dims(4),misdat)
c
           privat=ptrd
           if (privat.eq.0) return
c          assign units (m/s)
           data_units='m/s'
           data_display_units='m/s'
      endif ! wind

      if (isfunc(fld,'TRAJ_DIR',arg1,arg2,arg3)) then
c     show trajectories on the plot. read from file.
c     see t_dotraj about a description of the arguments.
         if (arg1.eq.'#') then
            print *,'first argument: directory name in: ',fld(1:iend)
            return
         endif
         call t_settraj_dir(arg1)
         privat=0
      endif ! traj_dir

      if (isfunc(fld,'TRAJ',arg1,arg2,arg3)) then
c     show trajectories on the plot. read from file.
c     see t_dotraj about a description of the arguments.
         if (arg1.eq.'#') then
            print *,'first argument: filename needed by ',fld(1:iend)
            return
         endif
         call t_dotraj(arg1,arg2,arg3)
         privat=0
      endif ! traj

      if (isfunc(fld,'DRAW',arg1,arg2,arg3)) then
      if ((arg1.eq.'#').or.(arg2.eq.'#').or.(arg3.eq.'#')) then
        print *,'illegal number of arguments in ',fld(1:iend)
        return
      endif
c     draw text on plot:
           call t_drawtext(arg1,arg2,arg3)
           privat=0
      endif ! draw

      if (isfunc(fld,'SYMB',arg1,arg2,arg3)) then
c     see the purpose of t_drawsymb about the meaning of argx.
         if ((arg1.eq.'#').or.(arg2.eq.'#').or.(arg3.eq.'#')) then
            print *,'illegal number of arguments in ',fld(1:iend)
            return
         endif
c     draw symbol on plot:
            call t_drawsymb(arg1,arg2,arg3)
            privat=0
      endif ! symb
      end
c
c
c--------------------------------------------------------------------------
c The following section contains all the subroutines called above.
c--------------------------------------------------------------------------
c
c
      subroutine t_drawtext(arg1,arg2,text)
c---------------------------------------------------------------------------
c     purpose:
c        this routine is used to put user-text on plots.
c        if there is a map set, the user passes lat and lon, if
c        not, values between 0 and 1 (left to right, bottom to top).
c     arguments:
c        arg1    char  input   the latitude of the text position. (as char)
c        arg2    char  input   the longitude of the text position. (as char)
c        text    char  input   the text itself.
c     history:
c        written by david n. bresch 931102
c-----------------------------------------------------------------------
c
c     argument declarations.
      real             userx,usery
      character *(80)  text,arg1,arg2
c
      logical          mapflg
c
c     local variable declarations.
c
      integer          ierror, line_index, linlog,  
     &                 qual, text_index
      real             vpl, vpr, vpb, vpt, wdl, wdr, wdb, wdt
      real             center, orient, size, szsf, 
     &                 xpos, ypos
      data  size / 0.019 /
c
      real             umin,umax,vmin,vmax,u,v
      integer          found
      character *(80)  command
c
c     external function declarations.
      integer         strbeg, strend
      real            str2nmb
c
c     internal vars:
      integer         ierr
c
c     convert characters to real
c
      userx=str2nmb(arg1(1:strend(arg1)),ierr)
      usery=str2nmb(arg2(1:strend(arg2)),ierr)
c
c     create the string which contains the command:
      command='field=draw['//arg1(strbeg(arg1):strend(arg1))//':'
     &       //arg2(strbeg(arg2):strend(arg2))//':'
     &       //text(strbeg(text):strend(text))//']'
c     &       //'; overlay_plot'
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
      szsf = size * (vpr - vpl)
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
c     test, if a map has been underlied:mapflg
      mapflg=1
      if (mapflg) then
c        convert coordinates:
c        userx and usery contain in this case longitude and latitude resp.
c        the foll. 2 lines are used if rotated coord. are given by the user:
c           call ph2ll(usery,userx,lat,lon)
c           call maptrn(lat,lon,u,v)
            call maptrn(usery,userx,u,v)
c        get the umin,umax ... values:
            call map_clip(umin,umax,vmin,vmax,found)
c        scale the position to screen-coordinates:
            xpos=(u-umin)/(umax-umin)
            ypos=(v-vmin)/(vmax-vmin)
         else
            xpos=userx
            ypos=usery
      endif
c
c     call set so we can use fractional coordinates.
c
      call set (vpl, vpr, vpb, vpt, 0.0, 1.0, 0.0, 1.0, 1)
c
c     place label just above top plot border.(folowing 2 lines)
c      xpos = cfux (vpl)
c      ypos = cfuy (vpt + 0.04)
      center = -1.0
      orient = 0.0
c      ibeg=strbeg(text)
c      iend=strend(text)
c      print *,text(ibeg:iend),' | the new text'
c      print *,'xpos: ',xpos,' ypos: ',ypos
c      print *,'orient: ',orient,' center: ',center
c      print *,'szsf: ',szsf
      call plchlq (xpos, ypos, text(strbeg(text):strend(text)), szsf,
     &            orient,center)
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
      call save_plot_command(0,1)

c     release all output.
c
      call plotit (0, 0, 0)
c
      call add_to_list(command(strbeg(command):strend(command)))
c
      return
      end ! of t_drawtext

c content of the draw-command list:
c  plot      command
c number
c ======     =======
c     1      transform=/home/david/ive/em_trafos/trans_em.o; file=/em2/w03/w03_chst_12; zloc=400; plat=47.5; plon=10.0; gridspacing=20; drawmap=on; maproj=st;c maplimits=cl; window; field=u; vector_comp=u,v; vector_int=2,2; input=0,0; ticme=12.000000; zloc=400.000000; time=12.000000; zloc=400.000000; time=12.000000c; new_plot; field=draw[47.:7.:ttt]; overlay_plot



      subroutine t_drawsymb(arg1,arg2,arg3)
c---------------------------------------------------------------------------
c     purpose:
c        this routine is used to draw symbols on plots.
c        if there is no map, the routine attempts to use coordinates
c        ranging from 0 to 1 in both directions.
c     arguments:
c        arg1   char  input   latitude
c        arg2   char  input   longitude
c        arg3   char  input   the color
c     history:
c        written by david n. bresch 931124
c-----------------------------------------------------------------------
c
c     argument declarations.
      character *(80)  arg1,arg2,arg3
c
c     local variable declarations.
c
      logical         mapflg
      real            lat,lon,size
      integer         color

c     internal vars:
      real            x1,x2,y1,y2,ierr

c     external function declarations.
c
      integer         strend,nint
      real            str2nmb
c
c     convert characters to real
c
      lat=str2nmb(arg1(1:strend(arg1)),ierr)
      lon=str2nmb(arg2(1:strend(arg2)),ierr)
      color=nint(str2nmb(arg3(1:strend(arg3)),ierr))
c
c     set the size:
      mapflg=1
      if (mapflg) then
         size=0.5
      else
         size=0.1
      endif
c
c     print a cross:
         x1=lon-size
         y1=lat-size
         x2=lon+size
         y2=lat+size
c     call the line-drawing routine:
c     hint:(lon,lat),(lon,lat),color
         call t_drawline(x1,y1,x2,y2,color)
         x1=lon-size
         y1=lat+size
         x2=lon+size
         y2=lat-size
c     hint:(lon,lat),(lon,lat),color
         call t_drawline(x1,y1,x2,y2,color)
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
c        x1     real  input   longitude start (or x start)
c        y1     real  input   latitude start (or ystart)
c        x2
c        y2
c        linecol integer inp  the color of the line:
c                             1:black,2:white,3:red,4:orange,5:yellow
c                             6:green,7:blue,8:violet,9:black,10:grey
c     history:
c        written by david n. bresch 931124
c-----------------------------------------------------------------------
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
      logical          error
      real             xpts(2),ypts(2)
c
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
c     test, if a map has been underlied:mapflg
      mapflg=1
      if (mapflg) then
c        get the umin,umax ... values:
            call map_clip(umin,umax,vmin,vmax,found)
c        convert coordinates:
c        userx and usery contain in this case longitude and latitude resp.
c        the foll. 2 lines are used if rotated coord. are given by the user:
c           call ph2ll(usery,userx,lat,lon)
c           call maptrn(y1,x1,u,v)
            call maptrn(y1,x1,u1,v1)
            call maptrn(y2,x2,u2,v2)
            call get_region(1,error)
c     test, if the points are in the plotted region:
            if ( (in_region(1,u1,v1).eq.0) .or. 
     &           (in_region(1,u2,v2).eq.0) ) then
c              print *,'point not in region'
               return
            endif
c        scale the position to screen-coordinates:
            xpts(1)=(u1-umin)/(umax-umin)
            ypts(1)=(v1-vmin)/(vmax-vmin)
            xpts(2)=(u2-umin)/(umax-umin)
            ypts(2)=(v2-vmin)/(vmax-vmin)
         else
            xpts(1)=x1
            xpts(2)=x2
            ypts(1)=y1
            ypts(2)=y2
      endif
c
c     call set so we can use fractional coordinates.
c
      call set (vpl, vpr, vpb, vpt, 0.0, 1.0, 0.0, 1.0, 1)
c
c     reset quality flag, and call set with original values.
c
      call getivar ('hicolor', lineclr1, error)
c     set user-linecolor, if defined
      if (linecol.ne.0) then
         lineclr1=linecol
      endif
      call gsplci (lineclr1)
      call getivar ('hclpat', linepat1, error)
      call dashdb (linepat1)
c     draw the line:
      call curved (xpts, ypts, 2)
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
c
c
      subroutine proj_wind(a,b,c,d,m,nx,ny,nz,nt,misdat)
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
      real       a(nx,ny,nz,nt),b(nx,ny,nz,nt),m(nx,ny,nz,nt),
     &           misdat,c
      integer    nx,ny,nz,nt,d

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
      end

      subroutine t_settraj_dir(arg1)
c---------------------------------------------------------------------------
c     purpose:
c        this routine sets the directory name of the one containing
c        trajectories.
c        Needed in t_dotraj
c     arguments:
c        arg1   char  input   the directory name.
c     history:
c        written by david n. bresch 931201
c-----------------------------------------------------------------------
c
      character*(80)  traj_dir,arg1
      common / traj_com / traj_dir
      integer         strbeg,strend
      data traj_dir / ' ' /
      
c
      traj_dir=arg1(strbeg(arg1):strend(arg1))
c
      return
      end ! of t_settraj_dir


      subroutine t_dotraj(arg1,arg2,arg3)
c---------------------------------------------------------------------------
c     purpose:
c        this routine governs the plotting of trajectories.
c        it calls t_plottraj to draw the trajectory.
c     arguments:
c        arg1   char  input   the filename of the netcdf-file containig
c                             trajectory-data.
c                             if arg1 contains no '/', it is just appended
c                             to the trajectory-directory stored in 
c                             traj_dir in the common block traj_com.
c                             if there are '/', arg1 is the full path.
c        arg2   char  input   the numbers of trajectories to be drawn.
c                             a string of the form '1/2/3/12/5' where the
c                             / separate numbers of the trajs to be drawn.
c                             or anly a '*' to draw all trajs.
c        arg3   char  input   the color of the trajectories.
c                             if not passed, the colors are varied cyclically.
c     history:
c        written by david n. bresch 931125
c-----------------------------------------------------------------------
c
c     argument declarations.
      character*(80)  arg1,arg2,arg3
c
      character*(80)  traj_dir
c
c     forward definition:
      integer         strbeg,strend,nint,ierr
      real            str2nmb
c
c     internal vars
      integer         begin,l,p,cyclic,ibeg,iend
      integer         color,cdfid,ierr,ntimes,ndim,nx,ny,nz
      integer         i,j,k,trajnumb,plot
      real            mdv
      integer         ntimax,nlevmax,maxparr  
c     ntimax: maximal number of timesteps
c             when changed: see also t_plottraj
c     nlevmax:maximal number of levels in vertical
c     maxparr:maximal number of entered trajectories selected    
      parameter(ntimax=97,nlevmax=32,maxparr=20)
c
c     internal dimensioned vars:
      real            tim(ntimax),stag(3),varmin(4),varmax(4)
      real            trax(ntimax),tray(ntimax),traz(ntimax)
      integer         vardim(4)
      integer         plotarr(maxparr)
c
c     common blocks:
      common / traj_com / traj_dir
c
      plot=0
c     do the conversion/verification of all arguments:
      if (arg2.ne.'#') then
c     extract the numbers of the lines to be plotted:
         if (arg2.eq.'*') then 
c     print all trajectories
            plot=2
            go to 99
         endif
         do p=1,maxparr
            plotarr(p)=0 ! clear plotarr
         enddo
c     search for the occurence of one '-' which defines a range:
         begin=1
         do l = strbeg(arg2), strend(arg2)
            if (arg2(l:l).eq.'-') then
               read(arg2(begin:l-1),'(i3)') plotarr(1)
               read(arg2(l+1:strend(arg2)),'(i3)') plotarr(2)
               do j=2,plotarr(2)-plotarr(1)+1
                  plotarr(j)=plotarr(1)+j-1
               enddo
               plot=1
c     overjump the search of '/'
               go to 99
            endif
         enddo
c
c     search for the occurence of '/' which separes numbers:
         p=1
         begin=1
         do l = strbeg(arg2), strend(arg2)
            if (arg2(l:l).eq.'/') then
               read(arg2(begin:l-1),'(i3)') plotarr(p)
               plot=1
               p=p+1
               begin=l+1
               if (p.gt.maxparr) then 
                  go to 99
               endif
            endif
         enddo
         read(arg2(begin:strend(arg2)),'(i3)') plotarr(p)
 99      continue
      else
         plot=0
      endif
      if (arg3.ne.'#') then
c     get the color-mode
         color=nint(str2nmb(arg3(1:strend(arg3)),ierr))
c         read(arg3(1:strend(arg3)),'(i5)') color
         cyclic=0
      else
         color=2
         cyclic=1
      endif
c
c     open netcdf-file:
      p=0
c     test, if arg1 contains any '/':
      do l = strbeg(arg1), strend(arg1)
         if (arg1(l:l).eq.'/') then
            p=1
         endif
      enddo
      if (p.ne.1) then
         if (traj_dir(strbeg(traj_dir):strend(traj_dir)).eq.' ') then
            print *,'use field=traj_dir[...] to set the directory'
            return
         else
         arg1=traj_dir(strbeg(traj_dir):strend(traj_dir))//'/'//
     &        arg1(strbeg(arg1):strend(arg1))
         endif
      endif
c
c     IVE passes only UPCASE strings, convert it:
      ibeg=strbeg(arg1)
      iend=strend(arg1)
      call locase(arg1(ibeg:iend),iend-ibeg+1)
      call cdfopn(arg1(1:strend(arg1)),cdfid,ierr)      
      if (ierr.ne.0) then
         return
      endif
c
c     get the timesteps:
      call gettimes(cdfid,tim,ntimes,ierr)
c
c     Get the data-dimensions etc.:
      call getdef(cdfid,'xpos',ndim,mdv,vardim,varmin,varmax,stag,ierr)
      nx=vardim(1)
      ny=vardim(2)
      nz=vardim(3)
c
c     get the trajectory-data:
      trajnumb=1
      if (plot.eq.0) then 
         print *,'Summary of all available trajectories:'
      endif
c     process every single trajectory:
      do k=1,nz
         do j=1,ny
             do i=1,nx
                call gettra(cdfid,'xpos',i,j,k,ntimes,trax,ierr)
                call gettra(cdfid,'ypos',i,j,k,ntimes,tray,ierr)
                call gettra(cdfid,'ppos',i,j,k,ntimes,traz,ierr)
                if (plot.eq.1) then
c     test, if this trajectory has to be plotted:
                   do p=1,maxparr
                      if (plotarr(p).eq.trajnumb) then
                         if (cyclic) then 
                             call t_cyclcol(color)
                         endif
                         call t_plottraj(trax,tray,traz,ntimes,color)
                      endif
                   enddo
                elseif (plot.eq.2) then
c     plot all trajectories
                   if (cyclic) then 
                      call t_cyclcol(color)
                   endif
                   call t_plottraj(trax,tray,traz,ntimes,color)
                else                                    
c     print the start-point out and number the trajectories:
                   print *,trajnumb,': ',trax(1),tray(1),traz(1),'to ',
     &                     trax(ntimes),tray(ntimes),traz(ntimes)
                endif
                trajnumb=trajnumb+1
             enddo
         enddo
      enddo
c
      return
      end ! of t_dotraj

      subroutine t_plottraj(trax,tray,traz,ntimes,color)
c---------------------------------------------------------------------------
c     Purpose:
c        This routine is used to draw trajectories on a map.
c     Arguments:
c        trax   real  input   array(ntimax), containing the x-positions.
c        tray   real  input   array(ntimax), containing the y-positions.
c        traz   real  input   array(ntimax), containing the p-positions.
c        ntimes integer inp   the number of points on the trajectory.
c        color  integer inp   the line color.
c     History:
c        written by David N. Bresch 931129
c-----------------------------------------------------------------------
c
c     Argument declarations.
      integer         ntimes,color
      integer         ntimax,nlevmax      
      parameter(ntimax=97,nlevmax=32)
      real            trax(ntimax),tray(ntimax),traz(ntimax)
c
c     Internal dimensioned vars:
      integer         i
      real            x1,y1,x2,y2
c
c     plot it:
      do i=1,ntimes-1 
         x1=trax(i)
         y1=tray(i)
         x2=trax(i+1)
         y2=tray(i+1)
         call t_drawline(x1,y1,x2,y2,color)
      enddo
c
      return
      end ! of t_plottraj   

      subroutine t_cyclcol(color)
c---------------------------------------------------------------------------
c     Purpose:
c        This routine selects the next color in the cycle.
c     Arguments:
c        color  integer in/out   the line color.
c               1:black,2:white,3:red,4:orange,5:yellow
c               6:green,7:blue,8:violet,9:black,10:grey
c     History:
c        written by David N. Bresch 931129
c-----------------------------------------------------------------------
      integer   color
c
c     Select the next color:
      color=color+1
c      print *,'color:',color
      if (color.gt.8) then
         color=3
      elseif (color.lt.3) then
         color=3
      endif
c
      return
      end ! of t_cyclcol
c
c
c------------------------------------------------------------------
c The following routine is from ~henry/lib/libcdfplus.f
c--------------------------------------------------------------------
c
      subroutine gettra(cdfid,varnam,ix,iy,iz,ntimes,array,ierr)
C------------------------------------------------------------------------
C
C     Reads the time-evolution for one grid-point of the variable
C     indicated by varnam.
C
C     cdfid     int     input   identifier for NetCDF file
C     varnam    char    input   name of variable
C     ix        int     input   x-index for values to read
C     iy        int     input   y-index for values to read
C     iz        int     input   z-index for values to read
C     ntimes    int     input   number of time-indices to read
C     array     real    output  array contains the readed values
C     ierr      int     output  error flag
C------------------------------------------------------------------------

C     Declaration of attributes

      integer   cdfid
      character*31 varnam
      integer   ix,iy,iz
      integer	ntimes
      real      array(ntimes)

C     Declaration of local variables

      integer   corner(4),edgeln(4)
      integer   idvar,ierr
      integer	ncvid

      corner(1)=ix
      corner(2)=iy
      corner(3)=iz
      corner(4)=1
      edgeln(1)=1
      edgeln(2)=1
      edgeln(3)=1
      edgeln(4)=ntimes

      idvar =ncvid(cdfid,varnam,ierr)
      call ncvgt(cdfid,idvar,corner,edgeln,array,ierr)
      if (ierr.ne.0) goto 991

      return
  991 stop 'Variable not found on NetCDF file in SR gettra'
      end
