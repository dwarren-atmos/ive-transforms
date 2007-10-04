      integer function fcalcfld (varnam, ndims, dims, stag, datmin, 
     &                           datmax, misdat, 
     &                           data_units, data_display_units) 
c-----------------------------------------------------------------------
c     $Id: fcalcfld.f,v 1.1 1994/11/14 22:37:06 warren Exp $
c     Purpose:
c        This is a user-written function that is used to calculate 
c        user-derived fields. 
c     Arguments:
c        varnam  char  input   the name of the variable to be derived.
c        ndims   int   output  the number of dimensions of the variable.
c                              Assumed to be 4, which is okay for any
c                              variable with not more than 4 dimensions.
c        dims    int   output  the number of data points along each 
c                              dimension. 
c                              THIS ARRAY MUST BE SUPPLIED BY THE USER.
c                              For example, if dims(1) = nx, 
c                              dims(2) = ny, dims(3), this is a 3D data
c                              set that is nx X ny X nz.
c        stag    real  output  the grid staggering along each dimension.
c                              THIS ARRAY MUST BE SUPPLIED BY THE USER.
c        datmin  real  output  the location of the origin of the data
c                              in physical space.
c                              THIS ARRAY MUST BE SUPPLIED BY THE USER.
c        datmax  real  output  the extent of the data in physical space
c                              along each dimension.
c                              THIS ARRAY MUST BE SUPPLIED BY THE USER.
c        misdat  real  output  the missing data value. Any point whose
c                              value is misdat will be ignored by the
c                              plotting routines.
c                              THIS VALUE MUST BE SUPPLIED BY THE USER.
c     Note : on output the value of the function is equal to either 0 or
c            some integer that points to memory where the derived field
c            is stored. If the function value is 0, it indicates to
c            the calling routine that there was some problem in 
c            calculating the field.
c     History:
c	$Log: fcalcfld.f,v $
c Revision 1.1  1994/11/14  22:37:06  warren
c Christoph Schaer's European Model transforms.
c
cRevision 1.40  1993/11/24  14:19:26  schaer
cSuche in INTER neu von oben nach unten.
c
cRevision 1.39  1993/11/19  15:36:03  schaer
cAdded function INTER[var1:var2:valvar2]
c
cRevision 1.38  1993/11/07  13:02:54  schaer
cCorrected bug in XMEAN (Adi Gamma).
c
cRevision 1.37  1993/10/25  10:12:56  schaer
cAenderung Andy Rossa in "initatt".
c
cRevision 1.36  1993/10/15  13:01:59  schaer
cAllow 4-dimensional arrays in GET.
c
cRevision 1.35  1993/10/07  13:36:18  schaer
cAdded new heading-routines from David Bresch.
c
cRevision 1.34  1993/08/31  08:12:08  schaer
cNeu organisierte Ableitung (subroutinene in File derivat.f), sowie
cmodifizierte Ableitungsfunktionen mit missing data checking (aehnlich
cwie frueher in den ecl-transformationen.
c
cRevision 1.33  1993/08/30  13:27:22  schaer
cKosmetik.
c
cRevision 1.31  1993/08/30  08:01:33  schaer
cNeu auch Theta-Koordinaten auf NETcdf-File erlaubt (datatype=3).
cAufnahme von ak und bk in /transcommon/.
c
cRevision 1.30  1993/08/27  12:39:06  schaer
cBug in 'gcdfvar' korrigiert.
c
cRevision 1.29  1993/08/27  12:00:01  schaer
cKleine Korrektur an horizontaler Ableitung.
c
cRevision 1.28  1993/08/27  11:40:42  schaer
cAb sofort ist derptype=1 der default.
c
cRevision 1.27  1993/08/27  09:05:56  schaer
cAdded 3-point vertical derivation (centered is default) from Andrea Rossa.
c
cRevision 1.26  1993/08/27  07:01:24  schaer
cNeue Funktion 'privat' fuer individuelle Aenderungen.
cNeue flags datatype und disptype zur Kennzeichnung der Daten (ersetzten
cdie alten flags prsrf und levtype, welche aber vorlaeufig noch
cvorhanden sind).
c
cRevision 1.25  1993/08/23  14:36:34  schaer
cNeue Philosophie mit PS: Falls nicht auf File, wird kein neuer
cButton generiert, sondern Verwaltung des Memories geschieht
cintern in ftr_read.
c
cRevision 1.24  1993/08/19  12:41:53  schaer
cGlobale Attribute der vertikalen Koordinate werden neu aus dem
cDatenfile gelesen, anstelle (1050,0) zu sein.
c
cRevision 1.23  1993/08/19  12:18:49  schaer
cVerwendung von freevar anstelle freemem.
c
cRevision 1.22  1993/08/18  06:14:37  schaer
cFehler in Memory-Return umschifft. Markiert mit c++. Vorlaeufig wird
cvom Benutzer kein Variablen-Memory zurueckgegeben (dies erfolgt erst
cverspaetet mit Garbage-Collect).
c
cRevision 1.21  1993/08/17  15:16:58  schaer
cNeue HELP-File Konvention.
c
cRevision 1.20  1993/08/06  06:41:04  schaer
cKorrektur eines Fehlers in der horizontalen Ableitung.
c(Andrea Rossa und Dani Luethi).
c
cRevision 1.19  1993/07/16  09:47:35  schaer
cEnhancement of MIMA-function.
c
cRevision 1.18  1993/07/13  11:43:45  schaer
cAdded function MIMA to ensure uniform data-range.
c
cRevision 1.17  1993/06/25  15:11:55  schaer
cFunktion TMAX zum Berechnen des Zeit-Maximums.
c
cRevision 1.16  1993/06/15  08:56:19  schaer
cImproved error checking on 0-returns from getvar.
c
cRevision 1.15  1993/06/03  14:31:37  schaer
cNeue Berechnung von aslay, etc in ftr_read.
c
cRevision 1.14  1993/06/02  09:26:18  schaer
cAdded code for Filter-Function.
c
cRevision 1.13  1993/06/01  12:56:48  schaer
cChanged definitions of COS, SIN, TAN to be associated with latitude
cin rotated coordinate system.
c
cRevision 1.12  1993/06/01  07:07:32  schaer
cRevision der Ableitung.
c
cRevision 1.10  1993/05/25  11:58:12  schaer
cChanges associated with read a var from a file.
c
cRevision 1.9  1993/05/14  12:04:22  schaer
cUmbenennen von U' > U, 2. Teil.
c
c---------------------------------------------------------------------- 



c     do compiler dependent definition of recursive variables
c#ifdef ultrix
c
c#else
      implicit automatic (a-z)
c#endif

c     include the common-block with the level information
      include 'constants.icl'

c     Argument declarations. MAXDIM is the maximum data dimension allowed,
c     LENVAR is the maximum length of the variable-name.
      integer   MAXDIM, LENVAR
      parameter (MAXDIM=4, LENVAR=80)
      character*(80)   data_units, data_display_units
      character*(LENVAR)  varnam
      integer          ndims, dims(MAXDIM)
      real             stag(MAXDIM), misdat,
     &                 datmin(MAXDIM), datmax(MAXDIM)

c     Local variable declarations
      character*(LENVAR) fld,fld1,fld2,arg1,arg2,arg3,dvar
      logical          err, local
      logical          errflg,localu,localv,local
      real             misdatu,misdatv,misdvel,misdatd,misdata
      real             delta1,delta2,str2nmb,afil,tstrt,tstop,umin,umax
      integer          ptru,ptrv,ptrt,ptrvel,ptrd,ptrrr,size,ptrdps,ierr
      integer          ptrh,i,k,ibeg,iend, ptr,ptr1,ptr2,dims1(MAXDIM)
      integer          ixsmin,ixsmax,iysmin,iysmax

      character*(80)   data1_units, data1_display_units
      character*(1)    del
      logical          local1
      integer          ndims1, dims1(MAXDIM), itime,ilevel,ired
      real             stag1(MAXDIM), misdat1,rred,valint,
     &                 datmin1(MAXDIM), datmax1(MAXDIM)
      real             datmin2(MAXDIM), datmax2(MAXDIM),
     &                 stag2(MAXDIM), misdat2
      integer          ndims2,dims2(MAXDIM)

c     External function declarations.
      integer        getmem, addvar, getvar, strbeg, strend, gcdfvar
      integer        andy, privat
      logical        isfunc

c     initially, set pointer 'fcalcfld' to unused, indicating a
c     failure in computing the requested field.
      fcalcfld=0

c     define field-name to be computed (removes leading and trailing blanks)
      ibeg = strbeg (varnam)
      iend = strend (varnam)
      fld=varnam(ibeg:iend)
      iend = strend (fld)

      print *,'Entering fcalcfld for ',varnam(ibeg:iend)

c     first attempt to provide the field from the user-written
c     function privat
      fcalcfld=privat(fld, ndims, dims, stag, datmin, 
     &                           datmax, misdat, 
     &                           data_units, data_display_units)
      if (fcalcfld.gt.0) then
         print *,'Field ',varnam(ibeg:iend),' is user-defined'
         return
      endif

c     start computation of fields

      if (fld.eq.'VEL') then
c       Compute the horizontal velocity (u^2+v^2)^0.5. In this
c       first example, this should be done by using the computational
c       package provided with IVE. Any variable which is on the
c       data-file (or which can be computed) can be assessed with
c       a call to getvar. This function provides the pointer to the
c       data as well as the attributes of the data-field. The variable
c       name (the first argument of getvar) must be defined according 
c       to the syntax for the 'field='command of IVE:
           fcalcfld = getvar('(U^2+V^2)^0.5', ndims, dims, stag, 
     &                       datmin, datmax, misdat,
     &                       data_units, data_display_units, local)
c       All intermediate operations (e.g. reading data from the data-file,
c       doing computations) are hidden to the user. If getvar returns
c       successfully (i.e. fcalcfld<>0), then a new button denoted 'vel'
c       will be installed after return from fcalcfld.

      else if (fld.eq.'THETA') then
c          get the temperature field
           ptrt  = getvar('T', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if (ptrt.eq.0) return

c          get memory
           do i=ndims+1,MAXDIM
             dims(i)=1
           enddo
           fcalcfld=getmem(dims(1)*dims(2)*dims(3)*dims(4))
           if (fcalcfld.eq.0) return

c          call routine to do the computation
           call pottemp(%val(fcalcfld),%val(ptrt),
     &                           dims(1),dims(2),dims(3),dims(4))

      else if (fld.eq.'THETAE') then
c          get the temperature field (PS is already available with)
           ptrt  = getvar('T', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           ptru  = getvar('QD', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if ((ptrt.eq.0).or.(ptru.eq.0)) return

c          get memory
           do i=ndims+1,MAXDIM
             dims(i)=1
           enddo
           fcalcfld=getmem(dims(1)*dims(2)*dims(3)*dims(4))
           if (fcalcfld.eq.0) return

c          call routine to do the computation
           call equpot(%val(fcalcfld),%val(ptrt),
     &                %val(ptru),dims(1),dims(2),dims(3),dims(4))

      else if (fld.eq.'PSNN') then
c          get the surface height
           ptru  = getvar('ZB', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
c          get the temperature field (PS is already available with)
           ptrt  = getvar('T', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if ((ptrt.eq.0).or.(ptru.eq.0)) return

c          get memory
           dims(3)=1
           do i=ndims+1,MAXDIM
             dims(i)=1
           enddo
           fcalcfld=getmem(dims(1)*dims(2)*dims(3)*dims(4))
           if (fcalcfld.eq.0) return

c          call routine to do the computation
           call pstonn(%val(fcalcfld),%val(ptrt),%val(ptru),
     &                        dims(1),dims(2),dims(3),dims(4))

      else if (fld.eq.'DIV') then

           fcalcfld=getvar('D[U:X]+D[V*COS:Y]/COS',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)


      else if (fld.eq.'VORT') then

           fcalcfld=getvar('D[V:X]-D[U*COS:Y]/COS',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if (fcalcfld.eq.0) return
c          assign units 
           data_units='1 / s'
           data_display_units='1.E-4 / s'
c          remove data from lowermost level
           misdat=1.
           call blnkbnd(%val(fcalcfld),misdat,1,
     &                             dims(1),dims(2),dims(3),dims(4) )


      else if (fld(1:5).eq.'VORT@') then
           arg2=fld(6:strend(fld))

           fcalcfld=getvar(
     &       'D[V@'//arg2(1:strend(arg2))//':XM]'//
     &       '-D[U@'//arg2(1:strend(arg2))//'*COS:YM]/COS',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if (fcalcfld.eq.0) return
c          assign units 
           data_units='1 / s'
           data_display_units='1.E-4 / s'


      else if (fld.eq.'F') then
           fcalcfld=getvar('0.000145444*SIN[PHI]', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if (fcalcfld.eq.0) return
c          assign units 
           data_units='1 / s'
           data_display_units='1.E-4 / s'
           

      else if (fld.eq.'PV') then

c          addvar und getvar sind in ~schaer/uwgap7/uars/getvar.F

           fcalcfld=getvar('-(VORT+F)*9.80665*D[THETA:P]',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if (fcalcfld.eq.0) return
c          assign units (PVU's)
           data_units='K m**2 / (kg s)'
           data_display_units='1.E-6  K m**2 / (kg s)'
c          remove data from lowermost level
           misdat=1.e-2
           call blnkbnd(%val(fcalcfld),misdat,1,
     &                             dims(1),dims(2),dims(3),dims(4) )


      else if (fld.eq.'RH') then
c          get the temperature field (PS is already available with)
           ptrt  = getvar('T', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           ptru  = getvar('QD', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           ptrv  = getvar('QW', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if ((ptrt.eq.0).or.(ptru.eq.0).or.(ptrv.eq.0)) return

c          get memory
           do i=ndims+1,MAXDIM
             dims(i)=1
           enddo
           fcalcfld=getmem(dims(1)*dims(2)*dims(3)*dims(4))
           if (fcalcfld.eq.0) return

c          call routine to do the computation
           call relhum(%val(fcalcfld),%val(ptrt),
     &       %val(ptru),%val(ptrv),dims(1),dims(2),dims(3),dims(4))

      else if (fld.eq.'P') then
           call initatt(ndims,dims, stag, datmin, datmax, misdat,
     &       data_units, data_display_units)
           stag(3)=-0.5

c          get memory
           fcalcfld=getmem(nx*ny*nz*nt)
           if (fcalcfld.eq.0) return

c          call routine to do the computation
           call pressure(%val(fcalcfld),stag(3),nx,ny,nz,nt)

      else if (fld.eq.'PLEV') then
           call initatt(ndims,dims, stag, datmin, datmax, misdat,
     &       data_units, data_display_units)
           stag(3)=0.

c          get memory
           fcalcfld=getmem(nx*ny*nz*nt)
           if (fcalcfld.eq.0) return

c          call routine to do the computation
           call pressure(%val(fcalcfld),stag(3),nx,ny,nz,nt)


      else if (fld.eq.'PSX') then
           fcalcfld=getvar('D[PS:XM]', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)

      else if (fld.eq.'PSY') then
           fcalcfld=getvar('D[PS:YM]', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)


      else if (isfunc(fld,'GET',arg1,arg2,arg3)) then
        if ((arg3.eq.'#').and.(arg2.ne.'#')) then
          fcalcfld=gcdfvar(arg1, arg2, ndims, dims, stag, 
     &                           datmin,datmax, misdat, 
     &                           data_units, data_display_units) 
        endif


      else if (isfunc(fld,'PUT',arg1,arg2,arg3)) then
c	 with PUT a variable arg2 is saved on the file arg1, the name
c	 of the variable can be changed to arg3. if arg3 is not defined
c	 the variable is saved under the name arg2.
c        the function 'isfunc(var,fun,arg1,arg2,arg3)' is a function
c        provided to the user by IVE. It checks whether the variable
c        'var' (string) is a function involving the function 'fun'
c        (string). The function PUT has three arguments, which are
c        returned in 'arg1',`arg2' and 'arg3'. The third argument is
c        optional.
         if ((arg1.eq.'#').or.(arg2.eq.'#')) then
           print *,'Illegal number of arguments in ',fld(1:iend)
           return
         endif

c        the integer function 'getvar' is a function provided to the
c        user by IVE. It gets the data-pointer to the selected variable
c        as well as its attributes. Here we have to get the variable
c        defined by 'arg2' (string).
         ptr = getvar(arg2, ndims, dims, stag,
     &                datmin, datmax, misdat,
     &                data_units, data_display_units, local)

c        if getvar returns with the pointer-value 0, then the variable was
c        not found. Return from subroutine
         if (ptr.eq.0) return

c        the variable 'arg2' was found and can be used
c	 if arg3 is not defined, use the same name to save the variable
         if (arg3.eq.'#') arg3=arg2

         call pcdfvar (arg1, arg3, ptr, ndims, dims,
     &                          stag, datmin, datmax, misdat, 
     &                          data_units, data_display_units) 

c        if the variable is a computed one, the memory is freed
         if (local) then
           call freevar(arg2)
         endif
         return

      else if ((fld.eq.'GZ0').or.(fld.eq.'ZB').or.(fld.eq.'BLA')
     &     .or.(fld.eq.'BTY')) then
        fcalcfld=gcdfvar(extfiln, fld, ndims, dims, stag, 
     &                           datmin,datmax, misdat, 
     &                           data_units, data_display_units) 
        return

      else if (fld.eq.'DERPTYPE1') then
c       definiere Typ der p-Ableitung (siehe constants.icl)
        derptype=1
        print *
        print *,'Switching to 3-point vertical dervative'
        print *
        return

      else if (fld.eq.'DERPTYPE2') then
c       definiere Typ der p-Ableitung (siehe constants.icl)
        derptype=2
        print *
        print *,'Switching to 2-point centered vertical dervative'
        print *
        return

      else if (isfunc(fld,'D',arg1,arg2,arg3)) then
c       Compute the derivative for a 3d field on pressure surfaces.
c       Consider the syntax D(<var>,<dir>) for the computation
c       of the derivative of the variable <var> into the direction
c       <dir>, where dir='X','Y','XM','YM','P'  or 'T'. Here 'XM' and
c       'YM' denote horizontal derivatives on model surfaces, while
c       'X' and 'Y' are derivatives on pressure surfaces (for datatype<3).
c       Vertical derivatives are obtained with dir='P' (for datatype<3) and
c       dir='T' (for datatype=3).
c       
c       Check Syntax with the user-callable logical  function 'isfunc' 
c       (see the function itself for full description). This function
c       tests for the function-syntax and returns the  arguments in 
c       arg1 .. arg3. A Maximum of 3 arguments is allowed for, 
c       unused arguments are returned as '#'.

c       First check for the proper number of arguments (here 2)
           if ((arg3.ne.'#').or.(arg2.eq.'#')) then
             print *,'Illegal number of arguments in ',fld(1:iend)
             return
           endif
c       Check if arg2 is consistent
           if (.not. (
     &         (arg2.eq.'X' ).or.(arg2.eq.'Y' ).or.
     &         (arg2.eq.'XM').or.(arg2.eq.'YM').or.
     &         ((arg2.eq.'P').and.(datatype.ne.3)).or.
     &         ((arg2.eq.'T').and.(datatype.eq.3)) ) ) then
               print *,'Illegal argument ',
     &            arg2(strbeg(arg2):strend(arg2)),' in ',fld(1:iend)
               return
           endif
c       Initialize the pointers 
           ptr=0
           ptrd=0
           ptrdps=0
c       get the field which the derivative is to be taken from 
           ptr = getvar(arg1, ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
           if (ptr.eq.0) goto 961
           misdata=misdat
c       get the memory for the result
           do i=ndims+1,4
             dims(i)=1
           enddo
           if ((dims(3).lt.2).and.(arg2.eq.'P')) goto 961
           ptrd  = getmem (dims(1)*dims(2)*dims(3)*dims(4))
           if (ptrd.eq.0) goto 961
c       call the routine to do the computation. The code depends on whether
c       the field is 2d or 3d, and on whether the data is on pressure
c       surfaces
           if ( (arg2.eq.'P').or.(arg2.eq.'T') ) then
c            -------------------
c            vertical derivation
c            -------------------
             call der1dm(%val(ptr),%val(ptrd),arg2,
     &           dims(1),dims(2),dims(3),dims(4),
     &           datmin,datmax,stag,misdata,misdat)
             fcalcfld=ptrd
           else if ((datatype.ne.1).or.(dims(3).eq.1).or.
     &       (arg2.eq.'XM').or.(arg2.eq.'YM')) then
c            -----------------------------------
c            derivation on data (model) surfaces
c            -----------------------------------
             call der2dm(%val(ptr),%val(ptrd),arg2,
     &           dims(1),dims(2),dims(3)*dims(4),
     &           datmin,datmax,misdata,misdat)
             fcalcfld=ptrd
           else
c            ---------------------------------------
c            derivation on p-sufaces (data on sigma)
c            ---------------------------------------
c            first get the respective derivative of the surface pressure
c            print *,'getting PSX or PSY
             dvar='PS'//arg2(1:1)
             ptrdps = getvar(dvar, ndims1, dims1, stag1, 
     &                  datmin1, datmax1, misdat1,
     &                  data1_units, data1_display_units, local1)
             if (ptrdps.eq.0) then
               call freemem(ptrd)
               fcalcfld=0
               return
             endif
             if ((dims1(1).ne.dims(1)).or.(dims1(2).ne.dims(2)).or.
     &         (dims1(4).ne.dims(4))) then
               print *,'The fields ',arg1(strbeg(arg1):strend(arg1)),
     &           ' and PS are incompatible.'
               call freemem(ptrd)
               if (local1) call freevar(dvar)
               fcalcfld=0
               return
             endif
             if (ptrdps.eq.0) goto 961

c            Call the subroutine to compute the derivative.
             call der3d(%val(ptr),%val(ptrd),arg2,
     &              dims(1),dims(2),dims(3),dims(4),
     &              %val(ptrps),%val(ptrdps),datmin,datmax,stag)
             fcalcfld=ptrd
           endif
c       return unused and temporary memory, return to calling subroutine
  961      continue
           if ((ptr.gt.0).and.(local)) call freevar(arg1)
           return

      else if ((fld.eq.'RR').or.(fld.eq.'RRK').or.(fld.eq.'RRS')) then
c        total precipitation rate
         if (fld.eq.'RRK') then
           arg1='RRKN+RSKN'
         else if (fld.eq.'RRS') then
           arg1='RRSN+RSSN'
         else
           arg1='RRSN+RSSN+RRKN+RSKN'
         endif
         ptrrr=getvar(arg1, ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, local)
         if (ptrrr.eq.0) goto 963
c        Get the memory for the result
         do i=ndims+1,4
           dims(i)=1
         enddo
         if (dims(4).lt.2) goto 963
         ptrd  = getmem (dims(1)*dims(2)*dims(3)*dims(4))
         if (ptrd.eq.0) goto 963
c        Call the subroutine to compute the precipitation rate.
         call dertime(%val(ptrrr),%val(ptrd),
     &                 dims(1),dims(2),dims(3),dims(4),ierr)
         if (ierr.eq.0) then
           fcalcfld=ptrd
         else
           call freemem(ptrd)
         endif
c        return unused and temporary memory, return to calling subroutine
 963     continue
         if ((ptrrr.gt.0).and.(local)) call freevar(arg1)
         return


      else if (isfunc(fld,'INTER',arg1,arg2,arg3)) then
c       interpolates the field arg1 onto the surface arg2=arg3, where
c       arg2 is another field, and arg3 a real number

c       Check for the proper number of arguments
        if ((arg1.eq.'#').or.(arg2.eq.'#').or.(arg3.eq.'#')) then
          print *,'Illegal number of arguments in ',fld(1:iend)
          return
        endif

c       Initialize the pointers and error-flags
        ptr1=0
        ptr2=0

c       get the field arg2
        ptr2 = getvar(arg2, ndims2, dims2, stag2, 
     &                datmin2, datmax2, misdat2,
     &                data_units, data_display_units, local)
        if (ptr2.eq.0) return

c       get the field arg1
        ptr1 = getvar(arg1, ndims, dims, stag, 
     &                datmin, datmax, misdat,
     &                data_units, data_display_units, local)
        if (ptr1.eq.0) return

c       check whether the fields are compatible
        do i = ndims2+1,MAXDIM
          dims2(i)=1
        enddo
        do i = ndims+1,MAXDIM
          dims(i)=1
        enddo
        if ( (dims(1).ne.dims2(1)).or.(dims(1).ne.dims2(1)).or.
     &       (dims(4).ne.dims2(4)).or.
     &       (datmin(1).ne.datmin2(1)).or.(datmax(1).ne.datmax2(1)).or.
     &       (datmin(2).ne.datmin2(2)).or.(datmax(2).ne.datmax2(2)).or.
     &       (stag(1).ne.stag2(1)).or.(stag(2).ne.stag2(2))) then
             print *,'The fields in ',fld(1:iend),' have incompatible',
     &               ' attributes.'
             return
        endif

c       get the contour of arg2 to be interpolated to
        valint=str2nmb(arg3,ierr) 
        if (ierr.eq.1) then
          print *,arg3(strbeg(arg3):strend(arg3)),' is not a number'
          return

        endif

c       allocate memory
        size=dims(1)*dims(2)*dims(4)
        fcalcfld  = getmem (size)
        if (fcalcfld.eq.0) then
          print *,'Could not allocate memory in fcalcfld'
          return
        endif

        call int2iso (%val(ptr1),%val(ptr2),valint,%val(fcalcfld),
     &                 dims(1),dims(2),dims(3),dims2(3),dims(4),
     &                 datmin ,datmax ,stag ,misdat ,
     &                 datmin2,datmax2,stag2,misdat2)
        dims(3)=1
        misdat=misdat2
        return

      else if (fld.eq.'TESTFUN') then
 111     continue
         print *
         print *,'enter expression'
         read *,fld1
         if ((fld1(1:3).eq.'END').or.(fld1(1:3).eq.'end')) return
         print *,'enter function-name'
         read *,fld2
         local=isfunc(fld1,fld2,arg1,arg2,arg3)
         if (local) then
           print *,'   is a function with args'
           print *,'   ',arg1(1:strend(arg1)),' ',
     &                   arg2(1:strend(arg2)),' ',
     &                   arg3(1:strend(arg3))
         else
           print *,'   is not a function'
         endif
         goto 111

      else if ((fld.eq.'GRIDX').or.(fld.eq.'GRIDY').or.
     &         (fld.eq.'RLAT').or.(fld.eq.'RLON').or.
     &         (fld.eq.'RLAM').or.(fld.eq.'RPHI').or.
     &         (fld.eq.'LAT').or.(fld.eq.'LON').or.
     &         (fld.eq.'LAM').or.(fld.eq.'PHI')) then
         
c        construct grid-array and attributes
         ndims=2
         dims(1)=nx
         dims(2)=ny
         dims(3)=1
         dims(4)=1
         data_units=' '
         data_display_units=' '
         datmin(1)=lonmin
         datmax(1)=lonmax
         datmin(2)=latmin
         datmax(2)=latmax
         stag(1)=0.
         stag(2)=0.
         misdat=0.0

c        allocate memory
         size=1
         do i = 1, ndims
            size = size * dims(i)
         enddo
         fcalcfld  = getmem (size)
         if (fcalcfld.eq.0) then
           print *,'Could not allocate memory in fcalcfld'
           return
         endif

         call cregrid (fld,%val(fcalcfld),dims(1),dims(2),datmin,datmax)


      else if (fld.eq.'TAN') then
c       compute trionometric function of latitude
        fcalcfld = getvar('SIN[RPHI]/COS[RPHI]',ndims, dims, stag, 
     &                    datmin, datmax, misdat,
     &                    data_units, data_display_units, 
     &                    local)

      else if (fld.eq.'SIN') then
c       compute trionometric function of latitude
        fcalcfld = getvar('SIN[RPHI]',ndims, dims, stag, 
     &                    datmin, datmax, misdat,
     &                    data_units, data_display_units, 
     &                    local)

      else if (fld.eq.'COS') then
c       compute trionometric function of latitude
        fcalcfld = getvar('COS[RPHI]',ndims, dims, stag, 
     &                    datmin, datmax, misdat,
     &                    data_units, data_display_units, 
     &                    local)

      else if (fld.eq.'COS_MD') then
c       compute trionometric function of latitude
        fcalcfld = getvar('COS[RPHI]',ndims, dims, stag, 
     &                    datmin, datmax, misdat,
     &                    data_units, data_display_units, 
     &                    local)
        call addmiss(%val(fcalcfld),dims(1)*dims(2),0.,misdat)

      else if (isfunc(fld,'TMAX',arg1,arg2,arg3)) then
c       Compute extremum in time

c       Check for the proper number of arguments, convert arguments
        if (arg1.eq.'#') then
          print *,'Illegal number of arguments in ',fld(1:iend)
          return
        endif
        if (arg2.ne.'#') then
          tstrt=str2nmb(arg2,ierr)
          if (ierr.ne.0) then
            print *,'Inconsistent 2nd argument ',arg2(1:strend(arg2)),
     &                 ' in ',fld(1:iend)
            return
          endif
        endif
        if (arg3.ne.'#') then
          tstop=str2nmb(arg3,ierr)
          if (ierr.ne.0) then
            print *,'Inconsistent 3rd argument ',arg3(1:strend(arg2)),
     &                 ' in ',fld(1:iend)
            return
          endif
        endif

c       Next get the field
        ptr = getvar(arg1, ndims, dims1, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, local)
        if (ptr.eq.0) return
        if ((ndims.lt.4).or.(dims1(4).eq.1)) then
            print *,fld(1:iend),' is not timedependent'
            return
        endif
      
c       create the tar-field
        dims(1)=dims1(1)
        dims(2)=dims1(2)
        dims(3)=dims1(3)
        dims(4)=1
        ndims=3
        ptrd  = getmem (dims(1)*dims(2)*dims(3))
        if (ptrd.eq.0) return

c       Call the subroutine to compute the max
        call tmax(%val(ptr),%val(ptrd),nint(tstrt),nint(tstop),
     &                  dims(1),dims(2),dims(3),dims1(4),misdat)
        fcalcfld=ptrd
        misdat=misdatd

      else if (isfunc(fld,'TMEAN',arg1,arg2,arg3)) then
c       Compute time mean

c       Check for the proper number of arguments
        if ((arg1.eq.'#').or.(arg2.ne.'#')) then
          print *,'Illegal number of arguments in ',fld(1:iend)
          return
        endif

c       Initialize the pointers and error-flags
        ptr=0
        ptrd=0

c       Next get the field which the time-mean is to be taken from 
        ptr = getvar(arg1, ndims, dims1, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, local)
        if (ptr.eq.0) return
        
c       Define the dimension of the mean and get the memory for the result
        dims(1)=dims1(1)
        dims(2)=dims1(2)
        dims(3)=dims1(3)
        dims(4)=1
        size=1
        do i = 1, ndims
          size = size * dims(i)
        enddo
        ptrd  = getmem (size)
        if (ptrd.eq.0) return

c       Call the subroutine to compute the time mean.
        call tmean(%val(ptr),%val(ptrd),
     &             dims1(1),dims1(2),dims1(3),dims1(4),misdat)
        fcalcfld=ptrd


      else if (isfunc(fld,'XMEAN',arg1,arg2,arg3)) then
c       Compute zonal mean

c       Check for the proper number of arguments
        if ((arg1.eq.'#').or.(arg2.ne.'#')) then
          print *,'Illegal number of arguments in ',fld(1:iend)
          return
        endif

c       Initialize the pointers and error-flags
        ptr=0
        ptrd=0

c       Next get the field which the derivative is to be taken from 
        ptr = getvar(arg1, ndims, dims1, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, local)
        if (ptr.eq.0) return
        
c       Define the dimension of the mean and get the memory for the result
        dims(1)=1
        dims(2)=dims1(2)
        dims(3)=dims1(3)
        dims(4)=dims1(4)
        size=1
        do i = 1, ndims
          size = size * dims(i)
        enddo
        ptrd  = getmem (size)
        if (ptrd.eq.0) return

c       Call the subroutine to compute the zonal mean.
        call xmean(%val(ptr),%val(ptrd),
     &             dims1(1),dims(2),dims(3),dims(4),misdat)
        fcalcfld=ptrd


      else if (isfunc(fld,'MIMA',arg1,arg2,arg3)) then
c       Overwrite upper-left and lower-right corner with the values
c       provided in arg2 and arg3.

c       Check for the proper number of arguments
        if ((arg1.eq.'#').or.(arg2.eq.'#')) then
          print *,'Illegal number of arguments in ',fld(1:iend)
          return
        endif

c       Initialize the pointers and error-flags
        ptr=0

c       Next get the field which the derivative is to be taken from 
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, local)
        if (ptr.eq.0) return
        
c       Overwrite these values
        umin=str2nmb(arg2,ierr) 
        if (ierr.eq.1) return
        umax=str2nmb(arg3,ierr) 
        if (ierr.eq.1) return

c       allocate memory
        size=1
        do i = 1, ndims
           size = size * dims(i)
        enddo
        fcalcfld  = getmem (size)
        if (fcalcfld.eq.0) then
          print *,'Could not allocate memory in fcalcfld'
          return
        endif

        call minmax(%val(ptr),%val(fcalcfld),umin,umax,
     &             dims(1),dims(2),dims(3),dims(4))
        return

      else if (isfunc(fld,'STAGX',arg1,arg2,arg3)) then
c       Compute staggered field

c       Check for the proper number of arguments
        if ((arg1.eq.'#').or.(arg2.ne.'#')) then
          print *,'Illegal number of arguments in ',fld(1:iend)
          return
        endif

c       Initialize the pointers and error-flags
        ptr=0
        ptrd=0

c       Next get the field which the derivative is to be taken from 
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, local)
        if (ptr.eq.0) return
        
c       Check whether the field is suited for staggering
        if (dims(1).eq.1) return

c       Define the attributes of the staggered field and get memory
        delta1=(datmax(1)-datmin(1))/real(dims(1)-1)
        dims(1)=dims(1)-1
        datmin(1)=datmin(1)+delta1/2.
        datmax(1)=datmax(1)-delta1/2.
        stag(1)=stag(1)+0.5
        size=1
        do i = 1, ndims
          size = size * dims(i)
        enddo
        ptrd  = getmem (size)
        if (ptrd.eq.0) return

c       Call the subroutine to do the staggering
        call dostagx(%val(ptr),%val(ptrd),
     &             dims(1),dims(2),dims(3),dims(4),misdat)
        fcalcfld=ptrd


      else if (isfunc(fld,'FLT',arg1,arg2,arg3)) then
c       apply diffusion operator onto the horizontal levels of
c       field 'arg1'.
        
c       check for consistent number of arguments
        if ((arg3.ne.'#').or.(arg1.eq.'#')) then
          print *,' Inconsistent number of arguments in ',
     &              fld(strbeg(fld):strend(fld)),'.'
          return
        endif     

c       get the filter-constant
        if (arg2.eq.'#') then
          afil=1.
        else
          afil=str2nmb(arg2,ierr)
        endif
        print *,'AFIL=',afil

c       get variable
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, local)
        if (ptr.eq.0) return
        do i=ndims+1,MAXDIM
          dims(i)=1
        enddo
        print *,'Got var'

c       get memory
        ptrd=getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (ptrd.eq.0) return
        print *,'Got mem'

c       call the routine to do the computation
        call filt4d (%val(ptr),%val(ptrd),
     &           dims(1),dims(2),dims(3),dims(4),afil,misdat,ierr)
        if (ierr.eq.1) then
          call freemem(ptrd)
          return
        else
          fcalcfld=ptrd
          return
        endif


      else if (isfunc(fld,'RED',arg1,arg2,arg3)) then
c       reduce the resolution of field 'arg1' by an integer-factor given by arg2
        
c       check for consistent number of arguments
        if ((arg3.ne.'#').or.(arg2.eq.'#')) then
          print *,' Inconsistent number of arguments in ',
     &              fld(strbeg(fld):strend(fld)),'.'
          return
        endif

c       get the reduction-factor
        rred=str2nmb(arg2,ierr)
        if (ierr.eq.1) then
          print *,' Inconsistent argument ',
     &              arg2(strbeg(arg2):strend(arg2)),' in ',
     &              fld(strbeg(fld):strend(fld)),'.'
          return
        else 
          ired=nint(rred)
          print *,'Reduction of resolution by factor ',ired
        endif

c       get variable
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, local)
        if (ptr.eq.0) return
        do i=ndims+1,MAXDIM
          dims(i)=1
        enddo

c       get memory
        ptrd=getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (ptrd.eq.0) return

c       call the routine to do the computation
        call red4d (%val(ptr),%val(ptrd),
     &           dims(1),dims(2),dims(3),dims(4),ired,misdat,ierr)
        fcalcfld=ptrd




      else if (isfunc(fld,'SUM',arg1,arg2,arg3)) then
c       take the average value over some region
        
c       check for consistent number of arguments
        if ((arg2.ne.'#').or.(arg1.eq.'#')) then
          print *,' Inconsistent number of arguments in ',
     &              fld(strbeg(fld):strend(fld)),'.'
          return
        endif

c       get variable
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, local)
        if (ptr.eq.0) return
        do i=ndims+1,MAXDIM
          dims(i)=1
        enddo

c       enter the domain
        print *,'Enter domain (ixmin,ixmax,iymin,iymax): '
        read  *, ixsmin,ixsmax,iysmin,iysmax

c       get memory
        ptrd=getmem(dims(3)*dims(4))
        if (ptrd.eq.0) return

c       call the routine to do the computation
        call sum4d (%val(ptr),%val(ptrd),
     &           dims(1),dims(2),dims(3),dims(4),
     &           ixsmin,ixsmax,iysmin,iysmax,misdat,ierr)
        fcalcfld=ptrd
        dims(1)=1
        dims(2)=1


      else if (isfunc(fld,'ANALYZE',arg1,arg2,arg3)) then
c        the function 'isfunc(var,fun,arg1,arg2,arg3)' is a function 
c        provided to the user by IVE. It checks whether the variable 
c        'var' (string) is a function involving the function 'fun' 
c        (string). The function ANALYZE has one single argument, which
c        is returned in 'arg1'. Check whether the specified number of
c        arguments is 1.
         if ((arg1.eq.'#').or.(arg2.ne.'#')) then
           print *,'Illegal number of arguments in ',fld(1:iend)
           return
         endif

c        the integer function 'getvar' is a function provided to the
c        user by IVE. It gets the data-pointer to the selected variable 
c        as well as its attributes. Here we have to get the variable 
c        defined by 'arg1' (string).
         ptr = getvar(arg1, ndims, dims, stag, 
     &                datmin, datmax, misdat,
     &                data_units, data_display_units, local)

c        if getvar returns with the pointer-value 0, then the variable was
c        not found. Return from subroutine
         if (ptr.eq.0) return

c        the variable 'arg1' was found and can be used
         print *
         print *  ,'ANALYSIS OF VARIABLE ',arg1(1:strend(arg1))
         print 810,'number of dimensions: ',ndims
         print 810,'dimensions            ',(dims(i),i=1,ndims)
         print 811,'staggering            ',(stag(i),i=1,ndims)
         print 811,'domain-min            ',(datmin(i),i=1,ndims)
         print 811,'domain-max            ',(datmax(i),i=1,ndims)
         print 811,'missing-data flag     ',misdat
         print *   ,'units                      ',
     &               data_units(1:strend(data_units))
         print *   ,'display units              ',
     &               data_display_units(1:strend(data_display_units))
         print *   ,'locally computed           ',local
  810    format(a,4i12)
  811    format(a,4f12.3)

c        compute the total dimension of the data-array
         size=1
         do i = 1, ndims
            size = size * dims(i)
         enddo

c        the variable 'ptr' is the pointer (type integer) to the data-
c        array. To access the data, it has to be transfer to a subroutine
c        using the %val-function, in order to make it appear as an
c        arbitrary array.
         call anadata(%val(ptr),size,misdat)

c        if the variable is a computed one, the memory is freed
         if (local) call freevar(arg1)
         return

      else if (fld.eq.'GRIDZ') then
         fld1='T'
         ptr = getvar(fld1,ndims, dims, stag, 
     &                    datmin, datmax, misdat,
     &                    data_units, data_display_units, 
     &                    local)
         if (ptr.eq.0) then
c          the field T was not found on the data-file, try another field
           fld1='M'
           ptr = getvar(fld1,ndims, dims, stag, 
     &                       datmin, datmax, misdat,
     &                       data_units, data_display_units, 
     &                       local)
           if (ptr.eq.0) return
         endif
         size=1
         do i = 1, ndims
            size = size * dims(i)
         enddo
         fcalcfld  = getmem (size)
         if (fcalcfld.eq.0) then
           print *,'Could not allocate memory in fcalcfld'
           return
         endif
c        misdat=0 switches off missing data checking
         misdat=0.0
         data_units=' '
         data_display_units=' '

         call levelar (%val(fcalcfld),dims(1),dims(2),dims(3),dims(4))
 
      elseif (fld.eq.'L_COORD') then
         levtype=0
         disptype=0
         print *
         print *,' NOW USING ARRAY-INDICES AS VERTICAL COORDINATE.'
         print *,' MAKE SURE TO RESELECT CURRENT BUTTON!.'
         print *

      elseif (fld.eq.'P_COORD') then
         if (ptrps.gt.0) then
           disptype=1
           levtype=1
           print *
           print *,' NOW USING PRESSURE AS VERTICAL COORDINATE.'
           print *,' MAKE SURE TO RESELECT CURRENT BUTTON!.'
           print *
         else
           print *
           print *,' WARNING: Unable to switch to pressure coordinates!'
           print *
         endif  

      elseif (fld.eq.'T_COORD') then
c       get pointer to theta-field
        ptrth=getvar ('THETA', ndims, dims, stag, 
     &                       datmin, datmax, misdat,
     &                       data_units, data_display_units, local)
        if ((ptrth.gt.0).and.(.not.(local))) then
          levtype=2
          disptype=2
          print *
          print *,' NOW USING THETA AS VERTICAL COORDINATE.'
          print *,' MAKE SURE TO RESELECT CURRENT BUTTON!.'
          print *
        else 
          if (ptrth.gt.0) then
            call freevar('THETA')
          endif
          print *
          print *,' WARNING: Unable to switch to theta coordinates!'
          print *,'          Compute first FIELD=THETA first!'
          print *
        endif  

      elseif ((fld.eq.'HELP').or.(fld.eq.'?')) then
c        Print a list of the diagnostic fields.
         call system('more /usr/local/lib/ive/help_em')

      else

c       teste ob ein Level oder eine Zeit gewuenscht wird
        call fnddel(fld,'%',arg1,arg2,del)
        if (del.eq.'%') then
           ptr = getvar(arg1,ndims, dims, stag, 
     &                       datmin, datmax, misdat,
     &                       data_units, data_display_units, 
     &                       local)
           if (ptr.eq.0) return

           if ((arg2(1:1).eq.'T').and.(strend(arg2).gt.1)) then
c            select time-index
             itime=0
             do i=2,strend(arg2)
               if ((ichar(arg2(i:i)).ge.ichar('0')).and.
     &             (ichar(arg2(i:i)).le.ichar('9'))) then
                 itime=itime*10+ichar(arg2(i:i))-ichar('0')
               else
                 return
               endif
             enddo
             if ((itime.lt.1).or.(itime.gt.dims(4))) then
               print *,'Illegal time-index specified in ',
     &                                fld(1:strend(fld))
               return
             endif
c            get memory
             size=dims(1)*dims(2)*dims(3)
             fcalcfld=getmem(size)
             if (fcalcfld.eq.0) return
c            extract time
             print *,'Field ',arg1(1:strend(arg1)),' at time IT=',itime
             call exttime(%val(fcalcfld),%val(ptr),
     &           dims(1),dims(2),dims(3),dims(4),itime)
c            define new attributes
             dims(4)=1
             ndims=3

           else if ((arg2(1:1).eq.'K').and.(strend(arg2).gt.1)) then
c            select level-index
             ilevel=0
             do i=2,strend(arg2)
               if ((ichar(arg2(i:i)).ge.ichar('0')).and.
     &             (ichar(arg2(i:i)).le.ichar('9'))) then
                 ilevel=ilevel*10+ichar(arg2(i:i))-ichar('0')
               else
                 return
               endif
               if (datatype.eq.1) ilevel=dims(3)-ilevel+1
               if ((ilevel.lt.1).or.(ilevel.gt.dims(3))) then
                 print *,'Illegal level specified in',
     &                                fld(1:strend(fld))
                 return
               endif
             enddo
c            get memory
             size=dims(1)*dims(2)*dims(4)
             fcalcfld=getmem(size)
             if (fcalcfld.eq.0) return
c            extract level
             print *,'Field ',arg1(1:strend(arg1)),' at level K=',ilevel
             call extlevel(%val(fcalcfld),%val(ptr),
     &           dims(1),dims(2),dims(3),dims(4),ilevel)
c            define new attributes
             dims(3)=1
             if (ndims.eq.3) ndims=2
           endif
        endif

      endif
    
c      if (fcalcfld.eq.0) then
c        print *,'Unable to compute ',varnam(ibeg:iend)
c      else
c        print *,'Leaving fcalcfld with ',varnam(ibeg:iend)
c        print *,'ndims,dims',ndims,(dims(i),i=1,ndims)
c      endif
      end


      subroutine pcdfvar (filnam, varnam, ptr, ndims, dims,
     &                          stag, datmin, datmax, misdat, 
     &                          data_units, data_display_units, error) 
c     ================================================================
c     puts the data associated with 'ptr' as variable 'varnam' to 
c     the cdf-file 'filnam'.

      include 'netcdf.inc'
      include 'constants.icl'

c     Argument declarations. MAXDIM is the maximum data dimension allowed,
c     LENVAR is the maximum length of the variable-name.
      integer   MAXDIM, LENVAR, MAXVAR
      parameter (MAXDIM=4, LENVAR=80, MAXVAR=80)
      character*(80)   data_units, data_display_units
      character*(*)    filnam,varnam
      integer          ndims, dims(MAXDIM), ptr, error
      real             stag(MAXDIM), misdat,
     &                 datmin(MAXDIM), datmax(MAXDIM)

c     local variables
      character*(80)   varlst(MAXVAR),cfn,char
      character*(1)    del
      integer          cdfid,cstid,ierr,i,idate(5),istdate(5)
      integer          nvar,datmin1(MAXDIM),datmax1(MAXDIM),ndims1
      real             diff

c     externals
      integer          strend


c     try to open the NetCDF file with name given by filnam
      call ncpopt(NCVERBOS)
      call opncdf(filnam(1:strend(filnam)),cdfid,datmin1,datmax1,ndims1,
     &               varlst,nvar,cfn,ierr)

      if (ierr.ne.0) then
c       if an error occured create the file and define the variable
        call fnddel(cstfiln(1:strend(cstfiln)),'/',char,cfn,del)
        call crecdf(filnam(1:strend(filnam)),cdfid,datmin,datmax,
     &             3,cfn(1:strend(cfn)),ierr)
        call putdef(cdfid,varnam,ndims,misdat,dims,datmin,datmax,
     &              stag,ierr)
        write(*,*)'defined variable ',varnam(1:strend(varnam))
        diff=0.
      else   ! check if variable is already defined
        do i=1,nvar
          if (varlst(i).eq.varnam) then
c           calculate time value for actual data for 3d-fields
            if (dims(4).eq.1) then
              call cdfopn(cfn(1:strend(cfn)),cstid,ierr)
              call getstart(cstid,istdate,ierr)
              call clscdf(cstid,ierr)
              idate(1)=starty
              idate(2)=startm
              idate(3)=startd
              idate(4)=starth
              idate(5)=0
              call timediff(idate,istdate,diff)
c              write(*,*)'timediff ist ',diff
            endif
            goto 201
          endif
        enddo
        call putdef(cdfid,varnam,ndims,misdat,dims,datmin,datmax,
     &              stag,ierr)
        write(*,*)'defined variable ',varnam(1:strend(varnam))
        diff=0.
      endif
  201 continue
c     put the field on the file
      call putit(cdfid,varnam,%val(ptr),dims,diff)
      call clscdf(cdfid,ierr)
      end


      subroutine putit(cdfid,vnam,ar,dims,tval)
c     =========================================
c     writes the array ar on the NetCDF file with identifier
c     cdfid by calling the subroutine putdat (from the library
c     libcdfio) for every time step.
c     For 4d-fields tval is a dummy variable, for 3d-fields tval
c     gives the time value for the subroutine putdat.

c     include the common-block with the time-array information
      include 'attributes.icl'

      integer   cdfid
      integer   dims(4)
      real      ar(dims(1),dims(2),dims(3),dims(4))
      real      tval
      character*(*) vnam

      integer   ierr,i,ntimes
      integer   strend

      if (dims(4).eq.1) then
        call putdat(cdfid,vnam,tval,0,ar(1,1,1,1),ierr)
        write(*,*)'variable ',vnam(1:strend(vnam)),' at time ',
     &             tval, ' written'
      else
        do i=1,dims(4)
          call putdat(cdfid,vnam,timeval(i),0,ar(1,1,1,i),ierr)
          write(*,*)'variable ',vnam(1:strend(vnam)),' at time ',
     &               timeval(i), ' written'
        enddo
      endif

      return
      end


      subroutine getstart(cdfid,idate,ierr)
C------------------------------------------------------------------------
C     Purpose:
C       Get start date for fields on specified NetCDF file
C     Arguments:
C       cdfid   int     input   identifier for NetCDF file
C       idate   int     output  array contains date (year,month,day,time,step)
C                               dimensioned as idate(5)
C       ierr    int     output  error flag
C------------------------------------------------------------------------

      integer   ierr

      integer   cdfid
      integer   ncvid

      integer   idiyear,idimonth,ididay,iditime,idistep
      integer   idate(5)

      idiyear   =ncvid(cdfid,'starty',ierr)
      if (ierr.ne.0) return
      idimonth  =ncvid(cdfid,'startm',ierr)
      if (ierr.ne.0) return
      ididay    =ncvid(cdfid,'startd',ierr)
      if (ierr.ne.0) return
      iditime   =ncvid(cdfid,'starth',ierr)
      if (ierr.ne.0) return
      idistep   =ncvid(cdfid,'starts',ierr)
      if (ierr.ne.0) return

      call ncvgt1(cdfid,idiyear, 1,idate(1),ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,idimonth,1,idate(2),ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,ididay,  1,idate(3),ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,iditime, 1,idate(4),ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,idistep, 1,idate(5),ierr)
      if (ierr.ne.0) return

      end


      subroutine cdfopn(filnam,cdfid,ierr)
C------------------------------------------------------------------------

C     Opens the NetCDF file 'filnam' and returns its identifier cdfid.

C     filnam    char    input   name of NetCDF file to open
C     cdfid     int     output  identifier of NetCDF file
C     ierr      int     output  error flag
C------------------------------------------------------------------------

      include 'netcdf.inc'

      integer   cdfid,ierr
      character*(*) filnam

      integer   ncopn,strend

      call ncpopt(NCVERBOS)
      cdfid=ncopn(filnam(1:strend(filnam)),NCNOWRIT,ierr)

      return
      end


      subroutine timediff(date1,date2,diff)
C     =====================================
C     Calculates the time difference in hours for the two dates specified
C     by the two arrays date1 and date2. They are expected to contain the
C     following date information:
C     year      month   day     time    step.
C
C     date1     array specifying the first date
C     date2     array specifying the second date
C     diff      time differenc between date1 and date2 in hours
C
C     Warning:  ihdiff is equal to 0 for date1/2 = 880203_12_00 and
C               880202_12_24 !!!

      integer   date1(5),date2(5)
      integer   idays(12)       ! array containing the days of the monthes
      real      diff
      integer   ixday,ihdiff,iddiff,j

      data idays/31,28,31,30,31,30,31,31,30,31,30,31/

C     Determine if the period between date1 and date2 contains a Feb.29

      ixday=0   ! extra day flag

      if (((mod(date1(1),4).eq.0).and.(date1(1).ne.0)).or.
     >    ((mod(date2(1),4).eq.0).and.(date2(1).ne.0))) then
        if (((date1(2).le.2).and.(date2(2).gt.2)).or.
     >      ((date1(2).gt.2).and.(date2(2).le.2))) then
          ixday=1       ! set extra day flag to 1
          idays(2)=29   ! february has 29 days
        endif
      endif

      ihdiff=0  ! diff. in hours between date1/date2
      iddiff=0  ! diff. in days  between date1/date2

      if (date1(1).gt.date2(1)) then            ! compare years
        do j=date2(1),date1(1)-1
          iddiff=iddiff+365+ixday
        enddo
      else if (date1(1).lt.date2(1)) then
        do j=date1(1),date2(1)-1
          iddiff=iddiff-365-ixday
        enddo
      endif

      if (date1(2).gt.date2(2)) then            ! compare monthes
        do j=date2(2),date1(2)-1
          iddiff=iddiff+idays(j)
        enddo
      else if (date1(2).lt.date2(2)) then
        do j=date1(2),date2(2)-1
          iddiff=iddiff-idays(j)
        enddo
      endif

      iddiff=iddiff+date1(3)-date2(3)
      ihdiff=iddiff*24+date1(4)-date2(4)+date1(5)-date2(5)

      diff=real(ihdiff)

      return
      end


      integer function gcdfvar (filnam, varnam, ndims, dims,
     &                          stag, datmin, datmax, misdat, 
     &                          data_units, data_display_units) 
c     ============================================================
c     reads the variable 'varnam' from the cdf-file 'filnam'.

      include 'attributes.icl'

c     Argument declarations. MAXDIM is the maximum data dimension allowed,
c     LENVAR is the maximum length of the variable-name.
      integer   MAXDIM, LENVAR, MAXVAR
      parameter (MAXDIM=4, LENVAR=80, MAXVAR=80)
      character*(80)   data_units, data_display_units
      character*(*)    filnam,varnam
      integer          ndims, dims(MAXDIM)
      real             stag(MAXDIM), misdat,
     &                 datmin(MAXDIM), datmax(MAXDIM)

c     local variables
      character*(80)   varlst(MAXVAR),varhlp,cfn
      integer          cdfid,error,i
      integer          timerr,ierr,nvar,nrtimes
      real             rtimeval(200)

c     externals
      integer          strbeg,strend,getmem

c     initialize routine
      gcdfvar=0
      cdfid=0

c     make sure the file name is long enough
      if (strend(filnam).lt.3) return

c     open external file
      call opncdf (filnam,cdfid,datmin,datmax,ndims,
     &                                        varlst,nvar,cfn,error)
      if (error.ne.0) then
        print *,'Unable to open external file ',
     &             filnam(strbeg(filnam):strend(filnam))
        return
      else
        print *,'Opened external file ',
     &             filnam(strbeg(filnam):strend(filnam)),
     &             ' with the variables'
        do i=1,nvar
          varhlp=varlst(i)
          print *,'   ',varhlp(1:strend(varhlp))
        enddo
      endif

c     try to access the variable 'fld'
      call getdef (cdfid, varnam, ndims, misdat, 
     &                          dims, datmin, datmax, stag, error)
      if (error.eq.1) then
        print *,'Unable to read variable ',
     &             varnam(strbeg(varnam):strend(varnam))
        return
      else
        print *,'Reading variable ',
     &             varnam(strbeg(varnam):strend(varnam))
      endif

cc     define attributes of t-direction
c**    funktioniert aus unverstaendlichen gruenden nicht
c**    gettimes kehrt mit fehlermeldung zurueck
c      if ((ndims.eq.4).and.(dims(4).gt.1)) then
cC       get times-array
c        call gettimes(cdfid,rtimeval,nrtimes,ierr)
c        if (ierr.ne.0) then
c          print *,'Unable to read time-array on file ',
c     &             filnam(strbeg(filnam):strend(filnam))
c          print *,'ierr=',ierr
c          return
c        endif
cc       compare time-levels
c        timerr=0
c        if (ntimes.ne.nrtimes) timerr=1
c        do i=1,min(nrtimes,ntimes)
c          if (abs(rtimeval(i)-timeval(i)).gt.0.001) timerr=1
c        enddo
c        if (timerr.eq.1) 
c     &    print *,'WARNING: data has incompatible time-array'
cc       define time-attributes
c        stag(4)=0.
c        datmin(4)=rtimeval(1)
c        datmax(4)=rtimeval(min(nrtimes,ntimes))
c      endif

c     provisorische loesung welche halbwegs funktioniert
      if ((ndims.eq.4).and.(dims(4).gt.1)) then
        if (ntimes.ne.dims(4))
     &    print *,'WARNING: data has incompatible time-array'
        stag(4)=0.
        datmin(4)=timeval(1)
        datmax(4)=timeval(min(ntimes,dims(4)))
      endif

c     get the memory
      do i=ndims+1,MAXDIM
        dims(i)=1
      enddo
      gcdfvar=getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (gcdfvar.eq.0) return

c     read the data (in a later version, a call to getdat will be used
c     instead)
      call getcdf (cdfid, varnam, ndims, misdat, 
     &                 dims, datmin, datmax, stag, %val(gcdfvar), error)
      if (error.ne.0) then
        print *,'An error occured while reading variable ',
     &             varnam(strbeg(varnam):strend(varnam))
        return
      endif

c     exit routine
 900  continue
      call clscdf(cdfid,error)
      end


      subroutine exttime(art,ar,nx,ny,nz,nt,it)
c     =========================================
      integer nx,ny,nz,nt,it,i,j,k
      real    ar(nx,ny,nz,nt),art(nx,ny,nz)
      do i=1,nx
      do j=1,ny
      do k=1,nz
        art(i,j,k)=ar(i,j,k,it)
      enddo
      enddo
      enddo
      end


      subroutine extlevel(arl,ar,nx,ny,nz,nt,ik)
c     ==========================================
      integer nx,ny,nz,nt,ik,i,j,l
      real    ar(nx,ny,nz,nt),arl(nx,ny,nt)
      do i=1,nx
      do j=1,ny
      do l=1,nt
        arl(i,j,l)=ar(i,j,ik,l)
      enddo
      enddo
      enddo
      end


      subroutine filt4d (a,af,nx,ny,nz,nt,fil,misdat,ierr)
c     ====================================================
c     this subroutine applies the diffusion-operator onto the horizontal
c     levels of a 4d-array

c     argument declarations
      integer  nx,ny,nz,nt,ierr
      real     a(nx,ny,nz*nt),af(nx,ny,nz*nt),fil,misdat

c     variable declarations
      integer  k,ptrf1,ptrf2,i,nfil
      integer  getmem
      real     rfil

      nfil=int(fil-0.00001)
      rfil=fil-float(nfil)

c     get work-memory
      ptrf1=getmem((nx+1)*ny)
      ptrf2=getmem(nx*(ny+1))
      if ((ptrf1.eq.0).or.(ptrf2.eq.0)) then
        ierr=1
        write (*,*) ' Cannot get memory in filt4d'
        return
      else
        ierr=0
      endif

c     call the 2d filter routine for every level
      do k=1,nz*nt
        print *,'Calling filt2d r for k=',k
        call filt2d(a(1,1,k),af(1,1,k),%val(ptrf1),%val(ptrf2),
     &                                          nx,ny,rfil,misdat)
        do i=1,nfil
          print *,'Calling filt2d i for k=',k
          call filt2d(af(1,1,k),af(1,1,k),%val(ptrf1),%val(ptrf2),
     &                                          nx,ny,1.  ,misdat)
        enddo
      enddo

c     return work-memory
      call freemem(ptrf1)
      call freemem(ptrf2)
      end


      subroutine filt2d (a,af,f1,f2,nx,ny,fil,misdat)
c     ================================================
c     Apply a conservative diffusion operator onto the 2d field a,
c     with full missing data checking.

c     fil   real   inp  filter-amplitude, 0<afil<=1. Maximum smoothing
c                        with afil=1.

c     argument declaration
      integer     nx,ny
      real        a(nx,ny),af(nx,ny),f1(nx+1,ny),f2(nx,ny+1),fil,misdat

c     local variable declaration
      integer     i,j
      real        fh

c     compute constant fh
      fh=0.25*fil

c     compute fluxes in x-direction
      if (misdat.eq.0.) then
        do j=1,ny
        do i=2,nx
          f1(i,j)=a(i-1,j)-a(i,j)
        enddo
        enddo
      else
        do j=1,ny
        do i=2,nx
          if ((a(i,j).eq.misdat).or.(a(i-1,j).eq.misdat)) then
            f1(i,j)=0.
          else
            f1(i,j)=a(i-1,j)-a(i,j)
          endif
        enddo
        enddo
      endif
c     set boundary-fluxes to zero
      do j=1,ny
        f1(1,j)=0.
        f1(nx+1,j)=0.
      enddo

c     compute fluxes in y-direction
      if (misdat.eq.0.) then
        do j=2,ny
        do i=1,nx
          f2(i,j)=a(i,j-1)-a(i,j)
        enddo
        enddo
      else
        do j=2,ny
        do i=1,nx
          if ((a(i,j).eq.misdat).or.(a(i,j-1).eq.misdat)) then
            f2(i,j)=0.
          else
            f2(i,j)=a(i,j-1)-a(i,j)
          endif
        enddo
        enddo
      endif
c     set boundary-fluxes to zero
      do i=1,nx
        f2(i,1)=0.
        f2(i,ny+1)=0.
      enddo

c     compute flux-convergence -> filter
      if (misdat.eq.0.) then
        do j=1,ny
        do i=1,nx
            af(i,j)=a(i,j)+fh*(f1(i,j)-f1(i+1,j)+f2(i,j)-f2(i,j+1))
        enddo
        enddo
      else
        do j=1,ny
        do i=1,nx
          if (a(i,j).eq.misdat) then
            af(i,j)=misdat
          else
            af(i,j)=a(i,j)+fh*(f1(i,j)-f1(i+1,j)+f2(i,j)-f2(i,j+1))
          endif
        enddo
        enddo
      endif
      end



      real function str2nmb(chrs,ierr)
c---------------------------------------------------------------
c     This function takes the string chrs and converts it into a
c     real number. If the string is not a number, the routine
c     returns with str2numb=0 and ierr=1
c---------------------------------------------------------------

      character*(*) chrs
      integer       ierr

      integer       i,ibeg,iend,cp,ce,cf,strbeg,strend

c     Check whether string is a number
      ierr=0
      ibeg=strbeg(chrs)
      iend=strend(chrs)
      if (ibeg.gt.iend) return
      cp=0
      ce=0
      cf=1
      do i=ibeg,iend
c       check for + and -
        if ((cf.eq.1).and.
     &    ((chrs(i:i).eq.'-').or.(chrs(i:i).eq.'+'))) then
          cf=0
        else
          cf=0        
          if (chrs(i:i).eq.'.') cp=cp+1
c         check for format with exponents
          if ((chrs(i:i).eq.'E').or.(chrs(i:i).eq.'e')) then
            ce=ce+1
            cf=1
            if (cp.le.1) cp=1
          else
            if (((chrs(i:i).lt.'0').or.(chrs(i:i).gt.'9')).and.
     &        (chrs(i:i).ne.'.')) return
          endif
        endif
      enddo

      if ((cp.gt.1).or.(ce.gt.1)) goto 100

c     convert the string into a number
      read (chrs,*,err=100) str2nmb
      return

c     error-exit
 100  continue
      str2nmb=0.
      ierr=1
      end


      subroutine initatt(ndims,dims, stag, datmin, datmax, misdat,
     &                               data_units, data_display_units)
c-------------------------------------------------------------------
c     Purpose: 
c        Initializes the attributes of a variable as required for
c        a thermodynamic variable. The information is taken from
c        the constants common block.
c-------------------------------------------------------------------
      include       'constants.icl'
      include       'attributes.icl'

c     argument declaration
      integer       MAXDIM
      parameter     (MAXDIM=4)
      integer       ndims, dims(MAXDIM)
      real          stag(MAXDIM), datmin(MAXDIM), datmax(MAXDIM),misdat
      character*80  data_units, data_display_units

      ndims=4
      dims(1)=nx
      dims(2)=ny
      dims(3)=nz
      dims(4)=max0(1,max0(ntimes,nt))
      stag(1)=0.
      stag(2)=0.
      stag(3)=0.
      misdat=0.
      datmin(1)=lonmin
      datmax(1)=lonmax
      datmin(2)=latmin
      datmax(2)=latmax
      datmin(3)=presssrf
      datmax(3)=presstop
      datmin(4)=timeval(1)
      datmax(4)=timeval(ntimes)
      data_units=" "
      data_display_units=" "
      end 


      subroutine dertime(fld,der,ie,je,ke,le,ierr)
c     ============================================
c     computes the time-derivative of fld, results int der
      include 'attributes.icl'

c     argument declaration
      integer  ie,je,ke,le,ierr
      real     fld(ie,je,ke,le),der(ie,je,ke,le)

c     variable declaration
      integer  i,j,k,l,lp1,lm1
      real     dtr

c     Check whether ntimes in attributes agrees with le
      print *,'ntimes,le',ntimes,le
      if ((ntimes.eq.le).and.(ntimes.gt.1)) then
        ierr=0
      else
        ierr=1
        return
      endif

c     computation of time-derivative (in units 1/h)
      do l=1,le
        lp1=min0(le,l+1)
        lm1=max0(1 ,l-1)
        print *,'dt=',(timeval(lp1)-timeval(lm1))
        dtr=1./(timeval(lp1)-timeval(lm1))
        if (dtr.le.0.) then
          ierr=1
          return
        endif
        do k=1,ke
          do j=1,je
            do i=1,ie
              der(i,j,k,l)=dtr*(fld(i,j,k,lp1)-fld(i,j,k,lm1))
            enddo
          enddo
        enddo
      enddo
      end

  
      subroutine pottemp(pt,t,ie,je,ke,le)
c     ====================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     pt(ie,je,ke,le),t(ie,je,ke,le),getps

c     variable declaration
      integer  i,j,k,l
      real     rdcp,tzero,psrf
      data     rdcp,tzero /0.286,273.15/

c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation of potential temperature
      do l=1,le
        do i=1,ie
          do j=1,je
            psrf=getps(i,j,l)
            do k=1,ke
              pt(i,j,k,l)=(t(i,j,k,l)+tzero)*( (1000./prlay(k))**rdcp )
            enddo
          enddo
        enddo
      enddo
      end

      subroutine equpot(ap,t,qd,ie,je,ke,le)
c     ======================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     ap(ie,je,ke,le),t(ie,je,ke,le),getps
      real     qd(ie,je,ke,le)
      
c     variable declaration
      integer  i,j,k,l
      real     rdcp,tzero,psrf
      data     rdcp,tzero /0.286,273.15/

c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation of potential temperature
      do l=1,le
        do i=1,ie
          do j=1,je
            psrf=getps(i,j,l)
            do k=1,ke
              ap(i,j,k,l) = (t(i,j,k,l)+tzero)*(1000./prlay(k))
     +           **(0.2854*(1.0-0.28*qd(i,j,k,l)))*exp(
     +           (3.376/(2840.0/(3.5*alog(t(i,j,k,l)+tzero)-alog(
     +           100.*prlay(k)*max(1.0E-10,qd(i,j,k,l))/(0.622+0.378*
     +           qd(i,j,k,l)))-0.1998)+55.0)-0.00254)*1.0E3*
     +           max(1.0E-10,qd(i,j,k,l))*(1.0+0.81*qd(i,j,k,l)))
            enddo
          enddo
        enddo
      enddo
      end

      subroutine relhum(rh,t,qd,qw,ie,je,ke,le)
c     =========================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     rh(ie,je,ke,le),t(ie,je,ke,le),getps
      real     qd(ie,je,ke,le),qw(ie,je,ke,le)

c     variable declaration
      integer  i,j,k,l
      real     rdcp,tzero,psrf
      real     b1,b2w,b3,b4w,r,rd,gqd,ge
      data     rdcp,tzero /0.286,273.15/
      data     b1,b2w,b3,b4w,r,rd /6.1078, 17.2693882, 273.16, 35.86,
     &                  287.05, 461.51/

     
c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation of potential temperature
      do l=1,le
        do i=1,ie
          do j=1,je
            psrf=getps(i,j,l)
            do k=1,ke
              ge = b1*exp(b2w*(t(i,j,k,l))/(t(i,j,k,l)+b3-b4w))
              gqd= r/rd*ge/(prlay(k)-(1.-r/rd)*ge)
              rh(i,j,k,l)=(qd(i,j,k,l)+qw(i,j,k,l))/gqd
            enddo
          enddo
        enddo
      enddo
      end

      subroutine pstonn(pnn,t,zb,ie,je,ke,le)
c     =======================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     pnn(ie,je,le),t(ie,je,ke,le),getps
      real     zb(ie,je)

c     variable declaration
      integer  i,j,l
      real     rdcp,tzero,psrf,r,g
      real     ztstar,zalpha,zt0
      data     rdcp,tzero,r,g /0.286,273.15,287.05,9.80665/
     
c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation of potential temperature
      do l=1,le
        do i=1,ie
          do j=1,je
            psrf=getps(i,j,l)
            if (zb(i,j).lt.1.) then
              pnn(i,j,l)=getps(i,j,l)
            else
              ztstar = (t(i,j,1,l)+tzero)*(1. + 0.0065*r/g*
     &                   (getps(i,j,l)/prlay(1) - 1.0))
              zalpha = 0.0065*r
              zt0    = ztstar + 0.0065*zb(i,j)
              if (zt0.gt.290.5) then
                if (ztstar.gt.290.5) then
                  zalpha = 0.0
                  ztstar = 0.5*(ztstar+290.5)
                else
                  zalpha = r*(290.5-ztstar)/zb(i,j)
                endif
              else if (ztstar.lt.255.) then
                ztstar = 0.5*(255.0+ztstar)
              endif
              pnn(i,j,l) = getps(i,j,l)* exp(g*zb(i,j)/(r*ztstar)*
     &               (1.0 - 0.5*(zalpha*zb(i,j)/(r*ztstar)) +
     &                    0.333*(zalpha*zb(i,j)/(r*ztstar))**2))
            endif
          enddo
        enddo
      enddo
      end

      subroutine pressure(pr,stag3,ie,je,ke,le)
c     =========================================
      include 'constants.icl'

c     argument declaration
      integer  ie,je,ke,le
      real     pr(ie,je,ke,le),getps,stag3

c     variable declaration
      integer  i,j,k,l
      real     psrf

c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     computation pressure
      do l=1,le
        do i=1,ie
          do j=1,je
            psrf=getps(i,j,l)
            do k=1,ke
              if (stag3.eq.0.) then
                pr(i,j,k,l)=prlev(k)
              else
                pr(i,j,k,l)=prlay(k)
              endif
            enddo
          enddo
        enddo
      enddo
      end

      subroutine tmax(a,m,it1,it2,nx,ny,nz,nt,misdat)
c     ===============================================
c     Computes the temporal max of a field'

c     declaration of parameters
      real       a(nx,ny,nz,nt),m(nx,ny,nz),misdat
      integer    nx,ny,nz,nt,it1,it2

c     declaration of variables
      integer    i,j,k,l
      real       defval
      data       defval/-1.e15/

c     initialize field
      do k=1,nz
      do j=1,ny
      do i=1,nx
        m(i,j,k)=defval
      enddo
      enddo
      enddo

c     start computation
      print *,'TMAX mit l=',it1,it2
      print *,nx,ny,nz,nt
      do l=it1,it2
        print *,l
        do k=1,nz
        do j=1,ny
        do i=1,nx
          if ((a(i,j,k,l).ne.misdat).or.(misdat.eq.0.)) then
            m(i,j,k)=amax1(m(i,j,k),a(i,j,k,l))
          endif
        enddo
        enddo
        enddo
      enddo

c     clean up
      if (misdat.ne.0.) then
        do k=1,nz
        do j=1,ny
        do i=1,nx
          if (m(i,j,k).eq.defval) m(i,j,k)=misdat
        enddo
        enddo
        enddo
      endif
      end


      subroutine tmean(a,m,nx,ny,nz,nt,misdat)
c     ========================================
c     Computes zonal mean. Missing values in 'a' will only yield a missing value
c     in 'd' if more than 5% of the values required for the mean are missing.

c     declaration of parameters
      real       a(nx,ny,nz,nt),m(nx,ny,nz),misdat
      integer    nx,ny,nz,nt

c     declaration of variables
      integer    i,j,k,l,cnt
      real       sum

      print *,'Computing TMEAN for an array dimensioned',nx,ny,nz,nt
      print *,'misdat=',misdat

c     start computation
      do i=1,nx
        do k=1,nz
          do j=1,ny
            sum=0.
            cnt=0
            do l=1,nt
              if ((a(i,j,k,l).ne.misdat).or.(misdat.eq.0.)) then
                sum=sum+a(i,j,k,l)
                cnt=cnt+1
              endif
            enddo
            if (float(cnt)/float(nt).ge.0.95) then
              m(i,j,k)=sum/float(cnt)
            else
              m(i,j,k)=misdat
            endif
          enddo
        enddo
      enddo
      print *,'misdat=',misdat
      end


      subroutine xmean(a,m,nx,ny,nz,nt,misdat)
c     ========================================
c     Computes zonal mean. Missing values in 'a' will only yield a missing value
c     in 'd' if more than 5% of the values required for the mean are missing.

c     declaration of parameters
      real       a(nx,ny,nz,nt),m(ny,nz,nt),misdat
      integer    nx,ny,nz,nt

c     declaration of variables
      integer    i,j,k,l,cnt
      real       sum

      print *,'Computing XMEAN for an array dimensioned',nx,ny,nz,nt

c     start computation
      do l=1,nt
        do k=1,nz
          do j=1,ny
            sum=0.
            cnt=0.
            do i=1,nx
              if ((a(i,j,k,l).ne.misdat).or.(misdat.eq.0.)) then
                sum=sum+a(i,j,k,l)
                cnt=cnt+1
              endif
            enddo
            if (float(cnt)/float(nx).ge.0.95) then
              m(j,k,l)=sum/float(cnt)
            else
              m(j,k,l)=misdat
            endif
          enddo
        enddo
      enddo
      end


      subroutine minmax(a,d,umin,umax,nx,ny,nz,nt)
c     ============================================
c     Overwrite upper-left and lower-right corner with the values
c     provided in umin and umax. In addition, the new value is bound into
c     the interval [umin,umax]

c     declaration of parameters
      real       a(nx,ny,nz,nt),d(nx,ny,nz,nt),umin,umax
      integer    nx,ny,nz,nt

c     declaration of variables
      integer    i,j,k,l

      print *,'entering minmax'

c     start computation
      do l=1,nt
        do k=1,nz
          do i=1,nx
            do j=1,ny
              d(i,j,k,l)=amax1(umin,amin1(umax,a(i,j,k,l)))
            enddo
          enddo
          d(1,ny,k,l)=umin
          d(nx,1,k,l)=umax
        enddo
      enddo
      end


      subroutine dostagx(a,s,nx,ny,nz,nt,misdat)
c     ==========================================
c     Computes an x-staggered version of array 'a'.

c     declaration of arguments
      real       a(nx+1,ny,nz,nt),s(nx,ny,nz,nt),misdat
      integer    nx,ny,nz,nt
c     declaration of variables
      integer    i,j,k,l

c     start computation
      do l=1,nt
        do k=1,nz
          do j=1,ny
            do i=1,nx
              if ( ((a(i,j,k,l).ne.misdat).and.(a(i+1,j,k,l).ne.misdat))
     &          .or. (misdat.eq.0.) ) then
                s(i,j,k,l)=0.5*(a(i,j,k,l)+a(i+1,j,k,l))
              else
                s(i,j,k,l)=misdat
              endif
            enddo
          enddo
        enddo
      enddo
      end
        
      
      subroutine blnkbnd(ar,misdat,nlev,nx,ny,nz,nt)
c     ==============================================
c     overwrites the lowevermost nlev levels with misdat
      integer   nx,ny,nz,nt,nlev,i,j,k,l
      real      ar(nx,ny,nz,nt),misdat
      do l=1,nt
      do k=1,nlev
      do j=1,ny
      do i=1,nx
        ar(i,j,k,l)=misdat
      enddo
      enddo
      enddo
      enddo
      end


      subroutine cregrid (fld,ar,nx,ny,datmin,datmax)
c     ===============================================

      include 'rotpol.icl'
c     creates grid-array
      integer       nx,ny ,i,j
      character*(*) fld
      real          ar(nx,ny),datmin(2),datmax(2),dx,dy
      real          lon,lat,lonr,latr

      dx=(datmax(1)-datmin(1))/real(nx-1)
      dy=(datmax(2)-datmin(2))/real(ny-1)

      do j=1,ny
      do i=1,nx
        if (fld.eq.'GRIDX') ar(i,j)=real(i)
        if (fld.eq.'GRIDY') ar(i,j)=real(j)
        if (fld.eq.'RLON')  ar(i,j)=datmin(1)+real(i-1)*dx
        if (fld.eq.'RLAT')  ar(i,j)=datmin(2)+real(j-1)*dy
        if (fld.eq.'RLAM')  ar(i,j)=(datmin(1)+real(i-1)*dx)*zpir18
        if (fld.eq.'RPHI')  ar(i,j)=(datmin(2)+real(j-1)*dy)*zpir18
        if ((fld.eq.'LAT').or.(fld.eq.'LON')
     &  .or.(fld.eq.'PHI').or.(fld.eq.'LAM')) then
          lonr=datmin(1)+real(i-1)*dx
          latr=datmin(2)+real(j-1)*dy
          call ph2ll(latr,lonr,lat,lon)
          if (fld.eq.'LON') ar(i,j)=lon
          if (fld.eq.'LAT') ar(i,j)=lat
          if (fld.eq.'LAM') ar(i,j)=lon*zpir18
          if (fld.eq.'PHI') ar(i,j)=lat*zpir18
        endif
      enddo
      enddo
      end


      subroutine addmiss(ar,n,val,misdat)
c-----------------------------------------------------------------------
c     ersetzt die Werte 'val' in ar(n) mit misdat.
c-----------------------------------------------------------------------
      integer     n,i
      real        ar(n),misd,val,misdat,eps

c     define a missing data value
      if (misdat.eq.0.) then
        misd=999.99
      else
        misd=misdat
      endif

c     find the extremas of the array
      eps=1.e-4
      do i=1,n
        if (abs(ar(i)-val).lt.eps) then
          ar(i)=misd
          misdat=misd
        endif
      enddo
      end


      subroutine levelar(ar,nx,ny,nz,nt)
c     ==================================
      integer  nx,ny,nz,nt,nt,i,j,k,l
      real     ar(nx,ny,nz,nt)
      do i=1,nx
       do j=1,ny
        do k=1,nz
         do l=1,nt
           ar(i,j,k,l)=float(k)
         end do
        end do
       end do
      end do
      end



      subroutine anadata (ar,nn,misdat)
c     =================================
      integer   nn
      real      ar(nn),misdat

      integer   powmin, powmax
      parameter (powmin=-10,powmax=10)
      integer   i,lrgpow,smlpow,power(powmin:powmax)
      integer   missing,undef,neg,pos,zero,oom,pmin,pmax
      real      absmin,absmax,minval,maxval,log10,exp10,x,ln10

c     define statement-function for logarithm and exponent with base 10
      log10(x)=log(x)/ln10
      exp10(x)=exp(ln10*x)

c     initialize constants
      ln10=log(10.)
      absmin=exp10(float(powmin))
      absmax=exp10(float(powmax+1))

c     initialize counters
      missing=0
      undef=0
      zero=0
      neg=0
      pos=0
      lrgpow=0
      smlpow=0
      do i=powmin,powmax
        power(i)=0
      enddo
      maxval=-1.e20
      minval= 1.e20

      do i=1,nn
        if ((ar(i).eq.misdat).and.(misdat.ne.0.)) then
          missing=missing+1
        else if ((ar(i).lt.1.).or.(ar(i).gt.-1.)) then
c         do evaluation of extremum
          minval=amin1(minval,ar(i))
          maxval=amax1(maxval,ar(i))
c         do zero/neg/pos classification
          if (ar(i).eq.0.) zero=zero+1
          if (ar(i).gt.0.) pos=pos+1
          if (ar(i).lt.0.) neg=neg+1
c         do order of magnitude statistics
          if (abs(ar(i)).lt.absmin) then
            smlpow=smlpow+1
          else if (abs(ar(i)).ge.absmax) then
            lrgpow=lrgpow+1
          else
            x=abs(ar(i))
            oom=int(log10(x)+1000.)-1000
            power(oom)=power(oom)+1
          endif
        else
c         it depends on the compiler and the compiling-options, whether an
c         undefined value is accounted for in this loop
          undef=undef+1
        endif
      enddo
           
c     printout
      print 800,'number of data-points:',nn
      print 800,'missing:              ',missing
      print 800,'undefined:            ',undef
      print 800,'negative values:      ',neg
      print 800,'zero     values:      ',zero
      print 800,'positive values:      ',pos
      print 801,'minimum:              ',minval
      print 801,'maximum:              ',maxval
      pmin=powmax
      pmax=powmin
      do i=powmin,powmax
        if (power(i).gt.0) then
          if (i.lt.pmin) pmin=i
          if (i.gt.pmax) pmax=i
        endif
      enddo
      if (pmax.ge.pmin) then
        print *,'order of magnitude:'
        pmin=max0(pmin-1,powmin)
        pmax=min0(pmax+1,powmax)
        if (pmin.eq.powmin) print 811,powmin,smlpow
        do i=pmin,pmax
          print 810,i,i+1,power(i)
        enddo
        if (pmax.eq.powmax) print 812,powmax+1,lrgpow
      endif
      
  800 format(a,i12)
  801 format(a,e12.3)
  810 format('1.E',i3,' <=|x|< 1.E',i3,': ',i12)
  811 format(   '         |x|< 1.E',i3,': ',i12)
  812 format('1.E',i3,' <=|x|:         '   ,i12)

      end


      subroutine int2iso (ar1,ar2,valint,res,nx,ny,nz1,nz2,nt,
     &        datmin1,datmax1,stag1,misdat1,
     &        datmin2,datmax2,stag2,misdat2)
c-----------------------------------------------------------------------
c     Purpose:
c-----------------------------------------------------------------------
      
c     arguments
      integer      nx,ny,nz1,nz2,nt
      real         ar1(nx,ny,nz1,nt),ar2(nx,ny,nz2,nt),res(nx,ny,nt)
      real         datmin1(4),datmax1(4),stag1(4),misdat1
      real         datmin2(4),datmax2(4),stag2(4),misdat2
      real         valint

c     local variables
      integer      i,j,k,l,iflagp(4),iflagc(4),kk
      logical      missing
      real         compt(4),phypt(4),fract,misdath,v1,v2,vv1,vv2,int4dm

c     initialize iflag
      data         iflagp/1,1,1,1/
      data         iflagc/0,0,1,0/

c     make sure the missing data flag is set reasonably
      if (misdat1.eq.0.) then
        misdath=-999.99
        missing=.false.
      else
        misdath=misdat1
        missing=.true.
      endif
        
      do l=1,nt
      do i=1,nx
      do j=1,ny
c       compute range of ar2 at this grid point
        v1=ar2(i,j,  1,l)
        v2=ar2(i,j,nz2,l)

c         Find location in the grid of the valint
c         First find the level-index which is just above valint. 
          kk=0
          do k=2,nz2-1
            vv2=ar2(i,j,k  ,l)
            vv1=ar2(i,j,k-1,l)
            if ((misdat2.eq.0).or.
     &        ((vv1.ne.misdat2).and.(vv2.ne.misdat2))) then
              if ( ((vv2.ge.valint).and.(vv1.le.valint)).or. 
     &             ((vv2.le.valint).and.(vv1.ge.valint)) ) then
                  kk=k
                  v1=vv1
                  v2=vv2
c                 wird folgende zeile aktiviert, so ist Suche von unten
c                 nach oben (funktioniert komischerweise nicht)
c                 goto 200
              endif
            endif
          enddo
 200      continue

          if (kk.eq.0) then
            res(i,j,l)=misdath
            missing=.true.
          else
c           compute location in computational space of ar2
            fract=(valint-v1)/(v2-v1)
            compt(1)=real(i)
            compt(2)=real(j)
            compt(3)=(1.-fract)*(kk-1)+fract*kk
            compt(4)=real(l)

c           compute location in physical space
            call settrcom (datmin2,datmax2,stag2,4,nx,ny,nz2,nt)
            call cp2ph(phypt,compt,iflagp,4,1)

c           compute location in computational space of ar1
            call settrcom (datmin1,datmax1,stag1,4,nx,ny,nz1,nt)
            call ph2cp(phypt,compt,iflagc,4,1)

            if ((compt(3).ge.1).and.(compt(3).le.real(nz1))) then
c             get value through interpolation
              compt(1)=real(i)
              compt(2)=real(j)
              compt(4)=real(l)
              res(i,j,l)=int4dm(ar1,nx,ny,nz1,nt,
     &             compt(1),compt(2),compt(3),compt(4),misdat1)
            else
              res(i,j,l)=misdath
              missing=.true.
            endif
          endif
      enddo
      enddo
      enddo

      if (.not.missing) then
        misdat2=0.
      else
        misdat2=misdath
      endif
      end





      subroutine sum4d (a,af,nx,ny,nz,nt,i1,i2,j1,j2,misdat,ierr)
c     ===========================================================
c     this subroutine reduces the horizontal resolution of the field a by some
c     factor ired.

c     argument declarations
      integer  nx,ny,nz,nt,ierr,i1,i2,j1,j2
      real     a(nx,ny,nz*nt),af(1,1,nz*nt),misdat

c     variable declarations
      integer  k,ii,jj
      real     cnti,cnt,sum

      print *,'i1,i2,j1,j2,',i1,i2,j1,j2

      do k=1,nz*nt
        cnt=0.
        sum=0.
        do ii=i1,i2
        do jj=j1,j2
          cnti=1.
          if ((ii.eq.i1).or.(ii.eq.i2)) cnti=cnti*0.5
          if ((jj.eq.j1).or.(jj.eq.j2)) cnti=cnti*0.5
          if ((misdat.eq.0.).or.(a(ii,jj,k).ne.misdat)) then
            cnt=cnt+cnti
            sum=sum+a(ii,jj,k)
          endif
        enddo
        enddo
        if (cnt.gt.0.) then
          af(1,1,k)=sum/cnt
          print *,'k,average',k,af(1,1,k)
        else
          af(1,1,k)=misdat
        endif

      enddo

      end


      subroutine red4d (a,af,nx,ny,nz,nt,ired,misdat,ierr)
c     ====================================================
c     this subroutine reduces the horizontal resolution of the field a by some
c     factor ired.

c     argument declarations
      integer  nx,ny,nz,nt,ierr,ired
      real     a(nx,ny,nz*nt),af(nx,ny,nz*nt),misdat

c     variable declarations
      integer  i,j,k,i1,i2,ii,j1,j2,jj,m
      real     cnti,cnt,sum

      if (ired.eq.2) then
        m=1
      else if (ired.eq.4) then
        m=2
      else
        print *,'not implemented'
      endif

      do k=1,nz*nt
      do j=1,ny
      do i=1,nx
        cnt=0.
        sum=0.
        i1=i-m
        i2=i+m
        j1=j-m
        j2=j+m
        do ii=max(1,i1),min(nx,i2)
        do jj=max(1,j1),min(ny,j2)
          cnti=1.
          if ((ii.eq.i1).or.(ii.eq.i2)) cnti=cnti*0.5
          if ((jj.eq.j1).or.(jj.eq.j2)) cnti=cnti*0.5
          if ((misdat.eq.0.).or.(a(ii,jj,k).ne.misdat)) then
            cnt=cnt+cnti
            sum=sum+a(ii,jj,k)
          endif
        enddo
        enddo
        if (cnt.gt.0.) then
          af(i,j,k)=sum/cnt
        else
          af(i,j,k)=misdat
        endif

      enddo
      enddo
      enddo

      end
