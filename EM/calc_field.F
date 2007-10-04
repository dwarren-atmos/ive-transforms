      Integer function calc_field (varnam, ndims, dims, stag, datmin, 
     &                           datmax, misdat, 
     &                           data_units, data_display_units, dimnam) 
c-----------------------------------------------------------------------
c     $Id: calc_field.f,v 1.1 1994/11/14 22:36:54 warren Exp $
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
c	$Log: calc_field.f,v $
c Revision 1.1  1994/11/14  22:36:54  warren
c Christoph Schaer's European Model transforms.
c
c Revision 1.27  1994/11/11  09:02:12  schaer
c Aenderung fuer IVE-3-2-beta
c Verbessertes Heading
c Startzeit neu auch auf Datenfile
c
c Revision 1.23  1994/07/06  08:21:41  dani
c removed misdat statement for function INTER (Dani Luethi)
c
c Revision 1.22  1994/06/17  12:06:08  schaer
c Neue Funktion STDDEV und neu 4 Argumente fuer PUT.
c
c Revision 1.21  1994/05/26  11:51:59  schaer
c Neue Masken-Funktion num1<field<num2.
c Erweiterete Funktion DMEAN[field:maske].
c Aufsplitten von calc_field.f in mehrere Files.
c
c Revision 1.20  1994/05/21  19:43:43  schaer
c Neue Plot-Befehle von David Bresch.
c
c Revision 1.19  1994/05/16  06:37:47  schaer
c Neue Funktion XYSHIFT.
c
c Revision 1.18  1994/05/13  12:13:27  schaer
c Bug in REDRES behoben.
c
c Revision 1.17  1994/04/29  08:34:49  schaer
c Fehler in Berechnung von THETA behoben (war dort von 22.Maerz bis
c 29.April 94).
c
c Revision 1.16  1994/04/29  08:06:34  dani
c added functios TVARI, PERSI and RMSDIFF
c
c Revision 1.15  1994/04/12  09:59:27  schaer
c Zwei kleinere Bugs behoben.
c
c Revision 1.14  1994/04/12  08:59:58  schaer
c Neue Funktionen SETDOM und DMEAN zur Berechnung des Mittels auf
c einem Teilgebiet.
c
c Revision 1.13  1994/03/22  15:49:40  schaer
c Neue Funktionen: UG,VG,UA,VA, DIV[x_comp:y_comp], SIGMA, M,B,IPV,
c    SETSTACK, INTSTACK, LOWERCASE, UPPERCASE, INTPOL[var], NOMISSDAT[var]
c Neu Theta-Daten auch auf Druck-Flaechen.
c
c Revision 1.11  1994/02/17  14:29:16  dani
c Changes for height as vertical coordinate incorporated (Dani Luethi)
c
c Revision 1.10  1994/01/07  09:29:55  schaer
c Verbesserte Version der %T und %K Funktion. Ab sofort kann in %T direkt
c die Zeit anstelle des Zeit-Indexes angegeben werden.
c
c Revision 1.9  1993/12/24  10:25:54  schaer
c Neue Funktion PV3 (Andrea Rossa), und neue Skalierung von PV und PV3.
c
c Revision 1.8  1993/12/17  14:30:54  schaer
c Alternative Berechnung von PSRED aus Z, sowie Verallgemeinerung der
c Ableitung fuer f-plane Simulationen (Dani Luethi).
c
c Revision 1.7  1993/12/17  07:58:30  schaer
c Bug in INTER corrected.
c
c Revision 1.6  1993/12/16  08:03:12  schaer
c Bug in PUT korrigiert (Heini Wernli).
c
c Revision 1.5  1993/12/07  14:16:19  schaer
c Kleine Korrektur an der Level-Funktion <var>%K<level>
c
c Revision 1.4  1993/12/07  13:32:02  schaer
c Neue Funktion PFLT fuer Filtern auf Druckflaechen.
c
c Revision 1.2  1993/12/06  12:36:35  schaer
c Neue Funktion RTOT.
c
c Revision 1.1  1993/12/02  13:15:32  schaer
c Initial revision
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
#ifdef sun
      implicit automatic (a-z)
#endif

c     include the common-block with the level information
      include 'constants.icl'
      include 'cdfconsts.icl'
      include '/usr/local/include/netcdf.inc'

c     Argument declarations. MAXDIM is the maximum data dimension allowed,
c     LENVAR is the maximum length of the variable-name.
      integer   MAXDIM, LENVAR
      parameter (MAXDIM=4, LENVAR=80)
      character*(80)   data_units, data_display_units
      character*(LENVAR)  varnam
      integer          ndims, dims(MAXDIM)
      real             stag(MAXDIM), misdat,
     &                 datmin(MAXDIM), datmax(MAXDIM)
      character*(80) dimnam(MAXDIM)

c     Local variable declarations
      character*(LENVAR) fld,fld1,fld2,arg0,arg1,arg2,arg3,arg4,dvar
      logical          err, local, there
      logical          errflg,localu,localv
      real             misdatu,misdatv,misdvel,misdatd,misdata
      real             delta1,delta2,afil,tstrt,tstop,umin,umax
      real             timinp,rlevel,hnumb,rnumb1,rnumb2,newtim
      integer          ptru,ptrv,ptrt,ptrvel,ptrd,ptrrr,size,ptrdps
      integer          ierr,ierr1,ierr2,icheck
      integer          ptrh,i,k,ibeg,iend, ptr,ptr1,ptr2,dims1(MAXDIM)
      integer          cdfid

      character*(80)   data1_units, data1_display_units
      character*(1)    del
      logical          local1
      integer          ndims1, itime,ilevel,ired
      real             stag1(MAXDIM), misdat1,rred,valint,
     &                 datmin1(MAXDIM), datmax1(MAXDIM)
      real             datmin2(MAXDIM), datmax2(MAXDIM),
     &                 stag2(MAXDIM), misdat2
      integer          ndims2,dims2(MAXDIM)

c     variables for subdomain-selection
      real             sdomx1,sdomx2,sdomy1,sdomy2
      integer          idomx1,idomx2,idomy1,idomy2
      common /subdom/  sdomx1,sdomx2,sdomy1,sdomy2
      data             sdomx1,sdomx2,sdomy1,sdomy2/0.,0.,0.,0./
      real             compt(4),phypt(4)
      integer          iflag(4)

c     added by david:
c     variables for wind:
      integer        wsize,d
      real           c

c     External function declarations.
      real           str2nmb
      integer        getmem, addvar, getvar, strbeg, strend, gcdfvar
      integer        privat, time2index, str2int
      logical        isfunc,isfunc4

c     data-statements
      data           cdfcstfiln /'*'/
      data           cdflowcase /0/
c      data           headtyp /1/

c     initially, set pointer 'calc_field' to unused, indicating a
c     failure in computing the requested field.
      calc_field=0

c     define field-name to be computed (removes leading and trailing blanks)
      ibeg = strbeg (varnam)
      iend = strend (varnam)
      fld=varnam(ibeg:iend)
      iend = strend (fld)

      print *,'Entering calc_field for ',varnam(ibeg:iend)

c     first attempt to provide the field from the user-written
c     function privat
      calc_field=privat(fld, ndims, dims, stag, datmin, 
     &                           datmax, misdat, 
     &                           data_units, data_display_units, dimnam)
      if (calc_field.gt.0) then
         print *,'Field ',varnam(ibeg:iend),' is user-defined'
         return
      endif

c     start computation of fields

c     -----------------------------------------------------------------------
      if (fld.eq.'VEL') then
c       Compute the horizontal velocity (u^2+v^2)^0.5. In this
c       first example, this should be done by using the computational
c       package provided with IVE. Any variable which is on the
c       data-file (or which can be computed) can be assessed with
c       a call to getvar. This function provides the pointer to the
c       data as well as the attributes of the data-field. The variable
c       name (the first argument of getvar) must be defined according 
c       to the syntax for the 'field='command of IVE:
           calc_field = getvar('(U^2+V^2)^0.5', ndims, dims, stag, 
     &             datmin, datmax, misdat,
     &             data_units, data_display_units, dimnam, local)
c       All intermediate operations (e.g. reading data from the data-file,
c       doing computations) are hidden to the user. If getvar returns
c       successfully (i.e. calc_field<>0), then a new button denoted 'vel'
c       will be installed after return from calc_field.
        return

c     -----------------------------------------------------------------------
      else if (fld.eq.'THETA') then
c          computes the potential temperature

           if (datatype.lt.3) then
c            pressure or hybrid-coordinates
c            get the temperature field
             ptrt  = getvar('T', ndims, dims, stag, 
     &                    datmin, datmax, misdat,
     &                    data_units, data_display_units, dimnam, local)
             if (ptrt.eq.0) return

c            get memory
             do i=ndims+1,MAXDIM
               dims(i)=1
             enddo
             calc_field=getmem(dims(1)*dims(2)*dims(3)*dims(4))
             if (calc_field.eq.0) return

c            call routine to do the computation
             call pottemp(%val(calc_field),%val(ptrt),
     &                      dims(1),dims(2),dims(3),dims(4),stag(3))
           else
c            theta-coordinates (computation of theta works as that of
c            pressure in sigma-coordinates)

             call initatt(ndims,dims, stag, datmin, datmax, misdat,
     &         data_units, data_display_units)
             stag(3)=0.

c            get memory
             calc_field=getmem(dims(1)*dims(2)*dims(3)*dims(4))
             if (calc_field.eq.0) return

c            call routine to do the computation
             call pressure(%val(calc_field),stag(3),
     &                     dims(1),dims(2),dims(3),dims(4))
           endif

c     -----------------------------------------------------------------------
      else if (fld.eq.'THETAE') then
c          get the temperature field (PS is already available with)
           ptrt  = getvar('T', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           ptru  = getvar('QD', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if ((ptrt.eq.0).or.(ptru.eq.0)) return

c          get memory
           do i=ndims+1,MAXDIM
             dims(i)=1
           enddo
           calc_field=getmem(dims(1)*dims(2)*dims(3)*dims(4))
           if (calc_field.eq.0) return

c          call routine to do the computation
           call equpot(%val(calc_field),%val(ptrt),
     &                %val(ptru),dims(1),dims(2),dims(3),dims(4))

c     -----------------------------------------------------------------------
      else if (fld.eq.'ZLAY') then
c        computes the height Z on layers (from the height on levels)
c        get the temperature field (PS is already available with)
	   ptru  = getvar('Z', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
	   ptrt  = getvar('T', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
	   if ((ptrt.eq.0).or.(ptru.eq.0)) return

c          get memory
	   do i=ndims+1,MAXDIM
	     dims(i)=1
	   enddo
	   calc_field=getmem(dims(1)*dims(2)*dims(3)*dims(4))
	   if (calc_field.eq.0) return

c          call routine to do the computation
	   call zlayer(%val(calc_field),%val(ptrt),
     &                %val(ptru),dims(1),dims(2),dims(3),dims(4))

c     -----------------------------------------------------------------------
      else if (fld.eq.'PSRED') then
c          compute the reduced sea-level pressure

c          get the surface height
           ptru  = getvar('ZB', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
c          get the temperature field (PS is already available with)
           ptrt  = getvar('T', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if ((ptrt.eq.0).or.(ptru.eq.0)) return

c          get memory
           dims(3)=1
           do i=ndims+1,MAXDIM
             dims(i)=1
           enddo
           calc_field=getmem(dims(1)*dims(2)*dims(3)*dims(4))
           if (calc_field.eq.0) return

c          call routine to do the computation
           call pstonn(%val(calc_field),%val(ptrt),%val(ptru),
     &                        dims(1),dims(2),dims(3),dims(4))

c     -----------------------------------------------------------------------
      else if (fld.eq.'ZBCOMP') then
c          computes the surface height from the lowermost level where the
c          height Z is available

c          get the geopotential field
           ptru  = getvar('Z', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
c          get the temperature field
           ptrt  = getvar('T', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if ((ptrt.eq.0).or.(ptru.eq.0)) return

c          get memory
           dims(3)=1
           dims(4)=1
           calc_field=getmem(dims(1)*dims(2)*dims(3)*dims(4))
           if (calc_field.eq.0) return

c          call routine to do the computation
           call zboden(%val(calc_field),%val(ptrt),%val(ptru),
     &                        dims(1),dims(2),dims(3),dims(4))

c     -----------------------------------------------------------------------
      else if (fld.eq.'DIV') then
c          computes the horizontal divergence of the wind-field

           calc_field=getvar('D[U:X]+D[V*COS:Y]/COS',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if (calc_field.eq.0) return
           call t_intpole(%val(calc_field),
     &           dims(1),dims(2),dims(3),dims(4),datmin,datmax,misdat)


c     -----------------------------------------------------------------------
      else if (isfunc(fld,'DIV',arg1,arg2,arg3)) then
c          horizontal divergence of a vector field (arg1:arg2)

           if ((arg1.eq.'#').or.(arg2.eq.'#')) return
           calc_field=
     &          getvar('D['//arg1(strbeg(arg1):strend(arg1))//':X]+'      
     &          //'D['//arg2(strbeg(arg2):strend(arg2))//'*COS:Y]/COS',
     &          ndims, dims, stag, 
     &          datmin, datmax, misdat,
     &          data_units, data_display_units, dimnam, local)
           if (calc_field.eq.0) return
           call t_intpole(%val(calc_field),
     &            dims(1),dims(2),dims(3),dims(4),datmin,datmax,misdat)


c     -----------------------------------------------------------------------
      else if (fld.eq.'UG') then
c          geostrophischer Wind

           if (datatype.eq.3) then
             calc_field=getvar('-D[M:Y]*(F/((ABS(F)+1.e-12)^2))',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           else
             calc_field=
     &         getvar('-D[ZLAY:Y]*(9810*F/((ABS(F)+1.e-12)^2))',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           endif


c     -----------------------------------------------------------------------
      else if (fld.eq.'VG') then
c          geostrophischer Wind

           if (datatype.eq.3) then
             calc_field=getvar('D[M:X]*(F/((ABS(F)+1.e-12)^2))',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           else
             calc_field=
     &         getvar('D[ZLAY:X]*(9810*F/((ABS(F)+1.e-12)^2))',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           endif

c     -----------------------------------------------------------------------
      else if (fld.eq.'UA') then
c          ageostrophischer Wind

           calc_field=getvar('U-UG',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)


c     -----------------------------------------------------------------------
      else if (fld.eq.'VA') then
c          ageostrophischer Wind

           calc_field=getvar('V-VG',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)


c     -----------------------------------------------------------------------
      else if (fld.eq.'M') then
c          Montgomery-Potential

           if (datatype.ne.3) return
           calc_field=getvar('1004.*T+9.81*Z',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if (calc_field.eq.0) return
c          assign units
           data_units='J/kg'
           data_display_units='J/kg'


c     -----------------------------------------------------------------------
      else if (fld.eq.'B') then
c          Bernoulli-Function

           if (datatype.ne.3) return
           calc_field=getvar('M+0.5*(u^2+v^2)',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if (calc_field.eq.0) return
c          assign units
           data_units='J/kg'
           data_display_units='J/kg'


c     -----------------------------------------------------------------------
      else if (fld.eq.'SIGMA') then
c          computes sigma=-(dp/dtheta)/g

           if (datatype.eq.3) then
             calc_field=getvar('-D[P:T]/0.0981',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           else
             calc_field=getvar('-1./(D[THETA:P]*9.81)',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           endif


c     -----------------------------------------------------------------------
      else if (fld.eq.'IPV') then
c          computes isentropic potential vorticity

           if (datatype.ne.3) return
c          calc_field=getvar('-(VORT+F)*9.80665*1.E6*D[THETA:P]',
           calc_field=getvar('(VORT+F)*1.E6/SIGMA',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
c          assign units (PVU's)
           data_units='1.E-6 K m**2 / (kg s)'
           data_display_units='1.E-6 K m**2 / (kg s)'

c     -----------------------------------------------------------------------
      else if (fld.eq.'VORT') then
c          computes the vorticity

           calc_field=getvar('D[V:X]-D[U*COS:Y]/COS',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if (calc_field.eq.0) return
c          assign units 
           data_units='1 / s'
           data_display_units='1.E-4 / s'
c          correct for missing pole-values
           call t_intpole(%val(calc_field),
     &           dims(1),dims(2),dims(3),dims(4),datmin,datmax,misdat)
c          remove data from lowermost level
           if (datatype.eq.1) then
             if (misdat.eq.0.) misdat=1.
             call blnkbnd(%val(calc_field),misdat,1,
     &                             dims(1),dims(2),dims(3),dims(4) )
           endif


c     -----------------------------------------------------------------------
      else if (fld(1:5).eq.'VORT@') then
c          comptues the vorticity on a certain level from the wind on this
c          level. The pressure-height or the potential temperature can
c          be specified with the syntax of the @-command

           arg2=fld(6:strend(fld))
           calc_field=getvar(
     &       'D[V@'//arg2(1:strend(arg2))//':XM]'//
     &       '-D[U@'//arg2(1:strend(arg2))//'*COS:YM]/COS',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if (calc_field.eq.0) return
c          assign units 
           data_units='1 / s'
           data_display_units='1.E-4 / s'


c     -----------------------------------------------------------------------
      else if (fld.eq.'F') then
c          computes the coriolis parameter

           calc_field=getvar('0.000145444*SIN[PHI]', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if (calc_field.eq.0) return
c          assign units 
           data_units='1 / s'
           data_display_units='1.E-4 / s'
           

c     -----------------------------------------------------------------------
      else if (fld.eq.'PV') then
c          potentiel vorticity

           calc_field=getvar('-(VORT+F)*9.80665*1.E6*D[THETA:P]',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if (calc_field.eq.0) return
c          assign units (PVU's)
           data_units='1.E-6 K m**2 / (kg s)'
           data_display_units='1.E-6 K m**2 / (kg s)'
c          remove data from lowermost level
           if (datatype.eq.1) then
             if (misdat.eq.0.) misdat=1.
             call blnkbnd(%val(calc_field),misdat,1,
     &                             dims(1),dims(2),dims(3),dims(4) )
           endif
           

c     -----------------------------------------------------------------------
      else if (fld.eq.'PV3') then
c          potential vorticity with horizontal vorticity-comp. included

           calc_field=getvar('PV-9.80665*1.E6*'//
     &     '{D[U:P]*D[THETA:Y]-D[V:P]*D[THETA:X]}',
     &                  ndims, dims, stag,
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if (calc_field.eq.0) return
c          assign units (PVU's)
           data_units='1.E-6 K m**2 / (kg s)'
           data_display_units='1.E-6 K m**2 / (kg s)'

c     -----------------------------------------------------------------------
      else if (fld.eq.'RH') then
c          computes the relative humidity

c          get the temperature field (PS is already available with)
           ptrt  = getvar('T', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           ptru  = getvar('QD', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           ptrv  = getvar('QW', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if ((ptrt.eq.0).or.(ptru.eq.0).or.(ptrv.eq.0)) return

c          get memory
           do i=ndims+1,MAXDIM
             dims(i)=1
           enddo
           calc_field=getmem(dims(1)*dims(2)*dims(3)*dims(4))
           if (calc_field.eq.0) return

c          call routine to do the computation
           call relhum(%val(calc_field),%val(ptrt),
     &       %val(ptru),%val(ptrv),dims(1),dims(2),dims(3),dims(4))

c     -----------------------------------------------------------------------
      else if (fld.eq.'P') then
c          computes the pressure on the layers

           call initatt(ndims,dims, stag, datmin, datmax, misdat,
     &       data_units, data_display_units)
           stag(3)=-0.5

c          get memory
           calc_field=getmem(dims(1)*dims(2)*dims(3)*dims(4))
           if (calc_field.eq.0) return

c          call routine to do the computation
           call pressure(%val(calc_field),stag(3),
     &                     dims(1),dims(2),dims(3),dims(4))

c     -----------------------------------------------------------------------
      else if (fld.eq.'PLEV') then
c          computes the pressure on the levels

           call initatt(ndims,dims, stag, datmin, datmax, misdat,
     &       data_units, data_display_units)
           stag(3)=0.

c          get memory
           calc_field=getmem(dims(1)*dims(2)*dims(3)*dims(4))
           if (calc_field.eq.0) return

c          call routine to do the computation
           call pressure(%val(calc_field),stag(3),
     &                     dims(1),dims(2),dims(3),dims(4))


c     -----------------------------------------------------------------------
      else if (fld.eq.'PSX') then
c          computes the horizontal derivative of PS
           calc_field=getvar('D[PS:XM]', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)

c     -----------------------------------------------------------------------
      else if (fld.eq.'PSY') then
c          computes the horizontal derivative of PS
           calc_field=getvar('D[PS:YM]', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)


c     -----------------------------------------------------------------------
      else if (isfunc(fld,'GET',arg1,arg2,arg3)) then
c       reads in the field ARG2 from the file ARG1
        if ((arg3.eq.'#').and.(arg2.ne.'#')) then
c         convert filename to lowercase
          if (cdflowcase.eq.1) 
     &         call locase(arg1(1:strend(arg1)),strend(arg1))
          calc_field=gcdfvar(arg1, arg2, ndims, dims, stag, 
     &                           datmin,datmax, misdat, 
     &                           data_units, data_display_units) 
        endif


c     -----------------------------------------------------------------------
      else if (isfunc4(fld,'PUT',arg1,arg2,arg3,arg4)) then
c      else if (isfunc(fld,'PUT',arg1,arg2,arg3)) then
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
c        convert filename to lowercase
         if (cdflowcase.eq.1)
     &          call locase(arg1(1:strend(arg1)),strend(arg1))

c        the integer function 'getvar' is a function provided to the
c        user by IVE. It gets the data-pointer to the selected variable
c        as well as its attributes. Here we have to get the variable
c        defined by 'arg2' (string).
         ptr = getvar(arg2, ndims, dims, stag,
     &                datmin, datmax, misdat,
     &                data_units, data_display_units, dimnam, local)

c        if getvar returns with the pointer-value 0, then the variable was
c        not found. Return from subroutine
         if (ptr.eq.0) return

c        the variable 'arg2' was found and can be used
c	 if arg3 is not defined, use the same name to save the variable
         if (arg3.eq.'#') arg3=arg2

c        if arg4 is not defined, keep time-location as is now. Otherwise
c        the field to be written to the file is assigned the new time arg4
         if (arg4.eq.'#') then
           newtim=-99.99
         else if ((dims(4).eq.1).or.(ndims.lt.4)) then
           newtim=str2nmb(arg4,ierr)
           if (ierr.ne.0) then
             print *,'Illegal 4th argument in ',fld(1:iend)
             return
           endif
         else
           print *,'Argument 4 in ',fld(1:iend),' only with dims(4)=1'
           return
         endif

c        make sure that a new file is opened with useful global
c        attributes. This is done here in order not to make pcdfvar
c        dependent upon the transforms.
         if (cdfcstfiln.ne.'*') then
           inquire (file = arg1(1:strend(arg1)), exist=there )
           if (.not.there) then
c            the file does not yet exist. Compute datmin1 and datmin2.
             do i=1,3
               datmin1(i)=datmin(i)
               datmax1(i)=datmax(i)
             enddo
             if (cdfdatatype.eq.2) then
               datmin1(3)=amax1(datmin(3),1000.)
               datmax1(3)=0.
             else if (cdfdatatype.eq.3) then
               datmin1(3)=0.
               datmax1(3)=amax1(datmax(3),1000.)
             endif  
c            create the file
             call ncpopt(NCVERBOS)
             call crecdf(arg1(1:strend(arg1)),cdfid,datmin1,datmax1,
     &             3,cdfcstfiln(1:strend(cdfcstfiln)),ierr)
             call clscdf(cdfid,ierr)
           endif
         endif

         call pcdfvar (arg1, arg3, ptr, newtim, ndims, dims,
     &                          stag, datmin, datmax, misdat, 
     &                          data_units, data_display_units) 

         return

c     -----------------------------------------------------------------------
      else if (isfunc(fld,'SETSTACK',arg1,arg2,arg3)) then
c        this command allows the definition of a constants-filename
c        later to be used with the commands PUT und INTSTACK
         if ((arg1.eq.'#').or.(arg2.ne.'#')) then
           print *,'Illegal number of arguments in ',fld(1:iend)
           return
         endif

c        convert constants-filename to lowercase
         if (cdflowcase.eq.1) 
     &          call locase(arg1(1:strend(arg1)),strend(arg1))
c        read constants-filename
         call setstack (arg1(1:strend(arg1)))
         return

c     -----------------------------------------------------------------------
      else if (isfunc(fld,'INTSTACK',arg1,arg2,arg3)) then
c        interpolates the field ARG1 to a stack of pressure of theta
c        levels, which are specified through a 
c        constants-file previously selected with a call to SETSTACK.
         if ((arg1.eq.'#').or.(arg2.ne.'#')) then
           print *,'Illegal number of arguments in ',fld(1:iend)
           return
         endif

c        Check whether the constants-filename is defined
         if (cdfcstfiln.eq.'*') then
           print *,'You have to specify the constants-filename with ',
     &             'SETSTACK first!'
           return
         endif
         if ((cdfdatatype.ne.3).and.(cdfdatatype.ne.2)) then
           print *,'INTSTACK is not implemented for the vertical ',
     &       'coordinates on ',cdfcstfiln(1:strend(cdfcstfiln)),'.'
           return 
         endif

c        get the field
         ptr = getvar(arg1(strbeg(arg1):strend(arg1)), 
     &                ndims, dims, stag, datmin, datmax, misdat,
     &                data_units, data_display_units, dimnam, local)
         if (ptr.eq.0) then
           print *,'Unable to get ',arg1(1:strend(arg1)),
     &             ' in ',fld(1:iend)
           return
         endif
         if ((ndims.lt.3).or.(dims(3).eq.1)) then
           print *,'The field ',arg1(1:strend(arg1)),
     &             ' is one-dimensional in the vertical.'
           return
         endif

c        redefine vertical coordinate
         dims(3)=cdfnz
         stag(3)=0.
         datmin(3)=cdfaklev(1)
         datmax(3)=cdfaklev(cdfnz)

c        compute the size of the new field and get memory
         size=1
         do i=1,ndims
           size=size*dims(i)
         enddo
         calc_field=getmem(size)
         if (calc_field.eq.0) then
           print *,'Unable to get Memory in INTSTACK'
           return
         endif

         call intstack (arg1(1:strend(arg1)),%val(calc_field),
     &                 dims(1),dims(2),dims(3),dims(4),misdat,ierr)
         if (ierr.ne.0) then
           call freemem(calc_field)
           calc_field=0
         endif
         return


c     -----------------------------------------------------------------------
      else if ((fld.eq.'GZ0').or.(fld.eq.'ZB').or.(fld.eq.'BLA')
     &     .or.(fld.eq.'BTY')) then
c       reads in external data
        calc_field=gcdfvar(extfiln, fld, ndims, dims, stag, 
     &                           datmin,datmax, misdat, 
     &                           data_units, data_display_units) 
        if((calc_field.eq.0).and.(fld.eq.'ZB')) then
          calc_field = getvar('ZBCOMP', ndims, dims, stag,
     &                datmin, datmax, misdat,
     &                data_units, data_display_units, dimnam, local)
        endif
        return

c     -----------------------------------------------------------------------
      else if (fld.eq.'LOWERCASE') then
c       definiere uppercase/lowercase flag
        cdflowcase=1
        print *
        print *,'Switching to LOWERCASE filenames.'
        print *
        return

c     -----------------------------------------------------------------------
      else if (fld.eq.'UPPERCASE') then
c       definiere uppercase/lowercase flag
        cdflowcase=0
        print *
        print *,'Switching to UPPERCASE filenames.'
        print *
        return

c     -----------------------------------------------------------------------
      else if (isfunc(fld,'HEADING',arg1,arg2,arg3)) then
c     set the appropriate heading
        calc_field=0
        if ((arg1.eq.'#').or.(arg2.ne.'#')) then
          print *,'FIELD=HEADING[] requires one argument'
          return
        endif
        call t_heading(arg1,ierr)
        if (ierr.gt.0) then
          print *,'Illegal argument in FIELD=',fld(1:iend)
        endif
        return

c     -----------------------------------------------------------------------
      else if (fld.eq.'DERPTYPE1') then
c       definiere Typ der p-Ableitung (siehe constants.icl)
        derptype=1
        print *
        print *,'Switching to 3-point vertical dervative'
        print *
        return

c     -----------------------------------------------------------------------
      else if (fld.eq.'DERPTYPE2') then
c       definiere Typ der p-Ableitung (siehe constants.icl)
        derptype=2
        print *
        print *,'Switching to 2-point centered vertical dervative'
        print *
        return

c     -----------------------------------------------------------------------
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
     &                  data_units, data_display_units, dimnam, local)
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
             calc_field=ptrd
           else if ((datatype.ne.1).or.(dims(3).eq.1).or.
     &       (arg2.eq.'XM').or.(arg2.eq.'YM')) then
c            -----------------------------------
c            derivation on data (model) surfaces
c            -----------------------------------
             ptru   = getvar('COS', ndims1, dims1, stag1, 
     &              datmin1, datmax1, misdat1,
     &              data1_units, data1_display_units, dimnam, local1)
             call der2dm(%val(ptr),%val(ptrd),%val(ptru),
     &           arg2,dims(1),dims(2),dims(3)*dims(4),
     &           datmin,datmax,misdata,misdat)
             calc_field=ptrd
             if (local1) call freemem(ptru)
           else
c            ---------------------------------------
c            derivation on p-sufaces (data on sigma)
c            ---------------------------------------
c            first get the respective derivative of the surface pressure
c            print *,'getting PSX or PSY
             dvar='PS'//arg2(1:1)
             ptru   = getvar('COS', ndims1, dims1, stag1, 
     &              datmin1, datmax1, misdat1,
     &              data1_units, data1_display_units, dimnam, local1)
             ptrdps = getvar(dvar, ndims1, dims1, stag1, 
     &              datmin1, datmax1, misdat1,
     &              data1_units, data1_display_units, dimnam, local1)
             if ((ptrdps.eq.0).or.(ptru.eq.0)) then
               call freemem(ptrd)
               calc_field=0
               return
             endif
             if ((dims1(1).ne.dims(1)).or.(dims1(2).ne.dims(2)).or.
     &         (dims1(4).ne.dims(4))) then
               print *,'The fields ',arg1(strbeg(arg1):strend(arg1)),
     &           ' and PS are incompatible.'
               call freemem(ptrd)
               if (local1) call freevar(dvar)
               calc_field=0
               return
             endif

c            Call the subroutine to compute the derivative.
             call der3d(%val(ptr),%val(ptrd),%val(ptru),
     &              arg2,dims(1),dims(2),dims(3),dims(4),
     &              %val(ptrps),%val(ptrdps),datmin,datmax,stag)
             calc_field=ptrd
           endif
c       return unused and temporary memory, return to calling subroutine
  961      continue
           if ((ptr.gt.0).and.(local)) call freevar(arg1)
           return

c     -----------------------------------------------------------------------
      else if (fld.eq.'RTOT') then
c        total precipitation sum
         calc_field=getvar('RRSN+RSSN+RRKN+RSKN', ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
         return


c     -----------------------------------------------------------------------
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
     &                  data_units, data_display_units, dimnam, local)
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
           calc_field=ptrd
         else
           call freemem(ptrd)
         endif
c        return unused and temporary memory, return to calling subroutine
 963     continue
         if ((ptrrr.gt.0).and.(local)) call freevar(arg1)
         return


c     -----------------------------------------------------------------------
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

c       get the contour of arg2 to be interpolated to
        valint=str2nmb(arg3,ierr) 
        if (ierr.eq.1) then
          print *,arg3(strbeg(arg3):strend(arg3)),' is not a number'
          return
        endif

c       get the field arg1
        ptr1 = getvar(arg1, ndims, dims, stag, 
     &                datmin, datmax, misdat,
     &                data_units, data_display_units, dimnam, local)
        if (ptr1.eq.0) return

c       get the field arg2
        if (arg2.ne.'P') then
          ptr2 = getvar(arg2, ndims2, dims2, stag2, 
     &                datmin2, datmax2, misdat2,
     &                data_units, data_display_units, dimnam, local)
          if (ptr2.eq.0) return

c         check whether the fields are compatible
          do i = ndims2+1,MAXDIM
            dims2(i)=1
          enddo
          do i = ndims+1,MAXDIM
            dims(i)=1
          enddo
          if ( (dims(1).ne.dims2(1)).or.(dims(1).ne.dims2(1)).or.
     &       (dims(4).ne.dims2(4)).or.
     &       (abs(datmin(1)-datmin2(1)).gt.1.e-2).or.
     &       (abs(datmin(2)-datmin2(2)).gt.1.e-2).or.
     &       (abs(datmin(3)-datmin2(3)).gt.1.e-2).or.
     &       (abs(datmax(1)-datmax2(1)).gt.1.e-2).or.
     &       (abs(datmax(2)-datmax2(2)).gt.1.e-2).or.
     &       (abs(datmax(3)-datmax2(3)).gt.1.e-2).or.
     &       (abs(stag(1)-stag2(1)).gt.1.e-2).or.
     &       (abs(stag(2)-stag2(2)).gt.1.e-2).or.
     &       (abs(stag(3)-stag2(3)).gt.1.e-2) ) then
             print *,'The fields in ',fld(1:iend),' have incompatible',
     &               ' attributes.'
             return
          endif
        endif

c       allocate memory
        size=dims(1)*dims(2)*dims(4)
        calc_field  = getmem (size)
        if (calc_field.eq.0) then
          print *,'Could not allocate memory in calc_field'
          return
        endif

        if (arg2.eq.'P') then
c         int2prs (nam,ar,prsval,res,nx,ny,nz,nt,datmin,datmax,stag,misdat)
          call int2prs (arg1(1:strend(arg1)),
     &                 %val(ptr1),valint,%val(calc_field),
     &                 dims(1),dims(2),dims(3),dims(4),
     &                 datmin ,datmax ,stag ,misdat)
        else
          call int2iso (arg1(1:strend(arg1)),arg2(1:strend(arg2)),
     &                 %val(ptr1),%val(ptr2),valint,%val(calc_field),
     &                 dims(1),dims(2),dims(3),dims2(3),dims(4),
     &                 datmin ,datmax ,stag ,misdat ,
     &                 datmin2,datmax2,stag2,misdat2)
        endif
        dims(3)=1
        return
        return

c     -----------------------------------------------------------------------
      else if (isfunc(fld,'XYSHIFT',arg1,arg2,arg3)) then
c        shift all the fields by (arg1,arg2) grid-points into the positive x
c        and y-direction.
         if ((arg3.ne.'#').or.(arg2.eq.'#')) then
           print *,'Illegal number of arguments in ',fld(1:strend(fld))
           return
         endif
         hnumb=str2nmb(arg1,ierr1)
         hnumb=str2nmb(arg2,ierr2)
         if ((ierr1.ne.0).or.(ierr2.ne.0)) then
           print *,'Illegal arguments in ',fld(1:strend(fld))
           return
         endif
         xyshift(1)=str2nmb(arg1,ierr)
         xyshift(2)=str2nmb(arg2,ierr)
         print *
         print *,'Fields will be shifted until new file is loaded.'
         print *,'Reselect current field before proceeding!'
         print *

c     -----------------------------------------------------------------------
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

c     -----------------------------------------------------------------------
      else if ((fld.eq.'GRIDX').or.(fld.eq.'GRIDY').or.
     &         (fld.eq.'RLAT').or.(fld.eq.'RLON').or.
     &         (fld.eq.'RLAM').or.(fld.eq.'RPHI').or.
     &         (fld.eq.'RPHIM').or.
     &         (fld.eq.'LAT').or.(fld.eq.'LON').or.
     &         (fld.eq.'LAM').or.(fld.eq.'PHI').or.
     &         (fld.eq.'COS').or.(fld.eq.'COSR').or.
     &         (fld.eq.'TAN').or.(fld.eq.'SIN')) then
         
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
         calc_field  = getmem (size)
         if (calc_field.eq.0) then
           print *,'Could not allocate memory in calc_field'
           return
         endif

         call cregrid (fld,%val(calc_field),dims(1),dims(2),
     &                datmin,datmax,misdat)




c     -----------------------------------------------------------------------
      else if (isfunc(fld,'INTPOL',arg1,arg2,arg3)) then
c       function to interpolate the pole values
        if ((arg2.eq.'#').and.(arg1.ne.'#')) then
          calc_field=getvar(arg1(strbeg(arg1):strend(arg1)),
     &              ndims, dims, stag, 
     &              datmin, datmax, misdat,
     &              data_units, data_display_units, dimnam, local)
          if (calc_field.eq.0) return
          if (.not.local) then
c           this is for a 'permanent' computed field, get memory
            size=1
            do i=1,ndims
              size=size*dims(i)
            enddo
            ptr=getmem(size)
            if (ptr.eq.0) return
            call copyar(%val(calc_field),%val(ptr),size)
            calc_field=ptr
          endif
c         interpolate pole-values
          call t_intpole(%val(calc_field),
     &       dims(1),dims(2),dims(3),dims(4),datmin,datmax,misdat)
        endif

c     -----------------------------------------------------------------------
      else if (isfunc(fld,'NOMISDAT',arg1,arg2,arg3)) then
c       interpolates and extrapolates all missing data in the vertical. No 
c       action is taken at those gridpoints (i,j,k) which have missing
c       data for all k.
        if ((arg1.eq.'#').or.(arg2.ne.'#')) return
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, dimnam, local)
        print *,'nomisdat',arg1(strbeg(arg1):strend(arg2)),ptr,misdat
        if ((ptr.eq.0).or.(dims(3).eq.1)) return

c       get memory
        size=dims(1)*dims(2)*dims(3)*dims(4)
        calc_field=getmem(size)
        if (calc_field.eq.0) then
          print *,'Unable to get memory in NOMISDAT'
          return
        endif

        call copyar(%val(ptr),%val(calc_field),size)
        if (arg1(strbeg(arg1):strend(arg1)).eq.'P') then
c            louzy bug fix
             misdat=-999.99
        endif
        print *,'nomisdat',arg1(strbeg(arg1):strend(arg2)),ptr,misdat

        if (misdat.ne.0.) call remmisdat(%val(calc_field),
     &                 dims(1),dims(2),dims(3),dims(4),misdat)
        misdat=0.


c     -----------------------------------------------------------------------
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
     &               data_units, data_display_units, dimnam, local)
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
        calc_field=ptrd
        misdat=misdatd

c     -----------------------------------------------------------------------
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
     &               data_units, data_display_units, dimnam, local)
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
        calc_field=ptrd


c     -----------------------------------------------------------------------
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
     &               data_units, data_display_units, dimnam, local)
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
        calc_field=ptrd


c     -----------------------------------------------------------------------
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
     &               data_units, data_display_units, dimnam, local)
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
        calc_field  = getmem (size)
        if (calc_field.eq.0) then
          print *,'Could not allocate memory in calc_field'
          return
        endif

        call minmax(%val(ptr),%val(calc_field),umin,umax,
     &             dims(1),dims(2),dims(3),dims(4))
        return

c     -----------------------------------------------------------------------
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

c       Next get the field
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, dimnam, local)
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
        calc_field=ptrd


c     -----------------------------------------------------------------------
      else if (isfunc(fld,'PFLT',arg1,arg2,arg3)) then
c       smooth the field <arg1> on the pressure-level <arg2> with the
c       filter-value <arg3>. Default for <arg3> is 1.
c       check for consistent number of arguments
        if ((arg1.eq.'#').or.(arg2.eq.'#')) then
          print *,' Inconsistent number of arguments in ',
     &              fld(strbeg(fld):strend(fld)),'.'
          return
        endif     
        if (arg3.eq.'#') arg3='1'

c       construct command for filter and interpolation operation
        fld1='FLT{INTER['//arg1(strbeg(arg1):strend(arg1))//
     &            ':P:'//arg2(strbeg(arg2):strend(arg2))//']:'//
     &            arg3(strbeg(arg3):strend(arg3))//'}'
        calc_field=getvar(fld1(strbeg(fld1):strend(fld1)),
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
        return


c     -----------------------------------------------------------------------
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
     &               data_units, data_display_units, dimnam, local)
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
        call filt4d (%val(ptr),%val(ptrd),dims(1),dims(2),
     &           dims(3),dims(4),datmin,datmax,afil,misdat,ierr)
        if (ierr.eq.1) then
          call freemem(ptrd)
          return
        else
          calc_field=ptrd
          return
        endif


c     -----------------------------------------------------------------------
      else if (isfunc(fld,'REDRES',arg1,arg2,arg3)) then
c       reduce the resolution of field ARG1 by an integer-factor given by ARG2
        
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
          if ((ired/2)*2-ired.ne.0) then
            print *,'Reducing-Factor must be a multiple of 2!'
            return
          else
            print *,'Reduction of resolution by factor ',ired
          endif
        endif

c       get variable
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, dimnam, local)
        if (ptr.eq.0) return
        do i=ndims+1,MAXDIM
          dims(i)=1
        enddo

c       get memory
        ptrd=getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (ptrd.eq.0) return

c       call the routine to do the computation
        call red4d (%val(ptr),%val(ptrd),
     &           dims(1),dims(2),dims(3),dims(4),ired,misdat)
        calc_field=ptrd


c     -----------------------------------------------------------------------
      else if (isfunc4(fld,'SETDOM',arg1,arg2,arg3,arg4)) then
c       define the subdomain for function AVRG
        if (arg4.eq.'#') then
          print *,'Illegal number of arguments in ',
     &                            fld(strbeg(fld):strend(fld))
          return
        endif
        sdomx1=str2nmb(arg1,ierr) 
        if (ierr.ne.0) goto 871
        sdomx2=str2nmb(arg2,ierr) 
        if (ierr.ne.0) goto 871
        sdomy1=str2nmb(arg3,ierr) 
        if (ierr.ne.0) goto 871
        sdomy2=str2nmb(arg4,ierr) 
        if (ierr.ne.0) goto 871
        print *,'Subdomain set to '
        print 872,'    ',sdomx1,' < x <',sdomx2
        print 872,'    ',sdomy1,' < y <',sdomy2
 872    format (a4,f8.3,a6,f8.3)
        return

c       error-exit
 871    continue
        print *,'Illegal arguments in ',fld(strbeg(fld):strend(fld))
        sdomx1=0.
        sdomx2=0.
        sdomy1=0.
        sdomy2=0.
        return

c     -----------------------------------------------------------------------
      else if (isfunc(fld,'DMEAN',arg1,arg2,arg3)) then
c       take the average value over a subdomain previously specified with
c       field=setdom[x1,x2,y1,y2]
        
c       check for consistent number of arguments
        if ((arg3.ne.'#').or.(arg1.eq.'#')) then
          print *,' Inconsistent number of arguments in ',
     &              fld(strbeg(fld):strend(fld)),'.'
          return
        endif

c       get variable
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, dimnam, local)
        if (ptr.eq.0) return
        do i=ndims+1,MAXDIM
          dims(i)=1
        enddo

c       check whether domain is defined
        if ((sdomx1.ge.sdomx2).or.(sdomy1.ge.sdomy2)) then
          print *,'Subdomain must first be defined with SETDOM.'
          return
        endif
c       check whether domain is a subdomain
        if ((sdomx1.lt.datmin(1)).or.(sdomx2.gt.datmax(1)).or.
     &      (sdomy1.lt.datmin(2)).or.(sdomy2.gt.datmax(2))) then
          print *,'Subdomain not inside full domain!'
          return
        endif

c       convert domain to grid-coordinates
        iflag(1)=1
        iflag(2)=1
        iflag(3)=0
        iflag(4)=0
        call new_field_trans(arg1,%val(ptr),dims(1),dims(2),1,1)
        phypt(1)=sdomx1
        phypt(2)=sdomy1
        call phys_2_index_trans(phypt,compt,iflag,4,1)
        idomx1=nint(compt(1))
        idomy1=nint(compt(2))
c       check whether the domain is integer-based
        if ((abs(real(idomx1)-compt(1)).gt.0.1).or.
     &      (abs(real(idomy1)-compt(2)).gt.0.1)) then
          print *,'Subdomain does not coincide with grid-structure'
          return
        endif
        phypt(1)=sdomx2
        phypt(2)=sdomy2
        call phys_2_index_trans(phypt,compt,iflag,4,1)
        idomx2=nint(compt(1))
        idomy2=nint(compt(2))
c       check whether the domain is integer-based
        if ((abs(real(idomx2)-compt(1)).gt.0.1).or.
     &      (abs(real(idomy2)-compt(2)).gt.0.1)) then
          print *,'Subdomain does not coincide with grid-structure'
          return
        endif

c       check for second argument
        icheck=0
        if (arg2.ne.'#') then
          ptrd = getvar(arg2, ndims1, dims1, stag1, 
     &            datmin1, datmax1, misdat1,
     &            data1_units, data1_display_units, dimnam, localu)
          if (ptrd.eq.0) return
          if ((dims1(1).ne.dims(1)).or.(dims1(2).ne.dims(2)).or.
     &        (datmin1(1).ne.datmin(1)).or.(datmin1(2).ne.datmin(2)).or.
     &        (datmax1(1).ne.datmax(1)).or.(datmax1(2).ne.datmax(2)).or.
     &        ((dims1(3).gt.1).and.(dims1(3).ne.dims(3))).or.
     &        ((dims1(4).gt.1).and.(dims1(4).ne.dims(4)))) then
              print *,'The fields ',arg1(1:strend(arg1)),' and ',
     &            arg2(1:strend(arg2)),' have incompatible attributes'
              return
          endif
          if (misdat1.ne.0) icheck=1
        else
          ptrd=ptr
          do i=1,4
            dims1(i)=dims(i)
          enddo
          misdat1=misdat
        endif
        print *,'icheck',icheck,dims1(1),dims1(2),misdat1

c       print out integer subdomain
        print *,'Computing ',fld(strbeg(fld):strend(fld)),' in'
        print 873,'    ',idomx1,' < ix <',idomx2
        print 873,'    ',idomy1,' < iy <',idomy2
 873    format (a4,i4,a7,i4)

c       get memory
        calc_field=getmem(dims(3)*dims(4))
        if (calc_field.eq.0) return

c       call the routine to do the computation
        call dmean (%val(ptr),%val(calc_field),
     &           dims(1),dims(2),dims(3),dims(4),
     &           idomx1,idomx2,idomy1,idomy2,
     &           datmin,datmax,misdat,
     &           icheck,
     &           %val(ptrd),dims1(1),dims1(2),dims1(3),dims1(4),misdat1)
        dims(1)=1
        dims(2)=1
        if ((ndims.lt.3).or.
     &      ((ndims.eq.3).and.(dims(3).eq.1)).or.
     &      ((ndims.eq.4).and.(dims(3).eq.1).and.(dims(4).eq.1))) 
     &      ndims=1
          

c     -----------------------------------------------------------------------
      else if (isfunc(fld,'TVARI',arg1,arg2,arg3)) then
c       take the average value over a subdomain previously specified with
c       field=setdom[x1,x2,y1,y2]
        
c       check for consistent number of arguments
        if ((arg2.ne.'#').or.(arg1.eq.'#')) then
          print *,' Inconsistent number of arguments in ',
     &              fld(strbeg(fld):strend(fld)),'.'
          return
        endif

c       get variable
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, dimnam, local)
        if (ptr.eq.0) return
        do i=ndims+1,MAXDIM
          dims(i)=1
        enddo

c       check whether domain is defined
        if ((sdomx1.ge.sdomx2).or.(sdomy1.ge.sdomy2)) then
          print *,'Subdomain must first be defined with SETDOM.'
          return
        endif
c       check whether domain is a subdomain
        if ((sdomx1.lt.datmin(1)).or.(sdomx2.gt.datmax(1)).or.
     &      (sdomy1.lt.datmin(2)).or.(sdomy2.gt.datmax(2))) then
          print *,'Subdomain not inside full domain!'
          return
        endif

c       convert domain to grid-coordinates
        iflag(1)=1
        iflag(2)=1
        iflag(3)=0
        iflag(4)=0
        call new_field_trans(arg1,%val(ptr),dims(1),dims(2),1,1)
        phypt(1)=sdomx1
        phypt(2)=sdomy1
        call phys_2_index_trans(phypt,compt,iflag,4,1)
        idomx1=nint(compt(1))
        idomy1=nint(compt(2))
c       check whether the domain is integer-based
        if ((abs(real(idomx1)-compt(1)).gt.0.1).or.
     &      (abs(real(idomy1)-compt(2)).gt.0.1)) then
          print *,'Subdomain does not coincide with grid-structure'
          return
        endif
        phypt(1)=sdomx2
        phypt(2)=sdomy2
        call phys_2_index_trans(phypt,compt,iflag,4,1)
        idomx2=nint(compt(1))
        idomy2=nint(compt(2))
c       check whether the domain is integer-based
        if ((abs(real(idomx2)-compt(1)).gt.0.1).or.
     &      (abs(real(idomy2)-compt(2)).gt.0.1)) then
          print *,'Subdomain does not coincide with grid-structure'
          return
        endif
c       print out integer subdomain
        print *,'Computing time variability for ',
     &        arg1(strbeg(arg1):strend(arg1)),' in'
        print 873,'    ',idomx1,' < ix <',idomx2
        print 873,'    ',idomy1,' < iy <',idomy2

c       get memory
        ptrd=getmem(dims(3)*dims(4))
        if (ptrd.eq.0) return

c       call the routine to do the computation
        call varblty (%val(ptr),%val(ptrd),
     &           dims(1),dims(2),dims(3),dims(4),
     &           idomx1,idomx2,idomy1,idomy2,misdat)
        calc_field=ptrd
        dims(1)=1
        dims(2)=1
        if ((ndims.lt.3).or.
     &      ((ndims.eq.3).and.(dims(3).eq.1)).or.
     &      ((ndims.eq.4).and.(dims(3).eq.1).and.(dims(4).eq.1))) 
     &      ndims=1
          

c     -----------------------------------------------------------------------
      else if (isfunc(fld,'PERSI',arg1,arg2,arg3)) then
c       take the average value over a subdomain previously specified with
c       field=setdom[x1,x2,y1,y2]
        
c       check for consistent number of arguments
        if ((arg2.ne.'#').or.(arg1.eq.'#')) then
          print *,' Inconsistent number of arguments in ',
     &              fld(strbeg(fld):strend(fld)),'.'
          return
        endif

c       get variable
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, dimnam, local)
        if (ptr.eq.0) return
        do i=ndims+1,MAXDIM
          dims(i)=1
        enddo

c       check whether domain is defined
        if ((sdomx1.ge.sdomx2).or.(sdomy1.ge.sdomy2)) then
          print *,'Subdomain must first be defined with SETDOM.'
          return
        endif
c       check whether domain is a subdomain
        if ((sdomx1.lt.datmin(1)).or.(sdomx2.gt.datmax(1)).or.
     &      (sdomy1.lt.datmin(2)).or.(sdomy2.gt.datmax(2))) then
          print *,'Subdomain not inside full domain!'
          return
        endif

c       convert domain to grid-coordinates
        iflag(1)=1
        iflag(2)=1
        iflag(3)=0
        iflag(4)=0
        call new_field_trans(arg1,%val(ptr),dims(1),dims(2),1,1)
        phypt(1)=sdomx1
        phypt(2)=sdomy1
        call phys_2_index_trans(phypt,compt,iflag,4,1)
        idomx1=nint(compt(1))
        idomy1=nint(compt(2))
c       check whether the domain is integer-based
        if ((abs(real(idomx1)-compt(1)).gt.0.1).or.
     &      (abs(real(idomy1)-compt(2)).gt.0.1)) then
          print *,'Subdomain does not coincide with grid-structure'
          return
        endif
        phypt(1)=sdomx2
        phypt(2)=sdomy2
        call phys_2_index_trans(phypt,compt,iflag,4,1)
        idomx2=nint(compt(1))
        idomy2=nint(compt(2))
c       check whether the domain is integer-based
        if ((abs(real(idomx2)-compt(1)).gt.0.1).or.
     &      (abs(real(idomy2)-compt(2)).gt.0.1)) then
          print *,'Subdomain does not coincide with grid-structure'
          return
        endif
c       print out integer subdomain
        print *,'Computing persistence for ',
     &        arg1(strbeg(arg1):strend(arg1)),' in'
        print 873,'    ',idomx1,' < ix <',idomx2
        print 873,'    ',idomy1,' < iy <',idomy2

c       get memory
        ptrd=getmem(dims(3)*dims(4))
        if (ptrd.eq.0) return

c       call the routine to do the computation
        call persdiff (%val(ptr),%val(ptrd),
     &           dims(1),dims(2),dims(3),dims(4),
     &           idomx1,idomx2,idomy1,idomy2,misdat)
        calc_field=ptrd
        dims(1)=1
        dims(2)=1
        if ((ndims.lt.3).or.
     &      ((ndims.eq.3).and.(dims(3).eq.1)).or.
     &      ((ndims.eq.4).and.(dims(3).eq.1).and.(dims(4).eq.1))) 
     &      ndims=1
          
          

c     -----------------------------------------------------------------------
      else if (isfunc(fld,'STDDEV',arg1,arg2,arg3)) then
c       take the standarddeviation of a field
c       check for consistent number of arguments
        if ((arg2.ne.'#').or.(arg1.eq.'#')) then
          print *,' Inconsistent number of arguments in ',
     &              fld(strbeg(fld):strend(fld)),'.'
          return
        endif

        calc_field=getvar('TMEAN{('//arg1(1:strend(arg1))//'-'//
     &                     'TMEAN['//arg1(1:strend(arg1))//'])^2}^0.5',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)


c     -----------------------------------------------------------------------
      else if (isfunc(fld,'RMSDIFF',arg1,arg2,arg3)) then
c       take the average value over a subdomain previously specified with
c       field=setdom[x1,x2,y1,y2]
        
c       check for consistent number of arguments
        if ((arg2.eq.'#').or.(arg1.eq.'#').or.(arg3.ne.'#')) then
          print *,' Inconsistent number of arguments in ',
     &              fld(strbeg(fld):strend(fld)),'.'
          return
        endif

c       get variable
        ptr = getvar(arg1, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, dimnam, local)
        if (ptr.eq.0) return
        ptru = getvar(arg2, ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, dimnam, local)
        if (ptru.eq.0) return
        do i=ndims+1,MAXDIM
          dims(i)=1
        enddo

c       check whether domain is defined
        if ((sdomx1.ge.sdomx2).or.(sdomy1.ge.sdomy2)) then
          print *,'Subdomain must first be defined with SETDOM.'
          return
        endif
c       check whether domain is a subdomain
        if ((sdomx1.lt.datmin(1)).or.(sdomx2.gt.datmax(1)).or.
     &      (sdomy1.lt.datmin(2)).or.(sdomy2.gt.datmax(2))) then
          print *,'Subdomain not inside full domain!'
          return
        endif

c       convert domain to grid-coordinates
        iflag(1)=1
        iflag(2)=1
        iflag(3)=0
        iflag(4)=0
        call new_field_trans(arg1,%val(ptr),dims(1),dims(2),1,1)
        call new_field_trans(arg1,%val(ptru),dims(1),dims(2),1,1)
        phypt(1)=sdomx1
        phypt(2)=sdomy1
        call phys_2_index_trans(phypt,compt,iflag,4,1)
        idomx1=nint(compt(1))
        idomy1=nint(compt(2))
c       check whether the domain is integer-based
        if ((abs(real(idomx1)-compt(1)).gt.0.1).or.
     &      (abs(real(idomy1)-compt(2)).gt.0.1)) then
          print *,'Subdomain does not coincide with grid-structure'
          return
        endif
        phypt(1)=sdomx2
        phypt(2)=sdomy2
        call phys_2_index_trans(phypt,compt,iflag,4,1)
        idomx2=nint(compt(1))
        idomy2=nint(compt(2))
c       check whether the domain is integer-based
        if ((abs(real(idomx2)-compt(1)).gt.0.1).or.
     &      (abs(real(idomy2)-compt(2)).gt.0.1)) then
          print *,'Subdomain does not coincide with grid-structure'
          return
        endif
c       print out integer subdomain
        print *,'Computing rmserr of ',arg1(strbeg(arg1):strend(arg1)),
     &        ' in'
        print 873,'    ',idomx1,' < ix <',idomx2
        print 873,'    ',idomy1,' < iy <',idomy2

c       get memory
        ptrd=getmem(dims(3)*dims(4))
        if (ptrd.eq.0) return

c       call the routine to do the computation
        call rmsdiff (%val(ptr),%val(ptru),%val(ptrd),
     &           dims(1),dims(2),dims(3),dims(4),
     &           idomx1,idomx2,idomy1,idomy2,misdat)
        calc_field=ptrd
        dims(1)=1
        dims(2)=1
        if ((ndims.lt.3).or.
     &      ((ndims.eq.3).and.(dims(3).eq.1)).or.
     &      ((ndims.eq.4).and.(dims(3).eq.1).and.(dims(4).eq.1))) 
     &      ndims=1
          

c     -----------------------------------------------------------------------
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
     &                data_units, data_display_units, dimnam, local)

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

c     -----------------------------------------------------------------------
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
     &                       dimnam, local)
           if (ptr.eq.0) return
         endif
         size=1
         do i = 1, ndims
            size = size * dims(i)
         enddo
         calc_field  = getmem (size)
         if (calc_field.eq.0) then
           print *,'Could not allocate memory in calc_field'
           return
         endif
c        misdat=0 switches off missing data checking
         misdat=0.0
         data_units=' '
         data_display_units=' '

         call levelar (%val(calc_field),dims(1),dims(2),dims(3),dims(4))
 
c     -----------------------------------------------------------------------
      elseif (fld.eq.'L_COORD') then
c        switch to level-koordinates
         disptype=0
         print *
         print *,' NOW USING ARRAY-INDICES AS VERTICAL COORDINATE.'
         print *,' MAKE SURE TO RESELECT CURRENT BUTTON!.'
         print *

c     -----------------------------------------------------------------------
      elseif (fld.eq.'P_COORD') then
c        switch to pressure-koordinates
         if ((ptrps.gt.0).or.(datatype.eq.2)) then
           disptype=1
           print *
           print *,' NOW USING PRESSURE AS VERTICAL COORDINATE.'
           print *,' MAKE SURE TO RESELECT CURRENT BUTTON!.'
           print *
         elseif (datatype.eq.3) then
c          get pointer to pressure-field
           ptrpre=getvar ('NOMISDAT[P]', ndims, dims, stag, 
     &               datmin, datmax, misdat,
     &               data_units, data_display_units, dimnam, local)
           if ((ptrpre.gt.0).and.(.not.(local))) then
             disptype=1
             print *
             print *,' NOW USING PRESSURE AS VERTICAL COORDINATE.'
             print *,' MAKE SURE TO RESELECT CURRENT BUTTON!.'
             print *
           else 
             print *
             print *,' WARNING: Unable to switch to P-coordinates!'
             print *,'          Compute first FIELD=NOMISDAT[P]'
             print *
           endif
           
         else
           print *
           print *,' WARNING: Unable to switch to pressure coordinates!'
           print *
         endif  

c     -----------------------------------------------------------------------
      elseif (fld.eq.'T_COORD') then
c       switch to theta-koordinates

c       get pointer to theta-field
        if (datatype.eq.3) then
          disptype=2
          print *
          print *,' NOW USING THETA AS VERTICAL COORDINATE.'
          print *,' MAKE SURE TO RESELECT CURRENT BUTTON!.'
          print *
        else
          ptrth=getvar ('THETA', ndims, dims, stag, 
     &             datmin, datmax, misdat,
     &             data_units, data_display_units, dimnam, local)
          if ((ptrth.gt.0).and.(.not.(local))) then
            disptype=2
            print *
            print *,' NOW USING THETA AS VERTICAL COORDINATE.'
            print *,' MAKE SURE TO RESELECT CURRENT BUTTON!.'
            print *
          else 
            print *
            print *,' WARNING: Unable to switch to theta coordinates!'
            print *,'          Compute first FIELD=THETA first!'
            print *
          endif
        endif  

c     -----------------------------------------------------------------------
      elseif (fld.eq.'Z_COORD') then
c        switch to z-koordinates

c       get pointer to height-field
        ptrzl=getvar ('ZLAY', ndims, dims, stag, 
     &                datmin, datmax, misdat,
     &                data_units, data_display_units, dimnam, local)
        if ((ptrzl.gt.0).and.(.not.(local))) then
          disptype=3
          print *
          print *,' NOW USING HEIGHT AS VERTICAL COORDINATE.'
          print *,' MAKE SURE TO RESELECT CURRENT BUTTON!.'
          print *
        else 
          if (ptrzl.gt.0) then
            call freevar('ZLAY')
          endif
          print *
          print *,' WARNING: Unable to switch to height-coordinates!'
          print *,'          Compute first FIELD=ZLAY first!'
          print *
        endif  

c     -----------------------------------------------------------------------
      elseif ((fld.eq.'HELP').or.(fld.eq.'?')) then
c        Print a list of the diagnostic fields.
         call system('more /usr/local/lib/ive/help_em')
c     -----------------------------------------------------------------------
      elseif (isfunc(fld,'WIND',arg1,arg2,arg3)) then
c        calculates the projection of the wind vector to the current
c        cross-section plane (if the passed argument is positive,
c        it is just used to number the field, if negative, it is
c        converted to positive and interpreted as the desired angle)
c        the second argument is not requested, but when set to 1, the
c        vectors are interpolated to a plane perpendicular to the
c        cross-section.
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
     &                  data_units, data_display_units, dimnam, local)
           if (ptr1.eq.0) return
c
c     get the v-component:
           ptr2=getvar('v',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if (ptr2.eq.0) return
c
c       define the dimension of the field and get the memory for it
        wsize=1
        do i = 1, ndims
          wsize = wsize * dims(i)
        enddo
        ptrd  = getmem (wsize)
        if (ptrd.eq.0) return

c       call the subroutine to compute the projected wind:
           call t_proj_wind(%val(ptr1),%val(ptr2),c,d,%val(ptrd),
     &             dims(1),dims(2),dims(3),dims(4),misdat)
c
           calc_field=ptrd
           if (calc_field.eq.0) return
c          assign units (m/s)
           data_units='m/s'
           data_display_units='m/s'
c-----------------------------------------------------------------------------
      elseif (isfunc(fld,'DRAWRANGE',arg1,arg2,arg3)) then
c     the first argument specifies, if coordinate ranges should be tested
c     1:they are tested,0:not
         if (arg1.eq.'#') then
            print *,'first argument needed in: ',fld(1:iend)
            return
         endif
         call t_set_rangecheck(arg1,arg2,arg3)
         calc_field=0
c-----------------------------------------------------------------------------
      elseif (isfunc(fld,'DRAWPARAM',arg1,arg2,arg3)) then
c     set parameters, see t_set_drawparam
         if (arg1.eq.'#') then
            print *,'first argument (=size) needed in: ',fld(1:iend)
            return
         endif
         call t_set_drawparam(arg1,arg2,arg3)
         calc_field=0
c-----------------------------------------------------------------------------
      elseif (isfunc(fld,'DRAW',arg1,arg2,arg3)) then
      if ((arg1.eq.'#').or.(arg2.eq.'#').or.(arg3.eq.'#')) then
        print *,'illegal number of arguments in ',fld(1:iend)
        return
      endif
c     draw text on plot:
c
c     following 2 calls to print to gks metafile:
         call setlvar('savflg',.true.,ierr)
         call makseg !allocate plot segment
c
         call t_drawtext(arg1,arg2,arg3,0)
c
c     following call to finish printing on cgs metafile:
         call gclsg() !push segment
c
         calc_field=0
c-----------------------------------------------------------------------------
      elseif (isfunc(fld,'DRAWFILE',arg1,arg2,arg3)) then
      if ((arg1.eq.'#').or.(arg2.eq.'#')) then
        print *,'illegal number of arguments in ',fld(1:iend)
        return
      endif
c
c     following 2 calls to print to gks metafile:
         call setlvar('savflg',.true.,ierr)
         call makseg !allocate plot segment
c
c     draw text on plot:
           call t_drawfile(arg1,arg2,arg3)
c
c     following call to finish printing on cgs metafile:
         call gclsg() !push segment
c
           calc_field=0
c-----------------------------------------------------------------------------
      elseif (isfunc(fld,'SYMB',arg1,arg2,arg3)) then
c     see the purpose of t_drawsymb about the meaning of argx.
         if ((arg1.eq.'#').or.(arg2.eq.'#').or.(arg3.eq.'#')) then
            print *,'illegal number of arguments in ',fld(1:iend)
            return
         endif
c     draw symbol on plot:
c
c     following 2 calls to print to gks metafile:
         call setlvar('savflg',.true.,ierr)
         call makseg !allocate plot segment
c
         call t_drawsymb(arg1,arg2,arg3)
c     following call to finish printing on cgs metafile:
         call gclsg() !push segment
c
         calc_field=0
c     -----------------------------------------------------------------------
      else

c       -----------------------------------------------------------------------
c       teste ob ein Level oder eine Zeit gewuenscht wird. Dieser Operator
c       hat hoechste Bindungskraft, weshalb alle anderen Operatoren vorgesetzt
c       werden muessen.
        call fnddel(fld,'+-*/^@<%',arg1,arg2,del)
        if (del.eq.'%') then
           ptr = getvar(arg1,ndims, dims, stag, 
     &                       datmin, datmax, misdat,
     &                       data_units, data_display_units, 
     &                       dimnam, local)
           if (ptr.eq.0) return

           if ((arg2(1:1).eq.'T').and.(strend(arg2).gt.1)) then
c            find time-index
             itime=0
             timinp=str2nmb(arg2(2:strend(arg2)),ierr)
             print *,'timinp ',arg2(2:strend(arg2)),timinp,ierr
             if (ierr.eq.0) itime=time2index(timinp)
             if ((itime.lt.1).or.(itime.gt.dims(4))) then
               print *,'Illegal time specified in ',fld(1:strend(fld))
               return
             endif
c            get memory
             size=dims(1)*dims(2)*dims(3)
             calc_field=getmem(size)
             if (calc_field.eq.0) return
c            extract time
             print *,'Take ',arg1(1:strend(arg1)),' at time IT=',itime
             call exttime(%val(calc_field),%val(ptr),
     &           dims(1),dims(2),dims(3),dims(4),itime)
c            define new attributes
             dims(4)=1
             ndims=3

           else if ((arg2(1:1).eq.'K').and.(strend(arg2).gt.1)) then
c            select level-index
             ilevel=str2int(arg2(2:strend(arg2)),ierr)
             if ((ilevel.lt.1).or.(ilevel.gt.dims(3)).or.
     &         (ierr.ne.0)) then
               print *,'Illegal level in ',fld(1:strend(fld))
               return
             endif
             if (datatype.eq.1) ilevel=dims(3)-ilevel+1
c            get memory
             size=dims(1)*dims(2)*dims(4)
             calc_field=getmem(size)
             if (calc_field.eq.0) return
c            extract level
             print *,'Field ',arg1(1:strend(arg1)),' at level K=',ilevel
             call extlevel(%val(calc_field),%val(ptr),
     &           dims(1),dims(2),dims(3),dims(4),ilevel)
c            define new attributes
             dims(3)=1
             if (ndims.eq.3) ndims=2
           else if ((arg2(1:1).eq.'L').and.(strend(arg2).gt.1)
     &       .and.(datatype.gt.1)) then
c            select level-value (for datatype=2 or datatype=3 only)
             rlevel=0.
             rlevel=str2nmb(arg2(2:strend(arg2)),ierr)
             ilevel=0
             do k=1,dims(3)
               if (abs(aklev(k)-rlevel).lt.0.001) ilevel=k
             enddo
             if ((ilevel.lt.1).or.(ilevel.gt.dims(3)).or.
     &         (ierr.ne.0)) then
               print *,'Illegal level in ',fld(1:strend(fld))
               return
             endif
c            get memory
             size=dims(1)*dims(2)*dims(4)
             calc_field=getmem(size)
             if (calc_field.eq.0) return
c            extract level
             print *,'Field ',arg1(1:strend(arg1)),' at level K=',ilevel
             call extlevel(%val(calc_field),%val(ptr),
     &           dims(1),dims(2),dims(3),dims(4),ilevel)
c            define new attributes
             dims(3)=1
             if (ndims.eq.3) ndims=2
           endif
        endif

c       -----------------------------------------------------------------------
c       teste auf die Maskensyntax numb<field<numb, field<numb, numb<field
        call fnddel(fld,'<',arg1,arg2,del)
        if (del.eq.'<') then
c         versuche arg2 als zahl zu interpretieren
          rnumb2=str2nmb(arg2,ierr)
          if (ierr.ne.0) then
            arg0=arg2
            rnumb1=str2nmb(arg1,ierr)
            if (ierr.ne.0) return
            rnumb2=1.e20
          else
            arg0=arg1
c           versuche arg1 abermals auseinanderzubrechen
            call fnddel(arg0,'<',arg1,arg2,del)
            if (del.eq.'<') then
c             versuche arg1 als zahl zu interpretieren
              rnumb1=str2nmb(arg1,ierr)
              if (ierr.ne.0) then
                rnumb1=-1.e20
              else
                arg0=arg2 
              endif
            else
              rnumb1=-1.e20
            endif
          endif

c         get the field
          ptr = getvar(arg0, ndims, dims, stag, 
     &             datmin, datmax, misdat,
     &             data_units, data_display_units, dimnam, local)
          if (ptr.eq.0) return

c         get the memory for the result
          size=1
          do i=1,4
            size=size*dims(i)
          enddo
          calc_field=getmem(size)
          if (calc_field.eq.0) return
c         do the computation
          call copyar(%val(ptr),%val(calc_field),size)
          call mask(%val(calc_field),size,rnumb1,rnumb2,misdat)
        endif

      endif
    
c      if (calc_field.eq.0) then
c        print *,'Unable to compute ',varnam(ibeg:iend)
c      else
c        print *,'Leaving calc_field with ',varnam(ibeg:iend)
c        print *,'ndims,dims',ndims,(dims(i),i=1,ndims)
c        print *,'misdat',misdat
c      endif
      end


      subroutine t_heading(arg1,ierr)
c---------------------------------------------------------------------------
c     purpose:
c        this routine sets the appropriate heading type
c     arguments:
c        arg1   char  input   the type as character:
c                             ive:default as in IVE
c                             date:date, automatically fc and ana decided.
c                             date_fc:forecast
c                             date_ana:analyse
c                             date_ana_day:ana with day
c                             date_fc_day:fc with day
c     history:
c        written by david n. bresch 941024
c-----------------------------------------------------------------------

c     include the common-block with headtyp-variable
      include 'constants.icl'

      character*(80)  arg1
      integer  ierr

c     decide which type:
c      print *,'arg1:|',arg1,'|'
      ierr=0
      if (arg1.eq.'IVE') then
        headtyp=0
      else if (arg1.eq.'DATE') then
        headtyp=1
      else if (arg1.eq.'DATE_ANA') then
        headtyp=2
      else if (arg1.eq.'DATE_FC') then
        headtyp=3
      else if (arg1.eq.'DATE_ANA_DAY') then
        headtyp=4
      else if (arg1.eq.'DATE_FC_DAY') then
        headtyp=5
      else 
        ierr=1
      endif
      end !of t_heading


      subroutine setstack (cfn)
c------------------------------------------------------------------------------
c     reads the specified constants-file, which becomes the default for
c     later use in the commands PUT und INTSTACK. 
c------------------------------------------------------------------------------

      include '/usr/local/include/netcdf.inc'
      include 'constants.icl'
      include 'cdfconsts.icl'
      include 'attributes.icl'

c     arguments
      character*(*)   cfn

c     local variables
      character*(80)  path,datfil,chars
      logical         there,error
      integer         cdfid, i, rcode, varid, k

c     functions
      integer         strbeg,strend

c     statement-functions for the computation of pressure
      real      prlev,prlay,psrf
      integer   is
      prlev(is)=cdfaklev(is)+cdfbklev(is)*psrf
      prlay(is)=cdfaklay(is)+cdfbklay(is)*psrf

c     with cfn='*', a previously selected constants-file is deselected
      if (cfn(1:strend(cfn)).eq.'*') then
        cdfcstfiln='*'
        return
      endif
 
c     prepend constants-filename with the path:
      call getavar('datfil',datfil,error) ! get current filename
c     Cut the filename from datfiln (-> get the path):
      path=' '
      do i=strend(datfil), strbeg(datfil), -1
        if ( (datfil(i:i) .eq. '/').or.(datfil(i:i) .eq. ']') )then
          path=datfil(strbeg(datfil):i)//path(strbeg(path):strend(path))
          go to 80
        endif
      enddo
c     the variable path includes the full path plus name of the selected cfn
 80   path=path(strbeg(path):strend(path))//cfn(strbeg(cfn):strend(cfn))

c     Open the specified netCDF file for reading.
      print *,'Reading constants-file ',
     &          path(strbeg(path):strend(path))
      inquire ( file=path(strbeg(path):strend(path)),exist=there )
      if (there) cdfid = ncopn (path(strbeg(path):strend(path)),
     &          NCNOWRIT,rcode)

      if ((.not.there).or.(rcode.ne.NCNOERR)) then
        write (*,100) path(strbeg(path):strend(path))
        return
      else
c       the file was successfully accessed
        print *,'Opened file ',path(strbeg(path):strend(path)),' '

c       get dimensions
        varid = ncdid (cdfid, 'nz', rcode)
        call ncdinq (cdfid, varid, chars, cdfnz, rcode)
        if (rcode.ne.NCNOERR) goto 950

c       get coordinate-coefficients for levels (Schichtgrenzen)
        varid = ncvid (cdfid, 'aklev', rcode)
        call ncvgt (cdfid, varid, 1, cdfnz, cdfaklev, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncvid (cdfid, 'bklev', rcode)
        call ncvgt (cdfid, varid, 1, cdfnz, cdfbklev, rcode)
        if (rcode.ne.NCNOERR) goto 950

c       get coordinate-coefficients for layers (Schichtmitten)
        varid = ncvid (cdfid, 'aklay', rcode)
        call ncvgt (cdfid, varid, 1, cdfnz, cdfaklay, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncvid (cdfid, 'bklay', rcode)
        call ncvgt (cdfid, varid, 1, cdfnz, cdfbklay, rcode)
        if (rcode.ne.NCNOERR) goto 950

c       define data-type
        rcode=0
        do k=1,cdfnz
          if (cdfbklev(k).ne.0.) rcode=1
          if (cdfbklay(k).ne.0.) rcode=1
        enddo
        if (cdfnz.eq.1) then
          cdfdatatype=0
        else if (cdfaklev(1).gt.cdfaklev(cdfnz)) then
          if (rcode.eq.0) then
            cdfdatatype=2
          else
            cdfdatatype=1
          endif
        else
          if (rcode.eq.0) then
            cdfdatatype=3
          else
            cdfdatatype=0
          endif
        endif
        print *,'datatype is ',cdfdatatype

c       define constants-name
        cdfcstfiln=cfn

c       compute vertical derivatives of ak's and bk's
        do k=2,cdfnz-1
          cdfaslay(k)=(cdfaklay(k-1)-cdfaklay(k+1))/2.
          cdfbslay(k)=(cdfbklay(k-1)-cdfbklay(k+1))/2.
          cdfaslev(k)=(cdfaklev(k-1)-cdfaklev(k+1))/2.
          cdfbslev(k)=(cdfbklev(k-1)-cdfbklev(k+1))/2.
        enddo
        cdfaslay(1 )=cdfaklay(1)-cdfaklay(2)
        cdfbslay(1 )=cdfbklay(1)-cdfbklay(2)
        cdfaslev(1 )=cdfaklev(1)-cdfaklev(2)
        cdfbslev(1 )=cdfbklev(1)-cdfbklev(2)
        cdfaslev(cdfnz)=cdfaklev(cdfnz-1)-cdfaklev(cdfnz)
        cdfbslev(cdfnz)=cdfbklev(cdfnz-1)-cdfbklev(cdfnz)
        cdfaslay(cdfnz)=cdfaklay(cdfnz-1)-cdfaklay(cdfnz)
        cdfbslay(cdfnz)=cdfbklay(cdfnz-1)-cdfbklay(cdfnz)
      endif

      write (6,*)
      write (6,*) '  k(ive)  k(em)      levels      layers',
     &                           '       aklev       bklev'
      psrf=1000.
      write (6,160),0,cdfnz+1,psrf
      do k=1,cdfnz
       write (6,160) k,cdfnz+1-k,prlev(k),prlay(k),
     &                                cdfaklev(k),cdfbklev(k)
      enddo
      return

c     error-exits
 950  continue
      write (6,110) path(strbeg(path):strend(path))
      return

c     formats
 100  format (1x,'Unable to open file ',a,'.')
 110  format (1x,'Unable to read level-information from file ',a,'.')
 160  format (1x,2i7,4f12.2)

      end


      subroutine intstack (fld,ar,ie,je,ke,le,misdat,ierr)
c------------------------------------------------------------------------------
c     Interpolates the field named FLD to a stack of levels.
c     The levels are specified through a constants-file
c     which is previously selected with a call to SETSTACK, and may
c     be pressure or theta-levels
c------------------------------------------------------------------------------
      include 'constants.icl'
      include 'cdfconsts.icl'

c     arguments
      character*(*)      fld
      integer            ie,je,ke,le,ierr
      real               ar(ie,je,ke,le), misdat

c     local variables    
      integer   MAXDIM
      parameter (MAXDIM=4)
      character*(80)   data_units, data_display_units,chars
      integer          ndims, dims(MAXDIM)
      real             stag(MAXDIM), misdati,
     &                 datmin(MAXDIM), datmax(MAXDIM)
      logical          local
      integer          k,ptr
      character*(80)   dimnam(MAXDIM)

c     functions
      integer          getvar,strbeg,strend

c     set error-flag
      ierr=0

c     do loop over all the levels
      do k=1,ke
c       get the variable FLD interpolated to the levels        
        write (chars,*) cdfaklev(k)
        if (cdfdatatype.eq.3) then
          chars='INTER[' // fld(1:strend(fld)) // ':THETA:' //
     &                    chars(strbeg(chars):strend(chars)) // ']'
        else if (cdfdatatype.eq.2) then
          chars='INTER[' // fld(1:strend(fld)) // ':P:' //
     &                    chars(strbeg(chars):strend(chars)) // ']'
        else
          ierr=1
          return
        endif
        ptr=0
        ptr=getvar(chars(1:strend(chars)), ndims, dims, stag, 
     &           datmin, datmax, misdat,
     &           data_units, data_display_units, dimnam, local)
        if (ptr.eq.0) then
          print *,'Unable to get ',chars(1:strend(chars)),
     &                                        ' in INTSTACK'
          ierr=1
          return
        endif
        call copylevel(%val(ptr),ar,ie,je,ke,le,k)
      enddo
      end




      subroutine pcdfvar (filnam, varnam, ptr, newtim, ndims, dims,
     &                          stag, datmin, datmax, misdat, 
     &                          data_units, data_display_units, error) 
c------------------------------------------------------------------------------
c     puts the data associated with 'ptr' as variable 'varnam' to 
c     the cdf-file 'filnam'. The constants-filename written to the file is
c     either that previously selected with a call to SETSTACK, or is identical
c     to that associated with the current file.
c------------------------------------------------------------------------------

      include '/usr/local/include/netcdf.inc'
      include 'constants.icl'
      include 'cdfconsts.icl'
c     include the common-block with the time-array information
      include 'attributes.icl'

c     Argument declarations. MAXDIM is the maximum data dimension allowed,
c     LENVAR is the maximum length of the variable-name.
      integer   MAXDIM, LENVAR, MAXVAR
      parameter (MAXDIM=4, LENVAR=80, MAXVAR=80)
      character*(80)   data_units, data_display_units
      character*(*)    filnam,varnam
      integer          ndims, dims(MAXDIM), ptr, error
      real             newtim, stag(MAXDIM), misdat,
     &                 datmin(MAXDIM), datmax(MAXDIM)

c     local variables
      character*(80)   varlst(MAXVAR),cfn,char
      character*(1)    del
      integer          cdfid,cstid,ierr,i,idate(5),istdate(5)
      integer          nvar,datmin1(MAXDIM),datmax1(MAXDIM),ndims1

c     externals
      integer          strend

c     try to open the NetCDF file with name given by filnam
      call ncpopt(NCVERBOS)
      call opncdf(filnam(1:strend(filnam)),cdfid,datmin1,datmax1,ndims1,
     &               varlst,nvar,cfn,ierr)

      if (ierr.ne.0) then
c       if an error occured create the file and define the variable. First
c       specify the constants-filename
        if (cdfcstfiln.ne.'*') then
c         take the constants-name specified in a previous call to setstack 
          cfn=cdfcstfiln(1:strend(cdfcstfiln))
        else
c         take the constants-name of the currently active file
          call fnddel(cstfiln(1:strend(cstfiln)),'/',char,cfn,del)
        endif
        call crecdf(filnam(1:strend(filnam)),cdfid,datmin,datmax,
     &             3,cfn(1:strend(cfn)),ierr)
        call putdef(cdfid,varnam,ndims,misdat,dims,datmin,datmax,
     &              stag,ierr)
        write(*,*)'defined variable ',varnam(1:strend(varnam))
      else   ! check if variable is already defined
        do i=1,nvar
          if (varlst(i).eq.varnam) goto 201
        enddo
        call putdef(cdfid,varnam,ndims,misdat,dims,datmin,datmax,
     &              stag,ierr)
        write(*,*)'defined variable ',varnam(1:strend(varnam))
      endif
  201 continue
c     put the field on the file
      call putit(cdfid,varnam,%val(ptr),dims(1),dims(2),dims(3),
     &           dims(4),newtim)
      call clscdf(cdfid,ierr)
      end


      subroutine putit(cdfid,vnam,ar,dims1,dims2,dims3,dims4,newtim)
c     ===========================================
c     writes the array ar on the NetCDF file with identifier
c     cdfid by calling the subroutine putdat (from the library
c     libcdfio) for every time step.

c     include the common-block with the time-array information
      include 'attributes.icl'

      integer   cdfid
      integer   dims1,dims2,dims3,dims4
      real      ar(dims1,dims2,dims3,dims4)
      real      newtim,newtim2
      character*(*) vnam
     
      integer   ierr,i
      integer   strend
      
      if (dims4.eq.1) then
c       define newtim
        if (newtim.ne.-99.99) then
          newtim2=newtim
        else
          newtim2=timeval(1)
        endif
        call putdat(cdfid,vnam,newtim2,0,ar(1,1,1,1),ierr)
        write(*,*)'variable ',vnam(1:strend(vnam)),' at time ',
     &             newtim2, ' written'
      else
        do i=1,dims4
          call putdat(cdfid,vnam,timeval(i),0,ar(1,1,1,i),ierr)
          write(*,*)'variable ',vnam(1:strend(vnam)),' at time ',
     &               timeval(i), ' written'
        enddo
      endif

      return
      end


      subroutine cdfopn(filnam,cdfid,ierr)
C------------------------------------------------------------------------

C     Opens the NetCDF file 'filnam' and returns its identifier cdfid.

C     filnam    char    input   name of NetCDF file to open
C     cdfid     int     output  identifier of NetCDF file
C     ierr      int     output  error flag
C------------------------------------------------------------------------

      include '/usr/local/include/netcdf.inc'

      integer   cdfid,ierr
      character*(*) filnam

      integer   strend

      call ncpopt(NCVERBOS)
      cdfid=ncopn(filnam(1:strend(filnam)),NCNOWRIT,ierr)

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
        
      
      subroutine t_intpole(ar,ie,je,ke,le,datmin,datmax,misdat)
c----------------------------------------------------------------------------
c     interpolates the pole-values
c----------------------------------------------------------------------------
      integer     ie,je,ke,le,i,k,l,cnt
      real        ar(ie,je,ke,le),datmin(2),datmax(2),misdat,sum

c     do north-pole (j=je)
      if ((datmax(2)-90.).lt.1.e-3) then
        do l=1,le
        do k=1,ke
          sum=0.
          cnt=0
          do i=1,ie
            if ((ar(i,je-1,k,l).ne.misdat).or.(misdat.eq.0.)) then
              cnt=cnt+1
              sum=sum+ar(i,je-1,k,l)
            endif
          enddo
          if (2*cnt.ge.ie) then
            sum=sum/real(cnt)
            do i=1,ie
              ar(i,je,k,l)=sum
            enddo
          endif
        enddo
        enddo
      endif

c     do south-pole (j=1)
      if ((datmin(2)+90.).lt.1.e-3) then
        do l=1,le
        do k=1,ke
          sum=0.
          cnt=0
          do i=1,ie
            if ((ar(i,2,k,l).ne.misdat).or.(misdat.eq.0.)) then
              cnt=cnt+1
              sum=sum+ar(i,2,k,l)
            endif
          enddo
          if (2*cnt.ge.ie) then
            sum=sum/real(cnt)
            do i=1,ie
              ar(i,1,k,l)=sum
            enddo
          endif
        enddo
        enddo
      endif
      end


      subroutine int2iso (nam1,nam2,ar1,ar2,valint,res,nx,ny,nz1,nz2,nt,
     &        datmin1,datmax1,stag1,misdat1,
     &        datmin2,datmax2,stag2,misdat2)
c-----------------------------------------------------------------------
c     Purpose:
c        Interpolates the field ar1 onto the isosurface ar2=valint.
c     Arguments:
c        nam*     char   inp    name of the field ar1 and ar2
c        ar*      char   inp    the fields
c        valint   real   inp    the value of the isosurface
c        res      real   out    the resulting field, dimensioned (nx,ny,1,nt)
c        nx,etc   int    inp    the dimensions
c        datmin,datmax,stag: no longer used, could be removed from list
c        misdat1  real inp/out  the missing data value of field 1. If misdat1=0
c                               on input, it can be overwritten with a new
c                               missing data value.
c        misdat2  real   out    the missing data value of field 2.
c-----------------------------------------------------------------------
      
c     arguments
      character*(*) nam1,nam2
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
c                 nach oben
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
            call new_field_trans (nam2,ar2,nx,ny,nz2,nt)
            call index_2_phys_trans(phypt,compt,iflagp,4,1)
  
c           compute location in computational space of ar1
            call new_field_trans (nam1,ar1,nx,ny,nz1,nt)
            call phys_2_index_trans(phypt,compt,iflagc,4,1)

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

c     overwrite previous missing data value
      if (missing) misdat1=misdath
      end



      subroutine int2prs (nam,ar,prsval,res,nx,ny,nz,nt,
     &                    datmin,datmax,stag,misdat)
c-----------------------------------------------------------------------
c     Purpose: interpolates the 4d field ar(nx,ny,nz,nt) onto the
c     pressure surface prsval, resulting in the array res(nx,ny,nz)
c-----------------------------------------------------------------------
      
c     arguments
      integer      nx,ny,nz,nt
      real         ar(nx,ny,nz,nt),res(nx,ny,nt)
      real         datmin(4),datmax(4),stag(4),misdat
      real         prsval
      character*(*) nam

c     local variables
      integer      i,j,l,iflagp(4),iflagc(4)
      logical      missing
      real         compt(4),phypt(4),misdath,int4dm

c     initialize iflag
      data         iflagp/1,1,0,1/
      data         iflagc/0,0,1,0/

      if (misdat.eq.0.) then
        misdath=-999.99
        missing=.false.
      else
        misdath=misdat
        missing=.true.
      endif

      do i=1,nx
      do j=1,ny
      do l=1,nt
            call new_field_trans(nam,ar,nx,ny,nz,nt)

c           compute location in physical space
            compt(1)=real(i)
            compt(2)=real(j)
            compt(3)=(nz+1)/2
            compt(4)=real(l)
            call index_2_phys_trans(phypt,compt,iflagp,4,1)

c           compute location in computational space
            phypt(3)=prsval
            call phys_2_index_trans(phypt,compt,iflagc,4,1)

            if ((compt(3).ge.1).and.(compt(3).le.real(nz))) then
c             get value through interpolation
              compt(1)=real(i)
              compt(2)=real(j)
              compt(4)=real(l)
              res(i,j,l)=int4dm(ar,nx,ny,nz,nt,
     &             compt(1),compt(2),compt(3),compt(4),misdat)
            else
              res(i,j,l)=misdath
              missing=.true.
            endif
      enddo
      enddo
      enddo

c     overwrite previous missing data value
      if (missing) misdat=misdath
      end


      subroutine red4d (a,af,nx,ny,nz,nt,ired,misdat)
c     ===============================================
c     this subroutine reduces the horizontal resolution of the field a by some
c     factor ired.

c     argument declarations
      integer  nx,ny,nz,nt,ired
      real     a(nx,ny,nz*nt),af(nx,ny,nz*nt),misdat

c     variable declarations
      integer  i,j,k,i1,i2,ii,j1,j2,jj,m
      real     cnti,cnt,sum

      m=ired/2
      if (m*2-ired.ne.0) then
        print *,'Reducing-Factor should be a multiple of 2!'
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
            sum=sum+cnti*a(ii,jj,k)
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







