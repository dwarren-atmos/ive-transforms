      integer function privat (varnam, ndims, dims, stag, datmin, 
     &                           datmax, misdat, 
     &                           data_units, data_display_units, dimnam) 
c-----------------------------------------------------------------------
c     $Id: privat.f,v 1.1 1994/11/14 22:37:22 warren Exp $
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
c	$Log: privat.f,v $
c Revision 1.1  1994/11/14  22:37:22  warren
c Christoph Schaer's European Model transforms.
c
c Revision 1.5  1994/11/11  09:02:12  schaer
c Aenderung fuer IVE-3-2-beta
c Verbessertes Heading
c Startzeit neu auch auf Datenfile
c
c Revision 1.4  1993/12/06  10:25:11  schaer
c Back to revision 1.2
c
c Revision 1.2  1993/09/02  07:25:49  schaer
c Kleinere Aenderungen damit es schoener aussieht.
c
c Revision 1.1  1993/08/27  07:01:24  schaer
c Initial revision
c
c---------------------------------------------------------------------- 



c     do compiler dependent definition of recursive variables

c     include the common-block with the level information
      include 'constants.icl'

c     Argument declarations. MAXDIM is the maximum data dimension allowed,
c     LENVAR is the maximum length of the variable-name.
      integer   MAXDIM, LENVAR
      parameter (MAXDIM=4, LENVAR=80)

      character*(80)   data_units, data_display_units
      character*(80)   dimnam(MAXDIM)
      character*(LENVAR)  varnam
      integer          ndims, dims(MAXDIM)
      real             stag(MAXDIM), misdat,
     &                 datmin(MAXDIM), datmax(MAXDIM)

c     Local variable declarations
      character*(LENVAR) fld1,fld2,arg1,arg2,arg3,dvar
      real             misdatu,misdatv,misdvel,misdatd
      integer          ptru,ptrv,ptrt,ptrvel,ptrd,ptrrr,size,ptrdps,ierr
      logical          errflg,localu,localv

c     Local variable declarations.
      character*(LENVAR) fld
      integer          i, ibeg, iend, ibeg1, iend1, ptr, dims1(MAXDIM)
      logical          err, local
      real             delta1,delta2

c     External function declarations.
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
           privat = getvar('(U^2+V^2)^0.5', ndims, dims, stag, 
     &             datmin, datmax, misdat,
     &             data_units, data_display_units, dimnam, local)
c       All intermediate operations (e.g. reading data from the data-file,
c       doing computations) are hidden to the user. If getvar returns
c       successfully (i.e. calc_field<>0), then a new button denoted 'vel'
c       will be installed after return from calc_field.
        return

c     -----------------------------------------------------------------------
      else if (fld.eq.'PVALL') then

c          addvar und getvar sind in ~schaer/uwgap7/uars/getvar.F

           privat=getvar('-(VORT+F)*9.80665*D[THETA:P]',
     &                  ndims, dims, stag, 
     &                  datmin, datmax, misdat,
     &                  data_units, data_display_units, dimnam, local)
           if (privat.eq.0) return
c          assign units (PVU's)
           data_units='K m**2 / (kg s)'
           data_display_units='1.E-6  K m**2 / (kg s)'

      endif

      end
