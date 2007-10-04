#include "pointer.icl"

      ive_ptr function calc_field(name, ndims, dims, stag, min,
     &   max, missing, data_units, data_display_units, dim_names)
c     
c     calc_field : This routine is used to calculate user-derived fields.
c     The return value of the routine is a pointer to the field values.
c     NULL return => cannot calculate.
c     
c     Arguments:
c     name	character	The name of the field to derive.
c     ndims	integer		Number of dimensions in field (output).
c     dims	integer		Number of points in the field in Fortran
c				order (x, y, z, t) (output).
c     stag	real		Grid staggering per dimension (output).
c     min	real		Physical space minimum per dimension (output).
c     max	real		Physical space maximum per dimension (output).
c     missing	real		Missing data value, zero => none (output).
c     data_units
c		character	Units for field (output).
c     data_display_units
c		character	Units to use to display field (output).
c     dim_names character	Names of the dimensions in the field (output).
c

c #include "default.icl"

      integer ndims, dims(4)
      real    stag(4), min(4), max(4), missing
      character*(*) name, data_units, data_display_units, dim_names(4)

      ive_ptr derivative, getvar
      external getvar
c
c     This is where you put the code to do your own calculation.
c
c     -----------------------------------------------------------------
c
c     Here is an example.  Say you have the temperature and pressure fields,
c     and you want to calculate the potential temperature:
c
c                         THETA = T(1000/P)**0.286
c
c     There are TWO ways to code the new field variable THETA in the
c     CALC_FIELD routine:
c
c     ----------------------- First Approach --------------------------
c
c       The first approach is to pass the IVE routine GETVAR a mathematical 
c     expression identical to that which could be entered directly on the
c     IVE command line using the FIELD command (see IVE Help document
c     entitled "Diagnostic Computations Using the Command Line).
c
      integer ibeg,iend,strbeg,strend
      logical flag

c      ibeg = strbeg(name)  ! These commands find the beginning and ending
c      iend = strend(name)  !  characters in the variable "name".
c      if (name(ibeg:iend) .eq. 'THETA') then
c
c   !  The following command returns a pointer to the memory location where
c   !  the field 'THETA' will be stored.  Note that it is necessary to use
c   !  IVE syntax in the mathematical expression for THETA.
c
c         print *,'in calc_field'
c         calc_field = getvar('TMK*(1000/PRS)^0.286', ndims, dims, stag,
c     &	             ,min,max, missing, data_units, data_display_units,
c     &			  dim_names, flag) ! This command returns a pointer
c					   !  to the memory location where
c					   !  the field 'THETA' is stored.
c         print *,'still in calc_field'
c         return
c      endif
c
c	In this example, "FIELD=T*(1000/P)^0.286, THETA" could have
c     been entered at the IVE command line to define the potential 
c     temperature field.  However, adding the above code to CALC_FIELD
c     would enable IVE to automatically recognize THETA as a diagnostic 
c     variable.  This allows you to streamline your IVE sessions and 
c     avoid typographical errors in frequently evaluated functions.
c
c     ---------------------- Second Approach --------------------------
c       
c	Some field calculations are too complex to code at the IVE command
c     line.  For example, performing a field variable integration would
c     require the user to loop through the array values.  Such manipu-
c     lations are very involved and are only feasible by grabbing the 
c     field variables and programming the calculations yourself.  A simple
c     example of this second approach is shown below:
c     (Again, we will calculate THETA and we will assume that T and P are
c     on the same grid.  If T and P were not on the same grid, it would be
c     necessary to interpolate the field values using the interpolation
c     operator '@' defined in the IVE Help file entitled "Diagnostic 
c     Computations at the Command Line".)
C
      ive_ptr t, p, getmem
      real pmissing,t1,p1,r_val
      integer num, i
C
      ibeg = strbeg(name)  ! These commands find the beginning and ending
      iend = strend(name)  !  characters in the variable "name".
      if (name(ibeg:iend) .eq. 'THETA') then
         p = getvar('PRS', ndims, dims, stag, min, max,
     &               pmissing, data_units, data_display_units,
     &		     dim_names, flag)  !  Here, getvar returns a pointer
				       !   to the memory location where
				       !   the field 'P' is stored.

         if (p .eq. 0) then  ! This 'if' loop flags user when 'P' missing  
c             make_help_widget('THETA - cannot get P')
             print *,'THETA - cannot get P'
             calc_field = 0
             return
         endif

         t = getvar('TMK', ndims, dims, stag, min, max,
     &	             missing, data_units, data_display_units, dim_names,
     & 		     flag)  ! Here, getvar returns a pointer to the memory
			    !  location where the field 'T' is stored.
         if (t .eq. 0) then  ! This 'if' loop flags user when 'T' missing  
c             make_help_widget('THETA - cannot get T')
             print *,'THETA - cannot get T'
             calc_field = 0
             return
         endif

         num = dims(1) * dims(2) * dims(3) * dims(4)
         calc_field = getmem(num)  ! This command allocates array memory.
         if (calc_field .eq. 0) then  ! This 'if' loop flags user when
				       !  the memory cannot be allocated.
c	       make_help_widget('THETA - cannot allocate memory')
             print *,'THETA - cannot allocate memory'
	       calc_field = 0
              return
         endif

         do i = 1, num
              t1 = r_val(%val(t), i)  ! These commands return the i-th 
              p1 = r_val(%val(p), i)  !  data value in the arrays p & t
              if (t1 .eq. missing .or. p1 .eq. missing) then
c         If the i-th value in either p or t is missing, the following
c         command stores the value 'missing' in the i-th value of the
c         array THETA
                 call s_val(%val(calc_field), i, missing)
              else
c         If the data is not missing, we calculate the i-th value of the
c         array THETA

                 call s_val(%val(calc_field), i, t1*(100000./p1)**0.286)
              endif
         enddo
         return

C code for ptot

      elseif (name(ibeg:iend) .eq. 'PTOT') then
         p = getvar('RTC', ndims, dims, stag, min, max,
     &               pmissing, data_units, data_display_units,
     &		     dim_names, flag)  !  Here, getvar returns a pointer
				       !   to the memory location where
				       !   the field 'P' is stored.

         if (p .eq. 0) then  ! This 'if' loop flags user when 'P' missing  
             print *,'PTOT - cannot get RTC'
             calc_field = 0
             return
         endif

         t = getvar('RTE', ndims, dims, stag, min, max,
     &	             missing, data_units, data_display_units, dim_names,
     & 		     flag)  ! Here, getvar returns a pointer to the memory
			    !  location where the field 'T' is stored.
         if (t .eq. 0) then  ! This 'if' loop flags user when 'T' missing  
             print *,'PTOT - cannot get RTE'
             calc_field = 0
             return
         endif

         num = dims(1) * dims(2) * dims(3) * dims(4)
         calc_field = getmem(num)  ! This command allocates array memory.
         if (calc_field .eq. 0) then  ! This 'if' loop flags user when
				       !  the memory cannot be allocated.
             print *,'PTOT - cannot allocate memory'
	       calc_field = 0
              return
         endif

         do i = 1, num
              t1 = r_val(%val(t), i)  ! These commands return the i-th 
              p1 = r_val(%val(p), i)  !  data value in the arrays p & t
              if (t1 .eq. missing .or. p1 .eq. missing) then
c         If the i-th value in either p or t is missing, the following
c         command stores the value 'missing' in the i-th value of the
c         array THETA
                 call s_val(%val(calc_field), i, missing)
              else
c         If the data is not missing, we calculate the i-th value of the
c         array THETA

                 call s_val(%val(calc_field), i, t1+p1)
              endif
         enddo
         return

C code for wind speed

      elseif (name(ibeg:iend) .eq. 'WSPD') then
         print *,'here'
         p = getvar('UUU', ndims, dims, stag, min, max,
     &               pmissing, data_units, data_display_units,
     &		     dim_names, flag)  !  Here, getvar returns a pointer
				       !   to the memory location where
				       !   the field 'P' is stored.

         if (p .eq. 0) then  ! This 'if' loop flags user when 'P' missing  
             print *,'WSPD - cannot get UUU'
             calc_field = 0
             return
         endif

         t = getvar('VVV', ndims, dims, stag, min, max,
     &	             missing, data_units, data_display_units, dim_names,
     & 		     flag)  ! Here, getvar returns a pointer to the memory
			    !  location where the field 'T' is stored.
         if (t .eq. 0) then  ! This 'if' loop flags user when 'T' missing  
             print *,'WSPD - cannot get VVV'
             calc_field = 0
             return
         endif

         num = dims(1) * dims(2) * dims(3) * dims(4)
         calc_field = getmem(num)  ! This command allocates array memory.
         if (calc_field .eq. 0) then  ! This 'if' loop flags user when
				       !  the memory cannot be allocated.
             print *,'WSPD - cannot allocate memory'
	       calc_field = 0
              return
         endif

         do i = 1, num
              t1 = r_val(%val(t), i)  ! These commands return the i-th 
              p1 = r_val(%val(p), i)  !  data value in the arrays p & t
              if (t1 .eq. missing .or. p1 .eq. missing) then
c         If the i-th value in either p or t is missing, the following
c         command stores the value 'missing' in the i-th value of the
c         array THETA
                 call s_val(%val(calc_field), i, missing)
              else
c         If the data is not missing, we calculate the i-th value of the
c         array THETA

                 call s_val(%val(calc_field), i, sqrt(t1**2 + p1**2))
              endif
         enddo
         return
      endif

c
c     -----------------------------------------------------------------
c
c     "derivative" is a routine in IVE that calculates the derivative
c     of a field if the syntax of "name" is D[field:dir:type]. If "name"
c     is not in that format, or if "derivative" cannot do the calculation,
c     it returns NULL (0).
c     
      calc_field = derivative(name, ndims, dims, stag, min, max,
     &      missing, data_units, data_display_units, dim_names)
      return
      end
c
