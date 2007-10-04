/*
 * $Id: nps.c,v 1.6 1997-12-18 20:32:55 harry Exp $
 *
 * Transform routines for the NPS grids.  The following three grid types
 * are supported:
 *
 * 1. Latitude-Longitude
 * 2. Polar stereographic with the pole point at the North or South Pole.
 *    The following global attributes must be in the netCDF file:
 *
 *	orientation	- lat/long of the pole point
 *	pole_point	- the grid point coordinates of the pole
 *	grid_increment	- the grid spacing in km at 60 degrees latitude
 *
 * 3. Lambert conformal tangent cone projection (i.e. Latin1 = Latin2)
 *    with equal x and y spacing of the grid. The netCDF file must have
 *    a "nav" variable with the following attributes:
 *
 *	grid_type	- Must be "Lambert conformal"
 *	La1		- Latitude of first grid point
 *	Lo1		- Longitude of first grid point
 *	Lov		- Orientation Longitude
 *	Dx		- X-direction grid length
 *	Dy		- Y-direction grid length (= Dx)
 *	Latin1		- First latitude cutting the cone
 *	Latin2		- Second latitude cutting the cone (= Latin1)
 *
 * Note: Math for Lambert conformal is from NMC programs w3fb11 & w3fb12
 * written by John Stackpole.
 *
 * $Log: nps.c,v $
 * Revision 1.6  1997-12-18 20:32:55  harry
 * Add field_delta.
 *
 * Revision 1.5  1997/09/04 21:30:06  harry
 * Move map routines out of the DEBUG ifdef.
 *
 * Revision 1.4  1995/12/15 23:50:26  harry
 * Add call to derivative in calc_field.
 *
 * Revision 1.3  1995/10/05  22:12:38  harry
 * Change nps to reflect changes in new_file.
 *
 * Revision 1.2  1995/09/05  20:35:59  harry
 * Modifications for new slicer.
 *
 * Revision 1.1  1994/08/30  17:18:20  harry
 * Moves to trans_nps directory
 *
 * Revision 1.13  1994/06/23  21:42:42  harry
 * Change index to strchr to avoid UCB library under Solaris 2.
 *
 * Revision 1.12  1994/06/14  05:24:56  harry
 * Fixes for Solaris 2.
 *
 * Revision 1.11  1994/02/24  19:23:51  harry
 * Fix index_2_phys for the case where the Z dimension in linear.
 *
 * Revision 1.10  1994/02/17  22:08:38  harry
 * Change read_var to getvar and add dim_names parameter.
 *
 * Revision 1.9  1993/12/22  00:17:40  harry
 * Add support for Lambert Conformal grids.
 *
 * Revision 1.8  1993/12/17  21:40:42  harry
 * Add support for Polar Stereographic grids.
 *
 * Revision 1.7  1993/11/24  19:19:45  harry
 * Correct typo (extact_times => exact_times).
 *
 * Revision 1.6  1993/11/24  19:14:48  harry
 * Add Dales comments to the README.trans file and the transform programs.
 *
 * Revision 1.5  1993/10/18  22:30:30  harry
 * Pass an int to malloc for the number of bytes to allocate.
 *
 * Revision 1.4  1993/09/23  23:09:10  harry
 * Change call to calc_field.
 *
 * Revision 1.3  1993/09/16  18:28:40  harry
 * Put log P values into separate array rather than writing over the
 * original values (perhaps more than one time).
 *
 * Revision 1.2  1993/08/13  23:08:40  harry
 * Change LATLON to LONLAT.
 *
 * Revision 1.1  1993/07/09  23:15:09  harry
 * Initial checkin of default and nps transforms, along with a script to
 * pull the default transforms from loadfunc.c in IVE.
 *
 *
 */

static char ident[] = "@(#)$Id: nps.c,v 1.6 1997-12-18 20:32:55 harry Exp $";

#include <netcdf.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define RDC (180./M_PI)
#define DRC (M_PI/180.)

/*
 * Common variables
 */

static int   frtime_size,	/* Size of frtime array */
             polar_stereo,	/* Is this a polar stereographic grid? */
             lamb_conf;		/* Is this a lambert conformal grid */

static float *frtime,		/* Forecast times */
             *frtime_sorted	/* Sorted frtime */
			=NULL,
             *logP=NULL;	/* Log P values */

static float *coord[4],		/* An array containing pointers to arrays
				   holding the "dimension variables".  The
				   use of dimension variables is optional,
				   but allows one to specify irregular grid
				   intervals along each spatial coordinate. */
             inmax[4],		/* Maximum array indices of the data. */
             phmin[4],		/* Minimum physical_space locations for
				   each dimension */
             phmax[4],		/* Maximum physical_space locations for
				   each dimension */
             delta[4];		/* If non-zero => physical distance between
				   any two array indicies is constant and delta
				   is the value.  If zero => physical distance
				   is not constant, transforms use coord to
				   do translations */
/*
 * Variables needed for Polar Stereographic.
 */
static float orientation[2],	/* Pole lat/lon for Polar Stereographic grid */
             pole_point[2],	/* Point for north/south pole */
             grid_increment;	/* Grid interval at 60N in km */
/*
 * Variables needed for Lambert Conformal.
 */
static float polex, poley,	/* Pole points */
             h,			/* 1 for NH, -1 for SH */
             Latin1,		/* Tangent latitude */
             Lov,		/* Orientation longitude in degrees */
             an,		/* h*sin(Latin1) */
             thing,		/* Partial calculation for phys_2_lonlat */
             rebydx,		/* radius of earth divided by Dx */
             cosltn;		/* cos(Latin1) */

static int fcompare(i, j) /* Needed for qsort */
float *i, *j;
{
    if (*i==*j) return(0);
    else if (*i>*j) return(1);
    else return(-1);
}

/*
 *  get_attribute_unit splits an attribute into a real value and its units.
 */

static double
get_attribute_unit(ncid, varid, attname, units)
int ncid, varid;
char *attname, *units;
{
    char *attribute, *p;
    int len;
    nc_type type;
    double result;

    units[0] = '\0';
    if (ncattinq(ncid, varid, attname, &type, &len) == -1) {
	return(0.0);
    }
    attribute = (char *) malloc(len+1);
    ncattget(ncid, varid, attname, attribute);
    attribute[len] = '\0';
    if (p=(char *)strchr(attribute, ' '))
	strcpy(units, p+1);
    result = atof(attribute);
    free(attribute);
    return(result);
}

/*
 * new_file_ : This routine is called every time a file is read in.
 * The routine can be used to read other information from the open
 * netCDF file.  It also can provide an array of coord values that limit the
 * values the window sliders can be set to.
 *
 * Arguments:
 * ncid	 	int *	  input	 Id of currently opened netCDF file.
 * exact_coord	float **  output Pointer to array of discrtet coord values
 *				 that each window slider should use.
 * n		int *	  output Number of values in exact_coord.  If this
 *				 is 0 => ignore values in exact_coord
 *				 => no restrictions on window slider.
 * coord_dep	int[4][4] output States the functional dependencies
 *                               between the coordinates, 1 => dependent,
 *				 0 => independent e.g. if z is a function of
 *				 x, y and z, then coord_dep[2][0],
 *				 coord_dep[2][1] and coord_dep[2][2]
 *				 will be one.  A coordinate is always 
 *				 assumed to be dependent on itself.
 */

void
new_file_(ncid, exact_coord, n, coord_dep)
int *ncid, n[4], coord_dep[4][4];
float *exact_coord[4];

{
    float *data,	/* Data returned by getvar_ */
          missing,	/* For getvar_ call */
          plmin[4],	/* Domain minimum values */
          plmax[4];	/* Domain maximum values */
    int ndims, dims[4], /* Dimension information from getvar_ */
        stagger[4], min[4], max[4],
        flag, i, j, m, error;	/* Scratch */
    char data_units[81], display_units[81],
         dim_names[4][81]; /* For getvar_ */
    /*
       Structure cdf_dimvar_ (common block cdf_dimvar) is set up
       by IVE to contain the values of a variable with the same name as the
       dimension of the largest variable in the netCDF file.  If the
       size of the domain is set by attributes, or if a dimension
       variable does not exist, size will be zero, and data will be NULL.
       */

    extern struct cdf_dimvar_ {
	float *data[4];
	int size[4];
    } cdf_dimvar_;

/*
  Constants for Fortran calls
  */
    int four=4;

    for (i=0; i < 4; ++i) 
	for (j=0; j < 4; ++j) coord_dep[i][j] = (i==j);
    /*
      Blank out unneeded variables.
      */
    set_button_name_("reftime", "", 7, 0);
    set_button_name_("grib_center", "", 11, 0);
    set_button_name_("grib_model", "", 10, 0);
    set_button_name_("nav", "", 3, 0);
    /*
      Get forecast times.
      */
    flag = 0;
    if ((frtime=(float* )getvar_("frtime", &ndims, dims, stagger,
				min, max, &missing, data_units,
				display_units, dim_names, &flag,
				6, 80, 80, 80)) == NULL) {
	make_help_widget_("Missing \"frtime\" variable");
	return;
    }
    frtime_size = dims[0];
    /*
      Get domain limits.
      */
    getrarr_("plmin", plmin, &four, &error, 5);
    getrarr_("plmax", plmax, &four, &error, 5);
    /*
      Look for missing level -> 3 D variable (lat, lon, time) or
      (x, y, time). If missing, set up proper values of plmin, plmax,
      and the domain units.
      */
    if (getvid_("level", 5) < 0) {
	char domain_units[4][81], domain_display_units[4][81];

	memset(domain_units, ' ', sizeof(domain_units));
	memset(domain_display_units, ' ', sizeof(domain_display_units));
	plmin[3] = frtime[0];
	plmax[3] = frtime[frtime_size-1];
	strncpy(domain_units[3], data_units, 81);
	strncpy(domain_display_units[3], display_units, 81);
	flag = 0;
	if ((data=(float *)getvar_("lon", &ndims, dims, stagger,
				   min, max, &missing, domain_units[0],
				   domain_display_units[0], dim_names, &flag,
				   3, 80, 80, 80)) != NULL ||
	    (data=(float *)getvar_("x", &ndims, dims, stagger,
				   min, max, &missing, domain_units[0],
				   domain_display_units[0], dim_names, &flag,
				   1, 80, 80, 80)) != NULL) {
	    plmin[0] = data[0];
	    plmax[0] = data[dims[0]-1];
	}
	flag = 0;
	if ((data=(float *)getvar_("lat", &ndims, dims, stagger,
				   min, max, &missing, domain_units[0],
				   domain_display_units[0], dim_names, &flag,
				   3, 80, 80, 80)) != NULL ||
	    (data=(float *)getvar_("y", &ndims, dims, stagger,
				   min, max, &missing, domain_units[0],
				   domain_display_units[0], dim_names, &flag,
				   1, 80, 80, 80)) != NULL) {
	    plmin[1] = data[0];
	    plmax[1] = data[dims[0]-1];
	}
	plmin[2] = 1.;
	plmax[2] = 1.;
	domain_units[2][0] = '\0';
	domain_display_units[2][0] = '\0';
	setaarr_("domain_units", domain_units, &four, &error, 12, 81);
	setaarr_("domain_display_units", domain_display_units, &four, &error,
		 20, 81);
	setaarr_("domain_display_units_orig", domain_display_units, &four,
		 &error, 25, 81);
    }
    for (i=0; i < 3; ++i) {
	exact_coord[i] = cdf_dimvar_.data[i];
	n[i] = cdf_dimvar_.size[i];
    }
    /*
      Get the forecast times and sort to get correct time limits.
      */
    if (frtime_sorted) frtime_sorted =
	(float *) realloc(frtime_sorted, frtime_size*sizeof(float));
    else frtime_sorted = (float *) malloc(frtime_size*sizeof(float));
    for (i=0; i < frtime_size; ++i) frtime_sorted[i] = frtime[i];
    qsort(frtime_sorted, frtime_size, sizeof(float), fcompare);
    plmin[3] = frtime_sorted[0];
    plmax[3] = frtime_sorted[frtime_size-1];
    if (memcmp(frtime, frtime_sorted, frtime_size*sizeof(float)) == 0) {
	free(frtime_sorted);
	frtime_sorted = NULL;
	exact_coord[3] = frtime;
    }
    else exact_coord[3] = frtime_sorted;
    n[3] = frtime_size;
    setrarr_("plmin", plmin, &four, &error, 5);
    setrarr_("plmax", plmax, &four, &error, 5);
    /*
       Get grid orientation for map transforms if this is a polar
       stereographic grid.
       */
    flag = ncopts;
    ncopts = 0;
    if (polar_stereo=(ncattget(*ncid, NC_GLOBAL,
				      "orientation", (void *) orientation)
	!= -1)) {
	ncattget(*ncid, NC_GLOBAL, "pole_point", (void *) pole_point);
	ncattget(*ncid, NC_GLOBAL, "grid_increment", (void *) &grid_increment);
    }
    /*
       Get grid orientation for map transforms if this is a 
       Lambert conformal grid.
       */
    if (lamb_conf=((m=getvid_("nav", 3)) >= 0
		   && ncattinq(*ncid, m, "grid_type", NULL, &i)
		   && (i == 17 || i == 18)
		   && ncattget(*ncid, m, "grid_type", (void *) display_units)
		   && strncmp(display_units, "Lambert conformal", 17)==0)) {
	
	float La1, Lo1, Dx, Dy, Latin2, rmll;
	float slope, intercept;
	char units[80];

	La1 = get_attribute_unit(*ncid, m, "La1", display_units);
	if (strcmp(display_units, "degrees_south") == 0)
	    La1 = -La1;
	Lo1 = get_attribute_unit(*ncid, m, "Lo1", display_units);
	if (strcmp(display_units, "degrees_west") == 0)
	    Lo1 = -Lo1;
	Lov = get_attribute_unit(*ncid, m, "Lov", display_units);
	if (strcmp(display_units, "degrees_west") == 0)
	    Lo1 = -Lo1;
	Dx = get_attribute_unit(*ncid, m, "Dx", display_units);
	if (convert_units_(display_units, "km", &slope, &intercept,
			   strlen(display_units), 2) == 0) {
	    Dx = Dx*slope + intercept;
	}
	Latin1 = get_attribute_unit(*ncid, m, "Latin1", display_units);
	if (strcmp(display_units, "degrees_south") == 0)
	    Latin1 = -Latin1;
        if (Latin1 > 0) h = 1;
	else h = -1;
	rebydx = 6371.2/Dx;
	cosltn = cos(Latin1*DRC);
	an = h * sin(Latin1*DRC);
	thing = pow(an/rebydx, 1./an)/(pow(cosltn, (1.-an)/an)*(1+an));
	rmll = rebydx*pow(cosltn, 1.-an)*pow(1.+an, an) *
			pow(cos(La1*DRC)/(1.+h*sin(La1*DRC)),an)/an;
	polex = 1. - h*rmll*sin(an*(Lo1-Lov)*DRC);
	poley = 1. + rmll*cos(an*(Lo1-Lov)*DRC);
    }
    ncopts = flag;
}

/*
 * new_field_ : This routine is called every time a new field is set.
 * The routine can be used to store necessary information to perform
 * the other transformations.  It also can use used to reorder the
 * data array.
 * 
 * Note: before IVE calls new_field it has already called
 * getvar to read in the data for the new field.  As a
 * result of this call, IVE will have also read in
 * ndims, dims, stagger, phmin, phmax, missing, data_units
 * and display_units---provided this information appears
 * in the netCDF.  The transform routines can access this
 * data by calling getrarr, getlvar or getaarr.

 * Arguments:
 * name		char *	input	The name of the current field.
 * field	float *	input	Values of "field".
 * inmax1	int *	input	Size of dimensions in "Fortran"
 * - inmax4			order (i.e., A(inmax1,...inmax4)
 *				will hold the data with the first
 *				index varying most rapidly)
 * len		int	input	Length of name in characters.
 */

void
new_field_(name, field, inmax1, inmax2, inmax3, inmax4, len)
char *name;
float **field;
int *inmax1, *inmax2, *inmax3, *inmax4, len;

{
    char dnames[4][MAX_NC_NAME],	/* Names of "field" dimensions */
         data_units[81], display_units[81], dim_names[4][81]; /* For getvar_ */
    float *s, *t, *temp,		/* Scratch */
          *missing, stagger[4];		/* For getvar_ call */
    int zero=0, four=4,			/* Constants for Fortran calls */
        i, j, k,			/* Scratch */
        min[4], max[4], dims[4], nd,	/* For getvar_ call */
        error,				/* Error flag for getrarr_ */
        read_from_nc,			/* Was this field read from the 
                                           netCDF file => may need to 
					   reorder the time dimension */
        size,				/* Size of time slabs to move */
        n,				/* Variable number for the dimension */
        ndims,				/* Number of dimensions in "field" */
        flag;				/* For getvar_ call */

    /*
       Save the physical and index limits for use in later calculations
       by other transform routines.
       */
    getrarr_("phmin", phmin, &four, &error, 5);
    getrarr_("phmax", phmax, &four, &error, 5);
    inmax[0] = *inmax1;
    inmax[1] = *inmax2;
    inmax[2] = *inmax3;
    inmax[3] = *inmax4;
    /*
       Get the names of the dimensions for this field.
       */
    getaarr_("data_dim_names", dnames, &four, &error, 14, MAX_NC_NAME);
    getivar_("num_dims", &ndims, &error, 8);
    getlvar_("read_from_nc", &read_from_nc, &error, 12);
    /*
       If this variable was just read from the netCDF file and the
       forecast times are out of order, then reorder the data array.
       */
    if (read_from_nc && frtime_sorted && strcmp(name, "frtime") != 0 &&
	strcmp(dnames[ndims-1], "frtime")==0) {
	size = inmax[0];
	for (i=1; i < ndims-1; ++i) size *= inmax[i];
	temp = (float *)malloc(frtime_size*sizeof(float));
	memcpy(temp, frtime, frtime_size*sizeof(float));
	for (i=0; i < frtime_size; ++i) {
	    if (temp[i] != frtime_sorted[i]) {
		/*
		  Exchange correct time with incorrect time.
		  */
		for (j=i+1; j < frtime_size; ++j) {
		    if (temp[j] == frtime_sorted[i]) break;
		}
		for (s=*field+size*i, t=*field+size*j, k=0; k < size; k++) {
		    float f=*s;
		    *s++ = *t;
		    *t++ = f;
		}
	        temp[j] = temp[i];
	    }
	}
	free(temp);
    }
    /*
       Fix parameters for fields that are missing the "Z" dimension.
       */
    if (ndims == 3 && strcmp(dnames[2], "frtime") == 0) {
	inmax[3] = inmax[2];
	inmax[2] = 1;
	*inmax4 = *inmax3;
	*inmax3 = 1;
	strcpy(dnames[3], dnames[2]);
	strcpy(dnames[2], "");
	phmin[3] = phmin[2];
	phmin[2] = 1;
	phmax[3] = phmax[2];
	phmax[2] = 1;
	getrarr_("stag", stagger, &four, &error, 4);
	stagger[3] = stagger[2];
	stagger[2] = 0;
	setrarr_("stag", stagger, &four, &error, 4);
	ndims++;
    }
    /*
       Fix and store the variable limits.
       */
    if (frtime_sorted) {
	phmin[3] = frtime_sorted[0];
	phmax[3] = frtime_sorted[frtime_size-1];
    }
    setrarr_("phmin", phmin, &four, &error, 5);
    setrarr_("phmax", phmax, &four, &error, 5);
    /*
       Set up P for logrithmic interpolation.
       */
    if (inmax[2] != 1) {
	phmin[2] = log10((double) phmin[2]);
	phmax[2] = log10((double) phmax[2]);
    }
    /*
       Get values (if they exist) in the 1D array associated with each
       dimension to see which dimensions are evenly spaced.
       delta is zero if they are not.
       */
    for (i=0; i < ndims; ++i) {
        /*
           Check that the data array depth exceeds one along this coordinate
           */
	if (inmax[i] != 1) {
            /*
               delta is the distance between grid point in this dimension if
               the points are evenly spaced.
               */
	    delta[i] = (phmax[i]-phmin[i])/(inmax[i]-1);
	    /*
	       Get the memory address of the values of the variable with
	       the same name as the dimension.
	       If such a variable does not exist, coord[i] == NULL.
                   */
	    flag = 0;
	    if (strcmp(dnames[i], "frtime") == 0 && frtime_sorted)
		coord[i] = frtime_sorted;
            else coord[i] =
		(float *)getvar_(dnames[i], &nd, dims, stagger,
				 min, max, &missing,
				 data_units, display_units, dim_names,
				 &flag, strlen(dnames[i]), 80, 80, 80);
	    if (coord[i]) {
		if (i == 2) {
		    if (logP) logP =
			(float *) realloc(logP, (int) inmax[2]*sizeof(float));
		    else logP = (float*) malloc((int) inmax[2]*sizeof(float));
		    for(j=0; j < inmax[2]; ++j)
			logP[j] = log10((double) coord[2][j]);
		    coord[2] = logP;
		}
                /*
                   Check to see if the dimension is evenly spaced by getting
                   the first value and then comparing the other values
                   to the value calculated using the delta computed above.
                   If the dimension is not evenly spaced, set delta to 0.
                   */
		for(j=1; j < inmax[i]; ++j) {
		    if (coord[i][0]+j*delta[i] != coord[i][j]) {
			delta[i] = 0.;
			break;
		    }
		}
	    }
	}
	setrarr_("field_delta", delta, &four, &error, 10);
    }
}

/*
 * index_2_phys_ : This routine translates array indices into physical
 * coordinates.
 *
 * Arguments:
 * phys		float *	output	Array of physical coordinates.
 * index	float *	input	Array of indicies.
 * flag		int *	input	Flags that indicate which coordinates
 *				to compute.  If iflag[i] = 1 index[j][i]
 *				will be converted to phys[j][i], otherwise
 *				the value of phys[j][i] is not changes.
 * ndims	int *	input	Number of dimensions in each point.
 * npts		int *	input	Number of points to convert.
 */

void
index_2_phys_(phys, index, flag, ndims, npts)
int *flag, *ndims, *npts;
float *phys, *index;

{
    int i, j, k, idx;

    for (idx=0, i=0; i < *npts; i++) {
	for(j=0; j<*ndims; j++, idx++) {
	    if(flag[j]) {
                /*
                   Trivial case where dimension does not exist (depth of
                   array is unity along this dimension)
                   */
		if (inmax[j] == 1) phys[idx] = phmin[j];
                /*
                   Case where interpolation must be performed between the grid
                   point locations given in a 1D "dimension variable" array.
                   */
		else {
		    if (delta[j] == 0.) {
			k = index[idx];
			if (k < inmax[j]) {
			    phys[idx] = coord[j][k-1] + 
				(index[idx]-k)*(coord[j][k]-coord[j][k-1]);
			}
			else
			    phys[idx] = phmax[j];
		    }
		    /*
		       Case of evenly spaced physical coord.
		       */
		    else phys[idx] = phmin[j] + (index[idx]-1)*delta[j];
		    /*
		       Covert pressure from log to normal.
		       */
		    if (j==2) phys[idx] = pow((double) 10.,
					      (double) phys[idx]);
		}
	    }
	}
    }
}

/*
 * phys_2_index_ : This routine translates physical coordinates into 
 * array indices.
 *
 * Arguments:
 * phys		float *	input	Array of physical coordinates.
 * index	float *	output	Array of indicies.  A negative index =>
 *				this point gets the missing data value
 *				(e.g. is under the terrain).
 * flag		int *	input	Flags that indicate which coordinates
 *				to compute.  If iflag[i] = 1 phys[j][i]
 *				will be converted to index[j][i], otherwise
 *				the value of index[j][i] is not changes.
 *	ndims	int *	input	Number of dimensions in each point.
 *	npts	int *	input	Number of points to convert.
 */

void phys_2_index_(phys, index, flag, ndims, npts)
int *flag, *ndims, *npts;
float *phys, *index;

{
    int i, j, k, idx;
    float p;

    for (idx=0, i=0; i < *npts; i++) {
	for(j=0; j<*ndims; j++, idx++) {
	    if(flag[j]) {
		/*
		  Interpolate log P.
		  */
		if (j == 2) p = log10((double) phys[idx]);
		else p = phys[idx];
		/*
                   Trivial case where dimension does not exist (depth of
                   array is unity along this dimension)
                   */
		if (inmax[j] == 1) index[idx] = 1.;
                /*
                   Case where interpolation must be performed between the grid
                   point locations given in a 1D "dimension variable" array.
                   */
		else if (delta[j] == 0.) {
		    /*
		      Physical coord. increasing with index.
		      */
		    if (phmin[j] < phmax[j]) {
			if (p <= phmin[j]) index[idx] = 1;
			else if (p >= phmax[j]) index[idx] = inmax[j];
			else {
			    for (k=1; k < inmax[j] && p > coord[j][k]; ++k);
			    index[idx] = (float) k - (p-coord[j][k])/
				(coord[j][k-1]-coord[j][k]) + 1.;
			}
		    }
		    /*
		      Physical coord. decreasing with index.
		      */
		    else {
			if (p >= phmin[j]) index[idx] = 1;
			else if (p <= phmax[j]) index[idx] = inmax[j];
			else {
			    for (k=1; k < inmax[j] && p < coord[j][k]; ++k);
			    index[idx] = (float) k - (p-coord[j][k])/
				(coord[j][k-1]-coord[j][k]) + 1.;
			}
		    }
		}
		/*
		  Case of evenly spaced physical coord.
		  */
		else index[idx] = (p-phmin[j])/delta[j] + 1.;
	    }
	}
    }
}

/*
 * phys_2_lonlat_ : This routine converts physical coordinates to longitude-
 * latitude coordinates.
 *
 * Arguments:
 *	x	float *		Physical x coordinate array.
 *	y	float *		Physical y coordinate array.
 *	lon	float *		Longitude array (output).
 *	lat	float *		Latitude array (output).
 *	npts	int *		Number of points to convert.
 */

void
phys_2_lonlat_(x, y, lon, lat, npts)
float *x, *y, *lon, *lat;
int *npts;
{
    int i, scale;
    double r;

    if (polar_stereo) {
	if (orientation[0] >= 0) scale = 1;
	else scale = -1;
	for (i=0; i < *npts; ++i) {
	    lon[i] = atan2(x[i]-pole_point[0],
			   -scale*(y[i]-pole_point[1])) * RDC + orientation[1];
	    r = sqrt((x[i]-pole_point[0])*(x[i]-pole_point[0]) +
		     (y[i]-pole_point[1])*(y[i]-pole_point[1]));
	    lat[i] = scale*90 - 2*scale*atan(r*grid_increment/11888.86) * RDC;
	}
    }
    else if (lamb_conf) {
	for (i=0; i < *npts; ++i) {
	    if((r=(x[i]-polex)*(x[i]-polex)+(poley-y[i])*(poley-y[i])) == 0) {
		lon[i] = Lov;
		lat[i] = h*90;
	    }
	    else {
		lon[i] = Lov + RDC*atan2(h*(x[i]-polex), poley-y[i])/an;
		lat[i] = h * (0.5*M_PI - 2.*atan(thing*pow(r, 1/(an+an))))*RDC;
	    }
	}
    }
    else {
	for (i=0; i < *npts; ++i) {
	    lon[i] = x[i];
	    lat[i] = y[i];
	}
    }
}

/*
 * lonlat_2_phys_ : This routine converts longitude-latitude coordinates
 * to physical coordinates.
 *
 * Arguments:
 *	x	float *		Physical x coordinate array (output).
 *	y	float *		Physical y coordinate array (output).
 *	lon	float *		Longitude array.
 *	lat	float *		Latitude array.
 *	npts	int *		Number of points to convert.
 */

void
lonlat_2_phys_(x, y, lon, lat, npts)
float *x, *y, *lon, *lat;
int *npts;
{
    int i, scale;
    double ang, r;

    if (polar_stereo) {
	if (orientation[0] >= 0) scale = 1;
	else scale = -1;
	for (i=0; i < *npts; ++i) {
	    r = 11888.86*tan((90-scale*lat[i])*DRC/2.)/grid_increment;
	    ang = (lon[i] - orientation[1])*DRC;
	    x[i] = r*sin(ang)+pole_point[0];
	    y[i] = -scale*r*cos(ang)+pole_point[1];
	}
    }
    else if (lamb_conf) {
	for (i=0; i < *npts; ++i) {
	    r = rebydx * pow(cosltn, 1.-an)*pow(1.+an, an) *
		pow(cos(lat[i]*DRC)/(1.+h*sin(lat[i]*DRC)), an)/an;
	    x[i] = polex + h*r*sin(an*(lon[i]-Lov)*DRC);
	    y[i] = poley - r*cos(an*(lon[i]-Lov)*DRC);
	}
    }
    else {
	for (i=0; i < *npts; ++i) {
	    x[i] = lon[i];
	    y[i] = lat[i];
	}
    }
}

/*
 * default_map_ : This routine returns the map settings when the default
 * map is specified.
 *
 * Arguments:
 * proj		char *	output	EZMAP map projection.
 * plon		float *	output	Projection longitude.
 * plat		float *	output	Projection latitude.
 * rota		float * output	Projection rotation.
 * limit	char *	output	limit and plm1-4 are the arguments for
 * plm1		float *	output	  the MAPSET call which specifies the
 * plm2		float *	output	  rectagular portion of the u/v plane
 * plm3		float *	output	  to be drawn.
 * plm4		float *	output
 * exact_fit	int *	output	If 1 => the contour plot will fit exactly
 *				  in this map with no transform.  Thus IVE
 *				  will just call conpack within the map
 *				  window without calling phys_2_lonlat.
 * len1		int	input	Maximum length of proj in characters.
 * len2		int	input	Maximum length of limit in characters.
 */

void
default_map_(proj, plon, plat, rota, limit, plm1, plm2, plm3, plm4,
                   exact_fit, len1, len2)
char *proj, *limit;
float *plon, *plat, *rota, *plm1, *plm2, *plm3, *plm4;
int *exact_fit;
int len1, len2;
{
    float plwmin[4], plwmax[4];
    int one=1, four=4, error;

    if (polar_stereo) {
	getrarr_("plwmin", plwmin, &four, &error, 6);
	getrarr_("plwmax", plwmax, &four, &error, 6);
	strcpy(proj, "ST");
	*plon = orientation[1];
	*plat = orientation[0];
	*rota = 0.;
	/*
	   11888.86 = 1.866*earth_radius
	   */
	strcpy(limit, "LI");
	*plm1 = (plwmin[0]-pole_point[0])*grid_increment/11888.86;
	*plm2 = (plwmax[0]-pole_point[0])*grid_increment/11888.86;
	*plm3 = (plwmin[1]-pole_point[1])*grid_increment/11888.86;
	*plm4 = (plwmax[1]-pole_point[1])*grid_increment/11888.86;
	*exact_fit = 1;
    }
    else if (lamb_conf) {
	getrarr_("plwmin", plwmin, &four, &error, 6);
	getrarr_("plwmax", plwmax, &four, &error, 6);
	strcpy(proj, "LC");
	*plon = Lov;
	*plat = *rota = Latin1;
	strcpy(limit, "CO");
	phys_2_lonlat_(plwmin, plwmin+1, plm2, plm1, &one);
	phys_2_lonlat_(plwmax, plwmax+1, plm4, plm3, &one);
	*exact_fit = 1;
    }
    else {
	/*
	   The default is the current window in a Cylindrical Equidistant
	   Projection.
	   */
	strcpy(proj, "CE");
	*plon = 0.;
	*plat = 0.;
	*rota = 0.;
	/*
	   Setting limit to CL tells IVE to clip to the current window using
	   phys_2_lonlat.
	   */
	strcpy(limit, "CL");
	*exact_fit = 0;
    }
}

#ifdef DEBUG

/*
  The rest of these are the same as the default transforms and are not
  necessary.
  */

/*
 * horiz_ter_ : This routine fills an array with terrain heights in the
 * windowed domain for horizontal cross sections.
 *
 * Arguments:
 *	topo	float *		Array containing the terrain heights (output).
 *	nxw	int *		Number of x points.
 *	nyw	int *		Number of y points.
 *	stagi	float *		Grid staggering along x dimension.
 *	stagj	float *		Grid staggering along y dimension.
 *	zero	float *		Minimum terrain height.
 *	error	int *		Error flag 0 => no errors (output).
 */

void
horiz_ter_(topo, nxw, nyw, stagi, stagj, zero, error)
int *nxw, *nyw, *error;
float *topo, *stagi, *stagj, *zero;
{
    int i;

    for (i=0; i < *nxw * *nyw; topo[i++] = *zero);
    *error = 0;
}

/*
 * vert_ter_ : This routine fills an array with terrain heights in the
 * windowed domain for vertical cross sections.
 *
 * Arguments:
 *	xter	float **	Pointer to array of x-locations for topo_ht
 *	yter	float **	Pointer to array of y-locations for topo_ht
 *	topo_ht float **	Pointer to array containing 
 *				the terrain heights.
 *	nter	float *		Number of points in topo_ht (output)
 *	pt1	float *		Minimum physical coordinates for slice.
 *	pt2	float *		Maximum physical coordinates for slice.
 *	zero	float *		Minimum terrain height.
 */

void
vert_ter_(xter, yter, topo_ht, nter, pt1, pt2, zero)
int *nter;
float **xter, **yter, **topo_ht, *pt1, *pt2, *zero;
{
    int i;

    *nter = 2;
    *xter = (float *)malloc(2*sizeof(float));
    *yter = (float *)malloc(2*sizeof(float));
    *topo_ht = (float *)malloc(2*sizeof(float));
    *xter[0] = pt1[0];
    *xter[1] = pt2[0];
    *yter[0] = pt1[1];
    *yter[1] = pt2[1];
    *topo_ht[0] = *topo_ht[1] = *zero;
}

/*
 * data_slicer_1d_ : This routine returns the computational points at
 * which a 1d slice should be taken.
 *
 * Arguments:
 * endpt	float *	input	2 endpoints of the line in physical space
 * ni		int *	output	Number of points
 * da		int *	input	axis coordinate of average
 * nda		int *	output	Number of points in the da-direction.
 */

float *
data_slicer_1d_(endpt, ni, da, nda)

float *endpt;
int *ni, *da, *nda;
{
#ifndef MAX
#define MAX(x, y) ((x) > (y)? (x):(y))
#endif
#ifndef MIN
#define MIN(x, y) ((x) < (y)? (x):(y))
#endif
    /*
       Let's do the obvious, just divide the slice evenly, using the
       number of computation points represented by the slice.
       */
    struct point4a {
	float v[4];
    } *pendpt = (struct point4a *) endpt, cendpt[2], *pslice, *slice, *p;
    int iflag[4], four=4, two=2, i, j, k, error;
    float dv[4], dda, cpmin[4], cmin, cmax;

    for (i=0; i < 4; i++) iflag[i] = 1;
    phys_2_index_(endpt, cendpt, iflag, &four, &two);
    *ni = 0;
    for (i=0; i < 4; i++) 
	if (i != *da-1) *ni += abs(cendpt[1].v[i]-cendpt[0].v[i]);
    for (i=0; i < 4; i++) dv[i] = (pendpt[1].v[i]-pendpt[0].v[i])/ (*ni - 1);
    if (*da == 0) {
	pslice = (struct point4a *) memalign(sizeof(struct point4a),
					   sizeof(struct point4a)**ni);
	slice = (struct point4a *) memalign(sizeof(struct point4a),
					    sizeof(struct point4a)**ni);
	p = pslice;
	for (i=0; i < *ni; i++, p++) {
	    for (j=0; j < 4; j++) p->v[j] = pendpt[0].v[j] + i*dv[j];
	}
	j = *ni;
    }
    else {
	dv[*da-1] = 0;
	getiarr_("cpmin",cpmin,&four,&error,5,4);
	if (cendpt[0].v[*da-1] < 0 ) cmin = cpmin[*da-1];
	else cmin = cendpt[0].v[*da-1];
	if (cendpt[1].v[*da-1] < 0 ) cmax = cpmin[*da-1];
	else cmax = cendpt[1].v[*da-1];
	if (cmax == cmin) *nda = 1;
	else {
	    *nda = ceil(cmax) - floor(cmin) + 1;
	    dda = (pendpt[1].v[*da-1] - pendpt[0].v[*da-1])/ (*nda - 1);
	}
	pslice = (struct point4a *) memalign(sizeof(struct point4a),
				     sizeof(struct point4a)*(*ni * *nda));
	slice = (struct point4a *) memalign(sizeof(struct point4a),
				    sizeof(struct point4a)*(*ni * *nda));
	p = pslice;
	for (i=0; i < *ni; i++) {
	    for (j=0; j < 4; j++) p->v[j] = pendpt[0].v[j] + i*dv[j];
	    p++;
	    for (k=1; k < *nda; k++, p++) {
		p->v[*da-1] = pendpt[0].v[*da-1] + k*dda;
	    }
	}
	j = *ni * *nda;
    }
    phys_2_index_(pslice, slice, iflag, &four, &j);
    free(pslice);
    return((float *) slice);
}

/*
 * data_slicer_2d_ : This routine returns the computational points at
 * which a 2d slice should be taken.
 *
 * Arguments:
 * corner	float *	input	4 corner points of slice in physical space.
 * ii		int *	input	axis coordinate of ni.
 * jj		int *	input	axis coordinate of nj.
 * ni		int *	output	Number of points in the x-direction.
 * nj		int *	output	Number of points in the y-direction.
 * da		int *	input	axis coordinate of average
 * nda		int *	output	Number of points in the da-direction.
 *
 * The corners will be:     2    3
 *                      jj
 *                          0    1
 *                            ii
 */

float *
data_slicer_2d_(corner, ii, jj, ni, nj, da, nda)

float *corner;
int *ii, *jj, *ni, *nj, *da, *nda;
{
    /*
       Use pslicer2d_
       */
    float *pslicer2d_();

    return(pslicer2d_(corner, ii, jj, ni, nj, da, nda));
}

/*
 * calc_field_ : This routine is used to calculate user-derived field.
 * The return value of the routine is a pointer to the field values.
 * NULL return => cannot calculate.
 *
 * Arguments:
 *	name	char *		The name of the field to derive.
 *	ndims	int *		Number of dimensions in field (output).
 *	dims	int *		Number of points in the field in Fortran
 *				order (x, y, z, t) (output).
 *	stag	real *		Grid staggering per dimension (output).
 *	min	real *		Physical space minimum per dimension (output).
 *	max	real *		Physical space maximum per dimension (output).
 *	missing	float *		Missing data value, zero => none (output).
 *	data_units
 *		char *		Units for field (output).
 *	data_display_units
 *		char *		Units to use to display field (output).
 *	dim_names
 *		char *		Names of the dimensions.
 *	len1	int		Number of characters in name.
 *	len2	int		Number of characters in data_units.
 *	len3	int		Number of characters in data_display_units.
 *	len4	int		Number of characters in dim_names.
 */

float *
calc_field_(name, ndims, dims, stag, min, max, missing, data_units,
	    data_display_units, dim_names, len1, len2, len3, len4)
int *ndims, *dims, len1, len2, len3, len4;
float *stag, *min, *max, *missing;
char *name, *data_units, *data_display_units, *dim_names;
{
    float *derivative_();

    return(derivative_(name, ndims, dims, stag, min, max, missing,
		       data_units, data_display_units, dim_names,
		       len1, len2, len3, len4));
}


/*
 * heading_ : This routine allows the user to change headings on plots
 *
 * Arguments:
 *      where   int *           Wherefrom called
 *      which   int *           Which heading (first (1) or second (2) line)
 *      line    char *          The line to be changed.
 *
 * written by david@atmos.umnw.ethz.ch
 */

void
heading_(where, which, line)
     int  *where;
     int  *which;
     char *line;
{
}

#endif
