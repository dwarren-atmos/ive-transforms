/*
 * Default transforms from IVE.
 */

#include <netcdf.h>
#include <stdio.h>
#include <math.h>

/*
 * Common variables
 */

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
    int i, j;

    for (i=0; i < 4; ++i) 
	for (j=0; j < 4; ++j) coord_dep[i][j] = (i==j);
    /*
       Structure cdf_dimvar_ (common block cdf_dimvar) is set up
       by IVE to contain the values of a variable with the same name as the
       dimension of the largest variable in the netCDF file.  If the
       size of the domain is set by attributes, or if a dimension
       variable does not exist, size will be zero, and data will be NULL.
       */
    {
    extern struct cdf_dimvar_ {
	float *data[4];
	int size[4];
    } cdf_dimvar_;

    for (i=0; i < 4; ++i) {
	exact_coord[i] = cdf_dimvar_.data[i];
	n[i] = cdf_dimvar_.size[i];
    }
    }
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
    float *missing, stagger[4];		/* For getvar_ call */
    int zero=0, four=4,			/* Constants for Fortran calls */
        i, j,				/* Scratch */
        min[4], max[4], nd,     	/* For getvar_ call */
        n,				/* Variable number for the dimension */
        ndims,				/* Number of dimensions in "field" */
        flag,				/* For getvar_ call */
	dims[4];			/* For getvar_ call */
    int error;   			/* Error flag for getrarr_ */

    /*
       Save the physical and index limits for use in later calculations
       by other transform routines.
       */
    getrarr_("phmin", phmin, &four, &error, 5);
    getrarr_("phmax", phmax, &four, &error, 5);
    getrarr_("grid_delta", delta, &four, &error, 10);
    inmax[0] = *inmax1;
    inmax[1] = *inmax2;
    inmax[2] = *inmax3;
    inmax[3] = *inmax4;
    /*
       Get the names of the dimensions for this field.
       */
    getaarr_("data_dim_names", dnames, &four, &error, 14, MAX_NC_NAME);
    getivar_("num_dims", &ndims, &error, 8);
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
	    if(delta[i] == 0.0) delta[i] = (phmax[i]-phmin[i])/(inmax[i]-1);
            /*
               Get the memory address of the values of the variable with
               the same name as the dimension.
               If such a variable does not exist, coord[i] == NULL.
               */
	    flag = 0;
            if ((coord[i] =
                 (float *)read_var_(dnames[i], &nd, dims, stagger,
				    min, max, &missing,
				    data_units, display_units, dim_names,
				    &flag, strlen(dnames[i]), 80, 80, 80))
                != NULL) {
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
    }
    (void) setrarr_("field_delta", delta, &four, &error, 11);
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
 *				If iflag[i] comes back as -1 => insufficient
 *				information to compute this coordinate.
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
		else if (delta[j] == 0.) {
		    k = index[idx];
		    if (k < inmax[j])
			phys[idx] = coord[j][k-1] + 
			    (index[idx]-k)*(coord[j][k]-coord[j][k-1]);
		    else
			phys[idx] = phmax[j];
		}
		/*
		  Case of evenly spaced physical coord.
		  */
		else phys[idx] = phmin[j] + (index[idx]-1)*delta[j];
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
 *				If iflag[i] comes back as -1 => insufficient
 *				information to compute this coordinate.
 *	ndims	int *	input	Number of dimensions in each point.
 *	npts	int *	input	Number of points to convert.
 */

void phys_2_index_(phys, index, flag, ndims, npts)
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
			if (phys[idx] <= phmin[j]) index[idx] = 1;
			else if (phys[idx] >= phmax[j])
			    index[idx] = inmax[j];
			else {
			    for (k=1; k < inmax[j] &&
				 phys[idx] > coord[j][k]; ++k);
			    index[idx] =(float)k - (phys[idx]-coord[j][k])/
				(coord[j][k-1]-coord[j][k]) + 1.;
			}
		    }
		    /*
		      Physical coord. decreasing with index.
		      */
		    else {
			if (phys[idx] >= phmin[j]) index[idx] = 1;
			else if (phys[idx] <= phmax[j])
			    index[idx] = inmax[j];
			else {
			    for (k=1; k < inmax[j] &&
				 phys[idx] < coord[j][k]; ++k);
			    index[idx] =(float)k - (phys[idx]-coord[j][k])/
				(coord[j][k-1]-coord[j][k]) + 1.;
			}
		    }
		}
		/*
		  Case of evenly spaced physical coord.
		  */
		else index[idx] = (phys[idx]-phmin[j])/delta[j] + 1.;
	    }
	}
    }
}

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
 *	error	int *	        Error flag 0 => no errors (output).
 */

void
horiz_ter_(topo, nxw, nyw, stagi, stagj, zero, error)
int *nxw, *nyw;
int *error;
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
 *	nter	int *		Number of points in topo_ht (output)
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
    (*xter)[0] = pt1[0];
    (*xter)[1] = pt2[0];
    (*yter)[0] = pt1[1];
    (*yter)[1] = pt2[1];
    (*topo_ht)[0] = (*topo_ht)[1] = *zero;
}

/*
 * phys_2_lonlat : This routine converts physical coordinates to longitude-
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
    int i;

    for (i=0; i < *npts; ++i) {
	lon[i] = x[i];
	lat[i] = y[i];
    }
}

/*
 * lonlat_2_phys : This routine converts longitude-latitude coordinates
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
    int i;

    for (i=0; i < *npts; ++i) {
	x[i] = lon[i];
	y[i] = lat[i];
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
    /*
       The default is the current window in a Cylindrical Equidistant
       Projection
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
 * ri		float *	input	Multiply normal # x-points by this
 * rj		float *	input	Multiply normal # y-points by this
 *
 * The corners will be:     2    3
 *                      jj
 *                          0    1
 *                            ii
 */

float *
data_slicer_2d_(corner, ii, jj, ni, nj, da, nda, ri, rj)

float *corner, *ri, *rj;
int *ii, *jj, *ni, *nj, *da, *nda;
{
    /*
       Use pslicer2d_
       */
    float *pslicer2d_();

    return(pslicer2d_(corner, ii, jj, ni, nj, da, nda, ri, rj));
}

/*
 * data_slicer_3d_ : This routine returns the computational points at
 * which a 3d slice should be taken.
 *
 * Arguments:
 * corner	float *	input	2 corner points of slice in physical space.
 * ii		int *	input	axis coordinate of ni.
 * jj		int *	input	axis coordinate of nj.
 * kk		int *	input	axis coordinate of nk.
 * ni		int *	output	Number of points in the x-direction.
 * nj		int *	output	Number of points in the y-direction.
 * nk		int *	output	Number of points in the z-direction.
 * da		int *	input	axis coordinate of average direction
 * nda		int *	output	Number of points in the da-direction
 *
 */

float *
data_slicer_3d_(corner, ii, jj, kk, ni, nj, nk, da, nda)

float *corner;
int *ii, *jj, *kk, *ni, *nj, *nk, *da, *nda;
{
    /*
       Use pslicer3d_
       */
    float *pslicer3d_();

    return(pslicer3d_(corner, ii, jj, kk, ni, nj, nk, da, nda));
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
int *ndims, len1, len2, len3, len4, *dims;
float *stag, *min, *max, *missing;
char *name, *data_units, *data_display_units, *dim_names;
{
    float *derivative_();

    /*
       This is where you put the code to do your own calculation.

       Here is an example.  Say you have temperature and pressure fields,
       and you want to calculate the potential temperature:

       THETA = T(1000/P)**0.277

       There are two ways to code this.  This first uses the ability of
       the IVE routine "getvar" to do the math for you:

       if (strcmp(name, "THETA")) {
           int flag;

           return(getvar_("T*(1000/P)^0.277", ndims, dims, stag, min,
	                  max, missing, data_units, data_display_units,
			  dim_names, &flag, 16, len2, len3, len4);
       }
       
       The second way is to grab the variables and do it yourself
       (this code assumes T and P are on the same grid):

       if (strcmp(name, "THETA")) {
           float *t, *p, *result, pmissing;
	   int flag, i, num;

           if ((p = getvar_("P", ndims, dims, stag, min, max,
	                    pmissing, data_units, data_display_units,
			    dim_names, &flag, 1, len2, len3, len4)) == NULL) {
	       make_help_widget_("THETA - cannot get P");
	       return NULL;
	   }
           if ((t = getvar_("T", ndims, dims, stag, min, max,
	                    missing, data_units, data_display_units, dim_names,
			    &flag, 1, len2, len3, len4)) == NULL) {
	       make_help_widget_("THETA - cannot get T");
	       return NULL;
	   }
	   num = dims[0] * dims[1] * dims[2] * dims[3];
	   if ((result=(float *)memalign(sizeof(float), sizeof(float)*num))
	       == NULL) {
	       make_help_widget_("THETA - cannot allocate memory");
	       return NULL;
	   }
	   for (i=0; i < num; i++) {
	       if (t[i] == missing || p[i] == pmiss) result[i] = missing;
	       else result[i] =
	           pow((double)t[i]*(1000./p[i]), (double) 0.277);
	   return(result);
       }
       */
    /*
       "derivative" is a routine in IVE that calculates the derivative
       of a field, if the syntax of "name" is D[field:dir:type].
       If "name" is not in that format, or if derivative cannot do the
       calculation, it returns NULL.
       */

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
heading_(where, which, line, len)
     int  *where;
     int  *which;
     char *line;
     int len; /*stupid fortran length int*/
{
}

/*
 * file_coordinate_: This routine allows the user to change the current
 * value of the file coordinate from what IVE thinks it should be.
 *
 * Arguments:
 *	value	float *		Current value of the file coordinate
 *	coord	int *		Index of file coordinate (0-3)
 *      ncid    int *           Id of currently opened netCDF file.
 */

float
file_coordinate_(value, coord, ncid)
    float *value;
    int *coord, *ncid;
{
    return(*value);
}

