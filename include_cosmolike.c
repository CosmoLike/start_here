#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>


#include "../cosmolike_public/theory/basics.c"
#include "../cosmolike_public/theory/structs.c"
#include "../cosmolike_public/theory/parameters.c"
#include "../cosmolike_public/emu17/P_cb/emu.c"
#include "../cosmolike_public/theory/recompute.c"
#include "../cosmolike_public/theory/cosmo3D_public.c"
#include "../cosmolike_public/theory/redshift_spline.c"
#include "../cosmolike_public/theory/halo.c"
#include "../cosmolike_public/theory/HOD.c"
#include "../cosmolike_public/theory/cosmo2D_fourier.c"
#include "../cosmolike_public/theory/IA.c"
#include "../cosmolike_public/theory/cluster.c"
#include "../cosmolike_public/theory/BAO.c"
#include "../cosmolike_public/theory/external_prior.c"
#include "../cosmolike_public/theory/covariances_3D.c"
#include "../cosmolike_public/theory/covariances_fourier.c"
#include "../cosmolike_public/theory/init.c"


