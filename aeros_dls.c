/***********************************************************/
/* AEROS_DLS ... Aerosol generator with dls.run            */
/* Dubovik, O., et al. (2006), Application of spheroid     */
/* models to account for aerosol particle nonsphericity in */
/* remote sensing of desert dust,                          */
/* J.Geophys.Res., 111, D11208                             */
/* Author: N.Manago                                        */
/* $Revision: 1117 $                                       */
/* $Date: 2015-09-04 15:48:04 +0900 (Fri, 04 Sep 2015) $   */
/***********************************************************/
#define	M_1_LNA			0.43429448190325176		// 1/Ln(10)
#define	PI4_3_LNA		1.8191684717891214		// 4*PI/3/Ln(10)
#define	MIE_RMIN		NAN				// R min in um
#define	MIE_RMAX		NAN				// R max in um
#define	MIE_LMIN		-4.0				// Min Log10(R)
#define	MIE_LMAX		+3.0				// Max Log10(R)
#define	MIE_WSGM		5.0				// Log10(R) width in sigma
#define	DLS_ROOT		"/usr/local/DLS"		// Root directory
#define	DLS_N_GRID		1001				// #Grid radii
#define	DLS_N_AXIS		1				// #Axis ratios
#define	DLS_N_ANGL		181				// #Angles
#define	DLS_MAXMODE		2				// Max number of modes
#define	DLS_MAXGRID		1001				// Max number of grid radii
#define	DLS_MAXAXIS		25				// Max number of axis ratios
#define	DLS_MAXANGL		181				// Max number of scattering angles
#define	DLS_KEY			2				// Key
#define	DLS_RMIN		0.6480964e-03			// Min radius
#define	DLS_RMAX		0.3388173e+02			// Max radius
#define	DLS_XMIN		0.01197679346504701		// Min size parameter
#define	DLS_XMAX		626.1329062288993		// Max size parameter
#define	DLS_WMIN		0.34				// Min wavelength
#define	DLS_WMAX		0.34				// Max wavelength
#define	DELTA			1.0e-6				// A small number
#define	ROUND_D(x,n)		(round(x*1.0e##n)*1.0e-##n)
#define	ROUND_E(x,n)		(round(x*pow(10.0,n-floor(log10(x))))*pow(0.1,n-floor(log10(x))))

#include "aeros_common.c"

int WriteInput(int ngrd,double* grid,double* dist);
int ReadOutput(int ngrd,double* grid,double* dist);

// Parameters for Mie calculation
double	mie_rmin		= MIE_RMIN;			// R min in um
double	mie_rmax		= MIE_RMAX;			// R max in um
double	mie_wsgm		= MIE_WSGM;			// Log10(R) width in sigma
double	mie_rmin_com[MIE_MAXCOMP];				// R min in um
double	mie_rmax_com[MIE_MAXCOMP];				// R min in um
// Parameters for DLS
int	dls_key			= DLS_KEY;
int	dls_key_RD		= 2; // surface area mixture of spheroids
int	dls_key_f11		= 0; // calculate scattering matrix
int	dls_key_f344		= 0; // calculate scattering matrix
int	dls_key_org		= 0; // read original kernels
int	dls_key_fx		= 0; // create fixed kernels (if dls_key == 1)
int	dls_key_SD		= 1; // calculate Size Distribution using Log Normal function
int	dls_ID			= 0; // dimension of d(...)/dlnR or d(...)/dR (0 -> number)
int	dls_n_axis		= DLS_N_AXIS;
double	dls_axis[DLS_MAXAXIS];
double	dls_axis_dist[DLS_MAXAXIS];
int	dls_n_grid			= DLS_N_GRID;
double	dls_grid[DLS_MAXGRID];
double	dls_dist[DLS_MAXGRID];
int	dls_n_mode			= 1;			// Number of modes (up to 2)
double	dls_mode_wcon[DLS_MAXMODE];				// Concentration
double	dls_mode_rmed[DLS_MAXMODE];				// Median radius in um
double	dls_mode_lsgm[DLS_MAXMODE];				// Standard deviation
double	dls_rmin			= DLS_RMIN;		// Min radius
double	dls_rmax			= DLS_RMAX;		// Max radius
double	dls_wmin			= DLS_WMIN;		// Min wavelength
double	dls_wmax			= DLS_WMAX;		// Max wavelength
double	dls_wlen			= NAN;
double	dls_refr			= NAN;
double	dls_refi			= NAN;
double	dls_aext;						// Extinction coefficient
double	dls_asca;						// Scattering coefficient
double	dls_asym;						// Asymmetry parameter
double	dls_omeg;						// Single scattering albedo
double	dls_nsum;						// Number sum
double	dls_vcom;						// Mean volume
double	*dls_mtrx_f11			= NULL;			// Stokes matrix (F11)
double	*dls_mtrx_f22			= NULL;			// Stokes matrix (F22)
double	*dls_mtrx_f33			= NULL;			// Stokes matrix (F33)
double	*dls_mtrx_f44			= NULL;			// Stokes matrix (F44)
double	*dls_mtrx_f12			= NULL;			// Stokes matrix (F12)
double	*dls_mtrx_f34			= NULL;			// Stokes matrix (F34)
int	dls_n_angl			= DLS_N_ANGL;
double	dls_angl[DLS_MAXANGL] =
{
    0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
   10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
   20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
   30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0,
   40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0,
   50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0,
   60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0,
   70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0,
   80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0,
   90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0,
  100.0,101.0,102.0,103.0,104.0,105.0,106.0,107.0,108.0,109.0,
  110.0,111.0,112.0,113.0,114.0,115.0,116.0,117.0,118.0,119.0,
  120.0,121.0,122.0,123.0,124.0,125.0,126.0,127.0,128.0,129.0,
  130.0,131.0,132.0,133.0,134.0,135.0,136.0,137.0,138.0,139.0,
  140.0,141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,
  150.0,151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,
  160.0,161.0,162.0,163.0,164.0,165.0,166.0,167.0,168.0,169.0,
  170.0,171.0,172.0,173.0,174.0,175.0,176.0,177.0,178.0,179.0,
  180.0,
};
char	dls_root[MAXLINE]		= DLS_ROOT;		// Root directory
char	dls_finp[MAXLINE]		= "dls_input.dat";	// Input  file
char	dls_fout[MAXLINE]		= "sca_mtrx.dat";	// Output file
char	dls_frat[MAXLINE]		= "";			// Axis Ratio file
// Parameters for control
int	cnt_norun			= 0;			// No run  mode
int	cnt_vol				= 0;			// Volume equivalent mode
int	cnt_full			= 0;			// Calculate full kernel

int MieCalc(void)
{
  int i,j,k,n;
  int ngrd;
  int imin,imax;
  int err;
  double r,l;
  double lmin,lmax,lstp;
  double cnst,fact,indx;
  double fact1,fact2;
  double nsum;
  double *grid = NULL;
  double *dist = NULL;
  double *f11 = NULL;
  double *f12 = NULL;
  char fnam[] = "MieCalc";

  // Initialization
  dls_n_mode = 1;
  dls_mtrx_f11 = (double *)malloc(dls_n_angl*sizeof(double));
  dls_mtrx_f22 = (double *)malloc(dls_n_angl*sizeof(double));
  dls_mtrx_f33 = (double *)malloc(dls_n_angl*sizeof(double));
  dls_mtrx_f44 = (double *)malloc(dls_n_angl*sizeof(double));
  dls_mtrx_f12 = (double *)malloc(dls_n_angl*sizeof(double));
  dls_mtrx_f34 = (double *)malloc(dls_n_angl*sizeof(double));
  f11 = (double *)malloc(mie_n_angl*sizeof(double));
  f12 = (double *)malloc(mie_n_angl*sizeof(double));
  if(dls_mtrx_f11==NULL || dls_mtrx_f22==NULL || dls_mtrx_f33==NULL ||
     dls_mtrx_f44==NULL || dls_mtrx_f12==NULL || dls_mtrx_f34==NULL ||
     f11==NULL || f12==NULL)
  {
    fprintf(stderr,"%s: failed in allocating memory\n",fnam);
    return -1;
  }
  for(n=0; n<mie_n_comp; n++)
  {
    if(mie_size_func[n] == MIE_FUNC_FILE)
    {
      mie_rmin_com[n] = mie_xval_com[n][0];
      mie_rmax_com[n] = mie_xval_com[n][mie_n_size[n]-1];
    }
    else
    {
      if(!isnan(mie_rmin))
      {
        mie_rmin_com[n] = mie_rmin;
      }
      else
      {
        l = log10(mie_rmod_com[n])-mie_lsgm_com[n]*mie_wsgm;
        if(l < MIE_LMIN)
        {
          l = MIE_LMIN;
        }
        mie_rmin_com[n] = pow(10.0,l);
      }
      if(!isnan(mie_rmax))
      {
        mie_rmax_com[n] = mie_rmax;
      }
      else
      {
        l = log10(mie_rmod_com[n])+mie_lsgm_com[n]*mie_wsgm;
        if(l > MIE_LMAX)
        {
          l = MIE_LMAX;
        }
        mie_rmax_com[n] = pow(10.0,l);
      }
    }
  }
  if(cnt_full)
  {
    dls_wmin = DLS_WMIN;
    dls_wmax = DLS_WMAX;
    dls_rmin = DLS_RMIN;
    dls_rmax = DLS_RMAX;
  }
  else
  {
    dls_wmin = mie_wlen_um[0];
    dls_wmax = mie_wlen_um[0];
    for(i=1; i<mie_n_wlen; i++)
    {
      if(mie_wlen_um[i] < dls_wmin)
      {
        dls_wmin = mie_wlen_um[i];
      }
      if(mie_wlen_um[i] > dls_wmax)
      {
        dls_wmax = mie_wlen_um[i];
      }
    }
    dls_rmin = mie_rmin_com[0];
    dls_rmax = mie_rmax_com[0];
    for(n=1; n<mie_n_comp; n++)
    {
      if(mie_rmin_com[n] < dls_rmin)
      {
        dls_rmin = mie_rmin_com[n];
      }
      if(mie_rmax_com[n] > dls_rmax)
      {
        dls_rmax = mie_rmax_com[n];
      }
    }
    if(dls_rmin*PI2/dls_wmax < DLS_XMIN)
    {
      dls_rmin = DLS_XMIN*0.5*M_1_PI*dls_wmax;
    }
    if(dls_rmax*PI2/dls_wmin > DLS_XMAX)
    {
      dls_rmax = DLS_XMAX*0.5*M_1_PI*dls_wmin;
    }
  }
  err = 0;
  for(n=0; n<mie_n_comp; n++)
  {
    switch(mie_size_func[n])
    {
      case MIE_FUNC_FILE:
        dls_key_SD = 0;
        break;
      case MIE_FUNC_LOGNORMAL:
      case MIE_FUNC_EFF_LOGNORMAL:
        dls_mode_wcon[0] = 1.0;
        dls_mode_rmed[0] = mie_rmod_com[n];
        dls_mode_lsgm[0] = mie_lsgm_com[n]*M_LNA;
        ngrd = dls_n_grid;
        grid = dls_grid;
        dist = NULL;
        dls_key_SD = 1;
        break;
      case MIE_FUNC_GAMMA:
        indx = 1.0/mie_lsgm_com[n]-2.0;
        cnst = PI4_3/(pow(mie_rmod_com[n]*mie_lsgm_com[n],indx)*tgamma(indx));
        fact = -1.0/(mie_rmod_com[n]*mie_lsgm_com[n]);
        ngrd = dls_n_grid;
        grid = dls_grid;
        dist = dls_dist;
        dls_key_SD = 0;
        break;
      default:
        fprintf(stderr,"%s: error, mie_size_func[%d]=%c\n",fnam,n,mie_size_func[n]);
        err = 1;
        break;
    }
    if(err) break;
    if(mie_size_shape[n] == MIE_SHAPE_SPHEROID)
    {
      if(mie_n_reps[n] > DLS_MAXAXIS)
      {
        fprintf(stderr,"%s: error, number of axis exceed the limit (%d) >>> %d\n",
                 fnam,DLS_MAXAXIS,mie_n_reps[n]);
        err = 1;
        break;
      }
      dls_n_axis = mie_n_reps[n];
      for(i=0; i<mie_n_reps[n]; i++)
      {
        dls_axis[i] = mie_xeps_com[n][i];
        dls_axis_dist[i] = mie_yeps_com[n][i];
      }
    }
    else
    {
      fprintf(stderr,"%s: error, mie_size_shape[%d]=%c\n",fnam,n,mie_size_shape[n]);
      err = 1;
      break;
    }
    nsum = -1.0e100;
    for(i=0; i<mie_n_wlen; i++)
    {
      dls_wlen = mie_wlen_um[i];
      dls_refr = mie_refr_com[n][i];
      dls_refi = mie_refi_com[n][i];
      if(mie_size_func[n] == MIE_FUNC_FILE)
      {
        fact1 = PI2/mie_wlen_um[i];
        for(j=0; j<mie_n_size[n]; j++)
        {
          if(mie_xval_com[n][j]*fact1 > DLS_XMIN) break;
        }
        imin = j;
        if(imin >= mie_n_size[n])
        {
          fprintf(stderr,"%s: error, imin=%d, xmax=%13.6e, DLS_XMIN=%13.6e\n",
                          fnam,imin,mie_xval_com[n][mie_n_size[n]-1]*fact1,DLS_XMIN);
          err = 1;
          break;
        }
        for(j=mie_n_size[n]-1; j>=0; j--)
        {
          if(mie_xval_com[n][j]*fact1 < DLS_XMAX) break;
        }
        imax = j;
        if(imax < 0)
        {
          fprintf(stderr,"%s: error, imax=%d, xmin=%13.6e, DLS_XMAX=%13.6e\n",
                          fnam,imax,mie_xval_com[n][0]*fact1,DLS_XMAX);
          err = 1;
          break;
        }
        ngrd = imax-imin+1;
        grid = &mie_xval_com[n][imin];
        dist = &mie_yval_com[n][imin];
      }
      else
      {
        fact1 = PI2/mie_wlen_um[i];
        fact2 = 0.5*M_1_PI*mie_wlen_um[i];
        if(mie_rmin_com[n]*fact1 < DLS_XMIN)
        {
          lmin = log10(DLS_XMIN*fact2);
        }
        else
        {
          lmin = log10(mie_rmin_com[n]);
        }
        if(mie_rmax_com[n]*fact1 > DLS_XMAX)
        {
          lmax = log10(DLS_XMAX*fact2);
        }
        else
        {
          lmax = log10(mie_rmax_com[n]);
        }
        lstp = (lmax-lmin)/(dls_n_grid-1);
        for(j=0; j<dls_n_grid; j++)
        {
          dls_grid[j] = pow(10.0,lmin+lstp*j);
        }
        switch(mie_size_func[n])
        {
          case MIE_FUNC_LOGNORMAL:
          case MIE_FUNC_EFF_LOGNORMAL:
            break;
          case MIE_FUNC_GAMMA:
            for(j=0; j<dls_n_grid; j++)
            {
              r = dls_grid[j];
              dls_dist[j] = cnst*pow(r,indx+3.0)*exp(fact*r); // dV/dLn(R)
            }
            break;
          default:
            fprintf(stderr,"%s: error, mie_size_func[%d]=%c\n",fnam,n,mie_size_func[n]);
            err = 1;
            break;
        }
        if(err) break;
      }
      if(!cnt_norun)
      {
        if(WriteInput(ngrd,grid,dist) < 0)
        {
          err = 1;
          break;
        }
        system("dls.run");
      }
      if(ReadOutput(ngrd,grid,dist) < 0)
      {
        err = 1;
        break;
      }
      if(dls_key == 1)
      {
        dls_key = 2;
      }
      mie_aext_com[n][i] = dls_aext;
      mie_asca_com[n][i] = dls_asca;
      mie_asym_com[n][i] = dls_asym;
      if(dls_nsum > nsum)
      {
        nsum = dls_nsum;
        mie_vcom_com[n] = dls_vcom;
      }
      SamplingE(dls_n_angl,dls_angl,dls_mtrx_f11,mie_n_angl,mie_angl,f11,1.0);
      SamplingE(dls_n_angl,dls_angl,dls_mtrx_f12,mie_n_angl,mie_angl,f12,1.0);
      for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
      {
        mie_phs1_com[n][k] = f11[j]-f12[j]; // perpendicular polarization
        mie_phs2_com[n][k] = f11[j]+f12[j]; // parallel polarization
        mie_phas_com[n][k] = f11[j];
      }
    }
    if(err) break;
  }
  // deallocate memory
  free(dls_mtrx_f11);
  free(dls_mtrx_f22);
  free(dls_mtrx_f33);
  free(dls_mtrx_f44);
  free(dls_mtrx_f12);
  free(dls_mtrx_f34);
  free(f11);
  free(f12);
  if(err)
  {
    return -1;
  }

  if(Normalize(0) < 0) return -1;

  return 0;
}

int WriteInput(int ngrd,double* grid,double* dist)
{
  int i;
  int err;
  FILE *fp;

  if((fp=fopen(dls_finp,"w")) == NULL)
  {
    fprintf(stderr,"Error, cannot open %s\n",dls_finp);
    return -1;
  }
  err = 0;
  do
  {
    fprintf(fp,"%d %d %d %d %d %d\n",
                dls_key,dls_key_RD,dls_key_f11,
                dls_key_f344,dls_key_org,dls_key_fx);
    fprintf(fp,"%13.6f %22.15e %22.15e %14.10f %14.10f %13.6f %13.6f\n",
                dls_wlen,dls_refr,dls_refi,
                dls_rmin,dls_rmax,dls_wmin,dls_wmax);
    fprintf(fp,"%d %d %d\n",dls_key_SD,dls_ID,dls_n_mode);
    if(dls_key_SD == 0)
    {
      fprintf(fp,"blank\n");
      fprintf(fp,"%d\n",ngrd); // KN
      for(i=0; i<ngrd; i++)
      {
        fprintf(fp,"%13.6e %22.15e\n",grid[i],dist[i]);
      }
    }
    else
    {
      for(i=0; i<dls_n_mode; i++)
      {
        fprintf(fp,"%22.15e %22.15e %22.15e%s",
                    dls_mode_wcon[i],dls_mode_lsgm[i],dls_mode_rmed[i],
                    (i==dls_n_mode-1?"\n":" "));
      }
      fprintf(fp,"%d\n",ngrd); // KN
      for(i=0; i<ngrd; i++)
      {
        fprintf(fp,"%13.6e\n",grid[i]);
      }
    }
    fprintf(fp,"KERNEL_n22_181\n");
    fprintf(fp,"KERNEL_n22_fix\n");
    fprintf(fp,"KERNEL_n22_new\n");
    fprintf(fp,"%d\n",dls_n_axis); // KR
    for(i=0; i<dls_n_axis; i++)
    {
      fprintf(fp,"%13.6f %22.15e\n",dls_axis[i],dls_axis_dist[i]);
    }
    fprintf(fp,"%d\n",dls_n_angl); // KM
    for(i=0; i<dls_n_angl; i++)
    {
      fprintf(fp,"%13.6f\n",dls_angl[i]); // ANGLE(KM)
    }
  }
  while(0);
  fclose(fp);
  if(err)
  {
    return -1;
  }

  return 0;
}

int ReadOutput(int ngrd,double* grid,double* dist)
{
  int i,j;
  int err;
  double v[8];
  double *l = NULL;
  double *dl = NULL;
  double *vinv = NULL;
  double dv,sn,sv;
  char line[MAXLINE];
  char str[8][MAXLINE];
  char *p;
  char fnam[] = "ReadOutput";
  FILE *fp;

  // Output example
  // Size distribution:
  //       r(mkm)             Sd(r)
  //     0.50000E-01     0.36367E-02
  //     0.65604E-01     0.12765E-01
  //     ...
  //
  // Axis ratio distribution:
  //          R                Rd(R)
  //     0.33490E+00     0.66185E-01
  //     0.36690E+00     0.65025E-01
  //     ...
  //
  //wavelength= 0.4400E+00  n= 0.1560E+01  k= 0.2900E-02
  //
  //         APPROXIMATED OPTICAL CHARACTERISTICS
  //                  (volume mixture)
  //
  //ext=  0.72684E+00  abs=  0.49962E-01  sca=  0.67688E+00  albedo=  0.93126E+00
  //
  //  ANGLE       F11        -F12/F11       F22/F11       F33/F11       F34/F11       F44/F11
  //   0.00   0.38325E+03  -0.00000E+00   0.99984E+00   0.99984E+00   0.00000E+00   0.99969E+00
  //   1.00   0.27348E+03   0.83746E-04   0.99978E+00   0.99978E+00  -0.39201E-03   0.99956E+00
  //   ...
  // 180.00   0.25674E+00  -0.00000E+00   0.71376E+00  -0.71240E+00   0.00000E+00  -0.44569E+00
  //asymmetry parameter=   0.6900E+00

  if(ngrd < 2)
  {
    fprintf(stderr,"%s: error, ngrd=%d\n",fnam,ngrd);
    return -1;
  }
  l = (double *)malloc(ngrd*sizeof(double));
  dl = (double *)malloc(ngrd*sizeof(double));
  vinv = (double *)malloc(ngrd*sizeof(double));
  if(l==NULL || dl==NULL || vinv==NULL)
  {
    fprintf(stderr,"%s: failed in allocating memory\n",fnam);
    return -1;
  }
  for(i=0; i<ngrd; i++)
  {
    l[i] = log(grid[i]);
    vinv[i] = 1.0/(PI4_3*grid[i]*grid[i]*grid[i]);
  }
  dl[0] = l[1]-l[0];
  for(i=1; i<ngrd-1; i++)
  {
    dl[i] = 0.5*(l[i+1]-l[i-1]);
  }
  dl[ngrd-1] = l[ngrd-1]-l[ngrd-2];

  if((fp=fopen(dls_fout,"r")) == NULL)
  {
    fprintf(stderr,"%s: cannot open %s\n",fnam,dls_fout);
    return -1;
  }
  err = 0;
  do
  {
    // Read size distribution
    err = 1;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      for(p=line; *p==' '; p++);
      if(strncasecmp("Size distribution:",p,18) == 0)
      {
        err = 0;
        break;
      }
    }
    if(err)
    {
      fprintf(stderr,"%s: failed in reading size distribution.\n",fnam);
      break;
    }
    if(fgets(line,MAXLINE,fp) == NULL) // skip one line
    {
      fprintf(stderr,"%s: error in reading size distribution.\n",fnam);
      err = 1;
      break;
    }
    sn = 0.0;
    sv = 0.0;
    for(i=0; i<ngrd; i++)
    {
      if(fgets(line,MAXLINE,fp) == NULL)
      {
        fprintf(stderr,"%s: error in reading size distribution, i=%d\n",fnam,i);
        err = 1;
        break;
      }
      if(sscanf(line,"%s%s",str[0],str[1]) != 2)
      {
        fprintf(stderr,"%s: error in reading size distribution, i=%d >>> %s\n",fnam,i,line);
        err = 1;
        break;
      }
      for(j=0; j<2; j++)
      {
        errno = 0;
        v[j] = strtod(str[j],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(err) break;
      if(fabs(v[0]/grid[i]-1.0) > DELTA)
      {
        fprintf(stderr,"%s: error, grid[%d]=%13.6e, v[0]=%13.6e (%13.6e)\n",fnam,i,grid[i],v[0],v[0]/grid[i]-1.0);
        err = 1;
        break;
      }
      if(dls_key_SD == 0)
      {
        if(fabs(dist[i]) < 1.0e-30)
        {
          if(fabs(v[1]-dist[i]) > 1.0e-30)
          {
            fprintf(stderr,"%s: error, dist[%d]=%13.6e, v[1]=%13.6e (%13.6e)\n",fnam,i,dist[i],v[1],v[1]-dist[i]);
            err = 1;
            break;
          }
        }
        else
        {
          if(fabs(v[1]/dist[i]-1.0) > DELTA)
          {
            fprintf(stderr,"%s: error, dist[%d]=%13.6e, v[1]=%13.6e (%13.6e)\n",fnam,i,dist[i],v[1],v[1]/dist[i]-1.0);
            err = 1;
            break;
          }
        }
        dv = dist[i]*dl[i];
        sn += vinv[i]*dv;
        sv += dv;
      }
      else
      {
        dv = v[1]*dl[i];
        sn += vinv[i]*dv;
        sv += dv;
      }
    }
    if(err) break;
    if(fabs(sn-1.0) > 1.0e-1)
    {
      fprintf(stderr,"%s: warning, sn=%13.6e\n",fnam,sn);
    }
    dls_nsum = sn;
    dls_vcom = sv/sn;
    // Read axis ratio distribution
    err = 1;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      for(p=line; *p==' '; p++);
      if(strncasecmp("Axis ratio distribution:",p,24) == 0)
      {
        err = 0;
        break;
      }
    }
    if(err)
    {
      fprintf(stderr,"%s: failed in reading axis ratio distribution.\n",fnam);
      break;
    }
    if(fgets(line,MAXLINE,fp) == NULL) // skip one line
    {
      fprintf(stderr,"%s: error in reading axis ratio distribution.\n",fnam);
      err = 1;
      break;
    }
    for(i=0; i<dls_n_axis; i++)
    {
      if(fgets(line,MAXLINE,fp) == NULL)
      {
        fprintf(stderr,"%s: error in reading axis ratio distribution, i=%d\n",fnam,i);
        err = 1;
        break;
      }
      if(sscanf(line,"%s%s",str[0],str[1]) != 2)
      {
        fprintf(stderr,"%s: error in reading axis ratio distribution, i=%d >>> %s\n",fnam,i,line);
        err = 1;
        break;
      }
      for(j=0; j<2; j++)
      {
        errno = 0;
        v[j] = strtod(str[j],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(err) break;
      if(fabs(v[0]/dls_axis[i]-1.0) > DELTA)
      {
        fprintf(stderr,"%s: error, dls_axis[%d]=%13.6e, v[0]=%13.6e (%13.6e)\n",fnam,i,dls_axis[i],v[0],v[0]/dls_axis[i]-1.0);
        err = 1;
        break;
      }
    }
    if(err) break;
    // Read wavelength etc.
    err = 1;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      for(p=line; *p==' '; p++);
      if(strncasecmp("wavelength=",p,11) == 0)
      {
        err = 0;
        break;
      }
    }
    if(err)
    {
      fprintf(stderr,"%s: failed in reading wavelength etc.\n",fnam);
      break;
    }
    if(sscanf(line,"%s%s%s%s%s%s",str[0],str[1],str[2],str[3],str[4],str[5]) != 6)
    {
      fprintf(stderr,"%s: error in reading wavelength etc >>> %s\n",fnam,line);
      err = 1;
      break;
    }
    if(strcasecmp(str[0],"wavelength=") !=0 ||
       strcasecmp(str[2],"n=") !=0 ||
       strcasecmp(str[4],"k=") !=0)
    {
      fprintf(stderr,"%s: error in reading keywords for wavelength etc, i=%d >>> %s\n",fnam,i,line);
      err = 1;
      break;
    }
    for(j=1; j<6; j+=2)
    {
      errno = 0;
      v[j] = strtod(str[j],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
        err = 1;
        break;
      }
    }
    if(err) break;
    if(fabs(v[1]/dls_wlen-1.0) > DELTA)
    {
      fprintf(stderr,"%s: error, dls_wlen=%13.6e, v[1]=%13.6e (%13.6e)\n",fnam,dls_wlen,v[1],v[1]/dls_wlen-1.0);
      err = 1;
      break;
    }
    if(fabs(v[3]/dls_refr-1.0) > 1.0e-2)
    {
      fprintf(stderr,"%s: error, dls_refr=%13.6e, v[3]=%13.6e (%13.6e)\n",fnam,dls_refr,v[3],v[3]/dls_refr-1.0);
      err = 1;
      break;
    }
    if(fabs(dls_refi/v[5]-1.0) > 1.0e-2)
    {
      fprintf(stderr,"%s: warning, dls_refi=%13.6e, v[5]=%13.6e (%13.6e)\n",fnam,dls_refi,v[5],dls_refi/v[5]-1.0);
    }
    // Read approximated optical characteristics
    err = 1;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      for(p=line; *p==' '; p++);
      if(strncasecmp("ext=",p,4) == 0)
      {
        err = 0;
        break;
      }
    }
    if(err)
    {
      fprintf(stderr,"%s: failed in reading approximated optical characteristics.\n",fnam);
      break;
    }
    if(sscanf(line,"%s%s%s%s%s%s%s%s",str[0],str[1],str[2],str[3],str[4],str[5],str[6],str[7]) != 8)
    {
      fprintf(stderr,"%s: error in reading approximated optical characteristics >>> %s\n",fnam,line);
      err = 1;
      break;
    }
    if(strcasecmp(str[0],"ext=")   !=0 ||
       strcasecmp(str[2],"abs=")   !=0 ||
       strcasecmp(str[4],"sca=")   !=0 ||
       strcasecmp(str[6],"albedo=")!=0)
    {
      fprintf(stderr,"%s: error in reading keywords for approximated optical characteristics, i=%d >>> %s\n",fnam,i,line);
      err = 1;
      break;
    }
    for(j=1; j<8; j+=2)
    {
      errno = 0;
      v[j] = strtod(str[j],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
        err = 1;
        break;
      }
    }
    if(err) break;
    dls_aext = v[1];
    dls_asca = v[5];
    dls_omeg = v[7];
    // Read matrix elements
    err = 1;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      for(p=line; *p==' '; p++);
      if(strncasecmp("ANGLE",p,5) == 0)
      {
        err = 0;
        break;
      }
    }
    if(err)
    {
      fprintf(stderr,"%s: failed in reading matrix elements.\n",fnam);
      break;
    }
    for(i=0; i<dls_n_angl; i++)
    {
      if(fgets(line,MAXLINE,fp) == NULL)
      {
        fprintf(stderr,"%s: error in reading matrix elements, i=%d\n",fnam,i);
        err = 1;
        break;
      }
      if(sscanf(line,"%s%s%s%s%s%s%s",str[0],str[1],str[2],str[3],str[4],str[5],str[6]) != 7)
      {
        fprintf(stderr,"%s: error in reading matrix elements, i=%d >>> %s\n",fnam,i,line);
        err = 1;
        break;
      }
      for(j=0; j<7; j++)
      {
        errno = 0;
        v[j] = strtod(str[j],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(err) break;
      if(fabs(v[0]/dls_angl[i]-1.0) > DELTA)
      {
        fprintf(stderr,"%s: error, dls_angl[%d]=%13.6e, v[0]=%13.6e (%13.6e)\n",fnam,i,dls_angl[i],v[0],v[0]/dls_angl[i]-1.0);
        err = 1;
        break;
      }
      dls_mtrx_f11[i] = v[1];
      dls_mtrx_f12[i] = v[1]*v[2]*(-1.0);
      dls_mtrx_f22[i] = v[1]*v[3];
      dls_mtrx_f33[i] = v[1]*v[4];
      dls_mtrx_f34[i] = v[1]*v[5];
      dls_mtrx_f44[i] = v[1]*v[6];
    }
    if(err) break;
    // Read asymmetry parameter
    err = 1;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      for(p=line; *p==' '; p++);
      if(strncasecmp("asymmetry parameter=",p,20) == 0)
      {
        err = 0;
        break;
      }
    }
    if(err)
    {
      fprintf(stderr,"%s: failed in reading asymmetry parameter.\n",fnam);
      break;
    }
    if(sscanf(line,"%s%s%s",str[0],str[1],str[2]) != 3)
    {
      fprintf(stderr,"%s: error in reading asymmetry parameter >>> %s\n",fnam,line);
      err = 1;
      break;
    }
    errno = 0;
    v[2] = strtod(str[2],&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
      err = 1;
      break;
    }
    dls_asym = v[2];
  }
  while(0);
  fclose(fp);
  free(l);
  free(dl);
  free(vinv);
  if(err)
  {
    return -1;
  }

  return 0;
}

int Init(void)
{
  int i,n;
  char line[MAXLINE];
  double fact;
  const char fnam[] = "Init";

  if(dls_key == 1)
  {
    if(access("KERNEL_n22_181",R_OK) < 0)
    {
      snprintf(line,MAXLINE,"%s/KERNEL_n22_181",dls_root);
      if(symlink(line,"KERNEL_n22_181") < 0)
      {
        fprintf(stderr,"%s: error, cannot create symlink to %s\n",fnam,line);
        return -1;
      }
    }
    if(access("KERNEL_n22_fix",R_OK) < 0)
    {
      if(mkdir("KERNEL_n22_fix",0755) < 0)
      {
        fprintf(stderr,"%s: error, cannot make directory %s\n",fnam,"KERNEL_n22_fix");
        return -1;
      }
    }
    if(access("name.dat",R_OK) < 0)
    {
      snprintf(line,MAXLINE,"%s/ROUTINE/name.dat",dls_root);
      if(symlink(line,"name.dat") < 0)
      {
        fprintf(stderr,"%s: error, cannot create symlink to %s\n",fnam,line);
        return -1;
      }
    }
  } else
  if(dls_key == 2)
  {
    if(access("KERNEL_n22_fix",R_OK) < 0)
    {
      snprintf(line,MAXLINE,"%s/KERNEL_n22_fix",dls_root);
      if(symlink(line,"KERNEL_n22_fix") < 0)
      {
        fprintf(stderr,"%s: error, cannot create symlink to %s\n",fnam,line);
        return -1;
      }
    }
  } else
  if(dls_key == 4)
  {
    if(access("KERNEL_n22_181",R_OK) < 0)
    {
      snprintf(line,MAXLINE,"%s/KERNEL_n22_181",dls_root);
      if(symlink(line,"KERNEL_n22_181") < 0)
      {
        fprintf(stderr,"%s: error, cannot create symlink to %s\n",fnam,line);
        return -1;
      }
    }
    if(access("name.dat",R_OK) < 0)
    {
      snprintf(line,MAXLINE,"%s/ROUTINE/name.dat",dls_root);
      if(symlink(line,"name.dat") < 0)
      {
        fprintf(stderr,"%s: error, cannot create symlink to %s\n",fnam,line);
        return -1;
      }
    }
  }
  if(dls_key_org == 1)
  {
    if(access("KERNEL_n22_new",R_OK) < 0)
    {
      if(mkdir("KERNEL_n22_new",0755) < 0)
      {
        fprintf(stderr,"%s: error, cannot make directory %s\n",fnam,"KERNEL_n22_new");
        return -1;
      }
    }
  }
  if(cnt_vol)
  {
    dls_key_RD = 1;
  }
  else
  {
    dls_key_RD = 2;
  }

  if(CommonInit() < 0)
  {
    return -1;
  }

  for(n=0; n<mie_n_comp; n++)
  {
    if(mie_size_func[n] == MIE_FUNC_FILE)
    {
      switch(mie_size_xtype[n])
      {
        case MIE_XTYPE_R:
          switch(mie_size_ytype[n])
          {
            case MIE_YTYPE_NUMBER:
              for(i=0; i<mie_n_size[n]; i++)
              {
                fact = mie_xval_com[n][i]; // r^1
                fact *= fact; // r^2
                fact *= fact; // r^4
                mie_yval_com[n][i] *= PI4_3*fact; // 4*PI/3*r^4 -> dV/dLn(R)
              }
              break;
            case MIE_YTYPE_VOLUME:
              for(i=0; i<mie_n_size[n]; i++)
              {
                fact = mie_xval_com[n][i]; // r^1
                mie_yval_com[n][i] *= fact; // r^1 -> dV/dLn(R)
              }
              break;
            default:
              fprintf(stderr,"%s: error, mie_size_ytype[%d]=%d\n",fnam,n,mie_size_ytype[n]);
              return -1;
              break;
          }
          break;
        case MIE_XTYPE_LOGR:
          switch(mie_size_ytype[n])
          {
            case MIE_YTYPE_NUMBER:
              for(i=0; i<mie_n_size[n]; i++)
              {
                mie_xval_com[n][i] = pow(10.0,mie_xval_com[n][i]); // Log10(R) -> R
                fact = mie_xval_com[n][i]; // r^1
                fact *= mie_xval_com[n][i]; // r^2
                fact *= mie_xval_com[n][i]; // r^3
                mie_yval_com[n][i] *= PI4_3_LNA*fact; // 4*PI/3/Ln(10)*r^3 -> dV/dLn(R)
              }
              break;
            case MIE_YTYPE_VOLUME:
              for(i=0; i<mie_n_size[n]; i++)
              {
                mie_xval_com[n][i] = pow(10.0,mie_xval_com[n][i]); // Log10(R) -> R
                fact = M_1_LNA;
                mie_yval_com[n][i] *= fact; // 1/Ln(10) -> dV/dLn(R)
              }
              break;
            default:
              fprintf(stderr,"%s: error, mie_size_ytype[%d]=%d\n",fnam,n,mie_size_ytype[n]);
              return -1;
              break;
          }
          break;
        case MIE_XTYPE_LNR:
          switch(mie_size_ytype[n])
          {
            case MIE_YTYPE_NUMBER:
              for(i=0; i<mie_n_size[n]; i++)
              {
                mie_xval_com[n][i] = exp(mie_xval_com[n][i]); // Ln(R) -> R
                fact = mie_xval_com[n][i]; // r^1
                fact *= mie_xval_com[n][i]; // r^2
                fact *= mie_xval_com[n][i]; // r^3
                mie_yval_com[n][i] *= PI4_3*fact; // 4*PI/3*r^3 -> dV/dLn(R)
              }
              break;
            case MIE_YTYPE_VOLUME:
              for(i=0; i<mie_n_size[n]; i++)
              {
                mie_xval_com[n][i] = exp(mie_xval_com[n][i]); // Ln(R) -> R
              }
              break;
            default:
              fprintf(stderr,"%s: error, mie_size_ytype[%d]=%d\n",fnam,n,mie_size_ytype[n]);
              return -1;
              break;
          }
          break;
        default:
          fprintf(stderr,"%s: error, mie_size_xtype[%d]=%d\n",fnam,n,mie_size_xtype[n]);
          return -1;
          break;
      }
    }
  }

  return 0;
}

void Finish(void)
{
  CommonFinish();
}

int GetOwnOpt(int argn,char **args)
{
  int c,rt;
  int option_index = 0;
  int this_option_optind;
  int ntmp;
  char *endp;
  double xtmp;
  struct option long_options[] =
  {
    {"mie_rmin",1,0,'x'},
    {"mie_rmax",1,0,'X'},
    {"mie_wsgm",1,0,'S'},
    {"dls_root",1,0,'D'},
    {"dls_n_grid",1,0,'n'},
    {"dls_key",1,0,'K'},
    {"cnt_norun",0,0,'N'},
    {"cnt_vol",0,0,'V'},
    {"cnt_full",0,0,'F'},
  };

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-x:X:S:D:n:K:NVF",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'x':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_rmin = xtmp;
        else
        {
          fprintf(stderr,"R min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'X':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_rmax = xtmp;
        else
        {
          fprintf(stderr,"R max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'S':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_wsgm = xtmp;
        else
        {
          fprintf(stderr,"Log10(R) width -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'D':
        strncpy(dls_root,optarg,MAXLINE);
        break;
      case 'n':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0 && ntmp<=DLS_MAXGRID) dls_n_grid = ntmp;
        else
        {
          fprintf(stderr,"#Grids -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'K':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=1 && ntmp<=4) dls_key = ntmp;
        else
        {
          fprintf(stderr,"Key -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'N':
        cnt_norun = 1;
        break;
      case 'V':
        cnt_vol = 1;
        break;
      case 'F':
        cnt_full = 1;
        break;
      case '?':
        if(cnt_optn >= MAXNOPT)
        {
          fprintf(stderr,"Error, #options exceed the limit >>> %d\n",cnt_optn);
          rt = -1;
        }
        else
        {
          if(optopt == '\0')
          {
            cnt_opts[cnt_optn] = args[this_option_optind];
            cnt_optn++;
          }
          else
          {
            if((cnt_opts[cnt_optn]=(char*)malloc(3)) == NULL)
            {
              fprintf(stderr,"Error, cannot allocate memory.");
              rt = -1;
            }
            else
            {
              snprintf(cnt_opts[cnt_optn],3,"-%c",optopt);
              cnt_optn++;
            }
          }
        }
        break;
      case 1:
        if(cnt_optn >= MAXNOPT)
        {
          fprintf(stderr,"Error, #options exceed the limit >>> %d\n",cnt_optn);
          rt = -1;
        }
        else
        {
          cnt_opts[cnt_optn] = optarg;
          cnt_optn++;
        }
        break;
      case ':':
        if(optopt == '\0')
        {
          fprintf(stderr,"Option requires an argument : %s\n",args[this_option_optind]);
        }
        else
        {
          fprintf(stderr,"Option requires an argument : %c\n",optopt);
        }
        rt = -1;
        break;
      default:
        fprintf(stderr,"?? getopt returned character code 0%o ??\n",c);
        rt = -1;
        break;
    }
  }
  optind = 1;

  if(cnt_hp) rt = 1;
  return rt;
}

int Usage(void)
{
  int  i;
  int  n = 16;
  char e[MAXLINE];
  char a[MAXLINE];
  char d[MAXLINE];
  char optstring[] = "fgwlLaWAPQiIruxXSDnKNVFmoOdvh";

  CommonUsage("aeros_dls",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'x':
        fprintf(stderr," x -mie_rmin      |%s|%s|%s| %e\n",As(e,"R Min",n),          As(a,"um",n),         Ae(d,MIE_RMIN,n),mie_rmin);
        break;
      case 'X':
        fprintf(stderr," X -mie_rmax      |%s|%s|%s| %e\n",As(e,"R Max",n),          As(a,"um",n),         Ae(d,MIE_RMAX,n),mie_rmax);
        break;
      case 'S':
        fprintf(stderr," S -mie_wsgm      |%s|%s|%s| %f\n",As(e,"Log10(R) width",n), As(a,"sigma",n),      Af(d,MIE_WSGM,n),mie_wsgm);
        break;
      case 'D':
        fprintf(stderr," D -dls_root      |%s|%s|%s| %s\n",As(e,"Root directory",n), As(a,"directory",n),  As(d,DLS_ROOT,n),dls_root);
        break;
      case 'n':
        fprintf(stderr," n -dls_n_grid    |%s|%s|%s| %d\n",As(e,"#Grids",n),         As(a,"#",n),          Ad(d,DLS_N_GRID,n),dls_n_grid);
        break;
      case 'K':
        fprintf(stderr," K -dls_key       |%s|%s|%s| %d\n",As(e,"Key",n),            As(a,"1-4",n),        Ad(d,DLS_KEY,n),dls_key);
        break;
      case 'N':
        fprintf(stderr," N -cnt_norun     |%s|%s|%s| %d\n",As(e,"No run  mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_norun);
        break;
      case 'V':
        fprintf(stderr," V -cnt_vol       |%s|%s|%s| %d\n",As(e,"Volume eq. mode",n),As(a,"nothing",n),    Ad(d,0,n),cnt_vol);
        break;
      case 'F':
        fprintf(stderr," F -cnt_full      |%s|%s|%s| %d\n",As(e,"Full kernel",n),    As(a,"nothing",n),    Ad(d,0,n),cnt_full);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);
  fprintf(stderr,"Key of dls.run\n");
  fprintf(stderr,"1 ... create fixed kernels (for fixed axis\n");
  fprintf(stderr,"      ratio distr.) and save them\n");
  fprintf(stderr,"      into 'Rke...fix...' files and calculate\n");
  fprintf(stderr,"      opt.characteristics\n");
  fprintf(stderr,"2 ... read fixed kernels from 'Rke...fix...'  files\n");
  fprintf(stderr,"3 ... create fixed kernels but don't save them\n");
  fprintf(stderr,"4 ... don't create fixed kernels, calculate\n");
  fprintf(stderr,"      opt.characteristics from original kernels\n");
  CommonUsage(NULL,2);
  CommonUsage(NULL,3);
  CommonUsage(NULL,4);

  return 0;
}
