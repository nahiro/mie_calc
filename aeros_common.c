/***********************************************************/
/* AEROS_COMMON ... Aerosol generator common file          */
/* Author: N.Manago                                        */
/* $Revision: 1108 $                                        */
/* $Date: 2015-09-02 13:31:28 +0900 (Wed, 02 Sep 2015) $   */
/***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <bits/nan.h>
#include "strutil.h"

#define	PI			3.141592653589793		// PI
#define	PI2			6.283185307179586		// 2*PI
#define	PI4_3			4.1887902047863905		// 4*PI/3
#define	M_LNA			2.3025850929940459		// Ln(10)
#define	M_1_LNA			0.43429448190325176		// 1/Ln(10)
#define	R_TO_D			57.29577951308232		// 180/PI rad -> deg
#define	D_TO_R			1.745329251994329e-02		// PI/180 deg -> rad
#define	EPSILON			1.0e-14				// A small number
#define	MIE_FLAG_SIZE		4
#define	MIE_FUNC_INPUT		'I'
#define	MIE_FUNC_FILE		'F'
#define	MIE_FUNC_LOGNORMAL	'L'
#define	MIE_FUNC_EFF_LOGNORMAL	'E'
#define	MIE_FUNC_GAMMA		'G'
#define	MIE_XTYPE_R		'R'
#define	MIE_XTYPE_LOGR		'L'
#define	MIE_XTYPE_LNR		'N'
#define	MIE_YTYPE_NUMBER	'N'
#define	MIE_YTYPE_VOLUME	'V'
#define	MIE_SHAPE_SPHERE	'S'
#define	MIE_SHAPE_SPHEROID	'T'
#define	MIE_N_WLEN		11				// #Wavelengths
#define	MIE_N_ANGL		121				// #Angles
#define	MIE_MAXDATA		5000				// Max #data
#define	MIE_DMAX		1000000				// Max #data for size dist
#define	MIE_IMIN		0				// Min line#
#define	MIE_IMAX		1000000000			// Min line#
#define	MIE_WLEN_NUM		0				// Out. wave  column#
#define	MIE_ANGL_NUM		0				// Out. angle column#
#define	MIE_REAL_NUM		1				// Refractive index (real) column#
#define	MIE_IMAG_NUM		2				// Refractive index (imag) column#
#define	MIE_WLEN_UNI		1.0				// Wavelength unit in nm
#define	MIE_WLEN_REF		550.0				// Reference wavelength in nm
#define	MIE_WLEN_MIN		250.0				// Min wavelength in nm
#define	MIE_WLEN_MAX		1500.0				// Max wavelength in nm
#define	MIE_ANGL_UNI		1.0				// Angle unit in degree
#define	MIE_ANGL_MIN		0.0				// Min angle in degree
#define	MIE_ANGL_MAX		180.0				// Max angle in degree
#define	MIE_MAXCOMP		10				// Max #components
#define	MIE_OUT0		"mie_out0.dat"			// Out. file 0
#define	MIE_OUT1		"mie_out1.dat"			// Out. file 1
#define	MIE_OUT2		"mie_out2.dat"			// Out. file 2
#define	MIE_WLEN_NAM		""				// Wavelen file
#define	MIE_ANGL_NAM		""				// Angle file
#define	NONAME			"NONE"				// No name
#define	MAXITEM			20				// Max #entries in a line
#define	MAXNOPT			100				// Max #options
#define	MAXLINE			256				// Max #chars in a line
#define	MIN(a,b)		((b)<(a)?(b):(a))
#define	MAX(a,b)		((b)>(a)?(b):(a))

int	mie_imin		= MIE_IMIN;			// Min line#
int	mie_imax		= MIE_IMAX;			// Max line#
int	mie_iref		= -1;				// Ref wave line#
int	mie_wlen_imin		= MIE_IMAX;			// Min wave line#
int	mie_wlen_imax		= MIE_IMIN;			// Max wave line#
int	mie_wlen_num		= MIE_WLEN_NUM;			// Out. wave  column#
int	mie_angl_num		= MIE_ANGL_NUM;			// Out. angle column#
int	mie_real_num		= MIE_REAL_NUM;			// Refractive index (real) column#
int	mie_imag_num		= MIE_IMAG_NUM;			// Refractive index (imag) column#
int	mie_n_wlen		= MIE_N_WLEN;			// #Wavelengths (output)
int	mie_n_angl		= MIE_N_ANGL;			// #Angles
int	mie_n_comp		= 0;				// #Components
int	mie_n_size[MIE_MAXCOMP];				// #Size parameters
int	mie_n_reps[MIE_MAXCOMP];				// #Aspect ratios
char	mie_size_func[MIE_MAXCOMP];				// Size distribution function
char	mie_size_xtype[MIE_MAXCOMP];				// Size distribution X type
char	mie_size_ytype[MIE_MAXCOMP];				// Size distribution Y type
char	mie_size_shape[MIE_MAXCOMP];				// Particle shape
double	mie_wlen_uni		= MIE_WLEN_UNI;			// Wavelength unit
double	mie_wlen_ref		= MIE_WLEN_REF;			// Reference wavelength
double	mie_wlen_min		= MIE_WLEN_MIN;			// Min wavelength in nm
double	mie_wlen_max		= MIE_WLEN_MAX;			// Max wavelength in nm
double	mie_wlen[MIE_MAXDATA]	=
{
  300.0, 337.1, 400.0,  488.0,  514.5, 550.0,
  632.8, 694.3, 860.0, 1060.0, 1300.0,
};
double	*mie_wlen_um		= NULL;
double	*mie_aext		= NULL;
double	*mie_asca		= NULL;
double	*mie_asym		= NULL;
double	mie_angl_uni		= MIE_ANGL_UNI;			// Angle unit in degree
double	mie_angl_min		= MIE_ANGL_MIN;			// Min angle in degree
double	mie_angl_max		= MIE_ANGL_MAX;			// Max angle in degree
double	mie_angl[MIE_MAXDATA]	=
{
    0.0,   0.5,   1.0,   1.5,   2.0,   2.5,   3.0,   3.5,   4.0,   4.5,   5.0,
    5.5,   6.0,   6.5,   7.0,   7.5,   8.0,   8.5,   9.0,   9.5,  10.0,  12.0,
   14.0,  16.0,  18.0,  20.0,  22.0,  24.0,  26.0,  28.0,  30.0,  32.0,  34.0,
   36.0,  38.0,  40.0,  42.0,  44.0,  46.0,  48.0,  50.0,  52.0,  54.0,  56.0,
   58.0,  60.0,  62.0,  64.0,  66.0,  68.0,  70.0,  72.0,  74.0,  76.0,  78.0,
   80.0,  82.0,  84.0,  86.0,  88.0,  90.0,  92.0,  94.0,  96.0,  98.0, 100.0,
  102.0, 104.0, 106.0, 108.0, 110.0, 112.0, 114.0, 116.0, 118.0, 120.0, 122.0,
  124.0, 126.0, 128.0, 130.0, 132.0, 134.0, 136.0, 138.0, 140.0, 142.0, 144.0,
  146.0, 148.0, 150.0, 152.0, 154.0, 156.0, 158.0, 160.0, 162.0, 164.0, 166.0,
  168.0, 170.0, 170.5, 171.0, 171.5, 172.0, 172.5, 173.0, 173.5, 174.0, 174.5,
  175.0, 175.5, 176.0, 176.5, 177.0, 177.5, 178.0, 178.5, 179.0, 179.5, 180.0,
};
double	*mie_angl_rad		= NULL;
double	*mie_angl_sin		= NULL;
double	*mie_angl_cos		= NULL;
double	*mie_angl_dif		= NULL;
double	*mie_phs1		= NULL;
double	*mie_phs2		= NULL;
double	*mie_phas		= NULL;
double	*mie_refr_com[MIE_MAXCOMP];				// Refractive index (real)
double	*mie_refi_com[MIE_MAXCOMP];				// Refractive index (imaginary)
double	*mie_aext_com[MIE_MAXCOMP];				// Extinction coefficient
double	*mie_asca_com[MIE_MAXCOMP];				// Scattering coefficient
double	*mie_asym_com[MIE_MAXCOMP];				// Asymmetry parameter
double	*mie_phs1_com[MIE_MAXCOMP];				// Phase function 1
double	*mie_phs2_com[MIE_MAXCOMP];				// Phase function 2
double	*mie_phas_com[MIE_MAXCOMP];				// Phase function 3
double	mie_mixr_com[MIE_MAXCOMP];				// Number mixing ratio
double	mie_wcom_com[MIE_MAXCOMP];				// Weight
double	mie_vcom_com[MIE_MAXCOMP];				// Volume
double	mie_rmod_com[MIE_MAXCOMP];				// Mode radius in um
double	mie_lsgm_com[MIE_MAXCOMP];				// Sigma Log10(R in um)
double	*mie_xeps_com[MIE_MAXCOMP];				// Aspect ratio
double	*mie_yeps_com[MIE_MAXCOMP];				// dN(epsilon)/d(epsilon)
double	*mie_xval_com[MIE_MAXCOMP];				// R, Log10(R), or Ln(R)
double	*mie_yval_com[MIE_MAXCOMP];				// dN(R) or dV(R)
char	mie_out0[MAXLINE]	= MIE_OUT0;			// Output file 0
char	mie_out1[MAXLINE]	= MIE_OUT1;			// Output file 1
char	mie_out2[MAXLINE]	= MIE_OUT2;			// Output file 2
char	*mie_inp1[MIE_MAXCOMP];					// Input file 1
char	*mie_inp2[MIE_MAXCOMP];					// Input file 2
char	*mie_inp3[MIE_MAXCOMP];					// Input file 3
char	mie_angl_nam[MAXLINE]	= MIE_ANGL_NAM;			// Out. angle file
char	mie_wlen_nam[MAXLINE]	= MIE_WLEN_NAM;			// Out. wave  file
int	cnt_cmod		= 0;				// Calc.   mode
int	cnt_xmod		= 0;				// Mixing  mode
int	cnt_output		= 0;				// Output mode
int	cnt_long		= 0;				// Long output mode
int	cnt_db			= 0;				// Debug   mode
int	cnt_vb			= 0;				// Verbose mode
int	cnt_hp			= 0;				// Help    mode
int	cnt_optn		= 1;
char	*cnt_opts[MAXNOPT];

int Init(void);
int CommonInit(void);
void Finish(void);
void CommonFinish(void);
int MieCalc(void);
int MixComp(void);
int Printout(void);
int Normalize(int c);
int Sampling(int ninp,const double *xinp,const double *yinp,
             int nout,const double *xout,double *yout,double yuni);
int SamplingE(int ninp,const double *xinp,const double *yinp,
              int nout,const double *xout,double *yout,double yuni);
int Read1A(char *s,int size,int cx,double ux,int (*func)(double*),double *x);
int Read2A(char *s,int size,int cx,int cy,double ux,double uy,int (*func)(double*),double *x,double *y);
int Read2P(char *s,int size,int cx,int cy,double ux,double uy,int (*func)(double*),double **x,double **y);
int ReadDC(char *s,int np,int size,const int *c,const double *u,int (*func)(double*),double **d);
int ReadDA(char *s,int np,int size,const int *c,const double *u,double **d);
int ReadComp(char *s,int size,int cr,int ci,int imin,int imax,double u,
             double *wcom,double *rmod,double *lsgm,double **refr,double **refi);
int WaveSelect(double *w);
int AnglSelect(double *a);
int GetOpt(int argn,char **args);
int GetOwnOpt(int argn,char **args);
int GetCommonOpt(int argn,char **args);
int Usage(void);
int CommonUsage(char *name,int c);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(cnt_hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  if(MieCalc() < 0) exit(-1);
  if(MixComp() < 0) exit(-1);

  Printout();
  Finish();

  return 0;
}

int CommonInit(void)
{
  int i,j,n;

  if(strcmp(mie_wlen_nam,"") != 0)
  {
    if((mie_n_wlen=Read1A(mie_wlen_nam,MIE_MAXDATA,mie_wlen_num,mie_wlen_uni,WaveSelect,mie_wlen)) < 1)
    {
      return -1;
    }
  }
  if(strcmp(mie_angl_nam,"") != 0)
  {
    if((mie_n_angl=Read1A(mie_angl_nam,MIE_MAXDATA,mie_angl_num,mie_angl_uni,AnglSelect,mie_angl)) < 1)
    {
      return -1;
    }
  }
  for(i=0; i<mie_n_wlen; i++)
  {
    if(fabs(mie_wlen[i]-mie_wlen_ref) < EPSILON)
    {
      mie_iref = i;
      break;
    }
  }
  if(mie_iref < 0)
  {
    fprintf(stderr,"Failed in finding reference wavelength %13.6e\n",mie_wlen_ref);
    return -1;
  }

  if((mie_n_comp=ReadComp("stdin",MIE_MAXCOMP,mie_real_num,mie_imag_num,mie_imin,mie_imax,mie_wlen_uni,
                          mie_wcom_com,mie_rmod_com,mie_lsgm_com,mie_refr_com,mie_refi_com)) < 1)
  {
    return -1;
  }

  mie_wlen_um = (double *)malloc(mie_n_wlen*sizeof(double));
  if(mie_wlen_um == NULL)
  {
    fprintf(stderr,"Failed in allocating memory\n");
    return -1;
  }
  for(i=0; i<mie_n_wlen; i++)
  {
    mie_wlen_um[i] = mie_wlen[i]*1.0e-3;
  }
  mie_aext = (double *)malloc(mie_n_wlen*sizeof(double));
  mie_asca = (double *)malloc(mie_n_wlen*sizeof(double));
  mie_asym = (double *)malloc(mie_n_wlen*sizeof(double));
  mie_phs1 = (double *)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
  mie_phs2 = (double *)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
  mie_phas = (double *)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
  if(mie_aext==NULL || mie_asca==NULL || mie_asym==NULL ||
     mie_phs1==NULL || mie_phs2==NULL || mie_phas==NULL)
  {
    fprintf(stderr,"Failed in allocating memory\n");
    return -1;
  }
  for(n=0; n<mie_n_comp; n++)
  {
    mie_aext_com[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_asca_com[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_asym_com[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_phs1_com[n] = (double *)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
    mie_phs2_com[n] = (double *)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
    mie_phas_com[n] = (double *)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
    if(mie_aext_com[n]==NULL || mie_asca_com[n]==NULL || mie_asym_com[n]==NULL ||
       mie_phs1_com[n]==NULL || mie_phs2_com[n]==NULL || mie_phas_com[n]==NULL)
    {
      fprintf(stderr,"Failed in allocating memory\n");
      return -1;
    }
  }

  mie_angl_rad = (double *)malloc(mie_n_angl*sizeof(double));
  mie_angl_sin = (double *)malloc(mie_n_angl*sizeof(double));
  mie_angl_cos = (double *)malloc(mie_n_angl*sizeof(double));
  mie_angl_dif = (double *)malloc(mie_n_angl*sizeof(double));
  if(mie_angl_rad==NULL || mie_angl_sin==NULL || mie_angl_cos==NULL || mie_angl_dif==NULL)
  {
    fprintf(stderr,"Failed in allocating memory\n");
    return -1;
  }
  for(j=0; j<mie_n_angl; j++)
  {
    mie_angl_rad[j] = mie_angl[j]*D_TO_R;
    mie_angl_sin[j] = sin(mie_angl_rad[j]);
    mie_angl_cos[j] = cos(mie_angl_rad[j]);
  }
  for(j=1; j<mie_n_angl; j++)
  {
    mie_angl_dif[j] = PI*fabs(mie_angl_rad[j]-mie_angl_rad[j-1]);
  }

  return 0;
}

void CommonFinish(void)
{
  int n;

  for(n=0; n<mie_n_comp; n++)
  {
    free(mie_refr_com[n]);
    free(mie_refi_com[n]);
    free(mie_aext_com[n]);
    free(mie_asca_com[n]);
    free(mie_asym_com[n]);
    free(mie_phs1_com[n]);
    free(mie_phs2_com[n]);
    free(mie_phas_com[n]);
    if(mie_size_func[n] == MIE_FUNC_INPUT)
    {
      free(mie_inp1[n]);
      free(mie_inp2[n]);
    } else
    if(mie_size_func[n] == MIE_FUNC_FILE)
    {
      free(mie_xval_com[n]);
      free(mie_yval_com[n]);
    }
    if(mie_n_reps[n] > 0)
    {
      free(mie_xeps_com[n]);
      free(mie_yeps_com[n]);
    }
  }
  free(mie_wlen_um);
  free(mie_aext);
  free(mie_asca);
  free(mie_asym);
  free(mie_angl_rad);
  free(mie_angl_sin);
  free(mie_angl_cos);
  free(mie_angl_dif);
  free(mie_phs1);
  free(mie_phs2);
  free(mie_phas);

  return;
}

int Printout(void)
{
  int i,j,k,n;
  char temp1[MAXLINE];
  char temp2[MAXLINE];
  char temp3[MAXLINE];
  char *p;
  FILE *fp;
  char fnam[] = "Printout";

  if(cnt_long)
  {
    if(cnt_output > 0)
    {
      if((fp=fopen(mie_out0,"w")) == NULL)
      {
        fprintf(stderr,"%s: cannot open %s\n",fnam,mie_out0);
        return -1;
      }
      fprintf(fp,"# %2d : mie_n_comp\n",mie_n_comp);
      for(n=0; n<mie_n_comp; n++)
      {
        fprintf(fp,"%2d %22.15e %22.15e %22.15e\n",n,mie_wcom_com[n],mie_mixr_com[n],mie_vcom_com[n]);
      }
      for(n=0; n<mie_n_comp; n++)
      {
        if(mie_size_func[n] == MIE_FUNC_INPUT)
        {
          continue;
        }
        fprintf(fp,"# %2d %4d : n mie_n_wlen\n",n,mie_n_wlen);
        for(i=0; i<mie_n_wlen; i++)
        {
          fprintf(fp,"%9.4f %22.15e %22.15e\n",mie_wlen[i],mie_refr_com[n][i],mie_refi_com[n][i]);
        }
      }
      fclose(fp);
    }
    if((fp=fopen(mie_out1,"w")) == NULL)
    {
      fprintf(stderr,"%s: cannot open %s\n",fnam,mie_out1);
      return -1;
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      fprintf(fp,"%9.4f %22.15e %22.15e %22.15e %22.15e\n",
                  mie_wlen[i],mie_aext[i],mie_aext[i]/mie_aext[mie_iref],
                  mie_asca[i]/mie_aext[i],mie_asym[i]);
    }
    fclose(fp);
    if((fp=fopen(mie_out2,"w")) == NULL)
    {
      fprintf(stderr,"%s: cannot open %s\n",fnam,mie_out2);
      return -1;
    }
    for(i=0,k=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++,k++)
      {
        fprintf(fp,"%9.4f %9.4f %22.15e %22.15e %22.15e\n",
                    mie_wlen[i],mie_angl[j],mie_phs1[k],mie_phs2[k],mie_phas[k]);
      }
    }
    fclose(fp);
    if(cnt_output > 1)
    {
      strncpy(temp1,mie_out1,MAXLINE);
      if((p=strrchr(temp1,'.')) != NULL)
      {
        *p = '\0';
      }
      strncpy(temp2,mie_out2,MAXLINE);
      if((p=strrchr(temp2,'.')) != NULL)
      {
        *p = '\0';
      }
      for(n=0; n<mie_n_comp; n++)
      {
        snprintf(temp3,MAXLINE,"%s_%02d.dat",temp1,n);
        if((fp=fopen(temp3,"w")) == NULL)
        {
          fprintf(stderr,"%s: cannot open %s\n",fnam,temp3);
          return -1;
        }
        for(i=0; i<mie_n_wlen; i++)
        {
          fprintf(fp,"%9.4f %22.15e %22.15e %22.15e %22.15e\n",
                      mie_wlen[i],mie_aext_com[n][i],mie_aext_com[n][i]/mie_aext_com[n][mie_iref],
                      mie_asca_com[n][i]/mie_aext_com[n][i],mie_asym_com[n][i]);
        }
        fclose(fp);
        snprintf(temp3,MAXLINE,"%s_%02d.dat",temp2,n);
        if((fp=fopen(temp3,"w")) == NULL)
        {
          fprintf(stderr,"%s: cannot open %s\n",fnam,temp3);
          return -1;
        }
        for(i=0,k=0; i<mie_n_wlen; i++)
        {
          for(j=0; j<mie_n_angl; j++,k++)
          {
            fprintf(fp,"%9.4f %9.4f %22.15e %22.15e %22.15e\n",
                        mie_wlen[i],mie_angl[j],mie_phs1_com[n][k],
                        mie_phs2_com[n][k],mie_phas_com[n][k]);
          }
        }
        fclose(fp);
      }
    }
  }
  else
  {
    if(cnt_output > 0)
    {
      if((fp=fopen(mie_out0,"w")) == NULL)
      {
        fprintf(stderr,"%s: cannot open %s\n",fnam,mie_out0);
        return -1;
      }
      fprintf(fp,"# %2d : mie_n_comp\n",mie_n_comp);
      for(n=0; n<mie_n_comp; n++)
      {
        fprintf(fp,"%2d %13.6e %13.6e %13.6e\n",n,mie_wcom_com[n],mie_mixr_com[n],mie_vcom_com[n]);
      }
      for(n=0; n<mie_n_comp; n++)
      {
        if(mie_size_func[n] == MIE_FUNC_INPUT)
        {
          continue;
        }
        fprintf(fp,"# %2d %4d : n mie_n_wlen\n",n,mie_n_wlen);
        for(i=0; i<mie_n_wlen; i++)
        {
          fprintf(fp,"%9.4f %13.6e %13.6e\n",mie_wlen[i],mie_refr_com[n][i],mie_refi_com[n][i]);
        }
      }
      fclose(fp);
    }
    if((fp=fopen(mie_out1,"w")) == NULL)
    {
      fprintf(stderr,"%s: cannot open %s\n",fnam,mie_out1);
      return -1;
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      fprintf(fp,"%9.4f %13.6e %13.6e %13.6e %13.6e\n",
                  mie_wlen[i],mie_aext[i],mie_aext[i]/mie_aext[mie_iref],
                  mie_asca[i]/mie_aext[i],mie_asym[i]);
    }
    fclose(fp);
    if((fp=fopen(mie_out2,"w")) == NULL)
    {
      fprintf(stderr,"%s: cannot open %s\n",fnam,mie_out2);
      return -1;
    }
    for(i=0,k=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++,k++)
      {
        fprintf(fp,"%9.4f %9.4f %13.6e %13.6e %13.6e\n",
                    mie_wlen[i],mie_angl[j],mie_phs1[k],mie_phs2[k],mie_phas[k]);
      }
    }
    fclose(fp);
    if(cnt_output > 1)
    {
      strncpy(temp1,mie_out1,MAXLINE);
      if((p=strrchr(temp1,'.')) != NULL)
      {
        *p = '\0';
      }
      strncpy(temp2,mie_out2,MAXLINE);
      if((p=strrchr(temp2,'.')) != NULL)
      {
        *p = '\0';
      }
      for(n=0; n<mie_n_comp; n++)
      {
        snprintf(temp3,MAXLINE,"%s_%02d.dat",temp1,n);
        if((fp=fopen(temp3,"w")) == NULL)
        {
          fprintf(stderr,"%s: cannot open %s\n",fnam,temp3);
          return -1;
        }
        for(i=0; i<mie_n_wlen; i++)
        {
          fprintf(fp,"%9.4f %13.6e %13.6e %13.6e %13.6e\n",
                      mie_wlen[i],mie_aext_com[n][i],mie_aext_com[n][i]/mie_aext_com[n][mie_iref],
                      mie_asca_com[n][i]/mie_aext_com[n][i],mie_asym_com[n][i]);
        }
        fclose(fp);
        snprintf(temp3,MAXLINE,"%s_%02d.dat",temp2,n);
        if((fp=fopen(temp3,"w")) == NULL)
        {
          fprintf(stderr,"%s: cannot open %s\n",fnam,temp3);
          return -1;
        }
        for(i=0,k=0; i<mie_n_wlen; i++)
        {
          for(j=0; j<mie_n_angl; j++,k++)
          {
            fprintf(fp,"%9.4f %9.4f %13.6e %13.6e %13.6e\n",
                        mie_wlen[i],mie_angl[j],mie_phs1_com[n][k],
                        mie_phs2_com[n][k],mie_phas_com[n][k]);
          }
        }
        fclose(fp);
      }
    }
  }

  return 0;
}

int MixComp(void)
{
  int i,j,k,n;
  double r,ws;
  double sw,se,ss,sa;
  double s1,s2,s3;

  // calculate number mixing ratio
  if(cnt_xmod == 2)
  {
    for(r=0.0,n=0; n<mie_n_comp; n++)
    {
      mie_mixr_com[n] = mie_wcom_com[n]/mie_vcom_com[n];
      r += mie_mixr_com[n];
    }
    for(n=0; n<mie_n_comp; n++)
    {
      mie_mixr_com[n] /= r;
      if(cnt_vb > 1)
      {
        message(stderr,"Component#   : %2d\n",n);
        message(stderr,"Mixing ratio : %13.6e (%13.6e)\n",mie_mixr_com[n],mie_wcom_com[n]);
      }
    }
  } else
  if(cnt_xmod == 1)
  {
    for(r=0.0,n=0; n<mie_n_comp; n++)
    {
      mie_mixr_com[n] = mie_wcom_com[n]/mie_aext_com[n][mie_iref];
      r += mie_mixr_com[n];
    }
    for(n=0; n<mie_n_comp; n++)
    {
      mie_mixr_com[n] /= r;
      if(cnt_vb > 1)
      {
        message(stderr,"Component#   : %2d\n",n);
        message(stderr,"Mixing ratio : %13.6e (%13.6e)\n",mie_mixr_com[n],mie_wcom_com[n]);
      }
    }
  }
  else
  {
    for(n=0; n<mie_n_comp; n++)
    {
      mie_mixr_com[n] = mie_wcom_com[n];
    }
  }

  for(i=0; i<mie_n_wlen; i++)
  {
    sw = se = ss = sa = 0.0;
    if(cnt_cmod == 1) // do not calculate mie_aext,mie_asca,mie_asym
    {
      for(n=0; n<mie_n_comp; n++)
      {
        ss += mie_mixr_com[n]*mie_asca_com[n][i];
      }
    }
    else
    {
      for(n=0; n<mie_n_comp; n++)
      {
        sw += mie_mixr_com[n];
        se += mie_mixr_com[n]*mie_aext_com[n][i];
        ws  = mie_mixr_com[n]*mie_asca_com[n][i];
        ss += ws;
        sa += ws*mie_asym_com[n][i];
      }
      mie_aext[i] = se/sw;
      mie_asca[i] = ss/sw;
      mie_asym[i] = sa/ss;
    }
    for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
    {
      s1 = s2 = s3 = 0.0;
      for(n=0; n<mie_n_comp; n++)
      {
        ws  = mie_mixr_com[n]*mie_asca_com[n][i];
        s1 += ws*mie_phs1_com[n][k];
        s2 += ws*mie_phs2_com[n][k];
        s3 += ws*mie_phas_com[n][k];
      }
      mie_phs1[k] = s1;
      mie_phs2[k] = s2;
      mie_phas[k] = s3;
    }
  }

  if(Normalize(1) < 0) return -1;

  return 0;
}

int Normalize(int c)
{
  int i,j,k,n;
  double sp;
  double sinv;

  if(c == 0)
  {
    for(n=0; n<mie_n_comp; n++)
    {
      for(i=0; i<mie_n_wlen; i++)
      {
        sp = 0.0;
        for(j=1,k=mie_n_angl*i+1; j<mie_n_angl; j++,k++)
        {
          sp += (mie_phas_com[n][k-1]*mie_angl_sin[j-1]+mie_phas_com[n][k]*mie_angl_sin[j])*mie_angl_dif[j];
        }
        sinv = 1.0/sp;
        for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
        {
          mie_phs1_com[n][k] *= sinv;
          mie_phs2_com[n][k] *= sinv;
          mie_phas_com[n][k] *= sinv;
        }
      }
    }
  }
  else
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      sp = 0.0;
      for(j=1,k=mie_n_angl*i+1; j<mie_n_angl; j++,k++)
      {
        sp += (mie_phas[k-1]*mie_angl_sin[j-1]+mie_phas[k]*mie_angl_sin[j])*mie_angl_dif[j];
      }
      sinv = 1.0/sp;
      for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
      {
        mie_phs1[k] *= sinv;
        mie_phs2[k] *= sinv;
        mie_phas[k] *= sinv;
      }
    }
  }

  return 0;
}

int Sampling(int ninp,const double *xinp,const double *yinp,
             int nout,const double *xout,double *yout,double yuni)
{
  int i,j;
  int i1,i2;
  int err;
  double x1,x2;
  double xstp;
  double epsilon = 1.0e-12;
  char fnam[] = "Sampling";

  // check input
  if(ninp < 2)
  {
    fprintf(stderr,"%s: error, ninp=%d\n",fnam,ninp);
    return -1;
  }
  // check whether data interval is constant or not
  xstp = xinp[1]-xinp[0];
  err = 0;
  for(i=1; i<ninp; i++)
  {
    if(fabs(xinp[i]-(xinp[0]+xstp*i)) > epsilon)
    {
      err = 1;
      break;
    }
  }
  // interpolate
  if(err)
  {
    if(cnt_vb > 2) fprintf(stderr,"%s: nonconstant data interval.\n",fnam);
    x1 = xinp[0];
    x2 = xinp[ninp-1];
    for(i=0; i<nout; i++)
    {
      if(xout[i] < x1)
      {
        i1 = 0;
        i2 = 1;
      } else
      if(xout[i] > x2)
      {
        i2 = ninp-1;
        i1 = i2-1;
      }
      else
      {
        i1 = -1;
        i2 = -1;
        for(j=1; j<ninp; j++)
        {
          if(xout[i]>=xinp[j-1] && xout[i]<=xinp[j])
          {
            i1 = j-1;
            i2 = j;
            break;
          }
        }
        if(i1<0 || i2<0)
        {
          fprintf(stderr,"%s: error, faild in finding region for %13.6e\n",fnam,xout[i]);
          return -1;
        }
      }
      yout[i] = (yinp[i1]+(yinp[i2]-yinp[i1])*(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1]))*yuni;
    }
  }
  else
  {
    if(cnt_vb > 2) fprintf(stderr,"%s: constant data interval.\n",fnam);
    for(i=0; i<nout; i++)
    {
      i1 = (int)((xout[i]-xinp[0])/xstp+epsilon);
      if(i1 < 0) i1 = 0;
      i2 = i1+1;
      if(i2 >= ninp)
      {
        i2 = ninp-1;
        i1 = i2-1;
      }
      yout[i] = (yinp[i1]+(yinp[i2]-yinp[i1])*(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1]))*yuni;
    }
  }

  return 0;
}

int SamplingE(int ninp,const double *xinp,const double *yinp,
              int nout,const double *xout,double *yout,double yuni)
{
  int i,j;
  int i1,i2;
  int err;
  double x1,x2;
  double xstp;
  double epsilon = 1.0e-12;
  double delta = 1.0e-100;
  char fnam[] = "SamplingE";

  // check input
  if(ninp < 2)
  {
    fprintf(stderr,"%s: error, ninp=%d\n",fnam,ninp);
    return -1;
  }
  // check whether data interval is constant or not
  xstp = xinp[1]-xinp[0];
  err = 0;
  for(i=1; i<ninp; i++)
  {
    if(fabs(xinp[i]-(xinp[0]+xstp*i)) > epsilon)
    {
      err = 1;
      break;
    }
  }
  // interpolate
  if(err)
  {
    if(cnt_vb > 2) fprintf(stderr,"%s: nonconstant data interval.\n",fnam);
    x1 = xinp[0];
    x2 = xinp[ninp-1];
    for(i=0; i<nout; i++)
    {
      if(xout[i] < x1)
      {
        i1 = 0;
        i2 = 1;
      } else
      if(xout[i] > x2)
      {
        i2 = ninp-1;
        i1 = i2-1;
      }
      else
      {
        i1 = -1;
        i2 = -1;
        for(j=1; j<ninp; j++)
        {
          if(xout[i]>=xinp[j-1] && xout[i]<=xinp[j])
          {
            i1 = j-1;
            i2 = j;
            break;
          }
        }
        if(i1<0 || i2<0)
        {
          fprintf(stderr,"%s: error, faild in finding region for %13.6e\n",fnam,xout[i]);
          return -1;
        }
      }
      if(yinp[i1]<delta || yinp[i2]<delta)
      {
        yout[i] = (yinp[i1]+(yinp[i2]-yinp[i1])*(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1]))*yuni;
      }
      else
      {
        yout[i] = (yinp[i1]*pow(yinp[i2]/yinp[i1],(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1])))*yuni;
      }
    }
  }
  else
  {
    if(cnt_vb > 2) fprintf(stderr,"%s: constant data interval.\n",fnam);
    for(i=0; i<nout; i++)
    {
      i1 = (int)((xout[i]-xinp[0])/xstp+epsilon);
      if(i1 < 0) i1 = 0;
      i2 = i1+1;
      if(i2 >= ninp)
      {
        i2 = ninp-1;
        i1 = i2-1;
      }
      if(yinp[i1]<delta || yinp[i2]<delta)
      {
        yout[i] = (yinp[i1]+(yinp[i2]-yinp[i1])*(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1]))*yuni;
      }
      else
      {
        yout[i] = (yinp[i1]*pow(yinp[i2]/yinp[i1],(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1])))*yuni;
      }
    }
  }

  return 0;
}

int Read1A(char *s,int size,int cx,double ux,int (*func)(double*),double *x)
{
  int nd;
  int np = 1;
  int c[1];
  double u[1];
  double *d[1];

  // fill array
  c[0] = cx; u[0] = ux; d[0] = x;
  // read input
  if(func == NULL)
  {
    nd = ReadDA(s,size,np,c,u,d);
  }
  else
  {
    nd = ReadDC(s,size,np,c,u,func,d);
  }

  return nd;
}

int Read2A(char *s,int size,int cx,int cy,double ux,double uy,int (*func)(double*),double *x,double *y)
{
  int nd;
  int np = 2;
  int c[2];
  double u[2];
  double *d[2];

  // fill array
  c[0] = cx; u[0] = ux; d[0] = x;
  c[1] = cy; u[1] = uy; d[1] = y;
  // read input
  if(func == NULL)
  {
    nd = ReadDA(s,size,np,c,u,d);
  }
  else
  {
    nd = ReadDC(s,size,np,c,u,func,d);
  }

  return nd;
}

int Read2P(char *s,int size,int cx,int cy,double ux,double uy,int (*func)(double*),double **x,double **y)
{
  int n;
  int nd;
  int err;
  int np = 2;
  int c[2];
  double u[2];
  double *d[2];
  double **p[2];
  char fnam[] = "Read2P";

  // fill array
  c[0] = cx; u[0] = ux; d[0] = *x; p[0] = x;
  c[1] = cy; u[1] = uy; d[1] = *y; p[1] = y;
  // allocate memory
  for(n=0; n<np; n++)
  {
    if(d[n] == NULL)
    {
      if((d[n]=(double*)malloc(size*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: failed in allocating memory (%s)\n",fnam,s);
        return -1;
      }
    }
  }
  // read input
  if(func == NULL)
  {
    nd = ReadDA(s,size,np,c,u,d);
  }
  else
  {
    nd = ReadDC(s,size,np,c,u,func,d);
  }
  if(nd < 0)
  {
    for(n=0; n<np; n++)
    {
      if(*p[n] == NULL) free(d[n]);
    }
    return -1;
  }
  // resize memory
  err = 0;
  for(n=0; n<np; n++)
  {
    if(*p[n] == NULL)
    {
      if((*p[n]=(double*)realloc(d[n],nd*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: failed in allocating memory (%s)\n",fnam,s);
        free(d[n]); // this is needed only when realloc failed
        err = 1;
      }
    }
  }
  if(err)
  {
    return -1;
  }

  return nd;
}

int ReadDC(char *s,int size,int np,const int *c,const double *u,int (*func)(double*),double **d)
{
  int i,j,n;
  int nc,nd;
  int num;
  int err;
  int flag;
  double v[3];
  char line[MAXLINE];
  char str[MAXITEM][MAXLINE];
  char fnam[] = "ReadDC";
  char *p;
  FILE *fp;

  // check input
  if(np<0 || np>3)
  {
    fprintf(stderr,"%s: error, np=%d\n",fnam,np);
    return -1;
  }
  for(n=0; n<np; n++)
  {
    if(c[n]<0 || c[n]>=MAXITEM)
    {
      fprintf(stderr,"%s: invalid column number >>> c[%d]=%d (%s)\n",fnam,n,c[n],s);
      return -1;
    }
  }
  // open input
  flag = 0;
  if(strcmp(s,"") == 0)
  {
    fp = stdin;
  }
  else
  {
    if((fp=fopen(s,"r")) == NULL)
    {
      fprintf(stderr,"%s: cannot open %s\n",fnam,s);
      return -1;
    }
    flag = 1;
  }
  // read input
  num = c[0];
  for(n=1; n<np; n++)
  {
    num = MAX(num,c[n]);
  }
  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    if(line[0] == '#')
    {
      continue;
    }
    for(j=0,p=line; j<=num; j++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[j],&nc) == EOF)
      {
        if(i >= size)
        {
          err = -1;
          break;
        }
        else
        {
          fprintf(stderr,"%s: read error >>> %s (%s)\n",fnam,line,s);
          err = 1;
          break;
        }
      }
    }
    if(err) break;
    for(n=0; n<np; n++)
    {
      errno = 0;
      v[n] = strtod(str[c[n]],&p)*u[n];
      if(errno==ERANGE || *p!='\0')
      {
        if(i >= size)
        {
          err = -1;
          break;
        }
        else
        {
          fprintf(stderr,"%s: convert error >>> %s (%s)\n",fnam,line,s);
          err = 1;
          break;
        }
      }
    }
    if(err) break;
    if(func(v))
    {
      continue;
    }
    if(i >= size)
    {
      fprintf(stderr,"%s: warning, #data exceed the limit %d (%s)\n",fnam,i,s);
      break;
    }
    for(n=0; n<np; n++)
    {
      d[n][i] = v[n];
    }
    i++;
  }
  if(flag)
  {
    fclose(fp);
  }
  nd = i;
  if(nd < 1)
  {
    fprintf(stderr,"%s: error, nd=%d (%s)\n",fnam,nd,s);
    err = 1;
  }
  if(err > 0)
  {
    return -1;
  }

  return nd;
}

int ReadDA(char *s,int size,int np,const int *c,const double *u,double **d)
{
  int i,j,n;
  int nc,nd;
  int num;
  int err;
  int flag;
  double v[3];
  char line[MAXLINE];
  char str[MAXITEM][MAXLINE];
  char fnam[] = "ReadDA";
  char *p;
  FILE *fp;

  // check input
  if(np<0 || np>3)
  {
    fprintf(stderr,"%s: error, np=%d\n",fnam,np);
    return -1;
  }
  for(n=0; n<np; n++)
  {
    if(c[n]<0 || c[n]>=MAXITEM)
    {
      fprintf(stderr,"%s: invalid column number >>> c[%d]=%d (%s)\n",fnam,n,c[n],s);
      return -1;
    }
  }
  // open input
  flag = 0;
  if(strcmp(s,"") == 0)
  {
    fp = stdin;
  }
  else
  {
    if((fp=fopen(s,"r")) == NULL)
    {
      fprintf(stderr,"%s: cannot open %s\n",fnam,s);
      return -1;
    }
    flag = 1;
  }
  // read input
  num = c[0];
  for(n=1; n<np; n++)
  {
    num = MAX(num,c[n]);
  }
  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    if(line[0] == '#')
    {
      continue;
    }
    for(j=0,p=line; j<=num; j++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[j],&nc) == EOF)
      {
        if(i >= size)
        {
          err = -1;
          break;
        }
        else
        {
          fprintf(stderr,"%s: read error >>> %s (%s)\n",fnam,line,s);
          err = 1;
          break;
        }
      }
    }
    if(err) break;
    for(n=0; n<np; n++)
    {
      errno = 0;
      v[n] = strtod(str[c[n]],&p)*u[n];
      if(errno==ERANGE || *p!='\0')
      {
        if(i >= size)
        {
          err = -1;
          break;
        }
        else
        {
          fprintf(stderr,"%s: convert error >>> %s (%s)\n",fnam,line,s);
          err = 1;
          break;
        }
      }
    }
    if(err) break;
    if(i >= size)
    {
      fprintf(stderr,"%s: warning, #data exceed the limit %d (%s)\n",fnam,i,s);
      break;
    }
    for(n=0; n<np; n++)
    {
      d[n][i] = v[n];
    }
    i++;
  }
  if(flag)
  {
    fclose(fp);
  }
  nd = i;
  if(nd < 1)
  {
    fprintf(stderr,"%s: error, nd=%d (%s)\n",fnam,nd,s);
    err = 1;
  }
  if(err > 0)
  {
    return -1;
  }

  return nd;
}

int ReadComp(char *s,int size,int cr,int ci,int imin,int imax,double u,
             double *wcom,double *rmod,double *lsgm,double **refr,double **refi)
{
  int i,j,n;
  int nc,ns;
  int n_comp = 0;
  int err;
  int flag;
  int indx;
  double v;
  double *x,*y;
  char line[MAXLINE];
  char str[MAXITEM][MAXLINE];
  char fnam[] = "ReadComp";
  char *p;
  FILE *fp;

  if(strcmp(s,"")!=0 && strcmp(s,NONAME)!=0)
  {
    if(strcmp(s,"stdin") == 0)
    {
      fp = stdin;
      flag = 0;
    }
    else
    {
      if((fp=fopen(s,"r")) == NULL)
      {
        fprintf(stderr,"%s: cannot open %s\n",fnam,s);
        return -1;
      }
      flag = 1;
    }
    i = n = 0;
    err = 0;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      if(i < imin)
      {
        strncpy(line,"",MAXLINE);
        i++;
        continue;
      } else
      if(i > imax)
      {
        break;
      }
      for(ns=nc=0,p=line; ns<MAXITEM; ns++,p+=nc)
      {
        if(sscanf(p,"%s%n",str[ns],&nc) == EOF) break;
      }
      if(ns == 0)
      {
        i++;
        continue;
      } else
      if(str[0][0] == '#')
      {
        i++;
        continue;
      }
      if(strlen(str[0]) != MIE_FLAG_SIZE)
      {
        fprintf(stderr,"%s: error in flag (%s) >>> %s\n",fnam,s,line);
        err = 1;
        break;
      }
      if(n >= size)
      {
        fprintf(stderr,"%s: warning, #data exceed the limit %d (%s)\n",fnam,n,s);
        break;
      }
      switch(str[0][0]) // Function
      {
        case MIE_FUNC_INPUT:
        case MIE_FUNC_FILE:
        case MIE_FUNC_LOGNORMAL:
        case MIE_FUNC_EFF_LOGNORMAL:
        case MIE_FUNC_GAMMA:
          mie_size_func[n] = str[0][0];
          break;
        case '!':
          mie_size_func[n] = MIE_FUNC_LOGNORMAL;
          break;
        default:
          fprintf(stderr,"%s: error in function (%s) >>> %s\n",fnam,s,line);
          err = 1;
          break;
      }
      if(err) break;
      switch(str[0][1]) // X type
      {
        case MIE_XTYPE_R:
        case MIE_XTYPE_LOGR:
        case MIE_XTYPE_LNR:
          mie_size_xtype[n] = str[0][1];
          break;
        case '!':
          if(mie_size_func[n] == MIE_FUNC_LOGNORMAL)
          {
            mie_size_xtype[n] = MIE_XTYPE_LOGR;
          }
          else
          {
            mie_size_xtype[n] = MIE_XTYPE_R;
          }
          break;
        default:
          fprintf(stderr,"%s: error in X type (%s) >>> %s\n",fnam,s,line);
          err = 1;
          break;
      }
      if(err) break;
      switch(str[0][2]) // Y type
      {
        case MIE_YTYPE_NUMBER:
        case MIE_YTYPE_VOLUME:
          mie_size_ytype[n] = str[0][2];
          break;
        case '!':
          mie_size_ytype[n] = MIE_YTYPE_NUMBER;
          break;
        default:
          fprintf(stderr,"%s: error in Y type (%s) >>> %s\n",fnam,s,line);
          err = 1;
          break;
      }
      if(err) break;
      switch(str[0][3]) // Shape
      {
        case MIE_SHAPE_SPHERE:
        case MIE_SHAPE_SPHEROID:
          mie_size_shape[n] = str[0][3];
          break;
        case '!':
          mie_size_shape[n] = MIE_SHAPE_SPHERE;
          break;
        default:
          fprintf(stderr,"%s: error in shape (%s) >>> %s\n",fnam,s,line);
          err = 1;
          break;
      }
      if(err) break;
      switch(mie_size_func[n])
      {
        case MIE_FUNC_INPUT:
          if(cnt_xmod == 2)
          {
            if(ns < 5)
            {
              fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
              fprintf(stderr,"%s: cnt_xmod=%d, ns=%d\n",fnam,cnt_xmod,ns);
              err = 1;
              break;
            }
            if(str[4][0] == '@')
            {
              v = strtod(&str[4][1],&p);
              if(errno==ERANGE || *p!='\0')
              {
                fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
                err = 1;
                break;
              }
              if(v <= 0.0)
              {
                fprintf(stderr,"%s: error in volume (%s) >>> %s\n",fnam,s,line);
                err = 1;
                break;
              }
              mie_vcom_com[n] = v;
            }
            else
            {
              if(ns > 5)
              {
                indx = strtol(str[5],&p,10);
                if(errno==ERANGE || *p!='\0')
                {
                  fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
                  err = 1;
                  break;
                }
              }
              else
              {
                indx = n;
              }
              x = NULL;
              y = NULL;
              if((j=Read2P(str[4],indx+1,0,3,1.0,1.0,NULL,&x,&y)) < 0)
              {
                err = 1;
                break;
              }
              do
              {
                if(j != indx+1)
                {
                  fprintf(stderr,"%s: error, j=%d, indx+1=%d\n",fnam,j,indx+1);
                  err = 1;
                  break;
                }
                if(fabs(x[indx]-(double)indx) > EPSILON)
                {
                  fprintf(stderr,"%s: error, x[%d]=%13.6e, indx=%d\n",fnam,indx,x[indx],indx);
                  err = 1;
                  break;
                }
                mie_vcom_com[n] = y[indx];
              }
              while(0);
              free(x);
              free(y);
              if(err) break;
            }
          }
          else
          {
            if(ns < 4)
            {
              fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
              fprintf(stderr,"%s: mie_size_func[%d]=%c, ns=%d\n",fnam,n,mie_size_func[n],ns);
              err = 1;
              break;
            }
          }
          // Mixing ratio
          wcom[n] = strtod(str[1],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(wcom[n] < 0.0)
          {
            fprintf(stderr,"%s: error in mixing ratio (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          // Optical parameters
          mie_inp1[n] = (char *)malloc((strlen(str[2])+1)*sizeof(char));
          mie_inp2[n] = (char *)malloc((strlen(str[3])+1)*sizeof(char));
          if(mie_inp1[n]==NULL || mie_inp2[n]==NULL)
          {
            fprintf(stderr,"%s: error in allocating memory.\n",fnam);
            err = 1;
            break;
          }
          strncpy(mie_inp1[n],str[2],strlen(str[2])+1);
          strncpy(mie_inp2[n],str[3],strlen(str[3])+1);
          break;
        case MIE_FUNC_FILE: // flag wcom fnam refr refi (reps) comment
          // Aspect ratio
          switch(mie_size_shape[n])
          {
            case MIE_SHAPE_SPHERE:
              if(ns < 5)
              {
                fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
                fprintf(stderr,"%s: mie_size_func[%d]=%c, ns=%d\n",fnam,n,mie_size_func[n],ns);
                err = 1;
                break;
              }
              break;
            case MIE_SHAPE_SPHEROID:
              if(ns < 6)
              {
                fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
                fprintf(stderr,"%s: mie_size_func[%d]=%c, ns=%d\n",fnam,n,mie_size_func[n],ns);
                err = 1;
                break;
              }
              if(str[5][0] == '@')
              {
                v = strtod(&str[5][1],&p);
                if(errno==ERANGE || *p!='\0')
                {
                  fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
                  err = 1;
                  break;
                }
                if(v <= 0.0)
                {
                  fprintf(stderr,"%s: error in aspect ratio (%s) >>> %s\n",fnam,s,line);
                  err = 1;
                  break;
                }
                mie_xeps_com[n] = (double *)malloc(sizeof(double));
                mie_yeps_com[n] = (double *)malloc(sizeof(double));
                if(mie_xeps_com[n]==NULL || mie_yeps_com[n]==NULL)
                {
                  fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                  err = 1;
                  break;
                }
                mie_xeps_com[n][0] = v;
                mie_yeps_com[n][0] = 1.0;
                mie_n_reps[n] = 1;
              }
              else
              {
                if((mie_n_reps[n]=Read2P(str[5],MIE_DMAX,0,1,1.0,1.0,NULL,&mie_xeps_com[n],&mie_yeps_com[n])) < 0)
                {
                  err = 1;
                  break;
                }
                if(mie_n_reps[n] < 1)
                {
                  fprintf(stderr,"%s: error, mie_n_reps[%d]=%d\n",fnam,n,mie_n_reps[n]);
                  free(mie_xeps_com[n]);
                  free(mie_yeps_com[n]);
                  err = 1;
                  break;
                }
              }
              if(cnt_vb > 2)
              {
                fprintf(stderr,"Axis ratio distribution of component %d:\n",n);
                for(j=0; j<mie_n_reps[n]; j++)
                {
                  fprintf(stderr,"%13.6e %13.6e\n",mie_xeps_com[n][j],mie_yeps_com[n][j]);
                }
              }
              break;
            default:
              fprintf(stderr,"%s: error in shape (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
          }
          if(err) break;
          // Mixing ratio
          wcom[n] = strtod(str[1],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(wcom[n] < 0.0)
          {
            fprintf(stderr,"%s: error in mixing ratio (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          // Size distribution
          if((mie_n_size[n]=Read2P(str[2],MIE_DMAX,0,1,1.0,1.0,NULL,&mie_xval_com[n],&mie_yval_com[n])) < 0)
          {
            err = 1;
            break;
          }
          if(mie_n_size[n] < 2)
          {
            fprintf(stderr,"%s: error, mie_n_size[%d]=%d\n",fnam,n,mie_n_size[n]);
            free(mie_xval_com[n]);
            free(mie_yval_com[n]);
            err = 1;
            break;
          }
          if(cnt_vb > 2)
          {
            fprintf(stderr,"Size distribution of component %d:\n",n);
            for(j=0; j<mie_n_size[n]; j++)
            {
              fprintf(stderr,"%13.6e %13.6e\n",mie_xval_com[n][j],mie_yval_com[n][j]);
            }
          }
          // Refractive index (real)
          if(str[3][0] == '@')
          {
            v = strtod(&str[3][1],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
            }
            refr[n] = NULL;
            if((refr[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
            {
              fprintf(stderr,"%s: error in allocating memory.\n",fnam);
              err = 1;
              break;
            }
            for(j=0; j<mie_n_wlen; j++)
            {
              refr[n][j] = v;
            }
          }
          else
          {
            x = NULL;
            y = NULL;
            if((j=Read2P(str[3],MIE_MAXDATA,0,cr,u,1.0,NULL,&x,&y)) < 0)
            {
              err = 1;
              break;
            }
            do
            {
              refr[n] = NULL;
              if((refr[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
              {
                fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                err = 1;
                break;
              }
              if(Sampling(j,x,y,mie_n_wlen,mie_wlen,refr[n],1.0) < 0)
              {
                err = 1;
                break;
              }
            }
            while(0);
            free(x);
            free(y);
            if(err) break;
          }
          // Refractive index (imag)
          if(str[4][0] == '@')
          {
            v = strtod(&str[4][1],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
            }
            refi[n] = NULL;
            if((refi[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
            {
              fprintf(stderr,"%s: error in allocating memory.\n",fnam);
              err = 1;
              break;
            }
            for(j=0; j<mie_n_wlen; j++)
            {
              refi[n][j] = v;
            }
          }
          else
          {
            x = NULL;
            y = NULL;
            if((j=Read2P(str[4],MIE_MAXDATA,0,ci,u,1.0,NULL,&x,&y)) < 0)
            {
              err = 1;
              break;
            }
            do
            {
              refi[n] = NULL;
              if((refi[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
              {
                fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                err = 1;
                break;
              }
              if(SamplingE(j,x,y,mie_n_wlen,mie_wlen,refi[n],1.0) < 0)
              {
                err = 1;
                break;
              }
            }
            while(0);
            free(x);
            free(y);
            if(err) break;
          }
          break;
        case MIE_FUNC_LOGNORMAL: // flag wcom rmod lsgm refr refi (reps) comment
          // Aspect ratio
          switch(mie_size_shape[n])
          {
            case MIE_SHAPE_SPHERE:
              if(ns < 6)
              {
                fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
                fprintf(stderr,"%s: mie_size_func[%d]=%c, ns=%d\n",fnam,n,mie_size_func[n],ns);
                err = 1;
                break;
              }
              break;
            case MIE_SHAPE_SPHEROID:
              if(ns < 7)
              {
                fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
                fprintf(stderr,"%s: mie_size_func[%d]=%c, ns=%d\n",fnam,n,mie_size_func[n],ns);
                err = 1;
                break;
              }
              if(str[6][0] == '@')
              {
                v = strtod(&str[6][1],&p);
                if(errno==ERANGE || *p!='\0')
                {
                  fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
                  err = 1;
                  break;
                }
                if(v <= 0.0)
                {
                  fprintf(stderr,"%s: error in aspect ratio (%s) >>> %s\n",fnam,s,line);
                  err = 1;
                  break;
                }
                mie_xeps_com[n] = (double *)malloc(sizeof(double));
                mie_yeps_com[n] = (double *)malloc(sizeof(double));
                if(mie_xeps_com[n]==NULL || mie_yeps_com[n]==NULL)
                {
                  fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                  err = 1;
                  break;
                }
                mie_xeps_com[n][0] = v;
                mie_yeps_com[n][0] = 1.0;
                mie_n_reps[n] = 1;
              }
              else
              {
                if((mie_n_reps[n]=Read2P(str[6],MIE_DMAX,0,1,1.0,1.0,NULL,&mie_xeps_com[n],&mie_yeps_com[n])) < 0)
                {
                  err = 1;
                  break;
                }
                if(mie_n_reps[n] < 1)
                {
                  fprintf(stderr,"%s: error, mie_n_reps[%d]=%d\n",fnam,n,mie_n_reps[n]);
                  free(mie_xeps_com[n]);
                  free(mie_yeps_com[n]);
                  err = 1;
                  break;
                }
              }
              if(cnt_vb > 2)
              {
                fprintf(stderr,"Axis ratio distribution of component %d:\n",n);
                for(j=0; j<mie_n_reps[n]; j++)
                {
                  fprintf(stderr,"%13.6e %13.6e\n",mie_xeps_com[n][j],mie_yeps_com[n][j]);
                }
              }
              break;
            default:
              fprintf(stderr,"%s: error in shape (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
          }
          if(err) break;
          // Mixing ratio
          wcom[n] = strtod(str[1],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(wcom[n] < 0.0)
          {
            fprintf(stderr,"%s: error in mixing ratio (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          // Mode radius
          rmod[n] = strtod(str[2],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(rmod[n] <= 0.0)
          {
            fprintf(stderr,"%s: error in mode radius (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          // Sigma
          lsgm[n] = strtod(str[3],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(lsgm[n] <= 0.0)
          {
            fprintf(stderr,"%s: error in sigma (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(mie_size_xtype[n] == MIE_XTYPE_LNR)
          {
            lsgm[n] *= M_1_LNA; // dY/dLnR -> dY/dLogR
          }
          // Refractive index (real)
          if(str[4][0] == '@')
          {
            v = strtod(&str[4][1],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
            }
            refr[n] = NULL;
            if((refr[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
            {
              fprintf(stderr,"%s: error in allocating memory.\n",fnam);
              err = 1;
              break;
            }
            for(j=0; j<mie_n_wlen; j++)
            {
              refr[n][j] = v;
            }
          }
          else
          {
            x = NULL;
            y = NULL;
            if((j=Read2P(str[4],MIE_MAXDATA,0,cr,u,1.0,NULL,&x,&y)) < 0)
            {
              err = 1;
              break;
            }
            do
            {
              refr[n] = NULL;
              if((refr[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
              {
                fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                err = 1;
                break;
              }
              if(Sampling(j,x,y,mie_n_wlen,mie_wlen,refr[n],1.0) < 0)
              {
                err = 1;
                break;
              }
            }
            while(0);
            free(x);
            free(y);
            if(err) break;
          }
          // Refractive index (imag)
          if(str[5][0] == '@')
          {
            v = strtod(&str[5][1],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
            }
            refi[n] = NULL;
            if((refi[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
            {
              fprintf(stderr,"%s: error in allocating memory.\n",fnam);
              err = 1;
              break;
            }
            for(j=0; j<mie_n_wlen; j++)
            {
              refi[n][j] = v;
            }
          }
          else
          {
            x = NULL;
            y = NULL;
            if((j=Read2P(str[5],MIE_MAXDATA,0,ci,u,1.0,NULL,&x,&y)) < 0)
            {
              err = 1;
              break;
            }
            do
            {
              refi[n] = NULL;
              if((refi[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
              {
                fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                err = 1;
                break;
              }
              if(SamplingE(j,x,y,mie_n_wlen,mie_wlen,refi[n],1.0) < 0)
              {
                err = 1;
                break;
              }
            }
            while(0);
            free(x);
            free(y);
            if(err) break;
          }
          break;
        case MIE_FUNC_EFF_LOGNORMAL: // flag wcom reff veff refr refi (reps) comment
          // Aspect ratio
          switch(mie_size_shape[n])
          {
            case MIE_SHAPE_SPHERE:
              if(ns < 6)
              {
                fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
                fprintf(stderr,"%s: mie_size_func[%d]=%c, ns=%d\n",fnam,n,mie_size_func[n],ns);
                err = 1;
                break;
              }
              break;
            case MIE_SHAPE_SPHEROID:
              if(ns < 7)
              {
                fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
                fprintf(stderr,"%s: mie_size_func[%d]=%c, ns=%d\n",fnam,n,mie_size_func[n],ns);
                err = 1;
                break;
              }
              if(str[6][0] == '@')
              {
                v = strtod(&str[6][1],&p);
                if(errno==ERANGE || *p!='\0')
                {
                  fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
                  err = 1;
                  break;
                }
                if(v <= 0.0)
                {
                  fprintf(stderr,"%s: error in aspect ratio (%s) >>> %s\n",fnam,s,line);
                  err = 1;
                  break;
                }
                mie_xeps_com[n] = (double *)malloc(sizeof(double));
                mie_yeps_com[n] = (double *)malloc(sizeof(double));
                if(mie_xeps_com[n]==NULL || mie_yeps_com[n]==NULL)
                {
                  fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                  err = 1;
                  break;
                }
                mie_xeps_com[n][0] = v;
                mie_yeps_com[n][0] = 1.0;
                mie_n_reps[n] = 1;
              }
              else
              {
                if((mie_n_reps[n]=Read2P(str[6],MIE_DMAX,0,1,1.0,1.0,NULL,&mie_xeps_com[n],&mie_yeps_com[n])) < 0)
                {
                  err = 1;
                  break;
                }
                if(mie_n_reps[n] < 1)
                {
                  fprintf(stderr,"%s: error, mie_n_reps[%d]=%d\n",fnam,n,mie_n_reps[n]);
                  free(mie_xeps_com[n]);
                  free(mie_yeps_com[n]);
                  err = 1;
                  break;
                }
              }
              if(cnt_vb > 2)
              {
                fprintf(stderr,"Axis ratio distribution of component %d:\n",n);
                for(j=0; j<mie_n_reps[n]; j++)
                {
                  fprintf(stderr,"%13.6e %13.6e\n",mie_xeps_com[n][j],mie_yeps_com[n][j]);
                }
              }
              break;
            default:
              fprintf(stderr,"%s: error in shape (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
          }
          if(err) break;
          // Mixing ratio
          wcom[n] = strtod(str[1],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(wcom[n] < 0.0)
          {
            fprintf(stderr,"%s: error in mixing ratio (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          // Mode radius
          rmod[n] = strtod(str[2],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(rmod[n] <= 0.0)
          {
            fprintf(stderr,"%s: error in effective radius (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          // Sigma
          lsgm[n] = strtod(str[3],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(lsgm[n] <= 0.0)
          {
            fprintf(stderr,"%s: error in effective variance (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          v = log(lsgm[n]+1.0);
          rmod[n] = rmod[n]*exp(-2.5*v);
          lsgm[n] = sqrt(v)*M_1_LNA; // dY/dLogR
          // Refractive index (real)
          if(str[4][0] == '@')
          {
            v = strtod(&str[4][1],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
            }
            refr[n] = NULL;
            if((refr[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
            {
              fprintf(stderr,"%s: error in allocating memory.\n",fnam);
              err = 1;
              break;
            }
            for(j=0; j<mie_n_wlen; j++)
            {
              refr[n][j] = v;
            }
          }
          else
          {
            x = NULL;
            y = NULL;
            if((j=Read2P(str[4],MIE_MAXDATA,0,cr,u,1.0,NULL,&x,&y)) < 0)
            {
              err = 1;
              break;
            }
            do
            {
              refr[n] = NULL;
              if((refr[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
              {
                fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                err = 1;
                break;
              }
              if(Sampling(j,x,y,mie_n_wlen,mie_wlen,refr[n],1.0) < 0)
              {
                err = 1;
                break;
              }
            }
            while(0);
            free(x);
            free(y);
            if(err) break;
          }
          // Refractive index (imag)
          if(str[5][0] == '@')
          {
            v = strtod(&str[5][1],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
            }
            refi[n] = NULL;
            if((refi[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
            {
              fprintf(stderr,"%s: error in allocating memory.\n",fnam);
              err = 1;
              break;
            }
            for(j=0; j<mie_n_wlen; j++)
            {
              refi[n][j] = v;
            }
          }
          else
          {
            x = NULL;
            y = NULL;
            if((j=Read2P(str[5],MIE_MAXDATA,0,ci,u,1.0,NULL,&x,&y)) < 0)
            {
              err = 1;
              break;
            }
            do
            {
              refi[n] = NULL;
              if((refi[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
              {
                fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                err = 1;
                break;
              }
              if(SamplingE(j,x,y,mie_n_wlen,mie_wlen,refi[n],1.0) < 0)
              {
                err = 1;
                break;
              }
            }
            while(0);
            free(x);
            free(y);
            if(err) break;
          }
          break;
        case MIE_FUNC_GAMMA: // flag wcom reff veff refr refi (reps) comment
          // Aspect ratio
          switch(mie_size_shape[n])
          {
            case MIE_SHAPE_SPHERE:
              if(ns < 6)
              {
                fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
                fprintf(stderr,"%s: mie_size_func[%d]=%c, ns=%d\n",fnam,n,mie_size_func[n],ns);
                err = 1;
                break;
              }
              break;
            case MIE_SHAPE_SPHEROID:
              if(ns < 7)
              {
                fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
                fprintf(stderr,"%s: mie_size_func[%d]=%c, ns=%d\n",fnam,n,mie_size_func[n],ns);
                err = 1;
                break;
              }
              if(str[6][0] == '@')
              {
                v = strtod(&str[6][1],&p);
                if(errno==ERANGE || *p!='\0')
                {
                  fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
                  err = 1;
                  break;
                }
                if(v <= 0.0)
                {
                  fprintf(stderr,"%s: error in aspect ratio (%s) >>> %s\n",fnam,s,line);
                  err = 1;
                  break;
                }
                mie_xeps_com[n] = (double *)malloc(sizeof(double));
                mie_yeps_com[n] = (double *)malloc(sizeof(double));
                if(mie_xeps_com[n]==NULL || mie_yeps_com[n]==NULL)
                {
                  fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                  err = 1;
                  break;
                }
                mie_xeps_com[n][0] = v;
                mie_yeps_com[n][0] = 1.0;
                mie_n_reps[n] = 1;
              }
              else
              {
                if((mie_n_reps[n]=Read2P(str[6],MIE_DMAX,0,1,1.0,1.0,NULL,&mie_xeps_com[n],&mie_yeps_com[n])) < 0)
                {
                  err = 1;
                  break;
                }
                if(mie_n_reps[n] < 1)
                {
                  fprintf(stderr,"%s: error, mie_n_reps[%d]=%d\n",fnam,n,mie_n_reps[n]);
                  free(mie_xeps_com[n]);
                  free(mie_yeps_com[n]);
                  err = 1;
                  break;
                }
              }
              if(cnt_vb > 2)
              {
                fprintf(stderr,"Axis ratio distribution of component %d:\n",n);
                for(j=0; j<mie_n_reps[n]; j++)
                {
                  fprintf(stderr,"%13.6e %13.6e\n",mie_xeps_com[n][j],mie_yeps_com[n][j]);
                }
              }
              break;
            default:
              fprintf(stderr,"%s: error in shape (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
          }
          if(err) break;
          // Mixing ratio
          wcom[n] = strtod(str[1],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(wcom[n] < 0.0)
          {
            fprintf(stderr,"%s: error in mixing ratio (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          // Effective radius
          rmod[n] = strtod(str[2],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(rmod[n] <= 0.0)
          {
            fprintf(stderr,"%s: error in effective radius (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          // Effective variance
          lsgm[n] = strtod(str[3],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          if(lsgm[n] <= 0.0)
          {
            fprintf(stderr,"%s: error in effective variance (%s) >>> %s\n",fnam,s,line);
            err = 1;
            break;
          }
          // Refractive index (real)
          if(str[4][0] == '@')
          {
            v = strtod(&str[4][1],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
            }
            refr[n] = NULL;
            if((refr[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
            {
              fprintf(stderr,"%s: error in allocating memory.\n",fnam);
              err = 1;
              break;
            }
            for(j=0; j<mie_n_wlen; j++)
            {
              refr[n][j] = v;
            }
          }
          else
          {
            x = NULL;
            y = NULL;
            if((j=Read2P(str[4],MIE_MAXDATA,0,cr,u,1.0,NULL,&x,&y)) < 0)
            {
              err = 1;
              break;
            }
            do
            {
              refr[n] = NULL;
              if((refr[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
              {
                fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                err = 1;
                break;
              }
              if(Sampling(j,x,y,mie_n_wlen,mie_wlen,refr[n],1.0) < 0)
              {
                err = 1;
                break;
              }
            }
            while(0);
            free(x);
            free(y);
            if(err) break;
          }
          // Refractive index (imag)
          if(str[5][0] == '@')
          {
            v = strtod(&str[5][1],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
              err = 1;
              break;
            }
            refi[n] = NULL;
            if((refi[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
            {
              fprintf(stderr,"%s: error in allocating memory.\n",fnam);
              err = 1;
              break;
            }
            for(j=0; j<mie_n_wlen; j++)
            {
              refi[n][j] = v;
            }
          }
          else
          {
            x = NULL;
            y = NULL;
            if((j=Read2P(str[5],MIE_MAXDATA,0,ci,u,1.0,NULL,&x,&y)) < 0)
            {
              err = 1;
              break;
            }
            do
            {
              refi[n] = NULL;
              if((refi[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
              {
                fprintf(stderr,"%s: error in allocating memory.\n",fnam);
                err = 1;
                break;
              }
              if(SamplingE(j,x,y,mie_n_wlen,mie_wlen,refi[n],1.0) < 0)
              {
                err = 1;
                break;
              }
            }
            while(0);
            free(x);
            free(y);
            if(err) break;
          }
          break;
        default:
          fprintf(stderr,"%s: error in function (%s) >>> %s\n",fnam,s,line);
          err = 1;
          break;
      }
      if(err) break;
      strncpy(line,"",MAXLINE);
      n++;
      i++;
    }
    if(flag)
    {
      fclose(fp);
    }
    if(err)
    {
      return -1;
    }
    n_comp = n;
    fprintf(stderr,"%s: %d data have been read (%s)\n",fnam,n_comp,s);
    if(n_comp < 1)
    {
      fprintf(stderr,"%s: error, n_comp=%d (%s)\n",fnam,n_comp,s);
      return -1;
    }
  }

  return n_comp;
}

int WaveSelect(double *w)
{
  if(w[0]<mie_wlen_min || w[0]>mie_wlen_max)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int AnglSelect(double *a)
{
  if(a[0]<mie_angl_min || a[0]>mie_angl_max)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int GetOpt(int argn,char **args)
{
  int rt;

  rt = GetOwnOpt(argn,args);
  if(cnt_optn > 1)
  {
    if(GetCommonOpt(cnt_optn,cnt_opts) < 0) rt = -1;
  }
  if(cnt_hp) rt = 1;

  return rt;
}

int GetCommonOpt(int argn,char **args)
{
  int c,rt;
  int option_index = 0;
  int this_option_optind;
  int ntmp;
  char *endp;
  double xtmp;
  struct option long_options[] =
  {
    {"mie_out1",1,0,'f'},
    {"mie_out2",1,0,'g'},
    {"mie_wlen_nam",1,0,'w'},
    {"mie_wlen_min",1,0,'l'},
    {"mie_wlen_max",1,0,'L'},
    {"mie_angl_nam",1,0,'a'},
    {"mie_wlen_num",1,0,'W'},
    {"mie_angl_num",1,0,'A'},
    {"mie_real_num",1,0,'P'},
    {"mie_imag_num",1,0,'Q'},
    {"mie_imin",1,0,'i'},
    {"mie_imax",1,0,'I'},
    {"mie_wlen_ref",1,0,'r'},
    {"mie_wlen_uni",1,0,'u'},
    {"cnt_xmod",1,0,'m'},
    {"output",0,0,'o'},
    {"long",0,0,'O'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":f:g:w:l:L:a:W:A:P:Q:i:I:r:u:m:oOdvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'f':
        strncpy(mie_out1,optarg,MAXLINE);
        break;
      case 'g':
        strncpy(mie_out2,optarg,MAXLINE);
        break;
      case 'w':
        strncpy(mie_wlen_nam,optarg,MAXLINE);
        break;
      case 'l':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_wlen_min = xtmp;
        else
        {
          fprintf(stderr,"Min wavelength -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'L':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_wlen_max = xtmp;
        else
        {
          fprintf(stderr,"Max wavelength -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'a':
        strncpy(mie_angl_nam,optarg,MAXLINE);
        break;
      case 'W':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<MAXITEM) mie_wlen_num = ntmp;
        else
        {
          fprintf(stderr,"Out. wave column# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'A':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<MAXITEM) mie_angl_num = ntmp;
        else
        {
          fprintf(stderr,"Out. angle column# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'P':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<MAXITEM) mie_real_num = ntmp;
        else
        {
          fprintf(stderr,"Refractive index (real) column# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'Q':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<MAXITEM) mie_imag_num = ntmp;
        else
        {
          fprintf(stderr,"Refractive index (imag) column# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'i':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0) mie_imin = ntmp;
        else
        {
          fprintf(stderr,"Minimum line# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'I':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0) mie_imax = ntmp;
        else
        {
          fprintf(stderr,"Maximum line# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'r':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_wlen_ref = xtmp;
        else
        {
          fprintf(stderr,"Reference wavelength -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'u':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_wlen_uni = xtmp;
        else
        {
          fprintf(stderr,"Wavelength unit -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'm':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<=2) cnt_xmod = ntmp;
        else
        {
          fprintf(stderr,"Mixing mode -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'o':
        cnt_output++;
        break;
      case 'O':
        cnt_long = 1;
        break;
      case 'd':
        cnt_db++;
        break;
      case 'v':
        cnt_vb++;
        break;
      case 'h':
        cnt_hp = 1;
        break;
      case '?':
        if(optopt == '\0')
        {
          fprintf(stderr,"Invalid option : %s\n",args[this_option_optind]);
        }
        else
        {
          fprintf(stderr,"Invalid option : %c\n",optopt);
        }
        rt = -1;
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

  if(optind < argn)
  {
    fprintf(stderr,"non-option ARGV-elements:\n");
    while(optind < argn)
    {
      fprintf(stderr,"%s\n",args[optind++]);
    }
    rt = -1;
  }

  if(cnt_hp) rt = 1;
  return rt;
}

int CommonUsage(char *name,int c)
{
  int  n = 16;
  char e[MAXLINE];
  char a[MAXLINE];
  char d[MAXLINE];

  switch(c)
  {
    case 0:
      fprintf(stderr,"Usage:\n");
      fprintf(stderr,"%s -[option] (argument) -[option] (argument) ...\n",name);
      fprintf(stderr,"------------------------------------------------------------------------------------\n");
      fprintf(stderr,"   option         |%s|%s|%s| current\n",As(e,"",n),          As(a,"argument",n),   As(d,"default",n));
      break;
    case 'f':
      fprintf(stderr," f -mie_out1      |%s|%s|%s| %s\n",As(e,"Out. file 1",n),    As(a,"name",n),       As(d,MIE_OUT1,n),mie_out1);
      break;
    case 'g':
      fprintf(stderr," g -mie_out2      |%s|%s|%s| %s\n",As(e,"Out. file 2",n),    As(a,"name",n),       As(d,MIE_OUT2,n),mie_out2);
      break;
    case 'w':
      fprintf(stderr," w -mie_wlen_nam  |%s|%s|%s| %s\n",As(e,"Out. wave file",n), As(a,"name",n),       As(d,MIE_WLEN_NAM,n),mie_wlen_nam);
      break;
    case 'l':
      fprintf(stderr," l -mie_wlen_min  |%s|%s|%s| %f\n",As(e,"Min wavelength",n), As(a,"nm",n),         Af(d,MIE_WLEN_MIN,n),mie_wlen_min);
      break;
    case 'L':
      fprintf(stderr," L -mie_wlen_max  |%s|%s|%s| %f\n",As(e,"Max wavelength",n), As(a,"nm",n),         Af(d,MIE_WLEN_MAX,n),mie_wlen_max);
      break;
    case 'a':
      fprintf(stderr," a -mie_angl_nam  |%s|%s|%s| %s\n",As(e,"Out. ang. file",n), As(a,"name",n),       As(d,MIE_ANGL_NAM,n),mie_angl_nam);
      break;
    case 'W':
      fprintf(stderr," W -mie_wlen_num  |%s|%s|%s| %d\n",As(e,"Out. wave col#",n), As(a,"#",n),          Ad(d,MIE_WLEN_NUM,n),mie_wlen_num);
      break;
    case 'A':
      fprintf(stderr," A -mie_angl_num  |%s|%s|%s| %d\n",As(e,"Out. ang. col#",n), As(a,"#",n),          Ad(d,MIE_ANGL_NUM,n),mie_angl_num);
      break;
    case 'P':
      fprintf(stderr," P -mie_real_num  |%s|%s|%s| %d\n",As(e,"Rindx (r) col#",n), As(a,"#",n),          Ad(d,MIE_REAL_NUM,n),mie_real_num);
      break;
    case 'Q':
      fprintf(stderr," Q -mie_imag_num  |%s|%s|%s| %d\n",As(e,"Rindx (i) col#",n), As(a,"#",n),          Ad(d,MIE_IMAG_NUM,n),mie_imag_num);
      break;
    case 'i':
      fprintf(stderr," i -mie_imin      |%s|%s|%s| %d\n",As(e,"Minimum line#",n),  As(a,"#",n),          Ad(d,MIE_IMIN,n),mie_imin);
      break;
    case 'I':
      fprintf(stderr," I -mie_imax      |%s|%s|%s| %d\n",As(e,"Maximum line#",n),  As(a,"#",n),          Ad(d,MIE_IMAX,n),mie_imax);
      break;
    case 'r':
      fprintf(stderr," r -mie_wlen_ref  |%s|%s|%s| %f\n",As(e,"Reference wave",n), As(a,"nm",n),         Af(d,MIE_WLEN_REF,n),mie_wlen_ref);
      break;
    case 'u':
      fprintf(stderr," u -mie_wlen_uni  |%s|%s|%s| %e\n",As(e,"Wave unit",n),      As(a,"nm",n),         Af(d,MIE_WLEN_UNI,n),mie_wlen_uni);
      break;
    case 'm':
      fprintf(stderr," m -cnt_xmod      |%s|%s|%s| %d\n",As(e,"Mixing mode",n),    As(a,"0|1|2",n),      Ad(d,0,n),cnt_xmod);
      break;
    case 'o':
      fprintf(stderr," o -output        |%s|%s|%s| %d\n",As(e,"Output  mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_output);
      break;
    case 'O':
      fprintf(stderr," O -long          |%s|%s|%s| %d\n",As(e,"Long    mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_long);
      break;
    case 'd':
      fprintf(stderr," d -debug         |%s|%s|%s| %d\n",As(e,"Debug   mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_db);
      break;
    case 'v':
      fprintf(stderr," v -verbose       |%s|%s|%s| %d\n",As(e,"Verbose mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_vb);
      break;
    case 'h':
      fprintf(stderr," h -help          |%s|%s|%s| %d\n",As(e,"Help    mode",n),   As(a,"nothing",n),    Ad(d,0,n),1);
      break;
    case 1:
      fprintf(stderr,"------------------------------------------------------------------------------------\n");
      break;
    case 2:
      fprintf(stderr,"Mixing mode:\n");
      fprintf(stderr,"  0 ... Number ratio\n");
      fprintf(stderr,"  1 ... Extinction cross-section ratio\n");
      fprintf(stderr,"  2 ... Volume ratio\n");
      break;
    case 3:
      fprintf(stderr,"Input data format (stdin): flag param1 param2 ...\n");
      fprintf(stderr,"A flag is composed of 4 characters.\n");
      fprintf(stderr,"  1st character ... Function (I=Input,F=File,L=Lognormal,E=Effective Lognormal,G=Gamma,!=Lognormal)\n");
      fprintf(stderr,"  2nd character ... X type (R=Radius,L=Log10(R),N=Ln(R),!=Log10(R) if Function==L else Radius)\n");
      fprintf(stderr,"  3rd character ... Y type (N=Number distribution,V=Volume distribution,!=Number)\n");
      fprintf(stderr,"  4th character ... Shape (S=Sphere,T=Spheroid,!=Sphere)\n");
      fprintf(stderr,"Input data format (stdin) for Input:               flag wcom inp1 inp2 (inp3) (indx) comment\n");
      fprintf(stderr,"Input data format (stdin) for File:                flag wcom fnam refr refi (reps) comment\n");
      fprintf(stderr,"Input data format (stdin) for Lognormal:           flag wcom rmod lsgm refr refi (reps) comment\n");
      fprintf(stderr,"Input data format (stdin) for Effective Lognormal: flag wcom reff veff refr refi (reps) comment\n");
      fprintf(stderr,"Input data format (stdin) for Gamma:               flag wcom reff veff refr refi (reps) comment\n");
      fprintf(stderr,"  flag ... Control flag composed of 4 characters\n");
      fprintf(stderr,"  wcom ... Weight of component [a.u.]\n");
      fprintf(stderr,"  rmod ... Mode radius [um]\n");
      fprintf(stderr,"  lsgm ... Sigma Log10(R [um]) or Sigma Ln(R [um])\n");
      fprintf(stderr,"  reff ... Effective radius [um]\n");
      fprintf(stderr,"  veff ... Effective variance [a.u.]\n");
      fprintf(stderr,"  inp1 ... File name of optical parameters 1 (w ea en o g)\n");
      fprintf(stderr,"  inp2 ... File name of optical parameters 2 (w a s1 s2 s3)\n");
      fprintf(stderr,"  inp3 ... File name of optical parameters 3\n");
      fprintf(stderr,"  indx ... Index of optical parameters 3\n");
      fprintf(stderr,"  fnam ... File name of size distribution\n");
      fprintf(stderr,"  refr ... File name or at mark (@) followed by single value of refractive index (real part)\n");
      fprintf(stderr,"  refi ... File name or at mark (@) followed by single value of refractive index (imaginary part)\n");
      fprintf(stderr,"  reps ... File name or at mark (@) followed by single value of aspect ratio of spheroid [a.u.]\n");
      fprintf(stderr,"  comment ... Name of component etc.\n");
      break;
    case 4:
      fprintf(stderr,"Output data format:\n");
      fprintf(stderr,"  Out. file 1: w ea en o g\n");
      fprintf(stderr,"    w  ... Wavelength [nm]\n");
      fprintf(stderr,"    ea ... Extinction coefficient (absolute) [a.u.]\n");
      fprintf(stderr,"    en ... Extinction coefficient (normalized at reference wavelength)\n");
      fprintf(stderr,"    o  ... Single scattering albedo\n");
      fprintf(stderr,"    g  ... Asymmetry parameter\n");
      fprintf(stderr,"  Note: If wcom is in particle/cm3, then ea is in 1/Mm.\n");
      fprintf(stderr,"        If wcom is in particle, then ea is in um2.\n");
      fprintf(stderr,"  Out. file 2: w a s1 s2 s3\n");
      fprintf(stderr,"    w  ... Wavelength [nm]\n");
      fprintf(stderr,"    a  ... Scattering angle [degree]\n");
      fprintf(stderr,"    s1 ... Phase function of s-polarized wave\n");
      fprintf(stderr,"    s2 ... Phase function of p-polarized wave\n");
      fprintf(stderr,"    s3 ... Phase function of averaged wave\n");
      break;
    default:
      return -1;
      break;
  }

  return 0;
}
