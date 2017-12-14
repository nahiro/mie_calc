/*********************************************************/
/* aeros_mix                                             */
/* Author: N.Manago Apr,30,2008                          */
/* $Revision: 1127 $                                      */
/* $Date: 2016-07-14 14:15:47 +0900 (Thu, 14 Jul 2016) $ */
/*********************************************************/
#include <cfortran.h>
#include <minuitfcn.h>
#include <minuit.h>
#include "aeros_comp.h"
#include "aeros_comp.c"

// Constants for Optimization
#define	OPT_NCOL			16
#define	OPT_NINT			3
#define	OPT_NDBL			13
#define	OPT_VMIE_MIN			1.0e+0
#define	OPT_VMIE_MAX			9.999999e+2
#define	OPT_VMIE_STP			1.0e+0
#define	OPT_WSCL_MIN			0.0
#define	OPT_WSCL_MAX			10.0
#define	OPT_WSCL_STP			1.0e-2
#define	OPT_NCMP			3.0
#define	OPT_NCMP_MIN			0.0
#define	OPT_NCMP_MAX			10.0
#define	OPT_NCMP_STP			1.0
#define	OPT_MIXR			0.0
#define	OPT_MIXR_MIN			-3.0
#define	OPT_MIXR_MAX			+3.0
#define	OPT_MIXR_STP			1.0e-2
#define	OPT_LMOD			0.3
#define	OPT_LMOD_MIN			-3.0
#define	OPT_LMOD_MAX			+2.0
#define	OPT_LMOD_STP			1.0e-2
#define	OPT_LSGM			0.3
#define	OPT_LSGM_MIN			0.1
#define	OPT_LSGM_MAX			0.6
#define	OPT_LSGM_STP			1.0e-2
#define	OPT_MAXVMIN			100000
#define	OPT_DSR_TOLERANCE		0.5
#define	OPT_SSR_TOLERANCE		0.1
#define	OPT_FLIM			1.0e10
#define	OPT_N_BUILTIN			3
#define	OPT_PREC			1.0e-6			// Precision
#define	OPT_OK				3			// Minuit status (CONVERGED)
#define	OPT_DONE			10000
#define	OPT_MODE_DSR_MODTRAN_0		0
#define	OPT_MODE_DSR_MODTRAN_1		1
#define	OPT_MODE_DSR_MODTRAN_2		2
#define	OPT_MODE_DSR_MODTRAN_3		3
#define	OPT_MODE_DSR_CUSTOM_0		4
#define	OPT_MODE_DSR_CUSTOM_1		5
#define	OPT_MODE_DSR_CUSTOM_2		6
#define	OPT_MODE_DSR_CUSTOM_3		7
#define	OPT_MODE_SSR_MODTRAN_S		100
#define	OPT_MODE_SSR_MODTRAN_B		101
#define	OPT_MODE_SSR_CUSTOM_S		102
#define	OPT_MODE_SSR_CUSTOM_B		103
#define	OPT_MODE_CRS_CUSTOM		1001
// Constants for Mixing
#define	MIX_LMOD_NSTP			4
#define	MIX_LSGM_NSTP			3
// Other constants
//#define	f2cFortran

// Parameters for Optimization
int	opt_mode			= -1;			// Optimization mode
int	opt_amod			= -1;			// The best fit model#
int	opt_vmie_skp			= 0;			// Skip vmie optimization
int	opt_wscl_skp			= 0;			// Skip wscl optimization
int	opt_amod_skp			= 0;			// Skip amod optimization
int	opt_ihaz[OPT_N_BUILTIN]		= {1,4,5};		// MODTRAN aerosol models
int	opt_n_builtin			= OPT_N_BUILTIN;	// #Models
int	opt_imin[OPT_MAXVMIN];					// Serial#
int	opt_flag[OPT_MAXVMIN];					// Flag
int	opt_n_pass			= 5000;			// #Values
double	opt_vmie_1st			= NAN;			// Visibility in km
double	opt_wscl_1st			= NAN;			// Water scale
double	opt_pval[OPT_NPAR];					// Optimized parameter value
double	opt_perr[OPT_NPAR];					// Optimized parameter error
double	opt_vmin[OPT_MAXVMIN];					// Minimum function value
char	*opt_command			= NULL;			// Minuit command
// Parameters for Mixing
int	mix_comp[OPT_N_COMP]		= {8,9,2};		// Mixture component#
int	mix_lmod_nstp			= MIX_LMOD_NSTP;	// Log10(R) #steps
int	mix_lsgm_nstp			= MIX_LSGM_NSTP;	// Sigma Log10(R) #steps
int	mix_lmd1_ngp			= 10;
int	mix_lmd2_ngp			= 10;
int	mix_lmd3_ngp			= 5;
int	mix_mxr2_ngp			= 10;
int	mix_mxr3_ngp			= 10;
double	mix_lmd1_min			= -1.4;
double	mix_lmd1_max			= +1.3;
double	mix_lmd2_min			= -3.0;
double	mix_lmd2_max			= +0.8;
double	mix_lmd3_min			= -2.0;
double	mix_lmd3_max			= -1.0;
double	mix_mxr2_min			= -2.0;
double	mix_mxr2_max			= +2.0;
double	mix_mxr3_min			= -2.0;
double	mix_mxr3_max			= +2.0;
double	mix_lsg1_val			= 0.31;
double	mix_lsg2_val			= 0.35;
double	mix_lsg3_val			= 0.34;

int Optimize(void);
int Compare(int *n1,int *n2);
void PrintTime(void);
int (*CRS_FCN_CUSTOM)(double *par,double *val);
int (*RES_FCN_CUSTOM)(double *par,double *val);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(cnt_hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  Optimize();
  MiePrintout();
  Printout();
  Finish();

  return 0;
}

int Init(void)
{
  int i,j,n;
  int rng;
  int rt;
  double pval[MIE_NPAR];
  double fval;

  cnt_aur_auto[1] = 1;
  if(ReadInput() < 0) return -1;
  if(cnt_aur[1] == 1)
  {
    if(dsr_n_data<1 || aur_n_data<1 || ssr_n_data<1)
    {
      fprintf(stderr,"Error, dsr_n_data=%d, aur_n_data=%d, ssr_n_data=%d\n",
                      dsr_n_data,aur_n_data,ssr_n_data);
      return -1;
    }
    CRS_FCN_CUSTOM = CRS_FCN_DSR_AUR_SSR;
    RES_FCN_CUSTOM = CRS_FCN_AUR_SSR;
  }
  else
  {
    if(dsr_n_data<1 || ssr_n_data<1)
    {
      fprintf(stderr,"Error, dsr_n_data=%d, ssr_n_data=%d\n",
                      dsr_n_data,ssr_n_data);
      return -1;
    }
    CRS_FCN_CUSTOM = CRS_FCN_DSR_SSR;
    RES_FCN_CUSTOM = SSR_FCN_CUSTOM_B1;
  }
  if(MieInit() < 0) return -1;
  if(tst_n_comp > 0)
  {
    tst_mode += 1;
  }
  if(ReadCRC() < 0) return -1;
  if((rt=ReadSimul(tst_fsim,1,1,1,0,1)) < 0) return -1;
  if(rt > 0)
  {
    tst_mode += 2;
  }
  if(upper_pfunc(50.0) < 0) return -1;

  if(tst_mode == 1)
  {
    for(n=0; n<tst_n_comp; n++)
    {
      for(i=0; i<mie_n_wlen; i++)
      {
        mie_refr_com[n][i] = tst_refr_com[n][i];
        mie_refi_com[n][i] = tst_refi_com[n][i];
      }
      if(MieTable(n,n) < 0) return -1;
      if(MieReff(n,n,tst_lsgm_com[n]) < 0) return -1;
    }
    pval[OPT_PAR_VMIE] = opt_vmie_init;
    pval[OPT_PAR_WSCL] = opt_wscl_init;
    pval[OPT_PAR_NCMP] = (double)tst_n_comp;
    pval[OPT_PAR_LMD1] = mod2eff(0,tst_lmod_com[0]);
    pval[OPT_PAR_LSG1] = tst_lsgm_com[0];
    for(i=OPT_PAR_MXR2,n=1; n<tst_n_comp; n++)
    {
      pval[i++] = log10(tst_wcom_com[0]/tst_wcom_com[n]);
      pval[i++] = mod2eff(n,tst_lmod_com[n]);
      pval[i++] = tst_lsgm_com[n];
    }
    sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
    if(cnt_dsr[0] == 1) DSR_FCN_CUSTOM_B0(pval,&fval);
    if(cnt_dsr[1] == 1) DSR_FCN_CUSTOM_B1(pval,&fval);
    if(cnt_dsr[2] == 1) DSR_FCN_CUSTOM_B2(pval,&fval);
    if(cnt_dsr[3] == 1) DSR_FCN_CUSTOM_B3(pval,&fval);
    if(cnt_aur[0] == 1) AUR_FCN_CUSTOM_S0(pval,&fval);
    if(cnt_aur[1] == 1) AUR_FCN_CUSTOM_S1(pval,&fval);
    if(cnt_ssr[0] == 1) SSR_FCN_CUSTOM_S0(pval,&fval);
    if(cnt_ssr[1] == 1) SSR_FCN_CUSTOM_S1(pval,&fval);
    if(cnt_ssr[2] == 1) SSR_FCN_CUSTOM_S2(pval,&fval);
    if(cnt_ssr[3] == 1) SSR_FCN_CUSTOM_S3(pval,&fval);
  }
  if(tst_mode) // replace data with simulation data
  {
    Printout();
    for(j=0; j<dsr_n_data; j++)
    {
      for(rng=0; rng<DSR_MAXRANG; rng++)
      {
        for(i=0; i<dsr_n_wlen[rng]; i++)
        {
          dsr_data[rng][j][i] = dsr_dsim[rng][j][i];
        }
      }
    }
    for(j=0; j<aur_n_data; j++)
    {
      for(rng=0; rng<AUR_MAXRANG; rng++)
      {
        for(i=0; i<aur_n_wlen[rng]; i++)
        {
          aur_data[rng][j][i] = aur_dsim[rng][j][i];
        }
      }
    }
    for(j=0; j<ssr_n_data; j++)
    {
      for(rng=0; rng<SSR_MAXRANG; rng++)
      {
        for(i=0; i<ssr_n_wlen[rng]; i++)
        {
          ssr_data[rng][j][i] = ssr_dsim[rng][j][i];
        }
      }
    }
  }

  opt_command = (char *)malloc(MAXLINE*sizeof(char));
  if(opt_command == NULL)
  {
    fprintf(stderr,"Error in allocating memory\n");
    return -1;
  }
  if(cnt_vb)
  {
    //MNINIT(5,6,7);
  }
  else
  {
    //MNINIT(2,2,2);
  }
  // set the print level
  if(cnt_db)
  {
    snprintf(opt_command,MAXLINE,"SET PRI 2");
  }
  else
  {
    snprintf(opt_command,MAXLINE,"SET PRI 0");
  }
  //MNCOMD(minuitfcn,opt_command,dummy,NULL);
  // set precision
  snprintf(opt_command,MAXLINE,"SET EPS %.4e",OPT_PREC);
  //MNCOMD(minuitfcn,opt_command,dummy,NULL);

  return 0;
}

int Optimize(void)
{
  int i,n;
  double fmin;
  double pval[MIE_NPAR];

  // initial guess
  if(!isnan(opt_vmie_1st))
  {
    pval[OPT_PAR_VMIE] = opt_vmie_1st;
  }
  else
  {
    pval[OPT_PAR_VMIE] = opt_vmie_init;
  }
  if(!isnan(opt_wscl_1st))
  {
    pval[OPT_PAR_WSCL] = opt_wscl_1st;
  }
  else
  {
    pval[OPT_PAR_WSCL] = opt_wscl_init;
  }
  pval[OPT_PAR_NCMP] = OPT_NCMP;
  pval[OPT_PAR_MXR2] = OPT_MIXR;
  pval[OPT_PAR_MXR3] = OPT_MIXR;
  pval[OPT_PAR_LMD1] = OPT_LMOD;
  pval[OPT_PAR_LMD2] = OPT_LMOD;
  pval[OPT_PAR_LMD3] = OPT_LMOD;
  pval[OPT_PAR_LSG1] = OPT_LSGM;
  pval[OPT_PAR_LSG2] = OPT_LSGM;
  pval[OPT_PAR_LSG3] = OPT_LSGM;
  if(crc_n_comp > 0)
  {
    for(n=0; n<crc_n_comp; n++)
    {
      for(i=0; i<mie_n_wlen; i++)
      {
        mie_refr_com[n][i] = crc_refr_com[n][i];
        mie_refi_com[n][i] = crc_refi_com[n][i];
      }
      if(MieTable(n,n) < 0) return -1;
      if(MieReff(n,n,crc_lsgm_com[n]) < 0) return -1;
    }
    pval[OPT_PAR_NCMP] = (double)crc_n_comp;
    pval[OPT_PAR_LMD1] = mod2eff(0,crc_lmod_com[0]);
    pval[OPT_PAR_LSG1] = crc_lsgm_com[0];
    for(i=OPT_PAR_MXR2,n=1; n<crc_n_comp; n++)
    {
      pval[i++] = log10(crc_wcom_com[0]/crc_wcom_com[n]);
      pval[i++] = mod2eff(n,crc_lmod_com[n]);
      pval[i++] = crc_lsgm_com[n];
    }
  }
  if(cnt_vb) PrintTime();

  // 3 components
  message(stderr,"3 components\n");
  Interp1D(cmp_wlen,cmp_refr[mix_comp[0]],cmp_n_wlen,mie_wlen_um,mie_refr_com[0],mie_n_wlen,1);
  Interp1D(cmp_wlen,cmp_refi[mix_comp[0]],cmp_n_wlen,mie_wlen_um,mie_refi_com[0],mie_n_wlen,0);
  Interp1D(cmp_wlen,cmp_refr[mix_comp[1]],cmp_n_wlen,mie_wlen_um,mie_refr_com[1],mie_n_wlen,0);
  Interp1D(cmp_wlen,cmp_refi[mix_comp[1]],cmp_n_wlen,mie_wlen_um,mie_refi_com[1],mie_n_wlen,0);
  Interp1D(cmp_wlen,cmp_refr[mix_comp[2]],cmp_n_wlen,mie_wlen_um,mie_refr_com[2],mie_n_wlen,0);
  Interp1D(cmp_wlen,cmp_refi[mix_comp[2]],cmp_n_wlen,mie_wlen_um,mie_refi_com[2],mie_n_wlen,0);
  for(n=0; n<OPT_N_COMP; n++)
  {
    if(MieTable(n,n) < 0) return -1;
    //if(MieReff(n,n,cmp_lsgm[mix_comp[n]]) < 0) return -1;
    if(MieReff(n,n,crc_lsgm_com[n]) < 0) return -1;
  }

  pval[ 0] =   5.79102539181652e+01;
  pval[ 1] =   1.07799360742185e+00;
  pval[ 2] =   3.00000000000000e+00;
  pval[ 3] =   1.70906715331877e-01;
  pval[ 4] =   3.10000000000000e-01;
  pval[ 5] =  -4.65334865520782e-01;
  pval[ 6] =  -1.10987930633307e+00;
  pval[ 7] =   3.50000000000000e-01;
  pval[ 8] =   1.24569944331153e-01;
  pval[ 9] =  -7.87460331767568e-01;
  pval[10] =   3.40000000000000e-01;
  pval[11] =   1.00000000000000e-06;
  pval[12] =  9.32762610054711e-262;
  pval[13] =  9.32771332461869e-262;
  pval[14] =  9.32762610054739e-262;
  pval[15] =  8.07332192848161e-316;
  pval[16] =   0.00000000000000e+00;
  pval[17] =   0.00000000000000e+00;
  pval[18] =   0.00000000000000e+00;
  pval[19] =  4.94065645841247e-323;
  pval[20] =  6.01946397098629e-270;
  pval[21] =  2.12914813663507e-313;
  pval[22] =  9.88131291682493e-324;
  pval[23] =   8.72735372383850e+43;
  pval[24] =  1.03977793756804e-312;
  pval[25] =  -4.85818609199095e-02;
  pval[26] =  6.02051007632726e-270;
  pval[27] =  1.58867574617256e-314;
  pval[28] =  6.01833391930700e-270;
  pval[29] =  4.46412403105209e-305;
  pval[30] =  -4.60933027506062e-39;

  // calculate DSR to set mie parameters
  DSR_FCN_CUSTOM_B1(pval,&fmin);

  fprintf(stderr,"pval[OPT_PAR_LMD1] = %22.14e\n",pval[OPT_PAR_LMD1]);
  fprintf(stderr,"lmod = %22.14e\n",pow(10.0,eff2mod(0,pval[OPT_PAR_LMD1])));
  fprintf(stderr,"pval[OPT_PAR_LMD2] = %22.14e\n",pval[OPT_PAR_LMD2]);
  fprintf(stderr,"lmod = %22.14e\n",pow(10.0,eff2mod(1,pval[OPT_PAR_LMD2])));
  fprintf(stderr,"pval[OPT_PAR_LMD3] = %22.14e\n",pval[OPT_PAR_LMD3]);
  fprintf(stderr,"lmod = %22.14e\n",pow(10.0,eff2mod(2,pval[OPT_PAR_LMD3])));

  return 0;
}

int Compare(int *n1,int *n2)
{
  double d;

  d = opt_vmin[*n1]-opt_vmin[*n2];

  return (d<0.0?-1:d>0.0?1:0);
}

void fcn(int npar,double* grad,double* fcnval,double* xval,int iflag,void* futil)
{
  switch(opt_mode)
  {
    case OPT_MODE_DSR_MODTRAN_0:
      DSR_FCN_MODTRAN_B0(xval,fcnval);
      break;
    case OPT_MODE_DSR_MODTRAN_1:
      DSR_FCN_MODTRAN_B1(xval,fcnval);
      break;
    case OPT_MODE_DSR_MODTRAN_2:
      DSR_FCN_MODTRAN_B2(xval,fcnval);
      break;
    case OPT_MODE_DSR_MODTRAN_3:
      DSR_FCN_MODTRAN_B3(xval,fcnval);
      break;
    case OPT_MODE_DSR_CUSTOM_0:
      DSR_FCN_CUSTOM_B0(xval,fcnval);
      break;
    case OPT_MODE_DSR_CUSTOM_1:
      DSR_FCN_CUSTOM_B1(xval,fcnval);
      break;
    case OPT_MODE_DSR_CUSTOM_2:
      DSR_FCN_CUSTOM_B2(xval,fcnval);
      break;
    case OPT_MODE_DSR_CUSTOM_3:
      DSR_FCN_CUSTOM_B3(xval,fcnval);
      break;
    case OPT_MODE_SSR_MODTRAN_S:
      SSR_FCN_MODTRAN_S1(xval,fcnval);
      break;
    case OPT_MODE_SSR_MODTRAN_B:
      SSR_FCN_MODTRAN_B1(xval,fcnval);
      break;
    case OPT_MODE_SSR_CUSTOM_S:
      SSR_FCN_CUSTOM_S1(xval,fcnval);
      break;
    case OPT_MODE_SSR_CUSTOM_B:
      SSR_FCN_CUSTOM_B1(xval,fcnval);
      break;
    case OPT_MODE_CRS_CUSTOM:
      CRS_FCN_CUSTOM(xval,fcnval);
      break;
    default:
      message(stderr,"Invalid optimizing mode >>> %d\n",opt_mode);
      break;
  }

  return;
}

int Printout(void)
{
  int i;

  CommonPrintout(stderr,1,1,1);
  for(i=0; i<OPT_NPAR; i++)
  {
    printf("%2d %18.10e %18.10e\n",i,opt_pval[i],opt_perr[i]);
  }
  fflush(stdout);
  fflush(stderr);

  return 0;
}

int ReadOwnConfig(char *s)
{
  int n,nc;
  int err;
  char temp[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char *p;

  err = 0;
  do
  {
    strncpy(temp,s,MAXLINE);
    if((p=strchr(temp,'#')) != NULL) *p = '\0';
    // read line
    for(n=nc=0,p=temp; n<MAXENTR; n++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
    }
    if(n < 1) break;
    if(strcasecmp(str[0],"opt_vmie_1st") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        opt_vmie_1st = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n","opt_vmie_1st",opt_vmie_1st);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_wscl_1st") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        opt_wscl_1st = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n","opt_wscl_1st",opt_wscl_1st);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_vmie_skp") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        opt_vmie_skp = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30d\n","opt_vmie_skp",opt_vmie_skp);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_wscl_skp") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        opt_wscl_skp = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30d\n","opt_wscl_skp",opt_wscl_skp);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_amod_skp") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        opt_amod_skp = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30d\n","opt_amod_skp",opt_amod_skp);
        cnt_n_cmnt++;
      }
    }
  }
  while(0);
  if(err)
  {
    return -1;
  }

  return 0;
}

void PrintTime(void)
{
  time_t t;
  char line[MAXLINE];

  time(&t);
  strftime(line,MAXLINE,"%Y-%m-%d %H:%M:%S",localtime(&t));
  message(stderr,"%s\n",line);
}

int GetOwnOpt(int argn,char **args)
{
  int c,rt;
  int ntmp;
  int option_index = 0;
  int this_option_optind;
  char *endp;
  struct option long_options[] =
  {
    {"mix_lmod_nstp",1,0,'l'},
    {"mix_lsgm_nstp",1,0,'L'},
  };

  fprintf(stderr,"%15s : $Revision: 1127 $ $Date: 2016-07-14 14:15:47 +0900 (Thu, 14 Jul 2016) $\n","aeros_mix.c");

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-l:L:",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'l':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) mix_lmod_nstp = ntmp;
        else
        {
          fprintf(stderr,"Log10(R) #steps for opt -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'L':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) mix_lsgm_nstp = ntmp;
        else
        {
          fprintf(stderr,"Sigma Log10(R) #steps for opt -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case '?':
        if(cnt_optn >= CNT_MAXNOPT)
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
        if(cnt_optn >= CNT_MAXNOPT)
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
  int  n = 15;
  char e[MAXLINE];
  char a[MAXLINE];
  char d[MAXLINE];
  char optstring[] = "DNteFlLVWAaGgwmMrxXnsSCcodvh";

  cnt_n_own_format = 6;
  sprintf(cnt_own_format[1],"opt_vmie_1st  vmie         | visibility\n");
  sprintf(cnt_own_format[2],"opt_wscl_1st  wscl         | water scaling factor\n");
  sprintf(cnt_own_format[3],"opt_vmie_skp  flag         | skip vmie optimization(%d)\n",0);
  sprintf(cnt_own_format[4],"opt_wscl_skp  flag         | skip wscl optimization(%d)\n",0);
  sprintf(cnt_own_format[5],"opt_amod_skp  flag         | skip amod optimization(%d)\n",0);

  fprintf(stderr,"aeros_mix ... fit DSR/SSR spectrum.\n");
  CommonUsage("aeros_mix",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'l':
        fprintf(stderr," l -mix_lmod_nstp |%s|%s|%s| %d\n",As(e,"Log10(R) #steps",n),As(a,"#",n),          Ad(d,MIX_LMOD_NSTP,n),mix_lmod_nstp);
        break;
      case 'L':
        fprintf(stderr," L -mix_lsgm_nstp |%s|%s|%s| %d\n",As(e,"Sigma    #steps",n),As(a,"#",n),          Ad(d,MIX_LSGM_NSTP,n),mix_lsgm_nstp);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);

  return 0;
}
