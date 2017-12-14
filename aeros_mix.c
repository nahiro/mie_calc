/*********************************************************/
/* aeros_mix                                             */
/* Author: N.Manago Apr,30,2008                          */
/* $Revision: 1127 $                                      */
/* $Date: 2016-07-14 14:15:47 +0900 (Thu, 14 Jul 2016) $ */
/*********************************************************/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cfortran.h>
#include <minuitfcn.h>
#include <minuit.h>
#include "aeros_comp.h"
#include "aeros_comp.c"

// Constants for Optimization
#define	OPT_NCOL			16
#define	OPT_NINT			3
#define	OPT_NDBL			13
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
#define	OPT_VTAU_MIN			0.0
#define	OPT_VTAU_MAX			10.0
#define	OPT_VTAU_STP			1.0e-4
#define	OPT_VMIE_MIN			1.0e+0
#define	OPT_VMIE_MAX			1.0e3
#define	OPT_VMIE_STP			1.0e+0
#define	OPT_WSCL_MIN			0.0
#define	OPT_WSCL_MAX			10.0
#define	OPT_WSCL_STP			1.0e-2
#define	OPT_OSCL_MIN			0.0
#define	OPT_OSCL_MAX			2.0
#define	OPT_OSCL_STP			1.0e-2
#define	OPT_XSCL_MIN			0.0
#define	OPT_XSCL_MAX			10.0
#define	OPT_XSCL_STP			1.0e-2
#define	OPT_YSCL_MIN			0.0
#define	OPT_YSCL_MAX			10.0
#define	OPT_YSCL_STP			1.0e-2
#define	OPT_ZSCL_MIN			0.0
#define	OPT_ZSCL_MAX			10.0
#define	OPT_ZSCL_STP			1.0e-2
#define	OPT_NCMP_MIN			0.0
#define	OPT_NCMP_MAX			10.0
#define	OPT_NCMP_STP			1.0
#define	OPT_MIXR_MIN			-3.0
#define	OPT_MIXR_MAX			+3.0
#define	OPT_MIXR_STP			1.0e-2
#define	OPT_LMOD_MIN			-3.0
#define	OPT_LMOD_MAX			+2.0
#define	OPT_LMOD_STP			1.0e-2
#define	OPT_LSGM_MIN			0.1
#define	OPT_LSGM_MAX			0.6
#define	OPT_LSGM_STP			1.0e-2
// Other constants
//#define	f2cFortran

// Parameters for Optimization
int	opt_test			= 0;			// Test mode
int	opt_mode			= -1;			// Optimization mode
int	opt_amod			= -1;			// The best fit model#
int	opt_amod_skp			= 0;			// Skip amod optimization
int	opt_ihaz[OPT_N_BUILTIN]		= {1,4,5};		// MODTRAN aerosol models
int	opt_n_builtin			= OPT_N_BUILTIN;	// #Models
int	opt_imin[OPT_MAXVMIN];					// Serial#
int	opt_flag[OPT_MAXVMIN];					// Flag
int	opt_n_pass			= 5000;			// #Values
int	opt_skip[OPT_NPAR];					// Skip optimization flag
int	opt_stat[OPT_NPAR];					// Optimization status
double	opt_vmin[OPT_MAXVMIN];					// Minimum function value
double	opt_pval[OPT_NPAR];					// Optimized parameter value
double	opt_perr[OPT_NPAR];					// Optimized parameter error
double	opt_fmin[OPT_NPAR];					// Optimized function value
double	opt_fedm[OPT_NPAR];					// Estimated vertical distance
char	*opt_command			= NULL;			// Minuit command
int	opt_comp[OPT_N_COMP]		= {8,3,2};		// Mixture component#
int	opt_lmod_ngp[OPT_N_COMP]	= {8,8,8};
int	opt_mixr_ngp[OPT_N_COMP]	= {0,6,6};
double	opt_lmod_min[OPT_N_COMP]	= {-2.0,-2.5,-2.5};
double	opt_lmod_max[OPT_N_COMP]	= { 0.8, 0.8, 0.0};
double	opt_mixr_min[OPT_N_COMP]	= { 0.0,-1.5,-0.5};
double	opt_mixr_max[OPT_N_COMP]	= { 0.0, 1.5, 2.5};
// Parameters for Test
double	tst_dsr_ernd			= NAN;
double	tst_dsr_esys			= NAN;
double	tst_aur_ernd			= NAN;
double	tst_aur_esys			= NAN;
double	tst_ssr_ernd			= NAN;
double	tst_ssr_esys			= NAN;
unsigned long tst_seed			= 0;

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
  int dummy;
  double pval[SIM_NPAR];
  double fval;
  struct timeval  tv;
  struct timezone tz;
  const gsl_rng_type *typ;
  gsl_rng *gen;

  switch(sim_mode)
  {
    case SIM_MODTRAN_VTAU:
      opt_skip[OPT_PAR_VMIE] = 1;
      break;
    case SIM_MODTRAN_VMIE:
      opt_skip[OPT_PAR_VTAU] = 1;
      break;
    default:
      fprintf(stderr,"Invalid simulation mode >>> %d\n",sim_mode);
      return -1;
      break;
  }
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
    for(i=0; i<OPT_PAR_NCMP; i++)
    {
      pval[i] = opt_init[i];
    }
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
  } else
  if(opt_test)
  {
    for(i=0; i<OPT_PAR_NCMP; i++)
    {
      pval[i] = opt_init[i];
    }
    pval[OPT_PAR_NCMP] = 0.0;
    sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
    if(cnt_dsr[0] == 1) DSR_FCN_MODTRAN_B0(pval,&fval);
    if(cnt_dsr[1] == 1) DSR_FCN_MODTRAN_B1(pval,&fval);
    if(cnt_dsr[2] == 1) DSR_FCN_MODTRAN_B2(pval,&fval);
    if(cnt_dsr[3] == 1) DSR_FCN_MODTRAN_B3(pval,&fval);
    if(cnt_aur[0] == 1) AUR_FCN_MODTRAN_S0(pval,&fval);
    if(cnt_aur[1] == 1) AUR_FCN_MODTRAN_S1(pval,&fval);
    if(cnt_ssr[0] == 1) SSR_FCN_MODTRAN_S0(pval,&fval);
    if(cnt_ssr[1] == 1) SSR_FCN_MODTRAN_S1(pval,&fval);
    if(cnt_ssr[2] == 1) SSR_FCN_MODTRAN_S2(pval,&fval);
    if(cnt_ssr[3] == 1) SSR_FCN_MODTRAN_S3(pval,&fval);
    tst_mode = 4;
  }
  if(tst_mode) // replace data with simulation data
  {
    if(cnt_vb > 0)
    {
      message(stderr,"Before replacement:\n");
    }
    Printout();
    if((!isnan(tst_dsr_ernd)) || (!isnan(tst_aur_ernd)) || (!isnan(tst_ssr_ernd)))
    {
      gsl_rng_env_setup();
      typ = gsl_rng_default;
      gen = gsl_rng_alloc(typ);
      if(tst_seed == 0)
      {
        gettimeofday(&tv,&tz);
        tst_seed = tv.tv_sec+tv.tv_usec;
      }
      if(cnt_vb > 0)
      {
        message(stderr,"tst_seed : %18ld\n",tst_seed);
      }
      gsl_rng_set(gen,tst_seed);
    }
    // DSR
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
    // add random noise
    if(!isnan(tst_dsr_ernd))
    {
      for(j=0; j<dsr_n_data; j++)
      {
        for(rng=0; rng<DSR_MAXRANG; rng++)
        {
          for(i=0; i<dsr_n_wlen[rng]; i++)
          {
            dsr_data[rng][j][i] += gsl_ran_gaussian(gen,tst_dsr_ernd*dsr_dsim[rng][j][i]);
          }
        }
      }
    }
    // add systematic noise
    if(!isnan(tst_dsr_esys))
    {
      for(j=0; j<dsr_n_data; j++)
      {
        for(rng=0; rng<DSR_MAXRANG; rng++)
        {
          for(i=0; i<dsr_n_wlen[rng]; i++)
          {
            dsr_data[rng][j][i] += tst_dsr_esys*dsr_dsim[rng][j][i];
          }
        }
      }
    }
    // AUR
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
    // add random noise
    if(!isnan(tst_aur_ernd))
    {
      for(j=0; j<aur_n_data; j++)
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          for(i=0; i<aur_n_wlen[rng]; i++)
          {
            aur_data[rng][j][i] += gsl_ran_gaussian(gen,tst_aur_ernd*aur_dsim[rng][j][i]);
          }
        }
      }
    }
    // add systematic noise
    if(!isnan(tst_aur_esys))
    {
      for(j=0; j<aur_n_data; j++)
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          for(i=0; i<aur_n_wlen[rng]; i++)
          {
            aur_data[rng][j][i] += tst_aur_esys*aur_dsim[rng][j][i];
          }
        }
      }
    }
    // SSR
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
    // add random noise
    if(!isnan(tst_ssr_ernd))
    {
      for(j=0; j<ssr_n_data; j++)
      {
        for(rng=0; rng<SSR_MAXRANG; rng++)
        {
          for(i=0; i<ssr_n_wlen[rng]; i++)
          {
            ssr_data[rng][j][i] += gsl_ran_gaussian(gen,tst_ssr_ernd*ssr_dsim[rng][j][i]);
          }
        }
      }
    }
    // add systematic noise
    if(!isnan(tst_ssr_esys))
    {
      for(j=0; j<ssr_n_data; j++)
      {
        for(rng=0; rng<SSR_MAXRANG; rng++)
        {
          for(i=0; i<ssr_n_wlen[rng]; i++)
          {
            ssr_data[rng][j][i] += tst_ssr_esys*ssr_dsim[rng][j][i];
          }
        }
      }
    }
    if((!isnan(tst_dsr_ernd)) || (!isnan(tst_aur_ernd)) || (!isnan(tst_ssr_ernd)))
    {
      gsl_rng_free(gen);
    }
    if(cnt_vb > 0)
    {
      message(stderr,"After replacement:\n");
      Printout();
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
    MNINIT(5,6,7);
  }
  else
  {
    MNINIT(2,2,2);
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
  MNCOMD(minuitfcn,opt_command,dummy,NULL);
  // set precision
  snprintf(opt_command,MAXLINE,"SET EPS %.4e",OPT_PREC);
  MNCOMD(minuitfcn,opt_command,dummy,NULL);

  return 0;
}

int Optimize(void)
{
  int i,j,k,n;
  int ngrd;
  int i_lmd1,i_lmd2,i_lmd3;
  int i_mxr2,i_mxr3;
  int flag;
  int dummy;
  int npari,nparx,istat;
  double opt_leff_min[OPT_N_COMP];
  double opt_leff_max[OPT_N_COMP];
  double opt_leff_stp[OPT_N_COMP];
  double opt_mixr_stp[OPT_N_COMP];
  double vmin;
  double fmin,fedm,errdef;
  double ptmp[SIM_NPAR];
  double pval[SIM_NPAR];
  double pstp[OPT_NPAR];
  double pmin[OPT_NPAR];
  double pmax[OPT_NPAR];
  char   pnam[OPT_NPAR][MAXLINE] = {"vtau","vmie","wscl","oscl",
                                    "xscl","yscl","zscl",
                                    "ncmp","lmd1","lsg1",
                                    "mxr2","lmd2","lsg2",
                                    "mxr3","lmd3","lsg3"};
  int    fitnum[OPT_NPAR];
  double fitval[OPT_NPAR];
  double fiterr[OPT_NPAR];
  double fitmin[OPT_NPAR];
  double fitmax[OPT_NPAR];
  char   fitnam[OPT_NPAR][MAXLINE];

  // Initialization
  // Minimum
  pmin[OPT_PAR_VTAU] = OPT_VTAU_MIN;
  pmin[OPT_PAR_VMIE] = OPT_VMIE_MIN;
  pmin[OPT_PAR_WSCL] = OPT_WSCL_MIN;
  pmin[OPT_PAR_OSCL] = OPT_OSCL_MIN;
  pmin[OPT_PAR_XSCL] = OPT_XSCL_MIN;
  pmin[OPT_PAR_YSCL] = OPT_YSCL_MIN;
  pmin[OPT_PAR_ZSCL] = OPT_ZSCL_MIN;
  pmin[OPT_PAR_NCMP] = OPT_NCMP_MIN;
  pmin[OPT_PAR_LMD1] = OPT_LMOD_MIN;
  pmin[OPT_PAR_LSG1] = OPT_LSGM_MIN;
  pmin[OPT_PAR_MXR2] = OPT_MIXR_MIN;
  pmin[OPT_PAR_LMD2] = OPT_LMOD_MIN;
  pmin[OPT_PAR_LSG2] = OPT_LSGM_MIN;
  pmin[OPT_PAR_MXR3] = OPT_MIXR_MIN;
  pmin[OPT_PAR_LMD3] = OPT_LMOD_MIN;
  pmin[OPT_PAR_LSG3] = OPT_LSGM_MIN;
  // Maximum
  pmax[OPT_PAR_VTAU] = OPT_VTAU_MAX;
  pmax[OPT_PAR_VMIE] = OPT_VMIE_MAX;
  pmax[OPT_PAR_WSCL] = OPT_WSCL_MAX;
  pmax[OPT_PAR_OSCL] = OPT_OSCL_MAX;
  pmax[OPT_PAR_XSCL] = OPT_XSCL_MAX;
  pmax[OPT_PAR_YSCL] = OPT_YSCL_MAX;
  pmax[OPT_PAR_ZSCL] = OPT_ZSCL_MAX;
  pmax[OPT_PAR_NCMP] = OPT_NCMP_MAX;
  pmax[OPT_PAR_LMD1] = OPT_LMOD_MAX;
  pmax[OPT_PAR_LSG1] = OPT_LSGM_MAX;
  pmax[OPT_PAR_MXR2] = OPT_MIXR_MAX;
  pmax[OPT_PAR_LMD2] = OPT_LMOD_MAX;
  pmax[OPT_PAR_LSG2] = OPT_LSGM_MAX;
  pmax[OPT_PAR_MXR3] = OPT_MIXR_MAX;
  pmax[OPT_PAR_LMD3] = OPT_LMOD_MAX;
  pmax[OPT_PAR_LSG3] = OPT_LSGM_MAX;
  // Step
  pstp[OPT_PAR_VTAU] = OPT_VTAU_STP;
  pstp[OPT_PAR_VMIE] = OPT_VMIE_STP;
  pstp[OPT_PAR_WSCL] = OPT_WSCL_STP;
  pstp[OPT_PAR_OSCL] = OPT_OSCL_STP;
  pstp[OPT_PAR_XSCL] = OPT_XSCL_STP;
  pstp[OPT_PAR_YSCL] = OPT_YSCL_STP;
  pstp[OPT_PAR_ZSCL] = OPT_ZSCL_STP;
  pstp[OPT_PAR_NCMP] = OPT_NCMP_STP;
  pstp[OPT_PAR_LMD1] = OPT_LMOD_STP;
  pstp[OPT_PAR_LSG1] = OPT_LSGM_STP;
  pstp[OPT_PAR_MXR2] = OPT_MIXR_STP;
  pstp[OPT_PAR_LMD2] = OPT_LMOD_STP;
  pstp[OPT_PAR_LSG2] = OPT_LSGM_STP;
  pstp[OPT_PAR_MXR3] = OPT_MIXR_STP;
  pstp[OPT_PAR_LMD3] = OPT_LMOD_STP;
  pstp[OPT_PAR_LSG3] = OPT_LSGM_STP;
  // Value
  for(i=0; i<SIM_NPAR; i++)
  {
    if(!isnan(crc_init[i]))
    {
      pval[i] = crc_init[i];
    }
    else
    {
      pval[i] = opt_init[i];
    }
  }
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
    flag = 1; // Custom model
  }
  else
  {
    pval[OPT_PAR_NCMP] = 0.0;
    flag = 0; // MODTRAN model
  }
  if(cnt_vb) PrintTime();

  // optimize aerosol optical depth
  if(opt_skip[OPT_PAR_VTAU] != 1)
  {
    message(stderr,"Optimizing aerosol optical depth\n");
    if(flag == 1)
    {
      opt_mode = OPT_MODE_DSR_CUSTOM_0;
    }
    else
    {
      opt_mode = OPT_MODE_DSR_MODTRAN_0;
    }
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_VTAU)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
    MNPOUT(OPT_PAR_VTAU+1,fitnam[OPT_PAR_VTAU],fitval[OPT_PAR_VTAU],fiterr[OPT_PAR_VTAU],
                          fitmin[OPT_PAR_VTAU],fitmax[OPT_PAR_VTAU],fitnum[OPT_PAR_VTAU]);
    pval[OPT_PAR_VTAU] = fitval[OPT_PAR_VTAU];
  }
  // optimize visibility
  if(opt_skip[OPT_PAR_VMIE] != 1)
  {
    message(stderr,"Optimizing visibility\n");
    if(flag == 1)
    {
      opt_mode = OPT_MODE_DSR_CUSTOM_0;
    }
    else
    {
      opt_mode = OPT_MODE_DSR_MODTRAN_0;
    }
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_VMIE)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
    MNPOUT(OPT_PAR_VMIE+1,fitnam[OPT_PAR_VMIE],fitval[OPT_PAR_VMIE],fiterr[OPT_PAR_VMIE],
                          fitmin[OPT_PAR_VMIE],fitmax[OPT_PAR_VMIE],fitnum[OPT_PAR_VMIE]);
    pval[OPT_PAR_VMIE] = fitval[OPT_PAR_VMIE];
  }
  // optimize water scaling factor
  if(opt_skip[OPT_PAR_WSCL] != 1)
  {
    message(stderr,"Optimizing water scaling factor\n");
    if(flag == 1)
    {
      opt_mode = OPT_MODE_DSR_CUSTOM_2;
    }
    else
    {
      opt_mode = OPT_MODE_DSR_MODTRAN_2;
    }
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_WSCL)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
    MNPOUT(OPT_PAR_WSCL+1,fitnam[OPT_PAR_WSCL],fitval[OPT_PAR_WSCL],fiterr[OPT_PAR_WSCL],
                          fitmin[OPT_PAR_WSCL],fitmax[OPT_PAR_WSCL],fitnum[OPT_PAR_WSCL]);
    pval[OPT_PAR_WSCL] = fitval[OPT_PAR_WSCL];
  }
  // optimize aerosol model
  if(opt_amod_skp != 1)
  {
    message(stderr,"Optimizing aerosol model\n");
    sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
    vmin = 1.0e100;
    opt_amod = -1;
    if(flag == 1)
    {
      opt_mode = OPT_MODE_DSR_CUSTOM_1;
      for(i=0; i<opt_n_builtin; i++)
      {
        sim_ihaz = opt_ihaz[i];
        DSR_FCN_CUSTOM_B1(pval,&opt_vmin[i]);
        if(opt_vmin[i] < vmin)
        {
          vmin = opt_vmin[i];
          opt_amod = i;
        }
      }
    }
    else
    {
      opt_mode = OPT_MODE_DSR_MODTRAN_1;
      for(i=0; i<opt_n_builtin; i++)
      {
        sim_ihaz = opt_ihaz[i];
        DSR_FCN_MODTRAN_B1(pval,&opt_vmin[i]);
        if(opt_vmin[i] < vmin)
        {
          vmin = opt_vmin[i];
          opt_amod = i;
        }
      }
    }
    message(stderr,"Model optimization summary\n");
    for(i=0; i<opt_n_builtin; i++)
    {
      message(stderr,"vmin[%2d]: %18.10e\n",i,opt_vmin[i]);
    }
    message(stderr,"opt_amod : %18d\n",opt_amod);
    sim_ihaz = opt_ihaz[opt_amod];
  }
  if(cnt_vb > 0)
  {
    message(stderr,"sim_ihaz : %18d\n",sim_ihaz);
  }
  // optimize aerosol optical depth
  if(opt_skip[OPT_PAR_VTAU] != 1)
  {
    message(stderr,"Optimizing aerosol optical depth\n");
    if(flag == 1)
    {
      opt_mode = OPT_MODE_DSR_CUSTOM_0;
    }
    else
    {
      opt_mode = OPT_MODE_DSR_MODTRAN_0;
    }
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_VTAU)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
    MNPOUT(OPT_PAR_VTAU+1,fitnam[OPT_PAR_VTAU],fitval[OPT_PAR_VTAU],fiterr[OPT_PAR_VTAU],
                          fitmin[OPT_PAR_VTAU],fitmax[OPT_PAR_VTAU],fitnum[OPT_PAR_VTAU]);
    pval[OPT_PAR_VTAU] = fitval[OPT_PAR_VTAU];
  }
  // optimize visibility
  if(opt_skip[OPT_PAR_VMIE] != 1)
  {
    message(stderr,"Optimizing visibility\n");
    if(flag == 1)
    {
      opt_mode = OPT_MODE_DSR_CUSTOM_0;
    }
    else
    {
      opt_mode = OPT_MODE_DSR_MODTRAN_0;
    }
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_VMIE)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
    MNPOUT(OPT_PAR_VMIE+1,fitnam[OPT_PAR_VMIE],fitval[OPT_PAR_VMIE],fiterr[OPT_PAR_VMIE],
                          fitmin[OPT_PAR_VMIE],fitmax[OPT_PAR_VMIE],fitnum[OPT_PAR_VMIE]);
    pval[OPT_PAR_VMIE] = fitval[OPT_PAR_VMIE];
  }
  // optimize water scaling factor
  if(opt_skip[OPT_PAR_WSCL] != 1)
  {
    message(stderr,"Optimizing water scaling factor\n");
    if(flag == 1)
    {
      opt_mode = OPT_MODE_DSR_CUSTOM_2;
    }
    else
    {
      opt_mode = OPT_MODE_DSR_MODTRAN_2;
    }
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_WSCL)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
    MNPOUT(OPT_PAR_WSCL+1,fitnam[OPT_PAR_WSCL],fitval[OPT_PAR_WSCL],fiterr[OPT_PAR_WSCL],
                          fitmin[OPT_PAR_WSCL],fitmax[OPT_PAR_WSCL],fitnum[OPT_PAR_WSCL]);
    pval[OPT_PAR_WSCL] = fitval[OPT_PAR_WSCL];
  }
  message(stderr,"Initial guess\n");
  if(opt_skip[OPT_PAR_VTAU] != 1)
  {
    message(stderr,"vtau: %18.10e\n",pval[OPT_PAR_VTAU]);
  }
  if(opt_skip[OPT_PAR_VMIE] != 1)
  {
    message(stderr,"vmie: %18.10e\n",pval[OPT_PAR_VMIE]);
  }
  if(opt_skip[OPT_PAR_WSCL] != 1)
  {
    message(stderr,"wscl: %18.10e\n",pval[OPT_PAR_WSCL]);
  }

  // calculate correction factors
  if(flag == 1)
  {
    if(CalcFactor(pval,1,1) < 0) return -1;
  }
  else
  {
    if(CalcFactor(pval,0,1) < 0) return -1;
  }
  sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
  sim_cflg = 1; // Apply DIS=f->t correction factor from here
  sim_dflg = 0; // Set DIS=f from here
  if(opt_test == 0)
  {
    if(flag == 1)
    {
      if(Background(pval,2) < 0) return -1;
      if(RingCenter(pval,2) < 0) return -1;
      if(Uniformity(pval,2) < 0) return -1;
    }
    else
    {
      if(Background(pval,0) < 0) return -1;
      if(RingCenter(pval,0) < 0) return -1;
      if(Uniformity(pval,0) < 0) return -1;
    }
  }

  // 3 components
  message(stderr,"3 components\n");
  if(opt_n_comp > 0)
  {
    for(n=0; n<opt_n_comp; n++)
    {
      for(i=0; i<mie_n_wlen; i++)
      {
        mie_refr_com[n][i] = opt_refr_com[n][i];
        mie_refi_com[n][i] = opt_refi_com[n][i];
      }
      if(MieTable(n,n) < 0) return -1;
      if(MieReff(n,n,opt_lsgm_com[n]) < 0) return -1;
    }
    pval[OPT_PAR_NCMP] = (double)opt_n_comp;
    pval[OPT_PAR_LMD1] = mod2eff(0,opt_lmod_com[0]);
    pval[OPT_PAR_LSG1] = opt_lsgm_com[0];
    for(i=OPT_PAR_MXR2,n=1; n<opt_n_comp; n++)
    {
      pval[i++] = log10(opt_wcom_com[0]/opt_wcom_com[n]);
      pval[i++] = mod2eff(n,opt_lmod_com[n]);
      pval[i++] = opt_lsgm_com[n];
    }
  }
  else
  {
    Interp1D(cmp_wlen,cmp_refr[opt_comp[0]],cmp_n_wlen,mie_wlen_um,mie_refr_com[0],mie_n_wlen,1);
    Interp1D(cmp_wlen,cmp_refi[opt_comp[0]],cmp_n_wlen,mie_wlen_um,mie_refi_com[0],mie_n_wlen,0);
    Interp1D(cmp_wlen,cmp_refr[opt_comp[1]],cmp_n_wlen,mie_wlen_um,mie_refr_com[1],mie_n_wlen,0);
    Interp1D(cmp_wlen,cmp_refi[opt_comp[1]],cmp_n_wlen,mie_wlen_um,mie_refi_com[1],mie_n_wlen,0);
    Interp1D(cmp_wlen,cmp_refr[opt_comp[2]],cmp_n_wlen,mie_wlen_um,mie_refr_com[2],mie_n_wlen,0);
    Interp1D(cmp_wlen,cmp_refi[opt_comp[2]],cmp_n_wlen,mie_wlen_um,mie_refi_com[2],mie_n_wlen,0);
    for(n=0; n<OPT_N_COMP; n++)
    {
      if(MieTable(n,n) < 0) return -1;
      if(MieReff(n,n,cmp_lsgm[opt_comp[n]]) < 0) return -1;
    }
    pval[OPT_PAR_NCMP] = OPT_NCMP;
    pval[OPT_PAR_LMD1] = mod2eff(0,cmp_lmod[opt_comp[0]]);
    pval[OPT_PAR_LSG1] = cmp_lsgm[opt_comp[0]];
    pval[OPT_PAR_MXR2] = OPT_MIXR; // to be modified for fixed case?
    pval[OPT_PAR_LMD2] = mod2eff(1,cmp_lmod[opt_comp[1]]);
    pval[OPT_PAR_LSG2] = cmp_lsgm[opt_comp[1]];
    pval[OPT_PAR_MXR3] = OPT_MIXR; // to be modified for fixed case?
    pval[OPT_PAR_LMD3] = mod2eff(2,cmp_lmod[opt_comp[2]]);
    pval[OPT_PAR_LSG3] = cmp_lsgm[opt_comp[2]];
  }

  // perform a grid seach to find the initial values
  if(opt_skip[OPT_PAR_LMD1] > 0)
  {
    opt_lmod_ngp[0] = 1;
    opt_lmod_min[0] = pval[OPT_PAR_LMD1];
    opt_lmod_max[0] = pval[OPT_PAR_LMD1];
  }
  if(opt_skip[OPT_PAR_LMD2] > 0)
  {
    opt_lmod_ngp[1] = 1;
    opt_lmod_min[1] = pval[OPT_PAR_LMD2];
    opt_lmod_max[1] = pval[OPT_PAR_LMD2];
  }
  if(opt_skip[OPT_PAR_LMD3] > 0)
  {
    opt_lmod_ngp[2] = 1;
    opt_lmod_min[2] = pval[OPT_PAR_LMD3];
    opt_lmod_max[2] = pval[OPT_PAR_LMD3];
  }
  if(opt_skip[OPT_PAR_MXR2] > 0)
  {
    opt_mixr_ngp[1] = 1;
    opt_mixr_min[1] = pval[OPT_PAR_MXR2];
    opt_mixr_max[1] = pval[OPT_PAR_MXR2];
  }
  if(opt_skip[OPT_PAR_MXR3] > 0)
  {
    opt_mixr_ngp[2] = 1;
    opt_mixr_min[2] = pval[OPT_PAR_MXR3];
    opt_mixr_max[2] = pval[OPT_PAR_MXR3];
  }
  for(i=0; i<OPT_N_COMP; i++)
  {
    opt_leff_min[i] = mod2eff(i,opt_lmod_min[i]);
    opt_leff_max[i] = mod2eff(i,opt_lmod_max[i]);
    opt_leff_stp[i] = (opt_lmod_ngp[i]==1?0.0:(opt_leff_max[i]-opt_leff_min[i])/(opt_lmod_ngp[i]-1));
    opt_mixr_stp[i] = (opt_mixr_ngp[i]==1?0.0:(opt_mixr_max[i]-opt_mixr_min[i])/(opt_mixr_ngp[i]-1));
  }
  // DSR
  ngrd = 0;
  for(i_lmd1=0; i_lmd1<opt_lmod_ngp[0]; i_lmd1++)
  {
    pval[OPT_PAR_LMD1] = opt_leff_min[0]+opt_leff_stp[0]*i_lmd1;
    for(i_lmd2=0; i_lmd2<opt_lmod_ngp[1]; i_lmd2++)
    {
      pval[OPT_PAR_LMD2] = opt_leff_min[1]+opt_leff_stp[1]*i_lmd2;
      for(i_lmd3=0; i_lmd3<opt_lmod_ngp[2]; i_lmd3++)
      {
        pval[OPT_PAR_LMD3] = opt_leff_min[2]+opt_leff_stp[2]*i_lmd3;
        for(i_mxr2=0; i_mxr2<opt_mixr_ngp[1]; i_mxr2++)
        {
          pval[OPT_PAR_MXR2] = opt_mixr_min[1]+opt_mixr_stp[1]*i_mxr2;
          for(i_mxr3=0; i_mxr3<opt_mixr_ngp[2]; i_mxr3++)
          {
            pval[OPT_PAR_MXR3] = opt_mixr_min[2]+opt_mixr_stp[2]*i_mxr3;
            if(ngrd >= OPT_MAXVMIN)
            {
              message(stderr,"Error: #grids exceed the limit >>> %d\n",ngrd);
              return -1;
            }
            DSR_FCN_CUSTOM_B1(pval,&fmin);
            opt_vmin[ngrd] = fmin;
            opt_imin[ngrd] = ngrd;
            ngrd++;
          }
        }
      }
    }
  }
  qsort(opt_imin,ngrd,sizeof(int),(int(*)(const void*,const void*))Compare);
  for(i=0; i<ngrd; i++)
  {
    opt_flag[i] = 0;
  }
  for(i=0; i<opt_n_pass; i++)
  {
    opt_flag[opt_imin[i]] = 1;
  }
  if(cnt_vb > 0)
  {
    message(stderr,"DSR grid search results:\n");
    for(i=0; i<opt_n_pass; i++)
    {
      message(stderr,"%5d %5d %13.4e%s",i,opt_imin[i],opt_vmin[opt_imin[i]],
                      (i==opt_n_pass-1?"\n":i%8==7?"\n":" "));
    }
  }
  // SSR
  ngrd = 0;
  vmin = 1.0e100;
  for(i_lmd1=0; i_lmd1<opt_lmod_ngp[0]; i_lmd1++)
  {
    pval[OPT_PAR_LMD1] = opt_leff_min[0]+opt_leff_stp[0]*i_lmd1;
    for(i_lmd2=0; i_lmd2<opt_lmod_ngp[1]; i_lmd2++)
    {
      pval[OPT_PAR_LMD2] = opt_leff_min[1]+opt_leff_stp[1]*i_lmd2;
      for(i_lmd3=0; i_lmd3<opt_lmod_ngp[2]; i_lmd3++)
      {
        pval[OPT_PAR_LMD3] = opt_leff_min[2]+opt_leff_stp[2]*i_lmd3;
        for(i_mxr2=0; i_mxr2<opt_mixr_ngp[1]; i_mxr2++)
        {
          pval[OPT_PAR_MXR2] = opt_mixr_min[1]+opt_mixr_stp[1]*i_mxr2;
          for(i_mxr3=0; i_mxr3<opt_mixr_ngp[2]; i_mxr3++)
          {
            pval[OPT_PAR_MXR3] = opt_mixr_min[2]+opt_mixr_stp[2]*i_mxr3;
            if(opt_flag[ngrd] == 1)
            {
              RES_FCN_CUSTOM(pval,&fmin);
              if(opt_vmin[ngrd] > 1.0)
              {
                fmin *= opt_vmin[ngrd];
              }
              if(fmin < vmin)
              {
                ptmp[OPT_PAR_LMD1] = pval[OPT_PAR_LMD1];
                ptmp[OPT_PAR_LMD2] = pval[OPT_PAR_LMD2];
                ptmp[OPT_PAR_LMD3] = pval[OPT_PAR_LMD3];
                ptmp[OPT_PAR_MXR2] = pval[OPT_PAR_MXR2];
                ptmp[OPT_PAR_MXR3] = pval[OPT_PAR_MXR3];
                vmin = fmin;
              }
            }
            ngrd++;
          }
        }
      }
    }
  }
  pval[OPT_PAR_LMD1] = ptmp[OPT_PAR_LMD1];
  pval[OPT_PAR_LMD2] = ptmp[OPT_PAR_LMD2];
  pval[OPT_PAR_LMD3] = ptmp[OPT_PAR_LMD3];
  pval[OPT_PAR_MXR2] = ptmp[OPT_PAR_MXR2];
  pval[OPT_PAR_MXR3] = ptmp[OPT_PAR_MXR3];
  //pmin[OPT_PAR_LMD1] = pval[OPT_PAR_LMD1]-opt_leff_stp[0];
  //pmin[OPT_PAR_LMD2] = pval[OPT_PAR_LMD2]-opt_leff_stp[1];
  //pmin[OPT_PAR_LMD3] = pval[OPT_PAR_LMD3]-opt_leff_stp[2];
  //pmin[OPT_PAR_MXR2] = pval[OPT_PAR_MXR2]-opt_mixr_stp[1];
  //pmin[OPT_PAR_MXR3] = pval[OPT_PAR_MXR3]-opt_mixr_stp[2];
  //pmax[OPT_PAR_LMD1] = pval[OPT_PAR_LMD1]+opt_leff_stp[0];
  //pmax[OPT_PAR_LMD2] = pval[OPT_PAR_LMD2]+opt_leff_stp[1];
  //pmax[OPT_PAR_LMD3] = pval[OPT_PAR_LMD3]+opt_leff_stp[2];
  //pmax[OPT_PAR_MXR2] = pval[OPT_PAR_MXR2]+opt_mixr_stp[1];
  //pmax[OPT_PAR_MXR3] = pval[OPT_PAR_MXR3]+opt_mixr_stp[2];
  if(cnt_vb > 0)
  {
    // Set parameters
    DSR_FCN_CUSTOM_B0(pval,&fmin);
    DSR_FCN_CUSTOM_B2(pval,&fmin);
    CRS_FCN_CUSTOM(pval,&fmin);
    message(stderr,"grid search results:\n");
    message(stderr,"vtau: %22.14e\n",pval[OPT_PAR_VTAU]);
    message(stderr,"vmie: %22.14e\n",pval[OPT_PAR_VMIE]);
    message(stderr,"wscl: %22.14e\n",pval[OPT_PAR_WSCL]);
    message(stderr,"lef1: %22.14e\n",pval[OPT_PAR_LMD1]);
    message(stderr,"lmd1: %22.14e\n",eff2mod(0,pval[OPT_PAR_LMD1]));
    message(stderr,"lsg1: %22.14e\n",pval[OPT_PAR_LSG1]);
    message(stderr,"mxr2: %22.14e\n",pval[OPT_PAR_MXR2]);
    message(stderr,"lef2: %22.14e\n",pval[OPT_PAR_LMD2]);
    message(stderr,"lmd2: %22.14e\n",eff2mod(1,pval[OPT_PAR_LMD2]));
    message(stderr,"lsg2: %22.14e\n",pval[OPT_PAR_LSG2]);
    message(stderr,"mxr3: %22.14e\n",pval[OPT_PAR_MXR3]);
    message(stderr,"lef3: %22.14e\n",pval[OPT_PAR_LMD3]);
    message(stderr,"lmd3: %22.14e\n",eff2mod(2,pval[OPT_PAR_LMD3]));
    message(stderr,"lsg3: %22.14e\n",pval[OPT_PAR_LSG3]);
    message(stderr,"vmin: %22.14e\n",vmin);
    for(n=0; n<mie_n_comp; n++)
    {
      message(stderr,"%22.14e %22.14e %22.14e\n",mie_wcom_com[n],pow(10.0,mie_lmod_com[n]),mie_lsgm_com[n]);
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      message(stderr,"%7.2f %13.4e %13.4e %13.4e %13.4e\n",mie_wlen[i],mie_aext[i],mie_aext[i]/mie_aext[mie_iref],
                                                           mie_asca[i]/mie_aext[i],mie_asym[i]);
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        message(stderr,"%7.2f %7.2f %13.4e\n",mie_wlen[i],mie_angl[j],mie_phas[k]);
      }
    }
    Printout();
  }

  // update correction factors
  if(CalcFactor(pval,1,0) < 0) return -1;
  sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
  sim_cflg = 1; // Apply DIS=f->t correction factor from here
  sim_dflg = 0; // Set DIS=f from here
  if(opt_test == 0)
  {
    if(Background(pval,2) < 0) return -1;
    if(RingCenter(pval,2) < 0) return -1;
    if(Uniformity(pval,2) < 0) return -1;
  }

  // optimize mode radius and mixing ratio
  if(opt_skip[OPT_PAR_LMD1]!=1 ||
     opt_skip[OPT_PAR_MXR2]!=1 || opt_skip[OPT_PAR_LMD2]!=1 ||
     opt_skip[OPT_PAR_MXR3]!=1 || opt_skip[OPT_PAR_LMD3]!=1)
  {
    message(stderr,"Optimizing mode radius and mixing ratio\n");
    opt_mode = OPT_MODE_CRS_CUSTOM;
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i!=OPT_PAR_LMD1 && i!=OPT_PAR_MXR2 && i!=OPT_PAR_LMD2 && i!=OPT_PAR_MXR3 && i!=OPT_PAR_LMD3)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      } else
      if(opt_skip[i] == 1)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
      MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPOUT(i+1,fitnam[i],fitval[i],fiterr[i],fitmin[i],fitmax[i],fitnum[i]);
      pval[i] = fitval[i];
    }
    if(opt_skip[OPT_PAR_LMD1] != 1)
    {
      //pmin[OPT_PAR_LMD1] = pval[OPT_PAR_LMD1]-opt_leff_stp[0];
      //pmax[OPT_PAR_LMD1] = pval[OPT_PAR_LMD1]+opt_leff_stp[0];
    }
    if(opt_skip[OPT_PAR_MXR2] != 1)
    {
      //pmin[OPT_PAR_MXR2] = pval[OPT_PAR_MXR2]-opt_mixr_stp[1];
      //pmax[OPT_PAR_MXR2] = pval[OPT_PAR_MXR2]+opt_mixr_stp[1];
    }
    if(opt_skip[OPT_PAR_LMD2] != 1)
    {
      //pmin[OPT_PAR_LMD2] = pval[OPT_PAR_LMD2]-opt_leff_stp[1];
      //pmax[OPT_PAR_LMD2] = pval[OPT_PAR_LMD2]+opt_leff_stp[1];
    }
    if(opt_skip[OPT_PAR_MXR3] != 1)
    {
      //pmin[OPT_PAR_MXR3] = pval[OPT_PAR_MXR3]-opt_mixr_stp[2];
      //pmax[OPT_PAR_MXR3] = pval[OPT_PAR_MXR3]+opt_mixr_stp[2];
    }
    if(opt_skip[OPT_PAR_LMD3] != 1)
    {
      //pmin[OPT_PAR_LMD3] = pval[OPT_PAR_LMD3]-opt_leff_stp[2];
      //pmax[OPT_PAR_LMD3] = pval[OPT_PAR_LMD3]+opt_leff_stp[2];
    }
    vmin = fmin;
  }
  // optimize aerosol optical depth
  if(opt_skip[OPT_PAR_VTAU] != 1)
  {
    message(stderr,"Optimizing aerosol optical depth\n");
    opt_mode = OPT_MODE_DSR_CUSTOM_0;
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_VTAU)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
    MNPOUT(OPT_PAR_VTAU+1,fitnam[OPT_PAR_VTAU],fitval[OPT_PAR_VTAU],fiterr[OPT_PAR_VTAU],
                          fitmin[OPT_PAR_VTAU],fitmax[OPT_PAR_VTAU],fitnum[OPT_PAR_VTAU]);
    pval[OPT_PAR_VTAU] = fitval[OPT_PAR_VTAU];
  }
  // optimize visibility
  if(opt_skip[OPT_PAR_VMIE] != 1)
  {
    message(stderr,"Optimizing visibility\n");
    opt_mode = OPT_MODE_DSR_CUSTOM_0;
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_VMIE)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
    MNPOUT(OPT_PAR_VMIE+1,fitnam[OPT_PAR_VMIE],fitval[OPT_PAR_VMIE],fiterr[OPT_PAR_VMIE],
                          fitmin[OPT_PAR_VMIE],fitmax[OPT_PAR_VMIE],fitnum[OPT_PAR_VMIE]);
    pval[OPT_PAR_VMIE] = fitval[OPT_PAR_VMIE];
  }
  // optimize water scaling factor
  if(opt_skip[OPT_PAR_WSCL] != 1)
  {
    message(stderr,"Optimizing water scaling factor\n");
    opt_mode = OPT_MODE_DSR_CUSTOM_2;
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_WSCL)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
    MNPOUT(OPT_PAR_WSCL+1,fitnam[OPT_PAR_WSCL],fitval[OPT_PAR_WSCL],fiterr[OPT_PAR_WSCL],
                          fitmin[OPT_PAR_WSCL],fitmax[OPT_PAR_WSCL],fitnum[OPT_PAR_WSCL]);
    pval[OPT_PAR_WSCL] = fitval[OPT_PAR_WSCL];
  }
  if(cnt_vb > 0)
  {
    // Set parameters
    DSR_FCN_CUSTOM_B0(pval,&fmin);
    DSR_FCN_CUSTOM_B2(pval,&fmin);
    CRS_FCN_CUSTOM(pval,&fmin);
    message(stderr,"1st MINUIT results:\n");
    message(stderr,"vtau: %22.14e\n",pval[OPT_PAR_VTAU]);
    message(stderr,"vmie: %22.14e\n",pval[OPT_PAR_VMIE]);
    message(stderr,"wscl: %22.14e\n",pval[OPT_PAR_WSCL]);
    message(stderr,"lef1: %22.14e\n",pval[OPT_PAR_LMD1]);
    message(stderr,"lmd1: %22.14e\n",eff2mod(0,pval[OPT_PAR_LMD1]));
    message(stderr,"lsg1: %22.14e\n",pval[OPT_PAR_LSG1]);
    message(stderr,"mxr2: %22.14e\n",pval[OPT_PAR_MXR2]);
    message(stderr,"lef2: %22.14e\n",pval[OPT_PAR_LMD2]);
    message(stderr,"lmd2: %22.14e\n",eff2mod(1,pval[OPT_PAR_LMD2]));
    message(stderr,"lsg2: %22.14e\n",pval[OPT_PAR_LSG2]);
    message(stderr,"mxr3: %22.14e\n",pval[OPT_PAR_MXR3]);
    message(stderr,"lef3: %22.14e\n",pval[OPT_PAR_LMD3]);
    message(stderr,"lmd3: %22.14e\n",eff2mod(2,pval[OPT_PAR_LMD3]));
    message(stderr,"lsg3: %22.14e\n",pval[OPT_PAR_LSG3]);
    message(stderr,"vmin: %22.14e\n",vmin);
    for(n=0; n<mie_n_comp; n++)
    {
      message(stderr,"%22.14e %22.14e %22.14e\n",mie_wcom_com[n],pow(10.0,mie_lmod_com[n]),mie_lsgm_com[n]);
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      message(stderr,"%7.2f %13.4e %13.4e %13.4e %13.4e\n",mie_wlen[i],mie_aext[i],mie_aext[i]/mie_aext[mie_iref],
                                                           mie_asca[i]/mie_aext[i],mie_asym[i]);
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        message(stderr,"%7.2f %7.2f %13.4e\n",mie_wlen[i],mie_angl[j],mie_phas[k]);
      }
    }
    Printout();
  }

  // update correction factors
  if(CalcFactor(pval,1,0) < 0) return -1;
  sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
  sim_cflg = 1; // Apply DIS=f->t correction factor from here
  sim_dflg = 0; // Set DIS=f from here
  if(opt_test == 0)
  {
    if(Background(pval,2) < 0) return -1;
    if(RingCenter(pval,2) < 0) return -1;
    if(Uniformity(pval,2) < 0) return -1;
  }

  // optimize mode radius and mixing ratio
  if(opt_skip[OPT_PAR_LMD1]!=1 ||
     opt_skip[OPT_PAR_MXR2]!=1 || opt_skip[OPT_PAR_LMD2]!=1 ||
     opt_skip[OPT_PAR_MXR3]!=1 || opt_skip[OPT_PAR_LMD3]!=1)
  {
    message(stderr,"Optimizing mode radius and mixing ratio\n");
    opt_mode = OPT_MODE_CRS_CUSTOM;
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i!=OPT_PAR_LMD1 && i!=OPT_PAR_MXR2 && i!=OPT_PAR_LMD2 && i!=OPT_PAR_MXR3 && i!=OPT_PAR_LMD3)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      } else
      if(opt_skip[i] == 1)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
      MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPOUT(i+1,fitnam[i],fitval[i],fiterr[i],fitmin[i],fitmax[i],fitnum[i]);
      pval[i] = fitval[i];
    }
    if(opt_skip[OPT_PAR_LMD1] != 1)
    {
      //pmin[OPT_PAR_LMD1] = pval[OPT_PAR_LMD1]-opt_leff_stp[0];
      //pmax[OPT_PAR_LMD1] = pval[OPT_PAR_LMD1]+opt_leff_stp[0];
    }
    if(opt_skip[OPT_PAR_MXR2] != 1)
    {
      //pmin[OPT_PAR_MXR2] = pval[OPT_PAR_MXR2]-opt_mixr_stp[1];
      //pmax[OPT_PAR_MXR2] = pval[OPT_PAR_MXR2]+opt_mixr_stp[1];
    }
    if(opt_skip[OPT_PAR_LMD2] != 1)
    {
      //pmin[OPT_PAR_LMD2] = pval[OPT_PAR_LMD2]-opt_leff_stp[1];
      //pmax[OPT_PAR_LMD2] = pval[OPT_PAR_LMD2]+opt_leff_stp[1];
    }
    if(opt_skip[OPT_PAR_MXR3] != 1)
    {
      //pmin[OPT_PAR_MXR3] = pval[OPT_PAR_MXR3]-opt_mixr_stp[2];
      //pmax[OPT_PAR_MXR3] = pval[OPT_PAR_MXR3]+opt_mixr_stp[2];
    }
    if(opt_skip[OPT_PAR_LMD3] != 1)
    {
      //pmin[OPT_PAR_LMD3] = pval[OPT_PAR_LMD3]-opt_leff_stp[2];
      //pmax[OPT_PAR_LMD3] = pval[OPT_PAR_LMD3]+opt_leff_stp[2];
    }
    vmin = fmin;
  }
  // optimize aerosol optical depth
  if(opt_skip[OPT_PAR_VTAU] != 1)
  {
    message(stderr,"Optimizing aerosol optical depth\n");
    opt_mode = OPT_MODE_DSR_CUSTOM_0;
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_VTAU)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
      MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    }
    MNPOUT(OPT_PAR_VTAU+1,fitnam[OPT_PAR_VTAU],fitval[OPT_PAR_VTAU],fiterr[OPT_PAR_VTAU],
                          fitmin[OPT_PAR_VTAU],fitmax[OPT_PAR_VTAU],fitnum[OPT_PAR_VTAU]);
    pval[OPT_PAR_VTAU] = fitval[OPT_PAR_VTAU];
    opt_pval[OPT_PAR_VTAU] = fitval[OPT_PAR_VTAU];
    opt_perr[OPT_PAR_VTAU] = fiterr[OPT_PAR_VTAU];
    opt_fmin[OPT_PAR_VTAU] = fmin;
    opt_fedm[OPT_PAR_VTAU] = fedm;
    opt_stat[OPT_PAR_VTAU] = istat;
  }
  // optimize visibility
  if(opt_skip[OPT_PAR_VMIE] != 1)
  {
    message(stderr,"Optimizing visibility\n");
    opt_mode = OPT_MODE_DSR_CUSTOM_0;
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_VMIE)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
      MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    }
    MNPOUT(OPT_PAR_VMIE+1,fitnam[OPT_PAR_VMIE],fitval[OPT_PAR_VMIE],fiterr[OPT_PAR_VMIE],
                          fitmin[OPT_PAR_VMIE],fitmax[OPT_PAR_VMIE],fitnum[OPT_PAR_VMIE]);
    pval[OPT_PAR_VMIE] = fitval[OPT_PAR_VMIE];
    opt_pval[OPT_PAR_VMIE] = fitval[OPT_PAR_VMIE];
    opt_perr[OPT_PAR_VMIE] = fiterr[OPT_PAR_VMIE];
    opt_fmin[OPT_PAR_VMIE] = fmin;
    opt_fedm[OPT_PAR_VMIE] = fedm;
    opt_stat[OPT_PAR_VMIE] = istat;
  }
  // optimize water scaling factor
  if(opt_skip[OPT_PAR_WSCL] != 1)
  {
    message(stderr,"Optimizing water scaling factor\n");
    opt_mode = OPT_MODE_DSR_CUSTOM_2;
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_WSCL)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
      MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    }
    MNPOUT(OPT_PAR_WSCL+1,fitnam[OPT_PAR_WSCL],fitval[OPT_PAR_WSCL],fiterr[OPT_PAR_WSCL],
                          fitmin[OPT_PAR_WSCL],fitmax[OPT_PAR_WSCL],fitnum[OPT_PAR_WSCL]);
    pval[OPT_PAR_WSCL] = fitval[OPT_PAR_WSCL];
    opt_pval[OPT_PAR_WSCL] = fitval[OPT_PAR_WSCL];
    opt_perr[OPT_PAR_WSCL] = fiterr[OPT_PAR_WSCL];
    opt_fmin[OPT_PAR_WSCL] = fmin;
    opt_fedm[OPT_PAR_WSCL] = fedm;
    opt_stat[OPT_PAR_WSCL] = istat;
  }
  if(cnt_vb > 0)
  {
    // Set parameters
    DSR_FCN_CUSTOM_B0(pval,&fmin);
    DSR_FCN_CUSTOM_B2(pval,&fmin);
    CRS_FCN_CUSTOM(pval,&fmin);
    message(stderr,"2nd MINUIT results:\n");
    message(stderr,"vtau: %22.14e\n",pval[OPT_PAR_VTAU]);
    message(stderr,"vmie: %22.14e\n",pval[OPT_PAR_VMIE]);
    message(stderr,"wscl: %22.14e\n",pval[OPT_PAR_WSCL]);
    message(stderr,"lef1: %22.14e\n",pval[OPT_PAR_LMD1]);
    message(stderr,"lmd1: %22.14e\n",eff2mod(0,pval[OPT_PAR_LMD1]));
    message(stderr,"lsg1: %22.14e\n",pval[OPT_PAR_LSG1]);
    message(stderr,"mxr2: %22.14e\n",pval[OPT_PAR_MXR2]);
    message(stderr,"lef2: %22.14e\n",pval[OPT_PAR_LMD2]);
    message(stderr,"lmd2: %22.14e\n",eff2mod(1,pval[OPT_PAR_LMD2]));
    message(stderr,"lsg2: %22.14e\n",pval[OPT_PAR_LSG2]);
    message(stderr,"mxr3: %22.14e\n",pval[OPT_PAR_MXR3]);
    message(stderr,"lef3: %22.14e\n",pval[OPT_PAR_LMD3]);
    message(stderr,"lmd3: %22.14e\n",eff2mod(2,pval[OPT_PAR_LMD3]));
    message(stderr,"lsg3: %22.14e\n",pval[OPT_PAR_LSG3]);
    message(stderr,"vmin: %22.14e\n",vmin);
    for(n=0; n<mie_n_comp; n++)
    {
      message(stderr,"%22.14e %22.14e %22.14e\n",mie_wcom_com[n],pow(10.0,mie_lmod_com[n]),mie_lsgm_com[n]);
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      message(stderr,"%7.2f %13.4e %13.4e %13.4e %13.4e\n",mie_wlen[i],mie_aext[i],mie_aext[i]/mie_aext[mie_iref],
                                                           mie_asca[i]/mie_aext[i],mie_asym[i]);
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        message(stderr,"%7.2f %7.2f %13.4e\n",mie_wlen[i],mie_angl[j],mie_phas[k]);
      }
    }
    Printout();
  }

  // update correction factors
  if(CalcFactor(pval,1,0) < 0) return -1;
  sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
  sim_cflg = 1; // Apply DIS=f->t correction factor from here
  sim_dflg = 0; // Set DIS=f from here
  if(opt_test == 0)
  {
    if(Background(pval,2) < 0) return -1;
    if(RingCenter(pval,2) < 0) return -1;
    if(Uniformity(pval,2) < 0) return -1;
  }

  // optimize mode radius and mixing ratio
  if(opt_skip[OPT_PAR_LMD1]!=1 ||
     opt_skip[OPT_PAR_MXR2]!=1 || opt_skip[OPT_PAR_LMD2]!=1 ||
     opt_skip[OPT_PAR_MXR3]!=1 || opt_skip[OPT_PAR_LMD3]!=1)
  {
    message(stderr,"Optimizing mode radius and mixing ratio\n");
    opt_mode = OPT_MODE_CRS_CUSTOM;
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i!=OPT_PAR_LMD1 && i!=OPT_PAR_MXR2 && i!=OPT_PAR_LMD2 && i!=OPT_PAR_MXR3 && i!=OPT_PAR_LMD3)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      } else
      if(opt_skip[i] == 1)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
      MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    }
    for(i=OPT_PAR_NCMP; i<OPT_NPAR; i++)
    {
      MNPOUT(i+1,fitnam[i],fitval[i],fiterr[i],fitmin[i],fitmax[i],fitnum[i]);
      pval[i] = fitval[i];
      opt_pval[i] = fitval[i];
      opt_perr[i] = fiterr[i];
      if(fitnum[i] > 0)
      {
        opt_fmin[i] = fmin;
        opt_fedm[i] = fedm;
        opt_stat[i] = istat;
      }
    }
    vmin = fmin;
  }
  // Set parameters
  DSR_FCN_CUSTOM_B0(pval,&fmin);
  DSR_FCN_CUSTOM_B2(pval,&fmin);
  CRS_FCN_CUSTOM(pval,&fmin);
  if(cnt_vb > 0)
  {
    message(stderr,"3rd MINUIT results:\n");
    message(stderr,"vtau: %22.14e\n",pval[OPT_PAR_VTAU]);
    message(stderr,"vmie: %22.14e\n",pval[OPT_PAR_VMIE]);
    message(stderr,"wscl: %22.14e\n",pval[OPT_PAR_WSCL]);
    message(stderr,"lef1: %22.14e\n",pval[OPT_PAR_LMD1]);
    message(stderr,"lmd1: %22.14e\n",eff2mod(0,pval[OPT_PAR_LMD1]));
    message(stderr,"lsg1: %22.14e\n",pval[OPT_PAR_LSG1]);
    message(stderr,"mxr2: %22.14e\n",pval[OPT_PAR_MXR2]);
    message(stderr,"lef2: %22.14e\n",pval[OPT_PAR_LMD2]);
    message(stderr,"lmd2: %22.14e\n",eff2mod(1,pval[OPT_PAR_LMD2]));
    message(stderr,"lsg2: %22.14e\n",pval[OPT_PAR_LSG2]);
    message(stderr,"mxr3: %22.14e\n",pval[OPT_PAR_MXR3]);
    message(stderr,"lef3: %22.14e\n",pval[OPT_PAR_LMD3]);
    message(stderr,"lmd3: %22.14e\n",eff2mod(2,pval[OPT_PAR_LMD3]));
    message(stderr,"lsg3: %22.14e\n",pval[OPT_PAR_LSG3]);
    message(stderr,"vmin: %22.14e\n",vmin);
    for(n=0; n<mie_n_comp; n++)
    {
      message(stderr,"%22.14e %22.14e %22.14e\n",mie_wcom_com[n],pow(10.0,mie_lmod_com[n]),mie_lsgm_com[n]);
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      message(stderr,"%7.2f %13.4e %13.4e %13.4e %13.4e\n",mie_wlen[i],mie_aext[i],mie_aext[i]/mie_aext[mie_iref],
                                                           mie_asca[i]/mie_aext[i],mie_asym[i]);
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        message(stderr,"%7.2f %7.2f %13.4e\n",mie_wlen[i],mie_angl[j],mie_phas[k]);
      }
    }
  }

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
    printf("%2d %18.10e %18.10e %18.10e %18.10e %d\n",
            i,opt_pval[i],opt_perr[i],opt_fmin[i],opt_fedm[i],opt_stat[i]);
  }
  fflush(stdout);
  fflush(stderr);

  return 0;
}

int OwnInit(void)
{
  int n;

  for(n=0; n<OPT_NPAR; n++)
  {
    opt_skip[n] = 0;
    opt_stat[n] = -1;
    opt_pval[n] = 0.0;
    opt_perr[n] = 0.0;
    opt_fmin[n] = 0.0;
    opt_fedm[n] = 0.0;
  }

  return 0;
}

int ReadOwnConfig(char *s)
{
  int n,nc;
  int idx;
  int err;
  char temp[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char fnam[] = "ReadOwnConfig";
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
    if(strcasecmp(str[0],"opt_skip") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=OPT_NPAR)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,s);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        opt_skip[idx] = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27d\n",str[0],idx,opt_skip[idx]);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30d\n",str[0],opt_amod_skp);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_comp") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=OPT_N_COMP)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,s);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        opt_comp[idx] = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27d\n",str[0],idx,opt_comp[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_lmod_ngp") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=OPT_N_COMP)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,s);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        opt_lmod_ngp[idx] = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27d\n",str[0],idx,opt_lmod_ngp[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_mixr_ngp") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(idx<1 || idx>=OPT_N_COMP)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,s);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        opt_mixr_ngp[idx] = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27d\n",str[0],idx,opt_mixr_ngp[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_lmod_min") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=OPT_N_COMP)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,s);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        opt_lmod_min[idx] = strtod(str[2],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27.4e\n",str[0],idx,opt_lmod_min[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_lmod_max") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=OPT_N_COMP)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,s);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        opt_lmod_max[idx] = strtod(str[2],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27.4e\n",str[0],idx,opt_lmod_max[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_mixr_min") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(idx<1 || idx>=OPT_N_COMP)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,s);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        opt_mixr_min[idx] = strtod(str[2],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27.4e\n",str[0],idx,opt_mixr_min[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_mixr_max") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(idx<1 || idx>=OPT_N_COMP)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,s);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        opt_mixr_max[idx] = strtod(str[2],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27.4e\n",str[0],idx,opt_mixr_max[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_dsr_ernd") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_dsr_ernd = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n",str[0],tst_dsr_ernd);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_dsr_esys") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_dsr_esys = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n",str[0],tst_dsr_esys);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_aur_ernd") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_aur_ernd = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n",str[0],tst_aur_ernd);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_aur_esys") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_aur_esys = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n",str[0],tst_aur_esys);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_ssr_ernd") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_ssr_ernd = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n",str[0],tst_ssr_ernd);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_ssr_esys") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_ssr_esys = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n",str[0],tst_ssr_esys);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_seed") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_seed = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30ld\n",str[0],tst_seed);
        cnt_n_cmnt++;
      }
    }
    else
    {
      return 1;
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
  int option_index = 0;
  int this_option_optind;
  struct option long_options[] =
  {
    {"test",0,0,'T'},
  };

  fprintf(stderr,"%15s : $Revision: 1127 $ $Date: 2016-07-14 14:15:47 +0900 (Thu, 14 Jul 2016) $\n","aeros_mix.c");

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-T",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'T':
        opt_test = 1;
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
  char optstring[] = "DNteFUVWOsAaGgwmMrxXnSCcoTdvh";

  cnt_n_own_format = 16;
  sprintf(cnt_own_format[ 0],"opt_skip      # flag       | parameter#(%d),skip optimization(%d)\n",0,0);
  sprintf(cnt_own_format[ 1],"opt_amod_skp  flag         | skip amod optimization(%d)\n",0);
  sprintf(cnt_own_format[ 2],"opt_comp      # #          | component#(%d),component ID(%d)\n",0,8);
  sprintf(cnt_own_format[ 3],"opt_lmod_ngp  # #          | component#(%d),lmod #grid(%d)\n",0,8);
  sprintf(cnt_own_format[ 4],"opt_mixr_ngp  # #          | component#(%d),mixr #grid(%d)\n",1,8);
  sprintf(cnt_own_format[ 5],"opt_lmod_min  # value      | component#(%d),Min lmod(%4.1f)\n",0,-1.4);
  sprintf(cnt_own_format[ 6],"opt_lmod_max  # value      | component#(%d),Max lmod(%4.1f)\n",0,1.3);
  sprintf(cnt_own_format[ 7],"opt_mixr_min  # value      | component#(%d),Min mixr(%4.1f)\n",1,-2.0);
  sprintf(cnt_own_format[ 8],"opt_mixr_max  # value      | component#(%d),Max mixr(%4.1f)\n",1,2.0);
  sprintf(cnt_own_format[ 9],"tst_dsr_ernd  value        | DSR random     error(%.1f)\n",0.0);
  sprintf(cnt_own_format[10],"tst_dsr_esys  value        | DSR systematic error(%.1f)\n",0.0);
  sprintf(cnt_own_format[11],"tst_aur_ernd  value        | AUR random     error(%.1f)\n",0.0);
  sprintf(cnt_own_format[12],"tst_aur_esys  value        | AUR systematic error(%.1f)\n",0.0);
  sprintf(cnt_own_format[13],"tst_ssr_ernd  value        | SSR random     error(%.1f)\n",0.0);
  sprintf(cnt_own_format[14],"tst_ssr_esys  value        | SSR systematic error(%.1f)\n",0.0);
  sprintf(cnt_own_format[15],"tst_seed      value        | Random seed(%d)\n",0);

  fprintf(stderr,"aeros_mix ... fit DSR/SSR spectrum.\n");
  CommonUsage("aeros_mix",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'T':
        fprintf(stderr," T -opt_test      |%s|%s|%s| %d\n",As(e,"Test    mode",n),   As(a,"nothing",n),    Ad(d,0,n),opt_test);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);

  return 0;
}
