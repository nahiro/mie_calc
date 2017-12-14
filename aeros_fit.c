/*********************************************************/
/* aeros_fit                                             */
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
// Other constants
//#define	f2cFortran

// Parameters for Optimization
int	opt_test			= 0;			// Test mode
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

int Optimize(void);
void PrintTime(void);

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

  cnt_dsr[0] = 1;
  cnt_dsr[1] = 1;
  cnt_dsr[2] = 1;
  cnt_dsr[3] = 1;
  cnt_aur[0] = 0;
  cnt_aur[1] = 0;
  cnt_ssr[0] = 0;
  cnt_ssr[1] = 0;
  cnt_ssr[2] = 0;
  cnt_ssr[3] = 0;
  if(ReadInput() < 0) return -1;
  if(dsr_n_data < 1)
  {
    fprintf(stderr,"Error, dsr_n_data=%d\n",dsr_n_data);
    return -1;
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
    pval[OPT_PAR_VMIE] = opt_init[OPT_PAR_VMIE];
    pval[OPT_PAR_WSCL] = opt_init[OPT_PAR_WSCL];
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
  int i,n;
  int flag;
  int dummy;
  int npari,nparx,istat;
  double vmin;
  double fmin,fedm,errdef;
  double pval[SIM_NPAR];
  double pstp[OPT_NPAR];
  double pmin[OPT_NPAR];
  double pmax[OPT_NPAR];
  char   pnam[OPT_NPAR][MAXLINE] = {"vmie","wscl",
                                    "ncmp","lmd1","lsg1",
                                    "mxr2","lmd2","lsg2",
                                    "mxr3","lmd3","lsg3"};
  int    fitnum[OPT_NPAR];
  double fitval[OPT_NPAR];
  double fiterr[OPT_NPAR];
  double fitmin[OPT_NPAR];
  double fitmax[OPT_NPAR];
  char   fitnam[OPT_NPAR][MAXLINE];

  // initial guess
  if(!isnan(opt_vmie_1st))
  {
    pval[OPT_PAR_VMIE] = opt_vmie_1st;
  }
  else
  {
    pval[OPT_PAR_VMIE] = opt_init[OPT_PAR_VMIE];
  }
  if(!isnan(opt_wscl_1st))
  {
    pval[OPT_PAR_WSCL] = opt_wscl_1st;
  }
  else
  {
    pval[OPT_PAR_WSCL] = opt_init[OPT_PAR_WSCL];
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
  pmin[OPT_PAR_VMIE] = OPT_VMIE_MIN;
  pmin[OPT_PAR_WSCL] = OPT_WSCL_MIN;
  pmin[OPT_PAR_NCMP] = OPT_NCMP_MIN;
  pmin[OPT_PAR_MXR2] = OPT_MIXR_MIN;
  pmin[OPT_PAR_MXR3] = OPT_MIXR_MIN;
  pmin[OPT_PAR_LMD1] = OPT_LMOD_MIN;
  pmin[OPT_PAR_LMD2] = OPT_LMOD_MIN;
  pmin[OPT_PAR_LMD3] = OPT_LMOD_MIN;
  pmin[OPT_PAR_LSG1] = OPT_LSGM_MIN;
  pmin[OPT_PAR_LSG2] = OPT_LSGM_MIN;
  pmin[OPT_PAR_LSG3] = OPT_LSGM_MIN;
  pmax[OPT_PAR_VMIE] = OPT_VMIE_MAX;
  pmax[OPT_PAR_WSCL] = OPT_WSCL_MAX;
  pmax[OPT_PAR_NCMP] = OPT_NCMP_MAX;
  pmax[OPT_PAR_MXR2] = OPT_MIXR_MAX;
  pmax[OPT_PAR_MXR3] = OPT_MIXR_MAX;
  pmax[OPT_PAR_LMD1] = OPT_LMOD_MAX;
  pmax[OPT_PAR_LMD2] = OPT_LMOD_MAX;
  pmax[OPT_PAR_LMD3] = OPT_LMOD_MAX;
  pmax[OPT_PAR_LSG1] = OPT_LSGM_MAX;
  pmax[OPT_PAR_LSG2] = OPT_LSGM_MAX;
  pmax[OPT_PAR_LSG3] = OPT_LSGM_MAX;
  pstp[OPT_PAR_VMIE] = OPT_VMIE_STP;
  pstp[OPT_PAR_WSCL] = OPT_WSCL_STP;
  pstp[OPT_PAR_NCMP] = OPT_NCMP_STP;
  pstp[OPT_PAR_MXR2] = OPT_MIXR_STP;
  pstp[OPT_PAR_MXR3] = OPT_MIXR_STP;
  pstp[OPT_PAR_LMD1] = OPT_LMOD_STP;
  pstp[OPT_PAR_LMD2] = OPT_LMOD_STP;
  pstp[OPT_PAR_LMD3] = OPT_LMOD_STP;
  pstp[OPT_PAR_LSG1] = OPT_LSGM_STP;
  pstp[OPT_PAR_LSG2] = OPT_LSGM_STP;
  pstp[OPT_PAR_LSG3] = OPT_LSGM_STP;
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
    flag = 0; // MODTRAN model
  }
  if(cnt_vb) PrintTime();

  // optimize visibility
  if(opt_vmie_skp != 1)
  {
    message(stderr,"Optimizing visibility\n");
    opt_mode = (flag==1?OPT_MODE_DSR_CUSTOM_0:OPT_MODE_DSR_MODTRAN_0);
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
  if(opt_wscl_skp != 1)
  {
    message(stderr,"Optimizing water scaling factor\n");
    opt_mode = (flag==1?OPT_MODE_DSR_CUSTOM_2:OPT_MODE_DSR_MODTRAN_2);
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
    opt_mode = (flag==1?OPT_MODE_DSR_CUSTOM_1:OPT_MODE_DSR_MODTRAN_1);
    vmin = 1.0e100;
    opt_amod = -1;
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
    message(stderr,"Model optimization summary\n");
    for(i=0; i<opt_n_builtin; i++)
    {
      message(stderr,"vmin[%2d]: %18.10e\n",i,opt_vmin[i]);
    }
    message(stderr,"opt_amod                     : %18d\n",opt_amod);
    sim_ihaz = opt_ihaz[opt_amod];
  }
  // optimize visibility
  if(opt_vmie_skp != 1)
  {
    message(stderr,"Optimizing visibility\n");
    opt_mode = (flag==1?OPT_MODE_DSR_CUSTOM_0:OPT_MODE_DSR_MODTRAN_0);
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
  if(opt_wscl_skp != 1)
  {
    message(stderr,"Optimizing water scaling factor\n");
    opt_mode = (flag==1?OPT_MODE_DSR_CUSTOM_2:OPT_MODE_DSR_MODTRAN_2);
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
  message(stderr,"vmie: %18.10e\n",pval[OPT_PAR_VMIE]);
  message(stderr,"wscl: %18.10e\n",pval[OPT_PAR_WSCL]);

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
  if(flag == 1)
  {
    if(Background(pval,2) < 0) return -1;
  }
  else
  {
    if(Background(pval,0) < 0) return -1;
  }

  // optimize visibility
  message(stderr,"Optimizing visibility\n");
  opt_mode = (flag==1?OPT_MODE_DSR_CUSTOM_0:OPT_MODE_DSR_MODTRAN_0);
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
  // optimize water scaling factor
  message(stderr,"Optimizing water scaling factor\n");
  opt_mode = (flag==1?OPT_MODE_DSR_CUSTOM_2:OPT_MODE_DSR_MODTRAN_2);
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

  // optimize visibility
  message(stderr,"Optimizing visibility\n");
  opt_mode = (flag==1?OPT_MODE_DSR_CUSTOM_0:OPT_MODE_DSR_MODTRAN_0);
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
  // optimize water scaling factor
  message(stderr,"Optimizing water scaling factor\n");
  opt_mode = (flag==1?OPT_MODE_DSR_CUSTOM_2:OPT_MODE_DSR_MODTRAN_2);
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
  opt_n_opt = 0;
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
  message(stderr,"wscl: %18.10e\n",pval[OPT_PAR_WSCL]);
  if(opt_test)
  {
    opt_n_opt = 0;
    snprintf(opt_command,MAXLINE,"SCAN %d 100 0.01 2.0",OPT_PAR_WSCL+1);
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    // to plot the procedure
    snprintf(opt_command,MAXLINE,"CLE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i == OPT_PAR_WSCL)
      {
        MNPARM(i+1,pnam[i],SIM_WSCL,pstp[i],pmin[i],pmax[i],dummy);
      }
      else
      {
        MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
      }
    }
    for(i=0; i<OPT_NPAR; i++)
    {
      if(i != OPT_PAR_WSCL)
      {
        snprintf(opt_command,MAXLINE,"FIX %d",i+1);
        MNCOMD(minuitfcn,opt_command,dummy,NULL);
      }
    }
    opt_n_opt = 0;
    snprintf(opt_command,MAXLINE,"MINI");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
    MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
    if(istat == OPT_OK)
    {
      snprintf(opt_command,MAXLINE,"IMPROVE");
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
  }
  sim_wscl = pval[OPT_PAR_WSCL];

  // optimize ozone scaling factor
  message(stderr,"Optimizing ozone scaling factor\n");
  opt_mode = (flag==1?OPT_MODE_DSR_CUSTOM_3:OPT_MODE_DSR_MODTRAN_3);
  pval[OPT_PAR_OSCL] = SIM_OSCL;
  snprintf(opt_command,MAXLINE,"CLE");
  MNCOMD(minuitfcn,opt_command,dummy,NULL);
  for(i=0; i<OPT_NPAR; i++)
  {
    MNPARM(i+1,pnam[i],pval[i],pstp[i],pmin[i],pmax[i],dummy);
  }
  for(i=0; i<OPT_NPAR; i++)
  {
    if(i != OPT_PAR_OSCL)
    {
      snprintf(opt_command,MAXLINE,"FIX %d",i+1);
      MNCOMD(minuitfcn,opt_command,dummy,NULL);
    }
  }
  opt_n_opt = 0;
  snprintf(opt_command,MAXLINE,"MINI");
  MNCOMD(minuitfcn,opt_command,dummy,NULL);
  MNSTAT(fmin,fedm,errdef,npari,nparx,istat);
  if(istat == OPT_OK)
  {
    snprintf(opt_command,MAXLINE,"IMPROVE");
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
  }
  MNPOUT(OPT_PAR_OSCL+1,fitnam[OPT_PAR_OSCL],fitval[OPT_PAR_OSCL],fiterr[OPT_PAR_OSCL],
                        fitmin[OPT_PAR_OSCL],fitmax[OPT_PAR_OSCL],fitnum[OPT_PAR_OSCL]);
  pval[OPT_PAR_OSCL] = fitval[OPT_PAR_OSCL];
  message(stderr,"oscl: %18.10e\n",pval[OPT_PAR_OSCL]);
  if(opt_test)
  {
    opt_n_opt = 0;
    snprintf(opt_command,MAXLINE,"SCAN %d 100 0.01 2.0",OPT_PAR_OSCL+1);
    MNCOMD(minuitfcn,opt_command,dummy,NULL);
  }

  return 0;
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

int OwnInit(void)
{
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
  int option_index = 0;
  int this_option_optind;
  struct option long_options[] =
  {
    {"test",0,0,'T'},
  };

  fprintf(stderr,"%15s : $Revision: 1127 $ $Date: 2016-07-14 14:15:47 +0900 (Thu, 14 Jul 2016) $\n","aeros_fit.c");

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
  char optstring[] = "DNteFVWAaGgwmMrxXnsSCcoTdvh";

  cnt_n_own_format = 6;
  sprintf(cnt_own_format[1],"opt_vmie_1st  vmie         | visibility\n");
  sprintf(cnt_own_format[2],"opt_wscl_1st  wscl         | water scaling factor\n");
  sprintf(cnt_own_format[3],"opt_vmie_skp  flag         | skip vmie optimization(%d)\n",0);
  sprintf(cnt_own_format[4],"opt_wscl_skp  flag         | skip wscl optimization(%d)\n",0);
  sprintf(cnt_own_format[5],"opt_amod_skp  flag         | skip amod optimization(%d)\n",0);

  fprintf(stderr,"aeros_fit ... fit DSR/SSR spectrum.\n");
  CommonUsage("aeros_fit",0);
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
