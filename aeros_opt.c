/*********************************************************/
/* aeros_opt                                             */
/* Author: N.Manago Apr,30,2008                          */
/* $Revision: 1.58 $                                     */
/* $Date: 2008-08-13 15:42:42 $ (UTC)                    */
/*********************************************************/
#include <cfortran.h>
#include <minuitfcn.h>
#include <minuit.h>
#include "aeros_comp.h"
#include "aeros_comp.c"

// Constants for Optimization
#define	OPT_NCOL	        	16
#define	OPT_NINT	        	3
#define	OPT_NDBL	        	13
#define	OPT_VMIE_MIN	        	1.0e+0
#define	OPT_VMIE_MAX	        	1.0e+3
#define	OPT_VMIE_STP	        	1.0e+0
#define	OPT_WSCL_MIN	        	0.0
#define	OPT_WSCL_MAX	        	10.0
#define	OPT_WSCL_STP	        	1.0e-2
#define	OPT_NCMP	        	3.0
#define	OPT_NCMP_MIN	        	0.0
#define	OPT_NCMP_MAX	        	3.0
#define	OPT_NCMP_STP	        	1.0
#define	OPT_MIXR	        	0.0
#define	OPT_MIXR_MIN	        	-5.5
#define	OPT_MIXR_MAX	        	+9.5
#define	OPT_MIXR_STP	        	1.0e-2
#define	OPT_LMOD	        	0.3
#define	OPT_LMOD_MIN	        	-3.0
#define	OPT_LMOD_MAX	        	+2.0
#define	OPT_LMOD_STP	        	1.0e-2
#define	OPT_LSGM	        	0.3
#define	OPT_LSGM_MIN	        	0.0
#define	OPT_LSGM_MAX	        	0.8
#define	OPT_LSGM_STP	        	1.0e-2
#define	OPT_DSR_TOLERANCE       	0.5
#define	OPT_SSR_TOLERANCE       	0.1
#define	OPT_FLIM	        	1.0e10
#define	OPT_N_BUILTIN	        	3
#define	OPT_MIXR_NSTP	        	12
#define	OPT_LMOD_NSTP	        	20
#define	OPT_LSGM_NSTP	        	15
#define	OPT_PREC	        	1.0e-12			// Precision
#define	OPT_OK		        	3			// Minuit status (CONVERGED)
#define	OPT_DONE	        	10000
#define	OPT_MODE_DSR_MODTRAN_0  	0
#define	OPT_MODE_DSR_MODTRAN_1  	1
#define	OPT_MODE_DSR_MODTRAN_2  	2
#define	OPT_MODE_DSR_MODTRAN_3  	3
#define	OPT_MODE_DSR_CUSTOM_0   	4
#define	OPT_MODE_DSR_CUSTOM_1   	5
#define	OPT_MODE_DSR_CUSTOM_2   	6
#define	OPT_MODE_DSR_CUSTOM_3   	7
#define	OPT_MODE_SSR_MODTRAN_S  	100
#define	OPT_MODE_SSR_MODTRAN_B  	101
#define	OPT_MODE_SSR_CUSTOM_S   	102
#define	OPT_MODE_SSR_CUSTOM_B   	103
// Other constants
//#define	f2cFortran

// Parameters for Optimization
int	opt_mode	        	= -1;			// Optimization mode
int	opt_amod	        	= -1;			// The best fit model#
int	opt_vmie_skp	        	= 0;			// Skip vmie optimization
int	opt_wscl_skp	        	= 0;			// Skip wscl optimization
int	opt_amod_skp	        	= 0;			// Skip amod optimization
int	opt_comp[OPT_N_COMP]    	= {-1,-1,-1};		// The best fit component#
int	opt_flag[OPT_N_COMP][CMP_N_COMP];	        	// Optimization flag
double	opt_vdsr[OPT_N_COMP][CMP_N_COMP];	        	// Optimized DSR function value
double	opt_vssr[OPT_N_COMP][CMP_N_COMP];	        	// Optimized SSR function value
double	opt_pval[OPT_N_COMP][CMP_N_COMP][OPT_NPAR];     	// Optimized parameter value
double	opt_vmie_1st	        	= NAN;			// Visibility in km
double	opt_wscl_1st	        	= NAN;			// Water scale
double	opt_dsr_tolerance       	= OPT_DSR_TOLERANCE;	// DSR tolerance
double	opt_ssr_tolerance       	= OPT_SSR_TOLERANCE;	// SSR tolerance
double	opt_flim	        	= OPT_FLIM;		// Maximum allowed function value
double	opt_vmin[OPT_N_BUILTIN];		        	// Minimum function value
double	opt_mixr[10] =
{          // wcom[0]  wcom[1]
  -1.2788, // 0.05000  0.95000
  -0.9542, // 0.10000  0.90000
  -0.4771, // 0.25000  0.75000
   0.0000, // 0.50000  0.50000
   0.4771, // 0.75000  0.25000
   0.9542, // 0.90000  0.10000
   1.2788, // 0.95000  0.05000
   1.5911, // 0.97500  0.02500
   1.9956, // 0.99000  0.01000
   2.9996, // 0.99900  0.00100
};
int	opt_ihaz[OPT_N_BUILTIN] 	= {1,4,5};		// MODTRAN aerosol models
int	opt_n_builtin	        	= OPT_N_BUILTIN;	// #Models
int	opt_mixr_nstp	        	= OPT_MIXR_NSTP;	// Mixing ratio #steps
int	opt_lmod_nstp	        	= OPT_LMOD_NSTP;	// Log10(R) #steps
int	opt_lsgm_nstp	        	= OPT_LSGM_NSTP;	// Sigma Log10(R) #steps
int	opt_dsr_best	        	= -1;			// The best fit component#
int	opt_ssr_best	        	= -1;			// The best fit component#
int	opt_n_candidate	        	= NODATA;		// #Candidates
int	opt_candidate[CMP_N_COMP];		        	// Candidates
char	*opt_command	        	= NULL;			// Minuit command
char	opt_fres[MAXLINE]       	= NONAME;		// Result file

int Optimize(void);
int ReadResult(void);
void PrintTime(void);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(cnt_hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  Optimize();
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

  if(ReadInput() < 0) return -1;
  if(dsr_n_data<1 || ssr_n_data<1)
  {
    fprintf(stderr,"Error, dsr_n_data=%d, ssr_n_data=%d\n",dsr_n_data,ssr_n_data);
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
  if(ReadResult() < 0) return -1;
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
    }
    pval[OPT_PAR_VMIE] = opt_init[OPT_PAR_VMIE];
    pval[OPT_PAR_WSCL] = opt_init[OPT_PAR_WSCL];
    pval[OPT_PAR_NCMP] = (double)tst_n_comp;
    pval[OPT_PAR_LMD1] = tst_lmod_com[0];
    pval[OPT_PAR_LSG1] = tst_lsgm_com[0];
    for(i=OPT_PAR_MXR2,n=1; n<tst_n_comp; n++)
    {
      pval[i++] = log10(tst_wcom_com[0]/tst_wcom_com[n]);
      pval[i++] = tst_lmod_com[n];
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
  int i,j,k,n;
  int dummy;
  double vmin,vmax;
  double pval[OPT_NPAR];
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
  double fmin;
  int candidate[1000];
  double lmin,lmax,lstp;
  double ltmp[1000];
  double mtmp[1000];
  double tmp_vmin[1000];

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
  if(cnt_vb) PrintTime();
  // optimize visibility
  if(opt_vmie_skp != 1)
  {
    message(stderr,"Optimizing visibility\n");
    opt_mode = OPT_MODE_DSR_MODTRAN_0;
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
    MNPOUT(OPT_PAR_VMIE+1,fitnam[OPT_PAR_VMIE],fitval[OPT_PAR_VMIE],fiterr[OPT_PAR_VMIE],
                          fitmin[OPT_PAR_VMIE],fitmax[OPT_PAR_VMIE],fitnum[OPT_PAR_VMIE]);
    pval[OPT_PAR_VMIE] = fitval[OPT_PAR_VMIE];
  }
  // optimize water scaling factor
  if(opt_wscl_skp != 1)
  {
    message(stderr,"Optimizing water scaling factor\n");
    opt_mode = OPT_MODE_DSR_MODTRAN_2;
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
    MNPOUT(OPT_PAR_WSCL+1,fitnam[OPT_PAR_WSCL],fitval[OPT_PAR_WSCL],fiterr[OPT_PAR_WSCL],
                          fitmin[OPT_PAR_WSCL],fitmax[OPT_PAR_WSCL],fitnum[OPT_PAR_WSCL]);
    pval[OPT_PAR_WSCL] = fitval[OPT_PAR_WSCL];
  }
  // optimize aerosol model
  if(opt_amod_skp != 1)
  {
    message(stderr,"Optimizing aerosol model\n");
    opt_mode = OPT_MODE_DSR_MODTRAN_1;
    vmin = 1.0e100;
    opt_amod = -1;
    for(n=0; n<opt_n_builtin; n++)
    {
      sim_ihaz = opt_ihaz[n];
      DSR_FCN_MODTRAN_B1(pval,&opt_vmin[n]);
      message(stderr,"opt_vmin[%2d]: %18.10e\n",n,opt_vmin[n]);
      if(opt_vmin[n] < vmin)
      {
        vmin = opt_vmin[n];
        opt_amod = n;
      }
    }
    message(stderr,"Model optimization summary\n");
    for(n=0; n<opt_n_builtin; n++)
    {
      message(stderr,"opt_vmin[%2d]: %18.10e\n",n,opt_vmin[n]);
    }
    message(stderr,"opt_amod                     : %18d\n",opt_amod);
    sim_ihaz = opt_ihaz[opt_amod];
  }
  // optimize visibility
  if(opt_vmie_skp != 1)
  {
    message(stderr,"Optimizing visibility\n");
    opt_mode = OPT_MODE_DSR_MODTRAN_0;
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
    MNPOUT(OPT_PAR_VMIE+1,fitnam[OPT_PAR_VMIE],fitval[OPT_PAR_VMIE],fiterr[OPT_PAR_VMIE],
                          fitmin[OPT_PAR_VMIE],fitmax[OPT_PAR_VMIE],fitnum[OPT_PAR_VMIE]);
    pval[OPT_PAR_VMIE] = fitval[OPT_PAR_VMIE];
  }
  // optimize water scaling factor
  if(opt_wscl_skp != 1)
  {
    message(stderr,"Optimizing water scaling factor\n");
    opt_mode = OPT_MODE_DSR_MODTRAN_2;
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
    MNPOUT(OPT_PAR_WSCL+1,fitnam[OPT_PAR_WSCL],fitval[OPT_PAR_WSCL],fiterr[OPT_PAR_WSCL],
                          fitmin[OPT_PAR_WSCL],fitmax[OPT_PAR_WSCL],fitnum[OPT_PAR_WSCL]);
    pval[OPT_PAR_WSCL] = fitval[OPT_PAR_WSCL];
  }
  message(stderr,"Initial guess\n");
  message(stderr,"vmie: %18.10e\n",pval[OPT_PAR_VMIE]);
  message(stderr,"wscl: %18.10e\n",pval[OPT_PAR_WSCL]);

  // calculate correction factor
  if(CalcFactor(pval,0,1) < 0) return -1;
  sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
  sim_cflg = 1; // Apply DIS=f->t correction factor from here
  sim_dflg = 0; // Set DIS=f from here

  // 1st component
  pval[OPT_PAR_NCMP] = 1.0;
  message(stderr,"1st component\n");
  for(n=0; n<cmp_n_comp; n++)
  {
    if(opt_flag[0][n] == OPT_DONE)
    {
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_VMIE]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_VMIE]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_WSCL]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_WSCL]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_NCMP]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_NCMP]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD1]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_LMD1]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG1]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_LSG1]);
      message(stderr,"opt_vdsr[%d][%2d]              : %18.10e\n",0,n,opt_vdsr[0][n]);
      message(stderr,"opt_vssr[%d][%2d]              : %18.10e\n",0,n,opt_vssr[0][n]);
      message(stderr,"opt_flag[%d][%2d]              : %18d\n",   0,n,opt_flag[0][n]);
      continue;
    }
    if(cnt_vb) PrintTime();
    Interp1D(cmp_wlen,cmp_refr[n],cmp_n_wlen,mie_wlen_um,mie_refr_com[0],mie_n_wlen,1);
    Interp1D(cmp_wlen,cmp_refi[n],cmp_n_wlen,mie_wlen_um,mie_refi_com[0],mie_n_wlen,0);
    if(MieTable(0,0) < 0) return -1;
    // set parameters
    pval[OPT_PAR_LSG1] = cmp_lsgm[n];
    // optimize mode radius for each candidate
    message(stderr,"Optimizing mode radius (1st, %2d/%2d)\n",n+1,cmp_n_comp);
    message(stderr,"DSR_FCN_CUSTOM_1\n");
    lmin = cmp_lmin[n];
    lmax = cmp_lmax[n];
    lstp = (lmax-lmin)/10;
    for(i=0; i<14; i++)
    {
      ltmp[i] = lmin+lstp*((double)(i-2)+0.5);
    }
    vmin = +1.0e100;
    vmax = -1.0e100;
    opt_dsr_best = 0;
    for(i=0; i<14; i++)
    {
      pval[OPT_PAR_LMD1] = ltmp[i];
      DSR_FCN_CUSTOM_B1(pval,&tmp_vmin[i]);
      if(tmp_vmin[i] < vmin)
      {
        vmin = tmp_vmin[i];
        opt_dsr_best = i;
      }
      if(tmp_vmin[i] > vmax)
      {
        vmax = tmp_vmin[i];
      }
    }
    message(stderr,"opt_dsr_best = %4d, vmin = %13.4e\n",opt_dsr_best,vmin);
    pval[OPT_PAR_LMD1] = ltmp[opt_dsr_best];
    if(opt_dsr_best==0 || opt_dsr_best==13 || vmin>0.1*vmax)
    {
      opt_pval[0][n][OPT_PAR_VMIE] = pval[OPT_PAR_VMIE];
      opt_pval[0][n][OPT_PAR_WSCL] = pval[OPT_PAR_WSCL];
      opt_pval[0][n][OPT_PAR_NCMP] = pval[OPT_PAR_NCMP];
      opt_pval[0][n][OPT_PAR_LMD1] = pval[OPT_PAR_LMD1];
      opt_pval[0][n][OPT_PAR_LSG1] = pval[OPT_PAR_LSG1];
      opt_vdsr[0][n]               = vmin;
      opt_vssr[0][n]               = HUGE;
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_VMIE]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_VMIE]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_WSCL]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_WSCL]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_NCMP]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_NCMP]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD1]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_LMD1]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG1]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_LSG1]);
      message(stderr,"opt_vdsr[%d][%2d]              : %18.10e\n",0,n,opt_vdsr[0][n]);
      message(stderr,"opt_vssr[%d][%2d]              : %18.10e\n",0,n,opt_vssr[0][n]);
      message(stderr,"opt_flag[%d][%2d]              : %18d\n",   0,n,opt_flag[0][n]);
      continue;
    }
    lmin = pval[OPT_PAR_LMD1]-lstp;
    lmax = pval[OPT_PAR_LMD1]+lstp;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LMD1];
    tmp_vmin[10] = vmin;
    opt_dsr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LMD1] = ltmp[i];
      DSR_FCN_CUSTOM_B1(pval,&tmp_vmin[i]);
      if(tmp_vmin[i] < vmin)
      {
        vmin = tmp_vmin[i];
        opt_dsr_best = i;
      }
    }
    message(stderr,"opt_dsr_best = %4d, vmin = %13.4e\n",opt_dsr_best,vmin);
    opt_n_candidate = 0;
    for(i=0; i<11; i++)
    {
      if(tmp_vmin[i]/vmin < 10.0)
      {
        candidate[opt_n_candidate] = i;
        opt_n_candidate++;
      }
    }
    for(j=0; j<opt_n_candidate; j++)
    {
      message(stderr,"candidate[%4d] = %4d\n",j,candidate[j]);
    }
    message(stderr,"SSR_FCN_CUSTOM_B\n");
    vmin = 1.0e100;
    for(j=0; j<opt_n_candidate; j++)
    {
      i = candidate[j];
      pval[OPT_PAR_LMD1] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&tmp_vmin[i]);
      if(tmp_vmin[i] < vmin)
      {
        vmin = tmp_vmin[i];
        opt_ssr_best = i;
      }
    }
    for(j=0; j<opt_n_candidate; j++)
    {
      i = candidate[j];
      message(stderr,"candidate[%4d] %4d %13.4e\n",j,i,tmp_vmin[i]);
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LMD1] = ltmp[opt_ssr_best];
    // optimize sigma for each candidate
    message(stderr,"Optimizing sigma (1st, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LSG1]*0.5;
    lmax = pval[OPT_PAR_LSG1]*1.5;
    lstp = (lmax-lmin)/5;
    ltmp[5] = pval[OPT_PAR_LSG1];
    opt_ssr_best = 5;
    for(i=0; i<5; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<5; i++)
    {
      pval[OPT_PAR_LSG1] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LSG1] = ltmp[opt_ssr_best];
    lmin = pval[OPT_PAR_LSG1]-lstp;
    lmax = pval[OPT_PAR_LSG1]+lstp;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LSG1];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LSG1] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    opt_vssr[0][n] = vmin;
    pval[OPT_PAR_LSG1] = ltmp[opt_ssr_best];
    DSR_FCN_CUSTOM_B1(pval,&opt_vdsr[0][n]);
    // store the results
    opt_pval[0][n][OPT_PAR_VMIE] = pval[OPT_PAR_VMIE];
    opt_pval[0][n][OPT_PAR_WSCL] = pval[OPT_PAR_WSCL];
    opt_pval[0][n][OPT_PAR_NCMP] = pval[OPT_PAR_NCMP];
    opt_pval[0][n][OPT_PAR_LMD1] = pval[OPT_PAR_LMD1];
    opt_pval[0][n][OPT_PAR_LSG1] = pval[OPT_PAR_LSG1];
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_VMIE]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_VMIE]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_WSCL]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_WSCL]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_NCMP]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_NCMP]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD1]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_LMD1]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG1]: %18.10e\n",0,n,opt_pval[0][n][OPT_PAR_LSG1]);
    message(stderr,"opt_vdsr[%d][%2d]              : %18.10e\n",0,n,opt_vdsr[0][n]);
    message(stderr,"opt_vssr[%d][%2d]              : %18.10e\n",0,n,opt_vssr[0][n]);
    message(stderr,"opt_flag[%d][%2d]              : %18d\n",   0,n,opt_flag[0][n]);
  }
  message(stderr,"1st component summary\n");
  for(n=0; n<cmp_n_comp; n++)
  {
    message(stderr,"%d %2d",0,n);
    for(i=0; i<=OPT_PAR_LSG1; i++)
    {
      message(stderr," %11.4e",opt_pval[0][n][i]);
    }
    message(stderr,"%11.4e %11.4e %3d\n",opt_vdsr[0][n],opt_vssr[0][n],opt_flag[0][n]);
  }
  // find the best fit component
  vmin = 1.0e100;
  opt_ssr_best = -1;
  for(n=0; n<cmp_n_comp; n++)
  {
    if(opt_vssr[0][n] < vmin)
    {
      vmin = opt_vssr[0][n];
      opt_ssr_best = n;
    }
  }
  message(stderr,"opt_ssr_best: %d\n",opt_ssr_best);
  opt_n_candidate = 0;
  for(n=0; n<cmp_n_comp; n++)
  {
    if((opt_vssr[0][n]-vmin)/vmin < opt_ssr_tolerance)
    {
      opt_candidate[opt_n_candidate] = n;
      opt_n_candidate++;
    }
  }
  vmin = 1.0e100;
  opt_dsr_best = -1;
  for(i=0; i<opt_n_candidate; i++)
  {
    n = opt_candidate[i];
    message(stderr,"opt_candidate[%2d]=%d, vdsr=%18.10e, vssr=%18.10e\n",i,n,opt_vdsr[0][n],opt_vssr[0][n]);
    if(opt_vdsr[0][n] < vmin)
    {
      vmin = opt_vdsr[0][n];
      opt_dsr_best = n;
    }
  }
  message(stderr,"opt_dsr_best: %d\n",opt_dsr_best);
  if((opt_vdsr[0][opt_ssr_best]-vmin)/vmin < opt_dsr_tolerance)
  {
    opt_comp[0] = opt_ssr_best;
  }
  else
  {
    opt_comp[0] = opt_dsr_best;
  }
  message(stderr,"opt_comp[%d]                  : %18d\n",0,opt_comp[0]);
  if(vmin > opt_flim)
  {
    message(stderr,"Error, failed in finding the first component >>> %18.10e\n",vmin);
    for(n=0; n<cmp_n_comp; n++)
    {
      message(stderr,"%2d %18.10e %18.10e %d\n",n,opt_vdsr[0][n],opt_vssr[0][n],opt_flag[0][n]);
    }
    pval[OPT_PAR_NCMP] = 0.0;
    return -1;
  }
  // set values from the results
  n = opt_comp[0];
  Interp1D(cmp_wlen,cmp_refr[n],cmp_n_wlen,mie_wlen_um,mie_refr_com[0],mie_n_wlen,1);
  Interp1D(cmp_wlen,cmp_refi[n],cmp_n_wlen,mie_wlen_um,mie_refi_com[0],mie_n_wlen,0);
  if(MieTable(0,0) < 0) return -1;
  pval[OPT_PAR_VMIE] = opt_pval[0][n][OPT_PAR_VMIE];
  pval[OPT_PAR_WSCL] = opt_pval[0][n][OPT_PAR_WSCL];
  pval[OPT_PAR_NCMP] = opt_pval[0][n][OPT_PAR_NCMP];
  pval[OPT_PAR_LMD1] = opt_pval[0][n][OPT_PAR_LMD1];
  pval[OPT_PAR_LSG1] = opt_pval[0][n][OPT_PAR_LSG1];
  // optimize visibility
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
  MNPOUT(OPT_PAR_VMIE+1,fitnam[OPT_PAR_VMIE],fitval[OPT_PAR_VMIE],fiterr[OPT_PAR_VMIE],
                        fitmin[OPT_PAR_VMIE],fitmax[OPT_PAR_VMIE],fitnum[OPT_PAR_VMIE]);
  pval[OPT_PAR_VMIE] = fitval[OPT_PAR_VMIE];
  // optimize water scaling factor
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
  MNPOUT(OPT_PAR_WSCL+1,fitnam[OPT_PAR_WSCL],fitval[OPT_PAR_WSCL],fiterr[OPT_PAR_WSCL],
                        fitmin[OPT_PAR_WSCL],fitmax[OPT_PAR_WSCL],fitnum[OPT_PAR_WSCL]);
  pval[OPT_PAR_WSCL] = fitval[OPT_PAR_WSCL];

  // 2nd component
  pval[OPT_PAR_NCMP] = 2.0;
  message(stderr,"2nd component\n");
  for(n=0; n<cmp_n_comp; n++)
  {
    if(opt_flag[1][n] == OPT_DONE)
    {
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_VMIE]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_VMIE]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_WSCL]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_WSCL]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_NCMP]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_NCMP]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD1]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_LMD1]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG1]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_LSG1]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_MXR2]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_MXR2]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD2]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_LMD2]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG2]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_LSG2]);
      message(stderr,"opt_vdsr[%d][%2d]              : %18.10e\n",1,n,opt_vdsr[1][n]);
      message(stderr,"opt_vssr[%d][%2d]              : %18.10e\n",1,n,opt_vssr[1][n]);
      message(stderr,"opt_flag[%d][%2d]              : %18d\n",   1,n,opt_flag[1][n]);
      continue;
    }
    if(cnt_vb) PrintTime();
    Interp1D(cmp_wlen,cmp_refr[n],cmp_n_wlen,mie_wlen_um,mie_refr_com[1],mie_n_wlen,1);
    Interp1D(cmp_wlen,cmp_refi[n],cmp_n_wlen,mie_wlen_um,mie_refi_com[1],mie_n_wlen,0);
    if(MieTable(1,1) < 0) return -1;
    // set parameters
    pval[OPT_PAR_LMD1] = opt_pval[0][opt_comp[0]][OPT_PAR_LMD1];
    pval[OPT_PAR_LSG1] = opt_pval[0][opt_comp[0]][OPT_PAR_LSG1];
    pval[OPT_PAR_LSG2] = cmp_lsgm[n];
    // optimize mode radius and mixing ratio for each candidate
    message(stderr,"Optimizing mode radius and mixing ratio (2nd, %2d/%2d)\n",n+1,cmp_n_comp);
    message(stderr,"DSR_FCN_CUSTOM_1\n");
    lmin = cmp_lmin[n];
    lmax = cmp_lmax[n];
    lstp = (lmax-lmin)/10;
    for(i=j=0; j<10; j++)
    {
      for(k=0; k<10; i++,k++)
      {
        ltmp[i] = lmin+lstp*((double)j+0.5);
        mtmp[i] = opt_mixr[k];
      }
    }
    vmin = +1.0e100;
    opt_dsr_best = 0;
    for(i=0; i<100; i++)
    {
      pval[OPT_PAR_MXR2] = mtmp[i];
      pval[OPT_PAR_LMD2] = ltmp[i];
      DSR_FCN_CUSTOM_B1(pval,&tmp_vmin[i]);
      if(tmp_vmin[i] < vmin)
      {
        vmin = tmp_vmin[i];
        opt_dsr_best = i;
      }
    }
    message(stderr,"opt_dsr_best = %4d, vmin = %13.4e\n",opt_dsr_best,vmin);
    opt_n_candidate = 0;
    for(i=0; i<100; i++)
    {
      if(tmp_vmin[i]/vmin < 2.0) // this factor must be optimized later
      {
        candidate[opt_n_candidate] = i;
        opt_n_candidate++;
      }
    }
    for(j=0; j<opt_n_candidate; j++)
    {
      message(stderr,"candidate[%4d] = %4d\n",j,candidate[j]);
    }
    message(stderr,"SSR_FCN_CUSTOM_B\n");
    vmin = 1.0e100;
    for(j=0; j<opt_n_candidate; j++)
    {
      i = candidate[j];
      pval[OPT_PAR_MXR2] = mtmp[i];
      pval[OPT_PAR_LMD2] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&tmp_vmin[i]);
      if(tmp_vmin[i] < vmin)
      {
        vmin = tmp_vmin[i];
        opt_ssr_best = i;
      }
    }
    for(j=0; j<opt_n_candidate; j++)
    {
      i = candidate[j];
      message(stderr,"candidate[%4d] %4d %13.4e\n",j,i,tmp_vmin[i]);
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_MXR2] = mtmp[opt_ssr_best];
    pval[OPT_PAR_LMD2] = ltmp[opt_ssr_best];
    // optimize mode radius of the 1st component
    message(stderr,"Optimizing mode radius of the 1st component (2nd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LMD1]-0.07;
    lmax = pval[OPT_PAR_LMD1]+0.07;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LMD1];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LMD1] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LMD1] = ltmp[opt_ssr_best];
    // optimize sigma of the 1st component
    message(stderr,"Optimizing sigma of the 1st component (2nd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LSG1]*0.85;
    lmax = pval[OPT_PAR_LSG1]*1.15;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LSG1];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LSG1] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LSG1] = ltmp[opt_ssr_best];
    // optimize mixing ratio of the 2nd component
    message(stderr,"Optimizing mixing ratio of the 2nd component (2nd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_MXR2]-0.5;
    lmax = pval[OPT_PAR_MXR2]+0.5;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_MXR2];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_MXR2] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_MXR2] = ltmp[opt_ssr_best];
    // optimize mode radius of the 2nd component
    message(stderr,"Optimizing mode radius of the 2nd component (2nd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LMD2]-0.07;
    lmax = pval[OPT_PAR_LMD2]+0.07;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LMD2];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LMD2] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LMD2] = ltmp[opt_ssr_best];
    // optimize sigma for each candidate
    message(stderr,"Optimizing sigma (2nd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LSG2]*0.5;
    lmax = pval[OPT_PAR_LSG2]*1.5;
    lstp = (lmax-lmin)/5;
    ltmp[5] = pval[OPT_PAR_LSG2];
    opt_ssr_best = 5;
    for(i=0; i<5; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<5; i++)
    {
      pval[OPT_PAR_LSG2] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LSG2] = ltmp[opt_ssr_best];
    lmin = pval[OPT_PAR_LSG2]-lstp;
    lmax = pval[OPT_PAR_LSG2]+lstp;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LSG2];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LSG2] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    opt_vssr[1][n] = vmin;
    pval[OPT_PAR_LSG2] = ltmp[opt_ssr_best];
    DSR_FCN_CUSTOM_B1(pval,&opt_vdsr[1][n]);
    // store the results
    opt_pval[1][n][OPT_PAR_VMIE] = pval[OPT_PAR_VMIE];
    opt_pval[1][n][OPT_PAR_WSCL] = pval[OPT_PAR_WSCL];
    opt_pval[1][n][OPT_PAR_NCMP] = pval[OPT_PAR_NCMP];
    opt_pval[1][n][OPT_PAR_LMD1] = pval[OPT_PAR_LMD1];
    opt_pval[1][n][OPT_PAR_LSG1] = pval[OPT_PAR_LSG1];
    opt_pval[1][n][OPT_PAR_MXR2] = pval[OPT_PAR_MXR2];
    opt_pval[1][n][OPT_PAR_LMD2] = pval[OPT_PAR_LMD2];
    opt_pval[1][n][OPT_PAR_LSG2] = pval[OPT_PAR_LSG2];
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_VMIE]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_VMIE]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_WSCL]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_WSCL]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_NCMP]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_NCMP]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD1]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_LMD1]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG1]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_LSG1]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_MXR2]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_MXR2]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD2]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_LMD2]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG2]: %18.10e\n",1,n,opt_pval[1][n][OPT_PAR_LSG2]);
    message(stderr,"opt_vdsr[%d][%2d]              : %18.10e\n",1,n,opt_vdsr[1][n]);
    message(stderr,"opt_vssr[%d][%2d]              : %18.10e\n",1,n,opt_vssr[1][n]);
    message(stderr,"opt_flag[%d][%2d]              : %18d\n",   1,n,opt_flag[1][n]);
  }
  message(stderr,"2nd component summary\n");
  for(n=0; n<cmp_n_comp; n++)
  {
    message(stderr,"%d %2d",1,n);
    for(i=0; i<=OPT_PAR_LSG2; i++)
    {
      message(stderr," %11.4e",opt_pval[1][n][i]);
    }
    message(stderr,"%11.4e %11.4e %3d\n",opt_vdsr[1][n],opt_vssr[1][n],opt_flag[1][n]);
  }
  // find the best fit component
  vmin = 1.0e100;
  opt_ssr_best = -1;
  for(n=0; n<cmp_n_comp; n++)
  {
    if(opt_vssr[1][n] < vmin)
    {
      vmin = opt_vssr[1][n];
      opt_ssr_best = n;
    }
  }
  message(stderr,"opt_ssr_best: %d\n",opt_ssr_best);
  opt_n_candidate = 0;
  for(n=0; n<cmp_n_comp; n++)
  {
    if((opt_vssr[1][n]-vmin)/vmin < opt_ssr_tolerance)
    {
      opt_candidate[opt_n_candidate] = n;
      opt_n_candidate++;
    }
  }
  vmin = 1.0e100;
  opt_dsr_best = -1;
  for(i=0; i<opt_n_candidate; i++)
  {
    n = opt_candidate[i];
    message(stderr,"opt_candidate[%2d]=%d, vdsr=%18.10e, vssr=%18.10e\n",i,n,opt_vdsr[1][n],opt_vssr[1][n]);
    if(opt_vdsr[1][n] < vmin)
    {
      vmin = opt_vdsr[1][n];
      opt_dsr_best = n;
    }
  }
  message(stderr,"opt_dsr_best: %d\n",opt_dsr_best);
  if((opt_vdsr[1][opt_ssr_best]-vmin)/vmin < opt_dsr_tolerance)
  {
    opt_comp[1] = opt_ssr_best;
  }
  else
  {
    opt_comp[1] = opt_dsr_best;
  }
  message(stderr,"opt_comp[%d]                  : %18d\n",1,opt_comp[1]);
  if(vmin > opt_flim)
  {
    message(stderr,"Warning, failed in finding the second component >>> %18.10e\n",vmin);
    for(n=0; n<cmp_n_comp; n++)
    {
      message(stderr,"%2d %18.10e %18.10e %d\n",n,opt_vdsr[1][n],opt_vssr[1][n],opt_flag[1][n]);
    }
    pval[OPT_PAR_NCMP] = 1.0;
    DSR_FCN_CUSTOM_B0(pval,&fmin);
    DSR_FCN_CUSTOM_B1(pval,&fmin);
    DSR_FCN_CUSTOM_B2(pval,&fmin);
    SSR_FCN_CUSTOM_B1(pval,&fmin);
    return 0;
  }
  // set values from the results
  n = opt_comp[1];
  Interp1D(cmp_wlen,cmp_refr[n],cmp_n_wlen,mie_wlen_um,mie_refr_com[1],mie_n_wlen,1);
  Interp1D(cmp_wlen,cmp_refi[n],cmp_n_wlen,mie_wlen_um,mie_refi_com[1],mie_n_wlen,0);
  if(MieTable(1,1) < 0) return -1;
  pval[OPT_PAR_VMIE] = opt_pval[1][n][OPT_PAR_VMIE];
  pval[OPT_PAR_WSCL] = opt_pval[1][n][OPT_PAR_WSCL];
  pval[OPT_PAR_NCMP] = opt_pval[1][n][OPT_PAR_NCMP];
  pval[OPT_PAR_LMD1] = opt_pval[1][n][OPT_PAR_LMD1];
  pval[OPT_PAR_LSG1] = opt_pval[1][n][OPT_PAR_LSG1];
  pval[OPT_PAR_MXR2] = opt_pval[1][n][OPT_PAR_MXR2];
  pval[OPT_PAR_LMD2] = opt_pval[1][n][OPT_PAR_LMD2];
  pval[OPT_PAR_LSG2] = opt_pval[1][n][OPT_PAR_LSG2];
  // optimize visibility
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
  MNPOUT(OPT_PAR_VMIE+1,fitnam[OPT_PAR_VMIE],fitval[OPT_PAR_VMIE],fiterr[OPT_PAR_VMIE],
                        fitmin[OPT_PAR_VMIE],fitmax[OPT_PAR_VMIE],fitnum[OPT_PAR_VMIE]);
  pval[OPT_PAR_VMIE] = fitval[OPT_PAR_VMIE];
  // optimize water scaling factor
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
  MNPOUT(OPT_PAR_WSCL+1,fitnam[OPT_PAR_WSCL],fitval[OPT_PAR_WSCL],fiterr[OPT_PAR_WSCL],
                        fitmin[OPT_PAR_WSCL],fitmax[OPT_PAR_WSCL],fitnum[OPT_PAR_WSCL]);
  pval[OPT_PAR_WSCL] = fitval[OPT_PAR_WSCL];

  // 3rd component
  pval[OPT_PAR_NCMP] = 3.0;
  message(stderr,"3rd component\n");
  for(n=0; n<cmp_n_comp; n++)
  {
    if(opt_flag[2][n] == OPT_DONE)
    {
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_VMIE]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_VMIE]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_WSCL]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_WSCL]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_NCMP]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_NCMP]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD1]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LMD1]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG1]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LSG1]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_MXR2]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_MXR2]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD2]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LMD2]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG2]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LSG2]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_MXR3]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_MXR3]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD3]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LMD3]);
      message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG3]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LSG3]);
      message(stderr,"opt_vdsr[%d][%2d]              : %18.10e\n",2,n,opt_vdsr[2][n]);
      message(stderr,"opt_vssr[%d][%2d]              : %18.10e\n",2,n,opt_vssr[2][n]);
      message(stderr,"opt_flag[%d][%2d]              : %18d\n",   2,n,opt_flag[2][n]);
      continue;
    }
    if(cnt_vb) PrintTime();
    Interp1D(cmp_wlen,cmp_refr[n],cmp_n_wlen,mie_wlen_um,mie_refr_com[2],mie_n_wlen,1);
    Interp1D(cmp_wlen,cmp_refi[n],cmp_n_wlen,mie_wlen_um,mie_refi_com[2],mie_n_wlen,0);
    if(MieTable(2,2) < 0) return -1;
    // set parameters
    pval[OPT_PAR_LMD1] = opt_pval[1][opt_comp[1]][OPT_PAR_LMD1];
    pval[OPT_PAR_LSG1] = opt_pval[1][opt_comp[1]][OPT_PAR_LSG1];
    pval[OPT_PAR_MXR2] = opt_pval[1][opt_comp[1]][OPT_PAR_MXR2];
    pval[OPT_PAR_LMD2] = opt_pval[1][opt_comp[1]][OPT_PAR_LMD2];
    pval[OPT_PAR_LSG2] = opt_pval[1][opt_comp[1]][OPT_PAR_LSG2];
    pval[OPT_PAR_LSG3] = cmp_lsgm[n];
    // optimize mode radius and mixing ratio for each candidate
    message(stderr,"Optimizing mode radius and mixing ratio (3rd, %2d/%2d)\n",n+1,cmp_n_comp);
    message(stderr,"DSR_FCN_CUSTOM_1\n");
    lmin = cmp_lmin[n];
    lmax = cmp_lmax[n];
    lstp = (lmax-lmin)/10;
    for(i=j=0; j<10; j++)
    {
      for(k=0; k<10; i++,k++)
      {
        ltmp[i] = lmin+lstp*((double)j+0.5);
        mtmp[i] = opt_mixr[k];
      }
    }
    vmin = +1.0e100;
    opt_dsr_best = 0;
    for(i=0; i<100; i++)
    {
      pval[OPT_PAR_MXR3] = mtmp[i];
      pval[OPT_PAR_LMD3] = ltmp[i];
      DSR_FCN_CUSTOM_B1(pval,&tmp_vmin[i]);
      if(tmp_vmin[i] < vmin)
      {
        vmin = tmp_vmin[i];
        opt_dsr_best = i;
      }
    }
    message(stderr,"opt_dsr_best = %4d, vmin = %13.4e\n",opt_dsr_best,vmin);
    opt_n_candidate = 0;
    for(i=0; i<100; i++)
    {
      if(tmp_vmin[i]/vmin < 2.0) // this factor must be optimized later
      {
        candidate[opt_n_candidate] = i;
        opt_n_candidate++;
      }
    }
    for(j=0; j<opt_n_candidate; j++)
    {
      message(stderr,"candidate[%4d] = %4d\n",j,candidate[j]);
    }
    message(stderr,"SSR_FCN_CUSTOM_B\n");
    vmin = 1.0e100;
    for(j=0; j<opt_n_candidate; j++)
    {
      i = candidate[j];
      pval[OPT_PAR_MXR3] = mtmp[i];
      pval[OPT_PAR_LMD3] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&tmp_vmin[i]);
      if(tmp_vmin[i] < vmin)
      {
        vmin = tmp_vmin[i];
        opt_ssr_best = i;
      }
    }
    for(j=0; j<opt_n_candidate; j++)
    {
      i = candidate[j];
      message(stderr,"candidate[%4d] %4d %13.4e\n",j,i,tmp_vmin[i]);
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_MXR3] = mtmp[opt_ssr_best];
    pval[OPT_PAR_LMD3] = ltmp[opt_ssr_best];
    // optimize mixing ratio of the 2nd component
    message(stderr,"Optimizing mixing ratio of the 2nd component (3rd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_MXR2]-0.5;
    lmax = pval[OPT_PAR_MXR2]+0.5;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_MXR2];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_MXR2] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_MXR2] = ltmp[opt_ssr_best];
    // optimize mode radius of the 1st component
    message(stderr,"Optimizing mode radius of the 1st component (3rd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LMD1]-0.07;
    lmax = pval[OPT_PAR_LMD1]+0.07;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LMD1];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LMD1] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LMD1] = ltmp[opt_ssr_best];
    // optimize sigma of the 1st component
    message(stderr,"Optimizing sigma of the 1st component (3rd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LSG1]*0.85;
    lmax = pval[OPT_PAR_LSG1]*1.15;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LSG1];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LSG1] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LSG1] = ltmp[opt_ssr_best];
    // optimize mode radius of the 2nd component
    message(stderr,"Optimizing mode radius of the 2nd component (3rd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LMD2]-0.07;
    lmax = pval[OPT_PAR_LMD2]+0.07;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LMD2];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LMD2] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LMD2] = ltmp[opt_ssr_best];
    // optimize sigma of the 2nd component
    message(stderr,"Optimizing sigma of the 2nd component (3rd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LSG2]*0.85;
    lmax = pval[OPT_PAR_LSG2]*1.15;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LSG2];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LSG2] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LSG2] = ltmp[opt_ssr_best];
    // optimize mixing ratio of the 3rd component
    message(stderr,"Optimizing mixing ratio of the 3rd component (3rd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_MXR3]-0.5;
    lmax = pval[OPT_PAR_MXR3]+0.5;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_MXR3];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_MXR3] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_MXR3] = ltmp[opt_ssr_best];
    // optimize mode radius of the 3rd component
    message(stderr,"Optimizing mode radius of the 3rd component (3rd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LMD3]-0.07;
    lmax = pval[OPT_PAR_LMD3]+0.07;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LMD3];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LMD3] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LMD3] = ltmp[opt_ssr_best];
    // optimize sigma for each candidate
    message(stderr,"Optimizing sigma (3rd, %2d/%2d)\n",n+1,cmp_n_comp);
    lmin = pval[OPT_PAR_LSG3]*0.5;
    lmax = pval[OPT_PAR_LSG3]*1.5;
    lstp = (lmax-lmin)/5;
    ltmp[5] = pval[OPT_PAR_LSG3];
    opt_ssr_best = 5;
    for(i=0; i<5; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<5; i++)
    {
      pval[OPT_PAR_LSG3] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    pval[OPT_PAR_LSG3] = ltmp[opt_ssr_best];
    lmin = pval[OPT_PAR_LSG3]-lstp;
    lmax = pval[OPT_PAR_LSG3]+lstp;
    lstp = (lmax-lmin)/10;
    ltmp[10] = pval[OPT_PAR_LSG3];
    opt_ssr_best = 10;
    for(i=0; i<10; i++)
    {
      ltmp[i] = lmin+lstp*((double)i+0.5);
    }
    for(i=0; i<10; i++)
    {
      pval[OPT_PAR_LSG3] = ltmp[i];
      SSR_FCN_CUSTOM_B1(pval,&fmin);
      if(fmin < vmin)
      {
        vmin = fmin;
        opt_ssr_best = i;
      }
    }
    message(stderr,"opt_ssr_best = %4d, vmin = %13.4e\n",opt_ssr_best,vmin);
    opt_vssr[2][n] = vmin;
    pval[OPT_PAR_LSG3] = ltmp[opt_ssr_best];
    DSR_FCN_CUSTOM_B1(pval,&opt_vdsr[2][n]);
    // store the results
    opt_pval[2][n][OPT_PAR_VMIE] = pval[OPT_PAR_VMIE];
    opt_pval[2][n][OPT_PAR_WSCL] = pval[OPT_PAR_WSCL];
    opt_pval[2][n][OPT_PAR_NCMP] = pval[OPT_PAR_NCMP];
    opt_pval[2][n][OPT_PAR_LMD1] = pval[OPT_PAR_LMD1];
    opt_pval[2][n][OPT_PAR_LSG1] = pval[OPT_PAR_LSG1];
    opt_pval[2][n][OPT_PAR_MXR2] = pval[OPT_PAR_MXR2];
    opt_pval[2][n][OPT_PAR_LMD2] = pval[OPT_PAR_LMD2];
    opt_pval[2][n][OPT_PAR_LSG2] = pval[OPT_PAR_LSG2];
    opt_pval[2][n][OPT_PAR_MXR3] = pval[OPT_PAR_MXR3];
    opt_pval[2][n][OPT_PAR_LMD3] = pval[OPT_PAR_LMD3];
    opt_pval[2][n][OPT_PAR_LSG3] = pval[OPT_PAR_LSG3];
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_VMIE]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_VMIE]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_WSCL]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_WSCL]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_NCMP]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_NCMP]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD1]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LMD1]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG1]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LSG1]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_MXR2]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_MXR2]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD2]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LMD2]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG2]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LSG2]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_MXR3]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_MXR3]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LMD3]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LMD3]);
    message(stderr,"opt_pval[%d][%2d][OPT_PAR_LSG3]: %18.10e\n",2,n,opt_pval[2][n][OPT_PAR_LSG3]);
    message(stderr,"opt_vdsr[%d][%2d]              : %18.10e\n",2,n,opt_vdsr[2][n]);
    message(stderr,"opt_vssr[%d][%2d]              : %18.10e\n",2,n,opt_vssr[2][n]);
    message(stderr,"opt_flag[%d][%2d]              : %18d\n",   2,n,opt_flag[2][n]);
  }
  message(stderr,"3rd component summary\n");
  for(n=0; n<cmp_n_comp; n++)
  {
    message(stderr,"%d %2d",2,n);
    for(i=0; i<OPT_NPAR; i++)
    {
      message(stderr," %11.4e",opt_pval[2][n][i]);
    }
    message(stderr,"%11.4e %11.4e %3d\n",opt_vdsr[2][n],opt_vssr[2][n],opt_flag[2][n]);
  }
  // find the best fit component
  vmin = 1.0e100;
  opt_ssr_best = -1;
  for(n=0; n<cmp_n_comp; n++)
  {
    if(opt_vssr[2][n] < vmin)
    {
      vmin = opt_vssr[2][n];
      opt_ssr_best = n;
    }
  }
  message(stderr,"opt_ssr_best: %d\n",opt_ssr_best);
  opt_n_candidate = 0;
  for(n=0; n<cmp_n_comp; n++)
  {
    if((opt_vssr[2][n]-vmin)/vmin < opt_ssr_tolerance)
    {
      opt_candidate[opt_n_candidate] = n;
      opt_n_candidate++;
    }
  }
  vmin = 1.0e100;
  opt_dsr_best = -1;
  for(i=0; i<opt_n_candidate; i++)
  {
    n = opt_candidate[i];
    message(stderr,"opt_candidate[%2d]=%d, vdsr=%18.10e, vssr=%18.10e\n",i,n,opt_vdsr[2][n],opt_vssr[2][n]);
    if(opt_vdsr[2][n] < vmin)
    {
      vmin = opt_vdsr[2][n];
      opt_dsr_best = n;
    }
  }
  message(stderr,"opt_dsr_best: %d\n",opt_dsr_best);
  if((opt_vdsr[2][opt_ssr_best]-vmin)/vmin < opt_dsr_tolerance)
  {
    opt_comp[2] = opt_ssr_best;
  }
  else
  {
    opt_comp[2] = opt_dsr_best;
  }
  message(stderr,"opt_comp[%d]                  : %18d\n",2,opt_comp[2]);
  if(vmin > opt_flim)
  {
    message(stderr,"Warning, failed in finding the 3rd component >>> %18.10e\n",vmin);
    for(n=0; n<cmp_n_comp; n++)
    {
      message(stderr,"%2d %18.10e %18.10e %d\n",n,opt_vdsr[2][n],opt_vssr[2][n],opt_flag[2][n]);
    }
    pval[OPT_PAR_NCMP] = 2.0;
    DSR_FCN_CUSTOM_B0(pval,&fmin);
    DSR_FCN_CUSTOM_B1(pval,&fmin);
    DSR_FCN_CUSTOM_B2(pval,&fmin);
    SSR_FCN_CUSTOM_B1(pval,&fmin);
    return 0;
  }
  // set values from the results
  n = opt_comp[2];
  Interp1D(cmp_wlen,cmp_refr[n],cmp_n_wlen,mie_wlen_um,mie_refr_com[2],mie_n_wlen,1);
  Interp1D(cmp_wlen,cmp_refi[n],cmp_n_wlen,mie_wlen_um,mie_refi_com[2],mie_n_wlen,0);
  if(MieTable(2,2) < 0) return -1;
  pval[OPT_PAR_VMIE] = opt_pval[2][n][OPT_PAR_VMIE];
  pval[OPT_PAR_WSCL] = opt_pval[2][n][OPT_PAR_WSCL];
  pval[OPT_PAR_NCMP] = opt_pval[2][n][OPT_PAR_NCMP];
  pval[OPT_PAR_LMD1] = opt_pval[2][n][OPT_PAR_LMD1];
  pval[OPT_PAR_LSG1] = opt_pval[2][n][OPT_PAR_LSG1];
  pval[OPT_PAR_MXR2] = opt_pval[2][n][OPT_PAR_MXR2];
  pval[OPT_PAR_LMD2] = opt_pval[2][n][OPT_PAR_LMD2];
  pval[OPT_PAR_LSG2] = opt_pval[2][n][OPT_PAR_LSG2];
  pval[OPT_PAR_MXR3] = opt_pval[2][n][OPT_PAR_MXR3];
  pval[OPT_PAR_LMD3] = opt_pval[2][n][OPT_PAR_LMD3];
  pval[OPT_PAR_LSG3] = opt_pval[2][n][OPT_PAR_LSG3];

  // simulate with the best fit parameters
  DSR_FCN_CUSTOM_B0(pval,&fmin);
  DSR_FCN_CUSTOM_B1(pval,&fmin);
  DSR_FCN_CUSTOM_B2(pval,&fmin);
  SSR_FCN_CUSTOM_B1(pval,&fmin);

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
    case OPT_MODE_DSR_CUSTOM_0:
      DSR_FCN_CUSTOM_B0(xval,fcnval);
      break;
    case OPT_MODE_DSR_CUSTOM_1:
      DSR_FCN_CUSTOM_B1(xval,fcnval);
      break;
    case OPT_MODE_DSR_CUSTOM_2:
      DSR_FCN_CUSTOM_B2(xval,fcnval);
      break;
    case OPT_MODE_SSR_MODTRAN_B:
      SSR_FCN_MODTRAN_B1(xval,fcnval);
      break;
    case OPT_MODE_SSR_CUSTOM_B:
      SSR_FCN_CUSTOM_B1(xval,fcnval);
      break;
  }

  return;
}

int Printout(void)
{
  int n;

  CommonPrintout(stderr,1,0,1);
  if(mie_n_comp>0 && opt_comp[0]>=0)
  {
    n = opt_comp[0];
    fprintf(stderr,"%2d %13.5e %13.5e %d\n",n,opt_vdsr[0][n],opt_vssr[0][n],opt_flag[0][n]);
    fprintf(stderr,"%13.5e\n",opt_pval[0][n][OPT_PAR_VMIE]);
    fprintf(stderr,"%13.5e\n",opt_pval[0][n][OPT_PAR_WSCL]);
    fprintf(stderr,"%13.5e\n",opt_pval[0][n][OPT_PAR_NCMP]);
    fprintf(stderr,"%13.5e\n",opt_pval[0][n][OPT_PAR_LMD1]);
    fprintf(stderr,"%13.5e\n",opt_pval[0][n][OPT_PAR_LSG1]);
  }
  if(mie_n_comp>1 && opt_comp[1]>=0)
  {
    n = opt_comp[1];
    fprintf(stderr,"%2d %13.5e %13.5e %d\n",n,opt_vdsr[1][n],opt_vssr[1][n],opt_flag[1][n]);
    fprintf(stderr,"%13.5e\n",opt_pval[1][n][OPT_PAR_VMIE]);
    fprintf(stderr,"%13.5e\n",opt_pval[1][n][OPT_PAR_WSCL]);
    fprintf(stderr,"%13.5e\n",opt_pval[1][n][OPT_PAR_NCMP]);
    fprintf(stderr,"%13.5e\n",opt_pval[1][n][OPT_PAR_LMD1]);
    fprintf(stderr,"%13.5e\n",opt_pval[1][n][OPT_PAR_LSG1]);
    fprintf(stderr,"%13.5e\n",opt_pval[1][n][OPT_PAR_MXR2]);
    fprintf(stderr,"%13.5e\n",opt_pval[1][n][OPT_PAR_LMD2]);
    fprintf(stderr,"%13.5e\n",opt_pval[1][n][OPT_PAR_LSG2]);
  }
  if(mie_n_comp>2 && opt_comp[2]>=0)
  {
    n = opt_comp[2];
    fprintf(stderr,"%2d %13.5e %13.5e %d\n",n,opt_vdsr[2][n],opt_vssr[2][n],opt_flag[2][n]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_VMIE]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_WSCL]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_NCMP]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_LMD1]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_LSG1]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_MXR2]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_LMD2]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_LSG2]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_MXR3]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_LMD3]);
    fprintf(stderr,"%13.5e\n",opt_pval[2][n][OPT_PAR_LSG3]);
  }
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
    if(strcasecmp(str[0],"opt_fres") == 0)
    {
      if(n > 1)
      {
        strncpy(opt_fres,str[1],MAXLINE);
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30s\n","opt_fres",opt_fres);
        cnt_n_cmnt++;
      }
    } else
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

int ReadResult(void)
{
  int i,j,k;
  int n,nc;
  int err;
  int v[OPT_NINT];
  double w[OPT_NDBL];
  char line[MAXLINE];
  char temp[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char *p;
  FILE *fp;

  if(strcmp(opt_fres,"")==0 || strcmp(opt_fres,NONAME)==0)
  {
    return 0;
  } else
  if((fp=fopen(opt_fres,"r")) == NULL)
  {
    fprintf(stderr,"ReadResult: error, cannot open %s\n",opt_fres);
    return -1;
  } else
  if(cnt_vb)
  {
    fprintf(stderr,"ReadResult: reading %s\n",opt_fres);
  }
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    strncpy(temp,line,MAXLINE);
    if((p=strchr(temp,'#')) != NULL) *p = '\0';
    // read line
    for(n=nc=0,p=temp; n<MAXENTR; n++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
    }
    if(n != OPT_NCOL) continue;
    // convert values
    for(n=0; n<OPT_NINT; n++)
    {
      errno = 0;
      v[n] = strtol(str[n],&p,10);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"ReadResult: convert error >>> %s\n",str[n]);
        err = 1;
        break;
      }
    }
    if(err) break;
    for(n=OPT_NINT; n<OPT_NCOL; n++)
    {
      errno = 0;
      w[n-OPT_NINT] = strtod(str[n],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"ReadResult: convert error >>> %s\n",str[n]);
        err = 1;
        break;
      }
    }
    if(err) break;
    if(v[0]<0 || v[0]>=OPT_N_COMP)
    {
      fprintf(stderr,"ReadResult: error, v[0]=%d, OPT_N_COMP=%d\n",v[0],OPT_N_COMP);
      err = 1;
      break;
    }
    if(v[1]<0 || v[1]>=CMP_N_COMP)
    {
      fprintf(stderr,"ReadResult: error, v[1]=%d, CMP_N_COMP=%d\n",v[1],CMP_N_COMP);
      err = 1;
      break;
    }
    i = v[0];
    j = v[1];
    for(k=0; k<OPT_NPAR; k++)
    {
      opt_pval[i][j][k] = w[k];
    }
    opt_vdsr[i][j] = w[k++];
    opt_vssr[i][j] = w[k];
    opt_flag[i][j] = OPT_DONE;
  }
  fclose(fp);
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
  double xtmp;
  char *endp;
  struct option long_options[] =
  {
    {"opt_flim",1,0,'l'},
    {"opt_lmod_nstp",1,0,'L'},
  };

  fprintf(stderr,"aeros_opt.c  : $Revision: 1.58 $  $Date: 2008-08-13 15:42:42 $ (UTC)\n");

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
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) opt_flim = xtmp;
        else
        {
          fprintf(stderr,"Limit value -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'L':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) opt_lmod_nstp = ntmp;
        else
        {
          fprintf(stderr,"Log10(R) #steps for opt -> out of range %s\n",optarg);
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
  sprintf(cnt_own_format[0],"opt_fres      name         | file name(%s)\n",As(a,NONAME,8));
  sprintf(cnt_own_format[1],"opt_vmie_1st  vmie         | visibility\n");
  sprintf(cnt_own_format[2],"opt_wscl_1st  wscl         | water scaling factor\n");
  sprintf(cnt_own_format[3],"opt_vmie_skp  flag         | skip vmie optimization(%d)\n",0);
  sprintf(cnt_own_format[4],"opt_wscl_skp  flag         | skip wscl optimization(%d)\n",0);
  sprintf(cnt_own_format[5],"opt_amod_skp  flag         | skip amod optimization(%d)\n",0);

  fprintf(stderr,"aeros_opt ... fit DSR/SSR spectrum.\n");
  CommonUsage("aeros_opt",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'l':
        fprintf(stderr," l -opt_flim      |%s|%s|%s| %e\n",As(e,"Limit value",n),    As(a,"value",n),      Ae(d,OPT_FLIM,n),opt_flim);
        break;
      case 'L':
        fprintf(stderr," L -opt_lmod_nstp |%s|%s|%s| %d\n",As(e,"Log10(R) #steps",n),As(a,"#",n),          Ad(d,OPT_LMOD_NSTP,n),opt_lmod_nstp);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);

  return 0;
}
