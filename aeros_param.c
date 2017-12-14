/*********************************************************/
/* aeros_param                                           */
/* Author: N.Manago Jan,31,2010                          */
/* $Revision: 711 $                                      */
/* $Date: 2010-03-23 11:40:06 +0900 (Tue, 23 Mar 2010) $ */
/*********************************************************/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "aeros_comp.h"
#include "aeros_comp.c"

// Constants for Optimization
#define	OPT_VTAU_MIN			0.0
#define	OPT_VTAU_MAX			10.0
#define	OPT_VMIE_MIN			1.0e+0
#define	OPT_VMIE_MAX			1.0e3
#define	OPT_WSCL_MIN			0.0
#define	OPT_WSCL_MAX			10.0
#define	OPT_OSCL_MIN			0.0
#define	OPT_OSCL_MAX			2.0
#define	OPT_XSCL_MIN			0.0
#define	OPT_XSCL_MAX			10.0
#define	OPT_YSCL_MIN			0.0
#define	OPT_YSCL_MAX			10.0
#define	OPT_ZSCL_MIN			0.0
#define	OPT_ZSCL_MAX			10.0
#define	OPT_NCMP_MIN			0.0
#define	OPT_NCMP_MAX			10.0
#define	OPT_MIXR_MIN			-3.0
#define	OPT_MIXR_MAX			+3.0
#define	OPT_LMOD_MIN			-3.0
#define	OPT_LMOD_MAX			+2.0
#define	OPT_LSGM_MIN			0.1
#define	OPT_LSGM_MAX			0.6
// Constants for Error estimation
#define	ERR_N_RAND			1000
#define	ERR_NCOL			6
#define	ERR_NINT			2
#define	ERR_NDBL			4
// Parameters for Optimization
int	opt_comp[OPT_N_COMP]		= {8,3,2};		// Mixture component#
// Parameters for Error estimation
int	err_n_rand			= ERR_N_RAND;
unsigned long err_seed			= 0;
double	*err_mixr_val			= NULL;
double	*err_wcom_val			= NULL;
double	*err_lmod_val			= NULL;
double	*err_lsgm_val			= NULL;
double	*err_mixr_rms			= NULL;
double	*err_wcom_rms			= NULL;
double	*err_lmod_rms			= NULL;
double	*err_lsgm_rms			= NULL;
double	*err_mixr_ofs			= NULL;
double	*err_wcom_ofs			= NULL;
double	*err_lmod_ofs			= NULL;
double	*err_lsgm_ofs			= NULL;
double	*err_aext_val			= NULL;
double	*err_omeg_val			= NULL;
double	*err_asym_val			= NULL;
double	*err_phas_val			= NULL;
double	*err_aext_rms			= NULL;
double	*err_omeg_rms			= NULL;
double	*err_asym_rms			= NULL;
double	*err_phas_rms			= NULL;
double	*err_aext_ofs			= NULL;
double	*err_omeg_ofs			= NULL;
double	*err_asym_ofs			= NULL;
double	*err_phas_ofs			= NULL;
double	err_pang_val			= 0.0;
double	err_pang_rms			= 0.0;
double	err_pang_ofs			= 0.0;
double	err_pval[OPT_NPAR];
double	err_perr[OPT_NPAR];
double	err_fmin[OPT_NPAR];
double	err_psgm[OPT_NPAR];

int GenEvents(void);
int End(void);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(cnt_hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  GenEvents();
  Printout();
  Finish();
  End();

  return 0;
}

int Init(void)
{
  int i,n;
  int nc;
  int err;
  int v[ERR_NINT];
  double w[ERR_NDBL];
  char line[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char *p;

  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,stdin) != NULL)
  {
    // read line
    for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
    }
    if(n != ERR_NCOL)
    {
      fprintf(stderr,"Error, expected %d columns >>> %d\n",ERR_NCOL,n);
      err = 1;
      break;
    }
    // convert values
    errno = 0;
    v[0] = strtol(str[0],&p,10);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",str[0]);
      err = 1;
      break;
    }
    errno = 0;
    v[1] = strtol(str[5],&p,10);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",str[5]);
      err = 1;
      break;
    }
    for(n=1; n<5; n++)
    {
      errno = 0;
      w[n-1] = strtod(str[n],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"Convert error >>> %s\n",str[n]);
        err = 1;
        break;
      }
    }
    if(err) break;
    if(i >= OPT_NPAR)
    {
      fprintf(stderr,"Error, #data exceed the limit >>> %d\n",i);
      err = 1;
      break;
    }
    err_pval[i] = w[0];
    err_perr[i] = w[1];
    err_fmin[i] = w[2];
    i++;
  }
  fprintf(stderr,"%d data have been read\n",i);
  if(i != OPT_NPAR)
  {
    fprintf(stderr,"Error, expected %d data\n",OPT_NPAR);
    err = 1;
  }
  if(err)
  {
    return -1;
  }
  for(i=0; i<OPT_NPAR; i++)
  {
    if(err_perr[i] > EPSILON)
    {
      err_psgm[i] = err_perr[i]*sqrt(err_fmin[i]);
    }
    else
    {
      err_psgm[i] = 0.0;
    }
  }

  if(MieInit() < 0) return -1;

  err_mixr_val = (double*)malloc(opt_n_comp*sizeof(double));
  err_wcom_val = (double*)malloc(opt_n_comp*sizeof(double));
  err_lmod_val = (double*)malloc(opt_n_comp*sizeof(double));
  err_lsgm_val = (double*)malloc(opt_n_comp*sizeof(double));
  err_mixr_rms = (double*)malloc(opt_n_comp*sizeof(double));
  err_wcom_rms = (double*)malloc(opt_n_comp*sizeof(double));
  err_lmod_rms = (double*)malloc(opt_n_comp*sizeof(double));
  err_lsgm_rms = (double*)malloc(opt_n_comp*sizeof(double));
  err_mixr_ofs = (double*)malloc(opt_n_comp*sizeof(double));
  err_wcom_ofs = (double*)malloc(opt_n_comp*sizeof(double));
  err_lmod_ofs = (double*)malloc(opt_n_comp*sizeof(double));
  err_lsgm_ofs = (double*)malloc(opt_n_comp*sizeof(double));
  if(err_mixr_val==NULL || err_wcom_val==NULL || err_lmod_val==NULL || err_lsgm_val==NULL ||
     err_mixr_rms==NULL || err_wcom_rms==NULL || err_lmod_rms==NULL || err_lsgm_rms==NULL ||
     err_mixr_ofs==NULL || err_wcom_ofs==NULL || err_lmod_ofs==NULL || err_lsgm_ofs==NULL)
  {
    fprintf(stderr,"Error in allocating memory.\n");
    return -1;
  }
  err_aext_val = (double*)malloc(mie_n_wlen*sizeof(double));
  err_omeg_val = (double*)malloc(mie_n_wlen*sizeof(double));
  err_asym_val = (double*)malloc(mie_n_wlen*sizeof(double));
  err_aext_rms = (double*)malloc(mie_n_wlen*sizeof(double));
  err_omeg_rms = (double*)malloc(mie_n_wlen*sizeof(double));
  err_asym_rms = (double*)malloc(mie_n_wlen*sizeof(double));
  err_aext_ofs = (double*)malloc(mie_n_wlen*sizeof(double));
  err_omeg_ofs = (double*)malloc(mie_n_wlen*sizeof(double));
  err_asym_ofs = (double*)malloc(mie_n_wlen*sizeof(double));
  if(err_aext_val==NULL || err_omeg_val==NULL || err_asym_val==NULL ||
     err_aext_rms==NULL || err_omeg_rms==NULL || err_asym_rms==NULL ||
     err_aext_ofs==NULL || err_omeg_ofs==NULL || err_asym_ofs==NULL)
  {
    fprintf(stderr,"Error in allocating memory.\n");
    return -1;
  }
  err_phas_val = (double*)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
  err_phas_rms = (double*)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
  err_phas_ofs = (double*)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
  if(err_phas_val==NULL || err_phas_rms==NULL || err_phas_ofs==NULL)
  {
    fprintf(stderr,"Error in allocating memory.\n");
    return -1;
  }

  return 0;
}

int GenEvents(void)
{
  int i,j,k,n;
  int err;
  int flag;
  int fpar[OPT_N_FPAR];
  double pval[SIM_NPAR];
  double ptmp[SIM_NPAR];
  double pmin[OPT_NPAR];
  double pmax[OPT_NPAR];
  double aext,omeg,pang;
  struct timeval  tv;
  struct timezone tz;
  const gsl_rng_type *typ;
  gsl_rng *gen;

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
    fprintf(stderr,"Error, opt_fcmp is not given.\n");
    exit(-1);
    /*
    opt_n_comp = OPT_N_COMP;
    Interp1D(cmp_wlen,cmp_refr[opt_comp[0]],cmp_n_wlen,mie_wlen_um,mie_refr_com[0],mie_n_wlen,1);
    Interp1D(cmp_wlen,cmp_refi[opt_comp[0]],cmp_n_wlen,mie_wlen_um,mie_refi_com[0],mie_n_wlen,0);
    Interp1D(cmp_wlen,cmp_refr[opt_comp[1]],cmp_n_wlen,mie_wlen_um,mie_refr_com[1],mie_n_wlen,0);
    Interp1D(cmp_wlen,cmp_refi[opt_comp[1]],cmp_n_wlen,mie_wlen_um,mie_refi_com[1],mie_n_wlen,0);
    Interp1D(cmp_wlen,cmp_refr[opt_comp[2]],cmp_n_wlen,mie_wlen_um,mie_refr_com[2],mie_n_wlen,0);
    Interp1D(cmp_wlen,cmp_refi[opt_comp[2]],cmp_n_wlen,mie_wlen_um,mie_refi_com[2],mie_n_wlen,0);
    for(n=0; n<opt_n_comp; n++)
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
    */
  }
  for(i=0; i<OPT_N_FPAR; i++)
  {
    fpar[i] = 1;
  }
  /*
  if(MieCalc(pval,fpar) < 0)
  {
    fprintf(stderr,"Error in MieCalc.\n");
    return -1;
  }
  if(MixComp(pval) < 0)
  {
    fprintf(stderr,"Error in MixComp.\n");
    return -1;
  }
  */

  // Initialize
  gsl_rng_env_setup();
  typ = gsl_rng_default;
  gen = gsl_rng_alloc(typ);
  if(err_seed == 0)
  {
    gettimeofday(&tv,&tz);
    err_seed = tv.tv_sec+tv.tv_usec;
  }
  gsl_rng_set(gen,err_seed);

  for(i=0; i<OPT_NPAR; i++)
  {
    ptmp[i] = err_pval[i];
  }
  for(i=OPT_NPAR; i<SIM_NPAR; i++)
  {
    ptmp[i] = pval[i];
  }
  if(MieCalc(ptmp,fpar) < 0)
  {
    fprintf(stderr,"Error in MieCalc.\n");
    return -1;
  }
  if(MixComp(ptmp) < 0)
  {
    fprintf(stderr,"Error in MixComp.\n");
    return -1;
  }
  printf("Original value:\n");
  for(i=0; i<opt_n_comp; i++)
  {
    printf("%d %13.6e %13.6e %13.6e %13.6e\n",i,
            mie_mixr_com[i],
            mie_wcom_com[i],
            mie_lmod_com[i],
            mie_lsgm_com[i]);
  }
  for(i=0; i<opt_n_comp; i++)
  {
    err_mixr_val[i] = 0.0;
    err_wcom_val[i] = 0.0;
    err_lmod_val[i] = 0.0;
    err_lsgm_val[i] = 0.0;
    err_mixr_rms[i] = 0.0;
    err_wcom_rms[i] = 0.0;
    err_lmod_rms[i] = 0.0;
    err_lsgm_rms[i] = 0.0;
    err_mixr_ofs[i] = mie_mixr_com[i]*mie_mixr_com[i];
    err_wcom_ofs[i] = mie_wcom_com[i]*mie_wcom_com[i];
    err_lmod_ofs[i] = mie_lmod_com[i]*mie_lmod_com[i];
    err_lsgm_ofs[i] = mie_lsgm_com[i]*mie_lsgm_com[i];
  }
  for(i=0; i<mie_n_wlen; i++)
  {
    err_aext_val[i] = 0.0;
    err_omeg_val[i] = 0.0;
    err_asym_val[i] = 0.0;
    err_aext_rms[i] = 0.0;
    err_omeg_rms[i] = 0.0;
    err_asym_rms[i] = 0.0;
    aext = mie_aext[i]/mie_aext[mie_iref];
    omeg = mie_asca[i]/mie_aext[i];
    err_aext_ofs[i] = aext*aext;
    err_omeg_ofs[i] = omeg*omeg;
    err_asym_ofs[i] = mie_asym[i]*mie_asym[i];
  }
  for(i=0; i<mie_n_wlen; i++)
  {
    for(j=0; j<mie_n_angl; j++)
    {
      k = mie_n_angl*i+j;
      err_phas_val[k] = 0.0;
      err_phas_rms[k] = 0.0;
      err_phas_ofs[k] = mie_phas[k]*mie_phas[k];
    }
  }
  err_pang_val = 0.0;
  err_pang_rms = 0.0;
  pang = -log(mie_aext[mie_iref-1]/mie_aext[mie_iref])/log(mie_wlen[mie_iref-1]/mie_wlen[mie_iref]);
  err_pang_ofs = pang*pang;

  // Generate events
  n = 0;
  err = 0;
  while(n < err_n_rand)
  {
    if(cnt_vb > 0)
    {
      message(stderr,"%d\n",n);
    }
    flag = 0;
    for(i=0; i<OPT_NPAR; i++)
    {
      if(err_psgm[i] > EPSILON)
      {
        ptmp[i] = err_pval[i]+gsl_ran_gaussian(gen,err_psgm[i]);
        if(ptmp[i]<=pmin[i] || ptmp[i]>=pmax[i])
        {
          flag = 1;
          break;
        }
      }
    }
    if(flag)
    {
      continue;
    }
    if(MieCalc(ptmp,fpar) < 0)
    {
      err = 1;
      break;
    }
    if(MixComp(ptmp) < 0)
    {
      err = 1;
      break;
    }
    for(i=0; i<opt_n_comp; i++)
    {
      err_mixr_val[i] += mie_mixr_com[i];
      err_wcom_val[i] += mie_wcom_com[i];
      err_lmod_val[i] += mie_lmod_com[i];
      err_lsgm_val[i] += mie_lsgm_com[i];
      err_mixr_rms[i] += (mie_mixr_com[i]*mie_mixr_com[i]-err_mixr_ofs[i]);
      err_wcom_rms[i] += (mie_wcom_com[i]*mie_wcom_com[i]-err_wcom_ofs[i]);
      err_lmod_rms[i] += (mie_lmod_com[i]*mie_lmod_com[i]-err_lmod_ofs[i]);
      err_lsgm_rms[i] += (mie_lsgm_com[i]*mie_lsgm_com[i]-err_lsgm_ofs[i]);
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      aext = mie_aext[i]/mie_aext[mie_iref];
      omeg = mie_asca[i]/mie_aext[i];
      err_aext_val[i] += aext;
      err_omeg_val[i] += omeg;
      err_asym_val[i] += mie_asym[i];
      err_aext_rms[i] += (aext*aext-err_aext_ofs[i]);
      err_omeg_rms[i] += (omeg*omeg-err_omeg_ofs[i]);
      err_asym_rms[i] += (mie_asym[i]*mie_asym[i]-err_asym_ofs[i]);
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        err_phas_val[k] += mie_phas[k];
        err_phas_rms[k] += (mie_phas[k]*mie_phas[k]-err_phas_ofs[k]);
      }
    }
    pang = -log(mie_aext[mie_iref-1]/mie_aext[mie_iref])/log(mie_wlen[mie_iref-1]/mie_wlen[mie_iref]);
    err_pang_val += pang;
    err_pang_rms += (pang*pang-err_pang_ofs);
    n++;
  }
  gsl_rng_free(gen);
  if(err)
  {
    return -1;
  }

  // Calculate statistics
  for(i=0; i<opt_n_comp; i++)
  {
    err_mixr_val[i] /= (double)err_n_rand;
    err_wcom_val[i] /= (double)err_n_rand;
    err_lmod_val[i] /= (double)err_n_rand;
    err_lsgm_val[i] /= (double)err_n_rand;
    err_mixr_rms[i] /= (double)err_n_rand;
    err_wcom_rms[i] /= (double)err_n_rand;
    err_lmod_rms[i] /= (double)err_n_rand;
    err_lsgm_rms[i] /= (double)err_n_rand;
    err_mixr_rms[i] = sqrt(err_mixr_rms[i]+(err_mixr_ofs[i]-err_mixr_val[i]*err_mixr_val[i]));
    err_wcom_rms[i] = sqrt(err_wcom_rms[i]+(err_wcom_ofs[i]-err_wcom_val[i]*err_wcom_val[i]));
    err_lmod_rms[i] = sqrt(err_lmod_rms[i]+(err_lmod_ofs[i]-err_lmod_val[i]*err_lmod_val[i]));
    err_lsgm_rms[i] = sqrt(err_lsgm_rms[i]+(err_lsgm_ofs[i]-err_lsgm_val[i]*err_lsgm_val[i]));
  }
  for(i=0; i<mie_n_wlen; i++)
  {
    err_aext_val[i] /= (double)err_n_rand;
    err_omeg_val[i] /= (double)err_n_rand;
    err_asym_val[i] /= (double)err_n_rand;
    err_aext_rms[i] /= (double)err_n_rand;
    err_omeg_rms[i] /= (double)err_n_rand;
    err_asym_rms[i] /= (double)err_n_rand;
    err_aext_rms[i] = sqrt(err_aext_rms[i]+(err_aext_ofs[i]-err_aext_val[i]*err_aext_val[i]));
    err_omeg_rms[i] = sqrt(err_omeg_rms[i]+(err_omeg_ofs[i]-err_omeg_val[i]*err_omeg_val[i]));
    err_asym_rms[i] = sqrt(err_asym_rms[i]+(err_asym_ofs[i]-err_asym_val[i]*err_asym_val[i]));
  }
  for(i=0; i<mie_n_wlen; i++)
  {
    for(j=0; j<mie_n_angl; j++)
    {
      k = mie_n_angl*i+j;
      err_phas_val[k] /= (double)err_n_rand;
      err_phas_rms[k] /= (double)err_n_rand;
      err_phas_rms[k] = sqrt(err_phas_rms[k]+(err_phas_ofs[k]-err_phas_val[k]*err_phas_val[k]));
    }
  }
  err_pang_val /= (double)err_n_rand;
  err_pang_rms /= (double)err_n_rand;
  err_pang_rms = sqrt(err_pang_rms+(err_pang_ofs-err_pang_val*err_pang_val));

  return 0;
}

int End(void)
{
  free(err_mixr_val);
  free(err_wcom_val);
  free(err_lmod_val);
  free(err_lsgm_val);
  free(err_mixr_rms);
  free(err_wcom_rms);
  free(err_lmod_rms);
  free(err_lsgm_rms);
  free(err_mixr_ofs);
  free(err_wcom_ofs);
  free(err_lmod_ofs);
  free(err_lsgm_ofs);
  free(err_aext_val);
  free(err_omeg_val);
  free(err_asym_val);
  free(err_aext_rms);
  free(err_omeg_rms);
  free(err_asym_rms);
  free(err_aext_ofs);
  free(err_omeg_ofs);
  free(err_asym_ofs);
  free(err_phas_val);
  free(err_phas_rms);
  free(err_phas_ofs);

  return 0;
}

int Printout(void)
{
  int i,j,k;

  for(i=0; i<OPT_NPAR; i++)
  {
    printf("%2d %13.6e %13.6e\n",i,err_pval[i],err_psgm[i]);
  }
  for(i=0; i<opt_n_comp; i++)
  {
    printf("%d %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",i,
            err_mixr_val[i],err_mixr_rms[i],
            err_wcom_val[i],err_wcom_rms[i],
            err_lmod_val[i],err_lmod_rms[i],
            err_lsgm_val[i],err_lsgm_rms[i]);
  }
  for(i=0; i<mie_n_wlen; i++)
  {
    printf("%9.4f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",mie_wlen[i],
            err_aext_val[i],err_aext_rms[i],
            err_omeg_val[i],err_omeg_rms[i],
            err_asym_val[i],err_asym_rms[i]);
  }
  for(i=0; i<mie_n_wlen; i++)
  {
    for(j=0; j<mie_n_angl; j++)
    {
      k = mie_n_angl*i+j;
      printf("%9.4f %9.4f %13.6e %13.6e\n",mie_wlen[i],mie_angl[j],
              err_phas_val[k],err_phas_rms[k]);
    }
  }
  printf("%13.6e %13.6e\n",err_pang_val,err_pang_rms);

  fflush(stdout);
  fflush(stderr);

  return 0;
}

int OwnInit(void)
{
  int n;

  for(n=0; n<OPT_NPAR; n++)
  {
    err_pval[n] = 0.0;
    err_perr[n] = 0.0;
    err_fmin[n] = 0.0;
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
    if(strcasecmp(str[0],"err_seed") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        err_seed = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30ld\n",str[0],err_seed);
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
  int ntmp;
  char *p;
  struct option long_options[] =
  {
    {"err_n_rand",1,0,'i'},
  };

  fprintf(stderr,"%15s : $Revision: 711 $ $Date: 2010-03-23 11:40:06 +0900 (Tue, 23 Mar 2010) $\n","aeros_param.c");

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-i:",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'i':
        errno = 0;
        ntmp = strtol(optarg,&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>0) err_n_rand = ntmp;
        else
        {
          fprintf(stderr,"#Events -> out of range %s\n",optarg);
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
  char optstring[] = "iDNteFUVWOsAaGgwmMrxXnSCcodvh";

  cnt_n_own_format = 3;
  sprintf(cnt_own_format[ 0],"opt_comp      # #          | component#(%d),component ID(%d)\n",0,8);
  sprintf(cnt_own_format[ 1],"err_seed      value        | Random seed(%d)\n",0);

  fprintf(stderr,"aeros_param ... calculate errors of optical parameters.\n");
  CommonUsage("aeros_param",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'i':
        fprintf(stderr," i -err_n_rand    |%s|%s|%s| %d\n",As(e,"#Events",n),        As(a,"#",n),          Ad(d,ERR_N_RAND,n),err_n_rand);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);

  return 0;
}
