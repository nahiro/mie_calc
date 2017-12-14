/*********************************************************/
/* aeros_simul                                           */
/* Author: N.Manago Apr,30,2008                          */
/* $Revision: 622 $                                      */
/* $Date: 2009-05-20 14:02:47 +0900 (Wed, 20 May 2009) $ */
/*********************************************************/
#include "aeros_comp.c"

// Constants for control
#define	CNT_DFLG			1			// DISORT flag

// Parameters for Test
int	tst_amod_1st			= INVALID;		// Aerosol model
// Parameters for control
int	cnt_dflg			= CNT_DFLG;		// DISORT flag

int main(int argc,char **argv)
{
  int i,n;
  double val;
  double pval[SIM_NPAR];

  if(GetOpt(argc,argv) < 0) exit(-1);
  if(cnt_hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  if(tst_n_comp > 0) // simulate with Custom aerosol model
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
    if(dsr_n_data > 0)
    {
      DSR_FCN_CUSTOM_B0(pval,&val);
      DSR_FCN_CUSTOM_B1(pval,&val);
      DSR_FCN_CUSTOM_B2(pval,&val);
    }
    if(aur_n_data > 0)
    {
      if(cnt_dflg == 0)
      {
        AUR_FCN_CUSTOM_B1(pval,&val);
      }
      else
      {
        AUR_FCN_CUSTOM_S1(pval,&val);
      }
    }
    if(ssr_n_data > 0)
    {
      if(cnt_dflg == 0)
      {
        SSR_FCN_CUSTOM_B1(pval,&val);
      }
      else
      {
        SSR_FCN_CUSTOM_S1(pval,&val);
      }
    }
    MiePrintout();
  }
  else // simulate with MODTRAN aerosol model
  {
    for(i=0; i<OPT_PAR_NCMP; i++)
    {
      pval[i] = opt_init[i];
    }
    pval[OPT_PAR_NCMP] = 0.0;
    if(dsr_n_data > 0)
    {
      DSR_FCN_MODTRAN_B0(pval,&val);
      DSR_FCN_MODTRAN_B1(pval,&val);
      DSR_FCN_MODTRAN_B2(pval,&val);
    }
    if(aur_n_data > 0)
    {
      if(cnt_dflg == 0)
      {
        AUR_FCN_MODTRAN_B1(pval,&val);
      }
      else
      {
        AUR_FCN_MODTRAN_S1(pval,&val);
      }
    }
    if(ssr_n_data > 0)
    {
      if(cnt_dflg == 0)
      {
        SSR_FCN_MODTRAN_B1(pval,&val);
      }
      else
      {
        SSR_FCN_MODTRAN_S1(pval,&val);
      }
    }
  }

  Printout();
  Finish();

  return 0;
}

int Init(void)
{
  int i,j,n;
  int rng;
  int rt;
  int flag;
  double ihaz;
  double pval[SIM_NPAR];

  cnt_aur_auto[1] = 1;
  if(ReadInput() < 0) return -1;
  if(dsr_n_data<1 && aur_n_data<1 && ssr_n_data<1)
  {
    fprintf(stderr,"Error, dsr_n_data=%d, aur_n_data=%d, ssr_n_data=%d\n",
                    dsr_n_data,aur_n_data,ssr_n_data);
    return -1;
  }
  if(dsr_n_data < 1)
  {
    for(rng=0; rng<DSR_MAXRANG; rng++)
    {
      cnt_dsr[rng] = 0;
    }
  }
  if(aur_n_data < 1)
  {
    for(rng=0; rng<AUR_MAXRANG; rng++)
    {
      cnt_aur[rng] = 0;
    }
  }
  if(ssr_n_data < 1)
  {
    for(rng=0; rng<SSR_MAXRANG; rng++)
    {
      cnt_ssr[rng] = 0;
    }
  }
  if(MieInit() < 0) return -1;
  if(tst_n_comp>0 || crc_n_comp>0) // simulate with Custom aerosol model
  {
    if(upper_pfunc(50.0) < 0) return -1;
  }
  if(ReadCRC() < 0) return -1;
  if((rt=ReadSimul(tst_fsim,1,1,1,1,1)) < 0) return -1;
  if(rt > 0) // replace data with simulation data
  {
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
    if(cnt_dflg == 1)
    {
      return 0;
    }
  }
  ihaz = sim_ihaz;
  if(tst_amod_1st >= 0)
  {
    sim_ihaz = tst_amod_1st;
  }
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
  if(rt > 0) // replace data with simulation data
  {
    sim_ihaz = ihaz;
    return 0;
  }
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
  sim_ihaz = ihaz;
  if(cnt_dflg == 1)
  {
    sim_cflg = 0;
    sim_dflg = 1;
  }

  return 0;
}

int Printout(void)
{
  int j;

  CommonPrintout(stderr,1,0,1);
  for(j=0; j<mie_n_comp; j++)
  {
    fprintf(stderr,"%13.5e %13.5e %13.5e\n",mie_wcom_com[j],mie_lmod_com[j],mie_lsgm_com[j]);
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
  char fnam[] = "ReadOwnConfig";

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
    if(strcasecmp(str[0],"tst_amod_1st") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_amod_1st = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,s);
          err = 1;
          break;
        }
        if(tst_amod_1st<0 || tst_amod_1st>10)
        {
          fprintf(stderr,"%s: aerosol model out of range >>> %s\n",fnam,s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30d\n",str[0],tst_amod_1st);
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

int GetOwnOpt(int argn,char **args)
{
  int c,rt;
  int option_index = 0;
  int this_option_optind;
  struct option long_options[] =
  {
  };

  fprintf(stderr,"%15s : $Revision: 622 $ $Date: 2009-05-20 14:02:47 +0900 (Wed, 20 May 2009) $\n","aeros_simul.c");

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-f",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'f':
        cnt_dflg = 0;
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
  char optstring[] = "DNteFUVWOsAaGgwmMrxXnSCcodfvh";

  cnt_n_own_format = 1;
  sprintf(cnt_own_format[0],"tst_amod_1st  ihaz         | aerosol model\n");

  fprintf(stderr,"aeros_simul ... simulate DSR/SSR spectrum.\n");
  CommonUsage("aeros_simul",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'f':
        fprintf(stderr," f -cnt_dflg      |%s|%s|%s| %d\n",As(e,"Set DISORT=f",n),   As(a,"nothing",n),    Ad(d,1-CNT_DFLG,n),1-cnt_dflg);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);

  return 0;
}
