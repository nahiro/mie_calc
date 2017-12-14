/*********************************************************/
/* ms720_correct_fov                                     */
/* Author: N.Manago Oct,28,2008                          */
/* $Revision: 679 $                                      */
/* $Date: 2009-09-22 22:01:39 +0900 (Tue, 22 Sep 2009) $ */
/*********************************************************/
#include "aeros_comp.c"

// Constants for control
#define	CNT_DFLG		1				// DISORT flag

// Parameters for Test
int		tst_amod_1st		= -1;			// Aerosol model
double		tst_vmie_1st		= NAN;			// Visibility in km
double		tst_wscl_1st		= NAN;			// Water scale
// Parameters for output
char		out_dnam[MAXLINE]	= ".";			// Output directory

int CorrectFOV(void);
int ReadInput2(void);
int WriteData(int rno,int n,const double *w,const double *d);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) return -1;
  if(cnt_hp) {Usage(); return 0;}
  if(Init() < 0) return -1;

  CorrectFOV();
  Finish();

  return 0;
}

int Init(void)
{
  if(ReadInput2() < 0) return -1;
  if(MieInit() < 0) return -1;
  if(ReadCRC() < 0) return -1;
  if(upper_pfunc(50.0) < 0) return -1;

  return 0;
}

int CorrectFOV(void)
{
  int i,j,n;
  int rng = 1;
  int flag_1,flag_2,flag_3;
  int tmp_type;
  time_t tmp_time;
  double tmp_th_los;
  double tmp_ph_los;
  double pval[SIM_NPAR];
  double etmp;

  // set parameters
  for(i=0; i<OPT_PAR_NCMP; i++)
  {
    pval[i] = opt_init[i];
  }
  if(!isnan(tst_vmie_1st))
  {
    pval[OPT_PAR_VMIE] = tst_vmie_1st;
  }
  if(!isnan(tst_wscl_1st))
  {
    pval[OPT_PAR_WSCL] = tst_wscl_1st;
  }
  if(tst_amod_1st >= 0)
  {
    sim_ihaz = tst_amod_1st;
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
    flag_1 = 1;
    flag_2 = 2;
  }
  else
  {
    flag_1 = 0;
    flag_2 = 0;
  }

  dsr_n_data = 0;
  aur_n_data = 0;
  ssr_n_data = 0;
  tmp_type = -1;
  tmp_time = -1;
  tmp_th_los = 1.0e100;
  tmp_ph_los = 1.0e100;
  for(i=0; i<inp_n_data; i++)
  {
    if(inp_type[i] == INP_TYPE_DSR) // DSR
    {
      if(tmp_type==inp_type[i] &&
         abs(tmp_time-inp_time[i])<inp_tdif)
      {
        flag_3 = 4;
      }
      else
      {
        flag_3 = 0;
      }
      dsr_n_data = 1;
      cnt_dsr[rng] = 1;
      dsr_time[0]   = inp_time[i]-inp_time[0];
      dsr_th_sun[0] = inp_th_sun[i];
      dsr_ph_sun[0] = inp_ph_sun[i];
      dsr_amas[0]   = inp_amas[i];
      for(j=0; j<dsr_n_wlen[rng]; j++)
      {
        dsr_data[rng][0][j] = inp_data[i][dsr_nw[rng][j]];
        dsr_dorg[rng][0][j] = inp_data[i][dsr_nw[rng][j]];
        dsr_derr[rng][0][j] = 1.0; // temporary
      }
      if(flag_3 == 4)
      {
        if(Background(pval,flag_3) < 0) return -1;
      }
      else
      {
        if(CalcFactor(pval,flag_1,1) < 0) return -1;
        sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
        sim_cflg = 1; // Apply DIS=f->t correction factor from here
        sim_dflg = 0; // Set DIS=f from here
        if(Background(pval,flag_2) < 0) return -1;
      }
      WriteData(inp_rno[i],dsr_n_wlen[rng],dsr_wlen[rng],dsr_data[rng][0]);
      dsr_n_data = 0;
      cnt_dsr[rng] = 0;
    } else
    if(inp_type[i] == INP_TYPE_AUR) // AUR
    {
      if(tmp_type==inp_type[i] &&
         abs(tmp_time-inp_time[i])<inp_tdif)
      {
        flag_3 = 4;
      }
      else
      {
        flag_3 = 0;
      }
      aur_n_data = 1;
      cnt_aur[rng] = 1;
      aur_time[0]   = inp_time[i]-inp_time[0];
      aur_th_sun[0] = inp_th_sun[i];
      aur_ph_sun[0] = inp_ph_sun[i];
      aur_amas[0]   = inp_amas[i];
      if(SepToAzi(aur_th_sun[0]*D_TO_R,aur_rfov*D_TO_R,&etmp) < 0)
      {
        fprintf(stderr,"Error in AUR geometry.\n");
        return -1;
      }
      aur_th_los[0] = aur_th_sun[0];
      aur_ph_los[0] = aur_ph_sun[0]-etmp*R_TO_D;
      for(j=0; j<aur_n_wlen[rng]; j++)
      {
        aur_data[rng][0][j] = inp_data[i][aur_nw[rng][j]]; // temporary (not correct if aur_type==AUR_TYPE_OLD)
        aur_dorg[rng][0][j] = inp_data[i][aur_nw[rng][j]]; // temporary (not correct if aur_type==AUR_TYPE_OLD)
        aur_derr[rng][0][j] = 1.0; // temporary
      }
      if(flag_3 == 4)
      {
        if(RingCenter(pval,flag_3) < 0) return -1;
      }
      else
      {
        if(CalcFactor(pval,flag_1,1) < 0) return -1;
        sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
        sim_cflg = 1; // Apply DIS=f->t correction factor from here
        sim_dflg = 0; // Set DIS=f from here
        if(RingCenter(pval,flag_2) < 0) return -1;
      }
      WriteData(inp_rno[i],aur_n_wlen[rng],aur_wlen[rng],aur_data[rng][0]);
      aur_n_data = 0;
      cnt_aur[rng] = 0;
    } else
    if(inp_type[i] == INP_TYPE_SSR) // SSR
    {
      if(tmp_type==inp_type[i] &&
         abs(tmp_time-inp_time[i])<inp_tdif &&
         fabs(tmp_th_los-inp_th_los[i])<DELTA &&
         fabs(tmp_ph_los-inp_ph_los[i])<DELTA)
      {
        flag_3 = 4;
      }
      else
      {
        flag_3 = 0;
      }
      ssr_n_data = 1;
      cnt_ssr[rng] = 1;
      ssr_time[0]   = inp_time[i]-inp_time[0];
      ssr_th_sun[0] = inp_th_sun[i];
      ssr_ph_sun[0] = inp_ph_sun[i];
      ssr_amas[0]   = inp_amas[i];
      ssr_th_los[0] = inp_th_los[i];
      ssr_ph_los[0] = inp_ph_los[i];
      for(j=0; j<ssr_n_wlen[rng]; j++)
      {
        ssr_data[rng][0][j] = inp_data[i][ssr_nw[rng][j]];
        ssr_dorg[rng][0][j] = inp_data[i][ssr_nw[rng][j]];
        ssr_derr[rng][0][j] = 1.0; // temporary
      }
      if(flag_3 == 4)
      {
        if(Uniformity(pval,flag_3) < 0) return -1;
      }
      else
      {
        if(CalcFactor(pval,flag_1,1) < 0) return -1;
        sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
        sim_cflg = 1; // Apply DIS=f->t correction factor from here
        sim_dflg = 0; // Set DIS=f from here
        if(Uniformity(pval,flag_2) < 0) return -1;
      }
      WriteData(inp_rno[i],ssr_n_wlen[rng],ssr_wlen[rng],ssr_data[rng][0]);
      ssr_n_data = 0;
      cnt_ssr[rng] = 0;
    }
    tmp_type   = inp_type[i];
    tmp_time   = inp_time[i];
    tmp_th_los = inp_th_los[i];
    tmp_ph_los = inp_ph_los[i];
  }

  return 0;
}

int ReadInput2(void)
{
  int i,j,k,m,n;
  int rng = 1;
  int nc;
  int err;
  int v[INP_NINT];
  double w[INP_NDBL];
  char line[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char fnam[MAXLINE];
  char format[MAXLINE];
  char *p;
  FILE *fp;
  struct tm t;

  // read stdin
  // format:
  // rno yr mo dy hr mi sc site_id inst_id cloud flag fov vazi valt par exposure temperature sazi salt airmass distance
  putenv("TZ=");
  t.tm_isdst = 0;
  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,stdin) != NULL)
  {
    // read line
    for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
    }
    if(n != INP_NCOL)
    {
      fprintf(stderr,"Error, expected %d columns >>> %d\n",INP_NCOL,n);
      err = 1;
      break;
    }
    // convert values
    for(n=0; n<INP_NINT; n++)
    {
      errno = 0;
      v[n] = strtol(str[n],&p,10);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"Convert error >>> %s\n",str[n]);
        err = 1;
        break;
      }
    }
    if(err) break;
    for(n=INP_NINT; n<INP_NCOL; n++)
    {
      errno = 0;
      w[n-INP_NINT] = strtod(str[n],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"Convert error >>> %s\n",str[n]);
        err = 1;
        break;
      }
    }
    if(err) break;
    if(i >= INP_MAXDATA)
    {
      fprintf(stderr,"Warning, #data exceed the limit >>> %d\n",i);
      break;
    }
    inp_rno[i]    = v[0];
    t.tm_year     = v[1]-1900;
    t.tm_mon      = v[2]-1;
    t.tm_mday     = v[3];
    t.tm_hour     = v[4];
    t.tm_min      = v[5];
    t.tm_sec      = v[6];
    t.tm_isdst    = 0;
    inp_time[i]   = mktime(&t);
    inp_ms720[i]  = v[8];
    inp_cloud[i]  = v[9];
    inp_flag[i]   = v[10];
    inp_ph_los[i] = w[1];
    inp_th_los[i] = 90.0-w[2];
    inp_temp[i]   = w[5];
    inp_ph_sun[i] = w[6];
    inp_th_sun[i] = 90.0-w[7];
    inp_amas[i]   = w[8];
    inp_dist[i]   = w[9];
    if(inp_flag[i]>=100 && inp_flag[i]<200) // DSR
    {
      inp_type[i] = INP_TYPE_DSR;
      if(dsr_ms720 < 0)
      {
        dsr_ms720 = inp_ms720[i];
      } else
      if(inp_ms720[i] != dsr_ms720)
      {
        fprintf(stderr,"Error, conflicting MS720ID for DSR.\n");
        err = 1;
        break;
      }
    } else
    if(inp_flag[i]>=300 && inp_flag[i]<400) // AUR (old method)
    {
      inp_type[i] = INP_TYPE_AUR;
      if(aur_ms720 < 0)
      {
        aur_ms720 = inp_ms720[i];
      } else
      if(inp_ms720[i] != aur_ms720)
      {
        fprintf(stderr,"Error, conflicting MS720ID for AUR.\n");
        err = 1;
        break;
      }
      if(aur_type < 0)
      {
        aur_type = AUR_TYPE_OLD;
      } else
      if(aur_type != AUR_TYPE_OLD)
      {
        fprintf(stderr,"Error, conflicting AUR type.\n");
        err = 1;
        break;
      }
    } else
    if(inp_flag[i]>=500 && inp_flag[i]<600) // AUR
    {
      inp_type[i] = INP_TYPE_AUR;
      if(aur_ms720 < 0)
      {
        aur_ms720 = inp_ms720[i];
      } else
      if(inp_ms720[i] != aur_ms720)
      {
        fprintf(stderr,"Error, conflicting MS720ID for AUR.\n");
        err = 1;
        break;
      }
      if(aur_type < 0)
      {
        aur_type = AUR_TYPE_NEW;
      } else
      if(aur_type != AUR_TYPE_NEW)
      {
        fprintf(stderr,"Error, conflicting aur type.\n");
        err = 1;
        break;
      }
    } else
    if(inp_flag[i]>=0 && inp_flag[i]<100)   // SSR
    {
      inp_type[i] = INP_TYPE_SSR;
      if(ssr_ms720 < 0)
      {
        ssr_ms720 = inp_ms720[i];
      } else
      if(inp_ms720[i] != ssr_ms720)
      {
        fprintf(stderr,"Error, conflicting MS720ID for SSR.\n");
        err = 1;
        break;
      }
    }
    else
    {
      inp_type[i] = INP_TYPE_UNKNOWN;
    }
    i++;
  }
  inp_n_data = i;
  fprintf(stderr,"%d data have been read\n",inp_n_data);
  if(inp_n_data < 1)
  {
    fprintf(stderr,"Error, #data is not enough.\n");
    err = 1;
  }
  if(err)
  {
    return -1;
  }
  if(sim_iday < 0)
  {
    localtime_r(&inp_time[0],&t);
    sim_iday = t.tm_yday+1;
  }

  // read input data
  for(j=0; j<INP_N_WLEN; j++)
  {
    inp_dsr_wlen[j] = 0.0;
    inp_aur_wlen[j] = 0.0;
    inp_ssr_wlen[j] = 0.0;
  }
  sprintf(format,"%%s/%%0%dd.dat",inp_ndig);
  for(i=0; i<inp_n_data; i++)
  {
    sprintf(fnam,format,inp_dnam,inp_rno[i]);
    if(cnt_vb > 0)
    {
      fprintf(stderr,"Reading %s\n",fnam);
    }
    if((fp=fopen(fnam,"r")) == NULL)
    {
      fprintf(stderr,"Error, cannot open %s\n",fnam);
      return -1;
    }
    err = 0;
    j = 0;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      if(sscanf(line,"%s%s",str[0],str[1]) != 2)
      {
        fprintf(stderr,"Read error >>> %s\n",line);
        err = 1;
        break;
      }
      errno = 0;
      for(n=0; n<2; n++)
      {
        w[n] = strtod(str[n],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"Convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(err) break;
      if(j >= INP_N_WLEN)
      {
        fprintf(stderr,"Error, #data exceed the limit >>> %d\n",j);
        err = 1;
        break;
      }
      if(inp_type[i] == INP_TYPE_DSR) // DSR
      {
        if(inp_dsr_wlen[j] < EPSILON)
        {
          inp_dsr_wlen[j] = w[0];
        }
        else
        {
          if(fabs(inp_dsr_wlen[j]-w[0]) > DELTA)
          {
            fprintf(stderr,"Error, inp_dsr_wlen[%d]=%13.4e, w[0]=%13.4e\n",j,inp_dsr_wlen[j],w[0]);
            err = 1;
            break;
          }
        }
      } else
      if(inp_type[i] == INP_TYPE_AUR) // AUR
      {
        if(inp_aur_wlen[j] < EPSILON)
        {
          inp_aur_wlen[j] = w[0];
        }
        else
        {
          if(fabs(inp_aur_wlen[j]-w[0]) > DELTA)
          {
            fprintf(stderr,"Error, inp_aur_wlen[%d]=%13.4e, w[0]=%13.4e\n",j,inp_aur_wlen[j],w[0]);
            err = 1;
            break;
          }
        }
      } else
      if(inp_type[i] == INP_TYPE_SSR) // SSR
      {
        if(inp_ssr_wlen[j] < EPSILON)
        {
          inp_ssr_wlen[j] = w[0];
        }
        else
        {
          if(fabs(inp_ssr_wlen[j]-w[0]) > DELTA)
          {
            fprintf(stderr,"Error, inp_ssr_wlen[%d]=%13.4e, w[0]=%13.4e\n",j,inp_ssr_wlen[j],w[0]);
            err = 1;
            break;
          }
        }
      }
      inp_data[i][j] = w[1];
      j++;
    }
    fclose(fp);
    if(err)
    {
      return -1;
    }
    if(j != INP_N_WLEN)
    {
      fprintf(stderr,"Error, i=%d j=%d\n",i,j);
      return -1;
    }
  }

  // set wavelengths
  for(i=0; i<DSR_MAXRANG; i++)
  {
    cnt_dsr[i] = 0;
    if(i != rng)
    {
      dsr_n_wlen[i] = 0;
    }
  }
  for(i=0; i<AUR_MAXRANG; i++)
  {
    cnt_aur[i] = 0;
    if(i != rng)
    {
      aur_n_wlen[i] = 0;
    }
  }
  for(i=0; i<SSR_MAXRANG; i++)
  {
    cnt_ssr[i] = 0;
    if(i != rng)
    {
      ssr_n_wlen[i] = 0;
    }
  }
  if(dsr_ms720 > 0)
  {
    if(dsr_n_wlen[rng] < 1)
    {
      dsr_n_wlen[rng] = INP_N_WLEN;
      for(j=0; j<dsr_n_wlen[rng]; j++)
      {
        dsr_nw[rng][j] = j;
        dsr_wlen[rng][j] = inp_dsr_wlen[j];
      }
    }
    else
    {
      for(j=0; j<dsr_n_wlen[rng]; j++)
      {
        dsr_nw[rng][j] = -1;
        for(k=0; k<INP_N_WLEN; k++)
        {
          if(fabs(dsr_wlen[rng][j]-inp_dsr_wlen[k]) < EPSILON)
          {
            dsr_nw[rng][j] = k;
            break;
          }
        }
        if(dsr_nw[rng][j] < 0)
        {
          fprintf(stderr,"Error, faild in finding DSR wavelength >>> %13.4f\n",dsr_wlen[rng][j]);
          return -1;
        }
      }
    }
  }
  else
  {
    dsr_n_wlen[rng] = 0;
  }
  if(aur_ms720 > 0)
  {
    if(aur_n_wlen[rng] < 1)
    {
      aur_n_wlen[rng] = INP_N_WLEN;
      for(j=0; j<aur_n_wlen[rng]; j++)
      {
        aur_nw[rng][j] = j;
        aur_wlen[rng][j] = inp_aur_wlen[j];
      }
    }
    else
    {
      for(j=0; j<aur_n_wlen[rng]; j++)
      {
        aur_nw[rng][j] = -1;
        for(k=0; k<INP_N_WLEN; k++)
        {
          if(fabs(aur_wlen[rng][j]-inp_aur_wlen[k]) < EPSILON)
          {
            aur_nw[rng][j] = k;
            break;
          }
        }
        if(aur_nw[rng][j] < 0)
        {
          fprintf(stderr,"Error, faild in finding AUR wavelength >>> %13.4f\n",aur_wlen[rng][j]);
          return -1;
        }
      }
    }
    if(aur_n_sbmp[rng] != aur_n_wlen[rng])
    {
      if(sim_sbmp_init < 0)
      {
        i = aur_ms720-1;
        if(i<0 || i>=INS_N_INST)
        {
          fprintf(stderr,"Error, invalid instrument ID for AUR >>> %d\n",aur_ms720);
          return -1;
        }
        for(j=0; j<aur_n_wlen[rng]; j++)
        {
          if(fabs(aur_wlen[rng][j]-ins_wlen[i][j]) > DELTA)
          {
            fprintf(stderr,"Error, aur_wlen[%d][%d]=%13.4e, ins_wlen[%d][%d]=%13.4e\n",
                            rng,j,aur_wlen[rng][j],i,j,ins_wlen[i][j]);
            return -1;
          }
          m = 0;
          for(n=1; n<rcr_n_band; n++)
          {
            if(rcr_rms[i][j][n] < rcr_rms[i][j][m])
            {
              m = n;
            }
          }
          aur_sbmp[rng][j] = rcr_band[m];
        }
      }
      else
      {
        for(j=0; j<aur_n_wlen[rng]; j++)
        {
          aur_sbmp[rng][j] = sim_sbmp_init;
        }
      }
      aur_n_sbmp[rng] = aur_n_wlen[rng];
    }
  }
  else
  {
    aur_n_wlen[rng] = 0;
  }
  if(ssr_ms720 > 0)
  {
    if(ssr_n_wlen[rng] < 1)
    {
      ssr_n_wlen[rng] = INP_N_WLEN;
      for(j=0; j<ssr_n_wlen[rng]; j++)
      {
        ssr_nw[rng][j] = j;
        ssr_wlen[rng][j] = inp_ssr_wlen[j];
      }
    }
    else
    {
      for(j=0; j<ssr_n_wlen[rng]; j++)
      {
        ssr_nw[rng][j] = -1;
        for(k=0; k<INP_N_WLEN; k++)
        {
          if(fabs(ssr_wlen[rng][j]-inp_ssr_wlen[k]) < EPSILON)
          {
            ssr_nw[rng][j] = k;
            break;
          }
        }
        if(ssr_nw[rng][j] < 0)
        {
          fprintf(stderr,"Error, faild in finding SSR wavelength >>> %13.4f\n",ssr_wlen[rng][j]);
          return -1;
        }
      }
    }
    if(ssr_n_sbmp[rng] != ssr_n_wlen[rng])
    {
      if(sim_sbmp_init < 0)
      {
        i = ssr_ms720-1;
        if(i<0 || i>=INS_N_INST)
        {
          fprintf(stderr,"Error, invalid instrument ID for SSR >>> %d\n",ssr_ms720);
          return -1;
        }
        for(j=0; j<ssr_n_wlen[rng]; j++)
        {
          if(fabs(ssr_wlen[rng][j]-ins_wlen[i][j]) > DELTA)
          {
            fprintf(stderr,"Error, ssr_wlen[%d][%d]=%13.4e, ins_wlen[%d][%d]=%13.4e\n",
                            rng,j,ssr_wlen[rng][j],i,j,ins_wlen[i][j]);
            return -1;
          }
          m = 0;
          for(n=1; n<rcr_n_band; n++)
          {
            if(rcr_rms[i][j][n] < rcr_rms[i][j][m])
            {
              m = n;
            }
          }
          ssr_sbmp[rng][j] = rcr_band[m];
        }
      }
      else
      {
        for(j=0; j<ssr_n_wlen[rng]; j++)
        {
          ssr_sbmp[rng][j] = sim_sbmp_init;
        }
      }
      ssr_n_sbmp[rng] = ssr_n_wlen[rng];
    }
  }
  else
  {
    ssr_n_wlen[rng] = 0;
  }

  return 0;
}

int WriteData(int rno,int n,const double *w,const double *d)
{
  int i;
  char fnam[MAXLINE];
  char format[MAXLINE];
  FILE *fp;

  sprintf(format,"%%s/%%0%dd.dat",inp_ndig);
  sprintf(fnam,format,out_dnam,rno);
  if((fp=fopen(fnam,"w")) == NULL)
  {
    fprintf(stderr,"WriteData: cannot open %s\n",fnam);
    return -1;
  }
  for(i=0; i<n; i++)
  {
    fprintf(fp,"%8.3f %12.5e\n",w[i],d[i]);
  }
  fclose(fp);

  return 0;
}

int Printout(void)
{
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
    if(strcasecmp(str[0],"tst_vmie_1st") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_vmie_1st = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadOwnConfig: convert error >>> %s\n",s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n","tst_vmie_1st",tst_vmie_1st);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_wscl_1st") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_wscl_1st = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadOwnConfig: convert error >>> %s\n",s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n","tst_wscl_1st",tst_wscl_1st);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_amod_1st") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tst_amod_1st = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadOwnConfig: convert error >>> %s\n",s);
          err = 1;
          break;
        }
        if(tst_amod_1st<0 || tst_amod_1st>10)
        {
          fprintf(stderr,"ReadOwnConfig: aerosol model out of range >>> %s\n",s);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30d\n","tst_amod_1st",tst_amod_1st);
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
    {"out_dnam",1,0,'O'},
  };

  fprintf(stderr,"%15s : $Revision: 679 $ $Date: 2009-09-22 22:01:39 +0900 (Tue, 22 Sep 2009) $\n","ms720_correct_fov.c");

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-O:",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'O':
        strncpy(out_dnam,optarg,MAXLINE);
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
  char optstring[] = "ODNteFVWAaGgwmMrxXnsSCcodvh";

  cnt_n_own_format = 3;
  sprintf(cnt_own_format[0],"tst_vmie_1st  vmie         | visibility\n");
  sprintf(cnt_own_format[1],"tst_wscl_1st  wscl         | water scaling factor\n");
  sprintf(cnt_own_format[2],"tst_amod_1st  ihaz         | aerosol model\n");

  fprintf(stderr,"ms720_correct_fov ... perform fov correction on MS720 data.\n");
  CommonUsage("ms720_correct_fov",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'O':
        fprintf(stderr," O -out_dnam      |%s|%s|%s| %s\n",As(e,"Output dir",n),     As(a,"name",n),       As(d,".",n),out_dnam);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);

  return 0;
}
