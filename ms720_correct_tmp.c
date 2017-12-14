/*********************************************************/
/* ms720_correct_tmp                                     */
/* Author: N.Manago Nov,21,2008                          */
/* $Revision: 673 $                                      */
/* $Date: 2009-07-21 17:39:43 +0900 (Tue, 21 Jul 2009) $ */
/*********************************************************/
#include "aeros_comp.c"

// Constants for temperature correction
#define	TCR_NCOL			6			// #Columns
#define	TCR_MAXWLEN			256			// Max #wavelengths
#define	TCR_N_WLEN			TCR_MAXWLEN		// #Wavelengths
#define	TCR_N_TEMP			5			// #Temperatures
#define	TCR_N_INST			2			// #Instruments
#define	TCR_TEMP_REF			15.0			// Reference temperature
// Constants for control
#define	CNT_IMOD_OLD			1
#define	CNT_IMOD_NEW			2
// Parameters for input data
double	inp_wlen[TCR_N_INST][INP_N_WLEN];			// Wavelength
// Parameters for temperature correction
int	tcr_n_wlen			= TCR_N_WLEN;		// #Wavelengths
int	tcr_n_temp			= TCR_N_TEMP;		// #Temperatures
int	tcr_n_inst			= TCR_N_INST;		// #Instruments
double	tcr_wlen[TCR_N_INST][TCR_MAXWLEN];			// Wavelength
double	tcr_fact[TCR_N_INST][TCR_MAXWLEN][TCR_N_TEMP];		// Temperature correction factor
double	tcr_slop[TCR_N_INST][TCR_MAXWLEN][TCR_N_TEMP];		// Temperature dependence in 1/dC
double	tcr_fact_ref[TCR_N_INST][TCR_MAXWLEN][TCR_N_TEMP];	// Temperature correction factor
double	tcr_slop_ref[TCR_N_INST][TCR_MAXWLEN][TCR_N_TEMP];	// Temperature dependence in 1/dC
double	tcr_data[TCR_MAXWLEN];					// Data
double	tcr_temp_ref			= TCR_TEMP_REF;		// Reference temperature
double	tcr_temp[TCR_N_TEMP]		=
{
  -5.0,5.0,15.0,25.0,35.0,
};
double	tcr_temp_min[TCR_N_TEMP]	=
{
  -10.0,0.0,10.0,20.0,30.0,
};
double	tcr_temp_max[TCR_N_TEMP]	=
{
  0.0,10.0,20.0,30.0,40.0,
};
char	tcr_dnam[MAXLINE]		= ".";			// Calib data directory
char	tcr_fnam[TCR_N_INST][MAXLINE]	=			// Calib file
{
  "tcal_no1.dat","tcal_no2.dat",
};
// Parameters for output
char	out_dnam[MAXLINE]		= ".";			// Output directory
// Parameters for control
int	cnt_inst[TCR_N_INST]		= {0,0};		// Instrument flag
int	cnt_imod			= CNT_IMOD_NEW;		// Input mode

int CorrectTemp(void);
int ConnectLines(double tref,const double *fact,const double *slop,double *fref,double *sref);
int ReadInput3(void);
int ReadInput4(void);
int ReadTCR(void);
int WriteData(int rno,int n,const double *w,const double *d);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) return -1;
  if(cnt_hp) {Usage(); return 0;}
  if(Init() < 0) return -1;

  CorrectTemp();

  return 0;
}

int Init(void)
{
  int i,j,n;

  if(cnt_imod == CNT_IMOD_OLD)
  {
    if(ReadInput3() < 0) return -1;
  }
  else
  {
    if(ReadInput4() < 0) return -1;
  }
  if(ReadTCR()    < 0) return -1;

  for(i=0; i<tcr_n_inst; i++)
  {
    if(cnt_inst[i] == 0)
    {
      continue;
    }
    for(j=0; j<tcr_n_wlen; j++)
    {
      for(n=0; n<tcr_n_temp; n++)
      {
        tcr_fact[i][j][n] = 1.0;
      }
      if(ConnectLines(tcr_temp_ref,tcr_fact[i][j],tcr_slop[i][j],
                      tcr_fact_ref[i][j],tcr_slop_ref[i][j]) < 0) return -1;
    }
  }
  if(cnt_vb > 0)
  {
    fprintf(stderr,"Temperature correction factor:\n");
    for(i=0; i<tcr_n_inst; i++)
    {
      if(cnt_inst[i] == 0)
      {
        continue;
      }
      for(j=0; j<tcr_n_wlen; j++)
      {
        for(n=0; n<tcr_n_temp; n++)
        {
          fprintf(stderr,"%2d %4d %2d %8.4f %6.2f %11.4e\n",i,j,n,tcr_wlen[i][j],
                          tcr_temp_min[n],tcr_fact_ref[i][j][n]+tcr_slop_ref[i][j][n]*(tcr_temp_min[n]-tcr_temp[n]));
          fprintf(stderr,"%2d %4d %2d %8.4f %6.2f %11.4e\n",i,j,n,tcr_wlen[i][j],
                          tcr_temp[n],tcr_fact_ref[i][j][n]);
          fprintf(stderr,"%2d %4d %2d %8.4f %6.2f %11.4e\n",i,j,n,tcr_wlen[i][j],
                          tcr_temp_max[n],tcr_fact_ref[i][j][n]+tcr_slop_ref[i][j][n]*(tcr_temp_max[n]-tcr_temp[n]));
        }
      }
    }
  }

  return 0;
}

int CorrectTemp(void)
{
  int i,j,k,m,n;
  double t;
  double r;

  for(i=0; i<inp_n_data; i++)
  {
    t = inp_temp[i];
    k = -1;
    if(t <= tcr_temp_min[0])
    {
      k = 0;
    } else
    if(t > tcr_temp_max[tcr_n_temp-1])
    {
      k = tcr_n_temp-1;
    }
    else
    {
      for(n=0; n<tcr_n_temp; n++)
      {
        if(t>tcr_temp_min[n] && t<=tcr_temp_max[n])
        {
          k = n;
          break;
        }
      }
    }
    if(k < 0)
    {
      fprintf(stderr,"Error in finding region >>> %13.4f\n",t);
      return -1;
    }
    m = inp_ms720[i]-1;
    for(j=0; j<tcr_n_wlen; j++)
    {
      r = tcr_fact_ref[m][j][k]+tcr_slop_ref[m][j][k]*(t-tcr_temp[k]);
      tcr_data[j] = inp_data[i][j]/r;
    }
    WriteData(inp_rno[i],tcr_n_wlen,tcr_wlen[m],tcr_data);
  }

  return 0;
}

int ConnectLines(double tref,const double *fact,const double *slop,double *fref,double *sref)
{
  int i;
  int iref;
  double xi,yi;
  double xr,yr;

  // find region
  iref = -1;
  if(tref <= tcr_temp_min[0])
  {
    iref = 0;
  } else
  if(tref > tcr_temp_max[tcr_n_temp-1])
  {
    iref = tcr_n_temp-1;
  }
  else
  {
    for(i=0; i<tcr_n_temp; i++)
    {
      if(tref>tcr_temp_min[i] && tref<=tcr_temp_max[i])
      {
        iref = i;
        break;
      }
    }
  }
  if(iref < 0)
  {
    fprintf(stderr,"ConnectLines: error in finding region >>> %13.4f\n",tref);
    return -1;
  }

  // connect lines
  xr = tref;
  yr = 1.0;
  for(i=iref; i<tcr_n_temp; i++)
  {
    xi = xr;
    yi = fact[i]+slop[i]*(xi-tcr_temp[i]);
    if(fabs(yi) < EPSILON)
    {
      fprintf(stderr,"ConnectLines: error, yi=%13.4f\n",yi);
      return -1;
    }
    fref[i] = fact[i]*yr/yi;
    sref[i] = slop[i]*yr/yi;
    xr = tcr_temp_max[i];
    yr = fref[i]+sref[i]*(xr-tcr_temp[i]);
  }
  xr = tcr_temp_min[iref];
  yr = fref[iref]+sref[iref]*(xr-tcr_temp[iref]);
  for(i=iref-1; i>=0; i--)
  {
    xi = xr;
    yi = fact[i]+slop[i]*(xi-tcr_temp[i]);
    if(fabs(yi) < EPSILON)
    {
      fprintf(stderr,"ConnectLines: error, yi=%13.4f\n",yi);
      return -1;
    }
    fref[i] = fact[i]*yr/yi;
    sref[i] = slop[i]*yr/yi;
    xr = tcr_temp_min[i];
    yr = fref[i]+sref[i]*(xr-tcr_temp[i]);
  }

  return 0;
}

int ReadInput3(void)
{
  int i,j,m,n;
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
  // year month day hour min sec site_id ms720_id flag rno salt sazi airmass temperature valt vazi fov exposure par
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
    t.tm_year = v[0]-1900;
    t.tm_mon  = v[1]-1;
    t.tm_mday = v[2];
    t.tm_hour = v[3];
    t.tm_min  = v[4];
    t.tm_sec  = v[5];
    t.tm_isdst = 0;
    inp_time[i] = mktime(&t);
    inp_ms720[i]  = v[7];
    inp_flag[i]   = v[8];
    inp_rno[i]    = v[9];
    inp_th_sun[i] = 90.0-w[0];
    inp_ph_sun[i] = w[1];
    inp_amas[i]   = w[2];
    inp_temp[i]   = w[3];
    inp_th_los[i] = 90.0-w[4];
    inp_ph_los[i] = w[5];
    m = inp_ms720[i]-1;
    if(m<0 || m>=TCR_N_INST)
    {
      fprintf(stderr,"Error, invalid instrument ID >>> %d\n",inp_ms720[i]);
      err = 1;
      break;
    }
    else
    {
      cnt_inst[m] = 1;
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

  // read input data
  for(j=0; j<INP_N_WLEN; j++)
  {
    for(m=0; m<TCR_N_INST; m++)
    {
      inp_wlen[m][j] = 0.0;
    }
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
    m = inp_ms720[i]-1;
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
      if(inp_wlen[m][j] < EPSILON)
      {
        inp_wlen[m][j] = w[0];
      }
      else
      {
        if(fabs(inp_wlen[m][j]-w[0]) > DELTA)
        {
          fprintf(stderr,"Error, inp_wlen[%d][%d]=%13.4e, w[0]=%13.4e\n",m,j,inp_wlen[m][j],w[0]);
          err = 1;
          break;
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

  return 0;
}

int ReadInput4(void)
{
  int i,j,m,n;
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
    m = inp_ms720[i]-1;
    if(m<0 || m>=TCR_N_INST)
    {
      fprintf(stderr,"Error, invalid instrument ID >>> %d\n",inp_ms720[i]);
      err = 1;
      break;
    }
    else
    {
      cnt_inst[m] = 1;
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

  // read input data
  for(j=0; j<INP_N_WLEN; j++)
  {
    for(m=0; m<TCR_N_INST; m++)
    {
      inp_wlen[m][j] = 0.0;
    }
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
    m = inp_ms720[i]-1;
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
      if(inp_wlen[m][j] < EPSILON)
      {
        inp_wlen[m][j] = w[0];
      }
      else
      {
        if(fabs(inp_wlen[m][j]-w[0]) > DELTA)
        {
          fprintf(stderr,"Error, inp_wlen[%d][%d]=%13.4e, w[0]=%13.4e\n",m,j,inp_wlen[m][j],w[0]);
          err = 1;
          break;
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

  return 0;
}

int ReadTCR(void)
{
  int i,j,n;
  int err;
  double w[TCR_NCOL];
  char fnam[MAXLINE];
  char line[MAXLINE];
  char str[TCR_NCOL][MAXLINE];
  char *p;
  FILE *fp;

  for(i=0; i<tcr_n_inst; i++)
  {
    if(cnt_inst[i] == 1)
    {
      snprintf(fnam,MAXLINE,"%s/%s",tcr_dnam,tcr_fnam[i]);
      if((fp=fopen(fnam,"r")) == NULL)
      {
        fprintf(stderr,"ReadTCR: error, cannot open %s\n",fnam);
        return -1;
      }
      j = 0;
      err = 0;
      while(fgets(line,MAXLINE,fp) != NULL)
      {
        if(sscanf(line,"%s%s%s%s%s%s",str[0],str[1],str[2],str[3],str[4],str[5]) != TCR_NCOL)
        {
          fprintf(stderr,"ReadTCR: read error >>> %s\n",line);
          err = 1;
          break;
        }
        errno = 0;
        for(n=0; n<TCR_NCOL; n++)
        {
          w[n] = strtod(str[n],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"ReadTCR: convert error >>> %s\n",line);
            err = 1;
            break;
          }
        }
        if(err) break;
        if(j >= tcr_n_wlen)
        {
          fprintf(stderr,"ReadTCR: error, #data exceed the limit >>> %d\n",j);
          err = 1;
          break;
        }
        tcr_wlen[i][j] = w[0];
        for(n=0; n<tcr_n_temp; n++)
        {
          tcr_slop[i][j][n] = w[n+1];
        }
        j++;
      }
      fclose(fp);
      if(j != tcr_n_wlen)
      {
        fprintf(stderr,"ReadTCR: error, j=%d, tcr_n_wlen=%d\n",j,tcr_n_wlen);
        err = 1;
      }
      if(err)
      {
        return -1;
      }
      for(j=0; j<tcr_n_wlen; j++)
      {
        if(fabs(inp_wlen[i][j]-tcr_wlen[i][j]) > DELTA)
        {
          fprintf(stderr,"ReadTCR: error, inp_wlen[%d][%d]=%13.4e, tcr_wlen[%d][%d]=%13.4e\n",
                          i,j,inp_wlen[i][j],i,j,tcr_wlen[i][j]);
          return -1;
        }
      }
    }
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
  double xtmp;
  char *endp;
  struct option long_options[] =
  {
    {"out_dnam",1,0,'O'},
    {"tcr_dnam",1,0,'T'},
    {"tcr_fnam_1",1,0,11},
    {"tcr_fnam_2",1,0,12},
    {"tcr_temp_ref",1,0,'R'},
  };

  fprintf(stderr,"%15s : $Revision: 673 $ $Date: 2009-07-21 17:39:43 +0900 (Tue, 21 Jul 2009) $\n","ms720_correct_tmp.c");

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-O:T:R:",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'O':
        strncpy(out_dnam,optarg,MAXLINE);
        break;
      case 'T':
        strncpy(tcr_dnam,optarg,MAXLINE);
        break;
      case 11:
        strncpy(tcr_fnam[0],optarg,MAXLINE);
        break;
      case 12:
        strncpy(tcr_fnam[1],optarg,MAXLINE);
        break;
      case 'R':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') tcr_temp_ref = xtmp;
        else
        {
          fprintf(stderr,"Reference temperature -> out of range %s\n",optarg);
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
  char optstring[] = "OT12RDNCdvh";

  fprintf(stderr,"ms720_correct_tmp ... perform temperature correction on MS720 data.\n");
  CommonUsage("ms720_correct_tmp",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'O':
        fprintf(stderr," O -out_dnam      |%s|%s|%s| %s\n",As(e,"Output dir   ",n),  As(a,"name",n),       As(d,".",n),out_dnam);
        break;
      case 'T':
        fprintf(stderr," T -tcr_dnam      |%s|%s|%s| %s\n",As(e,"Calib dir    ",n),  As(a,"name",n),       As(d,".",n),tcr_dnam);
        break;
      case '1':
        fprintf(stderr,"   -tcr_fnam_1    |%s|%s|%s| %s\n",As(e,"Calib file #1",n),  As(a,"name",n),       As(d,"tcal_no1.dat",n),tcr_fnam[0]);
        break;
      case '2':
        fprintf(stderr,"   -tcr_fnam_2    |%s|%s|%s| %s\n",As(e,"Calib file #2",n),  As(a,"name",n),       As(d,"tcal_no2.dat",n),tcr_fnam[1]);
        break;
      case 'R':
        fprintf(stderr," R -tcr_temp_ref  |%s|%s|%s| %f\n",As(e,"Ref temperature",n),As(a,"dC",n),         Af(d,TCR_TEMP_REF,n),tcr_temp_ref);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  fprintf(stderr,"------------------------------------------------------------------------------------\n");

  return 0;
}
