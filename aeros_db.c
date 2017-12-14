/***********************************************************/
/* AEROS_DB ... Aerosol generator with DB                  */
/* Author: N.Manago                                        */
/* $Revision: 1104 $                                       */
/* $Date: 2015-09-02 00:16:06 +0900 (Wed, 02 Sep 2015) $   */
/***********************************************************/
#include "aeros_common.c"

int	cnt_norun		= 0;				// No run  mode

int ReadTape6(void);
int ReadTape8(void);

int MieCalc(void)
{
  int i,n;
  FILE *fp;

  if(!cnt_norun)
  {
    if((fp=fopen("entry.dat","w")) == NULL)
    {
      fprintf(stderr,"Error, cannot open %s\n","entry.dat");
      return -1;
    }
    // User defined class
    fprintf(fp,"Type of class\n");
    fprintf(fp,"43\n\n");
    // No scaling
    fprintf(fp,"Scaling radii factor\n");
    fprintf(fp,"4\n\n");
    // Leave the mixing ration as they are in literature
    fprintf(fp,"User mixing ratio\n");
    fprintf(fp,"0\n\n");
    // External mixing
    fprintf(fp,"Type of mixing\n");
    fprintf(fp,"2\n\n");
    // RH 0%
    fprintf(fp,"Relative humidity level\n");
    fprintf(fp,"1\n\n");
    // Aerosol parameters
    fprintf(fp,"User parameters\n");
    fprintf(fp,"Original Model\n");
    fprintf(fp,"%d %d\n",mie_n_comp,mie_n_wlen);
    for(n=0; n<mie_n_comp; n++)
    {
      fprintf(fp,"%11.4e %11.4e %11.4e\n",mie_wcom_com[n],mie_rmod_com[n],mie_lsgm_com[n]);
      for(i=0; i<mie_n_wlen; i++)
      {
        fprintf(fp,"%11.4e %11.4e %11.4e\n",mie_wlen_um[i],mie_refr_com[n][i],fabs(mie_refi_com[n][i]));
      }
    }
    fprintf(fp,"\n");
    // Problem number
    fprintf(fp,"Problem number\n");
    fprintf(fp,"1\n\n");
    // Angle
    fprintf(fp,"Phase function tabulation\n");
    fprintf(fp,"3\n\n");
    fprintf(fp,"Number abscissae\n");
    fprintf(fp,"%d\n\n",mie_n_angl);
    fprintf(fp,"Abscissae\n");
    for(i=0; i<mie_n_angl; i++)
    {
      fprintf(fp," %9.5f%s",mie_angl[i],(i==mie_n_angl-1?"\n":i%8==7?"\n":""));
    }
    fprintf(fp,"\n");
    // Legendre polynomial
    fprintf(fp,"Legendre polynomial\n");
    fprintf(fp,"0\n\n");
    fprintf(fp,"Terms\n");
    fprintf(fp,"100\n\n");
    // Wavelength
    fprintf(fp,"Type of wavelenghts\n");
    fprintf(fp,"1\n\n");
    fprintf(fp,"Number wavelengths\n");
    fprintf(fp,"%d\n\n",mie_n_wlen);
    fprintf(fp,"Wavelength values\n");
    for(i=0; i<mie_n_wlen; i++)
    {
      fprintf(fp," %7.5f%s",mie_wlen_um[i],(i==mie_n_wlen-1?"\n":i%10==9?"\n":""));
    }
    fclose(fp);
    system("db");
  }

  // read results
  if(ReadTape6() < 0) return -1;
  if(ReadTape8() < 0) return -1;

  return 0;
}

int ReadTape6(void)
{
  int i,j,n;
  double v[4];
  char line[4][MAXLINE];
  char str[4][MAXLINE];
  char *p;
  FILE *fp;

  if((fp=fopen("mie6.out","r")) == NULL)
  {
    fprintf(stderr,"ReadTape6: cannot open %s\n","mie6.out");
    return -1;
  }

  // read components
  n = 0;
  i = 0;
  while(fgets(line[0],MAXLINE,fp) != NULL)
  {
    if(sscanf(line[0],"%s",str[0]) != 1) continue;
    if(strcmp(str[0],"EXT,SCAT,ABS,ASYM") == 0)
    {
      if(i == mie_n_wlen)
      {
        i = 0;
      }
      if(i == 0)
      {
        n++;
      }
      if((p=strchr(line[0],'=')) == NULL)
      {
        fprintf(stderr,"ReadTape6: error in finding \"=\" >>> %s\n",line[j]);
        fclose(fp);
        return -1;
      }
      if(sscanf(p+1,"%s%s%s%s",str[0],str[1],str[2],str[3]) != 4)
      {
        fprintf(stderr,"ReadTape6: scan error >>> %s\n",line[0]);
        fclose(fp);
        return -1;
      }
      for(j=0; j<4; j++)
      {
        errno = 0;
        v[j] = strtod(str[j],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadTape6: convert error >>> %s\n",line[0]);
        }
      }
      if(n>mie_n_comp || i>=mie_n_wlen)
      {
        fprintf(stderr,"ReadTape6: error, n=%d, i=%d\n",n,i);
        fclose(fp);
        return -1;
      }
      mie_aext_com[n-1][i] = v[0];
      mie_asca_com[n-1][i] = v[1];
      mie_asym_com[n-1][i] = v[3];
      i++;
    }
    if(n==mie_n_comp && i==mie_n_wlen) break;
  }
  if(n!=mie_n_comp || i!=mie_n_wlen)
  {
    fprintf(stderr,"ReadTape6: error, n=%d, mie_n_comp=%d, i=%d, mie_n_wlen=%d\n",n,mie_n_comp,i,mie_n_wlen);
    fclose(fp);
    return -1;
  }

  // read mixture
  i = 0;
  while(fgets(line[0],MAXLINE,fp) != NULL)
  {
    if(sscanf(line[0],"%s",str[0]) != 1) continue;
    if(strcmp(str[0],"WAVELENGTH") == 0)
    {
      if(sscanf(line[0],"%s%s%s%s",str[0],str[1],str[2],str[3]) > 3) continue;
      for(j=1; j<4; j++)
      {
        if(fgets(line[j],MAXLINE,fp) == NULL)
        {
          fprintf(stderr,"ReadTape6: read error (%d)\n",j);
          fclose(fp);
          return -1;
        }
      }
      for(j=0; j<4; j++)
      {
        if((p=strchr(line[j],'=')) == NULL)
        {
          fprintf(stderr,"ReadTape6: error in finding \"=\" (%d) >>> %s\n",j,line[j]);
          fclose(fp);
          return -1;
        }
        if(sscanf(p+1,"%s",str[0]) != 1)
        {
          fprintf(stderr,"ReadTape6: error in reading value (%d) >>> %s\n",j,line[j]);
          fclose(fp);
          return -1;
        }
        errno = 0;
        v[j] = strtod(str[0],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadTape6: convert error (%d) >>> %s\n",j,line[j]);
        }
      }
      if(i >= mie_n_wlen)
      {
        fprintf(stderr,"ReadTape6: error, i=%d\n",i);
        fclose(fp);
        return -1;
      }
      if(fabs(v[0]-mie_wlen_um[i]) > EPSILON)
      {
        fprintf(stderr,"ReadTape6: error, v[0]=%13.4e, mie_wlen_um[%d]=%13.4e\n",v[0],i,mie_wlen_um[i]);
        fclose(fp);
        return -1;
      }
      mie_aext[i] = v[1];
      mie_asca[i] = v[1]*v[2];
      mie_asym[i] = v[3];
      i++;
    }
  }
  fclose(fp);
  if(i != mie_n_wlen)
  {
    fprintf(stderr,"ReadTape6: error, i=%d, mie_n_wlen=%d\n",i,mie_n_wlen);
    return -1;
  }

  return 0;
}

int ReadTape8(void)
{
  int i,j,k,n;
  double v[5];
  char line[MAXLINE];
  char str[5][MAXLINE];
  char *p;
  FILE *fp;

  if((fp=fopen("mie8.out","r")) == NULL)
  {
    fprintf(stderr,"ReadTape8: cannot open %s\n","mie8.out");
    return -1;
  }
  for(n=0; n<mie_n_comp; n++)
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++)
      {
        if(fgets(line,MAXLINE,fp) == NULL)
        {
          fprintf(stderr,"ReadTape8: read error (%d,%d,%d)\n",n,i,j);
          fclose(fp);
          return -1;
        }
        if(sscanf(line,"%s%s%s%s",str[0],str[1],str[2],str[3]) != 4)
        {
          fprintf(stderr,"ReadTape8: scan error (%d,%d,%d) >>> %s\n",n,i,j,line);
          fclose(fp);
          return -1;
        }
        for(k=0; k<4; k++)
        {
          errno = 0;
          v[k] = strtod(str[k],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"ReadTape8: convert error (%d,%d,%d) >>> %s\n",n,i,j,line);
            fclose(fp);
            return -1;
          }
        }
        if(fabs(v[0]-mie_angl[j]) > EPSILON)
        {
          fprintf(stderr,"ReadTape8: error, v[0]=%13.4e, mie_angl[%d]=%13.4e\n",v[0],j,mie_angl[j]);
          fclose(fp);
          return -1;
        }
        k = mie_n_angl*i+j;
        mie_phs1_com[n][k] = v[1];
        mie_phs2_com[n][k] = v[2];
        mie_phas_com[n][k] = 0.5*(v[1]+v[2]);
      }
    }
  }

  return 0;
}

int Init(void)
{
  return CommonInit();
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
  struct option long_options[] =
  {
    {"cnt_cmod",0,0,'c'},
    {"cnt_norun",0,0,'N'},
  };

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-cN",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'c':
        cnt_cmod = 1;
        break;
      case 'N':
        cnt_norun = 1;
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
  char optstring[] = "fgwlLaWAPQiIruNmcdvh";

  CommonUsage("aeros_db",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'c':
        fprintf(stderr," c -cnt_cmod      |%s|%s|%s| %d\n",As(e,"Calc.   mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_cmod);
        break;
      case 'N':
        fprintf(stderr," N -cnt_norun     |%s|%s|%s| %d\n",As(e,"No run  mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_norun);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);

  return 0;
}
