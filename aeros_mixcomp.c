/***********************************************************/
/* AEROS_MIXCOMP ... Aerosol mixer                         */
/* Author: N.Manago                                        */
/* $Revision: 1105 $                                       */
/* $Date: 2015-09-02 12:35:25 +0900 (Wed, 02 Sep 2015) $   */
/***********************************************************/
#define	DELTA			1.0e-12				// A small number
#include "aeros_common.c"

int MieCalc(void)
{
  int i,j,k,n;
  int nc,nv,np;
  int err;
  double v[5];
  char line[MAXLINE];
  char str[MAXITEM][MAXLINE];
  char *p;
  FILE *fp;
  const char fnam[] = "MieCalc";

  for(n=0; n<mie_n_comp; n++)
  {
    if(mie_size_func[n] != MIE_FUNC_INPUT)
    {
      fprintf(stderr,"%s: error, mie_size_func[%d]=%c\n",fnam,n,mie_size_func[n]);
      return -1;
    }
    // Optical parameter 1 (w ea en o g)
    if((fp=fopen(mie_inp1[n],"r")) == NULL)
    {
      fprintf(stderr,"%s: cannot open %s\n",fnam,mie_inp1[n]);
      return -1;
    }
    np = 5;
    err = 0;
    for(i=0; i<mie_n_wlen; i++)
    {
      if(fgets(line,MAXLINE,fp) == NULL)
      {
        fprintf(stderr,"%s: read error (%s)\n",fnam,mie_inp1[n]);
        err = 1;
        break;
      }
      for(nv=0,p=line; nv<MAXITEM; nv++,p+=nc)
      {
        if(sscanf(p,"%s%n",str[nv],&nc) == EOF)
        {
          break;
        }
      }
      if(nv != np)
      {
        fprintf(stderr,"%s: read error >>> %s (%s)\n",fnam,line,mie_inp1[n]);
        err = 1;
        break;
      }
      for(nv=0; nv<np; nv++)
      {
        errno = 0;
        v[nv] = strtod(str[nv],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s (%s)\n",fnam,line,mie_inp1[n]);
          err = 1;
          break;
        }
      }
      if(err) break;
      if(fabs(v[0]/mie_wlen[i]-1.0) > DELTA)
      {
        fprintf(stderr,"%s: error, v[0]=%13.6e, mie_wlen[%d]=%13.6e (%s)\n",
                        fnam,v[0],i,mie_wlen[i],mie_inp1[n]);
        err = 1;
        break;
      }
      mie_aext_com[n][i] = v[1];
      mie_asca_com[n][i] = v[3]*v[1];
      mie_asym_com[n][i] = v[4];
    }
    if(!err)
    {
      if(fgets(line,MAXLINE,fp) != NULL)
      {
        fprintf(stderr,"%s: error, #data exceed the limit %d (%s)\n",fnam,i,mie_inp1[n]);
        err = 1;
      }
    }
    fclose(fp);
    if(err)
    {
      return -1;
    }
    // Optical parameters 2 (w a s1 s2 s3)
    if((fp=fopen(mie_inp2[n],"r")) == NULL)
    {
      fprintf(stderr,"%s: cannot open %s\n",fnam,mie_inp2[n]);
      return -1;
    }
    np = 5;
    err = 0;
    for(i=0,k=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++,k++)
      {
        if(fgets(line,MAXLINE,fp) == NULL)
        {
          fprintf(stderr,"%s: read error (%s)\n",fnam,mie_inp2[n]);
          err = 1;
          break;
        }
        for(nv=0,p=line; nv<MAXITEM; nv++,p+=nc)
        {
          if(sscanf(p,"%s%n",str[nv],&nc) == EOF)
          {
            break;
          }
        }
        if(nv != np)
        {
          fprintf(stderr,"%s: read error >>> %s (%s)\n",fnam,line,mie_inp2[n]);
          err = 1;
          break;
        }
        for(nv=0; nv<np; nv++)
        {
          errno = 0;
          v[nv] = strtod(str[nv],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error >>> %s (%s)\n",fnam,line,mie_inp2[n]);
            err = 1;
            break;
          }
        }
        if(err) break;
        if(fabs(v[0]/mie_wlen[i]-1.0) > DELTA)
        {
          fprintf(stderr,"%s: error, v[0]=%13.6e, mie_wlen[%d]=%13.6e (%s)\n",
                          fnam,v[0],i,mie_wlen[i],mie_inp2[n]);
          err = 1;
          break;
        }
        if(fabs(v[1]-mie_angl[j]) > DELTA)
        {
          fprintf(stderr,"%s: error, v[1]=%13.6e, mie_angl[%d]=%13.6e (%s)\n",
                          fnam,v[1],j,mie_angl[j],mie_inp2[n]);
          err = 1;
          break;
        }
        mie_phs1_com[n][k] = v[2];
        mie_phs2_com[n][k] = v[3];
        mie_phas_com[n][k] = v[4];
      }
      if(err) break;
    }
    if(!err)
    {
      if(fgets(line,MAXLINE,fp) != NULL)
      {
        fprintf(stderr,"%s: error, #data exceed the limit %d (%s)\n",fnam,k,mie_inp2[n]);
        err = 1;
      }
    }
    fclose(fp);
    if(err)
    {
      return -1;
    }
  }

  return 0;
}

int Init(void)
{
  int i,n;

  if(CommonInit() < 0)
  {
    return -1;
  }

  return 0;
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
  int ntmp;
  char *endp;
  double xtmp;
  struct option long_options[] =
  {
  };

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
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
  char optstring[] = "fgwlLWiIrumoOdvh";

  CommonUsage("aeros_mixcomp",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);
  CommonUsage(NULL,2);
  CommonUsage(NULL,3);
  CommonUsage(NULL,4);

  return 0;
}
