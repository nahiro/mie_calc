#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <bits/nan.h>
#include <complex.h>
#include <time.h>
#include <sys/time.h>
#include "strutil.h"
#include "aeros_comp.c"

#define		MAXNINP			100

int		ninp[MIE_MAXCOMP];
double		linp[MIE_MAXCOMP][MAXNINP];
char		finp[MAXLINE];
// Parameters for control
int		cnt_lg			= 0;			// Log     mode
int		cnt_ex			= 0;			// Exp     mode
int		cnt_rv			= 0;			// Reverse mode

int main(int argc,char **argv)
{
  int i,n;
  double lout;

  if(GetOpt(argc,argv) < 0) exit(-1);
  if(cnt_hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  for(n=0; n<tst_n_comp; n++)
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      mie_refr_com[n][i] = tst_refr_com[n][i];
      mie_refi_com[n][i] = tst_refi_com[n][i];
    }
    if(MieTable(n,n) < 0) return -1;
    if(MieReff(n,n,tst_lsgm_com[n]) < 0) return -1;
    for(i=0; i<ninp[n]; i++)
    {
      if(cnt_rv)
      {
        lout = mod2eff(n,linp[n][i]);
      }
      else
      {
        lout = eff2mod(n,linp[n][i]);
      }
      if(cnt_ex)
      {
        printf("%2d %13.4e %13.4e %13.4e %13.4e\n",n,linp[n][i],lout,pow(10.0,linp[n][i]),pow(10.0,lout));
      }
      else
      {
        printf("%2d %13.6f %13.6f %13.6f %13.6f\n",n,linp[n][i],lout,pow(10.0,linp[n][i]),pow(10.0,lout));
      }
    }
  }

  return 0;
}

int Init(void)
{
  int n;
  int err;
  int flag = 1;
  double v;
  char line[MAXLINE];
  char str1[MAXLINE];
  char str2[MAXLINE];
  char *p;
  FILE *fp;

  if(MieInit() < 0) return -1;
  if(tst_n_comp<1 || tst_n_comp>MIE_MAXCOMP)
  {
    fprintf(stderr,"Error, tst_n_comp=%d\n",tst_n_comp);
    return -1;
  }

  for(n=0; n<MIE_MAXCOMP; n++)
  {
    ninp[n] = 0;
  }
  if(strcmp(finp,"") != 0)
  {
    if((fp=fopen(finp,"r")) == NULL)
    {
      fprintf(stderr,"Error, cannot open %s\n",finp);
      return -1;
    }
  }
  else
  {
    fp = stdin;
    flag = 0;
  }
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    if(sscanf(line,"%s%s",str1,str2) != 2)
    {
      fprintf(stderr,"Read error >>> %s\n",line);
      err = 1;
      break;
    }
    errno = 0;
    n = strtol(str1,&p,10);
    if(errno==ERANGE || *p!='\0' || n<0 || n>=MIE_MAXCOMP)
    {
      fprintf(stderr,"Component# out of range >>> %s\n",line);
      err = 1;
      break;
    }
    errno = 0;
    v = strtod(str2,&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      err = 1;
      break;
    }
    if(ninp[n] >= MAXNINP)
    {
      fprintf(stderr,"Warning, #data exceed the limit (%d) >>> %d\n",n,ninp[n]);
      continue;
    }
    if(cnt_lg)
    {
      linp[n][ninp[n]] = log10(v);
    }
    else
    {
      linp[n][ninp[n]] = v;
    }
    ninp[n]++;
  }
  if(flag)
  {
    fclose(fp);
  }
  if(err)
  {
    return -1;
  }

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
  struct option long_options[] =
  {
    {"finp",1,0,'f'},
    {"log",0,0,'l'},
    {"exp",0,0,'e'},
    {"reverse",0,0,'R'},
  };

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-f:leR",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'f':
        strncpy(finp,optarg,MAXLINE);
        break;
      case 'l':
        cnt_lg = 1;
        break;
      case 'e':
        cnt_ex = 1;
        break;
      case 'R':
        cnt_rv = 1;
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
  char optstring[] = "srxXnSCfleRdvh";

  fprintf(stderr,"mie_eff2mod ... effective radius <-> mode radius conversion program.\n");
  CommonUsage("mie_eff2mod",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'f':
        fprintf(stderr," f -finp          |%s|%s|%s| %s\n",As(e,"Input file",n),     As(a,"name",n),       As(d,"",n),finp);
        break;
      case 'l':
        fprintf(stderr," l -log           |%s|%s|%s| %d\n",As(e,"Log     mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_lg);
        break;
      case 'e':
        fprintf(stderr," e -exp           |%s|%s|%s| %d\n",As(e,"Exp     mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_ex);
        break;
      case 'R':
        fprintf(stderr," R -reverse       |%s|%s|%s| %d\n",As(e,"Reverse mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_rv);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  fprintf(stderr,"------------------------------------------------------------------------------------\n");
  fprintf(stderr,"Config file (cnt_conf) format:\n");
  fprintf(stderr,"mie_wlen      name # u m M | file name(%s),wlen column#(%d),unit in nm(%.1f),min wlen(%.1f),max wlen(%.1f)\n",
                                               NONAME,0,1.0,MIE_WLEN_MIN,MIE_WLEN_MAX);
  fprintf(stderr,"tst_fcmp    name # # u # # | file name(%s),real column#(%d),imag column#(%d),wlen unit in nm(%.1f),"
                                               "min line#(%d),max line#(10^%.0f)\n",
                                               NONAME,MIE_REAL_NUM,MIE_IMAG_NUM,1.0,MIE_IMIN,log10((double)MIE_IMAX));
  fprintf(stderr,"Input data (finp or stdin) format: n v\n");
  fprintf(stderr,"  n ... Component#\n");
  fprintf(stderr,"  v ... log10(r) (default), r (Log mode)\n");

  return 0;
}
