#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include "aeros_comp.c"

// Parameters for Test
char	mie_dnam[MAXLINE]		= NONAME;		// directory name

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) return -1;
  if(cnt_hp) {Usage(); return 0;}
  if(Init() < 0) return -1;

  return 0;
}

int Init(void)
{
  if(MieInit() < 0) return -1;

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
  return 0;
}

int GetOwnOpt(int argn,char **args)
{
  int c,rt;
  int option_index = 0;
  int this_option_optind;
  struct option long_options[] =
  {
    {"mie_dnam",0,0,'f'},
  };

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-f:",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'f':
        strncpy(mie_dnam,optarg,MAXLINE);
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
  char optstring[] = "fDNteFUVWOsAaGgwmMrxXnSCcodvh";

  cnt_n_own_format = 1;
  sprintf(cnt_own_format[0],"tst_amod_1st  ihaz         | aerosol model\n");

  fprintf(stderr,"aeros_mie ... execute Mie calculation.\n");
  CommonUsage("aeros_mie",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'f':
        fprintf(stderr," f -mie_dnam      |%s|%s|%s| %s\n",As(e,"Data directory",n), As(a,"name",n),       As(d,NONAME,n),mie_dnam);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);

  return 0;
}
