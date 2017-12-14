/***********************************************************/
/* MIE_CALC ... Aerosol generator with MIEV0               */
/* Author: N.Manago                                        */
/* $Revision: 814 $                                        */
/* $Date: 2014-09-06 00:08:28 +0900 (Sat, 06 Sep 2014) $   */
/***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <bits/nan.h>
#include <complex.h>
#include "strutil.h"

#define	MAXLINE			256				// Max #chars in a line
#define	EPSILON			1.0e-14				// A small number
#define	DELTA			1.0e-13				// A small number
#define	MIE_TRUE		1				// True  value
#define	MIE_FALSE		0				// False value
#define	MIE_XMAX		2.0e4				// Max size parameter
#define	MIE_XVAL		1.0				// Size parameter
#define	MIE_REFR		1.53				// Real part of refractive index
#define	MIE_REFI		0.008				// Imaginal part of refractive index
#define	MIE_ICUT		1.0e-13				// Min imaginal part of refractive index
#define	MIE_NANG		361				// #Angles
#define	MIE_NMOM		0				// Highest Legendre moment
#define	MIE_DMOM		1				// First dimension of moment
#define	MIE_IPOL		0				// Legendre moment flags

double	mie_xval		= MIE_XVAL;
double	mie_refr		= MIE_REFR;
double	mie_refi		= MIE_REFI;
double	mie_icut		= MIE_ICUT;
int	mie_nang		= MIE_NANG;			// #Angles
int	mie_nmom		= MIE_NMOM;			// Highest Legendre moment
int	mie_dmom		= MIE_DMOM;			// First dimension of moment
int	mie_ipol		= MIE_IPOL;			// Legendre moment flags
int	cnt_bin			= 0;				// Binary  mode
int	cnt_db			= 0;				// Debug   mode
int	cnt_vb			= 0;				// Verbose mode
int	cnt_hp			= 0;				// Help    mode

int MieCalc(void);
int GetOpt(int argn,char **args);
int Usage(void);
extern void miev0_();

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(cnt_hp) {Usage(); return 0;}
  if(MieCalc() < 0) exit(-1);

  return 0;
}

int MieCalc(void)
{
  int i,j,imax;
  const double *XX = &mie_xval;
  double _Complex CREFIN = mie_refr-fabs(mie_refi)*I;
  int PERFCT = MIE_FALSE;
  const double *MIMCUT = &mie_icut;
  int ANYANG = MIE_FALSE;
  int *NUMANG = &mie_nang;
  double DMU;
  double *XMU;
  const int *NMOM = &mie_nmom;
  const int *IPOLZN = &mie_ipol;
  const int *MOMDIM = &mie_dmom;
  int PRNT[2];
  double QEXT,QSCA,GQSC;
  double _Complex *S1,*S2;
  double _Complex SFORW,SBACK;
  double _Complex TFORW[2],TBACK[2];
  double SPIKE;
  double *MOM;
  double *PMOM[4];

  // initialization
  XMU = (double *)malloc(mie_nang*sizeof(double));
  S1 = (double _Complex *)malloc(mie_nang*sizeof(double _Complex));
  S2 = (double _Complex *)malloc(mie_nang*sizeof(double _Complex));
  MOM = (double *)malloc((mie_dmom+1)*4*sizeof(double));
  if(XMU==NULL || S1==NULL || S2==NULL || MOM==NULL)
  {
    fprintf(stderr,"MieCalc: failed in allocating memory\n");
    return -1;
  }
  DMU = M_PI/(mie_nang-1);
  for(i=0; i<mie_nang; i++)
  {
    XMU[i] = cos(DMU*i);
  }
  if(mie_nmom > 0)
  {
    imax = (mie_dmom+1)*4;
    for(i=0; i<imax; i++)
    {
      MOM[i] = NAN;
    }
    for(i=0; i<4; i++)
    {
      PMOM[i] = &MOM[(mie_dmom+1)*i];
    }
  }
  PRNT[0] = (cnt_db>1?MIE_TRUE:MIE_FALSE);
  PRNT[1] = (cnt_db>0?MIE_TRUE:MIE_FALSE);

  // mie calculation
  miev0_(XX,&CREFIN,&PERFCT,MIMCUT,&ANYANG,NUMANG,XMU,
         NMOM,IPOLZN,MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
         MOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);

  // output
  if(cnt_bin)
  {
    fwrite(&mie_xval,sizeof(double),1,stdout);
    fwrite(&mie_refr,sizeof(double),1,stdout);
    fwrite(&mie_refi,sizeof(double),1,stdout);
    fwrite(&mie_icut,sizeof(double),1,stdout);
    fwrite(&PERFCT,sizeof(int),1,stdout);
    fwrite(&ANYANG,sizeof(int),1,stdout);
    fwrite(&QEXT,sizeof(double),1,stdout);
    fwrite(&QSCA,sizeof(double),1,stdout);
    fwrite(&GQSC,sizeof(double),1,stdout);
    fwrite(&SFORW,sizeof(double _Complex),1,stdout);
    fwrite(&SBACK,sizeof(double _Complex),1,stdout);
    fwrite(TFORW,sizeof(double _Complex),2,stdout);
    fwrite(TBACK,sizeof(double _Complex),2,stdout);
    fwrite(&SPIKE,sizeof(double),1,stdout);
    fwrite(&mie_nmom,sizeof(int),1,stdout);
    fwrite(&mie_dmom,sizeof(int),1,stdout);
    fwrite(&mie_ipol,sizeof(int),1,stdout);
    if(mie_nmom > 0)
    {
      fwrite(MOM,sizeof(double),(mie_nmom+1)*4,stdout);
    }
    fwrite(&mie_nang,sizeof(int),1,stdout);
    fwrite(S1,sizeof(double _Complex),mie_nang,stdout);
    fwrite(S2,sizeof(double _Complex),mie_nang,stdout);
  }
  else
  {
    printf("XVAL: %22.15e\n",mie_xval);
    printf("REFR: %22.15e\n",mie_refr);
    printf("REFI: %22.15e\n",mie_refi);
    printf("ICUT: %22.15e\n",mie_icut);
    printf("PFCT: %22d\n",PERFCT);
    printf("ANYA: %22d\n",ANYANG);
    printf("QEXT: %22.15e\n",QEXT);
    printf("QSCA: %22.15e\n",QSCA);
    printf("GQSC: %22.15e\n",GQSC);
    printf("SFRW: %22.15e %22.15e\n",creal(SFORW),cimag(SFORW));
    printf("SBCK: %22.15e %22.15e\n",creal(SBACK),cimag(SBACK));
    printf("TFW1: %22.15e %22.15e\n",creal(TFORW[0]),cimag(TFORW[0]));
    printf("TFW2: %22.15e %22.15e\n",creal(TFORW[1]),cimag(TFORW[1]));
    printf("TBK1: %22.15e %22.15e\n",creal(TBACK[0]),cimag(TBACK[0]));
    printf("TBK2: %22.15e %22.15e\n",creal(TBACK[1]),cimag(TBACK[1]));
    printf("SPIK: %22.15e\n",SPIKE);
    printf("NMOM: %22d\n",mie_nmom);
    printf("DMOM: %22d\n",mie_dmom);
    printf("IPOL: %22d\n",mie_ipol);
    if(mie_nmom > 0)
    {
      imax = mie_nmom;
      for(j=0; j<4; j++)
      {
        printf("PMOM%d:\n",j+1);
        for(i=0; i<=mie_nmom; i++)
        {
          printf("%22.15e%s",PMOM[j][i],(i%8==7?"\n":(i==imax?"\n":" ")));
        }
      }
    }
    printf("NANG: %22d\n",mie_nang);
    imax = mie_nang-1;
    printf("S1(real):\n");
    for(i=0; i<mie_nang; i++)
    {
      printf("%22.15e%s",creal(S1[i]),(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("S1(imag):\n");
    for(i=0; i<mie_nang; i++)
    {
      printf("%22.15e%s",cimag(S1[i]),(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("S2(real):\n");
    for(i=0; i<mie_nang; i++)
    {
      printf("%22.15e%s",creal(S2[i]),(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("S2(imag):\n");
    for(i=0; i<mie_nang; i++)
    {
      printf("%22.15e%s",cimag(S2[i]),(i%8==7?"\n":(i==imax?"\n":" ")));
    }
  }

  // cleanup
  free(XMU);
  free(S1);
  free(S2);
  free(MOM);

  return 0;
}

int GetOpt(int argn,char **args)
{
  int c,rt;
  int option_index = 0;
  int this_option_optind;
  int ntmp;
  double xtmp;
  char *endp;
  struct option long_options[] =
  {
    {"mie_xval",1,0,'x'},
    {"mie_refr",1,0,'R'},
    {"mie_refi",1,0,'I'},
    {"mie_icut",1,0,'i'},
    {"mie_nang",1,0,'M'},
    {"mie_nmom",1,0,'N'},
    {"mie_ipol",1,0,'p'},
    {"binary",0,0,'b'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"x:R:I:i:M:N:p:bdvh",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'x':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') mie_xval = xtmp;
        else
        {
          fprintf(stderr,"X value -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'R':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') mie_refr = xtmp;
        else
        {
          fprintf(stderr,"Refractive index (real) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'I':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') mie_refi = xtmp;
        else
        {
          fprintf(stderr,"Refractive index (imag) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'i':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') mie_icut = xtmp;
        else
        {
          fprintf(stderr,"Min ref index (imag) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'M':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) mie_nang = ntmp;
        else
        {
          fprintf(stderr,"#Angles -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'N':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0) mie_nmom = ntmp;
        else
        {
          fprintf(stderr,"#Moments -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'p':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && abs(ntmp)<=4444) mie_ipol = ntmp;
        else
        {
          fprintf(stderr,"Moment flags -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'b':
        cnt_bin = 1;
        break;
      case 'd':
        cnt_db++;
        break;
      case 'v':
        cnt_vb++;
        break;
      case 'h':
        cnt_hp = 1;
        break;
      case '?':
        if(optopt == '\0')
        {
          fprintf(stderr,"Invalid option : %s\n",args[this_option_optind]);
        }
        else
        {
          fprintf(stderr,"Invalid option : %c\n",optopt);
        }
        rt = -1;
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

  mie_dmom = (mie_nmom<1?1:mie_nmom);

  if(optind < argn)
  {
    fprintf(stderr,"non-option ARGV-elements:\n");
    while(optind < argn)
    {
      fprintf(stderr,"%s\n",args[optind++]);
    }
    rt = -1;
  }

  if(cnt_hp) rt = 1;
  return rt;
}

int Usage(void)
{
  int  n = 16;
  char e[MAXLINE];
  char a[MAXLINE];
  char d[MAXLINE];

  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"mie_calc -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"------------------------------------------------------------------------------------\n");
  fprintf(stderr,"   option         |%s|%s|%s| current\n",As(e,"",n),          As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," x -mie_xval      |%s|%s|%s| %f\n",As(e,"X Value",n),        As(a,"value",n),      Af(d,MIE_XVAL,n),mie_xval);
  fprintf(stderr," R -mie_refr      |%s|%s|%s| %f\n",As(e,"Re(Ref. index)",n), As(a,"value",n),      Af(d,MIE_REFR,n),mie_refr);
  fprintf(stderr," I -mie_refi      |%s|%s|%s| %e\n",As(e,"Im(Ref. index)",n), As(a,"value",n),      Ae(d,MIE_REFI,n),mie_refi);
  fprintf(stderr," i -mie_icut      |%s|%s|%s| %e\n",As(e,"Min ref. index",n), As(a,"value",n),      Ae(d,MIE_ICUT,n),mie_icut);
  fprintf(stderr," M -mie_nang      |%s|%s|%s| %d\n",As(e,"#Angles",n),        As(a,"#",n),          Ad(d,MIE_NANG,n),mie_nang);
  fprintf(stderr," N -mie_nmom      |%s|%s|%s| %d\n",As(e,"#Moments",n),       As(a,"#",n),          Ad(d,MIE_NMOM,n),mie_nmom);
  fprintf(stderr," p -mie_ipol      |%s|%s|%s| %d\n",As(e,"Moment flags",n),   As(a,"#",n),          Ad(d,MIE_IPOL,n),mie_ipol);
  fprintf(stderr," b -binary        |%s|%s|%s| %d\n",As(e,"Binary  mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_bin);
  fprintf(stderr," d -debug         |%s|%s|%s| %d\n",As(e,"Debug   mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_db);
  fprintf(stderr," v -verbose       |%s|%s|%s| %d\n",As(e,"Verbose mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_vb);
  fprintf(stderr," h -help          |%s|%s|%s| %d\n",As(e,"Help    mode",n),   As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"------------------------------------------------------------------------------------\n");

  return 0;
}
