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
#define	DELTA			1.0e-12				// A small number
#define	MIE_TRUE		1				// True  value
#define	MIE_FALSE		0				// False value
#define	MIE_LMIN		-4.0				// Log(Min size parameter)
#define	MIE_LMAX		4.0				// Log(Max size parameter)
#define	MIE_LSTP		0.0001				// Log(size parameter) step
#define	MIE_REFR		1.53				// Real part of refractive index
#define	MIE_REFI		0.008				// Imaginal part of refractive index
#define	MIE_ICUT		1.0e-13				// Min imaginal part of refractive index
#define	MIE_NANG		361				// #Angles
#define	MIE_NMOM		0				// Highest Legendre moment
#define	MIE_DMOM		1				// First dimension of moment
#define	MIE_IPOL		0				// Legendre moment flags

char	fnam_qext[MAXLINE]	= "qext.dat";
char	fnam_qsca[MAXLINE]	= "qsca.dat";
char	fnam_p1[MAXLINE]	= "p1.dat";
char	fnam_p2[MAXLINE]	= "p2.dat";
double	mie_lmin		= MIE_LMIN;
double	mie_lmax		= MIE_LMAX;
double	mie_lstp		= MIE_LSTP;
double	mie_refr		= MIE_REFR;
double	mie_refi		= MIE_REFI;
double	mie_icut		= MIE_ICUT;
int	mie_nang		= MIE_NANG;			// #Angles
int	mie_nmom		= MIE_NMOM;			// Highest Legendre moment
int	mie_dmom		= MIE_DMOM;			// First dimension of moment
int	mie_ipol		= MIE_IPOL;			// Legendre moment flags
int	cnt_db			= 0;				// Debug   mode
int	cnt_vb			= 0;				// Verbose mode
int	cnt_hp			= 0;				// Help    mode

int GetOpt(int argn,char **args);
int Usage(void);
int ReadMie(void);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(cnt_hp) {Usage(); return 0;}
  if(ReadMie() < 0) exit(-1);

  return 0;
}

int ReadMie(void)
{
  int i,j,imax;
  int err;
  double l;
  double mie_xval;
  double XX;
  double _Complex CREFIN;
  int PERFCT;
  double MIMCUT;
  int ANYANG;
  int NUMANG;
  int NMOM;
  int IPOLZN;
  int MOMDIM;
  double QEXT,QSCA,GQSC;
  double _Complex *S1,*S2;
  double _Complex SFORW,SBACK;
  double _Complex TFORW[2],TBACK[2];
  double SPIKE;
  double *MOM;
  double *PMOM[4];
  double *p1,*p2;
  FILE *fp_qext;
  FILE *fp_qsca;
  FILE *fp_p1;
  FILE *fp_p2;

  // allocate memory
  S1 = (double _Complex *)malloc(mie_nang*sizeof(double _Complex));
  S2 = (double _Complex *)malloc(mie_nang*sizeof(double _Complex));
  p1 = (double *)malloc(mie_nang*sizeof(double));
  p2 = (double *)malloc(mie_nang*sizeof(double));
  MOM = (double *)malloc((mie_dmom+1)*4*sizeof(double));
  if(S1==NULL || S2==NULL || p1==NULL || p2==NULL || MOM==NULL)
  {
    fprintf(stderr,"ReadMie: failed in allocating memory.\n");
    return -1;
  }
  if(mie_nmom > 0)
  {
    for(i=0; i<4; i++)
    {
      PMOM[i] = &MOM[(mie_dmom+1)*i];
    }
  }

  // open file
  fp_qext = fopen(fnam_qext,"w");
  fp_qsca = fopen(fnam_qsca,"w");
  fp_p1 = fopen(fnam_p1,"w");
  fp_p2 = fopen(fnam_p2,"w");
  if(fp_qext==NULL || fp_qsca==NULL || fp_p1==NULL || fp_p2==NULL)
  {
    fprintf(stderr,"ReadMie: failed in opening file.\n");
    free(S1);
    free(S2);
    free(p1);
    free(p2);
    free(MOM);
    return -1;
  }

  // read Mie
  for(l=mie_lmin; l<mie_lmax; l+=mie_lstp)
  {
    mie_xval = pow(10.0,l);
    fread(&XX,sizeof(double),1,stdin);
    fread(&CREFIN,sizeof(double _Complex),1,stdin);
    fread(&MIMCUT,sizeof(double),1,stdin);
    fread(&PERFCT,sizeof(int),1,stdin);
    fread(&ANYANG,sizeof(int),1,stdin);
    fread(&QEXT,sizeof(double),1,stdin);
    fread(&QSCA,sizeof(double),1,stdin);
    fread(&GQSC,sizeof(double),1,stdin);
    fread(&SFORW,sizeof(double _Complex),1,stdin);
    fread(&SBACK,sizeof(double _Complex),1,stdin);
    fread(TFORW,sizeof(double _Complex),2,stdin);
    fread(TBACK,sizeof(double _Complex),2,stdin);
    fread(&SPIKE,sizeof(double),1,stdin);
    fread(&NMOM,sizeof(int),1,stdin);
    fread(&MOMDIM,sizeof(int),1,stdin);
    fread(&IPOLZN,sizeof(int),1,stdin);
    if(mie_nmom > 0)
    {
      fread(MOM,sizeof(double),(mie_nmom+1)*4,stdin);
    }
    fread(&NUMANG,sizeof(int),1,stdin);
    fread(S1,sizeof(double _Complex),mie_nang,stdin);
    fread(S2,sizeof(double _Complex),mie_nang,stdin);
    err = 0;
    if(fabs(1.0-XX/mie_xval) > DELTA)
    {
      fprintf(stderr,"Error, 1-XX/mie_xval=%13.6e\n",1.0-XX/mie_xval);
      err = 1;
    }
    if(creal(CREFIN) != mie_refr || cimag(CREFIN) != mie_refi)
    {
      fprintf(stderr,"Error, CREFIN=%13.6e-%13.6ei, mie_refr=%13.6e, mie_refi=%13.6e\n",creal(CREFIN),cimag(CREFIN),mie_refr,mie_refi);
      err = 1;
    }
    if(MIMCUT != mie_icut)
    {
      fprintf(stderr,"Error, MIMCUT=%13.6e, mie_icut=%13.6e\n",MIMCUT,mie_icut);
      err = 1;
    }
    if(NUMANG != mie_nang)
    {
      fprintf(stderr,"Error, NUMANG=%d, mie_nang=%d\n",NUMANG,mie_nang);
      err = 1;
    }
    if(NMOM != mie_nmom)
    {
      fprintf(stderr,"Error, NMOM=%d, mie_nmom=%d\n",NMOM,mie_nmom);
      err = 1;
    }
    if(IPOLZN != mie_ipol)
    {
      fprintf(stderr,"Error, IPOLZN=%d, mie_ipol=%d\n",IPOLZN,mie_ipol);
      err = 1;
    }
    if(MOMDIM != mie_dmom)
    {
      fprintf(stderr,"Error, MOMDIM=%d, mie_dmom=%d\n",MOMDIM,mie_dmom);
      err = 1;
    }
    if(err != 0)
    {
      return -1;
    }
    if(cnt_vb)
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
    for(i=0; i<mie_nang; i++)
    {
      p1[i] = creal(S1[i]*conj(S1[i]));
      p2[i] = creal(S2[i]*conj(S2[i]));
    }
    fwrite(&QEXT,sizeof(double),1,fp_qext);
    fwrite(&QSCA,sizeof(double),1,fp_qsca);
    fwrite(p1,sizeof(double),mie_nang,fp_p1);
    fwrite(p2,sizeof(double),mie_nang,fp_p2);
  }

  // cleanup
  free(S1);
  free(S2);
  free(p1);
  free(p2);
  free(MOM);
  fclose(fp_qext);
  fclose(fp_qsca);
  fclose(fp_p1);
  fclose(fp_p2);
  
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
    {"mie_lmin",1,0,'l'},
    {"mie_lmax",1,0,'L'},
    {"mie_lstp",1,0,'s'},
    {"mie_refr",1,0,'R'},
    {"mie_refi",1,0,'I'},
    {"mie_icut",1,0,'i'},
    {"mie_nang",1,0,'M'},
    {"mie_nmom",1,0,'N'},
    {"mie_ipol",1,0,'p'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"l:L:s:R:I:i:M:N:p:bdvh",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'l':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') mie_lmin = xtmp;
        else
        {
          fprintf(stderr,"Log(X min) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'L':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') mie_lmax = xtmp;
        else
        {
          fprintf(stderr,"Log(X max) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 's':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') mie_lstp = xtmp;
        else
        {
          fprintf(stderr,"Log(X) step -> out of range %s\n",optarg);
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
  fprintf(stderr,"read_mie -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"------------------------------------------------------------------------------------\n");
  fprintf(stderr,"   option         |%s|%s|%s| current\n",As(e,"",n),          As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," l -mie_lmin      |%s|%s|%s| %f\n",As(e,"Log(X min)",n),     As(a,"value",n),      Af(d,MIE_LMIN,n),mie_lmin);
  fprintf(stderr," L -mie_lmax      |%s|%s|%s| %f\n",As(e,"Log(X max)",n),     As(a,"value",n),      Af(d,MIE_LMAX,n),mie_lmax);
  fprintf(stderr," s -mie_lstp      |%s|%s|%s| %f\n",As(e,"Log(X) step",n),    As(a,"value",n),      Af(d,MIE_LSTP,n),mie_lstp);
  fprintf(stderr," R -mie_refr      |%s|%s|%s| %f\n",As(e,"Re(Ref. index)",n), As(a,"value",n),      Af(d,MIE_REFR,n),mie_refr);
  fprintf(stderr," I -mie_refi      |%s|%s|%s| %e\n",As(e,"Im(Ref. index)",n), As(a,"value",n),      Ae(d,MIE_REFI,n),mie_refi);
  fprintf(stderr," i -mie_icut      |%s|%s|%s| %e\n",As(e,"Min ref. index",n), As(a,"value",n),      Ae(d,MIE_ICUT,n),mie_icut);
  fprintf(stderr," M -mie_nang      |%s|%s|%s| %d\n",As(e,"#Angles",n),        As(a,"#",n),          Ad(d,MIE_NANG,n),mie_nang);
  fprintf(stderr," N -mie_nmom      |%s|%s|%s| %d\n",As(e,"#Moments",n),       As(a,"#",n),          Ad(d,MIE_NMOM,n),mie_nmom);
  fprintf(stderr," p -mie_ipol      |%s|%s|%s| %d\n",As(e,"Moment flags",n),   As(a,"#",n),          Ad(d,MIE_IPOL,n),mie_ipol);
  fprintf(stderr," d -debug         |%s|%s|%s| %d\n",As(e,"Debug   mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_db);
  fprintf(stderr," v -verbose       |%s|%s|%s| %d\n",As(e,"Verbose mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_vb);
  fprintf(stderr," h -help          |%s|%s|%s| %d\n",As(e,"Help    mode",n),   As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"------------------------------------------------------------------------------------\n");

  return 0;
}
