#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include "strutil.h"

#define	NPN1			200
#define	NPNG1			(3*NPN1)
#define	NPNG2			(2*NPNG1)
#define	NPN2			(2*NPN1)
#define	NPL			(NPN2+1)
#define	NPN3			(NPN1+1)
#define	NPN4			(NPN1-25)
#define	NPN5			(2*NPN4)
#define	NPN6			(NPN4+1)
#define	NPL1			(NPN5+1)
#define	NAXMAX			10
#define	NAMAX			500
#define	NAFMAX			(NAMAX*6)
#define	NALMAX			(NAXMAX*NPL)
#define	NARMAX			(NAXMAX*NAFMAX)
#define	NANG			361
#define	REFR			1.53				// Refractive index (real)
#define	REFI			0.008				// Refractive index (imag)
#define	WLEN			0.550				// Wavelength in um
#define	RMIN			1.0e-3				// R min   in um
#define	RMAX			1.0e+2				// R max   in um
#define	RMOD			0.5				// R mode  in um
#define	LSGM			0.3				// R sigma in um (log)
#define	ARAT			1.000001			// Aspect ratio
#define	MAXLINE			256				// Max #chars in a line
#define	EPSILON			1.0e-14				// A small number

extern void tmatrix_();

int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);
int TmxCalc(void);

int	nang			= NANG;
int	db			= 0;
int	vb			= 0;
int	hp			= 0;
double	refr			= REFR;
double	refi			= REFI;
double	wlen			= WLEN;
double	rmin			= RMIN;
double	rmax			= RMAX;
double	rmod			= RMOD;
double	lsgm			= LSGM;
double	arat			= ARAT;

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) return -1;
  if(hp) {Usage(); return 0;}
  if(Init() < 0) return -1;

  if(TmxCalc() < 0) return -1;

  return 0;
}

int TmxCalc(void)
{
  int i,j,k;
  int iflg,ierr;
  int ndistr = 2;
  int npnax = 1;
  int nkmax = 5;
  int np = -1;
  int ndgs = 2;
  int lmax[NAXMAX];
  double rat = 0.5;
  double b = 0.1;
  double gam = 0.5;
  double ddelt = 0.001;
  double reff[NAXMAX];
  double veff[NAXMAX];
  double cext[NAXMAX];
  double csca[NAXMAX];
  double walb[NAXMAX];
  double asym[NAXMAX];
  double alp1[NALMAX];
  double alp2[NALMAX];
  double alp3[NALMAX];
  double alp4[NALMAX];
  double bet1[NALMAX];
  double bet2[NALMAX];
  double fmat[NARMAX];
  char fnam[] = "TmxCalc";

  if(lsgm < EPSILON)
  {
    iflg = vb;
    ierr = 0;
    nkmax = -1;
    ndistr = 4;
    b = 1.0e-1;
    rmin = 0.9999999*rmod;
    rmax = 1.0000001*rmod;
  }
  else
  {
    iflg = vb;
    ierr = 0;
    nkmax = 5;
    ndistr = 2;
    b = log(10.0)*lsgm;
    b = b*b;
  }
  tmatrix_(&iflg,&ierr,&rat,&ndistr,&rmod,&npnax,&rmin,&rmax,&b,&gam,
           &nkmax,&arat,&np,&wlen,&refr,&refi,&ddelt,&nang,&ndgs,
           reff,veff,cext,csca,walb,asym,
           lmax,alp1,alp2,alp3,alp4,bet1,bet2,fmat);
  if(ierr)
  {
    fprintf(stderr,"%s: error, ierr=%d\n",fnam,ierr);
    return -1;
  }

  for(i=0; i<npnax; i++)
  {
    printf("Wavelength      : %13.4f\n",wlen);
    printf("Effective radius: %13.4f\n",reff[i]);
    printf("Effective var   : %13.4f\n",veff[i]);
    printf("Extinction coeff: %13.6e\n",cext[i]);
    printf("Scattering coeff: %13.6e\n",csca[i]);
    printf("SSA             : %13.6e\n",walb[i]);
    printf("Asymmetry param : %13.6e\n",asym[i]);
    printf("LMAX            : %13d\n",lmax[i]);
    for(j=0; j<lmax[i]; j++)
    {
      k = i*NPL+j;
      printf("%4d %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n",j,alp1[k],alp2[k],alp3[k],alp4[k],bet1[k],bet2[k]);
    }
    for(j=0; j<nang; j++)
    {
      printf("%7.2f",180.0/(nang-1)*j);
      for(k=0; k<6; k++)
      {
        printf(" %10.4f",fmat[i*NAFMAX+j*6+k]);
      }
      printf("\n");
    }
  }

  return 0;
}

int Init(void)
{
  return 0;
}

int GetOpt(int argn,char **args)
{
  int c,rt;
  int ntmp;
  int option_index = 0;
  int this_option_optind;
  char *endp;
  double xtmp;
  struct option long_options[] =
  {
    {"nang",1,0,'n'},
    {"refr",1,0,'R'},
    {"refi",1,0,'I'},
    {"wlen",1,0,'w'},
    {"rmin",1,0,'x'},
    {"rmax",1,0,'X'},
    {"rmod",1,0,'r'},
    {"lsgm",1,0,'s'},
    {"arat",1,0,'e'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":n:R:I:w:x:X:r:s:e:dvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'n':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0) nang = ntmp;
        else
        {
          fprintf(stderr,"#Angles -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'R':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') refr = xtmp;
        else
        {
          fprintf(stderr,"Refractive index (real) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'I':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') refi = xtmp;
        else
        {
          fprintf(stderr,"Refractive index (imag) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'w':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) wlen = xtmp;
        else
        {
          fprintf(stderr,"Wavelength -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'x':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) rmin = xtmp;
        else
        {
          fprintf(stderr,"R min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'X':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) rmax = xtmp;
        else
        {
          fprintf(stderr,"R max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'r':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) rmod = xtmp;
        else
        {
          fprintf(stderr,"R mode -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 's':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>=0.0) lsgm = xtmp;
        else
        {
          fprintf(stderr,"R sigma -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'e':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) arat = xtmp;
        else
        {
          fprintf(stderr,"Aspect ratio -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'd':
        db++;
        break;
      case 'v':
        vb++;
        break;
      case 'h':
        hp = 1;
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

  if(optind < argn)
  {
    fprintf(stderr,"non-option ARGV-elements:\n");
    while(optind < argn)
    {
      fprintf(stderr,"%s\n",args[optind++]);
    }
    rt = -1;
  }

  if(hp) rt = 1;
  return rt;
}

int Usage(void)
{
  int  n = 15;
  char e[MAXLINE];
  char a[MAXLINE];
  char d[MAXLINE];

  fprintf(stderr,"tmx_test ... calculate mie scattering parameters.\n");
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"tmx_test -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," n -nang    |%s|%s|%s| %d\n",As(e,"#Angles",n),       As(a,"#",n),          Ad(d,NANG,n),nang);
  fprintf(stderr," R -refr    |%s|%s|%s| %f\n",As(e,"Re(Ref. index)",n),As(a,"value",n),      Af(d,REFR,n),refr);
  fprintf(stderr," I -refi    |%s|%s|%s| %f\n",As(e,"Im(Ref. index)",n),As(a,"value",n),      Af(d,REFI,n),refi);
  fprintf(stderr," w -wlen    |%s|%s|%s| %f\n",As(e,"Wavelength",n),    As(a,"um",n),         Af(d,WLEN,n),wlen);
  fprintf(stderr," x -rmin    |%s|%s|%s| %e\n",As(e,"R min",n),         As(a,"um",n),         Af(d,RMIN,n),rmin);
  fprintf(stderr," X -rmax    |%s|%s|%s| %e\n",As(e,"R max",n),         As(a,"um",n),         Af(d,RMAX,n),rmax);
  fprintf(stderr," r -rmod    |%s|%s|%s| %e\n",As(e,"R mode",n),        As(a,"um",n),         Af(d,RMOD,n),rmod);
  fprintf(stderr," s -lsgm    |%s|%s|%s| %e\n",As(e,"R sigma",n),       As(a,"value",n),      Af(d,LSGM,n),lsgm);
  fprintf(stderr," e -arat    |%s|%s|%s| %f\n",As(e,"Aspect ratio",n),  As(a,"value",n),      Af(d,ARAT,n),arat);
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"-----------------------------------------------------------------------------\n");

  return 0;
}
