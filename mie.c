/***********************************************************/
/* MIE ... calculate mie scattering parameters             */
/* Author: N.Manago                                        */
/* $Revision: 1.1.1.1 $                                       */
/* $Date: 2008-04-03 04:42:18 $ (UTC)                      */
/***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <complex.h>
#include "strutil.h"

#define	PI				3.141592653589793		// PI
#define	PI2				6.283185307179586		// 2*PI
#define	PI_2				1.570796326794896		// PI/2
#define	D_TO_R				1.745329251994329e-02		// PI/180 deg -> rad
#define	R_TO_D				57.29577951308232		// 180/PI rad -> deg
#define	TRUE				1
#define	FALSE				0

#define	NVAL				0				// Order
#define	REFR				1.53				// Real part of refractive index
#define	REFI				1.0e-9				// Imaginal part of refractive index
#define	LMIN				-2.0				// Log10(X) Min
#define	LMAX				3.0				// Log10(X) Max
#define	LSTP				1.0e-3				// Log10(X) step
#define	AMIN				0.0				// A Min
#define	AMAX				180.0				// A Max
#define	ASTP				1.0				// A step
#define	MAXNANG				2000				// Max #angles
#define	MAXLINE				256				// Max #chars in a line

int		nang			= 0;
int		vb			= 0;				// Verbose mode
int		db			= 0;				// Debug   mode
int		hp			= 0;				// Help    mode
double		lmin			= LMIN;
double		lmax			= LMAX;
double		lstp			= LSTP;
double		amin			= AMIN;
double		amax			= AMAX;
double		astp			= ASTP;
double		refr			= REFR;
double		refi			= REFI;

extern void miev0_();

int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);

int main(int argc,char **argv)
{
  int i;
  double l,x,a;
  double s1,s2;
  int ANYANG,PERFCT,PRNT[2];
  int IPOLZN,MOMDIM,NUMANG,NMOM;
  double GQSC,MIMCUT,PMOM[2][2],QEXT,QSCA,SPIKE,XMU[MAXNANG];
  double _Complex CREFIN,SFORW,SBACK,S1[MAXNANG],S2[MAXNANG],TFORW[2],TBACK[2];

  if(GetOpt(argc,argv) < 0) exit(-1);
  if(hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  CREFIN = refr-fabs(refi)*I;
  PERFCT = FALSE;
  MIMCUT = 1.0e-7;
  ANYANG = FALSE;
  NUMANG = nang;
  PRNT[0] = (db>0?TRUE:FALSE);
  PRNT[1] = (vb>0?TRUE:FALSE);
  NMOM = 0;
  IPOLZN = 0;
  MOMDIM = 1;

  for(l=lmin; l<=lmax; l+=lstp)
  {
    x = pow(10.0,l);
    for(i=0,a=amin; i<nang; i++,a+=astp)
    {
      XMU[i] = cos(a*D_TO_R);
    }
    miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,XMU,
           &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
           PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,
           &SPIKE);
    if(nang == 1)
    {
      printf("%13.6e %13.6e %13.6e\n",x,QEXT,QSCA);
    }
    else
    {
      for(i=0,a=amin; i<nang; i++,a+=astp)
      {
        s1 = cabs(S1[i]);
        s2 = cabs(S2[i]);
        printf("%13.6e %13.6f %13.6e\n",x,a,(s1*s1+s2*s2)/(PI2*x*x*QSCA));
      }
    }
  }

  return 0;
}

int Init(void)
{
  int i = 0;
  double a;

  for(a=amin; a<=amax; a+=astp)
  {
    i++;
  }
  nang = i;
  if(nang > MAXNANG)
  {
    fprintf(stderr,"Error, nang=%d >%d\n",nang,MAXNANG);
    return -1;
  }

  return 0;
}

int GetOpt(int argn,char **args)
{
  int c,rt;
  int option_index = 0;
  int this_option_optind;
  char *endp;
  double xtmp;
  struct option long_options[] =
  {
    {"refr",1,0,'R'},
    {"refi",1,0,'I'},
    {"lmin",1,0,'l'},
    {"lmax",1,0,'L'},
    {"lstp",1,0,'p'},
    {"aval",1,0,'a'},
    {"amax",1,0,'A'},
    {"astp",1,0,'q'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":n:R:I:l:L:p:a:A:q:dvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'R':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') refr = xtmp;
        else
        {
          fprintf(stderr,"Real part of refractive index -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'I':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') refi = xtmp;
        else
        {
          fprintf(stderr,"Imaginaly part of refractive index -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'l':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') lmin = xtmp;
        else
        {
          fprintf(stderr,"Log10(X) min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'L':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') lmax = xtmp;
        else
        {
          fprintf(stderr,"Log10(X) max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'p':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) lstp = xtmp;
        else
        {
          fprintf(stderr,"Log10(X) step -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'a':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp<90.0)
        {
          amin = xtmp;
          amax = xtmp;
        }
        else
        {
          fprintf(stderr,"A value -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'A':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') amax = xtmp;
        else
        {
          fprintf(stderr,"A max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'q':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) astp = xtmp;
        else
        {
          fprintf(stderr,"A step -> out of range %s\n",optarg);
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

  fprintf(stderr,"mie ... calculate mie scattering parameters.\n");
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"mie -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," R -refr    |%s|%s|%s| %f\n",As(e,"Re(Ref. index)",n),As(a,"value",n),      Af(d,REFR,n),refr);
  fprintf(stderr," I -refi    |%s|%s|%s| %e\n",As(e,"Im(Ref. index)",n),As(a,"value",n),      Ae(d,REFI,n),refi);
  fprintf(stderr," l -lmin    |%s|%s|%s| %f\n",As(e,"Log10(X) min",n),  As(a,"value",n),      Af(d,LMIN,n),lmin);
  fprintf(stderr," L -lmax    |%s|%s|%s| %f\n",As(e,"Log10(X) max",n),  As(a,"value",n),      Af(d,LMAX,n),lmax);
  fprintf(stderr," p -lstp    |%s|%s|%s| %f\n",As(e,"Log10(X) step",n), As(a,"value",n),      Af(d,LSTP,n),lstp);
  fprintf(stderr," a -aval    |%s|%s|%s| %f\n",As(e,"A Value",n),       As(a,"value",n),      Af(d,AMIN,n),amin);
  fprintf(stderr," A -amax    |%s|%s|%s| %f\n",As(e,"A Max",n),         As(a,"value",n),      Af(d,AMAX,n),amax);
  fprintf(stderr," q -astp    |%s|%s|%s| %f\n",As(e,"A Step",n),        As(a,"value",n),      Af(d,ASTP,n),astp);
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"-----------------------------------------------------------------------------\n");

  return 0;
}
