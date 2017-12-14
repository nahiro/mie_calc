/***********************************************************/
/* TMX_CALC ... Aerosol generator with T-Matrix            */
/* Author: N.Manago                                        */
/* $Revision: 790 $                                        */
/* $Date: 2014-09-01 04:11:02 +0900 (Mon, 01 Sep 2014) $   */
/***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <bits/nan.h>
#include "strutil.h"

#define	PI			3.141592653589793		// PI
#define	PI2			6.283185307179586		// 2*PI
#define	PI2_1			0.15915494309189535		// 1/(2*PI)
#define	TMX_NANG		361				// #Angles
#define	TMX_NDGS		2				// Control #divisions
#define	TMX_SHAP		-1				// Shape#
#define	TMX_REPS		1.000001			// Aspect ratio or deformation parameter
#define	TMX_DELT		0.001				// Accuracy of the computations
#define	EPSILON			1.0e-14				// A small number
#define	DELTA			1.0e-13				// A small number
#define	MAXLINE			256				// Max #chars in a line
#define	NPN1			300
#define	NDGSMAX			32
#define	NPNG1			((NDGSMAX+1)*NPN1)
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
#define	WLEN			0.550				// Wavelength in um
#define	XVAL			1.0				// Size parameter
#define	REFR			1.53				// Real part of refractive index
#define	REFI			0.008				// Imaginal part of refractive index

int	cnt_bin			= 0;				// Binary  mode
int	cnt_db			= 0;				// Debug   mode
int	cnt_vb			= 0;				// Verbose mode
int	cnt_hp			= 0;				// Help    mode
int	tmx_nang		= TMX_NANG;			// #Angles
int	tmx_ndgs		= TMX_NDGS;			// Control #divisions
int	tmx_shap		= TMX_SHAP;			// Shape#
double	tmx_reps		= TMX_REPS;
double	tmx_delt		= TMX_DELT;
double	wlen			= WLEN;
double	xval			= XVAL;
double	refr			= REFR;
double	refi			= REFI;

int MieCalc(void);
int GetOpt(int argn,char **args);
int Usage(void);
extern void tmatrix_();

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(cnt_hp) {Usage(); return 0;}
  if(MieCalc() < 0) exit(-1);

  return 0;
}

int MieCalc(void)
{
  int i,imax;
  int iflg,ierr;
  int ndistr,npnax,nkmax;
  int np,npna,ndgs;
  int lmax[NAXMAX];
  int dummy = 0.0;
  double rat,axmax,r,r1,r2;
  double b,gam,eps;
  double lam,mrr,mri;
  double ddelt;
  double reff[NAXMAX];
  double veff[NAXMAX];
  double cext[NAXMAX];
  double csca[NAXMAX];
  double walb[NAXMAX];
  double gsym[NAXMAX];
  double alp1[NALMAX];
  double alp2[NALMAX];
  double alp3[NALMAX];
  double alp4[NALMAX];
  double bet1[NALMAX];
  double bet2[NALMAX];
  double fmat[NARMAX];
  double f11[NAMAX];
  double f22[NAMAX];
  double f33[NAMAX];
  double f44[NAMAX];
  double f12[NAMAX];
  double f34[NAMAX];
  char fnam[] = "MieCalc";

  // mie calculation
  r = xval*wlen*PI2_1;
  iflg   = cnt_vb;
  ierr   = 0;
  rat    = 0.5; // equal-surface-area-sphere radius
  ndistr = 4;
  axmax  = r;
  npnax  = 1;
  r1     = 0.9999999*r;
  r2     = 1.0000001*r;
  b      = 1.0e-1;
  gam    = 0.5; // ignored
  nkmax  = -1;
  eps    = tmx_reps;
  np     = tmx_shap;
  lam    = wlen;
  mrr    = refr;
  mri    = refi;
  ddelt  = tmx_delt; // this will be modified by tmatrix_()!
  npna   = tmx_nang;
  ndgs   = tmx_ndgs;
  if(cnt_db)
  {
    fprintf(stderr,"%6s:%13d %6s:%13d %6s:%13.4f %6s:%13d %6s:%13.6e %6s:%13d\n",
                    "iflg",iflg,"ierr",ierr,"rat",rat,"ndistr",ndistr,"axmax",axmax,"npnax",npnax);
    fprintf(stderr,"%6s:%13.6e %6s:%13.6e %6s:%13.6e %6s:%13.6f %6s:%13d %6s:%13.6e %6s:%13d\n",
                    "r1",r1,"r2",r2,"b",b,"gam",gam,"nkmax",nkmax,"eps",eps,"np",np);
    fprintf(stderr,"%6s:%13.6f %6s:%13.6e %6s:%13.6e %6s:%13.6e %6s:%13d %6s:%13d\n",
                    "lam",lam,"mrr",mrr,"mri",mri,"ddelt",ddelt,"npna",npna,"ndgs",ndgs);
  }
  tmatrix_(&iflg,&ierr,&rat,&ndistr,&axmax,&npnax,&r1,&r2,&b,&gam,
           &nkmax,&eps,&np,&lam,&mrr,&mri,&ddelt,&npna,&ndgs,
           reff,veff,cext,csca,walb,gsym,
           lmax,alp1,alp2,alp3,alp4,bet1,bet2,fmat);
  if(ierr)
  {
    fprintf(stderr,"%s: warning, ierr=%d, r=%13.6e, lam=%13.6f, x=%13.6e, eps=%13.6e\n",
                    fnam,ierr,r,lam,xval,eps);
  }
  for(i=0; i<tmx_nang; i++)
  {
    f11[i] = fmat[i*6+0];
    f22[i] = fmat[i*6+1];
    f33[i] = fmat[i*6+2];
    f44[i] = fmat[i*6+3];
    f12[i] = fmat[i*6+4];
    f34[i] = fmat[i*6+5];
  }

  if(cnt_bin)
  {
    fwrite(&wlen,sizeof(double),1,stdout);
    fwrite(&xval,sizeof(double),1,stdout);
    fwrite(&refr,sizeof(double),1,stdout);
    fwrite(&refi,sizeof(double),1,stdout);
    fwrite(&tmx_reps,sizeof(double),1,stdout);
    fwrite(&tmx_delt,sizeof(double),1,stdout);
    fwrite(&tmx_ndgs,sizeof(int),1,stdout);
    fwrite(&tmx_shap,sizeof(int),1,stdout);
    fwrite(&ierr,sizeof(int),1,stdout);
    fwrite(&dummy,sizeof(int),1,stdout);
    fwrite(reff,sizeof(double),1,stdout);
    fwrite(veff,sizeof(double),1,stdout);
    fwrite(cext,sizeof(double),1,stdout);
    fwrite(csca,sizeof(double),1,stdout);
    fwrite(walb,sizeof(double),1,stdout);
    fwrite(gsym,sizeof(double),1,stdout);
    fwrite(lmax,sizeof(int),1,stdout);
    fwrite(alp1,sizeof(double),lmax[0],stdout);
    fwrite(alp2,sizeof(double),lmax[0],stdout);
    fwrite(alp3,sizeof(double),lmax[0],stdout);
    fwrite(alp4,sizeof(double),lmax[0],stdout);
    fwrite(bet1,sizeof(double),lmax[0],stdout);
    fwrite(bet2,sizeof(double),lmax[0],stdout);
    fwrite(&tmx_nang,sizeof(int),1,stdout);
    fwrite(f11,sizeof(double),tmx_nang,stdout);
    fwrite(f22,sizeof(double),tmx_nang,stdout);
    fwrite(f33,sizeof(double),tmx_nang,stdout);
    fwrite(f44,sizeof(double),tmx_nang,stdout);
    fwrite(f12,sizeof(double),tmx_nang,stdout);
    fwrite(f34,sizeof(double),tmx_nang,stdout);
  }
  else
  {
    printf("WLEN: %13.6e\n",wlen);
    printf("XVAL: %13.6e\n",xval);
    printf("REFR: %13.6e\n",refr);
    printf("REFI: %13.6e\n",refi);
    printf("REPS: %13.6e\n",tmx_reps);
    printf("DELT: %13.6e\n",tmx_delt);
    printf("NDGS: %13d\n",tmx_ndgs);
    printf("SHAP: %13d\n",tmx_shap);
    printf("IERR: %13d\n",ierr);
    printf("REFF: %13.6e\n",reff[0]);
    printf("VEFF: %13.6e\n",veff[0]);
    printf("CEXT: %13.6e\n",cext[0]);
    printf("CSCA: %13.6e\n",csca[0]);
    printf("WALB: %13.6e\n",walb[0]);
    printf("GSYM: %13.6e\n",gsym[0]);
    printf("LMAX: %13d\n",lmax[0]);
    imax = lmax[0]-1;
    printf("ALP1:\n");
    for(i=0; i<lmax[0]; i++)
    {
      printf("%13.6e%s",alp1[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("ALP2:\n");
    for(i=0; i<lmax[0]; i++)
    {
      printf("%13.6e%s",alp2[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("ALP3:\n");
    for(i=0; i<lmax[0]; i++)
    {
      printf("%13.6e%s",alp3[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("ALP4:\n");
    for(i=0; i<lmax[0]; i++)
    {
      printf("%13.6e%s",alp4[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("BET1:\n");
    for(i=0; i<lmax[0]; i++)
    {
      printf("%13.6e%s",bet1[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("BET2:\n");
    for(i=0; i<lmax[0]; i++)
    {
      printf("%13.6e%s",bet2[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("NANG: %13d\n",tmx_nang);
    imax = tmx_nang-1;
    printf("F11:\n");
    for(i=0; i<tmx_nang; i++)
    {
      printf("%13.6e%s",f11[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("F22:\n");
    for(i=0; i<tmx_nang; i++)
    {
      printf("%13.6e%s",f22[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("F33:\n");
    for(i=0; i<tmx_nang; i++)
    {
      printf("%13.6e%s",f33[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("F44:\n");
    for(i=0; i<tmx_nang; i++)
    {
      printf("%13.6e%s",f44[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("F12:\n");
    for(i=0; i<tmx_nang; i++)
    {
      printf("%13.6e%s",f12[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
    printf("F34:\n");
    for(i=0; i<tmx_nang; i++)
    {
      printf("%13.6e%s",f34[i],(i%8==7?"\n":(i==imax?"\n":" ")));
    }
  }

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
    {"wlen",1,0,'w'},
    {"xval",1,0,'x'},
    {"refr",1,0,'R'},
    {"refi",1,0,'I'},
    {"tmx_ndgs",1,0,'c'},
    {"tmx_nang",1,0,'M'},
    {"tmx_shap",1,0,'S'},
    {"tmx_reps",1,0,'e'},
    {"tmx_delt",1,0,'D'},
    {"binary",0,0,'b'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"w:x:R:I:c:M:S:e:D:bdvh",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
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
        if(errno!=ERANGE && *endp=='\0') xval = xtmp;
        else
        {
          fprintf(stderr,"X value -> out of range %s\n",optarg);
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
      case 'c':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) tmx_ndgs = ntmp;
        else
        {
          fprintf(stderr,"Control# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'M':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>1) tmx_nang = ntmp;
        else
        {
          fprintf(stderr,"#Angles -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'S':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=-2) tmx_shap = ntmp;
        else
        {
          fprintf(stderr,"Shape# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'e':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0)
        {
          tmx_reps = xtmp;
          if(fabs(tmx_reps-1.0) < 1.0e-6)
          {
            tmx_reps = 1.000001;
          }
        }
        else
        {
          fprintf(stderr,"Aspect ratio -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'D':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) tmx_delt = xtmp;
        else
        {
          fprintf(stderr,"Accuracy -> out of range %s\n",optarg);
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
  fprintf(stderr,"tmq_calc -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"------------------------------------------------------------------------------------\n");
  fprintf(stderr,"   option         |%s|%s|%s| current\n",As(e,"",n),          As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," w -wlen          |%s|%s|%s| %f\n",As(e,"Wavelength",n),     As(a,"um",n),         Af(d,WLEN,n),wlen);
  fprintf(stderr," x -xval          |%s|%s|%s| %f\n",As(e,"X Value",n),        As(a,"value",n),      Af(d,XVAL,n),xval);
  fprintf(stderr," R -refr          |%s|%s|%s| %f\n",As(e,"Re(Ref. index)",n), As(a,"value",n),      Af(d,REFR,n),refr);
  fprintf(stderr," I -refi          |%s|%s|%s| %e\n",As(e,"Im(Ref. index)",n), As(a,"value",n),      Ae(d,REFI,n),refi);
  fprintf(stderr," c -tmx_ndgs      |%s|%s|%s| %d\n",As(e,"Control#",n),       As(a,"#",n),          Ad(d,TMX_NDGS,n),tmx_ndgs);
  fprintf(stderr," M -tmx_nang      |%s|%s|%s| %d\n",As(e,"#Angles",n),        As(a,"#",n),          Ad(d,TMX_NANG,n),tmx_nang);
  fprintf(stderr," S -tmx_shap      |%s|%s|%s| %d\n",As(e,"Shape#",n),         As(a,"#",n),          Ad(d,TMX_SHAP,n),tmx_shap);
  fprintf(stderr," e -tmx_reps      |%s|%s|%s| %e\n",As(e,"Aspect ratio",n),   As(a,"value",n),      Af(d,TMX_REPS,n),tmx_reps);
  fprintf(stderr," D -tmx_delt      |%s|%s|%s| %e\n",As(e,"Accuracy",n),       As(a,"value",n),      Ae(d,TMX_DELT,n),tmx_delt);
  fprintf(stderr," b -binary        |%s|%s|%s| %d\n",As(e,"Binary  mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_bin);
  fprintf(stderr," d -debug         |%s|%s|%s| %d\n",As(e,"Debug   mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_db);
  fprintf(stderr," v -verbose       |%s|%s|%s| %d\n",As(e,"Verbose mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_vb);
  fprintf(stderr," h -help          |%s|%s|%s| %d\n",As(e,"Help    mode",n),   As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"------------------------------------------------------------------------------------\n");

  return 0;
}
