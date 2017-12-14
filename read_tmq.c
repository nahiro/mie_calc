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
#define	TMX_NANG		361				// #Angles
#define	TMX_NDGS		2				// Control #divisions
#define	TMX_SHAP		-1				// Shape#
#define	TMX_REPS		1.000001			// Aspect ratio or deformation parameter
#define	TMX_DELT		0.001				// Accuracy of the computations
#define	EPSILON			1.0e-14				// A small number
#define	DELTA			1.0e-11				// A small number
#define	MAXLINE			256				// Max #chars in a line
#define	NPN1			300
#define	NPN2			(2*NPN1)
#define	NPL			(NPN2+1)
#define	NAXMAX			10
#define	NAMAX			500
#define	NAFMAX			(NAMAX*6)
#define	NALMAX			(NAXMAX*NPL)
#define	NARMAX			(NAXMAX*NAFMAX)
#define	TMX_WLEN		0.550				// Wavelength in um
#define	TMX_LMIN		-3.0				// Log(Min size parameter)
#define	TMX_LMAX		2.0001				// Log(Max size parameter)
#define	TMX_LSTP		0.01				// Log(size parameter) step
#define	TMX_REFR		1.53				// Real part of refractive index
#define	TMX_REFI		0.008				// Imaginal part of refractive index

char	fnam_qext[MAXLINE]	= "qext.dat";
char	fnam_qsca[MAXLINE]	= "qsca.dat";
char	fnam_p1[MAXLINE]	= "p1.dat";
char	fnam_p2[MAXLINE]	= "p2.dat";
int	cnt_db			= 0;				// Debug   mode
int	cnt_vb			= 0;				// Verbose mode
int	cnt_hp			= 0;				// Help    mode
int	tmx_nang		= TMX_NANG;			// #Angles
int	tmx_ndgs		= TMX_NDGS;			// Control #divisions
int	tmx_shap		= TMX_SHAP;			// Shape#
double	tmx_reps		= TMX_REPS;
double	tmx_delt		= TMX_DELT;
double	tmx_wlen		= TMX_WLEN;
double	tmx_lmin		= TMX_LMIN;
double	tmx_lmax		= TMX_LMAX;
double	tmx_lstp		= TMX_LSTP;
double	tmx_refr		= TMX_REFR;
double	tmx_refi		= TMX_REFI;

int GetOpt(int argn,char **args);
int Usage(void);
int ReadTmq(void);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(cnt_hp) {Usage(); return 0;}
  if(ReadTmq() < 0) exit(-1);

  return 0;
}

int ReadTmq(void)
{
  int i,imax;
  int err;
  int ierr;
  int np,npna,ndgs;
  int lmax;
  int dummy = 0.0;
  double l;
  double eps;
  double lam,mrr,mri;
  double xval,tmx_xval;
  double ddelt;
  double reff;
  double veff;
  double cext;
  double csca;
  double walb;
  double gsym;
  double fact;
  double qext;
  double qsca;
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
  double p1[NAMAX];
  double p2[NAMAX];
  char fnam[] = "ReadTmq";
  FILE *fp_qext;
  FILE *fp_qsca;
  FILE *fp_p1;
  FILE *fp_p2;

  // open file
  fp_qext = fopen(fnam_qext,"w");
  fp_qsca = fopen(fnam_qsca,"w");
  fp_p1 = fopen(fnam_p1,"w");
  fp_p2 = fopen(fnam_p2,"w");
  if(fp_qext==NULL || fp_qsca==NULL || fp_p1==NULL || fp_p2==NULL)
  {
    fprintf(stderr,"%s: failed in opening file.\n",fnam);
    return -1;
  }

  // read T-matrix
  for(l=tmx_lmin; l<tmx_lmax; l+=tmx_lstp)
  {
    tmx_xval = pow(10.0,l);
    fread(&lam,sizeof(double),1,stdin);
    fread(&xval,sizeof(double),1,stdin);
    fread(&mrr,sizeof(double),1,stdin);
    fread(&mri,sizeof(double),1,stdin);
    fread(&eps,sizeof(double),1,stdin);
    fread(&ddelt,sizeof(double),1,stdin);
    fread(&ndgs,sizeof(int),1,stdin);
    fread(&np,sizeof(int),1,stdin);
    fread(&ierr,sizeof(int),1,stdin);
    fread(&dummy,sizeof(int),1,stdin);
    fread(&reff,sizeof(double),1,stdin);
    fread(&veff,sizeof(double),1,stdin);
    fread(&cext,sizeof(double),1,stdin);
    fread(&csca,sizeof(double),1,stdin);
    fread(&walb,sizeof(double),1,stdin);
    fread(&gsym,sizeof(double),1,stdin);
    fread(&lmax,sizeof(int),1,stdin);
    fread(alp1,sizeof(double),lmax,stdin);
    fread(alp2,sizeof(double),lmax,stdin);
    fread(alp3,sizeof(double),lmax,stdin);
    fread(alp4,sizeof(double),lmax,stdin);
    fread(bet1,sizeof(double),lmax,stdin);
    fread(bet2,sizeof(double),lmax,stdin);
    fread(&npna,sizeof(int),1,stdin);
    fread(f11,sizeof(double),npna,stdin);
    fread(f22,sizeof(double),npna,stdin);
    fread(f33,sizeof(double),npna,stdin);
    fread(f44,sizeof(double),npna,stdin);
    fread(f12,sizeof(double),npna,stdin);
    fread(f34,sizeof(double),npna,stdin);
    err = 0;
    if(fabs(1.0-lam/tmx_wlen) > DELTA)
    {
      fprintf(stderr,"Error, lam=%13.6e, tmx_wlen=%13.6e (x=%13.6e)\n",lam,tmx_wlen,xval);
      err = 1;
    }
    if(fabs(1.0-xval/tmx_xval) > DELTA)
    {
      fprintf(stderr,"Warning, xval=%13.6e, tmx_xval=%13.6e (x=%13.6e)\n",xval,tmx_xval,xval);
    }
    if(fabs(1.0-mrr/tmx_refr) > DELTA)
    {
      fprintf(stderr,"Error, mrr=%13.6e, tmx_refr=%13.6e (x=%13.6e)\n",mrr,tmx_refr,xval);
      err = 1;
    }
    if(fabs(1.0-mri/tmx_refi) > DELTA)
    {
      fprintf(stderr,"Error, mri=%13.6e, tmx_refi=%13.6e (x=%13.6e)\n",mri,tmx_refi,xval);
      err = 1;
    }
    if(fabs(1.0-eps/tmx_reps) > DELTA)
    {
      fprintf(stderr,"Error, eps=%13.6e, tmx_reps=%13.6e (x=%13.6e)\n",eps,tmx_reps,xval);
      err = 1;
    }
    if(fabs(1.0-ddelt/tmx_delt) > DELTA)
    {
      fprintf(stderr,"Error, ddelt=%13.6e, tmx_delt=%13.6e (x=%13.6e)\n",ddelt,tmx_delt,xval);
      err = 1;
    }
    if(ndgs != tmx_ndgs)
    {
      fprintf(stderr,"Error, ndgs=%d, tmx_ndgs=%d (x=%13.6e)\n",ndgs,tmx_ndgs,xval);
      err = 1;
    }
    if(np != tmx_shap)
    {
      fprintf(stderr,"Error, np=%d, tmx_shap=%d (x=%13.6e)\n",np,tmx_shap,xval);
      err = 1;
    }
    if(npna != tmx_nang)
    {
      fprintf(stderr,"Error, npna=%d, tmx_nang=%d (x=%13.6e)\n",npna,tmx_nang,xval);
      err = 1;
    }
    if(err != 0)
    {
      return -1;
    }
    if(cnt_vb)
    {
      printf("WLEN: %13.6e\n",tmx_wlen);
      printf("XVAL: %13.6e\n",tmx_xval);
      printf("REFR: %13.6e\n",tmx_refr);
      printf("REFI: %13.6e\n",tmx_refi);
      printf("REPS: %13.6e\n",tmx_reps);
      printf("DELT: %13.6e\n",tmx_delt);
      printf("NDGS: %13d\n",tmx_ndgs);
      printf("SHAP: %13d\n",tmx_shap);
      printf("IERR: %13d\n",ierr);
      printf("REFF: %13.6e\n",reff);
      printf("VEFF: %13.6e\n",veff);
      printf("CEXT: %13.6e\n",cext);
      printf("CSCA: %13.6e\n",csca);
      printf("WALB: %13.6e\n",walb);
      printf("GSYM: %13.6e\n",gsym);
      printf("LMAX: %13d\n",lmax);
      imax = lmax-1;
      printf("ALP1:\n");
      for(i=0; i<lmax; i++)
      {
        printf("%13.6e%s",alp1[i],(i%8==7?"\n":(i==imax?"\n":" ")));
      }
      printf("ALP2:\n");
      for(i=0; i<lmax; i++)
      {
        printf("%13.6e%s",alp2[i],(i%8==7?"\n":(i==imax?"\n":" ")));
      }
      printf("ALP3:\n");
      for(i=0; i<lmax; i++)
      {
        printf("%13.6e%s",alp3[i],(i%8==7?"\n":(i==imax?"\n":" ")));
      }
      printf("ALP4:\n");
      for(i=0; i<lmax; i++)
      {
        printf("%13.6e%s",alp4[i],(i%8==7?"\n":(i==imax?"\n":" ")));
      }
      printf("BET1:\n");
      for(i=0; i<lmax; i++)
      {
        printf("%13.6e%s",bet1[i],(i%8==7?"\n":(i==imax?"\n":" ")));
      }
      printf("BET2:\n");
      for(i=0; i<lmax; i++)
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
    fact = 1.0/(tmx_xval*tmx_wlen);
    fact *= fact;
    fact *= 4.0*M_PI;
    qext = cext*fact;
    qsca = csca*fact;
    for(i=0; i<tmx_nang; i++)
    {
      p1[i] = (f11[i]-f12[i])*csca;
      p2[i] = (f11[i]+f12[i])*csca;
    }
    fwrite(&qext,sizeof(double),1,fp_qext);
    fwrite(&qsca,sizeof(double),1,fp_qsca);
    fwrite(p1,sizeof(double),tmx_nang,fp_p1);
    fwrite(p2,sizeof(double),tmx_nang,fp_p2);
  }

  // cleanup
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
    {"tmx_wlen",1,0,'w'},
    {"tmx_lmin",1,0,'l'},
    {"tmx_lmax",1,0,'L'},
    {"tmx_lstp",1,0,'s'},
    {"tmx_refr",1,0,'R'},
    {"tmx_refi",1,0,'I'},
    {"tmx_ndgs",1,0,'c'},
    {"tmx_nang",1,0,'M'},
    {"tmx_shap",1,0,'S'},
    {"tmx_reps",1,0,'e'},
    {"tmx_delt",1,0,'D'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"w:l:L:s:R:I:c:M:S:e:D:dvh",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'w':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) tmx_wlen = xtmp;
        else
        {
          fprintf(stderr,"Wavelength -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'l':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') tmx_lmin = xtmp;
        else
        {
          fprintf(stderr,"Log(X min) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'L':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') tmx_lmax = xtmp;
        else
        {
          fprintf(stderr,"Log(X max) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 's':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') tmx_lstp = xtmp;
        else
        {
          fprintf(stderr,"Log(X) step -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'R':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') tmx_refr = xtmp;
        else
        {
          fprintf(stderr,"Refractive index (real) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'I':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') tmx_refi = xtmp;
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
  fprintf(stderr,"read_tmq -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"------------------------------------------------------------------------------------\n");
  fprintf(stderr,"   option         |%s|%s|%s| current\n",As(e,"",n),          As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," w -tmx_wlen      |%s|%s|%s| %f\n",As(e,"Wavelength",n),     As(a,"um",n),         Af(d,TMX_WLEN,n),tmx_wlen);
  fprintf(stderr," l -tmx_lmin      |%s|%s|%s| %f\n",As(e,"Log(X min)",n),     As(a,"value",n),      Af(d,TMX_LMIN,n),tmx_lmin);
  fprintf(stderr," L -tmx_lmax      |%s|%s|%s| %f\n",As(e,"Log(X max)",n),     As(a,"value",n),      Af(d,TMX_LMAX,n),tmx_lmax);
  fprintf(stderr," s -tmx_lstp      |%s|%s|%s| %f\n",As(e,"Log(X) step",n),    As(a,"value",n),      Af(d,TMX_LSTP,n),tmx_lstp);
  fprintf(stderr," R -tmx_refr      |%s|%s|%s| %f\n",As(e,"Re(Ref. index)",n), As(a,"value",n),      Af(d,TMX_REFR,n),tmx_refr);
  fprintf(stderr," I -tmx_refi      |%s|%s|%s| %e\n",As(e,"Im(Ref. index)",n), As(a,"value",n),      Ae(d,TMX_REFI,n),tmx_refi);
  fprintf(stderr," c -tmx_ndgs      |%s|%s|%s| %d\n",As(e,"Control#",n),       As(a,"#",n),          Ad(d,TMX_NDGS,n),tmx_ndgs);
  fprintf(stderr," M -tmx_nang      |%s|%s|%s| %d\n",As(e,"#Angles",n),        As(a,"#",n),          Ad(d,TMX_NANG,n),tmx_nang);
  fprintf(stderr," S -tmx_shap      |%s|%s|%s| %d\n",As(e,"Shape#",n),         As(a,"#",n),          Ad(d,TMX_SHAP,n),tmx_shap);
  fprintf(stderr," e -tmx_reps      |%s|%s|%s| %e\n",As(e,"Aspect ratio",n),   As(a,"value",n),      Af(d,TMX_REPS,n),tmx_reps);
  fprintf(stderr," D -tmx_delt      |%s|%s|%s| %e\n",As(e,"Accuracy",n),       As(a,"value",n),      Ae(d,TMX_DELT,n),tmx_delt);
  fprintf(stderr," d -debug         |%s|%s|%s| %d\n",As(e,"Debug   mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_db);
  fprintf(stderr," v -verbose       |%s|%s|%s| %d\n",As(e,"Verbose mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_vb);
  fprintf(stderr," h -help          |%s|%s|%s| %d\n",As(e,"Help    mode",n),   As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"------------------------------------------------------------------------------------\n");

  return 0;
}
