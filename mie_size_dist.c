/***********************************************************/
/* MIE_SIZE_DIST ... calculate mie scattering parameters   */
/* Author: N.Manago                                        */
/* $Revision: 1124 $                                        */
/* $Date: 2016-01-12 02:40:17 +0900 (Tue, 12 Jan 2016) $   */
/***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <bits/nan.h>
#include "strutil.h"
#include <complex.h>

#define	NAMAX			500
#define	NDIS			0				// Distribution#
#define	NANG			361				// #Angles
#define	REFR			1.53				// Refractive index (real)
#define	REFI			0.008				// Refractive index (imag)
#define	WLEN			0.550				// Wavelength in um
#define	RMIN			NAN				// R min in um
#define	RMAX			NAN				// R max in um
#define	LMIN			-3				// Log10(R min in um)
#define	LMAX			+2				// Log10(R max in um)
#define	NSTP			1000				// Log10(R in um) #steps
#define	LSTP			NAN				// Log10(R in um) step
#define	RMOD			0.5				// R mode  in um
#define	LSGM			0.0				// R sigma in um (log)
#define	WSGM			5.0				// R width in sigma
#define	XMAX			2.0e4				// Max size parameter
#define	TRUE			1				// True  value
#define	FALSE			0				// False value
#define	MAXLINE			256				// Max #chars in a line
#define	DELTA			1.0e-13				// A small number
#define	EPSILON			1.0e-14				// A small number
#define	PI			3.141592653589793		// PI
#define	PI2			6.283185307179586		// 2*PI
#define	M_LNA			2.3025850929940459		// Ln(10)
#define	R_TO_D			57.29577951308232		// 180/PI rad -> deg
#define	D_TO_R			1.745329251994329e-02		// PI/180 deg -> rad
#define	OMOD_OVERWR		1				// Overwrite mode
#define	OMOD_APPEND		2				// Append mode
#define	OUT1			"mie_out1.dat"			// Output file 1
#define	OUT2			"mie_out2.dat"			// Output file 2

extern void miev0_();

int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);
int MieCalc(void);
int Printout(void);

int	ndis			= NDIS;
int	nang			= NANG;
int	nstp			= NSTP;
int	omod			= OMOD_OVERWR;
int     gmod                    = 0;
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
double	wsgm			= WSGM;
double	lstp			= LSTP;
double	lmod,lmin,lmax;
double	aext,asca,omeg,asym;
double	angl[NAMAX];
double	phs1[NAMAX];
double	phs2[NAMAX];
double	phas[NAMAX];
char	out1[MAXLINE]		= OUT1;
char	out2[MAXLINE]		= OUT2;

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) return -1;
  if(hp) {Usage(); return 0;}
  if(Init() < 0) return -1;

  if(MieCalc() < 0) return -1;
  if(Printout() < 0) return -1;

  return 0;
}

int MieCalc(void)
{
  int j,m;
  double x,y,w,r,l;
  double sy,p1,p2,p3;
  double s1,s2,s3,s4;
  double norm,fact,indx;
  double epsilon = 1.0e-1;
  double rang[NAMAX];
  double sang[NAMAX];
  double cang[NAMAX];
  double dang[NAMAX];
  int ANYANG,PERFCT,PRNT[2];
  int IPOLZN,MOMDIM,NUMANG,NMOM;
  double GQSC,MIMCUT,PMOM[2][2],QEXT,QSCA,SPIKE;
  double _Complex CREFIN,SFORW,SBACK,S1[NAMAX],S2[NAMAX],TFORW[2],TBACK[2];

  // initialization
  aext = 0.0;
  asca = 0.0;
  omeg = 0.0;
  asym = 0.0;
  for(j=0; j<nang; j++)
  {
    phs1[j] = 0.0;
    phs2[j] = 0.0;
    phas[j] = 0.0;
    angl[j] = (180.0/(nang-1)*j);
    rang[j] = angl[j]*D_TO_R;
    sang[j] = sin(rang[j]);
    cang[j] = cos(rang[j]);
  }
  for(j=1; j<nang; j++)
  {
    dang[j] = PI*fabs(rang[j]-rang[j-1]);
  }
  PERFCT = FALSE;
  MIMCUT = DELTA;
  ANYANG = FALSE;
  NUMANG = nang;
  PRNT[0] = (db>1?TRUE:FALSE);
  PRNT[1] = (db>0?TRUE:FALSE);
  NMOM = 0;
  IPOLZN = 0;
  MOMDIM = 1;
  CREFIN = refr-fabs(refi)*I;

  // Mie calculation
  if(lsgm < EPSILON)
  {
    x = PI2*rmod/wlen;
    if(x >= XMAX)
    {
      fprintf(stderr,"MieCalc: error, x=%13.6e\n",x);
      return -1;
    }
    miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,cang,
           &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
           PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);
    aext = QEXT*PI*rmod*rmod;
    asca = QSCA*PI*rmod*rmod;
    omeg = QSCA/QEXT;
    for(j=0; j<nang; j++)
    {
      s1 = cabs(S1[j]);
      s2 = cabs(S2[j]);
      phs1[j] = s1*s1;
      phs2[j] = s2*s2;
      phas[j] = 0.5*(phs1[j]+phs2[j]);
    }
  }
  else
  {
    // integrate over size distribution
    if(isnan(rmin)) lmin = lmod-lsgm*wsgm;
    else            lmin = log10(rmin);
    if(lmin < LMIN) lmin = LMIN;
    if(isnan(rmax)) lmax = lmod+lsgm*wsgm;
    else            lmax = log10(rmax);
    if(lmax > LMAX) lmax = LMAX;
    if(isnan(lstp)) lstp = (lmax-lmin)/nstp;
    else            nstp = (int)((lmax-lmin)/lstp+0.5);
    if(nstp < 1)
    {
      fprintf(stderr,"MieCalc: error, invalid #steps (lmin=%13.4e lmax=%13.4e lstp=%13.4e nstp=%d\n",
                      lmin,lmax,lstp,nstp);
      return -1;
    }
    if(gmod)
    {
      indx = 1.0/lsgm-2.0;
      norm = lstp*M_LNA/(pow(rmod*lsgm,indx)*tgamma(indx));
      fact = -1.0/(rmod*lsgm);
    }
    else
    {
      norm = lstp/(sqrt(PI2)*lsgm);
      fact = -0.5/(lsgm*lsgm);
    }
    if(vb > 1)
    {
      // calculate averages
      s1 = s2 = s3 = s4 = 0.0;
      if(gmod)
      {
        for(m=0,sy=0.0; m<=nstp; m++)
        {
          l = lmin+lstp*m;
          r = pow(10.0,l);
          y = norm*pow(r,indx)*exp(fact*r);
          w = y*r;
          s1 += w; w *= r;
          s2 += w; w *= r;
          s3 += w; w *= r;
          s4 += w;
          sy += y;
        }
      }
      else
      {
        for(m=0,sy=0.0; m<=nstp; m++)
        {
          l = lmin+lstp*m;
          r = pow(10.0,l);
          y = norm*exp(fact*(l-lmod)*(l-lmod));
          w = y*r;
          s1 += w; w *= r;
          s2 += w; w *= r;
          s3 += w; w *= r;
          s4 += w;
          sy += y;
        }
      }
      if(fabs(sy-1.0) > epsilon)
      {
        fprintf(stderr,"Warning, sy=%13.4e\n",sy);
      }
      fprintf(stderr,"Min Log10(R) : %13.4e\n",lmin);
      fprintf(stderr,"Max Log10(R) : %13.4e\n",lmax);
      fprintf(stderr,"#Steps       : %13d\n",nstp);
      fprintf(stderr,"Integrated N : %13.6e\n",sy);
      fprintf(stderr,"Average R    : %13.6e\n",s1/sy);
      fprintf(stderr,"Average R^2  : %13.6e\n",s2/sy);
      fprintf(stderr,"Average R^3  : %13.6e\n",s3/sy);
      fprintf(stderr,"R^2 Weighted Average R: %13.6e\n",s3/s2);
      fprintf(stderr,"R^2 Weighted Variance R/Average R^2: %13.6e\n",s4/(s3*s3/s2)-1.0);
    }
    if(gmod)
    {
      for(m=0,sy=0.0; m<=nstp; m++)
      {
        l = lmin+lstp*m;
        r = pow(10.0,l);
        x = PI2*r/wlen;
        y = norm*pow(r,indx)*exp(fact*r);
        w = y*PI*r*r;
        if(x >= XMAX) break;
        miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,cang,
               &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
               PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);
        aext += QEXT*w;
        asca += QSCA*w;
        for(j=0; j<nang; j++)
        {
          s1 = cabs(S1[j]);
          s2 = cabs(S2[j]);
          phs1[j] += s1*s1*y;
          phs2[j] += s2*s2*y;
        }
        sy += y;
      }
    }
    else
    {
      for(m=0,sy=0.0; m<=nstp; m++)
      {
        l = lmin+lstp*m;
        r = pow(10.0,l);
        x = PI2*r/wlen;
        y = norm*exp(fact*(l-lmod)*(l-lmod));
        w = y*PI*r*r;
        if(x >= XMAX) break;
        miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,cang,
               &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
               PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);
        aext += QEXT*w;
        asca += QSCA*w;
        for(j=0; j<nang; j++)
        {
          s1 = cabs(S1[j]);
          s2 = cabs(S2[j]);
          phs1[j] += s1*s1*y;
          phs2[j] += s2*s2*y;
        }
        sy += y;
      }
    }
    if(fabs(sy-1.0) > epsilon)
    {
      fprintf(stderr,"Warning, sy=%13.4e\n",sy);
    }
    aext /= sy;
    asca /= sy;
    omeg = asca/aext;
    for(j=0; j<nang; j++)
    {
      phs1[j] /= sy;
      phs2[j] /= sy;
      phas[j] = 0.5*(phs1[j]+phs2[j]);
    }
  }

  // normalization
  p1 = p2 = p3 = 0.0;
  for(j=1; j<nang; j++)
  {
    p1 += (phs1[j-1]*sang[j-1]+phs1[j]*sang[j])*dang[j];
    p2 += (phs2[j-1]*sang[j-1]+phs2[j]*sang[j])*dang[j];
    p3 += (phas[j-1]*sang[j-1]+phas[j]*sang[j])*dang[j];
  }
  for(j=0; j<nang; j++)
  {
    phs1[j] /= p1;
    phs2[j] /= p2;
    phas[j] /= p3;
  }
  for(j=1; j<nang; j++)
  {
    asym += (phas[j-1]*sang[j-1]*cang[j-1]+
             phas[j]*sang[j]*cang[j])*dang[j];
  }

  return 0;
}

int Init(void)
{
  lmod = log10(rmod);

  return 0;
}

int Printout(void)
{
  int j;
  FILE *fp;

  fp = (omod==OMOD_APPEND?fopen(out1,"a"):fopen(out1,"w"));
  if(fp == NULL)
  {
    fprintf(stderr,"Error, cannot open %s\n",out1);
    return -1;
  }
  fprintf(fp,"%10.4f %22.15e %22.15e %22.15e %22.15e\n",wlen,aext,asca,omeg,asym);

  fp = (omod==OMOD_APPEND?fopen(out2,"a"):fopen(out2,"w"));
  if(fp == NULL)
  {
    fprintf(stderr,"Error, cannot open %s\n",out2);
    return -1;
  }
  for(j=0; j<nang; j++)
  {
    fprintf(fp,"%10.4f %10.4f %22.15e %22.15e %22.15e\n",wlen,angl[j],phs1[j],phs2[j],phas[j]);
  }
  fclose(fp);

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
    {"ndis",1,0,'N'},
    {"nang",1,0,'M'},
    {"refr",1,0,'R'},
    {"refi",1,0,'I'},
    {"wlen",1,0,'w'},
    {"rmin",1,0,'x'},
    {"rmax",1,0,'X'},
    {"rmod",1,0,'r'},
    {"lsgm",1,0,'s'},
    {"wsgm",1,0,'W'},
    {"nstp",1,0,'n'},
    {"lstp",1,0,'l'},
    {"out1",1,0,'o'},
    {"out2",1,0,'O'},
    {"append",0,0,'a'},
    {"gmod",0,0,'g'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":N:M:R:I:w:x:X:r:s:W:n:l:o:O:agdvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'N':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0) ndis = ntmp;
        else
        {
          fprintf(stderr,"Distribution# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'M':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<=NAMAX) nang = ntmp;
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
      case 'W':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) wsgm = xtmp;
        else
        {
          fprintf(stderr,"R width -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'n':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) nstp = ntmp;
        else
        {
          fprintf(stderr,"Log10(R) #steps -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'l':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) lstp = xtmp;
        else
        {
          fprintf(stderr,"Log10(R) step -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'o':
        strncpy(out1,optarg,MAXLINE);
        break;
      case 'O':
        strncpy(out2,optarg,MAXLINE);
        break;
      case 'a':
        omod = OMOD_APPEND;
        break;
      case 'g':
        gmod = 1;
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

  fprintf(stderr,"mie_size_dist ... calculate mie scattering parameters.\n");
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"mie_size_dist -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," N -ndis    |%s|%s|%s| %d\n",As(e,"Distribution#",n), As(a,"#",n),          Ad(d,NDIS,n),ndis);
  fprintf(stderr," M -nang    |%s|%s|%s| %d\n",As(e,"#Angles",n),       As(a,"#",n),          Ad(d,NANG,n),nang);
  fprintf(stderr," R -refr    |%s|%s|%s| %e\n",As(e,"Re(Ref. index)",n),As(a,"value",n),      Ae(d,REFR,n),refr);
  fprintf(stderr," I -refi    |%s|%s|%s| %e\n",As(e,"Im(Ref. index)",n),As(a,"value",n),      Ae(d,REFI,n),refi);
  fprintf(stderr," w -wlen    |%s|%s|%s| %f\n",As(e,"Wavelength",n),    As(a,"um",n),         Af(d,WLEN,n),wlen);
  fprintf(stderr," x -rmin    |%s|%s|%s| %e\n",As(e,"R min",n),         As(a,"um",n),         Af(d,RMIN,n),rmin);
  fprintf(stderr," X -rmax    |%s|%s|%s| %e\n",As(e,"R max",n),         As(a,"um",n),         Af(d,RMAX,n),rmax);
  fprintf(stderr," r -rmod    |%s|%s|%s| %e\n",As(e,"R mode",n),        As(a,"um",n),         Af(d,RMOD,n),rmod);
  fprintf(stderr," s -lsgm    |%s|%s|%s| %e\n",As(e,"R sigma",n),       As(a,"value",n),      Af(d,LSGM,n),lsgm);
  fprintf(stderr," W -wsgm    |%s|%s|%s| %e\n",As(e,"R width",n),       As(a,"sigma",n),      Af(d,WSGM,n),wsgm);
  fprintf(stderr," n -nstp    |%s|%s|%s| %d\n",As(e,"Log(R) #steps",n), As(a,"#",n),          Ad(d,NSTP,n),nstp);
  fprintf(stderr," l -lstp    |%s|%s|%s| %e\n",As(e,"Log(R) step",n),   As(a,"value",n),      Ae(d,LSTP,n),lstp);
  fprintf(stderr," o -out1    |%s|%s|%s| %s\n",As(e,"Output file 1",n), As(a,"name",n),       As(d,OUT1,n),out1);
  fprintf(stderr," O -out2    |%s|%s|%s| %s\n",As(e,"Output file 2",n), As(a,"name",n),       As(d,OUT2,n),out2);
  fprintf(stderr," a -append  |%s|%s|%s| %d\n",As(e,"Append  mode",n),  As(a,"nothing",n),    Ad(d,0,n),(omod==OMOD_APPEND?1:0));
  fprintf(stderr," g -gmod    |%s|%s|%s| %d\n",As(e,"Gamma   mode",n),  As(a,"nothing",n),    Ad(d,0,n),gmod);
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"-----------------------------------------------------------------------------\n");

  return 0;
}
