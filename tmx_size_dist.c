/***********************************************************/
/* TMX_SIZE_DIST ... calculate mie scattering parameters   */
/* Author: N.Manago                                        */
/* $Revision: 722 $                                        */
/* $Date: 2010-06-22 14:01:51 +0900 (Tue, 22 Jun 2010) $   */
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
#define	NDIS			0				// Distribution#
#define	NANG			361				// #Angles
#define	NCNT			2				// Control #divisions
#define	KMAX			100				// #Quadrature points
#define	SHAP			-1				// Shape#
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
#define	REPS			1.000001			// Aspect ratio
#define	DELT			0.001				// Accuracy of the computations
#define	XMAX			2.0e4				// Max size parameter
#define	MAXLINE			256				// Max #chars in a line
#define	EPSILON			1.0e-14				// A small number
#define	PI			3.141592653589793		// PI
#define	PI2			6.283185307179586		// 2*PI
#define	R_TO_D			57.29577951308232		// 180/PI rad -> deg
#define	D_TO_R			1.745329251994329e-02		// PI/180 deg -> rad
#define	OMOD_OVERWR		1				// Overwrite mode
#define	OMOD_APPEND		2				// Append mode
#define	OUT1			"tmx_out1.dat"			// Output file 1
#define	OUT2			"tmx_out2.dat"			// Output file 2

extern void tmatrix_();

int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);
int TmxCalc(void);
int Printout(void);

int	ndis			= NDIS;
int	nang			= NANG;
int	ncnt			= NCNT;
int	kmax			= KMAX;
int	shap			= SHAP;
int	nstp			= NSTP;
int	omod			= OMOD_OVERWR;
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
double	reps			= REPS;
double	delt			= DELT;
double	xmax			= XMAX;
double	lstp			= LSTP;
double	lmod,lmin,lmax;
double	aext,asca,omeg,asym,depl,depc;
double	angl[NAMAX];
double	phs1[NAMAX];
double	phs2[NAMAX];
double	phas[NAMAX];
double	mm11[NAMAX];
double	mm22[NAMAX];
double	mm33[NAMAX];
double	mm44[NAMAX];
double	mm12[NAMAX];
double	mm34[NAMAX];
char	out1[MAXLINE]		= OUT1;
char	out2[MAXLINE]		= OUT2;

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) return -1;
  if(hp) {Usage(); return 0;}
  if(Init() < 0) return -1;

  if(TmxCalc() < 0) return -1;
  if(Printout() < 0) return -1;

  return 0;
}

int TmxCalc(void)
{
  int j,m;
  double x,y,w,r,l;
  double sy,p1,p2,p3;
  double ss,s1,s2,s3;
  double norm,slop;
  double epsilon = 1.0e-1;
  int iflg,ierr;
  int ndistr,npnax,nkmax;
  int np,npna,ndgs;
  int l1max[NAXMAX];
  double rat,axmax,r1,r2;
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
  double rang[NAMAX];
  double sang[NAMAX];
  double cang[NAMAX];
  double dang[NAMAX];
  double f11[NAMAX];
  double f12[NAMAX];
  char fnam[] = "TmxCalc";

  // initialization
  aext = 0.0;
  asca = 0.0;
  omeg = 0.0;
  asym = 0.0;
  depl = 0.0;
  depc = 0.0;
  for(j=0; j<nang; j++)
  {
    phs1[j] = 0.0;
    phs2[j] = 0.0;
    phas[j] = 0.0;
    mm11[j] = 0.0;
    mm22[j] = 0.0;
    mm33[j] = 0.0;
    mm44[j] = 0.0;
    mm12[j] = 0.0;
    mm34[j] = 0.0;
    angl[j] = (180.0/(nang-1)*j);
    rang[j] = angl[j]*D_TO_R;
    sang[j] = sin(rang[j]);
    cang[j] = cos(rang[j]);
  }
  for(j=1; j<nang; j++)
  {
    dang[j] = PI*fabs(rang[j]-rang[j-1]);
  }

  // T-matrix calculation
  if(lsgm < EPSILON)
  {
    // single size
    iflg   = vb;
    ierr   = 0;
    rat    = 0.5;
    ndistr = 4;
    axmax  = rmod;
    npnax  = 1;
    r1     = 0.9999999*rmod;
    r2     = 1.0000001*rmod;
    b      = 1.0e-1;
    gam    = 0.5;
    nkmax  = -1;
    eps    = reps;
    np     = shap;
    lam    = wlen;
    mrr    = refr;
    mri    = refi;
    ddelt  = delt; // this will be modified by tmatrix_()!
    npna   = (nang<2?2:nang);
    ndgs   = ncnt;
    if(db)
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
             l1max,alp1,alp2,alp3,alp4,bet1,bet2,fmat);
    if(ierr)
    {
      fprintf(stderr,"%s: error, ierr=%d\n",fnam,ierr);
      return -1;
    }
    aext = cext[0];
    asca = csca[0];
    omeg = walb[0];
    asym = gsym[0];
    for(j=0; j<nang; j++)
    {
      f11[j] = fmat[j*6+0];
      f12[j] = fmat[j*6+4];
      phs1[j] = f11[j]-f12[j];
      phs2[j] = f11[j]+f12[j];
      phas[j] = f11[j];
      mm11[j] = fmat[j*6+0];
      mm22[j] = fmat[j*6+1];
      mm33[j] = fmat[j*6+2];
      mm44[j] = fmat[j*6+3];
      mm12[j] = fmat[j*6+4];
      mm34[j] = fmat[j*6+5];
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
      fprintf(stderr,"%s: error, invalid #steps (lmin=%13.4e lmax=%13.4e lstp=%13.4e nstp=%d\n",
                      fnam,lmin,lmax,lstp,nstp);
      return -1;
    }
    norm = lstp/(sqrt(PI2)*lsgm);
    slop = -0.5/(lsgm*lsgm);
    if(vb > 1)
    {
      fprintf(stderr,"Min Log10(R) : %13.4e\n",lmin);
      fprintf(stderr,"Max Log10(R) : %13.4e\n",lmax);
      fprintf(stderr,"#Steps       : %13d\n",nstp);
      // calculate averages
      s1 = s2 = s3 = 0.0;
      for(m=0,sy=0.0; m<=nstp; m++)
      {
        l = lmin+lstp*m;
        r = pow(10.0,l);
        y = norm*exp(slop*(l-lmod)*(l-lmod));
        w = y*r;
        s1 += w; w *= r;
        s2 += w; w *= r;
        s3 += w;
        sy += y;
      }
      if(fabs(sy-1.0) > epsilon)
      {
        fprintf(stderr,"%s: warning, sy=%13.4e\n",fnam,sy);
      }
      s1 /= sy;
      s2 /= sy;
      s3 /= sy;
      fprintf(stderr,"Average R    : %13.6e\n",s1);
      fprintf(stderr,"Average R^2  : %13.6e\n",s2);
      fprintf(stderr,"Average R^3  : %13.6e\n",s3);
    }
    if(ndis == 0)
    {
      for(m=0,sy=ss=0.0; m<=nstp; m++)
      {
        l = lmin+lstp*m;
        r = pow(10.0,l);
        x = PI2*r/wlen;
        y = norm*exp(slop*(l-lmod)*(l-lmod));
        if(x > xmax) break;
        iflg   = vb;
        ierr   = 0;
        rat    = 0.5;
        ndistr = 4;
        axmax  = r;
        npnax  = 1;
        r1     = 0.9999999*r;
        r2     = 1.0000001*r;
        b      = 1.0e-1;
        gam    = 0.5;
        nkmax  = -1;
        eps    = reps;
        np     = shap;
        lam    = wlen;
        mrr    = refr;
        mri    = refi;
        ddelt  = delt; // this will be modified by tmatrix_()!
        npna   = (nang<2?2:nang);
        ndgs   = ncnt;
        if(db)
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
                 l1max,alp1,alp2,alp3,alp4,bet1,bet2,fmat);
        if(ierr)
        {
          fprintf(stderr,"%s: warning, ierr=%d, r=%13.6e\n",fnam,ierr,r);
          break;
        }
        aext += cext[0]*y;
        asca += csca[0]*y;
        for(j=0; j<nang; j++)
        {
          f11[j] = fmat[j*6+0];
          f12[j] = fmat[j*6+4];
          phs1[j] += (f11[j]-f12[j])*csca[0]*y;
          phs2[j] += (f11[j]+f12[j])*csca[0]*y;
          mm11[j] += fmat[j*6+0]*csca[0]*y;
          mm22[j] += fmat[j*6+1]*csca[0]*y;
          mm33[j] += fmat[j*6+2]*csca[0]*y;
          mm44[j] += fmat[j*6+3]*csca[0]*y;
          mm12[j] += fmat[j*6+4]*csca[0]*y;
          mm34[j] += fmat[j*6+5]*csca[0]*y;
        }
        sy += y;
        ss += csca[0]*y;
      }
      if(fabs(sy-1.0) > epsilon)
      {
        fprintf(stderr,"%s: warning, sy=%13.4e\n",fnam,sy);
      }
      aext /= sy;
      asca /= sy;
      omeg = asca/aext;
      for(j=0; j<nang; j++)
      {
        phs1[j] /= ss;
        phs2[j] /= ss;
        phas[j] = 0.5*(phs1[j]+phs2[j]);
        mm11[j] /= ss;
        mm22[j] /= ss;
        mm33[j] /= ss;
        mm44[j] /= ss;
        mm12[j] /= ss;
        mm34[j] /= ss;
      }
    }
    else
    {
      iflg   = vb;
      ierr   = 0;
      rat    = 0.5;
      ndistr = 2;
      axmax  = rmod;
      npnax  = 1;
      r1     = pow(10.0,lmin);
      r2     = pow(10.0,lmax);
      x      = PI2*r2/wlen;
      if(x > xmax)
      {
        r2 = xmax*wlen/PI2;
        x = xmax;
      }
      b      = log(10.0)*lsgm;
      b      = b*b;
      gam    = 0.5;
      nkmax  = kmax;
      eps    = reps;
      np     = shap;
      lam    = wlen;
      mrr    = refr;
      mri    = refi;
      ddelt  = delt; // this will be modified by tmatrix_()!
      npna   = (nang<2?2:nang);
      ndgs   = ncnt;
      if(db)
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
               l1max,alp1,alp2,alp3,alp4,bet1,bet2,fmat);
      if(ierr)
      {
        fprintf(stderr,"%s: error, ierr=%d\n",fnam,ierr);
        return -1;
      }
      aext = cext[0];
      asca = csca[0];
      omeg = walb[0];
      asym = gsym[0];
      for(j=0; j<nang; j++)
      {
        f11[j] = fmat[j*6+0];
        f12[j] = fmat[j*6+4];
        phs1[j] = f11[j]-f12[j];
        phs2[j] = f11[j]+f12[j];
        phas[j] = f11[j];
        mm11[j] = fmat[j*6+0];
        mm22[j] = fmat[j*6+1];
        mm33[j] = fmat[j*6+2];
        mm44[j] = fmat[j*6+3];
        mm12[j] = fmat[j*6+4];
        mm34[j] = fmat[j*6+5];
      }
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
  if(asym == 0.0)
  {
    for(j=1; j<nang; j++)
    {
      asym += (phas[j-1]*sang[j-1]*cang[j-1]+
               phas[j]*sang[j]*cang[j])*dang[j];
    }
  }
  j = nang-1;
  depl = (mm11[j]-mm22[j])/(mm11[j]+mm22[j]+2.0*mm12[j]);
  depc = (mm11[j]+mm44[j])/(mm11[j]-mm44[j]);

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
  fprintf(fp,"%10.4f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",wlen,aext,asca,omeg,asym,depl,depc);

  fp = (omod==OMOD_APPEND?fopen(out2,"a"):fopen(out2,"w"));
  if(fp == NULL)
  {
    fprintf(stderr,"Error, cannot open %s\n",out2);
    return -1;
  }
  for(j=0; j<nang; j++)
  {
    fprintf(fp,"%10.4f %10.4f %13.6e %13.6e %13.6e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
                wlen,angl[j],phs1[j],phs2[j],phas[j],mm11[j],mm22[j],mm33[j],mm44[j],mm12[j],mm34[j]);
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
    {"ncnt",1,0,'c'},
    {"kmax",1,0,'k'},
    {"shap",1,0,'S'},
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
    {"xmax",1,0,'U'},
    {"reps",1,0,'e'},
    {"delt",1,0,'D'},
    {"out1",1,0,'o'},
    {"out2",1,0,'O'},
    {"append",0,0,'a'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":N:M:c:k:S:R:I:w:x:X:r:s:W:n:l:U:e:D:o:O:advh",long_options,&option_index);
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
      case 'c':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) ncnt = ntmp;
        else
        {
          fprintf(stderr,"Control# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'k':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>-2) kmax = ntmp;
        else
        {
          fprintf(stderr,"#Quadrature points -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'S':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0') shap = ntmp;
        else
        {
          fprintf(stderr,"Shape# -> out of range %s\n",optarg);
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
      case 'U':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) xmax = xtmp;
        else
        {
          fprintf(stderr,"X Max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'e':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) reps = xtmp;
        else
        {
          fprintf(stderr,"Aspect ratio -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'D':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) delt = xtmp;
        else
        {
          fprintf(stderr,"Accuracy -> out of range %s\n",optarg);
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

  fprintf(stderr,"tmx_size_dist ... calculate scattering parameters using T-matrix code.\n");
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"tmx_size_dist -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," N -ndis    |%s|%s|%s| %d\n",As(e,"Distribution#",n), As(a,"#",n),          Ad(d,NDIS,n),ndis);
  fprintf(stderr," M -nang    |%s|%s|%s| %d\n",As(e,"#Angles",n),       As(a,"#",n),          Ad(d,NANG,n),nang);
  fprintf(stderr," c -ncnt    |%s|%s|%s| %d\n",As(e,"Control#",n),      As(a,"#",n),          Ad(d,NCNT,n),ncnt);
  fprintf(stderr," k -kmax    |%s|%s|%s| %d\n",As(e,"#Quad points",n),  As(a,"#",n),          Ad(d,KMAX,n),kmax);
  fprintf(stderr," S -shap    |%s|%s|%s| %d\n",As(e,"Shape#",n),        As(a,"#",n),          Ad(d,SHAP,n),shap);
  fprintf(stderr," R -refr    |%s|%s|%s| %e\n",As(e,"Re(Ref. index)",n),As(a,"value",n),      Ae(d,REFR,n),refr);
  fprintf(stderr," I -refi    |%s|%s|%s| %e\n",As(e,"Im(Ref. index)",n),As(a,"value",n),      Ae(d,REFI,n),refi);
  fprintf(stderr," w -wlen    |%s|%s|%s| %f\n",As(e,"Wavelength",n),    As(a,"um",n),         Af(d,WLEN,n),wlen);
  fprintf(stderr," x -rmin    |%s|%s|%s| %e\n",As(e,"R min",n),         As(a,"um",n),         Af(d,RMIN,n),rmin);
  fprintf(stderr," X -rmax    |%s|%s|%s| %e\n",As(e,"R max",n),         As(a,"um",n),         Af(d,RMAX,n),rmax);
  fprintf(stderr," r -rmod    |%s|%s|%s| %e\n",As(e,"R mode",n),        As(a,"um",n),         Af(d,RMOD,n),rmod);
  fprintf(stderr," s -lsgm    |%s|%s|%s| %e\n",As(e,"Log(R) sigma",n),  As(a,"value",n),      Af(d,LSGM,n),lsgm);
  fprintf(stderr," W -wsgm    |%s|%s|%s| %e\n",As(e,"Log(R) width",n),  As(a,"sigma",n),      Af(d,WSGM,n),wsgm);
  fprintf(stderr," n -nstp    |%s|%s|%s| %d\n",As(e,"Log(R) #steps",n), As(a,"#",n),          Ad(d,NSTP,n),nstp);
  fprintf(stderr," l -lstp    |%s|%s|%s| %e\n",As(e,"Log(R) step",n),   As(a,"value",n),      Ae(d,LSTP,n),lstp);
  fprintf(stderr," U -xmax    |%s|%s|%s| %e\n",As(e,"X Max",n),         As(a,"value",n),      Ae(d,XMAX,n),xmax);
  fprintf(stderr," e -reps    |%s|%s|%s| %f\n",As(e,"Aspect ratio",n),  As(a,"value",n),      Af(d,REPS,n),reps);
  fprintf(stderr," D -delt    |%s|%s|%s| %e\n",As(e,"Accuracy",n),      As(a,"value",n),      Ae(d,DELT,n),delt);
  fprintf(stderr," o -out1    |%s|%s|%s| %s\n",As(e,"Output file 1",n), As(a,"name",n),       As(d,OUT1,n),out1);
  fprintf(stderr," O -out2    |%s|%s|%s| %s\n",As(e,"Output file 2",n), As(a,"name",n),       As(d,OUT2,n),out2);
  fprintf(stderr," a -append  |%s|%s|%s| %d\n",As(e,"Append  mode",n),  As(a,"nothing",n),    Ad(d,0,n),(omod==OMOD_APPEND?1:0));
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"-----------------------------------------------------------------------------\n");

  return 0;
}
