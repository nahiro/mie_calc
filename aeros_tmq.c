/***********************************************************/
/* AEROS_TMX ... Aerosol generator with T-Matrix           */
/* Author: N.Manago                                        */
/* $Revision: 1104 $                                       */
/* $Date: 2015-09-02 00:16:06 +0900 (Wed, 02 Sep 2015) $   */
/***********************************************************/
#define	MIE_TRUE		1				// True  value
#define	MIE_FALSE		0				// False value
#define	MIE_XMAX		2.0e4				// Max size parameter
#define	MIE_RMIN		NAN				// R min in um
#define	MIE_RMAX		NAN				// R max in um
#define	MIE_LMIN		-3.0				// Min Log10(R)
#define	MIE_LMAX		+2.0				// Max Log10(R)
#define	MIE_LSTP		1.0e-3				// Log10(R) step
#define	MIE_WSGM		5.0				// Log10(R) width in sigma
#define	MIE_NSTP		-1				// #Steps
#define	TMX_NDIS		0				// Distribution#
#define	TMX_NANG		361				// #Angles
#define	TMX_NDGS		2				// Control #divisions
#define	TMX_KMAX		100				// #Quad points
#define	TMX_SHAP		-1				// Shape#
#define	TMX_REPS		1.000001			// Aspect ratio or deformation parameter
#define	TMX_PGAM		0.5     			// Parameter for the modified gamma distribution
#define	TMX_DELT		0.001				// Accuracy of the computations
#define	TMX_XMAX		MIE_XMAX			// Max size parameter
#define	DELTA			1.0e-13				// A small number
#define	NPN1			300
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
#include "aeros_common.c"

int	mie_nstp		= MIE_NSTP;			// #Steps
int	tmx_ndis		= TMX_NDIS;			// Distribution#
int	tmx_nang		= TMX_NANG;			// #Angles
int	tmx_ndgs		= TMX_NDGS;			// Control #divisions
int	tmx_kmax		= TMX_KMAX;			// #Quad points
int	tmx_shap		= TMX_SHAP;			// Shape#
double	mie_rmin		= MIE_RMIN;			// R min in um
double	mie_rmax		= MIE_RMAX;			// R max in um
double	mie_lstp		= MIE_LSTP;			// Log10(R) step
double	mie_wsgm		= MIE_WSGM;			// Log10(R) width in sigma
double	mie_lmin;
double	mie_lmax;
double	mie_lmod_com[MIE_MAXCOMP];				// Log10(Mode radius in um)
double	tmx_reps		= TMX_REPS;
double	tmx_pgam		= TMX_PGAM;
double	tmx_delt		= TMX_DELT;
double	tmx_xmax		= TMX_XMAX;
int	cnt_gamma		= 0;				// Gamma   mode

extern void tmatrix_();

int MieCalc(void)
{
  int i,j,k,m,n;
  int n_step = 0;
  int err = 0;
  double x,y,w,r,l;
  double sy,sp,sinv;
  double s1,s2,s3,s4;
  double norm,fact,indx;
  double epsilon = 1.0e-1;
  int iflg,ierr;
  int ndistr,npnax,nkmax;
  int np,npna,ndgs;
  int lmax[NAXMAX];
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
  double angl[NAMAX];
  double f11[NAMAX];
  double f12[NAMAX];
  double a1[MIE_MAXDATA];
  double b1[MIE_MAXDATA];
  char fnam[] = "MieCalc";

  if(!isnan(mie_rmin))
  {
    mie_lmin = log10(mie_rmin);
    if(mie_lmin < MIE_LMIN)
    {
      mie_lmin = MIE_LMIN;
    }
  }
  if(!isnan(mie_rmax))
  {
    mie_lmax = log10(mie_rmax);
    if(mie_lmax > MIE_LMAX)
    {
      mie_lmax = MIE_LMAX;
    }
  }
  if(mie_nstp > 0)
  {
    n_step = mie_nstp;
  }
  if(!isnan(mie_rmin) && !isnan(mie_rmax))
  {
    if(mie_nstp > 0)
    {
      mie_lstp = (mie_lmax-mie_lmin)/mie_nstp;
    }
    else
    {
      n_step = (int)((mie_lmax-mie_lmin)/mie_lstp+0.5);
      if(n_step < 1)
      {
        fprintf(stderr,"%s: error, invalid #steps (lmin=%13.4e lmax=%13.4e lstp=%13.4e n_step=%d\n",
                        fnam,mie_lmin,mie_lmax,mie_lstp,n_step);
        return -1;
      }
    }
  }

  for(n=0; n<mie_n_comp; n++)
  {
    mie_lmod_com[n] = log10(mie_rmod_com[n]);
  }

  // initialization
  for(n=0; n<mie_n_comp; n++)
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      mie_aext_com[n][i] = 0.0;
      mie_asca_com[n][i] = 0.0;
      mie_asym_com[n][i] = 0.0;
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        mie_phs1_com[n][k] = 0.0;
        mie_phs2_com[n][k] = 0.0;
        mie_phas_com[n][k] = 0.0;
      }
    }
  }
  for(j=0; j<tmx_nang; j++)
  {
    angl[j] = (180.0/(tmx_nang-1)*j);
  }

  // mie calculation
  for(n=0; n<mie_n_comp; n++)
  {
    if(isnan(mie_rmin) || isnan(mie_rmax))
    {
      if(isnan(mie_rmin))
      {
        mie_lmin = mie_lmod_com[n]-mie_lsgm_com[n]*mie_wsgm;
        if(mie_lmin < MIE_LMIN)
        {
          mie_lmin = MIE_LMIN;
        }
      }
      if(isnan(mie_rmax))
      {
        mie_lmax = mie_lmod_com[n]+mie_lsgm_com[n]*mie_wsgm;
        if(mie_lmax > MIE_LMAX)
        {
          mie_lmax = MIE_LMAX;
        }
      }
      if(mie_nstp > 0)
      {
        mie_lstp = (mie_lmax-mie_lmin)/mie_nstp;
      }
      else
      {
        n_step = (int)((mie_lmax-mie_lmin)/mie_lstp+0.5);
        if(n_step < 1)
        {
          fprintf(stderr,"%s: error, invalid #steps (lmin=%13.4e lmax=%13.4e lstp=%13.4e n_step=%d\n",
                          fnam,mie_lmin,mie_lmax,mie_lstp,n_step);
          err = 1;
          break;
        }
      }
    }
    if(cnt_gamma)
    {
      indx = 1.0/mie_lsgm_com[n]-2.0;
      norm = mie_lstp*log(10.0)/(pow(mie_rmod_com[n]*mie_lsgm_com[n],indx)*tgamma(indx));
      fact = -1.0/(mie_rmod_com[n]*mie_lsgm_com[n]);
    }
    else
    {
      norm = mie_lstp/(sqrt(PI2)*mie_lsgm_com[n]);
      fact = -0.5/(mie_lsgm_com[n]*mie_lsgm_com[n]);
    }
    if(cnt_vb > 1 || cnt_xmod == 2)
    {
      if(cnt_vb > 1)
      {
        fprintf(stderr,"Component#   : %2d\n",n);
        fprintf(stderr,"Min Log10(R) : %13.4e\n",mie_lmin);
        fprintf(stderr,"Max Log10(R) : %13.4e\n",mie_lmax);
        fprintf(stderr,"#Steps       : %13d\n",n_step);
      }
      // calculate averages
      s1 = s2 = s3 = s4 = 0.0;
      if(cnt_gamma)
      {
        for(m=0,sy=0.0; m<=n_step; m++)
        {
          l = mie_lmin+mie_lstp*m;
          r = pow(10.0,l);
          y = norm*pow(r,indx)*exp(fact*r);
          if(cnt_vb > 2)
          {
            fprintf(stderr,"%13.6e %13.6e\n",l,y/mie_lstp);
          }
          w = y;
          sy += w; w *= r;
          s1 += w; w *= r;
          s2 += w; w *= r;
          s3 += w; w *= r;
          s4 += w;
        }
      }
      else
      {
        for(m=0,sy=0.0; m<=n_step; m++)
        {
          l = mie_lmin+mie_lstp*m;
          r = pow(10.0,l);
          y = norm*exp(fact*(l-mie_lmod_com[n])*(l-mie_lmod_com[n]));
          if(cnt_vb > 2)
          {
            fprintf(stderr,"%13.6e %13.6e\n",l,y/mie_lstp);
          }
          w = y;
          sy += w; w *= r;
          s1 += w; w *= r;
          s2 += w; w *= r;
          s3 += w; w *= r;
          s4 += w;
        }
      }
      mie_vcom_com[n] = s3/sy*PI4_3;
      if(cnt_vb > 1)
      {
        fprintf(stderr,"Integrated N : %13.6e\n",sy);
        fprintf(stderr,"Average R    : %13.6e\n",s1/sy);
        fprintf(stderr,"Average R^2  : %13.6e\n",s2/sy);
        fprintf(stderr,"Average R^3  : %13.6e\n",s3/sy);
        fprintf(stderr,"R^2 Weighted Average R: %13.6e\n",s3/s2);
        fprintf(stderr,"R^2 Weighted Variance R/Average R^2: %13.6e\n",s4/(s3*s3/s2)-1.0);
      }
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      if(tmx_ndis == 0)
      {
        if(cnt_gamma)
        {
          for(m=0,sy=0.0; m<=n_step; m++)
          {
            l = mie_lmin+mie_lstp*m;
            r = pow(10.0,l);
            x = PI2*r/(mie_wlen_um[i]);
            y = norm*pow(r,indx)*exp(fact*r);
            if(x > tmx_xmax) break;
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
            lam    = mie_wlen_um[i];
            mrr    = mie_refr_com[n][i];
            mri    = mie_refi_com[n][i];
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
                              fnam,ierr,r,lam,x,eps);
              break;
            }
            mie_aext_com[n][i] += cext[0]*y;
            mie_asca_com[n][i] += csca[0]*y;
            for(j=0; j<tmx_nang; j++)
            {
              f11[j] = fmat[j*6+0];
              f12[j] = fmat[j*6+4];
            }
            SamplingE(tmx_nang,angl,f11,mie_n_angl,mie_angl,a1,1.0);
            SamplingE(tmx_nang,angl,f12,mie_n_angl,mie_angl,b1,1.0);
            for(j=0; j<mie_n_angl; j++)
            {
              k = mie_n_angl*i+j;
              mie_phs1_com[n][k] += (a1[j]-b1[j])*csca[0]*y;
              mie_phs2_com[n][k] += (a1[j]+b1[j])*csca[0]*y;
            }
            sy += y;
          }
        }
        else
        {
          for(m=0,sy=0.0; m<=n_step; m++)
          {
            l = mie_lmin+mie_lstp*m;
            r = pow(10.0,l);
            x = PI2*r/(mie_wlen_um[i]);
            y = norm*exp(fact*(l-mie_lmod_com[n])*(l-mie_lmod_com[n]));
            if(x > tmx_xmax) break;
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
            lam    = mie_wlen_um[i];
            mrr    = mie_refr_com[n][i];
            mri    = mie_refi_com[n][i];
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
                              fnam,ierr,r,lam,x,eps);
              break;
            }
            mie_aext_com[n][i] += cext[0]*y;
            mie_asca_com[n][i] += csca[0]*y;
            for(j=0; j<tmx_nang; j++)
            {
              f11[j] = fmat[j*6+0];
              f12[j] = fmat[j*6+4];
            }
            SamplingE(tmx_nang,angl,f11,mie_n_angl,mie_angl,a1,1.0);
            SamplingE(tmx_nang,angl,f12,mie_n_angl,mie_angl,b1,1.0);
            for(j=0; j<mie_n_angl; j++)
            {
              k = mie_n_angl*i+j;
              mie_phs1_com[n][k] += (a1[j]-b1[j])*csca[0]*y;
              mie_phs2_com[n][k] += (a1[j]+b1[j])*csca[0]*y;
            }
            sy += y;
          }
        }
        if(fabs(sy-1.0) > epsilon)
        {
          fprintf(stderr,"%s: warning, sy=%13.4e\n",fnam,sy);
        }
        sinv = 1.0/sy;
        mie_aext_com[n][i] *= sinv;
        mie_asca_com[n][i] *= sinv;
        for(j=0; j<mie_n_angl; j++)
        {
          k = mie_n_angl*i+j;
          mie_phas_com[n][k] = 0.5*(mie_phs1_com[n][k]+mie_phs2_com[n][k]);
        }
      }
      else
      {
        iflg   = cnt_vb;
        ierr   = 0;
        rat    = 0.5; // equal-surface-area-sphere radius
        ndistr = tmx_ndis;
        axmax  = mie_rmod_com[n];
        npnax  = 1;
        r1     = pow(10.0,mie_lmin);
        r2     = pow(10.0,mie_lmax);
        x      = PI2*r2/(mie_wlen_um[i]);
        if(x > tmx_xmax)
        {
          r2 = tmx_xmax*mie_wlen_um[i]/PI2;
          x = tmx_xmax;
        }
        switch(tmx_ndis)
        {
          case 1: // modified gamma distribution
            b = mie_lsgm_com[n];
            break;
          case 3: // power law distribution
            b = mie_lsgm_com[n];
            break;
          case 4: // gamma distribution
            b = mie_lsgm_com[n];
            break;
          default: // log normal distribution
            b = log(10.0)*mie_lsgm_com[n];
            b = b*b;
            break;
        }
        gam    = tmx_pgam;
        nkmax  = tmx_kmax;
        eps    = tmx_reps;
        np     = tmx_shap;
        lam    = mie_wlen_um[i];
        mrr    = mie_refr_com[n][i];
        mri    = mie_refi_com[n][i];
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
          fprintf(stderr,"%s: error, ierr=%d, r1=%13.6e, r2=%13.6e, lam=%13.6f, x=%13.6e, eps=%13.6e\n",
                          fnam,ierr,r1,r2,lam,x,eps);
          return -1;
        }
        mie_aext_com[n][i] = cext[0];
        mie_asca_com[n][i] = csca[0];
        mie_asym_com[n][i] = gsym[0];
        for(j=0; j<tmx_nang; j++)
        {
          f11[j] = fmat[j*6+0];
          f12[j] = fmat[j*6+4];
        }
        SamplingE(tmx_nang,angl,f11,mie_n_angl,mie_angl,a1,1.0);
        SamplingE(tmx_nang,angl,f12,mie_n_angl,mie_angl,b1,1.0);
        for(j=0; j<mie_n_angl; j++)
        {
          k = mie_n_angl*i+j;
          mie_phs1_com[n][k] = a1[j]-b1[j];
          mie_phs2_com[n][k] = a1[j]+b1[j];
          mie_phas_com[n][k] = a1[j];
        }
      }
      sp = 0.0;
      for(j=1; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        sp += (mie_phas_com[n][k-1]*mie_angl_sin[j-1]+mie_phas_com[n][k]*mie_angl_sin[j])*mie_angl_dif[j];
      }
      sinv = 1.0/sp;
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        mie_phs1_com[n][k] *= sinv;
        mie_phs2_com[n][k] *= sinv;
        mie_phas_com[n][k] *= sinv;
      }
      if(tmx_ndis == 0)
      {
        for(j=1; j<mie_n_angl; j++)
        {
          k = mie_n_angl*i+j;
          mie_asym_com[n][i] += (mie_phas_com[n][k-1]*mie_angl_sin[j-1]*mie_angl_cos[j-1]+
                                 mie_phas_com[n][k]*mie_angl_sin[j]*mie_angl_cos[j])*mie_angl_dif[j];
        }
      }
    }
    if(err) break;
  }

  if(err)
  {
    return -1;
  }

  return 0;
}

int Init(void)
{
  return CommonInit();
}

void Finish(void)
{
  CommonFinish();
}

int GetOwnOpt(int argn,char **args)
{
  int c,rt;
  int option_index = 0;
  int this_option_optind;
  int ntmp;
  char *endp;
  double xtmp;
  struct option long_options[] =
  {
    {"mie_rmin",1,0,'x'},
    {"mie_rmax",1,0,'X'},
    {"mie_nstp",1,0,'n'},
    {"mie_lstp",1,0,'s'},
    {"mie_wsgm",1,0,'t'},
    {"tmx_ndis",1,0,'N'},
    {"tmx_ndgs",1,0,'c'},
    {"tmx_nang",1,0,'M'},
    {"tmx_kmax",1,0,'k'},
    {"tmx_shap",1,0,'S'},
    {"tmx_reps",1,0,'e'},
    {"tmx_pgam",1,0,'p'},
    {"tmx_delt",1,0,'D'},
    {"tmx_xmax",1,0,'U'},
    {"cnt_gamma",0,0,'G'},
  };

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-x:X:n:s:t:N:c:M:k:S:e:p:D:U:G",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'x':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_rmin = xtmp;
        else
        {
          fprintf(stderr,"R min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'X':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_rmax = xtmp;
        else
        {
          fprintf(stderr,"R max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'n':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) mie_nstp = ntmp;
        else
        {
          fprintf(stderr,"#Steps -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 's':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_lstp = xtmp;
        else
        {
          fprintf(stderr,"Log10(R) step -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 't':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_wsgm = xtmp;
        else
        {
          fprintf(stderr,"Log10(R) width -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'N':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0) tmx_ndis = ntmp;
        else
        {
          fprintf(stderr,"Distribution# -> out of range %s\n",optarg);
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
      case 'k':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>-2) tmx_kmax = ntmp;
        else
        {
          fprintf(stderr,"#Quad points -> out of range %s\n",optarg);
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
      case 'p':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') tmx_pgam = xtmp;
        else
        {
          fprintf(stderr,"Gamma parameter -> out of range %s\n",optarg);
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
      case 'U':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) tmx_xmax = xtmp;
        else
        {
          fprintf(stderr,"X Max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'G':
        cnt_gamma = 1;
        break;
      case '?':
        if(cnt_optn >= MAXNOPT)
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
        if(cnt_optn >= MAXNOPT)
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
  int  n = 16;
  char e[MAXLINE];
  char a[MAXLINE];
  char d[MAXLINE];
  char optstring[] = "fgwlLaWAPQiIruxXnstNcMkSepDUmGoOdvh";

  CommonUsage("aeros_tmq",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'x':
        fprintf(stderr," x -mie_rmin      |%s|%s|%s| %e\n",As(e,"R Min",n),          As(a,"um",n),         Ae(d,MIE_RMIN,n),mie_rmin);
        break;
      case 'X':
        fprintf(stderr," X -mie_rmax      |%s|%s|%s| %e\n",As(e,"R Max",n),          As(a,"um",n),         Ae(d,MIE_RMAX,n),mie_rmax);
        break;
      case 'n':
        fprintf(stderr," n -mie_nstp      |%s|%s|%s| %d\n",As(e,"#Steps",n),         As(a,"#",n),          Ad(d,MIE_NSTP,n),mie_nstp);
        break;
      case 's':
        fprintf(stderr," s -mie_lstp      |%s|%s|%s| %e\n",As(e,"Log10(R) step",n),  As(a,"value",n),      Ae(d,MIE_LSTP,n),mie_lstp);
        break;
      case 't':
        fprintf(stderr," t -mie_wsgm      |%s|%s|%s| %f\n",As(e,"Log10(R) width",n), As(a,"sigma",n),      Af(d,MIE_WSGM,n),mie_wsgm);
        break;
      case 'N':
        fprintf(stderr," N -tmx_ndis      |%s|%s|%s| %d\n",As(e,"Distribution#",n),  As(a,"#",n),          Ad(d,TMX_NDIS,n),tmx_ndis);
        break;
      case 'c':
        fprintf(stderr," c -tmx_ndgs      |%s|%s|%s| %d\n",As(e,"Control#",n),       As(a,"#",n),          Ad(d,TMX_NDGS,n),tmx_ndgs);
        break;
      case 'M':
        fprintf(stderr," M -tmx_nang      |%s|%s|%s| %d\n",As(e,"#Angles",n),        As(a,"#",n),          Ad(d,TMX_NANG,n),tmx_nang);
        break;
      case 'k':
        fprintf(stderr," k -tmx_kmax      |%s|%s|%s| %d\n",As(e,"#Quad points",n),   As(a,"#",n),          Ad(d,TMX_KMAX,n),tmx_kmax);
        break;
      case 'S':
        fprintf(stderr," S -tmx_shap      |%s|%s|%s| %d\n",As(e,"Shape#",n),         As(a,"#",n),          Ad(d,TMX_SHAP,n),tmx_shap);
        break;
      case 'e':
        fprintf(stderr," e -tmx_reps      |%s|%s|%s| %e\n",As(e,"Aspect ratio",n),   As(a,"value",n),      Af(d,TMX_REPS,n),tmx_reps);
        break;
      case 'p':
        fprintf(stderr," p -tmx_pgam      |%s|%s|%s| %e\n",As(e,"Gamma parameter",n),As(a,"value",n),      Af(d,TMX_PGAM,n),tmx_pgam);
        break;
      case 'D':
        fprintf(stderr," D -tmx_delt      |%s|%s|%s| %e\n",As(e,"Accuracy",n),       As(a,"value",n),      Ae(d,TMX_DELT,n),tmx_delt);
        break;
      case 'U':
        fprintf(stderr," U -tmx_xmax      |%s|%s|%s| %e\n",As(e,"X Max",n),          As(a,"value",n),      Ae(d,TMX_XMAX,n),tmx_xmax);
        break;
      case 'G':
        fprintf(stderr," G -gamma         |%s|%s|%s| %d\n",As(e,"Gamma   mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_gamma);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);
  fprintf(stderr,"Distribution#\n");
  fprintf(stderr,"  0 ... Log-normal distribution (monodisperse particles)\n");
  fprintf(stderr,"  1 ... Modified gamma distribution (rmod->alpha,rsgm->rc)\n");
  fprintf(stderr,"  2 ... Log-normal distribution\n");
  fprintf(stderr,"  3 ... Power law distribution (rmod->reff,rsgm->veff)\n");
  fprintf(stderr,"  4 ... Gamma distribution (rmod->reff,rsgm->veff)\n");
  CommonUsage(NULL,3);
  CommonUsage(NULL,4);
  CommonUsage(NULL,5);

  return 0;
}
