/***********************************************************/
/* AEROS_MIEV0 ... Aerosol generator with MIEV0            */
/* Author: N.Manago                                        */
/* $Revision: 1104 $                                       */
/* $Date: 2015-09-02 00:16:06 +0900 (Wed, 02 Sep 2015) $   */
/***********************************************************/
#define	MIE_TRUE		1				// True  value
#define	MIE_FALSE		0				// False value
#define	MIE_XMAX		2.0e4				// Max size parameter
#define	MIE_RMIN		NAN				// R min in um
#define	MIE_RMAX		NAN				// R max in um
#define	MIE_LMIN		-4.0				// Min Log10(R)
#define	MIE_LMAX		+3.0				// Max Log10(R)
#define	MIE_LSTP		1.0e-3				// Log10(R) step
#define	MIE_WSGM		5.0				// Log10(R) width in sigma
#define	MIE_NSTP		-1				// #Steps
#define	DELTA			1.0e-13				// A small number
#include <complex.h>
#include "aeros_common.c"

double	mie_rmin		= MIE_RMIN;			// R min in um
double	mie_rmax		= MIE_RMAX;			// R max in um
double	mie_lmin;						// Log10(R) min
double	mie_lmax;						// Log10(R) max
double	mie_lstp		= MIE_LSTP;			// Log10(R) step
double	mie_lmod_com[MIE_MAXCOMP];				// Log10(Mode radius in um)
double	mie_wsgm		= MIE_WSGM;			// Log10(R) width in sigma
int	mie_nstp		= MIE_NSTP;			// #Steps
int	cnt_anyang		= MIE_FALSE;			// Symmetric angle mode

extern void miev0_();

int MieCalc(void)
{
  int i,j,k,m,n;
  int n_step = 0;
  int err = 0;
  int ANYANG,PERFCT,PRNT[2];
  int IPOLZN,MOMDIM,NUMANG,NMOM;
  double v,w,x,y,r,l;
  double sw,sy,sp,sinv;
  double s1,s2,s3,s4;
  double norm,fact,indx;
  const double epsilon = 1.0e-1;
  double GQSC,MIMCUT,PMOM[2][2],QEXT,QSCA,SPIKE;
  double _Complex CREFIN,SFORW,SBACK,*S1,*S2,TFORW[2],TBACK[2];
  const char fnam[] = "MieCalc";

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

  // initialization
  S1 = (double _Complex *)malloc(mie_n_angl*sizeof(double _Complex));
  S2 = (double _Complex *)malloc(mie_n_angl*sizeof(double _Complex));
  if(S1==NULL || S2==NULL)
  {
    fprintf(stderr,"%s: failed in allocating memory\n",fnam);
    return -1;
  }
  PERFCT = MIE_FALSE;
  MIMCUT = DELTA;
  ANYANG = cnt_anyang;
  NUMANG = mie_n_angl;
  PRNT[0] = (cnt_db>1?MIE_TRUE:MIE_FALSE);
  PRNT[1] = (cnt_db>0?MIE_TRUE:MIE_FALSE);
  NMOM = 0;
  IPOLZN = 0;
  MOMDIM = 1;
  for(n=0; n<mie_n_comp; n++)
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      mie_aext_com[n][i] = 0.0;
      mie_asca_com[n][i] = 0.0;
      mie_asym_com[n][i] = 0.0;
      for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
      {
        mie_phs1_com[n][k] = 0.0;
        mie_phs2_com[n][k] = 0.0;
        mie_phas_com[n][k] = 0.0;
      }
    }
  }

  // mie calculation
  for(n=0; n<mie_n_comp; n++)
  {
    if(mie_size_func[n] == MIE_FUNC_FILE)
    {
      n_step = mie_n_size[n]-1;
      if(cnt_vb > 1 || cnt_xmod == 2)
      {
        mie_lmin = log10(mie_xval_com[n][0]);
        mie_lmax = log10(mie_xval_com[n][n_step]);
      }
    }
    else
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
      switch(mie_size_func[n])
      {
        case MIE_FUNC_LOGNORMAL:
        case MIE_FUNC_EFF_LOGNORMAL:
          norm = mie_lstp/(sqrt(PI2)*mie_lsgm_com[n]);
          fact = -0.5/(mie_lsgm_com[n]*mie_lsgm_com[n]);
          break;
        case MIE_FUNC_GAMMA:
          indx = 1.0/mie_lsgm_com[n]-2.0;
          norm = mie_lstp*M_LNA/(pow(mie_rmod_com[n]*mie_lsgm_com[n],indx)*tgamma(indx));
          fact = -1.0/(mie_rmod_com[n]*mie_lsgm_com[n]);
          break;
        default:
          fprintf(stderr,"%s: error, mie_size_func[%d]=%c\n",fnam,n,mie_size_func[n]);
          err = 1;
          break;
      }
      if(err) break;
    }
    if(cnt_vb > 1 || cnt_xmod == 2)
    {
      // calculate averages
      s1 = s2 = s3 = s4 = 0.0;
      switch(mie_size_func[n])
      {
        case MIE_FUNC_FILE:
          switch(mie_size_ytype[n])
          {
            case MIE_YTYPE_NUMBER: // Number distribution
              for(m=0,sy=0.0; m<=n_step; m++)
              {
                r = mie_xval_com[n][m];
                y = mie_yval_com[n][m];
                w = y;
                sy += w; w *= r;
                s1 += w; w *= r;
                s2 += w; w *= r;
                s3 += w; w *= r;
                s4 += w;
              }
              break;
            case MIE_YTYPE_VOLUME: // Volume distribution
              for(m=0,sy=0.0; m<=n_step; m++)
              {
                r = mie_xval_com[n][m];
                v = mie_yval_com[n][m];
                y = v/(PI4_3*r*r*r);
                w = y;
                sy += w; w *= r;
                s1 += w; w *= r;
                s2 += w; w *= r;
                s3 += w; w *= r;
                s4 += w;
              }
              break;
            default:
              fprintf(stderr,"%s: error, mie_size_ytype[%d]=%c\n",fnam,n,mie_size_ytype[n]);
              err = 1;
              break;
          }
          break;
        case MIE_FUNC_LOGNORMAL:
        case MIE_FUNC_EFF_LOGNORMAL:
          switch(mie_size_ytype[n])
          {
            case MIE_YTYPE_NUMBER: // Number distribution
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
              break;
            case MIE_YTYPE_VOLUME: // Volume distribution
              for(m=0,sy=0.0; m<=n_step; m++)
              {
                l = mie_lmin+mie_lstp*m;
                r = pow(10.0,l);
                v = norm*exp(fact*(l-mie_lmod_com[n])*(l-mie_lmod_com[n]));
                y = v/(PI4_3*r*r*r);
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
              break;
            default:
              fprintf(stderr,"%s: error, mie_size_ytype[%d]=%c\n",fnam,n,mie_size_ytype[n]);
              err = 1;
              break;
          }
          break;
        case MIE_FUNC_GAMMA:
          switch(mie_size_ytype[n])
          {
            case MIE_YTYPE_NUMBER: // Number distribution
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
              break;
            case MIE_YTYPE_VOLUME: // Volume distribution
              for(m=0,sy=0.0; m<=n_step; m++)
              {
                l = mie_lmin+mie_lstp*m;
                r = pow(10.0,l);
                v = norm*pow(r,indx)*exp(fact*r);
                y = v/(PI4_3*r*r*r);
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
              break;
            default:
              fprintf(stderr,"%s: error, mie_size_ytype[%d]=%c\n",fnam,n,mie_size_ytype[n]);
              err = 1;
              break;
          }
          break;
        default:
          fprintf(stderr,"%s: error, mie_size_func[%d]=%c\n",fnam,n,mie_size_func[n]);
          err = 1;
          break;
      }
      if(err) break;
      mie_vcom_com[n] = s3/sy*PI4_3;
      if(cnt_vb > 1)
      {
        fprintf(stderr,"Component#   : %2d\n",n);
        fprintf(stderr,"Min Log10(R) : %13.4e\n",mie_lmin);
        fprintf(stderr,"Max Log10(R) : %13.4e\n",mie_lmax);
        fprintf(stderr,"#Steps       : %13d\n",n_step);
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
      CREFIN = mie_refr_com[n][i]-fabs(mie_refi_com[n][i])*I;
      switch(mie_size_func[n])
      {
        case MIE_FUNC_FILE:
          switch(mie_size_ytype[n])
          {
            case MIE_YTYPE_NUMBER: // Number distribution
              for(m=0,sy=0.0; m<=n_step; m++)
              {
                r = mie_xval_com[n][m];
                x = PI2*r/(mie_wlen_um[i]);
                y = mie_yval_com[n][m];
                w = y*PI*r*r;
                if(x >= MIE_XMAX) break;
                miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,mie_angl_cos,
                       &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
                       PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);
                mie_aext_com[n][i] += QEXT*w;
                mie_asca_com[n][i] += QSCA*w;
                for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
                {
                  s1 = cabs(S1[j]);
                  s2 = cabs(S2[j]);
                  mie_phs1_com[n][k] += s1*s1*y;
                  mie_phs2_com[n][k] += s2*s2*y;
                }
                sy += y;
              }
              sw = sy;
              break;
            case MIE_YTYPE_VOLUME: // Volume distribution
              for(m=0,sw=0.0,sy=0.0; m<=n_step; m++)
              {
                r = mie_xval_com[n][m];
                x = PI2*r/(mie_wlen_um[i]);
                v = mie_yval_com[n][m];
                y = v/(PI4_3*r*r*r);
                w = y*PI*r*r;
                if(x >= MIE_XMAX) break;
                miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,mie_angl_cos,
                       &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
                       PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);
                mie_aext_com[n][i] += QEXT*w;
                mie_asca_com[n][i] += QSCA*w;
                for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
                {
                  s1 = cabs(S1[j]);
                  s2 = cabs(S2[j]);
                  mie_phs1_com[n][k] += s1*s1*y;
                  mie_phs2_com[n][k] += s2*s2*y;
                }
                sw += v;
                sy += y;
              }
              break;
            default:
              fprintf(stderr,"%s: error, mie_size_ytype[%d]=%c\n",fnam,n,mie_size_ytype[n]);
              err = 1;
              break;
          }
          break;
        case MIE_FUNC_LOGNORMAL:
        case MIE_FUNC_EFF_LOGNORMAL:
          switch(mie_size_ytype[n])
          {
            case MIE_YTYPE_NUMBER: // Number distribution
              for(m=0,sy=0.0; m<=n_step; m++)
              {
                l = mie_lmin+mie_lstp*m;
                r = pow(10.0,l);
                x = PI2*r/(mie_wlen_um[i]);
                y = norm*exp(fact*(l-mie_lmod_com[n])*(l-mie_lmod_com[n]));
                w = y*PI*r*r;
                if(x >= MIE_XMAX) break;
                miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,mie_angl_cos,
                       &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
                       PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);
                mie_aext_com[n][i] += QEXT*w;
                mie_asca_com[n][i] += QSCA*w;
                for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
                {
                  s1 = cabs(S1[j]);
                  s2 = cabs(S2[j]);
                  mie_phs1_com[n][k] += s1*s1*y;
                  mie_phs2_com[n][k] += s2*s2*y;
                }
                sy += y;
              }
              sw = sy;
              break;
            case MIE_YTYPE_VOLUME: // Volume distribution
              for(m=0,sw=0.0,sy=0.0; m<=n_step; m++)
              {
                l = mie_lmin+mie_lstp*m;
                r = pow(10.0,l);
                x = PI2*r/(mie_wlen_um[i]);
                v = norm*exp(fact*(l-mie_lmod_com[n])*(l-mie_lmod_com[n]));
                y = v/(PI4_3*r*r*r);
                w = y*PI*r*r;
                if(x >= MIE_XMAX) break;
                miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,mie_angl_cos,
                       &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
                       PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);
                mie_aext_com[n][i] += QEXT*w;
                mie_asca_com[n][i] += QSCA*w;
                for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
                {
                  s1 = cabs(S1[j]);
                  s2 = cabs(S2[j]);
                  mie_phs1_com[n][k] += s1*s1*y;
                  mie_phs2_com[n][k] += s2*s2*y;
                }
                sw += v;
                sy += y;
              }
              break;
            default:
              fprintf(stderr,"%s: error, mie_size_ytype[%d]=%c\n",fnam,n,mie_size_ytype[n]);
              err = 1;
              break;
          }
          break;
        case MIE_FUNC_GAMMA:
          switch(mie_size_ytype[n])
          {
            case MIE_YTYPE_NUMBER: // Number distribution
              for(m=0,sy=0.0; m<=n_step; m++)
              {
                l = mie_lmin+mie_lstp*m;
                r = pow(10.0,l);
                x = PI2*r/(mie_wlen_um[i]);
                y = norm*pow(r,indx)*exp(fact*r);
                w = y*PI*r*r;
                if(x >= MIE_XMAX) break;
                miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,mie_angl_cos,
                       &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
                       PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);
                mie_aext_com[n][i] += QEXT*w;
                mie_asca_com[n][i] += QSCA*w;
                for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
                {
                  s1 = cabs(S1[j]);
                  s2 = cabs(S2[j]);
                  mie_phs1_com[n][k] += s1*s1*y;
                  mie_phs2_com[n][k] += s2*s2*y;
                }
                sy += y;
              }
              sw = sy;
              break;
            case MIE_YTYPE_VOLUME: // Volume distribution
              for(m=0,sw=0.0,sy=0.0; m<=n_step; m++)
              {
                l = mie_lmin+mie_lstp*m;
                r = pow(10.0,l);
                x = PI2*r/(mie_wlen_um[i]);
                v = norm*pow(r,indx)*exp(fact*r);
                y = v/(PI4_3*r*r*r);
                w = y*PI*r*r;
                if(x >= MIE_XMAX) break;
                miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,mie_angl_cos,
                       &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
                       PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);
                mie_aext_com[n][i] += QEXT*w;
                mie_asca_com[n][i] += QSCA*w;
                for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
                {
                  s1 = cabs(S1[j]);
                  s2 = cabs(S2[j]);
                  mie_phs1_com[n][k] += s1*s1*y;
                  mie_phs2_com[n][k] += s2*s2*y;
                }
                sw += v;
                sy += y;
              }
              break;
            default:
              fprintf(stderr,"%s: error, mie_size_ytype[%d]=%c\n",fnam,n,mie_size_ytype[n]);
              err = 1;
              break;
          }
          break;
        default:
          fprintf(stderr,"%s: error, mie_size_func[%d]=%c\n",fnam,n,mie_size_func[n]);
          err = 1;
          break;
      }
      if(err) break;
      if(fabs(sw-1.0) > epsilon)
      {
        fprintf(stderr,"Warning, sw=%13.4e\n",sw);
      }
      sinv = 1.0/sy; // average per particle
      mie_aext_com[n][i] *= sinv;
      mie_asca_com[n][i] *= sinv;
      for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
      {
        mie_phas_com[n][k] = 0.5*(mie_phs1_com[n][k]+mie_phs2_com[n][k]);
      }
      sp = 0.0;
      for(j=1,k=mie_n_angl*i+1; j<mie_n_angl; j++,k++)
      {
        sp += (mie_phas_com[n][k-1]*mie_angl_sin[j-1]+mie_phas_com[n][k]*mie_angl_sin[j])*mie_angl_dif[j];
      }
      sinv = 1.0/sp; // normalized for natural light
      for(j=0,k=mie_n_angl*i; j<mie_n_angl; j++,k++)
      {
        mie_phs1_com[n][k] *= sinv;
        mie_phs2_com[n][k] *= sinv;
        mie_phas_com[n][k] *= sinv;
      }
      for(j=1,k=mie_n_angl*i+1; j<mie_n_angl; j++,k++)
      {
        mie_asym_com[n][i] += (mie_phas_com[n][k-1]*mie_angl_sin[j-1]*mie_angl_cos[j-1]+
                               mie_phas_com[n][k]*mie_angl_sin[j]*mie_angl_cos[j])*mie_angl_dif[j];
      }
    }
    if(err) break;
  }

  free(S1);
  free(S2);
  if(err)
  {
    return -1;
  }

  return 0;
}

int Init(void)
{
  int i,n;

  if(CommonInit() < 0)
  {
    return -1;
  }
  for(n=0; n<mie_n_comp; n++)
  {
    if(mie_size_func[n] == MIE_FUNC_FILE)
    {
      // dy/dx * dx -> dy
      mie_yval_com[n][0] *= (mie_xval_com[n][1]-mie_xval_com[n][0]);
      for(i=1; i<mie_n_size[n]-1; i++)
      {
        mie_yval_com[n][i] *= 0.5*(mie_xval_com[n][i+1]-mie_xval_com[n][i-1]);
      }
      i = mie_n_size[n]-1;
      mie_yval_com[n][i] *= (mie_xval_com[n][i]-mie_xval_com[n][i-1]);
      // Log10(R) -> R
      switch(mie_size_xtype[n])
      {
        case MIE_XTYPE_LOGR:
          for(i=0; i<mie_n_size[n]; i++)
          {
            mie_xval_com[n][i] = pow(10.0,mie_xval_com[n][i]);
          }
          break;
        case MIE_XTYPE_LNR:
          for(i=0; i<mie_n_size[n]; i++)
          {
            mie_xval_com[n][i] = exp(mie_xval_com[n][i]);
          }
          break;
        default:
          break;
      }
    }
    else
    {
      mie_lmod_com[n] = log10(mie_rmod_com[n]);
    }
  }
  n = mie_n_angl/2;
  if(mie_n_angl%2==1 && fabs(mie_angl[n]-90.0)>EPSILON)
  {
    cnt_anyang = MIE_TRUE;
  }
  else
  {
    for(i=0; i<n; i++)
    {
      if(fabs(mie_angl[mie_n_angl-1-i]+mie_angl[i]-180.0) > EPSILON)
      {
        cnt_anyang = MIE_TRUE;
        break;
      }
      if(mie_angl[i+1] <= mie_angl[i])
      {
        cnt_anyang = MIE_TRUE;
        break;
      }
    }
  }
  if(cnt_vb > 2)
  {
    fprintf(stderr,"cnt_anyang = %d\n",cnt_anyang);
  }

  return 0;
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
    {"mie_wsgm",1,0,'S'},
  };

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-x:X:n:s:S:",long_options,&option_index);
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
      case 'S':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_wsgm = xtmp;
        else
        {
          fprintf(stderr,"Log10(R) width -> out of range %s\n",optarg);
          rt = -1;
        }
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
  char optstring[] = "fgwlLaWAPQiIruxXnsSmoOdvh";

  CommonUsage("aeros_miev0",0);
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
      case 'S':
        fprintf(stderr," S -mie_wsgm      |%s|%s|%s| %f\n",As(e,"Log10(R) width",n), As(a,"sigma",n),      Af(d,MIE_WSGM,n),mie_wsgm);
        break;
      default:
        CommonUsage(NULL,optstring[i]);
        break;
    }
  }
  CommonUsage(NULL,1);
  CommonUsage(NULL,2);
  CommonUsage(NULL,3);
  CommonUsage(NULL,4);

  return 0;
}
