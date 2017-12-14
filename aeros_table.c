/***********************************************************/
/* AEROS_TABLE ... Aerosol generator with Table            */
/* Author: N.Manago (Created on Wed, 03 Sep 2014)          */
/* $Revision: 1104 $                                       */
/* $Date: 2015-09-02 00:16:06 +0900 (Wed, 02 Sep 2015) $   */
/***********************************************************/
#include "aeros_common.c"
// Common constants
#define	DELTA				1.0e-13			// A small number
#define	NODATA				0			// No data
#define	INVALID				-1			// Invalid value
#define	M_1_2SQRTPI			0.28209479177387814	// 1/(2sqrt(pi))
// Constants for table
#define	TAB_DATDIR			"data"
#define	TAB_N_ANGL			721
#define	TAB_X_LOG_MIN			-4.0
#define	TAB_X_LOG_MAX			4.0
#define	TAB_X_LOG_STP			0.001
#define	TAB_REFR_NDIG			8
#define	TAB_REFR_MULT			1.0e6
#define	TAB_REFR_MIN			1.3
#define	TAB_REFR_MAX			1.8
#define	TAB_REFR_STP			0.01
#define	TAB_REFI_NDIG			8
#define	TAB_REFI_MULT			1.0e6
#define	TAB_REFI_LOG_MIN		-5.0
#define	TAB_REFI_LOG_MAX		0.0
#define	TAB_REFI_LOG_STP		0.1
#define	TAB_E_NDIG 			6
#define	TAB_E_MULT 			1.0e5
#define	TAB_E_MIN			0.3
#define	TAB_E_MAX			3.0
#define	TAB_E_STP			0.1
// Constants for Mie calculation
#define	MIE_TRUE			1			// True  value
#define	MIE_FALSE			0			// False value
#define	MIE_XMAX			2.0e4			// Max size parameter
#define	MIE_RMIN			NAN			// R min in um
#define	MIE_RMAX			NAN			// R max in um
#define	MIE_LMIN			-4.0			// Min Log10(R)
#define	MIE_LMAX			+5.0			// Max Log10(R)
#define	MIE_LSTP			1.0e-3			// Log10(R) step
#define	MIE_WSGM			5.0			// Log10(R) width in sigma
// Constants for control
#define	CNT_CONF			NONAME
#define	CNT_MAXCMNT			30
#define	INTERP_NEAREST			0
#define	INTERP_LINEAR			1
#define	INTERP_SPLINE			2

char	tab_datdir[MAXLINE]		= TAB_DATDIR;
int	tab_n_angl			= TAB_N_ANGL;
int	tab_x_size			= 0;
double	tab_x_log_min			= TAB_X_LOG_MIN;
double	tab_x_log_max			= TAB_X_LOG_MAX;
double	tab_x_log_stp			= TAB_X_LOG_STP;
double	tab_x_log_stp_inv;
int	tab_refr_size			= 0;
int	tab_refr_ndig			= TAB_REFR_NDIG;
double	tab_refr_mult			= TAB_REFR_MULT;
double	tab_refr_min			= TAB_REFR_MIN;
double	tab_refr_max			= TAB_REFR_MAX;
double	tab_refr_stp			= TAB_REFR_STP;
double	tab_refr_stp_inv;
int	tab_refi_size			= 0;
int	tab_refi_ndig			= TAB_REFI_NDIG;
double	tab_refi_mult			= TAB_REFI_MULT;
double	tab_refi_log_min		= TAB_REFI_LOG_MIN;
double	tab_refi_log_max		= TAB_REFI_LOG_MAX;
double	tab_refi_log_stp		= TAB_REFI_LOG_STP;
double	tab_refi_log_stp_inv;
int	tab_e_size			= 0;
int	tab_e_ndig			= TAB_E_NDIG;
double	tab_e_mult			= TAB_E_MULT;
double	tab_e_min			= TAB_E_MIN;
double	tab_e_max			= TAB_E_MAX;
double	tab_e_stp			= TAB_E_STP;
double	tab_e_stp_inv;
double	*tab_x_val			= NULL;
double	*tab_x_log_val			= NULL;
double	*tab_refr_val			= NULL;
double	*tab_refi_log_val		= NULL;
double	*tab_e_val			= NULL;
double	*qval_init			= NULL;
double	*qval_refr			= NULL;
double	*qval_refi			= NULL;
double	*qval_e				= NULL;
double	*pval_init			= NULL;
double	*pval_refr			= NULL;
double	*pval_refi			= NULL;
double	*pval_e				= NULL;
double	*mie_aext_tab[MIE_MAXCOMP];				// Extinction coefficient
double	*mie_asca_tab[MIE_MAXCOMP];				// Scattering coefficient
double	*mie_phs1_tab[MIE_MAXCOMP];				// Phase function
double	*mie_phs2_tab[MIE_MAXCOMP];				// Phase function
double	*mie_aext_int[MIE_MAXCOMP];				// Extinction coefficient
double	*mie_asca_int[MIE_MAXCOMP];				// Scattering coefficient
double	*mie_phs1_int[MIE_MAXCOMP];				// Phase function
double	*mie_phs2_int[MIE_MAXCOMP];				// Phase function
double	*mie_lval_com[MIE_MAXCOMP];				// Log10(R)
double	*mie_refi_log_com[MIE_MAXCOMP];				// Refractive index (imaginary)
double	mie_rmin			= MIE_RMIN;		// R min in um
double	mie_rmax			= MIE_RMAX;		// R max in um
double	mie_lstp			= MIE_LSTP;		// Log10(R) step
double	mie_wsgm			= MIE_WSGM;		// Log10(R) width in sigma
double	mie_lmin;
double	mie_lmax;
double	mie_lmod_com[MIE_MAXCOMP];				// Log10(Mode radius in um)
double	*mie_rmax_com[MIE_MAXCOMP];				// R max in um
int	cnt_tmx				= 0;			// T-matrix mode
int	cnt_n_cmnt			= NODATA;
char	cnt_cmnt[CNT_MAXCMNT][MAXLINE];				// Comments
char	cnt_conf[MAXLINE]		= CNT_CONF;		// Configuration file

int MieTable(int n);
int TmxTable(int n);
int IntTable(void);
int ReadConfig(void);

int MieTable(int n)
{
  int i,j,k,p;
  int refr_skip;
  int refi_skip;
  int flag;
  int ndat;
  int refr_indx_1,refr_indx_2;
  int refi_indx_1,refi_indx_2;
  int nr_1,nr_2;
  int ni_1,ni_2;
  double refr_1,refr_2,refr_rat;
  double refi_1,refi_2;
  double refi_log_1,refi_log_2,refi_rat;
  double init_rat;
  double fact1,fact2;
  double lmin,lmax,lstp;
  const char tq_finp[2][MAXLINE] = {"_qext","_qsca"};
  const char tp_finp[2][MAXLINE] = {"_p1","_p2"};
  int tq_size = tab_x_size;
  int tp_size = tab_x_size*mie_n_angl;
  double *tq_para[2];
  double *tp_para[2];
  char fmt[MAXLINE];
  char finp[MAXLINE];
  char line[MAXLINE];
  FILE *fp;
  const char fnam[] = "MieTable";

  // Check the input range
  for(i=0; i<mie_n_wlen; i++)
  {
    if(mie_refr_com[n][i]<tab_refr_min-EPSILON ||
       mie_refr_com[n][i]>tab_refr_max+EPSILON)
    {
      fprintf(stderr,"%s: error, mie_refr_com[%d][%d]=%13.6e, out of range.\n",
                      fnam,n,i,mie_refr_com[n][i]);
      return -1;
    }
    if(mie_refi_log_com[n][i]<tab_refi_log_min-EPSILON ||
       mie_refi_log_com[n][i]>tab_refi_log_max+EPSILON)
    {
      fprintf(stderr,"%s: error, mie_refi_com[%d][%d]=%13.6e, out of range.\n",
                      fnam,n,i,mie_refi_com[n][i]);
      return -1;
    }
  }

  snprintf(fmt,MAXLINE,"%s/mie_refr%%0%dd_refi%%0%dd%%s.dat",tab_datdir,tab_refr_ndig,tab_refi_ndig);

  tq_para[0] = mie_aext_tab[n];
  tq_para[1] = mie_asca_tab[n];
  tp_para[0] = mie_phs1_tab[n];
  tp_para[1] = mie_phs2_tab[n];
  for(i=0; i<mie_n_wlen; i++)
  {
    fact1 = 0.5*M_1_PI*mie_wlen_um[i];
    fact2 = 0.25*M_1_PI*mie_wlen_um[i]*mie_wlen_um[i];
    init_rat = 1.0;
    // refr
    refr_indx_1 = (int)((mie_refr_com[n][i]-tab_refr_min)*tab_refr_stp_inv);
    refr_indx_2 = refr_indx_1+1;
    if(refr_indx_1 < 0)
    {
      refr_indx_1 = 0;
      refr_indx_2 = 1;
    } else
    if(refr_indx_2 >= tab_refr_size)
    {
      refr_indx_2 = tab_refr_size-1;
      refr_indx_1 = refr_indx_2-1;
    }
    refr_1 = tab_refr_val[refr_indx_1];
    refr_2 = tab_refr_val[refr_indx_2];
    refr_rat = (mie_refr_com[n][i]-refr_1)*tab_refr_stp_inv;
    if(fabs(refr_rat) < DELTA)
    {
      nr_1 = (int)(refr_1*tab_refr_mult+0.5);
      refr_skip = 1;
    } else
    if(fabs(1.0-refr_rat) < DELTA)
    {
      refr_1 = refr_2;
      nr_1 = (int)(refr_1*tab_refr_mult+0.5);
      refr_skip = 1;
    }
    else
    {
      nr_1 = (int)(refr_1*tab_refr_mult+0.5);
      nr_2 = (int)(refr_2*tab_refr_mult+0.5);
      init_rat -= refr_rat;
      refr_skip = 0;
    }
    // refi
    refi_indx_1 = (int)((mie_refi_log_com[n][i]-tab_refi_log_min)*tab_refi_log_stp_inv);
    refi_indx_2 = refi_indx_1+1;
    if(refi_indx_1 < 0)
    {
      refi_indx_1 = 0;
      refi_indx_2 = 1;
    } else
    if(refi_indx_2 >= tab_refi_size)
    {
      refi_indx_2 = tab_refi_size-1;
      refi_indx_1 = refi_indx_2-1;
    }
    refi_log_1 = tab_refi_log_val[refi_indx_1];
    refi_log_2 = tab_refi_log_val[refi_indx_2];
    refi_rat = (mie_refi_log_com[n][i]-refi_log_1)*tab_refi_log_stp_inv;
    if(fabs(refi_rat) < DELTA)
    {
      refi_1 = pow(10.0,refi_log_1);
      ni_1 = (int)(refi_1*tab_refi_mult+0.5);
      refi_skip = 1;
    } else
    if(fabs(1.0-refi_rat) < DELTA)
    {
      refi_log_1 = refi_log_2;
      refi_1 = pow(10.0,refi_log_1);
      ni_1 = (int)(refi_1*tab_refi_mult+0.5);
      refi_skip = 1;
    }
    else
    {
      refi_1 = pow(10.0,refi_log_1);
      refi_2 = pow(10.0,refi_log_2);
      ni_1 = (int)(refi_1*tab_refi_mult+0.5);
      ni_2 = (int)(refi_2*tab_refi_mult+0.5);
      init_rat -= refi_rat;
      refi_skip = 0;
    }
    flag = 0x0001*refr_skip+0x0010*refi_skip;
    // rmax
    snprintf(finp,MAXLINE,fmt,nr_1,ni_1,"");
    if((fp=fopen(finp,"r")) != NULL)
    {
      if(fgets(line,MAXLINE,fp) != NULL)
      {
        if(sscanf(line,"%lf%lf%lf",&lmin,&lmax,&lstp) == 3)
        {
          mie_rmax_com[n][i] = MIN(mie_rmax_com[n][i],pow(10.0,lmax)*fact1);
        }
      }
      fclose(fp);
    }
    // refr
    if(refr_skip == 0)
    {
      snprintf(finp,MAXLINE,fmt,nr_2,ni_1,"");
      if((fp=fopen(finp,"r")) != NULL)
      {
        if(fgets(line,MAXLINE,fp) != NULL)
        {
          if(sscanf(line,"%lf%lf%lf",&lmin,&lmax,&lstp) == 3)
          {
            mie_rmax_com[n][i] = MIN(mie_rmax_com[n][i],pow(10.0,lmax)*fact1);
          }
        }
        fclose(fp);
      }
    }
    // refi
    if(refi_skip == 0)
    {
      snprintf(finp,MAXLINE,fmt,nr_1,ni_2,"");
      if((fp=fopen(finp,"r")) != NULL)
      {
        if(fgets(line,MAXLINE,fp) != NULL)
        {
          if(sscanf(line,"%lf%lf%lf",&lmin,&lmax,&lstp) == 3)
          {
            mie_rmax_com[n][i] = MIN(mie_rmax_com[n][i],pow(10.0,lmax)*fact1);
          }
        }
        fclose(fp);
      }
    }
    if(cnt_vb > 1)
    {
      fprintf(stderr,"mie_rmax_com[%d][%d]: %13.6e\n",n,i,mie_rmax_com[n][i]);
    }
    // aext, asca
    for(k=0; k<2; k++)
    {
      // init
      snprintf(finp,MAXLINE,fmt,nr_1,ni_1,tq_finp[k]);
      if((fp=fopen(finp,"r")) == NULL)
      {
        fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
        return -1;
      }
      if((ndat=fread(qval_init,sizeof(double),tq_size+1,fp)) != tq_size)
      {
        fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                        ndat,tq_size,finp);
        return -1;
      }
      fclose(fp);
      // refr
      if(refr_skip == 0)
      {
        snprintf(finp,MAXLINE,fmt,nr_2,ni_1,tq_finp[k]);
        if((fp=fopen(finp,"r")) == NULL)
        {
          fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
          return -1;
        }
        if((ndat=fread(qval_refr,sizeof(double),tq_size+1,fp)) != tq_size)
        {
          fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                          ndat,tq_size,finp);
          return -1;
        }
        fclose(fp);
      }
      // refi
      if(refi_skip == 0)
      {
        snprintf(finp,MAXLINE,fmt,nr_1,ni_2,tq_finp[k]);
        if((fp=fopen(finp,"r")) == NULL)
        {
          fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
          return -1;
        }
        if((ndat=fread(qval_refi,sizeof(double),tq_size+1,fp)) != tq_size)
        {
          fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                          ndat,tq_size,finp);
          return -1;
        }
        fclose(fp);
      }
      switch(flag)
      {
        case 0x0011: // refi_skip==1 && refr_skip==1
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = qval_init[j]*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
        case 0x0001: // refi_skip==0 && refr_skip==1
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = (init_rat*qval_init[j]+refi_rat*qval_refi[j])*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
        case 0x0010: // refi_skip==1 && refr_skip==0
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = (init_rat*qval_init[j]+refr_rat*qval_refr[j])*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
        default: // refi_skip==0 && refr_skip==0
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = (init_rat*qval_init[j]+refr_rat*qval_refr[j]+refi_rat*qval_refi[j])*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
      }
    }
    // phs1, phs2
    for(k=0; k<2; k++)
    {
      // init
      snprintf(finp,MAXLINE,fmt,nr_1,ni_1,tp_finp[k]);
      if((fp=fopen(finp,"r")) == NULL)
      {
        fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
        return -1;
      }
      if((ndat=fread(pval_init,sizeof(double),tp_size+1,fp)) != tp_size)
      {
        fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                        ndat,tp_size,finp);
        return -1;
      }
      fclose(fp);
      // refr
      if(refr_skip == 0)
      {
        snprintf(finp,MAXLINE,fmt,nr_2,ni_1,tp_finp[k]);
        if((fp=fopen(finp,"r")) == NULL)
        {
          fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
          return -1;
        }
        if((ndat=fread(pval_refr,sizeof(double),tp_size+1,fp)) != tp_size)
        {
          fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                          ndat,tp_size,finp);
          return -1;
        }
        fclose(fp);
      }
      // refi
      if(refi_skip == 0)
      {
        snprintf(finp,MAXLINE,fmt,nr_1,ni_2,tp_finp[k]);
        if((fp=fopen(finp,"r")) == NULL)
        {
          fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
          return -1;
        }
        if((ndat=fread(pval_refi,sizeof(double),tp_size+1,fp)) != tp_size)
        {
          fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                          ndat,tp_size,finp);
          return -1;
        }
        fclose(fp);
      }
      switch(flag)
      {
        case 0x0011: // refi_skip==1 && refr_skip==1
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = pval_init[j];
          }
          break;
        case 0x0001: // refi_skip==0 && refr_skip==1
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = init_rat*pval_init[j]+refi_rat*pval_refi[j];
          }
          break;
        case 0x0010: // refi_skip==1 && refr_skip==0
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = init_rat*pval_init[j]+refr_rat*pval_refr[j];
          }
          break;
        default: // refi_skip==0 && refr_skip==0
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = init_rat*pval_init[j]+refr_rat*pval_refr[j]+refi_rat*pval_refi[j];
          }
          break;
      }
    }
  }

  return 0;
}

int TmxTable(int n)
{
  int i,j,k,p;
  int refr_skip;
  int refi_skip;
  int e_skip;
  int flag;
  int ndat;
  int refr_indx_1,refr_indx_2;
  int refi_indx_1,refi_indx_2;
  int e_indx_1,e_indx_2;
  int nr_1,nr_2;
  int ni_1,ni_2;
  int ne_1,ne_2;
  double refr_1,refr_2,refr_rat;
  double refi_1,refi_2;
  double refi_log_1,refi_log_2,refi_rat;
  double e_1,e_2,e_rat;
  double init_rat;
  double fact1,fact2;
  double lmin,lmax,lstp;
  const char tq_finp[2][MAXLINE] = {"_qext","_qsca"};
  const char tp_finp[2][MAXLINE] = {"_p1","_p2"};
  const int tq_size = tab_x_size;
  const int tp_size = tab_x_size*mie_n_angl;
  double *tq_para[2];
  double *tp_para[2];
  char fmt[MAXLINE];
  char finp[MAXLINE];
  char line[MAXLINE];
  FILE *fp;
  const char fnam[] = "TmxTable";

  // Check the input range
  for(i=0; i<mie_n_wlen; i++)
  {
    if(mie_refr_com[n][i]<tab_refr_min-EPSILON ||
       mie_refr_com[n][i]>tab_refr_max+EPSILON)
    {
      fprintf(stderr,"%s: error, mie_refr_com[%d][%d]=%13.6e, out of range.\n",
                      fnam,n,i,mie_refr_com[n][i]);
      return -1;
    }
    if(mie_refi_log_com[n][i]<tab_refi_log_min-EPSILON ||
       mie_refi_log_com[n][i]>tab_refi_log_max+EPSILON)
    {
      fprintf(stderr,"%s: error, mie_refi_com[%d][%d]=%13.6e, out of range.\n",
                      fnam,n,i,mie_refi_com[n][i]);
      return -1;
    }
  }

  snprintf(fmt,MAXLINE,"%s/tmq_e%%0%dd_refr%%0%dd_refi%%0%dd%%s.dat",tab_datdir,tab_e_ndig,tab_refr_ndig,tab_refi_ndig);

  tq_para[0] = mie_aext_tab[n];
  tq_para[1] = mie_asca_tab[n];
  tp_para[0] = mie_phs1_tab[n];
  tp_para[1] = mie_phs2_tab[n];
  for(i=0; i<mie_n_wlen; i++)
  {
    fact1 = 0.5*M_1_PI*mie_wlen_um[i];
    fact2 = 0.25*M_1_PI*mie_wlen_um[i]*mie_wlen_um[i];
    init_rat = 1.0;
    // refr
    refr_indx_1 = (int)((mie_refr_com[n][i]-tab_refr_min)*tab_refr_stp_inv);
    refr_indx_2 = refr_indx_1+1;
    if(refr_indx_1 < 0)
    {
      refr_indx_1 = 0;
      refr_indx_2 = 1;
    } else
    if(refr_indx_2 >= tab_refr_size)
    {
      refr_indx_2 = tab_refr_size-1;
      refr_indx_1 = refr_indx_2-1;
    }
    refr_1 = tab_refr_val[refr_indx_1];
    refr_2 = tab_refr_val[refr_indx_2];
    refr_rat = (mie_refr_com[n][i]-refr_1)*tab_refr_stp_inv;
    if(fabs(refr_rat) < DELTA)
    {
      nr_1 = (int)(refr_1*tab_refr_mult+0.5);
      refr_skip = 1;
    } else
    if(fabs(1.0-refr_rat) < DELTA)
    {
      refr_1 = refr_2;
      nr_1 = (int)(refr_1*tab_refr_mult+0.5);
      refr_skip = 1;
    }
    else
    {
      nr_1 = (int)(refr_1*tab_refr_mult+0.5);
      nr_2 = (int)(refr_2*tab_refr_mult+0.5);
      init_rat -= refr_rat;
      refr_skip = 0;
    }
    // refi
    refi_indx_1 = (int)((mie_refi_log_com[n][i]-tab_refi_log_min)*tab_refi_log_stp_inv);
    refi_indx_2 = refi_indx_1+1;
    if(refi_indx_1 < 0)
    {
      refi_indx_1 = 0;
      refi_indx_2 = 1;
    } else
    if(refi_indx_2 >= tab_refi_size)
    {
      refi_indx_2 = tab_refi_size-1;
      refi_indx_1 = refi_indx_2-1;
    }
    refi_log_1 = tab_refi_log_val[refi_indx_1];
    refi_log_2 = tab_refi_log_val[refi_indx_2];
    refi_rat = (mie_refi_log_com[n][i]-refi_log_1)*tab_refi_log_stp_inv;
    if(fabs(refi_rat) < DELTA)
    {
      refi_1 = pow(10.0,refi_log_1);
      ni_1 = (int)(refi_1*tab_refi_mult+0.5);
      refi_skip = 1;
    } else
    if(fabs(1.0-refi_rat) < DELTA)
    {
      refi_log_1 = refi_log_2;
      refi_1 = pow(10.0,refi_log_1);
      ni_1 = (int)(refi_1*tab_refi_mult+0.5);
      refi_skip = 1;
    }
    else
    {
      refi_1 = pow(10.0,refi_log_1);
      refi_2 = pow(10.0,refi_log_2);
      ni_1 = (int)(refi_1*tab_refi_mult+0.5);
      ni_2 = (int)(refi_2*tab_refi_mult+0.5);
      init_rat -= refi_rat;
      refi_skip = 0;
    }
    // e
    e_indx_1 = (int)((mie_xeps_com[n][0]-tab_e_min)*tab_e_stp_inv);
    e_indx_2 = e_indx_1+1;
    if(e_indx_2 >= tab_e_size)
    {
      e_indx_2 = tab_e_size-1;
      e_indx_1 = e_indx_2-1;
    }
    if(e_indx_1 < 0) // arrow tab_e_size == 1
    {
      e_indx_1 = 0;
      e_indx_2 = 1;
    }
    e_1 = tab_e_val[e_indx_1];
    e_2 = tab_e_val[e_indx_2];
    e_rat = (mie_xeps_com[n][0]-e_1)*tab_e_stp_inv;
    if(fabs(e_rat) < DELTA)
    {
      ne_1 = (int)(e_1*tab_e_mult+0.5);
      e_skip = 1;
    } else
    if(fabs(1.0-e_rat) < DELTA)
    {
      e_1 = e_2;
      ne_1 = (int)(e_1*tab_e_mult+0.5);
      e_skip = 1;
    }
    else
    {
      ne_1 = (int)(e_1*tab_e_mult+0.5);
      ne_2 = (int)(e_2*tab_e_mult+0.5);
      init_rat -= e_rat;
      e_skip = 0;
    }
    flag = 0x0001*refr_skip+0x0010*refi_skip+0x0100*e_skip;
    // rmax
    snprintf(finp,MAXLINE,fmt,ne_1,nr_1,ni_1,"");
    if((fp=fopen(finp,"r")) != NULL)
    {
      if(fgets(line,MAXLINE,fp) != NULL)
      {
        if(sscanf(line,"%lf%lf%lf",&lmin,&lmax,&lstp) == 3)
        {
          mie_rmax_com[n][i] = MIN(mie_rmax_com[n][i],pow(10.0,lmax)*fact1);
        }
      }
      fclose(fp);
    }
    // refr
    if(refr_skip == 0)
    {
      snprintf(finp,MAXLINE,fmt,ne_1,nr_2,ni_1,"");
      if((fp=fopen(finp,"r")) != NULL)
      {
        if(fgets(line,MAXLINE,fp) != NULL)
        {
          if(sscanf(line,"%lf%lf%lf",&lmin,&lmax,&lstp) == 3)
          {
            mie_rmax_com[n][i] = MIN(mie_rmax_com[n][i],pow(10.0,lmax)*fact1);
          }
        }
        fclose(fp);
      }
    }
    // refi
    if(refi_skip == 0)
    {
      snprintf(finp,MAXLINE,fmt,ne_1,nr_1,ni_2,"");
      if((fp=fopen(finp,"r")) != NULL)
      {
        if(fgets(line,MAXLINE,fp) != NULL)
        {
          if(sscanf(line,"%lf%lf%lf",&lmin,&lmax,&lstp) == 3)
          {
            mie_rmax_com[n][i] = MIN(mie_rmax_com[n][i],pow(10.0,lmax)*fact1);
          }
        }
        fclose(fp);
      }
    }
    // e
    if(e_skip == 0)
    {
      snprintf(finp,MAXLINE,fmt,ne_2,nr_1,ni_1,"");
      if((fp=fopen(finp,"r")) != NULL)
      {
        if(fgets(line,MAXLINE,fp) != NULL)
        {
          if(sscanf(line,"%lf%lf%lf",&lmin,&lmax,&lstp) == 3)
          {
            mie_rmax_com[n][i] = MIN(mie_rmax_com[n][i],pow(10.0,lmax)*fact1);
          }
        }
        fclose(fp);
      }
    }
    if(cnt_vb > 1)
    {
      fprintf(stderr,"mie_rmax_com[%d][%d]: %13.6e\n",n,i,mie_rmax_com[n][i]);
    }
    // aext, asca
    for(k=0; k<2; k++)
    {
      // init
      snprintf(finp,MAXLINE,fmt,ne_1,nr_1,ni_1,tq_finp[k]);
      if((fp=fopen(finp,"r")) == NULL)
      {
        fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
        return -1;
      }
      if((ndat=fread(qval_init,sizeof(double),tq_size+1,fp)) != tq_size)
      {
        fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                        ndat,tq_size,finp);
        return -1;
      }
      fclose(fp);
      // refr
      if(refr_skip == 0)
      {
        snprintf(finp,MAXLINE,fmt,ne_1,nr_2,ni_1,tq_finp[k]);
        if((fp=fopen(finp,"r")) == NULL)
        {
          fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
          return -1;
        }
        if((ndat=fread(qval_refr,sizeof(double),tq_size+1,fp)) != tq_size)
        {
          fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                          ndat,tq_size,finp);
          return -1;
        }
        fclose(fp);
      }
      // refi
      if(refi_skip == 0)
      {
        snprintf(finp,MAXLINE,fmt,ne_1,nr_1,ni_2,tq_finp[k]);
        if((fp=fopen(finp,"r")) == NULL)
        {
          fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
          return -1;
        }
        if((ndat=fread(qval_refi,sizeof(double),tq_size+1,fp)) != tq_size)
        {
          fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                          ndat,tq_size,finp);
          return -1;
        }
        fclose(fp);
      }
      // e
      if(e_skip == 0)
      {
        snprintf(finp,MAXLINE,fmt,ne_2,nr_1,ni_1,tq_finp[k]);
        if((fp=fopen(finp,"r")) == NULL)
        {
          fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
          return -1;
        }
        if((ndat=fread(qval_e,sizeof(double),tq_size+1,fp)) != tq_size)
        {
          fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                          ndat,tq_size,finp);
          return -1;
        }
        fclose(fp);
      }
      switch(flag)
      {
        case 0x0111: // e_skip==1 && refi_skip==1 && refr_skip==1
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = qval_init[j]*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
        case 0x0011: // e_skip==0 && refi_skip==1 && refr_skip==1
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = (init_rat*qval_init[j]+e_rat*qval_e[j])*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
        case 0x0101: // e_skip==1 && refi_skip==0 && refr_skip==1
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = (init_rat*qval_init[j]+refi_rat*qval_refi[j])*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
        case 0x0110: // e_skip==1 && refi_skip==1 && refr_skip==0
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = (init_rat*qval_init[j]+refr_rat*qval_refr[j])*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
        case 0x0001: // e_skip==0 && refi_skip==0 && refr_skip==1
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = (init_rat*qval_init[j]+refi_rat*qval_refi[j]+
                                                   e_rat*qval_e[j])*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
        case 0x0010: // e_skip==0 && refi_skip==1 && refr_skip==0
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = (init_rat*qval_init[j]+refr_rat*qval_refr[j]+
                                                   e_rat*qval_e[j])*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
        case 0x0100: // e_skip==1 && refi_skip==0 && refr_skip==0
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = (init_rat*qval_init[j]+refr_rat*qval_refr[j]+
                                                   refi_rat*qval_refi[j])*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
        default: // e_skip==0 && refi_skip==0 && refr_skip==0
          for(j=0,p=tq_size*i; j<tq_size; j++,p++)
          {
            tq_para[k][p] = (init_rat*qval_init[j]+refr_rat*qval_refr[j]+
                                                   refi_rat*qval_refi[j]+
                                                   e_rat*qval_e[j])*fact2*tab_x_val[j]*tab_x_val[j]; // Q*pi*r*r
          }
          break;
      }
    }
    // phs1, phs2
    for(k=0; k<2; k++)
    {
      // init
      snprintf(finp,MAXLINE,fmt,ne_1,nr_1,ni_1,tp_finp[k]);
      if((fp=fopen(finp,"r")) == NULL)
      {
        fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
        return -1;
      }
      if((ndat=fread(pval_init,sizeof(double),tp_size+1,fp)) != tp_size)
      {
        fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                        ndat,tp_size,finp);
        return -1;
      }
      fclose(fp);
      // refr
      if(refr_skip == 0)
      {
        snprintf(finp,MAXLINE,fmt,ne_1,nr_2,ni_1,tp_finp[k]);
        if((fp=fopen(finp,"r")) == NULL)
        {
          fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
          return -1;
        }
        if((ndat=fread(pval_refr,sizeof(double),tp_size+1,fp)) != tp_size)
        {
          fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                          ndat,tp_size,finp);
          return -1;
        }
        fclose(fp);
      }
      // refi
      if(refi_skip == 0)
      {
        snprintf(finp,MAXLINE,fmt,ne_1,nr_1,ni_2,tp_finp[k]);
        if((fp=fopen(finp,"r")) == NULL)
        {
          fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
          return -1;
        }
        if((ndat=fread(pval_refi,sizeof(double),tp_size+1,fp)) != tp_size)
        {
          fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                          ndat,tp_size,finp);
          return -1;
        }
        fclose(fp);
      }
      // e
      if(e_skip == 0)
      {
        snprintf(finp,MAXLINE,fmt,ne_2,nr_1,ni_1,tp_finp[k]);
        if((fp=fopen(finp,"r")) == NULL)
        {
          fprintf(stderr,"%s: error, cannot open %s\n",fnam,finp);
          return -1;
        }
        if((ndat=fread(pval_e,sizeof(double),tp_size+1,fp)) != tp_size)
        {
          fprintf(stderr,"%s: error, %d data have been read. Expected %d >>> %s\n",fnam,
                          ndat,tp_size,finp);
          return -1;
        }
        fclose(fp);
      }
      switch(flag)
      {
        case 0x0111: // e_skip==1 && refi_skip==1 && refr_skip==1
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = pval_init[j];
          }
          break;
        case 0x0011: // e_skip==0 && refi_skip==1 && refr_skip==1
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = init_rat*pval_init[j]+e_rat*pval_e[j];
          }
          break;
        case 0x0101: // e_skip==1 && refi_skip==0 && refr_skip==1
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = init_rat*pval_init[j]+refi_rat*pval_refi[j];
          }
          break;
        case 0x0110: // e_skip==1 && refi_skip==1 && refr_skip==0
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = init_rat*pval_init[j]+refr_rat*pval_refr[j];
          }
          break;
        case 0x0001: // e_skip==0 && refi_skip==0 && refr_skip==1
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = init_rat*pval_init[j]+refi_rat*pval_refi[j]+e_rat*pval_e[j];
          }
          break;
        case 0x0010: // e_skip==0 && refi_skip==1 && refr_skip==0
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = init_rat*pval_init[j]+refr_rat*pval_refr[j]+e_rat*pval_e[j];
          }
          break;
        case 0x0100: // e_skip==1 && refi_skip==0 && refr_skip==0
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = init_rat*pval_init[j]+refr_rat*pval_refr[j]+refi_rat*pval_refi[j];
          }
          break;
        default: // e_skip==0 && refi_skip==0 && refr_skip==0
          for(j=0,p=tp_size*i; j<tp_size; j++,p++)
          {
            tp_para[k][p] = init_rat*pval_init[j]+refr_rat*pval_refr[j]+refi_rat*pval_refi[j]+e_rat*pval_e[j];
          }
          break;
      }
    }
  }

  return 0;
}

int IntTable(void)
{
  int i,j,j1;
  int jmax = tab_x_size-2;
  int k,k1,k2,km,ktmp;
  int m,n;
  int err = 0;
  double x;
  double fact;
  double offs;
  double w1,w2;
  const double xmin = tab_x_val[0];
  const double xmax = tab_x_val[tab_x_size-2];
  const char fnam[] = "IntTable";

  for(n=0; n<mie_n_comp; n++)
  {
    if(mie_n_size[n] > 0)
    {
      mie_aext_int[n] = (double *)malloc(mie_n_wlen*mie_n_size[n]*sizeof(double));
      mie_asca_int[n] = (double *)malloc(mie_n_wlen*mie_n_size[n]*sizeof(double));
      mie_phs1_int[n] = (double *)malloc(mie_n_wlen*mie_n_size[n]*mie_n_angl*sizeof(double));
      mie_phs2_int[n] = (double *)malloc(mie_n_wlen*mie_n_size[n]*mie_n_angl*sizeof(double));
      if(mie_aext_int[n]==NULL || mie_asca_int[n]==NULL || mie_phs1_int[n]==NULL || mie_phs2_int[n]==NULL)
      {
        fprintf(stderr,"%s: error in allocating memory.\n",fnam);
        err = 1;
        break;
      }
      for(i=0; i<mie_n_wlen; i++)
      {
        fact = PI2/(mie_wlen_um[i]);
        offs = log10(fact); // log10(x) = log10(r) + offs
        j1 = 0;
        for(m=0,k=mie_n_size[n]*i,km=mie_n_size[n]*mie_n_angl*i; m<mie_n_size[n]; m++,k++,km+=mie_n_angl)
        {
          x = mie_xval_com[n][m]*fact;
          if(x <= xmin)
          {
            j1 = 0;
            if(cnt_vb > 1)
            {
              fprintf(stderr,"%s: warning, x=%13.6e <= xmin=%13.6e\n",fnam,x,xmin);
            }
          } else
          if(x >= xmax)
          {
            j1 = jmax;
            if(cnt_vb > 1)
            {
              fprintf(stderr,"%s: warning, x=%13.6e >= xmax=%13.6e\n",fnam,x,xmax);
            }
          }
          else
          {
            for(j=j1; j<jmax; j++)
            {
              if(x>=tab_x_val[j] && x<tab_x_val[j+1])
              {
                j1 = j;
                break;
              }
            }
            if(j == jmax)
            {
              fprintf(stderr,"%s: error, faild in finding range.\n",fnam);
              err = 1;
              break;
            }
          }
          w2 = (mie_lval_com[n][m]+offs-tab_x_log_val[j1])*tab_x_log_stp_inv;
          w1 = 1.0-w2;
          k1 = tab_x_size*i+j1;
          k2 = k1+1;
          mie_aext_int[n][k] = w1*mie_aext_tab[n][k1]+w2*mie_aext_tab[n][k2];
          mie_asca_int[n][k] = w1*mie_asca_tab[n][k1]+w2*mie_asca_tab[n][k2];
          for(j=0,ktmp=km,
              k1=tab_x_size*mie_n_angl*i+mie_n_angl*j1,
              k2=k1+mie_n_angl; j<mie_n_angl; j++,ktmp++,k1++,k2++)
          {
            mie_phs1_int[n][ktmp] = w1*mie_phs1_tab[n][k1]+w2*mie_phs1_tab[n][k2];
            mie_phs2_int[n][ktmp] = w1*mie_phs2_tab[n][k1]+w2*mie_phs2_tab[n][k2];
          }
        }
        if(err)
        {
          break;
        }
      }
      if(err)
      {
        free(mie_aext_int[n]);
        free(mie_asca_int[n]);
        free(mie_phs1_int[n]);
        free(mie_phs2_int[n]);
        break;
      }
    }
  }
  if(err)
  {
    return -1;
  }

  return 0;
}

int MieCalc(void)
{
  int i,j,k,m,n,p;
  int m1,m2;
  int err = 0;
  double v,w,y,r,l;
  double sw,sy,sp,sinv;
  double s1,s2,s3,s4;
  double cnst,norm,fact,indx,offs;
  double epsilon = 1.0e-1;
  const char fnam[] = "MieCalc";

  for(n=0; n<mie_n_comp; n++)
  {
    if(mie_size_shape[n] == MIE_SHAPE_SPHEROID)
    {
      if(TmxTable(n) < 0) return -1;
    }
    else
    {
      if(MieTable(n) < 0) return -1;
    }
  }
  if(IntTable() < 0) return -1;

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

  // initialization
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
      if(cnt_vb > 1 || cnt_xmod == 2)
      {
        mie_lmin = mie_lval_com[n][0];
        mie_lmax = mie_lval_com[n][mie_n_size[n]-1];
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
      }
      switch(mie_size_func[n])
      {
        case MIE_FUNC_LOGNORMAL:
        case MIE_FUNC_EFF_LOGNORMAL:
          cnst = 1.0/(sqrt(PI2)*mie_lsgm_com[n]);
          fact = -0.5/(mie_lsgm_com[n]*mie_lsgm_com[n]);
          break;
        case MIE_FUNC_GAMMA:
          indx = 1.0/mie_lsgm_com[n]-2.0;
          cnst = M_LNA/(pow(mie_rmod_com[n]*mie_lsgm_com[n],indx)*tgamma(indx));
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
              for(m=0,sy=0.0; m<mie_n_size[n]; m++)
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
              for(m=0,sy=0.0; m<mie_n_size[n]; m++)
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
              norm = cnst*mie_lstp;
              for(l=mie_lmin,sy=0.0; l<mie_lmax; l+=mie_lstp)
              {
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
              norm = cnst*mie_lstp;
              for(l=mie_lmin,sy=0.0; l<mie_lmax; l+=mie_lstp)
              {
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
              norm = cnst*mie_lstp;
              for(l=mie_lmin,sy=0.0; l<mie_lmax; l+=mie_lstp)
              {
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
              norm = cnst*mie_lstp;
              for(l=mie_lmin,sy=0.0; l<mie_lmax; l+=mie_lstp)
              {
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
      switch(mie_size_func[n])
      {
        case MIE_FUNC_FILE:
          switch(mie_size_ytype[n])
          {
            case MIE_YTYPE_NUMBER: // Number distribution
              for(m=0,sy=0.0; m<mie_n_size[n]; m++)
              {
                r = mie_xval_com[n][m];
                if(r > mie_rmax_com[n][i])
                {
                  break;
                }
                y = mie_yval_com[n][m];
                p = mie_n_size[n]*i+m;
                mie_aext_com[n][i] += mie_aext_int[n][p]*y;
                mie_asca_com[n][i] += mie_asca_int[n][p]*y;
                for(j=0,k=mie_n_angl*i,p=mie_n_size[n]*mie_n_angl*i+mie_n_angl*m; j<mie_n_angl; j++,k++,p++)
                {
                  mie_phs1_com[n][k] += mie_phs1_int[n][p]*y;
                  mie_phs2_com[n][k] += mie_phs1_int[n][p]*y;
                }
                sy += y;
              }
              sw = sy;
              break;
            case MIE_YTYPE_VOLUME: // Volume distribution
              for(m=0,sw=0.0,sy=0.0; m<mie_n_size[n]; m++)
              {
                r = mie_xval_com[n][m];
                if(r > mie_rmax_com[n][i])
                {
                  break;
                }
                v = mie_yval_com[n][m];
                y = v/(PI4_3*r*r*r);
                p = mie_n_size[n]*i+m;
                mie_aext_com[n][i] += mie_aext_int[n][p]*y;
                mie_asca_com[n][i] += mie_asca_int[n][p]*y;
                for(j=0,k=mie_n_angl*i,p=mie_n_size[n]*mie_n_angl*i+mie_n_angl*m; j<mie_n_angl; j++,k++,p++)
                {
                  mie_phs1_com[n][k] += mie_phs1_int[n][p]*y;
                  mie_phs2_com[n][k] += mie_phs2_int[n][p]*y;
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
              norm = cnst*tab_x_log_stp;
              offs = log10(PI2/mie_wlen_um[i]); // log10(x) = log10(r) + offs
              m1 = (int)((mie_lmin+offs-tab_x_log_min)*tab_x_log_stp_inv);
              m2 = (int)((mie_lmax+offs-tab_x_log_min)*tab_x_log_stp_inv);
              if(m1 < 0)
              {
                m1 = 0;
              }
              if(m2 >= tab_x_size)
              {
                m2 = tab_x_size-1;
              }
              for(m=m1,l=tab_x_log_min+tab_x_log_stp*m1-offs,sy=0.0; m<m2; m++,l+=tab_x_log_stp)
              {
                r = pow(10.0,l);
                if(r > mie_rmax_com[n][i])
                {
                  break;
                }
                y = norm*exp(fact*(l-mie_lmod_com[n])*(l-mie_lmod_com[n]));
                p = tab_x_size*i+m;
                mie_aext_com[n][i] += mie_aext_tab[n][p]*y;
                mie_asca_com[n][i] += mie_asca_tab[n][p]*y;
                for(j=0,k=mie_n_angl*i,p=tab_x_size*mie_n_angl*i+mie_n_angl*m; j<mie_n_angl; j++,k++,p++)
                {
                  mie_phs1_com[n][k] += mie_phs1_tab[n][p]*y;
                  mie_phs2_com[n][k] += mie_phs2_tab[n][p]*y;
                }
                sy += y;
              }
              sw = sy;
              break;
            case MIE_YTYPE_VOLUME: // Volume distribution
              norm = cnst*tab_x_log_stp;
              offs = log10(PI2/mie_wlen_um[i]); // log10(x) = log10(r) + offs
              m1 = (int)((mie_lmin+offs-tab_x_log_min)*tab_x_log_stp_inv);
              m2 = (int)((mie_lmax+offs-tab_x_log_min)*tab_x_log_stp_inv);
              if(m1 < 0)
              {
                m1 = 0;
              }
              if(m2 >= tab_x_size)
              {
                m2 = tab_x_size-1;
              }
              for(m=m1,l=tab_x_log_min+tab_x_log_stp*m1-offs,sw=0.0,sy=0.0; m<m2; m++,l+=tab_x_log_stp)
              {
                r = pow(10.0,l);
                if(r > mie_rmax_com[n][i])
                {
                  break;
                }
                v = norm*exp(fact*(l-mie_lmod_com[n])*(l-mie_lmod_com[n]));
                y = v/(PI4_3*r*r*r);
                p = tab_x_size*i+m;
                mie_aext_com[n][i] += mie_aext_tab[n][p]*y;
                mie_asca_com[n][i] += mie_asca_tab[n][p]*y;
                for(j=0,k=mie_n_angl*i,p=tab_x_size*mie_n_angl*i+mie_n_angl*m; j<mie_n_angl; j++,k++,p++)
                {
                  mie_phs1_com[n][k] += mie_phs1_tab[n][p]*y;
                  mie_phs2_com[n][k] += mie_phs2_tab[n][p]*y;
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
              norm = cnst*tab_x_log_stp;
              offs = log10(PI2/mie_wlen_um[i]); // log10(x) = log10(r) + offs
              m1 = (int)((mie_lmin+offs-tab_x_log_min)*tab_x_log_stp_inv);
              m2 = (int)((mie_lmax+offs-tab_x_log_min)*tab_x_log_stp_inv);
              if(m1 < 0)
              {
                m1 = 0;
              }
              if(m2 >= tab_x_size)
              {
                m2 = tab_x_size-1;
              }
              for(m=m1,l=tab_x_log_min+tab_x_log_stp*m1-offs,sy=0.0; m<m2; m++,l+=tab_x_log_stp)
              {
                r = pow(10.0,l);
                if(r > mie_rmax_com[n][i])
                {
                  break;
                }
                y = norm*pow(r,indx)*exp(fact*r);
                p = tab_x_size*i+m;
                mie_aext_com[n][i] += mie_aext_tab[n][p]*y;
                mie_asca_com[n][i] += mie_asca_tab[n][p]*y;
                for(j=0,k=mie_n_angl*i,p=tab_x_size*mie_n_angl*i+mie_n_angl*m; j<mie_n_angl; j++,k++,p++)
                {
                  mie_phs1_com[n][k] += mie_phs1_tab[n][p]*y;
                  mie_phs2_com[n][k] += mie_phs2_tab[n][p]*y;
                }
                sy += y;
              }
              sw = sy;
              break;
            case MIE_YTYPE_VOLUME: // Volume distribution
              norm = cnst*tab_x_log_stp;
              offs = log10(PI2/mie_wlen_um[i]); // log10(x) = log10(r) + offs
              m1 = (int)((mie_lmin+offs-tab_x_log_min)*tab_x_log_stp_inv);
              m2 = (int)((mie_lmax+offs-tab_x_log_min)*tab_x_log_stp_inv);
              if(m1 < 0)
              {
                m1 = 0;
              }
              if(m2 >= tab_x_size)
              {
                m2 = tab_x_size-1;
              }
              for(m=m1,l=tab_x_log_min+tab_x_log_stp*m1-offs,sw=0.0,sy=0.0; m<m2; m++,l+=tab_x_log_stp)
              {
                r = pow(10.0,l);
                if(r > mie_rmax_com[n][i])
                {
                  break;
                }
                v = norm*pow(r,indx)*exp(fact*r);
                y = v/(PI4_3*r*r*r);
                p = tab_x_size*i+m;
                mie_aext_com[n][i] += mie_aext_tab[n][p]*y;
                mie_asca_com[n][i] += mie_asca_tab[n][p]*y;
                for(j=0,k=mie_n_angl*i,p=tab_x_size*mie_n_angl*i+mie_n_angl*m; j<mie_n_angl; j++,k++,p++)
                {
                  mie_phs1_com[n][k] += mie_phs1_tab[n][p]*y;
                  mie_phs2_com[n][k] += mie_phs2_tab[n][p]*y;
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
        fprintf(stderr,"%s: warning, sw=%13.4e\n",fnam,sw);
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

  if(err)
  {
    return -1;
  }

  return 0;
}

int Init(void)
{
  int i,n;
  int err = 0;
  int flag_file = 0;
  double v;
  double dang;
  const char fnam[] = "Init";

  mie_n_angl = tab_n_angl;
  dang = 180.0/(mie_n_angl-1);
  for(i=0; i<mie_n_angl; i++)
  {
    mie_angl[i] = dang*i;
  }

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
      flag_file = 1;
      if((mie_lval_com[n]=(double *)malloc(mie_n_size[n]*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: error in allocating memory.\n",fnam);
        return -1;
      }
      switch(mie_size_xtype[n])
      {
        case MIE_XTYPE_R:
          for(i=0; i<mie_n_size[n]; i++)
          {
            mie_lval_com[n][i] = log10(mie_xval_com[n][i]);
          }
          break;
        case MIE_XTYPE_LOGR:
          for(i=0; i<mie_n_size[n]; i++)
          {
            mie_lval_com[n][i] = mie_xval_com[n][i];
            mie_xval_com[n][i] = pow(10.0,mie_xval_com[n][i]); // Log10(R) -> R
          }
          break;
        case MIE_XTYPE_LNR:
          for(i=0; i<mie_n_size[n]; i++)
          {
            mie_lval_com[n][i] = mie_xval_com[n][i]*M_1_LNA; // Ln(R) -> Log10(R)
            mie_xval_com[n][i] = exp(mie_xval_com[n][i]);
          }
          break;
        default:
          fprintf(stderr,"%s: error, mie_size_xtype[%d]=%d\n",fnam,n,mie_size_xtype[n]);
          return -1;
          break;
      }
    }
    else
    {
      mie_lmod_com[n] = log10(mie_rmod_com[n]);
    }
    if(mie_n_reps[n] > 0)
    {
      if(mie_n_reps[n] != 1)
      {
        fprintf(stderr,"%s: not implemented error, mie_n_reps[%d]=%d\n",
                        fnam,n,mie_n_reps[n]);
        return -1;
      }
      cnt_tmx = 1;
    }
    mie_refi_log_com[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_rmax_com[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    if(mie_refi_log_com[n]==NULL || mie_rmax_com[n]==NULL)
    {
      fprintf(stderr,"%s: error in allocating memory.\n",fnam);
      return -1;
    }
    for(i=0; i<mie_n_wlen; i++)
    {
      mie_refi_log_com[n][i] = log10(fabs(mie_refi_com[n][i]));
      mie_rmax_com[n][i] = pow(10.0,MIE_LMAX);
    }
  }

  tab_x_log_stp_inv = 1.0/tab_x_log_stp;
  tab_refr_stp_inv = 1.0/tab_refr_stp;
  tab_refi_log_stp_inv = 1.0/tab_refi_log_stp;
  tab_x_size = (int)((tab_x_log_max+1.1*tab_x_log_stp-tab_x_log_min)*tab_x_log_stp_inv);
  tab_refr_size = (int)((tab_refr_max+1.1*tab_refr_stp-tab_refr_min)*tab_refr_stp_inv);
  tab_refi_size = (int)((tab_refi_log_max+1.1*tab_refi_log_stp-tab_refi_log_min)*tab_refi_log_stp_inv);
  if(flag_file)
  {
    tab_x_val = (double *)malloc(tab_x_size*sizeof(double));
    tab_x_log_val = (double *)malloc(tab_x_size*sizeof(double));
    if(tab_x_val==NULL || tab_x_log_val==NULL)
    {
      fprintf(stderr,"%s: error in allocating memory.\n",fnam);
      return -1;
    }
    for(i=0,v=tab_x_log_min; i<tab_x_size; i++,v+=tab_x_log_stp)
    {
      tab_x_val[i] = pow(10.0,v);
      tab_x_log_val[i] = v;
    }
  }
  else
  {
    tab_x_val = (double *)malloc(tab_x_size*sizeof(double));
    if(tab_x_val == NULL)
    {
      fprintf(stderr,"%s: error in allocating memory.\n",fnam);
      return -1;
    }
    for(i=0,v=tab_x_log_min; i<tab_x_size; i++,v+=tab_x_log_stp)
    {
      tab_x_val[i] = pow(10.0,v);
    }
  }
  tab_refr_val = (double *)malloc(tab_refr_size*sizeof(double));
  tab_refi_log_val = (double *)malloc(tab_refi_size*sizeof(double));
  qval_init = (double *)malloc((tab_x_size+1)*sizeof(double));
  qval_refr = (double *)malloc((tab_x_size+1)*sizeof(double));
  qval_refi = (double *)malloc((tab_x_size+1)*sizeof(double));
  pval_init = (double *)malloc((tab_x_size*mie_n_angl+1)*sizeof(double));
  pval_refr = (double *)malloc((tab_x_size*mie_n_angl+1)*sizeof(double));
  pval_refi = (double *)malloc((tab_x_size*mie_n_angl+1)*sizeof(double));
  if(tab_refr_val==NULL || tab_refi_log_val==NULL ||
     qval_init==NULL || qval_refr==NULL || qval_refi==NULL ||
     pval_init==NULL || pval_refr==NULL || pval_refi==NULL)
  {
    fprintf(stderr,"%s: error in allocating memory.\n",fnam);
    return -1;
  }
  for(i=0,v=tab_refr_min; i<tab_refr_size; i++,v+=tab_refr_stp)
  {
    tab_refr_val[i] = v;
  }
  for(i=0,v=tab_refi_log_min; i<tab_refi_size; i++,v+=tab_refi_log_stp)
  {
    tab_refi_log_val[i] = v;
  }
  if(cnt_tmx)
  {
    tab_e_stp_inv = 1.0/tab_e_stp;
    tab_e_size = (int)((tab_e_max+1.1*tab_e_stp-tab_e_min)*tab_e_stp_inv);
    tab_e_val = (double *)malloc(tab_e_size*sizeof(double));
    qval_e = (double *)malloc((tab_x_size+1)*sizeof(double));
    pval_e = (double *)malloc((tab_x_size*mie_n_angl+1)*sizeof(double));
    if(tab_e_val==NULL || qval_e==NULL || pval_e==NULL)
    {
      fprintf(stderr,"%s: error in allocating memory.\n",fnam);
      return -1;
    }
    for(i=0,v=tab_e_min; i<tab_e_size; i++,v+=tab_e_stp)
    {
      tab_e_val[i] = v;
    }
  }

  for(n=0; n<mie_n_comp; n++)
  {
    mie_aext_tab[n] = (double *)malloc(mie_n_wlen*tab_x_size*sizeof(double));
    mie_asca_tab[n] = (double *)malloc(mie_n_wlen*tab_x_size*sizeof(double));
    mie_phs1_tab[n] = (double *)malloc(mie_n_wlen*tab_x_size*mie_n_angl*sizeof(double));
    mie_phs2_tab[n] = (double *)malloc(mie_n_wlen*tab_x_size*mie_n_angl*sizeof(double));
    if(mie_aext_tab[n]==NULL || mie_asca_tab[n]==NULL || mie_phs1_tab[n]==NULL || mie_phs2_tab[n]==NULL)
    {
      fprintf(stderr,"%s: error in allocating memory.\n",fnam);
      err = 1;
      break;
    }
  }
  if(err)
  {
    return -1;
  }

  return 0;
}

void Finish(void)
{
  int n;

  for(n=0; n<mie_n_comp; n++)
  {
    free(mie_aext_tab[n]);
    free(mie_asca_tab[n]);
    free(mie_phs1_tab[n]);
    free(mie_phs2_tab[n]);
    if(mie_n_size[n] > 0)
    {
      free(mie_lval_com[n]);
      free(mie_aext_int[n]);
      free(mie_asca_int[n]);
      free(mie_phs1_int[n]);
      free(mie_phs2_int[n]);
    }
    free(mie_refi_log_com[n]);
    free(mie_rmax_com[n]);
  }
  free(tab_x_val);
  free(tab_x_log_val);
  free(tab_refr_val);
  free(tab_refi_log_val);
  free(qval_init);
  free(qval_refr);
  free(qval_refi);
  free(pval_init);
  free(pval_refr);
  free(pval_refi);
  if(cnt_tmx)
  {
    free(tab_e_val);
    free(qval_e);
    free(pval_e);
  }

  CommonFinish();
}

int ReadConfig(void)
{
  int n,nc;
  int err;
  int ntmp;
  double xtmp;
  char line[MAXLINE];
  char temp[MAXLINE];
  char str[MAXITEM][MAXLINE];
  char *p;
  FILE *fp;
  char fnam[] = "ReadConfig";

  if(strcmp(cnt_conf,"")==0 || strcmp(cnt_conf,NONAME)==0)
  {
    return 0;
  }
  // Read parameters
  if((fp=fopen(cnt_conf,"r")) == NULL)
  {
    fprintf(stderr,"%s: error, cannot open %s\n",fnam,cnt_conf);
    return -1;
  }
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    strncpy(temp,line,MAXLINE);
    if((p=strchr(temp,'#')) != NULL) *p = '\0';
    // read line
    for(n=nc=0,p=temp; n<MAXITEM; n++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
    }
    if(n < 1) continue;
    if(strcasecmp(str[0],"tab_x_log_min") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_x_log_min = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_x_log_min);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_x_log_max") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_x_log_max = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_x_log_max);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_x_log_stp") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_x_log_stp = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_x_log_stp);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_refr_ndig") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        ntmp = strtol(str[1],&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>0) tab_refr_ndig = ntmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30d\n",str[0],tab_refr_ndig);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_refr_mult") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_refr_mult = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_refr_mult);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_refr_min") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_refr_min = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_refr_min);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_refr_max") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_refr_max = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_refr_max);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_refr_stp") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_refr_stp = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_refr_stp);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_refi_ndig") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        ntmp = strtol(str[1],&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>0) tab_refi_ndig = ntmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30d\n",str[0],tab_refi_ndig);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_refi_mult") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_refi_mult = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_refi_mult);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_refi_log_min") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_refi_log_min = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_refi_log_min);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_refi_log_max") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_refi_log_max = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_refi_log_max);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_refi_log_stp") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_refi_log_stp = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_refi_log_stp);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_e_ndig") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        ntmp = strtol(str[1],&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>0) tab_e_ndig = ntmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30d\n",str[0],tab_e_ndig);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_e_mult") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_e_mult = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_e_mult);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_e_min") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_e_min = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_e_min);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_e_max") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_e_max = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_e_max);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tab_e_stp") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        xtmp = strtod(str[1],&p);
        if(errno!=ERANGE && *p=='\0') tab_e_stp = xtmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.14e\n",str[0],tab_e_stp);
        cnt_n_cmnt++;
      }
    }
    else
    {
      fprintf(stderr,"%s: unrecognized option >>> %s\n",fnam,line);
      break;
    }
  }
  fclose(fp);
  if(err)
  {
    return -1;
  }

  return 0;
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
    {"tab_datdir",1,0,'D'},
    {"tab_n_angl",1,0,'N'},
    {"mie_rmin",1,0,'x'},
    {"mie_rmax",1,0,'X'},
    {"mie_lstp",1,0,'s'},
    {"mie_wsgm",1,0,'S'},
    {"cnt_conf",1,0,'C'},
    {"help",0,0,'h'}
  };

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-D:N:x:X:s:S:C:h",long_options,&option_index);
    if(c == -1) break;
    switch(c)
    {
      case 'D':
        strncpy(tab_datdir,optarg,MAXLINE);
        break;
      case 'N':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) tab_n_angl = ntmp;
        else
        {
          fprintf(stderr,"#Angles -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
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
      case 'C':
        strncpy(cnt_conf,optarg,MAXLINE);
        break;
      case 'h':
        cnt_hp = 1;
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

  if(ReadConfig() < 0) return -1;

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
  char optstring[] = "fgwlLaWAPQiIruDNxXsSmCoOdvh";

  CommonUsage("aeros_table",0);
  for(i=0; i<strlen(optstring); i++)
  {
    switch(optstring[i])
    {
      case 'D':
        fprintf(stderr," D -tab_datdir    |%s|%s|%s| %s\n",As(e,"Data directory",n), As(a,"name",n),       As(d,TAB_DATDIR,n),tab_datdir);
        break;
      case 'N':
        fprintf(stderr," N -tab_n_angl    |%s|%s|%s| %d\n",As(e,"#Angles",n),        As(a,"#",n),          Ad(d,TAB_N_ANGL,n),tab_n_angl);
        break;
      case 'x':
        fprintf(stderr," x -mie_rmin      |%s|%s|%s| %e\n",As(e,"R Min",n),          As(a,"um",n),         Ae(d,MIE_RMIN,n),mie_rmin);
        break;
      case 'X':
        fprintf(stderr," X -mie_rmax      |%s|%s|%s| %e\n",As(e,"R Max",n),          As(a,"um",n),         Ae(d,MIE_RMAX,n),mie_rmax);
        break;
      case 's':
        fprintf(stderr," s -mie_lstp      |%s|%s|%s| %e\n",As(e,"Log10(R) step",n),  As(a,"value",n),      Ae(d,MIE_LSTP,n),mie_lstp);
        break;
      case 'S':
        fprintf(stderr," S -mie_wsgm      |%s|%s|%s| %f\n",As(e,"Log10(R) width",n), As(a,"sigma",n),      Af(d,MIE_WSGM,n),mie_wsgm);
        break;
      case 'C':
        fprintf(stderr," C -cnt_conf      |%s|%s|%s| %s\n",As(e,"Config  file",n),   As(a,"name",n),       As(d,CNT_CONF,n),cnt_conf);
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
  fprintf(stderr,"Config file (cnt_conf) format:\n");
  fprintf(stderr,"tab_x_log_min    value     | Log(x) min. (%.1f)\n",TAB_X_LOG_MIN);
  fprintf(stderr,"tab_x_log_max    value     | Log(x) max. (%.1f)\n",TAB_X_LOG_MAX);
  fprintf(stderr,"tab_x_log_stp    value     | Log(x) step (%.3f)\n",TAB_X_LOG_STP);
  fprintf(stderr,"tab_refr_ndig    value     | refr #digits (%d)\n",TAB_REFR_NDIG);
  fprintf(stderr,"tab_refr_mult    value     | refr multiplying factor (%.1e)\n",TAB_REFR_MULT);
  fprintf(stderr,"tab_refr_min     value     | refr min. (%.1f)\n",TAB_REFR_MIN);
  fprintf(stderr,"tab_refr_max     value     | refr max. (%.1f)\n",TAB_REFR_MAX);
  fprintf(stderr,"tab_refr_stp     value     | refr step (%.2f)\n",TAB_REFR_STP);
  fprintf(stderr,"tab_refi_ndig    value     | refi #digits (%d)\n",TAB_REFI_NDIG);
  fprintf(stderr,"tab_refi_mult    value     | refi multiplying factor (%.1e)\n",TAB_REFI_MULT);
  fprintf(stderr,"tab_refi_log_min value     | Log(refi) min. (%.1f)\n",TAB_REFI_LOG_MIN);
  fprintf(stderr,"tab_refi_log_max value     | Log(refi) max. (%.1f)\n",TAB_REFI_LOG_MAX);
  fprintf(stderr,"tab_refi_log_stp value     | Log(refi) step (%.1f)\n",TAB_REFI_LOG_STP);
  fprintf(stderr,"tab_e_ndig       value     | e #digits (%d)\n",TAB_E_NDIG);
  fprintf(stderr,"tab_e_mult       value     | e multiplying factor (%.1e)\n",TAB_E_MULT);
  fprintf(stderr,"tab_e_min        value     | e min. (%.1f)\n",TAB_E_MIN);
  fprintf(stderr,"tab_e_max        value     | e max. (%.1f)\n",TAB_E_MAX);
  fprintf(stderr,"tab_e_stp        value     | e step (%.3f)\n",TAB_E_STP);
  if(cnt_n_cmnt > 0)
  {
    fprintf(stderr,"<Read from %s>\n",cnt_conf);
    for(n=0; n<cnt_n_cmnt; n++)
    {
      fprintf(stderr,"%s",cnt_cmnt[n]);
    }
  }

  return 0;
}
