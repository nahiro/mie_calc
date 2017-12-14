/***********************************************************/
/* AEROS_MIE2NEW ... Aerosol generator with mie2new        */
/* Author: N.Manago                                        */
/* $Revision: 1104 $                                       */
/* $Date: 2015-09-02 00:16:06 +0900 (Wed, 02 Sep 2015) $   */
/***********************************************************/
#define	MIE_RMIN		NAN				// R min in um
#define	MIE_RMAX		NAN				// R max in um
#define	MIE_LMIN		-3.0				// Min Log10(R)
#define	MIE_LMAX		+2.0				// Max Log10(R)
#define	MIE_WSGM		5.0				// Log10(R) width in sigma
#define	MIE_NSTP		799				// R step
#define	MAXSTEP			799				// Max #steps
#define	MAXSIZE			800				// Max #size
#define	DELTA			1.0e-7				// A small number

#include "aeros_common.c"

int	mie_nstp		= MIE_NSTP;
double	mie_rmin		= MIE_RMIN;			// R min in um
double	mie_rmax		= MIE_RMAX;			// R max in um
double	mie_wsgm		= MIE_WSGM;			// Log10(R) width in sigma
double	mie_rmin_com[MIE_MAXCOMP];				// R min in um
double	mie_rmax_com[MIE_MAXCOMP];				// R min in um
double	mie_lmod_com[MIE_MAXCOMP];				// Log10(Mode radius in um)
int	cnt_2peak		= 0;				// Double peak mode
int	cnt_norun		= 0;				// No run  mode

int ReadAngle(double *dthet,double *chang);
int ReadTape8(void);

extern void newrad_(double *r,double *y,int *n,int *m);

int MieCalc(void)
{
  const double SQRPI2 = 5.7717248988427263; // sqrt(2pi)*ln(10)
  int mie_nsiz;
  int i,j,n;
  int m = 2; // Moment of the size distribution
  double w,l;
  double sy;
  double s1,s2,s3,s4;
  double delt,norm,fact;
  double r[MAXSIZE];
  double d[MAXSIZE];
  double y[MAXSIZE];
  double dthet[5];
  double chang[4];
  FILE *fp;

  if(!isnan(mie_rmin))
  {
    for(n=0; n<mie_n_comp; n++)
    {
      mie_rmin_com[n] = mie_rmin;
    }
  }
  else
  {
    for(n=0; n<mie_n_comp; n++)
    {
      l = log10(mie_rmod_com[n])-mie_lsgm_com[n]*mie_wsgm;
      if(l < MIE_LMIN)
      {
        l = MIE_LMIN;
      }
      mie_rmin_com[n] = pow(10.0,l);
    }
  }
  if(!isnan(mie_rmax))
  {
    for(n=0; n<mie_n_comp; n++)
    {
      mie_rmax_com[n] = mie_rmax;
    }
  }
  else
  {
    for(n=0; n<mie_n_comp; n++)
    {
      l = log10(mie_rmod_com[n])+mie_lsgm_com[n]*mie_wsgm;
      if(l > MIE_LMAX)
      {
        l = MIE_LMAX;
      }
      mie_rmax_com[n] = pow(10.0,l);
    }
  }

  if(cnt_vb > 1 || cnt_xmod == 2)
  {
    mie_nsiz = mie_nstp+1;
    for(n=0; n<mie_n_comp; n++)
    {
      if(cnt_vb > 1)
      {
        fprintf(stderr,"Component#   : %2d\n",n);
        fprintf(stderr,"Min Log10(R) : %13.4e\n",log10(mie_rmin_com[n]));
        fprintf(stderr,"Max Log10(R) : %13.4e\n",log10(mie_rmax_com[n]));
        fprintf(stderr,"#Steps       : %13d\n",mie_nstp);
      }
      mie_lmod_com[n] = log10(mie_rmod_com[n]);
      delt = 2.0*(mie_rmax_com[n]-mie_rmin_com[n])/(mie_nstp*mie_nsiz);
      r[0] = mie_rmin_com[n];
      for(i=1; i<mie_nsiz; i++)
      {
        r[i] = r[i-1]+(double)i*delt;
      }
      norm = 1.0/(mie_lsgm_com[n]*SQRPI2);
      fact = -0.5/(mie_lsgm_com[n]*mie_lsgm_com[n]);
      for(i=0; i<mie_nsiz; i++)
      {
        l = log10(r[i]);
        y[i] = norm/r[i]*exp(fact*(l-mie_lmod_com[n])*(l-mie_lmod_com[n]));
      }
      for(j=0; j<2; j++)
      {
        newrad_(r,y,&mie_nsiz,&m);
        for(i=0; i<mie_nsiz; i++)
        {
          l = log10(r[i]);
          y[i] = norm/r[i]*exp(fact*(l-mie_lmod_com[n])*(l-mie_lmod_com[n]));
        }
      }
      d[0] = r[1]-r[0];
      for(i=1; i<mie_nstp; i++)
      {
        d[i] = 0.5*(r[i+1]-r[i-1]);
      }
      d[mie_nstp] = r[mie_nstp]-r[mie_nstp-1];
      // calculate averages
      s1 = s2 = s3 = s4 = 0.0;
      for(i=0,sy=0.0; i<mie_nsiz; i++)
      {
        w = y[i]*d[i];
        sy += w; w *= r[i];
        s1 += w; w *= r[i];
        s2 += w; w *= r[i];
        s3 += w; w *= r[i];
        s4 += w;
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
  }

  if(!cnt_norun)
  {
    if(ReadAngle(dthet,chang) < 0) return -1;
    if((fp=fopen("mie2new.tp5","w")) == NULL)
    {
      fprintf(stderr,"Error, cannot open %s\n","mie2new.tp5");
      return -1;
    }
    if(cnt_2peak)
    {
      mie_n_comp = 1;
      fprintf(fp,"%5d    1    0\n",mie_n_comp);
      fprintf(fp,"%5d%5d    0    0  component#%d\n",1,mie_n_wlen,0);
      fprintf(fp," %9.2f %9.2f %9.2f %9.2f %9.2f\n",dthet[0],chang[0],dthet[1],chang[1],dthet[2]);
      fprintf(fp," %9.2f %9.2f %9.2f %9.2f\n",chang[2],dthet[3],chang[3],dthet[4]);
      fprintf(fp,"%10.3e%10.3e   10000.0 %4d    2 %4d    3  component#%d\n",
                  MIN(mie_rmin_com[0],mie_rmin_com[1]),MAX(mie_rmax_com[0],mie_rmax_com[1]),mie_nstp,m,0);
      fprintf(fp,"%15.8e%15.8e%15.8e\n",mie_wcom_com[0],mie_rmod_com[0],mie_lsgm_com[0]); // 1st component
      fprintf(fp,"%15.8e%15.8e%15.8e\n",mie_wcom_com[1],mie_rmod_com[1],mie_lsgm_com[1]); // 2nd component
      for(i=0; i<mie_n_wlen; i++)
      {
        fprintf(fp,"%15.8e%15.8e%15.8e component#%d\n",mie_wlen_um[i],mie_refr_com[0][i],-fabs(mie_refi_com[0][i]),0);
      }
    }
    else
    {
      fprintf(fp,"%5d    1    0\n",mie_n_comp);
      for(n=0; n<mie_n_comp; n++)
      {
        fprintf(fp,"%5d%5d    0    0  component#%d\n",n+1,mie_n_wlen,n);
        fprintf(fp," %9.2f %9.2f %9.2f %9.2f %9.2f\n",dthet[0],chang[0],dthet[1],chang[1],dthet[2]);
        fprintf(fp," %9.2f %9.2f %9.2f %9.2f\n",chang[2],dthet[3],chang[3],dthet[4]);
        fprintf(fp,"%10.3e%10.3e   10000.0 %4d    2 %4d    3  component#%d\n",mie_rmin_com[n],mie_rmax_com[n],mie_nstp,m,n);
        fprintf(fp,"%15.8e%15.8e%15.8e\n",1.0,mie_rmod_com[n],mie_lsgm_com[n]); // 1st component
        fprintf(fp,"%15.8e%15.8e%15.8e\n",0.0,1.0,0.1); // 2nd component -> ignored
        for(i=0; i<mie_n_wlen; i++)
        {
          fprintf(fp,"%15.8e%15.8e%15.8e component#%d\n",mie_wlen_um[i],mie_refr_com[n][i],-fabs(mie_refi_com[n][i]),n);
        }
      }
    }
    fclose(fp);
    system("mie2new.exe");
  }

  // read results
  if(ReadTape8() < 0) return -1;

  return 0;
}

int ReadAngle(double *dthet,double *chang)
{
  // dthet must have size >= 5
  // chang must have size >= 4
  int i;
  int nchang;
  double *chang_array;
  double *dthet_array;
  double dtmp;
  char fnam[] = "ReadAngle";

  // Minimum #angles = 6
  // 0  0.0
  // 1  dthet1
  // 2  dthet1+dthet2 (chang1=mie_angl[1])
  // 3  dthet1+dthet2+dthet3 (chang2=mie_angl[2])
  // 4  dthet1+dthet2+dthet3+dthet4 (chang3=mie_angl[3])
  // 5  dthet1+dthet2+dthet3+dthet4+dthet5 (chang4=mie_angl[4])
  if(mie_n_angl < 6)
  {
    fprintf(stderr,"%s: error, mie_n_angl=%d < 6\n",fnam,mie_n_angl);
    return -1;
  }
  if(fabs(mie_angl[0]) > EPSILON)
  {
    fprintf(stderr,"%s: error, mie_angl[0]=%13.4e != 0.0\n",fnam,mie_angl[0]);
    return -1;
  }
  if(fabs(mie_angl[mie_n_angl-1]-180.0) > EPSILON)
  {
    fprintf(stderr,"%s: error, mie_angl[%d]=%13.4e != 180.0\n",fnam,mie_n_angl-1,mie_angl[mie_n_angl-1]);
    return -1;
  }
  chang_array = (double *)malloc(mie_n_angl*sizeof(double));
  dthet_array = (double *)malloc(mie_n_angl*sizeof(double));
  if(chang_array==NULL || dthet_array==NULL)
  {
    fprintf(stderr,"%s: failed in allocating memory\n",fnam);
    return -1;
  }

  nchang = 0;
  dthet_array[0] = mie_angl[1]-mie_angl[0];
  for(i=2; i<mie_n_angl; i++)
  {
    dtmp = mie_angl[i]-mie_angl[i-1];
    if(fabs(dtmp-dthet_array[nchang]) > EPSILON)
    {
      chang_array[nchang] = mie_angl[i-1];
      nchang++;
      dthet_array[nchang] = dtmp;
    }
  }

  if(nchang > 4)
  {
    fprintf(stderr,"%s: error, nchang=%d\n",fnam,nchang);
    return -1;
  } else
  if(nchang == 0)
  {
    dthet[0] = dthet_array[0];
    dthet[1] = dthet_array[0];
    dthet[2] = dthet_array[0];
    dthet[3] = dthet_array[0];
    dthet[4] = dthet_array[0];
    chang[0] = mie_angl[mie_n_angl-1];
    chang[1] = mie_angl[mie_n_angl-1];
    chang[2] = mie_angl[mie_n_angl-1];
    chang[3] = mie_angl[mie_n_angl-1];
  } else
  if(nchang == 1)
  {
    dthet[0] = dthet_array[0];
    dthet[1] = dthet_array[1];
    dthet[2] = dthet_array[1];
    dthet[3] = dthet_array[1];
    dthet[4] = dthet_array[1];
    chang[0] = chang_array[0];
    chang[1] = mie_angl[mie_n_angl-1];
    chang[2] = mie_angl[mie_n_angl-1];
    chang[3] = mie_angl[mie_n_angl-1];
  } else
  if(nchang == 2)
  {
    dthet[0] = dthet_array[0];
    dthet[1] = dthet_array[1];
    dthet[2] = dthet_array[2];
    dthet[3] = dthet_array[2];
    dthet[4] = dthet_array[2];
    chang[0] = chang_array[0];
    chang[1] = chang_array[1];
    chang[2] = mie_angl[mie_n_angl-1];
    chang[3] = mie_angl[mie_n_angl-1];

  } else
  if(nchang == 3)
  {
    dthet[0] = dthet_array[0];
    dthet[1] = dthet_array[1];
    dthet[2] = dthet_array[2];
    dthet[3] = dthet_array[3];
    dthet[4] = dthet_array[3];
    chang[0] = chang_array[0];
    chang[1] = chang_array[1];
    chang[2] = chang_array[2];
    chang[3] = mie_angl[mie_n_angl-1];
  }
  else
  {
    dthet[0] = dthet_array[0];
    dthet[1] = dthet_array[1];
    dthet[2] = dthet_array[2];
    dthet[3] = dthet_array[3];
    dthet[4] = dthet_array[4];
    chang[0] = chang_array[0];
    chang[1] = chang_array[1];
    chang[2] = chang_array[2];
    chang[3] = chang_array[3];
  }

  free(chang_array);
  free(dthet_array);

  return 0;
}

int ReadTape8(void)
{
  int i,j,k,l,m,n;
  int nnprob;
  int iprob,kprob;
  int na,ns;
  int v[10];
  double r1,r2;
  double wcom,rmod,lsgm;
  double wlen,refr,refi;
  double aext,asca,asym;
  double angl;
  double epsilon = 1.0e-3;
  double w[10];
  char line[MAXLINE];
  char str[10][MAXLINE];
  char *p;
  FILE *fp;
  char fnam[] = "ReadTape8";

// Tape8 example
//    2    1    0
//    1    6   WAVELENGTHS IN UM AND COEFFICIENTS IN KM-1
//  component#0
//      0.50     10.00      2.00    170.00      0.50  NO. OF ANGLES =  121
// 1.000E-02 2.000E-01 1.000E+04  799    1    2    3  component#0
// 0.1000000E+01  0.2700000E-01  0.3500000E+00  0.000E+00 1.000E+00 1.000E-01
// WAVELENGTH = 0.3000000E+00 UM, INDEX OF REF. =  0.1530000E+01 -0.8000000E-02  component#0
// EXT,SCAT,ABS,ASYM = 1.14656E-05  1.10206E-05  4.44995E-07 0.66823
// 0.00000E+00  7.56164E-01  7.56164E-01  7.56164E-01  0.00000E+00
// 5.00000E-01  7.55995E-01  7.55984E-01  7.55989E-01 -1.83537E-05
// 1.00000E+00  7.55487E-01  7.55443E-01  7.55465E-01 -7.33457E-05
// ...

  if((fp=fopen("mie2new.tp8","r")) == NULL)
  {
    fprintf(stderr,"%s: cannot open %s\n",fnam,"mie2new.tp8");
    return -1;
  }
  // line #1
  i = 0;
  if(fgets(line,MAXLINE,fp) == NULL)
  {
    fprintf(stderr,"%s: read error (line 1 %d)\n",fnam,i);
    fclose(fp);
    return -1;
  }
  if(sscanf(line,"%s",str[0]) != 1)
  {
    fprintf(stderr,"%s: scan error (line 1 %d) >>> %s\n",fnam,i,line);
    fclose(fp);
    return -1;
  }
  errno = 0;
  nnprob = strtol(str[0],&p,10);
  if(errno==ERANGE || *p!='\0')
  {
    fprintf(stderr,"%s: convert error (line 1 %d) >>> %s\n",fnam,i,line);
    fclose(fp);
    return -1;
  }
  if(nnprob != mie_n_comp)
  {
    fprintf(stderr,"%s: error (line 1 %d), nnprob=%d, mie_n_comp=%d\n",fnam,i,nnprob,mie_n_comp);
    fclose(fp);
    return -1;
  }
  for(n=0; n<mie_n_comp; n++)
  {
    // line #2
    i++;
    if(fgets(line,MAXLINE,fp) == NULL)
    {
      fprintf(stderr,"%s: read error (line 2 %d)\n",fnam,i);
      fclose(fp);
      return -1;
    }
    if(sscanf(line,"%s%s",str[0],str[1]) != 2)
    {
      fprintf(stderr,"%s: scan error (line 2 %d) >>> %s\n",fnam,i,line);
      fclose(fp);
      return -1;
    }
    for(j=0; j<2; j++)
    {
      errno = 0;
      v[j] = strtol(str[j],&p,10);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error (line 2 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
    }
    iprob = v[0];
    kprob = v[1];
    if(iprob != n+1)
    {
      fprintf(stderr,"%s: error (line 2 %d), iprob=%d, n=%d\n",fnam,i,iprob,n);
      fclose(fp);
      return -1;
    }
    if(kprob != mie_n_wlen)
    {
      fprintf(stderr,"%s: error (line 2 %d), kprob=%d, mie_n_wlen=%d\n",fnam,i,kprob,mie_n_wlen);
      fclose(fp);
      return -1;
    }
    // line #3
    i++;
    if(fgets(line,MAXLINE,fp) == NULL)
    {
      fprintf(stderr,"%s: read error (line 3 %d)\n",fnam,i);
      fclose(fp);
      return -1;
    }
    // line #4
    i++;
    if(fgets(line,MAXLINE,fp) == NULL)
    {
      fprintf(stderr,"%s: read error (line 4 %d)\n",fnam,i);
      fclose(fp);
      return -1;
    }
    for(j=0; j<MAXLINE; j++)
    {
      if(line[j] == '\0') break;
      if(line[j] == '=') line[j] = ' ';
    }
    if(sscanf(line,"%s%s%s%s%s%s%s%s%s",str[0],str[1],str[2],str[3],str[4],str[5],str[6],str[7],str[8]) != 9)
    {
      fprintf(stderr,"%s: scan error (line 4 %d) >>> %s\n",fnam,i,line);
      fclose(fp);
      return -1;
    }
    errno = 0;
    na = strtol(str[8],&p,10);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"%s: convert error (line 4 %d) >>> %s\n",fnam,i,line);
      fclose(fp);
      return -1;
    }
    if(na != mie_n_angl)
    {
      fprintf(stderr,"%s: error (line 4 %d), na=%d, mie_n_angl=%d\n",fnam,i,na,mie_n_angl);
      fclose(fp);
      return -1;
    }
    // line #5
    i++;
    if(fgets(line,MAXLINE,fp) == NULL)
    {
      fprintf(stderr,"%s: read error (line 5 %d)\n",fnam,i);
      fclose(fp);
      return -1;
    }
    if(sscanf(line,"%s%s%s%s",str[0],str[1],str[2],str[3]) != 4)
    {
      fprintf(stderr,"%s: scan error (line 5 %d) >>> %s\n",fnam,i,line);
      fclose(fp);
      return -1;
    }
    for(j=0; j<2; j++)
    {
      errno = 0;
      w[j] = strtod(str[j],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error (line 5 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
    }
    r1 = w[0];
    r2 = w[1];
    if(cnt_2peak)
    {
      if((fabs(r1/mie_rmin_com[0]-1.0)>epsilon && fabs(r1/mie_rmin_com[1]-1.0)>epsilon) ||
         (fabs(r2/mie_rmax_com[0]-1.0)>epsilon && fabs(r2/mie_rmax_com[1]-1.0)>epsilon))
      {
        fprintf(stderr,"%s: error (line 5 %d), r1=%13.4e, mie_rmin_com[%d]=%13.4e, r2=%13.4e, mie_rmax_com[%d]=%13.4e\n",fnam,
                        i,r1,n,mie_rmin_com[n],r2,n,mie_rmax_com[n]);
        fclose(fp);
        return -1;
      }
    }
    else
    {
      if(fabs(r1/mie_rmin_com[n]-1.0)>epsilon || fabs(r2/mie_rmax_com[n]-1.0)>epsilon)
      {
        fprintf(stderr,"%s: error (line 5 %d), r1=%13.4e, mie_rmin_com[%d]=%13.4e, r2=%13.4e, mie_rmax_com[%d]=%13.4e\n",fnam,
                        i,r1,n,mie_rmin_com[n],r2,n,mie_rmax_com[n]);
        fclose(fp);
        return -1;
      }
    }
    errno = 0;
    ns = strtol(str[3],&p,10);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"%s: convert error (line 5 %d) >>> %s\n",fnam,i,line);
      fclose(fp);
      return -1;
    }
    if(ns != mie_nstp)
    {
      fprintf(stderr,"%s: error (line 5 %d), ns=%d, mie_nstp=%d\n",fnam,i,ns,mie_nstp);
      fclose(fp);
      return -1;
    }
    // line #6
    i++;
    if(fgets(line,MAXLINE,fp) == NULL)
    {
      fprintf(stderr,"%s: read error (line 6 %d)\n",fnam,i);
      fclose(fp);
      return -1;
    }
    if(cnt_2peak)
    {
      if(sscanf(line,"%s%s%s%s%s%s",str[0],str[1],str[2],str[3],str[4],str[5]) != 6)
      {
        fprintf(stderr,"%s: scan error (line 6 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
      for(j=0; j<6; j++)
      {
        errno = 0;
        w[j] = strtod(str[j],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error (line 6 %d) >>> %s\n",fnam,i,line);
          fclose(fp);
          return -1;
        }
      }
      if(fabs(w[0]-mie_wcom_com[0])>DELTA ||
         fabs(w[1]-mie_rmod_com[0])>DELTA ||
         fabs(w[2]-mie_lsgm_com[0])>DELTA ||
         fabs(w[3]-mie_wcom_com[1])>DELTA ||
         fabs(w[4]-mie_rmod_com[1])>DELTA ||
         fabs(w[5]-mie_lsgm_com[1])>DELTA)
      {
        fprintf(stderr,"%s: error (line 6 %d)\n",fnam,i);
        fprintf(stderr,"  w[0]=%13.4e, mie_wcom_com[%d]=%13.4e, diff=%13.4e\n",w[0],0,mie_wcom_com[0],fabs(w[0]-mie_wcom_com[0]));
        fprintf(stderr,"  w[1]=%13.4e, mie_rmod_com[%d]=%13.4e, diff=%13.4e\n",w[1],0,mie_rmod_com[0],fabs(w[1]-mie_rmod_com[0]));
        fprintf(stderr,"  w[2]=%13.4e, mie_lsgm_com[%d]=%13.4e, diff=%13.4e\n",w[2],0,mie_lsgm_com[0],fabs(w[2]-mie_lsgm_com[0]));
        fprintf(stderr,"  w[3]=%13.4e, mie_wcom_com[%d]=%13.4e, diff=%13.4e\n",w[3],1,mie_wcom_com[1],fabs(w[3]-mie_wcom_com[1]));
        fprintf(stderr,"  w[4]=%13.4e, mie_rmod_com[%d]=%13.4e, diff=%13.4e\n",w[4],1,mie_rmod_com[1],fabs(w[4]-mie_rmod_com[1]));
        fprintf(stderr,"  w[5]=%13.4e, mie_lsgm_com[%d]=%13.4e, diff=%13.4e\n",w[5],1,mie_lsgm_com[1],fabs(w[5]-mie_lsgm_com[1]));
        fclose(fp);
        return -1;
      }
    }
    else
    {
      if(sscanf(line,"%s%s%s",str[0],str[1],str[2]) != 3)
      {
        fprintf(stderr,"%s: scan error (line 6 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
      for(j=0; j<3; j++)
      {
        errno = 0;
        w[j] = strtod(str[j],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error (line 6 %d) >>> %s\n",fnam,i,line);
          fclose(fp);
          return -1;
        }
      }
      wcom = w[0];
      rmod = w[1];
      lsgm = w[2];
      if(fabs(wcom-1.0)>DELTA || fabs(rmod-mie_rmod_com[n])>DELTA || fabs(lsgm-mie_lsgm_com[n])>DELTA)
      {
        fprintf(stderr,"%s: error (line 6 %d)\n",fnam,i);
        fprintf(stderr,"  wcom=%13.4e,        expected=%13.4e, diff=%13.4e\n",wcom,1.0,fabs(wcom-1.0));
        fprintf(stderr,"  rmod=%13.4e, mie_rmod_com[%d]=%13.4e, diff=%13.4e\n",rmod,n,mie_rmod_com[n],fabs(rmod-mie_rmod_com[n]));
        fprintf(stderr,"  lsgm=%13.4e, mie_lsgm_com[%d]=%13.4e, diff=%13.4e\n",lsgm,n,mie_lsgm_com[n],fabs(lsgm-mie_lsgm_com[n]));
        fclose(fp);
        return -1;
      }
    }
    for(k=0; k<mie_n_wlen; k++)
    {
      // line #7
      i++;
      if(fgets(line,MAXLINE,fp) == NULL)
      {
        fprintf(stderr,"%s: read error (line 7 %d)\n",fnam,i);
        fclose(fp);
        return -1;
      }
      for(j=0; j<MAXLINE; j++)
      {
        if(line[j] == '\0') break;
        if(line[j] == '=') line[j] = ' ';
      }
      if(sscanf(line,"%s%s%s%s%s%s%s%s",str[0],str[1],str[2],str[3],str[4],str[5],str[6],str[7]) != 8)
      {
        fprintf(stderr,"%s: scan error (line 7 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
      errno = 0;
      wlen = strtod(str[1],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error (line 7 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
      if(fabs(wlen-mie_wlen_um[k]) > DELTA)
      {
        fprintf(stderr,"%s: error (line 7 %d), wlen=%13.4e, mie_wlen_um[%d]=%13.4e\n",fnam,i,wlen,k,mie_wlen_um[k]);
        fclose(fp);
        return -1;
      }
      errno = 0;
      refr = strtod(str[6],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error (line 7 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
      if(fabs(refr-mie_refr_com[n][k]) > DELTA)
      {
        fprintf(stderr,"%s: error (line 7 %d), refr=%13.4e, mie_refr_com[%d][%d]=%13.4e\n",fnam,i,refr,n,k,mie_refr_com[n][k]);
        fclose(fp);
        return -1;
      }
      errno = 0;
      refi = strtod(str[7],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error (line 7 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
      if(fabs(fabs(refi)-fabs(mie_refi_com[n][k])) > DELTA)
      {
        fprintf(stderr,"%s: error (line 7 %d), refi=%13.4e, mie_refi_com[%d][%d]=%13.4e\n",fnam,i,refi,n,k,mie_refi_com[n][k]);
        fclose(fp);
        return -1;
      }
      // line #8
      i++;
      if(fgets(line,MAXLINE,fp) == NULL)
      {
        fprintf(stderr,"%s: read error (line 8 %d)\n",fnam,i);
        fclose(fp);
        return -1;
      }
      for(j=0; j<MAXLINE; j++)
      {
        if(line[j] == '\0') break;
        if(line[j] == '=') line[j] = ' ';
      }
      if(sscanf(line,"%s%s%s%s%s",str[0],str[1],str[2],str[3],str[4]) != 5)
      {
        fprintf(stderr,"%s: scan error (line 8 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
      errno = 0;
      aext = strtod(str[1],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error (line 8 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
      errno = 0;
      asca = strtod(str[2],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error (line 8 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
      errno = 0;
      asym = strtod(str[4],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error (line 8 %d) >>> %s\n",fnam,i,line);
        fclose(fp);
        return -1;
      }
      mie_aext_com[n][k] = aext*1.0e3; // 1/km -> 1/Mm
      mie_asca_com[n][k] = asca*1.0e3; // 1/km -> 1/Mm
      mie_asym_com[n][k] = asym;
      for(l=0; l<mie_n_angl; l++)
      {
        i++;
        if(fgets(line,MAXLINE,fp) == NULL)
        {
          fprintf(stderr,"%s: read error (line %d)\n",fnam,i);
          fclose(fp);
          return -1;
        }
        if(sscanf(line,"%s%s%s",str[0],str[1],str[2]) != 3)
        {
          fprintf(stderr,"%s: scan error (line %d) >>> %s\n",fnam,i,line);
          fclose(fp);
          return -1;
        }
        for(j=0; j<3; j++)
        {
          errno = 0;
          w[j] = strtod(str[j],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error (line %d) >>> %s\n",fnam,i,line);
            fclose(fp);
            return -1;
          }
        }
        angl = w[0];
        if(fabs(angl-mie_angl[l]) > DELTA)
        {
          fprintf(stderr,"%s: error (line %d), angl=%13.4e, mie_angl[%d]=%13.4e\n",fnam,i,angl,l,mie_angl[l]);
          fclose(fp);
          return -1;
        }
        m = mie_n_angl*k+l;
        mie_phs1_com[n][m] = w[1];
        mie_phs2_com[n][m] = w[2];
        mie_phas_com[n][m] = 0.5*(w[1]+w[2]);
      }
    }
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
    {"mie_wsgm",1,0,'S'},
    {"cnt_2peak",0,0,'D'},
    {"cnt_norun",0,0,'N'},
  };

  opterr = 0;
  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,"-x:X:n:S:DN",long_options,&option_index);
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
        if(errno!=ERANGE && *endp=='\0' && ntmp>0 && ntmp<=MAXSTEP) mie_nstp = ntmp;
        else
        {
          fprintf(stderr,"#Steps -> out of range %s\n",optarg);
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
      case 'D':
        cnt_2peak = 1;
        break;
      case 'N':
        cnt_norun = 1;
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
  char optstring[] = "fgwlLaWAPQiIruxXnSDNmoOdvh";

  CommonUsage("aeros_mie2new",0);
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
      case 'S':
        fprintf(stderr," S -mie_wsgm      |%s|%s|%s| %f\n",As(e,"Log10(R) width",n), As(a,"sigma",n),      Af(d,MIE_WSGM,n),mie_wsgm);
        break;
      case 'D':
        fprintf(stderr," D -cnt_2peak     |%s|%s|%s| %d\n",As(e,"2 peak mode",n),    As(a,"nothing",n),    Ad(d,0,n),cnt_2peak);
        break;
      case 'N':
        fprintf(stderr," N -cnt_norun     |%s|%s|%s| %d\n",As(e,"No run  mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_norun);
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
  CommonUsage(NULL,5);

  return 0;
}
