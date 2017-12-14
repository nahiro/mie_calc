#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include "instrument.h"

// Common constants
#define	NONAME				"NONE"			// No name
#define	NODATA				0			// No data
#define	INVALID				-1			// Invalid value
#define	MAXENTR				20			// Max #entries in a line
#define	MAXLINE				256			// Max #chars in a line
#define	EPSILON				1.0e-14			// A small number
#define	DELTA				1.0e-13			// A small number
// Constants for input data
#define	INP_MAXDATA			100000			// Max #input data
#define	INP_NCOL			2			// #Columns
// Constants for Simulation
#define	SIM_MAXDATA			20000			// #Data
#define	SIM_MAXCONV			1000			// #Data
#define	SIM_NSAM			4			// #Sampling
#define	SIM_XMGN			16.0			// X margin
#define	SIM_XSGM			4.0			// X sigma
#define	SIM_WSGM			20.0			// X width

int	inp_n_data			= NODATA;
int	inp_wnum_min			= INVALID;
int	inp_wnum_max			= INVALID;
int	inp_sbmp			= INVALID;
int	inp_wnum[INP_MAXDATA];
double	inp_wlen[INP_MAXDATA];
double	inp_data[INP_MAXDATA];
// Parameters for Simulation
int	sim_nsam			= SIM_NSAM;		// #Sampling
double	sim_xmgn			= SIM_XMGN;		// X margin
double	sim_xsgm			= SIM_XSGM;		// X sigma
double	sim_wsgm			= SIM_WSGM;		// X width in sigma
double	sim_wlen[SIM_MAXDATA];
double	sim_data[SIM_MAXDATA];
double	sim_wlen_cnv[SIM_MAXCONV];
double	sim_data_cnv[SIM_MAXCONV];

int Compare(void);
int Init(void);
int ReadInput(void);
int GetAverage(double wlen,double *v);
int GetGSample(double wlen,double *v);
int Convolute(double xmin,double xmax,double xstp,double xsgm,double wsgm,
              int ninp,const double *xinp,const double *yinp,
              int nout,double *xout,double *yout);
double Gauss(double g,double f0,double f);
int Sampling(int ninp,const double *xinp,const double *yinp,
             int nout,const double *xout,double *yout,double yuni);

int main(int argc,char **argv)
{
  if(Init() < 0) return -1;
  if(Compare() < 0) return -1;

  return 0;
}

int Compare(void)
{
  int i,j;
  double a,s;

  for(i=0; i<ins_n_inst; i++)
  {
    for(j=0; j<ins_n_wlen; j++)
    {
      if(GetAverage(ins_wlen[i][j],&a) < 0) return -1;
      if(GetGSample(ins_wlen[i][j],&s) < 0) return -1;
      printf("%d %13.4f %13.6e %13.6e %13.6e\n",i,ins_wlen[i][j],a,s,s/a);
    }
  }

  return 0;
}

int Init(void)
{
  if(ReadInput() < 0) return -1;

  return 0;
}

int ReadInput(void)
{
  int i,n;
  int err;
  int wnum;
  int sbmp;
  double wval;
  double w[INP_NCOL];
  double epsilon = 1.0e-1;
  char line[MAXLINE];
  char str[INP_NCOL][MAXLINE];
  char *p;

  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,stdin) != NULL)
  {
    if(sscanf(line,"%s%s",str[0],str[1]) != INP_NCOL)
    {
      fprintf(stderr,"Warning, read error >>> %s\n",line);
      continue;
    }
    for(n=0; n<INP_NCOL; n++)
    {
      errno = 0;
      w[n] = strtod(str[n],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"Convert error >>> %s\n",line);
        err = 1;
        break;
      }
    }
    if(err) break;
    wval = 1.0e7/w[0];
    wnum = (int)(wval+0.5);
    if(fabs(wval-(double)wnum) > epsilon)
    {
      fprintf(stderr,"Error, wval=%.4f wnum=%d\n",wval,wnum);
      err = 1;
      break;
    }
    if(i >= INP_MAXDATA)
    {
      fprintf(stderr,"Warning, #data exceed the limit >>> %d\n",i);
      break;
    }
    inp_wnum[i] = wnum;
    inp_wlen[i] = w[0];
    inp_data[i] = w[1];
    i++;
  }
  inp_n_data = i;
  fprintf(stderr,"%d data have been read.\n",inp_n_data);
  if(inp_n_data < 2)
  {
    fprintf(stderr,"#data is not enough.\n");
    err = 1;
  }
  if(err)
  {
    return -1;
  }

  inp_wnum_min = +RAND_MAX;
  inp_wnum_max = -RAND_MAX;
  for(i=0; i<inp_n_data; i++)
  {
    if(inp_wnum[i] < inp_wnum_min) inp_wnum_min = inp_wnum[i];
    if(inp_wnum[i] > inp_wnum_max) inp_wnum_max = inp_wnum[i];
  }
  for(i=1; i<inp_n_data; i++)
  {
    if(inp_sbmp < 0)
    {
      inp_sbmp = abs(inp_wnum[i]-inp_wnum[i-1]);
    }
    else
    {
      sbmp = abs(inp_wnum[i]-inp_wnum[i-1]);
      if(sbmp != inp_sbmp)
      {
        fprintf(stderr,"Error, sbmp=%d, inp_sbmp=%d\n",sbmp,inp_sbmp);
        return -1;
      }
    }
  }
  fprintf(stderr,"nmin: %10d\n",inp_wnum_min);
  fprintf(stderr,"nmax: %10d\n",inp_wnum_max);
  fprintf(stderr,"sbmp: %10d\n",inp_sbmp);

  return 0;
}

int GetAverage(double wlen,double *v)
{
  int i;
  int i1,i2;
  int n1,n2,ns;
  double nd;
  double s1;

  nd = (double)sim_nsam;
  n1 = inp_sbmp*((int)((1.0e7/(wlen*inp_sbmp)-0.5*(nd-1.0))+0.5));
  n2 = inp_sbmp*((int)(nd+0.5))+n1-inp_sbmp;

  i1 = (n1-inp_wnum_min)/inp_sbmp-1;
  i2 = (n2-inp_wnum_min)/inp_sbmp+1;
  if(i1 < 0) i1 = 0;
  if(i1 >= inp_n_data) i1 = inp_n_data-1;
  if(i2 < 0) i2 = 0;
  if(i2 >= inp_n_data) i2 = inp_n_data-1;
  ns = 0;
  s1 = 0.0;
  for(i=i1; i<=i2; i++)
  {
    if(inp_wnum[i]>=n1 && inp_wnum[i]<=n2)
    {
      s1 += inp_data[i];
      ns++;
    }
  }
  if(ns != sim_nsam)
  {
    fprintf(stderr,"GetAverage: error, sim_nsam=%d, ns=%d\n",sim_nsam,ns);
    return -1;
  }
  if(ns > 1)
  {
    s1 /= (double)ns;
  }
  *v = s1;

  return 0;
}

int GetGSample(double wlen,double *v)
{
  int n;
 
  if((n=Convolute(wlen-sim_xmgn,wlen+sim_xmgn,1.0,sim_xsgm,sim_wsgm,
                  inp_n_data,inp_wlen,inp_data,SIM_MAXCONV,sim_wlen_cnv,sim_data_cnv)) < 0)
  {
    fprintf(stderr,"GetGSample: error in Convolute().\n");
    return -1;
  }
  if(Sampling(n,sim_wlen_cnv,sim_data_cnv,1,&wlen,v,1.0) < 0)
  {
    fprintf(stderr,"GetGSample: error in Sampling().\n");
    return -1;
  }

  return 0;
}

int Convolute(double xmin,double xmax,double xstp,double xsgm,double wsgm,
              int ninp,const double *xinp,const double *yinp,
              int nout,double *xout,double *yout)
{
  int i,j,k;
  int nstp;
  int nwid;
  double x1,x2;
  double norm;
  double offset = 0.5000000000001;

  if(isnan(xmin)) xmin = xinp[0];
  if(isnan(xmax)) xmax = xinp[ninp-1];
  nstp = (int)((xmax-xmin)/xstp+0.5);
  if(nstp < 1)
  {
    fprintf(stderr,"Convolute: error, invalid #steps (xmin=%13.4e xmax=%13.4e xstp=%13.4e nstp=%d\n",
                    xmin,xmax,xstp,nstp);
    return -1;
  } else
  if(nstp >= nout)
  {
    fprintf(stderr,"Convolute: warning, nstp=%d -> changed to %d\n",nstp,nout-1);
    nstp = nout-1;
  }
  nwid = (int)(xsgm*wsgm/xstp+0.5);
  if(nwid < 1)
  {
    fprintf(stderr,"Error, nwid=%d, xsgm=%13.4e wsgm=%13.4e xstp=%13.4e\n",nwid,xsgm,wsgm,xstp);
    return -1;
  }

  x1 = xmin-xsgm*wsgm;
  x2 = xmax+xsgm*wsgm;
  for(i=0; i<=nstp; i++)
  {
    xout[i] = xmin+xstp*i;
    yout[i] = 0.0;
  }
  for(i=0; i<ninp; i++)
  {
    if(xinp[i] < x1) continue;
    if(xinp[i] > x2) continue;
    norm = yinp[i]*(i==ninp-1?fabs(xinp[i-1]-xinp[i]):fabs(xinp[i+1]-xinp[i]));
    j = (int)((xinp[i]-xmin)/xstp+offset);
    for(k=j-nwid; k<=j+nwid; k++)
    {
      if(k < 0)
      {
        continue;
      } else
      if(k > nstp)
      {
        break;
      }
      yout[k] += norm*Gauss(xsgm,xinp[i],xout[k]);
    }
  }

  return nstp+1;
}

double Gauss(double s,double x0,double x)
{
  return exp(-0.5*(x-x0)*(x-x0)/(s*s))/(s*sqrt(2.0*M_PI));
}

int Sampling(int ninp,const double *xinp,const double *yinp,
             int nout,const double *xout,double *yout,double yuni)
{
  int i,i1,i2;
  double xstp;
  double epsilon = 1.0e-12;

  // check input
  if(ninp < 2)
  {
    fprintf(stderr,"Sampling: error, ninp=%d\n",ninp);
    return -1;
  }
  xstp = xinp[1]-xinp[0];
  for(i=1; i<ninp; i++)
  {
    if(fabs(xinp[i]-(xinp[0]+xstp*i)) > epsilon)
    {
      fprintf(stderr,"Sampling: error, xinp[0]=%13.4e, xinp[%d]=%13.4e, xstp=%13.4e, diff=%13.4e\n",
                      xinp[0],i,xinp[i],xstp,xinp[i]-(xinp[0]+xstp*i));
      return -1;
    }
  }

  // interpolate
  for(i=0; i<nout; i++)
  {
    i1 = (int)((xout[i]-xinp[0])/xstp+epsilon);
    if(i1 < 0) i1 = 0;
    i2 = i1+1;
    if(i2 >= ninp)
    {
      i2 = ninp-1;
      i1 = i2-1;
    }
    yout[i] = (yinp[i1]+(yinp[i2]-yinp[i1])*(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1]))*yuni;
  }

  return 0;
}
