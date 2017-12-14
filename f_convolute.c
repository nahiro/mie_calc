/*********************************************************/
/* convolute                                             */
/* Author: N.Manago Apr,04,2008                          */
/* $Revision: 1.2 $                                      */
/* $Date: 2008-07-03 08:29:11 $ (UTC)                    */
/*********************************************************/
// required headers
// #include <stdio.h>
// #include <math.h>
// #include <bits/nan.h>

int Convolute(double xmin,double xmax,double xstp,double xsgm,double wsgm,
              int ninp,const double *xinp,const double *yinp,
              int nout,double *xout,double *yout);
double gauss(double g,double f0,double f);

int Convolute(double xmin,double xmax,double xstp,double xsgm,double wsgm,
              int ninp,const double *xinp,const double *yinp,
              int nout,double *xout,double *yout)
{
  int i,j,k;
  int nstp;
  int nwid;
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

  for(i=0; i<=nstp; i++)
  {
    xout[i] = xmin+xstp*i;
    yout[i] = 0.0;
  }
  for(i=0; i<ninp; i++)
  {
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
      yout[k] += norm*gauss(xsgm,xinp[i],xout[k]);
    }
  }

  return nstp+1;
}

double gauss(double s,double x0,double x)
{
  return exp(-0.5*(x-x0)*(x-x0)/(s*s))/(s*sqrt(2.0*M_PI));
}
