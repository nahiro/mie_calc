/*********************************************************/
/* sampling                                              */
/* Author: N.Manago Apr,04,2008                          */
/* $Revision: 1.2 $                                      */
/* $Date: 2008-07-03 17:06:46 $ (UTC)                    */
/*********************************************************/
// required headers
// #include <stdio.h>
// #include <math.h>
// #include <bits/nan.h>

int Sampling(int ninp,const double *xinp,const double *yinp,
             int nout,const double *xout,double *yout);
int Sampling(int ninp,const double *xinp,const double *yinp,
             int nout,const double *xout,double *yout)
{
  int i,i1,i2;
  double xstp;
  double epsilon = 1.0e-13;

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
    yout[i] = yinp[i1]+(yinp[i2]-yinp[i1])*(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1]);
  }

  return 0;
}
