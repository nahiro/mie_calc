/*********************************************************/
/* MS720_FOV                                             */
/* Author: N.Manago Dec,02,2008                          */
/* $Revision: 46 $                                       */
/* $Date: 2008-12-03 03:07:22 +0900 (Wed, 03 Dec 2008) $ */
/*********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TMarker.h>
#include <TMinuit.h>
#include <TPolyMarker.h>
#include <TPolyLine.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TPostScript.h>
#include <TStyle.h>
#include "strutil.h"

#define		PI			3.141592653589793	// PI
#define		PI2			6.283185307179586	// 2*PI
#define		PI_2			1.570796326794896	// PI/2
#define		R_TO_D			57.29577951308232	// 180/PI rad -> deg
#define		D_TO_R			1.745329251994329e-02	// PI/180 deg -> rad
#define		NCOL			8			// #Columns
#define		NCOL2			10			// #Columns
#define		NDVTH			4			// #Bins in Theta
#define		NDVPH			16			// #Bins in Phi
#define		XMIN			NAN			// X min
#define		XMAX			NAN			// X max
#define		YMIN			NAN			// Y min
#define		YMAX			NAN			// Y max
#define		ZMIN			NAN			// Z min
#define		ZMAX			NAN			// Z max
#define		ZSTP			100.0			// Z step
#define		THMIN			0.0			// Theta min
#define		THMAX			(10.0*D_TO_R)		// Theta max
#define		PHMIN			0.0			// Phi min
#define		PHMAX			PI2			// Phi max
#define		PHOFS			PI_2			// Phi offset
#define		THCTR			0.0			// Theta center
#define		PHCTR			0.0			// Phi   center
#define		XOFS			0.62			// Pad X offset
#define		YOFS			0.635			// Pad Y offset
#define		XMGN			0.3			// X Margin
#define		PWID			0.36			// Pad width
#define		CWID			800			// Canvas width
#define		MAXDATA			20000			// Max #points
#define		MAXENTR			20			// Max #columns
#define		MAXNPAR			10			// Max #Parameters
#define		MAXLINE			256			// Max #chars in a line
#define		EPSILON			1.0e-14			// A small number
#define		DELTA			1.0e-3			// A small number
#define		MAX(a,b)		((b)>(a)?(b):(a))

Int_t		ndvth			= NDVTH;		// #Bins
Int_t		ndvph			= NDVPH;		// #Bins
Double_t	xmin			= XMIN;			// X min
Double_t	xmax			= XMAX;			// X max
Double_t	ymin			= YMIN;			// Y min
Double_t	ymax			= YMAX;			// Y max
Double_t	zmin			= ZMIN;			// Z min
Double_t	zmax			= ZMAX;			// Z max
Double_t	zstp			= ZSTP;			// Z step
Double_t	thmin			= THMIN;		// Theta min
Double_t	thmax			= THMAX;		// Theta max
Double_t	phmin			= PHMIN;		// Phi min
Double_t	phmax			= PHMAX;		// Phi max
Double_t	phofs			= PHOFS;		// Phi offset
Double_t	thctr			= THCTR;		// Theta center
Double_t	phctr			= PHCTR;		// Phi   center
Double_t	xofs			= XOFS;			// Pad X offset
Double_t	yofs			= YOFS;			// Pad Y offset
Double_t	xmgn			= XMGN;			// X margin
Double_t	pwid			= PWID;			// Pad width
Int_t		cwid			= CWID;			// Canvas width
Int_t		cm			= 0;			// Circle  mode
Int_t		gm			= 0;			// Gray  mode
Int_t		bm			= 0;			// Batch   mode
Int_t		vb			= 0;			// Verbose mode
Int_t		db			= 0;			// Debug   mode
Int_t		hp			= 0;			// Help    mode
char		fout[MAXLINE]		= "";			// Output psfile name
char		title[MAXLINE]		= "title";		// Graph title
char		xtit[MAXLINE]		= "X";			// X axis   title
char		ytit[MAXLINE]		= "Y";			// Y axis   title
char		ztit[MAXLINE]		= "";			// Z axis   title
Int_t		ndat;
Double_t	xdat[MAXDATA];
Double_t	ydat[MAXDATA];
Double_t	zdat[MAXDATA];

TCanvas *canv;
TApplication *appl;

int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);
int DrawGraph(void);
double range(double a);
double range2(double a);
double SepCosine(double x1,double y1,double z1,double x2,double y2,double z2);
double SepCosine2(double th1,double ph1,double th2,double ph2);
double SepAngle(double x1,double y1,double z1,double x2,double y2,double z2);
double SepAngle2(double th1,double ph1,double th2,double ph2);
int SepToAzi(double th,double sep,double *dph);
int PolToXyz(double th,double ph,double *x,double *y,double *z);
int XyzToPol(double x,double y,double z,double *th,double *ph);
void IjkToXyz(double th,double ph,
              double xi,double yi,double zi,
              double *x,double *y,double *z);
void XyzToIjk(double th,double ph,
              double xi,double yi,double zi,
              double *x,double *y,double *z);
int RGB(double min,double max,double val);

int main(int argc,char *argv[])
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  DrawGraph();

  return 0;
}

int DrawGraph(void)
{
  Int_t i,j,k;
  Int_t m,n;
  Int_t nbin = 1000;
  Double_t th,dth;
  Double_t ph,dph;
  Double_t th0,th1,th2,dth0;
  Double_t ph0,ph1,ph2,dph0;
  Double_t ths,phs;
  Double_t x,y,z;
  Double_t xi,yi,zi;
  Double_t x1,x2;
  Double_t y1,y2;
  Double_t xwid,ywid;
  Double_t px[2010];
  Double_t py[2010];
  char temp[MAXLINE];
  Float_t c;
  TPolyLine *line;
  TBox *box;
  TLatex *text;
  TMarker *mark;

  if(isnan(xmin) || isnan(xmax) || isnan(ymin) || isnan(ymax))
  {
    x1 = +1.0e100;
    x2 = -1.0e100;
    y1 = +1.0e100;
    y2 = -1.0e100;
    dth = (thmax-thmin)/nbin;
    dph = (phmax-phmin)/nbin;
    for(m=0; m<=nbin; m++)
    {
      th = thmin+dth*m;
      for(n=0; n<=nbin; n++)
      {
        ph = phmin+dph*n;
        PolToXyz(th,ph+phofs,&xi,&yi,&zi);
        if(cm)
        {
          x = xi;
          y = yi;
          z = zi;
        }
        else
        {
          IjkToXyz(thctr,phctr,xi,yi,zi,&x,&y,&z);
        }
        if(x < x1) x1 = x;
        if(x > x2) x2 = x;
        if(y < y1) y1 = y;
        if(y > y2) y2 = y;
      }
    }
  }
  if(isnan(xmin) && isnan(xmax) && isnan(ymin) && isnan(ymax))
  {
    xwid = x2-x1;
    ywid = y2-y1;
    xwid = MAX(xwid,ywid);
    xwid *= 1.02;
    ywid = xwid;
    xmin = 0.5*(x1+x2-xwid);
    xmax = xmin+xwid;
    ymin = 0.5*(y1+y2-ywid);
    ymax = ymin+ywid;
  }
  else
  {
    if(isnan(xmin)) xmin = x1;
    if(isnan(xmax)) xmax = x2;
    if(isnan(ymin)) ymin = y1;
    if(isnan(ymax)) ymax = y2;
  }
  line = new TPolyLine();
  line->SetLineColor(1);
  canv->Draw();
  canv->Range(xmin,ymin,xmax+xwid*xmgn,ymax);

  dth0 = (thmax-thmin)/ndvth;
  dph0 = (phmax-phmin)/ndvph;
  // Fill segments
  for(i=0,k=1; i<ndvth; i++)
  {
    th1 = thmin+dth0*i;
    th2 = thmin+dth0*(i+1);
    for(j=0; j<ndvph; j++,k++)
    {
      ph1 = phmin+dph0*j;
      ph2 = phmin+dph0*(j+1);
      // Fill segment
      dph = (ph2-ph1)/nbin;
      m = 0;
      th = th1;
      for(n=0; n<=nbin; n++)
      {
        ph = ph1+dph*n;
        PolToXyz(th,ph+phofs,&xi,&yi,&zi);
        if(cm)
        {
          x = xi;
          y = yi;
          z = zi;
        }
        else
        {
          IjkToXyz(thctr,phctr,xi,yi,zi,&x,&y,&z);
        }
        px[m] = x;
        py[m] = y;
        m++;
      }
      th = th2;
      for(n=nbin; n>=0; n--)
      {
        ph = ph1+dph*n;
        PolToXyz(th,ph+phofs,&xi,&yi,&zi);
        if(cm)
        {
          x = xi;
          y = yi;
          z = zi;
        }
        else
        {
          IjkToXyz(thctr,phctr,xi,yi,zi,&x,&y,&z);
        }
        px[m] = x;
        py[m] = y;
        m++;
      }
      px[m] = px[0];
      py[m] = py[0];
      if(gm)
      {
        if(zdat[k] < zmin)
        {
          line->SetFillColor(1);
        } else
        if(zdat[k] > zmax)
        {
          line->SetFillColor(10);
        }
        else
        {
          c = (Float_t)(0.3+0.7*(zdat[k]-zmin)/(zmax-zmin));
          line->SetFillColor(TColor::GetColor(c,c,c));
        }
      }
      else
      {
        if(zdat[k] < zmin)
        {
          line->SetFillColor(TColor::GetColor(200,200,200));
        }
        else
        {
          line->SetFillColor(RGB(zmin,zmax,zdat[k]));
        }
      }
      line->DrawPolyLine(m,px,py,"F");
    }
  }
  // Draw lines
  for(i=0,k=1; i<ndvth; i++)
  {
    th1 = thmin+dth0*i;
    th2 = thmin+dth0*(i+1);
    th0 = 0.5*(th1+th2);
    for(j=0; j<ndvph; j++,k++)
    {
      ph1 = phmin+dph0*j;
      ph2 = phmin+dph0*(j+1);
      ph0 = 0.5*(ph1+ph2);
      PolToXyz(th0,ph0,&xi,&yi,&zi);
      IjkToXyz(thctr,phctr,xi,yi,zi,&x,&y,&z);
      XyzToPol(x,y,z,&ths,&phs);
      fprintf(stderr,"%13.4f %13.4f %13.4f %13.4f %13.4f %13.6e\n",th0*R_TO_D,ph0*R_TO_D,
                      ths*R_TO_D,phs*R_TO_D,SepAngle2(thctr,phctr,ths,phs)*R_TO_D,zdat[k]);
      if(SepAngle2(xdat[k],ydat[k],th0,ph0) > DELTA)
      {
        fprintf(stderr,"Warning, th[%3d]=%13.4e ph[%3d]=%13.4e th0=%13.4e ph0=%13.4e\n",
                        k,xdat[k]*R_TO_D,k,ydat[k]*R_TO_D,th0*R_TO_D,ph0*R_TO_D);
      }
      // Draw line
      dph = (ph2-ph1)/nbin;
      m = 0;
      th = th1;
      for(n=0; n<=nbin; n++)
      {
        ph = ph1+dph*n;
        PolToXyz(th,ph+phofs,&xi,&yi,&zi);
        if(cm)
        {
          x = xi;
          y = yi;
          z = zi;
        }
        else
        {
          IjkToXyz(thctr,phctr,xi,yi,zi,&x,&y,&z);
        }
        px[m] = x;
        py[m] = y;
        m++;
      }
      th = th2;
      for(n=nbin; n>=0; n--)
      {
        ph = ph1+dph*n;
        PolToXyz(th,ph+phofs,&xi,&yi,&zi);
        if(cm)
        {
          x = xi;
          y = yi;
          z = zi;
        }
        else
        {
          IjkToXyz(thctr,phctr,xi,yi,zi,&x,&y,&z);
        }
        px[m] = x;
        py[m] = y;
        m++;
      }
      px[m] = px[0];
      py[m] = py[0];
      line->DrawPolyLine(m,px,py,"");
    }
  }
  // Draw representative point
  th = xdat[0];
  ph = ydat[0];
  PolToXyz(th,ph+phofs,&xi,&yi,&zi);
  if(cm)
  {
    x = xi;
    y = yi;
    z = zi;
  }
  else
  {
    IjkToXyz(thctr,phctr,xi,yi,zi,&x,&y,&z);
  }
  mark = new TMarker();
  mark->SetMarkerStyle(20);
  mark->SetMarkerSize(3.0);
  mark->DrawMarker(x,y);
  // Draw color box
  box = new TBox();
  text = new TLatex();
  text->SetTextAlign(12);
  ymin = ymin+0.05*ywid;
  ymax = ymax-0.05*ywid;
  ywid = ymax-ymin;
  if(gm)
  {
    for(i=0; i<5000; i++)
    {
      c = 0.3+0.00014*i;
      box->SetFillColor(TColor::GetColor(c,c,c));
      box->DrawBox(xmax+0.05*xwid,ymin+0.0002*ywid*i,
                   xmax+0.10*xwid,ymin+0.0002*ywid*(i+1));
    }
  }
  else
  {
    for(i=0; i<5000; i++)
    {
      box->SetFillColor(RGB(0.0,5000.0,(double)i));
      box->DrawBox(xmax+0.05*xwid,ymin+0.0002*ywid*i,
                   xmax+0.10*xwid,ymin+0.0002*ywid*(i+1));
    }
  }
  for(z=zstp*((int)(zmin/zstp)); z<=zmax; z+=zstp)
  {
    if(z < zmin) continue;
    if(zmax-zmin < 1.0)
    {
      sprintf(temp,"%.2f",z);
    } else
    if(zmax-zmin < 10.0)
    {
      sprintf(temp,"%.1f",z);
    }
    else
    {
      sprintf(temp,"%.0f",z);
    }
    text->DrawLatex(xmax+0.11*xwid,ymin+ywid*(z-zmin)/(zmax-zmin),temp);
  }
  if(strcmp(ztit,"") != 0)
  {
    text->DrawLatex(xmax-0.03*xwid,ymax,ztit);
  }
  if(strcmp(fout,"") != 0)
  {
    TPostScript myps(fout,111);
    canv->Draw();
    myps.Close();
  }
  if(!bm)
  {
    fprintf(stderr,"Type Ctrl-c to exit\n");
    appl->Run();
    while(1);
  }

  return 0;
}

int Init(void)
{
  Int_t i,j,n;
  Int_t err;
  Double_t w,x,y,z;
  Double_t th,ph;
  Double_t z1,z2;
  char line[MAXLINE];
  char str[NCOL2][MAXLINE];
  char *p,*endp;

  z1 = +1.0e100;
  z2 = -1.0e100;
  thctr = 0.0;
  phctr = 0.0;
  err = 0;
  i = 0;
  // format: rng i th_sun ph_sun wlen dth dph th_los ph_los z
  //     or: rng i thc phc wlen dth dph z
  while(fgets(line,MAXLINE,stdin) != NULL)
  {
    for(j=0,p=line; j<MAXENTR; j++,p+=n)
    {
      if(sscanf(p,"%s%n",str[j],&n) == EOF) break;
    }
    if(j!=NCOL2 && j!=NCOL)
    {
      fprintf(stderr,"Error, expected %d (or %d) columns >>> %s\n",NCOL2,NCOL,line);
      err = 1;
      break;
    }
    errno = 0;
    th = strtod(str[2],&endp);
    if(errno==ERANGE || *endp!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      err = 1;
      break;
    }
    errno = 0;
    ph = strtod(str[3],&endp);
    if(errno==ERANGE || *endp!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      err = 1;
      break;
    }
    errno = 0;
    w = strtod(str[4],&endp);
    if(errno==ERANGE || *endp!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      err = 1;
      break;
    }
    errno = 0;
    x = strtod(str[5],&endp);
    if(errno==ERANGE || *endp!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      err = 1;
      break;
    }
    errno = 0;
    y = strtod(str[6],&endp);
    if(errno==ERANGE || *endp!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      err = 1;
      break;
    }
    errno = 0;
    z = strtod(str[j-1],&endp);
    if(errno==ERANGE || *endp!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      err = 1;
      break;
    }
    if(i >= MAXDATA)
    {
      fprintf(stderr,"Warning, #data exceed the limit >>> %d",i);
      break;
    }
    thctr += th;
    phctr += ph;
    xdat[i] = x*D_TO_R;
    ydat[i] = y*D_TO_R;
    zdat[i] = z;
    if(z < z1) z1 = z;
    if(z > z2) z2 = z;
    i++;
  }
  ndat = i;
  fprintf(stderr,"%d data have been read.\n",ndat);
  if(ndat != ndvth*ndvph+1)
  {
    fprintf(stderr,"Error, ndat=%d, ndvth=%d, ndvph=%d, expected=%d\n",ndat,ndvth,ndvph,ndvth*ndvph+1);
    err = 1;
  }
  if(err)
  {
    return -1;
  }
  if(isnan(zmin)) zmin = z1;
  if(isnan(zmax)) zmax = z2;
  thctr = thctr*D_TO_R/ndat;
  phctr = phctr*D_TO_R/ndat;
  fprintf(stderr,"thctr: %13.4f\n",thctr*R_TO_D);
  fprintf(stderr,"phctr: %13.4f\n",phctr*R_TO_D);
  fprintf(stderr,"zmin : %13.4e\n",zmin);
  fprintf(stderr,"zmax : %13.4e\n",zmax);

  if(!bm || strcmp(fout,"")!=0)
  {
    appl = new TApplication("App",NULL,NULL);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineColor(10);
    gStyle->SetLabelFont(62,"X,Y,Z");
    gStyle->SetLabelSize(0.04,"X,Y,Z");
    gStyle->SetTitleFont(62,"X,Y,Z");
    canv = new TCanvas("canv","",0,0,cwid*(1.0+xmgn),cwid+23);
    canv->ToggleEventStatus();
    canv->SetTopMargin(0.02);
    canv->SetLeftMargin(0.15);
    canv->SetRightMargin(0.15);
    canv->SetBottomMargin(0.10);
    canv->SetFrameFillColor(10);
    canv->SetFillColor(10);
    canv->SetGrid(1,1);
    if(bm) canv->SetBatch(kTRUE);
  }

  return 0;
}

double range(double a)
{
  double b;

  b = a/PI2;
  return a-PI2*floor(b);
}

double range2(double a)
{
  double b;

  b = (a+PI)/PI2;
  return a-PI2*floor(b);
}

double SepCosine(double x1,double y1,double z1,double x2,double y2,double z2)
{
  double r1,r2;
  double epsilon = 1.0e-14;

  r1 = sqrt(x1*x1+y1*y1+z1*z1);
  r2 = sqrt(x2*x2+y2*y2+z2*z2);
  if(r1<epsilon || r2<epsilon) return 0.0;

  return (x1*x2+y1*y2+z1*z2)/(r1*r2);
}

double SepCosine2(double th1,double ph1,double th2,double ph2)
{
  return (sin(th1)*sin(th2)*cos(ph1-ph2)+cos(th1)*cos(th2));
}

double SepAngle(double x1,double y1,double z1,double x2,double y2,double z2)
{
  double r1,r2;
  double epsilon = 1.0e-14;

  r1 = sqrt(x1*x1+y1*y1+z1*z1);
  r2 = sqrt(x2*x2+y2*y2+z2*z2);
  if(r1<epsilon || r2<epsilon) return 0.0;

  return acos((x1*x2+y1*y2+z1*z2)/(r1*r2));
}

double SepAngle2(double th1,double ph1,double th2,double ph2)
{
  return acos(sin(th1)*sin(th2)*cos(ph1-ph2)+cos(th1)*cos(th2));
}

int SepToAzi(double th,double sep,double *dph)
{
  double s1,s2;

  s1 = sin(th);
  s2 = sin(0.5*sep);
  if(s1 < s2)
  {
    fprintf(stderr,"SepToAzi: error, th=%13.4e, sep=%13.4e\n",th,sep);
    return -1;
  }
  *dph = range(2.0*asin(s2/s1));

  return 0;
}

int PolToXyz(double th,double ph,double *x,double *y,double *z)
{
  *x = sin(th)*cos(ph);
  *y = sin(th)*sin(ph);
  *z = cos(th);

  return 0;
}

int XyzToPol(double x,double y,double z,double *th,double *ph)
{
  double r;

  r   = sqrt(x*x+y*y);
  *th = atan2(r,z);
  *ph = range(atan2(y,x));

  return 0;
}

void IjkToXyz(double th,double ph,
              double xi,double yi,double zi,
              double *x,double *y,double *z)
{
  *x =  xi*cos(th)*cos(ph)-yi*sin(ph)+zi*sin(th)*cos(ph);
  *y =  xi*cos(th)*sin(ph)+yi*cos(ph)+zi*sin(th)*sin(ph);
  *z = -xi*sin(th)+zi*cos(th);
}

void XyzToIjk(double th,double ph,
              double xi,double yi,double zi,
              double *x,double *y,double *z)
{
  *x =  xi*cos(th)*cos(ph)+yi*cos(th)*sin(ph)-zi*sin(th);
  *y = -xi*sin(ph)+yi*cos(ph);
  *z =  xi*sin(th)*cos(ph)+yi*sin(th)*sin(ph)+zi*cos(th);
}

int RGB(double min,double max,double val)
{
  double x,r;
  int R,G,B;

  if(min >= max)
  {
    fprintf(stderr,"RGB: Error, min=%e max=%e\n",min,max);
    return -1;
  }

  r = 255.0/11.0;

  if(val < min)
  {
    x = 0.0;
  } else
  if(val > max)
  {
    x = 50.0;
  }
  else
  {
    x = 50.0*(val-min)/(max-min);
  }

  if(x < 6.0)
  {
    R = (int)rint(-r*(x-6.0));
    G = 0;
    B = 255;
  } else
  if(x < 17.0)
  {
    R = 0;
    G = (int)rint(r*(x-6.0));
    B = 255;
  } else
  if(x < 28.0)
  {
    R = 0;
    G = 255;
    B = (int)rint(-r*(x-28.0));
  } else
  if(x < 39.0)
  {
    R = (int)rint(r*(x-28.0));
    G = 255;
    B = 0;
  }
  else
  {
    R = 255;
    G = (int)rint(-r*(x-50.0));
    B = 0;
  }

  return TColor::GetColor(R,G,B);
}

int GetOpt(int argn,char **args)
{
  int c,rt;
  int ntmp;
  int option_index = 0;
  int this_option_optind;
  double xtmp;
  char *endp;
  struct option long_options[] =
  {
    {"fout",1,0,'F'},
    {"ndvth",1,0,'m'},
    {"ndvph",1,0,'n'},
    {"xmin",1,0,'x'},
    {"xmax",1,0,'X'},
    {"ymin",1,0,'y'},
    {"ymax",1,0,'Y'},
    {"zmin",1,0,'z'},
    {"zmax",1,0,'Z'},
    {"zstp",1,0,'s'},
    {"thmin",1,0,'l'},
    {"thmax",1,0,'u'},
    {"phmin",1,0,'L'},
    {"phmax",1,0,'U'},
    {"phofs",1,0,'O'},
    {"xofs",1,0,'p'},
    {"yofs",1,0,'q'},
    {"pwid",1,0,'w'},
    {"cwid",1,0,'W'},
    {"title",1,0,'T'},
    {"xtit",1,0,'A'},
    {"ytit",1,0,'B'},
    {"ztit",1,0,'C'},
    {"circle",0,0,'c'},
    {"gray",0,0,'g'},
    {"batch",0,0,'b'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":F:m:n:x:X:y:Y:z:Z:s:l:u:L:U:O:p:q:w:W:T:A:B:C:cgbdvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'F':
        strncpy(fout,optarg,MAXLINE);
        break;
      case 'm':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) ndvth = ntmp;
        else
        {
          fprintf(stderr,"#Bins in Theta -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'n':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) ndvph = ntmp;
        else
        {
          fprintf(stderr,"#Bins in Phi -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'x':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') xmin = xtmp;
        else
        {
          fprintf(stderr,"X min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'X':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') xmax = xtmp;
        else
        {
          fprintf(stderr,"X max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'y':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') ymin = xtmp;
        else
        {
          fprintf(stderr,"Y min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'Y':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') ymax = xtmp;
        else
        {
          fprintf(stderr,"Y max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'z':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') zmin = xtmp;
        else
        {
          fprintf(stderr,"Z min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'Z':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') zmax = xtmp;
        else
        {
          fprintf(stderr,"Z max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 's':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) zstp = xtmp;
        else
        {
          fprintf(stderr,"Z step -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'l':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') thmin = xtmp*D_TO_R;
        else
        {
          fprintf(stderr,"Theta min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'u':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') thmax = xtmp*D_TO_R;
        else
        {
          fprintf(stderr,"Theta max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'L':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') phmin = xtmp*D_TO_R;
        else
        {
          fprintf(stderr,"Phi min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'U':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') phmax = xtmp*D_TO_R;
        else
        {
          fprintf(stderr,"Phi max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'O':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') phofs = xtmp*D_TO_R;
        else
        {
          fprintf(stderr,"Phi offset -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'p':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') xofs = xtmp;
        else
        {
          fprintf(stderr,"Pad X offset -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'q':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') yofs = xtmp;
        else
        {
          fprintf(stderr,"Pad Y offset -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'w':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) pwid = xtmp;
        else
        {
          fprintf(stderr,"Pad width -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'W':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) cwid = ntmp;
        else
        {
          fprintf(stderr,"Canvas width -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'T':
        strncpy(title,optarg,MAXLINE);
        break;
      case 'A':
        strncpy(xtit,optarg,MAXLINE);
        break;
      case 'B':
        strncpy(ytit,optarg,MAXLINE);
        break;
      case 'C':
        strncpy(ztit,optarg,MAXLINE);
        break;
      case 'c':
        cm = 1;
        break;
      case 'g':
        gm = 1;
        break;
      case 'b':
        bm = 1;
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
  int n = 15;
  char e[MAXLINE];
  char a[MAXLINE];
  char d[MAXLINE];

  fprintf(stderr,"ms720_fov ... display simulation data in the FOV of MS-720.\n");
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"ms720_fov -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"------------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),          As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," F -fout    |%s|%s|%s| %s\n",As(e,"Output psfile",n),  As(a,"File name",n),  As(d,"none",n),fout);
  fprintf(stderr," m -ndvth   |%s|%s|%s| %d\n",As(e,"#Bins in Theta",n), As(a,"#",n),          Ad(d,NDVTH,n),ndvth);
  fprintf(stderr," n -ndvph   |%s|%s|%s| %d\n",As(e,"#Bins in Phi  ",n), As(a,"#",n),          Ad(d,NDVPH,n),ndvph);
  fprintf(stderr," x -xmin    |%s|%s|%s| %e\n",As(e,"X min",n),          As(a,"value",n),      Af(d,XMIN,n),xmin);
  fprintf(stderr," X -xmax    |%s|%s|%s| %e\n",As(e,"X max",n),          As(a,"value",n),      Af(d,XMAX,n),xmax);
  fprintf(stderr," y -ymin    |%s|%s|%s| %e\n",As(e,"Y min",n),          As(a,"value",n),      Af(d,YMIN,n),ymin);
  fprintf(stderr," Y -ymax    |%s|%s|%s| %e\n",As(e,"Y max",n),          As(a,"value",n),      Af(d,YMAX,n),ymax);
  fprintf(stderr," z -zmin    |%s|%s|%s| %e\n",As(e,"Z min",n),          As(a,"value",n),      Af(d,ZMIN,n),zmin);
  fprintf(stderr," Z -zmax    |%s|%s|%s| %e\n",As(e,"Z max",n),          As(a,"value",n),      Af(d,ZMAX,n),zmax);
  fprintf(stderr," s -zstp    |%s|%s|%s| %e\n",As(e,"Z step",n),         As(a,"value",n),      Af(d,ZSTP,n),zstp);
  fprintf(stderr," l -thmin   |%s|%s|%s| %e\n",As(e,"Theta min",n),      As(a,"degree",n),     Af(d,THMIN*R_TO_D,n),thmin*R_TO_D);
  fprintf(stderr," u -thmax   |%s|%s|%s| %e\n",As(e,"Theta max",n),      As(a,"degree",n),     Af(d,THMAX*R_TO_D,n),thmax*R_TO_D);
  fprintf(stderr," L -phmin   |%s|%s|%s| %e\n",As(e,"Phi   min",n),      As(a,"degree",n),     Af(d,PHMIN*R_TO_D,n),phmin*R_TO_D);
  fprintf(stderr," U -phmax   |%s|%s|%s| %e\n",As(e,"Phi   max",n),      As(a,"degree",n),     Af(d,PHMAX*R_TO_D,n),phmax*R_TO_D);
  fprintf(stderr," O -phofs   |%s|%s|%s| %e\n",As(e,"Phi offset",n),     As(a,"degree",n),     Af(d,PHOFS*R_TO_D,n),phofs*R_TO_D);
  fprintf(stderr," p -xofs    |%s|%s|%s| %f\n",As(e,"Pad X offset",n),   As(a,"value",n),      Af(d,XOFS,n),xofs);
  fprintf(stderr," q -yofs    |%s|%s|%s| %f\n",As(e,"Pad Y offset",n),   As(a,"value",n),      Af(d,YOFS,n),yofs);
  fprintf(stderr," w -pwid    |%s|%s|%s| %f\n",As(e,"Pad width",n),      As(a,"value",n),      Af(d,PWID,n),pwid);
  fprintf(stderr," W -cwid    |%s|%s|%s| %d\n",As(e,"Canvas width",n),   As(a,"value",n),      Ad(d,CWID,n),cwid);
  fprintf(stderr," T -title   |%s|%s|%s| %s\n",As(e,"Graph  title",n),   As(a,"text",n),       As(d,"title",n),title);
  fprintf(stderr," A -xtit    |%s|%s|%s| %s\n",As(e,"X axis title",n),   As(a,"text",n),       As(d,"X",n),xtit);
  fprintf(stderr," B -ytit    |%s|%s|%s| %s\n",As(e,"Y axis title",n),   As(a,"text",n),       As(d,"Y",n),ytit);
  fprintf(stderr," C -ztit    |%s|%s|%s| %s\n",As(e,"Z axis title",n),   As(a,"text",n),       As(d,"",n),ztit);
  fprintf(stderr," c -circle  |%s|%s|%s| %d\n",As(e,"Circle  mode",n),   As(a,"nothing",n),    Ad(d,0,n),cm);
  fprintf(stderr," g -gray    |%s|%s|%s| %d\n",As(e,"Gray    mode",n),   As(a,"nothing",n),    Ad(d,0,n),gm);
  fprintf(stderr," b -batch   |%s|%s|%s| %d\n",As(e,"Batch   mode",n),   As(a,"nothing",n),    Ad(d,0,n),bm);
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),   As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),   As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),   As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"------------------------------------------------------------------------------\n");
  fprintf(stderr,"Input data format: rng i th_sun ph_sun wlen dth dph th_los ph_los z\n");
  fprintf(stderr," alternate format: rng i thc phc wlen dth dph z\n");

  return 0;
}
