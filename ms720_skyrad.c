/*********************************************************/
/* MS720_SKYRAD                                          */
/* Author: N.Manago Oct,30,2008                          */
/* $Revision: 503 $                                      */
/* $Date: 2008-11-26 07:17:02 +0900 (Wed, 26 Nov 2008) $ */
/*********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <bits/nan.h>
#include "newkur.h"
#include "strutil.h"

// Mathematical constants
#define		PI			3.141592653589793	// PI
#define		PI2			6.283185307179586	// 2*PI
#define		PI_2			1.570796326794896	// PI/2
#define		R_TO_D			57.29577951308232	// 180/PI    rad  -> deg
#define		D_TO_R			1.745329251994329e-02	// PI/180    deg  -> rad
// Common constants
#define		NONAME			"NONE"			// No name
#define		NODATA			0			// No data
#define		INVALID			-1			// Invalid value
#define		MAXENTR			20			// Max #entries in a line
#define		MAXLINE			256			// Max #chars in a line
#define		EPSILON			1.0e-14			// A small number
#define		DELTA			1.0e-13			// A small number
// Constants for input data
#define		INP_MAXDATA		100			// Max #input data
#define		INP_NCOL		19			// #Columns
#define		INP_NINT		10			// #Ints
#define		INP_NDBL		9			// #Doubles
#define		INP_N_WLEN		256			// #Wavelengths
#define		INP_NDIG		4			// #Digits in run#
#define		INP_TDIF		600			// Time difference
#define		INP_DNAM		"../data"		// Data directory
#define		INP_TYPE_DSR		1			// DSR
#define		INP_TYPE_AUR		2			// AUR
#define		INP_TYPE_SSR		3			// SSR
#define		INP_TYPE_UNKNOWN	INVALID			// Unknwon
// Constants for Direct Solar radiation
#define		DSR_MAXDATA		30			// Max #DSRs
#define		DSR_MAXWLEN		256			// Max #wavelengths for DSR
#define		DSR_RFOV		2.5			// FOV radius (in degree) for DSR
#define		DSR_SFOV		5.98020017e-3		// FOV area (in sr) for DSR
// Constants for Scattered Solar Radiation
#define		SSR_MAXDATA		30			// Max #SSRs
#define		SSR_MAXWLEN		256			// Max #wavelengths for SSR
#define		SSR_TH_LOS_MIN		-90.0
#define		SSR_TH_LOS_MAX		90.0
#define		SSR_RFOV		10.0			// FOV radius (in degree) for SSR
#define		SSR_SFOV		9.54557030e-2		// FOV area (in sr) for SSR
// Constants for Simulation
#define		SIM_XMGN		16.0			// X margin
#define		SIM_XSGM		4.0			// X sigma
#define		SIM_WSGM		20.0			// X width
#define		SIM_MAXCONV		1000			// #Data
// Constants for skyrad.PACK
#define		SKY_TIMZ		900			// Timezone (hour*100)
#define		SKY_ITYP		30			// Instrument type
#define		SKY_DNUM		0			// Data number
#define		SKY_N_WLEN		7			// Number of wavelengths
#define		SKY_N_OZON		0			// Number of ozone bands
#define		SKY_N_WATV		0			// Number of water bands
#define		SKY_N_ANGL		1			// Number of angles
#define		SKY_N_SIZE		20			// Number of sizes
#define		SKY_N_RMOD		20			// Number of modes
#define		SKY_LPMX		5			// Max #iterations
#define		SKY_MAXWLEN		256			// Max #wavelengths
#define		SKY_MAXANGL		SSR_MAXDATA		// Max #angles
#define		SKY_MAXMODE		100			// Max #modes
#define		SKY_TLON		135.0			// Longitude[deg] for time standard (0.0 for GMT)
#define		SKY_SLON		140.104128		// Longitude[deg] of obs. site
#define		SKY_SLAT		35.624594		// Latitude[deg]  of obs. site
#define		SKY_GALT		30.0			// Altitude[m]    of obs. site
#define		SKY_PRES		1.0			// Pressure in atm
#define		SKY_OZON		0.3			// Ozone amount in cm,STP
#define		SKY_REFR		1.5			// Refractive index (real)
#define		SKY_REFI		5.0e-3			// Refractive index (imag)
#define		SKY_RIND_WUNI		1.0			// Wavelength unit in nm
#define		SKY_GRHO		0.1			// Ground albedo
#define		SKY_GRHO_WUNI		1.0			// Wavelength unit in nm
#define		SKY_WLEN_MIN		NAN			// Minimum wavelength
#define		SKY_WLEN_MAX		NAN			// Maximum wavelength
#define		SKY_ANGL_MIN		3.0			// Minimum scattering angle
#define		SKY_ANGL_MAX		180.0			// Maximum scattering angle
#define		SKY_RMIN		0.01e-4			// Minimum radius in cm
#define		SKY_RMAX		20.0e-4			// Maximum radius in cm
#define		SKY_EPSA		0.01			// Absolute criterion for convergence
#define		SKY_EPS1		0.001			// Relative criterion for convergence
#define		SKY_EPS2		0.2			// Criterion to give up convergence
#define		SKY_INUM_A		"S06079.05"		// Instrument S/N
#define		SKY_INUM_B		"S07032.06"		// Instrument S/N
#define		SKY_PNAM		""			// Project name
#define		SKY_SNAM		"CEReS"			// Site name
#define		SKY_CNAM		"JAPAN"			// Country name
#define		SKY_DTOP		"."			// Top directory
#define		SKY_MKER "/usr/local/skyrad.pack/V42/MIEKER"	// Mie kernel
#define		SKY_FMET		"METEO.DAT"		// Meteorological data
#define		SKY_FINS		"ins.para"		// Instrument parameters
#define		SKY_FOBS		"obs.para"		// Observation parameters
#define		SKY_FPAR		"sproc.par"		// Parameters for sproc4
#define		SKY_FDT4		"yymmddnn.DT4"		// Input data for sproc4
#define		SKY_FTAG		"yymmddnn.tag"		// Tag   data for sproc4
#define		SKY_FDAT		"yymmddnn.dat"		// Data names for sproc4
#define		SKY_FNAM		"fname"			// Data list  for sproc4
#define		SKY_WLEN_A_001		455.510
#define		SKY_WLEN_A_002		555.743
#define		SKY_WLEN_A_003		635.876
#define		SKY_WLEN_A_004		672.540
#define		SKY_WLEN_A_005		752.304
#define		SKY_WLEN_A_006		792.030
#define		SKY_WLEN_A_007		877.628
#define		SKY_WLEN_B_001		468.778
#define		SKY_WLEN_B_002		552.675
#define		SKY_WLEN_B_003		629.792
#define		SKY_WLEN_B_004		676.644
#define		SKY_WLEN_B_005		750.055
#define		SKY_WLEN_B_006		799.909
#define		SKY_WLEN_B_007		879.241
// Constants for control
#define		CNT_CONF		NONAME
#define		CNT_MAXCMNT		20
// Functions
#define		MAX(a,b)		((b)>(a)?(b):(a))
#define		INTERP(x,x1,x2,y1,y2)	(y1+(y2-y1)/(x2-x1)*(x-x1))

// Parameters for input data
int		inp_ndig		= INP_NDIG;		// #Digits in run#
int		inp_tdif		= INP_TDIF;		// Time difference
int		inp_n_data		= NODATA;		// #Input data
int		inp_ms720[INP_MAXDATA];				// MS720ID
int		inp_rno[INP_MAXDATA];				// Run#
int		inp_flag[INP_MAXDATA];				// Flag
int		inp_type[INP_MAXDATA];				// Type
time_t		inp_time[INP_MAXDATA];				// Time
double		inp_th_sun[INP_MAXDATA];			// Solar zenith
double		inp_ph_sun[INP_MAXDATA];			// Solar azimuth
double		inp_th_los[INP_MAXDATA];			// LOS zenith
double		inp_ph_los[INP_MAXDATA];			// LOS azimuth
double		inp_amas[INP_MAXDATA];				// Airmass
double		inp_temp[INP_MAXDATA];				// Temperature
double		inp_wlen[INP_N_WLEN];				// Wavelength
double		inp_data[INP_MAXDATA][INP_N_WLEN];		// Input data
char		inp_dnam[MAXLINE]	= INP_DNAM;		// Data directory
// Parameters for Direct Solar radiation
int		dsr_ms720		= INVALID;		// DSR MS720ID
int		dsr_n_data		= NODATA;		// #DSR data
int		dsr_time[DSR_MAXDATA];				// Time
int		dsr_n_wlen		= NODATA;		// #Wavelengths
double		dsr_data[DSR_MAXDATA][DSR_MAXWLEN];		// Input data
double		dsr_th_sun[DSR_MAXDATA];			// Solar zenith
double		dsr_ph_sun[DSR_MAXDATA];			// Solar azimuth
double		dsr_amas[DSR_MAXDATA];				// Airmass
double		dsr_rfov		= DSR_RFOV;		// FOV radius (in degree) for DSR
double		dsr_sfov		= DSR_SFOV;		// FOV area (in sr) for DSR
// Parameters for Scattered Solar Radiation
int		ssr_ms720		= INVALID;		// SSR MS720ID
int		ssr_n_data		= NODATA;		// #SSR data
int		ssr_time[SSR_MAXDATA];				// Time
int		ssr_n_wlen		= NODATA;		// #Wavelengths
double		ssr_data[SSR_MAXDATA][SSR_MAXWLEN];		// Input data
double		ssr_th_sun[SSR_MAXDATA];			// Solar zenith
double		ssr_ph_sun[SSR_MAXDATA];			// Solar azimuth
double		ssr_th_los[SSR_MAXDATA];			// LOS zenith
double		ssr_ph_los[SSR_MAXDATA];			// LOS azimuth
double		ssr_th_los_min		= SSR_TH_LOS_MIN;	// Min LOS zenith in degree
double		ssr_th_los_max		= SSR_TH_LOS_MAX;	// Max LOS zenith in degree
double		ssr_rfov		= SSR_RFOV;		// FOV radius (in degree) for SSR
double		ssr_sfov		= SSR_SFOV;		// FOV area (in sr) for SSR
// Parameters for Simulation
double		sim_xmgn		= SIM_XMGN;		// X margin
double		sim_xsgm		= SIM_XSGM;		// X sigma
double		sim_wsgm		= SIM_WSGM;		// X width in sigma
double		sim_wlen_cnv[SIM_MAXCONV];			// Wavelength
double		sim_data_cnv[SIM_MAXCONV];			// Data
// Parameters for skyrad.PACK
int		sky_timz		= SKY_TIMZ;		// TImezone (hour*100)
int		sky_ityp		= SKY_ITYP;		// Instrument type
int		sky_dnum		= SKY_DNUM;		// Data number
int		sky_ms720		= INVALID;		// MS720ID
int		sky_n_wlen		= NODATA;		// Number of wavelengths
int		sky_n_ozon		= SKY_N_OZON;		// Number of ozone bands
int		sky_n_watv		= SKY_N_WATV;		// Number of water bands
int		sky_n_angl		= SKY_N_ANGL;		// Number of angles
int		sky_n_size		= SKY_N_SIZE;		// Number of sizes
int		sky_n_rmod		= SKY_N_RMOD;		// Number of modes
int		sky_lpmx		= SKY_LPMX;		// Max #iterations
int		sky_ianl[SKY_MAXWLEN];				// Analysis flag
time_t		sky_utim;					// Unix time
double		sky_tlon		= SKY_TLON;		// Longitude[deg] for time standard (0.0 for GMT)
double		sky_slon		= SKY_SLON;		// Longitude[deg] of obs. site
double		sky_slat		= SKY_SLAT;		// Latitude[deg]  of obs. site
double		sky_galt		= SKY_GALT;		// Altitude[m]    of obs. site
double		sky_pres		= SKY_PRES;		// Pressure in atm
double		sky_ozon		= SKY_OZON;		// Ozone amount in cm,STP
double		sky_refr_init		= SKY_REFR;		// Refractive index (real)
double		sky_refi_init		= SKY_REFI;		// Refractive index (imag)
double		sky_rind_wuni		= SKY_RIND_WUNI;	// Wavelength unit in nm
double		sky_grho_init		= SKY_GRHO;		// Ground albedo
double		sky_grho_wuni		= SKY_GRHO_WUNI;	// Wavelength unit in nm
double		sky_wlen_min		= SKY_WLEN_MIN;		// Minimum wavelength
double		sky_wlen_max		= SKY_WLEN_MAX;		// Maximum wavelength
double		sky_angl_min		= SKY_ANGL_MIN;		// Minimum scattering angle
double		sky_angl_max		= SKY_ANGL_MAX;		// Maximum scattering angle
double		sky_rmin		= SKY_RMIN;		// Minimum radius in cm
double		sky_rmax		= SKY_RMAX;		// Maximum radius in cm
double		sky_epsa		= SKY_EPSA;		// Absolute criterion for convergence
double		sky_eps1		= SKY_EPS1;		// Relative criterion for convergence
double		sky_eps2		= SKY_EPS2;		// Criterion to give up convergence
double		sky_wlen[SKY_MAXWLEN];				// Wavelength in nm
double		sky_wlen_cm[SKY_MAXWLEN];			// Wavelength in cm
double		sky_sfov[SKY_MAXWLEN];				// FOV area in sr
double		sky_refr[SKY_MAXWLEN];				// Refractive index (real)
double		sky_refi[SKY_MAXWLEN];				// Refractive index (imag)
double		sky_grho[SKY_MAXWLEN];				// Ground albedo
double		sky_angl[SKY_MAXANGL];				// Angle in degree
double		sky_rmod[SKY_MAXMODE] =				// Mode in cm
{
  2.38E-06,3.36E-06,4.74E-06,6.70E-06,9.46E-06,1.34E-05,1.89E-05,2.67E-05,
  3.77E-05,5.32E-05,7.52E-05,1.06E-04,1.50E-04,2.12E-04,2.99E-04,4.23E-04,
  5.97E-04,8.43E-04,1.19E-03,1.68E-03,
};
double		sky_lsgm[SKY_MAXMODE] =				// Sigma
{
  0.4,     0.4,     0.4,     0.4,     0.4,     0.4,     0.4,     0.4,
  0.4,     0.4,     0.4,     0.4,     0.4,     0.4,     0.4,     0.4,
  0.4,     0.4,     0.4,     0.4,
};
double		sky_au_sun;					// Sun-Earth distance in au
double		sky_th_sun;					// Solar zenith  angle in degree
double		sky_ph_sun;					// Solar azimuth angle in degree
double		sky_amas;					// Airmass
double		sky_th_los[SKY_MAXANGL];			// LOS zenith  angle in degree
double		sky_ph_los[SKY_MAXANGL];			// LOS azimuth angle in degree
double		sky_ddsr[SKY_MAXWLEN];				// DSR data
double		sky_dssr[SKY_MAXANGL][SKY_MAXWLEN];		// SSR data
double		sky_dcal[SKY_MAXWLEN];				// Calibration constant
char		sky_inum[MAXLINE]	= NONAME;		// Instrument S/N
char		sky_pnam[MAXLINE]	= SKY_PNAM;		// Project name
char		sky_snam[MAXLINE]	= SKY_SNAM;		// Site name
char		sky_cnam[MAXLINE]	= SKY_CNAM;		// Country name
char		sky_dtop[MAXLINE]	= SKY_DTOP;		// Top directory
char		sky_mker[MAXLINE]	= SKY_MKER;		// Mie kernel
char		sky_fmet[MAXLINE]	= SKY_FMET;		// Meteorological data
char		sky_fins[MAXLINE]	= SKY_FINS;		// Instrument parameters
char		sky_fobs[MAXLINE]	= SKY_FOBS;		// Observation parameters
char		sky_fpar[MAXLINE]	= SKY_FPAR;		// Parameters for sproc4
char		sky_fdt4[MAXLINE]	= SKY_FDT4;		// Input data for sproc4
char		sky_ftag[MAXLINE]	= SKY_FTAG;		// Tag   data for sproc4
char		sky_fdat[MAXLINE]	= SKY_FDAT;		// Data names for sproc4
char		sky_fnam[MAXLINE]	= SKY_FNAM;		// Data list  for sproc4
char		sky_find[MAXLINE]	= NONAME;		// Refractive index
char		sky_frho[MAXLINE]	= NONAME;		// Ground albedo
// Parameters for control
int		cnt_vb			= 0;			// Verbose mode
int		cnt_db			= 0;			// Debug   mode
int		cnt_hp			= 0;			// Help    mode
int		cnt_n_cmnt		= NODATA;
char		cnt_conf[MAXLINE]	= CNT_CONF;		// Configuration file
char		cnt_cmnt[CNT_MAXCMNT][MAXLINE];			// Comments

int WriteMet(void);
int WriteIns(void);
int WriteObs(void);
int WritePar(void);
int WriteDT4(void);
int WriteTag(void);
int WriteNam(void);
int SetFiles(void);
int SetParams(void);
int CalcConst(void);
int ReadInput(void);
int ReadRind(void);
int ReadGrho(void);
int ReadConfig(void);
int ReadDouble(char *s,int c,double u,double *d,int size);
double SepAngle2(double th1,double ph1,double th2,double ph2);
int Interp1D(const double *x1,const double *z1,int nx1,
             const double *x2,      double *z2,int nx2,int f);
int Convolute(double xmin,double xmax,double xstp,double xsgm,double wsgm,
              int ninp,const double *xinp,const double *yinp,
              int nout,double *xout,double *yout);
double Gauss(double g,double f0,double f);
int Sampling(int ninp,const double *xinp,const double *yinp,
             int nout,const double *xout,double *yout,double yuni);
int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) return -1;
  if(cnt_hp) {Usage(); return 0;}
  if(Init() < 0) return -1;

  if(WriteMet() < 0) return -1;
  if(WriteIns() < 0) return -1;
  if(WriteObs() < 0) return -1;
  if(WritePar() < 0) return -1;
  if(WriteDT4() < 0) return -1;
  if(WriteTag() < 0) return -1;
  if(WriteNam() < 0) return -1;
  if(SetFiles() < 0) return -1;

  return 0;
}

int Init(void)
{
  if(ReadInput() < 0) return -1;
  if(ReadRind()  < 0) return -1;
  if(ReadGrho()  < 0) return -1;
  if(SetParams() < 0) return -1;
  return 0;
}

int WriteMet(void)
{
// METEO.DAT
// - meteorological data(pressure & ozone) -
// "worldwide default"  : comment [40chr]
// "worldwide default"  : observation site/country name [20chr]
// 0.0           : longitude[deg] for time standard (0.0 for GMT)
// 0.0  0.0  0.0 : longitude[deg] latitude[deg] altitude[m] of obs. site
// 1   : NDAY(number of meleor. data)/ YYYY MM DD HR(LST) PRS[atm] O3[cm,STP]
// 2000  1  1  0.0   1.0   0.3
  char name[MAXLINE];
  char line[MAXLINE];
  FILE *fp;
  struct tm *tp;

  snprintf(name,MAXLINE,"%s/%s",sky_dtop,sky_fmet);
  if((fp=fopen(name,"w")) == NULL)
  {
    fprintf(stderr,"WriteMet: cannot open %s\n",name);
    return -1;
  }
  fprintf(fp,"- meteorological data(pressure & ozone) -\n");
  snprintf(line,MAXLINE,"\"%s\"",sky_pnam);
  fprintf(fp,"%20s : comment [40chr]\n",line);
  snprintf(line,MAXLINE,"\"%s\"",sky_snam);
  fprintf(fp,"%20s : observation site/country name [20chr]\n",line);
  fprintf(fp,"%8.3f %7s %6s : longitude[deg] for time standard (0.0 for GMT)\n",sky_tlon,"","");
  fprintf(fp,"%8.3f %7.3f %6.1f : longitude[deg] latitude[deg] altitude[m] of obs. site\n",sky_slon,sky_slat,sky_galt);
  fprintf(fp,"%20d : NDAY(number of meleor. data)/ YYYY MM DD HR(LST) PRS[atm] O3[cm,STP]\n",1);
  tp = gmtime(&sky_utim);
  strftime(line,MAXLINE,"%Y %m %d",tp);
  fprintf(fp,"%s %5.2f %4.2f %4.2f\n",line,tp->tm_hour+tp->tm_min/60.0+tp->tm_sec/3600.0,sky_pres,sky_ozon);
  fclose(fp);

  return 0;
}

int WriteIns(void)
{
// ins.para
// - instrument parameters -
// "PS1202006"  : instrument S/N [20chr]
// 30  : instrument type (POM-01L:10,11/ POM-01MKII:20,21,22/ POM-01,POM-02:30)
//  7  : NW(number of wavelengths)/ WL[cm]/ SVA[sr]
// 0.315E-04 0.400E-04 0.500E-04 0.675E-04 0.870E-04 0.940E-04 1.020E-04
// 2.390E-04 2.339E-04 2.358E-04 2.349E-04 2.394E-04 2.385E-04 2.385E-04
// 1   : NDAY(number of calib. constant data)/ YYYY MM DD HR(GMT) F0(IW=1,NW)
// 2003  6  1  0.0   0.00E-4 1.428E-4 2.915E-4 3.554E-4 2.465E-4 0.000E-4 2.076E-4
  int i;
  char name[MAXLINE];
  char line[MAXLINE];
  FILE *fp;
  struct tm *tp;

  snprintf(name,MAXLINE,"%s/%s",sky_dtop,sky_fins);
  if((fp=fopen(name,"w")) == NULL)
  {
    fprintf(stderr,"WriteIns: cannot open %s\n",name);
    return -1;
  }
  fprintf(fp,"- instrument parameters -\n");
  snprintf(line,MAXLINE,"\"%s\"",sky_inum);
  fprintf(fp,"%20s : instrument S/N [20chr]\n",line);
  fprintf(fp,"%20d : instrument type (POM-01L:10,11/ POM-01MKII:20,21,22/ POM-01,POM-02:30)\n",sky_ityp);
  fprintf(fp,"%20d : NW(number of wavelengths)/ WL[cm]/ SVA[sr]\n",sky_n_wlen);
  for(i=0; i<sky_n_wlen; i++)
  {
    fprintf(fp,"%9.3e%s",sky_wlen_cm[i],(i==sky_n_wlen-1?"\n":i%8==7?"\n":" "));
  }
  for(i=0; i<sky_n_wlen; i++)
  {
    fprintf(fp,"%9.3e%s",sky_sfov[i],(i==sky_n_wlen-1?"\n":i%8==7?"\n":" "));
  }
  fprintf(fp,"%20d : NDAY(number of calib. constant data)/ YYYY MM DD HR(GMT) F0(IW=1,NW)\n",1);
  tp = gmtime(&sky_utim);
  strftime(line,MAXLINE,"%Y %m %d",tp);
  fprintf(fp,"%s %5.2f ",line,tp->tm_hour+tp->tm_min/60.0+tp->tm_sec/3600.0);
  for(i=0; i<sky_n_wlen; i++)
  {
    fprintf(fp,"%9.3e%s",sky_dcal[i],(i==sky_n_wlen-1?"\n":i%8==7?"\n":" "));
  }
  fclose(fp);

  return 0;
}

int WriteObs(void)
{
// obs.para
// - observation parameters -
// "Example1"           : project name (or comment) [40chr]
// "PREDE Co.,Ltd."     : observation site (or ship) name [20chr]
// "JAPAN"              : country (or nationality) name [20chr]
// 135.0                : longitude[deg] for time standard (0.0 for GMT)
// 139.323  35.752  0.0 : longitude[deg] latitude[deg] altitude[m] of obs. site
// 24                   : NA(number of angles)/ SA(scattering angles[deg])
//    0    2    3    4    5    7   10   15   20   25   30   40
//   50   60   70   80   90  100  110  120  130  140  150  160
  int i;
  char name[MAXLINE];
  char line[MAXLINE];
  FILE *fp;

  snprintf(name,MAXLINE,"%s/%s",sky_dtop,sky_fobs);
  if((fp=fopen(name,"w")) == NULL)
  {
    fprintf(stderr,"WriteObs: cannot open %s\n",name);
    return -1;
  }
  fprintf(fp,"- observation parameters -\n");
  snprintf(line,MAXLINE,"\"%s\"",sky_pnam);
  fprintf(fp,"%20s : project name (or comment) [40chr]\n",line);
  snprintf(line,MAXLINE,"\"%s\"",sky_snam);
  fprintf(fp,"%20s : observation site (or ship) name [20chr]\n",line);
  snprintf(line,MAXLINE,"\"%s\"",sky_cnam);
  fprintf(fp,"%20s : country (or nationality) name [20chr]\n",line);
  fprintf(fp,"%8.3f %7s %6s : longitude[deg] for time standard (0.0 for GMT)\n",sky_tlon,"","");
  fprintf(fp,"%8.3f %7.3f %6.1f : longitude[deg] latitude[deg] altitude[m] of obs. site\n",sky_slon,sky_slat,sky_galt);
  fprintf(fp,"%20d : NA(number of angles)/ SA(scattering angles[deg])\n",sky_n_angl+1);
  fprintf(fp,"%4.0f ",0.0);
  for(i=0; i<sky_n_angl; i++)
  {
    fprintf(fp,"%4.0f%s",sky_angl[i],(i==sky_n_angl-1?"\n":(i+1)%12==11?"\n":" "));
  }
  fclose(fp);

  return 0;
}

int WritePar(void)
{
// sproc.par
// 0 1 1 1 1 0 : IOUT IPAR IVOL IAUR IPHS IPF0(output option) - 1:create / 0:not
// "ins.para.POM-01"      : instrument parameter file name
// "METEO.DAT.default"    : meteorological data(pressure & ozone) file name
// C --- atmospheric model.
// 1       : NLN
// 1.0     : (CONCA(L),L=1,NLN)
// 1.0     : (CONCM(L),L=1,NLN)
// C --- data information.
// 7  1  6 : NW(number of wavelengths)  NO3(No.of O3 abs.)  NWV(No.of WV abs.)
// 0.315E-4 0.400E-4 0.500E-4 0.675E-4 0.870E-4 0.940E-4 1.020E-4 : WL[cm]
// 0        1        2        1        1        0        1        : IANL
// 1.5      1.5      1.5      1.5      1.5      1.5      1.5      : CR
// -0.005   -0.005   -0.005   -0.005   -0.005   -0.005   -0.005   : CI
// 0.1      0.1      0.1      0.1      0.1      0.1      0.1      : GA
// C --- processing options.
// 1  1  1  20  0  0     : JTAU IDCR IDCI INVM IPLC IPHS
// 3.0  180.0            : ANGMN[deg]  ANGMX[deg]
// 5   0.01  0.001 0.2   : LOOPMX EPSA EPSA1 EPSA2
// 20  0.01E-4  20.0E-4  : NSIZE  RMIN[cm]  RMAX[cm]
// 20                    : NMODE/ RMODE[cm] / SSL
// 2.38E-06 3.36E-06 4.74E-06 6.70E-06 9.46E-06 1.34E-05 1.89E-05 2.67E-05
// 3.77E-05 5.32E-05 7.52E-05 1.06E-04 1.50E-04 2.12E-04 2.99E-04 4.23E-04
// 5.97E-04 8.43E-04 1.19E-03 1.68E-03
// 0.4      0.4      0.4      0.4      0.4      0.4      0.4      0.4
// 0.4      0.4      0.4      0.4      0.4      0.4      0.4      0.4
// 0.4      0.4      0.4      0.4
  int i;
  char name[MAXLINE];
  char line[MAXLINE];
  FILE *fp;

  snprintf(name,MAXLINE,"%s/%s",sky_dtop,sky_fpar);
  if((fp=fopen(name,"w")) == NULL)
  {
    fprintf(stderr,"WritePar: cannot open %s\n",name);
    return -1;
  }
  fprintf(fp,"0 1 1 1 1 0 : IOUT IPAR IVOL IAUR IPHS IPF0(output option) - 1:create / 0:not\n");
  snprintf(line,MAXLINE,"\"%s\"",sky_fins);
  fprintf(fp,"%20s : instrument parameter file name\n",sky_fins);
  snprintf(line,MAXLINE,"\"%s\"",sky_fmet);
  fprintf(fp,"%20s : meteorological data(pressure & ozone) file name\n",sky_fmet);
  fprintf(fp,"C --- atmospheric model.\n");
  fprintf(fp,"1       : NLN\n");
  fprintf(fp,"1.0     : (CONCA(L),L=1,NLN)\n");
  fprintf(fp,"1.0     : (CONCM(L),L=1,NLN)\n");
  fprintf(fp,"C --- data information.\n");
  fprintf(fp,"%2d %2d %2d : NW(number of wavelengths)  NO3(No.of O3 abs.)  NWV(No.of WV abs.)\n",sky_n_wlen,sky_n_ozon,sky_n_watv);
  for(i=0; i<sky_n_wlen; i++)
  {
    fprintf(fp,"%9.3e%s",sky_wlen_cm[i],(i==sky_n_wlen-1?" : WL[cm]\n":i%8==7?"\n":" "));
  }
  for(i=0; i<sky_n_wlen; i++)
  {
    fprintf(fp,"%9d%s",sky_ianl[i],(i==sky_n_wlen-1?" : IANL\n":i%8==7?"\n":" "));
  }
  for(i=0; i<sky_n_wlen; i++)
  {
    fprintf(fp,"%9.3f%s",sky_refr[i],(i==sky_n_wlen-1?" : CR\n":i%8==7?"\n":" "));
  }
  for(i=0; i<sky_n_wlen; i++)
  {
    fprintf(fp,"%9.3f%s",sky_refi[i],(i==sky_n_wlen-1?" : CI\n":i%8==7?"\n":" "));
  }
  for(i=0; i<sky_n_wlen; i++)
  {
    fprintf(fp,"%9.3f%s",sky_grho[i],(i==sky_n_wlen-1?" : GA\n":i%8==7?"\n":" "));
  }
  fprintf(fp,"C --- processing options.\n");
  fprintf(fp,"1  1  1  20  0  0     : JTAU IDCR IDCI INVM IPLC IPHS\n");
  fprintf(fp,"%5.1f %5.1f : ANGMN[deg]  ANGMX[deg]\n",sky_angl_min,sky_angl_max);
  fprintf(fp,"%2d %4.2f %5.3f %3.1f : LOOPMX EPSA EPSA1 EPSA2\n",sky_lpmx,sky_epsa,sky_eps1,sky_eps2);
  fprintf(fp,"%3d %8.2e %8.2e : NSIZE  RMIN[cm]  RMAX[cm]\n",sky_n_size,sky_rmin,sky_rmax);
  fprintf(fp,"%20d : NMODE/ RMODE[cm] / SSL\n",sky_n_rmod);
  for(i=0; i<sky_n_rmod; i++)
  {
    fprintf(fp,"%8.2e%s",sky_rmod[i],(i==sky_n_rmod-1?"\n":i%8==7?"\n":" "));
  }
  for(i=0; i<sky_n_rmod; i++)
  {
    fprintf(fp,"%8.2f%s",sky_lsgm[i],(i==sky_n_rmod-1?"\n":i%8==7?"\n":" "));
  }
  fclose(fp);

  return 0;
}

int WriteDT4(void)
{
// 03053000.DT4
//     2003  5 30  5 20  2   5.3339 : IY IM ID IHH IMM ISS TM(hr)
//   135.00  139.32   35.75    0.00 : ALNGS ALNG ALAT ALT(m)
//       24    81.4  -110.5  1.0136 : NA TH0 FI0 DST
//    7 : NW/ WL
//                  3.150E-05                 4.000E-05                 5.000E-05
//                  6.750E-05                 8.700E-05                 9.400E-05
//                  1.020E-04
//    TH     FI    F (IF SCA=0)  or R=U/F/M/SOLID
//    81.4    0.0  1.3277E-09   81.4    0.0  1.5364E-06   81.4    0.0  1.9359E-05
//    81.4    0.0  7.8032E-05   81.4    0.0  9.4154E-05   81.4    0.0  8.3118E-06
//    81.4    0.0  9.4186E-05
//    81.4    2.0  1.1094E-01   81.4    2.0  1.9162E+00   81.4    2.0  1.2889E+00
//    81.4    2.0  8.7950E-01   81.4    2.0  8.3168E-01   81.4    2.0  1.1029E+00
//    81.4    2.0  9.6483E-01
//    ...
  int i,j;
  char dnam[MAXLINE];
  char name[MAXLINE];
  char line[MAXLINE];
  FILE *fp;
  struct tm *tp;
  struct stat buf;

  tp = gmtime(&sky_utim);
  snprintf(sky_fdt4,MAXLINE,"%02d%02d%02d%02d.DT4",tp->tm_year%100,tp->tm_mon+1,tp->tm_mday,sky_dnum);
  snprintf(sky_ftag,MAXLINE,"%02d%02d%02d%02d.tag",tp->tm_year%100,tp->tm_mon+1,tp->tm_mday,sky_dnum);
  snprintf(sky_fdat,MAXLINE,"%02d%02d%02d%02d.dat",tp->tm_year%100,tp->tm_mon+1,tp->tm_mday,sky_dnum);
  snprintf(dnam,MAXLINE,"%s/DT4",sky_dtop);
  if(stat(dnam,&buf)<0 || !S_ISDIR(buf.st_mode))
  {
    if(mkdir(dnam,S_IRUSR | S_IRGRP | S_IXUSR | S_IXGRP | S_IWUSR ) < 0)
    {
      fprintf(stderr,"WriteDT4: error in making directory >>> %s\n",dnam);
      return -1;
    }
  }
  snprintf(name,MAXLINE,"%s/%s",dnam,sky_fdt4);
  if((fp=fopen(name,"w")) == NULL)
  {
    fprintf(stderr,"WriteDT4: cannot open %s\n",name);
    return -1;
  }
  strftime(line,MAXLINE,"%Y %m %d %H %M %S",tp);
  fprintf(fp,"%s %7.4f : IY IM ID IHH IMM ISS TM(hr)\n",line,tp->tm_hour+tp->tm_min/60.0+tp->tm_sec/3600.0);
  fprintf(fp,"%8.3f %8.3f %7.3f %6.1f : ALNGS ALNG ALAT ALT(m)\n",sky_tlon,sky_slon,sky_slat,sky_galt);
  fprintf(fp,"%4d %6.2f %6.2f %7.4f : NA TH0 FI0 DST\n",sky_n_angl+1,sky_th_sun,sky_ph_sun,sky_au_sun);
  fprintf(fp,"%4d : NW/ WL\n",sky_n_wlen);
  for(i=0; i<sky_n_wlen; i++)
  {
    fprintf(fp,"%25.3e%s",sky_wlen_cm[i],(i==sky_n_wlen-1?"\n":i%3==2?"\n":" "));
  }
  fprintf(fp,"   TH     FI    F (IF SCA=0)  or R=U/F/M/SOLID\n");
  for(i=0; i<sky_n_wlen; i++)
  {
    fprintf(fp,"%6.2f %6.2f %11.4e%s",sky_th_sun,0.0,sky_ddsr[i],(i==sky_n_wlen-1?"\n":i%3==2?"\n":" "));
  }
  for(j=0; j<sky_n_angl; j++)
  {
    for(i=0; i<sky_n_wlen; i++)
    {
      fprintf(fp,"%6.2f %6.2f %11.4e%s",sky_th_los[j],sky_ph_los[j]-sky_ph_sun,sky_dssr[j][i],(i==sky_n_wlen-1?"\n":i%3==2?"\n":" "));
    }
  }
  fclose(fp);

  return 0;
}

int WriteTag(void)
{
// No yyyy mm dd Hour   Long    Lat    Hs    SA(max)
//  1 2003  5 30  5.33  139.32  35.75   8.61 160.0
//  2 2003  5 30  5.67  139.32  35.75  12.44 150.0
// ...
  char dnam[MAXLINE];
  char name[MAXLINE];
  char line[MAXLINE];
  FILE *fp;
  struct tm *tp;
  struct stat buf;

  snprintf(dnam,MAXLINE,"%s/Tag",sky_dtop);
  if(stat(dnam,&buf)<0 || !S_ISDIR(buf.st_mode))
  {
    if(mkdir(dnam,S_IRUSR | S_IRGRP | S_IXUSR | S_IXGRP | S_IWUSR ) < 0)
    {
      fprintf(stderr,"WriteTag: error in making directory >>> %s\n",dnam);
      return -1;
    }
  }
  snprintf(name,MAXLINE,"%s/%s",dnam,sky_ftag);
  if((fp=fopen(name,"w")) == NULL)
  {
    fprintf(stderr,"WriteTag: cannot open %s\n",name);
    return -1;
  }
  fprintf(fp,"No yyyy mm dd Hour   Long    Lat    Hs    SA(max)\n");
  tp = gmtime(&sky_utim);
  strftime(line,MAXLINE,"%Y %m %d",tp);
  fprintf(fp,"%2d %s %5.2f %8.3f %7.3f %6.2f %6.2f\n",
              1,line,tp->tm_hour+tp->tm_min/60.0+tp->tm_sec/3600.0,
              sky_slon,sky_slat,90.0-sky_th_sun,sky_angl_max);
  fclose(fp);

  return 0;
}

int WriteNam(void)
{
// 03052900.dat
// 03053000.dat
// ...
  char name[MAXLINE];
  FILE *fp;

  snprintf(name,MAXLINE,"%s/%s",sky_dtop,sky_fnam);
  if((fp=fopen(name,"w")) == NULL)
  {
    fprintf(stderr,"WriteNam: cannot open %s\n",name);
    return -1;
  }
  fprintf(fp,"%s\n",sky_fdat);
  fclose(fp);

  return 0;
}

int SetFiles(void)
{
  char name[MAXLINE];
  struct stat buf;

  snprintf(name,MAXLINE,"%s/Aur",sky_dtop);
  if(stat(name,&buf)<0 || !S_ISDIR(buf.st_mode))
  {
    if(mkdir(name,S_IRUSR | S_IRGRP | S_IXUSR | S_IXGRP | S_IWUSR ) < 0)
    {
      fprintf(stderr,"SetFiles: error in making directory >>> %s\n",name);
      return -1;
    }
  }
  snprintf(name,MAXLINE,"%s/Par",sky_dtop);
  if(stat(name,&buf)<0 || !S_ISDIR(buf.st_mode))
  {
    if(mkdir(name,S_IRUSR | S_IRGRP | S_IXUSR | S_IXGRP | S_IWUSR ) < 0)
    {
      fprintf(stderr,"SetFiles: error in making directory >>> %s\n",name);
      return -1;
    }
  }
  snprintf(name,MAXLINE,"%s/Phs",sky_dtop);
  if(stat(name,&buf)<0 || !S_ISDIR(buf.st_mode))
  {
    if(mkdir(name,S_IRUSR | S_IRGRP | S_IXUSR | S_IXGRP | S_IWUSR ) < 0)
    {
      fprintf(stderr,"SetFiles: error in making directory >>> %s\n",name);
      return -1;
    }
  }
  snprintf(name,MAXLINE,"%s/Vol",sky_dtop);
  if(stat(name,&buf)<0 || !S_ISDIR(buf.st_mode))
  {
    if(mkdir(name,S_IRUSR | S_IRGRP | S_IXUSR | S_IXGRP | S_IWUSR ) < 0)
    {
      fprintf(stderr,"SetFiles: error in making directory >>> %s\n",name);
      return -1;
    }
  }
  snprintf(name,MAXLINE,"%s/MIEKER",sky_dtop);
  if(stat(name,&buf)<0 || (!S_ISREG(buf.st_mode)&&!S_ISLNK(buf.st_mode)))
  {
    if(symlink(sky_mker,name) < 0)
    {
      fprintf(stderr,"SetFiles: error in linking file >>> %s\n",name);
      return -1;
    }
  }

  return 0;
}

int SetParams(void)
{
  int i,j,n;

  sky_utim   = 0;
  sky_th_sun = 0.0;
  sky_ph_sun = 0.0;
  sky_amas   = 0.0;
  for(i=0; i<sky_n_wlen; i++)
  {
    sky_ddsr[i] = 0.0;
  }
  for(j=0,n=0; j<dsr_n_data; j++)
  {
    sky_utim   += dsr_time[j];
    sky_th_sun += dsr_th_sun[j];
    sky_ph_sun += dsr_ph_sun[j];
    sky_amas   += dsr_amas[j];
    for(i=0; i<sky_n_wlen; i++)
    {
      sky_ddsr[i] += dsr_data[j][i];
    }
    n++;
  }
  if(n > 1)
  {
    sky_utim   /= n;
    sky_th_sun /= (double)n;
    sky_ph_sun /= (double)n;
    sky_amas   /= (double)n;
    for(i=0; i<sky_n_wlen; i++)
    {
      sky_ddsr[i] /= (double)n;
    }
  }
  sky_utim += inp_time[0];

  if(CalcConst() < 0) return -1;

  sky_n_angl = ssr_n_data;
  for(j=0; j<sky_n_angl; j++)
  {
    sky_th_los[j] = ssr_th_los[j];
    sky_ph_los[j] = ssr_ph_los[j];
    sky_angl[j] = SepAngle2(sky_th_sun*D_TO_R,sky_ph_sun*D_TO_R,
                            sky_th_los[j]*D_TO_R,sky_ph_los[j]*D_TO_R)*R_TO_D;
  }
  for(i=0; i<sky_n_wlen; i++)
  {
    sky_sfov[i] = ssr_sfov;
  }
  for(j=0; j<sky_n_angl; j++)
  {
    for(i=0; i<sky_n_wlen; i++)
    {
      sky_dssr[j][i] = ssr_data[j][i]/(sky_ddsr[i]*sky_amas*sky_sfov[i]);
    }
  }
  for(i=0; i<sky_n_wlen; i++)
  {
    sky_ddsr[i] *= (sky_au_sun*sky_au_sun);
  }

  for(i=0; i<sky_n_wlen; i++)
  {
    sky_ianl[i] = 1;
  }

  return 0;
}

int CalcConst(void)
{
  int n;
  double al,az;
  double th,ph,au;
  double epsilon = 1.0;
  char cmnd[MAXLINE];
  char line[MAXLINE];
  char str1[MAXLINE];
  char str2[MAXLINE];
  char str3[MAXLINE];
  char str4[MAXLINE];
  char str5[MAXLINE];
  char str6[MAXLINE];
  char str7[MAXLINE];
  char *p;
  FILE *fp;

  snprintf(cmnd,MAXLINE,"planet_pos -b 10 -U %ld",sky_utim-sky_timz*36);
  if((fp=popen(cmnd,"r")) == NULL)
  {
    fprintf(stderr,"CalcConst: failed in command >>> %s\n",cmnd);
    return -1;
  }
  while(1)
  {
    if(fgets(line,MAXLINE,fp) == NULL) break;
  }
  pclose(fp);
  if(sscanf(line,"%s%s%s%s%s%s%s",str1,str2,str3,str4,str5,str6,str7) != 7)
  {
    fprintf(stderr,"CalcConst: read error >>> %s\n",line);
    return -1;
  }
  errno = 0;
  al = strtod(str5,&p);
  if(errno==ERANGE || *p!='\0')
  {
    fprintf(stderr,"CalcConst: convert error >>> %s\n",line);
    return -1;
  }
  th = 90.0-al;
  errno = 0;
  az = strtod(str6,&p);
  if(errno==ERANGE || *p!='\0')
  {
    fprintf(stderr,"CalcConst: convert error >>> %s\n",line);
    return -1;
  }
  ph = az-180.0;
  errno = 0;
  au = strtod(str7,&p);
  if(errno==ERANGE || *p!='\0')
  {
    fprintf(stderr,"CalcConst: convert error >>> %s\n",line);
    return -1;
  }
  if(SepAngle2(th*D_TO_R,ph*D_TO_R,sky_th_sun*D_TO_R,sky_ph_sun*D_TO_R)*R_TO_D > epsilon)
  {
    fprintf(stderr,"CalcConst: error, th=%13.4f ph=%13.4f sky_th_sun=%13.4f sky_ph_sun=%13.4f\n",
                    th,ph,sky_th_sun,sky_ph_sun);
  }
  sky_au_sun = au;

  if((n=Convolute(sky_wlen_min-sim_xmgn,sky_wlen_max+sim_xmgn,1.0,sim_xsgm,sim_wsgm,
                  kur_n_wlen,kur_wlen,kur_data,SIM_MAXCONV,sim_wlen_cnv,sim_data_cnv)) < 0)
  {
    fprintf(stderr,"CalcConst: error in Convolute().\n");
    return -1;
  }
  if(Sampling(n,sim_wlen_cnv,sim_data_cnv,sky_n_wlen,sky_wlen,sky_dcal,1.0) < 0)
  {
    fprintf(stderr,"CalcConst: error in Sampling().\n");
    return -1;
  }

  return 0;
}

int ReadInput(void)
{
  int i,j,k,n;
  int nc;
  int err;
  int index[INP_N_WLEN];
  int v[INP_NINT];
  int dsr_repeat[DSR_MAXDATA];
  int ssr_repeat[SSR_MAXDATA];
  time_t dsr_time_tmp;
  time_t ssr_time_tmp;
  double w[INP_NDBL];
  double ssr_th_los_tmp;
  double ssr_ph_los_tmp;
  double dsr_inp_avg[DSR_MAXDATA][INP_N_WLEN];
  double ssr_inp_avg[SSR_MAXDATA][INP_N_WLEN];
  char line[MAXLINE];
  char item[MAXENTR][MAXLINE];
  char fnam[MAXLINE];
  char format[MAXLINE];
  char *p;
  FILE *fp;
  struct tm t;

  // read stdin
  // format:
  // year month day hour min sec site_id ms720_id flag rno salt sazi airmass temperature valt vazi fov exposure par
  putenv("TZ=");
  t.tm_isdst = 0;
  i = 0;
  while(fgets(line,MAXLINE,stdin) != NULL)
  {
    // read line
    for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
    {
      if(sscanf(p,"%s%n",item[n],&nc) == EOF) break;
    }
    if(n != INP_NCOL)
    {
      fprintf(stderr,"Error, expected %d columns >>> %d\n",INP_NCOL,n);
      return -1;
    }
    // convert values
    for(n=0; n<INP_NINT; n++)
    {
      errno = 0;
      v[n] = strtol(item[n],&p,10);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"Convert error >>> %s\n",item[n]);
        return -1;
      }
    }
    for(n=INP_NINT; n<INP_NCOL; n++)
    {
      errno = 0;
      w[n-INP_NINT] = strtod(item[n],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"Convert error >>> %s\n",item[n]);
        return -1;
      }
    }
    if(i >= INP_MAXDATA)
    {
      fprintf(stderr,"Warning, #data exceed the limit >>> %d\n",i);
      break;
    }
    t.tm_year = v[0]-1900;
    t.tm_mon  = v[1]-1;
    t.tm_mday = v[2];
    t.tm_hour = v[3];
    t.tm_min  = v[4];
    t.tm_sec  = v[5];
    t.tm_isdst = 0;
    inp_time[i] = mktime(&t);
    inp_ms720[i]  = v[7];
    inp_flag[i]   = v[8];
    inp_rno[i]    = v[9];
    inp_th_sun[i] = 90.0-w[0];
    inp_ph_sun[i] = w[1]-180.0;
    inp_amas[i]   = w[2];
    inp_temp[i]   = w[3];
    inp_th_los[i] = 90.0-w[4];
    inp_ph_los[i] = w[5]-180.0;
    if(inp_flag[i]>=100 && inp_flag[i]<200) // DSR
    {
      inp_type[i] = INP_TYPE_DSR;
    } else
    if(inp_flag[i]>=300 && inp_flag[i]<400) // AUR
    {
      inp_type[i] = INP_TYPE_AUR;
    } else
    if(inp_flag[i]>=0 && inp_flag[i]<100)   // SSR
    {
      inp_type[i] = INP_TYPE_SSR;
    }
    else
    {
      inp_type[i] = INP_TYPE_UNKNOWN;
    }
    i++;
  }
  inp_n_data = i;
  fprintf(stderr,"%d data have been read\n",inp_n_data);
  if(inp_n_data < 1)
  {
    fprintf(stderr,"Error, #data is not enough.\n");
    return -1;
  }

  // read input data
  for(j=0; j<INP_N_WLEN; j++)
  {
    inp_wlen[j] = 0.0;
  }
  sprintf(format,"%%s/%%0%dd.dat",inp_ndig);
  for(i=0; i<inp_n_data; i++)
  {
    sprintf(fnam,format,inp_dnam,inp_rno[i]);
    if(cnt_vb > 0)
    {
      fprintf(stderr,"Reading %s\n",fnam);
    }
    if((fp=fopen(fnam,"r")) == NULL)
    {
      fprintf(stderr,"Error, cannot open %s\n",fnam);
      return -1;
    }
    err = 0;
    j = 0;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      if(sscanf(line,"%s%s",item[0],item[1]) != 2)
      {
        fprintf(stderr,"Read error >>> %s\n",line);
        err = 1;
        break;
      }
      errno = 0;
      for(n=0; n<2; n++)
      {
        w[n] = strtod(item[n],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"Convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(err) break;
      if(j >= INP_N_WLEN)
      {
        fprintf(stderr,"Error, #data exceed the limit >>> %d\n",j);
        err = 1;
        break;
      }
      if(inp_wlen[j] < EPSILON)
      {
        inp_wlen[j] = w[0];
      }
      else
      {
        if(fabs(inp_wlen[j]-w[0]) > DELTA)
        {
          fprintf(stderr,"Error, inp_wlen[%d]=%13.4e, w[0]=%13.4e\n",j,inp_wlen[j],w[0]);
          err = 1;
          break;
        }
      }
      inp_data[i][j] = w[1];
      j++;
    }
    fclose(fp);
    if(err)
    {
      return -1;
    }
    if(j != INP_N_WLEN)
    {
      fprintf(stderr,"Error, i=%d j=%d\n",i,j);
      return -1;
    }
  }

  // SKY wavelength
  sky_ms720 = inp_ms720[0];
  switch(sky_ms720)
  {
    case 1:
      strncpy(sky_inum,SKY_INUM_A,MAXLINE);
      break;
    case 2:
      strncpy(sky_inum,SKY_INUM_B,MAXLINE);
      break;
    default:
      fprintf(stderr,"Error, invalid MS720ID >>> %d\n",sky_ms720);
      return -1;
      break;
  }
  if(sky_n_wlen == NODATA)
  {
    switch(sky_ms720)
    {
      case 1:
        sky_wlen[ 0] = SKY_WLEN_A_001;
        sky_wlen[ 1] = SKY_WLEN_A_002;
        sky_wlen[ 2] = SKY_WLEN_A_003;
        sky_wlen[ 3] = SKY_WLEN_A_004;
        sky_wlen[ 4] = SKY_WLEN_A_005;
        sky_wlen[ 5] = SKY_WLEN_A_006;
        sky_wlen[ 6] = SKY_WLEN_A_007;
        sky_n_wlen = SKY_N_WLEN;
        break;
      case 2:
        sky_wlen[ 0] = SKY_WLEN_B_001;
        sky_wlen[ 1] = SKY_WLEN_B_002;
        sky_wlen[ 2] = SKY_WLEN_B_003;
        sky_wlen[ 3] = SKY_WLEN_B_004;
        sky_wlen[ 4] = SKY_WLEN_B_005;
        sky_wlen[ 5] = SKY_WLEN_B_006;
        sky_wlen[ 6] = SKY_WLEN_B_007;
        sky_n_wlen = SKY_N_WLEN;
        break;
    }
  }
  sky_wlen_min = +1.0e100;
  sky_wlen_max = -1.0e100;
  for(i=0; i<sky_n_wlen; i++)
  {
    sky_wlen_cm[i] = sky_wlen[i]*1.0e-7;
    if(sky_wlen[i] < sky_wlen_min) sky_wlen_min = sky_wlen[i];
    if(sky_wlen[i] > sky_wlen_max) sky_wlen_max = sky_wlen[i];
  }
  for(j=0; j<sky_n_wlen; j++)
  {
    index[j] = -1;
    for(k=0; k<INP_N_WLEN; k++)
    {
      if(fabs(sky_wlen[j]-inp_wlen[k]) < EPSILON)
      {
        index[j] = k;
        break;
      }
    }
    if(index[j] < 0)
    {
      fprintf(stderr,"Error, failed in finding SKY wavelength >>> %13.4f\n",sky_wlen[j]);
      return -1;
    }
  }

  // DSR wavelength
  for(i=0; i<inp_n_data; i++)
  {
    if(inp_type[i] == INP_TYPE_DSR) // DSR
    {
      dsr_ms720 = inp_ms720[i];
      break;
    }
  }
  if(dsr_ms720 != sky_ms720)
  {
    fprintf(stderr,"SKY/DSR MS720ID mismatch >>> %d/%d\n",sky_ms720,dsr_ms720);
    return -1;
  }
  dsr_n_wlen = sky_n_wlen;

  // SSR wavelength
  for(i=0; i<inp_n_data; i++)
  {
    if(inp_type[i] == INP_TYPE_SSR) // SSR
    {
      ssr_ms720 = inp_ms720[i];
      break;
    }
  }
  if(ssr_ms720 != sky_ms720)
  {
    fprintf(stderr,"SKY/SSR MS720ID mismatch >>> %d/%d\n",sky_ms720,ssr_ms720);
    return -1;
  }
  ssr_n_wlen = sky_n_wlen;

  // copy data
  dsr_n_data = -1;
  dsr_time_tmp = -1;
  for(i=0; i<DSR_MAXDATA; i++)
  {
    dsr_time[i]   = 0;
    dsr_th_sun[i] = 0.0;
    dsr_ph_sun[i] = 0.0;
    dsr_amas[i]   = 0.0;
    dsr_repeat[i] = 0;
    for(j=0; j<INP_N_WLEN; j++)
    {
      dsr_inp_avg[i][j] = 0.0;
    }
  }
  ssr_n_data = -1;
  ssr_time_tmp = -1;
  ssr_th_los_tmp = -1.0;
  ssr_ph_los_tmp = -1.0;
  for(i=0; i<SSR_MAXDATA; i++)
  {
    ssr_time[i]   = 0;
    ssr_th_sun[i] = 0.0;
    ssr_ph_sun[i] = 0.0;
    ssr_repeat[i] = 0;
    for(j=0; j<INP_N_WLEN; j++)
    {
      ssr_inp_avg[i][j] = 0.0;
    }
  }
  for(i=0; i<inp_n_data; i++)
  {
    if(inp_type[i] == INP_TYPE_DSR) // DSR
    {
      if(abs(dsr_time_tmp-inp_time[i]) > inp_tdif)
      {
        dsr_time_tmp = inp_time[i];
        dsr_n_data++;
        if(dsr_n_data >= DSR_MAXDATA)
        {
          fprintf(stderr,"Warning, #DSR data exceed the limit >>> %d\n",dsr_n_data);
          dsr_n_data--;
          break;
        }
      }
      dsr_time[dsr_n_data]   += inp_time[i]-inp_time[0];
      dsr_th_sun[dsr_n_data] += inp_th_sun[i];
      dsr_ph_sun[dsr_n_data] += inp_ph_sun[i];
      dsr_amas[dsr_n_data]   += inp_amas[i];
      for(j=0; j<INP_N_WLEN; j++)
      {
        dsr_inp_avg[dsr_n_data][j] += inp_data[i][j];
      }
      dsr_repeat[dsr_n_data] += 1;
    } else
    if(inp_type[i] == INP_TYPE_SSR) // SSR
    {
      if(inp_th_los[i]<ssr_th_los_min || inp_th_los[i]>ssr_th_los_max) continue;
      if(abs(ssr_time_tmp-inp_time[i])>inp_tdif ||
         fabs(ssr_th_los_tmp-inp_th_los[i])>DELTA ||
         fabs(ssr_ph_los_tmp-inp_ph_los[i])>DELTA)
      {
        ssr_time_tmp   = inp_time[i];
        ssr_th_los_tmp = inp_th_los[i];
        ssr_ph_los_tmp = inp_ph_los[i];
        ssr_n_data++;
        if(ssr_n_data >= SSR_MAXDATA)
        {
          fprintf(stderr,"Warning, #SSR data exceed the limit >>> %d\n",ssr_n_data);
          ssr_n_data--;
          break;
        }
        ssr_th_los[ssr_n_data] = inp_th_los[i];
        ssr_ph_los[ssr_n_data] = inp_ph_los[i];
      }
      ssr_time[ssr_n_data]   += inp_time[i]-inp_time[0];
      ssr_th_sun[ssr_n_data] += inp_th_sun[i];
      ssr_ph_sun[ssr_n_data] += inp_ph_sun[i];
      for(j=0; j<INP_N_WLEN; j++)
      {
        ssr_inp_avg[ssr_n_data][j] += inp_data[i][j];
      }
      ssr_repeat[ssr_n_data] += 1;
    }
  }
  dsr_n_data++;
  ssr_n_data++;
  // store data for DSR
  for(i=0; i<dsr_n_data; i++)
  {
    if(dsr_repeat[i] > 1)
    {
      dsr_time[i]   /= dsr_repeat[i];
      dsr_th_sun[i] /= dsr_repeat[i];
      dsr_ph_sun[i] /= dsr_repeat[i];
      dsr_amas[i]   /= dsr_repeat[i];
      for(j=0; j<INP_N_WLEN; j++)
      {
        dsr_inp_avg[i][j] /= dsr_repeat[i];
      }
    }
    for(j=0; j<dsr_n_wlen; j++)
    {
      dsr_data[i][j] = dsr_inp_avg[i][index[j]];
    }
  }
  // store data for SSR
  for(i=0; i<ssr_n_data; i++)
  {
    if(ssr_repeat[i] > 1)
    {
      ssr_time[i]   /= ssr_repeat[i];
      ssr_th_sun[i] /= ssr_repeat[i];
      ssr_ph_sun[i] /= ssr_repeat[i];
      for(j=0; j<INP_N_WLEN; j++)
      {
        ssr_inp_avg[i][j] /= ssr_repeat[i];
      }
    }
    for(j=0; j<ssr_n_wlen; j++)
    {
      ssr_data[i][j] = ssr_inp_avg[i][index[j]];
    }
  }

  return 0;
}

int ReadRind(void)
{
  int i,j,n;
  int err;
  double v[3];
  double wlen[SKY_MAXWLEN];
  double refr[SKY_MAXWLEN];
  double refi[SKY_MAXWLEN];
  double temp[SKY_MAXWLEN];
  char line[MAXLINE];
  char str[3][MAXLINE];
  char *p;
  FILE *fp;

  if(strcmp(sky_find,NONAME) == 0)
  {
    for(i=0; i<sky_n_wlen; i++)
    {
      sky_refr[i] = sky_refr_init;
      sky_refi[i] = sky_refi_init;
    }
    return 0;
  } else
  if((fp=fopen(sky_find,"r")) == NULL)
  {
    fprintf(stderr,"ReadRind: cannot open %s\n",sky_find);
    return -1;
  }
  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    if(sscanf(line,"%s%s%s",str[0],str[1],str[2]) != 3)
    {
      fprintf(stderr,"ReadRind: read error >>> %s\n",line);
      err = 1;
      break;
    }
    for(j=0; j<3; j++)
    {
      errno = 0;
      v[i] = strtod(str[j],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"ReadRind: convert error >>> %s\n",line);
        err = 1;
        break;
      }
    }
    if(err) break;
    if(v[0]*sky_rind_wuni < sky_wlen_min-DELTA) continue;
    if(v[0]*sky_rind_wuni > sky_wlen_max+DELTA) break;
    if(i >= SKY_MAXWLEN)
    {
      fprintf(stderr,"ReadRind: warning, #data exceed the limit >>> %d\n",i);
      break;
    }
    wlen[i] = v[0]*sky_rind_wuni;
    refr[i] = v[1];
    refi[i] = v[2];
    i++;
  }
  n = i;
  fclose(fp);
  if(err)
  {
    return -1;
  }

  for(i=0; i<sky_n_wlen; i++)
  {
    err = 1;
    for(j=0; j<n; j++)
    {
      if(fabs(sky_wlen[i]-wlen[j]) < DELTA)
      {
        sky_refr[i] = refr[i];
        sky_refi[i] = refi[i];
        err = 0;
        break;
      }
    }
    if(err) break;
  }
  if(err)
  {
    Interp1D(wlen,refr,n,sky_wlen,sky_refr,sky_n_wlen,1);
    for(i=0; i<n; i++)
    {
      temp[i] = log10(refi[i]);
    }
    Interp1D(wlen,temp,n,sky_wlen,refi,sky_n_wlen,0);
    for(i=0; i<sky_n_wlen; i++)
    {
      sky_refi[i] = pow(10.0,refi[i]);
    }
  }

  return 0;
}

int ReadGrho(void)
{
  int i,j,n;
  int err;
  double v[2];
  double wlen[SKY_MAXWLEN];
  double grho[SKY_MAXWLEN];
  char line[MAXLINE];
  char str[2][MAXLINE];
  char *p;
  FILE *fp;

  if(strcmp(sky_frho,NONAME) == 0)
  {
    for(i=0; i<sky_n_wlen; i++)
    {
      sky_grho[i] = sky_grho_init;
    }
    return 0;
  } else
  if((fp=fopen(sky_frho,"r")) == NULL)
  {
    fprintf(stderr,"ReadGrho: cannot open %s\n",sky_frho);
    return -1;
  }
  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    if(sscanf(line,"%s%s",str[0],str[1]) != 2)
    {
      fprintf(stderr,"ReadGrho: read error >>> %s\n",line);
      err = 1;
      break;
    }
    for(j=0; j<2; j++)
    {
      errno = 0;
      v[i] = strtod(str[j],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"ReadGrho: convert error >>> %s\n",line);
        err = 1;
        break;
      }
    }
    if(err) break;
    if(v[0]*sky_grho_wuni < sky_wlen_min-DELTA) continue;
    if(v[0]*sky_grho_wuni > sky_wlen_max+DELTA) break;
    if(i >= SKY_MAXWLEN)
    {
      fprintf(stderr,"ReadGrho: warning, #data exceed the limit >>> %d\n",i);
      break;
    }
    wlen[i] = v[0]*sky_grho_wuni;
    grho[i] = v[1];
    i++;
  }
  n = i;
  fclose(fp);
  if(err)
  {
    return -1;
  }

  for(i=0; i<sky_n_wlen; i++)
  {
    err = 1;
    for(j=0; j<n; j++)
    {
      if(fabs(sky_wlen[i]-wlen[j]) < DELTA)
      {
        sky_grho[i] = grho[i];
        err = 0;
        break;
      }
    }
    if(err) break;
  }
  if(err)
  {
    Interp1D(wlen,grho,n,sky_wlen,sky_grho,sky_n_wlen,1);
  }

  return 0;
}

int ReadConfig(void)
{
  int n,nc;
  int num;
  int err;
  double uni;
  char line[MAXLINE];
  char temp[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char *p;
  FILE *fp;

  if(strcmp(cnt_conf,"")==0 || strcmp(cnt_conf,NONAME)==0)
  {
    return 0;
  } else
  if((fp=fopen(cnt_conf,"r")) == NULL)
  {
    fprintf(stderr,"ReadConfig: error, cannot open %s\n",cnt_conf);
    return -1;
  }
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    strncpy(temp,line,MAXLINE);
    if((p=strchr(temp,'#')) != NULL) *p = '\0';
    // read line
    for(n=nc=0,p=temp; n<MAXENTR; n++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
    }
    if(n < 1) continue;
    if(strcasecmp(str[0],"ssr_th_los_min") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        ssr_th_los_min = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n","ssr_th_los_min",ssr_th_los_min);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_th_los_max") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        ssr_th_los_max = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n","ssr_th_los_max",ssr_th_los_max);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sim_xsgm") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        sim_xsgm = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n","sim_xsgm",sim_xsgm);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sky_wlen") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(n > 3)
      {
        errno = 0;
        uni = strtod(str[3],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((sky_n_wlen=ReadDouble(str[1],num,uni,sky_wlen,SKY_MAXWLEN)) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30s %4d %13.4e\n","sky_wlen",str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sky_angl_min") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        sky_angl_min = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n","sky_angl_min",sky_angl_min);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sky_angl_max") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        sky_angl_max = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30.4e\n","sky_angl_max",sky_angl_max);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sky_mker") == 0)
    {
      if(n > 1)
      {
        strncpy(sky_mker,str[1],MAXLINE);
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30s\n","sky_mker",sky_mker);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sky_find") == 0)
    {
      if(n > 2)
      {
        errno = 0;
        sky_rind_wuni = strtod(str[2],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        strncpy(sky_find,str[1],MAXLINE);
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30s %13.4e\n","sky_find",sky_find,sky_rind_wuni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sky_frho") == 0)
    {
      if(n > 2)
      {
        errno = 0;
        sky_grho_wuni = strtod(str[2],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadConfig: convert error >>> %s\n",line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        strncpy(sky_frho,str[1],MAXLINE);
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30s %13.4e\n","sky_frho",sky_frho,sky_grho_wuni);
        cnt_n_cmnt++;
      }
    }
  }
  fclose(fp);
  if(err)
  {
    return -1;
  }

  return 0;
}

int ReadDouble(char *s,int c,double u,double *d,int size)
{
  int i,j,n;
  int nc,nd;
  int err;
  int flag;
  double vtmp;
  char line[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char *p;
  FILE *fp;

  // check input
  if(c<0 || c>=MAXENTR)
  {
    fprintf(stderr,"ReadDouble: invalid column number >>> %d (%s)\n",c,s);
    return -1;
  }
  // open input
  flag = 0;
  if(strcmp(s,"") == 0)
  {
    fp = stdin;
  }
  else
  {
    if((fp=fopen(s,"r")) == NULL)
    {
      fprintf(stderr,"ReadDouble: cannot open %s\n",s);
      return -1;
    }
    flag = 1;
  }
  // read input
  err = 0;
  i = n = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    for(j=0,p=line; j<=c; j++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[j],&nc) == EOF)
      {
        fprintf(stderr,"ReadDouble: read error >>> %s (%s)\n",line,s);
        err = 1;
        break;
      }
    }
    if(err) break;
    errno = 0;
    vtmp = strtod(str[c],&p)*u;
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"ReadDouble: convert error >>> %s (%s)\n",line,s);
      err = 1;
      break;
    }
    if(n >= size)
    {
      fprintf(stderr,"ReadDouble: warning, #data exceed the limit %d (%s)\n",n,s);
      break;
    }
    d[n] = vtmp;
    n++;
    i++;
  }
  if(flag)
  {
    fclose(fp);
  }
  nd = n;
  if(nd < 1)
  {
    fprintf(stderr,"ReadDouble: error, nd=%d (%s)\n",nd,s);
    err = 1;
  }
  if(err)
  {
    return -1;
  }

  return nd;
}

double SepAngle2(double th1,double ph1,double th2,double ph2)
{
  return acos(sin(th1)*sin(th2)*cos(ph1-ph2)+cos(th1)*cos(th2));
}

int Interp1D(const double *x1,const double *z1,int nx1,
             const double *x2,      double *z2,int nx2,int f)
{
  // x1[nx1],z1[nx1]
  // x2[nx2],z2[nx2]
  // x1 must be arranged in ascending order
  int i,j;
  static int flag = 1;
  static int *i_1 = NULL;
  static int *i_2 = NULL;

  if(flag==1 || f==1)
  {
    if(i_1 != NULL) free(i_1);
    if(i_2 != NULL) free(i_2);
    i_1 = (int*)malloc(nx2*sizeof(int));
    i_2 = (int*)malloc(nx2*sizeof(int));
    if(i_1==NULL || i_2==NULL)
    {
      fprintf(stderr,"Interp1D: error in allocating memory\n");
      return -1;
    }
    for(i=0; i<nx2; i++)
    {
      i_1[i] = i_2[i] = -1;
      for(j=nx1-1; j>=0; j--)
      {
        if(x1[j] <= x2[i])
        {
          i_1[i] = j;
          break;
        }
      }
      for(j=0; j<nx1; j++)
      {
        if(x1[j] >= x2[i])
        {
          i_2[i] = j;
          break;
        }
      }
      if(i_1[i] < 0)
      {
        if(nx1 < 2)
        {
          i_1[i] = 0;
          i_2[i] = 0;
        }
        else
        {
          i_1[i] = 0;
          i_2[i] = 1;
        }
      }
      if(i_2[i] < 0)
      {
        if(nx1 < 2)
        {
          i_1[i] = 0;
          i_2[i] = 0;
        }
        else
        {
          i_1[i] = nx1-2;
          i_2[i] = nx1-1;
        }
      }
      if(cnt_db)
      {
        fprintf(stderr,"Interp1D: i_1[%d]=%d, i_2[%d]=%d\n",i,i_1[i],i,i_2[i]);
      }
    }
    flag = 0;
  }

  for(i=0; i<nx2; i++)
  {
    if(fabs(x1[i_2[i]]-x1[i_1[i]]) < EPSILON)
    {
      // NO INTERPOLATION IS NECESSARY
      z2[i] = z1[i_1[i]];
    }
    else
    {
      z2[i] = INTERP(x2[i],x1[i_1[i]],x1[i_2[i]],z1[i_1[i]],z1[i_2[i]]);
    }
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
    yout[i] = (yinp[i1]+(yinp[i2]-yinp[i1])*(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1]))*yuni;
  }

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
    {"inp_dnam",1,0,'D'},
    {"inp_ndig",1,0,'N'},
    {"inp_tdif",1,0,'T'},
    {"sky_timz",1,0,'O'},
    {"sky_dnum",1,0,'n'},
    {"sky_dtop",1,0,'t'},
    {"sky_ozon",1,0,'o'},
    {"sky_pres",1,0,'p'},
    {"sky_lpmx",1,0,'i'},
    {"sky_epsa",1,0,'s'},
    {"sky_eps1",1,0,'e'},
    {"sky_eps2",1,0,'E'},
    {"cnt_conf",0,0,'C'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":D:N:T:O:n:t:o:p:i:s:e:E:C:dvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'D':
        strncpy(inp_dnam,optarg,MAXLINE);
        break;
      case 'N':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0 && ntmp<10) inp_ndig = ntmp;
        else
        {
          fprintf(stderr,"#Digits in run# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'T':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) inp_tdif = ntmp;
        else
        {
          fprintf(stderr,"Time difference -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'O':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=-2400 && ntmp<=2400) sky_timz = ntmp;
        else
        {
          fprintf(stderr,"Timezone -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'n':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0) sky_dnum = ntmp;
        else
        {
          fprintf(stderr,"Data number -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 't':
        strncpy(sky_dtop,optarg,MAXLINE);
        break;
      case 'o':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>=0.0 && xtmp<1.0) sky_ozon = xtmp;
        else
        {
          fprintf(stderr,"Ozone amount -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'p':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>=0.0 && xtmp<10.0) sky_pres = xtmp;
        else
        {
          fprintf(stderr,"Pressure -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'i':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) sky_lpmx = ntmp;
        else
        {
          fprintf(stderr,"Max #iteration -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 's':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>=0.0) sky_epsa = xtmp;
        else
        {
          fprintf(stderr,"EPSA -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'e':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>=0.0) sky_eps1 = xtmp;
        else
        {
          fprintf(stderr,"EPSA1 -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'E':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>=0.0) sky_eps2 = xtmp;
        else
        {
          fprintf(stderr,"EPSA2 -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'C':
        strncpy(cnt_conf,optarg,MAXLINE);
        break;
      case 'd':
        cnt_db++;
        break;
      case 'v':
        cnt_vb++;
        break;
      case 'h':
        cnt_hp = 1;
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

  if(ReadConfig() < 0)
  {
    rt = -1;
  }

  if(cnt_hp) rt = 1;
  return rt;
}

int Usage(void)
{
  int  n = 15;
  char e[MAXLINE];
  char a[MAXLINE];
  char d[MAXLINE];

  fprintf(stderr,"ms720_skyrad ... prepare files for processing ms720 data with skyrad.PACK_V42.\n");
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"ms720_skyrad -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"----------------------------------------------------------------------------------\n");
  fprintf(stderr,"      option     |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," D -inp_dnam     |%s|%s|%s| %s\n",As(e,"Data directory",n),As(a,"name",n),       As(d,INP_DNAM,n),inp_dnam);
  fprintf(stderr," N -inp_ndig     |%s|%s|%s| %d\n",As(e,"RNO #Digits",n),   As(a,"#",n),          Ad(d,INP_NDIG,n),inp_ndig);
  fprintf(stderr," T -inp_tdif     |%s|%s|%s| %d\n",As(e,"Time diff",n),     As(a,"sec",n),        Ad(d,INP_TDIF,n),inp_tdif);
  fprintf(stderr," O -sky_timz     |%s|%s|%s| %d\n",As(e,"Timezone",n),      As(a,"hour*100",n),   Ad(d,SKY_TIMZ,n),sky_timz);
  fprintf(stderr," n -sky_dnum     |%s|%s|%s| %d\n",As(e,"Data number",n),   As(a,"number",n),     Ad(d,SKY_DNUM,n),sky_dnum);
  fprintf(stderr," t -sky_dtop     |%s|%s|%s| %s\n",As(e,"Top directory",n), As(a,"name",n),       As(d,SKY_DTOP,n),sky_dtop);
  fprintf(stderr," o -sky_ozon     |%s|%s|%s| %e\n",As(e,"Ozone amount",n),  As(a,"cm,STP",n),     Af(d,SKY_OZON,n),sky_ozon);
  fprintf(stderr," p -sky_pres     |%s|%s|%s| %e\n",As(e,"Pressure",n),      As(a,"atm",n),        Af(d,SKY_PRES,n),sky_pres);
  fprintf(stderr," i -sky_lpmx     |%s|%s|%s| %d\n",As(e,"Max #iteration",n),As(a,"#",n),          Ad(d,SKY_LPMX,n),sky_lpmx);
  fprintf(stderr," s -sky_epsa     |%s|%s|%s| %e\n",As(e,"EPSA",n),          As(a,"value",n),      Ae(d,SKY_EPSA,n),sky_epsa);
  fprintf(stderr," e -sky_eps1     |%s|%s|%s| %e\n",As(e,"EPSA1",n),         As(a,"value",n),      Ae(d,SKY_EPS1,n),sky_eps1);
  fprintf(stderr," E -sky_eps2     |%s|%s|%s| %e\n",As(e,"EPSA2",n),         As(a,"value",n),      Ae(d,SKY_EPS2,n),sky_eps2);
  fprintf(stderr," C -cnt_conf     |%s|%s|%s| %s\n",As(e,"Config  file",n),  As(a,"name",n),       As(d,CNT_CONF,n),cnt_conf);
  fprintf(stderr," d -debug        |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),cnt_db);
  fprintf(stderr," v -verbose      |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),cnt_vb);
  fprintf(stderr," h -help         |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"----------------------------------------------------------------------------------\n");
  fprintf(stderr,"Config file (cnt_conf) format:\n");
  fprintf(stderr,"ssr_th_los_min angle       | min th(los) in degree(%.1f)\n",SSR_TH_LOS_MIN);
  fprintf(stderr,"ssr_th_los_max angle       | max th(los) in degree(%.1f)\n",SSR_TH_LOS_MAX);
  fprintf(stderr,"sim_xsgm       value       | resolution in nm(%.1f)\n",SIM_XSGM);
  fprintf(stderr,"sky_wlen       name # unit | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
  fprintf(stderr,"sky_angl_min   angle       | min scattering angle in degree(%.1f)\n",SKY_ANGL_MIN);
  fprintf(stderr,"sky_angl_max   angle       | max scattering angle in degree(%.1f)\n",SKY_ANGL_MAX);
  fprintf(stderr,"sky_mker       name        | Mie kernel(%s)\n",SKY_MKER);
  fprintf(stderr,"sky_find       name   unit | file name(%s),unit in nm(%.1f)\n",NONAME,SKY_RIND_WUNI);
  fprintf(stderr,"sky_frho       name   unit | file name(%s),unit in nm(%.1f)\n",NONAME,SKY_GRHO_WUNI);
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
