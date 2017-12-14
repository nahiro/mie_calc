/*********************************************************/
/* aeros_comp                                            */
/* Author: N.Manago Apr,30,2008                          */
/* $Revision: 788 $                                      */
/* $Date: 2014-09-01 01:36:09 +0900 (Mon, 01 Sep 2014) $ */
/*********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <bits/nan.h>
#include <complex.h>
#include <time.h>
#include <sys/time.h>
#include "resolution.h"
#include "aeros_profile.h"
#include "strutil.h"

// Mathematical constants
#define	PI				3.141592653589793	// PI
#define	PI2				6.283185307179586	// 2*PI
#define	R_TO_D				57.29577951308232	// 180/PI rad -> deg
#define	D_TO_R				1.745329251994329e-02	// PI/180 deg -> rad
// Physical constants
#define	SEC_DAY				86400
// Common constants
#define	NONAME				"NONE"			// No name
#define	NODATA				0			// No data
#define	INVALID				-1			// Invalid value
#define	MAXENTR				25			// Max #entries in a line
#define	MAXLINE				256			// Max #chars in a line
#define	EPSILON				1.0e-14			// A small number
#define	DELTA				1.0e-12			// A small number
// Constants for input data
#define	INP_MAXDATA			100			// Max #input data
#define	INP_NCOL			21			// #Columns
#define	INP_NINT			11			// #Ints
#define	INP_NDBL			10			// #Doubles
#define	INP_N_WLEN			256			// #Wavelengths
#define	INP_NDIG			4			// #Digits in run#
#define	INP_TDIF			600			// Time difference
#define	INP_DNAM			"../data"		// Data directory
#define	INP_TYPE_DSR			1			// DSR
#define	INP_TYPE_AUR			2			// AUR
#define	INP_TYPE_SSR			3			// SSR
#define	INP_TYPE_UNKNOWN		INVALID			// Unknwon
// Constants for error
#define	ERR_N_WLEN			17			// #Error range
// Constants for Direct Solar radiation
#define	DSR_MAXDATA			30			// Max #DSRs
#define	DSR_MAXRANG			4			// Max #ranges
#define	DSR_MAXWLEN			256			// Max #wavelengths for DSR
#define	DSR_MAXTDIV			10			// Max #divisions in theta
#define	DSR_N_TDIV			4			// #Divisions in theta
#define	DSR_N_PDIV			16			// #Divisions in phi
#define	DSR_N_WLEN_0			9
#define DSR_WLEN_0_A_001                535.694
#define DSR_WLEN_0_A_002                539.036
#define DSR_WLEN_0_A_003                542.378
#define DSR_WLEN_0_A_004                545.719
#define DSR_WLEN_0_A_005                549.061
#define DSR_WLEN_0_A_006                552.402
#define DSR_WLEN_0_A_007                555.743
#define DSR_WLEN_0_A_008                559.084
#define DSR_WLEN_0_A_009                562.425
#define	DSR_WLEN_0_B_001		535.897
#define	DSR_WLEN_0_B_002		539.253
#define	DSR_WLEN_0_B_003		542.609
#define	DSR_WLEN_0_B_004		545.964
#define	DSR_WLEN_0_B_005		549.320
#define	DSR_WLEN_0_B_006		552.675
#define	DSR_WLEN_0_B_007		556.031
#define	DSR_WLEN_0_B_008		559.386
#define	DSR_WLEN_0_B_009		562.741
#define	DSR_N_WLEN_1			7
#define	DSR_WLEN_1_A_001		478.890
#define	DSR_WLEN_1_A_002		535.694
#define	DSR_WLEN_1_A_003		579.128
#define	DSR_WLEN_1_A_004		669.209
#define	DSR_WLEN_1_A_005		748.988
#define	DSR_WLEN_1_A_006		775.492
#define	DSR_WLEN_1_A_007		874.349
#define	DSR_WLEN_1_B_001		478.844
#define	DSR_WLEN_1_B_002		535.897
#define	DSR_WLEN_1_B_003		579.513
#define	DSR_WLEN_1_B_004		669.956
#define	DSR_WLEN_1_B_005		746.725
#define	DSR_WLEN_1_B_006		779.990
#define	DSR_WLEN_1_B_007		875.948
#define	DSR_N_WLEN_2			12
#define	DSR_WLEN_2_A_001		709.144
#define	DSR_WLEN_2_A_002		712.468
#define	DSR_WLEN_2_A_003		715.791
#define	DSR_WLEN_2_A_004		719.114
#define	DSR_WLEN_2_A_005		722.436
#define	DSR_WLEN_2_A_006		725.758
#define	DSR_WLEN_2_A_007		729.079
#define	DSR_WLEN_2_A_008		732.399
#define	DSR_WLEN_2_A_009		735.718
#define	DSR_WLEN_2_A_010		739.037
#define	DSR_WLEN_2_A_011		742.355
#define	DSR_WLEN_2_A_012		745.672
#define	DSR_WLEN_2_B_001		710.050
#define	DSR_WLEN_2_B_002		713.388
#define	DSR_WLEN_2_B_003		716.724
#define	DSR_WLEN_2_B_004		720.060
#define	DSR_WLEN_2_B_005		723.396
#define	DSR_WLEN_2_B_006		726.731
#define	DSR_WLEN_2_B_007		730.065
#define	DSR_WLEN_2_B_008		733.398
#define	DSR_WLEN_2_B_009		736.731
#define	DSR_WLEN_2_B_010		740.063
#define	DSR_WLEN_2_B_011		743.394
#define	DSR_WLEN_2_B_012		746.725
#define	DSR_N_WLEN_3			7
#define	DSR_WLEN_3_A_001		508.961
#define	DSR_WLEN_3_A_002		552.402
#define	DSR_WLEN_3_A_003		582.469
#define	DSR_WLEN_3_A_004		605.843
#define	DSR_WLEN_3_A_005		635.876
#define	DSR_WLEN_3_A_006		669.209
#define	DSR_WLEN_3_A_007		709.144
#define	DSR_WLEN_3_B_001		509.049
#define	DSR_WLEN_3_B_002		552.675
#define	DSR_WLEN_3_B_003		582.867
#define	DSR_WLEN_3_B_004		602.984
#define	DSR_WLEN_3_B_005		633.141
#define	DSR_WLEN_3_B_006		669.956
#define	DSR_WLEN_3_B_007		706.712
#define	DSR_TMIN			0.0			// Min incident angle (in degree) for DSR
#define	DSR_TMAX			4.01255064		// Max incident angle (in degree) for DSR
#define	DSR_RFOV			2.5			// FOV radius (in degree) for DSR
#define	DSR_SFOV			5.98020017e-3		// FOV area (in sr) for DSR
#define	DSR_RERR			0.05
// Constants for Solar Aureole
#define	AUR_MAXDATA			30			// Max #AURs
#define	AUR_MAXRANG			2			// Max #ranges
#define	AUR_MAXWLEN			256			// Max #wavelengths for AUR
#define	AUR_MAXTDIV			10			// Max #divisions in theta
#define	AUR_N_TDIV			4			// #Divisions in theta
#define	AUR_N_PDIV			16			// #Divisions in phi
#define	AUR_TYPE_OLD			1			// AUR data by old method
#define	AUR_TYPE_NEW			2			// AUR data by new method
#define	AUR_N_WLEN_0			1
#define	AUR_WLEN_0_A			549.061
#define	AUR_WLEN_0_B			549.320
#define	AUR_N_WLEN_1			10
#define	AUR_WLEN_1_A_001		455.510
#define	AUR_WLEN_1_A_002		488.913
#define	AUR_WLEN_1_A_003		555.743
#define	AUR_WLEN_1_A_004		609.181
#define	AUR_WLEN_1_A_005		635.876
#define	AUR_WLEN_1_A_006		672.540
#define	AUR_WLEN_1_A_007		752.304
#define	AUR_WLEN_1_A_008		792.030
#define	AUR_WLEN_1_A_009		844.790
#define	AUR_WLEN_1_A_010		877.628
#define	AUR_WLEN_1_B_001		468.778
#define	AUR_WLEN_1_B_002		488.912
#define	AUR_WLEN_1_B_003		552.675
#define	AUR_WLEN_1_B_004		579.513
#define	AUR_WLEN_1_B_005		629.792
#define	AUR_WLEN_1_B_006		676.644
#define	AUR_WLEN_1_B_007		750.055
#define	AUR_WLEN_1_B_008		799.909
#define	AUR_WLEN_1_B_009		849.561
#define	AUR_WLEN_1_B_010		879.241
#define	AUR_SBMP			15
#define	AUR_SBMP_1_A_001		1
#define	AUR_SBMP_1_A_002		15
#define	AUR_SBMP_1_A_003		1
#define	AUR_SBMP_1_A_004		1
#define	AUR_SBMP_1_A_005		1
#define	AUR_SBMP_1_A_006		1
#define	AUR_SBMP_1_A_007		15
#define	AUR_SBMP_1_A_008		1
#define	AUR_SBMP_1_A_009		1
#define	AUR_SBMP_1_A_010		1
#define	AUR_SBMP_1_B_001		15
#define	AUR_SBMP_1_B_002		15
#define	AUR_SBMP_1_B_003		15
#define	AUR_SBMP_1_B_004		15
#define	AUR_SBMP_1_B_005		15
#define	AUR_SBMP_1_B_006		15
#define	AUR_SBMP_1_B_007		15
#define	AUR_SBMP_1_B_008		15
#define	AUR_SBMP_1_B_009		15
#define	AUR_SBMP_1_B_010		1
#define	AUR_TMIN			0.97474186		// Min incident angle (in degree) for AUR
#define	AUR_TMAX			12.06552397		// Max incident angle (in degree) for AUR
#define	AUR_RFOV			6.25			// FOV radius (in degree) for AUR
#define	AUR_SFOV			8.94755029e-2		// FOV area (in sr) for AUR
#define	AUR_RERR			0.05
// Constants for Scattered Solar Radiation
#define	SSR_MAXDATA			30			// Max #SSRs
#define	SSR_MAXRANG			4			// Max #ranges
#define	SSR_MAXWLEN			256			// Max #wavelengths for SSR
#define	SSR_MAXTDIV			10			// Max #divisions in theta
#define	SSR_N_TDIV			4			// #Divisions in theta
#define	SSR_N_PDIV			16			// #Divisions in phi
#define	SSR_N_WLEN_0			1
#define	SSR_WLEN_0_A			549.061
#define	SSR_WLEN_0_B			549.320
#define	SSR_N_WLEN_1			10
#define	SSR_WLEN_1_A_001		455.510
#define	SSR_WLEN_1_A_002		488.913
#define	SSR_WLEN_1_A_003		555.743
#define	SSR_WLEN_1_A_004		609.181
#define	SSR_WLEN_1_A_005		635.876
#define	SSR_WLEN_1_A_006		672.540
#define	SSR_WLEN_1_A_007		752.304
#define	SSR_WLEN_1_A_008		792.030
#define	SSR_WLEN_1_A_009		844.790
#define	SSR_WLEN_1_A_010		877.628
#define	SSR_WLEN_1_B_001		468.778
#define	SSR_WLEN_1_B_002		488.912
#define	SSR_WLEN_1_B_003		552.675
#define	SSR_WLEN_1_B_004		579.513
#define	SSR_WLEN_1_B_005		629.792
#define	SSR_WLEN_1_B_006		676.644
#define	SSR_WLEN_1_B_007		750.055
#define	SSR_WLEN_1_B_008		799.909
#define	SSR_WLEN_1_B_009		849.561
#define	SSR_WLEN_1_B_010		879.241
#define	SSR_N_WLEN_2			12
#define	SSR_WLEN_2_A_001		709.144
#define	SSR_WLEN_2_A_002		712.468
#define	SSR_WLEN_2_A_003		715.791
#define	SSR_WLEN_2_A_004		719.114
#define	SSR_WLEN_2_A_005		722.436
#define	SSR_WLEN_2_A_006		725.758
#define	SSR_WLEN_2_A_007		729.079
#define	SSR_WLEN_2_A_008		732.399
#define	SSR_WLEN_2_A_009		735.718
#define	SSR_WLEN_2_A_010		739.037
#define	SSR_WLEN_2_A_011		742.355
#define	SSR_WLEN_2_A_012		745.672
#define	SSR_WLEN_2_B_001		710.050
#define	SSR_WLEN_2_B_002		713.388
#define	SSR_WLEN_2_B_003		716.724
#define	SSR_WLEN_2_B_004		720.060
#define	SSR_WLEN_2_B_005		723.396
#define	SSR_WLEN_2_B_006		726.731
#define	SSR_WLEN_2_B_007		730.065
#define	SSR_WLEN_2_B_008		733.398
#define	SSR_WLEN_2_B_009		736.731
#define	SSR_WLEN_2_B_010		740.063
#define	SSR_WLEN_2_B_011		743.394
#define	SSR_WLEN_2_B_012		746.725
#define	SSR_N_WLEN_3			7
#define	SSR_WLEN_3_A_001		508.961
#define	SSR_WLEN_3_A_002		552.402
#define	SSR_WLEN_3_A_003		582.469
#define	SSR_WLEN_3_A_004		605.843
#define	SSR_WLEN_3_A_005		635.876
#define	SSR_WLEN_3_A_006		669.209
#define	SSR_WLEN_3_A_007		709.144
#define	SSR_WLEN_3_B_001		509.049
#define	SSR_WLEN_3_B_002		552.675
#define	SSR_WLEN_3_B_003		582.867
#define	SSR_WLEN_3_B_004		602.984
#define	SSR_WLEN_3_B_005		633.141
#define	SSR_WLEN_3_B_006		669.956
#define	SSR_WLEN_3_B_007		706.712
#define	SSR_SBMP			15
#define	SSR_SBMP_1_A_001		1
#define	SSR_SBMP_1_A_002		15
#define	SSR_SBMP_1_A_003		1
#define	SSR_SBMP_1_A_004		1
#define	SSR_SBMP_1_A_005		1
#define	SSR_SBMP_1_A_006		1
#define	SSR_SBMP_1_A_007		15
#define	SSR_SBMP_1_A_008		1
#define	SSR_SBMP_1_A_009		1
#define	SSR_SBMP_1_A_010		1
#define	SSR_SBMP_1_B_001		15
#define	SSR_SBMP_1_B_002		15
#define	SSR_SBMP_1_B_003		15
#define	SSR_SBMP_1_B_004		15
#define	SSR_SBMP_1_B_005		15
#define	SSR_SBMP_1_B_006		15
#define	SSR_SBMP_1_B_007		15
#define	SSR_SBMP_1_B_008		15
#define	SSR_SBMP_1_B_009		15
#define	SSR_SBMP_1_B_010		1
#define	SSR_TH_LOS_MIN			-90.0
#define	SSR_TH_LOS_MAX			90.0
#define	SSR_SEP_MIN			12.5
#define	SSR_SEP_MAX			180.0
#define	SSR_TMIN			0.0			// Min incident angle (in degree) for SSR
#define	SSR_TMAX			12.06552397		// Max incident angle (in degree) for SSR
#define	SSR_RFOV			10.0			// FOV radius (in degree) for SSR
#define	SSR_SFOV			9.54557030e-2		// FOV area (in sr) for SSR
#define	SSR_RERR			0.05
// Constants for function
#define	FCN_MODE_AER			0			// Aerosol
#define	FCN_MODE_H2O			1			// H2O
#define	FCN_MODE_OZN			2			// Ozone
// Constants for Simulation
#define	SIM_N_ANGL			50			// #Angles
#define	SIM_MAXDATA			50000			// #Data
#define	SIM_MAXCONV			1000			// #Data
#define	SIM_NCOL			8			// #Columns
#define	SIM_NINT			2			// #Ints
#define	SIM_NDBL			6			// #Doubles
#define	SIM_NSGL			6			// #Columns
#define	SIM_FLAG			2			// Simulation flag
#define	SIM_PRNT			1			// NOPRNT flag
#define	SIM_IATM			2			// Atmosphere modle number
#define	SIM_ISEA			1			// Season modle number
#define	SIM_IVUL			0			// Volcanic modle number
#define	SIM_MULT			1			// Multiple scattering flag
#define	SIM_SALB			53			// Spectral albedo#
#define	SIM_IDAY			-1			// Day of year
#define	SIM_IHAZ			1			// Aerosol model#
#define	SIM_DBMP			1			// DSR Band model number
#define	SIM_SBMP			-1			// SSR Band model number
#define	SIM_CFLG			0			// Correction flag
#define	SIM_DFLG			1			// DISORT flag
#define	SIM_WFLG			0			// H2OAER flag
#define	SIM_NSAM			4			// #Sampling
#define	SIM_IM				0			// CARD2C read flag
#define	SIM_2C				0			// CARD2C flag
#define	SIM_MIXR_CO2			365.0			// CO2 mixing ratio in ppmv
#define	SIM_NPAR			37			// #Parameters
#define	SIM_VTAU			0.36777			// Aerosol optical depth
#define	SIM_VMIE			20.0			// Visibility in km
#define	SIM_WSCL			1.0			// Water scale
#define	SIM_OSCL			1.0			// Ozone scale
#define	SIM_GRHO			0.3			// Ground albedo
#define	SIM_GALT			0.0			// Ground altitude
#define	SIM_XMGN			16.0			// X margin
#define	SIM_XSGM			4.0			// X sigma
#define	SIM_WSGM			20.0			// X width
#define	SIM_YUNI			10.0			// Y unit
#define	SIM_DMIN			1.0e-12			// Min parameter difference
#define	SIM_DMAX			10.0			// Max wavelength difference in nm
#define	SIM_MODE_DSR_MODTRAN_BO		1			// Simulation mode
#define	SIM_MODE_DSR_MODTRAN_BA		2			// Simulation mode
#define	SIM_MODE_DSR_CUSTOM_BO		3			// Simulation mode
#define	SIM_MODE_DSR_CUSTOM_BA		4			// Simulation mode
#define	SIM_MODE_AUR_MODTRAN_SA		102			// Simulation mode
#define	SIM_MODE_AUR_MODTRAN_BA		104			// Simulation mode
#define	SIM_MODE_AUR_CUSTOM_SA		106			// Simulation mode
#define	SIM_MODE_AUR_CUSTOM_BA		108			// Simulation mode
#define	SIM_MODE_SSR_MODTRAN_SO		201			// Simulation mode
#define	SIM_MODE_SSR_MODTRAN_SA		202			// Simulation mode
#define	SIM_MODE_SSR_MODTRAN_BO		203			// Simulation mode
#define	SIM_MODE_SSR_MODTRAN_BA		204			// Simulation mode
#define	SIM_MODE_SSR_CUSTOM_SO		205			// Simulation mode
#define	SIM_MODE_SSR_CUSTOM_SA		206			// Simulation mode
#define	SIM_MODE_SSR_CUSTOM_BO		207			// Simulation mode
#define	SIM_MODE_SSR_CUSTOM_BA		208			// Simulation mode
#define	SIM_MODTRAN_VTAU		1			// Simulation mode
#define	SIM_MODTRAN_VMIE		2			// Simulation mode
#define	SIM_PATH			"/usr/local/Mod4v2r1"	// Path
#define	SIM_FALB			"DATA/spec_alb.dat"	// Spectral albedo file
// Constants for Profile
#define	PRF_HAZ_NALT			34
#define	PRF_NALT			36
#define	PRF_MAXNALT			1000
#define	PRF_NMOL			12
#define	PRF_N_AMOL			7
#define	PRF_N_TRAC			5
#define	PRF_MOL_H2O			0
#define	PRF_MOL_CO2			1
#define	PRF_MOL_O3			2
#define	PRF_MOL_N2O			3
#define	PRF_MOL_CO			4
#define	PRF_MOL_CH4			5
#define	PRF_MOL_O2			6
#define	PRF_MOL_NO			7
#define	PRF_MOL_SO2			8
#define	PRF_MOL_NO2			9
#define	PRF_MOL_NH3			10
#define	PRF_MOL_HNO3			11
#define	PRF_NHAZ			4
#define	PRF_NHAZ_LOW			2
#define	PRF_NPAR			14
// Constants for Optimization
#define	OPT_N_COMP			3			// #Components
#define	OPT_NPAR			16			// #Parameters
#define	OPT_PAR_VTAU			0
#define	OPT_PAR_VMIE			1
#define	OPT_PAR_WSCL			2
#define	OPT_PAR_OSCL			3
#define	OPT_PAR_XSCL			4
#define	OPT_PAR_YSCL			5
#define	OPT_PAR_ZSCL			6
#define	OPT_PAR_NCMP			7
#define	OPT_PAR_LMD1			8
#define	OPT_PAR_LSG1			9
#define	OPT_PAR_MXR2			10
#define	OPT_PAR_LMD2			11
#define	OPT_PAR_LSG2			12
#define	OPT_PAR_MXR3			13
#define	OPT_PAR_LMD3			14
#define	OPT_PAR_LSG3			15
#define	OPT_VTAU			SIM_VTAU
#define	OPT_VMIE			SIM_VMIE
#define	OPT_WSCL			SIM_WSCL
#define	OPT_OSCL			SIM_OSCL
#define	OPT_XSCL			1.0
#define	OPT_YSCL			1.0
#define	OPT_ZSCL			1.0
#define	OPT_NCMP			3.0
#define	OPT_MIXR			0.0
#define	OPT_LMOD			0.3
#define	OPT_LSGM			0.3
#define	OPT_N_FPAR			14
#define	OPT_FPAR_PROF			0			// Vertical profiles
#define	OPT_FPAR_PMOD			1			// Modtran parameters
#define	OPT_FPAR_MIXR			2			// Mixing ratios or size distributions
#define	OPT_FPAR_DIST			3			// Size distributions
#define	OPT_FPAR_COMP			4			// Size distribution of each component (1-10)
// Constants for Mie calculation
#define	MIE_N_WLEN			11			// #Wavelengths
#define	MIE_N_ANGL			121			// #Angles
#define	MIE_MAXDATA			1000			// Max #data
#define	MIE_IMIN			0			// Min line#
#define	MIE_IMAX			1000000000		// Min line#
#define	MIE_REAL_NUM			1			// Ref. index (real) column#
#define	MIE_IMAG_NUM			2			// Ref. index (imag) column#
#define	MIE_WLEN_REF			550.0			// Reference wavelength in nm
#define	MIE_WLEN_MIN			250.0			// Min wavelength in nm
#define	MIE_WLEN_MAX			1500.0			// Max wavelength in nm
#define	MIE_ANGL_MIN			0.0			// Min angle in degree
#define	MIE_ANGL_MAX			180.0			// Max angle in degree
#define	MIE_MAXCOMP			10			// Max #components
#define	MIE_TRUE			1			// True  value
#define	MIE_FALSE			0			// False value
#define	MIE_XMAX			2.0e4			// Max size parameter
#define	MIE_NSTP			1000			// Log10(R) #steps
#define	MIE_NMAX			10000			// Max Log10(R) #steps
#define	MIE_RMIN			NAN			// R min in um
#define	MIE_RMAX			NAN			// R max in um
#define	MIE_LRAD_MIN			-4.0			// Min Log10(R)
#define	MIE_LRAD_MAX			+3.0			// Max Log10(R)
#define	MIE_LMOD_MIN			-3.5			// Min Log10(R)
#define	MIE_LMOD_MAX			+2.0			// Max Log10(R)
#define	MIE_WSGM			5.0			// Log10(R) width in sigma
#define	MIE_INP1			"mie_wlen.dat"		// Inp. file 1
#define	MIE_INP2			"mie_angl.dat"		// Inp. file 2
#define	MIE_INP3			"mie_refx.dat"		// Inp. file 3
#define	MIE_INP4			"sim_pval.dat"		// Inp. file 4
#define	MIE_INP5			"mie_comp.dat"		// Inp. file 5
#define	MIE_OUT1			"mie_out1.dat"		// Out. file 1
#define	MIE_OUT2			"mie_out2.dat"		// Out. file 2
// Constants for T-matrix calculation
#define	TMX_NPN1			200
#define	TMX_NPNG1			(3*TMX_NPN1)
#define	TMX_NPNG2			(2*TMX_NPNG1)
#define	TMX_NPN2			(2*TMX_NPN1)
#define	TMX_NPL				(TMX_NPN2+1)
#define	TMX_NPN3			(TMX_NPN1+1)
#define	TMX_NPN4			(TMX_NPN1-25)
#define	TMX_NPN5			(2*TMX_NPN4)
#define	TMX_NPN6			(TMX_NPN4+1)
#define	TMX_NPL1			(TMX_NPN5+1)
#define	TMX_NAXMAX			10
#define	TMX_NAMAX			500
#define	TMX_NAFMAX			(TMX_NAMAX*6)
#define	TMX_NALMAX			(TMX_NAXMAX*TMX_NPL)
#define	TMX_NARMAX			(TMX_NAXMAX*TMX_NAFMAX)
#define	TMX_NDIS			0			// Distribution#
#define	TMX_NANG			361			// #Angles
#define	TMX_NDGS			2			// Control #divisions
#define	TMX_KMAX			100			// #Quad points
#define	TMX_SHAP			-1			// Shape#
#define	TMX_REPS			1.000001		// Aspect ratio or deformation parameter
#define	TMX_DELT			0.001			// Accuracy of the computations
#define	TMX_RMAX			8.0			// Max radius in um
// Constants for upper atmosphere
#define	UPP_N_WLEN			6
#define	UPP_N_ANGL			34
#define	UPP_N_PHAS			204
#define	UPP_N_RHUM			4
// Constants for control
#define	CNT_CONF			NONAME
#define	CNT_MAXCMNT			30
#define	CNT_MAXFORM			30
#define	CNT_MAXNOPT			100			// Max #options
// Functions
#define	MAX(a,b)			((b)>(a)?(b):(a))
#define	INTERP(x,x1,x2,y1,y2)		(y1+(y2-y1)/(x2-x1)*(x-x1))

// Parameters for input data
int	inp_ndig			= INP_NDIG;		// #Digits in run#
int	inp_tdif			= INP_TDIF;		// Time difference
int	inp_n_data			= NODATA;		// #Input data
int	inp_rno[INP_MAXDATA];					// Run#
int	inp_ms720[INP_MAXDATA];					// MS720ID
int	inp_cloud[INP_MAXDATA];					// Cloud flag
int	inp_flag[INP_MAXDATA];					// Flag
int	inp_type[INP_MAXDATA];					// Type
time_t	inp_time[INP_MAXDATA];					// Time
double	inp_ph_los[INP_MAXDATA];				// LOS azimuth (N=0,CW)
double	inp_th_los[INP_MAXDATA];				// LOS zenith
double	inp_temp[INP_MAXDATA];					// Temperature
double	inp_ph_sun[INP_MAXDATA];				// Solar azimuth (N=0,CW)
double	inp_th_sun[INP_MAXDATA];				// Solar zenith
double	inp_amas[INP_MAXDATA];					// Airmass
double	inp_dist[INP_MAXDATA];					// Sun-Observer distance
double	inp_dsr_wlen[INP_N_WLEN];				// Wavelength
double	inp_aur_wlen[INP_N_WLEN];				// Wavelength
double	inp_ssr_wlen[INP_N_WLEN];				// Wavelength
double	inp_data[INP_MAXDATA][INP_N_WLEN];			// Input data
char	inp_dnam[MAXLINE]		= INP_DNAM;		// Data directory
// Parameters for error
double	err_wmin[ERR_N_WLEN] =
{
  250.0,  350.0,  400.0,  450.0,  500.0,  550.0,
  600.0,  650.0,  700.0,  750.0,  800.0,  850.0,
  900.0,  950.0, 1000.0, 1050.0, 1100.0,
};
double	err_wmax[ERR_N_WLEN] =
{
  350.0,  400.0,  450.0,  500.0,  550.0,  600.0,
  650.0,  700.0,  750.0,  800.0,  850.0,  900.0,
  950.0, 1000.0, 1050.0, 1100.0, 1200.0,
};
// Parameters for Direct Solar radiation
int	dsr_ms720			= INVALID;		// DSR MS720ID
int	dsr_n_data			= NODATA;		// #DSR data
int	dsr_n_tdiv			= DSR_N_TDIV;		// #Divisions in theta
int	dsr_n_pdiv			= DSR_N_PDIV;		// #Divisions in phi
int	dsr_time[DSR_MAXDATA];					// Time
int	dsr_n_wlen[DSR_MAXRANG];				// #Wavelengths
int	dsr_nw[DSR_MAXRANG][DSR_MAXWLEN];			// Wavelength number
double	dsr_wlen[DSR_MAXRANG][DSR_MAXWLEN];			// Wavelength
double	dsr_data[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];	// Input data
double	dsr_dorg[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];	// Input data (no correction)
double	dsr_derr[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];	// Input error
double	dsr_dsim[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];	// Simulation data
double	dsr_dsgl[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];	// Simulation data (single scattering)
double	dsr_dsrc[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];	// Source function
double	dsr_dfac[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];	// DIS=f->t correction factor
double	dsr_dctr[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];	// Center radiance
double	dsr_unif[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];	// Uniformity correction factor
double	dsr_th_sun[DSR_MAXDATA];				// Solar zenith
double	dsr_ph_sun[DSR_MAXDATA];				// Solar azimuth (N=0,CW)
double	dsr_amas[DSR_MAXDATA];					// Airmass
double	dsr_dist[DSR_MAXDATA];					// Sun-Observer distance
double	dsr_tmin			= DSR_TMIN;		// Min incident angle (in degree) for DSR
double	dsr_tmax			= DSR_TMAX;		// Max incident angle (in degree) for DSR
double	dsr_rfov			= DSR_RFOV;		// FOV radius (in degree) for DSR
double	dsr_sfov			= DSR_SFOV;		// FOV area (in sr) for DSR
double	dsr_rerr			= DSR_RERR;		// Relative error
double	dsr_romg[DSR_MAXTDIV]		=			// Correction factor of Omega
{
  9.999028e-01,
  8.165045e-01,
  4.212483e-01,
  9.264020e-02,
};

// Parameters for Solar Aureole
int	aur_ms720			= INVALID;		// AUR MS720ID
int	aur_n_data			= NODATA;		// #AUR data
int	aur_n_tdiv			= AUR_N_TDIV;		// #Divisions in theta
int	aur_n_pdiv			= AUR_N_PDIV;		// #Divisions in phi
int	aur_type			= INVALID;		// AUR data type
int	aur_time[AUR_MAXDATA];					// Time
int	aur_ndsr[AUR_MAXDATA];					// DSR data# to be subtracted
int	aur_n_wlen[AUR_MAXRANG];				// #Wavelengths
int	aur_n_sbmp[AUR_MAXRANG];				// #Models
int	aur_nw[AUR_MAXRANG][AUR_MAXWLEN];			// Wavelength number
int	aur_sbmp[AUR_MAXRANG][AUR_MAXWLEN];			// AUR Band model number
double	aur_wlen[AUR_MAXRANG][AUR_MAXWLEN];			// Wavelength
double	aur_data[AUR_MAXRANG][AUR_MAXDATA][AUR_MAXWLEN];	// Input data
double	aur_dorg[AUR_MAXRANG][AUR_MAXDATA][AUR_MAXWLEN];	// Input data (no correction)
double	aur_derr[AUR_MAXRANG][AUR_MAXDATA][AUR_MAXWLEN];	// Input error
double	aur_dsim[AUR_MAXRANG][AUR_MAXDATA][AUR_MAXWLEN];	// Simulation data
double	aur_dsgl[AUR_MAXRANG][AUR_MAXDATA][AUR_MAXWLEN];	// Simulation data (single scattering)
double	aur_dfac[AUR_MAXRANG][AUR_MAXDATA][AUR_MAXWLEN];	// DIS=f->t correction factor
double	aur_dctr[AUR_MAXRANG][AUR_MAXDATA][AUR_MAXWLEN];	// Center radiance
double	aur_unif[AUR_MAXRANG][AUR_MAXDATA][AUR_MAXWLEN];	// Uniformity correction factor
double	aur_th_sun[AUR_MAXDATA];				// Solar zenith
double	aur_ph_sun[AUR_MAXDATA];				// Solar azimuth (N=0,CW)
double	aur_amas[AUR_MAXDATA];					// Airmass
double	aur_dist[AUR_MAXDATA];					// Sun-Observer distance
double	aur_th_los[AUR_MAXDATA];				// LOS zenith
double	aur_ph_los[AUR_MAXDATA];				// LOS azimuth (N=0,CW)
double	aur_tmin			= AUR_TMIN;		// Min incident angle (in degree) for AUR
double	aur_tmax			= AUR_TMAX;		// Max incident angle (in degree) for AUR
double	aur_rfov			= AUR_RFOV;		// FOV radius (in degree) for AUR
double	aur_sfov			= AUR_SFOV;		// FOV area (in sr) for AUR
double	aur_rerr			= AUR_RERR;		// Relative error
double	aur_romg[AUR_MAXTDIV]		=			// Correction factor of Omega
{
  5.968121e-01,
  9.946576e-01,
  9.185929e-01,
  2.783299e-01,
};

// Parameters for Scattered Solar Radiation
int	ssr_ms720			= INVALID;		// SSR MS720ID
int	ssr_n_data			= NODATA;		// #SSR data
int	ssr_n_tdiv			= SSR_N_TDIV;		// #Divisions in theta
int	ssr_n_pdiv			= SSR_N_PDIV;		// #Divisions in phi
int	ssr_time[SSR_MAXDATA];					// Time
int	ssr_n_wlen[SSR_MAXRANG];				// #Wavelengths
int	ssr_n_sbmp[SSR_MAXRANG];				// #Models
int	ssr_nw[SSR_MAXRANG][SSR_MAXWLEN];			// Wavelength number
int	ssr_sbmp[SSR_MAXRANG][SSR_MAXWLEN];			// SSR Band model number
double	ssr_wlen[SSR_MAXRANG][SSR_MAXWLEN];			// Wavelength
double	ssr_data[SSR_MAXRANG][SSR_MAXDATA][SSR_MAXWLEN];	// Input data
double	ssr_dorg[SSR_MAXRANG][SSR_MAXDATA][SSR_MAXWLEN];	// Input data (no correction)
double	ssr_derr[SSR_MAXRANG][SSR_MAXDATA][SSR_MAXWLEN];	// Input error
double	ssr_dsim[SSR_MAXRANG][SSR_MAXDATA][SSR_MAXWLEN];	// Simulation data
double	ssr_dsgl[SSR_MAXRANG][SSR_MAXDATA][SSR_MAXWLEN];	// Simulation data (single scattering)
double	ssr_dfac[SSR_MAXRANG][SSR_MAXDATA][SSR_MAXWLEN];	// DIS=f->t correction factor
double	ssr_dctr[SSR_MAXRANG][SSR_MAXDATA][SSR_MAXWLEN];	// Center radiance
double	ssr_unif[SSR_MAXRANG][SSR_MAXDATA][SSR_MAXWLEN];	// Uniformity correction factor
double	ssr_th_sun[SSR_MAXDATA];				// Solar zenith
double	ssr_ph_sun[SSR_MAXDATA];				// Solar azimuth (N=0,CW)
double	ssr_amas[SSR_MAXDATA];					// Airmass
double	ssr_dist[SSR_MAXDATA];					// Sun-Observer distance
double	ssr_th_los[SSR_MAXDATA];				// LOS zenith
double	ssr_ph_los[SSR_MAXDATA];				// LOS azimuth (N=0,CW)
double	ssr_th_los_min			= SSR_TH_LOS_MIN;	// Min LOS zenith in degree
double	ssr_th_los_max			= SSR_TH_LOS_MAX;	// Max LOS zenith in degree
double	ssr_sep_min			= SSR_SEP_MIN;		// Min separation in degree
double	ssr_sep_max			= SSR_SEP_MAX;		// Max separation in degree
double	ssr_tmin			= SSR_TMIN;		// Min incident angle (in degree) for SSR
double	ssr_tmax			= SSR_TMAX;		// Max incident angle (in degree) for SSR
double	ssr_rfov			= SSR_RFOV;		// FOV radius (in degree) for SSR
double	ssr_sfov			= SSR_SFOV;		// FOV area (in sr) for SSR
double	ssr_rerr			= SSR_RERR;		// Relative error
double	ssr_romg[SSR_MAXTDIV]		=			// Correction factor of Omega
{
  9.993073e-01,
  9.965383e-01,
  9.484138e-01,
  3.092412e-01,
};
// Parameters for correction
int	crc_n_comp			= NODATA;		// #Components
int	crc_n_dsim			= NODATA;		// #Simulations
double	crc_dsr_dsim[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];
double	crc_dsr_dsgl[DSR_MAXRANG][DSR_MAXDATA][DSR_MAXWLEN];
double	crc_aur_dsim[AUR_MAXRANG][AUR_MAXDATA][AUR_MAXWLEN];
double	crc_aur_dsgl[AUR_MAXRANG][AUR_MAXDATA][AUR_MAXWLEN];
double	crc_ssr_dsim[SSR_MAXRANG][SSR_MAXDATA][SSR_MAXWLEN];
double	crc_ssr_dsgl[SSR_MAXRANG][SSR_MAXDATA][SSR_MAXWLEN];
double	crc_init[SIM_NPAR];					// Initial values
double	crc_wcom_com[MIE_MAXCOMP];				// Weight
double	crc_lmod_com[MIE_MAXCOMP];				// Mode radius in um
double	crc_lsgm_com[MIE_MAXCOMP];				// Sigma Log10(R in um)
double	*crc_refr_com[MIE_MAXCOMP];				// Refractive index (real)
double	*crc_refi_com[MIE_MAXCOMP];				// Refractive index (imaginary)
char	crc_fsim[MAXLINE]		= NONAME;		// Simulation file
char	crc_ffac[MAXLINE]		= NONAME;		// Factor file
// Parameters for Simulation
int	sim_mode			= SIM_MODTRAN_VMIE;	// Simulation mode
int	sim_n_angl			= SIM_N_ANGL;		// #Angles
int	sim_flag			= SIM_FLAG;		// Simulation flag
int	sim_prnt			= SIM_PRNT;		// NOPRNT flag
int	sim_iatm			= SIM_IATM;		// Atmosphere model number
int	sim_isea			= SIM_ISEA;		// Season model number
int	sim_ivul			= SIM_IVUL;		// Volcanic model number
int	sim_mult			= SIM_MULT;		// Multiple scattering flag
int	sim_salb			= SIM_SALB;		// Spectral albedo
int	sim_iday			= SIM_IDAY;		// Day of year
int	sim_ihaz			= SIM_IHAZ;		// Aerosol model#
int	sim_dbmp			= SIM_DBMP;		// DSR Band model number
int	sim_sbmp			= SIM_SBMP;		// SSR Band model number
int	sim_sbmp_init			= SIM_SBMP;		// SSR Band model number
int	sim_cflg			= SIM_CFLG;		// Correction flag
int	sim_dflg			= SIM_DFLG;		// DISORT flag
int	sim_wflg			= SIM_WFLG;		// H2OAER flag
int	sim_nsam			= SIM_NSAM;		// #Sampling
int	sim_im				= SIM_IM;		// CARD2C read flag
int	sim_2c				= SIM_2C;		// CARD2C flag
double	sim_mixr_co2			= SIM_MIXR_CO2;		// CO2 mixing ratio in ppmv
double	sim_vmie			= SIM_VMIE;		// Visibility in km
double	sim_wscl			= SIM_WSCL;		// Water scale
double	sim_oscl			= SIM_OSCL;		// Ozone scale
double	sim_grho			= SIM_GRHO;		// Ground albedo
double	sim_galt			= SIM_GALT;		// Ground altitude
double	sim_xmgn			= SIM_XMGN;		// X margin
double	sim_xsgm			= SIM_XSGM;		// X sigma
double	sim_wsgm			= SIM_WSGM;		// X width in sigma
double	sim_yuni			= SIM_YUNI;		// Y unit
double	sim_angl[SIM_N_ANGL] =
{
    0.0,   0.5,   1.0,   1.5,   2.0,   2.5,   3.0,   3.5,
    4.0,   4.5,   5.0,   6.0,   7.0,   8.0,   9.0,  10.0,
   12.0,  14.0,  16.0,  18.0,  20.0,  24.0,  28.0,  32.0,
   36.0,  40.0,  50.0,  60.0,  70.0,  80.0,  90.0,  95.0,
  100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0,
  140.0, 145.0, 150.0, 155.0, 160.0, 164.0, 168.0, 172.0,
  176.0, 180.0,
};
double	sim_pval[SIM_NPAR];
double	sim_wlen[SIM_MAXDATA];
double	sim_data[SIM_MAXDATA];
double	sim_wlen_cnv[SIM_MAXCONV];
double	sim_data_cnv[SIM_MAXCONV];
double	*sim_phas			= NULL;
double	*sim_tropo_phas			= NULL;
double	*sim_strat_phas			= NULL;
double	*sim_meteo_phas			= NULL;
int	sim_prf_pscl[PRF_NMOL];
double	sim_prf_pres[PRF_MAXNALT];
double	sim_prf_temp[PRF_MAXNALT];
double	sim_prf_wmol[PRF_NMOL][PRF_MAXNALT];
double	sim_prf_haze[PRF_NHAZ][PRF_MAXNALT];
double	sim_prf_vscl[PRF_NMOL];
char	sim_prf_jchar[PRF_NPAR];
char	sim_prf_jcharx;
char	sim_path[MAXLINE]		= SIM_PATH;
char	sim_falb[MAXLINE]		= SIM_FALB;		// Spectral albedo file
char	sim_form[OPT_NPAR][MAXLINE]	=
{
  " %9.5f", // VTAU
  " %9.5f", // VMIE
  " %8.6f", // WSCL
  " %8.6f", // OSCL
  " %8.6f", // XSCL
  " %8.6f", // YSCL
  " %8.6f", // ZSCL
  " %3.1f", // NCMP
  " %9.6f", // LMD1
  " %8.6f", // LSG1
  " %9.6f", // MXR2
  " %9.6f", // LMD2
  " %8.6f", // LSG2
  " %9.6f", // MXR3
  " %9.6f", // LMD3
  " %8.6f", // LSG3
};
// Parameters for Profile
int	prf_n_alt			= PRF_NALT;
double	prf_alt[PRF_MAXNALT]		=
{
   0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,
   8.0,  9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,
  16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0,
  24.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0,
  60.0, 70.0, 80.0, 100.0,
};
double	prf_haz_alt[PRF_HAZ_NALT]	=
{
    0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,
    8.0,  9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,
   16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0,
   24.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 70.0,
  100.0, 99999.0,
};
double	prf_pres[PRF_MAXNALT];
double	prf_temp[PRF_MAXNALT];
double	prf_wmol[PRF_NMOL][PRF_MAXNALT];
double	prf_haze[PRF_NHAZ][PRF_MAXNALT];
double	prf_pres_gnd			= NAN;
double	prf_temp_gnd			= NAN;
double	prf_aod				= 0.0;
double	prf_aod_lo			= 0.0;
double	prf_aod_hi			= 0.0;
char	prf_jchar[PRF_NPAR];
char	prf_wmol_nam[PRF_NMOL][MAXLINE]	=
{
  "H2O", "CO2", "O3", "N2O", "CO", "CH4",
  "O2", "NO", "SO2", "NO2", "NH3", "HNO3",
};
// Parameters for Optimization
int	opt_n_opt			= INVALID;		// #Optimization
int	opt_n_comp			= NODATA;		// #Components
double	opt_init[OPT_NPAR];					// Initial values
double	opt_wcom_com[OPT_N_COMP];				// Weight
double	opt_lmod_com[OPT_N_COMP];				// Mode radius in um
double	opt_lsgm_com[OPT_N_COMP];				// Sigma Log10(R in um)
double	*opt_refr_com[OPT_N_COMP];				// Refractive index (real)
double	*opt_refi_com[OPT_N_COMP];				// Refractive index (imaginary)
// Parameters for Mie calculation
int	mie_iref			= -1;			// Ref wave line#
int	mie_n_wlen			= MIE_N_WLEN;		// #Wavelengths (output)
int	mie_n_angl			= MIE_N_ANGL;		// #Angles
int	mie_n_comp			= NODATA;		// #Components
int	mie_n_step			= MIE_NSTP;		// Log10(R) #steps
double	mie_wlen_ref			= MIE_WLEN_REF;		// Reference wavelength
double	mie_wlen_min			= MIE_WLEN_MIN;		// Min wavelength in nm
double	mie_wlen_max			= MIE_WLEN_MAX;		// Max wavelength in nm
double	mie_wlen[MIE_MAXDATA]		=
{
  300.0, 337.1, 400.0,  488.0,  514.5, 550.0,
  632.8, 694.3, 860.0, 1060.0, 1300.0,
};
double	*mie_wlen_um			= NULL;
double	*mie_aext			= NULL;
double	*mie_asca			= NULL;
double	*mie_asym			= NULL;
double	mie_angl_min			= MIE_ANGL_MIN;		// Min angle in degree
double	mie_angl_max			= MIE_ANGL_MAX;		// Max angle in degree
double	mie_angl[MIE_MAXDATA]		=
{
    0.0,   0.5,   1.0,   1.5,   2.0,   2.5,   3.0,   3.5,   4.0,   4.5,   5.0,
    5.5,   6.0,   6.5,   7.0,   7.5,   8.0,   8.5,   9.0,   9.5,  10.0,  12.0,
   14.0,  16.0,  18.0,  20.0,  22.0,  24.0,  26.0,  28.0,  30.0,  32.0,  34.0,
   36.0,  38.0,  40.0,  42.0,  44.0,  46.0,  48.0,  50.0,  52.0,  54.0,  56.0,
   58.0,  60.0,  62.0,  64.0,  66.0,  68.0,  70.0,  72.0,  74.0,  76.0,  78.0,
   80.0,  82.0,  84.0,  86.0,  88.0,  90.0,  92.0,  94.0,  96.0,  98.0, 100.0,
  102.0, 104.0, 106.0, 108.0, 110.0, 112.0, 114.0, 116.0, 118.0, 120.0, 122.0,
  124.0, 126.0, 128.0, 130.0, 132.0, 134.0, 136.0, 138.0, 140.0, 142.0, 144.0,
  146.0, 148.0, 150.0, 152.0, 154.0, 156.0, 158.0, 160.0, 162.0, 164.0, 166.0,
  168.0, 170.0, 170.5, 171.0, 171.5, 172.0, 172.5, 173.0, 173.5, 174.0, 174.5,
  175.0, 175.5, 176.0, 176.5, 177.0, 177.5, 178.0, 178.5, 179.0, 179.5, 180.0,
};
double	*mie_angl_rad			= NULL;
double	*mie_angl_sin			= NULL;
double	*mie_angl_cos			= NULL;
double	*mie_angl_dif			= NULL;
double	*mie_phas			= NULL;
double	*mie_tropo_phas			= NULL;
double	*mie_strat_phas			= NULL;
double	*mie_meteo_phas			= NULL;
double	*mie_refr_com[MIE_MAXCOMP];				// Refractive index (real)
double	*mie_refi_com[MIE_MAXCOMP];				// Refractive index (imaginary)
double	*mie_aext_com[MIE_MAXCOMP];				// Extinction coefficient
double	*mie_asca_com[MIE_MAXCOMP];				// Scattering coefficient
double	*mie_asym_com[MIE_MAXCOMP];				// Asymmetry parameter
double	*mie_phs1_com[MIE_MAXCOMP];				// Phase function 1
double	*mie_phs2_com[MIE_MAXCOMP];				// Phase function 2
double	*mie_phas_com[MIE_MAXCOMP];				// Averaged phase function
double	*mie_aext_tab[MIE_MAXCOMP];				// Extinction coefficient
double	*mie_asca_tab[MIE_MAXCOMP];				// Scattering coefficient
double	*mie_phs1_tab[MIE_MAXCOMP];				// Phase function 1
double	*mie_phs2_tab[MIE_MAXCOMP];				// Phase function 2
double	*mie_lrad_lext_x[MIE_MAXCOMP];				// Log10(R)
double	*mie_lrad_lext_y[MIE_MAXCOMP];				// Log10(Qext*PI*r*r)
double	*mie_lext_lrad_x[MIE_MAXCOMP];				// Log10(Qext*PI*r*r)
double	*mie_lext_lrad_y[MIE_MAXCOMP];				// Log10(R)
double	*mie_lmod_leff_x[MIE_MAXCOMP];				// Log10(R)
double	*mie_lmod_leff_y[MIE_MAXCOMP];				// Log10(R)
double	*mie_leff_lmod_x[MIE_MAXCOMP];				// Log10(R)
double	*mie_leff_lmod_y[MIE_MAXCOMP];				// Log10(R)
double	mie_mixr_com[MIE_MAXCOMP];				// Number mixing ratio
double	mie_wcom_com[MIE_MAXCOMP];				// Weight
double	mie_lmod_com[MIE_MAXCOMP];				// Log10(R in um)
double	mie_lsgm_com[MIE_MAXCOMP];				// Sigma Log10(R in um)
double	*mie_lrad_min[MIE_MAXCOMP];				// Log10(R) min
double	*mie_lrad_max[MIE_MAXCOMP];				// Log10(R) max
double	*mie_lrad_stp[MIE_MAXCOMP];				// Log10(R) step
double	mie_lext_min[MIE_MAXCOMP];				// Log10(aext) min
double	mie_lext_max[MIE_MAXCOMP];				// Log10(aext) max
double	mie_lext_stp[MIE_MAXCOMP];				// Log10(aext) step
double	mie_lmod_min[MIE_MAXCOMP];				// Log10(R) min
double	mie_lmod_max[MIE_MAXCOMP];				// Log10(R) max
double	mie_lmod_stp[MIE_MAXCOMP];				// Log10(R) step
double	mie_leff_min[MIE_MAXCOMP];				// Log10(R) min
double	mie_leff_max[MIE_MAXCOMP];				// Log10(R) max
double	mie_leff_stp[MIE_MAXCOMP];				// Log10(R) step
double	mie_rmin			= MIE_RMIN;		// R min in um
double	mie_rmax			= MIE_RMAX;		// R max in um
double	mie_wsgm			= MIE_WSGM;		// Log10(R) width in sigma
char	mie_inp1[MAXLINE]		= MIE_INP1;		// Inp. file 1
char	mie_inp2[MAXLINE]		= MIE_INP2;		// Inp. file 2
char	mie_inp3[MAXLINE]		= MIE_INP3;		// Inp. file 3
char	mie_inp4[MAXLINE]		= MIE_INP4;		// Inp. file 4
char	mie_inp5[MAXLINE]		= MIE_INP5;		// Inp. file 5
char	mie_out1[MAXLINE]		= MIE_OUT1;		// Out. file 1
char	mie_out2[MAXLINE]		= MIE_OUT2;		// Out. file 2
// Parameters for T-matrix calculation
int	tmx_ndis			= TMX_NDIS;		// Distribution#
int	tmx_nang			= TMX_NANG;		// #Angles
int	tmx_ndgs			= TMX_NDGS;		// Control #divisions
int	tmx_kmax			= TMX_KMAX;		// #Quad points
int	tmx_shap[MIE_MAXCOMP];					// Shape#
double	tmx_reps[MIE_MAXCOMP];					// Aspect ratio or deformation parameter
double	tmx_delt			= TMX_DELT;		// Accuracy of the computations
// Parameters for Test
int	tst_n_comp			= NODATA;		// #Components
int	tst_mode			= 0;			// Test mode
double	tst_wcom_com[MIE_MAXCOMP];				// Weight
double	tst_lmod_com[MIE_MAXCOMP];				// Mode radius in um
double	tst_lsgm_com[MIE_MAXCOMP];				// Sigma Log10(R in um)
double	*tst_refr_com[MIE_MAXCOMP];				// Refractive index (real)
double	*tst_refi_com[MIE_MAXCOMP];				// Refractive index (imaginary)
char	tst_fsim[MAXLINE]		= NONAME;		// Simulation file
// Parameters for Upper Atmosphere
double	upp_wlen[UPP_N_WLEN] = {200.0, 300.0, 550.0, 694.3, 1060.0, 1536.0};
double	upp_angl[UPP_N_ANGL] =
{
    0.0,   2.0,   4.0,   6.0,   8.0,  10.0,  12.0,  16.0,
   20.0,  24.0,  28.0,  32.0,  36.0,  40.0,  50.0,  60.0,
   70.0,  80.0,  90.0, 100.0, 110.0, 120.0, 125.0, 130.0,
  135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0,
  175.0, 180.0,
};
double	upp_rhum[UPP_N_RHUM] = {0.0, 70.0, 80.0, 99.0};
double	upp_tropo_phas[UPP_N_RHUM][UPP_N_PHAS] =
{
  {// RH=0%
    // PHSFNC_29
    2.72707, 2.52017, 2.15725, 1.81463, 1.52555, 1.28665, 1.08932, 0.78963,
    0.58023, 0.43161, 0.32473, 0.24695, 0.18977, 0.14727, 0.08150, 0.04783,
    0.02980, 0.01981, 0.01414, 0.01091, 0.00911, 0.00819, 0.00796, 0.00784,
    0.00781, 0.00787, 0.00798, 0.00812, 0.00825, 0.00831, 0.00830, 0.00829,
    0.00847, 0.00868,
    // PHSFNC_16
    1.52473, 1.44645, 1.30800, 1.16736, 1.03716, 0.92088, 0.81802, 0.64725,
    0.51405, 0.40964, 0.32742, 0.26255, 0.21129, 0.17075, 0.10236, 0.06361,
    0.04130, 0.02819, 0.02041, 0.01583, 0.01330, 0.01222, 0.01213, 0.01235,
    0.01288, 0.01373, 0.01491, 0.01635, 0.01784, 0.01906, 0.02003, 0.02225,
    0.02715, 0.03082,
    // PHSFNC_32
    0.79503, 0.78714, 0.76605, 0.73601, 0.70017, 0.66092, 0.61994, 0.53736,
    0.45896, 0.38796, 0.32563, 0.27208, 0.22676, 0.18879, 0.11979, 0.07727,
    0.05125, 0.03541, 0.02586, 0.02025, 0.01717, 0.01576, 0.01552, 0.01551,
    0.01571, 0.01607, 0.01656, 0.01710, 0.01764, 0.01814, 0.01872, 0.01965,
    0.02093, 0.02164,
    // PHSFNC_32
    0.79503, 0.78714, 0.76605, 0.73601, 0.70017, 0.66092, 0.61994, 0.53736,
    0.45896, 0.38796, 0.32563, 0.27208, 0.22676, 0.18879, 0.11979, 0.07727,
    0.05125, 0.03541, 0.02586, 0.02025, 0.01717, 0.01576, 0.01552, 0.01551,
    0.01571, 0.01607, 0.01656, 0.01710, 0.01764, 0.01814, 0.01872, 0.01965,
    0.02093, 0.02164,
    // PHSFNC_36
    0.72999, 0.72086, 0.69632, 0.66227, 0.62369, 0.58366, 0.54390, 0.46850,
    0.40091, 0.34178, 0.29081, 0.24724, 0.21025, 0.17893, 0.12047, 0.08253,
    0.05810, 0.04253, 0.03278, 0.02688, 0.02352, 0.02183, 0.02142, 0.02122,
    0.02118, 0.02127, 0.02146, 0.02171, 0.02202, 0.02235, 0.02273, 0.02314,
    0.02351, 0.02367,
    // PHSFNC_36
    0.72999, 0.72086, 0.69632, 0.66227, 0.62369, 0.58366, 0.54390, 0.46850,
    0.40091, 0.34178, 0.29081, 0.24724, 0.21025, 0.17893, 0.12047, 0.08253,
    0.05810, 0.04253, 0.03278, 0.02688, 0.02352, 0.02183, 0.02142, 0.02122,
    0.02118, 0.02127, 0.02146, 0.02171, 0.02202, 0.02235, 0.02273, 0.02314,
    0.02351, 0.02367,
  },
  {// RH=70%
    // PHSFNC_29
    2.72707, 2.52017, 2.15725, 1.81463, 1.52555, 1.28665, 1.08932, 0.78963,
    0.58023, 0.43161, 0.32473, 0.24695, 0.18977, 0.14727, 0.08150, 0.04783,
    0.02980, 0.01981, 0.01414, 0.01091, 0.00911, 0.00819, 0.00796, 0.00784,
    0.00781, 0.00787, 0.00798, 0.00812, 0.00825, 0.00831, 0.00830, 0.00829,
    0.00847, 0.00868,
    // PHSFNC_16
    1.52473, 1.44645, 1.30800, 1.16736, 1.03716, 0.92088, 0.81802, 0.64725,
    0.51405, 0.40964, 0.32742, 0.26255, 0.21129, 0.17075, 0.10236, 0.06361,
    0.04130, 0.02819, 0.02041, 0.01583, 0.01330, 0.01222, 0.01213, 0.01235,
    0.01288, 0.01373, 0.01491, 0.01635, 0.01784, 0.01906, 0.02003, 0.02225,
    0.02715, 0.03082,
    // PHSFNC_32
    0.79503, 0.78714, 0.76605, 0.73601, 0.70017, 0.66092, 0.61994, 0.53736,
    0.45896, 0.38796, 0.32563, 0.27208, 0.22676, 0.18879, 0.11979, 0.07727,
    0.05125, 0.03541, 0.02586, 0.02025, 0.01717, 0.01576, 0.01552, 0.01551,
    0.01571, 0.01607, 0.01656, 0.01710, 0.01764, 0.01814, 0.01872, 0.01965,
    0.02093, 0.02164,
    // PHSFNC_32
    0.79503, 0.78714, 0.76605, 0.73601, 0.70017, 0.66092, 0.61994, 0.53736,
    0.45896, 0.38796, 0.32563, 0.27208, 0.22676, 0.18879, 0.11979, 0.07727,
    0.05125, 0.03541, 0.02586, 0.02025, 0.01717, 0.01576, 0.01552, 0.01551,
    0.01571, 0.01607, 0.01656, 0.01710, 0.01764, 0.01814, 0.01872, 0.01965,
    0.02093, 0.02164,
    // PHSFNC_32
    0.79503, 0.78714, 0.76605, 0.73601, 0.70017, 0.66092, 0.61994, 0.53736,
    0.45896, 0.38796, 0.32563, 0.27208, 0.22676, 0.18879, 0.11979, 0.07727,
    0.05125, 0.03541, 0.02586, 0.02025, 0.01717, 0.01576, 0.01552, 0.01551,
    0.01571, 0.01607, 0.01656, 0.01710, 0.01764, 0.01814, 0.01872, 0.01965,
    0.02093, 0.02164,
    // PHSFNC_36
    0.72999, 0.72086, 0.69632, 0.66227, 0.62369, 0.58366, 0.54390, 0.46850,
    0.40091, 0.34178, 0.29081, 0.24724, 0.21025, 0.17893, 0.12047, 0.08253,
    0.05810, 0.04253, 0.03278, 0.02688, 0.02352, 0.02183, 0.02142, 0.02122,
    0.02118, 0.02127, 0.02146, 0.02171, 0.02202, 0.02235, 0.02273, 0.02314,
    0.02351, 0.02367,
  },
  {// RH=80%
    // PHSFNC_29
    2.72707, 2.52017, 2.15725, 1.81463, 1.52555, 1.28665, 1.08932, 0.78963,
    0.58023, 0.43161, 0.32473, 0.24695, 0.18977, 0.14727, 0.08150, 0.04783,
    0.02980, 0.01981, 0.01414, 0.01091, 0.00911, 0.00819, 0.00796, 0.00784,
    0.00781, 0.00787, 0.00798, 0.00812, 0.00825, 0.00831, 0.00830, 0.00829,
    0.00847, 0.00868,
    // PHSFNC_28
    1.99262, 1.84586, 1.62814, 1.42586, 1.24733, 1.09157, 0.95572, 0.73435,
    0.56567, 0.43699, 0.33874, 0.26371, 0.20635, 0.16235, 0.09181, 0.05440,
    0.03397, 0.02253, 0.01601, 0.01230, 0.01031, 0.00946, 0.00938, 0.00951,
    0.00985, 0.01043, 0.01123, 0.01222, 0.01328, 0.01405, 0.01405, 0.01340,
    0.01391, 0.01521,
    // PHSFNC_37
    1.15457, 1.13481, 1.08637, 1.02208, 0.95035, 0.87581, 0.80158, 0.66139,
    0.53812, 0.43402, 0.34834, 0.27905, 0.22360, 0.17950, 0.10538, 0.06398,
    0.04056, 0.02712, 0.01935, 0.01490, 0.01250, 0.01143, 0.01126, 0.01130,
    0.01153, 0.01191, 0.01242, 0.01301, 0.01357, 0.01399, 0.01422, 0.01450,
    0.01526, 0.01585,
    // PHSFNC_37
    1.15457, 1.13481, 1.08637, 1.02208, 0.95035, 0.87581, 0.80158, 0.66139,
    0.53812, 0.43402, 0.34834, 0.27905, 0.22360, 0.17950, 0.10538, 0.06398,
    0.04056, 0.02712, 0.01935, 0.01490, 0.01250, 0.01143, 0.01126, 0.01130,
    0.01153, 0.01191, 0.01242, 0.01301, 0.01357, 0.01399, 0.01422, 0.01450,
    0.01526, 0.01585,
    // PHSFNC_32
    0.79503, 0.78714, 0.76605, 0.73601, 0.70017, 0.66092, 0.61994, 0.53736,
    0.45896, 0.38796, 0.32563, 0.27208, 0.22676, 0.18879, 0.11979, 0.07727,
    0.05125, 0.03541, 0.02586, 0.02025, 0.01717, 0.01576, 0.01552, 0.01551,
    0.01571, 0.01607, 0.01656, 0.01710, 0.01764, 0.01814, 0.01872, 0.01965,
    0.02093, 0.02164,
    // PHSFNC_36
    0.72999, 0.72086, 0.69632, 0.66227, 0.62369, 0.58366, 0.54390, 0.46850,
    0.40091, 0.34178, 0.29081, 0.24724, 0.21025, 0.17893, 0.12047, 0.08253,
    0.05810, 0.04253, 0.03278, 0.02688, 0.02352, 0.02183, 0.02142, 0.02122,
    0.02118, 0.02127, 0.02146, 0.02171, 0.02202, 0.02235, 0.02273, 0.02314,
    0.02351, 0.02367,
  },
  {// RH=99%
    // PHSFNC_15
    7.98100, 5.58389, 3.81067, 2.75011, 2.05900, 1.58178, 1.24011, 0.80424,
    0.55144, 0.39021, 0.28203, 0.20701, 0.15369, 0.11570, 0.05985, 0.03329,
    0.01999, 0.01295, 0.00920, 0.00726, 0.00645, 0.00646, 0.00683, 0.00743,
    0.00833, 0.00942, 0.01038, 0.01112, 0.01187, 0.01290, 0.01451, 0.01415,
    0.01176, 0.01646,
    // PHSFNC_26
    4.10720, 3.38127, 2.61113, 2.06147, 1.66320, 1.36400, 1.13140, 0.79725,
    0.57362, 0.41947, 0.31103, 0.23368, 0.17766, 0.13659, 0.07413, 0.04288,
    0.02647, 0.01753, 0.01255, 0.00978, 0.00837, 0.00791, 0.00800, 0.00831,
    0.00888, 0.00974, 0.01087, 0.01219, 0.01347, 0.01448, 0.01468, 0.01320,
    0.01237, 0.01467,
    // PHSFNC_28
    1.99262, 1.84586, 1.62814, 1.42586, 1.24733, 1.09157, 0.95572, 0.73435,
    0.56567, 0.43699, 0.33874, 0.26371, 0.20635, 0.16235, 0.09181, 0.05440,
    0.03397, 0.02253, 0.01601, 0.01230, 0.01031, 0.00946, 0.00938, 0.00951,
    0.00985, 0.01043, 0.01123, 0.01222, 0.01328, 0.01405, 0.01405, 0.01340,
    0.01391, 0.01521,
    // PHSFNC_28
    1.99262, 1.84586, 1.62814, 1.42586, 1.24733, 1.09157, 0.95572, 0.73435,
    0.56567, 0.43699, 0.33874, 0.26371, 0.20635, 0.16235, 0.09181, 0.05440,
    0.03397, 0.02253, 0.01601, 0.01230, 0.01031, 0.00946, 0.00938, 0.00951,
    0.00985, 0.01043, 0.01123, 0.01222, 0.01328, 0.01405, 0.01405, 0.01340,
    0.01391, 0.01521,
    // PHSFNC_37
    1.15457, 1.13481, 1.08637, 1.02208, 0.95035, 0.87581, 0.80158, 0.66139,
    0.53812, 0.43402, 0.34834, 0.27905, 0.22360, 0.17950, 0.10538, 0.06398,
    0.04056, 0.02712, 0.01935, 0.01490, 0.01250, 0.01143, 0.01126, 0.01130,
    0.01153, 0.01191, 0.01242, 0.01301, 0.01357, 0.01399, 0.01422, 0.01450,
    0.01526, 0.01585,
    // PHSFNC_37
    1.15457, 1.13481, 1.08637, 1.02208, 0.95035, 0.87581, 0.80158, 0.66139,
    0.53812, 0.43402, 0.34834, 0.27905, 0.22360, 0.17950, 0.10538, 0.06398,
    0.04056, 0.02712, 0.01935, 0.01490, 0.01250, 0.01143, 0.01126, 0.01130,
    0.01153, 0.01191, 0.01242, 0.01301, 0.01357, 0.01399, 0.01422, 0.01450,
    0.01526, 0.01585,
  },
};
double	upp_strat_phas[UPP_N_PHAS] =
{
  // PHSFNC_20
  3.35750, 2.77033, 2.20467, 1.77208, 1.43950, 1.18167, 0.98083, 0.69731,
  0.51368, 0.38864, 0.29989, 0.23472, 0.18552, 0.14785, 0.08680, 0.05339,
  0.03443, 0.02340, 0.01689, 0.01308, 0.01097, 0.01015, 0.01022, 0.01062,
  0.01145, 0.01289, 0.01512, 0.01829, 0.02228, 0.02580, 0.02615, 0.02523,
  0.02990, 0.03617,
  // PHSFNC_20
  3.35750, 2.77033, 2.20467, 1.77208, 1.43950, 1.18167, 0.98083, 0.69731,
  0.51368, 0.38864, 0.29989, 0.23472, 0.18552, 0.14785, 0.08680, 0.05339,
  0.03443, 0.02340, 0.01689, 0.01308, 0.01097, 0.01015, 0.01022, 0.01062,
  0.01145, 0.01289, 0.01512, 0.01829, 0.02228, 0.02580, 0.02615, 0.02523,
  0.02990, 0.03617,
  // PHSFNC_37
  1.15457, 1.13481, 1.08637, 1.02208, 0.95035, 0.87581, 0.80158, 0.66139,
  0.53812, 0.43402, 0.34834, 0.27905, 0.22360, 0.17950, 0.10538, 0.06398,
  0.04056, 0.02712, 0.01935, 0.01490, 0.01250, 0.01143, 0.01126, 0.01130,
  0.01153, 0.01191, 0.01242, 0.01301, 0.01357, 0.01399, 0.01422, 0.01450,
  0.01526, 0.01585,
  // PHSFNC_37
  1.15457, 1.13481, 1.08637, 1.02208, 0.95035, 0.87581, 0.80158, 0.66139,
  0.53812, 0.43402, 0.34834, 0.27905, 0.22360, 0.17950, 0.10538, 0.06398,
  0.04056, 0.02712, 0.01935, 0.01490, 0.01250, 0.01143, 0.01126, 0.01130,
  0.01153, 0.01191, 0.01242, 0.01301, 0.01357, 0.01399, 0.01422, 0.01450,
  0.01526, 0.01585,
  // PHSFNC_24
  0.51072, 0.50922, 0.50476, 0.49744, 0.48752, 0.47514, 0.46072, 0.42698,
  0.38890, 0.34888, 0.30902, 0.27072, 0.23504, 0.20254, 0.13652, 0.09070,
  0.06070, 0.04186, 0.03042, 0.02372, 0.02000, 0.01816, 0.01773, 0.01753,
  0.01753, 0.01767, 0.01792, 0.01826, 0.01866, 0.01907, 0.01948, 0.01984,
  0.02009, 0.02019,
  // PHSFNC_23
  0.45555, 0.45332, 0.44686, 0.43683, 0.42405, 0.40929, 0.39319, 0.35887,
  0.32388, 0.28985, 0.25771, 0.22804, 0.20101, 0.17671, 0.12722, 0.09162,
  0.06691, 0.05037, 0.03979, 0.03345, 0.03009, 0.02875, 0.02860, 0.02871,
  0.02902, 0.02946, 0.03000, 0.03059, 0.03120, 0.03179, 0.03232, 0.03276,
  0.03305, 0.03315,
};
double	upp_meteo_phas[UPP_N_PHAS] =
{
  // PHSFNC_59
  56.57000,10.86983,4.30450, 2.32817, 1.48817, 1.05608, 0.80615, 0.53590,
  0.39175, 0.30123, 0.23833, 0.19167, 0.15528, 0.12653, 0.07717, 0.04870,
  0.03126, 0.02128, 0.01493, 0.01091, 0.00856, 0.00765, 0.00767, 0.00787,
  0.00840, 0.00973, 0.01199, 0.01623, 0.02446, 0.03802, 0.04421, 0.04801,
  0.05991, 0.06067,
  // PHSFNC_59
  56.57000,10.86983,4.30450, 2.32817, 1.48817, 1.05608, 0.80615, 0.53590,
  0.39175, 0.30123, 0.23833, 0.19167, 0.15528, 0.12653, 0.07717, 0.04870,
  0.03126, 0.02128, 0.01493, 0.01091, 0.00856, 0.00765, 0.00767, 0.00787,
  0.00840, 0.00973, 0.01199, 0.01623, 0.02446, 0.03802, 0.04421, 0.04801,
  0.05991, 0.06067,
  // PHSFNC_11
  15.36571,6.08614, 3.15543, 1.98500, 1.39671, 1.05749, 0.84297, 0.58904,
  0.44120, 0.34216, 0.27140, 0.21837, 0.17669, 0.14371, 0.08723, 0.05498,
  0.03568, 0.02444, 0.01753, 0.01339, 0.01087, 0.00986, 0.00983, 0.01019,
  0.01094, 0.01243, 0.01473, 0.01871, 0.02423, 0.03204, 0.03592, 0.03792,
  0.04492, 0.05004,
  // PHSFNC_11
  15.36571,6.08614, 3.15543, 1.98500, 1.39671, 1.05749, 0.84297, 0.58904,
  0.44120, 0.34216, 0.27140, 0.21837, 0.17669, 0.14371, 0.08723, 0.05498,
  0.03568, 0.02444, 0.01753, 0.01339, 0.01087, 0.00986, 0.00983, 0.01019,
  0.01094, 0.01243, 0.01473, 0.01871, 0.02423, 0.03204, 0.03592, 0.03792,
  0.04492, 0.05004,
  // PHSFNC_20
  3.35750, 2.77033, 2.20467, 1.77208, 1.43950, 1.18167, 0.98083, 0.69731,
  0.51368, 0.38864, 0.29989, 0.23472, 0.18552, 0.14785, 0.08680, 0.05339,
  0.03443, 0.02340, 0.01689, 0.01308, 0.01097, 0.01015, 0.01022, 0.01062,
  0.01145, 0.01289, 0.01512, 0.01829, 0.02228, 0.02580, 0.02615, 0.02523,
  0.02990, 0.03617,
  // PHSFNC_20
  3.35750, 2.77033, 2.20467, 1.77208, 1.43950, 1.18167, 0.98083, 0.69731,
  0.51368, 0.38864, 0.29989, 0.23472, 0.18552, 0.14785, 0.08680, 0.05339,
  0.03443, 0.02340, 0.01689, 0.01308, 0.01097, 0.01015, 0.01022, 0.01062,
  0.01145, 0.01289, 0.01512, 0.01829, 0.02228, 0.02580, 0.02615, 0.02523,
  0.02990, 0.03617,
};
// Parameters for control
int	cnt_dsr[DSR_MAXRANG]		= {1,1,1,0};		// DSR     flag
int	cnt_aur[AUR_MAXRANG]		= {0,0};		// AUR     flag
int	cnt_aur_auto[AUR_MAXRANG]	= {0,0};		// AUR     flag
int	cnt_ssr[SSR_MAXRANG]		= {0,1,0,0};		// SSR     flag
int	cnt_enable_old_aur		= 0;			// Enable old AUR
int	cnt_cflg			= 0;			// Clear   flag
int	cnt_sflg			= 0;			// Single scattering flag
int	cnt_tmx[MIE_MAXCOMP]		= {0,0,0,0,0,0,0,0,0,0};// T-matrix flag
int	cnt_db				= 0;			// Debug   mode
int	cnt_vb				= 0;			// Verbose mode
int	cnt_hp				= 0;			// Help    mode
int	cnt_n_cmnt			= NODATA;
int	cnt_n_own_format		= NODATA;
int	cnt_optn			= 1;
char	cnt_conf[MAXLINE]		= CNT_CONF;		// Configuration file
char	cnt_cmnt[CNT_MAXCMNT][MAXLINE];				// Comments
char	cnt_own_format[CNT_MAXFORM][MAXLINE];			// Formats
char	*cnt_opts[CNT_MAXNOPT];

int Init(void);
void Finish(void);
int GetOpt(int argn,char **args);
int GetOwnOpt(int argn,char **args);
int GetCommonOpt(int argn,char **args);
int Usage(void);
int CommonUsage(char *name,int c);
int DSR_FCN_MODTRAN_B(int rng,double *par,double *val,int mode);
int AUR_FCN_MODTRAN_S(int rng,double *par,double *val,int mode);
int AUR_FCN_MODTRAN_B(int rng,double *par,double *val,int mode);
int SSR_FCN_MODTRAN_S(int rng,double *par,double *val,int mode);
int SSR_FCN_MODTRAN_B(int rng,double *par,double *val,int mode);
int DSR_FCN_CUSTOM_B(int rng,double *par,double *val,int mode);
int AUR_FCN_CUSTOM_S(int rng,double *par,double *val,int mode);
int AUR_FCN_CUSTOM_B(int rng,double *par,double *val,int mode);
int SSR_FCN_CUSTOM_S(int rng,double *par,double *val,int mode);
int SSR_FCN_CUSTOM_B(int rng,double *par,double *val,int mode);
int DSR_SSR_MODTRAN_S(int rng,double *par,int mode);
int DSR_SSR_MODTRAN_B(int rng,double *par,int mode);
int DSR_SSR_CUSTOM_S(int rng,double *par,int mode);
int DSR_SSR_CUSTOM_B(int rng,double *par,int mode);
int DSR_FCN_MODTRAN_B0(double *par,double *val);
int DSR_FCN_MODTRAN_B1(double *par,double *val);
int DSR_FCN_MODTRAN_B2(double *par,double *val);
int DSR_FCN_MODTRAN_B3(double *par,double *val);
int AUR_FCN_MODTRAN_S0(double *par,double *val);
int AUR_FCN_MODTRAN_S1(double *par,double *val);
int AUR_FCN_MODTRAN_B0(double *par,double *val);
int AUR_FCN_MODTRAN_B1(double *par,double *val);
int SSR_FCN_MODTRAN_S0(double *par,double *val);
int SSR_FCN_MODTRAN_S1(double *par,double *val);
int SSR_FCN_MODTRAN_S2(double *par,double *val);
int SSR_FCN_MODTRAN_S3(double *par,double *val);
int SSR_FCN_MODTRAN_B0(double *par,double *val);
int SSR_FCN_MODTRAN_B1(double *par,double *val);
int SSR_FCN_MODTRAN_B2(double *par,double *val);
int SSR_FCN_MODTRAN_B3(double *par,double *val);
int DSR_FCN_CUSTOM_B0(double *par,double *val);
int DSR_FCN_CUSTOM_B1(double *par,double *val);
int DSR_FCN_CUSTOM_B2(double *par,double *val);
int DSR_FCN_CUSTOM_B3(double *par,double *val);
int AUR_FCN_CUSTOM_S0(double *par,double *val);
int AUR_FCN_CUSTOM_S1(double *par,double *val);
int AUR_FCN_CUSTOM_B0(double *par,double *val);
int AUR_FCN_CUSTOM_B1(double *par,double *val);
int SSR_FCN_CUSTOM_S0(double *par,double *val);
int SSR_FCN_CUSTOM_S1(double *par,double *val);
int SSR_FCN_CUSTOM_S2(double *par,double *val);
int SSR_FCN_CUSTOM_S3(double *par,double *val);
int SSR_FCN_CUSTOM_B0(double *par,double *val);
int SSR_FCN_CUSTOM_B1(double *par,double *val);
int SSR_FCN_CUSTOM_B2(double *par,double *val);
int SSR_FCN_CUSTOM_B3(double *par,double *val);
int DSR_SSR_MODTRAN_S0(double *par);
int DSR_SSR_MODTRAN_S1(double *par);
int DSR_SSR_MODTRAN_S2(double *par);
int DSR_SSR_MODTRAN_S3(double *par);
int DSR_SSR_MODTRAN_B0(double *par);
int DSR_SSR_MODTRAN_B1(double *par);
int DSR_SSR_MODTRAN_B2(double *par);
int DSR_SSR_MODTRAN_B3(double *par);
int DSR_SSR_CUSTOM_S0(double *par);
int DSR_SSR_CUSTOM_S1(double *par);
int DSR_SSR_CUSTOM_S2(double *par);
int DSR_SSR_CUSTOM_S3(double *par);
int DSR_SSR_CUSTOM_B0(double *par);
int DSR_SSR_CUSTOM_B1(double *par);
int DSR_SSR_CUSTOM_B2(double *par);
int DSR_SSR_CUSTOM_B3(double *par);
int CRS_FCN_DSR_AUR(double *par,double *val);
int CRS_FCN_DSR_SSR(double *par,double *val);
int CRS_FCN_AUR_SSR(double *par,double *val);
int CRS_FCN_DSR_AUR_SSR(double *par,double *val);
int CheckParam(const double *par,int *fpar);
int PrintValue(const double *par,const double *val,int npar,...);
int MieInit(void);
int MieTable(int n1,int n2);
int MieReff(int n1,int n2,double lsgm);
int MieCalc(const double *par,const int *fpar);
int MixComp(const double *par);
int MieSmoothing(double *x,double *y,double xstp,double xmgn,double xwid,int nstp);
double mod2eff(int n,double lmod);
double eff2mod(int n,double leff);
int MiePrintout(void);
int ResFactor(int inst,int nw,int sbmp,double wlen,double amas,double *val);
int CalcFactor(double *par,int mode,int flag);
int Background(const double *par,int mode);
int RingCenter(const double *par,int mode);
int Uniformity(const double *par,int mode);
int PrfInitMolecule(int imod);
int PrfInitAerosol(int ihaz,int isea,int ivul,double vis);
int PrfGround(void);
int PrfScaleMolecule(void);
int PrfScaleAerosol(double t);
int PrfCalcAOD(void);
int upper_pfunc(double rh);
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
int Interp1D(const double *x1,const double *z1,int nx1,
             const double *x2,      double *z2,int nx2,int f);
int Interp2D(const double *x1,const double *y1,const double *z1,int nx1,int ny1,
             const double *x2,const double *y2,      double *z2,int nx2,int ny2,int f);
int WriteTape5(int rng,double *par,int mode);
int WriteCard1(FILE *fp,int mode);
int WriteCard2C(FILE *fp);
int WriteCard2D(FILE *fp);
int WriteCard3C(FILE *fp);
int RunModtran(void);
int ReadPlt_S(const double *par,double *y);
int ReadPlt_B(const double *par,int n,const double *x,double *y);
int ReadDSR_B(int rng,const double *par);
int ReadAUR_S(int rng,const double *par,double y[AUR_MAXDATA][AUR_MAXWLEN]);
int ReadAUR_B(int rng,const double *par,double y[AUR_MAXDATA][AUR_MAXWLEN]);
int ReadSSR_S(int rng,const double *par,double y[SSR_MAXDATA][SSR_MAXWLEN]);
int ReadSSR_B(int rng,const double *par,double y[SSR_MAXDATA][SSR_MAXWLEN]);
int Convolute(double xmin,double xmax,double xstp,double xsgm,double wsgm,
              int ninp,const double *xinp,const double *yinp,
              int nout,double *xout,double *yout);
double Gauss(double g,double f0,double f);
int Sampling(int ninp,const double *xinp,const double *yinp,
             int nout,const double *xout,double *yout,double yuni);
int SamplingE(int ninp,const double *xinp,const double *yinp,
              int nout,const double *xout,double *yout,double yuni);
int Printout(void);
int CommonPrintout(FILE *fp,int flag_dsr,int flag_aur,int flag_ssr);
int CommonInit(void);
int OwnInit(void);
int PreConfig(void);
int ReadConfig(void);
int PostConfig(void);
int ReadOwnConfig(char *s);
int ReadSimul(char *s,int flag_dsr,int flag_aur,int flag_ssr,int flag_dat,int flag_sim);
int ReadInput(void);
int ReadCRC(void);
int ReadFactor(void);
int ReadXA(char *s,int size,int cx,double ux,int (*func)(double*),double *x);
int ReadXYA(char *s,int size,int cx,int cy,double ux,double uy,int (*func)(double*),double *x,double *y);
int ReadXYZA(char *s,int size,int cx,int cy,int cz,double ux,double uy,double uz,int (*func)(double*),
             double *x,double *y,double *z);
int ReadXP(char *s,int size,int cx,double ux,int (*func)(double*),double **x);
int ReadXYP(char *s,int size,int cx,int cy,double ux,double uy,int (*func)(double*),double **x,double **y);
int ReadXYZP(char *s,int size,int cx,int cy,int cz,double ux,double uy,double uz,int (*func)(double*),
             double **x,double **y,double **z);
int ReadDC(char *s,int np,int size,const int *c,const double *u,int (*func)(double*),double **d);
int ReadDA(char *s,int np,int size,const int *c,const double *u,double **d);
int ReadNA(char *s,int size,int c,int *d);
int WaveSelect(double *w);
int AnglSelect(double *a);
int ReadComp(char *s,int size,int cr,int ci,int imin,int imax,double u,
             double *wcom,double *lmod,double *lsgm,double **refr,double **refi);

extern void miev0_();
extern void tmatrix_();
extern float aerprf_(int*,float*,int*,int*,int*);

void Finish(void)
{
  int n;
  int n_comp;

  n_comp = MAX(OPT_N_COMP,MAX(tst_n_comp,crc_n_comp));
  for(n=0; n<n_comp; n++)
  {
    free(mie_refr_com[n]);
    free(mie_refi_com[n]);
    free(mie_aext_com[n]);
    free(mie_asca_com[n]);
    free(mie_asym_com[n]);
    free(mie_phs1_com[n]);
    free(mie_phs2_com[n]);
    free(mie_phas_com[n]);
    free(mie_aext_tab[n]);
    free(mie_asca_tab[n]);
    free(mie_phs1_tab[n]);
    free(mie_phs2_tab[n]);
    free(mie_lrad_min[n]);
    free(mie_lrad_max[n]);
    free(mie_lrad_stp[n]);
    free(mie_lrad_lext_x[n]);
    free(mie_lrad_lext_y[n]);
    free(mie_lext_lrad_x[n]);
    free(mie_lext_lrad_y[n]);
    free(mie_lmod_leff_x[n]);
    free(mie_lmod_leff_y[n]);
    free(mie_leff_lmod_x[n]);
    free(mie_leff_lmod_y[n]);
  }
  for(n=0; n<tst_n_comp; n++)
  {
    free(tst_refr_com[n]);
    free(tst_refi_com[n]);
  }
  for(n=0; n<crc_n_comp; n++)
  {
    free(crc_refr_com[n]);
    free(crc_refi_com[n]);
  }
  for(n=0; n<opt_n_comp; n++)
  {
    free(opt_refr_com[n]);
    free(opt_refi_com[n]);
  }
  free(mie_wlen_um);
  free(mie_aext);
  free(mie_asca);
  free(mie_asym);
  free(mie_angl_rad);
  free(mie_angl_sin);
  free(mie_angl_cos);
  free(mie_angl_dif);
  free(mie_phas);
  free(mie_tropo_phas);
  free(mie_strat_phas);
  free(mie_meteo_phas);
  free(sim_phas);
  free(sim_tropo_phas);
  free(sim_strat_phas);
  free(sim_meteo_phas);

  if(cnt_cflg > 0)
  {
    system("rm -f modtran.7sc");
    system("rm -f modtran.7sr");
    system("rm -f modtran.mc");
    system("rm -f modtran.tp6");
    system("rm -f modtran.tp7");
    system("rm -f modtran.tp8");
  }
  if(cnt_cflg > 1)
  {
    system("rm -f DATA");
    system("rm -f modroot.in");
    system("rm -f modtran.tp5");
  }
  if(cnt_cflg > 2)
  {
    system("rm -f modtran.plt");
  }

  return;
}

// Direct Solar Radiation Fitting (MODTRAN aerosol model, wavelength band)
int DSR_FCN_MODTRAN_B(int rng,double *par,double *val,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j,k;
  int err;
  int fpar[OPT_N_FPAR];
  int nval;
  double p[2];
  double r,v;
  double d1,d2,dr;
  double s1,s2,sr;
  double w1,w2;

  if(CheckParam(par,fpar) < 0)
  {
    *val = HUGE;
    return -1;
  }

  k = dsr_n_wlen[rng]-1;
  w1 = dsr_wlen[rng][0];
  w2 = dsr_wlen[rng][k];
  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    do
    {
      // individual setup
      switch(sim_mode)
      {
        case SIM_MODTRAN_VTAU:
          if(fpar[OPT_FPAR_PROF] > 0)
          {
            if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
            {
              err = 1;
            }
          }
          break;
        case SIM_MODTRAN_VMIE:
          sim_vmie = par[OPT_PAR_VMIE];
          break;
      }
      // common setup
      sim_wscl = par[OPT_PAR_WSCL];
      sim_oscl = par[OPT_PAR_OSCL];
      p[0] = w1-sim_xmgn;
      p[1] = w2+sim_xmgn;
      if(err) break;
      if(WriteTape5(rng,p,SIM_MODE_DSR_MODTRAN_BA) < 0)
      {
        err = 1;
        break;
      }
      if(RunModtran() < 0)
      {
        err = 1;
        break;
      }
      if(ReadDSR_B(rng,p) < 0)
      {
        err = 1;
        break;
      }
    }
    while(0);
    if(err)
    {
      *val = HUGE;
      return -1;
    }
  }

  nval = 0;
  *val = 0.0;
  if(mode > 0)
  {
    for(i=0; i<dsr_n_data; i++)
    {
      d1 = dsr_data[rng][i][0];
      d2 = dsr_data[rng][i][k];
      dr = (d2-d1)/(w2-w1);
      s1 = dsr_dsim[rng][i][0];
      s2 = dsr_dsim[rng][i][k];
      sr = (s2-s1)/(w2-w1);
      for(j=1; j<k; j++)
      {
        r = (dr*(dsr_wlen[rng][j]-w1)+d1)/(sr*(dsr_wlen[rng][j]-w1)+s1);
        v = dsr_data[rng][i][j]-r*dsr_dsim[rng][i][j];
        *val += v*v/(dsr_derr[rng][i][j]*dsr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  else
  {
    for(i=0; i<dsr_n_data; i++)
    {
      for(j=0; j<=k; j++)
      {
        v = dsr_data[rng][i][j]-dsr_dsim[rng][i][j];
        *val += v*v/(dsr_derr[rng][i][j]*dsr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  if(nval > 0)
  {
    *val /= nval;
  }

  if(cnt_vb)
  {
    PrintValue(par,val,2,(sim_mode==SIM_MODTRAN_VMIE?OPT_PAR_VMIE:OPT_PAR_VTAU),(mode==FCN_MODE_OZN?OPT_PAR_OSCL:OPT_PAR_WSCL));
  }

  return 0;
}

// Solar Aureole Fitting (MODTRAN aerosol model, single wavelength)
int AUR_FCN_MODTRAN_S(int rng,double *par,double *val,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j,k;
  int err;
  int fpar[OPT_N_FPAR];
  int nval;
  double p[1];
  double r,v;
  double d1,d2,dr;
  double s1,s2,sr;
  double w1,w2;

  if(CheckParam(par,fpar) < 0)
  {
    *val = HUGE;
    return -1;
  }

  k = aur_n_wlen[rng]-1;
  w1 = aur_wlen[rng][0];
  w2 = aur_wlen[rng][k];
  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    do
    {
      // individual setup
      switch(sim_mode)
      {
        case SIM_MODTRAN_VTAU:
          if(fpar[OPT_FPAR_PROF] > 0)
          {
            if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
            {
              err = 1;
            }
          }
          break;
        case SIM_MODTRAN_VMIE:
          sim_vmie = par[OPT_PAR_VMIE];
          break;
      }
      // common setup
      sim_wscl = par[OPT_PAR_WSCL];
      sim_oscl = par[OPT_PAR_OSCL];
      p[0] = (double)sim_nsam;
      if(err) break;
      if(WriteTape5(rng,p,SIM_MODE_AUR_MODTRAN_SA) < 0)
      {
        err = 1;
        break;
      }
      if(RunModtran() < 0)
      {
        err = 1;
        break;
      }
      if(ReadAUR_S(rng,p,aur_dsim[rng]) < 0)
      {
        err = 1;
        break;
      }
    }
    while(0);
    if(err)
    {
      *val = HUGE;
      return -1;
    }
    if(cnt_sflg == 1) // calculate single scattering
    {
      sim_mult = 0;
      err = 0;
      do
      {
        if(WriteTape5(rng,p,SIM_MODE_AUR_MODTRAN_SA) < 0)
        {
          err = 1;
          break;
        }
        if(RunModtran() < 0)
        {
          err = 1;
          break;
        }
        if(ReadAUR_S(rng,p,aur_dsgl[rng]) < 0)
        {
          err = 1;
          break;
        }
      }
      while(0);
      sim_mult = 1;
      if(err)
      {
        *val = HUGE;
        return -1;
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    // warning: correction factors are calculated for wavelength band
    // This function should be called with DIS==t
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<aur_n_data; i++)
        {
          for(j=0; j<aur_n_wlen[rng]; j++)
          {
            aur_dsim[rng][i][j] = aur_dsgl[rng][i][j]+(aur_dsim[rng][i][j]-aur_dsgl[rng][i][j])*aur_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<aur_n_data; i++)
        {
          for(j=0; j<aur_n_wlen[rng]; j++)
          {
            aur_dsim[rng][i][j] *= aur_dfac[rng][i][j];
          }
        }
      }
    }
  }

  nval = 0;
  *val = 0.0;
  if(mode > 0)
  {
    for(i=0; i<aur_n_data; i++)
    {
      d1 = aur_data[rng][i][0];
      d2 = aur_data[rng][i][k];
      dr = (d2-d1)/(w2-w1);
      s1 = aur_dsim[rng][i][0];
      s2 = aur_dsim[rng][i][k];
      sr = (s2-s1)/(w2-w1);
      for(j=1; j<k; j++)
      {
        r = (dr*(aur_wlen[rng][j]-w1)+d1)/(sr*(aur_wlen[rng][j]-w1)+s1);
        v = aur_data[rng][i][j]-r*aur_dsim[rng][i][j];
        *val += v*v/(aur_derr[rng][i][j]*aur_derr[rng][i][j]);
        nval++;
      }
    }
  }
  else
  {
    for(i=0; i<aur_n_data; i++)
    {
      for(j=0; j<=k; j++)
      {
        v = aur_data[rng][i][j]-aur_dsim[rng][i][j];
        *val += v*v/(aur_derr[rng][i][j]*aur_derr[rng][i][j]);
        nval++;
      }
    }
  }
  if(nval > 0)
  {
    *val /= nval;
  }

  if(cnt_vb)
  {
    PrintValue(par,val,2,(sim_mode==SIM_MODTRAN_VMIE?OPT_PAR_VMIE:OPT_PAR_VTAU),(mode==FCN_MODE_OZN?OPT_PAR_OSCL:OPT_PAR_WSCL));
  }

  return 0;
}

// Solar Aureole Fitting (MODTRAN aerosol model, wavelength band)
int AUR_FCN_MODTRAN_B(int rng,double *par,double *val,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j,k;
  int err;
  int fpar[OPT_N_FPAR];
  int nval;
  double p[1];
  double r,v;
  double d1,d2,dr;
  double s1,s2,sr;
  double w1,w2;

  if(CheckParam(par,fpar) < 0)
  {
    *val = HUGE;
    return -1;
  }

  k = aur_n_wlen[rng]-1;
  w1 = aur_wlen[rng][0];
  w2 = aur_wlen[rng][k];
  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    do
    {
      // individual setup
      switch(sim_mode)
      {
        case SIM_MODTRAN_VTAU:
          if(fpar[OPT_FPAR_PROF] > 0)
          {
            if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
            {
              err = 1;
            }
          }
          break;
        case SIM_MODTRAN_VMIE:
          sim_vmie = par[OPT_PAR_VMIE];
          break;
      }
      // common setup
      sim_wscl = par[OPT_PAR_WSCL];
      sim_oscl = par[OPT_PAR_OSCL];
      sim_sbmp = AUR_SBMP;
      p[0] = sim_xmgn;
      if(err) break;
      if(WriteTape5(rng,p,SIM_MODE_AUR_MODTRAN_BA) < 0)
      {
        err = 1;
        break;
      }
      if(RunModtran() < 0)
      {
        err = 1;
        break;
      }
      if(ReadAUR_B(rng,p,aur_dsim[rng]) < 0)
      {
        err = 1;
        break;
      }
    }
    while(0);
    if(err)
    {
      *val = HUGE;
      return -1;
    }
    if(cnt_sflg == 1) // calculate single scattering
    {
      sim_mult = 0;
      err = 0;
      do
      {
        if(WriteTape5(rng,p,SIM_MODE_AUR_MODTRAN_BA) < 0)
        {
          err = 1;
          break;
        }
        if(RunModtran() < 0)
        {
          err = 1;
          break;
        }
        if(ReadAUR_B(rng,p,aur_dsgl[rng]) < 0)
        {
          err = 1;
          break;
        }
      }
      while(0);
      sim_mult = 1;
      if(err)
      {
        *val = HUGE;
        return -1;
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<aur_n_data; i++)
        {
          for(j=0; j<aur_n_wlen[rng]; j++)
          {
            aur_dsim[rng][i][j] = aur_dsgl[rng][i][j]+(aur_dsim[rng][i][j]-aur_dsgl[rng][i][j])*aur_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<aur_n_data; i++)
        {
          for(j=0; j<aur_n_wlen[rng]; j++)
          {
            aur_dsim[rng][i][j] *= aur_dfac[rng][i][j];
          }
        }
      }
    }
  }

  nval = 0;
  *val = 0.0;
  if(mode > 0)
  {
    for(i=0; i<aur_n_data; i++)
    {
      d1 = aur_data[rng][i][0];
      d2 = aur_data[rng][i][k];
      dr = (d2-d1)/(w2-w1);
      s1 = aur_dsim[rng][i][0];
      s2 = aur_dsim[rng][i][k];
      sr = (s2-s1)/(w2-w1);
      for(j=1; j<k; j++)
      {
        r = (dr*(aur_wlen[rng][j]-w1)+d1)/(sr*(aur_wlen[rng][j]-w1)+s1);
        v = aur_data[rng][i][j]-r*aur_dsim[rng][i][j];
        *val += v*v/(aur_derr[rng][i][j]*aur_derr[rng][i][j]);
        nval++;
      }
    }
  }
  else
  {
    for(i=0; i<aur_n_data; i++)
    {
      for(j=0; j<=k; j++)
      {
        v = aur_data[rng][i][j]-aur_dsim[rng][i][j];
        *val += v*v/(aur_derr[rng][i][j]*aur_derr[rng][i][j]);
        nval++;
      }
    }
  }
  if(nval > 0)
  {
    *val /= nval;
  }

  if(cnt_vb)
  {
    PrintValue(par,val,2,(sim_mode==SIM_MODTRAN_VMIE?OPT_PAR_VMIE:OPT_PAR_VTAU),(mode==FCN_MODE_OZN?OPT_PAR_OSCL:OPT_PAR_WSCL));
  }

  return 0;
}

// Scattered Solar Radiation Fitting (MODTRAN aerosol model, single wavelength)
int SSR_FCN_MODTRAN_S(int rng,double *par,double *val,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j,k;
  int err;
  int fpar[OPT_N_FPAR];
  int nval;
  double p[1];
  double r,v;
  double d1,d2,dr;
  double s1,s2,sr;
  double w1,w2;

  if(CheckParam(par,fpar) < 0)
  {
    *val = HUGE;
    return -1;
  }

  k = ssr_n_wlen[rng]-1;
  w1 = ssr_wlen[rng][0];
  w2 = ssr_wlen[rng][k];
  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    do
    {
      // individual setup
      switch(sim_mode)
      {
        case SIM_MODTRAN_VTAU:
          if(fpar[OPT_FPAR_PROF] > 0)
          {
            if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
            {
              err = 1;
            }
          }
          break;
        case SIM_MODTRAN_VMIE:
          sim_vmie = par[OPT_PAR_VMIE];
          break;
      }
      // common setup
      sim_wscl = par[OPT_PAR_WSCL];
      sim_oscl = par[OPT_PAR_OSCL];
      p[0] = (double)sim_nsam;
      if(err) break;
      if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_SA) < 0)
      {
        err = 1;
        break;
      }
      if(RunModtran() < 0)
      {
        err = 1;
        break;
      }
      if(ReadSSR_S(rng,p,ssr_dsim[rng]) < 0)
      {
        err = 1;
        break;
      }
    }
    while(0);
    if(err)
    {
      *val = HUGE;
      return -1;
    }
    if(cnt_sflg == 1) // calculate single scattering
    {
      sim_mult = 0;
      err = 0;
      do
      {
        if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_SA) < 0)
        {
          err = 1;
          break;
        }
        if(RunModtran() < 0)
        {
          err = 1;
          break;
        }
        if(ReadSSR_S(rng,p,ssr_dsgl[rng]) < 0)
        {
          err = 1;
          break;
        }
      }
      while(0);
      sim_mult = 1;
      if(err)
      {
        *val = HUGE;
        return -1;
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    // warning: correction factors are calculated for wavelength band
    // This function should be called with DIS==t
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<ssr_n_data; i++)
        {
          for(j=0; j<ssr_n_wlen[rng]; j++)
          {
            ssr_dsim[rng][i][j] = ssr_dsgl[rng][i][j]+(ssr_dsim[rng][i][j]-ssr_dsgl[rng][i][j])*ssr_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<ssr_n_data; i++)
        {
          for(j=0; j<ssr_n_wlen[rng]; j++)
          {
            ssr_dsim[rng][i][j] *= ssr_dfac[rng][i][j];
          }
        }
      }
    }
  }

  nval = 0;
  *val = 0.0;
  if(mode > 0)
  {
    for(i=0; i<ssr_n_data; i++)
    {
      d1 = ssr_data[rng][i][0];
      d2 = ssr_data[rng][i][k];
      dr = (d2-d1)/(w2-w1);
      s1 = ssr_dsim[rng][i][0];
      s2 = ssr_dsim[rng][i][k];
      sr = (s2-s1)/(w2-w1);
      for(j=1; j<k; j++)
      {
        r = (dr*(ssr_wlen[rng][j]-w1)+d1)/(sr*(ssr_wlen[rng][j]-w1)+s1);
        v = ssr_data[rng][i][j]-r*ssr_dsim[rng][i][j];
        *val += v*v/(ssr_derr[rng][i][j]*ssr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  else
  {
    for(i=0; i<ssr_n_data; i++)
    {
      for(j=0; j<=k; j++)
      {
        v = ssr_data[rng][i][j]-ssr_dsim[rng][i][j];
        *val += v*v/(ssr_derr[rng][i][j]*ssr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  if(nval > 0)
  {
    *val /= nval;
  }

  if(cnt_vb)
  {
    PrintValue(par,val,2,(sim_mode==SIM_MODTRAN_VMIE?OPT_PAR_VMIE:OPT_PAR_VTAU),(mode==FCN_MODE_OZN?OPT_PAR_OSCL:OPT_PAR_WSCL));
  }

  return 0;
}

// Scattered Solar Radiation Fitting (MODTRAN aerosol model, wavelength band)
int SSR_FCN_MODTRAN_B(int rng,double *par,double *val,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j,k;
  int err;
  int fpar[OPT_N_FPAR];
  int nval;
  double p[1];
  double r,v;
  double d1,d2,dr;
  double s1,s2,sr;
  double w1,w2;

  if(CheckParam(par,fpar) < 0)
  {
    *val = HUGE;
    return -1;
  }

  k = ssr_n_wlen[rng]-1;
  w1 = ssr_wlen[rng][0];
  w2 = ssr_wlen[rng][k];
  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    do
    {
      // individual setup
      switch(sim_mode)
      {
        case SIM_MODTRAN_VTAU:
          if(fpar[OPT_FPAR_PROF] > 0)
          {
            if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
            {
              err = 1;
            }
          }
          break;
        case SIM_MODTRAN_VMIE:
          sim_vmie = par[OPT_PAR_VMIE];
          break;
      }
      // common setup
      sim_wscl = par[OPT_PAR_WSCL];
      sim_oscl = par[OPT_PAR_OSCL];
      sim_sbmp = SSR_SBMP;
      p[0] = sim_xmgn;
      if(err) break;
      if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_BA) < 0)
      {
        err = 1;
        break;
      }
      if(RunModtran() < 0)
      {
        err = 1;
        break;
      }
      if(ReadSSR_B(rng,p,ssr_dsim[rng]) < 0)
      {
        err = 1;
        break;
      }
    }
    while(0);
    if(err)
    {
      *val = HUGE;
      return -1;
    }
    if(cnt_sflg == 1) // calculate single scattering
    {
      sim_mult = 0;
      err = 0;
      do
      {
        if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_BA) < 0)
        {
          err = 1;
          break;
        }
        if(RunModtran() < 0)
        {
          err = 1;
          break;
        }
        if(ReadSSR_B(rng,p,ssr_dsgl[rng]) < 0)
        {
          err = 1;
          break;
        }
      }
      while(0);
      sim_mult = 1;
      if(err)
      {
        *val = HUGE;
        return -1;
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<ssr_n_data; i++)
        {
          for(j=0; j<ssr_n_wlen[rng]; j++)
          {
            ssr_dsim[rng][i][j] = ssr_dsgl[rng][i][j]+(ssr_dsim[rng][i][j]-ssr_dsgl[rng][i][j])*ssr_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<ssr_n_data; i++)
        {
          for(j=0; j<ssr_n_wlen[rng]; j++)
          {
            ssr_dsim[rng][i][j] *= ssr_dfac[rng][i][j];
          }
        }
      }
    }
  }

  nval = 0;
  *val = 0.0;
  if(mode > 0)
  {
    for(i=0; i<ssr_n_data; i++)
    {
      d1 = ssr_data[rng][i][0];
      d2 = ssr_data[rng][i][k];
      dr = (d2-d1)/(w2-w1);
      s1 = ssr_dsim[rng][i][0];
      s2 = ssr_dsim[rng][i][k];
      sr = (s2-s1)/(w2-w1);
      for(j=1; j<k; j++)
      {
        r = (dr*(ssr_wlen[rng][j]-w1)+d1)/(sr*(ssr_wlen[rng][j]-w1)+s1);
        v = ssr_data[rng][i][j]-r*ssr_dsim[rng][i][j];
        *val += v*v/(ssr_derr[rng][i][j]*ssr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  else
  {
    for(i=0; i<ssr_n_data; i++)
    {
      for(j=0; j<=k; j++)
      {
        v = ssr_data[rng][i][j]-ssr_dsim[rng][i][j];
        *val += v*v/(ssr_derr[rng][i][j]*ssr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  if(nval > 0)
  {
    *val /= nval;
  }

  if(cnt_vb)
  {
    PrintValue(par,val,2,(sim_mode==SIM_MODTRAN_VMIE?OPT_PAR_VMIE:OPT_PAR_VTAU),(mode==FCN_MODE_OZN?OPT_PAR_OSCL:OPT_PAR_WSCL));
  }

  return 0;
}

// Direct Solar Radiation Fitting (Custom aerosol model, wavelength band)
int DSR_FCN_CUSTOM_B(int rng,double *par,double *val,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j,k;
  int err;
  int fpar[OPT_N_FPAR];
  int nval;
  double p[2];
  double r,v;
  double d1,d2,dr;
  double s1,s2,sr;
  double w1,w2;
  double dif;

  err = 0;
  do
  {
    if(CheckParam(par,fpar) < 0)
    {
      err = 1;
      break;
    }
    if(fpar[OPT_FPAR_DIST] > 0)
    {
      if(MieCalc(par,fpar) < 0)
      {
        err = 1;
        break;
      }
    }
    if(fpar[OPT_FPAR_MIXR] > 0)
    {
      if(MixComp(par) < 0)
      {
        err = 1;
        break;
      }
      if(cnt_vb > 2) MiePrintout();
    }
  }
  while(0);
  if(err)
  {
    *val = HUGE;
    return -1;
  }

  k = dsr_n_wlen[rng]-1;
  w1 = dsr_wlen[rng][0];
  w2 = dsr_wlen[rng][k];
  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    do
    {
      // individual setup
      switch(sim_mode)
      {
        case SIM_MODTRAN_VTAU:
          if(fpar[OPT_FPAR_PROF] > 0)
          {
            if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
            {
              err = 1;
            }
          }
          break;
        case SIM_MODTRAN_VMIE:
          sim_vmie = par[OPT_PAR_VMIE];
          break;
      }
      // common setup
      sim_wscl = par[OPT_PAR_WSCL];
      sim_oscl = par[OPT_PAR_OSCL];
      p[0] = w1-sim_xmgn;
      p[1] = w2+sim_xmgn;
      if(err) break;
      if(WriteTape5(rng,p,SIM_MODE_DSR_CUSTOM_BA) < 0)
      {
        err = 1;
        break;
      }
      if(RunModtran() < 0)
      {
        err = 1;
        break;
      }
      if(ReadDSR_B(rng,p) < 0)
      {
        err = 1;
        break;
      }
    }
    while(0);
    if(err)
    {
      *val = HUGE;
      return -1;
    }
  }

  nval = 0;
  *val = 0.0;
  if(mode > 0)
  {
    for(i=0; i<dsr_n_data; i++)
    {
      d1 = dsr_data[rng][i][0];
      d2 = dsr_data[rng][i][k];
      dr = (d2-d1)/(w2-w1);
      s1 = dsr_dsim[rng][i][0];
      s2 = dsr_dsim[rng][i][k];
      sr = (s2-s1)/(w2-w1);
      for(j=1; j<k; j++)
      {
        r = (dr*(dsr_wlen[rng][j]-w1)+d1)/(sr*(dsr_wlen[rng][j]-w1)+s1);
        v = dsr_data[rng][i][j]-r*dsr_dsim[rng][i][j];
        *val += v*v/(dsr_derr[rng][i][j]*dsr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  else
  {
    for(i=0; i<dsr_n_data; i++)
    {
      for(j=0; j<=k; j++)
      {
        v = dsr_data[rng][i][j]-dsr_dsim[rng][i][j];
        *val += v*v/(dsr_derr[rng][i][j]*dsr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  if(nval > 0)
  {
    *val /= nval;
  }

  if(mode > 0)
  {
    if(opt_n_opt >= 0)
    {
      for(i=0; i<dsr_n_data; i++)
      {
        d1 = dsr_data[rng][i][0];
        d2 = dsr_data[rng][i][k];
        dr = (d2-d1)/(w2-w1);
        s1 = dsr_dsim[rng][i][0];
        s2 = dsr_dsim[rng][i][k];
        sr = (s2-s1)/(w2-w1);
        dif = 0.0;
        for(j=1; j<k; j++)
        {
          r = (dr*(dsr_wlen[rng][j]-w1)+d1)/(sr*(dsr_wlen[rng][j]-w1)+s1);
          v = dsr_data[rng][i][j]-r*dsr_dsim[rng][i][j];
          dif += v*v/(dsr_derr[rng][i][j]*dsr_derr[rng][i][j]);
        }
        for(j=0; j<=k; j++)
        {
          r = (dr*(dsr_wlen[rng][j]-w1)+d1)/(sr*(dsr_wlen[rng][j]-w1)+s1);
          message(stderr,"%5d %2d %2d %13.4f %13.6e %13.6e %18.10e %18.10e %18.10e\n",
                          opt_n_opt,i,j,dsr_wlen[rng][j],dsr_data[rng][i][j],r*dsr_dsim[rng][i][j],
                          par[OPT_PAR_WSCL],dif,*val);
        }
      }
      opt_n_opt++;
    }
  }

  if(cnt_vb)
  {
    PrintValue(par,val,11,(sim_mode==SIM_MODTRAN_VMIE?OPT_PAR_VMIE:OPT_PAR_VTAU),(mode==FCN_MODE_OZN?OPT_PAR_OSCL:OPT_PAR_WSCL),
               OPT_PAR_NCMP,OPT_PAR_LMD1,OPT_PAR_LSG1,
               OPT_PAR_MXR2,OPT_PAR_LMD2,OPT_PAR_LSG2,
               OPT_PAR_MXR3,OPT_PAR_LMD3,OPT_PAR_LSG3);
  }

  return 0;
}

// Solar Aureole Fitting (Custom aerosol model, single wavelength)
int AUR_FCN_CUSTOM_S(int rng,double *par,double *val,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j,k;
  int err;
  int fpar[OPT_N_FPAR];
  int nval;
  double p[1];
  double r,v;
  double d1,d2,dr;
  double s1,s2,sr;
  double w1,w2;

  err = 0;
  do
  {
    if(CheckParam(par,fpar) < 0)
    {
      err = 1;
      break;
    }
    if(fpar[OPT_FPAR_DIST] > 0)
    {
      if(MieCalc(par,fpar) < 0)
      {
        err = 1;
        break;
      }
    }
    if(fpar[OPT_FPAR_MIXR] > 0)
    {
      if(MixComp(par) < 0)
      {
        err = 1;
        break;
      }
      if(cnt_vb > 2) MiePrintout();
    }
  }
  while(0);
  if(err)
  {
    *val = HUGE;
    return -1;
  }

  k = aur_n_wlen[rng]-1;
  w1 = aur_wlen[rng][0];
  w2 = aur_wlen[rng][k];
  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    do
    {
      // individual setup
      switch(sim_mode)
      {
        case SIM_MODTRAN_VTAU:
          if(fpar[OPT_FPAR_PROF] > 0)
          {
            if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
            {
              err = 1;
            }
          }
          break;
        case SIM_MODTRAN_VMIE:
          sim_vmie = par[OPT_PAR_VMIE];
          break;
      }
      // common setup
      sim_wscl = par[OPT_PAR_WSCL];
      sim_oscl = par[OPT_PAR_OSCL];
      p[0] = (double)sim_nsam;
      if(err) break;
      if(WriteTape5(rng,p,SIM_MODE_AUR_CUSTOM_SA) < 0)
      {
        err = 1;
        break;
      }
      if(RunModtran() < 0)
      {
        err = 1;
        break;
      }
      if(ReadAUR_S(rng,p,aur_dsim[rng]) < 0)
      {
        err = 1;
        break;
      }
    }
    while(0);
    if(err)
    {
      *val = HUGE;
      return -1;
    }
    if(cnt_sflg == 1) // calculate single scattering
    {
      sim_mult = 0;
      err = 0;
      do
      {
        if(WriteTape5(rng,p,SIM_MODE_AUR_CUSTOM_SA) < 0)
        {
          err = 1;
          break;
        }
        if(RunModtran() < 0)
        {
          err = 1;
          break;
        }
        if(ReadAUR_S(rng,p,aur_dsgl[rng]) < 0)
        {
          err = 1;
          break;
        }
      }
      while(0);
      sim_mult = 1;
      if(err)
      {
        *val = HUGE;
        return -1;
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    // warning: correction factors are calculated for wavelength band
    // This function should be called with DIS==t
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<aur_n_data; i++)
        {
          for(j=0; j<aur_n_wlen[rng]; j++)
          {
            aur_dsim[rng][i][j] = aur_dsgl[rng][i][j]+(aur_dsim[rng][i][j]-aur_dsgl[rng][i][j])*aur_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<aur_n_data; i++)
        {
          for(j=0; j<aur_n_wlen[rng]; j++)
          {
            aur_dsim[rng][i][j] *= aur_dfac[rng][i][j];
          }
        }
      }
    }
  }

  nval = 0;
  *val = 0.0;
  if(mode > 0)
  {
    for(i=0; i<aur_n_data; i++)
    {
      d1 = aur_data[rng][i][0];
      d2 = aur_data[rng][i][k];
      dr = (d2-d1)/(w2-w1);
      s1 = aur_dsim[rng][i][0];
      s2 = aur_dsim[rng][i][k];
      sr = (s2-s1)/(w2-w1);
      for(j=1; j<k; j++)
      {
        r = (dr*(aur_wlen[rng][j]-w1)+d1)/(sr*(aur_wlen[rng][j]-w1)+s1);
        v = aur_data[rng][i][j]-r*aur_dsim[rng][i][j];
        *val += v*v/(aur_derr[rng][i][j]*aur_derr[rng][i][j]);
        nval++;
      }
    }
  }
  else
  {
    for(i=0; i<aur_n_data; i++)
    {
      for(j=0; j<=k; j++)
      {
        v = aur_data[rng][i][j]-aur_dsim[rng][i][j];
        *val += v*v/(aur_derr[rng][i][j]*aur_derr[rng][i][j]);
        nval++;
      }
    }
  }
  if(nval > 0)
  {
    *val /= nval;
  }

  if(cnt_vb)
  {
    PrintValue(par,val,11,(sim_mode==SIM_MODTRAN_VMIE?OPT_PAR_VMIE:OPT_PAR_VTAU),(mode==FCN_MODE_OZN?OPT_PAR_OSCL:OPT_PAR_WSCL),
               OPT_PAR_NCMP,OPT_PAR_LMD1,OPT_PAR_LSG1,
               OPT_PAR_MXR2,OPT_PAR_LMD2,OPT_PAR_LSG2,
               OPT_PAR_MXR3,OPT_PAR_LMD3,OPT_PAR_LSG3);
  }

  return 0;
}

// Solar Aureole Fitting (Custom aerosol model, wavelength band)
int AUR_FCN_CUSTOM_B(int rng,double *par,double *val,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j,k;
  int err;
  int fpar[OPT_N_FPAR];
  int nval;
  double p[1];
  double r,v;
  double d1,d2,dr;
  double s1,s2,sr;
  double w1,w2;

  err = 0;
  do
  {
    if(CheckParam(par,fpar) < 0)
    {
      err = 1;
      break;
    }
    if(fpar[OPT_FPAR_DIST] > 0)
    {
      if(MieCalc(par,fpar) < 0)
      {
        err = 1;
        break;
      }
    }
    if(fpar[OPT_FPAR_MIXR] > 0)
    {
      if(MixComp(par) < 0)
      {
        err = 1;
        break;
      }
      if(cnt_vb > 2) MiePrintout();
    }
  }
  while(0);
  if(err)
  {
    *val = HUGE;
    return -1;
  }

  k = aur_n_wlen[rng]-1;
  w1 = aur_wlen[rng][0];
  w2 = aur_wlen[rng][k];
  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    do
    {
      // individual setup
      switch(sim_mode)
      {
        case SIM_MODTRAN_VTAU:
          if(fpar[OPT_FPAR_PROF] > 0)
          {
            if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
            {
              err = 1;
            }
          }
          break;
        case SIM_MODTRAN_VMIE:
          sim_vmie = par[OPT_PAR_VMIE];
          break;
      }
      // common setup
      sim_wscl = par[OPT_PAR_WSCL];
      sim_oscl = par[OPT_PAR_OSCL];
      sim_sbmp = AUR_SBMP;
      p[0] = sim_xmgn;
      if(err) break;
      if(WriteTape5(rng,p,SIM_MODE_AUR_CUSTOM_BA) < 0)
      {
        err = 1;
        break;
      }
      if(RunModtran() < 0)
      {
        err = 1;
        break;
      }
      if(ReadAUR_B(rng,p,aur_dsim[rng]) < 0)
      {
        err = 1;
        break;
      }
    }
    while(0);
    if(err)
    {
      *val = HUGE;
      return -1;
    }
    if(cnt_sflg == 1) // calculate single scattering
    {
      sim_mult = 0;
      err = 0;
      do
      {
        if(WriteTape5(rng,p,SIM_MODE_AUR_CUSTOM_BA) < 0)
        {
          err = 1;
          break;
        }
        if(RunModtran() < 0)
        {
          err = 1;
          break;
        }
        if(ReadAUR_B(rng,p,aur_dsgl[rng]) < 0)
        {
          err = 1;
          break;
        }
      }
      while(0);
      sim_mult = 1;
      if(err)
      {
        *val = HUGE;
        return -1;
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<aur_n_data; i++)
        {
          for(j=0; j<aur_n_wlen[rng]; j++)
          {
            aur_dsim[rng][i][j] = aur_dsgl[rng][i][j]+(aur_dsim[rng][i][j]-aur_dsgl[rng][i][j])*aur_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<aur_n_data; i++)
        {
          for(j=0; j<aur_n_wlen[rng]; j++)
          {
            aur_dsim[rng][i][j] *= aur_dfac[rng][i][j];
          }
        }
      }
    }
  }

  nval = 0;
  *val = 0.0;
  if(mode > 0)
  {
    for(i=0; i<aur_n_data; i++)
    {
      d1 = aur_data[rng][i][0];
      d2 = aur_data[rng][i][k];
      dr = (d2-d1)/(w2-w1);
      s1 = aur_dsim[rng][i][0];
      s2 = aur_dsim[rng][i][k];
      sr = (s2-s1)/(w2-w1);
      for(j=1; j<k; j++)
      {
        r = (dr*(aur_wlen[rng][j]-w1)+d1)/(sr*(aur_wlen[rng][j]-w1)+s1);
        v = aur_data[rng][i][j]-r*aur_dsim[rng][i][j];
        *val += v*v/(aur_derr[rng][i][j]*aur_derr[rng][i][j]);
        nval++;
      }
    }
  }
  else
  {
    for(i=0; i<aur_n_data; i++)
    {
      for(j=0; j<=k; j++)
      {
        v = aur_data[rng][i][j]-aur_dsim[rng][i][j];
        *val += v*v/(aur_derr[rng][i][j]*aur_derr[rng][i][j]);
        nval++;
      }
    }
  }
  if(nval > 0)
  {
    *val /= nval;
  }

  if(cnt_vb)
  {
    PrintValue(par,val,11,(sim_mode==SIM_MODTRAN_VMIE?OPT_PAR_VMIE:OPT_PAR_VTAU),(mode==FCN_MODE_OZN?OPT_PAR_OSCL:OPT_PAR_WSCL),
               OPT_PAR_NCMP,OPT_PAR_LMD1,OPT_PAR_LSG1,
               OPT_PAR_MXR2,OPT_PAR_LMD2,OPT_PAR_LSG2,
               OPT_PAR_MXR3,OPT_PAR_LMD3,OPT_PAR_LSG3);
  }

  return 0;
}

// Scattered Solar Radiation Fitting (Custom aerosol model, single wavelength)
int SSR_FCN_CUSTOM_S(int rng,double *par,double *val,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j,k;
  int err;
  int fpar[OPT_N_FPAR];
  int nval;
  double p[1];
  double r,v;
  double d1,d2,dr;
  double s1,s2,sr;
  double w1,w2;

  err = 0;
  do
  {
    if(CheckParam(par,fpar) < 0)
    {
      err = 1;
      break;
    }
    if(fpar[OPT_FPAR_DIST] > 0)
    {
      if(MieCalc(par,fpar) < 0)
      {
        err = 1;
        break;
      }
    }
    if(fpar[OPT_FPAR_MIXR] > 0)
    {
      if(MixComp(par) < 0)
      {
        err = 1;
        break;
      }
      if(cnt_vb > 2) MiePrintout();
    }
  }
  while(0);
  if(err)
  {
    *val = HUGE;
    return -1;
  }

  k = ssr_n_wlen[rng]-1;
  w1 = ssr_wlen[rng][0];
  w2 = ssr_wlen[rng][k];
  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    do
    {
      // individual setup
      switch(sim_mode)
      {
        case SIM_MODTRAN_VTAU:
          if(fpar[OPT_FPAR_PROF] > 0)
          {
            if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
            {
              err = 1;
            }
          }
          break;
        case SIM_MODTRAN_VMIE:
          sim_vmie = par[OPT_PAR_VMIE];
          break;
      }
      // common setup
      sim_wscl = par[OPT_PAR_WSCL];
      sim_oscl = par[OPT_PAR_OSCL];
      p[0] = (double)sim_nsam;
      if(err) break;
      if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_SA) < 0)
      {
        err = 1;
        break;
      }
      if(RunModtran() < 0)
      {
        err = 1;
        break;
      }
      if(ReadSSR_S(rng,p,ssr_dsim[rng]) < 0)
      {
        err = 1;
        break;
      }
    }
    while(0);
    if(err)
    {
      *val = HUGE;
      return -1;
    }
    if(cnt_sflg == 1) // calculate single scattering
    {
      sim_mult = 0;
      err = 0;
      do
      {
        if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_SA) < 0)
        {
          err = 1;
          break;
        }
        if(RunModtran() < 0)
        {
          err = 1;
          break;
        }
        if(ReadSSR_S(rng,p,ssr_dsgl[rng]) < 0)
        {
          err = 1;
          break;
        }
      }
      while(0);
      sim_mult = 1;
      if(err)
      {
        *val = HUGE;
        return -1;
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    // warning: correction factors are calculated for wavelength band
    // This function should be called with DIS==t
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<ssr_n_data; i++)
        {
          for(j=0; j<ssr_n_wlen[rng]; j++)
          {
            ssr_dsim[rng][i][j] = ssr_dsgl[rng][i][j]+(ssr_dsim[rng][i][j]-ssr_dsgl[rng][i][j])*ssr_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<ssr_n_data; i++)
        {
          for(j=0; j<ssr_n_wlen[rng]; j++)
          {
            ssr_dsim[rng][i][j] *= ssr_dfac[rng][i][j];
          }
        }
      }
    }
  }

  nval = 0;
  *val = 0.0;
  if(mode > 0)
  {
    for(i=0; i<ssr_n_data; i++)
    {
      d1 = ssr_data[rng][i][0];
      d2 = ssr_data[rng][i][k];
      dr = (d2-d1)/(w2-w1);
      s1 = ssr_dsim[rng][i][0];
      s2 = ssr_dsim[rng][i][k];
      sr = (s2-s1)/(w2-w1);
      for(j=1; j<k; j++)
      {
        r = (dr*(ssr_wlen[rng][j]-w1)+d1)/(sr*(ssr_wlen[rng][j]-w1)+s1);
        v = ssr_data[rng][i][j]-r*ssr_dsim[rng][i][j];
        *val += v*v/(ssr_derr[rng][i][j]*ssr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  else
  {
    for(i=0; i<ssr_n_data; i++)
    {
      for(j=0; j<=k; j++)
      {
        v = ssr_data[rng][i][j]-ssr_dsim[rng][i][j];
        *val += v*v/(ssr_derr[rng][i][j]*ssr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  if(nval > 0)
  {
    *val /= nval;
  }

  if(cnt_vb)
  {
    PrintValue(par,val,11,(sim_mode==SIM_MODTRAN_VMIE?OPT_PAR_VMIE:OPT_PAR_VTAU),(mode==FCN_MODE_OZN?OPT_PAR_OSCL:OPT_PAR_WSCL),
               OPT_PAR_NCMP,OPT_PAR_LMD1,OPT_PAR_LSG1,
               OPT_PAR_MXR2,OPT_PAR_LMD2,OPT_PAR_LSG2,
               OPT_PAR_MXR3,OPT_PAR_LMD3,OPT_PAR_LSG3);
  }

  return 0;
}

// Scattered Solar Radiation Fitting (Custom aerosol model, wavelength band)
int SSR_FCN_CUSTOM_B(int rng,double *par,double *val,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j,k;
  int err;
  int fpar[OPT_N_FPAR];
  int nval;
  double p[1];
  double r,v;
  double d1,d2,dr;
  double s1,s2,sr;
  double w1,w2;

  err = 0;
  do
  {
    if(CheckParam(par,fpar) < 0)
    {
      err = 1;
      break;
    }
    if(fpar[OPT_FPAR_DIST] > 0)
    {
      if(MieCalc(par,fpar) < 0)
      {
        err = 1;
        break;
      }
    }
    if(fpar[OPT_FPAR_MIXR] > 0)
    {
      if(MixComp(par) < 0)
      {
        err = 1;
        break;
      }
      if(cnt_vb > 2) MiePrintout();
    }
  }
  while(0);
  if(err)
  {
    *val = HUGE;
    return -1;
  }

  k = ssr_n_wlen[rng]-1;
  w1 = ssr_wlen[rng][0];
  w2 = ssr_wlen[rng][k];
  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    do
    {
      // individual setup
      switch(sim_mode)
      {
        case SIM_MODTRAN_VTAU:
          if(fpar[OPT_FPAR_PROF] > 0)
          {
            if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
            {
              err = 1;
            }
          }
          break;
        case SIM_MODTRAN_VMIE:
          sim_vmie = par[OPT_PAR_VMIE];
          break;
      }
      // common setup
      sim_wscl = par[OPT_PAR_WSCL];
      sim_oscl = par[OPT_PAR_OSCL];
      sim_sbmp = SSR_SBMP;
      p[0] = sim_xmgn;
      if(err) break;
      if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_BA) < 0)
      {
        err = 1;
        break;
      }
      if(RunModtran() < 0)
      {
        err = 1;
        break;
      }
      if(ReadSSR_B(rng,p,ssr_dsim[rng]) < 0)
      {
        err = 1;
        break;
      }
    }
    while(0);
    if(err)
    {
      *val = HUGE;
      return -1;
    }
    if(cnt_sflg == 1) // calculate single scattering
    {
      sim_mult = 0;
      err = 0;
      do
      {
        if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_BA) < 0)
        {
          err = 1;
          break;
        }
        if(RunModtran() < 0)
        {
          err = 1;
          break;
        }
        if(ReadSSR_B(rng,p,ssr_dsgl[rng]) < 0)
        {
          err = 1;
          break;
        }
      }
      while(0);
      sim_mult = 1;
      if(err)
      {
        *val = HUGE;
        return -1;
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<ssr_n_data; i++)
        {
          for(j=0; j<ssr_n_wlen[rng]; j++)
          {
            ssr_dsim[rng][i][j] = ssr_dsgl[rng][i][j]+(ssr_dsim[rng][i][j]-ssr_dsgl[rng][i][j])*ssr_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<ssr_n_data; i++)
        {
          for(j=0; j<ssr_n_wlen[rng]; j++)
          {
            ssr_dsim[rng][i][j] *= ssr_dfac[rng][i][j];
          }
        }
      }
    }
  }

  nval = 0;
  *val = 0.0;
  if(mode > 0)
  {
    for(i=0; i<ssr_n_data; i++)
    {
      d1 = ssr_data[rng][i][0];
      d2 = ssr_data[rng][i][k];
      dr = (d2-d1)/(w2-w1);
      s1 = ssr_dsim[rng][i][0];
      s2 = ssr_dsim[rng][i][k];
      sr = (s2-s1)/(w2-w1);
      for(j=1; j<k; j++)
      {
        r = (dr*(ssr_wlen[rng][j]-w1)+d1)/(sr*(ssr_wlen[rng][j]-w1)+s1);
        v = ssr_data[rng][i][j]-r*ssr_dsim[rng][i][j];
        *val += v*v/(ssr_derr[rng][i][j]*ssr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  else
  {
    for(i=0; i<ssr_n_data; i++)
    {
      for(j=0; j<=k; j++)
      {
        v = ssr_data[rng][i][j]-ssr_dsim[rng][i][j];
        *val += v*v/(ssr_derr[rng][i][j]*ssr_derr[rng][i][j]);
        nval++;
      }
    }
  }
  if(nval > 0)
  {
    *val /= nval;
  }

  if(cnt_vb)
  {
    PrintValue(par,val,11,(sim_mode==SIM_MODTRAN_VMIE?OPT_PAR_VMIE:OPT_PAR_VTAU),(mode==FCN_MODE_OZN?OPT_PAR_OSCL:OPT_PAR_WSCL),
               OPT_PAR_NCMP,OPT_PAR_LMD1,OPT_PAR_LSG1,
               OPT_PAR_MXR2,OPT_PAR_LMD2,OPT_PAR_LSG2,
               OPT_PAR_MXR3,OPT_PAR_LMD3,OPT_PAR_LSG3);
  }

  return 0;
}

// Scattered Solar Radiation at DSR wavelengths (MODTRAN aerosol model, single wavelength)
int DSR_SSR_MODTRAN_S(int rng,double *par,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j;
  int err;
  int fpar[OPT_N_FPAR];
  double p[6];
  double r;

  if(CheckParam(par,fpar) < 0)
  {
    return -1;
  }

  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    // individual setup
    switch(sim_mode)
    {
      case SIM_MODTRAN_VTAU:
        if(fpar[OPT_FPAR_PROF] > 0)
        {
          if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
          {
            err = 1;
          }
        }
        break;
      case SIM_MODTRAN_VMIE:
        sim_vmie = par[OPT_PAR_VMIE];
        break;
    }
    // common setup
    sim_wscl = par[OPT_PAR_WSCL];
    sim_oscl = par[OPT_PAR_OSCL];
    sim_sbmp = SSR_SBMP;
    if(err)
    {
      return -1;
    }
    for(i=0; i<dsr_n_data; i++)
    {
      p[2] = dsr_th_sun[i];
      p[3] = dsr_ph_sun[i];
      p[4] = dsr_th_sun[i];
      p[5] = dsr_ph_sun[i];
      for(j=0; j<dsr_n_wlen[rng]; j++)
      {
        p[0] = dsr_wlen[rng][j];
        p[1] = (double)sim_nsam;
        err = 0;
        do
        {
          if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_SO) < 0)
          {
            err = 1;
            break;
          }
          if(RunModtran() < 0)
          {
            err = 1;
            break;
          }
          if(ReadPlt_S(p,&dsr_dsim[rng][i][j]) < 0)
          {
            err = 1;
            break;
          }
          if(ResFactor(dsr_ms720,dsr_nw[rng][j],sim_sbmp,dsr_wlen[rng][j],dsr_amas[i],&r) < 0)
          {
            err = 1;
            break;
          }
          dsr_dsim[rng][i][j] *= r;
        }
        while(0);
        if(err)
        {
          return -1;
        }
        if(cnt_sflg == 1) // calculate single scattering
        {
          sim_mult = 0;
          err = 0;
          do
          {
            if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_SO) < 0)
            {
              err = 1;
              break;
            }
            if(RunModtran() < 0)
            {
              err = 1;
              break;
            }
            if(ReadPlt_S(p,&dsr_dsgl[rng][i][j]) < 0)
            {
              err = 1;
              break;
            }
            if(ResFactor(dsr_ms720,dsr_nw[rng][j],sim_sbmp,dsr_wlen[rng][j],dsr_amas[i],&r) < 0)
            {
              err = 1;
              break;
            }
            dsr_dsgl[rng][i][j] *= r;
          }
          while(0);
          sim_mult = 1;
          if(err)
          {
            return -1;
          }
        }
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    // warning: correction factors are calculated for wavelength band
    // This function should be called with DIS==t
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<dsr_n_data; i++)
        {
          for(j=0; j<dsr_n_wlen[rng]; j++)
          {
            dsr_dsim[rng][i][j] = dsr_dsgl[rng][i][j]+(dsr_dsim[rng][i][j]-dsr_dsgl[rng][i][j])*dsr_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<dsr_n_data; i++)
        {
          for(j=0; j<dsr_n_wlen[rng]; j++)
          {
            dsr_dsim[rng][i][j] *= dsr_dfac[rng][i][j];
          }
        }
      }
    }
  }

  return 0;
}

// Scattered Solar Radiation at DSR wavelengths (MODTRAN aerosol model, wavelength band)
int DSR_SSR_MODTRAN_B(int rng,double *par,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j;
  int err;
  int fpar[OPT_N_FPAR];
  double p[6];

  if(CheckParam(par,fpar) < 0)
  {
    return -1;
  }

  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    // individual setup
    switch(sim_mode)
    {
      case SIM_MODTRAN_VTAU:
        if(fpar[OPT_FPAR_PROF] > 0)
        {
          if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
          {
            err = 1;
          }
        }
        break;
      case SIM_MODTRAN_VMIE:
        sim_vmie = par[OPT_PAR_VMIE];
        break;
    }
    // common setup
    sim_wscl = par[OPT_PAR_WSCL];
    sim_oscl = par[OPT_PAR_OSCL];
    sim_sbmp = SSR_SBMP;
    if(err)
    {
      return -1;
    }
    for(i=0; i<dsr_n_data; i++)
    {
      p[2] = dsr_th_sun[i];
      p[3] = dsr_ph_sun[i];
      p[4] = dsr_th_sun[i];
      p[5] = dsr_ph_sun[i];
      for(j=0; j<dsr_n_wlen[rng]; j++)
      {
        p[0] = dsr_wlen[rng][j]-sim_xmgn;
        p[1] = dsr_wlen[rng][j]+sim_xmgn;
        err = 0;
        do
        {
          if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_BO) < 0)
          {
            err = 1;
            break;
          }
          if(RunModtran() < 0)
          {
            err = 1;
            break;
          }
          if(ReadPlt_B(p,1,&dsr_wlen[rng][j],&dsr_dsim[rng][i][j]) < 0)
          {
            err = 1;
            break;
          }
        }
        while(0);
        if(err)
        {
          return -1;
        }
        if(cnt_sflg == 1) // calculate single scattering
        {
          sim_mult = 0;
          err = 0;
          do
          {
            if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_BO) < 0)
            {
              err = 1;
              break;
            }
            if(RunModtran() < 0)
            {
              err = 1;
              break;
            }
            if(ReadPlt_B(p,1,&dsr_wlen[rng][j],&dsr_dsgl[rng][i][j]) < 0)
            {
              err = 1;
              break;
            }
          }
          while(0);
          sim_mult = 1;
          if(err)
          {
            return -1;
          }
        }
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<dsr_n_data; i++)
        {
          for(j=0; j<dsr_n_wlen[rng]; j++)
          {
            dsr_dsim[rng][i][j] = dsr_dsgl[rng][i][j]+(dsr_dsim[rng][i][j]-dsr_dsgl[rng][i][j])*dsr_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<dsr_n_data; i++)
        {
          for(j=0; j<dsr_n_wlen[rng]; j++)
          {
            dsr_dsim[rng][i][j] *= dsr_dfac[rng][i][j];
          }
        }
      }
    }
  }

  return 0;
}

// Scattered Solar Radiation at DSR wavelengths (Custom aerosol model, single wavelength)
int DSR_SSR_CUSTOM_S(int rng,double *par,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j;
  int err;
  int fpar[OPT_N_FPAR];
  double p[6];
  double r;

  err = 0;
  do
  {
    if(CheckParam(par,fpar) < 0)
    {
      err = 1;
      break;
    }
    if(fpar[OPT_FPAR_DIST] > 0)
    {
      if(MieCalc(par,fpar) < 0)
      {
        err = 1;
        break;
      }
    }
    if(fpar[OPT_FPAR_MIXR] > 0)
    {
      if(MixComp(par) < 0)
      {
        err = 1;
        break;
      }
      if(cnt_vb > 2) MiePrintout();
    }
  }
  while(0);
  if(err)
  {
    return -1;
  }

  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    // individual setup
    switch(sim_mode)
    {
      case SIM_MODTRAN_VTAU:
        if(fpar[OPT_FPAR_PROF] > 0)
        {
          if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
          {
            err = 1;
          }
        }
        break;
      case SIM_MODTRAN_VMIE:
        sim_vmie = par[OPT_PAR_VMIE];
        break;
    }
    // common setup
    sim_wscl = par[OPT_PAR_WSCL];
    sim_oscl = par[OPT_PAR_OSCL];
    sim_sbmp = SSR_SBMP;
    if(err)
    {
      return -1;
    }
    for(i=0; i<dsr_n_data; i++)
    {
      p[2] = dsr_th_sun[i];
      p[3] = dsr_ph_sun[i];
      p[4] = dsr_th_sun[i];
      p[5] = dsr_ph_sun[i];
      for(j=0; j<dsr_n_wlen[rng]; j++)
      {
        p[0] = dsr_wlen[rng][j];
        p[1] = (double)sim_nsam;
        err = 0;
        do
        {
          if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_SO) < 0)
          {
            err = 1;
            break;
          }
          if(RunModtran() < 0)
          {
            err = 1;
            break;
          }
          if(ReadPlt_S(p,&dsr_dsim[rng][i][j]) < 0)
          {
            err = 1;
            break;
          }
          if(ResFactor(dsr_ms720,dsr_nw[rng][j],sim_sbmp,dsr_wlen[rng][j],dsr_amas[i],&r) < 0)
          {
            err = 1;
            break;
          }
          dsr_dsim[rng][i][j] *= r;
        }
        while(0);
        if(err)
        {
          return -1;
        }
        if(cnt_sflg == 1) // calculate single scattering
        {
          sim_mult = 0;
          err = 0;
          do
          {
            if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_SO) < 0)
            {
              err = 1;
              break;
            }
            if(RunModtran() < 0)
            {
              err = 1;
              break;
            }
            if(ReadPlt_S(p,&dsr_dsgl[rng][i][j]) < 0)
            {
              err = 1;
              break;
            }
            if(ResFactor(dsr_ms720,dsr_nw[rng][j],sim_sbmp,dsr_wlen[rng][j],dsr_amas[i],&r) < 0)
            {
              err = 1;
              break;
            }
            dsr_dsgl[rng][i][j] *= r;
          }
          while(0);
          sim_mult = 1;
          if(err)
          {
            return -1;
          }
        }
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    // warning: correction factors are calculated for wavelength band
    // This function should be called with DIS==t
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<dsr_n_data; i++)
        {
          for(j=0; j<dsr_n_wlen[rng]; j++)
          {
            dsr_dsim[rng][i][j] = dsr_dsgl[rng][i][j]+(dsr_dsim[rng][i][j]-dsr_dsgl[rng][i][j])*dsr_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<dsr_n_data; i++)
        {
          for(j=0; j<dsr_n_wlen[rng]; j++)
          {
            dsr_dsim[rng][i][j] *= dsr_dfac[rng][i][j];
          }
        }
      }
    }
  }

  return 0;
}

// Scattered Solar Radiation at DSR wavelengths (Custom aerosol model, wavelength band)
int DSR_SSR_CUSTOM_B(int rng,double *par,int mode)
{
  // mode: 0=non-molecular abs. band, 1=H2O abs. band, 2=O3 abs. band
  // par[0] = vtau;
  // par[1] = vmie;
  // par[2] = wscl;
  // par[3] = oscl;
  // ...
  int i,j;
  int err;
  int fpar[OPT_N_FPAR];
  double p[6];

  err = 0;
  do
  {
    if(CheckParam(par,fpar) < 0)
    {
      err = 1;
      break;
    }
    if(fpar[OPT_FPAR_DIST] > 0)
    {
      if(MieCalc(par,fpar) < 0)
      {
        err = 1;
        break;
      }
    }
    if(fpar[OPT_FPAR_MIXR] > 0)
    {
      if(MixComp(par) < 0)
      {
        err = 1;
        break;
      }
      if(cnt_vb > 2) MiePrintout();
    }
  }
  while(0);
  if(err)
  {
    return -1;
  }

  if(fpar[OPT_FPAR_PMOD] > 0)
  {
    err = 0;
    // individual setup
    switch(sim_mode)
    {
      case SIM_MODTRAN_VTAU:
        if(fpar[OPT_FPAR_PROF] > 0)
        {
          if(PrfScaleAerosol(par[OPT_PAR_VTAU]) < 0)
          {
            err = 1;
          }
        }
        break;
      case SIM_MODTRAN_VMIE:
        sim_vmie = par[OPT_PAR_VMIE];
        break;
    }
    // common setup
    sim_wscl = par[OPT_PAR_WSCL];
    sim_oscl = par[OPT_PAR_OSCL];
    sim_sbmp = SSR_SBMP;
    if(err)
    {
      return -1;
    }
    for(i=0; i<dsr_n_data; i++)
    {
      p[2] = dsr_th_sun[i];
      p[3] = dsr_ph_sun[i];
      p[4] = dsr_th_sun[i];
      p[5] = dsr_ph_sun[i];
      for(j=0; j<dsr_n_wlen[rng]; j++)
      {
        p[0] = dsr_wlen[rng][j]-sim_xmgn;
        p[1] = dsr_wlen[rng][j]+sim_xmgn;
        err = 0;
        do
        {
          if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_BO) < 0)
          {
            err = 1;
            break;
          }
          if(RunModtran() < 0)
          {
            err = 1;
            break;
          }
          if(ReadPlt_B(p,1,&dsr_wlen[rng][j],&dsr_dsim[rng][i][j]) < 0)
          {
            err = 1;
            break;
          }
        }
        while(0);
        if(err)
        {
          return -1;
        }
        if(cnt_sflg == 1) // calculate single scattering
        {
          sim_mult = 0;
          err = 0;
          do
          {
            if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_BO) < 0)
            {
              err = 1;
              break;
            }
            if(RunModtran() < 0)
            {
              err = 1;
              break;
            }
            if(ReadPlt_B(p,1,&dsr_wlen[rng][j],&dsr_dsgl[rng][i][j]) < 0)
            {
              err = 1;
              break;
            }
          }
          while(0);
          sim_mult = 1;
          if(err)
          {
            return -1;
          }
        }
      }
    }
    if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
    {
      if(cnt_sflg == 1) // calculate single scattering
      {
        for(i=0; i<dsr_n_data; i++)
        {
          for(j=0; j<dsr_n_wlen[rng]; j++)
          {
            dsr_dsim[rng][i][j] = dsr_dsgl[rng][i][j]+(dsr_dsim[rng][i][j]-dsr_dsgl[rng][i][j])*dsr_dfac[rng][i][j];
          }
        }
      }
      else
      {
        for(i=0; i<dsr_n_data; i++)
        {
          for(j=0; j<dsr_n_wlen[rng]; j++)
          {
            dsr_dsim[rng][i][j] *= dsr_dfac[rng][i][j];
          }
        }
      }
    }
  }

  return 0;
}

// Direct Solar Radiation Fitting (Reference wavelength, MODTRAN aerosol model)
int DSR_FCN_MODTRAN_B0(double *par,double *val)
{
  return DSR_FCN_MODTRAN_B(0,par,val,FCN_MODE_AER);
}

// Direct Solar Radiation Fitting (non-molecular abs. band, MODTRAN aerosol model)
int DSR_FCN_MODTRAN_B1(double *par,double *val)
{
  return DSR_FCN_MODTRAN_B(1,par,val,FCN_MODE_AER);
}

// Direct Solar Radiation Fitting (H2O abs. band, MODTRAN aerosol model)
int DSR_FCN_MODTRAN_B2(double *par,double *val)
{
  return DSR_FCN_MODTRAN_B(2,par,val,FCN_MODE_H2O);
}

// Direct Solar Radiation Fitting (O3 abs. band, MODTRAN aerosol model)
int DSR_FCN_MODTRAN_B3(double *par,double *val)
{
  return DSR_FCN_MODTRAN_B(3,par,val,FCN_MODE_OZN);
}

// Solar Aureole Fitting (Reference wavelength, MODTRAN aerosol model)
int AUR_FCN_MODTRAN_S0(double *par,double *val)
{
  return AUR_FCN_MODTRAN_S(0,par,val,FCN_MODE_AER);
}

// Solar Aureole Fitting (non-molecular abs. band, MODTRAN aerosol model)
int AUR_FCN_MODTRAN_S1(double *par,double *val)
{
  return AUR_FCN_MODTRAN_S(1,par,val,FCN_MODE_AER);
}

// Solar Aureole Fitting (Reference wavelength, MODTRAN aerosol model)
int AUR_FCN_MODTRAN_B0(double *par,double *val)
{
  return AUR_FCN_MODTRAN_B(0,par,val,FCN_MODE_AER);
}

// Solar Aureole Fitting (non-molecular abs. band, MODTRAN aerosol model)
int AUR_FCN_MODTRAN_B1(double *par,double *val)
{
  return AUR_FCN_MODTRAN_B(1,par,val,FCN_MODE_AER);
}

// Scattered Solar Radiation Fitting (Reference wavelength, MODTRAN aerosol model)
int SSR_FCN_MODTRAN_S0(double *par,double *val)
{
  return SSR_FCN_MODTRAN_S(0,par,val,FCN_MODE_AER);
}

// Scattered Solar Radiation Fitting (non-molecular abs. band, MODTRAN aerosol model)
int SSR_FCN_MODTRAN_S1(double *par,double *val)
{
  return SSR_FCN_MODTRAN_S(1,par,val,FCN_MODE_AER);
}

// Scattered Solar Radiation Fitting (H2O abs. band, MODTRAN aerosol model)
int SSR_FCN_MODTRAN_S2(double *par,double *val)
{
  return SSR_FCN_MODTRAN_S(2,par,val,FCN_MODE_H2O);
}

// Scattered Solar Radiation Fitting (O3 abs. band, MODTRAN aerosol model)
int SSR_FCN_MODTRAN_S3(double *par,double *val)
{
  return SSR_FCN_MODTRAN_S(3,par,val,FCN_MODE_OZN);
}

// Scattered Solar Radiation Fitting (Reference wavelength, MODTRAN aerosol model)
int SSR_FCN_MODTRAN_B0(double *par,double *val)
{
  return SSR_FCN_MODTRAN_B(0,par,val,FCN_MODE_AER);
}

// Scattered Solar Radiation Fitting (non-molecular abs. band, MODTRAN aerosol model)
int SSR_FCN_MODTRAN_B1(double *par,double *val)
{
  return SSR_FCN_MODTRAN_B(1,par,val,FCN_MODE_AER);
}

// Scattered Solar Radiation Fitting (H2O abs. band, MODTRAN aerosol model)
int SSR_FCN_MODTRAN_B2(double *par,double *val)
{
  return SSR_FCN_MODTRAN_B(2,par,val,FCN_MODE_H2O);
}

// Scattered Solar Radiation Fitting (O3 abs. band, MODTRAN aerosol model)
int SSR_FCN_MODTRAN_B3(double *par,double *val)
{
  return SSR_FCN_MODTRAN_B(3,par,val,FCN_MODE_OZN);
}

// Direct Solar Radiation Fitting (Reference wavelength, Custom aerosol model)
int DSR_FCN_CUSTOM_B0(double *par,double *val)
{
  return DSR_FCN_CUSTOM_B(0,par,val,FCN_MODE_AER);
}

// Direct Solar Radiation Fitting (non-molecular abs. band, Custom aerosol model)
int DSR_FCN_CUSTOM_B1(double *par,double *val)
{
  return DSR_FCN_CUSTOM_B(1,par,val,FCN_MODE_AER);
}

// Direct Solar Radiation Fitting (H2O abs. band, Custom aerosol model)
int DSR_FCN_CUSTOM_B2(double *par,double *val)
{
  return DSR_FCN_CUSTOM_B(2,par,val,FCN_MODE_H2O);
}

// Direct Solar Radiation Fitting (O3 abs. band, Custom aerosol model)
int DSR_FCN_CUSTOM_B3(double *par,double *val)
{
  return DSR_FCN_CUSTOM_B(3,par,val,FCN_MODE_OZN);
}

// Solar Aureole Fitting (Reference wavelength, Custom aerosol model)
int AUR_FCN_CUSTOM_S0(double *par,double *val)
{
  return AUR_FCN_CUSTOM_S(0,par,val,FCN_MODE_AER);
}

// Solar Aureole Fitting (non-molecular abs. band, Custom aerosol model)
int AUR_FCN_CUSTOM_S1(double *par,double *val)
{
  return AUR_FCN_CUSTOM_S(1,par,val,FCN_MODE_AER);
}

// Solar Aureole Fitting (Reference wavelength, Custom aerosol model)
int AUR_FCN_CUSTOM_B0(double *par,double *val)
{
  return AUR_FCN_CUSTOM_B(0,par,val,FCN_MODE_AER);
}

// Solar Aureole Fitting (non-molecular abs. band, Custom aerosol model)
int AUR_FCN_CUSTOM_B1(double *par,double *val)
{
  return AUR_FCN_CUSTOM_B(1,par,val,FCN_MODE_AER);
}

// Scattered Solar Radiation Fitting (Reference wavelength, Custom aerosol model)
int SSR_FCN_CUSTOM_S0(double *par,double *val)
{
  return SSR_FCN_CUSTOM_S(0,par,val,FCN_MODE_AER);
}

// Scattered Solar Radiation Fitting (non-molecular abs. band, Custom aerosol model)
int SSR_FCN_CUSTOM_S1(double *par,double *val)
{
  return SSR_FCN_CUSTOM_S(1,par,val,FCN_MODE_AER);
}

// Scattered Solar Radiation Fitting (H2O abs. band, Custom aerosol model)
int SSR_FCN_CUSTOM_S2(double *par,double *val)
{
  return SSR_FCN_CUSTOM_S(2,par,val,FCN_MODE_H2O);
}

// Scattered Solar Radiation Fitting (O3 abs. band, Custom aerosol model)
int SSR_FCN_CUSTOM_S3(double *par,double *val)
{
  return SSR_FCN_CUSTOM_S(3,par,val,FCN_MODE_OZN);
}

// Scattered Solar Radiation Fitting (Reference wavelength, Custom aerosol model)
int SSR_FCN_CUSTOM_B0(double *par,double *val)
{
  return SSR_FCN_CUSTOM_B(0,par,val,FCN_MODE_AER);
}

// Scattered Solar Radiation Fitting (non-molecular abs. band, Custom aerosol model)
int SSR_FCN_CUSTOM_B1(double *par,double *val)
{
  return SSR_FCN_CUSTOM_B(1,par,val,FCN_MODE_AER);
}

// Scattered Solar Radiation Fitting (H2O abs. band, Custom aerosol model)
int SSR_FCN_CUSTOM_B2(double *par,double *val)
{
  return SSR_FCN_CUSTOM_B(2,par,val,FCN_MODE_H2O);
}

// Scattered Solar Radiation Fitting (O3 abs. band, Custom aerosol model)
int SSR_FCN_CUSTOM_B3(double *par,double *val)
{
  return SSR_FCN_CUSTOM_B(3,par,val,FCN_MODE_OZN);
}

// Scattered Solar Radiation at DSR wavelength (Reference wavelength, MODTRAN aerosol model)
int DSR_SSR_MODTRAN_S0(double *par)
{
  return DSR_SSR_MODTRAN_S(0,par,FCN_MODE_AER);
}

// Scattered Solar Radiation at DSR wavelength (non-molecular abs. band, MODTRAN aerosol model)
int DSR_SSR_MODTRAN_S1(double *par)
{
  return DSR_SSR_MODTRAN_S(1,par,FCN_MODE_AER);
}

// Scattered Solar Radiation at DSR wavelength (H2O abs. band, MODTRAN aerosol model)
int DSR_SSR_MODTRAN_S2(double *par)
{
  return DSR_SSR_MODTRAN_S(2,par,FCN_MODE_H2O);
}

// Scattered Solar Radiation at DSR wavelength (O3 abs. band, MODTRAN aerosol model)
int DSR_SSR_MODTRAN_S3(double *par)
{
  return DSR_SSR_MODTRAN_S(3,par,FCN_MODE_OZN);
}

// Scattered Solar Radiation at DSR wavelength (Reference wavelength, MODTRAN aerosol model)
int DSR_SSR_MODTRAN_B0(double *par)
{
  return DSR_SSR_MODTRAN_B(0,par,FCN_MODE_AER);
}

// Scattered Solar Radiation at DSR wavelength (non-molecular abs. band, MODTRAN aerosol model)
int DSR_SSR_MODTRAN_B1(double *par)
{
  return DSR_SSR_MODTRAN_B(1,par,FCN_MODE_AER);
}

// Scattered Solar Radiation at DSR wavelength (H2O abs. band, MODTRAN aerosol model)
int DSR_SSR_MODTRAN_B2(double *par)
{
  return DSR_SSR_MODTRAN_B(2,par,FCN_MODE_H2O);
}

// Scattered Solar Radiation at DSR wavelength (O3 abs. band, MODTRAN aerosol model)
int DSR_SSR_MODTRAN_B3(double *par)
{
  return DSR_SSR_MODTRAN_B(3,par,FCN_MODE_OZN);
}

// Scattered Solar Radiation at DSR wavelength (Reference wavelength, Custom aerosol model)
int DSR_SSR_CUSTOM_S0(double *par)
{
  return DSR_SSR_CUSTOM_S(0,par,FCN_MODE_AER);
}

// Scattered Solar Radiation at DSR wavelength (non-molecular abs. band, Custom aerosol model)
int DSR_SSR_CUSTOM_S1(double *par)
{
  return DSR_SSR_CUSTOM_S(1,par,FCN_MODE_AER);
}

// Scattered Solar Radiation at DSR wavelength (H2O abs. band, Custom aerosol model)
int DSR_SSR_CUSTOM_S2(double *par)
{
  return DSR_SSR_CUSTOM_S(2,par,FCN_MODE_H2O);
}

// Scattered Solar Radiation at DSR wavelength (O3 abs. band, Custom aerosol model)
int DSR_SSR_CUSTOM_S3(double *par)
{
  return DSR_SSR_CUSTOM_S(3,par,FCN_MODE_OZN);
}

// Scattered Solar Radiation at DSR wavelength (Reference wavelength, Custom aerosol model)
int DSR_SSR_CUSTOM_B0(double *par)
{
  return DSR_SSR_CUSTOM_B(0,par,FCN_MODE_AER);
}

// Scattered Solar Radiation at DSR wavelength (non-molecular abs. band, Custom aerosol model)
int DSR_SSR_CUSTOM_B1(double *par)
{
  return DSR_SSR_CUSTOM_B(1,par,FCN_MODE_AER);
}

// Scattered Solar Radiation at DSR wavelength (H2O abs. band, Custom aerosol model)
int DSR_SSR_CUSTOM_B2(double *par)
{
  return DSR_SSR_CUSTOM_B(2,par,FCN_MODE_H2O);
}

// Scattered Solar Radiation at DSR wavelength (O3 abs. band, Custom aerosol model)
int DSR_SSR_CUSTOM_B3(double *par)
{
  return DSR_SSR_CUSTOM_B(3,par,FCN_MODE_OZN);
}

int CRS_FCN_DSR_AUR(double *par,double *val)
{
  double vdsr,vaur;

  if(DSR_FCN_CUSTOM_B(1,par,&vdsr,FCN_MODE_AER) < 0) return -1;
  if(AUR_FCN_CUSTOM_B(1,par,&vaur,FCN_MODE_AER) < 0) return -1;
  if(vdsr > 1.0)
  {
    *val = vdsr*vaur;
  }
  else
  {
    *val = vaur;
  }

  return 0;
}

int CRS_FCN_DSR_SSR(double *par,double *val)
{
  double vdsr,vssr;

  if(DSR_FCN_CUSTOM_B(1,par,&vdsr,FCN_MODE_AER) < 0) return -1;
  if(SSR_FCN_CUSTOM_B(1,par,&vssr,FCN_MODE_AER) < 0) return -1;
  if(vdsr > 1.0)
  {
    *val = vdsr*vssr;
  }
  else
  {
    *val = vssr;
  }

  return 0;
}

int CRS_FCN_AUR_SSR(double *par,double *val)
{
  double vaur,vssr;

  if(AUR_FCN_CUSTOM_B(1,par,&vaur,FCN_MODE_AER) < 0) return -1;
  if(SSR_FCN_CUSTOM_B(1,par,&vssr,FCN_MODE_AER) < 0) return -1;
  *val = vaur*vssr;

  return 0;
}

int CRS_FCN_DSR_AUR_SSR(double *par,double *val)
{
  double vdsr,vaur,vssr;

  if(DSR_FCN_CUSTOM_B(1,par,&vdsr,FCN_MODE_AER) < 0) return -1;
  if(AUR_FCN_CUSTOM_B(1,par,&vaur,FCN_MODE_AER) < 0) return -1;
  if(SSR_FCN_CUSTOM_B(1,par,&vssr,FCN_MODE_AER) < 0) return -1;
  if(vdsr > 1.0)
  {
    *val = vdsr*vaur*vssr;
  }
  else
  {
    *val = vaur*vssr;
  }

  return 0;
}

int CheckParam(const double *par,int *fpar)
{
  // par[0]  = vtau;
  // par[1]  = vmie;
  // par[2]  = wscl;
  // par[3]  = oscl;
  // par[4]  = xscl;
  // par[5]  = yscl;
  // par[6]  = zscl;
  // par[7]  = n_comp;
  // par[8]  = lmod_com[0];
  // par[9]  = lsgm_com[0];
  // par[10] = log10(wcom_com[0]/wcom_com[1]);
  // par[11] = lmod_com[1];
  // par[12] = lsgm_com[1];
  // ...
  //  fpar[0] : Vertical profiles
  //  fpar[1] : MODTRAN parameters
  //  fpar[2] : Mixing ratios or size distributions
  //  fpar[3] : Size distributions
  //  fpar[4] : Comp1 size distribution
  //  fpar[5] : Comp2 size distribution
  // ...
  int i,j;
  int n_comp;
  char fnam[] = "CheckParam";

  for(i=0; i<OPT_N_FPAR; i++)
  {
    fpar[i] = 0;
  }
  if(sim_flag > 0) fpar[OPT_FPAR_PMOD] = 1;
  if(sim_flag > 1) fpar[OPT_FPAR_MIXR] = 1;
  if(sim_flag > 2) fpar[OPT_FPAR_DIST] = 1;
  // check parameters
  if(par[OPT_PAR_VTAU] != sim_pval[OPT_PAR_VTAU])
  {
    fpar[OPT_FPAR_PROF] = 1;
  }
  for(i=0; i<OPT_N_FPAR; i++)
  {
    if(par[i] != sim_pval[i])
    {
      fpar[OPT_FPAR_PMOD] = 1;
      break;
    }
  }
  n_comp = (int)(par[OPT_PAR_NCMP]+0.5);
  if(n_comp<0 || n_comp>MIE_MAXCOMP || fabs(par[OPT_PAR_NCMP]-(double)n_comp)>DELTA)
  {
    message(stderr,"%s: error, par[%d]=%18.10e, n_comp=%d\n",fnam,
                    OPT_PAR_NCMP,par[OPT_PAR_NCMP],n_comp);
    return -1;
  }
  if(par[OPT_PAR_NCMP]!=sim_pval[OPT_PAR_NCMP] || sim_flag>2)
  {
    fpar[OPT_FPAR_MIXR] = 1;
    fpar[OPT_FPAR_DIST] = 1;
    for(i=0; i<n_comp; i++)
    {
      fpar[OPT_FPAR_COMP+i] = 1;
    }
  }
  else
  {
    if(n_comp > 0)
    {
      if(par[OPT_PAR_LMD1] != sim_pval[OPT_PAR_LMD1]) // LMOD
      {
        fpar[OPT_FPAR_MIXR] = 1;
        fpar[OPT_FPAR_DIST] = 1;
        fpar[OPT_FPAR_COMP] = 1;
      }
      if(par[OPT_PAR_LSG1] != sim_pval[OPT_PAR_LSG1]) // LSGM
      {
        fpar[OPT_FPAR_MIXR] = 1;
        fpar[OPT_FPAR_DIST] = 1;
        fpar[OPT_FPAR_COMP] = 1;
      }
    }
    for(i=OPT_PAR_MXR2,j=1; j<n_comp; j++)
    {
      if(par[i] != sim_pval[i]) // MIXR
      {
        fpar[OPT_FPAR_MIXR] = 1;
      }
      i++;
      if(par[i] != sim_pval[i]) // LMOD
      {
        fpar[OPT_FPAR_MIXR] = 1;
        fpar[OPT_FPAR_DIST] = 1;
        fpar[OPT_FPAR_COMP+j] = 1;
      }
      i++;
      if(par[i] != sim_pval[i]) // LSGM
      {
        fpar[OPT_FPAR_MIXR] = 1;
        fpar[OPT_FPAR_DIST] = 1;
        fpar[OPT_FPAR_COMP+j] = 1;
      }
      i++;
    }
  }
  for(i=0; i<SIM_NPAR; i++)
  {
    sim_pval[i] = par[i];
  }

  if(cnt_db > 1)
  {
    for(i=0; i<OPT_N_FPAR; i++)
    {
      message(stderr,"%2d %d\n",i,fpar[i]);
    }
  }

  return 0;
}

int PrintValue(const double *par,const double *val,int npar,...)
{
  int i,n;
  int list[OPT_NPAR];
  va_list ap;
  char fnam[] = "PrintValue";

  if(npar<1 || npar>OPT_NPAR)
  {
    fprintf(stderr,"%s: #parameter out of range >>> %d\n",fnam,npar);
    return -1;
  }
  va_start(ap,npar);
  for(i=0; i<npar; i++)
  {
    n = va_arg(ap,int);
    if(n<0 || n>=OPT_NPAR)
    {
      fprintf(stderr,"%s: parameter# out of range >>> %d\n",fnam,n);
      return -1;
    }
    list[i] = n;
  }
  va_end(ap);
  message(stderr,"par:",n);
  for(i=0; i<npar; i++)
  {
    n = list[i];
    message(stderr,sim_form[n],par[n]);
  }
  message(stderr," val: %11.4e\n",*val);

  return 0;
}

int MieInit(void)
{
  int i,j,n;
  int n_comp;
  char fnam[] = "MieInit";

  mie_iref = -1;
  for(i=0; i<mie_n_wlen; i++)
  {
    if(fabs(mie_wlen[i]-mie_wlen_ref) < EPSILON)
    {
      mie_iref = i;
      break;
    }
  }
  if(mie_iref < 0)
  {
    fprintf(stderr,"%s: failed in finding reference wavelength %13.6e\n",fnam,mie_wlen_ref);
    return -1;
  }

  mie_wlen_um = (double *)malloc(mie_n_wlen*sizeof(double));
  if(mie_wlen_um == NULL)
  {
    fprintf(stderr,"%s: failed in allocating memory\n",fnam);
    return -1;
  }
  for(i=0; i<mie_n_wlen; i++)
  {
    mie_wlen_um[i] = mie_wlen[i]*1.0e-3;
  }
  mie_aext = (double *)malloc(mie_n_wlen*sizeof(double));
  mie_asca = (double *)malloc(mie_n_wlen*sizeof(double));
  mie_asym = (double *)malloc(mie_n_wlen*sizeof(double));
  mie_phas = (double *)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
  sim_phas = (double *)malloc(mie_n_wlen*sim_n_angl*sizeof(double));
  if(mie_aext==NULL || mie_asca==NULL || mie_asym==NULL || mie_phas==NULL || sim_phas==NULL)
  {
    fprintf(stderr,"%s: failed in allocating memory\n",fnam);
    return -1;
  }
  n_comp = MAX(OPT_N_COMP,MAX(tst_n_comp,crc_n_comp));
  for(n=0; n<n_comp; n++)
  {
    mie_refr_com[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_refi_com[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_aext_com[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_asca_com[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_asym_com[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_phs1_com[n] = (double *)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
    mie_phs2_com[n] = (double *)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
    mie_phas_com[n] = (double *)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
    if(mie_refr_com[n]==NULL || mie_refi_com[n]==NULL ||
       mie_aext_com[n]==NULL || mie_asca_com[n]==NULL || mie_asym_com[n]==NULL ||
       mie_phs1_com[n]==NULL || mie_phs2_com[n]==NULL || mie_phas_com[n]==NULL)
    {
      fprintf(stderr,"%s: failed in allocating memory\n",fnam);
      return -1;
    }
    mie_aext_tab[n] = NULL;
    mie_asca_tab[n] = NULL;
    mie_phs1_tab[n] = NULL;
    mie_phs2_tab[n] = NULL;
    mie_lrad_min[n] = NULL;
    mie_lrad_max[n] = NULL;
    mie_lrad_stp[n] = NULL;
    mie_lrad_lext_x[n] = NULL;
    mie_lrad_lext_y[n] = NULL;
    mie_lext_lrad_x[n] = NULL;
    mie_lext_lrad_y[n] = NULL;
    mie_lmod_leff_x[n] = NULL;
    mie_lmod_leff_y[n] = NULL;
    mie_leff_lmod_x[n] = NULL;
    mie_leff_lmod_y[n] = NULL;
  }

  mie_angl_rad = (double *)malloc(mie_n_angl*sizeof(double));
  mie_angl_sin = (double *)malloc(mie_n_angl*sizeof(double));
  mie_angl_cos = (double *)malloc(mie_n_angl*sizeof(double));
  mie_angl_dif = (double *)malloc(mie_n_angl*sizeof(double));
  if(mie_angl_rad==NULL || mie_angl_sin==NULL || mie_angl_cos==NULL || mie_angl_dif==NULL)
  {
    fprintf(stderr,"%s: failed in allocating memory\n",fnam);
    return -1;
  }
  for(j=0; j<mie_n_angl; j++)
  {
    mie_angl_rad[j] = mie_angl[j]*D_TO_R;
    mie_angl_sin[j] = sin(mie_angl_rad[j]);
    mie_angl_cos[j] = cos(mie_angl_rad[j]);
  }
  for(j=1; j<mie_n_angl; j++)
  {
    mie_angl_dif[j] = PI*fabs(mie_angl_rad[j]-mie_angl_rad[j-1]);
  }

  return 0;
}

int MieTable(int n1,int n2)
{
  // IMPORTANT: This function MUST be called before MieCalc().
  //            (after setting mie_refr_com & mie_refi_com.)
  int i,j,k,m,n;
  int m1,m2,mt;
  int err = 0;
  int fmie,ftmx;
  double x,r,l;
  double s1,s2;
  int ANYANG,PERFCT,PRNT[2];
  int IPOLZN,MOMDIM,NUMANG,NMOM;
  double GQSC,MIMCUT,PMOM[2][2],QEXT,QSCA,SPIKE;
  double _Complex CREFIN,SFORW,SBACK,*S1,*S2,TFORW[2],TBACK[2];
  int iflg,ierr,ndistr,npnax,nkmax;
  int np,npna,ndgs,l1max[TMX_NAXMAX];
  double rat,axmax,r1,r2,b,gam,eps;
  double lam,mrr,mri,ddelt;
  double *reff,*veff,*cext,*csca,*walb,*gsym;
  double *alp1,*alp2,*alp3,*alp4,*bet1,*bet2;
  double *fmat,*angl,*f11,*f12,*a1,*b1;
  char fnam[] = "MieTable";

  // initialization
  fmie = ftmx = 0;
  for(n=n1; n<=n2; n++)
  {
    if(cnt_tmx[n])
    {
      ftmx = 1;
    }
    else
    {
      fmie = 1;
    }
  }
  if(fmie)
  {
    S1 = (double _Complex *)malloc(mie_n_angl*sizeof(double _Complex));
    S2 = (double _Complex *)malloc(mie_n_angl*sizeof(double _Complex));
    if(S1==NULL || S2==NULL)
    {
      fprintf(stderr,"%s: failed in allocating memory\n",fnam);
      return -1;
    }
    PERFCT = MIE_FALSE;
    MIMCUT = DELTA;
    ANYANG = MIE_FALSE;
    NUMANG = mie_n_angl;
    PRNT[0] = (cnt_db>1?MIE_TRUE:MIE_FALSE);
    PRNT[1] = (cnt_db>0?MIE_TRUE:MIE_FALSE);
    NMOM = 0;
    IPOLZN = 0;
    MOMDIM = 1;
  }
  if(ftmx)
  {
    reff = (double*)malloc(TMX_NAXMAX*sizeof(double));
    veff = (double*)malloc(TMX_NAXMAX*sizeof(double));
    cext = (double*)malloc(TMX_NAXMAX*sizeof(double));
    csca = (double*)malloc(TMX_NAXMAX*sizeof(double));
    walb = (double*)malloc(TMX_NAXMAX*sizeof(double));
    gsym = (double*)malloc(TMX_NAXMAX*sizeof(double));
    alp1 = (double*)malloc(TMX_NALMAX*sizeof(double));
    alp2 = (double*)malloc(TMX_NALMAX*sizeof(double));
    alp3 = (double*)malloc(TMX_NALMAX*sizeof(double));
    alp4 = (double*)malloc(TMX_NALMAX*sizeof(double));
    bet1 = (double*)malloc(TMX_NALMAX*sizeof(double));
    bet2 = (double*)malloc(TMX_NALMAX*sizeof(double));
    fmat = (double*)malloc(TMX_NARMAX*sizeof(double));
    angl = (double*)malloc(TMX_NAMAX*sizeof(double));
    f11  = (double*)malloc(TMX_NAMAX*sizeof(double));
    f12  = (double*)malloc(TMX_NAMAX*sizeof(double));
    a1   = (double*)malloc(MIE_MAXDATA*sizeof(double));
    b1   = (double*)malloc(MIE_MAXDATA*sizeof(double));
    if(reff==NULL || veff==NULL || cext==NULL || csca==NULL || walb==NULL || gsym==NULL ||
       alp1==NULL || alp2==NULL || alp3==NULL || alp4==NULL || bet1==NULL || bet2==NULL ||
       fmat==NULL || angl==NULL ||  f11==NULL ||  f12==NULL ||   a1==NULL ||   b1==NULL)
    {
      fprintf(stderr,"%s: failed in allocating memory\n",fnam);
      return -1;
    }
    for(j=0; j<tmx_nang; j++)
    {
      angl[j] = (180.0/(tmx_nang-1)*j);
    }
  }

  for(n=n1; n<=n2; n++)
  {
    // allocate memory
    if(mie_aext_tab[n] != NULL) free(mie_aext_tab[n]);
    if(mie_asca_tab[n] != NULL) free(mie_asca_tab[n]);
    if(mie_phs1_tab[n] != NULL) free(mie_phs1_tab[n]);
    if(mie_phs2_tab[n] != NULL) free(mie_phs2_tab[n]);
    mie_aext_tab[n] = (double *)calloc(mie_n_wlen*mie_n_step,sizeof(double));
    mie_asca_tab[n] = (double *)calloc(mie_n_wlen*mie_n_step,sizeof(double));
    mie_phs1_tab[n] = (double *)calloc(mie_n_wlen*mie_n_angl*mie_n_step,sizeof(double));
    mie_phs2_tab[n] = (double *)calloc(mie_n_wlen*mie_n_angl*mie_n_step,sizeof(double));
    if(mie_aext_tab[n]==NULL || mie_asca_tab[n]==NULL ||
       mie_phs1_tab[n]==NULL || mie_phs2_tab[n]==NULL)
    {
      fprintf(stderr,"%s: error in allocating memory\n",fnam);
      err = 1;
      break;
    }
    if(mie_lrad_min[n] != NULL) free(mie_lrad_min[n]);
    if(mie_lrad_max[n] != NULL) free(mie_lrad_max[n]);
    if(mie_lrad_stp[n] != NULL) free(mie_lrad_stp[n]);
    mie_lrad_min[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_lrad_max[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    mie_lrad_stp[n] = (double *)malloc(mie_n_wlen*sizeof(double));
    if(mie_lrad_min[n]==NULL || mie_lrad_max[n]==NULL || mie_lrad_stp[n]==NULL)
    {
      fprintf(stderr,"%s: error in allocating memory\n",fnam);
      err = 1;
      break;
    }
    if(mie_lrad_lext_x[n] != NULL) free(mie_lrad_lext_x[n]);
    if(mie_lrad_lext_y[n] != NULL) free(mie_lrad_lext_y[n]);
    if(mie_lext_lrad_x[n] != NULL) free(mie_lext_lrad_x[n]);
    if(mie_lext_lrad_y[n] != NULL) free(mie_lext_lrad_y[n]);
    mie_lrad_lext_x[n] = (double *)malloc(mie_n_step*sizeof(double));
    mie_lrad_lext_y[n] = (double *)malloc(mie_n_step*sizeof(double));
    mie_lext_lrad_x[n] = (double *)malloc(mie_n_step*sizeof(double));
    mie_lext_lrad_y[n] = (double *)malloc(mie_n_step*sizeof(double));
    if(mie_lrad_lext_x[n]==NULL || mie_lrad_lext_y[n]==NULL ||
       mie_lext_lrad_x[n]==NULL || mie_lext_lrad_y[n]==NULL)
    {
      fprintf(stderr,"%s: error in allocating memory\n",fnam);
      err = 1;
      break;
    }

    // Mie calculation
    for(i=0; i<mie_n_wlen; i++)
    {
      if(isnan(mie_rmin))
      {
        mie_lrad_min[n][i] = MIE_LRAD_MIN;
      }
      else
      {
        mie_lrad_min[n][i] = log10(mie_rmin);
        if(mie_lrad_min[n][i] < MIE_LRAD_MIN)
        {
          mie_lrad_min[n][i] = MIE_LRAD_MIN;
        }
      }
      if(isnan(mie_rmax))
      {
        mie_lrad_max[n][i] = MIE_LRAD_MAX;
      }
      else
      {
        mie_lrad_max[n][i] = log10(mie_rmax);
        if(mie_lrad_max[n][i] > MIE_LRAD_MAX)
        {
          mie_lrad_max[n][i] = MIE_LRAD_MAX;
        }
      }
      if(cnt_tmx[n])
      {
        l = log10(TMX_RMAX);
      }
      else
      {
        l = log10(MIE_XMAX*mie_wlen_um[i]/PI2);
      }
      if(mie_lrad_max[n][i] > l)
      {
        mie_lrad_max[n][i] = l;
      }
      mie_lrad_stp[n][i] = (mie_lrad_max[n][i]-mie_lrad_min[n][i])/mie_n_step;
      for(m=0; m<mie_n_step; m++)
      {
        l = mie_lrad_min[n][i]+mie_lrad_stp[n][i]*m;
        r = pow(10.0,l);
        x = PI2*r/(mie_wlen_um[i]);
        if(cnt_tmx[n])
        {
          iflg   = cnt_vb;
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
          eps    = tmx_reps[n];
          np     = tmx_shap[n];
          lam    = mie_wlen_um[i];
          mrr    = mie_refr_com[n][i];
          mri    = mie_refi_com[n][i];
          ddelt  = tmx_delt; // this will be modified by tmatrix_()!
          npna   = tmx_nang;
          ndgs   = tmx_ndgs;
          tmatrix_(&iflg,&ierr,&rat,&ndistr,&axmax,&npnax,&r1,&r2,&b,&gam,
                   &nkmax,&eps,&np,&lam,&mrr,&mri,&ddelt,&npna,&ndgs,
                   reff,veff,cext,csca,walb,gsym,
                   l1max,alp1,alp2,alp3,alp4,bet1,bet2,fmat);
          if(ierr)
          {
            fprintf(stderr,"%s: warning, ierr=%d, r=%13.6e, lam=%13.6f, x=%13.6e, eps=%13.6e\n",
                            fnam,ierr,r,lam,x,eps);
            break;
          }
          k = mie_n_step*i+m;
          mie_aext_tab[n][k] = cext[0];
          mie_asca_tab[n][k] = csca[0];
          if(i == mie_iref)
          {
            mie_lrad_lext_x[n][m] = l;
            mie_lrad_lext_y[n][m] = log10(cext[0]);
          }
          for(j=0; j<tmx_nang; j++)
          {
            f11[j] = fmat[j*6+0];
            f12[j] = fmat[j*6+4];
          }
          SamplingE(tmx_nang,angl,f11,mie_n_angl,mie_angl,a1,1.0);
          SamplingE(tmx_nang,angl,f12,mie_n_angl,mie_angl,b1,1.0);
          for(j=0; j<mie_n_angl; j++)
          {
            k = (mie_n_angl*mie_n_step)*i+mie_n_step*j+m;
            mie_phs1_tab[n][k] = (a1[j]-b1[j])*csca[0];
            mie_phs2_tab[n][k] = (a1[j]+b1[j])*csca[0];
          }
        }
        else
        {
          CREFIN = mie_refr_com[n][i]-fabs(mie_refi_com[n][i])*I;
          miev0_(&x,&CREFIN,&PERFCT,&MIMCUT,&ANYANG,&NUMANG,mie_angl_cos,
                 &NMOM,&IPOLZN,&MOMDIM,PRNT,&QEXT,&QSCA,&GQSC,
                 PMOM,&SFORW,&SBACK,S1,S2,TFORW,TBACK,&SPIKE);
          k = mie_n_step*i+m;
          mie_aext_tab[n][k] = QEXT*PI*r*r;
          mie_asca_tab[n][k] = QSCA*PI*r*r;
          if(i == mie_iref)
          {
            mie_lrad_lext_x[n][m] = l;
            mie_lrad_lext_y[n][m] = log10(QEXT*PI*r*r);
          }
          for(j=0; j<mie_n_angl; j++)
          {
            k = (mie_n_angl*mie_n_step)*i+mie_n_step*j+m;
            s1 = cabs(S1[j]);
            s2 = cabs(S2[j]);
            mie_phs1_tab[n][k] = s1*s1;
            mie_phs2_tab[n][k] = s2*s2;
          }
        }
      }
    }
    if(err) break;
  }
  if(fmie)
  {
    free(S1);
    free(S2);
  }
  if(ftmx)
  {
    free(reff);
    free(veff);
    free(cext);
    free(csca);
    free(walb);
    free(gsym);
    free(alp1);
    free(alp2);
    free(alp3);
    free(alp4);
    free(bet1);
    free(bet2);
    free(fmat);
    free(angl);
    free(f11);
    free(f12);
    free(a1);
    free(b1);
  }
  if(err)
  {
    return -1;
  }
  if(cnt_db > 1)
  {
    for(n=n1; n<=n2; n++)
    {
      for(m=0; m<mie_n_step; m++)
      {
        fprintf(stderr,"%13.6e %13.6e %d %d\n",mie_lrad_lext_x[n][m],mie_lrad_lext_y[n][m],n,1);
      }
    }
  }

  // prepair inversion table
  for(n=n1; n<=n2; n++)
  {
    if(MieSmoothing(mie_lrad_lext_x[n],mie_lrad_lext_y[n],mie_lrad_stp[n][mie_iref],0.10,0.05,mie_n_step) < 0) return -1;
    mie_lext_min[n] = mie_lrad_lext_y[n][0];
    mie_lext_max[n] = mie_lrad_lext_y[n][mie_n_step-1];
    mie_lext_stp[n] = (mie_lext_max[n]-mie_lext_min[n])/mie_n_step;
    for(m=0; m<mie_n_step; m++)
    {
      mie_lext_lrad_x[n][m] = mie_lext_min[n]+mie_lext_stp[n]*m;
      m1 = m2 = -1;
      for(mt=1; mt<mie_n_step; mt++)
      {
        if(mie_lrad_lext_y[n][mt] > mie_lext_lrad_x[n][m])
        {
          m1 = mt-1;
          break;
        }
      }
      for(mt=mie_n_step-2; mt>=0; mt--)
      {
        if(mie_lrad_lext_y[n][mt] < mie_lext_lrad_x[n][m])
        {
          m2 = mt+1;
          break;
        }
      }
      if(m1 < 0)
      {
        m1 = mie_n_step-2;
        m2 = mie_n_step-1;
      } else
      if(m2 < 0)
      {
        m1 = 0;
        m2 = 1;
      }
      mie_lext_lrad_y[n][m] = INTERP(mie_lext_lrad_x[n][m],mie_lrad_lext_y[n][m1],mie_lrad_lext_y[n][m2],
                                                           mie_lrad_lext_x[n][m1],mie_lrad_lext_x[n][m2]);
    }
  }
  if(cnt_db > 1)
  {
    for(n=n1; n<=n2; n++)
    {
      for(m=0; m<mie_n_step; m++)
      {
        fprintf(stderr,"%13.6e %13.6e %d %d\n",mie_lrad_lext_x[n][m],mie_lrad_lext_y[n][m],n,2);
        fprintf(stderr,"%13.6e %13.6e %d %d\n",mie_lext_lrad_y[n][m],mie_lext_lrad_x[n][m],n,3);
      }
    }
  }

  return 0;
}

int MieReff(int n1,int n2,double lsgm)
{
  int j,k,m,n;
  int m1,m2,mt;
  int err = 0;
  double x,y,r,l;
  double sy;
  double norm,slop;
  double lmin,lmax;
  double aext,lext;
  double lmod;
  double epsilon = 2.0e-1;
  char fnam[] = "MieReff";

  for(n=n1; n<=n2; n++)
  {
    // allocate memory
    if(mie_lmod_leff_x[n] != NULL) free(mie_lmod_leff_x[n]);
    if(mie_lmod_leff_y[n] != NULL) free(mie_lmod_leff_y[n]);
    if(mie_leff_lmod_x[n] != NULL) free(mie_leff_lmod_x[n]);
    if(mie_leff_lmod_y[n] != NULL) free(mie_leff_lmod_y[n]);
    mie_lmod_leff_x[n] = (double *)malloc(mie_n_step*sizeof(double));
    mie_lmod_leff_y[n] = (double *)malloc(mie_n_step*sizeof(double));
    mie_leff_lmod_x[n] = (double *)malloc(mie_n_step*sizeof(double));
    mie_leff_lmod_y[n] = (double *)malloc(mie_n_step*sizeof(double));
    if(mie_lmod_leff_x[n]==NULL || mie_lmod_leff_y[n]==NULL ||
       mie_leff_lmod_x[n]==NULL || mie_leff_lmod_y[n]==NULL)
    {
      fprintf(stderr,"%s: error in allocating memory\n",fnam);
      err = 1;
      break;
    }
    // integrate over size
    slop = -0.5/(lsgm*lsgm);
    mie_lmod_min[n] = MIE_LMOD_MIN;
    mie_lmod_max[n] = MIE_LMOD_MAX;
    mie_lmod_stp[n] = (mie_lmod_max[n]-mie_lmod_min[n])/mie_n_step;
    norm = mie_lrad_stp[n][mie_iref]/(sqrt(PI2)*lsgm);
    for(j=0; j<mie_n_step; j++)
    {
      lmod = mie_lmod_min[n]+mie_lmod_stp[n]*j;
      lmin = lmod-lsgm*mie_wsgm;
      lmax = lmod+lsgm*mie_wsgm;
      aext = 0.0;
      for(m=0,sy=0.0; m<mie_n_step; m++)
      {
        l = mie_lrad_min[n][mie_iref]+mie_lrad_stp[n][mie_iref]*m;
        if(l < lmin) continue;
        if(l > lmax) break;
        r = pow(10.0,l);
        x = PI2*r/(mie_wlen_um[mie_iref]);
        y = norm*exp(slop*(l-lmod)*(l-lmod));
        k = mie_n_step*mie_iref+m;
        aext += mie_aext_tab[n][k]*y;
        sy += y;
      }
      if(fabs(sy-1.0) > epsilon)
      {
        fprintf(stderr,"%s: warning, sy=%13.6e\n",fnam,sy);
      }
      aext /= sy;
      lext = log10(aext);
      mie_lmod_leff_x[n][j] = lmod;
      m1 = (int)((lext-mie_lext_min[n])/mie_lext_stp[n]);
      if(m1 < 0) m1 = 0;
      m2 = m1+1;
      if(m2 >= mie_n_step)
      {
        m2 = mie_n_step-1;
        m1 = m2-1;
      }
      mie_lmod_leff_y[n][j] = INTERP(lext,mie_lext_lrad_x[n][m1],mie_lext_lrad_x[n][m2],
                                          mie_lext_lrad_y[n][m1],mie_lext_lrad_y[n][m2]);
    }
  }
  if(err)
  {
    return -1;
  }
  if(cnt_db > 1)
  {
    for(n=n1; n<=n2; n++)
    {
      for(m=0; m<mie_n_step; m++)
      {
        fprintf(stderr,"%13.6e %13.6e %d %d\n",mie_lmod_leff_x[n][m],mie_lmod_leff_y[n][m],n,4);
      }
    }
  }

  // prepair inversion table
  for(n=n1; n<=n2; n++)
  {
    if(MieSmoothing(mie_lmod_leff_x[n],mie_lmod_leff_y[n],mie_lmod_stp[n],0.10,0.05,mie_n_step) < 0) return -1;
    mie_leff_min[n] = mie_lmod_leff_y[n][0];
    mie_leff_max[n] = mie_lmod_leff_y[n][mie_n_step-1];
    mie_leff_stp[n] = (mie_leff_max[n]-mie_leff_min[n])/mie_n_step;
    for(m=0; m<mie_n_step; m++)
    {
      mie_leff_lmod_x[n][m] = mie_leff_min[n]+mie_leff_stp[n]*m;
      m1 = m2 = -1;
      for(mt=1; mt<mie_n_step; mt++)
      {
        if(mie_lmod_leff_y[n][mt] > mie_leff_lmod_x[n][m])
        {
          m1 = mt-1;
          break;
        }
      }
      for(mt=mie_n_step-2; mt>=0; mt--)
      {
        if(mie_lmod_leff_y[n][mt] < mie_leff_lmod_x[n][m])
        {
          m2 = mt+1;
          break;
        }
      }
      if(m1 < 0)
      {
        m1 = mie_n_step-2;
        m2 = mie_n_step-1;
      } else
      if(m2 < 0)
      {
        m1 = 0;
        m2 = 1;
      }
      mie_leff_lmod_y[n][m] = INTERP(mie_leff_lmod_x[n][m],mie_lmod_leff_y[n][m1],mie_lmod_leff_y[n][m2],
                                                           mie_lmod_leff_x[n][m1],mie_lmod_leff_x[n][m2]);
    }
  }
  if(cnt_db > 1)
  {
    for(n=n1; n<=n2; n++)
    {
      for(m=0; m<mie_n_step; m++)
      {
        fprintf(stderr,"%13.6e %13.6e %d %d\n",mie_lmod_leff_x[n][m],mie_lmod_leff_y[n][m],n,5);
        fprintf(stderr,"%13.6e %13.6e %d %d\n",mie_leff_lmod_y[n][m],mie_leff_lmod_x[n][m],n,6);
      }
    }
  }

  return 0;
}

int MieCalc(const double *par,const int *fpar)
{
  // par[7]  = n_comp;
  // par[8]  = lmod_com[0];
  // par[9]  = lsgm_com[0];
  // par[10] = log10(wcom_com[0]/wcom_com[1]);
  // par[11] = lmod_com[1];
  // par[12] = lsgm_com[1];
  // ...
  int i,j,k,m,n,p;
  double x,y,r,l;
  double sy,sp;
  double norm,slop;
  double lmin,lmax;
  double epsilon = 2.0e-1;
  char fnam[] = "MieCalc";

  // set parameters
  mie_n_comp = (int)(par[OPT_PAR_NCMP]+0.5);
  mie_lmod_com[0] = eff2mod(0,par[OPT_PAR_LMD1]);
  mie_lsgm_com[0] = par[OPT_PAR_LSG1];
  for(i=OPT_PAR_LMD2,j=1; j<mie_n_comp; i+=3,j++)
  {
    mie_lmod_com[j] = eff2mod(j,par[i]);
    mie_lsgm_com[j] = par[i+1];
  }

  // initialization
  for(n=0; n<mie_n_comp; n++)
  {
    if(fpar[OPT_FPAR_COMP+n] == 0) continue;
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

  // integrate over size
  for(n=0; n<mie_n_comp; n++)
  {
    if(fpar[OPT_FPAR_COMP+n] == 0) continue;
    slop = -0.5/(mie_lsgm_com[n]*mie_lsgm_com[n]);
    lmin = mie_lmod_com[n]-mie_lsgm_com[n]*mie_wsgm;
    lmax = mie_lmod_com[n]+mie_lsgm_com[n]*mie_wsgm;
    for(i=0; i<mie_n_wlen; i++)
    {
      norm = mie_lrad_stp[n][i]/(sqrt(PI2)*mie_lsgm_com[n]);
      for(m=0,sy=0.0; m<mie_n_step; m++)
      {
        l = mie_lrad_min[n][i]+mie_lrad_stp[n][i]*m;
        if(l < lmin) continue;
        if(l > lmax) break;
        r = pow(10.0,l);
        x = PI2*r/(mie_wlen_um[i]);
        y = norm*exp(slop*(l-mie_lmod_com[n])*(l-mie_lmod_com[n]));
        k = mie_n_step*i+m;
        mie_aext_com[n][i] += mie_aext_tab[n][k]*y;
        mie_asca_com[n][i] += mie_asca_tab[n][k]*y;
        for(j=0; j<mie_n_angl; j++)
        {
          k = mie_n_angl*i+j;
          p = (mie_n_angl*mie_n_step)*i+mie_n_step*j+m;
          mie_phs1_com[n][k] += mie_phs1_tab[n][p]*y;
          mie_phs2_com[n][k] += mie_phs2_tab[n][p]*y;
        }
        sy += y;
      }
      if(fabs(sy-1.0) > epsilon)
      {
        fprintf(stderr,"%s: warning, sy=%13.6e\n",fnam,sy);
      }
      mie_aext_com[n][i] /= sy;
      mie_asca_com[n][i] /= sy;
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        mie_phs1_com[n][k] /= sy;
        mie_phs2_com[n][k] /= sy;
        mie_phas_com[n][k] = 0.5*(mie_phs1_com[n][k]+mie_phs2_com[n][k]);
      }
      sp = 0.0;
      for(j=1; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        sp += (mie_phas_com[n][k-1]*mie_angl_sin[j-1]+mie_phas_com[n][k]*mie_angl_sin[j])*mie_angl_dif[j];
      }
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        mie_phas_com[n][k] /= sp;
      }
      for(j=1; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        mie_asym_com[n][i] += (mie_phas_com[n][k-1]*mie_angl_sin[j-1]*mie_angl_cos[j-1]+
                               mie_phas_com[n][k]*mie_angl_sin[j]*mie_angl_cos[j])*mie_angl_dif[j];
      }
    }
  }

  return 0;
}

int MixComp(const double *par)
{
  // par[7]  = n_comp;
  // par[8]  = lmod_com[0];
  // par[9]  = lsgm_com[0];
  // par[10] = log10(wcom_com[0]/wcom_com[1]);
  // par[11] = lmod_com[1];
  // par[12] = lsgm_com[1];
  // ...
  int i,j,k,n;
  double r;
  double we,ws;
  double sw,se,ss,sa,sp;

  // set parameters
  mie_n_comp = (int)(par[OPT_PAR_NCMP]+0.5);
  for(r=1.0,i=OPT_PAR_MXR2,j=1; j<mie_n_comp; i+=3,j++)
  {
    mie_wcom_com[j] = 1.0/pow(10.0,par[i]);
    r += mie_wcom_com[j];
  }
  mie_wcom_com[0] = 1.0/r;
  for(j=1; j<mie_n_comp; j++)
  {
    mie_wcom_com[j] *= mie_wcom_com[0];
  }

  // calculate number mixing ratio
  for(r=0.0,n=0; n<mie_n_comp; n++)
  {
    mie_mixr_com[n] = mie_wcom_com[n]/mie_aext_com[n][mie_iref];
    r += mie_mixr_com[n];
  }
  for(n=0; n<mie_n_comp; n++)
  {
    mie_mixr_com[n] /= r;
    if(cnt_vb > 2)
    {
      message(stderr,"Component#   : %2d\n",n);
      message(stderr,"Mixing ratio : %22.14e (%22.14e)\n",mie_mixr_com[n],mie_wcom_com[n]);
    }
  }

  for(i=0; i<mie_n_wlen; i++)
  {
    sw = se = ss = sa = 0.0;
    for(n=0; n<mie_n_comp; n++)
    {
      we  = mie_mixr_com[n]*mie_aext_com[n][i];
      ws  = mie_mixr_com[n]*mie_asca_com[n][i];
      sw += mie_mixr_com[n];
      se += we;
      ss += ws;
      sa += ws*mie_asym_com[n][i];
    }
    mie_aext[i] = se/sw;
    mie_asca[i] = ss/sw;
    mie_asym[i] = sa/ss;
    for(j=0; j<mie_n_angl; j++)
    {
      k = mie_n_angl*i+j;
      sp = 0.0;
      for(n=0; n<mie_n_comp; n++)
      {
        ws  = mie_mixr_com[n]*mie_asca_com[n][i];
        sp += ws*mie_phas_com[n][k];
      }
      mie_phas[k] = sp/ss;
    }
  }

  for(i=0; i<mie_n_wlen; i++)
  {
    sp = 0.0;
    for(j=1; j<mie_n_angl; j++)
    {
      k = mie_n_angl*i+j;
      sp += (mie_phas[k-1]*mie_angl_sin[j-1]+mie_phas[k]*mie_angl_sin[j])*mie_angl_dif[j];
    }
    for(j=0; j<mie_n_angl; j++)
    {
      k = mie_n_angl*i+j;
      mie_phas[k] /= sp;
    }
  }

  if(Interp2D(mie_wlen,mie_angl,mie_phas,mie_n_wlen,mie_n_angl,
              mie_wlen,sim_angl,sim_phas,mie_n_wlen,sim_n_angl,0) < 0)
  {
    return -1;
  }

  return 0;
}

int MieSmoothing(double *x,double *y,double xstp,double xmgn,double xwid,int nstp)
{
  int m,n;
  int m1,m2,m3,m4,m5,m6;
  int n1,n2;
  int nmgn,nwid;
  double d1,d2;
  double vmin,vmax;
  double yorg[MIE_NMAX];
  double ynew[MIE_NMAX];
  char fnam[] = "MieSmoothing";

  if(nstp<=1 || nstp>MIE_NMAX)
  {
    fprintf(stderr,"%s: error, nstp=%d\n",fnam,nstp);
    return -1;
  }
  else
  {
    nmgn = (int)(xmgn/xstp);
    nwid = (int)(xwid/xstp);
    if(nmgn<=1 || nwid<=1)
    {
      fprintf(stderr,"%s: error, nmgn=%d, nwid=%d\n",fnam,nmgn,nwid);
      return -1;
    }
  }

  // initialize
  for(m=0; m<nstp; m++)
  {
    ynew[m] = y[m];
  }
  if(cnt_db)
  {
    fprintf(stderr,"%s: before smoothing (%d)\n",fnam,nstp);
    for(m=0; m<nstp; m++)
    {
      fprintf(stderr,"%13.6e %13.6e\n",x[m],y[m]);
    }
    fprintf(stderr,"-----\n");
  }

  // smooth the function and make it increase monotonically
  while(1)
  {
    for(m=0; m<nstp; m++)
    {
      yorg[m] = ynew[m];
    }
    m1 = m2 = -1;
    // find inflection points
    vmax = -1.0e100;
    for(m=0; m<nstp; m++)
    {
      if(yorg[m] < vmax)
      {
        m1 = m;
        break;
      }
      else
      {
        vmax = yorg[m];
      }
    }
    vmin = 1.0e100;
    for(m=nstp-1; m>=0; m--)
    {
      if(yorg[m] > vmin)
      {
        m2 = m;
        break;
      }
      else
      {
        vmin = yorg[m];
      }
    }
    if(m1<0 && m2<0)
    {
      break;
    } else
    if(m1<0 || m2<0)
    {
      fprintf(stderr,"%s: error, m1=%d, m2=%d\n",fnam,m1,m2);
      return -1;
    } else
    if(cnt_db > 1)
    {
      fprintf(stderr,"%s: m1=%d(%13.6e) m2=%d(%13.6e)\n",fnam,m1,x[m1],m2,x[m2]);
    }
    // make monotonically increasing functions
    m3 = m1-nmgn;
    m4 = m2+nmgn;
    if(m3 < 0) m3 = 0;
    if(m4 >= nstp) m4 = nstp-1;
    for(m=m3; m<=m4; m++)
    {
      n1 = m-nwid;
      n2 = m+nwid;
      if(n1 < 0) n1 = 0;
      if(n2 >= nstp) n2 = nstp-1;
      if(n2-n1 < 0)
      {
        fprintf(stderr,"%s: error, n1=%d, n2=%d\n",fnam,n1,n2);
        return -1;
      }
      // calculate moving averages
      ynew[m] = 0.0;
      for(n=n1; n<=n2; n++)
      {
        ynew[m] += yorg[n];
      }
      ynew[m] /= (n2-n1+1);
    }
    // find the connecting points
    m5 = m6 = -1;
    d1 = yorg[m1]-ynew[m1];
    for(m=m1; m>=m3; m--)
    {
      d2 = yorg[m]-ynew[m];
      if(d1*d2 <= 0.0)
      {
        m5 = m;
        break;
      }
    }
    d1 = yorg[m2]-ynew[m2];
    for(m=m2; m<=m4; m++)
    {
      d2 = yorg[m]-ynew[m];
      if(d1*d2 <= 0.0)
      {
        m6 = m;
        break;
      }
    }
    if(m5 >= 0)
    {
      for(m=m3; m<=m5; m++)
      {
        ynew[m] = yorg[m];
      }
    }
    if(m6 >= 0)
    {
      for(m=m6; m<=m4; m++)
      {
        ynew[m] = yorg[m];
      }
    }
  }

  // reflect the results
  for(m=0; m<nstp; m++)
  {
    y[m] = ynew[m];
  }
  if(cnt_db)
  {
    fprintf(stderr,"%s: after  smoothing (%d)\n",fnam,nstp);
    for(m=0; m<nstp; m++)
    {
      fprintf(stderr,"%13.6e %13.6e\n",x[m],y[m]);
    }
    fprintf(stderr,"-----\n");
  }

  return 0;
}

double mod2eff(int n,double lmod)
{
  int m1,m2;
  double leff;

  m1 = (int)((lmod-mie_lmod_min[n])/mie_lmod_stp[n]);
  m2 = m1+1;
  if(m1 < 0)
  {
    m1 = 0;
    m2 = 1;
  }
  if(m2 >= mie_n_step)
  {
    m2 = mie_n_step-1;
    m1 = m2-1;
  }
  leff = INTERP(lmod,mie_lmod_leff_x[n][m1],mie_lmod_leff_x[n][m2],
                     mie_lmod_leff_y[n][m1],mie_lmod_leff_y[n][m2]);

  return leff;
}

double eff2mod(int n,double leff)
{
  int m1,m2;
  double lmod;

  m1 = (int)((leff-mie_leff_min[n])/mie_leff_stp[n]);
  m2 = m1+1;
  if(m1 < 0)
  {
    m1 = 0;
    m2 = 1;
  }
  if(m2 >= mie_n_step)
  {
    m2 = mie_n_step-1;
    m1 = m2-1;
  }
  lmod = INTERP(leff,mie_leff_lmod_x[n][m1],mie_leff_lmod_x[n][m2],
                     mie_leff_lmod_y[n][m1],mie_leff_lmod_y[n][m2]);

  return lmod;
}

int MiePrintout(void)
{
  int i,j,k,n;
  char fnam[] = "MiePrintout";
  FILE *fp;

  if((fp=fopen(mie_inp1,"w")) == NULL)
  {
    fprintf(stderr,"%s: cannot open %s\n",fnam,mie_inp1);
  }
  else
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      fprintf(fp,"%9.4f\n",mie_wlen[i]);
    }
    fclose(fp);
  }

  if((fp=fopen(mie_inp2,"w")) == NULL)
  {
    fprintf(stderr,"%s: cannot open %s\n",fnam,mie_inp2);
  }
  else
  {
    for(j=0; j<mie_n_angl; j++)
    {
      fprintf(fp,"%9.4f\n",mie_angl[j]);
    }
    fclose(fp);
  }

  for(n=0; n<mie_n_comp; n++)
  {
    snprintf(mie_inp3,MAXLINE,"mie_refl_%02d.dat",n);
    if((fp=fopen(mie_inp3,"w")) == NULL)
    {
      fprintf(stderr,"%s: cannot open %s\n",fnam,mie_inp3);
    }
    else
    {
      for(i=0; i<mie_n_wlen; i++)
      {
        fprintf(fp,"%9.4f %13.5e %13.5e\n",mie_wlen[i],mie_refr_com[n][i],mie_refi_com[n][i]);
      }
      fclose(fp);
    }
  }

  if((fp=fopen(mie_inp4,"w")) == NULL)
  {
    fprintf(stderr,"%s: cannot open %s\n",fnam,mie_inp4);
  }
  else
  {
    for(i=0; i<SIM_NPAR; i++)
    {
      fprintf(fp,"%2d %22.14e\n",i,sim_pval[i]);
    }
    fclose(fp);
  }

  if((fp=fopen(mie_inp5,"w")) == NULL)
  {
    fprintf(stderr,"%s: cannot open %s\n",fnam,mie_inp5);
  }
  else
  {
    for(n=0; n<mie_n_comp; n++)
    {
      fprintf(fp,"%22.14e %22.14e %22.14e mie_refl_%02d.dat mie_refl_%02d.dat\n",
                  mie_wcom_com[n],pow(10.0,mie_lmod_com[n]),mie_lsgm_com[n],n,n);
    }
    fclose(fp);
  }

  if((fp=fopen(mie_out1,"w")) == NULL)
  {
    fprintf(stderr,"%s: cannot open %s\n",fnam,mie_out1);
  }
  else
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      fprintf(fp,"%7.2f %13.6e %13.6e %13.6e %13.6e\n",mie_wlen[i],mie_aext[i],mie_aext[i]/mie_aext[mie_iref],
                                                       mie_asca[i]/mie_aext[i],mie_asym[i]);
    }
    fclose(fp);
  }

  if((fp=fopen(mie_out2,"w")) == NULL)
  {
    fprintf(stderr,"%s: cannot open %s\n",fnam,mie_out2);
  }
  else
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        fprintf(fp,"%7.2f %7.2f %13.6e\n",mie_wlen[i],mie_angl[j],mie_phas[k]);
      }
    }
    fclose(fp);
  }

  return 0;
}

int ResFactor(int inst,int nw,int sbmp,double wlen,double amas,double *val)
{
  int m,n;
  double r;
  char fnam[] = "ResFactor";

  m = inst-1;
  if(m<0 || m>=INS_N_INST)
  {
    fprintf(stderr,"%s: invalid instrument ID >>> %d\n",fnam,inst);
    return -1;
  }
  if(nw<0 || nw>=INS_N_WLEN)
  {
    fprintf(stderr,"%s: invalid wavelength number >>> %d\n",fnam,nw);
    return -1;
  }
  n = (sbmp==1?0:(sbmp==5?1:(sbmp==15?2:-1)));
  if(n < 0)
  {
    fprintf(stderr,"%s: invalid band model >>> %d\n",fnam,sbmp);
    return -1;
  }
  if(fabs(wlen-ins_wlen[m][nw]) > DELTA)
  {
    fprintf(stderr,"%s: error, wlen=%13.6e, ins_wlen[%d][%d]=%13.6e\n",fnam,
                    wlen,m,nw,ins_wlen[m][nw]);
    return -1;
  }
  if(amas < EPSILON)
  {
    r = rcr_fact[m][nw][n];
  }
  else
  {
    r = rcr_fact[m][nw][n]+rcr_slop[m][nw][n]*(amas-rcr_amas);
  }
  *val = 1.0/r;

  return 0;
}

int CalcFactor(double *par,int mode,int flag)
{
  // Notice: The following flags (global variables) are modified in this function
  // and they MUST be set appropriately AFTER calling this function.
  // sim_flag, sim_cflg, sim_dflg
  // mode: 0=MODTRAN model, 1=Custom model
  // flag: 0=Use given parameters 1=Use crc parameters if possible
  int i,n;
  int rng;
  int rid;
  double fmin;

  // calculate correction factor
  if(flag==1 && strcmp(crc_ffac,NONAME)!=0)
  {
    if(ReadFactor() < 0) return -1;
  }
  else
  {
    sim_cflg = 0; // Do NOT apply correction
    sim_dflg = 1; // DIS=t
    if(flag==1 && crc_n_dsim>0)
    {
      for(n=0; n<dsr_n_data; n++)
      {
        for(rng=0; rng<DSR_MAXRANG; rng++)
        {
          for(i=0; i<dsr_n_wlen[rng]; i++)
          {
            dsr_dsim[rng][n][i] = crc_dsr_dsim[rng][n][i];
          }
          if(cnt_sflg == 1)
          {
            for(i=0; i<dsr_n_wlen[rng]; i++)
            {
              dsr_dsgl[rng][n][i] = crc_dsr_dsgl[rng][n][i];
            }
          }
        }
      }
      for(n=0; n<aur_n_data; n++)
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          for(i=0; i<aur_n_wlen[rng]; i++)
          {
            aur_dsim[rng][n][i] = crc_aur_dsim[rng][n][i];
          }
          if(cnt_sflg == 1)
          {
            for(i=0; i<aur_n_wlen[rng]; i++)
            {
              aur_dsgl[rng][n][i] = crc_aur_dsgl[rng][n][i];
            }
          }
        }
      }
      for(n=0; n<ssr_n_data; n++)
      {
        for(rng=0; rng<SSR_MAXRANG; rng++)
        {
          for(i=0; i<ssr_n_wlen[rng]; i++)
          {
            ssr_dsim[rng][n][i] = crc_ssr_dsim[rng][n][i];
          }
          if(cnt_sflg == 1)
          {
            for(i=0; i<ssr_n_wlen[rng]; i++)
            {
              ssr_dsgl[rng][n][i] = crc_ssr_dsgl[rng][n][i];
            }
          }
        }
      }
    }
    else
    {
      if(mode == 1)
      {
        sim_flag = 1; // Mie calculation is omissible, but MODTRAN calculation is NOT omissible
        if(cnt_dsr[0] == 1) DSR_SSR_CUSTOM_S0(par);
        if(cnt_dsr[1] == 1) DSR_SSR_CUSTOM_S1(par);
        if(cnt_dsr[2] == 1) DSR_SSR_CUSTOM_S2(par);
        if(cnt_dsr[3] == 1) DSR_SSR_CUSTOM_S3(par);
        if(cnt_aur[0] == 1) AUR_FCN_CUSTOM_S0(par,&fmin);
        if(cnt_aur[1] == 1) AUR_FCN_CUSTOM_S1(par,&fmin);
        if(cnt_ssr[0] == 1) SSR_FCN_CUSTOM_S0(par,&fmin);
        if(cnt_ssr[1] == 1) SSR_FCN_CUSTOM_S1(par,&fmin);
        if(cnt_ssr[2] == 1) SSR_FCN_CUSTOM_S2(par,&fmin);
        if(cnt_ssr[3] == 1) SSR_FCN_CUSTOM_S3(par,&fmin);
      } else
      if(mode == 0)
      {
        if(cnt_dsr[0] == 1) DSR_SSR_MODTRAN_S0(par);
        if(cnt_dsr[1] == 1) DSR_SSR_MODTRAN_S1(par);
        if(cnt_dsr[2] == 1) DSR_SSR_MODTRAN_S2(par);
        if(cnt_dsr[3] == 1) DSR_SSR_MODTRAN_S3(par);
        if(cnt_aur[0] == 1) AUR_FCN_MODTRAN_S0(par,&fmin);
        if(cnt_aur[1] == 1) AUR_FCN_MODTRAN_S1(par,&fmin);
        if(cnt_ssr[0] == 1) SSR_FCN_MODTRAN_S0(par,&fmin);
        if(cnt_ssr[1] == 1) SSR_FCN_MODTRAN_S1(par,&fmin);
        if(cnt_ssr[2] == 1) SSR_FCN_MODTRAN_S2(par,&fmin);
        if(cnt_ssr[3] == 1) SSR_FCN_MODTRAN_S3(par,&fmin);
      }
      else
      {
        fprintf(stderr,"CalcFactor: invalid mode >>> %d\n",mode);
        return -1;
      }
    }
    if(cnt_sflg == 1)
    {
      for(n=0; n<dsr_n_data; n++)
      {
        for(rng=0; rng<DSR_MAXRANG; rng++)
        {
          for(i=0; i<dsr_n_wlen[rng]; i++)
          {
            dsr_dfac[rng][n][i] = dsr_dsim[rng][n][i]-dsr_dsgl[rng][n][i];
          }
        }
      }
      for(n=0; n<aur_n_data; n++)
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          for(i=0; i<aur_n_wlen[rng]; i++)
          {
            aur_dfac[rng][n][i] = aur_dsim[rng][n][i]-aur_dsgl[rng][n][i];
          }
        }
      }
      for(n=0; n<ssr_n_data; n++)
      {
        for(rng=0; rng<SSR_MAXRANG; rng++)
        {
          for(i=0; i<ssr_n_wlen[rng]; i++)
          {
            ssr_dfac[rng][n][i] = ssr_dsim[rng][n][i]-ssr_dsgl[rng][n][i];
          }
        }
      }
    }
    else
    {
      for(n=0; n<dsr_n_data; n++)
      {
        for(rng=0; rng<DSR_MAXRANG; rng++)
        {
          for(i=0; i<dsr_n_wlen[rng]; i++)
          {
            dsr_dfac[rng][n][i] = dsr_dsim[rng][n][i];
          }
        }
      }
      for(n=0; n<aur_n_data; n++)
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          for(i=0; i<aur_n_wlen[rng]; i++)
          {
            aur_dfac[rng][n][i] = aur_dsim[rng][n][i];
          }
        }
      }
      for(n=0; n<ssr_n_data; n++)
      {
        for(rng=0; rng<SSR_MAXRANG; rng++)
        {
          for(i=0; i<ssr_n_wlen[rng]; i++)
          {
            ssr_dfac[rng][n][i] = ssr_dsim[rng][n][i];
          }
        }
      }
    }
    sim_dflg = 0; // DIS=f
    if(mode == 1)
    {
      if(cnt_dsr[0] == 1) DSR_SSR_CUSTOM_B0(par);
      if(cnt_dsr[1] == 1) DSR_SSR_CUSTOM_B1(par);
      if(cnt_dsr[2] == 1) DSR_SSR_CUSTOM_B2(par);
      if(cnt_dsr[3] == 1) DSR_SSR_CUSTOM_B3(par);
      if(cnt_aur[0] == 1) AUR_FCN_CUSTOM_B0(par,&fmin);
      if(cnt_aur[1] == 1) AUR_FCN_CUSTOM_B1(par,&fmin);
      if(cnt_ssr[0] == 1) SSR_FCN_CUSTOM_B0(par,&fmin);
      if(cnt_ssr[1] == 1) SSR_FCN_CUSTOM_B1(par,&fmin);
      if(cnt_ssr[2] == 1) SSR_FCN_CUSTOM_B2(par,&fmin);
      if(cnt_ssr[3] == 1) SSR_FCN_CUSTOM_B3(par,&fmin);
    } else
    if(mode == 0)
    {
      if(cnt_dsr[0] == 1) DSR_SSR_MODTRAN_B0(par);
      if(cnt_dsr[1] == 1) DSR_SSR_MODTRAN_B1(par);
      if(cnt_dsr[2] == 1) DSR_SSR_MODTRAN_B2(par);
      if(cnt_dsr[3] == 1) DSR_SSR_MODTRAN_B3(par);
      if(cnt_aur[0] == 1) AUR_FCN_MODTRAN_B0(par,&fmin);
      if(cnt_aur[1] == 1) AUR_FCN_MODTRAN_B1(par,&fmin);
      if(cnt_ssr[0] == 1) SSR_FCN_MODTRAN_B0(par,&fmin);
      if(cnt_ssr[1] == 1) SSR_FCN_MODTRAN_B1(par,&fmin);
      if(cnt_ssr[2] == 1) SSR_FCN_MODTRAN_B2(par,&fmin);
      if(cnt_ssr[3] == 1) SSR_FCN_MODTRAN_B3(par,&fmin);
    }
    else
    {
      fprintf(stderr,"CalcFactor: invalid mode >>> %d\n",mode);
      return -1;
    }
    if(cnt_sflg == 1)
    {
      for(n=0; n<dsr_n_data; n++)
      {
        for(rng=0; rng<DSR_MAXRANG; rng++)
        {
          for(i=0; i<dsr_n_wlen[rng]; i++)
          {
            dsr_dfac[rng][n][i] /= (dsr_dsim[rng][n][i]-dsr_dsgl[rng][n][i]);
          }
        }
      }
      for(n=0; n<aur_n_data; n++)
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          for(i=0; i<aur_n_wlen[rng]; i++)
          {
            aur_dfac[rng][n][i] /= (aur_dsim[rng][n][i]-aur_dsgl[rng][n][i]);
          }
        }
      }
      for(n=0; n<ssr_n_data; n++)
      {
        for(rng=0; rng<SSR_MAXRANG; rng++)
        {
          for(i=0; i<ssr_n_wlen[rng]; i++)
          {
            ssr_dfac[rng][n][i] /= (ssr_dsim[rng][n][i]-ssr_dsgl[rng][n][i]);
          }
        }
      }
    }
    else
    {
      for(n=0; n<dsr_n_data; n++)
      {
        for(rng=0; rng<DSR_MAXRANG; rng++)
        {
          for(i=0; i<dsr_n_wlen[rng]; i++)
          {
            dsr_dfac[rng][n][i] /= dsr_dsim[rng][n][i];
          }
        }
      }
      for(n=0; n<aur_n_data; n++)
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          for(i=0; i<aur_n_wlen[rng]; i++)
          {
            aur_dfac[rng][n][i] /= aur_dsim[rng][n][i];
          }
        }
      }
      for(n=0; n<ssr_n_data; n++)
      {
        for(rng=0; rng<SSR_MAXRANG; rng++)
        {
          for(i=0; i<ssr_n_wlen[rng]; i++)
          {
            ssr_dfac[rng][n][i] /= ssr_dsim[rng][n][i];
          }
        }
      }
    }
  }
  if(cnt_vb)
  {
    fprintf(stderr,"DIS=f->t correction factor:\n");
    for(n=0; n<dsr_n_data; n++)
    {
      for(rng=0; rng<DSR_MAXRANG; rng++)
      {
        rid = rng;
        for(i=0; i<dsr_n_wlen[rng]; i++)
        {
          fprintf(stderr,"%2d %4d %9.4f %22.14e\n",rid,n,dsr_wlen[rng][i],dsr_dfac[rng][n][i]);
        }
      }
    }
    for(n=0; n<aur_n_data; n++)
    {
      for(rng=0; rng<AUR_MAXRANG; rng++)
      {
        rid = rng+DSR_MAXRANG*2;
        for(i=0; i<aur_n_wlen[rng]; i++)
        {
          fprintf(stderr,"%2d %4d %9.4f %22.14e\n",rid,n,aur_wlen[rng][i],aur_dfac[rng][n][i]);
        }
      }
    }
    for(n=0; n<ssr_n_data; n++)
    {
      for(rng=0; rng<SSR_MAXRANG; rng++)
      {
        rid = rng+DSR_MAXRANG*2+AUR_MAXRANG*2;
        for(i=0; i<ssr_n_wlen[rng]; i++)
        {
          fprintf(stderr,"%2d %4d %9.4f %22.14e\n",rid,n,ssr_wlen[rng][i],ssr_dfac[rng][n][i]);
        }
      }
    }
  }

  return 0;
}

int Background(const double *par,int mode)
{
  // mode: 0=MODTRAN model, single
  //       1=MODTRAN model, band
  //       2=Custom  model, single
  //       3=Custom  model, band
  //       4=Use the same factor
  int i,j;
  int rng;
  int ith,iph;
  int err;
  int fpar[OPT_N_FPAR];
  double r,s,x,y,z;
  double th,ph;
  double xi,yi,zi;
  double thi,phi;
  double th0,th1,th2;
  double dth,dph,ds;
  double p[6];

  if(cnt_vb)
  {
    fprintf(stderr,"Background:\n");
  }
  if(mode < 4)
  {
    if(mode==2 || mode==3)
    {
      err = 0;
      do
      {
        if(CheckParam(par,fpar) < 0)
        {
          err = 1;
          break;
        }
        if(fpar[OPT_FPAR_DIST] > 0)
        {
          if(MieCalc(par,fpar) < 0)
          {
            err = 1;
            break;
          }
        }
        if(fpar[OPT_FPAR_MIXR] > 0)
        {
          if(MixComp(par) < 0)
          {
            return -1;
          }
        }
      }
      while(0);
      if(err)
      {
        return -1;
      }
    }

    // Calculate solar irradiance
    sim_vmie = par[OPT_PAR_VMIE];
    sim_wscl = par[OPT_PAR_WSCL];
    for(rng=0; rng<DSR_MAXRANG; rng++)
    {
      if(cnt_dsr[rng] == 0) continue;
      p[0] = dsr_wlen[rng][0]-sim_xmgn;
      p[1] = dsr_wlen[rng][dsr_n_wlen[rng]-1]+sim_xmgn;
      err = 0;
      do
      {
        switch(mode)
        {
          case 0:
          case 1:
            if(WriteTape5(rng,p,SIM_MODE_DSR_MODTRAN_BA) < 0) err = 1;
            break;
          case 2:
          case 3:
            if(WriteTape5(rng,p,SIM_MODE_DSR_CUSTOM_BA) < 0) err = 1;
            break;
          default:
            fprintf(stderr,"Background: invalid mode >>> %d\n",mode);
            err = 1;
            break;
        }
        if(err) break;
        if(RunModtran() < 0)
        {
          err = 1;
          break;
        }
        if(ReadDSR_B(rng,p) < 0)
        {
          err = 1;
          break;
        }
      }
      while(0);
      if(err)
      {
        return -1;
      }
    }

    // Calculate sky radiance
    sim_sbmp = SSR_SBMP;
    th0 = dsr_tmin*D_TO_R;
    dth = (dsr_tmax-dsr_tmin)*D_TO_R/dsr_n_tdiv;
    dph = PI2/dsr_n_pdiv;
    for(i=0; i<dsr_n_data; i++)
    {
      p[2] = dsr_th_sun[i];
      p[3] = dsr_ph_sun[i];
      for(rng=0; rng<DSR_MAXRANG; rng++)
      {
        for(j=0; j<dsr_n_wlen[rng]; j++)
        {
          if(mode==1 || mode==3)
          {
            p[0] = dsr_wlen[rng][j]-sim_xmgn;
            p[1] = dsr_wlen[rng][j]+sim_xmgn;
          }
          else
          {
            p[0] = dsr_wlen[rng][j];
            p[1] = (double)sim_nsam;
          }
          dsr_dctr[rng][i][j] = 0.0;
          dsr_dsrc[rng][i][j] = 0.0;
          dsr_unif[rng][i][j] = 0.0;
          for(ith=-1; ith<dsr_n_tdiv; ith++)
          {
            if(ith < 0)
            {
              thi = 0.0;
            }
            else
            {
              th1 = th0+dth*ith;
              th2 = th0+dth*(ith+1);
              thi = 0.5*(th1+th2);
              ds = dsr_romg[ith]*(cos(th1)-cos(th2))*dph;
            }
            for(iph=0; iph<dsr_n_pdiv; iph++)
            {
              phi = (ith<0?0.0:dph*((double)iph+0.5));
              PolToXyz(thi,phi,&xi,&yi,&zi);
              IjkToXyz(dsr_th_sun[i]*D_TO_R,dsr_ph_sun[i]*D_TO_R,xi,yi,zi,&x,&y,&z);
              XyzToPol(x,y,z,&th,&ph);
              p[4] = th*R_TO_D;
              p[5] = ph*R_TO_D;
              err = 0;
              do
              {
                switch(mode)
                {
                  case 0:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_SO) < 0) err = 1;
                    break;
                  case 1:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_BO) < 0) err = 1;
                    break;
                  case 2:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_SO) < 0) err = 1;
                    break;
                  case 3:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_BO) < 0) err = 1;
                    break;
                  default:
                    fprintf(stderr,"Background: invalid mode >>> %d\n",mode);
                    err = 1;
                    break;
                }
                if(err) break;
                if(RunModtran() < 0)
                {
                  err = 1;
                  break;
                }
                if(mode==1 || mode==3)
                {
                  if(ReadPlt_B(p,1,&dsr_wlen[rng][j],&y) < 0)
                  {
                    err = 1;
                    break;
                  }
                }
                else
                {
                  if(ReadPlt_S(p,&y) < 0)
                  {
                    err = 1;
                    break;
                  }
                  if(ResFactor(dsr_ms720,dsr_nw[rng][j],sim_sbmp,dsr_wlen[rng][j],dsr_amas[i],&r) < 0)
                  {
                    err = 1;
                    break;
                  }
                  y *= r;
                }
              }
              while(0);
              if(err)
              {
                return -1;
              }
              if(mode==1 || mode==3)
              {
                if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
                {
                  if(cnt_sflg == 1) // calculate single scattering
                  {
                    sim_mult = 0;
                    err = 0;
                    do
                    {
                      switch(mode)
                      {
                        case 1:
                          if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_BO) < 0) err = 1;
                          break;
                        case 3:
                          if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_BO) < 0) err = 1;
                          break;
                      }
                      if(err) break;
                      if(RunModtran() < 0)
                      {
                        err = 1;
                        break;
                      }
                      if(ReadPlt_B(p,1,&dsr_wlen[rng][j],&s) < 0)
                      {
                        err = 1;
                        break;
                      }
                    }
                    while(0);
                    sim_mult = 1;
                    if(err)
                    {
                      return -1;
                    }
                    y = s+(y-s)*dsr_dfac[rng][i][j];
                  }
                  else
                  {
                    y *= dsr_dfac[rng][i][j];
                  }
                }
              }
              if(cnt_vb > 1)
              {
                fprintf(stderr,"%2d %4d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %22.14e\n",
                                rng,i,p[2],p[3],dsr_wlen[rng][j],
                                thi*R_TO_D,phi*R_TO_D,p[4],p[5],y);
              }
              if(ith < 0)
              {
                dsr_dctr[rng][i][j] = y;
                break;
              }
              else
              {
                dsr_unif[rng][i][j] += y*ds;
              }
            }
          }
          dsr_dsrc[rng][i][j] = dsr_unif[rng][i][j]/dsr_dsim[rng][i][j];
          dsr_unif[rng][i][j] = dsr_dctr[rng][i][j]*dsr_sfov/dsr_unif[rng][i][j];
        }
      }
    }
  } else
  if(cnt_vb)
  {
    fprintf(stderr,"Using the same factor.\n");
  }

  // correct for sky background
  for(i=0; i<dsr_n_data; i++)
  {
    for(rng=0; rng<DSR_MAXRANG; rng++)
    {
      for(j=0; j<dsr_n_wlen[rng]; j++)
      {
        dsr_data[rng][i][j] = dsr_dorg[rng][i][j]*(1.0-dsr_dsrc[rng][i][j]);
      }
    }
  }

  if(cnt_vb)
  {
    for(i=0; i<dsr_n_data; i++)
    {
      for(rng=0; rng<DSR_MAXRANG; rng++)
      {
        for(j=0; j<dsr_n_wlen[rng]; j++)
        {
          fprintf(stderr,"%2d %4d %9.4f %9.4f %9.4f %22.14e %22.14e %22.14e\n",
                          rng,i,dsr_th_sun[i],dsr_ph_sun[i],dsr_wlen[rng][j],
                          dsr_dctr[rng][i][j],dsr_unif[rng][i][j],dsr_dsrc[rng][i][j]);
        }
      }
    }
  }

  return 0;
}

int RingCenter(const double *par,int mode)
{
  // mode: 0=MODTRAN model, single
  //       1=MODTRAN model, band
  //       2=Custom  model, single
  //       3=Custom  model, band
  //       4=Use the same factor
  int i,j;
  int rng;
  int ith,iph;
  int err;
  int fpar[OPT_N_FPAR];
  double r,s,x,y,z;
  double th,ph;
  double xi,yi,zi;
  double thi,phi;
  double th0,th1,th2;
  double dth,dph,ds;
  double p[6];

  if(cnt_vb)
  {
    fprintf(stderr,"RingCenter:\n");
  }
  if(mode < 4)
  {
    if(mode==2 || mode==3)
    {
      err = 0;
      do
      {
        if(CheckParam(par,fpar) < 0)
        {
          err = 1;
          break;
        }
        if(fpar[OPT_FPAR_DIST] > 0)
        {
          if(MieCalc(par,fpar) < 0)
          {
            err = 1;
            break;
          }
        }
        if(fpar[OPT_FPAR_MIXR] > 0)
        {
          if(MixComp(par) < 0)
          {
            err = 1;
            break;
          }
        }
      }
      while(0);
      if(err)
      {
        return -1;
      }
    }

    // Calculate sky radiance
    sim_vmie = par[OPT_PAR_VMIE];
    sim_wscl = par[OPT_PAR_WSCL];
    th0 = aur_tmin*D_TO_R;
    dth = (aur_tmax-aur_tmin)*D_TO_R/aur_n_tdiv;
    dph = PI2/aur_n_pdiv;
    for(i=0; i<aur_n_data; i++)
    {
      p[2] = aur_th_sun[i];
      p[3] = aur_ph_sun[i];
      for(rng=0; rng<AUR_MAXRANG; rng++)
      {
        for(j=0; j<aur_n_wlen[rng]; j++)
        {
          if(mode==1 || mode==3)
          {
            sim_sbmp = AUR_SBMP;
            p[0] = aur_wlen[rng][j]-sim_xmgn;
            p[1] = aur_wlen[rng][j]+sim_xmgn;
          }
          else
          {
            sim_sbmp = aur_sbmp[rng][j];
            p[0] = aur_wlen[rng][j];
            p[1] = (double)sim_nsam;
          }
          aur_dctr[rng][i][j] = 0.0;
          aur_unif[rng][i][j] = 0.0;
          for(ith=-1; ith<aur_n_tdiv; ith++)
          {
            if(ith >= 0)
            {
              th1 = th0+dth*ith;
              th2 = th0+dth*(ith+1);
              thi = 0.5*(th1+th2);
              ds = aur_romg[ith]*(cos(th1)-cos(th2))*dph;
            }
            for(iph=0; iph<aur_n_pdiv; iph++)
            {
              if(ith < 0)
              {
                p[4] = aur_th_los[i];
                p[5] = aur_ph_los[i];
                PolToXyz(aur_th_los[i]*D_TO_R,aur_ph_los[i]*D_TO_R,&x,&y,&z);
                XyzToIjk(aur_th_sun[i]*D_TO_R,aur_ph_sun[i]*D_TO_R,x,y,z,&xi,&yi,&zi);
                XyzToPol(xi,yi,zi,&thi,&phi);
              }
              else
              {
                phi = dph*((double)iph+0.5);
                PolToXyz(thi,phi,&xi,&yi,&zi);
                IjkToXyz(aur_th_sun[i]*D_TO_R,aur_ph_sun[i]*D_TO_R,xi,yi,zi,&x,&y,&z);
                XyzToPol(x,y,z,&th,&ph);
                p[4] = th*R_TO_D;
                p[5] = ph*R_TO_D;
              }
              err = 0;
              do
              {
                switch(mode)
                {
                  case 0:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_SO) < 0) err = 1;
                    break;
                  case 1:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_BO) < 0) err = 1;
                    break;
                  case 2:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_SO) < 0) err = 1;
                    break;
                  case 3:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_BO) < 0) err = 1;
                    break;
                  default:
                    fprintf(stderr,"RingCenter: invalid mode >>> %d\n",mode);
                    err = 1;
                    break;
                }
                if(err) break;
                if(RunModtran() < 0)
                {
                  err = 1;
                  break;
                }
                if(mode==1 || mode==3)
                {
                  if(ReadPlt_B(p,1,&aur_wlen[rng][j],&y) < 0)
                  {
                    err = 1;
                    break;
                  }
                }
                else
                {
                  if(ReadPlt_S(p,&y) < 0)
                  {
                    err = 1;
                    break;
                  }
                  if(ResFactor(aur_ms720,aur_nw[rng][j],sim_sbmp,aur_wlen[rng][j],aur_amas[i],&r) < 0)
                  {
                    err = 1;
                    break;
                  }
                  y *= r;
                }
              }
              while(0);
              if(err)
              {
                return -1;
              }
              if(mode==1 || mode==3)
              {
                if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
                {
                  if(cnt_sflg == 1) // calculate single scattering
                  {
                    sim_mult = 0;
                    err = 0;
                    do
                    {
                      switch(mode)
                      {
                        case 1:
                          if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_BO) < 0) err = 1;
                          break;
                        case 3:
                          if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_BO) < 0) err = 1;
                          break;
                      }
                      if(err) break;
                      if(RunModtran() < 0)
                      {
                        err = 1;
                        break;
                      }
                      if(ReadPlt_B(p,1,&aur_wlen[rng][j],&s) < 0)
                      {
                        err = 1;
                        break;
                      }
                    }
                    while(0);
                    sim_mult = 1;
                    if(err)
                    {
                      return -1;
                    }
                    y = s+(y-s)*aur_dfac[rng][i][j];
                  }
                  else
                  {
                    y *= aur_dfac[rng][i][j];
                  }
                }
              }
              if(cnt_vb > 1)
              {
                fprintf(stderr,"%2d %4d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %22.14e\n",
                                rng,i,p[2],p[3],aur_wlen[rng][j],
                                thi*R_TO_D,phi*R_TO_D,p[4],p[5],y);
              }
              if(ith < 0)
              {
                aur_dctr[rng][i][j] = y;
                break;
              }
              else
              {
                aur_unif[rng][i][j] += y*ds;
              }
            }
          }
          aur_unif[rng][i][j] = aur_dctr[rng][i][j]*aur_sfov/aur_unif[rng][i][j];
        }
      }
    }
  } else
  if(cnt_vb)
  {
    fprintf(stderr,"Using the same factor.\n");
  }

  // correct for ring center
  for(i=0; i<aur_n_data; i++)
  {
    for(rng=0; rng<AUR_MAXRANG; rng++)
    {
      for(j=0; j<aur_n_wlen[rng]; j++)
      {
        aur_data[rng][i][j] = aur_dorg[rng][i][j]*aur_unif[rng][i][j];
      }
    }
  }

  if(cnt_vb)
  {
    for(i=0; i<aur_n_data; i++)
    {
      for(rng=0; rng<AUR_MAXRANG; rng++)
      {
        for(j=0; j<aur_n_wlen[rng]; j++)
        {
          fprintf(stderr,"%2d %4d %9.4f %9.4f %9.4f %22.14e %22.14e\n",
                          rng,i,aur_th_sun[i],aur_ph_sun[i],aur_wlen[rng][j],
                          aur_dctr[rng][i][j],aur_unif[rng][i][j]);
        }
      }
    }
  }

  return 0;
}

int Uniformity(const double *par,int mode)
{
  // mode: 0=MODTRAN model, single
  //       1=MODTRAN model, band
  //       2=Custom  model, single
  //       3=Custom  model, band
  //       4=Use the same factor
  int i,j;
  int rng;
  int ith,iph;
  int err;
  int fpar[OPT_N_FPAR];
  double r,s,x,y,z;
  double th,ph;
  double xi,yi,zi;
  double thi,phi;
  double th0,th1,th2;
  double dth,dph,ds;
  double p[6];

  if(cnt_vb)
  {
    fprintf(stderr,"Uniformity:\n");
  }
  if(mode < 4)
  {
    if(mode==2 || mode==3)
    {
      err = 0;
      do
      {
        if(CheckParam(par,fpar) < 0)
        {
          err = 1;
          break;
        }
        if(fpar[OPT_FPAR_DIST] > 0)
        {
          if(MieCalc(par,fpar) < 0)
          {
            err = 1;
            break;
          }
        }
        if(fpar[OPT_FPAR_MIXR] > 0)
        {
          if(MixComp(par) < 0)
          {
            err = 1;
            break;
          }
        }
      }
      while(0);
      if(err)
      {
        return -1;
      }
    }

    // Calculate sky radiance
    sim_vmie = par[OPT_PAR_VMIE];
    sim_wscl = par[OPT_PAR_WSCL];
    th0 = ssr_tmin*D_TO_R;
    dth = (ssr_tmax-ssr_tmin)*D_TO_R/ssr_n_tdiv;
    dph = PI2/ssr_n_pdiv;
    for(i=0; i<ssr_n_data; i++)
    {
      p[2] = ssr_th_sun[i];
      p[3] = ssr_ph_sun[i];
      for(rng=0; rng<SSR_MAXRANG; rng++)
      {
        for(j=0; j<ssr_n_wlen[rng]; j++)
        {
          if(mode==1 || mode==3)
          {
            sim_sbmp = SSR_SBMP;
            p[0] = ssr_wlen[rng][j]-sim_xmgn;
            p[1] = ssr_wlen[rng][j]+sim_xmgn;
          }
          else
          {
            sim_sbmp = ssr_sbmp[rng][j];
            p[0] = ssr_wlen[rng][j];
            p[1] = (double)sim_nsam;
          }
          ssr_dctr[rng][i][j] = 0.0;
          ssr_unif[rng][i][j] = 0.0;
          for(ith=-1; ith<ssr_n_tdiv; ith++)
          {
            if(ith < 0)
            {
              thi = 0.0;
            }
            else
            {
              th1 = th0+dth*ith;
              th2 = th0+dth*(ith+1);
              thi = 0.5*(th1+th2);
              ds = ssr_romg[ith]*(cos(th1)-cos(th2))*dph;
            }
            for(iph=0; iph<ssr_n_pdiv; iph++)
            {
              phi = (ith<0?0.0:dph*((double)iph+0.5));
              PolToXyz(thi,phi,&xi,&yi,&zi);
              IjkToXyz(ssr_th_los[i]*D_TO_R,ssr_ph_los[i]*D_TO_R,xi,yi,zi,&x,&y,&z);
              XyzToPol(x,y,z,&th,&ph);
              p[4] = th*R_TO_D;
              p[5] = ph*R_TO_D;
              err = 0;
              do
              {
                switch(mode)
                {
                  case 0:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_SO) < 0) err = 1;
                    break;
                  case 1:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_BO) < 0) err = 1;
                    break;
                  case 2:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_SO) < 0) err = 1;
                    break;
                  case 3:
                    if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_BO) < 0) err = 1;
                    break;
                  default:
                    fprintf(stderr,"Uniformity: invalid mode >>> %d\n",mode);
                    err = 1;
                    break;
                }
                if(err) break;
                if(RunModtran() < 0)
                {
                  err = 1;
                  break;
                }
                if(mode==1 || mode==3)
                {
                  if(ReadPlt_B(p,1,&ssr_wlen[rng][j],&y) < 0)
                  {
                    err = 1;
                    break;
                  }
                }
                else
                {
                  if(ReadPlt_S(p,&y) < 0)
                  {
                    err = 1;
                    break;
                  }
                  if(ResFactor(ssr_ms720,ssr_nw[rng][j],sim_sbmp,ssr_wlen[rng][j],ssr_amas[i],&r) < 0)
                  {
                    err = 1;
                    break;
                  }
                  y *= r;
                }
              }
              while(0);
              if(err)
              {
                return -1;
              }
              if(mode==1 || mode==3)
              {
                if(sim_cflg==1 && sim_dflg==0) // apply correction && DIS==f
                {
                  if(cnt_sflg == 1) // calculate single scattering
                  {
                    sim_mult = 0;
                    err = 0;
                    do
                    {
                      switch(mode)
                      {
                        case 1:
                          if(WriteTape5(rng,p,SIM_MODE_SSR_MODTRAN_BO) < 0) err = 1;
                          break;
                        case 3:
                          if(WriteTape5(rng,p,SIM_MODE_SSR_CUSTOM_BO) < 0) err = 1;
                          break;
                      }
                      if(err) break;
                      if(RunModtran() < 0)
                      {
                        err = 1;
                        break;
                      }
                      if(ReadPlt_B(p,1,&ssr_wlen[rng][j],&s) < 0)
                      {
                        err = 1;
                        break;
                      }
                    }
                    while(0);
                    sim_mult = 1;
                    if(err)
                    {
                      return -1;
                    }
                    y = s+(y-s)*ssr_dfac[rng][i][j];
                  }
                  else
                  {
                    y *= ssr_dfac[rng][i][j];
                  }
                }
              }
              if(cnt_vb > 1)
              {
                fprintf(stderr,"%2d %4d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %22.14e\n",
                                rng,i,p[2],p[3],ssr_wlen[rng][j],
                                thi*R_TO_D,phi*R_TO_D,p[4],p[5],y);
              }
              if(ith < 0)
              {
                ssr_dctr[rng][i][j] = y;
                break;
              }
              else
              {
                ssr_unif[rng][i][j] += y*ds;
              }
            }
          }
          ssr_unif[rng][i][j] = ssr_dctr[rng][i][j]*ssr_sfov/ssr_unif[rng][i][j];
        }
      }
    }
  } else
  if(cnt_vb)
  {
    fprintf(stderr,"Using the same factor.\n");
  }

  // correct for uniformity
  for(i=0; i<ssr_n_data; i++)
  {
    for(rng=0; rng<SSR_MAXRANG; rng++)
    {
      for(j=0; j<ssr_n_wlen[rng]; j++)
      {
        ssr_data[rng][i][j] = ssr_dorg[rng][i][j]*ssr_unif[rng][i][j];
      }
    }
  }

  if(cnt_vb)
  {
    for(i=0; i<ssr_n_data; i++)
    {
      for(rng=0; rng<SSR_MAXRANG; rng++)
      {
        for(j=0; j<ssr_n_wlen[rng]; j++)
        {
          fprintf(stderr,"%2d %4d %9.4f %9.4f %9.4f %22.14e %22.14e\n",
                          rng,i,ssr_th_los[i],ssr_ph_los[i],ssr_wlen[rng][j],
                          ssr_dctr[rng][i][j],ssr_unif[rng][i][j]);
        }
      }
    }
  }

  return 0;
}

int PrfInitMolecule(int imod)
{
  int i,j;
  int nmod;

  if(imod<=0 || imod>PRF_DEF_NMOD)
  {
    fprintf(stderr,"PrfInitMolecule: invalid model number >>> %d\n",imod);
    return -1;
  }
  else
  {
    nmod = imod-1;
  }

  if(Sampling(PRF_DEF_NALT,prf_def_alt,prf_def_pres[nmod],prf_n_alt,prf_alt,prf_pres,1.0) < 0) return -1;
  if(Sampling(PRF_DEF_NALT,prf_def_alt,prf_def_temp[nmod],prf_n_alt,prf_alt,prf_temp,1.0) < 0) return -1;
  for(i=0; i<PRF_N_AMOL; i++)
  {
    if(Sampling(PRF_DEF_NALT,prf_def_alt,prf_def_amol[i][nmod],prf_n_alt,prf_alt,prf_wmol[i],1.0) < 0) return -1;
  }
  for(j=0; j<PRF_N_TRAC; i++,j++)
  {
    if(Sampling(PRF_DEF_NALT,prf_def_alt,prf_def_trac[j],prf_n_alt,prf_alt,prf_wmol[i],1.0) < 0) return -1;
  }

  for(j=0; j<prf_n_alt; j++)
  {
    sim_prf_pres[j] = 0.0;
    sim_prf_temp[j] = 0.0;
  }
  for(i=0; i<PRF_NMOL; i++)
  {
    sim_prf_pscl[i] = 0;
    sim_prf_vscl[i] = NAN;
    for(j=0; j<prf_n_alt; j++)
    {
      sim_prf_wmol[i][j] = 0.0;
    }
  }

  for(i=0; i<PRF_NPAR; i++)
  {
    prf_jchar[i] = 'A';
    sim_prf_jchar[i] = (char)('0'+imod);
  }
  sim_prf_jcharx = (char)('0'+imod);

  sim_2c = 0;

  return 0;
}

int PrfInitAerosol(int ihaz,int isea,int ivul,double vis)
{
  int i,j;
  float v;
  double p[PRF_HAZ_NALT];
  double q[PRF_MAXNALT];

  v = (float)vis;
  for(j=0; j<PRF_HAZ_NALT; j++)
  {
    i = j+1;
    p[j] = (double)aerprf_(&i,&v,&ihaz,&isea,&ivul);
  }
  if(SamplingE(PRF_HAZ_NALT,prf_haz_alt,p,prf_n_alt,prf_alt,q,1.0) < 0) return -1;
  for(i=0; i<PRF_NHAZ; i++)
  {
    for(j=0; j<prf_n_alt; j++)
    {
      prf_haze[i][j] = 0.0;
    }
  }
  if(ihaz != 0)
  {
    for(j=0; j<prf_n_alt; j++)
    {
      if(prf_alt[j] < 2.0+EPSILON)
      {
        prf_haze[0][j] = q[j];
      } else
      if(prf_alt[j] < 10.0+EPSILON)
      {
        prf_haze[1][j] = q[j];
      } else
      if(prf_alt[j] < 30.0+EPSILON)
      {
        prf_haze[2][j] = q[j];
      }
      else
      {
        prf_haze[3][j] = q[j];
      }
    }
  }

  for(i=0; i<PRF_NHAZ; i++)
  {
    for(j=0; j<prf_n_alt; j++)
    {
      sim_prf_haze[i][j] = prf_haze[i][j];
    }
  }

  return 0;
}

int PrfGround(void)
{
  int i,j;
  double f;

  if(!isnan(prf_pres_gnd))
  {
    f = prf_pres_gnd/prf_pres[0];
    for(j=0; j<prf_n_alt; j++)
    {
      sim_prf_pres[j] = prf_pres[j]*f;
    }
    sim_prf_jchar[0] = prf_jchar[0];
    for(i=0; i<PRF_NMOL; i++)
    {
      if(sim_prf_pscl[i] != 0)
      {
        for(j=0; j<prf_n_alt; j++)
        {
          sim_prf_wmol[i][j] = prf_wmol[i][j]*f;
        }
        sim_prf_jchar[i+2] = prf_jchar[i+2];
        if(i > 2) sim_2c |= 0x0010;
      }
    }
  }

  if(!isnan(prf_temp_gnd))
  {
    f = prf_temp_gnd/prf_temp[0];
    for(j=0; j<prf_n_alt; j++)
    {
      sim_prf_temp[j] = prf_temp[j]*f;
    }
    sim_prf_jchar[1] = prf_jchar[1];
  }

  sim_2c |= 0x0001;

  return 0;
}

int PrfScaleMolecule(void)
{
  int i,j;

  for(i=0; i<PRF_NMOL; i++)
  {
    if(!isnan(sim_prf_vscl[i]))
    {
      for(j=0; j<prf_n_alt; j++)
      {
        sim_prf_wmol[i][j] = prf_wmol[i][j]*sim_prf_vscl[i];
      }
      sim_prf_jchar[i+2] = prf_jchar[i+2];
      if(i > 2) sim_2c |= 0x0010;
    }
  }

  sim_2c |= 0x0001;

  return 0;
}

int PrfScaleAerosol(double t)
{
  int i,j;
  double f;

  f = (t-prf_aod_hi)/prf_aod_lo;
  for(i=0; i<PRF_NHAZ_LOW; i++)
  {
    for(j=0; j<prf_n_alt; j++)
    {
      sim_prf_haze[i][j] = prf_haze[i][j]*f;
    }
  }

  sim_2c |= 0x1000;

  return 0;
}

int PrfCalcAOD(void)
{
  int i,j;
  double a1,a2;
  double dz,dt;
  double r;
  double t[PRF_NHAZ];

  for(i=0; i<PRF_NHAZ; i++)
  {
    t[i] = 0.0;
    for(j=0; j<prf_n_alt-1; j++)
    {
      a1 = prf_haze[i][j];
      a2 = prf_haze[i][j+1];
      dz = prf_alt[j+1]-prf_alt[j];
      if(a1>0.0 && a2>0.0 && (fabs(a1-a2)>a2*EPSILON))
      {
        r = dz/log(a1/a2);
        if(dz/r < EPSILON)
        {
          r = 0.0;
        }
      }
      else
      {
        r = 0.0;
      }
      if(r > 0.0)
      {
        dt = a1*r*(1.0-exp(-dz/r));
      }
      else
      {
        dt = 0.5*(a1+a2)*dz;
      }
      t[i] += dt;
    }
  }

  prf_aod_lo = 0.0;
  for(i=0; i<PRF_NHAZ_LOW; i++)
  {
    prf_aod_lo += t[i];
  }
  prf_aod_hi = 0.0;
  for(i=PRF_NHAZ_LOW; i<PRF_NHAZ; i++)
  {
    prf_aod_hi += t[i];
  }
  prf_aod = prf_aod_lo+prf_aod_hi;

  return 0;
}

int upper_pfunc(double rh)
{
  int i,j,k;
  int i_1,i_2;
  double p1,p2,p3;
  double *tropo_phas_tmp1 = NULL;
  double *tropo_phas_tmp2 = NULL;

  mie_tropo_phas = (double*)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
  mie_strat_phas = (double*)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
  mie_meteo_phas = (double*)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
  sim_tropo_phas = (double*)malloc(mie_n_wlen*sim_n_angl*sizeof(double));
  sim_strat_phas = (double*)malloc(mie_n_wlen*sim_n_angl*sizeof(double));
  sim_meteo_phas = (double*)malloc(mie_n_wlen*sim_n_angl*sizeof(double));
  if(mie_tropo_phas==NULL || mie_strat_phas==NULL || mie_meteo_phas==NULL ||
     sim_tropo_phas==NULL || sim_strat_phas==NULL || sim_meteo_phas==NULL)
  {
    fprintf(stderr,"upper_pfunc: error in allocating memory\n");
    return -1;
  }

  Interp2D(upp_wlen,upp_angl,upp_strat_phas,UPP_N_WLEN,UPP_N_ANGL,mie_wlen,mie_angl,mie_strat_phas,mie_n_wlen,mie_n_angl,1);
  Interp2D(upp_wlen,upp_angl,upp_meteo_phas,UPP_N_WLEN,UPP_N_ANGL,mie_wlen,mie_angl,mie_meteo_phas,mie_n_wlen,mie_n_angl,0);
  i_1 = i_2 = -1;
  for(i=UPP_N_RHUM-1; i>=0; i--)
  {
    if(upp_rhum[i] <= rh)
    {
      i_1 = i;
      break;
    }
  }
  if(i_1 < 0) i_1 = 0;
  for(i=0; i<UPP_N_RHUM; i++)
  {
    if(upp_rhum[i] >= rh)
    {
      i_2 = i;
      break;
    }
  }
  if(i_2 < 0) i_2 = UPP_N_RHUM-1;
  if(i_1 != i_2)
  {
    tropo_phas_tmp1 = (double*)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
    tropo_phas_tmp2 = (double*)malloc(mie_n_wlen*mie_n_angl*sizeof(double));
    Interp2D(upp_wlen,upp_angl,upp_tropo_phas[i_1],UPP_N_WLEN,UPP_N_ANGL,mie_wlen,mie_angl,tropo_phas_tmp1,mie_n_wlen,mie_n_angl,0);
    Interp2D(upp_wlen,upp_angl,upp_tropo_phas[i_2],UPP_N_WLEN,UPP_N_ANGL,mie_wlen,mie_angl,tropo_phas_tmp2,mie_n_wlen,mie_n_angl,0);
    for(i=0; i<mie_n_wlen; i++)
    {
      for(j=0; j<mie_n_angl; j++)
      {
        k = mie_n_angl*i+j;
        mie_tropo_phas[k] = INTERP(rh,upp_rhum[i_1],upp_rhum[i_2],tropo_phas_tmp1[k],tropo_phas_tmp2[k]);
      }
    }
    free(tropo_phas_tmp1);
    free(tropo_phas_tmp2);
  }
  else
  {
    Interp2D(upp_wlen,upp_angl,upp_tropo_phas[i_1],UPP_N_WLEN,UPP_N_ANGL,mie_wlen,mie_angl,mie_tropo_phas,mie_n_wlen,mie_n_angl,0);
  }

  for(i=0; i<mie_n_wlen; i++)
  {
    p1 = p2 = p3 = 0.0;
    for(j=1; j<mie_n_angl; j++)
    {
      k = mie_n_angl*i+j;
      p1 += (mie_tropo_phas[k-1]*mie_angl_sin[j-1]+mie_tropo_phas[k]*mie_angl_sin[j])*mie_angl_dif[j];
      p2 += (mie_strat_phas[k-1]*mie_angl_sin[j-1]+mie_strat_phas[k]*mie_angl_sin[j])*mie_angl_dif[j];
      p3 += (mie_meteo_phas[k-1]*mie_angl_sin[j-1]+mie_meteo_phas[k]*mie_angl_sin[j])*mie_angl_dif[j];
    }
    for(j=0; j<mie_n_angl; j++)
    {
      k = mie_n_angl*i+j;
      mie_tropo_phas[k] /= p1;
      mie_strat_phas[k] /= p2;
      mie_meteo_phas[k] /= p3;
    }
  }

  Interp2D(mie_wlen,mie_angl,mie_tropo_phas,mie_n_wlen,mie_n_angl,mie_wlen,sim_angl,sim_tropo_phas,mie_n_wlen,sim_n_angl,1);
  Interp2D(mie_wlen,mie_angl,mie_strat_phas,mie_n_wlen,mie_n_angl,mie_wlen,sim_angl,sim_strat_phas,mie_n_wlen,sim_n_angl,0);
  Interp2D(mie_wlen,mie_angl,mie_meteo_phas,mie_n_wlen,mie_n_angl,mie_wlen,sim_angl,sim_meteo_phas,mie_n_wlen,sim_n_angl,0);

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
    fprintf(stderr,"SepToAzi: error, th=%13.6e, sep=%13.6e\n",th,sep);
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

int Interp2D(const double *x1,const double *y1,const double *z1,int nx1,int ny1,
             const double *x2,const double *y2,      double *z2,int nx2,int ny2,int f)
{
  // x1[nx1],y1[ny1],z1[nx1*ny1]
  // x2[nx2],y2[ny2],z2[nx2*ny2]
  // x1,y1 must be arranged in ascending order
  int i,j,k;
  static int flag = 1;
  static int *i_1 = NULL;
  static int *i_2 = NULL;
  static int *j_1 = NULL;
  static int *j_2 = NULL;
  int k_1,k_2;
  double p1,p2;

  if(flag==1 || f==1)
  {
    if(i_1 != NULL) free(i_1);
    if(i_2 != NULL) free(i_2);
    if(j_1 != NULL) free(j_1);
    if(j_2 != NULL) free(j_2);
    i_1 = (int*)malloc(nx2*sizeof(int));
    i_2 = (int*)malloc(nx2*sizeof(int));
    j_1 = (int*)malloc(ny2*sizeof(int));
    j_2 = (int*)malloc(ny2*sizeof(int));
    if(i_1==NULL || i_2==NULL || j_1==NULL || j_2==NULL)
    {
      fprintf(stderr,"Interp2D: error in allocating memory\n");
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
        fprintf(stderr,"Interp2D: i_1[%d]=%d, i_2[%d]=%d\n",i,i_1[i],i,i_2[i]);
      }
    }
    for(i=0; i<ny2; i++)
    {
      j_1[i] = j_2[i] = -1;
      for(j=ny1-1; j>=0; j--)
      {
        if(y1[j] <= y2[i])
        {
          j_1[i] = j;
          break;
        }
      }
      for(j=0; j<ny1; j++)
      {
        if(y1[j] >= y2[i])
        {
          j_2[i] = j;
          break;
        }
      }
      if(j_1[i] < 0)
      {
        if(ny1 < 2)
        {
          j_1[i] = 0;
          j_2[i] = 0;
        }
        else
        {
          j_1[i] = 0;
          j_2[i] = 1;
        }
      }
      if(j_2[i] < 0)
      {
        if(ny1 < 2)
        {
          j_1[i] = 0;
          j_2[i] = 0;
        }
        else
        {
          j_1[i] = ny1-2;
          j_2[i] = ny1-1;
        }
      }
      if(cnt_db)
      {
        fprintf(stderr,"Interp2D: j_1[%d]=%d, j_2[%d]=%d\n",i,j_1[i],i,j_2[i]);
      }
    }
    flag = 0;
  }

  for(i=0; i<nx2; i++)
  {
    if(fabs(x1[i_2[i]]-x1[i_1[i]]) < EPSILON)
    {
      for(j=0; j<ny2; j++)
      {
        k = ny2*i+j;
        if(fabs(y1[j_2[j]]-y1[j_1[j]]) < EPSILON)
        {
          // NO INTERPOLATION IS NECESSARY
          k_1 = ny1*i_1[i]+j_1[j];
          z2[k] = z1[k_1];
        }
        else
        {
          // ONLY Y INTERPOLATION IS NECESSARY
          k_1 = ny1*i_1[i]+j_1[j];
          k_2 = ny1*i_1[i]+j_2[j];
          z2[k] = INTERP(y2[j],y1[j_1[j]],y1[j_2[j]],z1[k_1],z1[k_2]);
        }
      }
    }
    else
    {
      for(j=0; j<ny2; j++)
      {
        k = ny2*i+j;
        if(fabs(y1[j_2[j]]-y1[j_1[j]]) < EPSILON)
        {
          // ONLY X INTERPOLATION IS NECESSARY
          k_1 = ny1*i_1[i]+j_1[j];
          k_2 = ny1*i_2[i]+j_1[j];
          z2[k] = INTERP(x2[i],x1[i_1[i]],x1[i_2[i]],z1[k_1],z1[k_2]);
        }
        else
        {
          // INTERPOLATE IN X THEN Y
          k_1 = ny1*i_1[i]+j_1[j];
          k_2 = ny1*i_2[i]+j_1[j];
          p1 = INTERP(x2[i],x1[i_1[i]],x1[i_2[i]],z1[k_1],z1[k_2]);
          k_1 = ny1*i_1[i]+j_2[j];
          k_2 = ny1*i_2[i]+j_2[j];
          p2 = INTERP(x2[i],x1[i_1[i]],x1[i_2[i]],z1[k_1],z1[k_2]);
          z2[k] = INTERP(y2[j],y1[j_1[j]],y1[j_2[j]],p1,p2);
        }
      }
    }
  }

  return 0;
}

int WriteTape5(int rng,double *par,int mode)
{
  int m,n;
  int m1,m2,m3;
  int n1,n2;
  int err;
  double ph;
  FILE *fp;
  char fnam[] = "WriteTape5";

  if((fp=fopen("modtran.tp5","w")) == NULL)
  {
    fprintf(stderr,"Error, cannot open %s\n","modtran.tp5");
    return -1;
  }

  sim_im = (sim_2c>0?1:0);
  err = 0;
  switch(mode)
  {
    case SIM_MODE_DSR_MODTRAN_BO: // DSR simulation with MODTRAN aerosol model
      // par[0] = min wavelength
      // par[1] = max wavelength
      // par[2] = solar zenith
      // CARD 1
      WriteCard1(fp,3);
      // CARD 2
      cprintf(fp,"  %3d    0   %2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      // CARD 3
      fprintf(fp,"%10.4f     0.000%10.6f%5d           .000    0      .000\n",sim_galt,par[2],sim_iday);
      // CARD 4
      fprintf(fp,"%10.0f%10.0f%10.1f%10.1frn\n",1.0e7/par[1],1.0e7/par[0],1.0*sim_dbmp,2.0*sim_dbmp);
      if(sim_salb > 0)
      {
        // CARD 4A
        fprintf(fp,"1      0.0\n");
        // CARD 4L1
        fprintf(fp,"%s\n",sim_falb);
        // CARD 4L2
        fprintf(fp,"%d\n",sim_salb);
      }
      // CARD 5
      fprintf(fp,"    0\n");
      break;
    case SIM_MODE_DSR_MODTRAN_BA: // DSR simulation with MODTRAN aerosol model
      // par[0] = min wavelength
      // par[1] = max wavelength
      // CARD 1
      WriteCard1(fp,3);
      // CARD 2
      cprintf(fp,"  %3d    0   %2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      for(n=0; n<dsr_n_data; n++)
      {
        // CARD 3
        fprintf(fp,"%10.4f     0.000%10.6f%5d           .000    0      .000\n",sim_galt,dsr_th_sun[n],sim_iday);
        if(n == 0)
        {
          // CARD 4
          fprintf(fp,"%10.0f%10.0f%10.1f%10.1frn\n",1.0e7/par[1],1.0e7/par[0],1.0*sim_dbmp,2.0*sim_dbmp);
        }
        if(sim_salb > 0)
        {
          // CARD 4A
          fprintf(fp,"1      0.0\n");
          // CARD 4L1
          fprintf(fp,"%s\n",sim_falb);
          // CARD 4L2
          fprintf(fp,"%d\n",sim_salb);
        }
        // CARD 5
        if(n == dsr_n_data-1)
        {
          fprintf(fp,"    0\n");
        }
        else
        {
          fprintf(fp,"    3\n");
        }
      }
      break;
    case SIM_MODE_DSR_CUSTOM_BO: // DSR simulation with custom aerosol model
      // par[0] = min wavelength
      // par[1] = max wavelength
      // par[2] = solar zenith
      // CARD 1
      WriteCard1(fp,3);
      // CARD 2
      cprintf(fp,"  %3d    0USS%2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      // CARD 2D
      WriteCard2D(fp);
      // CARD 3
      fprintf(fp,"%10.4f     0.000%10.6f%5d           .000    0      .000\n",sim_galt,par[2],sim_iday);
      // CARD 4
      fprintf(fp,"%10.0f%10.0f%10.1f%10.1frn\n",1.0e7/par[1],1.0e7/par[0],1.0*sim_dbmp,2.0*sim_dbmp);
      if(sim_salb > 0)
      {
        // CARD 4A
        fprintf(fp,"1      0.0\n");
        // CARD 4L1
        fprintf(fp,"%s\n",sim_falb);
        // CARD 4L2
        fprintf(fp,"%d\n",sim_salb);
      }
      // CARD 5
      fprintf(fp,"    0\n");
      break;
    case SIM_MODE_DSR_CUSTOM_BA: // DSR simulation with custom aerosol model
      // par[0] = min wavelength
      // par[1] = max wavelength
      // CARD 1
      WriteCard1(fp,3);
      // CARD 2
      cprintf(fp,"  %3d    0USS%2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      // CARD 2D
      WriteCard2D(fp);
      for(n=0; n<dsr_n_data; n++)
      {
        // CARD 3
        fprintf(fp,"%10.4f     0.000%10.6f%5d           .000    0      .000\n",sim_galt,dsr_th_sun[n],sim_iday);
        if(n == 0)
        {
          // CARD 4
          fprintf(fp,"%10.0f%10.0f%10.1f%10.1frn\n",1.0e7/par[1],1.0e7/par[0],1.0*sim_dbmp,2.0*sim_dbmp);
        }
        if(sim_salb > 0)
        {
          // CARD 4A
          fprintf(fp,"1      0.0\n");
          // CARD 4L1
          fprintf(fp,"%s\n",sim_falb);
          // CARD 4L2
          fprintf(fp,"%d\n",sim_salb);
        }
        // CARD 5
        if(n == dsr_n_data-1)
        {
          fprintf(fp,"    0\n");
        }
        else
        {
          fprintf(fp,"    3\n");
        }
      }
      break;
    case SIM_MODE_AUR_MODTRAN_SA: // AUR simulation with MODTRAN aerosol model
      // par[0] = #data in modtran.plt
      for(m=0; m<aur_n_wlen[rng]; m++)
      {
        sim_sbmp = aur_sbmp[rng][m];
        // CARD 1
        WriteCard1(fp,2);
        // CARD 2
        cprintf(fp,"  %3d    0   %2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
        if(sim_2c>0 && sim_im>0)
        {
          // CARD 2C
          WriteCard2C(fp);
          sim_im = 0;
        }
        for(n=0; n<aur_n_data; n++)
        {
          ph = aur_ph_sun[n]-aur_ph_los[n];
          if(ph < -180.0) ph += 360.0;
          // CARD 3
          fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,aur_th_los[n]);
          // CARD 3A1
          fprintf(fp,"    2    2%5d    0\n",sim_iday);
          // CARD 3A2
          fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,aur_th_sun[n]);
          if(n == 0)
          {
            // CARD 4
            n1 = sim_sbmp*((int)((1.0e7/(aur_wlen[rng][m]*sim_sbmp)-0.5*(par[0]-1.0))+0.5));
            n2 = sim_sbmp*((int)(par[0]+0.5))+n1-sim_sbmp;
            fprintf(fp,"%10d%10d%10.1f%10.1frn\n",n1,n2,1.0*sim_sbmp,2.0*sim_sbmp);
          }
          if(sim_salb > 0)
          {
            // CARD 4A
            fprintf(fp,"1      0.0\n");
            // CARD 4L1
            fprintf(fp,"%s\n",sim_falb);
            // CARD 4L2
            fprintf(fp,"%d\n",sim_salb);
          }
          // CARD 5
          if(m==aur_n_wlen[rng]-1 && n==aur_n_data-1)
          {
            fprintf(fp,"    0\n");
          } else
          if(n == aur_n_data-1)
          {
            fprintf(fp,"    1\n");
          }
          else
          {
            fprintf(fp,"    3\n");
          }
        }
      }
      break;
    case SIM_MODE_AUR_MODTRAN_BA: // AUR simulation with MODTRAN aerosol model
      // par[0] = wavelength margin
      // CARD 1
      WriteCard1(fp,2);
      // CARD 2
      cprintf(fp,"  %3d    0   %2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      for(n=0; n<aur_n_data; n++)
      {
        ph = aur_ph_sun[n]-aur_ph_los[n];
        if(ph < -180.0) ph += 360.0;
        if(n%2 == 0)
        {
          m1 = 0;
          m2 = aur_n_wlen[rng];
          m3 = 1;
        }
        else
        {
          m1 = aur_n_wlen[rng]-1;
          m2 = -1;
          m3 = -1;
        }
        for(m=m1; m!=m2; m+=m3)
        {
          if(m == m1)
          {
            // CARD 3
            fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,aur_th_los[n]);
            // CARD 3A1
            fprintf(fp,"    2    2%5d    0\n",sim_iday);
            // CARD 3A2
            fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,aur_th_sun[n]);
          }
          if(n==0 || m!=m1)
          {
            // CARD 4
            fprintf(fp,"%10.0f%10.0f%10.1f%10.1frn\n",1.0e7/(aur_wlen[rng][m]+par[0]),
                                                      1.0e7/(aur_wlen[rng][m]-par[0]),1.0*sim_sbmp,2.0*sim_sbmp);
          }
          if(sim_salb > 0)
          {
            // CARD 4A
            fprintf(fp,"1      0.0\n");
            // CARD 4L1
            fprintf(fp,"%s\n",sim_falb);
            // CARD 4L2
            fprintf(fp,"%d\n",sim_salb);
          }
          // CARD 5
          if(m==m2-m3 && n==aur_n_data-1)
          {
            fprintf(fp,"    0\n");
          } else
          if(m == m2-m3)
          {
            fprintf(fp,"    3\n");
          }
          else
          {
            fprintf(fp,"    4\n");
          }
        }
      }
      break;
    case SIM_MODE_AUR_CUSTOM_SA: // AUR simulation with custom aerosol model
      // par[0] = #data in modtran.plt
      for(m=0; m<aur_n_wlen[rng]; m++)
      {
        sim_sbmp = aur_sbmp[rng][m];
        // CARD 1
        WriteCard1(fp,2);
        // CARD 2
        cprintf(fp,"  %3d    0USS%2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
        if(sim_2c>0 && sim_im>0)
        {
          // CARD 2C
          WriteCard2C(fp);
          sim_im = 0;
        }
        // CARD 2D
        WriteCard2D(fp);
        for(n=0; n<aur_n_data; n++)
        {
          ph = aur_ph_sun[n]-aur_ph_los[n];
          if(ph < -180.0) ph += 360.0;
          // CARD 3
          fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,aur_th_los[n]);
          // CARD 3A1
          fprintf(fp,"    2    1%5d    0\n",sim_iday);
          // CARD 3A2
          fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,aur_th_sun[n]);
          // CARD 3B1
          fprintf(fp,"%5d%5d\n",sim_n_angl,mie_n_wlen);
          // CARD 3C
          WriteCard3C(fp);
          if(n == 0)
          {
            // CARD 4
            n1 = sim_sbmp*((int)((1.0e7/(aur_wlen[rng][m]*sim_sbmp)-0.5*(par[0]-1.0))+0.5));
            n2 = sim_sbmp*((int)(par[0]+0.5))+n1-sim_sbmp;
            fprintf(fp,"%10d%10d%10.1f%10.1frn\n",n1,n2,1.0*sim_sbmp,2.0*sim_sbmp);
          }
          if(sim_salb > 0)
          {
            // CARD 4A
            fprintf(fp,"1      0.0\n");
            // CARD 4L1
            fprintf(fp,"%s\n",sim_falb);
            // CARD 4L2
            fprintf(fp,"%d\n",sim_salb);
          }
          // CARD 5
          if(m==aur_n_wlen[rng]-1 && n==aur_n_data-1)
          {
            fprintf(fp,"    0\n");
          } else
          if(n == aur_n_data-1)
          {
            fprintf(fp,"    1\n");
          }
          else
          {
            fprintf(fp,"    3\n");
          }
        }
      }
      break;
    case SIM_MODE_AUR_CUSTOM_BA: // AUR simulation with custom aerosol model
      // par[0] = wavelength margin
      // CARD 1
      WriteCard1(fp,2);
      // CARD 2
      cprintf(fp,"  %3d    0USS%2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      // CARD 2D
      WriteCard2D(fp);
      for(n=0; n<aur_n_data; n++)
      {
        ph = aur_ph_sun[n]-aur_ph_los[n];
        if(ph < -180.0) ph += 360.0;
        if(n%2 == 0)
        {
          m1 = 0;
          m2 = aur_n_wlen[rng];
          m3 = 1;
        }
        else
        {
          m1 = aur_n_wlen[rng]-1;
          m2 = -1;
          m3 = -1;
        }
        for(m=m1; m!=m2; m+=m3)
        {
          if(m == m1)
          {
            // CARD 3
            fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,aur_th_los[n]);
            // CARD 3A1
            fprintf(fp,"    2    1%5d    0\n",sim_iday);
            // CARD 3A2
            fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,aur_th_sun[n]);
            // CARD 3B1
            fprintf(fp,"%5d%5d\n",sim_n_angl,mie_n_wlen);
            // CARD 3C
            WriteCard3C(fp);
          }
          if(n==0 || m!=m1)
          {
            // CARD 4
            fprintf(fp,"%10.0f%10.0f%10.1f%10.1frn\n",1.0e7/(aur_wlen[rng][m]+par[0]),
                                                      1.0e7/(aur_wlen[rng][m]-par[0]),1.0*sim_sbmp,2.0*sim_sbmp);
          }
          if(sim_salb > 0)
          {
            // CARD 4A
            fprintf(fp,"1      0.0\n");
            // CARD 4L1
            fprintf(fp,"%s\n",sim_falb);
            // CARD 4L2
            fprintf(fp,"%d\n",sim_salb);
          }
          // CARD 5
          if(m==m2-m3 && n==aur_n_data-1)
          {
            fprintf(fp,"    0\n");
          } else
          if(m == m2-m3)
          {
            fprintf(fp,"    3\n");
          }
          else
          {
            fprintf(fp,"    4\n");
          }
        }
      }
      break;
    case SIM_MODE_SSR_MODTRAN_SO: // SSR simulation with MODTRAN aerosol model
      // par[0] = wavelength
      // par[1] = #data in modtran.plt
      // par[2] = solar zenith
      // par[3] = solar azimuth
      // par[4] = los   zenith
      // par[5] = los   azimuth
      // CARD 1
      WriteCard1(fp,2);
      // CARD 2
      cprintf(fp,"  %3d    0   %2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      ph = par[3]-par[5];
      if(ph < -180.0) ph += 360.0;
      // CARD 3
      fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,par[4]);
      // CARD 3A1
      fprintf(fp,"    2    2%5d    0\n",sim_iday);
      // CARD 3A2
      fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,par[2]);
      // CARD 4
      n1 = sim_sbmp*((int)((1.0e7/(par[0]*sim_sbmp)-0.5*(par[1]-1.0))+0.5));
      n2 = sim_sbmp*((int)(par[1]+0.5))+n1-sim_sbmp;
      fprintf(fp,"%10d%10d%10.1f%10.1frn\n",n1,n2,1.0*sim_sbmp,2.0*sim_sbmp);
      if(sim_salb > 0)
      {
        // CARD 4A
        fprintf(fp,"1      0.0\n");
        // CARD 4L1
        fprintf(fp,"%s\n",sim_falb);
        // CARD 4L2
        fprintf(fp,"%d\n",sim_salb);
      }
      // CARD 5
      fprintf(fp,"    0\n");
      break;
    case SIM_MODE_SSR_MODTRAN_SA: // SSR simulation with MODTRAN aerosol model
      // par[0] = #data in modtran.plt
      for(m=0; m<ssr_n_wlen[rng]; m++)
      {
        sim_sbmp = ssr_sbmp[rng][m];
        // CARD 1
        WriteCard1(fp,2);
        // CARD 2
        cprintf(fp,"  %3d    0   %2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
        if(sim_2c>0 && sim_im>0)
        {
          // CARD 2C
          WriteCard2C(fp);
          sim_im = 0;
        }
        for(n=0; n<ssr_n_data; n++)
        {
          ph = ssr_ph_sun[n]-ssr_ph_los[n];
          if(ph < -180.0) ph += 360.0;
          // CARD 3
          fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,ssr_th_los[n]);
          // CARD 3A1
          fprintf(fp,"    2    2%5d    0\n",sim_iday);
          // CARD 3A2
          fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,ssr_th_sun[n]);
          if(n == 0)
          {
            // CARD 4
            n1 = sim_sbmp*((int)((1.0e7/(ssr_wlen[rng][m]*sim_sbmp)-0.5*(par[0]-1.0))+0.5));
            n2 = sim_sbmp*((int)(par[0]+0.5))+n1-sim_sbmp;
            fprintf(fp,"%10d%10d%10.1f%10.1frn\n",n1,n2,1.0*sim_sbmp,2.0*sim_sbmp);
          }
          if(sim_salb > 0)
          {
            // CARD 4A
            fprintf(fp,"1      0.0\n");
            // CARD 4L1
            fprintf(fp,"%s\n",sim_falb);
            // CARD 4L2
            fprintf(fp,"%d\n",sim_salb);
          }
          // CARD 5
          if(m==ssr_n_wlen[rng]-1 && n==ssr_n_data-1)
          {
            fprintf(fp,"    0\n");
          } else
          if(n == ssr_n_data-1)
          {
            fprintf(fp,"    1\n");
          }
          else
          {
            fprintf(fp,"    3\n");
          }
        }
      }
      break;
    case SIM_MODE_SSR_MODTRAN_BO: // SSR simulation with MODTRAN aerosol model
      // par[0] = min wavelength
      // par[1] = max wavelength
      // par[2] = solar zenith
      // par[3] = solar azimuth
      // par[4] = los   zenith
      // par[5] = los   azimuth
      // CARD 1
      WriteCard1(fp,2);
      // CARD 2
      cprintf(fp,"  %3d    0   %2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      ph = par[3]-par[5];
      if(ph < -180.0) ph += 360.0;
      // CARD 3
      fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,par[4]);
      // CARD 3A1
      fprintf(fp,"    2    2%5d    0\n",sim_iday);
      // CARD 3A2
      fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,par[2]);
      // CARD 4
      fprintf(fp,"%10.0f%10.0f%10.1f%10.1frn\n",1.0e7/par[1],1.0e7/par[0],1.0*sim_sbmp,2.0*sim_sbmp);
      if(sim_salb > 0)
      {
        // CARD 4A
        fprintf(fp,"1      0.0\n");
        // CARD 4L1
        fprintf(fp,"%s\n",sim_falb);
        // CARD 4L2
        fprintf(fp,"%d\n",sim_salb);
      }
      // CARD 5
      fprintf(fp,"    0\n");
      break;
    case SIM_MODE_SSR_MODTRAN_BA: // SSR simulation with MODTRAN aerosol model
      // par[0] = wavelength margin
      // CARD 1
      WriteCard1(fp,2);
      // CARD 2
      cprintf(fp,"  %3d    0   %2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      for(n=0; n<ssr_n_data; n++)
      {
        ph = ssr_ph_sun[n]-ssr_ph_los[n];
        if(ph < -180.0) ph += 360.0;
        if(n%2 == 0)
        {
          m1 = 0;
          m2 = ssr_n_wlen[rng];
          m3 = 1;
        }
        else
        {
          m1 = ssr_n_wlen[rng]-1;
          m2 = -1;
          m3 = -1;
        }
        for(m=m1; m!=m2; m+=m3)
        {
          if(m == m1)
          {
            // CARD 3
            fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,ssr_th_los[n]);
            // CARD 3A1
            fprintf(fp,"    2    2%5d    0\n",sim_iday);
            // CARD 3A2
            fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,ssr_th_sun[n]);
          }
          if(n==0 || m!=m1)
          {
            // CARD 4
            fprintf(fp,"%10.0f%10.0f%10.1f%10.1frn\n",1.0e7/(ssr_wlen[rng][m]+par[0]),
                                                      1.0e7/(ssr_wlen[rng][m]-par[0]),1.0*sim_sbmp,2.0*sim_sbmp);
          }
          if(sim_salb > 0)
          {
            // CARD 4A
            fprintf(fp,"1      0.0\n");
            // CARD 4L1
            fprintf(fp,"%s\n",sim_falb);
            // CARD 4L2
            fprintf(fp,"%d\n",sim_salb);
          }
          // CARD 5
          if(m==m2-m3 && n==ssr_n_data-1)
          {
            fprintf(fp,"    0\n");
          } else
          if(m == m2-m3)
          {
            fprintf(fp,"    3\n");
          }
          else
          {
            fprintf(fp,"    4\n");
          }
        }
      }
      break;
    case SIM_MODE_SSR_CUSTOM_SO: // SSR simulation with custom aerosol model
      // par[0] = wavelength
      // par[1] = #data in modtran.plt
      // par[2] = solar zenith
      // par[3] = solar azimuth
      // par[4] = los   zenith
      // par[5] = los   azimuth
      // CARD 1
      WriteCard1(fp,2);
      // CARD 2
      cprintf(fp,"  %3d    0USS%2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      // CARD 2D
      WriteCard2D(fp);
      ph = par[3]-par[5];
      if(ph < -180.0) ph += 360.0;
      // CARD 3
      fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,par[4]);
      // CARD 3A1
      fprintf(fp,"    2    1%5d    0\n",sim_iday);
      // CARD 3A2
      fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,par[2]);
      // CARD 3B1
      fprintf(fp,"%5d%5d\n",sim_n_angl,mie_n_wlen);
      // CARD 3C
      WriteCard3C(fp);
      // CARD 4
      n1 = sim_sbmp*((int)((1.0e7/(par[0]*sim_sbmp)-0.5*(par[1]-1.0))+0.5));
      n2 = sim_sbmp*((int)(par[1]+0.5))+n1-sim_sbmp;
      fprintf(fp,"%10d%10d%10.1f%10.1frn\n",n1,n2,1.0*sim_sbmp,2.0*sim_sbmp);
      if(sim_salb > 0)
      {
        // CARD 4A
        fprintf(fp,"1      0.0\n");
        // CARD 4L1
        fprintf(fp,"%s\n",sim_falb);
        // CARD 4L2
        fprintf(fp,"%d\n",sim_salb);
      }
      // CARD 5
      fprintf(fp,"    0\n");
      break;
    case SIM_MODE_SSR_CUSTOM_SA: // SSR simulation with custom aerosol model
      // par[0] = #data in modtran.plt
      for(m=0; m<ssr_n_wlen[rng]; m++)
      {
        sim_sbmp = ssr_sbmp[rng][m];
        // CARD 1
        WriteCard1(fp,2);
        // CARD 2
        cprintf(fp,"  %3d    0USS%2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
        if(sim_2c>0 && sim_im>0)
        {
          // CARD 2C
          WriteCard2C(fp);
          sim_im = 0;
        }
        // CARD 2D
        WriteCard2D(fp);
        for(n=0; n<ssr_n_data; n++)
        {
          ph = ssr_ph_sun[n]-ssr_ph_los[n];
          if(ph < -180.0) ph += 360.0;
          // CARD 3
          fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,ssr_th_los[n]);
          // CARD 3A1
          fprintf(fp,"    2    1%5d    0\n",sim_iday);
          // CARD 3A2
          fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,ssr_th_sun[n]);
          // CARD 3B1
          fprintf(fp,"%5d%5d\n",sim_n_angl,mie_n_wlen);
          // CARD 3C
          WriteCard3C(fp);
          if(n == 0)
          {
            // CARD 4
            n1 = sim_sbmp*((int)((1.0e7/(ssr_wlen[rng][m]*sim_sbmp)-0.5*(par[0]-1.0))+0.5));
            n2 = sim_sbmp*((int)(par[0]+0.5))+n1-sim_sbmp;
            fprintf(fp,"%10d%10d%10.1f%10.1frn\n",n1,n2,1.0*sim_sbmp,2.0*sim_sbmp);
          }
          if(sim_salb > 0)
          {
            // CARD 4A
            fprintf(fp,"1      0.0\n");
            // CARD 4L1
            fprintf(fp,"%s\n",sim_falb);
            // CARD 4L2
            fprintf(fp,"%d\n",sim_salb);
          }
          // CARD 5
          if(m==ssr_n_wlen[rng]-1 && n==ssr_n_data-1)
          {
            fprintf(fp,"    0\n");
          } else
          if(n == ssr_n_data-1)
          {
            fprintf(fp,"    1\n");
          }
          else
          {
            fprintf(fp,"    3\n");
          }
        }
      }
      break;
    case SIM_MODE_SSR_CUSTOM_BO: // SSR simulation with custom aerosol model
      // par[0] = min wavelength
      // par[1] = max wavelength
      // par[2] = solar zenith
      // par[3] = solar azimuth
      // par[4] = los   zenith
      // par[5] = los   azimuth
      // CARD 1
      WriteCard1(fp,2);
      // CARD 2
      cprintf(fp,"  %3d    0USS%2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      // CARD 2D
      WriteCard2D(fp);
      ph = par[3]-par[5];
      if(ph < -180.0) ph += 360.0;
      // CARD 3
      fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,par[4]);
      // CARD 3A1
      fprintf(fp,"    2    1%5d    0\n",sim_iday);
      // CARD 3A2
      fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,par[2]);
      // CARD 3B1
      fprintf(fp,"%5d%5d\n",sim_n_angl,mie_n_wlen);
      // CARD 3C
      WriteCard3C(fp);
      // CARD 4
      fprintf(fp,"%10.0f%10.0f%10.1f%10.1frn\n",1.0e7/par[1],1.0e7/par[0],1.0*sim_sbmp,2.0*sim_sbmp);
      if(sim_salb > 0)
      {
        // CARD 4A
        fprintf(fp,"1      0.0\n");
        // CARD 4L1
        fprintf(fp,"%s\n",sim_falb);
        // CARD 4L2
        fprintf(fp,"%d\n",sim_salb);
      }
      // CARD 5
      fprintf(fp,"    0\n");
      break;
    case SIM_MODE_SSR_CUSTOM_BA: // SSR simulation with custom aerosol model
      // par[0] = wavelength margin
      // CARD 1
      WriteCard1(fp,2);
      // CARD 2
      cprintf(fp,"  %3d    0USS%2d    0    0    0%10.5f      .000      .000      .000      .000\n",sim_ihaz,sim_ivul,sim_vmie);
      if(sim_2c>0 && sim_im>0)
      {
        // CARD 2C
        WriteCard2C(fp);
      }
      // CARD 2D
      WriteCard2D(fp);
      for(n=0; n<ssr_n_data; n++)
      {
        ph = ssr_ph_sun[n]-ssr_ph_los[n];
        if(ph < -180.0) ph += 360.0;
        if(n%2 == 0)
        {
          m1 = 0;
          m2 = ssr_n_wlen[rng];
          m3 = 1;
        }
        else
        {
          m1 = ssr_n_wlen[rng]-1;
          m2 = -1;
          m3 = -1;
        }
        for(m=m1; m!=m2; m+=m3)
        {
          if(m == m1)
          {
            // CARD 3
            fprintf(fp,"%10.4f     0.000%10.6f      .000      .000      .000    0        0.00000\n",sim_galt,ssr_th_los[n]);
            // CARD 3A1
            fprintf(fp,"    2    1%5d    0\n",sim_iday);
            // CARD 3A2
            fprintf(fp,"%10.5f%10.6f      .000      .000      .000      .000      .000      .000\n",ph,ssr_th_sun[n]);
            // CARD 3B1
            fprintf(fp,"%5d%5d\n",sim_n_angl,mie_n_wlen);
            // CARD 3C
            WriteCard3C(fp);
          }
          if(n==0 || m!=m1)
          {
            // CARD 4
            fprintf(fp,"%10.0f%10.0f%10.1f%10.1frn\n",1.0e7/(ssr_wlen[rng][m]+par[0]),
                                                      1.0e7/(ssr_wlen[rng][m]-par[0]),1.0*sim_sbmp,2.0*sim_sbmp);
          }
          if(sim_salb > 0)
          {
            // CARD 4A
            fprintf(fp,"1      0.0\n");
            // CARD 4L1
            fprintf(fp,"%s\n",sim_falb);
            // CARD 4L2
            fprintf(fp,"%d\n",sim_salb);
          }
          // CARD 5
          if(m==m2-m3 && n==ssr_n_data-1)
          {
            fprintf(fp,"    0\n");
          } else
          if(m == m2-m3)
          {
            fprintf(fp,"    3\n");
          }
          else
          {
            fprintf(fp,"    4\n");
          }
        }
      }
      break;
    default:
      fprintf(stderr,"%s: invalid mode >>> %d\n",fnam,mode);
      err = 1;
      break;
  }
  fclose(fp);
  if(err)
  {
    return -1;
  }

  return 0;
}

int WriteCard1(FILE *fp,int mode)
{
  if(mode == 2)
  {
    // CARD 1
    if(sim_salb > 0)
    {
      fprintf(fp,"MM%3d    3    2%5d    0    0    0    0    0    0    0%5d%5d    .000 LAMBER\n",
                  (sim_2c>0?7:sim_iatm),sim_mult,sim_im,sim_prnt);
    }
    else
    {
      fprintf(fp,"MM%3d    3    2%5d    0    0    0    0    0    0    0%5d%5d    .000%7.4f\n",
                  (sim_2c>0?7:sim_iatm),sim_mult,sim_im,sim_prnt,sim_grho);
    }
    // CARD 1A
    cprintf(fp,"%st 16f   0%10.3f%10.0f%10.0f f t f %s\n",sim_dflg==0?"f":"t",
                sim_mixr_co2,sim_wscl,sim_oscl,sim_wflg==1?"t":"f");
    // CARD 1A2
    fprintf(fp,"DATA/B2001_%02d.BIN\n",sim_sbmp);
  } else
  if(mode == 3)
  {
    // CARD 1
    if(sim_salb > 0)
    {
      fprintf(fp,"MM%3d    3    3%5d    0    0    0    0    0    0    0%5d%5d    .000 LAMBER\n",
                  (sim_2c>0?7:sim_iatm),sim_mult,sim_im,sim_prnt);
    }
    else
    {
      fprintf(fp,"MM%3d    3    3%5d    0    0    0    0    0    0    0%5d%5d    .000%7.4f\n",
                  (sim_2c>0?7:sim_iatm),sim_mult,sim_im,sim_prnt,sim_grho);
    }
    // CARD 1A
    cprintf(fp,"tt  8f   0%10.3f%10.0f%10.0f f t f %s\n",
                sim_mixr_co2,sim_wscl,sim_oscl,sim_wflg==1?"t":"f");
    // CARD 1A2
    fprintf(fp,"DATA/B2001_%02d.BIN\n",sim_dbmp);
  }
  else
  {
    fprintf(stderr,"WriteCard1: error, mode=%d\n",mode);
    return -1;
  }

  return 0;
}

int WriteCard2C(FILE *fp)
{
  int i;
  int f2,f3;
  char jchar[15] = "";
  char jcharx[2] = " ";

  if(sim_2c <= 0)
  {
    return 0;
  }
  else
  {
    f2 = (sim_2c&0x0010);
    f3 = (sim_2c&0x1000);
  }

  for(i=0; i<PRF_NPAR; i++)
  {
    jchar[i] = sim_prf_jchar[i];
  }
  jchar[i] = '\0';
  jcharx[0] = sim_prf_jcharx;

  // CARD 2C
  fprintf(fp,"%5d%5d%5d%20s%10.0f\n",prf_n_alt,(f2>0?1:0),(f3>0?2:0),"Custom",0.0);
  for(i=0; i<prf_n_alt; i++)
  {
    // CARD 2C1
    fprintf(fp,"%10.3f%10.3e%10.3e%10.3e%10.3e%10.3e%14s %1s\n",
                prf_alt[i],sim_prf_pres[i],sim_prf_temp[i],
                sim_prf_wmol[0][i],sim_prf_wmol[1][i],sim_prf_wmol[2][i],jchar,jcharx);
    if(f2 > 0)
    {
      // CARD 2C2
      fprintf(fp,"%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e\n",
                  sim_prf_wmol[3][i],sim_prf_wmol[4][i],sim_prf_wmol[5][i],sim_prf_wmol[6][i],
                  sim_prf_wmol[7][i],sim_prf_wmol[8][i],sim_prf_wmol[9][i],sim_prf_wmol[10][i]);
      fprintf(fp,"%10.3e\n",sim_prf_wmol[11][i]);
    }
    if(f3 > 0)
    {
      // CARD 2C3
      cprintf(fp,"%10s%10.0f%10s%10.0f%10.0f%10.0f%10.0f\n","",sim_prf_haze[0][i],"",
                  0.0,sim_prf_haze[1][i],sim_prf_haze[2][i],sim_prf_haze[3][i]);
    }
  }

  return 0;
}

int WriteCard2D(FILE *fp)
{
  int i;

  // CARD 2D
  fprintf(fp,"%5d%5d%5d%5d\n",mie_n_wlen,0,0,0);
  // CARD 2D1
  fprintf(fp,"%10.3e%70s\n",0.0,"title");
  // CARD 2D2
  for(i=0; i<mie_n_wlen; i++)
  {
    cprintf(fp,"%6.2f%7.5f%7.5f%6.4f%s",mie_wlen_um[i],mie_aext[i]/mie_aext[mie_iref],
                                        (mie_aext[i]-mie_asca[i])/mie_aext[mie_iref],mie_asym[i],
                                        i==mie_n_wlen-1?"\n":(i%3)==2?"\n":"");
  }

  return 0;
}

int WriteCard3C(FILE *fp)
{
  int i,j,k;

  // CARD 3C1
  for(i=0; i<sim_n_angl; i++)
  {
    fprintf(fp,"%10.5f%s",sim_angl[i],i==sim_n_angl-1?"\n":i%8==7?"\n":"");
  }
  // CARD 3C2
  for(i=0; i<mie_n_wlen; i++)
  {
    fprintf(fp,"%10.6f%s",mie_wlen_um[i],i==mie_n_wlen-1?"\n":i%8==7?"\n":"");
  }
  // CARD 3C3
  for(j=0; j<sim_n_angl; j++)
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      k = sim_n_angl*i+j;
      cprintf(fp," %9.3e%s",sim_phas[k],i==mie_n_wlen-1?"\n":i%8==7?"\n":"");
    }
  }
  // CARD 3C4
  for(j=0; j<sim_n_angl; j++)
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      k = sim_n_angl*i+j;
      fprintf(fp,"%10.3e%s",sim_tropo_phas[k],i==mie_n_wlen-1?"\n":i%8==7?"\n":"");
    }
  }
  // CARD 3C5
  for(j=0; j<sim_n_angl; j++)
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      k = sim_n_angl*i+j;
      fprintf(fp,"%10.3e%s",sim_strat_phas[k],i==mie_n_wlen-1?"\n":i%8==7?"\n":"");
    }
  }
  // CARD 3C6
  for(j=0; j<sim_n_angl; j++)
  {
    for(i=0; i<mie_n_wlen; i++)
    {
      k = sim_n_angl*i+j;
      fprintf(fp,"%10.3e%s",sim_meteo_phas[k],i==mie_n_wlen-1?"\n":i%8==7?"\n":"");
    }
  }

  return 0;
}

int RunModtran(void)
{
  char line[MAXLINE];
  FILE *fp;

  if((fp=fopen("modroot.in","w")) == NULL)
  {
    fprintf(stderr,"Error, cannot open %s\n","modroot.in");
    return -1;
  }
  fprintf(fp,"modtran\n");
  fclose(fp);
  if(access("DATA",R_OK) < 0)
  {
    snprintf(line,MAXLINE,"%s/DATA",sim_path);
    if(symlink(line,"DATA") < 0)
    {
      fprintf(stderr,"Error, cannot create symlink to %s\n",line);
      return -1;
    }
  }
  snprintf(line,MAXLINE,"%s/Mod4v2r1_F90.exe",sim_path);
  system(line);

  return 0;
}

int ReadPlt_S(const double *par,double *y)
{
  // par[0] = wavelength
  // par[1] = #data in modtran.plt
  int i;
  int n1,n2;
  int err;
  int flag;
  double s;
  double v[2];
  char line[MAXLINE];
  char str1[MAXLINE];
  char str2[MAXLINE];
  char fnam[] = "ReadPlt_S";
  char *p;
  FILE *fp;

  if((fp=fopen("modtran.plt","r")) == NULL)
  {
    fprintf(stderr,"%s: error, cannot open %s\n",fnam,"modtran.plt");
    return -1;
  }
  n1 = (int)(par[1]+0.5);
  i = 0;
  s = 0.0;
  err = 0;
  flag = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    if(sscanf(line,"%s%s",str1,str2) != 2)
    {
      flag = 1;
      break;
    }
    errno = 0;
    v[0] = strtod(str1,&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
      err = 1;
      break;
    }
    errno = 0;
    v[1] = strtod(str2,&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
      err = 1;
      break;
    }
    if(fabs(v[0]-par[0]) < SIM_DMAX)
    {
      s += v[1];
      i += 1;
    }
  }
  n2 = i;
  if(flag == 0)
  {
    fprintf(stderr,"%s: error, no blank line\n",fnam);
    err = 1;
  }
  if(n2 != n1)
  {
    fprintf(stderr,"%s: error, n1=%d, n2=%d\n",fnam,n1,n2);
    err = 1;
  }
  *y = s*sim_yuni/n1;
  fclose(fp);
  if(err)
  {
    return -1;
  }

  return 0;
}

int ReadPlt_B(const double *par,int n,const double *x,double *y)
{
  // par[0] = min wavelength
  // par[1] = max wavelength
  int i;
  int n1,n2;
  int err;
  int flag;
  double v[2];
  char line[MAXLINE];
  char str1[MAXLINE];
  char str2[MAXLINE];
  char fnam[] = "ReadPlt_B";
  char *p;
  FILE *fp;

  if((fp=fopen("modtran.plt","r")) == NULL)
  {
    fprintf(stderr,"%s: error, cannot open %s\n",fnam,"modtran.plt");
    return -1;
  }
  i = 0;
  err = 0;
  flag = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    if(sscanf(line,"%s%s",str1,str2) != 2)
    {
      flag = 1;
      break;
    }
    errno = 0;
    v[0] = strtod(str1,&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
      err = 1;
      break;
    }
    errno = 0;
    v[1] = strtod(str2,&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
      err = 1;
      break;
    }
    if(i >= SIM_MAXDATA)
    {
      fprintf(stderr,"%s: error, #data exceed the limit >>> %d\n",fnam,i);
      err = 1;
      break;
    }
    sim_wlen[i] = v[0];
    sim_data[i] = v[1];
    i++;
  }
  n1 = i;
  if(flag == 0)
  {
    fprintf(stderr,"%s: error, no blank line\n",fnam);
    err = 1;
  }
  if(n1 < 1)
  {
    fprintf(stderr,"%s: no data\n",fnam);
    err = 1;
  }
  if((n2=Convolute(par[0],par[1],1.0,sim_xsgm,sim_wsgm,
                   n1,sim_wlen,sim_data,SIM_MAXCONV,sim_wlen_cnv,sim_data_cnv)) < 0)
  {
    err = 1;
  }
  if(Sampling(n2,sim_wlen_cnv,sim_data_cnv,n,x,y,sim_yuni) < 0)
  {
    err = 1;
  }
  fclose(fp);
  if(err)
  {
    return -1;
  }

  return 0;
}

int ReadDSR_B(int rng,const double *par)
{
  // par[0] = min wavelength
  // par[1] = max wavelength
  int i,n;
  int n1,n2;
  int err;
  int flag;
  double v[2];
  char line[MAXLINE];
  char str1[MAXLINE];
  char str2[MAXLINE];
  char fnam[] = "ReadDSR_B";
  char *p;
  FILE *fp;

  if((fp=fopen("modtran.plt","r")) == NULL)
  {
    fprintf(stderr,"%s: error, cannot open %s\n",fnam,"modtran.plt");
    return -1;
  }
  for(n=0; n<dsr_n_data; n++)
  {
    i = 0;
    err = 0;
    flag = 0;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      if(sscanf(line,"%s%s",str1,str2) != 2)
      {
        flag = 1;
        break;
      }
      errno = 0;
      v[0] = strtod(str1,&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      errno = 0;
      v[1] = strtod(str2,&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      if(i >= SIM_MAXDATA)
      {
        fprintf(stderr,"%s: error, #data exceed the limit >>> %d\n",fnam,i);
        err = 1;
        break;
      }
      sim_wlen[i] = v[0];
      sim_data[i] = v[1];
      i++;
    }
    n1 = i;
    if(flag == 0)
    {
      fprintf(stderr,"%s: error, no blank line for %d\n",fnam,n);
      err = 1;
    }
    if(n1 < 1)
    {
      fprintf(stderr,"%s: no data\n",fnam);
      err = 1;
    }
    if(err) break;
    if((n2=Convolute(par[0],par[1],1.0,sim_xsgm,sim_wsgm,
                     n1,sim_wlen,sim_data,SIM_MAXCONV,sim_wlen_cnv,sim_data_cnv)) < 0)
    {
      err = 1;
      break;
    }
    if(Sampling(n2,sim_wlen_cnv,sim_data_cnv,dsr_n_wlen[rng],dsr_wlen[rng],dsr_dsim[rng][n],sim_yuni) < 0)
    {
      err = 1;
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

int ReadAUR_S(int rng,const double *par,double y[AUR_MAXDATA][AUR_MAXWLEN])
{
  // par[0] = #data in modtran.plt
  int i,m,n;
  int n1,n2;
  int err;
  int flag;
  double r,s;
  double v[2];
  char line[MAXLINE];
  char str1[MAXLINE];
  char str2[MAXLINE];
  char fnam[] = "ReadAUR_S";
  char *p;
  FILE *fp;

  if((fp=fopen("modtran.plt","r")) == NULL)
  {
    fprintf(stderr,"%s: error, cannot open %s\n",fnam,"modtran.plt");
    return -1;
  }
  n1 = (int)(par[0]+0.5);
  for(m=0; m<aur_n_wlen[rng]; m++)
  {
    for(n=0; n<aur_n_data; n++)
    {
      i = 0;
      s = 0.0;
      err = 0;
      flag = 0;
      while(fgets(line,MAXLINE,fp) != NULL)
      {
        if(sscanf(line,"%s%s",str1,str2) != 2)
        {
          flag = 1;
          break;
        }
        errno = 0;
        v[0] = strtod(str1,&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        errno = 0;
        v[1] = strtod(str2,&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        if(fabs(v[0]-aur_wlen[rng][m]) < SIM_DMAX)
        {
          s += v[1];
          i += 1;
        }
      }
      n2 = i;
      if(flag == 0)
      {
        fprintf(stderr,"%s: error, no blank line for (%d,%d)\n",fnam,m,n);
        err = 1;
      }
      if(n2 != n1)
      {
        fprintf(stderr,"%s: error, n1=%d, n2=%d\n",fnam,n1,n2);
        err = 1;
      }
      if(ResFactor(aur_ms720,aur_nw[rng][m],aur_sbmp[rng][m],aur_wlen[rng][m],aur_amas[n],&r) < 0)
      {
        err = 1;
      }
      if(err) break;
      y[n][m] = s*r*sim_yuni/n1;
    }
    if(err) break;
  }
  fclose(fp);
  if(err)
  {
    return -1;
  }

  return 0;
}

int ReadAUR_B(int rng,const double *par,double y[AUR_MAXDATA][AUR_MAXWLEN])
{
  // par[0] = wavelength margin
  int i,m,n;
  int m1,m2,m3;
  int n1,n2;
  int err;
  int flag;
  double v[2];
  char line[MAXLINE];
  char str1[MAXLINE];
  char str2[MAXLINE];
  char fnam[] = "ReadAUR_B";
  char *p;
  FILE *fp;

  if((fp=fopen("modtran.plt","r")) == NULL)
  {
    fprintf(stderr,"%s: error, cannot open %s\n",fnam,"modtran.plt");
    return -1;
  }
  for(n=0; n<aur_n_data; n++)
  {
    if(n%2 == 0)
    {
      m1 = 0;
      m2 = aur_n_wlen[rng];
      m3 = 1;
    }
    else
    {
      m1 = aur_n_wlen[rng]-1;
      m2 = -1;
      m3 = -1;
    }
    for(m=m1; m!=m2; m+=m3)
    {
      i = 0;
      err = 0;
      flag = 0;
      while(fgets(line,MAXLINE,fp) != NULL)
      {
        if(sscanf(line,"%s%s",str1,str2) != 2)
        {
          flag = 1;
          break;
        }
        errno = 0;
        v[0] = strtod(str1,&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        errno = 0;
        v[1] = strtod(str2,&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        if(i >= SIM_MAXDATA)
        {
          fprintf(stderr,"%s: error, #data exceed the limit >>> %d\n",fnam,i);
          err = 1;
          break;
        }
        sim_wlen[i] = v[0];
        sim_data[i] = v[1];
        i++;
      }
      n1 = i;
      if(flag == 0)
      {
        fprintf(stderr,"%s: error, no blank line for (%d,%d)\n",fnam,m,n);
        err = 1;
      }
      if(n1 < 1)
      {
        fprintf(stderr,"%s: no data\n",fnam);
        err = 1;
      }
      if(err) break;
      if((n2=Convolute(aur_wlen[rng][m]-par[0],aur_wlen[rng][m]+par[0],1.0,sim_xsgm,sim_wsgm,
                       n1,sim_wlen,sim_data,SIM_MAXCONV,sim_wlen_cnv,sim_data_cnv)) < 0)
      {
        err = 1;
        break;
      }
      if(Sampling(n2,sim_wlen_cnv,sim_data_cnv,1,&aur_wlen[rng][m],&y[n][m],sim_yuni) < 0)
      {
        err = 1;
        break;
      }
    }
    if(err) break;
  }
  fclose(fp);
  if(err)
  {
    return -1;
  }

  return 0;
}

int ReadSSR_S(int rng,const double *par,double y[SSR_MAXDATA][SSR_MAXWLEN])
{
  // par[0] = #data in modtran.plt
  int i,m,n;
  int n1,n2;
  int err;
  int flag;
  double r,s;
  double v[2];
  char line[MAXLINE];
  char str1[MAXLINE];
  char str2[MAXLINE];
  char fnam[] = "ReadSSR_S";
  char *p;
  FILE *fp;

  if((fp=fopen("modtran.plt","r")) == NULL)
  {
    fprintf(stderr,"%s: error, cannot open %s\n",fnam,"modtran.plt");
    return -1;
  }
  n1 = (int)(par[0]+0.5);
  for(m=0; m<ssr_n_wlen[rng]; m++)
  {
    for(n=0; n<ssr_n_data; n++)
    {
      i = 0;
      s = 0.0;
      err = 0;
      flag = 0;
      while(fgets(line,MAXLINE,fp) != NULL)
      {
        if(sscanf(line,"%s%s",str1,str2) != 2)
        {
          flag = 1;
          break;
        }
        errno = 0;
        v[0] = strtod(str1,&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        errno = 0;
        v[1] = strtod(str2,&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        if(fabs(v[0]-ssr_wlen[rng][m]) < SIM_DMAX)
        {
          s += v[1];
          i += 1;
        }
      }
      n2 = i;
      if(flag == 0)
      {
        fprintf(stderr,"%s: error, no blank line for (%d,%d)\n",fnam,m,n);
        err = 1;
      }
      if(n2 != n1)
      {
        fprintf(stderr,"%s: error, n1=%d, n2=%d\n",fnam,n1,n2);
        err = 1;
      }
      if(ResFactor(ssr_ms720,ssr_nw[rng][m],ssr_sbmp[rng][m],ssr_wlen[rng][m],ssr_amas[n],&r) < 0)
      {
        err = 1;
      }
      if(err) break;
      y[n][m] = s*r*sim_yuni/n1;
    }
    if(err) break;
  }
  fclose(fp);
  if(err)
  {
    return -1;
  }

  return 0;
}

int ReadSSR_B(int rng,const double *par,double y[SSR_MAXDATA][SSR_MAXWLEN])
{
  // par[0] = wavelength margin
  int i,m,n;
  int m1,m2,m3;
  int n1,n2;
  int err;
  int flag;
  double v[2];
  char line[MAXLINE];
  char str1[MAXLINE];
  char str2[MAXLINE];
  char fnam[] = "ReadSSR_B";
  char *p;
  FILE *fp;

  if((fp=fopen("modtran.plt","r")) == NULL)
  {
    fprintf(stderr,"%s: error, cannot open %s\n",fnam,"modtran.plt");
    return -1;
  }
  for(n=0; n<ssr_n_data; n++)
  {
    if(n%2 == 0)
    {
      m1 = 0;
      m2 = ssr_n_wlen[rng];
      m3 = 1;
    }
    else
    {
      m1 = ssr_n_wlen[rng]-1;
      m2 = -1;
      m3 = -1;
    }
    for(m=m1; m!=m2; m+=m3)
    {
      i = 0;
      err = 0;
      flag = 0;
      while(fgets(line,MAXLINE,fp) != NULL)
      {
        if(sscanf(line,"%s%s",str1,str2) != 2)
        {
          flag = 1;
          break;
        }
        errno = 0;
        v[0] = strtod(str1,&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        errno = 0;
        v[1] = strtod(str2,&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        if(i >= SIM_MAXDATA)
        {
          fprintf(stderr,"%s: error, #data exceed the limit >>> %d\n",fnam,i);
          err = 1;
          break;
        }
        sim_wlen[i] = v[0];
        sim_data[i] = v[1];
        i++;
      }
      n1 = i;
      if(flag == 0)
      {
        fprintf(stderr,"%s: error, no blank line for (%d,%d)\n",fnam,m,n);
        err = 1;
      }
      if(n1 < 1)
      {
        fprintf(stderr,"%s: no data\n",fnam);
        err = 1;
      }
      if(err) break;
      if((n2=Convolute(ssr_wlen[rng][m]-par[0],ssr_wlen[rng][m]+par[0],1.0,sim_xsgm,sim_wsgm,
                       n1,sim_wlen,sim_data,SIM_MAXCONV,sim_wlen_cnv,sim_data_cnv)) < 0)
      {
        err = 1;
        break;
      }
      if(Sampling(n2,sim_wlen_cnv,sim_data_cnv,1,&ssr_wlen[rng][m],&y[n][m],sim_yuni) < 0)
      {
        err = 1;
        break;
      }
    }
    if(err) break;
  }
  fclose(fp);
  if(err)
  {
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
  double norm;
  double offset = 0.5000000000001;
  char fnam[] = "Convolute";

  if(isnan(xmin)) xmin = xinp[0];
  if(isnan(xmax)) xmax = xinp[ninp-1];
  nstp = (int)((xmax-xmin)/xstp+0.5);
  if(nstp < 1)
  {
    fprintf(stderr,"%s: error, invalid #steps (xmin=%13.6e xmax=%13.6e xstp=%13.6e nstp=%d\n",
                    fnam,xmin,xmax,xstp,nstp);
    return -1;
  } else
  if(nstp >= nout)
  {
    fprintf(stderr,"%s: warning, nstp=%d -> changed to %d\n",fnam,nstp,nout-1);
    nstp = nout-1;
  }
  nwid = (int)(xsgm*wsgm/xstp+0.5);
  if(nwid < 1)
  {
    fprintf(stderr,"%s: error, nwid=%d, xsgm=%13.6e wsgm=%13.6e xstp=%13.6e\n",fnam,nwid,xsgm,wsgm,xstp);
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
  int i,j;
  int i1,i2;
  int err;
  double x1,x2;
  double xstp;
  double epsilon = 1.0e-12;
  char fnam[] = "Sampling";

  // check input
  if(ninp < 2)
  {
    fprintf(stderr,"%s: error, ninp=%d\n",fnam,ninp);
    return -1;
  }
  // check whether data interval is constant or not
  xstp = xinp[1]-xinp[0];
  err = 0;
  for(i=1; i<ninp; i++)
  {
    if(fabs(xinp[i]-(xinp[0]+xstp*i)) > epsilon)
    {
      err = 1;
      break;
    }
  }
  // interpolate
  if(err)
  {
    if(cnt_vb > 2) fprintf(stderr,"%s: nonconstant data interval.\n",fnam);
    x1 = xinp[0];
    x2 = xinp[ninp-1];
    for(i=0; i<nout; i++)
    {
      if(xout[i] < x1)
      {
        i1 = 0;
        i2 = 1;
      } else
      if(xout[i] > x2)
      {
        i2 = ninp-1;
        i1 = i2-1;
      }
      else
      {
        i1 = -1;
        i2 = -1;
        for(j=1; j<ninp; j++)
        {
          if(xout[i]>=xinp[j-1] && xout[i]<=xinp[j])
          {
            i1 = j-1;
            i2 = j;
            break;
          }
        }
        if(i1<0 || i2<0)
        {
          fprintf(stderr,"%s: error, faild in finding region for %13.6e\n",fnam,xout[i]);
          return -1;
        }
      }
      yout[i] = (yinp[i1]+(yinp[i2]-yinp[i1])*(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1]))*yuni;
    }
  }
  else
  {
    if(cnt_vb > 2) fprintf(stderr,"%s: constant data interval.\n",fnam);
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
  }

  return 0;
}

int SamplingE(int ninp,const double *xinp,const double *yinp,
              int nout,const double *xout,double *yout,double yuni)
{
  int i,j;
  int i1,i2;
  int err;
  double x1,x2;
  double xstp;
  double epsilon = 1.0e-12;
  double delta = 1.0e-100;
  char fnam[] = "SamplingE";

  // check input
  if(ninp < 2)
  {
    fprintf(stderr,"%s: error, ninp=%d\n",fnam,ninp);
    return -1;
  }
  // check whether data interval is constant or not
  xstp = xinp[1]-xinp[0];
  err = 0;
  for(i=1; i<ninp; i++)
  {
    if(fabs(xinp[i]-(xinp[0]+xstp*i)) > epsilon)
    {
      err = 1;
      break;
    }
  }
  // interpolate
  if(err)
  {
    if(cnt_vb > 2) fprintf(stderr,"%s: nonconstant data interval.\n",fnam);
    x1 = xinp[0];
    x2 = xinp[ninp-1];
    for(i=0; i<nout; i++)
    {
      if(xout[i] < x1)
      {
        i1 = 0;
        i2 = 1;
      } else
      if(xout[i] > x2)
      {
        i2 = ninp-1;
        i1 = i2-1;
      }
      else
      {
        i1 = -1;
        i2 = -1;
        for(j=1; j<ninp; j++)
        {
          if(xout[i]>=xinp[j-1] && xout[i]<=xinp[j])
          {
            i1 = j-1;
            i2 = j;
            break;
          }
        }
        if(i1<0 || i2<0)
        {
          fprintf(stderr,"%s: error, faild in finding region for %13.6e\n",fnam,xout[i]);
          return -1;
        }
      }
      if(yinp[i1]<delta || yinp[i2]<delta)
      {
        yout[i] = (yinp[i1]+(yinp[i2]-yinp[i1])*(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1]))*yuni;
      }
      else
      {
        yout[i] = (yinp[i1]*pow(yinp[i2]/yinp[i1],(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1])))*yuni;
      }
    }
  }
  else
  {
    if(cnt_vb > 2) fprintf(stderr,"%s: constant data interval.\n",fnam);
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
      if(yinp[i1]<delta || yinp[i2]<delta)
      {
        yout[i] = (yinp[i1]+(yinp[i2]-yinp[i1])*(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1]))*yuni;
      }
      else
      {
        yout[i] = (yinp[i1]*pow(yinp[i2]/yinp[i1],(xout[i]-xinp[i1])/(xinp[i2]-xinp[i1])))*yuni;
      }
    }
  }

  return 0;
}

int CommonPrintout(FILE *fp,int flag_dsr,int flag_aur,int flag_ssr)
{
  int i,j;
  int rng;
  int rid;

  fflush(fp);

  // Write DSR data
  if(flag_dsr != 0)
  {
    for(j=0; j<dsr_n_data; j++)
    {
      for(rng=0; rng<DSR_MAXRANG; rng++)
      {
        rid = rng;
        for(i=0; i<dsr_n_wlen[rng]; i++)
        {
          fprintf(fp,"%2d %4d %9.4f %9.4f %9.4f %22.14e %22.14e %22.14e\n",rid,j,dsr_th_sun[j],dsr_ph_sun[j],
                      dsr_wlen[rng][i],dsr_data[rng][j][i],dsr_dsim[rng][j][i],dsr_derr[rng][j][i]);
        }
      }
    }
    if(cnt_sflg==1 && flag_dsr>1)
    {
      for(j=0; j<dsr_n_data; j++)
      {
        for(rng=0; rng<DSR_MAXRANG; rng++)
        {
          rid = rng+DSR_MAXRANG;
          for(i=0; i<dsr_n_wlen[rng]; i++)
          {
            fprintf(fp,"%2d %4d %9.4f %9.4f %9.4f %22.14e\n",rid,j,dsr_th_sun[j],dsr_ph_sun[j],dsr_wlen[rng][i],dsr_dsgl[rng][j][i]);
          }
        }
      }
    }
  }
  // Write AUR data
  if(flag_aur != 0)
  {
    for(j=0; j<aur_n_data; j++)
    {
      for(rng=0; rng<AUR_MAXRANG; rng++)
      {
        rid = rng+DSR_MAXRANG*2;
        for(i=0; i<aur_n_wlen[rng]; i++)
        {
          fprintf(fp,"%2d %4d %9.4f %9.4f %9.4f %22.14e %22.14e %22.14e\n",rid,j,aur_th_sun[j],aur_ph_sun[j],
                      aur_wlen[rng][i],aur_data[rng][j][i],aur_dsim[rng][j][i],aur_derr[rng][j][i]);
        }
      }
    }
    if(cnt_sflg == 1)
    {
      for(j=0; j<aur_n_data; j++)
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          rid = rng+DSR_MAXRANG*2+AUR_MAXRANG;
          for(i=0; i<aur_n_wlen[rng]; i++)
          {
            fprintf(fp,"%2d %4d %9.4f %9.4f %9.4f %22.14e\n",rid,j,aur_th_sun[j],aur_ph_sun[j],aur_wlen[rng][i],aur_dsgl[rng][j][i]);
          }
        }
      }
    }
  }
  // Write SSR data
  if(flag_ssr != 0)
  {
    for(j=0; j<ssr_n_data; j++)
    {
      for(rng=0; rng<SSR_MAXRANG; rng++)
      {
        rid = rng+DSR_MAXRANG*2+AUR_MAXRANG*2;
        for(i=0; i<ssr_n_wlen[rng]; i++)
        {
          fprintf(fp,"%2d %4d %9.4f %9.4f %9.4f %22.14e %22.14e %22.14e\n",rid,j,ssr_th_los[j],ssr_ph_los[j],
                      ssr_wlen[rng][i],ssr_data[rng][j][i],ssr_dsim[rng][j][i],ssr_derr[rng][j][i]);
        }
      }
    }
    if(cnt_sflg == 1)
    {
      for(j=0; j<ssr_n_data; j++)
      {
        for(rng=0; rng<SSR_MAXRANG; rng++)
        {
          rid = rng+DSR_MAXRANG*2+AUR_MAXRANG*2+SSR_MAXRANG;
          for(i=0; i<ssr_n_wlen[rng]; i++)
          {
            fprintf(fp,"%2d %4d %9.4f %9.4f %9.4f %22.14e\n",rid,j,ssr_th_los[j],ssr_ph_los[j],ssr_wlen[rng][i],ssr_dsgl[rng][j][i]);
          }
        }
      }
    }
  }

  fflush(fp);

  return 0;
}

int CommonInit(void)
{
  int n;

  // initialization
  for(n=0; n<DSR_MAXRANG; n++)
  {
    dsr_n_wlen[n] = NODATA;
  }
  for(n=0; n<AUR_MAXRANG; n++)
  {
    aur_n_wlen[n] = NODATA;
    aur_n_sbmp[n] = NODATA;
  }
  for(n=0; n<SSR_MAXRANG; n++)
  {
    ssr_n_wlen[n] = NODATA;
    ssr_n_sbmp[n] = NODATA;
  }

  for(n=0; n<OPT_NPAR; n++)
  {
    sim_pval[n] = NAN;
  }
  for(n=OPT_NPAR; n<SIM_NPAR; n++)
  {
    sim_pval[n] = 0.0;
  }

  for(n=0; n<SIM_NPAR; n++)
  {
    crc_init[n] = NAN;
  }
  opt_init[OPT_PAR_VTAU] = OPT_VTAU;
  opt_init[OPT_PAR_VMIE] = OPT_VMIE;
  opt_init[OPT_PAR_WSCL] = OPT_WSCL;
  opt_init[OPT_PAR_OSCL] = OPT_OSCL;
  opt_init[OPT_PAR_XSCL] = OPT_XSCL;
  opt_init[OPT_PAR_YSCL] = OPT_YSCL;
  opt_init[OPT_PAR_ZSCL] = OPT_ZSCL;
  opt_init[OPT_PAR_NCMP] = OPT_NCMP;
  opt_init[OPT_PAR_LMD1] = OPT_LMOD;
  opt_init[OPT_PAR_LSG1] = OPT_LSGM;
  opt_init[OPT_PAR_MXR2] = OPT_MIXR;
  opt_init[OPT_PAR_LMD2] = OPT_LMOD;
  opt_init[OPT_PAR_LSG2] = OPT_LSGM;
  opt_init[OPT_PAR_MXR3] = OPT_MIXR;
  opt_init[OPT_PAR_LMD3] = OPT_LMOD;
  opt_init[OPT_PAR_LSG3] = OPT_LSGM;

  for(n=0; n<MIE_MAXCOMP; n++)
  {
    cnt_tmx[n] = 0;
    tmx_shap[n] = TMX_SHAP;
    tmx_reps[n] = TMX_REPS;
  }

  return 0;
}

int PreConfig(void)
{
  switch(sim_iatm)
  {
    case 0:
    case 1:
    case 2:
    case 4:
    case 6:
    case 7:
      sim_isea = 1;
      break;
    case 3:
    case 5:
      sim_isea = 2;
      break;
    default:
      sim_isea = 1;
      break;
  }
  if(sim_mode == SIM_MODTRAN_VTAU)
  {
    if(PrfInitMolecule(sim_iatm) < 0) return -1;
    if(PrfInitAerosol(sim_ihaz,sim_isea,sim_ivul,sim_vmie) < 0) return -1;
  }

  return 0;
}

int ReadConfig(void)
{
  int i,n,nc;
  int num;
  int cx,cy;
  int imin,imax;
  int idx;
  int err;
  double uni;
  double *x,*y;
  char line[MAXLINE];
  char temp[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char fnam[] = "ReadConfig";
  char *p;
  FILE *fp;

  if(strcmp(cnt_conf,"")==0 || strcmp(cnt_conf,NONAME)==0)
  {
    return 0;
  }
  // Read parameters which may have affects on others
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
    for(n=nc=0,p=temp; n<MAXENTR; n++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
    }
    if(n < 1) continue;
    if(strcasecmp(str[0],"mie_angl") == 0)
    {
      num = 0;
      uni = 1.0;
      mie_angl_min = MIE_ANGL_MIN;
      mie_angl_max = MIE_ANGL_MAX;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 3)
      {
        if(strncasecmp(str[3],"deg",3) == 0)
        {
          uni = 1.0;
        } else
        if(strncasecmp(str[3],"rad",3) == 0)
        {
          uni = R_TO_D;
        }
        else
        {
          fprintf(stderr,"%s: parameter error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 4)
      {
        errno = 0;
        mie_angl_min = strtod(str[4],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 5)
      {
        errno = 0;
        mie_angl_max = strtod(str[5],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((mie_n_angl=ReadXA(str[1],MIE_MAXDATA,num,uni,AnglSelect,mie_angl)) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %s %13.4f %13.4f\n",str[0],str[1],
                                               num,fabs(uni-1.0)>DELTA?"rad":"deg",
                                               mie_angl_min,mie_angl_max);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"mie_wlen") == 0)
    {
      num = 0;
      uni = 1.0;
      mie_wlen_min = MIE_WLEN_MIN;
      mie_wlen_max = MIE_WLEN_MAX;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 4)
      {
        errno = 0;
        mie_wlen_min = strtod(str[4],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 5)
      {
        errno = 0;
        mie_wlen_max = strtod(str[5],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((mie_n_wlen=ReadXA(str[1],MIE_MAXDATA,num,uni,WaveSelect,mie_wlen)) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e %13.4f %13.4f\n",
                                               str[0],str[1],num,uni,mie_wlen_min,mie_wlen_max);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"prf_falt") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((prf_n_alt=ReadXA(str[1],PRF_MAXNALT,num,uni,NULL,prf_alt)) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    }
  }
  fclose(fp);
  if(err)
  {
    return -1;
  }
  if(PreConfig() < 0) return -1;
  // Read others
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
    for(n=nc=0,p=temp; n<MAXENTR; n++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
    }
    if(n < 1) continue;
    if(strcasecmp(str[0],"dsr_wlen_0") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((dsr_n_wlen[0]=ReadXA(str[1],DSR_MAXWLEN,num,uni,NULL,dsr_wlen[0])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"dsr_wlen_1") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((dsr_n_wlen[1]=ReadXA(str[1],DSR_MAXWLEN,num,uni,NULL,dsr_wlen[1])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"dsr_wlen_2") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((dsr_n_wlen[2]=ReadXA(str[1],DSR_MAXWLEN,num,uni,NULL,dsr_wlen[2])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"dsr_wlen_3") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((dsr_n_wlen[3]=ReadXA(str[1],DSR_MAXWLEN,num,uni,NULL,dsr_wlen[3])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"dsr_fomg") == 0)
    {
      num = 1;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          x = NULL;
          y = NULL;
          if((nc=ReadXYP(str[1],DSR_MAXTDIV+1,0,num,1.0,1.0,NULL,&x,&y)) < 0)
          {
            err = 1;
            break;
          }
          if(nc < 2)
          {
            fprintf(stderr,"%s: error, nc=%d >>> %s\n",fnam,nc,line);
            err = 1;
          }
          else
          {
            uni = y[1]-y[0];
            for(i=2; i<nc; i++)
            {
              if(fabs(y[i]-y[i-1]-uni) > 0.1)
              {
                fprintf(stderr,"%s: error, nonconstant interval >>> %s\n",fnam,line);
                err = 1;
              }
            }
          }
          if(err == 0)
          {
            dsr_n_tdiv = nc-1;
            dsr_tmin = x[0];
            dsr_tmax = x[dsr_n_tdiv];
            for(i=0; i<dsr_n_tdiv; i++)
            {
              dsr_romg[i] = y[i];
            }
          }
          free(x);
          free(y);
          if(err)
          {
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %2d\n",str[0],str[1],num);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"aur_wlen_0") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((aur_n_wlen[0]=ReadXA(str[1],AUR_MAXWLEN,num,uni,NULL,aur_wlen[0])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"aur_wlen_1") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((aur_n_wlen[1]=ReadXA(str[1],AUR_MAXWLEN,num,uni,NULL,aur_wlen[1])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"aur_sbmp_0") == 0)
    {
      num = 0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((aur_n_sbmp[0]=ReadNA(str[1],AUR_MAXWLEN,num,aur_sbmp[0])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d\n",str[0],str[1],num);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"aur_sbmp_1") == 0)
    {
      num = 0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((aur_n_sbmp[1]=ReadNA(str[1],AUR_MAXWLEN,num,aur_sbmp[1])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d\n",str[0],str[1],num);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"aur_fomg") == 0)
    {
      num = 1;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          x = NULL;
          y = NULL;
          if((nc=ReadXYP(str[1],DSR_MAXTDIV+1,0,num,1.0,1.0,NULL,&x,&y)) < 0)
          {
            err = 1;
            break;
          }
          if(nc < 2)
          {
            fprintf(stderr,"%s: error, nc=%d >>> %s\n",fnam,nc,line);
            err = 1;
          }
          else
          {
            uni = y[1]-y[0];
            for(i=2; i<nc; i++)
            {
              if(fabs(y[i]-y[i-1]-uni) > 0.1)
              {
                fprintf(stderr,"%s: error, nonconstant interval >>> %s\n",fnam,line);
                err = 1;
              }
            }
          }
          if(err == 0)
          {
            aur_n_tdiv = nc-1;
            aur_tmin = x[0];
            aur_tmax = x[aur_n_tdiv];
            for(i=0; i<aur_n_tdiv; i++)
            {
              aur_romg[i] = y[i];
            }
          }
          free(x);
          free(y);
          if(err)
          {
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %2d\n",str[0],str[1],num);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_wlen_0") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((ssr_n_wlen[0]=ReadXA(str[1],SSR_MAXWLEN,num,uni,NULL,ssr_wlen[0])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_wlen_1") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((ssr_n_wlen[1]=ReadXA(str[1],SSR_MAXWLEN,num,uni,NULL,ssr_wlen[1])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_wlen_2") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((ssr_n_wlen[2]=ReadXA(str[1],SSR_MAXWLEN,num,uni,NULL,ssr_wlen[2])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_wlen_3") == 0)
    {
      num = 0;
      uni = 1.0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((ssr_n_wlen[3]=ReadXA(str[1],SSR_MAXWLEN,num,uni,NULL,ssr_wlen[3])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %13.6e\n",str[0],str[1],num,uni);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_sbmp_0") == 0)
    {
      num = 0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((ssr_n_sbmp[0]=ReadNA(str[1],SSR_MAXWLEN,num,ssr_sbmp[0])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d\n",str[0],str[1],num);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_sbmp_1") == 0)
    {
      num = 0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((ssr_n_sbmp[1]=ReadNA(str[1],SSR_MAXWLEN,num,ssr_sbmp[1])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d\n",str[0],str[1],num);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_sbmp_2") == 0)
    {
      num = 0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((ssr_n_sbmp[2]=ReadNA(str[1],SSR_MAXWLEN,num,ssr_sbmp[2])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d\n",str[0],str[1],num);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_sbmp_3") == 0)
    {
      num = 0;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          if((ssr_n_sbmp[3]=ReadNA(str[1],SSR_MAXWLEN,num,ssr_sbmp[3])) < 1)
          {
            err = 1;
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d\n",str[0],str[1],num);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_fomg") == 0)
    {
      num = 1;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          x = NULL;
          y = NULL;
          if((nc=ReadXYP(str[1],DSR_MAXTDIV+1,0,num,1.0,1.0,NULL,&x,&y)) < 0)
          {
            err = 1;
            break;
          }
          if(nc < 2)
          {
            fprintf(stderr,"%s: error, nc=%d >>> %s\n",fnam,nc,line);
            err = 1;
          }
          else
          {
            uni = y[1]-y[0];
            for(i=2; i<nc; i++)
            {
              if(fabs(y[i]-y[i-1]-uni) > 0.1)
              {
                fprintf(stderr,"%s: error, nonconstant interval >>> %s\n",fnam,line);
                err = 1;
              }
            }
          }
          if(err == 0)
          {
            ssr_n_tdiv = nc-1;
            ssr_tmin = x[0];
            ssr_tmax = x[ssr_n_tdiv];
            for(i=0; i<ssr_n_tdiv; i++)
            {
              ssr_romg[i] = y[i];
            }
          }
          free(x);
          free(y);
          if(err)
          {
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %2d\n",str[0],str[1],num);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_th_los_min") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        ssr_th_los_min = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.4e\n",str[0],ssr_th_los_min);
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
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.4e\n",str[0],ssr_th_los_max);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_sep_min") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        ssr_sep_min = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.4e\n",str[0],ssr_sep_min);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"ssr_sep_max") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        ssr_sep_max = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.4e\n",str[0],ssr_sep_max);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"crc_ffac") == 0)
    {
      if(n > 1)
      {
        strncpy(crc_ffac,str[1],MAXLINE);
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s\n",str[0],crc_ffac);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"crc_fsim") == 0)
    {
      if(n > 1)
      {
        strncpy(crc_fsim,str[1],MAXLINE);
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s\n",str[0],crc_fsim);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"crc_init") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=SIM_NPAR)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        crc_init[idx] = strtod(str[2],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27.14e\n",str[0],idx,crc_init[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"crc_fcmp") == 0)
    {
      cx = MIE_REAL_NUM;
      cy = MIE_IMAG_NUM;
      uni = 1.0;
      imin = MIE_IMIN;
      imax = MIE_IMAX;
      if(n > 2)
      {
        errno = 0;
        cx = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 3)
      {
        errno = 0;
        cy = strtol(str[3],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 4)
      {
        errno = 0;
        uni = strtod(str[4],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 5)
      {
        errno = 0;
        imin = strtol(str[5],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 6)
      {
        errno = 0;
        imax = strtol(str[6],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if((crc_n_comp=ReadComp(str[1],MIE_MAXCOMP,cx,cy,imin,imax,uni,
                              crc_wcom_com,crc_lmod_com,crc_lsgm_com,crc_refr_com,crc_refi_com)) < 0)
      {
        err = 1;
        break;
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %4d %13.6e %4d %4d\n",str[0],str[1],
                                               cx,cy,uni,imin,imax);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_fsim") == 0)
    {
      if(n > 1)
      {
        strncpy(tst_fsim,str[1],MAXLINE);
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s\n",str[0],tst_fsim);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tst_fcmp") == 0)
    {
      cx = MIE_REAL_NUM;
      cy = MIE_IMAG_NUM;
      uni = 1.0;
      imin = MIE_IMIN;
      imax = MIE_IMAX;
      if(n > 2)
      {
        errno = 0;
        cx = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 3)
      {
        errno = 0;
        cy = strtol(str[3],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 4)
      {
        errno = 0;
        uni = strtod(str[4],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 5)
      {
        errno = 0;
        imin = strtol(str[5],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 6)
      {
        errno = 0;
        imax = strtol(str[6],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if((tst_n_comp=ReadComp(str[1],MIE_MAXCOMP,cx,cy,imin,imax,uni,
                              tst_wcom_com,tst_lmod_com,tst_lsgm_com,tst_refr_com,tst_refi_com)) < 0)
      {
        err = 1;
        break;
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %4d %13.6e %4d %4d\n",str[0],str[1],
                                               cx,cy,uni,imin,imax);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_init") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=OPT_NPAR)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        opt_init[idx] = strtod(str[2],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27.14e\n",str[0],idx,opt_init[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"opt_fcmp") == 0)
    {
      cx = MIE_REAL_NUM;
      cy = MIE_IMAG_NUM;
      uni = 1.0;
      imin = MIE_IMIN;
      imax = MIE_IMAX;
      if(n > 2)
      {
        errno = 0;
        cx = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 3)
      {
        errno = 0;
        cy = strtol(str[3],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 4)
      {
        errno = 0;
        uni = strtod(str[4],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 5)
      {
        errno = 0;
        imin = strtol(str[5],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 6)
      {
        errno = 0;
        imax = strtol(str[6],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if((opt_n_comp=ReadComp(str[1],MIE_MAXCOMP,cx,cy,imin,imax,uni,
                              opt_wcom_com,opt_lmod_com,opt_lsgm_com,opt_refr_com,opt_refi_com)) < 0)
      {
        err = 1;
        break;
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %4d %13.6e %4d %4d\n",str[0],str[1],
                                               cx,cy,uni,imin,imax);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sim_path") == 0)
    {
      if(n > 1)
      {
        strncpy(sim_path,str[1],MAXLINE);
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s\n",str[0],sim_path);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sim_falb") == 0)
    {
      if(n > 1)
      {
        strncpy(sim_falb,str[1],MAXLINE);
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s\n",str[0],sim_falb);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sim_prf_vscl") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=PRF_NMOL)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        sim_prf_vscl[idx] = strtod(str[2],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27.14e\n",str[0],idx,sim_prf_vscl[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sim_prf_pscl") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=PRF_NMOL)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        sim_prf_pscl[idx] = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27d\n",str[0],idx,sim_prf_pscl[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"sim_prf_jcharx") == 0)
    {
      if(n > 1)
      {
        sim_prf_jcharx = str[1][0];
      }
      else
      {
        sim_prf_jcharx = ' ';
      }
      if(cnt_hp && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %30c\n",str[0],sim_prf_jcharx);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"prf_fpre") == 0)
    {
      num = 1;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 3)
      {
        prf_jchar[0] = str[3][0];
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          x = NULL;
          y = NULL;
          if((nc=ReadXYP(str[1],PRF_MAXNALT,0,num,1.0,1.0,NULL,&x,&y)) < 0)
          {
            err = 1;
            break;
          }
          if(Sampling(nc,x,y,prf_n_alt,prf_alt,prf_pres,1.0) < 0)
          {
            err = 1;
          }
          free(x);
          free(y);
          if(err)
          {
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %c\n",str[0],str[1],num,prf_jchar[0]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"prf_ftmp") == 0)
    {
      num = 1;
      if(n > 2)
      {
        errno = 0;
        num = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 3)
      {
        prf_jchar[1] = str[3][0];
      }
      if(n > 1)
      {
        if(strcmp(str[1],NONAME) != 0)
        {
          x = NULL;
          y = NULL;
          if((nc=ReadXYP(str[1],PRF_MAXNALT,0,num,1.0,1.0,NULL,&x,&y)) < 0)
          {
            err = 1;
            break;
          }
          if(Sampling(nc,x,y,prf_n_alt,prf_alt,prf_temp,1.0) < 0)
          {
            err = 1;
          }
          free(x);
          free(y);
          if(err)
          {
            break;
          }
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30s %4d %c\n",str[0],str[1],num,prf_jchar[1]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"prf_fmol") == 0)
    {
      idx = 0;
      num = 1;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=PRF_NMOL)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      if(n > 3)
      {
        errno = 0;
        num = strtol(str[3],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 4)
      {
        prf_jchar[idx+2] = str[4][0];
      }
      if(n > 2)
      {
        if(strcmp(str[2],NONAME) != 0)
        {
          x = NULL;
          y = NULL;
          if((nc=ReadXYP(str[2],PRF_MAXNALT,0,num,1.0,1.0,NULL,&x,&y)) < 0)
          {
            err = 1;
            break;
          }
          if(Sampling(nc,x,y,prf_n_alt,prf_alt,prf_wmol[idx],1.0) < 0)
          {
            err = 1;
          }
          free(x);
          free(y);
          if(err)
          {
            break;
          }
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %4d %25s %4d %c\n",str[0],idx,str[2],num,prf_jchar[idx+2]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"prf_fhaz") == 0)
    {
      idx = 0;
      num = 1;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=PRF_NHAZ)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      if(n > 3)
      {
        errno = 0;
        num = strtol(str[3],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(n > 2)
      {
        if(strcmp(str[2],NONAME) != 0)
        {
          x = NULL;
          y = NULL;
          if((nc=ReadXYP(str[2],PRF_MAXNALT,0,num,1.0,1.0,NULL,&x,&y)) < 0)
          {
            err = 1;
            break;
          }
          if(SamplingE(nc,x,y,prf_n_alt,prf_alt,prf_haze[idx],1.0) < 0)
          {
            err = 1;
          }
          free(x);
          free(y);
          if(err)
          {
            break;
          }
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %4d %25s %4d\n",str[0],idx,str[2],num);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"prf_pres_gnd") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        prf_pres_gnd = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.4e\n",str[0],prf_pres_gnd);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"prf_temp_gnd") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        prf_temp_gnd = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.4e\n",str[0],prf_temp_gnd);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"cnt_tmx") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        /* TO BE MODIFIED
        ntmp = strtol(str[1],&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>=0 && ntmp<MIE_MAXCOMP) idx = ntmp;
        else
        {
          fprintf(stderr,"%s: index out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        */
      }
      if(idx<0 || idx>=MIE_MAXCOMP)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        cnt_tmx[idx] = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        /* TO BE MODIFIED
        ntmp = strtol(str[2],&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>=0 && ntmp<MIE_MAXCOMP) cnt_tmx[idx] = ntmp;
        else
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
        */
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27d\n",str[0],idx,cnt_tmx[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tmx_nang") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tmx_nang = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0' || tmx_nang<2 || tmx_nang>=TMX_NAMAX)
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30d\n",str[0],tmx_nang);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tmx_ndgs") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tmx_ndgs = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0' || tmx_ndgs<0)
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30d\n",str[0],tmx_ndgs);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tmx_kmax") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tmx_kmax = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0' || tmx_kmax<-2)
        {
          fprintf(stderr,"%s: out of range >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30d\n",str[0],tmx_kmax);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tmx_delt") == 0)
    {
      if(n > 1)
      {
        errno = 0;
        tmx_delt = strtod(str[1],&p);
        if(errno==ERANGE || *p!='\0' || tmx_delt<=0.0)
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>1 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        snprintf(cnt_cmnt[cnt_n_cmnt],MAXLINE,"%-14s: %30.4e\n",str[0],tmx_delt);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tmx_shap") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=MIE_MAXCOMP)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        tmx_shap[idx] = strtol(str[2],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27d\n",str[0],idx,tmx_shap[idx]);
        cnt_n_cmnt++;
      }
    } else
    if(strcasecmp(str[0],"tmx_reps") == 0)
    {
      idx = 0;
      if(n > 1)
      {
        errno = 0;
        idx = strtol(str[1],&p,10);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(idx<0 || idx>=MIE_MAXCOMP)
      {
        fprintf(stderr,"%s: index out of range >>> %s\n",fnam,line);
        err = 1;
        break;
      }
      if(n > 2)
      {
        errno = 0;
        tmx_reps[idx] = strtod(str[2],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error >>> %s\n",fnam,line);
          err = 1;
          break;
        }
      }
      if(cnt_hp && n>2 && cnt_n_cmnt<CNT_MAXCMNT)
      {
        sprintf(cnt_cmnt[cnt_n_cmnt],"%-14s: %2d %27.14e\n",str[0],idx,tmx_reps[idx]);
        cnt_n_cmnt++;
      }
    }
    else
    {
      err = ReadOwnConfig(line);
      if(err > 0)
      {
        fprintf(stderr,"%s: unrecognized option >>> %s\n",fnam,line);
        break;
      } else
      if(err < 0)
      {
        err = 1;
        break;
      }
    }
  }
  fclose(fp);
  if(err)
  {
    return -1;
  }
  if(PostConfig() < 0) return -1;

  return 0;
}

int PostConfig(void)
{
  int i,j;
  int flag;

  // post configuration
  if(sim_mode == SIM_MODTRAN_VTAU)
  {
    if(PrfGround() < 0) return -1;
    if(PrfScaleMolecule() < 0) return -1;
    if(PrfCalcAOD() < 0) return -1;
    if(cnt_vb > 1)
    {
      flag = 0;
      fprintf(stderr,"Molecular profile:\n");
      fprintf(stderr,"%6s %9s %9s","Alt","Pres","Temp");
      if(!isnan(prf_pres_gnd)) flag = 1;
      if(!isnan(prf_temp_gnd)) flag = 1;
      for(j=0; j<PRF_NMOL; j++)
      {
        fprintf(stderr," %9s",prf_wmol_nam[j]);
        if(sim_prf_pscl[j]!=0 || (!isnan(sim_prf_vscl[j]))) flag = 1;
      }
      fprintf(stderr,"\n");
      for(i=0; i<prf_n_alt; i++)
      {
        fprintf(stderr,"%6.2f %9.2e %9.2e",prf_alt[i],prf_pres[i],prf_temp[i]);
        for(j=0; j<PRF_NMOL; j++)
        {
          fprintf(stderr," %9.2e",prf_wmol[j][i]);
        }
        fprintf(stderr,"\n");
      }
      if(flag)
      {
        fprintf(stderr,"After scaling:\n");
        fprintf(stderr,"%6s","Alt");
        if(!isnan(prf_pres_gnd))
        {
          fprintf(stderr," %9s","Pres");
        }
        if(!isnan(prf_temp_gnd))
        {
          fprintf(stderr," %9s","Temp");
        }
        for(j=0; j<PRF_NMOL; j++)
        {
          if(sim_prf_pscl[j]!=0 || (!isnan(sim_prf_vscl[j])))
          {
            fprintf(stderr," %9s",prf_wmol_nam[j]);
          }
        }
        fprintf(stderr,"\n");
        for(i=0; i<prf_n_alt; i++)
        {
          fprintf(stderr,"%6.2f",prf_alt[i]);
          if(!isnan(prf_pres_gnd))
          {
            fprintf(stderr," %9.2e",sim_prf_pres[i]);
          }
          if(!isnan(prf_temp_gnd))
          {
            fprintf(stderr," %9.2e",sim_prf_temp[i]);
          }
          for(j=0; j<PRF_NMOL; j++)
          {
            if(sim_prf_pscl[j]!=0 || (!isnan(sim_prf_vscl[j])))
            {
              fprintf(stderr," %9.2e",sim_prf_wmol[j][i]);
            }
          }
          fprintf(stderr,"\n");
        }
      }
      fprintf(stderr,"Aerosol profile:\n");
      for(i=0; i<prf_n_alt; i++)
      {
        fprintf(stderr,"%6.2f %13.10f %13.10f %13.10f %13.10f\n",prf_alt[i],prf_haze[0][i],prf_haze[1][i],prf_haze[2][i],prf_haze[3][i]);
      }
    }
  }

  return 0;
}

int ReadSimul(char *s,int flag_dsr,int flag_aur,int flag_ssr,int flag_dat,int flag_sim)
{
  int i,j;
  int rng;
  int rid;
  int n,nc;
  int err;
  int v[SIM_NINT];
  double w[SIM_NDBL];
  double epsilon;
  char line[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char fnam[] = "ReadSimul";
  char *p;
  FILE *fp;

  if(strcmp(s,"")==0 || strcmp(s,NONAME)==0)
  {
    return 0;
  } else
  if((fp=fopen(s,"r")) == NULL)
  {
    fprintf(stderr,"%s: error, cannot open %s\n",fnam,s);
    return -1;
  } else
  if(cnt_vb)
  {
    fprintf(stderr,"%s: reading %s\n",fnam,s);
  }
  epsilon = 360.0*inp_tdif/SEC_DAY;
  err = 0;
  // Read DSR data
  if(flag_dsr != 0)
  {
    for(j=0; j<dsr_n_data; j++)
    {
      for(rng=0; rng<DSR_MAXRANG; rng++)
      {
        rid = rng;
        for(i=0; i<dsr_n_wlen[rng]; i++)
        {
          if(fgets(line,MAXLINE,fp) == NULL)
          {
            fprintf(stderr,"%s: read error 1 (i=%d,j=%d,rng=%d)\n",fnam,i,j,rng);
            err = 1;
            break;
          }
          // read line
          for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
          {
            if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
          }
          if(n != SIM_NCOL)
          {
            fprintf(stderr,"%s: error, expected %d columns >>> %d\n",fnam,SIM_NCOL,n);
            err = 1;
            break;
          }
          // convert values
          for(n=0; n<SIM_NINT; n++)
          {
            errno = 0;
            v[n] = strtol(str[n],&p,10);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
              err = 1;
              break;
            }
          }
          if(err) break;
          for(n=SIM_NINT; n<SIM_NCOL; n++)
          {
            errno = 0;
            w[n-SIM_NINT] = strtod(str[n],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
              err = 1;
              break;
            }
          }
          if(err) break;
          if(v[0] != rid)
          {
            fprintf(stderr,"%s: error, v[0]=%d, rid=%d\n",fnam,v[0],rid);
            err = 1;
            break;
          } else
          if(v[1] != j)
          {
            fprintf(stderr,"%s: error, v[1]=%d, j=%d\n",fnam,v[1],j);
            err = 1;
            break;
          } else
          if(fabs(w[0]-dsr_th_sun[j]) > epsilon)
          {
            fprintf(stderr,"%s: error, w[0]=%13.6e, dsr_th_sun[%d]=%13.6e\n",fnam,w[0],j,dsr_th_sun[j]);
            err = 1;
            break;
          } else
          if(fabs(w[1]-dsr_ph_sun[j]) > epsilon)
          {
            fprintf(stderr,"%s: error, w[1]=%13.6e, dsr_ph_sun[%d]=%13.6e\n",fnam,w[1],j,dsr_ph_sun[j]);
            err = 1;
            break;
          } else
          if(fabs(w[2]-dsr_wlen[rng][i]) > DELTA)
          {
            fprintf(stderr,"%s: error, w[2]=%13.6e, dsr_wlen[%d][%d]=%13.6e\n",fnam,w[2],rng,i,dsr_wlen[rng][i]);
            err = 1;
            break;
          }
          else
          {
            if(flag_dat != 0) dsr_data[rng][j][i] = w[3];
            if(flag_sim != 0) dsr_dsim[rng][j][i] = w[4];
          }
        }
        if(err) break;
      }
      if(err) break;
    }
    if(err)
    {
      fclose(fp);
      return -1;
    }
    if(cnt_sflg==1 && flag_dsr>1)
    {
      for(j=0; j<dsr_n_data; j++)
      {
        for(rng=0; rng<DSR_MAXRANG; rng++)
        {
          rid = rng+DSR_MAXRANG;
          for(i=0; i<dsr_n_wlen[rng]; i++)
          {
            if(fgets(line,MAXLINE,fp) == NULL)
            {
              fprintf(stderr,"%s: read error 2 (i=%d,j=%d,rng=%d)\n",fnam,i,j,rng);
              err = 1;
              break;
            }
            // read line
            for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
            {
              if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
            }
            if(n != SIM_NSGL)
            {
              fprintf(stderr,"%s: error, expected %d columns >>> %d\n",fnam,SIM_NSGL,n);
              err = 1;
              break;
            }
            // convert values
            for(n=0; n<SIM_NINT; n++)
            {
              errno = 0;
              v[n] = strtol(str[n],&p,10);
              if(errno==ERANGE || *p!='\0')
              {
                fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
                err = 1;
                break;
              }
            }
            if(err) break;
            for(n=SIM_NINT; n<SIM_NSGL; n++)
            {
              errno = 0;
              w[n-SIM_NINT] = strtod(str[n],&p);
              if(errno==ERANGE || *p!='\0')
              {
                fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
                err = 1;
                break;
              }
            }
            if(err) break;
            if(v[0] != rid)
            {
              fprintf(stderr,"%s: error, v[0]=%d, rid=%d\n",fnam,v[0],rid);
              err = 1;
              break;
            } else
            if(v[1] != j)
            {
              fprintf(stderr,"%s: error, v[1]=%d, j=%d\n",fnam,v[1],j);
              err = 1;
              break;
            } else
            if(fabs(w[0]-dsr_th_sun[j]) > epsilon)
            {
              fprintf(stderr,"%s: error, w[0]=%13.6e, dsr_th_sun[%d]=%13.6e\n",fnam,w[0],j,dsr_th_sun[j]);
              err = 1;
              break;
            } else
            if(fabs(w[1]-dsr_ph_sun[j]) > epsilon)
            {
              fprintf(stderr,"%s: error, w[1]=%13.6e, dsr_ph_sun[%d]=%13.6e\n",fnam,w[1],j,dsr_ph_sun[j]);
              err = 1;
              break;
            } else
            if(fabs(w[2]-dsr_wlen[rng][i]) > DELTA)
            {
              fprintf(stderr,"%s: error, w[2]=%13.6e, dsr_wlen[%d][%d]=%13.6e\n",fnam,w[2],rng,i,dsr_wlen[rng][i]);
              err = 1;
              break;
            }
            else
            {
              if(flag_sim != 0) dsr_dsgl[rng][j][i] = w[3];
            }
          }
          if(err) break;
        }
        if(err) break;
      }
      if(err)
      {
        fclose(fp);
        return -1;
      }
    }
  }
  // Read AUR data
  if(flag_aur != 0)
  {
    for(j=0; j<aur_n_data; j++)
    {
      for(rng=0; rng<AUR_MAXRANG; rng++)
      {
        rid = rng+DSR_MAXRANG*2;
        for(i=0; i<aur_n_wlen[rng]; i++)
        {
          if(fgets(line,MAXLINE,fp) == NULL)
          {
            fprintf(stderr,"%s: read error 3 (i=%d,j=%d,rng=%d)\n",fnam,i,j,rng);
            err = 1;
            break;
          }
          // read line
          for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
          {
            if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
          }
          if(n != SIM_NCOL)
          {
            fprintf(stderr,"%s: error, expected %d columns >>> %d\n",fnam,SIM_NCOL,n);
            err = 1;
            break;
          }
          // convert values
          for(n=0; n<SIM_NINT; n++)
          {
            errno = 0;
            v[n] = strtol(str[n],&p,10);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
              err = 1;
              break;
            }
          }
          if(err) break;
          for(n=SIM_NINT; n<SIM_NCOL; n++)
          {
            errno = 0;
            w[n-SIM_NINT] = strtod(str[n],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
              err = 1;
              break;
            }
          }
          if(err) break;
          if(v[0] != rid)
          {
            fprintf(stderr,"%s: error, v[0]=%d, rid=%d\n",fnam,v[0],rid);
            err = 1;
            break;
          } else
          if(v[1] != j)
          {
            fprintf(stderr,"%s: error, v[1]=%d, j=%d\n",fnam,v[1],j);
            err = 1;
            break;
          } else
          if(fabs(w[0]-aur_th_sun[j]) > epsilon)
          {
            fprintf(stderr,"%s: error, w[0]=%13.6e, aur_th_sun[%d]=%13.6e\n",fnam,w[0],j,aur_th_sun[j]);
            err = 1;
            break;
          } else
          if(fabs(w[1]-aur_ph_sun[j]) > epsilon)
          {
            fprintf(stderr,"%s: error, w[1]=%13.6e, aur_ph_sun[%d]=%13.6e\n",fnam,w[1],j,aur_ph_sun[j]);
            err = 1;
            break;
          } else
          if(fabs(w[2]-aur_wlen[rng][i]) > DELTA)
          {
            fprintf(stderr,"%s: error, w[2]=%13.6e, aur_wlen[%d][%d]=%13.6e\n",fnam,w[2],rng,i,aur_wlen[rng][i]);
            err = 1;
            break;
          }
          else
          {
            if(flag_dat != 0) aur_data[rng][j][i] = w[3];
            if(flag_sim != 0) aur_dsim[rng][j][i] = w[4];
          }
        }
        if(err) break;
      }
      if(err) break;
    }
    if(err)
    {
      fclose(fp);
      return -1;
    }
    if(cnt_sflg == 1)
    {
      for(j=0; j<aur_n_data; j++)
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          rid = rng+DSR_MAXRANG*2+AUR_MAXRANG;
          for(i=0; i<aur_n_wlen[rng]; i++)
          {
            if(fgets(line,MAXLINE,fp) == NULL)
            {
              fprintf(stderr,"%s: read error 4 (i=%d,j=%d,rng=%d)\n",fnam,i,j,rng);
              err = 1;
              break;
            }
            // read line
            for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
            {
              if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
            }
            if(n != SIM_NSGL)
            {
              fprintf(stderr,"%s: error, expected %d columns >>> %d\n",fnam,SIM_NSGL,n);
              err = 1;
              break;
            }
            // convert values
            for(n=0; n<SIM_NINT; n++)
            {
              errno = 0;
              v[n] = strtol(str[n],&p,10);
              if(errno==ERANGE || *p!='\0')
              {
                fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
                err = 1;
                break;
              }
            }
            if(err) break;
            for(n=SIM_NINT; n<SIM_NSGL; n++)
            {
              errno = 0;
              w[n-SIM_NINT] = strtod(str[n],&p);
              if(errno==ERANGE || *p!='\0')
              {
                fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
                err = 1;
                break;
              }
            }
            if(err) break;
            if(v[0] != rid)
            {
              fprintf(stderr,"%s: error, v[0]=%d, rid=%d\n",fnam,v[0],rid);
              err = 1;
              break;
            } else
            if(v[1] != j)
            {
              fprintf(stderr,"%s: error, v[1]=%d, j=%d\n",fnam,v[1],j);
              err = 1;
              break;
            } else
            if(fabs(w[0]-aur_th_sun[j]) > epsilon)
            {
              fprintf(stderr,"%s: error, w[0]=%13.6e, aur_th_sun[%d]=%13.6e\n",fnam,w[0],j,aur_th_sun[j]);
              err = 1;
              break;
            } else
            if(fabs(w[1]-aur_ph_sun[j]) > epsilon)
            {
              fprintf(stderr,"%s: error, w[1]=%13.6e, aur_ph_sun[%d]=%13.6e\n",fnam,w[1],j,aur_ph_sun[j]);
              err = 1;
              break;
            } else
            if(fabs(w[2]-aur_wlen[rng][i]) > DELTA)
            {
              fprintf(stderr,"%s: error, w[2]=%13.6e, aur_wlen[%d][%d]=%13.6e\n",fnam,w[2],rng,i,aur_wlen[rng][i]);
              err = 1;
              break;
            }
            else
            {
              if(flag_sim != 0) aur_dsgl[rng][j][i] = w[3];
            }
          }
          if(err) break;
        }
        if(err) break;
      }
      if(err)
      {
        fclose(fp);
        return -1;
      }
    }
  }
  // Read SSR data
  if(flag_ssr != 0)
  {
    for(j=0; j<ssr_n_data; j++)
    {
      for(rng=0; rng<SSR_MAXRANG; rng++)
      {
        rid = rng+DSR_MAXRANG*2+AUR_MAXRANG*2;
        for(i=0; i<ssr_n_wlen[rng]; i++)
        {
          if(fgets(line,MAXLINE,fp) == NULL)
          {
            fprintf(stderr,"%s: read error 5 (i=%d,j=%d,rng=%d)\n",fnam,i,j,rng);
            err = 1;
            break;
          }
          // read line
          for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
          {
            if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
          }
          if(n != SIM_NCOL)
          {
            fprintf(stderr,"%s: error, expected %d columns >>> %d\n",fnam,SIM_NCOL,n);
            err = 1;
            break;
          }
          // convert values
          for(n=0; n<SIM_NINT; n++)
          {
            errno = 0;
            v[n] = strtol(str[n],&p,10);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
              err = 1;
              break;
            }
          }
          if(err) break;
          for(n=SIM_NINT; n<SIM_NCOL; n++)
          {
            errno = 0;
            w[n-SIM_NINT] = strtod(str[n],&p);
            if(errno==ERANGE || *p!='\0')
            {
              fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
              err = 1;
              break;
            }
          }
          if(err) break;
          if(v[0] != rid)
          {
            fprintf(stderr,"%s: error, v[0]=%d, rid=%d\n",fnam,v[0],rid);
            err = 1;
            break;
          } else
          if(v[1] != j)
          {
            fprintf(stderr,"%s: error, v[1]=%d, j=%d\n",fnam,v[1],j);
            err = 1;
            break;
          } else
          if(fabs(w[0]-ssr_th_los[j]) > DELTA)
          {
            fprintf(stderr,"%s: error, w[0]=%13.6e, ssr_th_los[%d]=%13.6e\n",fnam,w[0],j,ssr_th_los[j]);
            err = 1;
            break;
          } else
          if(fabs(w[1]-ssr_ph_los[j]) > DELTA)
          {
            fprintf(stderr,"%s: error, w[1]=%13.6e, ssr_ph_los[%d]=%13.6e\n",fnam,w[1],j,ssr_ph_los[j]);
            err = 1;
            break;
          } else
          if(fabs(w[2]-ssr_wlen[rng][i]) > DELTA)
          {
            fprintf(stderr,"%s: error, w[2]=%13.6e, ssr_wlen[%d][%d]=%13.6e\n",fnam,w[2],rng,i,ssr_wlen[rng][i]);
            err = 1;
            break;
          }
          else
          {
            if(flag_dat != 0) ssr_data[rng][j][i] = w[3];
            if(flag_sim != 0) ssr_dsim[rng][j][i] = w[4];
          }
        }
        if(err) break;
      }
      if(err) break;
    }
    if(err)
    {
      fclose(fp);
      return -1;
    }
    if(cnt_sflg == 1)
    {
      for(j=0; j<ssr_n_data; j++)
      {
        for(rng=0; rng<SSR_MAXRANG; rng++)
        {
          rid = rng+DSR_MAXRANG*2+AUR_MAXRANG*2+SSR_MAXRANG;
          for(i=0; i<ssr_n_wlen[rng]; i++)
          {
            if(fgets(line,MAXLINE,fp) == NULL)
            {
              fprintf(stderr,"%s: read error 6 (i=%d,j=%d,rng=%d)\n",fnam,i,j,rng);
              err = 1;
              break;
            }
            // read line
            for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
            {
              if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
            }
            if(n != SIM_NSGL)
            {
              fprintf(stderr,"%s: error, expected %d columns >>> %d\n",fnam,SIM_NSGL,n);
              err = 1;
              break;
            }
            // convert values
            for(n=0; n<SIM_NINT; n++)
            {
              errno = 0;
              v[n] = strtol(str[n],&p,10);
              if(errno==ERANGE || *p!='\0')
              {
                fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
                err = 1;
                break;
              }
            }
            if(err) break;
            for(n=SIM_NINT; n<SIM_NSGL; n++)
            {
              errno = 0;
              w[n-SIM_NINT] = strtod(str[n],&p);
              if(errno==ERANGE || *p!='\0')
              {
                fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
                err = 1;
                break;
              }
            }
            if(err) break;
            if(v[0] != rid)
            {
              fprintf(stderr,"%s: error, v[0]=%d, rid=%d\n",fnam,v[0],rid);
              err = 1;
              break;
            } else
            if(v[1] != j)
            {
              fprintf(stderr,"%s: error, v[1]=%d, j=%d\n",fnam,v[1],j);
              err = 1;
              break;
            } else
            if(fabs(w[0]-ssr_th_los[j]) > DELTA)
            {
              fprintf(stderr,"%s: error, w[0]=%13.6e, ssr_th_los[%d]=%13.6e\n",fnam,w[0],j,ssr_th_los[j]);
              err = 1;
              break;
            } else
            if(fabs(w[1]-ssr_ph_los[j]) > DELTA)
            {
              fprintf(stderr,"%s: error, w[1]=%13.6e, ssr_ph_los[%d]=%13.6e\n",fnam,w[1],j,ssr_ph_los[j]);
              err = 1;
              break;
            } else
            if(fabs(w[2]-ssr_wlen[rng][i]) > DELTA)
            {
              fprintf(stderr,"%s: error, w[2]=%13.6e, ssr_wlen[%d][%d]=%13.6e\n",fnam,w[2],rng,i,ssr_wlen[rng][i]);
              err = 1;
              break;
            }
            else
            {
              if(flag_sim != 0) ssr_dsgl[rng][j][i] = w[3];
            }
          }
          if(err) break;
        }
        if(err) break;
      }
      if(err)
      {
        fclose(fp);
        return -1;
      }
    }
  }
  fclose(fp);

  return 1;
}

int ReadInput(void)
{
  int i,j,k,n;
  int rng;
  int nc;
  int err;
  int v[INP_NINT];
  int err_num[ERR_N_WLEN];
  int dsr_repeat[DSR_MAXDATA];
  int aur_repeat[AUR_MAXDATA];
  int ssr_repeat[SSR_MAXDATA];
  int err_dsr_nw[INP_N_WLEN];
  int err_aur_nw[INP_N_WLEN];
  int err_ssr_nw[INP_N_WLEN];
  time_t dsr_time_tmp;
  time_t aur_time_tmp;
  time_t ssr_time_tmp;
  double w[INP_NDBL];
  double err_avg[ERR_N_WLEN];
  double ssr_th_los_tmp;
  double ssr_ph_los_tmp;
  double etmp;
  double epsilon = 1.0e-9;
  double err_tmp[INP_N_WLEN];
  double dsr_inp_avg[DSR_MAXDATA][INP_N_WLEN];
  double dsr_inp_sqr[DSR_MAXDATA][INP_N_WLEN];
  double aur_inp_avg[AUR_MAXDATA][INP_N_WLEN];
  double aur_inp_sqr[AUR_MAXDATA][INP_N_WLEN];
  double ssr_inp_avg[SSR_MAXDATA][INP_N_WLEN];
  double ssr_inp_sqr[SSR_MAXDATA][INP_N_WLEN];
  double dsr_aur_avg[INP_N_WLEN];
  double dsr_aur_sqr[INP_N_WLEN];
  char line[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char fnam[MAXLINE];
  char format[MAXLINE];
  char *p;
  FILE *fp;
  struct tm t;

  // read stdin
  // format:
  // rno yr mo dy hr mi sc site_id inst_id cloud flag fov vazi valt par exposure temperature sazi salt airmass distance
  putenv("TZ=");
  t.tm_isdst = 0;
  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,stdin) != NULL)
  {
    // read line
    for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
    }
    if(n != INP_NCOL)
    {
      fprintf(stderr,"Error, expected %d columns >>> %d\n",INP_NCOL,n);
      err = 1;
      break;
    }
    // convert values
    for(n=0; n<INP_NINT; n++)
    {
      errno = 0;
      v[n] = strtol(str[n],&p,10);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"Convert error >>> %s\n",str[n]);
        err = 1;
        break;
      }
    }
    if(err) break;
    for(n=INP_NINT; n<INP_NCOL; n++)
    {
      errno = 0;
      w[n-INP_NINT] = strtod(str[n],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"Convert error >>> %s\n",str[n]);
        err = 1;
        break;
      }
    }
    if(err) break;
    if(i >= INP_MAXDATA)
    {
      fprintf(stderr,"Warning, #data exceed the limit >>> %d\n",i);
      break;
    }
    inp_rno[i]    = v[0];
    t.tm_year     = v[1]-1900;
    t.tm_mon      = v[2]-1;
    t.tm_mday     = v[3];
    t.tm_hour     = v[4];
    t.tm_min      = v[5];
    t.tm_sec      = v[6];
    t.tm_isdst    = 0;
    inp_time[i]   = mktime(&t);
    inp_ms720[i]  = v[8];
    inp_cloud[i]  = v[9];
    inp_flag[i]   = v[10];
    inp_ph_los[i] = w[1];
    inp_th_los[i] = 90.0-w[2];
    inp_temp[i]   = w[5];
    inp_ph_sun[i] = w[6];
    inp_th_sun[i] = 90.0-w[7];
    inp_amas[i]   = w[8];
    inp_dist[i]   = w[9];
    if(inp_flag[i]>=100 && inp_flag[i]<200) // DSR
    {
      inp_type[i] = INP_TYPE_DSR;
      if(dsr_ms720 < 0)
      {
        dsr_ms720 = inp_ms720[i];
      } else
      if(inp_ms720[i] != dsr_ms720)
      {
        fprintf(stderr,"Error, conflicting MS720ID for DSR.\n");
        err = 1;
        break;
      }
    } else
    if(inp_flag[i]>=300 && inp_flag[i]<400 && cnt_enable_old_aur) // AUR (old method)
    {
      inp_type[i] = INP_TYPE_AUR;
      if(aur_ms720 < 0)
      {
        aur_ms720 = inp_ms720[i];
      } else
      if(inp_ms720[i] != aur_ms720)
      {
        fprintf(stderr,"Error, conflicting MS720ID for AUR.\n");
        err = 1;
        break;
      }
      if(aur_type < 0)
      {
        aur_type = AUR_TYPE_OLD;
      } else
      if(aur_type != AUR_TYPE_OLD)
      {
        fprintf(stderr,"Error, conflicting AUR type.\n");
        err = 1;
        break;
      }
    } else
    if(inp_flag[i]>=500 && inp_flag[i]<600) // AUR
    {
      inp_type[i] = INP_TYPE_AUR;
      if(aur_ms720 < 0)
      {
        aur_ms720 = inp_ms720[i];
      } else
      if(inp_ms720[i] != aur_ms720)
      {
        fprintf(stderr,"Error, conflicting MS720ID for AUR.\n");
        err = 1;
        break;
      }
      if(aur_type < 0)
      {
        aur_type = AUR_TYPE_NEW;
      } else
      if(aur_type != AUR_TYPE_NEW)
      {
        fprintf(stderr,"Error, conflicting aur type.\n");
        err = 1;
        break;
      }
    } else
    if(inp_flag[i]>=0 && inp_flag[i]<100)   // SSR
    {
      inp_type[i] = INP_TYPE_SSR;
      if(ssr_ms720 < 0)
      {
        ssr_ms720 = inp_ms720[i];
      } else
      if(inp_ms720[i] != ssr_ms720)
      {
        fprintf(stderr,"Error, conflicting MS720ID for SSR.\n");
        err = 1;
        break;
      }
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
    err = 1;
  }
  if(err)
  {
    return -1;
  }
  if(sim_iday < 0)
  {
    localtime_r(&inp_time[0],&t);
    sim_iday = t.tm_yday+1;
  }

  // read input data
  for(j=0; j<INP_N_WLEN; j++)
  {
    inp_dsr_wlen[j] = 0.0;
    inp_aur_wlen[j] = 0.0;
    inp_ssr_wlen[j] = 0.0;
  }
  snprintf(format,MAXLINE,"%%s/%%0%dd.dat",inp_ndig);
  for(i=0; i<inp_n_data; i++)
  {
    snprintf(fnam,MAXLINE,format,inp_dnam,inp_rno[i]);
    if(cnt_vb > 0)
    {
      fprintf(stderr,"Reading %s\n",fnam);
    }
    if((fp=fopen(fnam,"r")) == NULL)
    {
      fprintf(stderr,"Error, cannot open %s\n",fnam);
      return -1;
    }
    j = 0;
    err = 0;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      if(sscanf(line,"%s%s",str[0],str[1]) != 2)
      {
        fprintf(stderr,"Read error >>> %s\n",line);
        err = 1;
        break;
      }
      errno = 0;
      for(n=0; n<2; n++)
      {
        w[n] = strtod(str[n],&p);
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
      if(inp_type[i] == INP_TYPE_DSR) // DSR
      {
        if(inp_dsr_wlen[j] < EPSILON)
        {
          inp_dsr_wlen[j] = w[0];
        }
        else
        {
          if(fabs(inp_dsr_wlen[j]-w[0]) > DELTA)
          {
            fprintf(stderr,"Error, inp_dsr_wlen[%d]=%13.6e, w[0]=%13.6e\n",j,inp_dsr_wlen[j],w[0]);
            err = 1;
            break;
          }
        }
      } else
      if(inp_type[i] == INP_TYPE_AUR) // AUR
      {
        if(inp_aur_wlen[j] < EPSILON)
        {
          inp_aur_wlen[j] = w[0];
        }
        else
        {
          if(fabs(inp_aur_wlen[j]-w[0]) > DELTA)
          {
            fprintf(stderr,"Error, inp_aur_wlen[%d]=%13.6e, w[0]=%13.6e\n",j,inp_aur_wlen[j],w[0]);
            err = 1;
            break;
          }
        }
      } else
      if(inp_type[i] == INP_TYPE_SSR) // SSR
      {
        if(inp_ssr_wlen[j] < EPSILON)
        {
          inp_ssr_wlen[j] = w[0];
        }
        else
        {
          if(fabs(inp_ssr_wlen[j]-w[0]) > DELTA)
          {
            fprintf(stderr,"Error, inp_ssr_wlen[%d]=%13.6e, w[0]=%13.6e\n",j,inp_ssr_wlen[j],w[0]);
            err = 1;
            break;
          }
        }
      }
      inp_data[i][j] = w[1];
      j++;
    }
    fclose(fp);
    if(j != INP_N_WLEN)
    {
      fprintf(stderr,"Error, i=%d j=%d\n",i,j);
      err = 1;
    }
    if(err)
    {
      return -1;
    }
  }

  // DSR wavelength
  if(dsr_ms720 >= 0)
  {
    if(dsr_ms720<1 || dsr_ms720>2)
    {
      fprintf(stderr,"Invalid DSR MS720ID >>> %d\n",dsr_ms720);
      return -1;
    }
  }
  if(dsr_n_wlen[0]==0 && cnt_dsr[0]==1)
  {
    switch(dsr_ms720)
    {
      case 1:
        dsr_wlen[0][ 0] = DSR_WLEN_0_A_001;
        dsr_wlen[0][ 1] = DSR_WLEN_0_A_002;
        dsr_wlen[0][ 2] = DSR_WLEN_0_A_003;
        dsr_wlen[0][ 3] = DSR_WLEN_0_A_004;
        dsr_wlen[0][ 4] = DSR_WLEN_0_A_005;
        dsr_wlen[0][ 5] = DSR_WLEN_0_A_006;
        dsr_wlen[0][ 6] = DSR_WLEN_0_A_007;
        dsr_wlen[0][ 7] = DSR_WLEN_0_A_008;
        dsr_wlen[0][ 8] = DSR_WLEN_0_A_009;
        dsr_n_wlen[0] = DSR_N_WLEN_0;
        break;
      case 2:
        dsr_wlen[0][ 0] = DSR_WLEN_0_B_001;
        dsr_wlen[0][ 1] = DSR_WLEN_0_B_002;
        dsr_wlen[0][ 2] = DSR_WLEN_0_B_003;
        dsr_wlen[0][ 3] = DSR_WLEN_0_B_004;
        dsr_wlen[0][ 4] = DSR_WLEN_0_B_005;
        dsr_wlen[0][ 5] = DSR_WLEN_0_B_006;
        dsr_wlen[0][ 6] = DSR_WLEN_0_B_007;
        dsr_wlen[0][ 7] = DSR_WLEN_0_B_008;
        dsr_wlen[0][ 8] = DSR_WLEN_0_B_009;
        dsr_n_wlen[0] = DSR_N_WLEN_0;
        break;
    }
  }
  if(dsr_n_wlen[1]==0 && cnt_dsr[1]==1)
  {
    switch(dsr_ms720)
    {
      case 1:
        dsr_wlen[1][ 0] = DSR_WLEN_1_A_001;
        dsr_wlen[1][ 1] = DSR_WLEN_1_A_002;
        dsr_wlen[1][ 2] = DSR_WLEN_1_A_003;
        dsr_wlen[1][ 3] = DSR_WLEN_1_A_004;
        dsr_wlen[1][ 4] = DSR_WLEN_1_A_005;
        dsr_wlen[1][ 5] = DSR_WLEN_1_A_006;
        dsr_wlen[1][ 6] = DSR_WLEN_1_A_007;
        dsr_n_wlen[1] = DSR_N_WLEN_1;
        break;
      case 2:
        dsr_wlen[1][ 0] = DSR_WLEN_1_B_001;
        dsr_wlen[1][ 1] = DSR_WLEN_1_B_002;
        dsr_wlen[1][ 2] = DSR_WLEN_1_B_003;
        dsr_wlen[1][ 3] = DSR_WLEN_1_B_004;
        dsr_wlen[1][ 4] = DSR_WLEN_1_B_005;
        dsr_wlen[1][ 5] = DSR_WLEN_1_B_006;
        dsr_wlen[1][ 6] = DSR_WLEN_1_B_007;
        dsr_n_wlen[1] = DSR_N_WLEN_1;
        break;
    }
  }
  if(dsr_n_wlen[2]==0 && cnt_dsr[2]==1)
  {
    switch(dsr_ms720)
    {
      case 1:
        dsr_wlen[2][ 0] = DSR_WLEN_2_A_001;
        dsr_wlen[2][ 1] = DSR_WLEN_2_A_002;
        dsr_wlen[2][ 2] = DSR_WLEN_2_A_003;
        dsr_wlen[2][ 3] = DSR_WLEN_2_A_004;
        dsr_wlen[2][ 4] = DSR_WLEN_2_A_005;
        dsr_wlen[2][ 5] = DSR_WLEN_2_A_006;
        dsr_wlen[2][ 6] = DSR_WLEN_2_A_007;
        dsr_wlen[2][ 7] = DSR_WLEN_2_A_008;
        dsr_wlen[2][ 8] = DSR_WLEN_2_A_009;
        dsr_wlen[2][ 9] = DSR_WLEN_2_A_010;
        dsr_wlen[2][10] = DSR_WLEN_2_A_011;
        dsr_wlen[2][11] = DSR_WLEN_2_A_012;
        dsr_n_wlen[2] = DSR_N_WLEN_2;
        break;
      case 2:
        dsr_wlen[2][ 0] = DSR_WLEN_2_B_001;
        dsr_wlen[2][ 1] = DSR_WLEN_2_B_002;
        dsr_wlen[2][ 2] = DSR_WLEN_2_B_003;
        dsr_wlen[2][ 3] = DSR_WLEN_2_B_004;
        dsr_wlen[2][ 4] = DSR_WLEN_2_B_005;
        dsr_wlen[2][ 5] = DSR_WLEN_2_B_006;
        dsr_wlen[2][ 6] = DSR_WLEN_2_B_007;
        dsr_wlen[2][ 7] = DSR_WLEN_2_B_008;
        dsr_wlen[2][ 8] = DSR_WLEN_2_B_009;
        dsr_wlen[2][ 9] = DSR_WLEN_2_B_010;
        dsr_wlen[2][10] = DSR_WLEN_2_B_011;
        dsr_wlen[2][11] = DSR_WLEN_2_B_012;
        dsr_n_wlen[2] = DSR_N_WLEN_2;
        break;
    }
  }
  if(dsr_n_wlen[3]==0 && cnt_dsr[3]==1)
  {
    switch(dsr_ms720)
    {
      case 1:
        dsr_wlen[3][ 0] = DSR_WLEN_3_A_001;
        dsr_wlen[3][ 1] = DSR_WLEN_3_A_002;
        dsr_wlen[3][ 2] = DSR_WLEN_3_A_003;
        dsr_wlen[3][ 3] = DSR_WLEN_3_A_004;
        dsr_wlen[3][ 4] = DSR_WLEN_3_A_005;
        dsr_wlen[3][ 5] = DSR_WLEN_3_A_006;
        dsr_wlen[3][ 6] = DSR_WLEN_3_A_007;
        dsr_n_wlen[3] = DSR_N_WLEN_3;
        break;
      case 2:
        dsr_wlen[3][ 0] = DSR_WLEN_3_B_001;
        dsr_wlen[3][ 1] = DSR_WLEN_3_B_002;
        dsr_wlen[3][ 2] = DSR_WLEN_3_B_003;
        dsr_wlen[3][ 3] = DSR_WLEN_3_B_004;
        dsr_wlen[3][ 4] = DSR_WLEN_3_B_005;
        dsr_wlen[3][ 5] = DSR_WLEN_3_B_006;
        dsr_wlen[3][ 6] = DSR_WLEN_3_B_007;
        dsr_n_wlen[3] = DSR_N_WLEN_3;
        break;
    }
  }
  for(rng=0; rng<DSR_MAXRANG; rng++)
  {
    for(j=0; j<dsr_n_wlen[rng]; j++)
    {
      dsr_nw[rng][j] = -1;
      for(k=0; k<INP_N_WLEN; k++)
      {
        if(fabs(dsr_wlen[rng][j]-inp_dsr_wlen[k]) < EPSILON)
        {
          dsr_nw[rng][j] = k;
          break;
        }
      }
      if(dsr_nw[rng][j] < 0)
      {
        fprintf(stderr,"Error, faild in finding DSR wavelength >>> %13.4f\n",dsr_wlen[rng][j]);
        return -1;
      }
    }
  }

  // AUR wavelength
  if(aur_ms720 >= 0)
  {
    if(aur_ms720<1 || aur_ms720>2)
    {
      fprintf(stderr,"Invalid AUR MS720ID >>> %d\n",aur_ms720);
      return -1;
    }
    for(rng=0; rng<AUR_MAXRANG; rng++)
    {
      if(cnt_aur_auto[rng] == 1)
      {
        cnt_aur[rng] = 1;
        if(cnt_vb > 0)
        {
          fprintf(stderr,"Activated AUR (%d)\n",rng);
        }
      }
    }
  }
  if(aur_n_wlen[0]==0 && cnt_aur[0]==1)
  {
    switch(aur_ms720)
    {
      case 1:
        aur_wlen[0][0] = AUR_WLEN_0_A;
        aur_sbmp[0][0] = AUR_SBMP;
        aur_n_wlen[0] = AUR_N_WLEN_0;
        aur_n_sbmp[0] = AUR_N_WLEN_0;
        break;
      case 2:
        aur_wlen[0][0] = AUR_WLEN_0_B;
        aur_sbmp[0][0] = AUR_SBMP;
        aur_n_wlen[0] = AUR_N_WLEN_0;
        aur_n_sbmp[0] = AUR_N_WLEN_0;
        break;
    }
  } else
  if(aur_n_sbmp[0] != aur_n_wlen[0])
  {
    if(sim_sbmp_init < 0)
    {
      fprintf(stderr,"Error, aur_n_sbmp[%d]=%d, aur_n_wlen[%d]=%d\n",0,aur_n_sbmp[0],0,aur_n_wlen[0]);
      return -1;
    }
  }
  if(aur_n_wlen[1]==0 && cnt_aur[1]==1)
  {
    switch(aur_ms720)
    {
      case 1:
        aur_wlen[1][ 0] = AUR_WLEN_1_A_001;
        aur_wlen[1][ 1] = AUR_WLEN_1_A_002;
        aur_wlen[1][ 2] = AUR_WLEN_1_A_003;
        aur_wlen[1][ 3] = AUR_WLEN_1_A_004;
        aur_wlen[1][ 4] = AUR_WLEN_1_A_005;
        aur_wlen[1][ 5] = AUR_WLEN_1_A_006;
        aur_wlen[1][ 6] = AUR_WLEN_1_A_007;
        aur_wlen[1][ 7] = AUR_WLEN_1_A_008;
        aur_wlen[1][ 8] = AUR_WLEN_1_A_009;
        aur_wlen[1][ 9] = AUR_WLEN_1_A_010;
        aur_sbmp[1][ 0] = AUR_SBMP_1_A_001;
        aur_sbmp[1][ 1] = AUR_SBMP_1_A_002;
        aur_sbmp[1][ 2] = AUR_SBMP_1_A_003;
        aur_sbmp[1][ 3] = AUR_SBMP_1_A_004;
        aur_sbmp[1][ 4] = AUR_SBMP_1_A_005;
        aur_sbmp[1][ 5] = AUR_SBMP_1_A_006;
        aur_sbmp[1][ 6] = AUR_SBMP_1_A_007;
        aur_sbmp[1][ 7] = AUR_SBMP_1_A_008;
        aur_sbmp[1][ 8] = AUR_SBMP_1_A_009;
        aur_sbmp[1][ 9] = AUR_SBMP_1_A_010;
        aur_n_wlen[1] = AUR_N_WLEN_1;
        aur_n_sbmp[1] = AUR_N_WLEN_1;
        break;
      case 2:
        aur_wlen[1][ 0] = AUR_WLEN_1_B_001;
        aur_wlen[1][ 1] = AUR_WLEN_1_B_002;
        aur_wlen[1][ 2] = AUR_WLEN_1_B_003;
        aur_wlen[1][ 3] = AUR_WLEN_1_B_004;
        aur_wlen[1][ 4] = AUR_WLEN_1_B_005;
        aur_wlen[1][ 5] = AUR_WLEN_1_B_006;
        aur_wlen[1][ 6] = AUR_WLEN_1_B_007;
        aur_wlen[1][ 7] = AUR_WLEN_1_B_008;
        aur_wlen[1][ 8] = AUR_WLEN_1_B_009;
        aur_wlen[1][ 9] = AUR_WLEN_1_B_010;
        aur_sbmp[1][ 0] = AUR_SBMP_1_B_001;
        aur_sbmp[1][ 1] = AUR_SBMP_1_B_002;
        aur_sbmp[1][ 2] = AUR_SBMP_1_B_003;
        aur_sbmp[1][ 3] = AUR_SBMP_1_B_004;
        aur_sbmp[1][ 4] = AUR_SBMP_1_B_005;
        aur_sbmp[1][ 5] = AUR_SBMP_1_B_006;
        aur_sbmp[1][ 6] = AUR_SBMP_1_B_007;
        aur_sbmp[1][ 7] = AUR_SBMP_1_B_008;
        aur_sbmp[1][ 8] = AUR_SBMP_1_B_009;
        aur_sbmp[1][ 9] = AUR_SBMP_1_B_010;
        aur_n_wlen[1] = AUR_N_WLEN_1;
        aur_n_sbmp[1] = AUR_N_WLEN_1;
        break;
    }
  } else
  if(aur_n_sbmp[1] != aur_n_wlen[1])
  {
    if(sim_sbmp_init < 0)
    {
      fprintf(stderr,"Error, aur_n_sbmp[%d]=%d, aur_n_wlen[%d]=%d\n",1,aur_n_sbmp[1],1,aur_n_wlen[1]);
      return -1;
    }
  }
  for(rng=0; rng<AUR_MAXRANG; rng++)
  {
    if(sim_sbmp_init > 0)
    {
      for(j=0; j<aur_n_wlen[rng]; j++)
      {
        aur_sbmp[rng][j] = sim_sbmp_init;
      }
      if(cnt_aur[rng] == 0) continue;
      aur_n_sbmp[rng] = aur_n_wlen[rng];
    }
    for(j=0; j<aur_n_wlen[rng]; j++)
    {
      aur_nw[rng][j] = -1;
      for(k=0; k<INP_N_WLEN; k++)
      {
        if(fabs(aur_wlen[rng][j]-inp_aur_wlen[k]) < EPSILON)
        {
          aur_nw[rng][j] = k;
          break;
        }
      }
      if(aur_nw[rng][j] < 0)
      {
        fprintf(stderr,"Error, faild in finding AUR wavelength >>> %13.4f\n",aur_wlen[rng][j]);
        return -1;
      }
    }
  }

  // SSR wavelength
  if(ssr_ms720 >= 0)
  {
    if(ssr_ms720<1 || ssr_ms720>2)
    {
      fprintf(stderr,"Invalid SSR MS720ID >>> %d\n",ssr_ms720);
      return -1;
    }
  }
  if(ssr_n_wlen[0]==0 && cnt_ssr[0]==1)
  {
    switch(ssr_ms720)
    {
      case 1:
        ssr_wlen[0][0] = SSR_WLEN_0_A;
        ssr_sbmp[0][0] = SSR_SBMP;
        ssr_n_wlen[0] = SSR_N_WLEN_0;
        ssr_n_sbmp[0] = SSR_N_WLEN_0;
        break;
      case 2:
        ssr_wlen[0][0] = SSR_WLEN_0_B;
        ssr_sbmp[0][0] = SSR_SBMP;
        ssr_n_wlen[0] = SSR_N_WLEN_0;
        ssr_n_sbmp[0] = SSR_N_WLEN_0;
        break;
    }
  } else
  if(ssr_n_sbmp[0] != ssr_n_wlen[0])
  {
    if(sim_sbmp_init < 0)
    {
      fprintf(stderr,"Error, ssr_n_sbmp[%d]=%d, ssr_n_wlen[%d]=%d\n",0,ssr_n_sbmp[0],0,ssr_n_wlen[0]);
      return -1;
    }
  }
  if(ssr_n_wlen[1]==0 && cnt_ssr[1]==1)
  {
    switch(ssr_ms720)
    {
      case 1:
        ssr_wlen[1][ 0] = SSR_WLEN_1_A_001;
        ssr_wlen[1][ 1] = SSR_WLEN_1_A_002;
        ssr_wlen[1][ 2] = SSR_WLEN_1_A_003;
        ssr_wlen[1][ 3] = SSR_WLEN_1_A_004;
        ssr_wlen[1][ 4] = SSR_WLEN_1_A_005;
        ssr_wlen[1][ 5] = SSR_WLEN_1_A_006;
        ssr_wlen[1][ 6] = SSR_WLEN_1_A_007;
        ssr_wlen[1][ 7] = SSR_WLEN_1_A_008;
        ssr_wlen[1][ 8] = SSR_WLEN_1_A_009;
        ssr_wlen[1][ 9] = SSR_WLEN_1_A_010;
        ssr_sbmp[1][ 0] = SSR_SBMP_1_A_001;
        ssr_sbmp[1][ 1] = SSR_SBMP_1_A_002;
        ssr_sbmp[1][ 2] = SSR_SBMP_1_A_003;
        ssr_sbmp[1][ 3] = SSR_SBMP_1_A_004;
        ssr_sbmp[1][ 4] = SSR_SBMP_1_A_005;
        ssr_sbmp[1][ 5] = SSR_SBMP_1_A_006;
        ssr_sbmp[1][ 6] = SSR_SBMP_1_A_007;
        ssr_sbmp[1][ 7] = SSR_SBMP_1_A_008;
        ssr_sbmp[1][ 8] = SSR_SBMP_1_A_009;
        ssr_sbmp[1][ 9] = SSR_SBMP_1_A_010;
        ssr_n_wlen[1] = SSR_N_WLEN_1;
        ssr_n_sbmp[1] = SSR_N_WLEN_1;
        break;
      case 2:
        ssr_wlen[1][ 0] = SSR_WLEN_1_B_001;
        ssr_wlen[1][ 1] = SSR_WLEN_1_B_002;
        ssr_wlen[1][ 2] = SSR_WLEN_1_B_003;
        ssr_wlen[1][ 3] = SSR_WLEN_1_B_004;
        ssr_wlen[1][ 4] = SSR_WLEN_1_B_005;
        ssr_wlen[1][ 5] = SSR_WLEN_1_B_006;
        ssr_wlen[1][ 6] = SSR_WLEN_1_B_007;
        ssr_wlen[1][ 7] = SSR_WLEN_1_B_008;
        ssr_wlen[1][ 8] = SSR_WLEN_1_B_009;
        ssr_wlen[1][ 9] = SSR_WLEN_1_B_010;
        ssr_sbmp[1][ 0] = SSR_SBMP_1_B_001;
        ssr_sbmp[1][ 1] = SSR_SBMP_1_B_002;
        ssr_sbmp[1][ 2] = SSR_SBMP_1_B_003;
        ssr_sbmp[1][ 3] = SSR_SBMP_1_B_004;
        ssr_sbmp[1][ 4] = SSR_SBMP_1_B_005;
        ssr_sbmp[1][ 5] = SSR_SBMP_1_B_006;
        ssr_sbmp[1][ 6] = SSR_SBMP_1_B_007;
        ssr_sbmp[1][ 7] = SSR_SBMP_1_B_008;
        ssr_sbmp[1][ 8] = SSR_SBMP_1_B_009;
        ssr_sbmp[1][ 9] = SSR_SBMP_1_B_010;
        ssr_n_wlen[1] = SSR_N_WLEN_1;
        ssr_n_sbmp[1] = SSR_N_WLEN_1;
        break;
    }
  } else
  if(ssr_n_sbmp[1] != ssr_n_wlen[1])
  {
    if(sim_sbmp_init < 0)
    {
      fprintf(stderr,"Error, ssr_n_sbmp[%d]=%d, ssr_n_wlen[%d]=%d\n",1,ssr_n_sbmp[1],1,ssr_n_wlen[1]);
      return -1;
    }
  }
  if(ssr_n_wlen[2]==0 && cnt_ssr[2]==1)
  {
    switch(ssr_ms720)
    {
      case 1:
        ssr_wlen[2][ 0] = SSR_WLEN_2_A_001;
        ssr_wlen[2][ 1] = SSR_WLEN_2_A_002;
        ssr_wlen[2][ 2] = SSR_WLEN_2_A_003;
        ssr_wlen[2][ 3] = SSR_WLEN_2_A_004;
        ssr_wlen[2][ 4] = SSR_WLEN_2_A_005;
        ssr_wlen[2][ 5] = SSR_WLEN_2_A_006;
        ssr_wlen[2][ 6] = SSR_WLEN_2_A_007;
        ssr_wlen[2][ 7] = SSR_WLEN_2_A_008;
        ssr_wlen[2][ 8] = SSR_WLEN_2_A_009;
        ssr_wlen[2][ 9] = SSR_WLEN_2_A_010;
        ssr_wlen[2][10] = SSR_WLEN_2_A_011;
        ssr_wlen[2][11] = SSR_WLEN_2_A_012;
        ssr_sbmp[2][ 0] = SSR_SBMP;
        ssr_sbmp[2][ 1] = SSR_SBMP;
        ssr_sbmp[2][ 2] = SSR_SBMP;
        ssr_sbmp[2][ 3] = SSR_SBMP;
        ssr_sbmp[2][ 4] = SSR_SBMP;
        ssr_sbmp[2][ 5] = SSR_SBMP;
        ssr_sbmp[2][ 6] = SSR_SBMP;
        ssr_sbmp[2][ 7] = SSR_SBMP;
        ssr_sbmp[2][ 8] = SSR_SBMP;
        ssr_sbmp[2][ 9] = SSR_SBMP;
        ssr_sbmp[2][10] = SSR_SBMP;
        ssr_sbmp[2][11] = SSR_SBMP;
        ssr_n_wlen[2] = SSR_N_WLEN_2;
        ssr_n_sbmp[2] = SSR_N_WLEN_2;
        break;
      case 2:
        ssr_wlen[2][ 0] = SSR_WLEN_2_B_001;
        ssr_wlen[2][ 1] = SSR_WLEN_2_B_002;
        ssr_wlen[2][ 2] = SSR_WLEN_2_B_003;
        ssr_wlen[2][ 3] = SSR_WLEN_2_B_004;
        ssr_wlen[2][ 4] = SSR_WLEN_2_B_005;
        ssr_wlen[2][ 5] = SSR_WLEN_2_B_006;
        ssr_wlen[2][ 6] = SSR_WLEN_2_B_007;
        ssr_wlen[2][ 7] = SSR_WLEN_2_B_008;
        ssr_wlen[2][ 8] = SSR_WLEN_2_B_009;
        ssr_wlen[2][ 9] = SSR_WLEN_2_B_010;
        ssr_wlen[2][10] = SSR_WLEN_2_B_011;
        ssr_wlen[2][11] = SSR_WLEN_2_B_012;
        ssr_sbmp[2][ 0] = SSR_SBMP;
        ssr_sbmp[2][ 1] = SSR_SBMP;
        ssr_sbmp[2][ 2] = SSR_SBMP;
        ssr_sbmp[2][ 3] = SSR_SBMP;
        ssr_sbmp[2][ 4] = SSR_SBMP;
        ssr_sbmp[2][ 5] = SSR_SBMP;
        ssr_sbmp[2][ 6] = SSR_SBMP;
        ssr_sbmp[2][ 7] = SSR_SBMP;
        ssr_sbmp[2][ 8] = SSR_SBMP;
        ssr_sbmp[2][ 9] = SSR_SBMP;
        ssr_sbmp[2][10] = SSR_SBMP;
        ssr_sbmp[2][11] = SSR_SBMP;
        ssr_n_wlen[2] = SSR_N_WLEN_2;
        ssr_n_sbmp[2] = SSR_N_WLEN_2;
        break;
    }
  } else
  if(ssr_n_sbmp[2] != ssr_n_wlen[2])
  {
    if(sim_sbmp_init < 0)
    {
      fprintf(stderr,"Error, ssr_n_sbmp[%d]=%d, ssr_n_wlen[%d]=%d\n",2,ssr_n_sbmp[2],2,ssr_n_wlen[2]);
      return -1;
    }
  }
  if(ssr_n_wlen[3]==0 && cnt_ssr[3]==1)
  {
    switch(ssr_ms720)
    {
      case 1:
        ssr_wlen[3][ 0] = SSR_WLEN_3_A_001;
        ssr_wlen[3][ 1] = SSR_WLEN_3_A_002;
        ssr_wlen[3][ 2] = SSR_WLEN_3_A_003;
        ssr_wlen[3][ 3] = SSR_WLEN_3_A_004;
        ssr_wlen[3][ 4] = SSR_WLEN_3_A_005;
        ssr_wlen[3][ 5] = SSR_WLEN_3_A_006;
        ssr_wlen[3][ 6] = SSR_WLEN_3_A_007;
        ssr_sbmp[3][ 0] = SSR_SBMP;
        ssr_sbmp[3][ 1] = SSR_SBMP;
        ssr_sbmp[3][ 2] = SSR_SBMP;
        ssr_sbmp[3][ 3] = SSR_SBMP;
        ssr_sbmp[3][ 4] = SSR_SBMP;
        ssr_sbmp[3][ 5] = SSR_SBMP;
        ssr_sbmp[3][ 6] = SSR_SBMP;
        ssr_n_wlen[3] = SSR_N_WLEN_3;
        ssr_n_sbmp[3] = SSR_N_WLEN_3;
        break;
      case 2:
        ssr_wlen[3][ 0] = SSR_WLEN_3_B_001;
        ssr_wlen[3][ 1] = SSR_WLEN_3_B_002;
        ssr_wlen[3][ 2] = SSR_WLEN_3_B_003;
        ssr_wlen[3][ 3] = SSR_WLEN_3_B_004;
        ssr_wlen[3][ 4] = SSR_WLEN_3_B_005;
        ssr_wlen[3][ 5] = SSR_WLEN_3_B_006;
        ssr_wlen[3][ 6] = SSR_WLEN_3_B_007;
        ssr_sbmp[3][ 0] = SSR_SBMP;
        ssr_sbmp[3][ 1] = SSR_SBMP;
        ssr_sbmp[3][ 2] = SSR_SBMP;
        ssr_sbmp[3][ 3] = SSR_SBMP;
        ssr_sbmp[3][ 4] = SSR_SBMP;
        ssr_sbmp[3][ 5] = SSR_SBMP;
        ssr_sbmp[3][ 6] = SSR_SBMP;
        ssr_n_wlen[3] = SSR_N_WLEN_3;
        ssr_n_sbmp[3] = SSR_N_WLEN_3;
        break;
    }
  } else
  if(ssr_n_sbmp[3] != ssr_n_wlen[3])
  {
    if(sim_sbmp_init < 0)
    {
      fprintf(stderr,"Error, ssr_n_sbmp[%d]=%d, ssr_n_wlen[%d]=%d\n",3,ssr_n_sbmp[3],3,ssr_n_wlen[3]);
      return -1;
    }
  }
  for(rng=0; rng<SSR_MAXRANG; rng++)
  {
    if(sim_sbmp_init > 0)
    {
      for(j=0; j<ssr_n_wlen[rng]; j++)
      {
        ssr_sbmp[rng][j] = sim_sbmp_init;
      }
      if(cnt_ssr[rng] == 0) continue;
      ssr_n_sbmp[rng] = ssr_n_wlen[rng];
    }
    for(j=0; j<ssr_n_wlen[rng]; j++)
    {
      ssr_nw[rng][j] = -1;
      for(k=0; k<INP_N_WLEN; k++)
      {
        if(fabs(ssr_wlen[rng][j]-inp_ssr_wlen[k]) < EPSILON)
        {
          ssr_nw[rng][j] = k;
          break;
        }
      }
      if(ssr_nw[rng][j] < 0)
      {
        fprintf(stderr,"Error, faild in finding SSR wavelength >>> %13.4f\n",ssr_wlen[rng][j]);
        return -1;
      }
    }
  }

  // Error wavelength
  if(dsr_ms720 > 0)
  {
    i = 0;
    for(j=0; j<INP_N_WLEN; j++)
    {
      err_dsr_nw[j] = -1;
      while(i < ERR_N_WLEN)
      {
        if(inp_dsr_wlen[j]>=err_wmin[i] && inp_dsr_wlen[j]<err_wmax[i])
        {
          err_dsr_nw[j] = i;
          break;
        }
        i++;
      }
      if(err_dsr_nw[j] < 0)
      {
        fprintf(stderr,"Error, failed in finding error for %13.4f\n",inp_dsr_wlen[j]);
        return -1;
      }
    }
  }
  if(aur_ms720 > 0)
  {
    i = 0;
    for(j=0; j<INP_N_WLEN; j++)
    {
      err_aur_nw[j] = -1;
      while(i < ERR_N_WLEN)
      {
        if(inp_aur_wlen[j]>=err_wmin[i] && inp_aur_wlen[j]<err_wmax[i])
        {
          err_aur_nw[j] = i;
          break;
        }
        i++;
      }
      if(err_aur_nw[j] < 0)
      {
        fprintf(stderr,"Error, failed in finding error for %13.4f\n",inp_aur_wlen[j]);
        return -1;
      }
    }
  }
  if(ssr_ms720 > 0)
  {
    i = 0;
    for(j=0; j<INP_N_WLEN; j++)
    {
      err_ssr_nw[j] = -1;
      while(i < ERR_N_WLEN)
      {
        if(inp_ssr_wlen[j]>=err_wmin[i] && inp_ssr_wlen[j]<err_wmax[i])
        {
          err_ssr_nw[j] = i;
          break;
        }
        i++;
      }
      if(err_ssr_nw[j] < 0)
      {
        fprintf(stderr,"Error, failed in finding error for %13.4f\n",inp_ssr_wlen[j]);
        return -1;
      }
    }
  }

  // copy data
  dsr_n_data = -1;
  dsr_time_tmp = -1;
  for(i=0; i<DSR_MAXDATA; i++)
  {
    dsr_time[i]   = 0;
    dsr_th_sun[i] = 0.0;
    dsr_ph_sun[i] = 0.0;
    dsr_amas[i]   = 0.0;
    dsr_dist[i]   = 0.0;
    dsr_repeat[i] = 0;
    for(j=0; j<INP_N_WLEN; j++)
    {
      dsr_inp_avg[i][j] = 0.0;
      dsr_inp_sqr[i][j] = 0.0;
    }
  }
  aur_n_data = -1;
  aur_time_tmp = -1;
  for(i=0; i<AUR_MAXDATA; i++)
  {
    aur_time[i]   = 0;
    aur_th_sun[i] = 0.0;
    aur_ph_sun[i] = 0.0;
    aur_amas[i]   = 0.0;
    aur_dist[i]   = 0.0;
    aur_repeat[i] = 0;
    for(j=0; j<INP_N_WLEN; j++)
    {
      aur_inp_avg[i][j] = 0.0;
      aur_inp_sqr[i][j] = 0.0;
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
    ssr_amas[i]   = 0.0;
    ssr_dist[i]   = 0.0;
    ssr_repeat[i] = 0;
    for(j=0; j<INP_N_WLEN; j++)
    {
      ssr_inp_avg[i][j] = 0.0;
      ssr_inp_sqr[i][j] = 0.0;
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
      dsr_time[dsr_n_data]   += (inp_time[i]-inp_time[0]);
      dsr_th_sun[dsr_n_data] += inp_th_sun[i];
      dsr_ph_sun[dsr_n_data] += inp_ph_sun[i];
      dsr_amas[dsr_n_data]   += inp_amas[i];
      dsr_dist[dsr_n_data]   += inp_dist[i];
      for(j=0; j<INP_N_WLEN; j++)
      {
        dsr_inp_avg[dsr_n_data][j] += inp_data[i][j];
        dsr_inp_sqr[dsr_n_data][j] += inp_data[i][j]*inp_data[i][j];
      }
      dsr_repeat[dsr_n_data] += 1;
    } else
    if(inp_type[i] == INP_TYPE_AUR) // AUR
    {
      if(abs(aur_time_tmp-inp_time[i]) > inp_tdif)
      {
        aur_time_tmp = inp_time[i];
        aur_n_data++;
        if(aur_n_data >= AUR_MAXDATA)
        {
          fprintf(stderr,"Warning, #AUR data exceed the limit >>> %d\n",aur_n_data);
          aur_n_data--;
          break;
        }
      }
      aur_time[aur_n_data]   += (inp_time[i]-inp_time[0]);
      aur_th_sun[aur_n_data] += inp_th_sun[i];
      aur_ph_sun[aur_n_data] += inp_ph_sun[i];
      aur_amas[aur_n_data]   += inp_amas[i];
      aur_dist[aur_n_data]   += inp_dist[i];
      for(j=0; j<INP_N_WLEN; j++)
      {
        aur_inp_avg[aur_n_data][j] += inp_data[i][j];
        aur_inp_sqr[aur_n_data][j] += inp_data[i][j]*inp_data[i][j];
      }
      aur_repeat[aur_n_data] += 1;
    } else
    if(inp_type[i] == INP_TYPE_SSR) // SSR
    {
      if(inp_th_los[i]<ssr_th_los_min || inp_th_los[i]>ssr_th_los_max) continue;
      etmp = SepAngle2(inp_th_sun[i]*D_TO_R,inp_ph_sun[i]*D_TO_R,inp_th_los[i]*D_TO_R,inp_ph_los[i]*D_TO_R)*R_TO_D;
      if(etmp<ssr_sep_min || etmp>ssr_sep_max) continue;
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
      ssr_time[ssr_n_data]   += (inp_time[i]-inp_time[0]);
      ssr_th_sun[ssr_n_data] += inp_th_sun[i];
      ssr_ph_sun[ssr_n_data] += inp_ph_sun[i];
      ssr_amas[ssr_n_data]   += inp_amas[i];
      ssr_dist[ssr_n_data]   += inp_dist[i];
      for(j=0; j<INP_N_WLEN; j++)
      {
        ssr_inp_avg[ssr_n_data][j] += inp_data[i][j];
        ssr_inp_sqr[ssr_n_data][j] += inp_data[i][j]*inp_data[i][j];
      }
      ssr_repeat[ssr_n_data] += 1;
    }
  }
  dsr_n_data++;
  aur_n_data++;
  ssr_n_data++;
  // store data&err for DSR
  for(i=0; i<dsr_n_data; i++)
  {
    if(dsr_repeat[i] > 1)
    {
      dsr_time[i]   /= dsr_repeat[i];
      dsr_th_sun[i] /= dsr_repeat[i];
      dsr_ph_sun[i] /= dsr_repeat[i];
      dsr_amas[i]   /= dsr_repeat[i];
      dsr_dist[i]   /= dsr_repeat[i];
      for(j=0; j<INP_N_WLEN; j++)
      {
        dsr_inp_avg[i][j] /= dsr_repeat[i];
        dsr_inp_sqr[i][j] /= dsr_repeat[i];
      }
    }
    for(j=0; j<ERR_N_WLEN; j++)
    {
      err_avg[j] = 0.0;
      err_num[j] = 0;
    }
    for(j=0; j<INP_N_WLEN; j++)
    {
      if((etmp=dsr_inp_sqr[i][j]-dsr_inp_avg[i][j]*dsr_inp_avg[i][j]) < 0.0)
      {
        if(fabs(etmp) > epsilon)
        {
          fprintf(stderr,"Error, dsr_inp_sqr[%d][%d]=%13.6e, dsr_inp_avg[%d][%d]=%13.6e, etmp=%13.6e\n",
                          i,j,dsr_inp_sqr[i][j],i,j,dsr_inp_avg[i][j],etmp);
          return -1;
        }
        etmp = 0.0;
      }
      err_avg[err_dsr_nw[j]] += etmp;
      err_num[err_dsr_nw[j]] += 1;
    }
    for(j=0; j<ERR_N_WLEN; j++)
    {
      if(err_num[j] > 1)
      {
        err_avg[j] /= (double)err_num[j];
      }
    }
    for(j=0; j<INP_N_WLEN; j++)
    {
      err_tmp[j] = sqrt(err_avg[err_dsr_nw[j]]);
      if(cnt_db)
      {
        fprintf(stderr,"%d %3d %9.4f %22.14e %22.14e\n",0,i,inp_dsr_wlen[j],dsr_inp_avg[i][j],err_tmp[j]);
      }
    }
    for(rng=0; rng<DSR_MAXRANG; rng++)
    {
      for(j=0; j<dsr_n_wlen[rng]; j++)
      {
        dsr_data[rng][i][j] = dsr_inp_avg[i][dsr_nw[rng][j]];
        dsr_derr[rng][i][j] = sqrt(err_tmp[dsr_nw[rng][j]]*err_tmp[dsr_nw[rng][j]]+
                                  (dsr_rerr*dsr_data[rng][i][0])*(dsr_rerr*dsr_data[rng][i][0])); // independent of wlen
      }
    }
  }
  // store data&err for AUR
  for(i=0; i<aur_n_data; i++)
  {
    if(aur_repeat[i] > 1)
    {
      aur_time[i]   /= aur_repeat[i];
      aur_th_sun[i] /= aur_repeat[i];
      aur_ph_sun[i] /= aur_repeat[i];
      aur_amas[i]   /= aur_repeat[i];
      aur_dist[i]   /= aur_repeat[i];
      for(j=0; j<INP_N_WLEN; j++)
      {
        aur_inp_avg[i][j] /= aur_repeat[i];
        aur_inp_sqr[i][j] /= aur_repeat[i];
      }
    }
    if(SepToAzi(aur_th_sun[i]*D_TO_R,aur_rfov*D_TO_R,&etmp) < 0)
    {
      fprintf(stderr,"Error in AUR geometry.\n");
      return -1;
    }
    aur_th_los[i] = aur_th_sun[i];
    aur_ph_los[i] = aur_ph_sun[i]-etmp*R_TO_D;
    for(j=0; j<ERR_N_WLEN; j++)
    {
      err_avg[j] = 0.0;
      err_num[j] = 0;
    }
    for(j=0; j<INP_N_WLEN; j++)
    {
      if((etmp=aur_inp_sqr[i][j]-aur_inp_avg[i][j]*aur_inp_avg[i][j]) < 0.0)
      {
        if(fabs(etmp) > epsilon)
        {
          fprintf(stderr,"Error, aur_inp_sqr[%d][%d]=%13.6e, aur_inp_avg[%d][%d]=%13.6e, etmp=%13.6e\n",
                          i,j,aur_inp_sqr[i][j],i,j,aur_inp_avg[i][j],etmp);
          return -1;
        }
        etmp = 0.0;
      }
      err_avg[err_aur_nw[j]] += etmp;
      err_num[err_aur_nw[j]] += 1;
    }
    for(j=0; j<ERR_N_WLEN; j++)
    {
      if(err_num[j] > 1)
      {
        err_avg[j] /= (double)err_num[j];
      }
    }
    for(j=0; j<INP_N_WLEN; j++)
    {
      err_tmp[j] = sqrt(err_avg[err_aur_nw[j]]);
      if(cnt_db)
      {
        fprintf(stderr,"%d %3d %9.4f %22.14e %22.14e\n",1,i,inp_aur_wlen[j],aur_inp_avg[i][j],err_tmp[j]);
      }
    }
    for(rng=0; rng<AUR_MAXRANG; rng++)
    {
      for(j=0; j<aur_n_wlen[rng]; j++)
      {
        aur_data[rng][i][j] = aur_inp_avg[i][aur_nw[rng][j]];
        aur_derr[rng][i][j] = sqrt(err_tmp[aur_nw[rng][j]]*err_tmp[aur_nw[rng][j]]+
                                  (aur_rerr*aur_data[rng][i][0])*(aur_rerr*aur_data[rng][i][0])); // independent of wlen
      }
    }
  }
  if(aur_type == AUR_TYPE_OLD)
  {
    // subtract DSR
    for(i=0; i<aur_n_data; i++)
    {
      aur_ndsr[i] = -1;
      for(n=0; n<dsr_n_data; n++)
      {
        if(abs(aur_time[i]-dsr_time[n]) < inp_tdif)
        {
          aur_ndsr[i] = n;
          break;
        }
      }
      if(aur_ndsr[i] < 0)
      {
        fprintf(stderr,"Error, failed in finding DSR data for i=%d\n",i);
        return -1;
      }
      if(aur_ms720 != dsr_ms720)
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          Interp1D(inp_dsr_wlen,dsr_inp_avg[aur_ndsr[i]],INP_N_WLEN,inp_aur_wlen,dsr_aur_avg,INP_N_WLEN,1);
          Interp1D(inp_dsr_wlen,dsr_inp_sqr[aur_ndsr[i]],INP_N_WLEN,inp_aur_wlen,dsr_aur_sqr,INP_N_WLEN,0);
          for(j=0; j<aur_n_wlen[rng]; j++)
          {
            aur_data[rng][i][j] -= dsr_aur_avg[aur_nw[rng][j]];
            aur_derr[rng][i][j] = sqrt(aur_derr[rng][i][j]*aur_derr[rng][i][j]+
                                      dsr_aur_sqr[aur_nw[rng][j]]-
                                      dsr_aur_avg[aur_nw[rng][j]]*dsr_aur_avg[aur_nw[rng][j]]+
                                      (dsr_rerr*dsr_aur_avg[aur_nw[rng][j]])*
                                      (dsr_rerr*dsr_aur_avg[aur_nw[rng][j]]));
          }
        }
      }
      else
      {
        for(rng=0; rng<AUR_MAXRANG; rng++)
        {
          for(j=0; j<aur_n_wlen[rng]; j++)
          {
            if(fabs(inp_aur_wlen[aur_nw[rng][j]]-inp_dsr_wlen[aur_nw[rng][j]]) > DELTA)
            {
              fprintf(stderr,"Error, inp_aur_wlen[%d]=%13.4f, inp_dsr_wlen[%d]=%13.4f\n",
                              aur_nw[rng][j],inp_aur_wlen[aur_nw[rng][j]],
                              aur_nw[rng][j],inp_dsr_wlen[aur_nw[rng][j]]);
              return -1;
            }
            aur_data[rng][i][j] -= dsr_inp_avg[aur_ndsr[i]][aur_nw[rng][j]];
            aur_derr[rng][i][j] = sqrt(aur_derr[rng][i][j]*aur_derr[rng][i][j]+
                                      dsr_inp_sqr[aur_ndsr[i]][aur_nw[rng][j]]-
                                      dsr_inp_avg[aur_ndsr[i]][aur_nw[rng][j]]*dsr_inp_avg[aur_ndsr[i]][aur_nw[rng][j]]+
                                      (dsr_rerr*dsr_inp_avg[aur_ndsr[i]][aur_nw[rng][j]])*
                                      (dsr_rerr*dsr_inp_avg[aur_ndsr[i]][aur_nw[rng][j]]));
          }
        }
      }
    }
  }
  // Irradiance -> Radiance
  for(i=0; i<aur_n_data; i++)
  {
    for(rng=0; rng<AUR_MAXRANG; rng++)
    {
      for(j=0; j<aur_n_wlen[rng]; j++)
      {
        aur_data[rng][i][j] /= aur_sfov;
        aur_derr[rng][i][j] /= aur_sfov;
      }
    }
  }
  // store data&err for SSR
  for(i=0; i<ssr_n_data; i++)
  {
    if(ssr_repeat[i] > 1)
    {
      ssr_time[i]   /= ssr_repeat[i];
      ssr_th_sun[i] /= ssr_repeat[i];
      ssr_ph_sun[i] /= ssr_repeat[i];
      ssr_amas[i]   /= ssr_repeat[i];
      ssr_dist[i]   /= ssr_repeat[i];
      for(j=0; j<INP_N_WLEN; j++)
      {
        ssr_inp_avg[i][j] /= ssr_repeat[i];
        ssr_inp_sqr[i][j] /= ssr_repeat[i];
      }
    }
    for(j=0; j<ERR_N_WLEN; j++)
    {
      err_avg[j] = 0.0;
      err_num[j] = 0;
    }
    for(j=0; j<INP_N_WLEN; j++)
    {
      if((etmp=ssr_inp_sqr[i][j]-ssr_inp_avg[i][j]*ssr_inp_avg[i][j]) < 0.0)
      {
        if(fabs(etmp) > epsilon)
        {
          fprintf(stderr,"Error, ssr_inp_sqr[%d][%d]=%13.6e, ssr_inp_avg[%d][%d]=%13.6e, etmp=%13.6e\n",
                          i,j,ssr_inp_sqr[i][j],i,j,ssr_inp_avg[i][j],etmp);
          return -1;
        }
        etmp = 0.0;
      }
      err_avg[err_ssr_nw[j]] += etmp;
      err_num[err_ssr_nw[j]] += 1;
    }
    for(j=0; j<ERR_N_WLEN; j++)
    {
      if(err_num[j] > 1)
      {
        err_avg[j] /= (double)err_num[j];
      }
    }
    for(j=0; j<INP_N_WLEN; j++)
    {
      err_tmp[j] = sqrt(err_avg[err_ssr_nw[j]]);
      if(cnt_db)
      {
        fprintf(stderr,"%d %3d %9.4f %22.14e %22.14e\n",2,i,inp_ssr_wlen[j],ssr_inp_avg[i][j],err_tmp[j]);
      }
    }
    for(rng=0; rng<SSR_MAXRANG; rng++)
    {
      for(j=0; j<ssr_n_wlen[rng]; j++)
      {
        ssr_data[rng][i][j] = ssr_inp_avg[i][ssr_nw[rng][j]];
        ssr_derr[rng][i][j] = sqrt(err_tmp[ssr_nw[rng][j]]*err_tmp[ssr_nw[rng][j]]+
                                  (ssr_rerr*ssr_data[rng][i][0])*(ssr_rerr*ssr_data[rng][i][0])); // independent of wlen
      }
    }
  }
  // Irradiance -> Radiance
  // (this procedure was incorrectly inserted in the j-loop above until Feb. 26, 2012)
  for(i=0; i<ssr_n_data; i++)
  {
    for(rng=0; rng<SSR_MAXRANG; rng++)
    {
      for(j=0; j<ssr_n_wlen[rng]; j++)
      {
        ssr_data[rng][i][j] /= ssr_sfov;
        ssr_derr[rng][i][j] /= ssr_sfov;
      }
    }
  }

  // save the original data
  for(i=0; i<dsr_n_data; i++)
  {
    for(rng=0; rng<DSR_MAXRANG; rng++)
    {
      for(j=0; j<dsr_n_wlen[rng]; j++)
      {
        dsr_dorg[rng][i][j] = dsr_data[rng][i][j];
      }
    }
  }
  for(i=0; i<aur_n_data; i++)
  {
    for(rng=0; rng<AUR_MAXRANG; rng++)
    {
      for(j=0; j<aur_n_wlen[rng]; j++)
      {
        aur_dorg[rng][i][j] = aur_data[rng][i][j];
      }
    }
  }
  for(i=0; i<ssr_n_data; i++)
  {
    for(rng=0; rng<SSR_MAXRANG; rng++)
    {
      for(j=0; j<ssr_n_wlen[rng]; j++)
      {
        ssr_dorg[rng][i][j] = ssr_data[rng][i][j];
      }
    }
  }

  return 0;
}

int ReadCRC(void)
{
  int i,n;
  int rng;

  // read simulation data
  if(strcmp(crc_fsim,NONAME) != 0)
  {
    if(ReadSimul(crc_fsim,2,1,1,0,1) < 0) return -1;
    crc_n_dsim = ssr_n_data;
    for(n=0; n<dsr_n_data; n++)
    {
      for(rng=0; rng<DSR_MAXRANG; rng++)
      {
        for(i=0; i<dsr_n_wlen[rng]; i++)
        {
          crc_dsr_dsim[rng][n][i] = dsr_dsim[rng][n][i];
        }
        if(cnt_sflg == 1)
        {
          for(i=0; i<dsr_n_wlen[rng]; i++)
          {
            crc_dsr_dsgl[rng][n][i] = dsr_dsgl[rng][n][i];
          }
        }
      }
    }
    for(n=0; n<aur_n_data; n++)
    {
      for(rng=0; rng<AUR_MAXRANG; rng++)
      {
        for(i=0; i<aur_n_wlen[rng]; i++)
        {
          crc_aur_dsim[rng][n][i] = aur_dsim[rng][n][i];
        }
        if(cnt_sflg == 1)
        {
          for(i=0; i<aur_n_wlen[rng]; i++)
          {
            crc_aur_dsgl[rng][n][i] = aur_dsgl[rng][n][i];
          }
        }
      }
    }
    for(n=0; n<ssr_n_data; n++)
    {
      for(rng=0; rng<SSR_MAXRANG; rng++)
      {
        for(i=0; i<ssr_n_wlen[rng]; i++)
        {
          crc_ssr_dsim[rng][n][i] = ssr_dsim[rng][n][i];
        }
        if(cnt_sflg == 1)
        {
          for(i=0; i<ssr_n_wlen[rng]; i++)
          {
            crc_ssr_dsgl[rng][n][i] = ssr_dsgl[rng][n][i];
          }
        }
      }
    }
  }

  return 0;
}

int ReadFactor(void)
{
  int i,j;
  int rng;
  int rid;
  int n,nc;
  int err;
  int v[2];
  double w[2];
  char line[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char fnam[] = "ReadFactor";
  char *p;
  FILE *fp;

  if(strcmp(crc_ffac,"")==0 || strcmp(crc_ffac,NONAME)==0)
  {
    return 0;
  } else
  if((fp=fopen(crc_ffac,"r")) == NULL)
  {
    fprintf(stderr,"%s: error, cannot open %s\n",fnam,crc_ffac);
    return -1;
  } else
  if(cnt_vb)
  {
    fprintf(stderr,"%s: reading %s\n",fnam,crc_ffac);
  }
  err = 0;
  // Read factors for DSR
  for(j=0; j<dsr_n_data; j++)
  {
    for(rng=0; rng<DSR_MAXRANG; rng++)
    {
      rid = rng;
      for(i=0; i<dsr_n_wlen[rng]; i++)
      {
        if(fgets(line,MAXLINE,fp) == NULL)
        {
          fprintf(stderr,"%s: read error (i=%d,j=%d,rng=%d)\n",fnam,i,j,rng);
          err = 1;
          break;
        }
        // read line
        for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
        {
          if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
        }
        if(n != 4)
        {
          fprintf(stderr,"%s: error, expected %d columns >>> %d\n",fnam,4,n);
          err = 1;
          break;
        }
        // convert values
        for(n=0; n<2; n++)
        {
          errno = 0;
          v[n] = strtol(str[n],&p,10);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
            err = 1;
            break;
          }
        }
        if(err) break;
        for(n=2; n<4; n++)
        {
          errno = 0;
          w[n-2] = strtod(str[n],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
            err = 1;
            break;
          }
        }
        if(err) break;
        if(v[0] != rid)
        {
          fprintf(stderr,"%s: error, v[0]=%d, rid=%d\n",fnam,v[0],rid);
          err = 1;
          break;
        } else
        if(v[1] != j)
        {
          fprintf(stderr,"%s: error, v[1]=%d, j=%d\n",fnam,v[1],j);
          err = 1;
          break;
        } else
        if(fabs(w[0]-dsr_wlen[rng][i]) > DELTA)
        {
          fprintf(stderr,"%s: error, w[0]=%13.6e, dsr_wlen[%d][%d]=%13.6e\n",fnam,w[0],rng,i,dsr_wlen[rng][i]);
          err = 1;
          break;
        }
        else
        {
          dsr_dfac[rng][j][i] = w[1];
        }
      }
      if(err) break;
    }
    if(err) break;
  }
  if(err)
  {
    fclose(fp);
    return -1;
  }
  // Read factors for AUR
  for(j=0; j<aur_n_data; j++)
  {
    for(rng=0; rng<AUR_MAXRANG; rng++)
    {
      rid = rng+DSR_MAXRANG*2;
      for(i=0; i<aur_n_wlen[rng]; i++)
      {
        if(fgets(line,MAXLINE,fp) == NULL)
        {
          fprintf(stderr,"%s: read error (i=%d,j=%d,rng=%d)\n",fnam,i,j,rng);
          err = 1;
          break;
        }
        // read line
        for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
        {
          if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
        }
        if(n != 4)
        {
          fprintf(stderr,"%s: error, expected %d columns >>> %d\n",fnam,4,n);
          err = 1;
          break;
        }
        // convert values
        for(n=0; n<2; n++)
        {
          errno = 0;
          v[n] = strtol(str[n],&p,10);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
            err = 1;
            break;
          }
        }
        if(err) break;
        for(n=2; n<4; n++)
        {
          errno = 0;
          w[n-2] = strtod(str[n],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
            err = 1;
            break;
          }
        }
        if(err) break;
        if(v[0] != rid)
        {
          fprintf(stderr,"%s: error, v[0]=%d, rid=%d\n",fnam,v[0],rid);
          err = 1;
          break;
        } else
        if(v[1] != j)
        {
          fprintf(stderr,"%s: error, v[1]=%d, j=%d\n",fnam,v[1],j);
          err = 1;
          break;
        } else
        if(fabs(w[0]-aur_wlen[rng][i]) > DELTA)
        {
          fprintf(stderr,"%s: error, w[0]=%13.6e, aur_wlen[%d][%d]=%13.6e\n",fnam,w[0],rng,i,aur_wlen[rng][i]);
          err = 1;
          break;
        }
        else
        {
          aur_dfac[rng][j][i] = w[1];
        }
      }
      if(err) break;
    }
    if(err) break;
  }
  if(err)
  {
    fclose(fp);
    return -1;
  }
  // Read factors for SSR
  for(j=0; j<ssr_n_data; j++)
  {
    for(rng=0; rng<SSR_MAXRANG; rng++)
    {
      rid = rng+DSR_MAXRANG*2+AUR_MAXRANG*2;
      for(i=0; i<ssr_n_wlen[rng]; i++)
      {
        if(fgets(line,MAXLINE,fp) == NULL)
        {
          fprintf(stderr,"%s: read error (i=%d,j=%d,rng=%d)\n",fnam,i,j,rng);
          err = 1;
          break;
        }
        // read line
        for(n=nc=0,p=line; n<MAXENTR; n++,p+=nc)
        {
          if(sscanf(p,"%s%n",str[n],&nc) == EOF) break;
        }
        if(n != 4)
        {
          fprintf(stderr,"%s: error, expected %d columns >>> %d\n",fnam,4,n);
          err = 1;
          break;
        }
        // convert values
        for(n=0; n<2; n++)
        {
          errno = 0;
          v[n] = strtol(str[n],&p,10);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
            err = 1;
            break;
          }
        }
        if(err) break;
        for(n=2; n<4; n++)
        {
          errno = 0;
          w[n-2] = strtod(str[n],&p);
          if(errno==ERANGE || *p!='\0')
          {
            fprintf(stderr,"%s: convert error >>> %s\n",fnam,str[n]);
            err = 1;
            break;
          }
        }
        if(err) break;
        if(v[0] != rid)
        {
          fprintf(stderr,"%s: error, v[0]=%d, rid=%d\n",fnam,v[0],rid);
          err = 1;
          break;
        } else
        if(v[1] != j)
        {
          fprintf(stderr,"%s: error, v[1]=%d, j=%d\n",fnam,v[1],j);
          err = 1;
          break;
        } else
        if(fabs(w[0]-ssr_wlen[rng][i]) > DELTA)
        {
          fprintf(stderr,"%s: error, w[0]=%13.6e, ssr_wlen[%d][%d]=%13.6e\n",fnam,w[0],rng,i,ssr_wlen[rng][i]);
          err = 1;
          break;
        }
        else
        {
          ssr_dfac[rng][j][i] = w[1];
        }
      }
      if(err) break;
    }
    if(err) break;
  }
  if(err)
  {
    fclose(fp);
    return -1;
  }
  fclose(fp);

  return 1;
}

int ReadXA(char *s,int size,int cx,double ux,int (*func)(double*),double *x)
{
  int nd;
  int np = 1;
  int c[1];
  double u[1];
  double *d[1];

  // fill array
  c[0] = cx; u[0] = ux; d[0] = x;
  // read input
  if(func == NULL)
  {
    nd = ReadDA(s,size,np,c,u,d);
  }
  else
  {
    nd = ReadDC(s,size,np,c,u,func,d);
  }

  return nd;
}

int ReadXYA(char *s,int size,int cx,int cy,double ux,double uy,int (*func)(double*),double *x,double *y)
{
  int nd;
  int np = 2;
  int c[2];
  double u[2];
  double *d[2];

  // fill array
  c[0] = cx; u[0] = ux; d[0] = x;
  c[1] = cy; u[1] = uy; d[1] = y;
  // read input
  if(func == NULL)
  {
    nd = ReadDA(s,size,np,c,u,d);
  }
  else
  {
    nd = ReadDC(s,size,np,c,u,func,d);
  }

  return nd;
}

int ReadXYZA(char *s,int size,int cx,int cy,int cz,double ux,double uy,double uz,int (*func)(double*),
             double *x,double *y,double *z)
{
  int nd;
  int np = 3;
  int c[3];
  double u[3];
  double *d[3];

  // fill array
  c[0] = cx; u[0] = ux; d[0] = x;
  c[1] = cy; u[1] = uy; d[1] = y;
  c[2] = cz; u[2] = uz; d[2] = z;
  // read input
  if(func == NULL)
  {
    nd = ReadDA(s,size,np,c,u,d);
  }
  else
  {
    nd = ReadDC(s,size,np,c,u,func,d);
  }

  return nd;
}

int ReadXP(char *s,int size,int cx,double ux,int (*func)(double*),double **x)
{
  int n;
  int nd;
  int err;
  int np = 1;
  int c[1];
  double u[1];
  double *d[1];
  double **p[1];
  char fnam[] = "ReadXP";

  // fill array
  c[0] = cx; u[0] = ux; d[0] = *x; p[0] = x;
  // allocate memory
  for(n=0; n<np; n++)
  {
    if(d[n] == NULL)
    {
      if((d[n]=(double*)malloc(size*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: failed in allocating memory (%s)\n",fnam,s);
        return -1;
      }
    }
  }
  // read input
  if(func == NULL)
  {
    nd = ReadDA(s,size,np,c,u,d);
  }
  else
  {
    nd = ReadDC(s,size,np,c,u,func,d);
  }
  if(nd < 0)
  {
    for(n=0; n<np; n++)
    {
      if(*p[n] == NULL) free(d[n]);
    }
    return -1;
  }
  // resize memory
  err = 0;
  for(n=0; n<np; n++)
  {
    if(*p[n] == NULL)
    {
      if((*p[n]=(double*)realloc(d[n],nd*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: failed in allocating memory (%s)\n",fnam,s);
        free(d[n]); // this is needed only when realloc failed
        err = 1;
      }
    }
  }
  if(err)
  {
    return -1;
  }

  return nd;
}

int ReadXYP(char *s,int size,int cx,int cy,double ux,double uy,int (*func)(double*),double **x,double **y)
{
  int n;
  int nd;
  int err;
  int np = 2;
  int c[2];
  double u[2];
  double *d[2];
  double **p[2];
  char fnam[] = "ReadXYP";

  // fill array
  c[0] = cx; u[0] = ux; d[0] = *x; p[0] = x;
  c[1] = cy; u[1] = uy; d[1] = *y; p[1] = y;
  // allocate memory
  for(n=0; n<np; n++)
  {
    if(d[n] == NULL)
    {
      if((d[n]=(double*)malloc(size*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: failed in allocating memory (%s)\n",fnam,s);
        return -1;
      }
    }
  }
  // read input
  if(func == NULL)
  {
    nd = ReadDA(s,size,np,c,u,d);
  }
  else
  {
    nd = ReadDC(s,size,np,c,u,func,d);
  }
  if(nd < 0)
  {
    for(n=0; n<np; n++)
    {
      if(*p[n] == NULL) free(d[n]);
    }
    return -1;
  }
  // resize memory
  err = 0;
  for(n=0; n<np; n++)
  {
    if(*p[n] == NULL)
    {
      if((*p[n]=(double*)realloc(d[n],nd*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: failed in allocating memory (%s)\n",fnam,s);
        free(d[n]); // this is needed only when realloc failed
        err = 1;
      }
    }
  }
  if(err)
  {
    return -1;
  }

  return nd;
}

int ReadXYZP(char *s,int size,int cx,int cy,int cz,double ux,double uy,double uz,int (*func)(double*),
             double **x,double **y,double **z)
{
  int n;
  int nd;
  int err;
  int np = 3;
  int c[3];
  double u[3];
  double *d[3];
  double **p[3];
  char fnam[] = "ReadXYZP";

  // fill array
  c[0] = cx; u[0] = ux; d[0] = *x; p[0] = x;
  c[1] = cy; u[1] = uy; d[1] = *y; p[1] = y;
  c[2] = cz; u[2] = uz; d[2] = *z; p[2] = z;
  // allocate memory
  for(n=0; n<np; n++)
  {
    if(d[n] == NULL)
    {
      if((d[n]=(double*)malloc(size*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: failed in allocating memory (%s)\n",fnam,s);
        return -1;
      }
    }
  }
  // read input
  if(func == NULL)
  {
    nd = ReadDA(s,size,np,c,u,d);
  }
  else
  {
    nd = ReadDC(s,size,np,c,u,func,d);
  }
  if(nd < 0)
  {
    for(n=0; n<np; n++)
    {
      if(*p[n] == NULL) free(d[n]);
    }
    return -1;
  }
  // resize memory
  err = 0;
  for(n=0; n<np; n++)
  {
    if(*p[n] == NULL)
    {
      if((*p[n]=(double*)realloc(d[n],nd*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: failed in allocating memory (%s)\n",fnam,s);
        free(d[n]); // this is needed only when realloc failed
        err = 1;
      }
    }
  }
  if(err)
  {
    return -1;
  }

  return nd;
}

int ReadDC(char *s,int size,int np,const int *c,const double *u,int (*func)(double*),double **d)
{
  int i,j,n;
  int nc,nd;
  int num;
  int err;
  int flag;
  double v[3];
  char line[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char fnam[] = "ReadDC";
  char *p;
  FILE *fp;

  // check input
  if(np<0 || np>3)
  {
    fprintf(stderr,"%s: error, np=%d\n",fnam,np);
    return -1;
  }
  for(n=0; n<np; n++)
  {
    if(c[n]<0 || c[n]>=MAXENTR)
    {
      fprintf(stderr,"%s: invalid column number >>> c[%d]=%d (%s)\n",fnam,n,c[n],s);
      return -1;
    }
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
      fprintf(stderr,"%s: cannot open %s\n",fnam,s);
      return -1;
    }
    flag = 1;
  }
  // read input
  num = c[0];
  for(n=1; n<np; n++)
  {
    num = MAX(num,c[n]);
  }
  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    for(j=0,p=line; j<=num; j++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[j],&nc) == EOF)
      {
        fprintf(stderr,"%s: read error >>> %s (%s)\n",fnam,line,s);
        err = 1;
        break;
      }
    }
    if(err) break;
    for(n=0; n<np; n++)
    {
      errno = 0;
      v[n] = strtod(str[c[n]],&p)*u[n];
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error >>> %s (%s)\n",fnam,line,s);
        err = 1;
        break;
      }
    }
    if(func(v))
    {
      continue;
    }
    if(i >= size)
    {
      fprintf(stderr,"%s: warning, #data exceed the limit %d (%s)\n",fnam,i,s);
      break;
    }
    for(n=0; n<np; n++)
    {
      d[n][i] = v[n];
    }
    i++;
  }
  if(flag)
  {
    fclose(fp);
  }
  nd = i;
  if(nd < 1)
  {
    fprintf(stderr,"%s: error, nd=%d (%s)\n",fnam,nd,s);
    err = 1;
  }
  if(err)
  {
    return -1;
  }

  return nd;
}

int ReadDA(char *s,int size,int np,const int *c,const double *u,double **d)
{
  int i,j,n;
  int nc,nd;
  int num;
  int err;
  int flag;
  double v[3];
  char line[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char fnam[] = "ReadDA";
  char *p;
  FILE *fp;

  // check input
  if(np<0 || np>3)
  {
    fprintf(stderr,"%s: error, np=%d\n",fnam,np);
    return -1;
  }
  for(n=0; n<np; n++)
  {
    if(c[n]<0 || c[n]>=MAXENTR)
    {
      fprintf(stderr,"%s: invalid column number >>> c[%d]=%d (%s)\n",fnam,n,c[n],s);
      return -1;
    }
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
      fprintf(stderr,"%s: cannot open %s\n",fnam,s);
      return -1;
    }
    flag = 1;
  }
  // read input
  num = c[0];
  for(n=1; n<np; n++)
  {
    num = MAX(num,c[n]);
  }
  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    for(j=0,p=line; j<=num; j++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[j],&nc) == EOF)
      {
        fprintf(stderr,"%s: read error >>> %s (%s)\n",fnam,line,s);
        err = 1;
        break;
      }
    }
    if(err) break;
    for(n=0; n<np; n++)
    {
      errno = 0;
      v[n] = strtod(str[c[n]],&p)*u[n];
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"%s: convert error >>> %s (%s)\n",fnam,line,s);
        err = 1;
        break;
      }
    }
    if(i >= size)
    {
      fprintf(stderr,"%s: warning, #data exceed the limit %d (%s)\n",fnam,i,s);
      break;
    }
    for(n=0; n<np; n++)
    {
      d[n][i] = v[n];
    }
    i++;
  }
  if(flag)
  {
    fclose(fp);
  }
  nd = i;
  if(nd < 1)
  {
    fprintf(stderr,"%s: error, nd=%d (%s)\n",fnam,nd,s);
    err = 1;
  }
  if(err)
  {
    return -1;
  }

  return nd;
}

int ReadNA(char *s,int size,int c,int *d)
{
  int i,j;
  int nc,nd;
  int err;
  int flag;
  int ntmp;
  char line[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char fnam[] = "ReadNA";
  char *p;
  FILE *fp;

  // check input
  if(c<0 || c>=MAXENTR)
  {
    fprintf(stderr,"%s: invalid column number >>> %d (%s)\n",fnam,c,s);
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
      fprintf(stderr,"%s: cannot open %s\n",fnam,s);
      return -1;
    }
    flag = 1;
  }
  // read input
  i = 0;
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    for(j=0,p=line; j<=c; j++,p+=nc)
    {
      if(sscanf(p,"%s%n",str[j],&nc) == EOF)
      {
        fprintf(stderr,"%s: read error >>> %s (%s)\n",fnam,line,s);
        err = 1;
        break;
      }
    }
    if(err) break;
    errno = 0;
    ntmp = strtol(str[c],&p,10);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"%s: convert error >>> %s (%s)\n",fnam,line,s);
      err = 1;
      break;
    }
    if(i >= size)
    {
      fprintf(stderr,"%s: warning, #data exceed the limit %d (%s)\n",fnam,i,s);
      break;
    }
    d[i] = ntmp;
    i++;
  }
  if(flag)
  {
    fclose(fp);
  }
  nd = i;
  if(nd < 1)
  {
    fprintf(stderr,"%s: error, nd=%d (%s)\n",fnam,nd,s);
    err = 1;
  }
  if(err)
  {
    return -1;
  }

  return nd;
}

int WaveSelect(double *w)
{
  if(w[0]<mie_wlen_min || w[0]>mie_wlen_max)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int AnglSelect(double *a)
{
  if(a[0]<mie_angl_min || a[0]>mie_angl_max)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int ReadComp(char *s,int size,int cr,int ci,int imin,int imax,double u,
             double *wcom,double *lmod,double *lsgm,double **refr,double **refi)
{
  int i,j,n;
  int nc;
  int n_comp = 0;
  int err;
  double v[3];
  double *x,*y;
  char line[MAXLINE];
  //char cmnt[MIE_MAXCOMP][MAXLINE];
  char str[6][MAXLINE];
  char fnam[] = "ReadComp";
  char *p;
  FILE *fp;

  if(strcmp(s,"")!=0 && strcmp(s,NONAME)!=0)
  {
    if((fp=fopen(s,"r")) == NULL)
    {
      fprintf(stderr,"%s: cannot open %s\n",fnam,s);
      return -1;
    }
    i = n = 0;
    err = 0;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      if(i < imin)
      {
        strncpy(line,"",MAXLINE);
        i++;
        continue;
      } else
      if(i > imax)
      {
        break;
      }
      if(sscanf(line,"%s%s%s%s%s%n",str[0],str[1],str[2],str[3],str[4],&nc) == EOF)
      {
        fprintf(stderr,"%s: read error (%s) >>> %s\n",fnam,s,line);
        err = 1;
        break;
      }
      strcpy(str[5],line+nc);
      errno = 0;
      for(j=0; j<3; j++)
      {
        v[j] = strtod(str[j],&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"%s: convert error (%s) >>> %s\n",fnam,s,line);
          err = 1;
          break;
        }
      }
      if(err) break;
      if(n >= size)
      {
        fprintf(stderr,"%s: warning, #data exceed the limit %d (%s)\n",fnam,n,s);
        break;
      }
      // Mixing ratio
      wcom[n] = v[0];
      // Mode radius
      if(v[1] < 0.0)
      {
        fprintf(stderr,"%s: error in mode radius (%s) >>> %13.6e\n",fnam,s,v[1]);
        err = 1;
        break;
      }
      else
      {
        lmod[n] = log10(v[1]);
      }
      // Sigma
      lsgm[n] = v[2];
      // Refractive index (real)
      x = NULL;
      y = NULL;
      if((j=ReadXYP(str[3],MIE_MAXDATA,0,cr,u,1.0,NULL,&x,&y)) < 0)
      {
        return -1;
      }
      refr[n] = NULL;
      if((refr[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: error in allocating memory.\n",fnam);
        return -1;
      }
      if(Sampling(j,x,y,mie_n_wlen,mie_wlen,refr[n],1.0) < 0)
      {
        return -1;
      }
      free(x);
      free(y);
      // Refractive index (imag)
      x = NULL;
      y = NULL;
      if((j=ReadXYP(str[4],MIE_MAXDATA,0,ci,u,1.0,NULL,&x,&y)) < 0)
      {
        return -1;
      }
      refi[n] = NULL;
      if((refi[n]=(double *)malloc(mie_n_wlen*sizeof(double))) == NULL)
      {
        fprintf(stderr,"%s: error in allocating memory.\n",fnam);
        return -1;
      }
      if(SamplingE(j,x,y,mie_n_wlen,mie_wlen,refi[n],1.0) < 0)
      {
        return -1;
      }
      free(x);
      free(y);
      // Comment
      //for(j=0; j<MAXLINE-1&&str[5][j]==' '; j++)
      //strncpy(cmnt[n],str[5]+j,MAXLINE);
      strncpy(line,"",MAXLINE);
      n++;
      i++;
    }
    fclose(fp);
    if(err)
    {
      return -1;
    }
    n_comp = n;
    fprintf(stderr,"%s: %d data have been read (%s)\n",fnam,n_comp,s);
    if(n_comp < 1)
    {
      fprintf(stderr,"%s: error, n_comp=%d (%s)\n",fnam,n_comp,s);
      return -1;
    }
  }

  return n_comp;
}

int GetOpt(int argn,char **args)
{
  int rt;

  fprintf(stderr,"%15s : $Revision: 788 $ $Date: 2014-09-01 01:36:09 +0900 (Mon, 01 Sep 2014) $\n","aeros_comp.c");

  if(CommonInit() < 0) return -1;
  if(OwnInit() < 0) return -1;
  rt = GetOwnOpt(argn,args);
  if(cnt_optn > 1)
  {
    if(GetCommonOpt(cnt_optn,cnt_opts) < 0) rt = -1;
  }
  if(ReadConfig() < 0) rt = -1;
  if(cnt_hp) rt = 1;

  return rt;
}

int GetCommonOpt(int argn,char **args)
{
  int c,rt;
  int ntmp;
  int option_index = 0;
  int this_option_optind;
  double xtmp;
  char *endp;
  struct option long_options[] =
  {
    {"inp_dnam",1,0,'D'},
    {"inp_ndig",1,0,'N'},
    {"inp_tdif",1,0,'t'},
    {"dsr_rerr",1,0,'e'},
    {"ssr_sfov",1,0,'F'},
    {"opt_vtau_init",1,0,'U'},
    {"opt_vmie_init",1,0,'V'},
    {"opt_wscl_init",1,0,'W'},
    {"opt_oscl_init",1,0,'O'},
    {"sim_mode",1,0,'s'},
    {"sim_iatm",1,0,'A'},
    {"sim_ihaz",1,0,'a'},
    {"sim_grho",1,0,'G'},
    {"sim_galt",1,0,'g'},
    {"sim_xsgm",1,0,'w'},
    {"sim_dbmp",1,0,'m'},
    {"sim_sbmp_init",1,0,'M'},
    {"mie_wlen_ref",1,0,'r'},
    {"mie_rmin",1,0,'x'},
    {"mie_rmax",1,0,'X'},
    {"mie_nstp",1,0,'n'},
    {"mie_wsgm",1,0,'S'},
    {"cnt_conf",0,0,'C'},
    {"cnt_cflg",0,0,'c'},
    {"cnt_sflg",0,0,'o'},
    {"cnt_db",0,0,'d'},
    {"cnt_vb",0,0,'v'},
    {"cnt_hp",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":D:N:t:e:F:U:V:W:O:s:A:a:G:g:w:m:M:r:x:X:n:S:C:codvh",long_options,&option_index);
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
      case 't':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) inp_tdif = ntmp;
        else
        {
          fprintf(stderr,"Time difference -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'e':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>=0.0) dsr_rerr = xtmp;
        else
        {
          fprintf(stderr,"DSR relative error -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'F':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) ssr_sfov = xtmp;
        else
        {
          fprintf(stderr,"SSR FOV area -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'U':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) opt_init[OPT_PAR_VTAU] = xtmp;
        else
        {
          fprintf(stderr,"Aerosol optical depth -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'V':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) opt_init[OPT_PAR_VMIE] = xtmp;
        else
        {
          fprintf(stderr,"Visibility -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'W':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>=0.0) opt_init[OPT_PAR_WSCL] = xtmp;
        else
        {
          fprintf(stderr,"Water scale -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'O':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>=0.0) opt_init[OPT_PAR_OSCL] = xtmp;
        else
        {
          fprintf(stderr,"Ozone scale -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 's':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=1 && ntmp<=2) sim_mode = ntmp;
        else
        {
          fprintf(stderr,"Simulation mode -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'A':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<=8) sim_iatm = ntmp;
        else
        {
          fprintf(stderr,"Atmosphere model# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'a':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<=10) sim_ihaz = ntmp;
        else
        {
          fprintf(stderr,"Aerosol model# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'G':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>=0.0)
        {
          sim_grho = xtmp;
          sim_salb = -1;
        }
        else
        {
          fprintf(stderr,"Ground albedo -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'g':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0') sim_galt = xtmp;
        else
        {
          fprintf(stderr,"Ground altitude -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'w':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) sim_xsgm = xtmp;
        else
        {
          fprintf(stderr,"X sigma -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'm':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && (ntmp==1||ntmp==5||ntmp==15)) sim_dbmp = ntmp;
        else
        {
          fprintf(stderr,"DSR Band model number -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'M':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && (ntmp==1||ntmp==5||ntmp==15)) sim_sbmp_init = ntmp;
        else
        {
          fprintf(stderr,"SSR Band model number -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'r':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) mie_wlen_ref = xtmp;
        else
        {
          fprintf(stderr,"Reference wavelength -> out of range %s\n",optarg);
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
      case 'n':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>0) mie_n_step = ntmp;
        else
        {
          fprintf(stderr,"Log10(R) #steps for mie -> out of range %s\n",optarg);
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
      case 'c':
        cnt_cflg++;
        break;
      case 'o':
        cnt_sflg = 1;
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

  if(cnt_hp) rt = 1;
  return rt;
}

int CommonUsage(char *name,int c)
{
  int n = 15;
  char e[MAXLINE];
  char a[MAXLINE];
  char d[MAXLINE];

  switch(c)
  {
    case 0:
      fprintf(stderr,"Usage:\n");
      fprintf(stderr,"%s -[option] (argument) -[option] (argument) ...\n",name);
      fprintf(stderr,"------------------------------------------------------------------------------------\n");
      fprintf(stderr,"   option         |%s|%s|%s| current\n",As(e,"",n),          As(a,"argument",n),   As(d,"default",n));
      break;
    case 'D':
      fprintf(stderr," D -inp_dnam      |%s|%s|%s| %s\n",As(e,"Data directory",n), As(a,"name",n),       As(d,INP_DNAM,n),inp_dnam);
      break;
    case 'N':
      fprintf(stderr," N -inp_ndig      |%s|%s|%s| %d\n",As(e,"RNO #Digits",n),    As(a,"#",n),          Ad(d,INP_NDIG,n),inp_ndig);
      break;
    case 't':
      fprintf(stderr," t -inp_tdif      |%s|%s|%s| %d\n",As(e,"Time diff",n),      As(a,"sec",n),        Ad(d,INP_TDIF,n),inp_tdif);
      break;
    case 'e':
      fprintf(stderr," e -dsr_rerr      |%s|%s|%s| %e\n",As(e,"DSR error",n),      As(a,"ratio",n),      Af(d,DSR_RERR,n),dsr_rerr);
      break;
    case 'F':
      fprintf(stderr," F -ssr_sfov      |%s|%s|%s| %e\n",As(e,"SSR FOV",n),        As(a,"sr",n),         Ae(d,SSR_SFOV,n),ssr_sfov);
      break;
    case 'U':
      fprintf(stderr," U -opt_vtau_init |%s|%s|%s| %e\n",As(e,"Optical depth",n),  As(a,"value",n),      Af(d,OPT_VTAU,n),opt_init[OPT_PAR_VTAU]);
      break;
    case 'V':
      fprintf(stderr," V -opt_vmie_init |%s|%s|%s| %e\n",As(e,"Visibility",n),     As(a,"km",n),         Af(d,OPT_VMIE,n),opt_init[OPT_PAR_VMIE]);
      break;
    case 'W':
      fprintf(stderr," W -opt_wscl_init |%s|%s|%s| %e\n",As(e,"Water scale",n),    As(a,"value",n),      Af(d,OPT_WSCL,n),opt_init[OPT_PAR_WSCL]);
      break;
    case 'O':
      fprintf(stderr," O -opt_oscl_init |%s|%s|%s| %e\n",As(e,"Ozone scale",n),    As(a,"value",n),      Af(d,OPT_OSCL,n),opt_init[OPT_PAR_OSCL]);
      break;
    case 's':
      fprintf(stderr," s -sim_mode      |%s|%s|%s| %d\n",As(e,"Simulation mode",n),As(a,"1|2",n),        Ad(d,SIM_MODTRAN_VMIE,n),sim_mode);
      break;
    case 'A':
      fprintf(stderr," A -sim_iatm      |%s|%s|%s| %d\n",As(e,"Atmosphere",n),     As(a,"#",n),          Ad(d,SIM_IATM,n),sim_iatm);
      break;
    case 'a':
      fprintf(stderr," a -sim_ihaz      |%s|%s|%s| %d\n",As(e,"Aerosol",n),        As(a,"#",n),          Ad(d,SIM_IHAZ,n),sim_ihaz);
      break;
    case 'G':
      fprintf(stderr," G -sim_grho      |%s|%s|%s| %f\n",As(e,"Ground albedo  ",n),As(a,"value",n),      Af(d,SIM_GRHO,n),sim_grho);
      break;
    case 'g':
      fprintf(stderr," g -sim_galt      |%s|%s|%s| %f\n",As(e,"Ground altitude",n),As(a,"km",n),         Af(d,SIM_GALT,n),sim_galt);
      break;
    case 'w':
      fprintf(stderr," w -sim_xsgm      |%s|%s|%s| %e\n",As(e,"X sigma",n),        As(a,"value",n),      Af(d,SIM_XSGM,n),sim_xsgm);
      break;
    case 'm':
      fprintf(stderr," m -sim_dbmp      |%s|%s|%s| %d\n",As(e,"DSR Band model",n), As(a,"1|5|15",n),     Ad(d,SIM_DBMP,n),sim_dbmp);
      break;
    case 'M':
      fprintf(stderr," M -sim_sbmp_init |%s|%s|%s| %d\n",As(e,"SSR Band model",n), As(a,"1|5|15",n),     Ad(d,SIM_SBMP,n),sim_sbmp_init);
      break;
    case 'r':
      fprintf(stderr," r -mie_wlen_ref  |%s|%s|%s| %f\n",As(e,"Reference wlen",n), As(a,"nm",n),         Af(d,MIE_WLEN_REF,n),mie_wlen_ref);
      break;
    case 'x':
      fprintf(stderr," x -mie_rmin      |%s|%s|%s| %e\n",As(e,"R Min",n),          As(a,"um",n),         Ae(d,MIE_RMIN,n),mie_rmin);
      break;
    case 'X':
      fprintf(stderr," X -mie_rmax      |%s|%s|%s| %e\n",As(e,"R Max",n),          As(a,"um",n),         Ae(d,MIE_RMAX,n),mie_rmax);
      break;
    case 'n':
      fprintf(stderr," n -mie_nstp      |%s|%s|%s| %d\n",As(e,"Log10(R) #steps",n),As(a,"#",n),          Ad(d,MIE_NSTP,n),mie_n_step);
      break;
    case 'S':
      fprintf(stderr," S -mie_wsgm      |%s|%s|%s| %f\n",As(e,"Log10(R) width",n), As(a,"sigma",n),      Af(d,MIE_WSGM,n),mie_wsgm);
      break;
    case 'C':
      fprintf(stderr," C -cnt_conf      |%s|%s|%s| %s\n",As(e,"Config  file",n),   As(a,"name",n),       As(d,CNT_CONF,n),cnt_conf);
      break;
    case 'c':
      fprintf(stderr," c -cnt_cflg      |%s|%s|%s| %d\n",As(e,"Clear   flag",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_cflg);
      break;
    case 'o':
      fprintf(stderr," o -cnt_sflg      |%s|%s|%s| %d\n",As(e,"Single  flag",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_sflg);
      break;
    case 'd':
      fprintf(stderr," d -cnt_db        |%s|%s|%s| %d\n",As(e,"Debug   mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_db);
      break;
    case 'v':
      fprintf(stderr," v -cnt_vb        |%s|%s|%s| %d\n",As(e,"Verbose mode",n),   As(a,"nothing",n),    Ad(d,0,n),cnt_vb);
      break;
    case 'h':
      fprintf(stderr," h -cnt_hp        |%s|%s|%s| %d\n",As(e,"Help    mode",n),   As(a,"nothing",n),    Ad(d,0,n),1);
      break;
    case 1:
      fprintf(stderr,"------------------------------------------------------------------------------------\n");
      fprintf(stderr,"Input data (stdin) format:\n");
      fprintf(stderr,"rno yr mo dy hr mi sc site_id inst_id cloud flag fov vazi valt par exposure temperature sazi salt airmass distance\n");
      fprintf(stderr,"site_id: 1:CEReS-1 2:CEReS-2 3:CEReS-3 4:CEReS-4 11:MRD-1 12:MRD-2 101:SUBARU 102:HALEPOHAKU\n");
      fprintf(stderr,"inst_id: 1:S06079.05 2:S07032.06\n");
      fprintf(stderr,"flag: 0-99:SSR 100-199:DSR 300-399:AUR(old) 500-599:AUR(new)\n");
      fprintf(stderr,"------------------------------------------------------------------------------------\n");
      fprintf(stderr,"Config file (cnt_conf) format:\n");
      fprintf(stderr,"mie_angl      name # u m M | file name(%s),angl column#(%d),deg or rad(%s),min angl(%5.1f),max angl(%6.1f)\n",
                                                   NONAME,0,"deg",MIE_ANGL_MIN,MIE_ANGL_MAX);
      fprintf(stderr,"mie_wlen      name # u m M | file name(%s),wlen column#(%d),unit in nm(%.1f),min wlen(%.1f),max wlen(%.1f)\n",
                                                   NONAME,0,1.0,MIE_WLEN_MIN,MIE_WLEN_MAX);
      fprintf(stderr,"dsr_wlen_0    name # unit  | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"dsr_wlen_1    name # unit  | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"dsr_wlen_2    name # unit  | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"dsr_wlen_3    name # unit  | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"dsr_fomg      name #       | file name(%s),romg column#(%d)\n",NONAME,1);
      fprintf(stderr,"aur_wlen_0    name # unit  | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"aur_wlen_1    name # unit  | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"aur_sbmp_0    name #       | file name(%s),sbmp column#(%d)\n",NONAME,0);
      fprintf(stderr,"aur_sbmp_1    name #       | file name(%s),sbmp column#(%d)\n",NONAME,0);
      fprintf(stderr,"aur_fomg      name #       | file name(%s),romg column#(%d)\n",NONAME,1);
      fprintf(stderr,"ssr_wlen_0    name # unit  | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"ssr_wlen_1    name # unit  | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"ssr_wlen_2    name # unit  | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"ssr_wlen_3    name # unit  | file name(%s),wlen column#(%d),unit in nm(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"ssr_sbmp_0    name #       | file name(%s),sbmp column#(%d)\n",NONAME,0);
      fprintf(stderr,"ssr_sbmp_1    name #       | file name(%s),sbmp column#(%d)\n",NONAME,0);
      fprintf(stderr,"ssr_sbmp_2    name #       | file name(%s),sbmp column#(%d)\n",NONAME,0);
      fprintf(stderr,"ssr_sbmp_3    name #       | file name(%s),sbmp column#(%d)\n",NONAME,0);
      fprintf(stderr,"ssr_fomg      name #       | file name(%s),romg column#(%d)\n",NONAME,1);
      fprintf(stderr,"ssr_th_los_min angle       | min th(los) in degree(%.1f)\n",SSR_TH_LOS_MIN);
      fprintf(stderr,"ssr_th_los_max angle       | max th(los) in degree(%.1f)\n",SSR_TH_LOS_MAX);
      fprintf(stderr,"ssr_sep_min    angle       | min separation in degree(%.1f)\n",SSR_SEP_MIN);
      fprintf(stderr,"ssr_sep_max    angle       | max separation in degree(%.1f)\n",SSR_SEP_MAX);
      fprintf(stderr,"crc_ffac      name         | file name(%s)\n",NONAME);
      fprintf(stderr,"crc_fsim      name         | file name(%s)\n",NONAME);
      fprintf(stderr,"crc_init      # value      | parameter#(%d),initial value(%f)\n",0,NAN);
      fprintf(stderr,"crc_fcmp    name # # u # # | file name(%s),real column#(%d),imag column#(%d),wlen unit in nm(%.1f),"
                                                   "min line#(%d),max line#(10^%.0f)\n",
                                                   NONAME,MIE_REAL_NUM,MIE_IMAG_NUM,1.0,MIE_IMIN,log10((double)MIE_IMAX));
      fprintf(stderr,"tst_fsim      name         | file name(%s)\n",NONAME);
      fprintf(stderr,"tst_fcmp    name # # u # # | file name(%s),real column#(%d),imag column#(%d),wlen unit in nm(%.1f),"
                                                   "min line#(%d),max line#(10^%.0f)\n",
                                                   NONAME,MIE_REAL_NUM,MIE_IMAG_NUM,1.0,MIE_IMIN,log10((double)MIE_IMAX));
      fprintf(stderr,"opt_init      # value      | parameter#(%d),initial value(%.5f)\n",0,OPT_VTAU);
      fprintf(stderr,"opt_fcmp    name # # u # # | file name(%s),real column#(%d),imag column#(%d),wlen unit in nm(%.1f),"
                                                   "min line#(%d),max line#(10^%.0f)\n",
                                                   NONAME,MIE_REAL_NUM,MIE_IMAG_NUM,1.0,MIE_IMIN,log10((double)MIE_IMAX));
      fprintf(stderr,"sim_path      name         | path name(%s)\n",SIM_PATH);
      fprintf(stderr,"sim_falb      name         | file name(%s)\n",SIM_FALB);
      fprintf(stderr,"sim_prf_vscl  # value      | mol#(%d),scale factor (%f)\n",0,NAN);
      fprintf(stderr,"sim_prf_pscl  # flag       | mol#(%d),scale flag(%d)\n",0,0);
      fprintf(stderr,"sim_prf_jcharx value       | jcharx(%d)\n",SIM_IATM);
      fprintf(stderr,"prf_falt      name # unit  | file name(%s),column#(%d),unit in km(%.1f)\n",NONAME,0,1.0);
      fprintf(stderr,"prf_fpre      name # unit  | file name(%s),column#(%d),unit(%c)\n",NONAME,1,'A');
      fprintf(stderr,"prf_ftmp      name # unit  | file name(%s),column#(%d),unit(%c)\n",NONAME,1,'A');
      fprintf(stderr,"prf_fmol    # name # unit  | mol#(%d),file name(%s),column#(%d),unit(%c)\n",0,NONAME,1,'A');
      fprintf(stderr,"prf_fhaz    # name #       | haz#(%d),file name(%s),column#(%d)\n",0,NONAME,1);
      fprintf(stderr,"prf_pres_gnd  value        | Ground pressure(%.1f)\n",NAN);
      fprintf(stderr,"prf_temp_gnd  value        | Ground temperature(%.1f)\n",NAN);
      fprintf(stderr,"cnt_tmx       # value      | comp#(%d),tmx flag(%d)\n",0,0);
      fprintf(stderr,"tmx_nang      #            | #angles(%d)\n",TMX_NANG);
      fprintf(stderr,"tmx_ndgs      #            | control #divisions(%d)\n",TMX_NDGS);
      fprintf(stderr,"tmx_kmax      #            | #quad points(%d)\n",TMX_KMAX);
      fprintf(stderr,"tmx_delt      value        | accuracy(%.8f)\n",TMX_DELT);
      fprintf(stderr,"tmx_shap      # #          | comp#(%d),shape#(%d)\n",0,TMX_SHAP);
      fprintf(stderr,"tmx_reps      # value      | comp#(%d),aspect ratio(%.8f)\n",0,TMX_REPS);
      for(n=0; n<cnt_n_own_format; n++)
      {
        fprintf(stderr,"%s",cnt_own_format[n]);
      }
      if(cnt_n_cmnt > 0)
      {
        fprintf(stderr,"<Read from %s>\n",cnt_conf);
        for(n=0; n<cnt_n_cmnt; n++)
        {
          fprintf(stderr,"%s",cnt_cmnt[n]);
        }
      }
      fprintf(stderr,"------------------------------------------------------------------------------------\n");
      fprintf(stderr,"Component file (crc_fcmp,tst_fcmp,opt_fcmp) format: wcom rmod lsgm refr refi comment\n");
      fprintf(stderr,"  wcom ... Weight of component [a.u.]\n");
      fprintf(stderr,"  rmod ... Mode radius [um]\n");
      fprintf(stderr,"  lsgm ... Sigma Log10(R [um])\n");
      fprintf(stderr,"  refr ... File name for refractive index (real part)\n");
      fprintf(stderr,"  refi ... File name for refractive index (imaginary part)\n");
      fprintf(stderr,"  comment ... Name of component etc.\n");
      fprintf(stderr,"Factor file (crc_ffac) format: dnum wlen dfac\n");
      fprintf(stderr,"  dnum ... Data number\n");
      fprintf(stderr,"  wlen ... Wavelength [nm]\n");
      fprintf(stderr,"  dfac ... Factor (DIS=t/DIS=f)\n");
      fprintf(stderr,"Profile file (prf_fpre,prf_ftmp,prf_fmol,prf_fhaz) format: z v\n");
      fprintf(stderr,"    z  ... Altitude [km]\n");
      fprintf(stderr,"    v  ... Value\n");
      fprintf(stderr,"Omega correction factor file (dsr_fomg,aur_fomg,ssr_fomg) format: th r\n");
      fprintf(stderr,"    th ... Incident angle [degree]\n");
      fprintf(stderr,"    r  ... Correction factor (relative sensitivity)\n");
      fprintf(stderr,"    Incident angle is low edge of each bin.\n");
      fprintf(stderr,"    The 1st/last incident angle must be the min/max incident angle.\n");
      fprintf(stderr,"    Angle interval must be a constant value.\n");
      fprintf(stderr,"    The last correction factor is supposed to be 0 (not used)\n");
      fprintf(stderr,"Output data format (debug):\n");
      fprintf(stderr,"  %s: w ea en o g\n",mie_out1);
      fprintf(stderr,"    w  ... Wavelength [nm]\n");
      fprintf(stderr,"    ea ... Extinction coefficient (absolute) [a.u.]\n");
      fprintf(stderr,"    en ... Extinction coefficient (normalized at reference wavelength)\n");
      fprintf(stderr,"    o  ... Single scattering albedo\n");
      fprintf(stderr,"    g  ... Asymmetry parameter\n");
      fprintf(stderr,"  %s: w a s\n",mie_out2);
      fprintf(stderr,"    w  ... Wavelength [nm]\n");
      fprintf(stderr,"    a  ... Scattering angle [degree]\n");
      fprintf(stderr,"    s  ... Phase function of averaged wave\n");
      break;
    default:
      return -1;
      break;
  }

  return 0;
}
