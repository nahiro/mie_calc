#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <bits/nan.h>
#include "strutil.h"

#define	R_TO_D			57.29577951308232		// 180/PI rad -> deg
#define	D_TO_R			1.745329251994329e-02		// PI/180 deg -> rad
#define	ANUM			0				// Out. angle column#
#define	WNUM			0				// Out. wave  column#
#define	PNUM			0				// Phase fun. (real) column#
#define	QNUM			0				// Phase fun. (imag) column#
#define	IMIN			0				// Min line#
#define	IMAX			1000000000			// Min line#
#define	AMIN			0.0				// Min angle in degree
#define	AMAX			180.0				// Max angle in degree
#define	WMIN			250.0				// Min wavelength in nm
#define	WMAX			1100.0				// Max wavelength in nm
#define	WREF			550.0				// Reference wavelength in nm
#define	WUNI			1.0				// Wavelength unit in nm
#define	RMIN			1.0e-3				// R min in um
#define	RMAX			1.0e+1				// R max in um
#define	LSTP			1.0e-3				// Log10(R) step
#define MAXENTR			20				// Max #entries in a line
#define MAXLINE			256				// Max #chars in a line
#define MAXDATA			150000				// Max #data
#define MAXCOMP			10				// Max #components
#define MAXSTEP			799				// Max #steps
#define EPSILON			1.0e-14				// A small number
#define FNAM			"out1.dat"			// Out. file 1
#define GNAM			"out2.dat"			// Out. file 2
#define ANAM			"ANG.dat"			// Out. angle file
#define WNAM			"WAVE.dat"			// Out. wave  file
#define	MAX(a,b)		((b)>(a)?(b):(a))

int	anum			= ANUM;				// Out. angle column#
int	wnum			= WNUM;				// Out. wave  column#
int	pnum			= PNUM;				// Phase fun. (real) column#
int	qnum			= QNUM;				// Phase fun. (imag) column#
int	imin			= IMIN;				// Min line#
int	imax			= IMAX;				// Max line#
int	db			= 0;				// Debug   mode
int	vb			= 0;				// Verbose mode
int	hp			= 0;				// Help    mode
int	n_angl_out		= 0;				// #Angles
int	n_wave_out		= 0;				// #Wavelengths (output)
int	n_com			= 0;				// #Components
int	n_stp			= 0;				// #R steps
int	iref			= -1;
int	i_wave_min		= IMAX;				// Min wave line#
int	i_wave_max		= IMIN;				// Max wave line#
char	fnam[MAXLINE]		= FNAM;				// Out. file 1
char	gnam[MAXLINE]		= GNAM;				// Out. file 2
char	anam[MAXLINE]		= ANAM;				// Out. angle file
char	wnam[MAXLINE]		= WNAM;				// Out. wave  file
double	amin			= AMIN;				// Min angle in degree
double	amax			= AMAX;				// Max angle in degree
double	wmin			= WMIN;				// Min wavelength in nm
double	wmax			= WMAX;				// Max wavelength in nm
double	wref			= WREF;				// Reference wavelength
double	wuni			= WUNI;				// Wavelength unit
double	rmin			= RMIN;				// R min in um
double	rmax			= RMAX;				// R max in um
double	lstp			= LSTP;				// Log10(R) step
double	lmin;
double	lmax;
double	*angl_out;
double	*wave_out;
double	*aext_out;
double	*asca_out;
double	*asym_out;
double	*phs1_out;
double	*phs2_out;
double	*phs3_out;
double	*refr_com[MAXCOMP];					// Refractive index (real)
double	*refi_com[MAXCOMP];					// Refractive index (imaginary)
double	*aext_com[MAXCOMP];					// Extinction coefficient
double	*asca_com[MAXCOMP];					// Scattering coefficient
double	*asym_com[MAXCOMP];					// Asymmetry parameter
double	*phs1_com[MAXCOMP];					// Phase function 1
double	*phs2_com[MAXCOMP];					// Phase function 2
double	*phs3_com[MAXCOMP];					// Phase function 3
double	wcom_com[MAXCOMP];					// Weight
double	lmod_com[MAXCOMP];					// Log10(Mode radius in um)
double	lsgm_com[MAXCOMP];					// Log10(sigma in um)
char	fnam_refr_com[MAXCOMP][MAXLINE];			// Refractive index (real)
char	fnam_refi_com[MAXCOMP][MAXLINE];			// Refractive index (imaginary)
char	cmnt_com[MAXCOMP][MAXLINE];				// Comment

int Init(void);
void Finish(void);
int MieCalc(void);
int MixComp(void);
int Printout(void);
int ReadColumn(char *s,int c,double u,int (*func)(int a,double b),double **d,int size);
int AnglSelect(int a,double b);
int WaveSelect(int a,double b);
int PfunSelect(int a,double b);
int GetOpt(int argn,char **args);
int Usage(void);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  if(MieCalc() < 0) exit(-1);
  if(MixComp() < 0) exit(-1);

  Printout();
  Finish();

  return 0;
}

int Init(void)
{
  int i,j,n;
  int nc;
  double v[3];
  char line[MAXLINE];
  char str[6][MAXLINE];
  char *p;

  if((n_angl_out=ReadColumn(anam,anum,1.00,AnglSelect,&angl_out,MAXDATA)) < 0) return -1;
  if((n_wave_out=ReadColumn(wnam,wnum,wuni,WaveSelect,&wave_out,MAXDATA)) < 0) return -1;
  for(i=0; i<n_wave_out; i++)
  {
    if(fabs(wave_out[i]-wref) < EPSILON)
    {
      iref = i;
      break;
    }
  }
  if(iref < 0)
  {
    fprintf(stderr,"Failed in finding reference wavelength %13.4e\n",wref);
    return -1;
  }

  i = n = 0;
  while(fgets(line,MAXLINE,stdin) != NULL)
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
      fprintf(stderr,"Read error >>> %s\n",line);
      return -1;
    }
    strcpy(str[5],line+nc);
    errno = 0;
    for(j=0; j<3; j++)
    {
      v[j] = strtod(str[j],&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"Convert error >>> %s\n",line);
        return -1;
      }
    }
    if(n >= MAXCOMP)
    {
      fprintf(stderr,"Warning, #data exceed the limit %d\n",n);
      break;
    }
    wcom_com[n] = v[0];
    lmod_com[n] = log10(v[1]);
    lsgm_com[n] = v[2];
    strncpy(fnam_refr_com[n],str[3],MAXLINE);
    strncpy(fnam_refi_com[n],str[4],MAXLINE);
    for(j=0; j<MAXLINE-1&&str[5][j]==' '; j++)
    strncpy(cmnt_com[n],str[5]+j,MAXLINE);
    strncpy(line,"",MAXLINE);
    n++;
    i++;
  }
  n_com = n;
  fprintf(stderr,"%d data have been read.\n",n_com);
  if(n_com < 1)
  {
    fprintf(stderr,"Error, n_com=%d\n",n_com);
    return -1;
  }

  for(n=0; n<n_com; n++)
  {
    if((i=ReadColumn(fnam_refr_com[n],pnum,1.0,PfunSelect,&refr_com[n],MAXDATA)) != n_wave_out)
    {
      fprintf(stderr,"Error, i=%d, n_wave_out=%d\n",i,n_wave_out);
      return -1;
    }
    if((i=ReadColumn(fnam_refi_com[n],qnum,1.0,PfunSelect,&refi_com[n],MAXDATA)) != n_wave_out)
    {
      fprintf(stderr,"Error, i=%d, n_wave_out=%d\n",i,n_wave_out);
      return -1;
    }
  }

  for(n=0; n<n_com; n++)
  {
    aext_com[n] = (double *)malloc(n_wave_out*sizeof(double));
    asca_com[n] = (double *)malloc(n_wave_out*sizeof(double));
    asym_com[n] = (double *)malloc(n_wave_out*sizeof(double));
    phs1_com[n] = (double *)malloc(n_wave_out*n_angl_out*sizeof(double));
    phs2_com[n] = (double *)malloc(n_wave_out*n_angl_out*sizeof(double));
    phs3_com[n] = (double *)malloc(n_wave_out*n_angl_out*sizeof(double));
    if(aext_com[n]==NULL || asca_com[n]==NULL || asym_com[n]==NULL ||
       phs1_com[n]==NULL || phs2_com[n]==NULL || phs3_com[n]==NULL)
    {
      fprintf(stderr,"Failed in allocating memory\n");
      return -1;
    }
  }
  aext_out = (double *)malloc(n_wave_out*sizeof(double));
  asca_out = (double *)malloc(n_wave_out*sizeof(double));
  asym_out = (double *)malloc(n_wave_out*sizeof(double));
  phs1_out = (double *)malloc(n_wave_out*n_angl_out*sizeof(double));
  phs2_out = (double *)malloc(n_wave_out*n_angl_out*sizeof(double));
  phs3_out = (double *)malloc(n_wave_out*n_angl_out*sizeof(double));
  if(aext_out==NULL || asca_out==NULL || asym_out==NULL ||
     phs1_out==NULL || phs2_out==NULL || phs3_out==NULL)
  {
    fprintf(stderr,"Failed in allocating memory\n");
    return -1;
  }

  lmin = log10(rmin);
  lmax = log10(rmax);
  if((n_stp=(lmax-lmin)/lstp) > MAXSTEP)
  {
    fprintf(stderr,"Error, n_stp=%d, MAXSTEP=%d\n",n_stp,MAXSTEP);
    return -1;
  }

  return 0;
}

void Finish(void)
{
  int n;

  for(n=0; n<n_com; n++)
  {
    free(refr_com[n]);
    free(refi_com[n]);
    free(aext_com[n]);
    free(asca_com[n]);
    free(asym_com[n]);
    free(phs1_com[n]);
    free(phs2_com[n]);
    free(phs3_com[n]);
  }
  free(angl_out);
  free(wave_out);
  free(aext_out);
  free(asca_out);
  free(asym_out);
  free(phs1_out);
  free(phs2_out);
  free(phs3_out);

  return;
}

int MieCalc(void)
{
  int i,j,k,n;

  for(n=0; n<n_com; n++)
  {
    for(i=0; i<n_wave_out; i++)
    {
      aext_com[n][i] = 0.0;
      asca_com[n][i] = 0.0;
      asym_com[n][i] = 0.0;
      for(j=0; j<n_angl_out; j++)
      {
        k = n_angl_out*i+j;
        phs1_com[n][k] = 0.0;
        phs2_com[n][k] = 0.0;
        phs3_com[n][k] = 0.0;
      }
    }
  }
  for(n=0; n<n_com; n++)
  {
    for(i=0; i<n_wave_out; i++)
    {
      for(j=0; j<n_angl_out; j++)
      {
      }
    }
  }

  return 0;
}

int MixComp(void)
{
  int i,j,k,n;
  double ws;
  double sw,se,ss,sa;
  double s1,s2,s3;

  for(i=0; i<n_wave_out; i++)
  {
    se = ss = sa = 0.0;
    for(n=0; n<n_com; n++)
    {
      sw += wcom_com[n];
      se += wcom_com[n]*aext_com[n][i];
      ws  = wcom_com[n]*asca_com[n][i];
      ss += ws;
      sa += ws*asym_com[n][i];
    }
    aext_out[i] = se/sw;
    asca_out[i] = ss/sw;
    asym_out[i] = sa/ss;
    for(j=0; j<n_angl_out; j++)
    {
      k = n_angl_out*i+j;
      s1 = s2 = s3 = 0.0;
      for(n=0; n<n_com; n++)
      {
        ws  = wcom_com[n]*asca_com[n][i];
        s1 += ws*phs1_com[n][k];
        s2 += ws*phs2_com[n][k];
        s3 += ws*phs3_com[n][k];
      }
      phs1_out[k] = s1/ss;
      phs2_out[k] = s2/ss;
      phs3_out[k] = s3/ss;
    }
  }

  return 0;
}

int Printout(void)
{
  int i,j,k;
  FILE *fp;

  if((fp=fopen(fnam,"w")) == NULL)
  {
    fprintf(stderr,"Printout: cannot open %s\n",fnam);
    return -1;
  }
  for(i=0; i<n_wave_out; i++)
  {
    fprintf(fp,"%7.2f %13.4e %13.4e %13.4e %13.4e\n",wave_out[i],aext_out[i],aext_out[i]/aext_out[iref],asca_out[i]/aext_out[i],asym_out[i]);
  }
  fclose(fp);

  if((fp=fopen(gnam,"w")) == NULL)
  {
    fprintf(stderr,"Printout: cannot open %s\n",fnam);
    return -1;
  }
  for(i=0; i<n_wave_out; i++)
  {
    for(j=0; j<n_angl_out; j++)
    {
      k = n_angl_out*i+j;
      fprintf(fp,"%7.2f %7.2f %13.4e %13.4e %13.4e\n",wave_out[i],angl_out[j],phs1_out[k],phs2_out[k],phs3_out[k]);
    }
  }
  fclose(fp);

  return 0;
}

int ReadColumn(char *s,int c,double u,int (*func)(int a,double b),double **d,int size)
{
  int i,j,n;
  int nc,nd;
  int err;
  int flag;
  char line[MAXLINE];
  char str[MAXENTR][MAXLINE];
  char *p;
  double vtmp;
  double *v;
  FILE *fp;

  // check input
  if(c<0 || c>=MAXENTR)
  {
    fprintf(stderr,"ReadColumn: invalid column number >>> %d (%s)\n",c,s);
    return -1;
  }
  // allocate memory
  if((v=(double*)malloc(size*sizeof(double))) == NULL)
  {
    fprintf(stderr,"ReadColumn: failed in allocating memory (%s)\n",s);
    return -1;
  }
  // open input
  if(strcmp(s,"") == 0)
  {
    flag = 1;
    fp = stdin;
  }
  else
  {
    flag = 0;
    if((fp=fopen(s,"r")) == NULL)
    {
      fprintf(stderr,"ReadColumn: cannot open %s\n",s);
      free(v);
      return -1;
    }
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
        fprintf(stderr,"ReadColumn: read error >>> %s (%s)\n",line,s);
        err = 1;
        break;
      }
    }
    if(j != c+1)
    {
      break;
    } else
    errno = 0;
    vtmp = strtod(str[c],&p)*u;
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"ReadColumn: convert error >>> %s (%s)\n",line,s);
      err = 1;
      break;
    }
    if(func(i,vtmp))
    {
      i++;
      continue;
    }
    if(n >= size)
    {
      fprintf(stderr,"Warning, #data exceed the limit %d (%s)\n",n,s);
      break;
    }
    v[n] = vtmp;
    n++;
    i++;
  }
  if(flag)
  {
    fclose(fp);
  }
  if(err)
  {
    free(v);
    return -1;
  }
  nd = n;
  if(nd < 1)
  {
    fprintf(stderr,"ReadColumn: error, nd=%d (%s)\n",nd,s);
    free(v);
    return -1;
  }

  // resize memory
  if((*d=(double*)realloc(v,nd*sizeof(double))) == NULL)
  {
    fprintf(stderr,"ReadColumn: failed in allocating memory (%s)\n",s);
    free(v);
    return -1;
  }

  return nd;
}

int AnglSelect(int a,double b)
{

  if(b<amin || b>amax) return 1;

  return 0;
}

int WaveSelect(int a,double b)
{
  if(b<wmin || b>wmax) return 1;
  if(a < i_wave_min) i_wave_min = a;
  if(a > i_wave_max) i_wave_max = a;

  return 0;
}

int PfunSelect(int a,double b)
{
  if(a<i_wave_min || a>i_wave_max) return 1;

  return 0;
}

int GetOpt(int argn,char **args)
{
  int c,rt;
  int option_index = 0;
  int this_option_optind;
  int ntmp;
  char *endp;
  double xtmp;
  struct option long_options[] =
  {
    {"fnam",1,0,'f'},
    {"gnam",1,0,'g'},
    {"anam",1,0,'a'},
    {"wnam",1,0,'w'},
    {"anum",1,0,'A'},
    {"wnum",1,0,'W'},
    {"pnum",1,0,'P'},
    {"qnum",1,0,'Q'},
    {"imin",1,0,'i'},
    {"imax",1,0,'I'},
    {"wref",1,0,'r'},
    {"wuni",1,0,'u'},
    {"rmin",1,0,'x'},
    {"rmax",1,0,'X'},
    {"lstp",1,0,'s'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":f:g:a:w:A:W:P:Q:i:I:r:u:x:X:s:dvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'f':
        strncpy(fnam,optarg,MAXLINE);
        break;
      case 'g':
        strncpy(gnam,optarg,MAXLINE);
        break;
      case 'a':
        strncpy(anam,optarg,MAXLINE);
        break;
      case 'w':
        strncpy(wnam,optarg,MAXLINE);
        break;
      case 'A':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<MAXENTR) anum = ntmp;
        else
        {
          fprintf(stderr,"Out. angle column# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'W':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<MAXENTR) wnum = ntmp;
        else
        {
          fprintf(stderr,"Out. wave column# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'P':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<MAXENTR) pnum = ntmp;
        else
        {
          fprintf(stderr,"Phase function column# (real) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'Q':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0 && ntmp<MAXENTR) qnum = ntmp;
        else
        {
          fprintf(stderr,"Phase function column# (imag) -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'i':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0) imin = ntmp;
        else
        {
          fprintf(stderr,"Minimum line# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'I':
        errno = 0;
        ntmp = strtol(optarg,&endp,10);
        if(errno!=ERANGE && *endp=='\0' && ntmp>=0) imax = ntmp;
        else
        {
          fprintf(stderr,"Maximum line# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'r':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) wref = xtmp;
        else
        {
          fprintf(stderr,"Reference wavelength -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'u':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) wuni = xtmp;
        else
        {
          fprintf(stderr,"Wavelength unit -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'x':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) rmin = xtmp;
        else
        {
          fprintf(stderr,"R min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'X':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) rmax = xtmp;
        else
        {
          fprintf(stderr,"R max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 's':
        errno = 0;
        xtmp = strtod(optarg,&endp);
        if(errno!=ERANGE && *endp=='\0' && xtmp>0.0) lstp = xtmp;
        else
        {
          fprintf(stderr,"Log10(R) step -> out of range %s\n",optarg);
          rt = -1;
        }
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
  int  n = 16;
  char e[MAXLINE];
  char a[MAXLINE];
  char d[MAXLINE];

  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"common -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," f -fnam    |%s|%s|%s| %s\n",As(e,"Out. file 1",n),   As(a,"name",n),       As(d,FNAM,n),fnam);
  fprintf(stderr," g -gnam    |%s|%s|%s| %s\n",As(e,"Out. file 2",n),   As(a,"name",n),       As(d,GNAM,n),gnam);
  fprintf(stderr," a -anam    |%s|%s|%s| %s\n",As(e,"Out. ang. file",n),As(a,"name",n),       As(d,ANAM,n),anam);
  fprintf(stderr," w -wnam    |%s|%s|%s| %s\n",As(e,"Out. wave file",n),As(a,"name",n),       As(d,WNAM,n),wnam);
  fprintf(stderr," A -anum    |%s|%s|%s| %d\n",As(e,"Out. ang. col#",n),As(a,"#",n),          Ad(d,ANUM,n),anum);
  fprintf(stderr," W -wnum    |%s|%s|%s| %d\n",As(e,"Out. wave col#",n),As(a,"#",n),          Ad(d,WNUM,n),wnum);
  fprintf(stderr," P -pnum    |%s|%s|%s| %d\n",As(e,"Pfun. (r) col#",n),As(a,"#",n),          Ad(d,PNUM,n),pnum);
  fprintf(stderr," Q -qnum    |%s|%s|%s| %d\n",As(e,"Pfun. (i) col#",n),As(a,"#",n),          Ad(d,QNUM,n),qnum);
  fprintf(stderr," i -imin    |%s|%s|%s| %d\n",As(e,"Minimum line#",n), As(a,"#",n),          Ad(d,IMIN,n),imin);
  fprintf(stderr," I -imax    |%s|%s|%s| %d\n",As(e,"Maximum line#",n), As(a,"#",n),          Ad(d,IMAX,n),imax);
  fprintf(stderr," r -wref    |%s|%s|%s| %f\n",As(e,"Reference wave",n),As(a,"nm",n),         Af(d,WREF,n),wref);
  fprintf(stderr," u -wuni    |%s|%s|%s| %e\n",As(e,"Wave unit",n),     As(a,"nm",n),         Af(d,WUNI,n),wuni);
  fprintf(stderr," x -rmin    |%s|%s|%s| %e\n",As(e,"R Min",n),         As(a,"um",n),         Ae(d,RMIN,n),rmin);
  fprintf(stderr," X -rmax    |%s|%s|%s| %e\n",As(e,"R Max",n),         As(a,"um",n),         Ae(d,RMAX,n),rmax);
  fprintf(stderr," s -lstp    |%s|%s|%s| %e\n",As(e,"Log10(R) step",n), As(a,"value",n),      Ae(d,LSTP,n),lstp);
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"Input data format (stdin): wcom rmod rsgm refr refi comment\n");
  fprintf(stderr,"  wcom ... Weight of component [a.u.]\n");
  fprintf(stderr,"  rmod ... Mode radius [um]\n");
  fprintf(stderr,"  rsgm ... Sigma Log10(R [um])\n");
  fprintf(stderr,"  refr ... File name for refractive index (real part)\n");
  fprintf(stderr,"  refi ... File name for refractive index (imaginary part)\n");
  fprintf(stderr,"  comment ... Name of component etc.\n");
  fprintf(stderr,"Output data format:\n");
  fprintf(stderr,"  Out. file 1: w ea en o g\n");
  fprintf(stderr,"    w  ... Wavelength [nm]\n");
  fprintf(stderr,"    ea ... Extinction coefficient (absolute) [a.u.]\n");
  fprintf(stderr,"    en ... Extinction coefficient (normalized at reference wavelength)\n");
  fprintf(stderr,"    o  ... Single scattering albedo\n");
  fprintf(stderr,"    g  ... Asymmetry parameter\n");
  fprintf(stderr,"  Out. file 2: w a s1 s2 s3\n");
  fprintf(stderr,"    w  ... Wavelength [nm]\n");
  fprintf(stderr,"    a  ... Scattering angle [degree]\n");
  fprintf(stderr,"    s1 ... Phase function of s-polarized wave\n");
  fprintf(stderr,"    s2 ... Phase function of p-polarized wave\n");
  fprintf(stderr,"    s3 ... Phase function of averaged wave\n");

  return 0;
}
