#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>

#define	PRF_NALT		50
#define	PRF_NMOD		6
// 0 Tropical Atmosphere (15° North Latitude)
// 1 Mid-Latitude Summer (45° North Latitude)
// 2 Mid-Latitude Winter (45° North Latitude)
// 3 Sub-Arctic Summer (60° North Latitude)
// 4 Sub-Arctic Winter (60° North Latitude)
// 5 1976 US Standard Atmosphere
#define	PRF_NMOL		8
//  0:H2O    1:CO2    2:OZONE  3:N2O
//  4:CO     5:CH4    6:O2     7:DENSITY
#define	PRF_NTRC		21
//  0:NO     1:SO2    2:NO2    3:NH3
//  4:HNO3   5:OH     6:HF     7:HCL
//  8:HBR    9:HI    10:CLO   11:OCS
// 12:H2CO  13:HOCL  14:N2    15:HCN
// 16:CH3CL 17:H2O2  18:C2H2  19:C2H6
// 20:PH3
#define	MAXLINE			256

int Init(void);
int ReadAlt(void);
int ReadPres(void);
int ReadTemp(void);
int ReadAmol(void);
int ReadTrac(void);

double	prf_alt[PRF_NALT];
double	prf_pres[PRF_NMOD][PRF_NALT];
double	prf_temp[PRF_NMOD][PRF_NALT];
double	prf_amol[PRF_NMOL][PRF_NMOD][PRF_NALT];
double	prf_trac[PRF_NTRC][PRF_NALT];
char	prf_ddir[MAXLINE]	= "/home/naohiro/work/MODTRAN4/mie_comp/data/profile";
char	prf_mod_name[PRF_NMOD][MAXLINE] =
{
  "Tropical Atmosphere",
  "Mid-Latitude Summer",
  "Mid-Latitude Winter",
  "Sub-Arctic Summer",
  "Sub-Arctic Winter",
  "1976 US Standard Atmosphere",
};
char	prf_amol_name[PRF_NMOL][MAXLINE] =
{
  "H2O","CO2","OZONE","N2O",
  "CO","CH4","O2","DENSITY",
};
char	prf_trac_name[PRF_NTRC][MAXLINE] =
{
  "NO","SO2","NO2","NH3",
  "HNO3","OH","HF","HCL",
  "HBR","HI","CLO","OCS",
  "H2CO","HOCL","N2","HCN",
  "CH3CL","H2O2","C2H2","C2H6",
  "PH3",
};

int main(void)
{
  int i,j,n;

  if(Init() < 0) return -1;

  // Defines
  printf("#define\tPRF_NALT\t\t50\n");
  printf("#define\tPRF_NMOD\t\t6\n");
  printf("#define\tPRF_NMOL\t\t8\n");
  printf("#define\tPRF_NTRC\t\t21\n");
  printf("\n");
  // Altitude
  printf("double prf_alt[PRF_NALT] =\n");
  printf("{\n");
  for(n=0; n<PRF_NALT; n++)
  {
    printf("%s%5.1f,%s",(n%8==0?"  ":""),prf_alt[n],(n==PRF_NALT-1?"\n":(n%8==7?"\n":" ")));
  }
  printf("};\n");
  printf("\n");
  // Pressure
  printf("double prf_pres[PRF_NMOD][PRF_NALT] =\n");
  printf("{\n");
  for(i=0; i<PRF_NMOD; i++)
  {
    printf("  // %d. %s\n",i+1,prf_mod_name[i]);
    printf("  {\n");
    for(n=0; n<PRF_NALT; n++)
    {
      printf("%s%14.8e,%s",(n%8==0?"    ":""),exp(prf_pres[i][n]),(n==PRF_NALT-1?"\n":(n%8==7?"\n":"")));
    }
    printf("  },\n");
  }
  printf("};\n");
  printf("\n");
  // Temperature
  printf("double prf_temp[PRF_NMOD][PRF_NALT] =\n");
  printf("{\n");
  for(i=0; i<PRF_NMOD; i++)
  {
    printf("  // %d. %s\n",i+1,prf_mod_name[i]);
    printf("  {\n");
    for(n=0; n<PRF_NALT; n++)
    {
      printf("%s%7.3f,%s",(n%8==0?"    ":""),prf_temp[i][n],(n==PRF_NALT-1?"\n":(n%8==7?"\n":" ")));
    }
    printf("  },\n");
  }
  printf("};\n");
  printf("\n");
  // Molecules
  printf("double prf_amol[PRF_NMOL][PRF_NMOD][PRF_NALT] =\n");
  printf("{\n");
  for(j=0; j<PRF_NMOL; j++)
  {
    printf("  // %2d. %s\n",j+1,prf_amol_name[j]);
    printf("  {\n");
    for(i=0; i<PRF_NMOD; i++)
    {
      printf("    // %d. %s\n",i+1,prf_mod_name[i]);
      printf("    {\n");
      for(n=0; n<PRF_NALT; n++)
      {
        printf("%s%10.3e,%s",(n%8==0?"      ":""),prf_amol[j][i][n],(n==PRF_NALT-1?"\n":(n%8==7?"\n":" ")));
      }
      printf("    },\n");
    }
    printf("  },\n");
  }
  printf("};\n");
  printf("\n");
  // Trace gases
  printf("double prf_trac[PRF_NTRC][PRF_NALT] =\n");
  printf("{\n");
  for(i=0; i<PRF_NTRC; i++)
  {
    printf("  // %2d. %s\n",i+1,prf_trac_name[i]);
    printf("  {\n");
    for(n=0; n<PRF_NALT; n++)
    {
      printf("%s%9.2e,%s",(n%8==0?"    ":""),prf_trac[i][n],(n==PRF_NALT-1?"\n":(n%8==7?"\n":" ")));
    }
    printf("  },\n");
  }
  printf("};\n");

  return 0;
}

int Init(void)
{
  if(ReadAlt()  < 0) return -1;
  if(ReadPres() < 0) return -1;
  if(ReadTemp() < 0) return -1;
  if(ReadAmol() < 0) return -1;
  if(ReadTrac() < 0) return -1;
  return 0;
}

int ReadAlt(void)
{
  int n;
  int err;
  double v;
  char line[MAXLINE];
  char str1[MAXLINE];
  char fnam[MAXLINE];
  char *p;
  FILE *fp;

  sprintf(fnam,"%s/ALT.dat",prf_ddir);
  if((fp=fopen(fnam,"r")) == NULL)
  {
    fprintf(stderr,"ReadAlt: cannot open %s\n",fnam);
    return -1;
  }
  n = 0;
  err = 0;
  while(fgets(line,MAXLINE,fp) != NULL)
  {
    if(sscanf(line,"%s",str1) != 1)
    {
      fprintf(stderr,"ReadAlt: read error >>> %s\n",line);
      err = 1;
      break;
    }
    errno = 0;
    v = strtod(str1,&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"ReadAlt: convert error >>> %s\n",line);
      err = 1;
      break;
    }
    if(n >= PRF_NALT)
    {
      fprintf(stderr,"ReadAlt: error, #data exceed the limit >>> %d\n",n);
      err = 1;
      break;
    }
    prf_alt[n] = v;
    n++;
  }
  fclose(fp);
  if(n != PRF_NALT)
  {
    fprintf(stderr,"ReadAlt: error, n=%d, expected=%d\n",n,PRF_NALT);
    err = 1;
  }
  if(err)
  {
    return -1;
  }

  return 0;
}

int ReadPres(void)
{
  int i,n;
  int err;
  double v;
  char line[MAXLINE];
  char str1[MAXLINE];
  char fnam[MAXLINE];
  char *p;
  FILE *fp;

  for(i=0; i<PRF_NMOD; i++)
  {
    sprintf(fnam,"%s/PMLOG_IALT_%d.dat",prf_ddir,i+1);
    if((fp=fopen(fnam,"r")) == NULL)
    {
      fprintf(stderr,"ReadPres: cannot open %s\n",fnam);
      return -1;
    }
    n = 0;
    err = 0;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      if(sscanf(line,"%s",str1) != 1)
      {
        fprintf(stderr,"ReadPres: read error >>> %s\n",line);
        err = 1;
        break;
      }
      errno = 0;
      v = strtod(str1,&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"ReadPres: convert error >>> %s\n",line);
        err = 1;
        break;
      }
      if(n >= PRF_NALT)
      {
        fprintf(stderr,"ReadPres: error, #data exceed the limit >>> %d\n",n);
        err = 1;
        break;
      }
      prf_pres[i][n] = v;
      n++;
    }
    fclose(fp);
    if(n != PRF_NALT)
    {
      fprintf(stderr,"ReadPres: error, n=%d, expected=%d\n",n,PRF_NALT);
      err = 1;
    }
    if(err)
    {
      return -1;
    }
  }

  return 0;
}

int ReadTemp(void)
{
  int i,n;
  int err;
  double v;
  char line[MAXLINE];
  char str1[MAXLINE];
  char fnam[MAXLINE];
  char *p;
  FILE *fp;

  for(i=0; i<PRF_NMOD; i++)
  {
    sprintf(fnam,"%s/TMATM_IALT_%d.dat",prf_ddir,i+1);
    if((fp=fopen(fnam,"r")) == NULL)
    {
      fprintf(stderr,"ReadTemp: cannot open %s\n",fnam);
      return -1;
    }
    n = 0;
    err = 0;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      if(sscanf(line,"%s",str1) != 1)
      {
        fprintf(stderr,"ReadTemp: read error >>> %s\n",line);
        err = 1;
        break;
      }
      errno = 0;
      v = strtod(str1,&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"ReadTemp: convert error >>> %s\n",line);
        err = 1;
        break;
      }
      if(n >= PRF_NALT)
      {
        fprintf(stderr,"ReadTemp: error, #data exceed the limit >>> %d\n",n);
        err = 1;
        break;
      }
      prf_temp[i][n] = v;
      n++;
    }
    fclose(fp);
    if(n != PRF_NALT)
    {
      fprintf(stderr,"ReadTemp: error, n=%d, expected=%d\n",n,PRF_NALT);
      err = 1;
    }
    if(err)
    {
      return -1;
    }
  }

  return 0;
}

int ReadAmol(void)
{
  int i,j,n;
  int err;
  double v;
  char line[MAXLINE];
  char str1[MAXLINE];
  char fnam[MAXLINE];
  char *p;
  FILE *fp;

  for(i=0; i<PRF_NMOD; i++)
  {
    for(j=0; j<PRF_NMOL; j++)
    {
      sprintf(fnam,"%s/AMOL_IALT_%d_%d.dat",prf_ddir,i+1,j+1);
      if((fp=fopen(fnam,"r")) == NULL)
      {
        fprintf(stderr,"ReadAmol: cannot open %s\n",fnam);
        return -1;
      }
      n = 0;
      err = 0;
      while(fgets(line,MAXLINE,fp) != NULL)
      {
        if(sscanf(line,"%s",str1) != 1)
        {
          fprintf(stderr,"ReadAmol: read error >>> %s\n",line);
          err = 1;
          break;
        }
        errno = 0;
        v = strtod(str1,&p);
        if(errno==ERANGE || *p!='\0')
        {
          fprintf(stderr,"ReadAmol: convert error >>> %s\n",line);
          err = 1;
          break;
        }
        if(n >= PRF_NALT)
        {
          fprintf(stderr,"ReadAmol: error, #data exceed the limit >>> %d\n",n);
          err = 1;
          break;
        }
        prf_amol[j][i][n] = v;
        n++;
      }
      fclose(fp);
      if(n != PRF_NALT)
      {
        fprintf(stderr,"ReadAmol: error, n=%d, expected=%d\n",n,PRF_NALT);
        err = 1;
      }
      if(err)
      {
        return -1;
      }
    }
  }

  return 0;
}

int ReadTrac(void)
{
  int i,n;
  int err;
  double v;
  char line[MAXLINE];
  char str1[MAXLINE];
  char fnam[MAXLINE];
  char *p;
  FILE *fp;

  for(i=0; i<PRF_NTRC; i++)
  {
    sprintf(fnam,"%s/TRAC_IALT_%d.dat",prf_ddir,i+1);
    if((fp=fopen(fnam,"r")) == NULL)
    {
      fprintf(stderr,"ReadTrac: cannot open %s\n",fnam);
      return -1;
    }
    n = 0;
    err = 0;
    while(fgets(line,MAXLINE,fp) != NULL)
    {
      if(sscanf(line,"%s",str1) != 1)
      {
        fprintf(stderr,"ReadTrac: read error >>> %s\n",line);
        err = 1;
        break;
      }
      errno = 0;
      v = strtod(str1,&p);
      if(errno==ERANGE || *p!='\0')
      {
        fprintf(stderr,"ReadTrac: convert error >>> %s\n",line);
        err = 1;
        break;
      }
      if(n >= PRF_NALT)
      {
        fprintf(stderr,"ReadTrac: error, #data exceed the limit >>> %d\n",n);
        err = 1;
        break;
      }
      prf_trac[i][n] = v;
      n++;
    }
    fclose(fp);
    if(n != PRF_NALT)
    {
      fprintf(stderr,"ReadTrac: error, n=%d, expected=%d\n",n,PRF_NALT);
      err = 1;
    }
    if(err)
    {
      return -1;
    }
  }

  return 0;
}
