#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#define	PI			3.141592653589793		// PI
#define	PI2			6.283185307179586		// 2*PI
#define	R_TO_D			57.29577951308232		// 180/PI rad -> deg
#define	D_TO_R			1.745329251994329e-02		// PI/180 deg -> rad
#define	MAXNANG			1000				// Max #angles
#define	MAXLINE			256				// Max #chars in a line
#define	EPSILON			1.0e-14				// A small number
#define	DELTA			1.0e-13				// A small number

double angl[MAXNANG];
double cang[MAXNANG];
double dang[MAXNANG];
double sang[MAXNANG];
double phas[MAXNANG];

int main(void)
{
  int i;
  int nang;
  double x,y;
  double sp;
  double asym;
  char line[MAXLINE];
  char str1[MAXLINE];
  char str2[MAXLINE];
  char *p;

  i = 0;
  while(fgets(line,MAXLINE,stdin) != NULL)
  {
    if(sscanf(line,"%s%s",str1,str2) != 2)
    {
      fprintf(stderr,"Read error >>> %s\n",line);
      return -1;
    }
    errno = 0;
    x = strtod(str1,&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      return -1;
    }
    errno = 0;
    y = strtod(str2,&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      return -1;
    }
    if(i >= MAXNANG)
    {
      fprintf(stderr,"Warning, #data exceed the limit >>> %d\n",i);
      break;
    }
    angl[i] = x*D_TO_R;
    phas[i] = y;
    i++;
  }
  nang = i;

  for(i=0; i<nang; i++)
  {
    sang[i] = sin(angl[i]);
    cang[i] = cos(angl[i]);
  }
  for(i=1; i<nang; i++)
  {
    dang[i] = PI*fabs(angl[i]-angl[i-1]);
  }
  sp = 0.0;
  for(i=1; i<nang; i++)
  {
    sp += (phas[i-1]*sang[i-1]+phas[i]*sang[i])*dang[i];
  }
  for(i=0; i<nang; i++)
  {
    phas[i] /= sp;
  }
  asym = 0.0;
  for(i=1; i<nang; i++)
  {
    asym += (phas[i-1]*sang[i-1]*cang[i-1]+phas[i]*sang[i]*cang[i])*dang[i];
  }
  printf("%13.9f %13.10f\n",sp,asym);

  return 0;
}
