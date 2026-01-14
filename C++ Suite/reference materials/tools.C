#include "tools.h"
#include <ctype.h>

//!Convert a string to all lower case.
void LowCase(char *str)
{
  for (char *iter=&str[0]; *iter != '\0'; ++iter)
      *iter = tolower(*iter);
}

//! double factorial function
int dfac(int i)
{
  //return ((i<=2) ? i : i*fac(i-2));
  int dfaci=1;
  /*
  if (i<=0)
    dfaci = 1;
  else
    for (int j=i; j>0; j-=2)
      dfaci *= j;
  */

  for(int j=i; j>0; j-=2)
    dfaci *= j;

  return dfaci;
  
}

//! factorial function
int fac(int i)
{
  //return ((i<=1) ? i : i*fac(i-1));
  int faci = 1;
  
  /*
  if (i==0)
    faci = 1;
  else
    for (int j=1; j<=i; j++)
    faci *= j; */

  for (int j=1; j<=i; j++)
    faci *= j; 
  
  return faci;
  
}


void error(const char *const msg)
{
  //perror(msg);
  //AIK: send all error messages to stdout, not stderr
  printf(msg);
  exit(1);
}


double Min (double a, double b)
{
  double c;
  c = (a < b) ? a : b;
  return c;
}

double Max(double a, double b)
{
  double c;
  c = (a > b) ? a : b;
  return c;
}

