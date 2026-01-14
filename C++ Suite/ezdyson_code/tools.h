#ifndef _tools_h
#define _tools_h

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <time.h>
#include <stddef.h>
#include <string.h>
#include <strings.h>

#include <string>
#include <iomanip>
#include <iostream>

//! The dimension of Cartesian space
#define XYZ 3
typedef enum {X,Y,Z} Carts;

//! Type of radial part (Spherical, Coulomb wave)

//Max and Min f-s
//! A macro that returns the maximum of a and b
#define MAX(a,b) ( (a>=b) ? (a) : (b) )
//! A macro that returns the minimum of a and b
#define MIN(a,b) ( (a<b) ? (a) : (b) ) 
//! A macro t hat returns the sign of as a (double)+/-1.
#define SIGN(a)  ( (a<0.) ? (-1.) : (1.) )
//! Normalized random numbers
#define RANDMAX ((double)pow(2.,31.)-1.)
//! random number in [0,1] 
#define RANDOM ((double)random()/RANDMAX) 

//! Error handling (print msg, exit(1))
inline void error(std::string msg) { error(msg.c_str());}
void error(const char *const msg);
//! if expr -> error(msg);
inline void check(bool expr, const char *const msg) {if (expr) error(msg);}
inline void check(bool expr, std::string msg) { check(expr,msg.c_str()); }

//! Randomizer (use system time)
void Randomize(void );
//! Swaping two integer numbers
inline void Swap(int& i, int& j)
{
  static int tmp;
  tmp=i;
  i=j;
  j=tmp;
}

#define BOHR_TO_ANGS 0.529177249
#define ANGS_TO_BOHR (1./BOHR_TO_ANGS)
#define HAR_TO_EV 27.2107
#define EV_TO_HAR (1./HAR_TO_EV)
#define MASS_ELEC 9.1093826E-31
#define HBAR 1.05457168E-34
#define EV_TO_J 1.60210E-19
#define C_AU 137.037

#define MEGABITE (1024*1024)

double const pi = M_PI;

//!Convert a string to all lower case.
void LowCase(char *str);
//!double factorial function
int dfac(int i);
//!factorial function
int fac(int i);

double Min(double a, double b);
double Max(double a, double b);

namespace WignerSymbols {

template <typename T>
double sgn(T val)
{
    int sgn = (T(0) < val) - (val < T(0));
    if (sgn == 0)
        return 1.0;
    else
        return (double)sgn;
}

} // namespace WignerSymbols

#endif

