#ifndef _complexno_h
#define _complexno_h

#include "tools.h"
#include <math.h>

//!Complex numbers
class Complex
{
  double real;
  double imag;

 public: 

  Complex(double re=0.0, double im=0.0) : real(re), imag(im) {}
  Complex(const Complex& other) : real(other.real), imag(other.imag) {}

  inline double& Re() { return real; }
  inline double& Im() { return imag; }
  inline double GetRe() const { return real; }
  inline double GetIm() const { return imag; }

  inline double Norm() const { return real*real+imag*imag; }
  inline Complex Conj() const { return Complex(real,-imag); }

  Complex& operator=(const Complex& other)
    {
      if(this!=&other)
	{
	  real=other.real;
	  imag=other.imag;
	}
      return *this; 
    }

  inline Complex& operator+=(const Complex& other)
    {
      real+=other.real;
      imag+=other.imag;
      return *this;
    }
  
  inline Complex& operator+=(double db)
    {
      real += db;
      return *this;
    }

  inline Complex& operator-=(const Complex& other)
    {
      real -= other.real;
      imag -= other.imag;
      return *this;
    }

  Complex& operator-=(double db)
    {
      real -= db;
      return *this;
    }

  inline Complex& operator*=(const Complex& other)
    {
      double tmp_real = real*other.real-imag*other.imag;
      double tmp_imag = real*other.imag+imag*other.real;
      real=tmp_real;
      imag=tmp_imag;
      return *this;
    }

  inline Complex& operator*=(double db)
    {
      real *= db;
      imag *= db;
      return *this;
    }

  inline Complex& operator/=(double db)
    {
      real /= db;
      imag /= db;
      return *this;
    }
  
  inline friend Complex operator+(const Complex& cx1, const Complex& cx2)
    {
      return Complex(cx1.real+cx2.real,cx1.imag+cx2.imag);
    }

  inline friend Complex operator-(const Complex& cx1, const Complex& cx2)
    {
      return Complex(cx1.real-cx2.real,cx1.imag-cx2.imag);
    }

  inline friend Complex operator*(const Complex& cx1, const Complex& cx2)
    {
      return Complex(cx1.real*cx2.real-cx1.imag*cx2.imag,
		     cx1.real*cx2.imag+cx1.imag*cx2.real);
      /*      
	      Complex tmp(cx1);
	      tmp *= cx2;
	      return tmp;
      */
    }

  inline friend Complex operator*(const Complex& cx1, double db)
    {
      return Complex(cx1.real*db,cx1.imag*db);
    }

  inline friend Complex operator*(double db, const Complex& cx1)
    {
      return Complex(cx1.real*db,cx1.imag*db);
    }
  
  inline friend Complex operator/(const Complex& cx1, double db)
    {
      return Complex(cx1.real/db,cx1.imag/db);
    }

  inline friend Complex CPower(const Complex& cx1, int n)
    {
      Complex cx2(1.0);
      // if(n>0) the loop will not execute, thus do not need this check...
      for(int i=0; i<n; i++)
	cx2 *= cx1;
      
      return cx2;
    }

  inline friend Complex CPower(double r, const Complex& cx1)
    {
      double tmparg=cx1.imag*log(r);
      double powr=pow(r,cx1.real);
      return  Complex(powr*cos(tmparg), powr*sin(tmparg));
    }
};



#endif


