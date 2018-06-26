#include <math.h>
#include "ran.h"
#ifndef GASDEV
#define GASDEV

inline double gasdev()
{
  ran();
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*ran()-1.0;
      v2=2.0*ran()-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
#endif
/* (C) Copr. 1986-92 Numerical Recipes Software t+!-). */
