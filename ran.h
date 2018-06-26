#ifndef RAN
#define RAN
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

#include<time.h>

inline double ran()
{
  static int initial=0;
  static long idum;
  time_t m_time;
  struct tm *jikoku;
  long k;
  double ans;

  if(initial==0){
    initial=1;
    time(&m_time);
    jikoku=localtime(&m_time);
    idum=jikoku[0].tm_sec;
    //printf("ran begins.\n");
  }    
  
  idum ^= MASK;
  k=(idum)/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum < 0) idum += IM;
  ans=AM*(idum);
  idum ^= MASK;
  return ans;
}

#endif
/* (C) Copr. 1986-92 Numerical Recipes Software t+!-). */
