/***********************************************************
 * generate a random vector on a unit sphere
 * Marsaglia 1972
 * page 349 on Allen and Tildesley 
 ***********************************************************/

#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

VECTOR RandomSphere(double norm)
{
 VECTOR r;

 double xx1,xx2; 
 double yy1,yy2;
 double yysq;

 do
 {
 xx1 = RandomNumber();
 xx2 = RandomNumber();
 yy1 = 1.-2.*xx1;
 yy2 = 1.-2.*xx2;
 yysq = yy1*yy1+yy2*yy2;
 }
 while(yysq > 1.);

 r.x = 2.*yy1*sqrt(1.-yysq);
 r.y = 2.*yy2*sqrt(1.-yysq);
 r.z = 1.-2.*yysq;

 r.x = r.x*norm;
 r.y = r.y*norm;
 r.z = r.z*norm;
 return r;
}
