/*****************************************
 * Calculate radial distribution function
 * g(r), gAA(r), gBB(r)
 *****************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

void RadialDis(void)
{
 int i,j;
 double rij,rij2;

 for(i=0;i<NumberOfParticles-1;i++)
  for(j=i+1;j<NumberOfParticles;j++)
  {
   rij2 = Distance(i,j);
   rij = sqrt(rij2);

   if(rij < L/2.)
   {
    g[(int)(rij/dradial)] += 2.;
    if(identity[i] == identity[j])
    {
     if(identity[i] == 1) // AA
      gAA[(int)(rij/dradial)] += 2.;
     else //BB
      gBB[(int)(rij/dradial)] += 2.;
    }// endif same particle
    else //AB
      gAB[(int)(rij/dradial)] += 2.;
   }//endif r<L/2
  }
return;

}
