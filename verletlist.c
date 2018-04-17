/***************************************
 * make a new Verlet list
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void VerletCheck(void)
{
 int i;
 double drmax;

 drmax = 0.;
  
 for(i=0;i<NumberOfParticles;i++)
 {
  // use original coordinates RX RY RZ, which don't put back to box
  drmax = MAX(fabs(RX[i]-RX0[i]),drmax);
  drmax = MAX(fabs(RY[i]-RY0[i]),drmax);
  drmax = MAX(fabs(RZ[i]-RZ0[i]),drmax);
 }
 
 drmax = 2.*sqrt(3.*SQR(drmax));

 if(drmax > rv-rc) 
  VerletCheckIndex = 1;
 else
  VerletCheckIndex = 0;
  
 return;
}
