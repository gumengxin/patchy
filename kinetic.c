/***************************************************
* Kinetic
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
//#include "ran_uniform.h"

void Kinetic(void)
{
 int i;

 Ktrans = 0.;
 Krot = 0.;
 PdotP = 0.;
 WdotL = 0.;
 for(i=0;i<NumberOfParticles;i++)
 {
   Ktrans += mass[i]*(SQR(RX1[i]) + SQR(RY1[i]) + SQR(RZ1[i]));
   Krot += IXX[i]*SQR(WX[i]) + IYY[i]*SQR(WY[i]) + IZZ[i]*SQR(WZ[i]);
   PdotP += (SQR(PX[i]) + SQR(PY[i]) + SQR(PZ[i]))/mass[i];
   WdotL += IXX[i]*SQR(WX[i]) + IYY[i]*SQR(WY[i]) + IZZ[i]*SQR(WZ[i]);
 }
 Ktrans *= 0.5; // instantenous translational kinetic energy
 Ttrans = 2.0*Ktrans/Nftrans/kB;
 Krot = 0.5*Krot;
 Trot = 2.0*Krot/Nfrot/kB;

 Kinstant = Ktrans+Krot;
 Tinstant = 2.0*Kinstant/Nf/kB;

 return;
}
