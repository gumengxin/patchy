/***************************************************
* Scale velocity to make isokinetic
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

void ScaleVelocity(void) //  need to use after Kinetic()
{
 int i;
 double scale_trans,scale_rot;

 /**************************************************
 Kinstant and Tinstant need to be done by Kinetic()
 ***************************************************/
 scale_trans = sqrt(T/Ttrans);
 scale_rot = sqrt(T/Trot);

 for(i=0;i<NumberOfParticles;i++) 
 {
  //p
  PX[i] *= scale_trans;
  PY[i] *= scale_trans;
  PZ[i] *= scale_trans;
  //r'=v
  RX1[i] *= scale_trans;
  RY1[i] *= scale_trans;
  RZ1[i] *= scale_trans;

  WX[i] *= scale_rot;
  WY[i] *= scale_rot;
  WZ[i] *= scale_rot;
 }
 Ttrans = T;
 Ktrans = Nftrans*kB/2.0*Ttrans;
 Trot = T;
 Krot = Nfrot*kB/2.0*Trot;
 Tinstant = T;
 Kinstant = Nf*kB/2.0*Tinstant;
 
 return;
}
