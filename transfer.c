/***********************************************************
 * molecular coordinates -> patch coordinates
 * patch coordinates -> molecular coordinates
 ***********************************************************/

#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void MolAtm(void)
{
 int i,j;
 double AXX,AXY,AXZ,AYX,AYY,AYZ,AZX,AZY,AZZ;// rotation matrix
 
 for(i=0;i<NumberOfParticles;i++)
 {
  AXX = QW[i]*QW[i] + QX[i]*QX[i] - QY[i]*QY[i] - QZ[i]*QZ[i];
  AXY = 2.0*(QX[i]*QY[i]+QW[i]*QZ[i]);
  AXZ = 2.0*(QX[i]*QZ[i]-QW[i]*QY[i]);
  AYX = 2.0*(QX[i]*QY[i]-QW[i]*QZ[i]);
  AYY = QW[i]*QW[i] - QX[i]*QX[i] + QY[i]*QY[i] - QZ[i]*QZ[i];
  AYZ = 2.0*(QY[i]*QZ[i]+QW[i]*QX[i]);
  AZX = 2.0*(QX[i]*QZ[i]+QW[i]*QY[i]);
  AZY = 2.0*(QY[i]*QZ[i]-QW[i]*QX[i]);
  AZZ = QW[i]*QW[i] - QX[i]*QX[i] - QY[i]*QY[i] + QZ[i]*QZ[i];

  for(j=0;j<zcoordinate[i];j++)
  {
  RSX[i][j] = RX[i] + AXX*dRSX[i][j] + AYX*dRSY[i][j] + AZX*dRSZ[i][j];
  RSY[i][j] = RY[i] + AXY*dRSX[i][j] + AYY*dRSY[i][j] + AZY*dRSZ[i][j];
  RSZ[i][j] = RZ[i] + AXZ*dRSX[i][j] + AYZ*dRSY[i][j] + AZZ*dRSZ[i][j];
  }
 }

 return;
}

void AtmMol(void)
{
 int i,j;
 double AXX,AXY,AXZ,AYX,AYY,AYZ,AZX,AZY,AZZ;// rotation matrix
  //for(j=0;j<zcoordinate[i];j++)

 double tx,ty,tz;
 
 PdotF = 0.;
 TdotL = 0.;
 for(i=0;i<NumberOfParticles;i++)
 {
  AXX = QW[i]*QW[i] + QX[i]*QX[i] - QY[i]*QY[i] - QZ[i]*QZ[i];
  AXY = 2.0*(QX[i]*QY[i]+QW[i]*QZ[i]);
  AXZ = 2.0*(QX[i]*QZ[i]-QW[i]*QY[i]);
  AYX = 2.0*(QX[i]*QY[i]-QW[i]*QZ[i]);
  AYY = QW[i]*QW[i] - QX[i]*QX[i] + QY[i]*QY[i] - QZ[i]*QZ[i];
  AYZ = 2.0*(QY[i]*QZ[i]+QW[i]*QX[i]);
  AZX = 2.0*(QX[i]*QZ[i]+QW[i]*QY[i]);
  AZY = 2.0*(QY[i]*QZ[i]-QW[i]*QX[i]);
  AZZ = QW[i]*QW[i] - QX[i]*QX[i] - QY[i]*QY[i] + QZ[i]*QZ[i];

  tx = torque[i].x; ty=torque[i].y; tz=torque[i].z;
  torque[i].x = AXX*tx+AXY*ty+AXZ*tz;
  torque[i].y = AYX*tx+AYY*ty+AYZ*tz;
  torque[i].z = AZX*tx+AZY*ty+AZZ*tz;

  TdotL += torque[i].x*WX[i] + torque[i].y*WY[i] + torque[i].z*WZ[i];
  PdotF += (PX[i]*force[i].x + PY[i]*force[i].y + PZ[i]*force[i].z)/mass[i]; 
 }
 return;
}

