/********************************************
 * NVE ensemble 
 * Gear Predictor-Corrector
 * Predict Step
 * 4-value, 1st order ODE for translation
 * 5-value, 1st order ODE for rotation
 ********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void PredictNVE2(double DT)
{
 int i;
 double c1,c2,c3,c4;
 double norm;

 c1 = DT;
 c2 = c1*DT/2.0;
 c3 = c2*DT/3.0;
 c4 = c3*DT/4.0;

 for(i=0;i<NumberOfParticles;i++)
 {
  // original r
  RX[i]  = RX[i]  + c1*RX1[i] + c2*RX2[i] + c3*RX3[i];
  RY[i]  = RY[i]  + c1*RY1[i] + c2*RY2[i] + c3*RY3[i];
  RZ[i]  = RZ[i]  + c1*RZ1[i] + c2*RZ2[i] + c3*RZ3[i];
  // v = r'
  RX1[i] = RX1[i] + c1*RX2[i] + c2*RX3[i];
  RY1[i] = RY1[i] + c1*RY2[i] + c2*RY3[i];
  RZ1[i] = RZ1[i] + c1*RZ2[i] + c2*RZ3[i];
  // v' = r''
  RX2[i] = RX2[i] + c1*RX3[i];
  RY2[i] = RY2[i] + c1*RY3[i];
  RZ2[i] = RZ2[i] + c1*RZ3[i];
  
  // momentum p
  PX[i]  = PX[i]  + c1*PX1[i] + c2*PX2[i] + c3*PX3[i];
  PY[i]  = PY[i]  + c1*PY1[i] + c2*PY2[i] + c3*PY3[i];
  PZ[i]  = PZ[i]  + c1*PZ1[i] + c2*PZ2[i] + c3*PZ3[i];
  // p' = f
  PX1[i] = PX1[i] + c1*PX2[i] + c2*PX3[i];
  PY1[i] = PY1[i] + c1*PY2[i] + c2*PY3[i];
  PZ1[i] = PZ1[i] + c1*PZ2[i] + c2*PZ3[i];
  // p''
  PX2[i] = PX2[i] + c1*PX3[i];
  PY2[i] = PY2[i] + c1*PY3[i];
  PZ2[i] = PZ2[i] + c1*PZ3[i];

  QW[i] = QW[i] + c1*QW1[i] + c2*QW2[i] + c3*QW3[i] + c4*QW4[i];
  QX[i] = QX[i] + c1*QX1[i] + c2*QX2[i] + c3*QX3[i] + c4*QX4[i];
  QY[i] = QY[i] + c1*QY1[i] + c2*QY2[i] + c3*QY3[i] + c4*QY4[i];
  QZ[i] = QZ[i] + c1*QZ1[i] + c2*QZ2[i] + c3*QZ3[i] + c4*QZ4[i];
  QW1[i] = QW1[i] + c1*QW2[i] + c2*QW3[i] + c3*QW4[i];
  QX1[i] = QX1[i] + c1*QX2[i] + c2*QX3[i] + c3*QX4[i];
  QY1[i] = QY1[i] + c1*QY2[i] + c2*QY3[i] + c3*QY4[i];
  QZ1[i] = QZ1[i] + c1*QZ2[i] + c2*QZ3[i] + c3*QZ4[i];
  QW2[i] = QW2[i] + c1*QW3[i] + c2*QW4[i];
  QX2[i] = QX2[i] + c1*QX3[i] + c2*QX4[i];
  QY2[i] = QY2[i] + c1*QY3[i] + c2*QY4[i];
  QZ2[i] = QZ2[i] + c1*QZ3[i] + c2*QZ4[i];
  QW3[i] = QW3[i] + c1*QW4[i];
  QX3[i] = QX3[i] + c1*QX4[i];
  QY3[i] = QY3[i] + c1*QY4[i];
  QZ3[i] = QZ3[i] + c1*QZ4[i];

 // printf("%d %lf %lf %lf %lf = %lf\n",i,QW[i],QX[i],QY[i],QZ[i],SQR(QW[i])+SQR(QZ[i])+SQR(QY[i])+SQR(QZ[i]));

  WX[i] = WX[i] + c1*WX1[i] + c2*WX2[i] + c3*WX3[i] + c4*WX4[i];
  WY[i] = WY[i] + c1*WY1[i] + c2*WY2[i] + c3*WY3[i] + c4*WY4[i];
  WZ[i] = WZ[i] + c1*WZ1[i] + c2*WZ2[i] + c3*WZ3[i] + c4*WZ4[i];
  WX1[i] = WX1[i] + c1*WX2[i] + c2*WX3[i] + c3*WX4[i];
  WY1[i] = WY1[i] + c1*WY2[i] + c2*WY3[i] + c3*WY4[i];
  WZ1[i] = WZ1[i] + c1*WZ2[i] + c2*WZ3[i] + c3*WZ4[i];
  WX2[i] = WX2[i] + c1*WX3[i] + c2*WX4[i];
  WY2[i] = WY2[i] + c1*WY3[i] + c2*WY4[i];
  WZ2[i] = WZ2[i] + c1*WZ3[i] + c2*WZ4[i];
  WX3[i] = WX3[i] + c1*WX4[i];
  WY3[i] = WY3[i] + c1*WY4[i];
  WZ3[i] = WZ3[i] + c1*WZ4[i];

 // norm = SQR(QW[i])+SQR(QX[i])+SQR(QY[i])+SQR(QZ[i]);
 // norm = sqrt(norm);
 // QW[i] = QW[i]/norm;
 // QX[i] = QX[i]/norm;
 // QY[i] = QY[i]/norm;
 // QZ[i] = QZ[i]/norm;
 }
 return;
}

/********************************************
 * NVE ensemble
 * Gear Predictor-Corrector
 * Correct Step
 * 4-value, 1st order ODE for translation
 * 5-value, 1st order ODE for rotation
 ********************************************/

void CorrectNVE2(double DT)
{
 int i;
 double c1,c2,c3,c4;
 double GEART0,GEART1,GEART2,GEART3; // translation
 double GEARR0,GEARR1,GEARR2,GEARR3,GEARR4; // rotation

 double ctran0,ctran1,ctran2,ctran3; // translation
 double crot0,crot1,crot2,crot3,crot4; // rotation

 double corrw,corrx, corry,corrz;
 //double RX2I,RY2I,RZ2I;
 double RX1I,RY1I,RZ1I;
 double PX1I,PY1I,PZ1I;
 double QW1I,QX1I,QY1I,QZ1I;
 double WX1I,WY1I,WZ1I;

 double norm;

 // 2nd order coefficient
 //GEART0 = 1.0/6.0;
 //GEART1 = 5.0/6.0; 
 //GEART2 = 1.0; // not used
 //GEART3 = 1.0/3.0;

 GEART0 = 3.0/8.0;
 GEART1 = 1.0; // not used
 GEART2 = 3.0/4.0;
 GEART3 = 1.0/6.0;
 
 GEARR0 = 251.0/720.0;
 GEARR1 = 1.0; // not used
 GEARR2 = 11.0/12.0; 
 GEARR3 = 1.0/3.0;
 GEARR4 = 1.0/24.0;
 
 c1 = DT;
 c2 = c1*DT/2.0;
 c3 = c2*DT/3.0;
 c4 = c3*DT/4.0;

 // 2nd order coefficient
 //ctran0 = GEART0*c2;
 //ctran1 = GEART1*c2/c1;
 //ctran2 = 1.0; // not used
 //ctran3 = GEART3*c2/c3;
 
 ctran0 = GEART0*c1;
 ctran1 = 1.0; // not used
 ctran2 = GEART2*c1/c2; 
 ctran3 = GEART3*c1/c3;

 crot0 = GEARR0*c1;
 crot1 = 1.0; // not used
 crot2 = GEARR2*c1/c2;
 crot3 = GEARR3*c1/c3;
 crot4 = GEARR4*c1/c4;

 for(i=0;i<NumberOfParticles;i++)
 {
  // r'' = f/m
  //RX2I = force[i].x/mass[i];
  //RY2I = force[i].y/mass[i];
  //RZ2I = force[i].z/mass[i];
  //corrx = RX2I - RX2[i];
  //corry = RY2I - RY2[i];
  //corrz = RZ2I - RZ2[i];
  // r
  //RX[i] = RX[i] + ctran0*corrx;
  //RY[i] = RY[i] + ctran0*corry;
  //RZ[i] = RZ[i] + ctran0*corrz;
  // r'
  //RX1[i] = RX1[i] + ctran1*corrx; 
  //RY1[i] = RY1[i] + ctran1*corry; 
  //RZ1[i] = RZ1[i] + ctran1*corrz; 
  // r''
  //RX2[i] = RX2I;
  //RY2[i] = RY2I;
  //RZ2[i] = RZ2I;
  // r'''
  //RX3[i] = RX3[i] + ctran3*corrx;
  //RY3[i] = RY3[i] + ctran3*corry;
  //RZ3[i] = RZ3[i] + ctran3*corrz;

  // r' = p/m
  RX1I = PX[i]/mass[i];
  RY1I = PY[i]/mass[i];
  RZ1I = PZ[i]/mass[i];
  corrx = RX1I - RX1[i];
  corry = RY1I - RY1[i];
  corrz = RZ1I - RZ1[i];
  RX[i] = RX[i] + ctran0*corrx;
  RY[i] = RY[i] + ctran0*corry;
  RZ[i] = RZ[i] + ctran0*corrz;
  RX1[i] = RX1I; 
  RY1[i] = RY1I; 
  RZ1[i] = RZ1I; 
  RX2[i] = RX2[i] + ctran2*corrx;
  RY2[i] = RY2[i] + ctran2*corry;
  RZ2[i] = RZ2[i] + ctran2*corrz;
  RX3[i] = RX3[i] + ctran3*corrx;
  RY3[i] = RY3[i] + ctran3*corry;
  RZ3[i] = RZ3[i] + ctran3*corrz;

  //p' = f
  PX1I = force[i].x; // - XI*PX[i];
  PY1I = force[i].y; // - XI*PY[i];
  PZ1I = force[i].z; // - XI*PZ[i];
  corrx = PX1I - PX1[i];
  corry = PY1I - PY1[i];
  corrz = PZ1I - PZ1[i];
  PX[i] = PX[i] + ctran0*corrx;
  PY[i] = PY[i] + ctran0*corry;
  PZ[i] = PZ[i] + ctran0*corrz;
  PX1[i] = PX1I;
  PY1[i] = PY1I;
  PZ1[i] = PZ1I;
  PX2[i] = PX2[i] + ctran2*corrx;
  PY2[i] = PY2[i] + ctran2*corry;
  PZ2[i] = PZ2[i] + ctran2*corrz;
  PX3[i] = PX3[i] + ctran3*corrx;
  PY3[i] = PY3[i] + ctran3*corry;
  PZ3[i] = PZ3[i] + ctran3*corrz;

  QW1I = (-QX[i]*WX[i]-QY[i]*WY[i]-QZ[i]*WZ[i])*0.5;
  QX1I = ( QW[i]*WX[i]-QZ[i]*WY[i]+QY[i]*WZ[i])*0.5;
  QY1I = ( QZ[i]*WX[i]+QW[i]*WY[i]-QX[i]*WZ[i])*0.5;
  QZ1I = (-QY[i]*WX[i]+QX[i]*WY[i]+QW[i]*WZ[i])*0.5;
  corrw = QW1I - QW1[i];
  corrx = QX1I - QX1[i];
  corry = QY1I - QY1[i];
  corrz = QZ1I - QZ1[i];
  QW[i] = QW[i] + crot0*corrw;
  QX[i] = QX[i] + crot0*corrx;
  QY[i] = QY[i] + crot0*corry;
  QZ[i] = QZ[i] + crot0*corrz;
  QW1[i] = QW1I;
  QX1[i] = QX1I;
  QY1[i] = QY1I;
  QZ1[i] = QZ1I;
  QW2[i] = QW2[i] + crot2*corrw;
  QX2[i] = QX2[i] + crot2*corrx;
  QY2[i] = QY2[i] + crot2*corry;
  QZ2[i] = QZ2[i] + crot2*corrz;
  QW3[i] = QW3[i] + crot3*corrw;
  QX3[i] = QX3[i] + crot3*corrx;
  QY3[i] = QY3[i] + crot3*corry;
  QZ3[i] = QZ3[i] + crot3*corrz;
  QW4[i] = QW4[i] + crot4*corrw;
  QX4[i] = QX4[i] + crot4*corrx;
  QY4[i] = QY4[i] + crot4*corry;
  QZ4[i] = QZ4[i] + crot4*corrz;

  // Wx' = T/IXX + (IYY-IZZ)/IXX*Wy*Wz
 // WX1I = (torque[i].x+WY[i]*WZ[i]*(IYY[i]-IZZ[i]))/IXX[i];
 // WY1I = (torque[i].y+WZ[i]*WX[i]*(IZZ[i]-IXX[i]))/IYY[i];
 // WZ1I = (torque[i].z+WX[i]*WY[i]*(IXX[i]-IYY[i]))/IZZ[i];
  WX1I = torque[i].x/IXX[i];
  WY1I = torque[i].y/IYY[i];
  WZ1I = torque[i].z/IZZ[i];
  corrx = WX1I-WX1[i];
  corry = WY1I-WY1[i];
  corrz = WZ1I-WZ1[i];
  WX[i] = WX[i] + crot0*corrx;
  WY[i] = WY[i] + crot0*corry;
  WZ[i] = WZ[i] + crot0*corrz;
  WX1[i] = WX1I;
  WY1[i] = WY1I;
  WZ1[i] = WZ1I;
  WX2[i] = WX2[i] + crot2*corrx;
  WY2[i] = WY2[i] + crot2*corry;
  WZ2[i] = WZ2[i] + crot2*corrz;
  WX3[i] = WX3[i] + crot3*corrx;
  WY3[i] = WY3[i] + crot3*corry;
  WZ3[i] = WZ3[i] + crot3*corrz;
  WX4[i] = WX4[i] + crot4*corrx;
  WY4[i] = WY4[i] + crot4*corry;
  WZ4[i] = WZ4[i] + crot4*corrz;
 
  norm = SQR(QW[i])+SQR(QX[i])+SQR(QY[i])+SQR(QZ[i]);
  norm = sqrt(norm);
  QW[i] = QW[i]/norm;
  QX[i] = QX[i]/norm;
  QY[i] = QY[i]/norm;
  QZ[i] = QZ[i]/norm;
  
 }

 return;
}
