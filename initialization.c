/***************************************************
* initialize particle parameters
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void Initialization(int job)
{
 char filename[20];
 int i,j;
 int nx,ny,nz;
 double scale;
 double particlepercell;
 //double Kinstant,scale,Tinstant;
 double Mcom;
 VECTOR Vcom;
 VECTOR Omega;

 FILE *fp;


 for(i=0;i<NumberOfParticles;i++)
 {
  if(i<NA)
  {
  sigma[i] = sigmaA;
  epsilon[i] = epsilonA;
  mass[i] = massA;
  identity[i] = 1;
  zcoordinate[i] = ZA;
  }
  else
  {
  sigma[i] = sigmaB;
  epsilon[i] = epsilonB;
  mass[i] = massB;
  identity[i] = 2;
  zcoordinate[i] = ZB;
  }
  IXX[i] = 2./5.*mass[i]*SQR(sigma[i]/2.);
  IYY[i] = 2./5.*mass[i]*SQR(sigma[i]/2.);
  IZZ[i] = 2./5.*mass[i]*SQR(sigma[i]/2.);
 }
printf("I = (%lf %lf %lf)\n",IXX[0],IYY[0],IZZ[0]);

/****************  position **********/
 if(rInitialType == 0) // read from file
 {
  sprintf(filename,"position_%d",job);
  fp = fopen(filename,"r");
  for(i=0;i<NumberOfParticles;i++)
  fscanf(fp,"%lf\t%lf\t%lf\n",&RX[i],&RY[i],&RZ[i]);
  fclose(fp);
  
  sprintf(filename,"quaternion_%d",job);
  fp = fopen(filename,"r");
  for(i=0;i<NumberOfParticles;i++)
  fscanf(fp,"%lf\t%lf\t%lf\t%lf\n",&QW[i],&QX[i],&QY[i],&QZ[i]);
  fclose(fp);
 }// end read position from file
 else if(rInitialType == 1) // on lattice 
 {
  /***********************place particle on lattice **********************************/
  if(PackType == 1) // simple cubic
  {
   //latticeconst = pow(4.0*V/NumberOfLatticeSites,1.0/3.0);
   particlepercell = 1.;
   latticeconst = pow(1.0/rho,1.0/3.0);

   printf("SC with Nc = %d and a = %lf\n",NumberOfLatticeSites,latticeconst);
   i=0;
   do
   {
    for(nx=0;nx<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ),1.0/3.0)*AspectRatioX;nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ),1.0/3.0)*AspectRatioY;ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ),1.0/3.0)*AspectRatioZ;nz++) 	
    {
     RX[i] = nx*latticeconst;
     RY[i] = ny*latticeconst;
     RZ[i] = nz*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }
   }//end do
   while(i<NumberOfParticles);
  }// end if sc
  if(PackType == 2) // bcc
  {
   particlepercell = 2.;
   latticeconst = pow(2.0/rho,1.0/3.0);
   printf("BCC with Nc = %d and a = %lf\n",NumberOfLatticeSites,latticeconst);
   i=0;	
   do
   {
    for(nx=0;nx<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/2,1.0/3.0)*AspectRatioX;nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/2,1.0/3.0)*AspectRatioY;ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/2,1.0/3.0)*AspectRatioZ;nz++) 	
    {
     RX[i] = nx*latticeconst;
     RY[i] = ny*latticeconst;
     RZ[i] = nz*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }
    for(nx=0;nx<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/2,1.0/3.0)*AspectRatioX;nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/2,1.0/3.0)*AspectRatioY;ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/2,1.0/3.0)*AspectRatioZ;nz++) 	
    {
     RX[i] = nx*latticeconst+0.5*latticeconst;
     RY[i] = ny*latticeconst+0.5*latticeconst;
     RZ[i] = nz*latticeconst+0.5*latticeconst;
     i++;
     if(i >= NumberOfParticles) break;
    }
   } //end do
   while(i<NumberOfParticles);
  }//end if bcc
  if(PackType == 3) // fcc
  {
   //latticeconst = pow(4.0*V/NumberOfLatticeSites,1.0/3.0);
   particlepercell = 4.;
   latticeconst = pow(4.0/rho,1.0/3.0);

   printf("FCC with Nc = %d and a = %lf L=(%f %f %f)\n",NumberOfLatticeSites,latticeconst,pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioX,pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioY,pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioZ);
   i=0;
   do
   {
    for(nx=0;nx<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioX;nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioY;ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioZ;nz++) 	
    {
     RX[i] = nx*latticeconst;
     RY[i] = ny*latticeconst;
     RZ[i] = nz*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }
    for(nx=0;nx<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioX;nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioY;ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioZ;nz++) 	
    {
     RX[i] = nx*latticeconst+0.5*latticeconst;
     RY[i] = ny*latticeconst+0.5*latticeconst;
     RZ[i] = nz*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }
 
    for(nx=0;nx<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioX;nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioY;ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioZ;nz++) 	
    {
     RX[i] = nx*latticeconst+0.5*latticeconst;
     RY[i] = ny*latticeconst;
     RZ[i] = nz*latticeconst+0.5*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }

    for(nx=0;nx<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioX;nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioY;ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/(AspectRatioX*AspectRatioY*AspectRatioZ)/4,1.0/3.0)*AspectRatioZ;nz++) 	
    {
     RX[i] = nx*latticeconst;
     RY[i] = ny*latticeconst+0.5*latticeconst;
     RZ[i] = nz*latticeconst+0.5*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }
   }//end do
   while(i<NumberOfParticles);
  }// end if fcc

 for(i=0;i<NumberOfParticles;i++)
 {
  QW[i] = 1.;
  QX[i] = 0.;
  QY[i] = 0.;
  QZ[i] = 0.;
 }
 }// end if on lattice 
/******************************* end position***************************************/

   
/******************************* patch position***************************************/
for(i=0;i<NumberOfParticles;i++)
{
 if(zcoordinate[i] == 6)
 {
  dRSX[i][0] = sigma[i]/2.;
  dRSY[i][0] = 0.;
  dRSZ[i][0] = 0.;
  
  dRSX[i][1] = -sigma[i]/2.;
  dRSY[i][1] = 0.;
  dRSZ[i][1] = 0.;
  
  dRSX[i][2] = 0.;
  dRSY[i][2] = sigma[i]/2.;
  dRSZ[i][2] = 0.;
  
  dRSX[i][3] = 0.;
  dRSY[i][3] = -sigma[i]/2.;
  dRSZ[i][3] = 0.;
  
  dRSX[i][4] = 0.;
  dRSY[i][4] = 0.;
  dRSZ[i][4] = sigma[i]/2.;
  
  dRSX[i][5] = 0.;
  dRSY[i][5] = 0.;
  dRSZ[i][5] = -sigma[i]/2.;
 }// SC
 
 if(zcoordinate[i] == 8)
 {
  dRSX[i][0] = sigma[i]/2./sqrt(3.);
  dRSY[i][0] = sigma[i]/2./sqrt(3.);
  dRSZ[i][0] = sigma[i]/2./sqrt(3.);
  
  dRSX[i][1] = sigma[i]/2./sqrt(3.);
  dRSY[i][1] = -sigma[i]/2./sqrt(3.);
  dRSZ[i][1] = sigma[i]/2./sqrt(3.);
  
  dRSX[i][2] = sigma[i]/2./sqrt(3.);
  dRSY[i][2] = sigma[i]/2./sqrt(3.);
  dRSZ[i][2] = -sigma[i]/2./sqrt(3.);
  
  dRSX[i][3] = -sigma[i]/2./sqrt(3.);
  dRSY[i][3] = sigma[i]/2./sqrt(3.);
  dRSZ[i][3] = sigma[i]/2./sqrt(3.);
  
  dRSX[i][4] = -sigma[i]/2./sqrt(3.);
  dRSY[i][4] = -sigma[i]/2./sqrt(3.);
  dRSZ[i][4] = sigma[i]/2./sqrt(3.);
  
  dRSX[i][5] = -sigma[i]/2./sqrt(3.);
  dRSY[i][5] = sigma[i]/2./sqrt(3.);
  dRSZ[i][5] = -sigma[i]/2./sqrt(3.);
  
  dRSX[i][6] = sigma[i]/2./sqrt(3.);
  dRSY[i][6] = -sigma[i]/2./sqrt(3.);
  dRSZ[i][6] = -sigma[i]/2./sqrt(3.);
  
  dRSX[i][7] = -sigma[i]/2./sqrt(3.);
  dRSY[i][7] = -sigma[i]/2./sqrt(3.);
  dRSZ[i][7] = -sigma[i]/2./sqrt(3.);
 }// end bcc

 if(zcoordinate[i] == 12)
 {
  dRSX[i][0] = sigma[i]/2./sqrt(2.);
  dRSY[i][0] = sigma[i]/2./sqrt(2.);
  dRSZ[i][0] = 0.;
  dRSX[i][1] = sigma[i]/2./sqrt(2.);
  dRSY[i][1] = -sigma[i]/2./sqrt(2.);
  dRSZ[i][1] = 0.;
  dRSX[i][2] = -sigma[i]/2./sqrt(2.);
  dRSY[i][2] = sigma[i]/2./sqrt(2.);
  dRSZ[i][2] = 0.;
  dRSX[i][3] = -sigma[i]/2./sqrt(2.);
  dRSY[i][3] = -sigma[i]/2./sqrt(2.);
  dRSZ[i][3] = 0.;
  
  dRSX[i][4] = sigma[i]/2./sqrt(2.);
  dRSY[i][4] = 0.;
  dRSZ[i][4] = sigma[i]/2./sqrt(2.);
  dRSX[i][5] = sigma[i]/2./sqrt(2.);
  dRSY[i][5] = 0.;
  dRSZ[i][5] = -sigma[i]/2./sqrt(2.);
  dRSX[i][6] = -sigma[i]/2./sqrt(2.);
  dRSY[i][6] = 0.;
  dRSZ[i][6] = sigma[i]/2./sqrt(2.);
  dRSX[i][7] = -sigma[i]/2./sqrt(2.);
  dRSY[i][7] = 0.;
  dRSZ[i][7] = -sigma[i]/2./sqrt(2.);

  dRSX[i][8] = 0.;
  dRSY[i][8] = sigma[i]/2./sqrt(2.);
  dRSZ[i][8] = sigma[i]/2./sqrt(2.);
  dRSX[i][9] = 0.;
  dRSY[i][9] = sigma[i]/2./sqrt(2.);
  dRSZ[i][9] = -sigma[i]/2./sqrt(2.);
  dRSX[i][10] = 0.;
  dRSY[i][10] = -sigma[i]/2./sqrt(2.);
  dRSZ[i][10] = sigma[i]/2./sqrt(2.);
  dRSX[i][11] = 0.;
  dRSY[i][11] = -sigma[i]/2./sqrt(2.);
  dRSZ[i][11] = -sigma[i]/2./sqrt(2.);
 }// end fcc

 for(j=0;j<zcoordinate[i];j++)
 {
 RSX[i][j] = RX[i] + dRSX[i][j];
 RSY[i][j] = RY[i] + dRSY[i][j];
 RSZ[i][j] = RZ[i] + dRSZ[i][j];
 }// end loop j
}//end loop i
/******************************* end position***************************************/
  
/********************************** velocity ********************************************/

if(vInitialType == 0) // read from file
{
  sprintf(filename,"velocity_%d",job);
  fp = fopen(filename,"r");
  for(i=0;i<NumberOfParticles;i++)
  {
  fscanf(fp,"%lf\t%lf\t%lf\n",&RX1[i],&RY1[i],&RZ1[i]);
  PX[i] = RX1[i]*mass[i];
  PY[i] = RY1[i]*mass[i];
  PZ[i] = RZ1[i]*mass[i];
  }
  fclose(fp);
  
  sprintf(filename,"angularvelocity_%d",job);
  fp = fopen(filename,"r");
  for(i=0;i<NumberOfParticles;i++)
  fscanf(fp,"%lf\t%lf\t%lf\n",&WX[i],&WY[i],&WZ[i]);
  fclose(fp);
}
else if(vInitialType == 1) // random velocity
{
  /********************random velocity**************************************/
  Mcom = 0.;
  Vcom.x = 0.;
  Vcom.y = 0.;
  Vcom.z = 0.;

  for(i=0;i<NumberOfParticles;i++)
  {
   RX1[i] = BoxMuller(0.,1.);
   RY1[i] = BoxMuller(0.,1.);
   RZ1[i] = BoxMuller(0.,1.);

   Vcom.x += mass[i]*RX1[i];
   Vcom.y += mass[i]*RY1[i];
   Vcom.z += mass[i]*RZ1[i];

   Mcom += mass[i];
  }

   Vcom.x /= Mcom;
   Vcom.y /= Mcom;
   Vcom.z /= Mcom;

  Ktrans = 0.;
  for(i=0;i<NumberOfParticles;i++) 
  {
   RX1[i] -= Vcom.x;
   RY1[i] -= Vcom.y;
   RZ1[i] -= Vcom.z;

   PX[i] = RX1[i]*mass[i];
   PY[i] = RY1[i]*mass[i];
   PZ[i] = RZ1[i]*mass[i];

   Ktrans += mass[i]*(SQR(RX1[i]) + SQR(RY1[i]) + SQR(RZ1[i]));
  }
 
   Ktrans *= 0.5; // instantenous translational kinetic energy
   // 0.5*kT*Nf = K = 0.5* sum_mv^2, Nf = 3N-3
   Ttrans = 2.0*Ktrans/Nftrans/kB;
 //  Tinstant = 2.0*Kinstant/Nf/kB;
   scale = sqrt(T/Ttrans);
   for(i=0;i<NumberOfParticles;i++) 
   {
    RX1[i] *= scale;
    RY1[i] *= scale;
    RZ1[i] *= scale;
    
    PX[i] *= scale;
    PY[i] *= scale;
    PZ[i] *= scale;
   }
   
  Krot = 0.;
  for(i=0;i<NumberOfParticles;i++) 
  {
   if(i%2==0)
   {
   Omega = RandomSphere(sqrt(kB*T*Nfrot/NumberOfParticles));
   WX[i] = Omega.x/sqrt(IXX[i]);
   WY[i] = Omega.y/sqrt(IYY[i]);
   WZ[i] = Omega.z/sqrt(IZZ[i]);
   }
   else
   {
   WX[i] = -WX[i-1];
   WY[i] = -WY[i-1];
   WZ[i] = -WZ[i-1];
   }


  WX[i] = 0.;
  WY[i] = 0.;
  WZ[i] = 0.;

   Krot += IXX[i]*SQR(WX[i]) + IYY[i]*SQR(WY[i]) + IZZ[i]*SQR(WZ[i]);
  }

  Krot = 0.5*Krot;
  Trot = 2.0*Krot/Nfrot/kB;
}//end if maxwell velocity
 
// check velocity
  sprintf(filename,"vmaxwell0_%d.dat",JobIndex);
  fp=fopen(filename,"w");
  COMcheck(fp);
  fclose(fp);
  
  printf("Ktrans = %lf Ttrans = %lf Krot = %lf Trot = %lf\n",Ktrans,Ttrans,Krot,Trot);
 /****************** end velocity ***************************************/
 
 for(i=0;i<NumberOfParticles;i++)
 {
  QW1[i] = (-QX[i]*WX[i]-QY[i]*WY[i]-QZ[i]*WZ[i])*0.5;
  QX1[i] = ( QW[i]*WX[i]-QZ[i]*WY[i]+QY[i]*WZ[i])*0.5;
  QY1[i] = ( QZ[i]*WX[i]+QW[i]*WY[i]-QX[i]*WZ[i])*0.5;
  QZ1[i] = (-QY[i]*WX[i]+QX[i]*WY[i]+QW[i]*WZ[i])*0.5;
 }


return;
}

//check center of mass velocity
void COMcheck(FILE *fp)
{
 int i;
 double Mcom; // total mass
 VECTOR Vcom;

  Mcom = 0.;
  Vcom.x = 0.;
  Vcom.y = 0.;
  Vcom.z = 0.;
  Ktrans = 0.;
  for(i=0;i<NumberOfParticles;i++) 
  {
   Mcom += mass[i];
   Vcom.x += mass[i]*RX1[i];
   Vcom.y += mass[i]*RY1[i];
   Vcom.z += mass[i]*RZ1[i];
   Ktrans += mass[i]*(SQR(RX1[i]) + SQR(RY1[i]) + SQR(RZ1[i]));

   fprintf(fp,"%lf\n",sqrt(SQR(RX1[i]) + SQR(RY1[i]) + SQR(RZ1[i])));
  }

   Vcom.x /= Mcom;
   Vcom.y /= Mcom;
   Vcom.z /= Mcom;
   Ktrans *= 0.5;
   Ttrans = 2.0*Ktrans/Nftrans/kB;
   printf("center of mass velosity (%lf, %lf, %lf) and translational Ttrans = %lf\n",Vcom.x,Vcom.y,Vcom.z,Ttrans);
/********************************************************************/
 return;
}
 
