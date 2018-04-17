/***************************************
 * calculate of force fx,fy,fz of 
 * each particle i at corrent position t
 * calculate potential energy Upotential
 * and virial = 1/d <sum f*r> as well
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void Force(void) 
{
  int i,j;
  int jbegin,jend,jlist;
  double rvij,rvij2; // Verlet list radius
  double rcij,rcij2; // rcutoff
  double sigmaij,sigmaij2;
  double rij,rij2,rij6; // avoid sqrt() if possible
  double xij,yij,zij;
 // double dudr; // du/dr
  double fij;
  double fxij,fyij,fzij;
  double txij,tyij,tzij; // torque
  double txji,tyji,tzji; // torque
  double ftxij,ftyij,ftzij; // force contribute to torque
  double ftxji,ftyji,ftzji; // force contribute to torque
  
  double pxij,pyij,pzij;
  double wij,uij,Xij,rijpij,pijfij;
   
  int ip,jp; // patch index
  int ip0,jp0; // patch index
  double thetai,thetaj;
  double prod,prodmaxi,prodmaxj;

  double dcosdx,dcosdy,dcosdz;
  double urA;

 // Virial = 0.;
 // HyperVirial = 0.;
  PdotF = 0.;
  TdotL = 0.;
 // RdotPdotX = 0.;
  Upotential = 0.;

  for(i=0;i<NumberOfParticles;i++)
  {
    force[i].x = 0.;
    force[i].y = 0.;
    force[i].z = 0.;
    torque[i].x = 0.;
    torque[i].y = 0.;
    torque[i].z = 0.;
  }
  
  if(VerletCheckIndex == 1)   
  {
   VerletUpdate++;

   for(i=0;i<NumberOfParticles;i++)
   {
    RX0[i] = RX[i];
    RY0[i] = RY[i];
    RZ0[i] = RZ[i];
   } // end store old position loop

  nlist = -1; // Verlet list index
  for(i=0;i<NumberOfParticles-1;i++)
  {
   point[i] = nlist + 1; //starting index of particle i in the list
   for(j=i+1;j<NumberOfParticles;j++)
   {
     sigmaij = SigmaIJ(identity[i],identity[j]);
     rvij = rv/sigmaA*sigmaij;
     rvij2 = SQR(rvij);
     rcij = rc/sigmaA*sigmaij;
     rcij2 = SQR(rcij);
     // rij = ri - rj from j to i 
     xij = RX[i] - RX[j];
     yij = RY[i] - RY[j];
     zij = RZ[i] - RZ[j];
     xij = xij - LX*round(xij/LX);
     yij = yij - LY*round(yij/LY);
     zij = zij - LZ*round(zij/LZ);
     rij2 = SQR(xij)+SQR(yij)+SQR(zij);

    if(rij2 < rvij2) // if within Verlet neighbor cell
    {
      nlist++;
      VerletList[nlist] = j;
     if(rij2 < rcij2) // try avoid sqrt for r>rc
     {
       rij = sqrt(rij2); // avoid sqrt if possible

       prodmaxi = -1.;
       for(ip=0;ip<zcoordinate[i];ip++)
       {
        // dRS.rji
        prod = (RSX[i][ip]-RX[i])*(-xij)+(RSY[i][ip]-RY[i])*(-yij)+(RSZ[i][ip]-RZ[i])*(-zij);
        if(prod>prodmaxi)
        {
         prodmaxi = prod; 
         thetai = (prod/rij/(sigma[i]/2.)); // cos(theta)
         ip0 = ip;
        }
       }
       
       prodmaxj = -1.;
       for(ip=0;ip<zcoordinate[j];ip++)
       {
        // dRS.rij
        prod = (RSX[j][ip]-RX[j])*(xij)+(RSY[j][ip]-RY[j])*(yij)+(RSZ[j][ip]-RZ[j])*(zij);
        if(prod>prodmaxj)
        {
         prodmaxj = prod; 
         thetaj = (prod/rij/(sigma[j]/2.)); // cos(theta)
         jp0 = ip;
        }
       }

       uij =  Potential(identity[i],identity[j],rij,thetai,thetaj);
       Upotential += uij;
       // force on center of mass, translational accelaration
       fij =  -Potential_dr(identity[i],identity[j],rij,thetai,thetaj); // along rij
       // along inter particle distance rij
       fxij = fij*xij/rij;
       fyij = fij*yij/rij;
       fzij = fij*zij/rij;
       //du/dthetai and du/dthetaj
       urA= Potential_Attract(identity[i],identity[j],rij,thetai,thetaj);
       dcosdx = (-(RSX[i][ip0]-RX[i])/rij-prodmaxi*xij/CUBIC(rij))/(sigma[i]/2.);
       dcosdy = (-(RSY[i][ip0]-RY[i])/rij-prodmaxi*yij/CUBIC(rij))/(sigma[i]/2.);
       dcosdz = (-(RSZ[i][ip0]-RZ[i])/rij-prodmaxi*zij/CUBIC(rij))/(sigma[i]/2.);
       fxij += -urA/SQR(patchsize)*2*(1.-thetai)*dcosdx;
       fyij += -urA/SQR(patchsize)*2*(1.-thetai)*dcosdy;
       fzij += -urA/SQR(patchsize)*2*(1.-thetai)*dcosdz;
       dcosdx = ((RSX[j][jp0]-RX[j])/rij-prodmaxj*xij/CUBIC(rij))/(sigma[j]/2.);
       dcosdy = ((RSY[j][jp0]-RY[j])/rij-prodmaxj*yij/CUBIC(rij))/(sigma[j]/2.);
       dcosdz = ((RSZ[j][jp0]-RZ[j])/rij-prodmaxj*zij/CUBIC(rij))/(sigma[j]/2.);
       fxij += -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdx;
       fyij += -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdy;
       fzij += -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdz;
       
       force[i].x += fxij;
       force[i].y += fyij;
       force[i].z += fzij;

       force[j].x += -fxij;
       force[j].y += -fyij;
       force[j].z += -fzij;
       
       // dcos/dsi
       dcosdx = (-xij/(sigma[i]/2.)-prodmaxi*(RSX[i][ip0]-RX[i])/CUBIC(sigma[i]/2.))/rij;
       dcosdy = (-yij/(sigma[i]/2.)-prodmaxi*(RSY[i][ip0]-RY[i])/CUBIC(sigma[i]/2.))/rij;
       dcosdz = (-zij/(sigma[i]/2.)-prodmaxi*(RSZ[i][ip0]-RZ[i])/CUBIC(sigma[i]/2.))/rij;
       ftxij = -urA/SQR(patchsize)*2*(1.-thetai)*dcosdx;
       ftyij = -urA/SQR(patchsize)*2*(1.-thetai)*dcosdy;
       ftzij = -urA/SQR(patchsize)*2*(1.-thetai)*dcosdz;
       // dcos/dsj 
       dcosdx = (xij/(sigma[j]/2.)-prodmaxj*(RSX[j][jp0]-RX[j])/CUBIC(sigma[j]/2.))/rij;
       dcosdy = (yij/(sigma[j]/2.)-prodmaxj*(RSY[j][jp0]-RY[j])/CUBIC(sigma[j]/2.))/rij;
       dcosdz = (zij/(sigma[j]/2.)-prodmaxj*(RSZ[j][jp0]-RZ[j])/CUBIC(sigma[j]/2.))/rij;
       ftxji = -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdx;
       ftyji = -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdy;
       ftzji = -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdz;
       // torque on i due to j 
       txij = (RSY[i][ip0]-RY[i])*ftzij - (RSZ[i][ip0]-RZ[i])*ftyij;
       tyij = (RSZ[i][ip0]-RZ[i])*ftxij - (RSX[i][ip0]-RX[i])*ftzij;
       tzij = (RSX[i][ip0]-RX[i])*ftyij - (RSY[i][ip0]-RY[i])*ftxij;
       // torque on j due to i 
       txji = (RSY[j][jp0]-RY[j])*ftzji - (RSZ[j][jp0]-RZ[j])*ftyji;
       tyji = (RSZ[j][jp0]-RZ[j])*ftxji - (RSX[j][jp0]-RX[j])*ftzji;
       tzji = (RSX[j][jp0]-RX[j])*ftyji - (RSY[j][jp0]-RY[j])*ftxji;
       
       torque[i].x += txij; 
       torque[i].y += tyij; 
       torque[i].z += tzij; 
       
       torque[j].x += txji; 
       torque[j].y += tyji; 
       torque[j].z += tzji; 
     }//end if rij < rc
    }//end if rij < rv
   }//end loop j
  }// end loop i
   point[NumberOfParticles-1] = nlist +1;
 } // endif VerletIndexCheck = 1 and update
 else if(VerletCheckIndex == 0)   // no need to update, Verlet list is used  
 {
  for(i=0;i<NumberOfParticles-1;i++)
  {
   jbegin = point[i];
   jend = point[i+1]-1;
   if(jbegin<=jend) // if list is not empty
   {
    for(jlist=jbegin;jlist<=jend;jlist++)
    {
     j = VerletList[jlist];  
     sigmaij = SigmaIJ(identity[i],identity[j]);
     rcij = rc/sigmaA*sigmaij;
     rcij2 = SQR(rcij);
     // rij = ri - rj from j to i 
     xij = RX[i] - RX[j];
     yij = RY[i] - RY[j];
     zij = RZ[i] - RZ[j];
     xij = xij - LX*round(xij/LX);
     yij = yij - LY*round(yij/LY);
     zij = zij - LZ*round(zij/LZ);
     rij2 = SQR(xij)+SQR(yij)+SQR(zij);

     if(rij2 < rcij2) // try avoid sqrt for r>rc
     {
       rij = sqrt(rij2); // avoid sqrt if possible

       prodmaxi = -1.;
       for(ip=0;ip<zcoordinate[i];ip++)
       {
        // dRS.rji
        prod = (RSX[i][ip]-RX[i])*(-xij)+(RSY[i][ip]-RY[i])*(-yij)+(RSZ[i][ip]-RZ[i])*(-zij);
        if(prod>prodmaxi)
        {
         prodmaxi = prod; 
         thetai = (prod/rij/(sigma[i]/2.)); // cos(theta)
         ip0 = ip;
        }
       }
       
       prodmaxj = -1.;
       for(ip=0;ip<zcoordinate[j];ip++)
       {
        // dRS.rij
        prod = (RSX[j][ip]-RX[j])*(xij)+(RSY[j][ip]-RY[j])*(yij)+(RSZ[j][ip]-RZ[j])*(zij);
        if(prod>prodmaxj)
        {
         prodmaxj = prod; 
         thetaj = (prod/rij/(sigma[j]/2.)); // cos(theta)
         jp0 = ip;
        }
       }

       uij =  Potential(identity[i],identity[j],rij,thetai,thetaj);
       Upotential += uij;
       // force on center of mass, translational accelaration
       fij =  -Potential_dr(identity[i],identity[j],rij,thetai,thetaj); // along rij
       // along inter particle distance rij
       fxij = fij*xij/rij;
       fyij = fij*yij/rij;
       fzij = fij*zij/rij;
       //du/dthetai and du/dthetaj
       urA= Potential_Attract(identity[i],identity[j],rij,thetai,thetaj);
       dcosdx = (-(RSX[i][ip0]-RX[i])/rij-prodmaxi*xij/CUBIC(rij))/(sigma[i]/2.);
       dcosdy = (-(RSY[i][ip0]-RY[i])/rij-prodmaxi*yij/CUBIC(rij))/(sigma[i]/2.);
       dcosdz = (-(RSZ[i][ip0]-RZ[i])/rij-prodmaxi*zij/CUBIC(rij))/(sigma[i]/2.);
       fxij += -urA/SQR(patchsize)*2*(1.-thetai)*dcosdx;
       fyij += -urA/SQR(patchsize)*2*(1.-thetai)*dcosdy;
       fzij += -urA/SQR(patchsize)*2*(1.-thetai)*dcosdz;
       dcosdx = ((RSX[j][jp0]-RX[j])/rij-prodmaxj*xij/CUBIC(rij))/(sigma[j]/2.);
       dcosdy = ((RSY[j][jp0]-RY[j])/rij-prodmaxj*yij/CUBIC(rij))/(sigma[j]/2.);
       dcosdz = ((RSZ[j][jp0]-RZ[j])/rij-prodmaxj*zij/CUBIC(rij))/(sigma[j]/2.);
       fxij += -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdx;
       fyij += -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdy;
       fzij += -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdz;
       
       force[i].x += fxij;
       force[i].y += fyij;
       force[i].z += fzij;

       force[j].x += -fxij;
       force[j].y += -fyij;
       force[j].z += -fzij;
       
       // dcos/dsi
       dcosdx = (-xij/(sigma[i]/2.)-prodmaxi*(RSX[i][ip0]-RX[i])/CUBIC(sigma[i]/2.))/rij;
       dcosdy = (-yij/(sigma[i]/2.)-prodmaxi*(RSY[i][ip0]-RY[i])/CUBIC(sigma[i]/2.))/rij;
       dcosdz = (-zij/(sigma[i]/2.)-prodmaxi*(RSZ[i][ip0]-RZ[i])/CUBIC(sigma[i]/2.))/rij;
       ftxij = -urA/SQR(patchsize)*2*(1.-thetai)*dcosdx;
       ftyij = -urA/SQR(patchsize)*2*(1.-thetai)*dcosdy;
       ftzij = -urA/SQR(patchsize)*2*(1.-thetai)*dcosdz;
       // dcos/dsj 
       dcosdx = (xij/(sigma[j]/2.)-prodmaxj*(RSX[j][jp0]-RX[j])/CUBIC(sigma[j]/2.))/rij;
       dcosdy = (yij/(sigma[j]/2.)-prodmaxj*(RSY[j][jp0]-RY[j])/CUBIC(sigma[j]/2.))/rij;
       dcosdz = (zij/(sigma[j]/2.)-prodmaxj*(RSZ[j][jp0]-RZ[j])/CUBIC(sigma[j]/2.))/rij;
       ftxji = -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdx;
       ftyji = -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdy;
       ftzji = -urA/SQR(patchsize)*2*(1.-thetaj)*dcosdz;
       // torque on i due to j 
       txij = (RSY[i][ip0]-RY[i])*ftzij - (RSZ[i][ip0]-RZ[i])*ftyij;
       tyij = (RSZ[i][ip0]-RZ[i])*ftxij - (RSX[i][ip0]-RX[i])*ftzij;
       tzij = (RSX[i][ip0]-RX[i])*ftyij - (RSY[i][ip0]-RY[i])*ftxij;
       // torque on j due to i 
       txji = (RSY[j][jp0]-RY[j])*ftzji - (RSZ[j][jp0]-RZ[j])*ftyji;
       tyji = (RSZ[j][jp0]-RZ[j])*ftxji - (RSX[j][jp0]-RX[j])*ftzji;
       tzji = (RSX[j][jp0]-RX[j])*ftyji - (RSY[j][jp0]-RY[j])*ftxji;
       
       torque[i].x += txij; 
       torque[i].y += tyij; 
       torque[i].z += tzij; 
       
       torque[j].x += txji; 
       torque[j].y += tyji; 
       torque[j].z += tzji; 
     }//end if rij < rc
    }//end loop jlist
   }//endif list is not empty
  }// end loop i
  } // endif VerletIndexCheck = 0 
  
return;
} // end function Force()

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
