/***************************************************
* return distance r, u(r), du/dr, u_tail, w_tail
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

double PBC(double xx,double box) // put back to box
{
  if(xx < 0.0)
  {
  do
   xx += box;
  while(xx < 0.0);
  }
  else if(xx >= box)
  {
  do
   xx -= box; 
  while(xx >= box);
  }

 return xx;
}

double Distance(int i,int j) // return rij^2
{
  double r2,dx,dy,dz;
  
  // miminum image convention under periodic boundary condition
  dx = RX[i]-RX[j];
  if(MinimumImage[0]==0)
  dx = dx - LX*round(dx/LX);

  dy = RY[i]-RY[j];
  if(MinimumImage[1]==0)
  dy = dy - LY*round(dy/LY);

  dz = RZ[i]-RZ[j];
  if(MinimumImage[2]==0)
  dz = dz - LZ*round(dz/LZ);
 
  r2 = SQR(dx)+SQR(dy)+SQR(dz); // reduced unit
  //r = sqrt(SQR(dx)+SQR(dy)+SQR(dz));

  return r2;
}

double SigmaIJ(int iID,int jID)
{
 double sigmaij;


 if(iID == jID)
 {
  if(iID == 1) // A-A
  sigmaij = sigmaA;
  else
  sigmaij = sigmaB;
 }
 else
  sigmaij = sigmaAB;

 return sigmaij;
}

// A: ID = 1 B:ID = 2
double Potential(int iID,int jID,double r,double thetai,double thetaj) // u(rij)
{
  double u;
  double ang;
 
  if(PotentialType == 0) // L-J
  {
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r < rcA)
      u = 4.0*epsilonA*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.));
      else
      u = 0.;
     else // B-B
      if(r < rcB)
      u = 4.0*epsilonB*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.));
      else
      u = 0.;
    }
    else // A-B
      if(r < rcAB)
      u = 4.0*epsilonAB*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.));
      else
      u = 0.;
  } //endif L-J
  else if(PotentialType == 1) // shifted L-J
  {
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r < rcA)
      u = 4.0*epsilonA*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.)) - ucA;
      else
      u = 0.;
     else // B-B
      if(r < rcB)
      u = 4.0*epsilonB*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.)) - ucB;
      else
      u = 0.;
    }
    else // A-B
      if(r < rcAB)
      u = 4.0*epsilonAB*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.)) - ucAB;
      else
      u = 0.;
  } //endif shifted L-J
  else if(PotentialType == 2) // shifted-force L-J  1/r^12 - 1/r^6
  {
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r < rcA)
      u = 4.0*epsilonA*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.)) - ucA - ducdrA*(r-rcA);
      else
      u = 0.;
     else // B-B
      if(r < rcB)
      u = 4.0*epsilonB*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.)) - ucB - ducdrB*(r-rcB);
      else
      u = 0.;
    }
    else // A-B
      if(r < rcAB)
      u = 4.0*epsilonAB*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.)) - ucAB - ducdrAB*(r-rcAB);
      else
      u = 0.;
  } // endif shifted-force L-J
  else if(PotentialType == 3) // patchy interaction shifted L-J
  {
    //ang=exp(-SQR(thetai/patchsize))*exp(-SQR(thetaj/patchsize));
    ang=exp(-SQR((1.0-thetai)/patchsize))*exp(-SQR((1.0-thetaj)/patchsize));
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r <= rminA)
      u = 4.0*epsilonA*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.)) - ucA - uminA + uminA*ang;
      else if (r > rminA && r <= rcA)
      u = (4.0*epsilonA*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.)) - ucA )*ang;
      else
      u = 0.;
     else // B-B
      if(r <= rminB)
      u = 4.0*epsilonB*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.)) - ucB - uminB + uminB*ang;
      else if (r > rminB && r <= rcB)
      u = (4.0*epsilonB*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.)) - ucB )*ang;
      else
      u = 0.;
    }
    else // A-B
      if(r <= rminAB)
      u = 4.0*epsilonAB*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.)) - ucAB - uminAB + uminAB*ang;
      else if (r > rminAB && r <= rcAB)
      u = (4.0*epsilonAB*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.)) - ucAB )*ang;
      else
      u = 0.;
  }

  return u;
}

double Potential_Attract(int iID,int jID,double r,double thetai,double thetaj) // u(rij)
{
 double u;
 double ang;
    //ang=exp(-SQR(thetai/patchsize))*exp(-SQR(thetaj/patchsize));
    ang=exp(-SQR((1.0-thetai)/patchsize))*exp(-SQR((1.0-thetaj)/patchsize));
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r <= rminA)
      u = uminA*ang;
      else if (r > rminA && r <= rcA)
      u = (4.0*epsilonA*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.)) - ucA )*ang;
      else
      u = 0.;
     else // B-B
      if(r <= rminB)
      u = uminB*ang;
      else if (r > rminB && r <= rcB)
      u = (4.0*epsilonB*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.)) - ucB )*ang;
      else
      u = 0.;
    }
    else // A-B
      if(r <= rminAB)
      u = uminAB*ang;
      else if (r > rminAB && r <= rcAB)
      u = (4.0*epsilonAB*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.)) - ucAB )*ang;
      else
      u = 0.;

 return u;
}


double Potential_dr(int iID,int jID,double r,double thetai,double thetaj) // du(rij)/dr
{
  double dudr;
  double ang;
 
  if(PotentialType == 0) // L-J
  {
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r < rcA)
      dudr = -48*epsilonA/r*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.)/2.);
      else
      dudr = 0.;
     else // B-B
      if(r < rcB)
      dudr = -48*epsilonB/r*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.)/2.);
      else
      dudr = 0.;
    }
    else // A-B
      if(r < rcAB)
      dudr = -48*epsilonAB/r*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.)/2.);
      else
      dudr = 0.;
  } //endif L-J
  else if(PotentialType == 1) // shifted L-J
  {
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r < rcA)
      dudr = -48*epsilonA/r*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.)/2.);
      else
      dudr = 0.;
     else // B-B
      if(r < rcB)
      dudr = -48*epsilonB/r*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.)/2.);
      else
      dudr = 0.;
    }
    else // A-B
      if(r < rcAB)
      dudr = -48*epsilonAB/r*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.)/2.);
      else
      dudr = 0.;
  } //endif shifted L-J
  else if(PotentialType == 2) // shifted-force L-J  1/r^12 - 1/r^6
  {
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r < rcA)
      dudr = -48*epsilonA/r*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.)/2.) - ducdrA;
      else
      dudr = 0.;
     else // B-B
      if(r < rcB)
      dudr = -48*epsilonB/r*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.)/2.) - ducdrB;
      else
      dudr = 0.;
    }
    else // A-B
      if(r < rcAB)
      dudr = -48*epsilonAB/r*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.)/2.) - ducdrAB;
      else
      dudr = 0.;
  } // endif shifted-force L-J
  else if(PotentialType == 3) // patchy interaction shifted L-J
  {
    //ang=exp(-SQR(thetai/patchsize))*exp(-SQR(thetaj/patchsize));
    ang=exp(-SQR((1.0-thetai)/patchsize))*exp(-SQR((1.0-thetaj)/patchsize));
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r <= rminA)
      dudr = -48*epsilonA/r*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.)/2.);
      else if (r > rminA && r <= rcA)
      dudr = -48*epsilonA/r*(pow(sigmaA/r,12.)-pow(sigmaA/r,6.)/2.)*ang;
      else
      dudr = 0.;
     else // B-B
      if(r <= rminB)
      dudr = -48*epsilonB/r*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.)/2.);
      else if (r > rminB && r <= rcB)
      dudr = -48*epsilonB/r*(pow(sigmaB/r,12.)-pow(sigmaB/r,6.)/2.)*ang;
      else
      dudr = 0.;
    }
    else // A-B
      if(r <= rminAB)
      dudr = -48*epsilonAB/r*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.)/2.);
      else if (r > rminAB && r <= rcAB)
      dudr = -48*epsilonAB/r*(pow(sigmaAB/r,12.)-pow(sigmaAB/r,6.)/2.)*ang;
      else
      dudr = 0.;
  }// endif patchy

  return dudr;
}

double Potential_dr2(int iID,int jID,double r) // d2u(rij)/dr2
{
 double dudr2; // d2u/dr2
  
 if(PotentialType == 0) // L-J
 {
   if(iID == jID)
   {
    if(iID == 1) // A-A
     if(r < rcA)
     dudr2 = 4.*epsilonA/SQR(r)*(12.*13.*pow(sigmaA/r,12.)-6.*7.*pow(sigmaA/r,6.));
     else
     dudr2 = 0.;
    else // B-B
     if(r < rcB)
     dudr2 = 4.*epsilonB/SQR(r)*(12.*13.*pow(sigmaB/r,12.)-6.*7.*pow(sigmaB/r,6.));
     else
     dudr2 = 0.;
   }
   else // A-B
     if(r < rcAB)
     dudr2 = 4.*epsilonAB/SQR(r)*(12.*13.*pow(sigmaAB/r,12.)-6.*7.*pow(sigmaAB/r,6.));
     else
     dudr2 = 0.;
 } //endif L-J
 else if(PotentialType == 1) // shifted L-J
 {
   if(iID == jID)
   {
    if(iID == 1) // A-A
     if(r < rcA)
     dudr2 = 4.*epsilonA/SQR(r)*(12.*13.*pow(sigmaA/r,12.)-6.*7.*pow(sigmaA/r,6.));
     else
     dudr2 = 0.;
    else // B-B
     if(r < rcB)
     dudr2 = 4.*epsilonB/SQR(r)*(12.*13.*pow(sigmaB/r,12.)-6.*7.*pow(sigmaB/r,6.));
     else
     dudr2 = 0.;
   }
   else // A-B
     if(r < rcAB)
     dudr2 = 4.*epsilonAB/SQR(r)*(12.*13.*pow(sigmaAB/r,12.)-6.*7.*pow(sigmaAB/r,6.));
     else
     dudr2 = 0.;
 }// endif shifted L-J
 else if(PotentialType == 2) // shifted force  L-J
 {
   if(iID == jID)
   {
    if(iID == 1) // A-A
     if(r < rcA)
     dudr2 = 4.*epsilonA/SQR(r)*(12.*13.*pow(sigmaA/r,12.)-6.*7.*pow(sigmaA/r,6.));
     else
     dudr2 = 0.;
    else // B-B
     if(r < rcB)
     dudr2 = 4.*epsilonB/SQR(r)*(12.*13.*pow(sigmaB/r,12.)-6.*7.*pow(sigmaB/r,6.));
     else
     dudr2 = 0.;
   }
   else // A-B
     if(r < rcAB)
     dudr2 = 4.*epsilonAB/SQR(r)*(12.*13.*pow(sigmaAB/r,12.)-6.*7.*pow(sigmaAB/r,6.));
     else
     dudr2 = 0.;
 } // endif shifted force LJ

 return dudr2;
}
