/***************************************************
* 3D Molecular Dynamics (MD) simulation of
* binary mixtures: particle A and B
* with continuous potential
* patchy particles
* Lennard-Jones (LJ) 	            0
* shifted Lennard-Jones (LJ)        1
* shifted-force Lennard-Jones (LJ)  2
* Kai Zhang, Yale University, 2014
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "system.h"
#include "ran_uniform.h"

//int main(void)
int main(int argc, char *argv[])
{
 int step;
 int i,j;
 char filename[20];
 FILE *fp;
 FILE *fpmovie;
 FILE *fpvelocity;
 FILE *fpsample;
 FILE *fpr; // position
 FILE *fpv; // velocity
 FILE *fpf; // force

 double Sum_K,Sum_U,Sum_E,Sum_W,Sum_T,Sum_V,Sum_rho,Sum_P,Sum_H;
 double Count;

 sscanf(argv[1],"%d",&JobIndex);
 
 printf("**************** 3D Molecular Dynamics simulation ********************");
 printf("\n");

 ReadInput();

 if(randomseed == 0.)  randomseed = (double)(JobIndex);
 printf("randomseed = %lf\n",randomseed);
 InitializeRandomNumberGenerator(randomseed); //time(0L)

 Initialization(JobIndex); // particle positions and velocities
 fflush(stdout);
 
 Sum_K = 0.;
 Sum_U = 0.;
 Sum_E = 0.;
 Sum_W = 0.;
 Sum_T = 0.;
 Sum_V = 0.;
 Sum_P = 0.;
 Sum_H = 0.;
 Sum_rho = 0.;
 Count = 0.;
 VerletUpdate = 0;

 for(i=0;i<drBins;i++)
 {
  g[i] = 0.;
  gAA[i] = 0.;
  gBB[i] = 0.;
  gAB[i] = 0.;
 }

 Kinetic();
 MolAtm();
 VerletCheckIndex = 1; // make list in the beginning
 Force();

 
 sprintf(filename,"movie_%d.xyz",JobIndex);
 fpmovie=fopen(filename,"w");
 //sprintf(filename,"sample_%d.dat",JobIndex);
// fpsample=fopen(filename,"w");

 for(step=0;step<NumberOfSteps;step++)
 {
  if(step%MovieMultiplier==0) Writemovie(fpmovie);

  if(step%SampleMultiplier==0)
  {
  printf("step = %d...\n",step);

 sprintf(filename,"sample_%d.dat",JobIndex);
 fpsample=fopen(filename,"a+");
  fprintf(fpsample,"t = %lf K = %lf U = %lf E = %lf Ttrans = %lf Trot = %lf T = %lf\n", \
 step*dt,Kinstant/NumberOfParticles,Upotential/NumberOfParticles,(Kinstant+Upotential)/NumberOfParticles,Ttrans,Trot,Tinstant);
 fclose(fpsample);

  /********************/
  if(step>=NumberOfInitialSteps && rate < 0.0000001) 
  {
   Sum_K += Kinstant;
   Sum_U += Upotential;
   Sum_E += Kinstant + Upotential;
   Sum_T += Tinstant;
   Sum_V += V;
   Sum_rho += rho;
   Count += 1.;
   RadialDis(); // g(r)
  }//endif after equilibration
 /********************/
  }//endif every ? steps 

  if(EnsembleType == 0) // NVE
  {
   PredictNVE(dt);
   MolAtm();
   VerletCheck();
   Force();   
   Kinetic();
   AtmMol();
   CorrectNVE(dt);
  }// end NVE
  else if(EnsembleType == 1) // NVK isokinetic scale
  {
   PredictNVE(dt);
   MolAtm();
   VerletCheck();
   Force();   
   Kinetic();
   AtmMol();
   CorrectNVE(dt);
   if(step>=NumberOfInitialSteps) ScaleVelocity();
  }// end NVK
  else if(EnsembleType == 2)// NVT with constraint method
  {
   if(step < NumberOfInitialSteps/2)
   {
   PredictNVE2(dt);
   MolAtm();
   VerletCheck();
   Force();
   Kinetic();
   AtmMol();
   CorrectNVE2(dt);
   }
   else
   {
   PredictNVT(dt);
   MolAtm();
   VerletCheck();
   Force();
   Kinetic();
   AtmMol();
  // Pinstant = rho*Tinstant + Virial/V;
   xi_trans = PdotF/PdotP;
   xi_rot = TdotL/WdotL;
   CorrectNVT(dt,xi_trans,xi_rot);
   if(step%SampleMultiplier==0) ScaleVelocity(); //avoid long time drifting
   }
  }// endif NVT

  
  if(step>=NumberOfInitialSteps)
   T = T0*exp(-rate*(step-NumberOfInitialSteps)*dt);
   //T = T0*exp(-rate*step*dt);

  if(T <= Tf)
  {
   printf("finial temperature %lf is reached in %d step\n",T,step);
   Sum_K += Kinstant;
   Sum_U += Upotential;
   Sum_E += Kinstant + Upotential;
   Sum_T += Tinstant;
   Sum_V += V;
   Sum_rho += rho;
   Count += 1.;
   RadialDis(); // g(r)
   break;
  }
 fflush(stdout);
 }// end MD loop

 //fclose(fpsample);
 
  Writemovie(fpmovie);
  fclose(fpmovie);


/*******************/
 Sum_K /= Count;
 Sum_U /= Count;
 Sum_E /= Count;
 Sum_H /= Count;
 Sum_W /= Count;
 Sum_T /= Count;
 Sum_V /= Count;
 Sum_rho /= Count;
 Sum_P /= Count;
/*******************/
 
printf("T = %lf\trho = %lf\tP = %lf\t<K>/N = %lf\t<U>/N = %lf\t<E>/N = %lf\t<H>/N = %lf\t<T> = %lf\t<V> = %lf\t<P> = %lf <rho> = %lf\n", \
 T, rho,P,Sum_K/NumberOfParticles,Sum_U/NumberOfParticles,Sum_E/NumberOfParticles,Sum_H/NumberOfParticles,Sum_T,Sum_V,Sum_P,Sum_rho);

 sprintf(filename,"gr_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  g[i] /= Count;
  g[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial)*Sum_rho;
  fprintf(fp,"r = %lf\tg(r) = %lf\n",(i+0.5)*dradial,g[i]/NumberOfParticles);
 }
 fclose(fp);
 
 sprintf(filename,"grAA_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  gAA[i] /= Count;
  gAA[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial)*Sum_rho*fA;
  fprintf(fp,"r = %lf\tgAA(r) = %lf\n",(i+0.5)*dradial,gAA[i]/NA);
 }
 fclose(fp);

 sprintf(filename,"grBB_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  gBB[i] /= Count;
  gBB[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial)*Sum_rho*fB;
  fprintf(fp,"r = %lf\tgBB(r) = %lf\n",(i+0.5)*dradial,gBB[i]/NB);
 }
 fclose(fp);

 sprintf(filename,"grAB_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  gAB[i] /= Count;
  gAB[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial);
  fprintf(fp,"r = %lf\tgAB(r) = %lf\n",(i+0.5)*dradial,gAB[i]*Sum_V/(NA*NB)/2.);
 }
 fclose(fp);

  sprintf(filename,"position1_%d.dat",JobIndex);
  fp = fopen(filename,"w");
  for(i=0;i<NumberOfParticles;i++)
  fprintf(fp,"%lf\t%lf\t%lf\n",RX[i],RY[i],RZ[i]);
  fclose(fp);
  
  sprintf(filename,"quaternion1_%d.dat",JobIndex);
  fp = fopen(filename,"w");
  for(i=0;i<NumberOfParticles;i++)
  fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",QW[i],QX[i],QY[i],QZ[i]);
  fclose(fp);
  
  sprintf(filename,"velocity1_%d.dat",JobIndex);
  fp = fopen(filename,"w");
  for(i=0;i<NumberOfParticles;i++)
  fprintf(fp,"%lf\t%lf\t%lf\n",RX1[i],RY1[i],RZ1[i]);
  fclose(fp);
  
  sprintf(filename,"angularvelocity1_%d.dat",JobIndex);
  fp = fopen(filename,"w");
  for(i=0;i<NumberOfParticles;i++)
  fprintf(fp,"%lf\t%lf\t%lf\n",WX[i],WY[i],WZ[i]);
  fclose(fp);

  printf("VerletUpdate = %d\n",VerletUpdate);

 printf("****************************** the end *******************************");
 return 0;
}
