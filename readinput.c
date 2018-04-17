/*********************************************************************
 * input simulation parameters from file "input"
 *********************************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void ReadInput(void)
{
 int i;
 FILE *fp;
 
 fp=fopen("input","r"); 
 fscanf(fp,"%d",&EnsembleType);
 fscanf(fp,"%d %d %d %d",&rInitialType,&PackType,&NumberOfLatticeSites,&vInitialType);
 fscanf(fp,"%d",&PotentialType);
 fscanf(fp,"%d %d %d",&NumberOfParticles,&NA,&NB);
 fscanf(fp,"%lf %lf",&massA,&massB);
 fscanf(fp,"%lf %lf %lf",&sigmaA,&sigmaB,&sigmaAB);
 fscanf(fp,"%lf %lf %lf",&epsilonA,&epsilonB,&epsilonAB);
 fscanf(fp,"%lf %lf %lf",&T0,&Tf,&P);
 fscanf(fp,"%lf",&V);
 fscanf(fp,"%lf %lf",&rc,&rv);
 fscanf(fp,"%d %d",&NumberOfSteps,&NumberOfInitialSteps);
 fscanf(fp,"%d %d",&SampleMultiplier,&MovieMultiplier);
 fscanf(fp,"%lf",&dt);
 fscanf(fp,"%d",&drBins);
 fscanf(fp,"%lf",&randomseed);
 fscanf(fp,"%lf %lf %lf",&AspectRatioX,&AspectRatioY,&AspectRatioZ);
 fscanf(fp,"%d %d %d",&BCType[0],&BCType[1],&BCType[2]);
 fscanf(fp,"%d %d %d",&MinimumImage[0],&MinimumImage[1],&MinimumImage[2]);
 fscanf(fp,"%d %d",&ZA,&ZB);
 fscanf(fp,"%lf",&patchsize);
 fscanf(fp,"%lf",&rate);
 fclose(fp);

 // input parameters are all in reduced units
 // rc and rv are multiplier of sigmaA

 fA = 1.*NA/NumberOfParticles;
 fB = 1.*NB/NumberOfParticles;

 Nftrans = 3*NumberOfParticles-3;
 Nfrot = 3*NumberOfParticles;
 Nf = Nftrans+Nfrot;
 
 printf("Nftrans = %d Nfrot = %d Nf = %d\n",Nftrans,Nfrot,Nf);

// Nf = 3*NumberOfParticles-3;
// Nf = 3*NumberOfParticles;

 //InitializeRandomNumberGenerator(time(0l)); //time(0L)
 //InitializeRandomNumberGenerator(randomseed); //time(0L)


 rcA = rc*sigmaA;
 rcB = rc*sigmaB;
 rcAB = rc*sigmaAB;

 rc = rc*sigmaA;
 rv = rv*sigmaA;


 printf("\n");
 if(PotentialType == 0) printf("Lennard-Jones potential\n");
 if(PotentialType == 1) printf("shifted Lennard-Jones potential u(rc) = 0\n");
 if(PotentialType == 2) printf("shifted-force Lennard-Jones potential f(rc) = 0\n");
 printf("\n");
 if(EnsembleType == 0) printf("NVE ensemble\n");
 if(EnsembleType == 1) printf("isokinetic ensemble by velocity rescaling\n");
 if(EnsembleType == 2) printf("NVT ensemble via constraint\n");
 if(EnsembleType == 3) printf("NPH ensemble via constraint\n");
 if(EnsembleType == 4) printf("NPT ensemble via constraint\n");
 printf("\n");

 if(PotentialType == 0 || PotentialType == 1 || PotentialType == 2 || PotentialType ==3) // L-J
 {
  ucA = 4.0*epsilonA*(pow(sigmaA/rcA,12.)-pow(sigmaA/rcA,6.));
  ucB = 4.0*epsilonB*(pow(sigmaB/rcB,12.)-pow(sigmaB/rcB,6.));
  ucAB = 4.0*epsilonAB*(pow(sigmaAB/rcAB,12.)-pow(sigmaAB/rcAB,6.));

  ducdrA = -48.*epsilonA/rcA*(pow(sigmaA/rcA,12.) - pow(sigmaA/rcA,6.)/2.);
  ducdrB = -48.*epsilonB/rcB*(pow(sigmaB/rcB,12.) - pow(sigmaB/rcB,6.)/2.);
  ducdrAB = -48.*epsilonAB/rcAB*(pow(sigmaAB/rcAB,12.) - pow(sigmaAB/rcAB,6.)/2.);
 }
 
 rminA = pow(2.,1./6.)*sigmaA;
 rminB = pow(2.,1./6.)*sigmaB;
 rminAB = pow(2.,1./6.)*sigmaAB;

 if(PotentialType == 0)
 {
  uminA = -epsilonA;
  uminB = -epsilonB;
  uminAB = -epsilonAB;
 }
 
 if(PotentialType == 1 || PotentialType == 3)
 {
  uminA = -epsilonA-ucA;
  uminB = -epsilonB-ucB;
  uminAB = -epsilonAB-ucAB;
 }

 printf("rmin = %lf umin = %lf\n",rminA,uminA);

 T = T0; // initial temperature
 kB = 1.;
 beta = 1.0 / T / kB;

 //T = T + dT*JobIndex;
 //rho = rho + dT*JobIndex;

 rho = NumberOfParticles/V*CUBIC(sigmaA);
 rhoA = NA/V*CUBIC(sigmaA);
 rhoB = NB/V*CUBIC(sigmaA);
 //L = pow(V,1.0/3.0); 
 LX = pow(V/(AspectRatioX*AspectRatioY*AspectRatioZ),1.0/3.0)*AspectRatioX; 
 LY = pow(V/(AspectRatioX*AspectRatioY*AspectRatioZ),1.0/3.0)*AspectRatioY; 
 LZ = pow(V/(AspectRatioX*AspectRatioY*AspectRatioZ),1.0/3.0)*AspectRatioZ; 
 L = MIN(MIN(LX,LY),LZ); 
 dradial = L/2./drBins; // PBC is used, L/2

 //packingfraction = M_PI/6.*rho*CUBIC(sigmaA);
 packingfraction = M_PI/6.*rho*(fA*CUBIC(sigmaA)+fB*CUBIC(sigmaB));

 printf("\n");
 printf("Number of MD steps: %d\n",NumberOfSteps);
 printf("Number of equilibrium MD steps: %d\n",NumberOfInitialSteps);
 printf("time increment dt = %lf\n",dt);
 printf("total time = %lf\n",dt*NumberOfSteps);
 printf("Sample multiplier (sampling frequency): %d\n",SampleMultiplier);
 printf("Movie multiplier (draw snapshot frequency): %d\n",MovieMultiplier);
 printf("\n");
 printf("\n");

 printf("T0 = %lf\tTf = %lf\n",T0,Tf);
 printf("T = %lf\tbeta = %lf\n",T,beta);
 printf("rho = %lf\n",rho);
 printf("external pressure P = %lf\n",P);
 printf("packing fraction = %lf\n",packingfraction);
 printf("V = %lf\n",V);
 printf("LX LY LZ = %lf %lf %lf\n",LX,LY,LZ);
 printf("rcutoff = %lf\n",rc);
 printf("rv = %lf\n",rv);
 
 printf("\n");
 printf("N = %d\tNA:NB = %d:%d\tfB = NB/(NA+NB)= %lf\n",NumberOfParticles,NA,NB,fB);
 printf("massA = %lf\tmassB = %lf\n",massA,massB);
 printf("sigmaA = %lf\tsigmaB = %lf\tsigmaAB = %lf\tsizeratio = %lf\n",sigmaA,sigmaB,sigmaAB,sigmaB/sigmaA);
 printf("epsilonA = %lf\tepsilonB = %lf\tepsilonAB = %lf\n",epsilonA,epsilonB,epsilonAB);
 printf("ucA = %lf\tucB = %lf\tucAB = %lf\n",ucA,ucB,ucAB);
 printf("\n");


 fp=fopen("uAA0.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tu(r) = %lf\n",i*0.003,Potential(1,1,i*0.003,1.0,1.0));
 fclose(fp);
 fp=fopen("uAA1.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tu(r) = %lf\n",i*0.003,Potential(1,1,i*0.003,1.0-0.5*patchsize,1.0-0.5*patchsize));
 fclose(fp);
 fp=fopen("uAA2.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tu(r) = %lf\n",i*0.003,Potential(1,1,i*0.003,1.-patchsize,1.-patchsize));
 fclose(fp);
 fp=fopen("uAA3.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tu(r) = %lf\n",i*0.003,Potential(1,1,i*0.003,1.-3.*patchsize,1.-3.*patchsize));
 fclose(fp);
 fp=fopen("uBB.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tu(r) = %lf\n",i*0.003,Potential(2,2,i*0.003,0.,0.));
 fclose(fp);
 fp=fopen("uAB.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tu(r) = %lf\n",i*0.003,Potential(1,2,i*0.003,0.,0.));
 fclose(fp);
 
 fp=fopen("fAA0.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tf(r) = %lf\n",i*0.003,-Potential_dr(1,1,i*0.003,1.,1.));
 fclose(fp);
 fp=fopen("fAA1.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tf(r) = %lf\n",i*0.003,-Potential_dr(1,1,i*0.003,1.-0.5*patchsize,1.0-0.5*patchsize));
 fclose(fp);
 fp=fopen("fAA2.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tf(r) = %lf\n",i*0.003,-Potential_dr(1,1,i*0.003,1.-patchsize,1.-patchsize));
 fclose(fp);
 fp=fopen("fAA3.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tf(r) = %lf\n",i*0.003,-Potential_dr(1,1,i*0.003,1.-3.*patchsize,1.-3.*patchsize));
 fclose(fp);
 fp=fopen("fBB.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tf(r) = %lf\n",i*0.003,-Potential_dr(2,2,i*0.003,0.,0.));
 fclose(fp);
 fp=fopen("fAB.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\tf(r) = %lf\n",i*0.003,-Potential_dr(1,2,i*0.003,0.,0.));
 fclose(fp);
 
 fp=fopen("d2uAA.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\td2u(r) = %lf\n",i*0.003,Potential_dr2(1,1,i*0.003));
 fclose(fp);
 fp=fopen("d2uBB.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\td2u(r) = %lf\n",i*0.003,Potential_dr2(2,2,i*0.003));
 fclose(fp);
 fp=fopen("d2uAB.dat","w"); 
 for(i=0;i<1500;i++)
  fprintf(fp,"r = %lf\td2u(r) = %lf\n",i*0.003,Potential_dr2(1,2,i*0.003));
 fclose(fp);
 
 fp=fopen("vtheta.dat","w"); 
 for(i=0;i<=2000;i++)
  fprintf(fp,"theta = %lf\tv(theta) = %lf\n",(-M_PI+i*(M_PI/1000))/M_PI,exp(-SQR((1.0-cos(-M_PI+i*(M_PI/1000)))/patchsize)));
 fclose(fp);

 return;
}
