#include "system.h"

int VerletUpdate;

int JobIndex;
double randomseed;

int PotentialType; // 
int EnsembleType; // 
int rInitialType,vInitialType;
int PackType; // SC, BCC or FCC

int NumberOfParticles; //N
int Nftrans,Nfrot,Nf; // translational and rotational degrees of freedom; -3 if velocity center of mass fixed
int NumberOfLatticeSites;
int NA,NB; // N = NA+NB
double fA,fB;// fraction
int ZA,ZB; // coordination number or number of patches per A(or B) particle: SC(6),BCC(8),FCC(12)

int NumberOfSteps;
int NumberOfInitialSteps;
int SampleMultiplier,MovieMultiplier;

double rc; // cutoff distance of pair potential u(r)=0 for r > rcutoff
double rcA,rcB,rcAB;
double ucA,ucB,ucAB; // u(r=rc) at cutoff
double ducdrA,ducdrB,ducdrAB; // du(r=rc)/dr at cutoff
double rminA,rminB,rminAB; // at which f(r)=0
double uminA,uminB,uminAB; // at which f(r)=0
double rv; // radius of Verlet neighbor shell

/**************************************************/
double dt; // time increment
/**************************************************/

/**************************************************/
double kB; // Boltzmann constant, set to 1
double T,T0,Tf; //external, initial, final temperature
double dT; // temperature increment
double Ttrans,Trot,Tinstant; //translational,rotational, instant temperature
double Ktrans,Krot,Kinstant,Upotential,Virial,HyperVirial,PdotF,RdotPdotX,PdotP,Enthalpy,TdotL,WdotL;
double beta; // 1/(kB*T)
double rho,rhoA,rhoB; // number density rho = N/V
double packingfraction; // phi = pi/6 rho in 3D
double V;  // Volume
double L; // min LXLYLZ
double LX,LX1,LX2,LX3; //simulation box length V = LX LY LZ
double LY,LY1,LY2,LY3; //simulation box length V = LX LY LZ
double LZ,LZ1,LZ2,LZ3; //simulation box length V = LX LY LZ
double AspectRatioX,AspectRatioY,AspectRatioZ;
double latticeconst; // lattice constant
double P; // pressure
double Pinstant;
int BCType[3]; // boundary condition
int MinimumImage[3]; // Minimum Image distance switch
double rate; // cooling rate
/**************************************************/

double sigma[MAX_NumberOfParticles]; // hardcore diameter  sigma_i
double epsilon[MAX_NumberOfParticles]; // attraction well depth epsilon_i
double mass[MAX_NumberOfParticles]; // particle mass  m_i
int  zcoordinate[MAX_NumberOfParticles]; // particle coordinate number (patch number) 
double IXX[MAX_NumberOfParticles],IYY[MAX_NumberOfParticles],IZZ[MAX_NumberOfParticles]; // principal inertias
double identity[MAX_NumberOfParticles]; // particle identity 1:A 2:B
double epsilonA,epsilonB,epsilonAB; // attraction well depth for binary mixture
double sigmaA,sigmaB,sigmaAB;
double massA,massB;
double g[MAX_drBins],gAA[MAX_drBins],gBB[MAX_drBins],gAB[MAX_drBins];
double dradial;
int drBins;
double patchsize; // exp(-theta^2/patchsize^2)

//VECTOR position[MAX_NumberOfParticles]; // t
//VECTOR position0[MAX_NumberOfParticles]; // position since last Verlet list update
//VECTOR position_old[MAX_NumberOfParticles]; // t-dt
//VECTOR position_new[MAX_NumberOfParticles]; // t+dt
//VECTOR velocity[MAX_NumberOfParticles];
//VECTOR velocity_old[MAX_NumberOfParticles];
//VECTOR velocity_new[MAX_NumberOfParticles];
VECTOR force[MAX_NumberOfParticles]; // total force on particle i
VECTOR torque[MAX_NumberOfParticles]; // total torque on particle i

/************ Gear-Predictor-Corrector**************************/
// position r
double RX[MAX_NumberOfParticles],RY[MAX_NumberOfParticles],RZ[MAX_NumberOfParticles]; // position not put particle back to box
double RX1[MAX_NumberOfParticles],RY1[MAX_NumberOfParticles],RZ1[MAX_NumberOfParticles]; // r'=v
double RX2[MAX_NumberOfParticles],RY2[MAX_NumberOfParticles],RZ2[MAX_NumberOfParticles]; // r''=a
double RX3[MAX_NumberOfParticles],RY3[MAX_NumberOfParticles],RZ3[MAX_NumberOfParticles]; // r'''
//double RX4[MAX_NumberOfParticles],RY4[MAX_NumberOfParticles],RZ4[MAX_NumberOfParticles]; // r''''

//momentum p
double PX[MAX_NumberOfParticles],PY[MAX_NumberOfParticles],PZ[MAX_NumberOfParticles]; // p=mv
double PX1[MAX_NumberOfParticles],PY1[MAX_NumberOfParticles],PZ1[MAX_NumberOfParticles]; // p'
double PX2[MAX_NumberOfParticles],PY2[MAX_NumberOfParticles],PZ2[MAX_NumberOfParticles]; // p''
double PX3[MAX_NumberOfParticles],PY3[MAX_NumberOfParticles],PZ3[MAX_NumberOfParticles]; // p'''
//double PX4[MAX_NumberOfParticles],PY4[MAX_NumberOfParticles],PZ4[MAX_NumberOfParticles]; // p''''

// quaternion Q = (QW,QX,QY,QZ)
double QW[MAX_NumberOfParticles],QX[MAX_NumberOfParticles],QY[MAX_NumberOfParticles],QZ[MAX_NumberOfParticles]; // q
double QW1[MAX_NumberOfParticles],QX1[MAX_NumberOfParticles],QY1[MAX_NumberOfParticles],QZ1[MAX_NumberOfParticles]; // q'
double QW2[MAX_NumberOfParticles],QX2[MAX_NumberOfParticles],QY2[MAX_NumberOfParticles],QZ2[MAX_NumberOfParticles]; // q''
double QW3[MAX_NumberOfParticles],QX3[MAX_NumberOfParticles],QY3[MAX_NumberOfParticles],QZ3[MAX_NumberOfParticles]; // q'''
double QW4[MAX_NumberOfParticles],QX4[MAX_NumberOfParticles],QY4[MAX_NumberOfParticles],QZ4[MAX_NumberOfParticles]; // q''''

//angular velocity W: omega
double WX[MAX_NumberOfParticles],WY[MAX_NumberOfParticles],WZ[MAX_NumberOfParticles]; // W
double WX1[MAX_NumberOfParticles],WY1[MAX_NumberOfParticles],WZ1[MAX_NumberOfParticles]; // W' 
double WX2[MAX_NumberOfParticles],WY2[MAX_NumberOfParticles],WZ2[MAX_NumberOfParticles]; // W''
double WX3[MAX_NumberOfParticles],WY3[MAX_NumberOfParticles],WZ3[MAX_NumberOfParticles]; // W'''
double WX4[MAX_NumberOfParticles],WY4[MAX_NumberOfParticles],WZ4[MAX_NumberOfParticles]; // W''''

double RSX[MAX_NumberOfParticles][MAX_NumberOfPatches],RSY[MAX_NumberOfParticles][MAX_NumberOfPatches],RSZ[MAX_NumberOfParticles][MAX_NumberOfPatches]; // absolute position of a patch
double dRSX[MAX_NumberOfParticles][MAX_NumberOfPatches],dRSY[MAX_NumberOfParticles][MAX_NumberOfPatches],dRSZ[MAX_NumberOfParticles][MAX_NumberOfPatches]; // relative position of a patch
double FSX[MAX_NumberOfParticles][MAX_NumberOfPatches],FSY[MAX_NumberOfParticles][MAX_NumberOfPatches],FSZ[MAX_NumberOfParticles][MAX_NumberOfPatches]; // total force on a patch

double xi_trans,xi_rot,chi,xi,chipxi; // chi and chi+xi
/******************************************************************/

/************************* Verlet list *********************************/
double RX0[MAX_NumberOfParticles],RY0[MAX_NumberOfParticles],RZ0[MAX_NumberOfParticles]; // store old position for updating Verlet list
int VerletCheckIndex; // 0: no need to update  1: need to update
int point[MAX_NumberOfParticles]; // point to the index number of where the verlet list of particle i starts from
int nlist; // verlet list index
int VerletList[MAX_NumberOfParticles*MAX_NumberOfNeighbors]; // verlet list stores all the neighbors of all particles
/******************************************************************/

