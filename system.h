/******************************************
 * headfile containing most global variables
 ******************************************/
#include <stdio.h>

#define MAX_NumberOfParticles 20000
#define MAX_NumberOfPatches 20 // Z=6 for SC, 8 for BCC, 12 for FCC
#define MAX_NumberOfNeighbors 500 // max number of neighboring particles in the Verlet list
#define MAX_drBins 10000 // for g(r)
//#define MAX_Dimension 3
#define SQR(x) ((x)*(x))
#define CUBIC(x) ((x)*(x)*(x))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define MIN(x,y) ((x)<(y) ? (x) : (y))
typedef struct
{
	double x;
	double y;
	double z;
} VECTOR;


extern int VerletUpdate;

extern int JobIndex;
extern double randomseed;

extern int PotentialType; // 
extern int EnsembleType; // 
extern int rInitialType,vInitialType;
extern int PackType; // SC, BCC or FCC

extern int NumberOfParticles; //N
extern int Nftrans,Nfrot,Nf; // translational and rotational degrees of freedom; -3 if velocity center of mass fixed
extern int NumberOfLatticeSites;
extern int NA,NB; // N = NA+NB
extern double fA,fB;// fraction
extern int ZA,ZB; // coordination number or number of patches per A(or B) particle: SC(6),BCC(8),FCC(12)

extern int NumberOfSteps;
extern int NumberOfInitialSteps;
extern int SampleMultiplier,MovieMultiplier;

extern double rc; // cutoff distance of pair potential u(r)=0 for r > rcutoff
extern double rcA,rcB,rcAB;
extern double ucA,ucB,ucAB; // u(r=rc) at cutoff
extern double ducdrA,ducdrB,ducdrAB; // du(r=rc)/dr at cutoff
extern double rminA,rminB,rminAB; // at which f(r)=0
extern double uminA,uminB,uminAB; // at which f(r)=0
extern double rv; // radius of Verlet neighbor shell

/**************************************************/
extern double dt; // time increment
/**************************************************/

/**************************************************/
extern double kB; // Boltzmann constant, set to 1
extern double T,T0,Tf; //external, initial, final temperature
extern double dT; // temperature increment
extern double Ttrans,Trot,Tinstant; //translational,rotational, instant temperature
extern double Ktrans,Krot,Kinstant,Upotential,Virial,HyperVirial,PdotF,RdotPdotX,PdotP,Enthalpy,TdotL,WdotL;
extern double beta; // 1/(kB*T)
extern double rho,rhoA,rhoB; // number density rho = N/V
extern double packingfraction; // phi = pi/6 rho in 3D
extern double V;  // Volume
extern double L; // min LXLYLZ
extern double LX,LX1,LX2,LX3; //simulation box length V = LX LY LZ
extern double LY,LY1,LY2,LY3; //simulation box length V = LX LY LZ
extern double LZ,LZ1,LZ2,LZ3; //simulation box length V = LX LY LZ
extern double AspectRatioX,AspectRatioY,AspectRatioZ;
extern double latticeconst; // lattice constant
extern double P; // pressure
extern double Pinstant;
extern int BCType[3]; // boundary condition
extern int MinimumImage[3]; // Minimum Image distance switch
extern double rate; // cooling rate
/**************************************************/

extern double sigma[MAX_NumberOfParticles]; // hardcore diameter  sigma_i
extern double epsilon[MAX_NumberOfParticles]; // attraction well depth epsilon_i
extern double mass[MAX_NumberOfParticles]; // particle mass  m_i
extern int  zcoordinate[MAX_NumberOfParticles]; // particle coordinate number (patch number) 
extern double IXX[MAX_NumberOfParticles],IYY[MAX_NumberOfParticles],IZZ[MAX_NumberOfParticles]; // principal inertias
extern double identity[MAX_NumberOfParticles]; // particle identity 1:A 2:B
extern double epsilonA,epsilonB,epsilonAB; // attraction well depth for binary mixture
extern double sigmaA,sigmaB,sigmaAB;
extern double massA,massB;
extern double g[MAX_drBins],gAA[MAX_drBins],gBB[MAX_drBins],gAB[MAX_drBins];
extern double dradial;
extern int drBins;
extern double patchsize; // exp(-theta^2/patchsize^2)

//extern VECTOR position[MAX_NumberOfParticles]; // t
//extern VECTOR position0[MAX_NumberOfParticles]; // position since last Verlet list update
//extern VECTOR position_old[MAX_NumberOfParticles]; // t-dt
//extern VECTOR position_new[MAX_NumberOfParticles]; // t+dt
//extern VECTOR velocity[MAX_NumberOfParticles];
//extern VECTOR velocity_old[MAX_NumberOfParticles];
//extern VECTOR velocity_new[MAX_NumberOfParticles];
extern VECTOR force[MAX_NumberOfParticles]; // total force on particle i
extern VECTOR torque[MAX_NumberOfParticles]; // total torque on particle i

/************ Gear-Predictor-Corrector**************************/
// position r
extern double RX[MAX_NumberOfParticles],RY[MAX_NumberOfParticles],RZ[MAX_NumberOfParticles]; // position not put particle back to box
extern double RX1[MAX_NumberOfParticles],RY1[MAX_NumberOfParticles],RZ1[MAX_NumberOfParticles]; // r'=v
extern double RX2[MAX_NumberOfParticles],RY2[MAX_NumberOfParticles],RZ2[MAX_NumberOfParticles]; // r''=a
extern double RX3[MAX_NumberOfParticles],RY3[MAX_NumberOfParticles],RZ3[MAX_NumberOfParticles]; // r'''
//extern double RX4[MAX_NumberOfParticles],RY4[MAX_NumberOfParticles],RZ4[MAX_NumberOfParticles]; // r''''

//momentum p
extern double PX[MAX_NumberOfParticles],PY[MAX_NumberOfParticles],PZ[MAX_NumberOfParticles]; // p=mv
extern double PX1[MAX_NumberOfParticles],PY1[MAX_NumberOfParticles],PZ1[MAX_NumberOfParticles]; // p'
extern double PX2[MAX_NumberOfParticles],PY2[MAX_NumberOfParticles],PZ2[MAX_NumberOfParticles]; // p''
extern double PX3[MAX_NumberOfParticles],PY3[MAX_NumberOfParticles],PZ3[MAX_NumberOfParticles]; // p'''
//extern double PX4[MAX_NumberOfParticles],PY4[MAX_NumberOfParticles],PZ4[MAX_NumberOfParticles]; // p''''

//quaternion Q = (QW,QX,QY,QZ)
extern double QW[MAX_NumberOfParticles],QX[MAX_NumberOfParticles],QY[MAX_NumberOfParticles],QZ[MAX_NumberOfParticles]; // q
extern double QW1[MAX_NumberOfParticles],QX1[MAX_NumberOfParticles],QY1[MAX_NumberOfParticles],QZ1[MAX_NumberOfParticles]; // q'
extern double QW2[MAX_NumberOfParticles],QX2[MAX_NumberOfParticles],QY2[MAX_NumberOfParticles],QZ2[MAX_NumberOfParticles]; // q''
extern double QW3[MAX_NumberOfParticles],QX3[MAX_NumberOfParticles],QY3[MAX_NumberOfParticles],QZ3[MAX_NumberOfParticles]; // q'''
extern double QW4[MAX_NumberOfParticles],QX4[MAX_NumberOfParticles],QY4[MAX_NumberOfParticles],QZ4[MAX_NumberOfParticles]; // q''''

// angular velocity W: omega
extern double WX[MAX_NumberOfParticles],WY[MAX_NumberOfParticles],WZ[MAX_NumberOfParticles]; // W
extern double WX1[MAX_NumberOfParticles],WY1[MAX_NumberOfParticles],WZ1[MAX_NumberOfParticles]; // W' 
extern double WX2[MAX_NumberOfParticles],WY2[MAX_NumberOfParticles],WZ2[MAX_NumberOfParticles]; // W''
extern double WX3[MAX_NumberOfParticles],WY3[MAX_NumberOfParticles],WZ3[MAX_NumberOfParticles]; // W'''
extern double WX4[MAX_NumberOfParticles],WY4[MAX_NumberOfParticles],WZ4[MAX_NumberOfParticles]; // W''''

extern double RSX[MAX_NumberOfParticles][MAX_NumberOfPatches],RSY[MAX_NumberOfParticles][MAX_NumberOfPatches],RSZ[MAX_NumberOfParticles][MAX_NumberOfPatches]; // absolute position of a patch
extern double dRSX[MAX_NumberOfParticles][MAX_NumberOfPatches],dRSY[MAX_NumberOfParticles][MAX_NumberOfPatches],dRSZ[MAX_NumberOfParticles][MAX_NumberOfPatches]; // relative position of a patch
extern double FSX[MAX_NumberOfParticles][MAX_NumberOfPatches],FSY[MAX_NumberOfParticles][MAX_NumberOfPatches],FSZ[MAX_NumberOfParticles][MAX_NumberOfPatches]; // total force on a patch

extern double xi_trans,xi_rot,chi,xi,chipxi; // chi and chi+xi
/******************************************************************/

/************************* Verlet list *********************************/
extern double RX0[MAX_NumberOfParticles],RY0[MAX_NumberOfParticles],RZ0[MAX_NumberOfParticles]; // store old position for updating Verlet list
extern int VerletCheckIndex; // 0: no need to update  1: need to update
extern int point[MAX_NumberOfParticles]; // point to the index number of where the verlet list of particle i starts from
extern int nlist; // verlet list index
extern int VerletList[MAX_NumberOfParticles*MAX_NumberOfNeighbors]; // verlet list stores all the neighbors of all particles
/******************************************************************/

void ReadInput(void);
void Initialization(int job);
double BoxMuller(double mm, double ss);
VECTOR RandomSphere(double norm);
void Writemovie(FILE *FilePtr);
void COMcheck(FILE *fp);


double Distance(int i,int j); // return rij^2
double SigmaIJ(int iID,int jID);
double Potential(int iID,int jID,double r,double thetai,double thetaj); // u(rij)
double Potential_dr(int iID,int jID,double r,double thetai,double thetaj); // du(rij)/dr
double Potential_Attract(int iID,int jID,double r,double thetai,double thetaj); // uA(rij)
double Potential_dr2(int iID,int jID,double r); // d2u(rij)/dr2
double PBC(double xx,double box);

void PredictNVE(double DT);
void CorrectNVE(double DT);
void PredictNVE2(double DT);
void CorrectNVE2(double DT);
void PredictNVT(double DT);
void CorrectNVT(double DT,double XI_TRANS,double XI_ROT);
void MolAtm(void);
void AtmMol(void);
void Force(void); 
void Kinetic(void);
void RadialDis(void);
void ScaleVelocity(void); //  need to use after Kinetic()
void VerletCheck(void);
