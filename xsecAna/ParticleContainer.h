#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include "AnaHelper.h"

namespace xsec_ana {

class ParticleContainer {

private:


public:

ParticleContainer();

virtual ~ParticleContainer();

//-----------------------
//Setter functions
/------------------------

void SetIndex(int);

void SetOrigin(int);
void SetmcPdgCode(int);
void SetmcNuPdgCode(int);
void SetmcParentPdgCode(int);
void SetmcIsNeutrino(bool);
void SetmcIsPrimary(bool);
void SetpfpPdgCode(int);
void SetpfpNuPdgCode(int);
void SetpfpParentPdgCode(int);
void SetpfpIsNeutrino(bool);

void SetmcVtxX  (double);
void SetmcVtxY  (double);
void SetmcVtxZ  (double);
void SetpfpVtxX (double);
void SetpfpVtxY (double);
void SetpfpVtxZ (double);

void SetmcDirX  (double);
void SetmcDirY  (double);
void SetmcDirZ  (double);
void SetpfpDirX (double);
void SetpfpDirY (double);
void SetpfpDirZ (double);

void SetmcTheta (double);
void SetmcPhi (double);
void SetpfpTheta (double);
void SetpfpPhi (double);

void SetmcLength(double);
void SetpfpLength(double);

void SetmcEnergy(double);
void SetmcMomentum(double);
void SetpfpMomentum(double);

void SetCompleteness(double);
void SetPurity(double);

void SetNumMCHits   (int);
void SetNumMCHitsU  (int);
void SetNumMCHitsV  (int);
void SetNumMCHitsW  (int);
void SetNumPFPHits  (int);
void SetNumPFPHitsU (int);
void SetNumPFPHitsV (int);
void SetNumPFPHitsW (int);

void SetNumMatchedHits  (int);
void SetNumMatchedHitsU (int);
void SetNumMatchedHitsV (int);
void SetNumMatchedHitsW (int);

void SetmcOpenAngle(double);
void SetpfpOpenAngle(double);

int origin;//this is the particle origin

int mcPdg;
int mcNuPdg;
int mcNuIndex;
int mcParentPdg;
bool mcIsNeutirno;
bool mcIsPrimary;
int pfpPdg;
int pfpNuPdg;
int pfpNuIndex;
int pfpIndex;
int pfpParentPdg;
bool pfpIsNeutrino;

double mcVtxX;
double mcVtxY;
double mcVtxZ;
double pfpVtxX;
double pfpVtxY;
double pfpVtxZ;

double mcDirX;
double mcDirY;
double mcDirZ;
double pfpDirX;
double pfpDirY;
double pfpDirZ;

double mcTheta;
double mcPhi;
double pfpTheta;
double pfpPhi;

double mcLength;
double pfpLength;

double mcEnergy;
double mcMomentum;
//double pfpEnergy;
double pfpMomentum;

double completeness;
double purity;

int nMCHits;
int nMCHitsU;
int nMCHitsV;
int nMCHitsY;
int nPFPHits;
int nPFPHitsU;
int nPFPHitsV;
int nPFPHitsY;
int nMatchedHits;
int nMatchedHitsU;
int nMatchedHitsV;
int nMatchedHitsY;

double mcOpenAngle;
double pfpOpenAngle;

};//end class

}//end namespace

#endif
