
#include <Pythia8/Pythia.h>
	R__LOAD_LIBRARY(libPHPythia8.so)
	R__LOAD_LIBRARY(EvtGen/lib/libEvtGenExternal.so)
	R__LOAD_LIBRARY(EvtGen/lib/libEvtGen.so)
	




#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TRandom3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"
//#include "phys_constants.h"
//#include "SystemOfUnits.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TTimer.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <TStyle.h>
#include <TROOT.h>


//Event Gen//
#include "EvtGen/include/EvtGen/EvtGen.hh"
#include "EvtGen/include/EvtGenBase/EvtParticle.hh"
#include "EvtGen/include/EvtGenBase/EvtParticleFactory.hh"
#include "EvtGen/include/EvtGenBase/EvtPatches.hh"
#include "EvtGen/include/EvtGenBase/EvtPDL.hh"
#include "EvtGen/include/EvtGenBase/EvtRandom.hh"
#include "EvtGen/include/EvtGenBase/EvtReport.hh"
#include "EvtGen/include/EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGen/include/EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGen/include/EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGen/include/EvtGenBase/EvtDecayBase.hh"
#include "EvtGen/include/EvtGenExternal/EvtExternalGenList.hh"

#include "PHEvtGenDecayer.h"
#include "PHEvtGenDecayer.cxx"



	//#include <hfmltrigger/HFMLTriggerHepMCTrigger.h>
	//R__LOAD_LIBRARY(libHFMLTrigger.so)
	using namespace std;

	bool RunBackground = false;
	bool RunSignal = false;
	bool RunEVTGEN = true;

	//PHPythia8 *pythia8 = NULL;



	//StarPythia8 *pythia8 = NULL;
	PHEvtGenDecayer* myEvtGenDecayer = NULL;

	double TrackPTMin = 0;
	double TrackPTMax = 20;

	int NEvents = 200;

	bool UseSTAR = true;

	double PhiMass = 1.019461;
	double DsMass = 1.96847;
	double KKMassWindow = 0.10;
	double PiKKMassWindow = 0.10;

	double minPtCut = 0.0;
	double etaCut = 1.1;

	float ThetaCut = -1;
	double DCACut = 200;
	double DCATrackCut = 0;

	bool mUseHist = false;
	bool mUseTree = true;

	const Int_t nPtBinsDca = 17;
	const Double_t ptEdgeDca[nPtBinsDca + 1] = {0 ,  0.1 , 0.5 , 0.6 , 0.7,  0.8 , 0.9 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2. , 2.4 , 3.0 , 4. , 6., 50. };


int DrawInput = 0;

int nPlus;
int nMinus;


bool NoCut = true;
bool doCut = false;

int MaxPlus = 10;
int MaxMinus = 10;
double PassPre;


int nPiPlus;
int nKPlus;
int nPiMinus;
int nKMinus;


float pthat;

Pythia8::Pythia pythia8;

float const sigmaPos0 = 15.2;
float const pxlLayer1Thickness = 0.00486;

const Int_t nParticles = 3;  //0--pi, 1--k, 2--p

TH2F* h2Dca[nParticles][nPtBinsDca];    //0--pi, 1--k, 2--p

int const nEtasDca = 1;
float const EtaEdgeDca[nEtasDca + 1] = { -1.1, 1.1};

int const nVzsDca = 1;
float const VzEdgeDca[nVzsDca + 1] = { -6.e4, 6.e4};



float BsMassPDG = 5.36684;

TFile * outFile = new TFile("output.root","RECREATE");



TTree * nt_sig = new TTree("BsDsKKPiPi","BsDsKKPiPi");

//Bs Kinematics
float BsPt;
float BsMass;
float BsY;


//Daughters Kinematics

float kpPt;
float kmPt;
float pipPt;
float pimPt;

float kpEta;
float kmEta;
float pipEta;
float pimEta;

float kpDca;
float kmDca;
float pipDca;
float pimDca;


float Bvtxx;
float Bvtxy;
float Bvtxz;

float dcadaughter;
float dcaToPv_Bs;
float decayLength_Bs;
float cosTheta_Bs;
float weight;
float primaryPt;
float primaryY;
float primaryPhi;

float ptWg; 

float dcapippim;
float dcakppip;
float dcakppim;
float dcakmpip;
float dcakmpim;
float dcakpkm;


int FromBsPiP;
int FromBsPiM;
int FromBsKP;
int FromBsKM;

int IsSigPiP;
int IsSigPiM;
int IsSigKP;
int IsSigKM;

int PairNow = 0;


float KKMass;
float PiKKMass;




float const acceptanceRapidity = 1.1; 


//TTree * nt_gen = new TTree("BsGenTree","BsGenTree");

int gEvent;
float gBvtxx;
float gBvtxy;
float gBvtxz;
float gBpx;
float gBpy;
float gBpz;
float gBpt;


int gInEvent;
float gInBvtxx;
float gInBvtxy;
float gInBvtxz;
float gInBpx;
float gInBpy;
float gInBpz;
float gInBpt;

TF1 * fKaonMomResolution = new TF1("fKaonMomResolution","0.0175 + 0.0011666667 * x",TrackPTMin,TrackPTMax);
TF1 * fPionMomResolution = new TF1("fPionMomResolution","0.0175 + 0.0011666667 * x",TrackPTMin,TrackPTMax);
TF1 * fProtonMomResolution = new TF1("fProtonMomResolution","0.0175 + 0.0011666667 * x",TrackPTMin,TrackPTMax);


//TF1 * hkPtWg = new TF1("hkPtWg","x*x*TMath::Exp(-x)",0,20);
//TF1 * hpiPtWg = new TF1("hpiPtWg","x*x*TMath::Exp(-x/2)",0,20);
TH1F * hkPtWg;
TH1F * hpiPtWg;
TH1F * hpPtWg;


TFile * fDca2D;
TFile * InputSpectra;

//int MaxPion = 20;
//int MaxKaon = 20;



vector<int> vecTrackID_km; // k minus
vector<float> vecPt_km;
vector<float> vecEta_km;
vector<float> vecPhi_km;
vector<TVector3> vecPos_km;
vector<int> vecType_km;
vector<bool> vecClean_km;
vector<bool> vecHybrid_km;
vector<bool> vecCleanIdeal_km;
vector<bool> vecHybridIdeal_km;
vector<float> vecDCA_km;
vector<int> vecBottom_km;


vector<int> vecTrackID_kp; 
vector<float> vecPt_kp; // k plus
vector<float> vecEta_kp;
vector<float> vecPhi_kp;
vector<TVector3> vecPos_kp;
vector<int> vecType_kp;
vector<bool> vecClean_kp;
vector<bool> vecHybrid_kp;
vector<bool> vecCleanIdeal_kp;
vector<bool> vecHybridIdeal_kp;
vector<float> vecDCA_kp;
vector<int> vecBottom_kp;



vector<int> vecTrackID_pip; 
vector<float> vecPt_pip; // pi plus
vector<float> vecEta_pip;
vector<float> vecPhi_pip;
vector<TVector3> vecPos_pip;
vector<int> vecType_pip;
vector<bool> vecClean_pip;
vector<bool> vecHybrid_pip;
vector<bool> vecCleanIdeal_pip;
vector<bool> vecHybridIdeal_pip;
vector<float> vecDCA_pip;
vector<int> vecBottom_pip;

vector<int> vecTrackID_pim; 
vector<float> vecPt_pim; // pi minus
vector<float> vecEta_pim;
vector<float> vecPhi_pim;
vector<TVector3> vecPos_pim;
vector<int> vecType_pim;
vector<bool> vecClean_pim;
vector<bool> vecHybrid_pim;
vector<bool> vecCleanIdeal_pim;
vector<bool> vecHybridIdeal_pim;
vector<float> vecDCA_pim;
vector<int> vecBottom_pim;


//void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters);
//void Fill(TLorentzVector const& kpRMom,  TLorentzVector const& kmRMom, TLorentzVector const& pipMom, TLorentzVector const& pimMom, TVector3 kpPos, TVector3 kmPos,TVector3 pipPos, TVector3 pimPos, TVector3 v0_Bs);
//void fillDouble(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TVector3 v00, TVector3 kPos, TVector3 pPos);

void GenerateBackground( int & EventID);
//void GetSignal(int & EventID);
//void GetGen(int & EventID);
//void SigEvttoEntry();



void getKinematics(TLorentzVector& b, double const mass);
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos);
TVector3 smearPosData(int iParticleIndex, double vz, int cent, TLorentzVector const& rMom, TVector3 const& pos);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
void initEvtGen();
void initPYTHIA8();
void GenerateDecay(int & EventID);
//void GenerateSignal(int & EventID);
void GenerateEvtGen(int & EventID);

int getPtIndexDca(double);
int getEtaIndexDca(double);
int getVzIndexDca(double);


void init();
void clean();
void end();
void BackgroundSimulations(int & EventID);


int Event;
float EventSignal;

double PI = 3.1415926538;

double M_KAON_MINUS=0.493677;
double M_PION_MINUS=0.13957018;

double M_KAON_PLUS=0.493677;
double M_PION_PLUS=0.13957018;


float gpx;
float gpy;
float gpz;
float gflavor;
float gvx;
float gvy;
float gvz;


float EventNow;



int EventSig;
int TotalCand;

std::vector<int> EvtToSigEntry;

int EntryMin;
int EntryMax;


TFile * SignalFile;
TTree * ntp_track;
TTree * ntp_gtrack;

const float ptCut = 0;
TVector3 const vertex(0,0,0);
TVector3 Pos(0,0,0);
TLorentzVector rcFourMom;
TVector3 rcPos;
TLorentzVector FourMom;


