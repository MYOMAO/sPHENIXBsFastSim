Int_t EventSignal;

//Bs Kinematics
float BsPt;
float BsMass;
float BsY;
float pthat;
float BsP;
float BsLifeTime;

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



//Ds Vertex//

float dcadaughter_Ds;
float dcaToPv_Ds;
float decayLength_Ds;
float cosTheta_Ds;


float Dvtxx;
float Dvtxy;
float Dvtxz;



//Bs Vertex//

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
float BCandEffWg;

float Eff_PIP;
float Eff_KP;
float Eff_PIM;
float Eff_KM;


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



float KKMass;
float PiKKMass;

int PIDPIP;
int PIDKP;
int PIDPIM;
int PIDKM;


int MomPIDPIP;
int MomPIDKP;
int MomPIDPIM;
int MomPIDKM;

int MotherPIDPIP;
int MotherPIDKP;
int MotherPIDPIM;
int MotherPIDKM;

vector<int> vecTrackID_km; // k minus
vector<int> vecPID_km;
vector<int> vecMomPID_km; 
vector<int> vecMotherPID_km; 
vector<float> vecPt_km;
vector<float> vecEta_km;
vector<float> vecPhi_km;
vector<float> vecEff_km;
vector<TVector3> vecPos_km;
vector<int> vecType_km;
vector<bool> vecClean_km;
vector<bool> vecHybrid_km;
vector<bool> vecCleanIdeal_km;
vector<bool> vecHybridIdeal_km;
vector<float> vecDCA_km;
vector<int> vecBottom_km;



vector<int> vecTrackID_kp; 
vector<int> vecPID_kp; 
vector<int> vecMomPID_kp; 
vector<int> vecMotherPID_kp; 
vector<float> vecPt_kp; // k plus
vector<float> vecEta_kp;
vector<float> vecPhi_kp;
vector<float> vecEff_kp;
vector<TVector3> vecPos_kp;
vector<int> vecType_kp;
vector<bool> vecClean_kp;
vector<bool> vecHybrid_kp;
vector<bool> vecCleanIdeal_kp;
vector<bool> vecHybridIdeal_kp;
vector<float> vecDCA_kp;
vector<int> vecBottom_kp;



vector<int> vecTrackID_pip; 
vector<int> vecPID_pip; 
vector<int> vecMomPID_pip; 
vector<int> vecMotherPID_pip; 
vector<float> vecPt_pip; // pi plus
vector<float> vecEta_pip;
vector<float> vecPhi_pip;
vector<float> vecEff_pip;
vector<TVector3> vecPos_pip;
vector<int> vecType_pip;
vector<bool> vecClean_pip;
vector<bool> vecHybrid_pip;
vector<bool> vecCleanIdeal_pip;
vector<bool> vecHybridIdeal_pip;
vector<float> vecDCA_pip;
vector<int> vecBottom_pip;

vector<int> vecTrackID_pim; 
vector<int> vecPID_pim; 
vector<int> vecMomPID_pim; 
vector<int> vecMotherPID_pim; 
vector<float> vecPt_pim; // pi minus
vector<float> vecEta_pim;
vector<float> vecPhi_pim;
vector<float> vecEff_pim;
vector<TVector3> vecPos_pim;
vector<int> vecType_pim;
vector<bool> vecClean_pim;
vector<bool> vecHybrid_pim;
vector<bool> vecCleanIdeal_pim;
vector<bool> vecHybridIdeal_pim;
vector<float> vecDCA_pim;
vector<int> vecBottom_pim;


/*
vector<vector<int>> vecMomList_pip;
vector<vector<int>> vecMomList_pim;
vector<vector<int>> vecMomList_kp;
vector<vector<int>> vecMomList_km;



vector<int> PiPMotherList;
vector<int> PiMMotherList;
vector<int> KPMotherList;
vector<int> KMMotherList;

*/

Int_t EventGen;

//Bs Kinematics
float BsPtGen;
float BsMassGen;
float BsYGen;
float pthatGen;
float BsPGen;
float BsLifeTimeGen;

//Daughters Kinematics

float kpPtGen;
float kmPtGen;
float pipPtGen;
float pimPtGen;

float kpEtaGen;
float kmEtaGen;
float pipEtaGen;
float pimEtaGen;

float kpDcaGen;
float kmDcaGen;
float pipDcaGen;
float pimDcaGen;



//Ds Vertex//

float dcadaughter_DsGen;
float dcaToPv_DsGen;
float decayLength_DsGen;
float cosTheta_DsGen;


float DvtxxGen;
float DvtxyGen;
float DvtxzGen;

//Bs Vertex//

float BvtxxGen;
float BvtxyGen;
float BvtxzGen;

float dcadaughterGen;
float dcaToPv_BsGen;
float decayLength_BsGen;
float cosTheta_BsGen;
float weightGen;
float primaryPtGen;
float primaryYGen;
float primaryPhiGen;

float ptWgGen; 

float dcapippimGen;
float dcakppipGen;
float dcakppimGen;
float dcakmpipGen;
float dcakmpimGen;
float dcakpkmGen;


int FromBsPiPGen;
int FromBsPiMGen;
int FromBsKPGen;
int FromBsKMGen;

int IsSigPiPGen;
int IsSigPiMGen;
int IsSigKPGen;
int IsSigKMGen;


float KKMassGen;
float PiKKMassGen;

int PIDPIPGen;
int PIDKPGen;
int PIDPIMGen;
int PIDKMGen;


int MomPIDPIPGen;
int MomPIDKPGen;
int MomPIDPIMGen;
int MomPIDKMGen;

int MotherPIDPIPGen;
int MotherPIDKPGen;
int MotherPIDPIMGen;
int MotherPIDKMGen;



vector<int> vecTrackID_kmGen; // k minus
vector<float> vecPt_kmGen;
vector<float> vecPID_kmGen;
vector<float> vecMomPID_kmGen;
vector<float> vecMotherPID_kmGen;
vector<float> vecEta_kmGen;
vector<float> vecPhi_kmGen;
vector<TVector3> vecPos_kmGen;
vector<int> vecType_kmGen;
vector<bool> vecClean_kmGen;
vector<bool> vecHybrid_kmGen;
vector<bool> vecCleanIdeal_kmGen;
vector<bool> vecHybridIdeal_kmGen;
vector<float> vecDCA_kmGen;
vector<int> vecBottom_kmGen;


vector<int> vecTrackID_kpGen; 
vector<float> vecPID_kpGen;
vector<float> vecMomPID_kpGen;
vector<float> vecMotherPID_kpGen;
vector<float> vecPt_kpGen; // k plus
vector<float> vecEta_kpGen;
vector<float> vecPhi_kpGen;
vector<TVector3> vecPos_kpGen;
vector<int> vecType_kpGen;
vector<bool> vecClean_kpGen;
vector<bool> vecHybrid_kpGen;
vector<bool> vecCleanIdeal_kpGen;
vector<bool> vecHybridIdeal_kpGen;
vector<float> vecDCA_kpGen;
vector<int> vecBottom_kpGen;



vector<int> vecTrackID_pipGen; 
vector<float> vecPID_pipGen;
vector<float> vecMomPID_pipGen;
vector<float> vecMotherPID_pipGen;
vector<float> vecPt_pipGen; // pi plus
vector<float> vecEta_pipGen;
vector<float> vecPhi_pipGen;
vector<TVector3> vecPos_pipGen;
vector<int> vecType_pipGen;
vector<bool> vecClean_pipGen;
vector<bool> vecHybrid_pipGen;
vector<bool> vecCleanIdeal_pipGen;
vector<bool> vecHybridIdeal_pipGen;
vector<float> vecDCA_pipGen;
vector<int> vecBottom_pipGen;


vector<int> vecTrackID_pimGen;
vector<float> vecPID_pimGen;
vector<float> vecMomPID_pimGen;
vector<float> vecMotherPID_pimGen;
vector<float> vecPt_pimGen; // pi minus
vector<float> vecEta_pimGen;
vector<float> vecPhi_pimGen;
vector<TVector3> vecPos_pimGen;
vector<int> vecType_pimGen;
vector<bool> vecClean_pimGen;
vector<bool> vecHybrid_pimGen;
vector<bool> vecCleanIdeal_pimGen;
vector<bool> vecHybridIdeal_pimGen;
vector<float> vecDCA_pimGen;
vector<int> vecBottom_pimGen;


