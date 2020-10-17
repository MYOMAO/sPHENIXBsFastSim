Int_t EventSignal;

//Bs Kinematics
float BsPt;
float BsMass;
float BsY;
float pthat;

//Daughters Kinematics

float kpPt;
float kmPt;
float epPt;
float emPt;

float kpEta;
float kmEta;
float epEta;
float emEta;

float kpDca;
float kmDca;
float epDca;
float emDca;


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

float dcaepem;
float dcakpep;
float dcakpem;
float dcakmep;
float dcakmem;
float dcakpkm;


int FromBsEP;
int FromBsEM;
int FromBsKP;
int FromBsKM;

int IsSigEP;
int IsSigEM;
int IsSigKP;
int IsSigKM;


float KKMass;
float eeMass;



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



vector<int> vecTrackID_ep; 
vector<float> vecPt_ep; // pi plus
vector<float> vecEta_ep;
vector<float> vecPhi_ep;
vector<TVector3> vecPos_ep;
vector<int> vecType_ep;
vector<bool> vecClean_ep;
vector<bool> vecHybrid_ep;
vector<bool> vecCleanIdeal_ep;
vector<bool> vecHybridIdeal_ep;
vector<float> vecDCA_ep;
vector<int> vecBottom_ep;

vector<int> vecTrackID_em; 
vector<float> vecPt_em; // pi minus
vector<float> vecEta_em;
vector<float> vecPhi_em;
vector<TVector3> vecPos_em;
vector<int> vecType_em;
vector<bool> vecClean_em;
vector<bool> vecHybrid_em;
vector<bool> vecCleanIdeal_em;
vector<bool> vecHybridIdeal_em;
vector<float> vecDCA_em;
vector<int> vecBottom_em;


Int_t EventGen;

//Bs Kinematics
float BsPtGen;
float BsMassGen;
float BsYGen;
float pthatGen;

//Daughters Kinematics

float kpPtGen;
float kmPtGen;
float epPtGen;
float emPtGen;

float kpEtaGen;
float kmEtaGen;
float epEtaGen;
float emEtaGen;

float kpDcaGen;
float kmDcaGen;
float epDcaGen;
float emDcaGen;


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

float dcaepemGen;
float dcakpepGen;
float dcakpemGen;
float dcakmepGen;
float dcakmemGen;
float dcakpkmGen;


int FromBsEPGen;
int FromBsEMGen;
int FromBsKPGen;
int FromBsKMGen;

int IsSigEPGen;
int IsSigEMGen;
int IsSigKPGen;
int IsSigKMGen;


float KKMassGen;
float eeMassGen;



vector<int> vecTrackID_kmGen; // k minus
vector<float> vecPt_kmGen;
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



vector<int> vecTrackID_epGen; 
vector<float> vecPt_epGen; // pi plus
vector<float> vecEta_epGen;
vector<float> vecPhi_epGen;
vector<TVector3> vecPos_epGen;
vector<int> vecType_epGen;
vector<bool> vecClean_epGen;
vector<bool> vecHybrid_epGen;
vector<bool> vecCleanIdeal_epGen;
vector<bool> vecHybridIdeal_epGen;
vector<float> vecDCA_epGen;
vector<int> vecBottom_epGen;

vector<int> vecTrackID_emGen; 
vector<float> vecPt_emGen; // pi minus
vector<float> vecEta_emGen;
vector<float> vecPhi_emGen;
vector<TVector3> vecPos_emGen;
vector<int> vecType_emGen;
vector<bool> vecClean_emGen;
vector<bool> vecHybrid_emGen;
vector<bool> vecCleanIdeal_emGen;
vector<bool> vecHybridIdeal_emGen;
vector<float> vecDCA_emGen;
vector<int> vecBottom_emGen;


