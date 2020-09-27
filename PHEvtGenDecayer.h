/**
   \class PHEvtGenDecayer

   \brief PH wrapper for EvtGen Decayer

   Authors: Xiaozhi Bai (xiaozhi@uic.edu),
            Mustafa Mustafa (mmustafa@lbl.gov)
			Zhaozhong Shi for Modification
*/

#ifndef PHEvtGenDecayer__h
#define PHEvtGenDecayer__h

#include <cstddef>

#include "TVirtualMCDecayer.h"
#include "TString.h"
//#include "EvtGen/EvtGen.hh"
#include  "EvtGen/include/EvtGen/EvtGen.hh"
#include <TLorentzVector.h>
#include "EvtGen/include/EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGen/include/EvtGenBase/EvtRandomEngine.hh"


class TLorentzVector;
class TClonesArray;
class EvtStdlibRandomEngine;
class EvtParticle;
class EvtRandomEngine;

class PHEvtGenDecayer : public TVirtualMCDecayer
{
  public:
   PHEvtGenDecayer(EvtGen* evtGen = NULL);
   virtual ~PHEvtGenDecayer();

   virtual void Init();
   virtual void Decay(Int_t pdgId, TLorentzVector* p);
   virtual Int_t ImportParticles(TClonesArray* particles);
   virtual void SetForceDecay(Int_t type);
   virtual void ForceDecay();
   virtual Float_t GetPartialBranchingRatio(Int_t ipart);
   virtual Float_t GetLifetime(Int_t pdgid);
   virtual void ReadDecayTable();
	
   void setVertex(TLorentzVector* r);
   void setDecayTable(TString decayTable);
   void ClearEvent();
   void AppendParticle(Int_t pdg, TLorentzVector* _p);	
   void SetDecayTable(TString decayTable);
   void SetVertex(TLorentzVector* r);	


  private:
  // EvtStdlibRandomEngine* mEvtGenRandomEngine;
   EvtRandomEngine  * mEvtGenRandomEngine;
   EvtGen* mEvtGen;
   EvtParticle* mParticle;
   EvtVector4R r_init;
   bool   mOwner;
   bool   mDebug = false; 
   EvtVector4R * mVertex;
};

inline void PHEvtGenDecayer::setDecayTable(TString decayTable) { mEvtGen->readUDecay(decayTable); }
inline void PHEvtGenDecayer::setVertex(TLorentzVector* r) {r_init.set(r->T(),r->X(),r->Y(),r->Z()); }
#endif
