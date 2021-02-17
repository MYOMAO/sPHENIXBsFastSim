#include "FastSim.h"
#include "VecBranch.h"
//R__LOAD_LIBRARY(EvtGen/lib/libEvtGen.so)
//R__LOAD_LIBRARY(EvtGen/lib/libEvtGenExternal.so)
#include "TF2.h"
#include "TRandom.h"
#include "TH2.h"

R__LOAD_LIBRARY(libEvtGenExternal.so)

	void FastSim(int NEvents){

		TotalEvents = NEvents;

		//declare vectors//

		Event = 0;

		init();


		//if(RunSignal)	SigEvttoEntry();



		if(RunEVTGEN){


			initPYTHIA8();

			cout << "Now Initiate EventGen" << endl;
			initEvtGen();

		}


		for(int i =0; i < NEvents; i++){

			if(i%10000 == 0)cout << "Now Working on Event = " << i << endl;
			Event = i;
			if(RunUncorrBkgd)	GenerateBackground(Event);
			//		if(RunSignal)	GetSignal(Event);
			//		if(RunSignal) GetGen(Event);
			if(RunCorrBkgd) GenerateDecay(Event);
			if(RunEVTGEN) GenerateEvtGen(Event);
			BackgroundSimulations(Event);
			
			if(RunEVTGEN == false) nTrk = 0;
			Evt = Event;
			Multi = nPiPlus + nKPlus + nPPlus + nPiMinus + nKMinus + nPMinus + nTrk;

			nt_evt->Fill();
			
		
			Event = Event + 1;
			


		}

		cout << "Now End" << endl;
		end();
	}


void BackgroundSimulations(int &EventID){





	TLorentzVector * KP =  new TLorentzVector; 
	TLorentzVector * KM =  new TLorentzVector; 
	TLorentzVector * PIP =  new TLorentzVector; 
	TLorentzVector * PIM =  new TLorentzVector; 


	//cout << "Pass Here" << endl;

	/*
	   if( vecPt_kp.size()  > 4 && vecPt_km.size() > 4){

	   cout << "nPlus = " << nPlus << "   vecPt_kp.size() = " << vecPt_kp.size() << "   vecPt_kpGen.size() = " << vecPt_kpGen.size() << endl;
	   cout << "nMinus = " << nMinus << "   vecPt_km.size() = " << vecPt_km.size() << "   vecPt_kmGen.size() = " << vecPt_kmGen.size() << endl;

	   }
	   */

	//cout << "nPiPlusGen = " << nPiPlusGen << "  nPiPlusComb =  " << nPiPlusComb << "  nPiPlusSig = " << nPiPlusSig << "   nPiPlusDecay = " << nPiPlusDecay << "   Sum of Shits = " << nPiPlusComb + nPiPlusSig + nPiPlusDecay << " vecPt_pipGen.size() = " << vecPt_pipGen.size() << " vecPt_kpGen.size() = " << vecPt_kpGen.size() << "   vecPt_kp = " <<   vecPt_kp.size() << "   vecPt_kpGen = " <<   vecPt_kpGen.size() << endl;



	outFile->cd();
	TVector3 v0_Bs(0,0,0);


	for(int i = 0; i < vecPt_kp.size(); i++){
		for(int j = 0; j <  vecPt_km.size(); j++){

			for(int k = 0; k <  vecPt_pip.size() ; k++){


				for(int l = 0; l < vecPt_pim.size(); l++){

					//	if(PairNow%100000==0) cout << "Now Working on Combination Pair = " << PairNow << endl;

					KP->SetPtEtaPhiM(vecPt_kp[i],vecEta_kp[i],vecPhi_kp[i],M_KAON_PLUS);
					KM->SetPtEtaPhiM(vecPt_km[j],vecEta_km[j],vecPhi_km[j],M_KAON_MINUS);
					PIP->SetPtEtaPhiM(vecPt_pip[k],vecEta_pip[k],vecPhi_pip[k],M_PION_PLUS);
					PIM->SetPtEtaPhiM(vecPt_pim[l],vecEta_pim[l],vecPhi_pim[l],M_PION_MINUS);

				
					if(i !=k && j != l){  

						 //Guaranteeing no signal track combine itself 

						if(vecTrackID_kp[i] == vecTrackID_pip[k])  continue;   
						if(vecTrackID_km[i] == vecTrackID_pim[k])  continue;

						EventSignal = EventID;
						FromBsKP = vecBottom_kp[i];
						FromBsKM = vecBottom_km[j];
						FromBsPiP = vecBottom_pip[k];
						FromBsPiM = vecBottom_pim[l];

						IsSigKP = vecType_kp[i];
						IsSigKM = vecType_km[j];
						IsSigPiP = vecType_pip[k];
						IsSigPiM = vecType_pim[l];

						PIDKP = vecPID_kp[i];
						PIDKM = vecPID_km[j];
						PIDPIP = vecPID_pip[k];
						PIDPIM = vecPID_pim[l];

						MomPIDKP = vecMomPID_kp[i];
						MomPIDKM = vecMomPID_km[j];
						MomPIDPIP = vecMomPID_pip[k];
						MomPIDPIM = vecMomPID_pim[l];

						MotherPIDKP = vecMotherPID_kp[i];
						MotherPIDKM = vecMotherPID_km[j];
						MotherPIDPIP = vecMotherPID_pip[k];
						MotherPIDPIM = vecMotherPID_pim[l];


						Eff_KP = vecEff_kp[i];
						Eff_KM = vecEff_km[j];
						Eff_PIP = vecEff_pip[k];
						Eff_PIM = vecEff_pim[l];
	
						BCandEffWg = Eff_KP * Eff_KM * Eff_PIP * Eff_PIM;

						TLorentzVector const BsMom = *KP + *KM + *PIP + *PIM;
						TLorentzVector const PhiMom = *KP + *KM;
					//	TLorentzVector const DsMom = *KP + *KM + *PIM;
						//TLorentzVector * DsMom = new TLorentzVector;
					//	if(PiPFromDs == 0) *DsMom = *KP + *KM + *PIM;
					//	if(PiPFromDs == 1) *DsMom = *KP + *KM + *PIP;
					//	if(PiPFromDs == -1) *DsMom = *KP + *KM + *PIM;
						*DsMom = *KP + *KM + *PIM;


						TVector3 v0pippim,v0kppip,v0kppim,v0kmpip,v0kmpim,v0kpkm;
						dcapippim = dca1To2(PIP->Vect(), vecPos_pip[k], PIM->Vect(), vecPos_pim[l], v0pippim);
						dcakppip = dca1To2(KP->Vect(), vecPos_kp[i], PIP->Vect(), vecPos_pip[k], v0kppip);
						dcakppim = dca1To2(KP->Vect(), vecPos_kp[i], PIM->Vect(), vecPos_pim[l], v0kppim);
						dcakmpip = dca1To2(KM->Vect(), vecPos_km[j], PIP->Vect(), vecPos_pip[k], v0kmpip);
						dcakmpim = dca1To2(KM->Vect(), vecPos_km[j], PIM->Vect(), vecPos_pim[l], v0kmpim);
						dcakpkm = dca1To2(KP->Vect(), vecPos_kp[i], KM->Vect(), vecPos_km[j], v0kpkm);

						TVector3 tmp = v0pippim+v0kppip+v0kppim+v0kmpip+v0kmpim+v0kpkm;
						TVector3 v0_Bs(tmp.X()/6.,tmp.Y()/6.,tmp.Z()/6.);



						decayLength_Bs = (v0_Bs - vertex).Mag();
						dcaToPv_Bs = dca(BsMom.Vect(), v0_Bs, vertex);
						cosTheta_Bs = (v0_Bs - vertex).Unit().Dot(BsMom.Vect().Unit());


						Bvtxx = v0_Bs.X();
						Bvtxy = v0_Bs.Y();
						Bvtxz = v0_Bs.Z();

						kpDca = vecDCA_kp[i];
						kmDca = vecDCA_km[j];
						pipDca = vecDCA_pip[k];
						pimDca = vecDCA_pim[l];




						kpPt = KP->Perp();
						kmPt = KM->Perp();
						pipPt = PIP->Perp();
						pimPt = PIM->Perp();

						kpEta = KP->Rapidity();
						kmEta = KM->Rapidity();
						pipEta = PIP->Rapidity();
						pimEta = PIM->Rapidity();

					//	if((RunUncorrBkgd == false && RunCorrBkgd == false && RunEVTGEN == true) && !(IsSigKP == 1 && PIDPIP == 211 && PIDPIM == -211 && PIDKP == 321 && PIDKM == -321))  continue;
						if((RunUncorrBkgd == false && RunCorrBkgd == false && RunEVTGEN == true) && !(IsSigPiP*IsSigPiM*IsSigKP*IsSigKM==1)) continue;

						bool afterAcc = false; 
						if(kpPt > minPtCut && kmPt >  minPtCut && pipPt > minPtCut && pimPt > minPtCut && fabs(kpEta) < etaCut && fabs(kmEta) < etaCut && fabs(pipEta) < etaCut && fabs(pimEta) < etaCut) afterAcc = true;
						if(!afterAcc) continue;


						dcadaughter = dcapippim;						
						if(dcadaughter > dcakppip) dcadaughter = dcakppip;
						if(dcadaughter > dcakppim) dcadaughter = dcakppim;
						if(dcadaughter > dcakmpip) dcadaughter = dcakmpip;
						if(dcadaughter > dcakmpim) dcadaughter = dcakmpim;
						if(dcadaughter > dcakpkm) dcadaughter = dcakpkm;




						bool PassDCA = false;
						if(dcadaughter < DCACut) PassDCA = true;

						if(!PassDCA) continue;



						BsPt = BsMom.Perp();
						BsMass = BsMom.M();
						BsY = BsMom.Rapidity();
						KKMass = PhiMom.M();
						PiKKMass = DsMom->M();


						if(fabs(KKMass - PhiMass) > KKMassWindow) continue;
						if(fabs(PiKKMass - DsMass) > PiKKMassWindow) continue;

						ptWg = 1;
						weight = 1;

						if(fabs(BsMom.Rapidity()) > BetaCut) continue;
						if(cosTheta_Bs < ThetaCut) continue;

						if(BsMass > BsMassUp || BsMass < BsMassDown) continue;
						
						//Now Ds//


						TVector3 tmp_Ds = v0kppim+v0kmpim+v0kpkm;
						TVector3 v0_Ds(tmp_Ds.X()/3.,tmp_Ds.Y()/3.,tmp_Ds.Z()/3.);


						decayLength_Ds = (v0_Ds - vertex).Mag();
						dcaToPv_Ds = dca(DsMom->Vect(), v0_Bs, vertex);
						cosTheta_Ds = (v0_Ds - vertex).Unit().Dot(DsMom->Vect().Unit());


						Dvtxx = v0_Ds.X();
						Dvtxy = v0_Ds.Y();
						Dvtxz = v0_Ds.Z();

						//Life Time Calculations//
						BsP = sqrt(BsMom.Perp() * BsMom.Perp() + BsMom.Pz() *  BsMom.Pz());
						BsLifeTime = decayLength_Ds * BsMass/BsP;


					
						nt_sig->Fill();

						EffWeightCounts = EffWeightCounts + BCandEffWg;
						//Fill(*KP,*KM,*PIP,*PIM,  vecPos_kp[i], vecPos_km[j], vecPos_pip[k],vecPos_pim[l], v0_Bs);
					}
					//	nt_sig->Fill();

					PairNow = PairNow + 1;
				}

			}

		}

	}





	for(int i = 0; i < vecPt_kpGen.size(); i++){
		for(int j = 0; j <  vecPt_kmGen.size(); j++){

			for(int k = 0; k <  vecPt_pipGen.size() ; k++){


				for(int l = 0; l < vecPt_pimGen.size(); l++){

					//	if(PairNow%100000==0) cout << "Now Working on Combination Pair = " << PairNow << endl;

					KP->SetPtEtaPhiM(vecPt_kpGen[i],vecEta_kpGen[i],vecPhi_kpGen[i],M_KAON_PLUS);
					KM->SetPtEtaPhiM(vecPt_kmGen[j],vecEta_kmGen[j],vecPhi_kmGen[j],M_KAON_MINUS);
					PIP->SetPtEtaPhiM(vecPt_pipGen[k],vecEta_pipGen[k],vecPhi_pipGen[k],M_PION_PLUS);
					PIM->SetPtEtaPhiM(vecPt_pimGen[l],vecEta_pimGen[l],vecPhi_pimGen[l],M_PION_MINUS);

					int rejectionflag = 0;
					
					if(AllowGenPID == 0){
						if(i !=k && j != l){  

							//Guaranteeing no signal track combine itself 

							if(vecTrackID_kp[i] == vecTrackID_pip[k])  rejectionflag = 1;   
							if(vecTrackID_km[i] == vecTrackID_pim[k])  rejectionflag = 1;

						}
					}

					if(AllowGenPID == 1) rejectionflag = 0;

					if(rejectionflag == 0){


						EventGen = EventID;
						FromBsKPGen  = vecBottom_kpGen[i];
						FromBsKMGen  = vecBottom_kmGen[j];				
						FromBsPiPGen  = vecBottom_pipGen[k];
						FromBsPiMGen  = vecBottom_pimGen[l];


						/*

						   PiPMotherList = vecMomList_pip[i];
						   PiMMotherList = vecMomList_pim[i];
						   KPMotherList = vecMomList_kp[i];
						   KMMotherList = vecMomList_km[i];
						   */

						IsSigKPGen  = vecType_kpGen[i];
						IsSigKMGen  = vecType_kmGen[j];
						IsSigPiPGen  = vecType_pipGen[k];
						IsSigPiMGen  = vecType_pimGen[l];



						PIDKPGen = vecPID_kpGen[i];
						PIDKMGen = vecPID_kmGen[j];
						PIDPIPGen = vecPID_pipGen[k];
						PIDPIMGen = vecPID_pimGen[l];

						MomPIDKPGen = vecMomPID_kpGen[i];
						MomPIDKMGen = vecMomPID_kmGen[j];
						MomPIDPIPGen = vecMomPID_pipGen[k];
						MomPIDPIMGen = vecMomPID_pimGen[l];

						MotherPIDKPGen = vecMotherPID_kpGen[i];
						MotherPIDKMGen = vecMotherPID_kmGen[j];
						MotherPIDPIPGen = vecMotherPID_pipGen[k];
						MotherPIDPIMGen = vecMotherPID_pimGen[l];


						TLorentzVector const BsMom = *KP + *KM + *PIP + *PIM;
						TLorentzVector const PhiMom = *KP + *KM;
						TLorentzVector * DsMom = new TLorentzVector;
						if(PiPFromDs == 0) *DsMom = *KP + *KM + *PIM;
						if(PiPFromDs == 1) *DsMom = *KP + *KM + *PIP;
						if(PiPFromDs == -1) *DsMom = *KP + *KM + *PIM;
	
						TVector3 v0pippim,v0kppip,v0kppim,v0kmpip,v0kmpim,v0kpkm;
						dcapippim = dca1To2(PIP->Vect(), vecPos_pipGen[k], PIM->Vect(), vecPos_pimGen[l], v0pippim);
						dcakppip = dca1To2(KP->Vect(), vecPos_kpGen[i], PIP->Vect(), vecPos_pipGen[k], v0kppip);
						dcakppim = dca1To2(KP->Vect(), vecPos_kpGen[i], PIM->Vect(), vecPos_pimGen[l], v0kppim);
						dcakmpip = dca1To2(KM->Vect(), vecPos_kmGen[j], PIP->Vect(), vecPos_pipGen[k], v0kmpip);
						dcakmpim = dca1To2(KM->Vect(), vecPos_kmGen[j], PIM->Vect(), vecPos_pimGen[l], v0kmpim);
						dcakpkm = dca1To2(KP->Vect(), vecPos_kpGen[i], KM->Vect(), vecPos_kmGen[j], v0kpkm);


						TVector3 tmp = v0pippim+v0kppip+v0kppim+v0kmpip+v0kmpim+v0kpkm;
						TVector3 v0_Bs(tmp.X()/6.,tmp.Y()/6.,tmp.Z()/6.);

						decayLength_BsGen = (v0_Bs - vertex).Mag();
						dcaToPv_BsGen = dca(BsMom.Vect(), v0_Bs, vertex);
						cosTheta_BsGen = (v0_Bs - vertex).Unit().Dot(BsMom.Vect().Unit());


						BvtxxGen = v0_Bs.X();
						BvtxyGen = v0_Bs.Y();
						BvtxzGen = v0_Bs.Z();

						kpDcaGen = vecDCA_kpGen[i];
						kmDcaGen = vecDCA_kmGen[j];
						pipDcaGen = vecDCA_pipGen[k];
						pimDcaGen = vecDCA_pimGen[l];



						kpPtGen = KP->Perp();
						kmPtGen = KM->Perp();
						pipPtGen = PIP->Perp();
						pimPtGen = PIM->Perp();

						kpEtaGen = KP->Rapidity();
						kmEtaGen = KM->Rapidity();
						pipEtaGen = PIP->Rapidity();
						pimEtaGen = PIM->Rapidity();


						bool afterAcc = false; 
					//	if(kpPtGen > minPtCut && kmPtGen >  minPtCut && pipPtGen > minPtCut && pimPtGen > minPtCut && fabs(kpEtaGen) < etaCut && fabs(kmEtaGen) < etaCut && fabs(pipEtaGen) < etaCut && fabs(pimEtaGen) < etaCut) afterAcc = true;
						if(!(IsSigPiPGen*IsSigPiMGen*IsSigKPGen*IsSigKMGen==1)) continue;
						
						if(kpPtGen > minPtCutGen && kmPtGen >  minPtCutGen && pipPtGen > minPtCutGen && pimPtGen > minPtCutGen && fabs(kpEtaGen) < etaCutGen && fabs(kmEtaGen) < etaCutGen && fabs(pipEtaGen) < etaCutGen && fabs(pimEtaGen) < etaCutGen) afterAcc = true;
						if(!afterAcc) continue;



						dcadaughterGen = dcapippim;
						if(dcadaughterGen > dcakppip) dcadaughterGen = dcakppip;
						if(dcadaughterGen > dcakppim) dcadaughterGen = dcakppim;
						if(dcadaughterGen > dcakmpip) dcadaughterGen = dcakmpip;
						if(dcadaughterGen > dcakmpim) dcadaughterGen = dcakmpim;
						if(dcadaughterGen > dcakpkm) dcadaughterGen = dcakpkm;

						bool PassDCA = false;
						if(dcadaughterGen < DCACut) PassDCA = true;

					//	if(!PassDCA) continue;


						BsPtGen = BsMom.Perp();
						BsMassGen = BsMom.M();
						BsYGen = BsMom.Rapidity();
						KKMassGen = PhiMom.M();
						PiKKMassGen = DsMom->M();


						//if(fabs(KKMassGen - PhiMass) > KKMassWindow) continue;
						//if(fabs(PiKKMassGen - DsMass) > PiKKMassWindow) continue;

						ptWgGen = 1;
						weightGen = 1;

						//if(fabs(BsMom.Rapidity()) > 1.1) continue;
						//if(cosTheta_BsGen < ThetaCut) continue;
					

						//Now Ds//


						TVector3 tmp_Ds = v0kppim+v0kmpim+v0kpkm;
						TVector3 v0_Ds(tmp_Ds.X()/3.,tmp_Ds.Y()/3.,tmp_Ds.Z()/3.);


						decayLength_DsGen = (v0_Ds - vertex).Mag();
						dcaToPv_DsGen = dca(DsMom->Vect(), v0_Bs, vertex);
						cosTheta_DsGen = (v0_Ds - vertex).Unit().Dot(DsMom->Vect().Unit());


						DvtxxGen = v0_Ds.X();
						DvtxyGen = v0_Ds.Y();
						DvtxzGen = v0_Ds.Z();

						//Life Time Calculations//
						BsPGen = sqrt(BsMom.Perp() * BsMom.Perp() + BsMom.Pz() *  BsMom.Pz());
						BsLifeTimeGen = decayLength_DsGen * BsMassGen/BsPGen;



						nt_Gen->Fill();

						//Fill(*KP,*KM,*PIP,*PIM,  vecPos_kp[i], vecPos_km[j], vecPos_pip[k],vecPos_pim[l], v0_Bs);
					}
					//	nt_sig->Fill();

					PairNow = PairNow + 1;
				}

			}

		}

	}

	//	cout << "Pass 3" << endl;

	nPlus = vecPt_kp.size();
	nMinus = vecPt_km.size();


	double TotalEvents = vecPt_kp.size() * vecPt_pip.size() * vecPt_km.size() * vecPt_pim.size();
	double PassEvents = nt_sig->GetEntries() - PassPre;
	double Efficiency = 0;	
	double EfficiencyWeighted = 0;	
	
	if(TotalEvents > 0) Efficiency = PassEvents/TotalEvents;
//	if(TotalEvents > 0) EfficiencyWeighted = EffWeightCounts/TotalEvents;


	PassPre = nt_sig->GetEntries();
//	if(Efficiency > 0) cout << "Total Pairs = " << TotalEvents << "   Passed Pairs = " << PassEvents << "   Efficiency =  " << Efficiency << endl;
//	if(Efficiency > 0) cout << "Efficiency Weighted:  Total Pairs = " << TotalEvents << "   Passed Pairs = " << EffWeightCounts << "   Efficiency =  " << EfficiencyWeighted << endl;

	clean();

}



void GenerateBackground( int &EventID){


	gRandom->SetSeed();


	//Initial Histogram




	//	  Double_t nPions_tmp, nKaons_tmp;

	nPiPlus = gRandom->Integer(MaxPlus) + BaseCharge;
	nKPlus = gRandom->Integer(MaxMinus) + BaseCharge;

	nPPlus = gRandom->Integer(MaxProton) + BaseProton;

	nPiPlusGen = nPiPlus;


	nPMinus = nPPlus;
	nPiMinus = nPiPlus;
	nKMinus = nKPlus;


	nPlus = nPiPlus + nKPlus + nPPlus;
	nMinus = nPiMinus + nKMinus + nPMinus;
	nProtons = nPPlus + nPMinus;

	cout << "nPlus = " << nPlus << "    nMinus = " << nMinus << "  nProtons = " << nProtons << endl;
	Pos.SetXYZ(0,0,0);

	nPiPlusComb = 0;

	/*
	   vector<int> MomList;
	   MomList.clear();
	   MomList.push_back(0);
	   */

	//Assuming Pi+ = Pi-//

	//Pi+ Loop
	for(int ipi=0; ipi<nPiPlus; ipi++) {
		//	cout << "hpiPtWg->Integral() = " << hpiPtWg->Integral() << endl;	
		float pt = hpiPtWg->GetRandom();//gRandom->Uniform(0.6,20);
		//	cout << "Pass Inside" << endl;
		float eta = gRandom->Uniform(-1,1);
		float phi = gRandom->Uniform(-PI,PI);
		FourMom.SetPtEtaPhiM(pt, eta , phi, M_PION_PLUS);
		float rcdca = dca(FourMom.Vect(), Pos, vertex);

		//		if(abs(FourMom.Eta()) > TrackEtaCut) continue;
		//		if(abs(FourMom.Perp()) < minPtCut) continue;
		//		if(rcdca < DCATrackCut) continue;

		/*
		   vecMomList_pip.push_back(MomList);
		   vecMomList_kp.push_back(MomList);
		   */

		nPiPlusComb = nPiPlusComb + 1;

		if(SaveCombGen == 1){
			vecTrackID_pipGen.push_back(ipi);
			vecDCA_pipGen.push_back(rcdca);   	
			vecPt_pipGen.push_back(pt);
			vecEta_pipGen.push_back(eta);
			vecPhi_pipGen.push_back(phi);
			vecPos_pipGen.push_back(Pos);
			vecType_pipGen.push_back(0);
			vecBottom_pipGen.push_back(0); 


			vecPID_pipGen.push_back(211);
			vecMomPID_pipGen.push_back(-99999);
			vecMotherPID_pipGen.push_back(-99999);


			if(AllowGenPID == 0){
				vecTrackID_kpGen.push_back(ipi);
				vecDCA_kpGen.push_back(rcdca);               
				vecPt_kpGen.push_back(pt);
				vecEta_kpGen.push_back(eta);
				vecPhi_kpGen.push_back(phi);
				vecPos_kpGen.push_back(Pos);
				vecType_kpGen.push_back(0);
				vecBottom_kpGen.push_back(0); 

				vecPID_kpGen.push_back(211);
				vecMomPID_kpGen.push_back(-99999);
				vecMotherPID_kpGen.push_back(-99999);
			}
		}
		TrackEff = fTpcPi->Eval(pt);
		rcFourMom = smearMom(FourMom,fPionMomResolution);
		if (rcFourMom.Perp()<minPtCut) continue;
		rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
		if(rcdca < DCATrackCut) continue;


		//Identify as Pi+//
		vecTrackID_pip.push_back(ipi);
		vecDCA_pip.push_back(rcdca);               
		vecPt_pip.push_back(rcFourMom.Perp());
		vecEta_pip.push_back(rcFourMom.Eta());
		vecPhi_pip.push_back(rcFourMom.Phi());
		vecPos_pip.push_back(rcPos);
		vecType_pip.push_back(0);
		vecBottom_pip.push_back(0); 

		vecPID_pip.push_back(211);
		vecMomPID_pip.push_back(-99999);
		vecMotherPID_pip.push_back(-99999);
		vecEff_pip.push_back(TrackEff);

		//Identify as K+//
		vecTrackID_kp.push_back(ipi);
		vecDCA_kp.push_back(rcdca);               
		vecPt_kp.push_back(rcFourMom.Perp());
		vecEta_kp.push_back(rcFourMom.Eta());
		vecPhi_kp.push_back(rcFourMom.Phi());
		vecPos_kp.push_back(rcPos);
		vecType_kp.push_back(0);
		vecBottom_kp.push_back(0); 
		vecEff_kp.push_back(TrackEff);

		vecPID_kp.push_back(211);
		vecMomPID_kp.push_back(-99999);
		vecMotherPID_kp.push_back(-99999);

	}


	//Pi- Loop
	for(int ipi=0; ipi<nPiMinus; ipi++) {
		//	cout << "hpiPtWg->Integral() = " << hpiPtWg->Integral() << endl;	
		float pt = hpiPtWg->GetRandom();//gRandom->Uniform(0.6,20);
		//	cout << "Pass Inside" << endl;
		float eta = gRandom->Uniform(-1,1);
		float phi = gRandom->Uniform(-PI,PI);
		FourMom.SetPtEtaPhiM(pt, eta , phi, M_PION_MINUS);
		float rcdca = dca(FourMom.Vect(), Pos, vertex);


		/*
		   vecMomList_pim.push_back(MomList);
		   vecMomList_km.push_back(MomList);
		   */

		//		if(abs(FourMom.Eta()) > TrackEtaCut) continue;
		//		if(abs(FourMom.Perp()) < minPtCut) continue;
		//		if(rcdca < DCATrackCut) continue;


		if(SaveCombGen == 1){

			vecTrackID_pimGen.push_back(ipi);
			vecDCA_pimGen.push_back(rcdca);   	
			vecPt_pimGen.push_back(pt);
			vecEta_pimGen.push_back(eta);
			vecPhi_pimGen.push_back(phi);
			vecPos_pimGen.push_back(Pos);
			vecType_pimGen.push_back(0);
			vecBottom_pimGen.push_back(0); 

			vecPID_pimGen.push_back(-211);
			vecMomPID_pimGen.push_back(-99999);
			vecMotherPID_pimGen.push_back(-99999);


			if(AllowGenPID == 0){

				vecTrackID_kmGen.push_back(ipi);
				vecDCA_kmGen.push_back(rcdca);               
				vecPt_kmGen.push_back(pt);
				vecEta_kmGen.push_back(eta);
				vecPhi_kmGen.push_back(phi);
				vecPos_kmGen.push_back(Pos);
				vecType_kmGen.push_back(0);
				vecBottom_kmGen.push_back(0); 

				vecPID_kmGen.push_back(-211);
				vecMomPID_kmGen.push_back(-99999);
				vecMotherPID_kmGen.push_back(-99999);

			}

		}
		TrackEff = fTpcPi->Eval(pt);
		rcFourMom = smearMom(FourMom,fPionMomResolution);
		if (rcFourMom.Perp()<minPtCut) continue;
		rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
		if(rcdca < DCATrackCut) continue;


		//Identify as Pi-//
		vecTrackID_pim.push_back(ipi);
		vecDCA_pim.push_back(rcdca);               
		vecPt_pim.push_back(rcFourMom.Perp());
		vecEta_pim.push_back(rcFourMom.Eta());
		vecPhi_pim.push_back(rcFourMom.Phi());
		vecPos_pim.push_back(rcPos);
		vecType_pim.push_back(0);
		vecBottom_pim.push_back(0); 
		vecEff_pim.push_back(TrackEff);

		vecPID_pim.push_back(-211);
		vecMomPID_pim.push_back(-99999);
		vecMotherPID_pim.push_back(-99999);



		//Identify as K-//
		vecTrackID_km.push_back(ipi);
		vecDCA_km.push_back(rcdca);               
		vecPt_km.push_back(rcFourMom.Perp());
		vecEta_km.push_back(rcFourMom.Eta());
		vecPhi_km.push_back(rcFourMom.Phi());
		vecPos_km.push_back(rcPos);
		vecType_km.push_back(0);
		vecBottom_km.push_back(0); 
		vecEff_km.push_back(TrackEff);


		vecPID_km.push_back(-211);
		vecMomPID_km.push_back(-99999);
		vecMotherPID_km.push_back(-99999);


	}




	/*
	   for(int ik=0; ik<nPions; ik++) {
	   float pt = hkPtWg->GetRandom();//gRandom->Uniform(0.6,20);
	   float eta = gRandom->Uniform(-1,1);
	   float phi = gRandom->Uniform(-PI,PI);
	   FourMom.SetPtEtaPhiM(pt, eta , phi, M_KAON_PLUS);
	   rcFourMom = smearMom(FourMom,fKaonMomResolution);
	   if (rcFourMom.Perp()<minPtCut) continue;
	   rcPos = smearPosData(1, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
	   float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
	   if(rcdca < DCATrackCut) continue;
	   vecDCA_kp.push_back(rcdca);               
	   vecPt_kp.push_back(rcFourMom.Perp());
	   vecEta_kp.push_back(rcFourMom.Eta());
	   vecPhi_kp.push_back(rcFourMom.Phi());
	   vecPos_kp.push_back(rcPos);
	   vecType_kp.push_back(1);
	   vecBottom_kp.push_back(0); 

	   }
	   */



	//K+ Loop
	for(int ik=0; ik<nKPlus; ik++) {
		float pt = hkPtWg->GetRandom();//gRandom->Uniform(0.6,20);
		float eta = gRandom->Uniform(-1,1);
		float phi = gRandom->Uniform(-PI,PI);
		FourMom.SetPtEtaPhiM(pt, eta , phi, M_KAON_PLUS);	
		float rcdca = dca(FourMom.Vect(), Pos, vertex);

		/*
		   vecMomList_kp.push_back(MomList);
		   vecMomList_pip.push_back(MomList);
		   */

		//		if(abs(FourMom.Eta()) > TrackEtaCut) continue;
		//		if(abs(FourMom.Perp()) < minPtCut) continue;
		//		if(rcdca < DCATrackCut) continue;


		//	nPiPlusComb = nPiPlusComb + 1;
	
				if(SaveCombGen == 1){
		if(AllowGenPID == 0){

			vecTrackID_pipGen.push_back(ik+nPiPlus);
			vecDCA_pipGen.push_back(rcdca);   	
			vecPt_pipGen.push_back(pt);
			vecEta_pipGen.push_back(eta);
			vecPhi_pipGen.push_back(phi);
			vecPos_pipGen.push_back(Pos);
			vecType_pipGen.push_back(0);
			vecBottom_pipGen.push_back(0); 

			vecPID_pipGen.push_back(321);
			vecMomPID_pipGen.push_back(-99999);
			vecMotherPID_pipGen.push_back(-99999);

		}


		vecTrackID_kpGen.push_back(ik+nPiPlus);
		vecDCA_kpGen.push_back(rcdca);               
		vecPt_kpGen.push_back(pt);
		vecEta_kpGen.push_back(eta);
		vecPhi_kpGen.push_back(phi);
		vecPos_kpGen.push_back(Pos);
		vecType_kpGen.push_back(0);
		vecBottom_kpGen.push_back(0); 


		vecPID_kpGen.push_back(321);
		vecMomPID_kpGen.push_back(-99999);
		vecMotherPID_kpGen.push_back(-99999);

				}

		TrackEff = fTpcK->Eval(pt);
		rcFourMom = smearMom(FourMom,fKaonMomResolution);
		if (rcFourMom.Perp()<minPtCut) continue;
		rcPos = smearPosData(1, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
		if(rcdca < DCATrackCut) continue;

		//Identify as Pi+//
		vecTrackID_pip.push_back(ik+nPiPlus);
		vecDCA_pip.push_back(rcdca);               
		vecPt_pip.push_back(rcFourMom.Perp());
		vecEta_pip.push_back(rcFourMom.Eta());
		vecPhi_pip.push_back(rcFourMom.Phi());
		vecPos_pip.push_back(rcPos);
		vecType_pip.push_back(0);
		vecBottom_pip.push_back(0); 
		vecEff_pip.push_back(TrackEff);

		vecPID_pip.push_back(321);
		vecMomPID_pip.push_back(-99999);
		vecMotherPID_pip.push_back(-99999);


		//Identify as K+//

		vecTrackID_kp.push_back(ik+nPiPlus);
		vecDCA_kp.push_back(rcdca);               
		vecPt_kp.push_back(rcFourMom.Perp());
		vecEta_kp.push_back(rcFourMom.Eta());
		vecPhi_kp.push_back(rcFourMom.Phi());
		vecPos_kp.push_back(rcPos);
		vecType_kp.push_back(0);
		vecBottom_kp.push_back(0); 
		vecEff_kp.push_back(TrackEff);

		vecPID_kp.push_back(321);
		vecMomPID_kp.push_back(-99999);
		vecMotherPID_kp.push_back(-99999);

	}


	//K- Loop
	for(int ik=0; ik<nKMinus; ik++) {
		float pt = hkPtWg->GetRandom();//gRandom->Uniform(0.6,20);
		float eta = gRandom->Uniform(-1,1);
		float phi = gRandom->Uniform(-PI,PI);
		FourMom.SetPtEtaPhiM(pt, eta , phi, M_KAON_MINUS);
		float rcdca = dca(FourMom.Vect(), Pos, vertex);

		/*
		   vecMomList_km.push_back(MomList);
		   vecMomList_pim.push_back(MomList);
		   */

		//		if(abs(FourMom.Eta()) > TrackEtaCut) continue;
		//		if(abs(FourMom.Perp()) < minPtCut) continue;
		//		if(rcdca < DCATrackCut) continue;

		if(SaveCombGen == 1){


			if(AllowGenPID == 0){


				vecTrackID_pimGen.push_back(ik+nPiMinus);
				vecDCA_pimGen.push_back(rcdca);   	
				vecPt_pimGen.push_back(pt);
				vecEta_pimGen.push_back(eta);
				vecPhi_pimGen.push_back(phi);
				vecPos_pimGen.push_back(Pos);
				vecType_pimGen.push_back(0);
				vecBottom_pimGen.push_back(0); 

				vecPID_pimGen.push_back(-321);
				vecMomPID_pimGen.push_back(-99999);
				vecMotherPID_pimGen.push_back(-99999);

			}


			vecTrackID_kmGen.push_back(ik+nPiMinus);
			vecDCA_kmGen.push_back(rcdca);               
			vecPt_kmGen.push_back(pt);
			vecEta_kmGen.push_back(eta);
			vecPhi_kmGen.push_back(phi);
			vecPos_kmGen.push_back(Pos);
			vecType_kmGen.push_back(0);
			vecBottom_kmGen.push_back(0); 

			vecPID_kmGen.push_back(-321);
			vecMomPID_kmGen.push_back(-99999);
			vecMotherPID_kmGen.push_back(-99999);

		}
		TrackEff = fTpcK->Eval(pt);
		rcFourMom = smearMom(FourMom,fKaonMomResolution);
		if (rcFourMom.Perp()<minPtCut) continue;
		rcPos = smearPosData(1, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
		if(rcdca < DCATrackCut) continue;

		//Identify as Pi-//
		vecTrackID_pim.push_back(ik+nPiMinus);
		vecDCA_pim.push_back(rcdca);               
		vecPt_pim.push_back(rcFourMom.Perp());
		vecEta_pim.push_back(rcFourMom.Eta());
		vecPhi_pim.push_back(rcFourMom.Phi());
		vecPos_pim.push_back(rcPos);
		vecType_pim.push_back(0);
		vecBottom_pim.push_back(0); 
		vecEff_pim.push_back(TrackEff);


		vecPID_pim.push_back(-321);
		vecMomPID_pim.push_back(-99999);
		vecMotherPID_pim.push_back(-99999);


		//Identify as K-//

		vecTrackID_km.push_back(ik+nPiMinus);
		vecDCA_km.push_back(rcdca);               
		vecPt_km.push_back(rcFourMom.Perp());
		vecEta_km.push_back(rcFourMom.Eta());
		vecPhi_km.push_back(rcFourMom.Phi());
		vecPos_km.push_back(rcPos);
		vecType_km.push_back(0);
		vecBottom_km.push_back(0); 
		vecEff_km.push_back(TrackEff);


		vecPID_km.push_back(-321);
		vecMomPID_km.push_back(-99999);
		vecMotherPID_km.push_back(-99999);

	}



	//Proton+ Loop
	for(int ip=0; ip<nPPlus; ip++) {
		//	cout << "hpiPtWg->Integral() = " << hpiPtWg->Integral() << endl;	
		float pt = hpPtWg->GetRandom();//gRandom->Uniform(0.6,20);
		//	cout << "Pass Inside" << endl;
		float eta = gRandom->Uniform(-1,1);
		float phi = gRandom->Uniform(-PI,PI);
		FourMom.SetPtEtaPhiM(pt, eta , phi, M_PROTON_PLUS);
		float rcdca = dca(FourMom.Vect(), Pos, vertex);


		/*
		   vecMomList_kp.push_back(MomList);
		   vecMomList_pip.push_back(MomList);
		   */

		//		if(abs(FourMom.Eta()) > TrackEtaCut) continue;
		//		if(abs(FourMom.Perp()) < minPtCut) continue;
		//		if(rcdca < DCATrackCut) continue;


		if(SaveCombGen == 1){

			if(AllowGenPID == 0){


				vecTrackID_pipGen.push_back(ip+nPiPlus+nPiPlus);
				vecDCA_pipGen.push_back(rcdca);   	
				vecPt_pipGen.push_back(pt);
				vecEta_pipGen.push_back(eta);
				vecPhi_pipGen.push_back(phi);
				vecPos_pipGen.push_back(Pos);
				vecType_pipGen.push_back(0);
				vecBottom_pipGen.push_back(0); 


				vecPID_pipGen.push_back(2212);
				vecMomPID_pipGen.push_back(-99999);
				vecMotherPID_pipGen.push_back(-99999);

				//		nPiPlusComb = nPiPlusComb + 1;

				vecTrackID_kpGen.push_back(ip+nPiPlus+nPiPlus);
				vecDCA_kpGen.push_back(rcdca);               
				vecPt_kpGen.push_back(pt);
				vecEta_kpGen.push_back(eta);
				vecPhi_kpGen.push_back(phi);
				vecPos_kpGen.push_back(Pos);
				vecType_kpGen.push_back(0);
				vecBottom_kpGen.push_back(0); 


				vecPID_kpGen.push_back(2212);
				vecMomPID_kpGen.push_back(-99999);
				vecMotherPID_kpGen.push_back(-99999);

			}
		}
		TrackEff = fTpcP->Eval(pt);		
		rcFourMom = smearMom(FourMom,fProtonMomResolution);
		if (rcFourMom.Perp()<minPtCut) continue;
		rcPos = smearPosData(2, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
		//	if(rcdca < DCATrackCut) continue;


		//Identify as Pi+//
		vecTrackID_pip.push_back(ip+nPiPlus+nPiPlus);
		vecDCA_pip.push_back(rcdca);               
		vecPt_pip.push_back(rcFourMom.Perp());
		vecEta_pip.push_back(rcFourMom.Eta());
		vecPhi_pip.push_back(rcFourMom.Phi());
		vecPos_pip.push_back(rcPos);
		vecType_pip.push_back(0);
		vecBottom_pip.push_back(0); 
		vecEff_pip.push_back(TrackEff);

		vecPID_pip.push_back(2212);
		vecMomPID_pip.push_back(-99999);
		vecMotherPID_pip.push_back(-99999);

		//Identify as K+//
		vecTrackID_kp.push_back(ip+nPiPlus+nPiPlus);
		vecDCA_kp.push_back(rcdca);               
		vecPt_kp.push_back(rcFourMom.Perp());
		vecEta_kp.push_back(rcFourMom.Eta());
		vecPhi_kp.push_back(rcFourMom.Phi());
		vecPos_kp.push_back(rcPos);
		vecType_kp.push_back(0);
		vecBottom_kp.push_back(0); 
		vecEff_kp.push_back(TrackEff);

		vecPID_kp.push_back(2212);
		vecMomPID_kp.push_back(-99999);
		vecMotherPID_kp.push_back(-99999);

	}


	//Proton- Loop
	for(int ip=0; ip<nPMinus; ip++) {
		//	cout << "hpiPtWg->Integral() = " << hpiPtWg->Integral() << endl;	
		float pt = hpPtWg->GetRandom();//gRandom->Uniform(0.6,20);
		//	cout << "Pass Inside" << endl;
		float eta = gRandom->Uniform(-1,1);
		float phi = gRandom->Uniform(-PI,PI);
		FourMom.SetPtEtaPhiM(pt, eta , phi, M_PROTON_MINUS);
		float rcdca = dca(FourMom.Vect(), Pos, vertex);

		//		if(abs(FourMom.Eta()) > TrackEtaCut) continue;
		//		if(abs(FourMom.Perp()) < minPtCut) continue;
		//if(rcdca < DCATrackCut) continue;

		/*

		   vecMomList_km.push_back(MomList);
		   vecMomList_pim.push_back(MomList);
		   */

		if(SaveCombGen == 1){

			if(AllowGenPID == 0){

				vecTrackID_pimGen.push_back(ip+nPiMinus+nKMinus);
				vecDCA_pimGen.push_back(rcdca);   	
				vecPt_pimGen.push_back(pt);
				vecEta_pimGen.push_back(eta);
				vecPhi_pimGen.push_back(phi);
				vecPos_pimGen.push_back(Pos);
				vecType_pimGen.push_back(0);
				vecBottom_pimGen.push_back(0); 


				vecPID_pimGen.push_back(-2212);
				vecMomPID_pimGen.push_back(-99999);
				vecMotherPID_pimGen.push_back(-99999);


				vecTrackID_kmGen.push_back(ip+nPiMinus+nKMinus);
				vecDCA_kmGen.push_back(rcdca);               
				vecPt_kmGen.push_back(pt);
				vecEta_kmGen.push_back(eta);
				vecPhi_kmGen.push_back(phi);
				vecPos_kmGen.push_back(Pos);
				vecType_kmGen.push_back(0);
				vecBottom_kmGen.push_back(0); 

				vecPID_kmGen.push_back(-2212);
				vecMomPID_kmGen.push_back(-99999);
				vecMotherPID_kmGen.push_back(-99999);

			}
		}
		TrackEff = fTpcP->Eval(pt);				
		rcFourMom = smearMom(FourMom,fProtonMomResolution);
		if (rcFourMom.Perp()<minPtCut) continue;
		rcPos = smearPosData(2, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
		if(rcdca < DCATrackCut) continue;


		//Identify as Pi-//
		vecTrackID_pim.push_back(ip+nPiMinus+nKMinus);
		vecDCA_pim.push_back(rcdca);               
		vecPt_pim.push_back(rcFourMom.Perp());
		vecEta_pim.push_back(rcFourMom.Eta());
		vecPhi_pim.push_back(rcFourMom.Phi());
		vecPos_pim.push_back(rcPos);
		vecType_pim.push_back(0);
		vecBottom_pim.push_back(0); 
		vecEff_pim.push_back(TrackEff);

		vecPID_pim.push_back(-2212);
		vecMomPID_pim.push_back(-99999);
		vecMotherPID_pim.push_back(-99999);

		//Identify as K-//
		vecTrackID_km.push_back(ip+nPiMinus+nKMinus);
		vecDCA_km.push_back(rcdca);               
		vecPt_km.push_back(rcFourMom.Perp());
		vecEta_km.push_back(rcFourMom.Eta());
		vecPhi_km.push_back(rcFourMom.Phi());
		vecPos_km.push_back(rcPos);
		vecType_km.push_back(0);
		vecBottom_km.push_back(0); 
		vecEff_km.push_back(TrackEff);

		vecPID_km.push_back(-2212);
		vecMomPID_km.push_back(-99999);
		vecMotherPID_km.push_back(-99999);


	}




	/*
	   for(int ik=0; ik<nKaons; ik++) {
	   float pt = hkPtWg->GetRandom();//gRandom->Uniform(0.6,20);
	   float eta = gRandom->Uniform(-1,1);
	   float phi = gRandom->Uniform(-PI,PI);
	   FourMom.SetPtEtaPhiM(pt, eta , phi, M_KAON_MINUS);
	   rcFourMom = smearMom(FourMom,fKaonMomResolution);
	   if (rcFourMom.Perp()<minPtCut) continue;
	   rcPos = smearPosData(1, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
	   float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
	   if(rcdca < DCATrackCut) continue;
	   vecDCA_km.push_back(rcdca);               
	   vecPt_km.push_back(rcFourMom.Perp());
	   vecEta_km.push_back(rcFourMom.Eta());
	   vecPhi_km.push_back(rcFourMom.Phi());
	   vecPos_km.push_back(rcPos);
	   vecType_km.push_back(1);
	   vecBottom_km.push_back(0); 

	   }

*/




}


void GenerateDecay(int & EventID){

	//cout << "BRO It WORK BRO" << endl;




	pythia8.next();

	pthat = pythia8.info.pTHat();
	//cout << "Event = " << EventID << "   Total Particles = " << pythia8.event.nFinal() << endl; 

	int totalpar = pythia8.event.nFinal();
	int EventSize =  pythia8.event.size();


	//cout << "Inside PYTHIA: " << "   Pvx  = " <<  vertex.X()<< "   Pvy  = " <<  vertex.Y()<< "   Pvz  = " <<  vertex.Z() << endl;


	nPiPlusDecay = 0;

	for(int i = 0; i < EventSize; i++){
		if(!pythia8.event[i].isFinal()){
			//cout << "SUCK BRO" << endl;
			continue;

		}

		int id = pythia8.event[i].id();
		float pt = pythia8.event[i].pT();
		float eta = pythia8.event[i].eta();
		float phi = pythia8.event[i].phi();



		float vx  = pythia8.event[i].xProd();
		float vy  = pythia8.event[i].yProd();
		float vz  = pythia8.event[i].zProd();

		Pos.SetXYZ(vx*1000,vy*1000,vz*1000);

		int motheridx = pythia8.event[i].mother1();
		int momidx = pythia8.event[i].mother2();

		int motherid = pythia8.event[motheridx].id();
		int momid = pythia8.event[momidx].id();

		/*
		   if(abs(id) == 323|| abs(id) ==  313 || abs(id) == 113 || abs(id) == 213) continue;
		   if(abs(id) == 310) continue; //Reject Kshort
		   if(abs(id) == 221) continue; //Reject eta
		   if(abs(id) == 331) continue; //Reject eta' meson
		   if(abs(id) < 6) continue; //Reject quarks
		   */


		//	cout << "particle i = " << i << "   PDG ID = " << TMath::Abs(pythia8.event[i].id()) << endl;


		//		vector<int> MomList;
		//		MomList = pythia8.event[i].motherList();

		//if (abs(id)!=321 && abs(id)!=211) continue;
		if(id == 211 || id == 321 || id == 2212) {
			/*
			   if(id == 321){
			   cout << "  --------------------------------------------------------------------------------------------  " << endl;
			   cout << "Decay from First (Grand)Mother  = " << mom << "  To Last Mother =  " << mother << endl;
			   cout << "  --------------------------------------------------------------------------------------------  " << endl;

			   }
			   */		
			//cout << "K+/Pi+ Recorded!!" << endl;


			FourMom.SetPtEtaPhiM(pt,eta,phi, M_PION_PLUS);
			float rcdca = dca(FourMom.Vect(), Pos, vertex);
			//		if(rcdca < DCATrackCut) continue;
			//			if(abs(FourMom.Eta()) > TrackEtaCut) continue;
			//			if(abs(FourMom.Perp()) < minPtCut) continue;
			//			if(rcdca < DCATrackCut) continue;
			if(id == 211) nPiPlusDecay = nPiPlusDecay + 1;

			/*
			   vecMomList_kp.push_back(MomList);
			   vecMomList_pip.push_back(MomList);

			   MomList.clear();
			   */


			if(	SaveDecayGen == 1){
				int flagPi = 0;
				int flagK = 0;

				if(AllowGenPID == 0){
					flagPi = 1;
					flagK = 1;
				}


				if(AllowGenPID == 1){
					if(id == 211) flagPi = 1;
					if(id == 321) flagK = 1;
				}


				if(flagPi == 1){
					vecTrackID_pipGen.push_back(nPlus);
					vecDCA_pipGen.push_back(rcdca);               
					vecPt_pipGen.push_back(FourMom.Perp());
					vecEta_pipGen.push_back(FourMom.Eta());
					vecPhi_pipGen.push_back(FourMom.Phi());
					vecPos_pipGen.push_back(Pos);
					if(id == 211 )	vecType_pipGen.push_back(3);
					if(id == 321 )	vecType_pipGen.push_back(4);
					if(id == 2212 )	vecType_pipGen.push_back(4);

					vecBottom_pipGen.push_back(1); 
	

					vecPID_pipGen.push_back(id);
					vecMomPID_pipGen.push_back(momid);
					vecMotherPID_pipGen.push_back(motherid);
				}


				//Identify as K+//
				if(flagK == 1){
					vecTrackID_kpGen.push_back(nPlus);
					vecDCA_kpGen.push_back(rcdca);               
					vecPt_kpGen.push_back(FourMom.Perp());
					vecEta_kpGen.push_back(FourMom.Eta());
					vecPhi_kpGen.push_back(FourMom.Phi());
					vecPos_kpGen.push_back(Pos);
					if(id == 211 )	vecType_kpGen.push_back(4);
					if(id == 321 )	vecType_kpGen.push_back(3);
					if(id == 2212 )	vecType_kpGen.push_back(4);

					vecBottom_kpGen.push_back(1); 

					vecPID_kpGen.push_back(id);
					vecMomPID_kpGen.push_back(momid);
					vecMotherPID_kpGen.push_back(motherid);

				}
			}

			if(id == 211){
				FourMom.SetPtEtaPhiM(pt,eta,phi, M_PION_PLUS);
				TrackEff = fTpcPi->Eval(pt);							
				rcFourMom = smearMom(FourMom,fPionMomResolution);
				rcPos = smearPosData(0, 0, 7, rcFourMom, Pos);

			}

			if(id == 321){
				FourMom.SetPtEtaPhiM(pt,eta,phi, M_KAON_PLUS);
				//	FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
				TrackEff = fTpcK->Eval(pt);				
				rcFourMom = smearMom(FourMom,fKaonMomResolution);
				rcPos = smearPosData(1, 0, 7, rcFourMom, Pos);
			}

			if(id == 2212){
				FourMom.SetPtEtaPhiM(pt,eta,phi, M_PROTON_PLUS);
				//	FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
				TrackEff = fTpcP->Eval(pt);			
				rcFourMom = smearMom(FourMom,fProtonMomResolution);
				rcPos = smearPosData(2, 0, 7, rcFourMom, Pos);
			}


			if (rcFourMom.Perp()<minPtCut) continue;
			//rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
			//			TVector3 rcPosSig(gvx * 10000,gvy * 10000,gvz * 10000);

			rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
			if(rcdca < DCATrackCut) continue;
			//	cout << "Pass 4" << endl;


			//Identify as Pi+//

			vecTrackID_pip.push_back(nPlus);
			vecDCA_pip.push_back(rcdca);               
			vecPt_pip.push_back(rcFourMom.Perp());
			vecEta_pip.push_back(rcFourMom.Eta());
			vecPhi_pip.push_back(rcFourMom.Phi());
			vecPos_pip.push_back(rcPos);
			if(id == 211 )	vecType_pip.push_back(3);
			if(id == 321 )	vecType_pip.push_back(4);
			if(id == 2212 )	vecType_pip.push_back(4);	
			vecBottom_pip.push_back(1); 
			vecEff_pip.push_back(TrackEff);


			vecPID_pip.push_back(id);
			vecMomPID_pip.push_back(momid);
			vecMotherPID_pip.push_back(motherid);



			//Identify as K+//
			vecTrackID_kp.push_back(nPlus);
			vecDCA_kp.push_back(rcdca);               
			vecPt_kp.push_back(rcFourMom.Perp());
			vecEta_kp.push_back(rcFourMom.Eta());
			vecPhi_kp.push_back(rcFourMom.Phi());
			vecPos_kp.push_back(rcPos);
			if(id == 211 )	vecType_kp.push_back(4);
			if(id == 321 )	vecType_kp.push_back(3);
			if(id == 2212 )	vecType_kp.push_back(4);		
			vecBottom_kp.push_back(1); 
			vecEff_kp.push_back(TrackEff);

			vecPID_kp.push_back(id);
			vecMomPID_kp.push_back(momid);
			vecMotherPID_kp.push_back(motherid);


			nPlus = nPlus + 1;
		}



		if(id == -211 || id == -321 || id == -2212) {

			//		cout << "K-/Pi- Recorded!!" << endl;

			FourMom.SetPtEtaPhiM(pt,eta,phi, M_PION_MINUS);
			float rcdca = dca(FourMom.Vect(), Pos, vertex);
			//		if(rcdca < DCATrackCut) continue;
			//			if(abs(FourMom.Eta()) > TrackEtaCut) continue;
			//			if(abs(FourMom.Perp()) < minPtCut) continue;
			//			if(rcdca < DCATrackCut) continue;

			/*
			   vecMomList_km.push_back(MomList);
			   vecMomList_pim.push_back(MomList);

			   MomList.clear();
			   */
			if(	SaveDecayGen == 1){

				int flagPi = 0;
				int flagK = 0;

				if(AllowGenPID == 0){
					flagPi = 1;
					flagK = 1;
				}


				if(AllowGenPID == 1){
					if(id == -211) flagPi = 1;
					if(id == -321) flagK = 1;
				}



				if(flagPi == 1){

					vecTrackID_pimGen.push_back(nMinus);
					vecDCA_pimGen.push_back(rcdca);               
					vecPt_pimGen.push_back(FourMom.Perp());
					vecEta_pimGen.push_back(FourMom.Eta());
					vecPhi_pimGen.push_back(FourMom.Phi());
					vecPos_pimGen.push_back(Pos);
					if(id == -211 )	vecType_pimGen.push_back(3);
					if(id == -321 )	vecType_pimGen.push_back(4);
					if(id == -2212 ) vecType_pimGen.push_back(4);		
					vecBottom_pimGen.push_back(1); 
					//			cout << "Pass 2 " << endl;

					vecPID_pimGen.push_back(id);
					vecMomPID_pimGen.push_back(momid);
					vecMotherPID_pimGen.push_back(motherid);
				}

				//Identify as K-//


				if(flagK == 1){

					vecTrackID_kmGen.push_back(nMinus);
					vecDCA_kmGen.push_back(rcdca);               
					vecPt_kmGen.push_back(FourMom.Perp());
					vecEta_kmGen.push_back(FourMom.Eta());
					vecPhi_kmGen.push_back(FourMom.Phi());
					vecPos_kmGen.push_back(Pos);
					if(id == -211 )	vecType_kmGen.push_back(4);
					if(id == -321 )	vecType_kmGen.push_back(3);
					if(id == -2212 )	vecType_kmGen.push_back(4);

					vecBottom_kmGen.push_back(1); 

					vecPID_kmGen.push_back(id);
					vecMomPID_kmGen.push_back(momid);
					vecMotherPID_kmGen.push_back(motherid);

				}

			}
			if(id == -211){

				FourMom.SetPtEtaPhiM(pt,eta,phi, M_PION_MINUS);
				TrackEff = fTpcPi->Eval(pt);											
				rcFourMom = smearMom(FourMom,fPionMomResolution);
				rcPos = smearPosData(0, 0, 7, rcFourMom, Pos);
				//			cout << "Pass Pos Pi" << endl;

			}

			if(id == -321){

				FourMom.SetPtEtaPhiM(pt,eta,phi, M_KAON_MINUS);
				TrackEff = fTpcK->Eval(pt);							
				rcFourMom = smearMom(FourMom,fKaonMomResolution);
				rcPos = smearPosData(1, 0, 7, rcFourMom, Pos);

			}

			if(id == -2212){

				FourMom.SetPtEtaPhiM(pt,eta,phi, M_PROTON_MINUS);
				TrackEff = fTpcP->Eval(pt);											
				rcFourMom = smearMom(FourMom,fProtonMomResolution);
				rcPos = smearPosData(2, 0, 7, rcFourMom, Pos);

			}

			//			cout << "Pass 0 " << endl;

			if (rcFourMom.Perp()<minPtCut) continue;
			//rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
			//		TVector3 rcPosSig(gvx * 10000,gvy * 10000,gvz * 10000);
			//		rcPos = smearPosData(0, gvz, 8, rcFourMom, rcPosSig);

			rcdca = dca(rcFourMom.Vect(), rcPos, vertex);


			//			cout << "Pass 1 " << endl;
			//Identify as Pi-//

			vecTrackID_pim.push_back(nMinus);
			vecDCA_pim.push_back(rcdca);               
			vecPt_pim.push_back(rcFourMom.Perp());
			vecEta_pim.push_back(rcFourMom.Eta());
			vecPhi_pim.push_back(rcFourMom.Phi());
			vecPos_pim.push_back(rcPos);
			if(id == -211 )	vecType_pim.push_back(3);
			if(id == -321 )	vecType_pim.push_back(4);
			if(id == -2212 )	vecType_pim.push_back(4);

			vecBottom_pim.push_back(1); 
			//			cout << "Pass 2 " << endl;
			vecEff_pim.push_back(TrackEff);

			vecPID_pim.push_back(id);
			vecMomPID_pim.push_back(momid);
			vecMotherPID_pim.push_back(motherid);
			//Identify as K-//

			vecTrackID_km.push_back(nMinus);
			vecDCA_km.push_back(rcdca);               
			vecPt_km.push_back(rcFourMom.Perp());
			vecEta_km.push_back(rcFourMom.Eta());
			vecPhi_km.push_back(rcFourMom.Phi());
			vecPos_km.push_back(rcPos);
			if(id == -211 )	vecType_km.push_back(4);
			if(id == -321 )	vecType_km.push_back(3);
			if(id == -2212 )	vecType_km.push_back(4);

			vecBottom_km.push_back(1); 
			vecEff_km.push_back(TrackEff);
			
			//		cout << "Pass 3 " << endl;
			vecPID_km.push_back(id);
			vecMomPID_km.push_back(momid);
			vecMotherPID_km.push_back(motherid);

			nMinus = nMinus + 1;
		}





	}





	//	return;


}


/*

   void GenerateDecay(  int & EventID){


   if(EventID%10000 == 0)	cout << "Now Inside Fucking EvtGen for Event = " << EventID  << endl;

   gRandom->SetSeed();


//pythia8 = new PHPythia8();

//cout << "Pass 1" << endl;

pythia8->Clear();
pythia8->Make();
int npar = pythia8->Event()->GetNumberOfParticles();
StarGenEvent* evt = pythia8->Event();	
//	pythia8->init();
//	pythia8->next();


//int npar = pythia8->Event()->nFinal();

//	evt = pythia8->Event();	

//		cout << "Pass 2" << endl;

//cout << "NParticle = " << npar << endl;




for (int ipar=0;ipar<npar;ipar++){

StarGenParticle* p =  ((*evt)[ipar]);
int id = p->GetId();

int momidx = p->GetFirstMother();
int momid = ((*evt)[momidx])->GetId();

//	int mammyidx =  p->GetSecondMother();
//	int mammyid = ((*evt)[mammyidx])->GetId();

int motheridx = p->GetLastMother();
int motherpdg = ((*evt)[motheridx])->GetId();

int mammyidx = ((*evt)[motheridx])->GetFirstDaughter();
int mammyid = ((*evt)[mammyidx])->GetId();


//		if(abs(motherpdg) != 531) continue;

if(doCut){


if (abs(id)!=321 && abs(id)!=211) continue;
if (fabs(p->momentum().PseudoRapidity())>=1.1) continue;

if(abs(id) == 321 && abs(momid) != 333) continue;
if(abs(id) == 321 && abs(mammyidx) != 431) continue;
if(abs(id) == 221 && abs(momid) != 531 && abs(momid) != 431 ) continue;

}

Pos.SetXYZ(p->GetVx()*1000,p->GetVy()*1000,p->GetVz()*1000);

//	cout << "Pass 3" << endl;

if(id == 211 || id == 321) {

//cout << "K+/Pi+ Recorded!!" << endl;

FourMom = p->momentum();

if(id == 211){
	//FourMom.SetXYZM(gpx, gpy , gpz, M_PION_MINUS);
	rcFourMom = smearMom(FourMom,fPionMomResolution);
	rcPos = smearPosData(0, 0, 7, rcFourMom, Pos);

}

if(id == 321){
	//	FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
	rcFourMom = smearMom(FourMom,fKaonMomResolution);
	rcPos = smearPosData(1, 0, 7, rcFourMom, Pos);
}

if (rcFourMom.Perp()<minPtCut) continue;
//rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
//			TVector3 rcPosSig(gvx * 10000,gvy * 10000,gvz * 10000);

float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
if(rcdca < DCATrackCut) continue;
//	cout << "Pass 4" << endl;


//Identify as Pi+//

vecTrackID_pip.push_back(nPlus);
vecDCA_pip.push_back(rcdca);               
vecPt_pip.push_back(rcFourMom.Perp());
vecEta_pip.push_back(rcFourMom.Eta());
vecPhi_pip.push_back(rcFourMom.Phi());
vecPos_pip.push_back(rcPos);
if(id == 211 )	vecType_pip.push_back(3);
if(id == 321 )	vecType_pip.push_back(4);
vecBottom_pip.push_back(1); 


//Identify as K+//
vecTrackID_kp.push_back(nPlus);
vecDCA_kp.push_back(rcdca);               
vecPt_kp.push_back(rcFourMom.Perp());
vecEta_kp.push_back(rcFourMom.Eta());
vecPhi_kp.push_back(rcFourMom.Phi());
vecPos_kp.push_back(rcPos);
if(id == 211 )	vecType_kp.push_back(4);
if(id == 321 )	vecType_kp.push_back(3);
vecBottom_kp.push_back(1); 


nPlus = nPlus + 1;
}




//	cout << "Pass 5" << endl;

if(id == -211 || id == -321) {

	//	cout << "K-/Pi- Recorded!!" << endl;


	FourMom = p->momentum();

	//cout << "Pass 4 Momentum " << endl;

	if(id == -211){
		//	FourMom.SetXYZM(gpx, gpy , gpz, M_PION_MINUS);
		rcFourMom = smearMom(FourMom,fPionMomResolution);
		//			cout << "Pass Reso Pi" << endl;

		//			cout << "rcFourMom.Pt() = " << rcFourMom.Pt() << "   Pos = " << Pos.X() << endl;

		rcPos = smearPosData(0, 0, 7, rcFourMom, Pos);
		//			cout << "Pass Pos Pi" << endl;

	}

	if(id == -321){
		//	FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
		rcFourMom = smearMom(FourMom,fKaonMomResolution);
		//		cout << "Pass Reso K" << endl;
		rcPos = smearPosData(1, 0, 7, rcFourMom, Pos);
		//		cout << "Pass Pos K" << endl;

	}

	//			cout << "Pass 0 " << endl;

	if (rcFourMom.Perp()<minPtCut) continue;
	//rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
	//		TVector3 rcPosSig(gvx * 10000,gvy * 10000,gvz * 10000);
	//		rcPos = smearPosData(0, gvz, 8, rcFourMom, rcPosSig);

	float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
	if(rcdca < DCATrackCut) continue;

	//			cout << "Pass 1 " << endl;
	//Identify as Pi-//

	vecTrackID_pim.push_back(nMinus);
	vecDCA_pim.push_back(rcdca);               
	vecPt_pim.push_back(rcFourMom.Perp());
	vecEta_pim.push_back(rcFourMom.Eta());
	vecPhi_pim.push_back(rcFourMom.Phi());
	vecPos_pim.push_back(rcPos);
	if(id == -211 )	vecType_pim.push_back(3);
	if(id == -321 )	vecType_pim.push_back(4);
	vecBottom_pim.push_back(1); 
	//			cout << "Pass 2 " << endl;


	//Identify as K-//

	vecTrackID_km.push_back(nMinus);
	vecDCA_km.push_back(rcdca);               
	vecPt_km.push_back(rcFourMom.Perp());
	vecEta_km.push_back(rcFourMom.Eta());
	vecPhi_km.push_back(rcFourMom.Phi());
	vecPos_km.push_back(rcPos);
	if(id == -211 )	vecType_km.push_back(4);
	if(id == -321 )	vecType_km.push_back(3);
	vecBottom_km.push_back(1); 
	//		cout << "Pass 3 " << endl;


	nMinus = nMinus + 1;
}



//	cout << "Pass 6" << endl;


}




}

*/

void GenerateEvtGen(int & EventID){


	cout << "Now Doing FUCKING EVT GEN BRO on FUCKING EVENT = " << EventID  << endl;

	if(!RunUncorrBkgd){
		nPlus = 0;
		nMinus = 0;
	}


	TClonesArray daughters("TParticle", 10);
	TLorentzVector* BsVecSig = new TLorentzVector;

	getKinematics(*BsVecSig, BsMassPDG);
		
	myEvtGenDecayer->Decay(531, BsVecSig);
	myEvtGenDecayer->ImportParticles(&daughters);

	nTrk = daughters.GetEntriesFast();

	//	vector<int> MomList;

	nPiPlusSig = 0;

	
	
	for (int iTrk = 0; iTrk < nTrk; ++iTrk)
	{

		TParticle* ptl0 = (TParticle*)daughters.At(iTrk);
		int id = ptl0->GetPdgCode();
		//cout << "id = " << id << endl;
		//cout << "iTrk = " << iTrk << "    GetMother = " << ptl0->GetMother(2) << endl;
		if(id == 431) break;
		if(id == 431) PiPFromDs = 1;
		
		
		if(id == -431) PiPFromDs = 0;
	
		ptl0->Momentum(FourMom);
		TVector3 rcPosSig(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.);

	
		if(id == 211 || id == 321 || id == 2212) {

			if(id == 211 )nPiPlusSig = nPiPlusSig + 1;


			/*
			   if(id == 321){

			   MomList.push_back(531);
			   MomList.push_back(431);
			   MomList.push_back(333);

			   }

			   if(id == 211){

			   MomList.push_back(531);
			   MomList.push_back(431);

			   }

			   if(id == 2212){
			   MomList.push_back(0);
			   MomList.push_back(0);
			   MomList.push_back(0);

			   }

			   MomList.clear();
			   */

			//cout << "K+/Pi+ Recorded!!" << endl;

			//	FourMom = p->momentum();
			float rcdca = dca(FourMom.Vect(), rcPosSig, vertex);

			//Apply Track cuts//

			//			if(abs(FourMom.Eta()) > TrackEtaCut) continue;
			//			if(abs(FourMom.Perp()) < minPtCut) continue;
			//			if(rcdca < DCATrackCut) continue;

			/*
			   vecMomList_kp.push_back(MomList);
			   vecMomList_pip.push_back(MomList);
			   */

			int flagPi = 0;
			int flagK = 0;

			if(AllowGenPID == 0){
				flagPi = 1;
				flagK = 1;
			}


			if(AllowGenPID == 1){
				if(id == 211) flagPi = 1;
				if(id == 321) flagK = 1;
			}

			if(	SaveSigGen == 1){


				if(flagPi == 1){
					vecTrackID_pipGen.push_back(nPlus);
					vecDCA_pipGen.push_back(rcdca);               
					vecPt_pipGen.push_back(FourMom.Perp());
					vecEta_pipGen.push_back(FourMom.Eta());
					vecPhi_pipGen.push_back(FourMom.Phi());
					vecPos_pipGen.push_back(rcPosSig);
					if(id == 211)	vecType_pipGen.push_back(1);
					if(id == 321 )	vecType_pipGen.push_back(2);
					if(id == 2212 )	vecType_pipGen.push_back(2);

					vecPID_pipGen.push_back(id);
					vecMomPID_pipGen.push_back(99999);
					vecMotherPID_pipGen.push_back(99999);

					vecBottom_pipGen.push_back(1); 

				}


				if(flagK == 1){

					vecTrackID_kpGen.push_back(nPlus);
					vecDCA_kpGen.push_back(rcdca);               
					vecPt_kpGen.push_back(FourMom.Perp());
					vecEta_kpGen.push_back(FourMom.Eta());
					vecPhi_kpGen.push_back(FourMom.Phi());
					vecPos_kpGen.push_back(rcPosSig);
					if(id == 211 )	vecType_kpGen.push_back(2);
					if(id == 321 )	vecType_kpGen.push_back(1);
					if(id == 2212 )	vecType_kpGen.push_back(2);

					vecPID_kpGen.push_back(id);
					vecMomPID_kpGen.push_back(99999);
					vecMotherPID_kpGen.push_back(99999);


					vecBottom_kpGen.push_back(1); 

				}
			}

			if(id == 211){
				//FourMom.SetXYZM(gpx, gpy , gpz, M_PION_MINUS);
				TrackEff = fTpcPi->Eval(FourMom.Perp());											
				rcFourMom = smearMom(FourMom,fPionMomResolution);
				rcPos = smearPosData(0, 0, 7, rcFourMom, rcPosSig);

			}

			if(id == 321){
				//	FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
				TrackEff = fTpcK->Eval(FourMom.Perp());															
				rcFourMom = smearMom(FourMom,fKaonMomResolution);
				rcPos = smearPosData(1, 0, 7, rcFourMom, rcPosSig);
			}

			if(id == 2212){
				//	FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
				TrackEff = fTpcP->Eval(FourMom.Perp());																		
				rcFourMom = smearMom(FourMom,fProtonMomResolution);
				rcPos = smearPosData(2, 0, 7, rcFourMom, rcPosSig);
			}


			if (rcFourMom.Perp()<minPtCut) continue;
			//rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
			//			TVector3 rcPosSig(gvx * 10000,gvy * 10000,gvz * 10000);

			rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
				
			if(rcdca < DCATrackCut) continue;
			//	cout << "Pass 4" << endl;


			//Identify as Pi+//

			vecTrackID_pip.push_back(nPlus);
			vecDCA_pip.push_back(rcdca);               
			vecPt_pip.push_back(rcFourMom.Perp());
			vecEta_pip.push_back(rcFourMom.Eta());
			vecPhi_pip.push_back(rcFourMom.Phi());
			vecPos_pip.push_back(rcPos);
			if(id == 211)	vecType_pip.push_back(1);
			if(id == 321 )	vecType_pip.push_back(2);
			if(id == 2212 )	vecType_pip.push_back(2);

			vecBottom_pip.push_back(1); 
			vecEff_pip.push_back(TrackEff);

			vecPID_pip.push_back(id);
			vecMomPID_pip.push_back(99999);
			vecMotherPID_pip.push_back(99999);


			//Identify as K+//
			vecTrackID_kp.push_back(nPlus);
			vecDCA_kp.push_back(rcdca);               
			vecPt_kp.push_back(rcFourMom.Perp());
			vecEta_kp.push_back(rcFourMom.Eta());
			vecPhi_kp.push_back(rcFourMom.Phi());
			vecPos_kp.push_back(rcPos);
			if(id == 211 )	vecType_kp.push_back(2);
			if(id == 321 )	vecType_kp.push_back(1);
			if(id == 2212 )	vecType_kp.push_back(2);	
			vecBottom_kp.push_back(1); 
			vecEff_kp.push_back(TrackEff);

			vecPID_kp.push_back(id);
			vecMomPID_kp.push_back(99999);
			vecMotherPID_kp.push_back(99999);



			nPlus = nPlus + 1;
		}




		//	cout << "Pass 5" << endl;

		if(id == -211 || id == -321 || id == -2212) {

			//	cout << "K-/Pi- Recorded!!" << endl;


			/*
			   if(id == -321){

			   MomList.push_back(531);
			   MomList.push_back(431);
			   MomList.push_back(333);

			   }

			   if(id == -211){

			   MomList.push_back(531);
			   MomList.push_back(431);

			   }

			   if(id == -2212){
			   MomList.push_back(0);
			   MomList.push_back(0);
			   MomList.push_back(0);

			   }



			   MomList.clear();
			   */


			float rcdca = dca(FourMom.Vect(), rcPosSig, vertex);

			//			if(abs(FourMom.Eta()) > TrackEtaCut) continue;
			//			if(abs(FourMom.Perp()) < minPtCut) continue;
			//			if(rcdca < DCATrackCut) continue;

			/*
			   vecMomList_km.push_back(MomList);
			   vecMomList_pim.push_back(MomList);
			   */

			if(	SaveSigGen == 1){

				int flagPi = 0;
				int flagK = 0;

				if(AllowGenPID == 0){
					flagPi = 1;
					flagK = 1;
				}


				if(AllowGenPID == 1){
					if(id == -211) flagPi = 1;
					if(id == -321) flagK = 1;
				}


				if(flagPi == 1){
					vecTrackID_pimGen.push_back(nMinus);
					vecDCA_pimGen.push_back(rcdca);               
					vecPt_pimGen.push_back(FourMom.Perp());
					vecEta_pimGen.push_back(FourMom.Eta());
					vecPhi_pimGen.push_back(FourMom.Phi());
					vecPos_pimGen.push_back(rcPosSig);
					if(id == -211)	vecType_pimGen.push_back(1);
					if(id == -321 )	vecType_pimGen.push_back(2);
					if(id == -2212 )	vecType_pimGen.push_back(2);

					vecPID_pimGen.push_back(id);
					vecMomPID_pimGen.push_back(99999);
					vecMotherPID_pimGen.push_back(99999);

					vecBottom_pimGen.push_back(1); 
				}



				if(flagK == 1){

					vecTrackID_kmGen.push_back(nMinus);
					vecDCA_kmGen.push_back(rcdca);               
					vecPt_kmGen.push_back(FourMom.Perp());
					vecEta_kmGen.push_back(FourMom.Eta());
					vecPhi_kmGen.push_back(FourMom.Phi());
					vecPos_kmGen.push_back(rcPosSig);
					if(id == -211 )	vecType_kmGen.push_back(2);
					if(id == -321 )	vecType_kmGen.push_back(1);
					if(id == -2212 ) vecType_kmGen.push_back(2);	
					vecBottom_kmGen.push_back(1); 

					vecPID_kmGen.push_back(id);
					vecMomPID_kmGen.push_back(99999);
					vecMotherPID_kmGen.push_back(99999);
				}

			}
			//	FourMom = p->momentum();

			//cout << "Pass 4 Momentum " << endl;

			if(id == -211){
				//	FourMom.SetXYZM(gpx, gpy , gpz, M_PION_MINUS);
			
				TrackEff = fTpcPi->Eval(FourMom.Perp());											
				rcFourMom = smearMom(FourMom,fPionMomResolution);

				rcPos = smearPosData(0, 0, 7, rcFourMom, rcPosSig);
				//			cout << "Pass Pos Pi" << endl;
				
			}

			if(id == -321){
				//	FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
				TrackEff = fTpcK->Eval(FourMom.Perp());															
				rcFourMom = smearMom(FourMom,fKaonMomResolution);
				//		cout << "Pass Reso K" << endl;
				rcPos = smearPosData(1, 0, 7, rcFourMom, rcPosSig);
				//		cout << "Pass Pos K" << endl;

			}

			if(id == -2212){
				//	FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
				TrackEff = fTpcP->Eval(FourMom.Perp());															
				rcFourMom = smearMom(FourMom,fProtonMomResolution);
				//		cout << "Pass Reso K" << endl;
				rcPos = smearPosData(2, 0, 7, rcFourMom, rcPosSig);
				//		cout << "Pass Pos K" << endl;

			}


			if (rcFourMom.Perp()<minPtCut) continue;



			rcdca = dca(rcFourMom.Vect(), rcPos, vertex);


			if(rcdca < DCATrackCut) continue;


			vecTrackID_pim.push_back(nMinus);
			vecDCA_pim.push_back(rcdca);               
			vecPt_pim.push_back(rcFourMom.Perp());
			vecEta_pim.push_back(rcFourMom.Eta());
			vecPhi_pim.push_back(rcFourMom.Phi());
			vecPos_pim.push_back(rcPos);
			if(id == -211 )	vecType_pim.push_back(1);
			if(id == -321 )	vecType_pim.push_back(2);
			if(id == -2212 )	vecType_pim.push_back(2);	
			vecBottom_pim.push_back(1); 
			vecEff_pim.push_back(TrackEff);

			vecPID_pim.push_back(id);
			vecMomPID_pim.push_back(99999);
			vecMotherPID_pim.push_back(99999);


			//Identify as K-//

			vecTrackID_km.push_back(nMinus);
			vecDCA_km.push_back(rcdca);               
			vecPt_km.push_back(rcFourMom.Perp());
			vecEta_km.push_back(rcFourMom.Eta());
			vecPhi_km.push_back(rcFourMom.Phi());
			vecPos_km.push_back(rcPos);
			if(id == -211 )	vecType_km.push_back(2);
			if(id == -321 )	vecType_km.push_back(1);
			if(id == -2212 )	vecType_km.push_back(2);	
			vecBottom_km.push_back(1); 
			vecEff_km.push_back(TrackEff);


			vecPID_km.push_back(id);
			vecMomPID_km.push_back(99999);
			vecMotherPID_km.push_back(99999);


			nMinus = nMinus + 1;
		}


	}


	daughters.Clear();


}


/*


   void SigEvttoEntry(){


   cout << "Now Fucking Deciding the Entries" << endl;

   TFile * SignalFile = new TFile("G4sPHENIX.root_g4svtx_eval.root");
   SignalFile->cd();

   TTree * ntp_track = (TTree *) SignalFile->Get("ntp_track");
   ntp_track->SetBranchAddress("event",&EventNow);


   TotalCand = ntp_track->GetEntries();

//	EventComp = int(EventNow);

int EventPre = -1;

cout << "TotalCand = " << TotalCand << endl;

for(int i = 0 ; i < TotalCand; i++){

ntp_track->GetEntry(i);

//cout <<  "EventNow = " << EventNow << endl;

if(EventPre != EventNow)  EvtToSigEntry.push_back(i);

EventPre = EventNow;

}


for(int i = 0; i < EvtToSigEntry.size(); i ++){

//	cout << "EntryBoundaries " << EvtToSigEntry[i] << endl;
}


}


void GetGen(int & EventID){


ntp_gtrack->GetEntry(EventID);

gEvent = gInEvent;
gBvtxx = gInBvtxx;
gBvtxy = gInBvtxy;
gBvtxz = gInBvtxz;

gBpx = gInBpx;
gBpy = gInBpy;
gBpz = gInBpz;
gBpt = gInBpt;

nt_gen->Fill();

}


void GetSignal(int & EventID){



if(!RunBackground){
nPlus = 0;
nMinus = 0;
}

EntryMin = EvtToSigEntry[EventID];
EntryMax = EvtToSigEntry[EventID+1];

cout << "Signal Event = " << EventID << "    EntryMin = " << EntryMin << "   EntryMax = " << EntryMax << endl;

cout << "Now Getting the Fucking Signal From Existing Codes" << endl;


for(int i = EntryMin; i < EntryMax; i++){

	ntp_track->GetEntry(i);





	if(gflavor == 211 || gflavor == 321) {

		if(gflavor == 211){
			FourMom.SetXYZM(gpx, gpy , gpz, M_PION_MINUS);
			rcFourMom = smearMom(FourMom,fPionMomResolution);
		}

		if(gflavor == 321){
			FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
			rcFourMom = smearMom(FourMom,fKaonMomResolution);
		}

		if (rcFourMom.Perp()<minPtCut) continue;
		//rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		TVector3 rcPosSig(gvx,gvy ,gvz );
		rcPos = smearPosData(0, gvz, 8, rcFourMom, rcPosSig);

		float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
		if(rcdca < DCATrackCut) continue;


		//Identify as Pi+//

		vecTrackID_pip.push_back(nPlus);
		vecDCA_pip.push_back(rcdca);               
		vecPt_pip.push_back(rcFourMom.Perp());
		vecEta_pip.push_back(rcFourMom.Eta());
		vecPhi_pip.push_back(rcFourMom.Phi());
		vecPos_pip.push_back(rcPos);
		if(gflavor == 211 )	vecType_pip.push_back(1);
		if(gflavor == 321 )	vecType_pip.push_back(0);
		vecBottom_pip.push_back(1); 


		//Identify as K+//
		vecTrackID_kp.push_back(nPlus);
		vecDCA_kp.push_back(rcdca);               
		vecPt_kp.push_back(rcFourMom.Perp());
		vecEta_kp.push_back(rcFourMom.Eta());
		vecPhi_kp.push_back(rcFourMom.Phi());
		vecPos_kp.push_back(rcPos);
		if(gflavor == 211 )	vecType_kp.push_back(0);
		if(gflavor == 321 )	vecType_kp.push_back(1);
		vecBottom_kp.push_back(1); 


		nPlus = nPlus + 1;
	}





	if(gflavor == -211 || gflavor == -321) {

		if(gflavor == -211){
			FourMom.SetXYZM(gpx, gpy , gpz, M_PION_MINUS);
			rcFourMom = smearMom(FourMom,fPionMomResolution);
		}

		if(gflavor == -321){
			FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
			rcFourMom = smearMom(FourMom,fKaonMomResolution);
		}

		if (rcFourMom.Perp()<minPtCut) continue;
		//rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		TVector3 rcPosSig(gvx ,gvy ,gvz );
		rcPos = smearPosData(0, gvz, 8, rcFourMom, rcPosSig);

		float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
		if(rcdca < DCATrackCut) continue;


		//Identify as Pi-//

		vecTrackID_pim.push_back(nMinus);
		vecDCA_pim.push_back(rcdca);               
		vecPt_pim.push_back(rcFourMom.Perp());
		vecEta_pim.push_back(rcFourMom.Eta());
		vecPhi_pim.push_back(rcFourMom.Phi());
		vecPos_pim.push_back(rcPos);
		if(gflavor == -211 )	vecType_pim.push_back(1);
		if(gflavor == -321 )	vecType_pim.push_back(0);
		vecBottom_pim.push_back(1); 


		//Identify as K-//

		vecTrackID_km.push_back(nMinus);
		vecDCA_km.push_back(rcdca);               
		vecPt_km.push_back(rcFourMom.Perp());
		vecEta_km.push_back(rcFourMom.Eta());
		vecPhi_km.push_back(rcFourMom.Phi());
		vecPos_km.push_back(rcPos);
		if(gflavor == -211 )	vecType_km.push_back(0);
		if(gflavor == -321 )	vecType_km.push_back(1);
		vecBottom_km.push_back(1); 


		nMinus = nMinus + 1;
	}





}


cout << "DONE SIGNAL" << endl;




}

*/

void initPYTHIA8(){


	pythia8.readString("HardQCD:all = on");
	pythia8.readString("Random:setSeed = on");
	pythia8.readString("Random:seed = 0");
	pythia8.readFile("phpythia8.cfg");
	pythia8.readString("HardQCD:qqbar2bbbar=on");
	pythia8.readString("HardQCD:gg2bbbar=on");
	pythia8.readString("HardQCD:qqbar2ccbar=off");
	pythia8.readString("HardQCD:gg2ccbar=off");


	pythia8.init();


	cout << "PYTHIA8 INIT" << endl;

	/*
	   pythia8 = new StarPythia8();
	   cout << "Pass 4" << endl;

	   pythia8->SetFrame("CMS", 200.0);
	   pythia8->SetBlue("proton");
	   cout << "Pass 5" << endl;

	   pythia8->SetYell("proton");            
	// pythia8->ReadFile("star_hf_tune_v1.1.cmnd");
	pythia8->Set("HardQCD:qqbar2bbbar=on");
	pythia8->Set("HardQCD:gg2bbbar=on");
	pythia8->Set("HardQCD:qqbar2ccbar=off");
	pythia8->Set("HardQCD:gg2ccbar=off");
	pythia8->Init();
	cout << "Pass 6" << endl;
	*/

	//gSystem->Load("libPHPythia8.so");

	/*
	   char *charPath = getenv("PYTHIA8");
	   if (!charPath)
	   {
	   cout << "PHPythia8::Could not find $PYTHIA8 path!" << endl;
	   return;
	   }

	   std::string thePath(charPath);
	   thePath += "/xmldoc/";
	//Pythia8::Pythia m_Pythia8 = new Pythia8::Pythia(thePath.c_str());
	*/





}


void initEvtGen()
{
	cout << "Now Loading the Fucking EVT GEN with PHENIX SOFTWARE" << endl;


	gSystem->Load("EvtGen/lib/libEvtGenExternal.so");
	gSystem->Load("EvtGen/lib/libEvtGen.so");


	cout << "DONE LOADING EVTGEN::: GO	!!!!" << endl;	
	EvtRandomEngine* eng = 0;
	eng = new EvtSimpleRandomEngine();
	EvtRandom::setRandomEngine((EvtRandomEngine*)eng);
	EvtAbsRadCorr* radCorrEngine = 0;
	std::list<EvtDecayBase*> extraModels;
	cout << "Pass 8" << endl;

	EvtExternalGenList genList;
	cout << "Pass 9" << endl;
	radCorrEngine = genList.getPhotosModel();
	cout << "Pass 10" << endl;

	extraModels = genList.getListOfModels();

	TString Decay_DEC = "InputDECAYFiles/DECAY.DEC";
	//TString Decay_DEC = "Bs.Mypythia0.DEC";

	TString Evt_pdl = "InputDECAYFiles/evt.pdl";
	EvtGen *myGenerator = new EvtGen(Decay_DEC, Evt_pdl, (EvtRandomEngine*)eng, radCorrEngine, &extraModels);
	myEvtGenDecayer = new PHEvtGenDecayer(myGenerator);
	cout << "Setting Decay Table" << endl;
	myEvtGenDecayer->SetDecayTable("InputDECAYFiles/Bs.KKPiPi.DEC");
	cout << "Pass 10" << endl;



}



void getKinematics(TLorentzVector& b, double const mass)
{
	//float const pt = gRandom->Uniform(0, 20);
	//float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
	
	double pt = 0;
	double y = 0;
	
	FONLL2DMap->GetRandom2(pt,y);
	

	float const phi = TMath::TwoPi() * gRandom->Rndm();

	primaryPt = pt;
	primaryY = y;
	primaryPhi = phi;

	float const mT = sqrt(mass * mass + pt * pt);
	float const pz = mT * sinh(y);
	float const E = mT * cosh(y);

	b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
	TVector3 posDiff = pos - vertex;
	return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
	//return fabs(p.Cross(posDiff)/p.Mag());
}

float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
	TVector3 posDiff = pos - vertex;
	float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;

	return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}

float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
	TVector3 newPos(pos);
	newPos.SetZ(0);

	TVector3 newP(p);
	newP.SetZ(0);

	TVector3 newVertex(vertex);
	newVertex.SetZ(0);

	TVector3 posDiff = newPos - newVertex;
	float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
	return sign * newP.Cross(posDiff.Cross(newP)).Unit().Dot(posDiff);
}

float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
	TVector3 posDiff = pos - vertex;
	if (sin(p.Theta()) == 0) return 0;
	else return (-posDiff.x() * cos(p.Phi()) - posDiff.y() * sin(p.Phi())) * cos(p.Theta()) / sin(p.Theta()) + posDiff.z();
}

float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0)
{
	TVector3 posDiff = pos2 - pos1;
	TVector3 pu1 = p1.Unit();
	TVector3 pu2 = p2.Unit();
	double pu1Pu2 = pu1.Dot(pu2);
	double g = posDiff.Dot(pu1);
	double k = posDiff.Dot(pu2);
	double s2 = (k - pu1Pu2 * g) / (pu1Pu2 * pu1Pu2 - 1.);
	double s1 = g + s2 * pu1Pu2;
	TVector3 posDca1 = pos1 + pu1 * s1;
	TVector3 posDca2 = pos2 + pu2 * s2;
	v0 = 0.5 * (posDca1 + posDca2);
	return (posDca1 - posDca2).Mag();
}

TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution)
{
	float const pt = b.Perp();
	float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt));

	TLorentzVector sMom;
	sMom.SetXYZM(sPt * cos(b.Phi()), sPt * sin(b.Phi()), sPt * sinh(b.PseudoRapidity()), b.M());
	return sMom;
}

TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos)
{
	float thetaMCS = 13.6 / mom.Beta() / rMom.P() / 1000 * sqrt(pxlLayer1Thickness / fabs(sin(mom.Theta())));
	float sigmaMCS = thetaMCS * 28000 / fabs(sin(mom.Theta()));
	float sigmaPos = sqrt(pow(sigmaMCS, 2) + pow(sigmaPos0, 2));

	return TVector3(gRandom->Gaus(pos.X(), sigmaPos), gRandom->Gaus(pos.Y(), sigmaPos), gRandom->Gaus(pos.Z(), sigmaPos));
}

TVector3 smearPosData(int const iParticleIndex, double const vz, int cent, TLorentzVector const& rMom, TVector3 const& pos)
{
	int const iEtaIndex = getEtaIndexDca(rMom.PseudoRapidity());
	int const iVzIndex = getVzIndexDca(vz);
	// int const iPhiIndex = getPhiIndexDca(rMom.Phi());
	int const iPtIndex = getPtIndexDca(rMom.Perp());

	//	cout << "iPtIndex = " << iPtIndex << endl;
	//	cout << "iParticleIndex = " << iParticleIndex << endl;
	double sigmaPosZ = 0;
	double sigmaPosXY = 0;

	//	cout << "Pass Inside 1 " << endl;

	if (cent == 8) cent = 7;
	//All the centrality position smear was based on 0-10% centrality input, so here the cent==0
	//	cout << "Pass Inside 1.5 " << endl;

	// h2Dca[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom2(sigmaPosXY,sigmaPosZ);
	h2Dca[iParticleIndex][iPtIndex]->GetRandom2(sigmaPosXY, sigmaPosZ);
	

	//	cout << "Pass Inside 1.6 " << endl;
	sigmaPosZ *= 1.e4;
	sigmaPosXY *= 1.e4;

	//	cout << "Pass Inside 2 " << endl;


	/*if (h1DcaZ1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->ComputeIntegral())
	  {
	  do sigmaPosZ = h1DcaZ1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom() * 1e4;
	  while (fabs(sigmaPosZ) > 1.e3);
	  }

	  if (h1DcaXY1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->ComputeIntegral())
	  {
	  do sigmaPosXY = h1DcaXY1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom() * 1e4;
	  while (fabs(sigmaPosXY) > 1.e3);
	  }
	  */

	TVector3 newPos(pos);
	newPos.SetZ(0); 
	TVector3 momPerp(-rMom.Vect().Y(), rMom.Vect().X(), 0.0);
	newPos -= momPerp.Unit() * sigmaPosXY;

	DCATrackCut = DCAR[iParticleIndex][iPtIndex] * NSigmaDCA;

	//cout << "Pass Inside 3 " << endl;

	//ActualSmearSigma = sqrt(sigmaPosXY * sigmaPosXY + sigmaPosZ * sigmaPosZ);

	return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

int getEtaIndexDca(double Eta)
{
	for (int i = 0; i < nEtasDca; i++)
	{
		if ((Eta >= EtaEdgeDca[i]) && (Eta < EtaEdgeDca[i + 1]))
			return i;
	}
	return nEtasDca - 1 ;
}

int getPtIndexDca(double pT)
{
	for (int i = 0; i < nPtBinsDca; i++)
	{
		if ((pT >= ptEdgeDca[i]) && (pT < ptEdgeDca[i + 1]))
			return i;
	}
	return nPtBinsDca - 1 ;
}

int getVzIndexDca(double Vz)
{
	for (int i = 0; i < nVzsDca; i++)
	{
		if ((Vz >= VzEdgeDca[i]) && (Vz < VzEdgeDca[i + 1]))
			return i;
	}
	return nVzsDca - 1 ;
}



void init(){
	gSystem->ListLibraries();

	gStyle->SetOptStat(0);
	gRandom = new TRandom3(0);

	Multi = 0;
	
	PiPFromDs = -1;

	TrackEff = 1;

	EffWeightCounts = 0;	

	if(doCut == 0){
		PiKKMassWindow = 100;
		KKMassWindow = 100;
		ThetaCut = -1;
		DCACut = 999999;
		//	DCATrackCut = -1;
	//	DCATrackCut = 200;
		DCATrackCut = 0;
		minPtCut=0.0;
		TrackEtaCut = 9999;
		
	
		etaCut = 10000;
		BetaCut = 10000;
		DCACut = 1000;
		BsMassDown = -1;
		BsMassUp = 999;
		NSigmaDCA = 0.0;
	//	DCACut = 25;
		
	//	ThetaCut = 0.8;
	//	DCATrackCut = 100;
	//	DCACut = 50;

		
		/*
		//New Cuts
		BsMassDown = -1;
		BsMassUp = 999;
		etaCut = 10000;
		BetaCut = 10000;
		minPtCut = 0.1;
		DCATrackCut = 0;
		DCACut = 99999;
		ThetaCut = -1;
		PiKKMassWindow = 100;
		KKMassWindow = 100;
		*/
		
			
	}

	if(doCut == 1){
		PiKKMassWindow = 100;
		KKMassWindow = 100;
		ThetaCut = -1;
		DCACut = 99999;
		//	DCATrackCut = -1;
	//	DCATrackCut = 200;
		DCATrackCut = 200;
		minPtCut=0.1;
		TrackEtaCut = 1.1;
	}



	if(doCut == 2){
		
		PiKKMassWindow = 0.2;
		KKMassWindow = 0.1;
		ThetaCut = 0.9;
		DCACut = 100;
		DCATrackCut = 200;
		minPtCut=0.1;
		TrackEtaCut = 1.1;
	}

	if(doCut == 3){

		PiKKMassWindow = 0.2;
		KKMassWindow = 0.1;
		ThetaCut = 0.80;
		DCACut = 1000;
		DCATrackCut = 100;
		minPtCut=0.1;
		TrackEtaCut = 1.1;
	}


	if(doCut == 4){

		PiKKMassWindow = 0.2;
		KKMassWindow = 0.1;
		ThetaCut = 0.80;
		DCACut = 25;
		//DCATrackCut = 50;
		NSigmaDCA = 2.0;
		minPtCut=0.1;
		TrackEtaCut = 1.1;
	}

	if(UseSTAR){

		fPionMomResolution = new TF1("fPionMomResolution","0.00332099 -0.000868694/x -2.04427e-05*x*x + 0.0015437 * x +  0.000544441/(x*x)",TrackPTMin,TrackPTMax);
		fKaonMomResolution = new TF1("fKaonMomResolution","0.00271293 -0.000125306/x  -3.26199e-05*x*x + 0.00171111 * x +  0.000674163/(x*x)",TrackPTMin,TrackPTMax);
		fProtonMomResolution = new TF1("fProtonMomResolution","0.00346318 -0.00152079/x -2.8538e-05*x*x + 0.00161474 * x +  0.00203283/(x*x)",TrackPTMin,TrackPTMax);
			
	    fTpcPi = new TF1("fTpcPi","0.940158*exp(-pow(x/0.229646,-2.02656))",0.2,20.);
	    fTpcK = new TF1("fTpcK","0.947194*exp(-pow(x/0.259668,-1.46521))",0.2,20.);
		fTpcP = new TF1("fTpcP","0.924992*exp(-pow(x/0.562471,-3.37468))",0.2,20.);

	}
	cout << "Pass 0" << endl;
	cout << "DONE LOADING LIBRARIES!!!" << endl;


	nt_sig->Branch("Event",&EventSignal,"Event/I");
	nt_sig->Branch("pthat",&pthat,"pthat/F");

	nt_sig->Branch("BsPt",&BsPt,"BsPt/F");
	nt_sig->Branch("BsMass",&BsMass,"BsMass/F");
	nt_sig->Branch("BsY",&BsY,"BsY/F");
	nt_sig->Branch("BsP",&BsP,"BsP/F");
	nt_sig->Branch("BsLifeTime",&BsLifeTime,"BsLifeTime/F");
	nt_sig->Branch("KKMass",&KKMass,"KKMass/F");
	nt_sig->Branch("PiKKMass",&PiKKMass,"PiKKMass/F");

	nt_sig->Branch("Eff_PIP",&Eff_PIP,"Eff_PIP/F");
	nt_sig->Branch("Eff_PIM",&Eff_PIM,"Eff_PIM/F");
	nt_sig->Branch("Eff_KP",&Eff_KP,"Eff_KP/F");
	nt_sig->Branch("Eff_KM",&Eff_KM,"Eff_KM/F");
	nt_sig->Branch("BCandEffWg",&BCandEffWg,"BCandEffWg/F");

	nt_sig->Branch("Bvtxx",&Bvtxx,"Bvtxx/F");
	nt_sig->Branch("Bvtxy",&Bvtxy,"Bvtxy/F");
	nt_sig->Branch("Bvtxz",&Bvtxz,"Bvtxz/F");

//From Bs//
/*
	nt_sig->Branch("FromBsPiP",&FromBsPiP,"FromBsPiP/I");
	nt_sig->Branch("FromBsPiM",&FromBsPiM,"FromBsPiM/I");
	nt_sig->Branch("FromBsKP",&FromBsKP,"FromBsKP/I");
	nt_sig->Branch("FromBsKM",&FromBsKM,"FromBsKM/I");
	nt_sig->Branch("FromBsKM",&FromBsKM,"FromBsKM/I");
*/



	nt_sig->Branch("IsSigPiP",&IsSigPiP,"IsSigPiP/I");
	nt_sig->Branch("IsSigPiM",&IsSigPiM,"IsSigPiM/I");
	nt_sig->Branch("IsSigKP",&IsSigKP,"IsSigKP/I");
	nt_sig->Branch("IsSigKM",&IsSigKM,"IsSigKM/I");

	cout << "Pass 2" << endl;


	nt_sig->Branch("kpPt",&kpPt,"kpPt/F");
	nt_sig->Branch("kmPt",&kmPt,"kmPt/F");
	nt_sig->Branch("pipPt",&pipPt,"pipPt/F");
	nt_sig->Branch("pimPt",&pimPt,"pimPt/F");


	nt_sig->Branch("kpEta",&kpEta,"kpEta/F");
	nt_sig->Branch("kmEta",&kmEta,"kmEta/F");
	nt_sig->Branch("pipEta",&pipEta,"pipEta/F");
	nt_sig->Branch("pimEta",&pimEta,"pimEta/F");


	nt_sig->Branch("kpDca",&kpDca,"kpDca/F");
	nt_sig->Branch("kmDca",&kmDca,"kmDca/F");
	nt_sig->Branch("pipDca",&pipDca,"pipDca/F");
	nt_sig->Branch("pimDca",&pimDca,"pimDca/F");


	nt_sig->Branch("dcadaughter",&dcadaughter,"dcadaughter/F");
	nt_sig->Branch("dcaToPv_Bs",&dcaToPv_Bs,"dcaToPv_Bs/F");
	nt_sig->Branch("decayLength_Bs",&decayLength_Bs,"decayLength_Bs/F");
	nt_sig->Branch("cosTheta_Bs",&cosTheta_Bs,"cosTheta_Bs/F");
	nt_sig->Branch("weight",&weight,"weight/F");
	nt_sig->Branch("primaryPt",&primaryPt,"primaryPt/F");
	nt_sig->Branch("ptWg",&ptWg,"ptWg/F");


	nt_sig->Branch("dcapippim",&dcapippim,"dcapippim/F");
	nt_sig->Branch("dcakppip",&dcakppip,"dcakppip/F");
	nt_sig->Branch("dcakppim",&dcakppim,"dcakppim/F");
	nt_sig->Branch("dcakmpip",&dcakmpip,"dcakmpip/F");
	nt_sig->Branch("dcakmpim",&dcakmpim,"dcakmpim/F");
	nt_sig->Branch("dcakpkm",&dcakpkm,"dcakpkm/F");

	cout << "Pass 3" << endl;



	//Ds Vertex//
	nt_sig->Branch("dcadaughter_Ds",&dcadaughter_Ds,"dcadaughter_Ds/F");
	nt_sig->Branch("dcaToPv_Ds",&dcaToPv_Ds,"dcaToPv_Ds/F");
	nt_sig->Branch("decayLength_Ds",&decayLength_Ds,"decayLength_Ds/F");
	nt_sig->Branch("cosTheta_Ds",&cosTheta_Ds,"cosTheta_Ds/F");


	nt_sig->Branch("Dvtxx",&Dvtxx,"Dvtxx/F");
	nt_sig->Branch("Dvtxy",&Dvtxy,"Dvtxy/F");
	nt_sig->Branch("Dvtxz",&Dvtxz,"Dvtxz/F");

	
	
	nt_sig->Branch("PIDPIP",&PIDPIP,"PIDPIP/I");
	nt_sig->Branch("PIDPIM",&PIDPIM,"PIDPIM/I");
	nt_sig->Branch("PIDKP",&PIDKP,"PIDKP/I");
	nt_sig->Branch("PIDKM",&PIDKM,"PIDKM/I");

	//Mother PID Info
	/*



	nt_sig->Branch("MomPIDPIP",&MomPIDPIP,"MomPIDPIP/I");
	nt_sig->Branch("MomPIDPIM",&MomPIDPIM,"MomPIDPIM/I");
	nt_sig->Branch("MomPIDKP",&MomPIDKP,"MomPIDKP/I");
	nt_sig->Branch("MomPIDKM",&MomPIDKM,"MomPIDKM/I");

	nt_sig->Branch("MotherPIDPIP",&MotherPIDPIP,"MotherPIDPIP/I");
	nt_sig->Branch("MotherPIDPIM",&MotherPIDPIM,"MotherPIDPIM/I");
	nt_sig->Branch("MotherPIDKP",&MotherPIDKP,"MotherPIDKP/I");
	nt_sig->Branch("MotherPIDKM",&MotherPIDKM,"MotherPIDKM/I");
	*/


	cout << "DONE BRANCHING NT SIG" << endl;

	/*
	   nt_Gen->Branch("PiPMotherList",&PiPMotherList,"PiPMotherList/I");
	   nt_Gen->Branch("PiMMotherList",&PiMMotherList,"PiMMotherList/I");
	   nt_Gen->Branch("KPMotherList",&KPMotherList,"KPMotherList/I");
	   nt_Gen->Branch("KMMotherList",&KMMotherList,"KMMotherList/I");
	   */
	nt_Gen->Branch("EventGen",&EventGen,"EventGen/I");
	nt_Gen->Branch("pthatGen",&pthatGen,"pthatGen/F");
	nt_Gen->Branch("BsPtGen",&BsPtGen,"BsPtGen/F");
	nt_Gen->Branch("BsMassGen",&BsMassGen,"BsMassGen/F");
	nt_Gen->Branch("BsYGen",&BsYGen,"BsYGen/F");
	nt_Gen->Branch("BsPGen",&BsPGen,"BsPGen/F");
	nt_Gen->Branch("BsLifeTimeGen",&BsLifeTimeGen,"BsLifeTimeGen/F");

	nt_Gen->Branch("KKMassGen",&KKMassGen,"KKMassGen/F");
	nt_Gen->Branch("PiKKMassGen",&PiKKMassGen,"PiKKMassGen/F");


	nt_Gen->Branch("BvtxxGen",&BvtxxGen,"BvtxxGen/F");
	nt_Gen->Branch("BvtxyGen",&BvtxyGen,"BvtxyGen/F");
	nt_Gen->Branch("BvtxzGen",&BvtxzGen,"BvtxzGen/F");



	nt_Gen->Branch("FromBsPiPGen",&FromBsPiPGen,"FromBsPiPGen/I");
	nt_Gen->Branch("FromBsPiMGen",&FromBsPiMGen,"FromBsPiMGen/I");
	nt_Gen->Branch("FromBsKPGen",&FromBsKPGen,"FromBsKPGen/I");
	nt_Gen->Branch("FromBsKMGen",&FromBsKMGen,"FromBsKMGen/I");
	nt_Gen->Branch("FromBsKMGen",&FromBsKMGen,"FromBsKMGen/I");



	nt_Gen->Branch("IsSigPiPGen",&IsSigPiPGen,"IsSigPiPGen/I");
	nt_Gen->Branch("IsSigPiMGen",&IsSigPiMGen,"IsSigPiMGen/I");
	nt_Gen->Branch("IsSigKPGen",&IsSigKPGen,"IsSigKPGen/I");
	nt_Gen->Branch("IsSigKMGen",&IsSigKMGen,"IsSigKMGen/I");



	nt_Gen->Branch("kpPtGen",&kpPtGen,"kpPtGen/F");
	nt_Gen->Branch("kmPtGen",&kmPtGen,"kmPtGen/F");
	nt_Gen->Branch("pipPtGen",&pipPtGen,"pipPtGen/F");
	nt_Gen->Branch("pimPtGen",&pimPtGen,"pimPtGen/F");


	nt_Gen->Branch("kpEtaGen",&kpEtaGen,"kpEtaGen/F");
	nt_Gen->Branch("kmEtaGen",&kmEtaGen,"kmEtaGen/F");
	nt_Gen->Branch("pipEtaGen",&pipEtaGen,"pipEtaGen/F");
	nt_Gen->Branch("pimEtaGen",&pimEtaGen,"pimEtaGen/F");


	nt_Gen->Branch("kpDcaGen",&kpDcaGen,"kpDcaGen/F");
	nt_Gen->Branch("kmDcaGen",&kmDcaGen,"kmDcaGen/F");
	nt_Gen->Branch("pipDcaGen",&pipDcaGen,"pipDcaGen/F");
	nt_Gen->Branch("pimDcaGen",&pimDcaGen,"pimDcaGen/F");


	nt_Gen->Branch("dcadaughterGen",&dcadaughterGen,"dcadaughterGen/F");
	nt_Gen->Branch("dcaToPv_BsGen",&dcaToPv_BsGen,"dcaToPv_BsGen/F");
	nt_Gen->Branch("decayLength_BsGen",&decayLength_BsGen,"decayLength_BsGen/F");
	nt_Gen->Branch("cosTheta_BsGen",&cosTheta_BsGen,"cosTheta_BsGen/F");
	nt_Gen->Branch("weightGen",&weightGen,"weightGen/F");
	nt_Gen->Branch("primaryPtGen",&primaryPtGen,"primaryPtGen/F");
	nt_Gen->Branch("ptWgGen",&ptWgGen,"ptWgGen/F");


	nt_Gen->Branch("dcapippimGen",&dcapippimGen,"dcapippimGen/F");
	nt_Gen->Branch("dcakppipGen",&dcakppipGen,"dcakppipGen/F");
	nt_Gen->Branch("dcakppimGen",&dcakppimGen,"dcakppimGen/F");
	nt_Gen->Branch("dcakmpipGen",&dcakmpipGen,"dcakmpipGen/F");
	nt_Gen->Branch("dcakmpimGen",&dcakmpimGen,"dcakmpimGen/F");
	nt_Gen->Branch("dcakpkmGen",&dcakpkmGen,"dcakpkmGen/F");


	//Ds Vertex//
	nt_Gen->Branch("dcadaughter_DsGen",&dcadaughter_DsGen,"dcadaughter_DsGen/F");
	nt_Gen->Branch("dcaToPv_DsGen",&dcaToPv_DsGen,"dcaToPv_DsGen/F");
	nt_Gen->Branch("decayLength_DsGen",&decayLength_DsGen,"decayLength_DsGen/F");
	nt_Gen->Branch("cosTheta_DsGen",&cosTheta_DsGen,"cosTheta_DsGen/F");


	nt_Gen->Branch("DvtxxGen",&DvtxxGen,"DvtxxGen/F");
	nt_Gen->Branch("DvtxyGen",&DvtxyGen,"DvtxyGen/F");
	nt_Gen->Branch("DvtxzGen",&DvtxzGen,"DvtxzGen/F");

	nt_Gen->Branch("PIDPIPGen",&PIDPIPGen,"PIDPIPGen/I");
	nt_Gen->Branch("PIDPIMGen",&PIDPIMGen,"PIDPIMGen/I");
	nt_Gen->Branch("PIDKPGen",&PIDKPGen,"PIDKPGen/I");
	nt_Gen->Branch("PIDKMGen",&PIDKMGen,"PIDKMGen/I");


	nt_Gen->Branch("MomPIDPIPGen",&MomPIDPIPGen,"MomPIDPIPGen/I");
	nt_Gen->Branch("MomPIDPIMGen",&MomPIDPIMGen,"MomPIDPIMGen/I");
	nt_Gen->Branch("MomPIDKPGen",&MomPIDKPGen,"MomPIDKPGen/I");
	nt_Gen->Branch("MomPIDKMGen",&MomPIDKMGen,"MomPIDKMGen/I");

	nt_Gen->Branch("MotherPIDPIPGen",&MotherPIDPIPGen,"MotherPIDPIPGen/I");
	nt_Gen->Branch("MotherPIDPIMGen",&MotherPIDPIMGen,"MotherPIDPIMGen/I");
	nt_Gen->Branch("MotherPIDKPGen",&MotherPIDKPGen,"MotherPIDKPGen/I");
	nt_Gen->Branch("MotherPIDKMGen",&MotherPIDKMGen,"MotherPIDKMGen/I");


	cout << "DONE BRANCHING NT SIG" << endl;

	nt_evt->Branch("Evt",&Evt,"Evt/I");
	nt_evt->Branch("Multi",&Multi,"Multi/I");

	cout << "DONE BRANCHING NT TRACK" << endl;

	/*
	//Small Sample Fast//

	nt_sig->Branch("IsSigPiP",&IsSigPiP,"IsSigPiP/I");
	nt_sig->Branch("IsSigPiM",&IsSigPiM,"IsSigPiM/I");
	nt_sig->Branch("IsSigKP",&IsSigKP,"IsSigKP/I");
	nt_sig->Branch("IsSigKM",&IsSigKM,"IsSigKM/I");



	nt_sig->Branch("MomPIDPIP",&MomPIDPIP,"MomPIDPIP/I");
	nt_sig->Branch("MomPIDPIM",&MomPIDPIM,"MomPIDPIM/I");
	nt_sig->Branch("MomPIDKP",&MomPIDKP,"MomPIDKP/I");
	nt_sig->Branch("MomPIDKM",&MomPIDKM,"MomPIDKM/I");

	nt_sig->Branch("MotherPIDPIP",&MotherPIDPIP,"MotherPIDPIP/I");
	nt_sig->Branch("MotherPIDPIM",&MotherPIDPIM,"MotherPIDPIM/I");
	nt_sig->Branch("MotherPIDKP",&MotherPIDKP,"MotherPIDKP/I");
	nt_sig->Branch("MotherPIDKM",&MotherPIDKM,"MotherPIDKM/I");
	nt_sig->Branch("decayLength_Bs",&decayLength_Bs,"decayLength_Bs/F");


	nt_Gen->Branch("IsSigPiPGen",&IsSigPiPGen,"IsSigPiPGen/I");
	nt_Gen->Branch("IsSigPiMGen",&IsSigPiMGen,"IsSigPiMGen/I");
	nt_Gen->Branch("IsSigKPGen",&IsSigKPGen,"IsSigKPGen/I");
	nt_Gen->Branch("IsSigKMGen",&IsSigKMGen,"IsSigKMGen/I");

	nt_Gen->Branch("PIDPIPGen",&PIDPIPGen,"PIDPIPGen/I");
	nt_Gen->Branch("PIDPIMGen",&PIDPIMGen,"PIDPIMGen/I");
	nt_Gen->Branch("PIDKPGen",&PIDKPGen,"PIDKPGen/I");
	nt_Gen->Branch("PIDKMGen",&PIDKMGen,"PIDKMGen/I");


	nt_Gen->Branch("MomPIDPIPGen",&MomPIDPIPGen,"MomPIDPIPGen/I");
	nt_Gen->Branch("MomPIDPIMGen",&MomPIDPIMGen,"MomPIDPIMGen/I");
	nt_Gen->Branch("MomPIDKPGen",&MomPIDKPGen,"MomPIDKPGen/I");
	nt_Gen->Branch("MomPIDKMGen",&MomPIDKMGen,"MomPIDKMGen/I");

	nt_Gen->Branch("MotherPIDPIPGen",&MotherPIDPIPGen,"MotherPIDPIPGen/I");
	nt_Gen->Branch("MotherPIDPIMGen",&MotherPIDPIMGen,"MotherPIDPIMGen/I");
	nt_Gen->Branch("MotherPIDKPGen",&MotherPIDKPGen,"MotherPIDKPGen/I");
	nt_Gen->Branch("MotherPIDKMGen",&MotherPIDKMGen,"MotherPIDKMGen/I");
	nt_Gen->Branch("decayLength_BsGen",&decayLength_BsGen,"decayLength_BsGen/F");
	*/

	//DONE FAST//

	/*

	   SignalFile = new TFile("InputROOT/G4sPHENIX.root_g4svtx_eval.root");
	   SignalFile->cd();

	   ntp_track = (TTree *) SignalFile->Get("ntp_track");




	   TotalCand = ntp_track->GetEntries();

	   ntp_track->SetBranchAddress("gpx",&gpx);
	   ntp_track->SetBranchAddress("gpy",&gpy);
	   ntp_track->SetBranchAddress("gpz",&gpz);
	   ntp_track->SetBranchAddress("gflavor",&gflavor);
	//ntp_track->SetBranchAddress("event",&EventSig);


	ntp_track->SetBranchAddress("gvx",&gvx);
	ntp_track->SetBranchAddress("gvy",&gvy);
	ntp_track->SetBranchAddress("gvz",&gvz);



	PassPre = 0;




	nt_gen->Branch("gBvtxx",&gBvtxx,"gBvtxx/F");
	nt_gen->Branch("gBvtxy",&gBvtxy,"gBvtxy/F");
	nt_gen->Branch("gBvtxz",&gBvtxz,"gBvtxz/F");

	nt_gen->Branch("gBpx",&gBpx,"gBpx/F");
	nt_gen->Branch("gBpy",&gBpy,"gBpy/F");
	nt_gen->Branch("gBpz",&gBpz,"gBpz/F");
	nt_gen->Branch("gBpt",&gBpt,"gBpt/F");

	ntp_gtrack = (TTree *) SignalFile->Get("ntp_gtrack");

	//ntp_gtrack->SetBranchAddress("event",&gInEvent);
	ntp_gtrack->SetBranchAddress("gvx",&gInBvtxx);
	ntp_gtrack->SetBranchAddress("gvy",&gInBvtxy);
	ntp_gtrack->SetBranchAddress("gvz",&gInBvtxz);

	ntp_gtrack->SetBranchAddress("gpx",&gInBpx);
	ntp_gtrack->SetBranchAddress("gpy",&gInBpy);
	ntp_gtrack->SetBranchAddress("gpz",&gInBpz);
	ntp_gtrack->SetBranchAddress("gpt",&gInBpt);
	*/
	TString DCA2DFile = "InputROOT/2DProjection_DcaXyZ_sPHENIX_40_80.root";

	fDca2D = new TFile(DCA2DFile.Data());
	fDca2D->cd();

	for (int iParticle = 0; iParticle < nParticles; ++iParticle)
	{
		for (int iPt = 0; iPt < nPtBinsDca; ++iPt)
		{
			h2Dca[iParticle][iPt] = (TH2F*) fDca2D->Get(Form("mh2DcaPtPart_%i_%i", iParticle, iPt)) ;
			//cout << "SMEAR Track RMS = " << h2Dca[iParticle][iPt]->GetRMS() << endl;
			h2DcaXY[iParticle][iPt] = (TH1F *) h2Dca[iParticle][iPt]->ProjectionX(); 
			h2DcaZ[iParticle][iPt] = (TH1F *) h2Dca[iParticle][iPt]->ProjectionY(); 
			
			DCAXY[iParticle][iPt] = h2DcaXY[iParticle][iPt]->GetRMS() * 1000;
			DCAZ[iParticle][iPt] = h2DcaZ[iParticle][iPt]->GetRMS() * 1000;
			DCAR[iParticle][iPt] = sqrt(DCAXY[iParticle][iPt] * DCAXY[iParticle][iPt]  + DCAZ[iParticle][iPt] * DCAZ[iParticle][iPt]);
			//cout << "DCAR[iParticle][iPt] = " << DCAR[iParticle][iPt] << endl;
			//h2Dca[iParticle][iPt]->SetDirectory(0);
		}
	}

	cout << "DONE Loading DCA Bro" << endl;

	TString FONLLFile = "InputROOT/BsFONLL2DMap.root"; 

	FONLLInput = new TFile(FONLLFile.Data());
	FONLLInput->cd();


	FONLL2DMap = (TH2D * ) FONLLInput->Get("BsFONLL2DMap");





	TString SpectraFiles = "InputROOT/input_DaughterPtWg.root"; 

	InputSpectra = new TFile(SpectraFiles.Data());

	hpiPtWg = (TH1F *) InputSpectra->Get("hpiPtWg");
	hkPtWg = (TH1F *) InputSpectra->Get("hkPtWg");
	hpPtWg = (TH1F *) InputSpectra->Get("hpPtWg");


	if(DrawInput == 1){

		TCanvas * c = new TCanvas("c","c",600,600);
		c->cd();
		
		FONLL2DMap->Draw("COLZ");

		c->SaveAs("Plots/FONLL2DMaps.png");

		c->SetLogy();



		hpiPtWg->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		hpiPtWg->GetYaxis()->SetTitle("Counts");

		hpiPtWg->SetMarkerStyle(20);
		hkPtWg->SetMarkerStyle(21);
		hpPtWg->SetMarkerStyle(22);


		hpiPtWg->SetMarkerSize(0.5);
		hkPtWg->SetMarkerSize(0.5);
		hpPtWg->SetMarkerSize(0.5);

		hpiPtWg->SetLineColor(kBlack);
		hkPtWg->SetLineColor(kRed);
		hpPtWg->SetLineColor(kBlue);




		hpiPtWg->SetMarkerColor(kBlack);
		hkPtWg->SetMarkerColor(kRed);
		hpPtWg->SetMarkerColor(kBlue);

		hpiPtWg->Draw("p");
		hkPtWg->Draw("pSAME");
		hpPtWg->Draw("pSAME");



		TLegend *leg = new TLegend(0.36,0.60,0.75,0.85,NULL,"brNDC");
		leg->SetBorderSize(0);
		leg->SetTextSize(0.040);
		leg->SetTextFont(42);
		leg->SetFillStyle(0);
		leg->SetLineWidth(3);

		leg->AddEntry(hpiPtWg,"Pion p_{T} Spectrum","l");
		leg->AddEntry(hkPtWg,"Kaon p_{T} Spectrum","l");
		leg->AddEntry(hpPtWg,"Proton p_{T} Spectrum","l");

		leg->Draw("SAME");



		//c->SetLogy();
		c->SaveAs("Plots/InputSpectra.png");


		TCanvas * c2 = new TCanvas("c2","c2",600,600);
		c2->cd();

		c2->SetLogy();

		fPionMomResolution->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		fPionMomResolution->GetYaxis()->SetTitle("#sigma p_{T}/p_{T}");
		fPionMomResolution->SetLineColor(kBlack);
		fPionMomResolution->SetLineWidth(2);

		fKaonMomResolution->SetLineColor(kRed);
		fKaonMomResolution->SetLineWidth(2);

		fProtonMomResolution->SetLineColor(kBlue);
		fProtonMomResolution->SetLineWidth(2);

		fPionMomResolution->GetYaxis()->SetTitleOffset(1.0);

		fPionMomResolution->Draw("R");
		fKaonMomResolution->Draw("SAME");
		fProtonMomResolution->Draw("SAME");

		


		TLegend *leg2 = new TLegend(0.36,0.60,0.75,0.85,NULL,"brNDC");
		leg2->SetBorderSize(0);
		leg2->SetTextSize(0.040);
		leg2->SetTextFont(42);
		leg2->SetFillStyle(0);
		leg2->SetLineWidth(3);

		leg2->AddEntry(fPionMomResolution,"Pion p_{T} Resolution","l");
		leg2->AddEntry(fKaonMomResolution,"Kaon p_{T} Resolution","l");
		leg2->AddEntry(fProtonMomResolution,"Proton p_{T} Resolution","l");

		leg2->Draw("SAME");


		c2->SaveAs("Plots/InputPtReso.png");


		//Tracking Efficiency//

		TCanvas * c3 = new TCanvas("c3","c3",600,600);
		c3->cd();

		c3->SetLogy();

		fTpcPi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		fTpcPi->GetYaxis()->SetTitle("Tracking Efficiency");
		fTpcPi->SetLineColor(kBlack);
		fTpcPi->SetLineWidth(2);

		fTpcK->SetLineColor(kRed);
		fTpcK->SetLineWidth(2);

		fTpcP->SetLineColor(kBlue);
		fTpcP->SetLineWidth(2);

		fPionMomResolution->GetYaxis()->SetTitleOffset(1.0);

		fTpcPi->Draw("R");
		fTpcK->Draw("SAME");
		fTpcP->Draw("SAME");

	


		TLegend *leg4 = new TLegend(0.36,0.60,0.75,0.85,NULL,"brNDC");
		leg4->SetBorderSize(0);
		leg4->SetTextSize(0.040);
		leg4->SetTextFont(42);
		leg4->SetFillStyle(0);
		leg4->SetLineWidth(3);

		leg4->AddEntry(fTpcPi,"Pion Tracking Efficiency","l");
		leg4->AddEntry(fTpcK,"Kaon Tracking Efficiency","l");
		leg4->AddEntry(fTpcP,"Proton Tracking Efficiency","l");

		leg4->Draw("SAME");


		c3->SaveAs("Plots/TrackingEff.png");

		//DOne Eff//




		TCanvas * c1 = new TCanvas("c1","c1",600,600);
		c1->cd();

		for(int i = 0; i < 	nParticles; i++){

			for(int j = 0; j < nPtBinsDca; j ++){

				h2Dca[i][j]->GetXaxis()->SetTitle("DCA XY (cm)");
				h2Dca[i][j]->GetXaxis()->SetTitle("DCA Z (cm)");

				h2Dca[i][j]->Draw("COLZ");

				c1->SaveAs(Form("Plots/DCAInput/DCADis_%d_%d.png",i,j));

			}


		}

	}

	cout << "DONE INIT" << endl;

}



void clean(){

	//Special debug variables

	/*
	   vecMomList_pip.clear();
	   vecMomList_pim.clear();
	   vecMomList_kp.clear();
	   vecMomList_km.clear();
	   */


	//DONE

	vecPt_km.clear(); // k minus
	vecPID_km.clear();
	vecMomPID_km.clear();
	vecMotherPID_km.clear();	
	vecEta_km.clear();
	vecPhi_km.clear();
	vecPos_km.clear();
	vecType_km.clear();
	vecClean_km.clear();
	vecHybrid_km.clear();
	vecCleanIdeal_km.clear();
	vecHybridIdeal_km.clear();
	vecDCA_km.clear();
	vecBottom_km.clear();
	vecTrackID_km.clear();
	vecEff_km.clear();


	vecPt_kp.clear(); // k plus
	vecPID_kp.clear();
	vecMomPID_kp.clear();
	vecMotherPID_kp.clear();	
	vecEta_kp.clear();
	vecPhi_kp.clear();
	vecPos_kp.clear();
	vecType_kp.clear();
	vecClean_kp.clear();
	vecHybrid_kp.clear();
	vecCleanIdeal_kp.clear();
	vecHybridIdeal_kp.clear();
	vecDCA_kp.clear();
	vecBottom_kp.clear();
	vecTrackID_kp.clear();
	vecEff_kp.clear();


	vecPt_pip.clear(); // pi plus
	vecPID_pip.clear();
	vecMomPID_pip.clear();
	vecMotherPID_pip.clear();		
	vecEta_pip.clear();
	vecPhi_pip.clear();
	vecPos_pip.clear();
	vecType_pip.clear();
	vecClean_pip.clear();
	vecHybrid_pip.clear();
	vecCleanIdeal_pip.clear();
	vecHybridIdeal_pip.clear();
	vecDCA_pip.clear();
	vecBottom_pip.clear();
	vecTrackID_pip.clear();
	vecEff_pip.clear();


	vecPt_pim.clear(); // pi minus
	vecPID_pim.clear();
	vecMomPID_pim.clear();
	vecMotherPID_pim.clear();		
	vecEta_pim.clear();
	vecPhi_pim.clear();
	vecPos_pim.clear();
	vecType_pim.clear();
	vecClean_pim.clear();
	vecHybrid_pim.clear();
	vecCleanIdeal_pim.clear();
	vecHybridIdeal_pim.clear();
	vecDCA_pim.clear();
	vecBottom_pim.clear();
	vecTrackID_pim.clear();
	vecEff_pim.clear();

		



	//Gen Vector//



	vecPt_kmGen.clear(); // k minus
	vecPID_kmGen.clear();	
	vecMomPID_kmGen.clear();
	vecMotherPID_kmGen.clear();		
	vecEta_kmGen.clear();
	vecPhi_kmGen.clear();
	vecPos_kmGen.clear();
	vecType_kmGen.clear();
	vecClean_kmGen.clear();
	vecHybrid_kmGen.clear();
	vecCleanIdeal_kmGen.clear();
	vecHybridIdeal_kmGen.clear();
	vecDCA_kmGen.clear();
	vecBottom_kmGen.clear();
	vecTrackID_kmGen.clear();


	vecPt_kpGen.clear(); // k plus
	vecPID_kpGen.clear();	
	vecMomPID_kpGen.clear();
	vecMotherPID_kpGen.clear();			
	vecEta_kpGen.clear();
	vecPhi_kpGen.clear();
	vecPos_kpGen.clear();
	vecType_kpGen.clear();
	vecClean_kpGen.clear();
	vecHybrid_kpGen.clear();
	vecCleanIdeal_kpGen.clear();
	vecHybridIdeal_kpGen.clear();
	vecDCA_kpGen.clear();
	vecBottom_kpGen.clear();
	vecTrackID_kpGen.clear();


	vecPt_pipGen.clear(); // pi plus
	vecPID_pipGen.clear();		
	vecMomPID_pipGen.clear();
	vecMotherPID_pipGen.clear();			
	vecEta_pipGen.clear();
	vecPhi_pipGen.clear();
	vecPos_pipGen.clear();
	vecType_pipGen.clear();
	vecClean_pipGen.clear();
	vecHybrid_pipGen.clear();
	vecCleanIdeal_pipGen.clear();
	vecHybridIdeal_pipGen.clear();
	vecDCA_pipGen.clear();
	vecBottom_pipGen.clear();
	vecTrackID_pipGen.clear();


	vecPt_pimGen.clear(); // pi minus
	vecPID_pimGen.clear();			
	vecMomPID_pimGen.clear();
	vecMotherPID_pimGen.clear();	
	vecEta_pimGen.clear();
	vecPhi_pimGen.clear();
	vecPos_pimGen.clear();
	vecType_pimGen.clear();
	vecClean_pimGen.clear();
	vecHybrid_pimGen.clear();
	vecCleanIdeal_pimGen.clear();
	vecHybridIdeal_pimGen.clear();
	vecDCA_pimGen.clear();
	vecBottom_pimGen.clear();
	vecTrackID_pimGen.clear();



}


void end(){

	cout << "-------------------------------------- Summary --------------------------------------" << endl;
	cout << "Total Candidates Input = " << TotalEvents << "   Total Candidates Generated " << nt_Gen->GetEntries()   << "   Total Candidates Left " << nt_sig->GetEntries() << "  Pre Selection Efficiency = " << float(nt_sig->GetEntries())/nt_Gen->GetEntries() *100 << " %" << endl;
	cout << "Efficiency Weighted: Total Candates Input = " << TotalEvents << "   Total Candidates Generated " << nt_Gen->GetEntries()   << "   Total Candidates Left " << EffWeightCounts << "  Tracking + Pre Selections Efficiency = " << float(EffWeightCounts)/nt_Gen->GetEntries() *100 << " %" << endl;	
	cout << "-------------------------------------- DONE SUMMARY --------------------------------------" << endl;

	outFile->Write();
	outFile->Close();
	fDca2D->Close();
	//SignalFile->Close();
	InputSpectra->Close();
	

}
