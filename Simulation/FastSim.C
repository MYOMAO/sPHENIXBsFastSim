#include "FastSim.h"

void FastSim(){


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
		if(RunBackground)	GenerateBackground(Event);
//		if(RunSignal)	GetSignal(Event);
//		if(RunSignal) GetGen(Event);
		if(RunBackground) GenerateDecay(Event);
		if(RunEVTGEN) GenerateEvtGen(Event);
		BackgroundSimulations(Event);
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

	if( vecPt_kp.size()  > 4 && vecPt_km.size() > 4){

		cout << "nPlus = " << nPlus << "   vecPt_kp.size() = " << vecPt_kp.size() << endl;
		cout << "nMinus = " << nMinus << "   vecPt_km.size() = " << vecPt_km.size() << endl;

	}

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

						EventSignal = EventID;
						FromBsKP = vecBottom_kp[i];
						FromBsKM = vecBottom_km[j];
						FromBsPiP = vecBottom_pip[k];
						FromBsPiM = vecBottom_pim[l];

						IsSigKP = vecType_kp[i];
						IsSigKM = vecType_km[j];
						IsSigPiP = vecType_pip[k];
						IsSigPiM = vecType_pim[l];

						TLorentzVector const BsMom = *KP + *KM + *PIP + *PIM;
						TLorentzVector const PhiMom = *KP + *KM;
						TLorentzVector const DsMom = *KP + *KM + *PIM;

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
						PiKKMass = DsMom.M();


						if(fabs(KKMass - PhiMass) > KKMassWindow) continue;
						if(fabs(PiKKMass - DsMass) > PiKKMassWindow) continue;

						ptWg = 1;
						weight = 1;

						if(fabs(BsMom.Rapidity()) > 1.1) continue;
						if(cosTheta_Bs < ThetaCut) continue;

			
						nt_sig->Fill();

						//Fill(*KP,*KM,*PIP,*PIM,  vecPos_kp[i], vecPos_km[j], vecPos_pip[k],vecPos_pim[l], v0_Bs);
					}
					//	nt_sig->Fill();

					PairNow = PairNow + 1;
				}

			}

		}

	}


	double TotalEvents = nPlus *nPlus * nMinus * nMinus;
	double PassEvents = nt_sig->GetEntries() - PassPre;
	double Efficiency = PassEvents/TotalEvents;

	PassPre = nt_sig->GetEntries();
	if(Efficiency > 0) cout << "Total Pairs = " << TotalEvents << "   Passed Pairs = " << PassEvents << "   Efficiency =  " << Efficiency << endl;

	clean();

}



void GenerateBackground( int &EventID){


	gRandom->SetSeed();


	//Initial Histogram




	//	  Double_t nPions_tmp, nKaons_tmp;

	nPiPlus = gRandom->Integer(MaxPlus);
	nKPlus = gRandom->Integer(MaxMinus);





	nPiMinus = nPiPlus;
	nKMinus = nKPlus;

	nPlus = nPiPlus + nKPlus;
	nMinus = nPiMinus + nKMinus;

	cout << "nPlus = " << nPlus << "    nKaons = " << nMinus << endl;




	//Assuming Pi+ = Pi-//

	//Pi+ Loop
	for(int ipi=0; ipi<nPiPlus; ipi++) {
		//	cout << "hpiPtWg->Integral() = " << hpiPtWg->Integral() << endl;	
		float pt = hpiPtWg->GetRandom();//gRandom->Uniform(0.6,20);
		//	cout << "Pass Inside" << endl;
		float eta = gRandom->Uniform(-1,1);
		float phi = gRandom->Uniform(-PI,PI);
		FourMom.SetPtEtaPhiM(pt, eta , phi, M_PION_PLUS);
		rcFourMom = smearMom(FourMom,fPionMomResolution);
		if (rcFourMom.Perp()<minPtCut) continue;
		rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
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

		//Identify as K+//
		vecTrackID_kp.push_back(ipi);
		vecDCA_kp.push_back(rcdca);               
		vecPt_kp.push_back(rcFourMom.Perp());
		vecEta_kp.push_back(rcFourMom.Eta());
		vecPhi_kp.push_back(rcFourMom.Phi());
		vecPos_kp.push_back(rcPos);
		vecType_kp.push_back(0);
		vecBottom_kp.push_back(0); 

	}


	//Pi- Loop
	for(int ipi=0; ipi<nPiMinus; ipi++) {
		//	cout << "hpiPtWg->Integral() = " << hpiPtWg->Integral() << endl;	
		float pt = hpiPtWg->GetRandom();//gRandom->Uniform(0.6,20);
		//	cout << "Pass Inside" << endl;
		float eta = gRandom->Uniform(-1,1);
		float phi = gRandom->Uniform(-PI,PI);
		FourMom.SetPtEtaPhiM(pt, eta , phi, M_PION_MINUS);
		rcFourMom = smearMom(FourMom,fPionMomResolution);
		if (rcFourMom.Perp()<minPtCut) continue;
		rcPos = smearPosData(0, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
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

		//Identify as K-//
		vecTrackID_km.push_back(ipi);
		vecDCA_km.push_back(rcdca);               
		vecPt_km.push_back(rcFourMom.Perp());
		vecEta_km.push_back(rcFourMom.Eta());
		vecPhi_km.push_back(rcFourMom.Phi());
		vecPos_km.push_back(rcPos);
		vecType_km.push_back(0);
		vecBottom_km.push_back(0); 

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
		rcFourMom = smearMom(FourMom,fKaonMomResolution);
		if (rcFourMom.Perp()<minPtCut) continue;
		rcPos = smearPosData(1, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
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


		//Identify as K+//

		vecTrackID_kp.push_back(ik+nPiPlus);
		vecDCA_kp.push_back(rcdca);               
		vecPt_kp.push_back(rcFourMom.Perp());
		vecEta_kp.push_back(rcFourMom.Eta());
		vecPhi_kp.push_back(rcFourMom.Phi());
		vecPos_kp.push_back(rcPos);
		vecType_kp.push_back(0);
		vecBottom_kp.push_back(0); 


	}


	//K- Loop
	for(int ik=0; ik<nKMinus; ik++) {
		float pt = hkPtWg->GetRandom();//gRandom->Uniform(0.6,20);
		float eta = gRandom->Uniform(-1,1);
		float phi = gRandom->Uniform(-PI,PI);
		FourMom.SetPtEtaPhiM(pt, eta , phi, M_KAON_MINUS);
		rcFourMom = smearMom(FourMom,fKaonMomResolution);
		if (rcFourMom.Perp()<minPtCut) continue;
		rcPos = smearPosData(1, 0, 7, rcFourMom, Pos); //1--kaon 0--pi
		float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);
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


		//Identify as K-//

		vecTrackID_km.push_back(ik+nPiMinus);
		vecDCA_km.push_back(rcdca);               
		vecPt_km.push_back(rcFourMom.Perp());
		vecEta_km.push_back(rcFourMom.Eta());
		vecPhi_km.push_back(rcFourMom.Phi());
		vecPos_km.push_back(rcPos);
		vecType_km.push_back(0);
		vecBottom_km.push_back(0); 


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



		//	cout << "particle i = " << i << "   PDG ID = " << TMath::Abs(pythia8.event[i].id()) << endl;


		if (abs(id)!=321 && abs(id)!=211) continue;
		if(id == 211 || id == 321) {

			//cout << "K+/Pi+ Recorded!!" << endl;


			if(id == 211){
				FourMom.SetPtEtaPhiM(pt,eta,phi, M_PION_PLUS);
				rcFourMom = smearMom(FourMom,fPionMomResolution);
				rcPos = smearPosData(0, 0, 7, rcFourMom, Pos);

			}

			if(id == 321){
				FourMom.SetPtEtaPhiM(pt,eta,phi, M_KAON_PLUS);
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



		if(id == -211 || id == -321) {

			//		cout << "K-/Pi- Recorded!!" << endl;


			if(id == -211){

				FourMom.SetPtEtaPhiM(pt,eta,phi, M_PION_MINUS);
				rcFourMom = smearMom(FourMom,fPionMomResolution);
				rcPos = smearPosData(0, 0, 7, rcFourMom, Pos);
				//			cout << "Pass Pos Pi" << endl;

			}

			if(id == -321){

				FourMom.SetPtEtaPhiM(pt,eta,phi, M_KAON_MINUS);
				rcFourMom = smearMom(FourMom,fKaonMomResolution);
				rcPos = smearPosData(1, 0, 7, rcFourMom, Pos);

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

	if(!RunBackground){
		nPlus = 0;
		nMinus = 0;
	}

	TClonesArray daughters("TParticle", 10);
	TLorentzVector* BsVecSig = new TLorentzVector;

	getKinematics(*BsVecSig, BsMassPDG);

	myEvtGenDecayer->Decay(531, BsVecSig);
	myEvtGenDecayer->ImportParticles(&daughters);

	int nTrk = daughters.GetEntriesFast();




	for (int iTrk = 0; iTrk < nTrk; ++iTrk)
	{

		TParticle* ptl0 = (TParticle*)daughters.At(iTrk);
		int id = ptl0->GetPdgCode();
		ptl0->Momentum(FourMom);
		TVector3 rcPosSig(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.);


		if(id == 211 || id == 321) {

			//cout << "K+/Pi+ Recorded!!" << endl;

			//	FourMom = p->momentum();

			if(id == 211){
				//FourMom.SetXYZM(gpx, gpy , gpz, M_PION_MINUS);
				rcFourMom = smearMom(FourMom,fPionMomResolution);
				rcPos = smearPosData(0, 0, 7, rcFourMom, rcPosSig);

			}

			if(id == 321){
				//	FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
				rcFourMom = smearMom(FourMom,fKaonMomResolution);
				rcPos = smearPosData(1, 0, 7, rcFourMom, rcPosSig);
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
			if(id == 211 )	vecType_pip.push_back(1);
			if(id == 321 )	vecType_pip.push_back(2);
			vecBottom_pip.push_back(1); 


			//Identify as K+//
			vecTrackID_kp.push_back(nPlus);
			vecDCA_kp.push_back(rcdca);               
			vecPt_kp.push_back(rcFourMom.Perp());
			vecEta_kp.push_back(rcFourMom.Eta());
			vecPhi_kp.push_back(rcFourMom.Phi());
			vecPos_kp.push_back(rcPos);
			if(id == 211 )	vecType_kp.push_back(2);
			if(id == 321 )	vecType_kp.push_back(1);
			vecBottom_kp.push_back(1); 


			nPlus = nPlus + 1;
		}




		//	cout << "Pass 5" << endl;

		if(id == -211 || id == -321) {

			//	cout << "K-/Pi- Recorded!!" << endl;


			//	FourMom = p->momentum();

			//cout << "Pass 4 Momentum " << endl;

			if(id == -211){
				//	FourMom.SetXYZM(gpx, gpy , gpz, M_PION_MINUS);
				rcFourMom = smearMom(FourMom,fPionMomResolution);

				rcPos = smearPosData(0, 0, 7, rcFourMom, rcPosSig);
				//			cout << "Pass Pos Pi" << endl;

			}

			if(id == -321){
				//	FourMom.SetXYZM(gpx, gpy , gpz, M_KAON_MINUS);
				rcFourMom = smearMom(FourMom,fKaonMomResolution);
				//		cout << "Pass Reso K" << endl;
				rcPos = smearPosData(1, 0, 7, rcFourMom, rcPosSig);
				//		cout << "Pass Pos K" << endl;

			}


			if (rcFourMom.Perp()<minPtCut) continue;

			float rcdca = dca(rcFourMom.Vect(), rcPos, vertex);

			vecTrackID_pim.push_back(nMinus);
			vecDCA_pim.push_back(rcdca);               
			vecPt_pim.push_back(rcFourMom.Perp());
			vecEta_pim.push_back(rcFourMom.Eta());
			vecPhi_pim.push_back(rcFourMom.Phi());
			vecPos_pim.push_back(rcPos);
			if(id == -211 )	vecType_pim.push_back(1);
			if(id == -321 )	vecType_pim.push_back(2);
			vecBottom_pim.push_back(1); 


			//Identify as K-//

			vecTrackID_km.push_back(nMinus);
			vecDCA_km.push_back(rcdca);               
			vecPt_km.push_back(rcFourMom.Perp());
			vecEta_km.push_back(rcFourMom.Eta());
			vecPhi_km.push_back(rcFourMom.Phi());
			vecPos_km.push_back(rcPos);
			if(id == -211 )	vecType_km.push_back(2);
			if(id == -321 )	vecType_km.push_back(1);
			vecBottom_km.push_back(1); 


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
	myEvtGenDecayer->SetDecayTable("InputDECAYFiles/Bs.Mypythia0.DEC");
	cout << "Pass 10" << endl;



}



void getKinematics(TLorentzVector& b, double const mass)
{
	float const pt = gRandom->Uniform(0, 20);
	float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
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

	//cout << "Pass Inside 3 " << endl;

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

	gStyle->SetOptStat(0);
	if(NoCut){
		PiKKMassWindow = 100;
		KKMassWindow = 100;
		ThetaCut = -1;
		DCACut = 99999;
		DCATrackCut = 200;
		minPtCut=0;
	}

	if(NoCut == false){
		PiKKMassWindow = 0.2;
		KKMassWindow = 0.1;
		ThetaCut = 0.0;
		DCACut = 1000;
		DCATrackCut = 0;
		minPtCut=2;
	}

	if(UseSTAR){

		fPionMomResolution = new TF1("fPionMomResolution","0.00332099 -0.000868694/x -2.04427e-05*x*x + 0.0015437 * x +  0.000544441/(x*x)",TrackPTMin,TrackPTMax);
		fKaonMomResolution = new TF1("fKaonMomResolution","0.00271293 -0.000125306/x  -3.26199e-05*x*x + 0.00171111 * x +  0.000674163/(x*x)",TrackPTMin,TrackPTMax);
		fProtonMomResolution = new TF1("fProtonMomResolution","0.00346318 -0.00152079/x -2.8538e-05*x*x + 0.00161474 * x +  0.00203283/(x*x)",TrackPTMin,TrackPTMax);

	}

	cout << "DONE LOADING LIBRARIES!!!" << endl;

	nt_sig->Branch("Event",&EventSignal,"Event/I");
	nt_sig->Branch("pthat",&pthat,"pthat/F");
	nt_sig->Branch("BsPt",&BsPt,"BsPt/F");
	nt_sig->Branch("BsMass",&BsMass,"BsMass/F");
	nt_sig->Branch("BsY",&BsY,"BsY/F");
	nt_sig->Branch("KKMass",&KKMass,"KKMass/F");
	nt_sig->Branch("PiKKMass",&PiKKMass,"PiKKMass/F");


	nt_sig->Branch("Bvtxx",&Bvtxx,"Bvtxx/F");
	nt_sig->Branch("Bvtxy",&Bvtxy,"Bvtxy/F");
	nt_sig->Branch("Bvtxz",&Bvtxz,"Bvtxz/F");



	nt_sig->Branch("FromBsPiP",&FromBsPiP,"FromBsPiP/I");
	nt_sig->Branch("FromBsPiM",&FromBsPiM,"FromBsPiM/I");
	nt_sig->Branch("FromBsKP",&FromBsKP,"FromBsKP/I");
	nt_sig->Branch("FromBsKM",&FromBsKM,"FromBsKM/I");
	nt_sig->Branch("FromBsKM",&FromBsKM,"FromBsKM/I");



	nt_sig->Branch("IsSigPiP",&IsSigPiP,"IsSigPiP/I");
	nt_sig->Branch("IsSigPiM",&IsSigPiM,"IsSigPiM/I");
	nt_sig->Branch("IsSigKP",&IsSigKP,"IsSigKP/I");
	nt_sig->Branch("IsSigKM",&IsSigKM,"IsSigKM/I");



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

	cout << "DONE nt_sig Branching BRO" << endl;

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
			//h2Dca[iParticle][iPt]->SetDirectory(0);
		}
	}

	cout << "DONE Loading DCA Bro" << endl;


	TString SpectraFiles = "InputROOT/input_DaughterPtWg.root"; 

	InputSpectra = new TFile(SpectraFiles.Data());

	hpiPtWg = (TH1F *) InputSpectra->Get("hpiPtWg");
	hkPtWg = (TH1F *) InputSpectra->Get("hkPtWg");
	hpPtWg = (TH1F *) InputSpectra->Get("hpPtWg");

	if(DrawInput == 1){

		TCanvas * c = new TCanvas("c","c",600,600);
		c->cd();
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




	vecPt_km.clear(); // k minus
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


	vecPt_kp.clear(); // k plus
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


	vecPt_pip.clear(); // pi plus
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


	vecPt_pim.clear(); // pi minus
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


}


void end(){

	outFile->Write();
	outFile->Close();
	fDca2D->Close();
	//SignalFile->Close();
	InputSpectra->Close();


}
