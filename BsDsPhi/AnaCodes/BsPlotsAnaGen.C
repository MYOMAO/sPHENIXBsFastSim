#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>


#include "sPhenixStyle.C"
#include "sPhenixStyle.h"

using namespace std;

using std::cout;
using std::endl;


void BsPlotsAnaGen(){

	int UseTrackptCut = 1;

	double	TrackptCutValue =  0.2;

	TString TrackptCut = Form("(pipPtGen > %f && pimPtGen > %f && kpPtGen > %f && kmPtGen > %f)",TrackptCutValue,TrackptCutValue,TrackptCutValue,TrackptCutValue);


	gStyle->SetOptStat(0);
	SetsPhenixStyle();

	gStyle->SetPadRightMargin(1.2);

	//Cuts//

	TString BsSigCut = "(IsSigPiPGen*IsSigPiMGen*IsSigKPGen*IsSigKMGen==1)";
	TString BsSelCut = "(abs(KKMassGen - 1.019461) < 0.1 && abs(PiKKMassGen - 1.96847) < 0.10)";

	TString BsDecayCut = "(IsSigPiPGen>2 && IsSigPiMGen>2 && IsSigKPGen>2 && IsSigKMGen>2)";
	TString BsCombCut = "(IsSigPiPGen==0 && IsSigPiMGen==0 && IsSigKPGen==0 && IsSigKMGen==0)";

	TString BsDecayIDCut = "(IsSigPiPGen>2 && IsSigPiMGen>2 && IsSigKPGen>2 && IsSigKMGen>2) && (PIDPIPGen == 211 && PIDPIMGen == -211 && PIDKPGen == 321 && PIDKMGen == -321)";



	if(UseTrackptCut == 1){

		BsSigCut = Form("(IsSigPiPGen*IsSigPiMGen*IsSigKPGen*IsSigKMGen==1) && %s",TrackptCut.Data());
		BsSelCut = Form("(abs(KKMassGen - 1.019461) < 0.1 && abs(PiKKMassGen - 1.96847) < 0.10) && %s",TrackptCut.Data());

	}

	TString infile = "output.root";

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TTree * BsDsKKPiPi = (TTree *) fin->Get("BsDsKKPiPiGen");

	/*
	const int NVar = 27;


	TString XAxisName[NVar] = {"B_{s}^{0} (K^{+} K^{-} #pi^{+}  #pi^{-}) (GeV/c^{2})","DCA Daughters (#mu m)","SV to PV (Decay) Length (#mu m)","cos(#theta)","|m_{KK} - m_{#phi}| (GeV/c^2)","|m_{#piKK} - m_{D_{s}}| (GeV/c^2)","Postive Track p_{T} (GeV/c)","Negative Track p_{T} (GeV/c)","Postive Track #eta","Negative Track #eta","Postive Track DCA (#mu m)","Negative Track DCA (#mu m)","B_{s}^{0} (K^{+} K^{-} #pi^{+}  #pi^{-}) (GeV/c^{2})","m_{KK}","m_{#piKK}","B_{s}^{0} Vertex X (#mu m)","B_{s}^{0} Vertex Y (#mu m)","B_{s}^{0} Vertex Z (#mu m)","D_{s}^{+} Vertex X (#mu m)","D_{s}^{+} Vertex Y (#mu m)","D_{s}^{+} Vertex Z (#mu m)","D_{s}^{+} DCA Daughters (#mu m)","D_{s}^{+} SV to PV (Decay) Length (#mu m)","cos(#theta_{D})","B_{s}^{0} Proper Lifetime (#mu m/c)","BsPt (GeV/c)","BsP (GeV/c)"};
	double XMin[NVar]={0,0,0,-1,0,0,0,0,-1.2,-1.2,0,0,5,1.00,1.85,-4000,-4000,-4000,-4000,-4000,-4000,0,0,-1,0,0,0};
	double XMax[NVar]={10,280,10000 * 5,1,0.1,0.2,8,8,1.2,1.2,4000,4000,6,1.04,2.05,4000,4000,4000,4000,4000,4000,300,6000 * 5,1,2500,25,30};

	TString SaveName[NVar] = {"BsMass","DCADaughters","SVPVDis","Angle","PhiMassDis","DsMassDis","PiPPt","PiMPt","PiPEta","PiMEta","PIPDCA","PIMDCA","BsMassZoom","KKMass","PiKKMass","BVX","BVY","BVZ","DVX","DVY","DVZ","DCADaughters_Ds","SVPVDis_Ds","DsAngle","BsLifeTime","BsPt","BsP"};
	*/

	const int NVar = 35;

	TString XAxisName[NVar] = {"B_{s}^{0} (K^{+} K^{-} #pi^{+}  #pi^{-}) (GeV/c^{2})","DCA Daughters (#mu m)","SV to PV (Decay) Length (#mu m)","cos(#theta)","|m_{KK} - m_{#phi}| (GeV/c^2)","|m_{#piKK} - m_{D_{s}}| (GeV/c^2)","Postive Track p_{T} (GeV/c)","Negative Track p_{T} (GeV/c)","Postive Track #eta","Negative Track #eta","Postive Track DCA (#mu m)","Negative Track DCA (#mu m)","B_{s}^{0} (K^{+} K^{-} #pi^{+}  #pi^{-}) (GeV/c^{2})","m_{KK}","m_{#piKK}","B_{s}^{0} Vertex X (#mu m)","B_{s}^{0} Vertex Y (#mu m)","B_{s}^{0} Vertex Z (#mu m)","D_{s}^{+} Vertex X (#mu m)","D_{s}^{+} Vertex Y (#mu m)","D_{s}^{+} Vertex Z (#mu m)","D_{s}^{+} DCA Daughters (#mu m)","D_{s}^{+} SV to PV (Decay) Length (#mu m)","cos(#theta_{D})","B_{s}^{0} Proper Lifetime (#mu m/c)","BsPt (GeV/c)","BsP (GeV/c)","First Mother Of #pi^{+}","First Mother Of #pi^{-}","First Mother Of K^{+}","First Mother Of K^{-}","GRAND Mother Of #pi^{+}","GRAND Mother Of #pi^{-}","GRAND Mother Of K^{+}","GRAND Mother Of K^{-}"};
	double XMin[NVar]={0,0,0,-1,0,0,0,0,-1.2,-1.2,0,0,5,1.00,1.85,-4000,-4000,-4000,-4000,-4000,-4000,0,0,-1,0,0,0,0,0,0,0,0,0,0,0};
	double XMax[NVar]={10,280,1000 * 5,1,0.1,0.2,8,8,1.2,1.2,4000,4000,6,1.04,2.05,4000,4000,4000,4000,4000,4000,300,1000,1,2500,25,30,560,560,560,560,800,800,800,800};

	TString SaveName[NVar] = {"BsMass","DCADaughters","SVPVDis","Angle","PhiMassDis","DsMassDis","PiPPt","PiMPt","PiPEta","PiMEta","PIPDCA","PIMDCA","BsMassZoom","KKMass","PiKKMass","BVX","BVY","BVZ","DVX","DVY","DVZ","DCADaughters_Ds","SVPVDis_Ds","DsAngle","BsLifeTime","BsPt","BsP","MomPIDPIP","MomPIDPIM","MomPIDKP","MomPIDKM","MotherPIDPIP","MotherPIDPIM","MotherPIDKP","MotherPIDKM"};


	TString VarName[NVar] = {"BsMassGen","dcadaughterGen","decayLength_BsGen","cosTheta_BsGen","abs(KKMassGen - 1.019461)","abs(PiKKMassGen - 1.96847)","pipPtGen","pimPtGen","pipEtaGen","pimEtaGen","pipDcaGen","pimDcaGen","BsMassGen","KKMassGen","PiKKMassGen","BvtxxGen","BvtxyGen","BvtxzGen","DvtxxGen","DvtxyGen","DvtxzGen","dcadaughter_DsGen","decayLength_DsGen","cosTheta_DsGen","BsLifeTimeGen","BsPtGen","BsPGen","MomPIDPIPGen","MomPIDPIMGen","MomPIDKPGen","MomPIDKMGen","MotherPIDPIPGen","MotherPIDPIMGen","MotherPIDKPGen","MotherPIDKMGen"};


	TLegend * leg[NVar];

	//All Files 

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	
	TString OutName;

	TString AllName;
	TString SigName;
	TString SelName;
	TString CombName;
	TString DecayName;
	TString DecayIDName;

	TFile * fout = new TFile("AnalyzedFileGen.root","RECREATE");
	fout->cd();

	TLatex *latE = new TLatex();
	latE->SetNDC();
	latE->SetTextSize(0.05);
	latE->SetTextColor(kBlack);


	for(int i = 2; i < 3; i++){



		cout << "Now Working on Variable:  " << SaveName[i].Data() << endl;

		AllName = Form("%sAll",SaveName[i].Data());
		SigName = Form("%sSig",SaveName[i].Data());
		SelName = Form("%sSel",SaveName[i].Data());
		CombName = Form("%sComb",SaveName[i].Data());
		DecayName = Form("%sDecay",SaveName[i].Data());
		DecayIDName = Form("%sDecayID",SaveName[i].Data());

		//Inclusive Stuffs//

		TH1D * BsHisAll = new TH1D(AllName.Data(),"",100,XMin[i],XMax[i]); 
		BsHisAll->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisAll->GetYaxis()->SetTitle("Normalized Counts");
		BsHisAll->GetXaxis()->CenterTitle();
		BsHisAll->GetYaxis()->CenterTitle();
		BsHisAll->SetLineWidth(2);
		BsHisAll->SetLineColor(kBlack);
		



		BsDsKKPiPi->Project(AllName.Data(),VarName[i].Data());
		BsHisAll->Write();

		BsHisAll->Sumw2();
		BsHisAll->Scale(1.0/BsHisAll->Integral());


		//Combinatorial Background//
		
		TH1D * BsHisComb = new TH1D(CombName.Data(),"",100,XMin[i],XMax[i]); 
		BsHisComb->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisComb->GetYaxis()->SetTitle("Normalized Counts");
		BsHisComb->GetXaxis()->CenterTitle();
		BsHisComb->GetYaxis()->CenterTitle();
		BsHisComb->SetLineWidth(2);
		BsHisComb->SetLineColor(kOrange);



		BsDsKKPiPi->Project(CombName.Data(),VarName[i].Data(),BsCombCut.Data());
		
		BsHisComb->Write();

		BsHisComb->Sumw2();
		BsHisComb->Scale(1.0/BsHisComb->Integral());


		//Bs Decay Background//
		
		TH1D * BsHisDecay = new TH1D(DecayName.Data(),"",100,XMin[i],XMax[i]); 
		BsHisDecay->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisDecay->GetYaxis()->SetTitle("Normalized Counts");
		BsHisDecay->GetXaxis()->CenterTitle();
		BsHisDecay->GetYaxis()->CenterTitle();
		BsHisDecay->SetLineWidth(2);
		BsHisDecay->SetLineColor(kGreen);


	
		BsDsKKPiPi->Project(DecayName.Data(),VarName[i].Data(),BsDecayCut.Data());

		BsHisDecay->Write();

		BsHisDecay->Sumw2();
		BsHisDecay->Scale(1.0/BsHisDecay->Integral());


		//Bs Decay ID Matched Required//
	

		TH1D * BsHisDecayID = new TH1D(DecayIDName.Data(),"",100,XMin[i],XMax[i]); 
		BsHisDecayID->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisDecayID->GetYaxis()->SetTitle("Normalized Counts");
		BsHisDecayID->GetXaxis()->CenterTitle();
		BsHisDecayID->GetYaxis()->CenterTitle();
		BsHisDecayID->SetLineWidth(2);
		BsHisDecayID->SetLineColor(kPink);


	
		BsDsKKPiPi->Project(DecayIDName.Data(),VarName[i].Data(),BsDecayIDCut.Data());

		BsHisDecayID->Write();

		BsHisDecayID->Sumw2();
		BsHisDecayID->Scale(1.0/BsHisDecayID->Integral());

		//Bs Signal//


		TH1D * BsHisSig = new TH1D(SigName.Data(),"",100,XMin[i],XMax[i]); 
		BsHisSig->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisSig->GetYaxis()->SetTitle("Normalized Counts");
		BsHisSig->GetXaxis()->CenterTitle();
		BsHisSig->GetYaxis()->CenterTitle();
		BsHisSig->SetLineWidth(2);
		BsHisSig->SetLineColor(kRed);

		BsDsKKPiPi->Project(SigName.Data(),VarName[i].Data(),BsSigCut.Data());

		BsHisSig->Write();


		cout << "Total Signal = " << BsHisSig->Integral() << endl;

		BsHisSig->Sumw2();
		BsHisSig->Scale(1.0/BsHisSig->Integral());

		


		TH1D * BsHisSel = new TH1D(SelName.Data(),"",100,XMin[i],XMax[i]); 
		BsHisSel->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisSel->GetYaxis()->SetTitle("Normalized Counts");
		BsHisSel->GetXaxis()->CenterTitle();
		BsHisSel->GetYaxis()->CenterTitle();
		BsHisSel->SetLineWidth(2);
		BsHisSel->SetLineColor(kBlue);

		BsDsKKPiPi->Project(SelName.Data(),VarName[i].Data(),BsSelCut.Data());

		BsHisSel->Write();

		BsHisSel->Sumw2();
		BsHisSel->Scale(1.0/BsHisSel->Integral());


		if(i != 2){
			BsHisSig->Draw("hist");	
			BsHisAll->Draw("histSAME");
			BsHisSel->Draw("histSAME");
			BsHisComb->Draw("histSAME");
			BsHisDecay->Draw("histSAME");
		}

		if(i == 2){
			c->SetLogy();	
			BsHisAll->SetMaximum(0.5);
			BsHisAll->GetXaxis()->SetMaxDigits(3);
			
			BsHisAll->Draw("hist");
			BsHisSig->Draw("histSAME");	
			BsHisSel->Draw("histSAME");
			BsHisComb->Draw("histSAME");
			BsHisDecay->Draw("histSAME");

		}




		leg[i] = new TLegend(0.30,0.65,0.75,0.85,NULL,"brNDC");
		leg[i]->SetBorderSize(0);
		leg[i]->SetTextSize(0.040);
		leg[i]->SetTextFont(42);
		leg[i]->SetFillStyle(0);
		leg[i]->SetLineWidth(3);

		leg[i]->AddEntry(BsHisAll,"Signal + Background","l");
		leg[i]->AddEntry(BsHisSig,"B_{s}^{0} Signal Decay","l");
		leg[i]->AddEntry(BsHisSel,"Mass Window Selection Applied","l");
		leg[i]->AddEntry(BsHisComb,"Random Background from PV","l");
		leg[i]->AddEntry(BsHisDecay,"PYTHIA 8 pp #rightarrow b #bar{b} Decay","l");


		leg[i]->Draw("SAME");



		latE->DrawLatex(0.42,0.87,"#it{#bf{GEN (Truth)}}");



		OutName = Form("PlotsInclusive/Gen/%s.png",SaveName[i].Data());

		c->SaveAs(OutName.Data());


	}

	//Resolution Plots//

	BsDsKKPiPi->AddFriend("BsDsKKPiPi");
	const int ResoNVar = 6;

	TString ResoXAxisName[ResoNVar] = {"B_{s}^{0} Vertex X: RECO - GEN (#mu m)","B_{s}^{0} Vertex Y: RECO - GEN (#mu m)","B_{s}^{0} Vertex Z: RECO - GEN (#mu m)","D_{s}^{+} Vertex X: RECO - GEN (#mu m)","D_{s}^{+} Vertex Y: RECO - GEN (#mu m)","D_{s}^{+} Vertex Z: RECO - GEN (#mu m)"};
	double ResoXMin[ResoNVar]={-1500,-1500,-1500,-1500,-1500,-1500};
	double ResoXMax[ResoNVar]={1500,1500,1500,1500,1500,1500};

	TString ResoSaveName[ResoNVar] = {"BVXRESO","BVYRESO","BVZRESO","DVXRESO","DVYRESO","DVZRESO"};

	TString ResoVarName[ResoNVar] = {"Bvtxx - BvtxxGen","Bvtxy - BvtxyGen","Bvtxz - BvtxzGen","Dvtxx - DvtxxGen","Dvtxy - DvtxyGen","Dvtxz - DvtxzGen"};
	TLegend * Resoleg[ResoNVar];



	

	for(int i = 0; i < ResoNVar; i++){



		cout << "Now Working on Variable:  " << ResoSaveName[i].Data() << endl;

		AllName = Form("%sAll",ResoSaveName[i].Data());
		SigName = Form("%sSig",ResoSaveName[i].Data());
		SelName = Form("%sSel",ResoSaveName[i].Data());
		CombName = Form("%sComb",ResoSaveName[i].Data());
		DecayName = Form("%sDecay",ResoSaveName[i].Data());

		//Inclusive Stuffs//

		TH1D * BsHisAll = new TH1D(AllName.Data(),"",100,ResoXMin[i],ResoXMax[i]); 
		BsHisAll->GetXaxis()->SetTitle(ResoXAxisName[i].Data());
		BsHisAll->GetYaxis()->SetTitle("Normalized Counts");
		BsHisAll->GetXaxis()->CenterTitle();
		BsHisAll->GetYaxis()->CenterTitle();
		BsHisAll->SetLineWidth(2);
		BsHisAll->SetLineColor(kBlack);
		



		BsDsKKPiPi->Project(AllName.Data(),ResoVarName[i].Data());
		BsHisAll->Write();

		BsHisAll->Sumw2();
		BsHisAll->Scale(1.0/BsHisAll->Integral());


		//Combinatorial Background//
		
		TH1D * BsHisComb = new TH1D(CombName.Data(),"",100,ResoXMin[i],ResoXMax[i]); 
		BsHisComb->GetXaxis()->SetTitle(ResoXAxisName[i].Data());
		BsHisComb->GetYaxis()->SetTitle("Normalized Counts");
		BsHisComb->GetXaxis()->CenterTitle();
		BsHisComb->GetYaxis()->CenterTitle();
		BsHisComb->SetLineWidth(2);
		BsHisComb->SetLineColor(kOrange);



		BsDsKKPiPi->Project(CombName.Data(),ResoVarName[i].Data(),BsCombCut.Data());
		
		BsHisComb->Write();

		BsHisComb->Sumw2();
		BsHisComb->Scale(1.0/BsHisComb->Integral());


		//Bs Decay Background//
		
		TH1D * BsHisDecay = new TH1D(DecayName.Data(),"",100,ResoXMin[i],ResoXMax[i]); 
		BsHisDecay->GetXaxis()->SetTitle(ResoXAxisName[i].Data());
		BsHisDecay->GetYaxis()->SetTitle("Normalized Counts");
		BsHisDecay->GetXaxis()->CenterTitle();
		BsHisDecay->GetYaxis()->CenterTitle();
		BsHisDecay->SetLineWidth(2);
		BsHisDecay->SetLineColor(kGreen);


	
		BsDsKKPiPi->Project(DecayName.Data(),ResoVarName[i].Data(),BsDecayCut.Data());

		BsHisDecay->Write();

		BsHisDecay->Sumw2();
		BsHisDecay->Scale(1.0/BsHisDecay->Integral());

		//Bs Signal//


		TH1D * BsHisSig = new TH1D(SigName.Data(),"",100,ResoXMin[i],ResoXMax[i]); 
		BsHisSig->GetXaxis()->SetTitle(ResoXAxisName[i].Data());
		BsHisSig->GetYaxis()->SetTitle("Normalized Counts");
		BsHisSig->GetXaxis()->CenterTitle();
		BsHisSig->GetYaxis()->CenterTitle();
		BsHisSig->SetLineWidth(2);
		BsHisSig->SetLineColor(kRed);

		BsDsKKPiPi->Project(SigName.Data(),ResoVarName[i].Data(),BsSigCut.Data());

		BsHisSig->Write();


		cout << "Total Signal = " << BsHisSig->Integral() << endl;

		BsHisSig->Sumw2();
		BsHisSig->Scale(1.0/BsHisSig->Integral());

		


		TH1D * BsHisSel = new TH1D(SelName.Data(),"",100,ResoXMin[i],ResoXMax[i]); 
		BsHisSel->GetXaxis()->SetTitle(ResoXAxisName[i].Data());
		BsHisSel->GetYaxis()->SetTitle("Normalized Counts");
		BsHisSel->GetXaxis()->CenterTitle();
		BsHisSel->GetYaxis()->CenterTitle();
		BsHisSel->SetLineWidth(2);
		BsHisSel->SetLineColor(kBlue);

		BsDsKKPiPi->Project(SelName.Data(),ResoVarName[i].Data(),BsSelCut.Data());

		BsHisSel->Write();

		BsHisSel->Sumw2();
		BsHisSel->Scale(1.0/BsHisSel->Integral());


		if(i != 2){
			BsHisSig->Draw("hist");	
			BsHisAll->Draw("histSAME");
			BsHisSel->Draw("histSAME");
			BsHisComb->Draw("histSAME");
			BsHisDecay->Draw("histSAME");
		}

		if(i == 2){
			BsHisAll->Draw("hist");
			BsHisSig->Draw("histSAME");	
			BsHisSel->Draw("histSAME");
			BsHisComb->Draw("histSAME");
			BsHisDecay->Draw("histSAME");

		}




		Resoleg[i] = new TLegend(0.36,0.65,0.75,0.85,NULL,"brNDC");
		Resoleg[i]->SetBorderSize(0);
		Resoleg[i]->SetTextSize(0.040);
		Resoleg[i]->SetTextFont(42);
		Resoleg[i]->SetFillStyle(0);
		Resoleg[i]->SetLineWidth(3);

		Resoleg[i]->AddEntry(BsHisAll,"Signal + Background","l");
		Resoleg[i]->AddEntry(BsHisSig,"B_{s}^{0} Signal Decay","l");
		Resoleg[i]->AddEntry(BsHisSel,"Mass Window Selection Applied","l");
		Resoleg[i]->AddEntry(BsHisComb,"Random Background from PV","l");
		Resoleg[i]->AddEntry(BsHisDecay,"PYTHIA 8 pp #rightarrow b #bar{b} Decay","l");


		Resoleg[i]->Draw("SAME");


		OutName = Form("PlotsInclusive/RESOPlots/%s.png",ResoSaveName[i].Data());

		c->SaveAs(OutName.Data());


	}


	fout->Close();




}
