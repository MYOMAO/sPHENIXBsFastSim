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


void BsPlotsAna(){

	int UseTrackptCut = 1;

	double	TrackptCutValue =  0.2;

	TString TrackptCut = Form("(pipPt > %f && pimPt > %f && kpPt > %f && kmPt > %f)",TrackptCutValue,TrackptCutValue,TrackptCutValue,TrackptCutValue);


	gStyle->SetOptStat(0);
	SetsPhenixStyle();

	gStyle->SetPadRightMargin(1.2);

	//Cuts//

	TString BsSigCut = "(IsSigPiP*IsSigPiM*IsSigKP*IsSigKM==1)";
	TString BsSelCut = "(abs(KKMass - 1.019461) < 0.1 && abs(PiKKMass - 1.96847) < 0.10)";

	TString BsDecayCut = "(IsSigPiP>2 && IsSigPiM>2 && IsSigKP>2 && IsSigKM>2)";
	TString BsCombCut = "(IsSigPiP==0 && IsSigPiM==0 && IsSigKP==0 && IsSigKM==0)";



	if(UseTrackptCut == 1){

		BsSigCut = Form("(IsSigPiP*IsSigPiM*IsSigKP*IsSigKM==1) && %s",TrackptCut.Data());
		BsSelCut = Form("(abs(KKMass - 1.019461) < 0.1 && abs(PiKKMass - 1.96847) < 0.10) && %s",TrackptCut.Data());

	}

	TString infile = "output.root";

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TTree * BsDsKKPiPi = (TTree *) fin->Get("BsDsKKPiPi");


	const int NVar = 24;

	TString XAxisName[NVar] = {"B_{s}^{0} (K^{+} K^{-} #pi^{+}  #pi^{-}) (GeV/c^{2})","DCA Daughters (#mu m)","SV to PV (Decay) Length (#mu m)","cos(#theta)","|m_{KK} - m_{#phi}| (GeV/c^2)","|m_{#piKK} - m_{D_{s}}| (GeV/c^2)","Postive Track p_{T} (GeV/c)","Negative Track p_{T} (GeV/c)","Postive Track #eta","Negative Track #eta","Postive Track DCA (#mu m)","Negative Track DCA (#mu m)","B_{s}^{0} (K^{+} K^{-} #pi^{+}  #pi^{-}) (GeV/c^{2})","m_{KK}","m_{#piKK}","B_{s}^{0} Vertex X (#mu m)","B_{s}^{0} Vertex Y (#mu m)","B_{s}^{0} Vertex Z (#mu m)","D_{s}^{0} Vertex X (#mu m)","D_{s}^{0} Vertex Y (#mu m)","D_{s}^{0} Vertex Z (#mu m)"};
	double XMin[NVar]={0,0,0,-1,0,0,0,0,-1.2,-1.2,0,0,5,1.00,1.85,-4000,-4000,-4000,-4000,-4000,-4000,0,0,-1};
	double XMax[NVar]={10,280,10000,1,0.1,0.2,8,8,1.2,1.2,4000,4000,6,1.04,2.05,4000,4000,4000,4000,4000,4000,300,6000,1};

	TString SaveName[NVar] = {"BsMass","DCADaughters","SVPVDis","Angle","PhiMassDis","DsMassDis","PiPPt","PiMPt","PiPEta","PiMEta","PIPDCA","PIMDCA","BsMassZoom","KKMass","PiKKMass","BVX","BVY","BVZ","DVX","DVY","DVZ","DCADaughters_Ds","SVPVDis_Ds","DsAngle"};
	TString VarName[NVar] = {"BsMass","dcadaughter","decayLength_Bs","cosTheta_Bs","abs(KKMass - 1.019461)","abs(PiKKMass - 1.96847)","pipPt","pimPt","pipEta","pimEta","pipDca","pimDca","BsMass","KKMass","PiKKMass","Bvtxx","Bvtxy","Bvtxz","Dvtxx","Dvtxy","Dvtxz","dcadaughter_Ds","decayLength_Ds","cosTheta_Ds"};

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


	TFile * fout = new TFile("AnalyzedFile.root","RECREATE");
	fout->cd();

	for(int i = 0; i < NVar; i++){


		cout << "Now Working on Variable:  " << SaveName[i].Data() << endl;

		AllName = Form("%sAll",SaveName[i].Data());
		SigName = Form("%sSig",SaveName[i].Data());
		SelName = Form("%sSel",SaveName[i].Data());
		CombName = Form("%sComb",SaveName[i].Data());
		DecayName = Form("%sDecay",SaveName[i].Data());

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
			BsHisAll->Draw("hist");
			BsHisSig->Draw("histSAME");	
			BsHisSel->Draw("histSAME");
			BsHisComb->Draw("histSAME");
			BsHisDecay->Draw("histSAME");

		}




		leg[i] = new TLegend(0.36,0.65,0.75,0.85,NULL,"brNDC");
		leg[i]->SetBorderSize(0);
		leg[i]->SetTextSize(0.040);
		leg[i]->SetTextFont(42);
		leg[i]->SetFillStyle(0);
		leg[i]->SetLineWidth(3);

		leg[i]->AddEntry(BsHisAll,"Signal + Background","l");
		leg[i]->AddEntry(BsHisSig,"B_{s}^{0} Signal Decay","l");
		leg[i]->AddEntry(BsHisSel,"Mass Window Selection Applied","l");
		leg[i]->AddEntry(BsHisComb,"Random Background from PV","l");
		leg[i]->AddEntry(BsHisDecay,"Inclusive B_{s}^{0} Decay","l");


		leg[i]->Draw("SAME");


		OutName = Form("PlotsInclusive/RECO/%s.png",SaveName[i].Data());

		c->SaveAs(OutName.Data());


	}


	fout->Close();




}
