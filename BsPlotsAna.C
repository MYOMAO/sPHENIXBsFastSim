#ifndef __CINT__
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
#endif


void BsPlotsAna(){


	gStyle->SetOptStat(0);
	SetsPhenixStyle();

	gStyle->SetPadRightMargin(1.2);

	TString BsSigCut = "(IsSigPiP*IsSigPiM*IsSigKP*IsSigKM==1)";
	TString BsSelCut = "(abs(KKMass - 1.019461) < 0.1 && abs(PiKKMass - 1.96847) < 0.10)";


	TString infile = "/sphenix/user/zshi/sPHENIXHF/Bs/Simulations/tutorials/AnaTutorial/src/condor/directory/macros/output.root";

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TTree * BsDsKKPiPi = (TTree *) fin->Get("BsDsKKPiPi");

	const int NVar = 8;

	TString XAxisName[NVar] = {"B_{s}^{0} (K^{+} K^{-} #pi^{+}  #pi^{-}) (GeV/c^{2})","DCA Daughters (#mu m)","SV to PV (Decay) Length (#mu m)","cos(#theta)","|m_{KK} - m_{#phi}| (GeV/c^2)","|m_{#piKK} - m_{D_{s}}| (GeV/c^2)","Positive Track p_{T} (GeV/c)","Negative Track p_{T} (GeV/c)"};
	double XMin[NVar]={0,0,0,-1,0,0,0,0};
	double XMax[NVar]={10,1000,1000,1,2,2,20,20};


	TString SaveName[NVar] = {"BsMass","DCADaughters","SVPVDis","Angle","PhiMassDis","DsMassDis","PlusPt","MinusPt"};
	

	TString VarName[NVar] = {"BsMass","dcadaughter","decayLength_Bs","cosTheta_Bs","abs(KKMass - 1.019461)","abs(PiKKMass - 1.96847)","pipPt","pimPt"};

	TLegend * leg[NVar];

	//All Files 

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	TString OutName;
		



	for(int i = 0; i < NVar; i++){
	

		cout << "Now Working on Variable:  " << SaveName[i].Data() << endl;

		TH1D * BsHisAll = new TH1D("BsHisAll","",100,XMin[i],XMax[i]); 
		BsHisAll->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisAll->GetYaxis()->SetTitle("Normalized Counts");
		BsHisAll->GetXaxis()->CenterTitle();
		BsHisAll->GetYaxis()->CenterTitle();
		BsHisAll->SetLineWidth(1);
		BsHisAll->SetLineColor(kBlack);



		BsDsKKPiPi->Project("BsHisAll",VarName[i].Data());

		BsHisAll->Sumw2();
		BsHisAll->Scale(1.0/BsHisAll->Integral());


		TH1D * BsHisSig = new TH1D("BsHisSig","",100,XMin[i],XMax[i]); 
		BsHisSig->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisSig->GetYaxis()->SetTitle("Normalized Counts");
		BsHisSig->GetXaxis()->CenterTitle();
		BsHisSig->GetYaxis()->CenterTitle();
		BsHisSig->SetLineWidth(1);
		BsHisSig->SetLineColor(kRed);

		BsDsKKPiPi->Project("BsHisSig",VarName[i].Data(),BsSigCut.Data());


		BsHisSig->Sumw2();
		BsHisSig->Scale(1.0/BsHisSig->Integral());

		
		TH1D * BsHisSel = new TH1D("BsHisSel","",100,XMin[i],XMax[i]); 
		BsHisSel->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisSel->GetYaxis()->SetTitle("Normalized Counts");
		BsHisSel->GetXaxis()->CenterTitle();
		BsHisSel->GetYaxis()->CenterTitle();
		BsHisSel->SetLineWidth(1);
		BsHisSel->SetLineColor(kRed);

		BsDsKKPiPi->Project("BsHisSel",VarName[i].Data(),BsSelCut.Data());


		BsHisSel->Sumw2();
		BsHisSel->Scale(1.0/BsHisSel->Integral());



		BsHisSig->Draw("hist");
		BsHisAll->Draw("histSAME");
		BsHisSel->Draw("histSAME");

		
		


		leg[i] = new TLegend(0.36,0.65,0.75,0.85,NULL,"brNDC");
		leg[i]->SetBorderSize(0);
		leg[i]->SetTextSize(0.040);
		leg[i]->SetTextFont(42);
		leg[i]->SetFillStyle(0);
		leg[i]->SetLineWidth(3);

		leg[i]->AddEntry(BsHisAll,"Signal + Background","l");
		leg[i]->AddEntry(BsHisSel,"Mass Window Selection Applied","l");
		leg[i]->AddEntry(BsHisSig,"B_{s}^{0} Signal Decay","l");

		leg[i]->Draw("SAME");
		
		
		OutName = Form("Plots/%s.png",SaveName[i].Data());

		c->SaveAs(OutName.Data());


	}







}
