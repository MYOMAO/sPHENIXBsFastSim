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


void PostAna(){

	gStyle->SetOptStat(0);
	SetsPhenixStyle();

	TString infile = "AnalyzedFileAll.root";

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TString infileGen = "AnalyzedFileGenAll.root";
	TFile * finGen = new TFile(infileGen.Data());


	const int NVar = 24;
	TString XAxisName[NVar] = {"B_{s}^{0} (K^{+} K^{-} e^{+}  e^{-}) (GeV/c^{2})","DCA Daughters (#mu m)","SV to PV (Decay) Length (#mu m)","cos(#theta)","|m_{KK} - m_{#phi}| (GeV/c^2)","|m_{ee} - m_{J/#psi}}| (GeV/c^2)","K^{+} p_{T} (GeV/c)","K^{-} p_{T} (GeV/c)","K^{+} #eta","K^{-} #eta","K^{+} DCA (#mu m)","K^{-} DCA (#mu m)","e^{+} p_{T} (GeV/c)","e^{-} p_{T} (GeV/c)","e^{+} #eta","e^{-} #eta","e^{+} DCA (#mu m)","e^{-} DCA (#mu m)","B_{s}^{0} (K^{+} K^{-} e^{+}  e^{-}) (GeV/c^{2})","m_{KK}","m_{ee}","B_{s}^{0} Vertex X (#mu m)","B_{s}^{0} Vertex Y (#mu m)","B_{s}^{0} Vertex Z (#mu m)"};
	double XMin[NVar]={0,0,0,-1,0,0,0,0,-1.2,-1.2,0,0,0,0,-1.2,-1.2,0,0,5,1.00,2.90,-4000,-4000,-4000};
	double XMax[NVar]={10,280,10000,1,0.1,0.2,8,8,1.2,1.2,4000,4000,8,8,1.2,1.2,4000,4000,6,1.04,3.30,4000,4000,4000};

	TString SaveName[NVar] = {"BsMass","DCADaughters","SVPVDis","Angle","PhiMassDis","eeMassDis","KPPt","KMPt","KPEta","KMEta","KPDCA","KMDCA","EPPt","EMPt","EPEta","EMEta","EPDCA","EMDCA","BsMassZoom","KKMass","eeMass","BVX","BVY","BVZ"};


	//	TString VarName[NVar] = {"BsMassGen","dcadaughterGen","decayLength_BsGen","cosTheta_BsGen","abs(KKMassGen - 1.019461)","abs(eeMassGen - 1.96847)","epPtGen","emPtGen","epEtaGen","emEtaGen","epDcaGen","emDcaGen","BsMassGen","KKMassGen","eeMassGen","BvtxxGen","BvtxyGen","BvtxzGen"};
	
	TLegend * leg[NVar];
	TString OutName;

	TString AllName;
	TString SigName;
	TString SelName;
	TString CombName;
	TString DecayName;

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();


	
	for(int i = 0; i < NVar; i++){


		cout << "Now Working on Variable:  " << SaveName[i].Data() << endl;

		AllName = Form("%sAll",SaveName[i].Data());
		SigName = Form("%sSig",SaveName[i].Data());
		SelName = Form("%sSel",SaveName[i].Data());
		CombName = Form("%sComb",SaveName[i].Data());
		DecayName = Form("%sDecay",SaveName[i].Data());

		
		TH1D * BsHisAll = (TH1D *)  fin->Get(AllName.Data());
		TH1D * BsHisSig = (TH1D *)  fin->Get(SigName.Data());
		TH1D * BsHisSel = (TH1D *)  fin->Get(SelName.Data());
		TH1D * BsHisComb = (TH1D *)  fin->Get(CombName.Data());
		TH1D * BsHisDecay = (TH1D *)  fin->Get(DecayName.Data());

	
		BsHisAll->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisAll->GetYaxis()->SetTitle("Normalized Counts");
		BsHisAll->GetXaxis()->CenterTitle();
		BsHisAll->GetYaxis()->CenterTitle();
		BsHisAll->GetXaxis()->SetTitleOffset(1.2);
		BsHisAll->GetYaxis()->SetTitleOffset(1.2);

		BsHisAll->Scale(1.0/BsHisAll->Integral());
		BsHisAll->SetLineWidth(2);
		BsHisAll->SetLineColor(kBlack);
		
		BsHisSig->Scale(1.0/BsHisSig->Integral());
		BsHisSig->SetLineWidth(2);
		BsHisSig->SetLineColor(kRed);


		BsHisComb->Scale(1.0/BsHisComb->Integral());
		BsHisComb->SetLineWidth(2);
		BsHisComb->SetLineColor(kOrange);


		BsHisDecay->Scale(1.0/BsHisDecay->Integral());
		BsHisDecay->SetLineWidth(2);
		BsHisDecay->SetLineColor(kGreen);


		BsHisSel->Scale(1.0/BsHisSel->Integral());
		BsHisSel->SetLineWidth(2);
		BsHisSel->SetLineColor(kBlue);

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


	TLegend * legGen[NVar];

	for(int i = 0; i < NVar; i++){


		cout << "Now Working on Variable:  " << SaveName[i].Data() << endl;

		AllName = Form("%sAll",SaveName[i].Data());
		SigName = Form("%sSig",SaveName[i].Data());
		SelName = Form("%sSel",SaveName[i].Data());
		CombName = Form("%sComb",SaveName[i].Data());
		DecayName = Form("%sDecay",SaveName[i].Data());

	
		TH1D * BsHisAll = (TH1D *)  finGen->Get(AllName.Data());
		TH1D * BsHisSig = (TH1D *)  finGen->Get(SigName.Data());
		TH1D * BsHisSel = (TH1D *)  finGen->Get(SelName.Data());
		TH1D * BsHisComb = (TH1D *)  finGen->Get(CombName.Data());
		TH1D * BsHisDecay = (TH1D *)  finGen->Get(DecayName.Data());

	
		BsHisAll->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisAll->GetYaxis()->SetTitle("Normalized Counts");
		BsHisAll->GetXaxis()->CenterTitle();
		BsHisAll->GetYaxis()->CenterTitle();
		BsHisAll->GetXaxis()->SetTitleOffset(1.2);
		BsHisAll->GetYaxis()->SetTitleOffset(1.2);

		BsHisAll->Scale(1.0/BsHisAll->Integral());
		BsHisAll->SetLineWidth(2);
		BsHisAll->SetLineColor(kBlack);
		
		BsHisSig->Scale(1.0/BsHisSig->Integral());
		BsHisSig->SetLineWidth(2);
		BsHisSig->SetLineColor(kRed);


		BsHisComb->Scale(1.0/BsHisComb->Integral());
		BsHisComb->SetLineWidth(2);
		BsHisComb->SetLineColor(kOrange);


		BsHisDecay->Scale(1.0/BsHisDecay->Integral());
		BsHisDecay->SetLineWidth(2);
		BsHisDecay->SetLineColor(kGreen);


		BsHisSel->Scale(1.0/BsHisSel->Integral());
		BsHisSel->SetLineWidth(2);
		BsHisSel->SetLineColor(kBlue);

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


		legGen[i] = new TLegend(0.36,0.65,0.75,0.85,NULL,"brNDC");
		legGen[i]->SetBorderSize(0);
		legGen[i]->SetTextSize(0.040);
		legGen[i]->SetTextFont(42);
		legGen[i]->SetFillStyle(0);
		legGen[i]->SetLineWidth(3);

		legGen[i]->AddEntry(BsHisAll,"Signal + Background","l");
		legGen[i]->AddEntry(BsHisSig,"B_{s}^{0} Signal Decay","l");
		legGen[i]->AddEntry(BsHisSel,"Mass Window Selection Applied","l");
		legGen[i]->AddEntry(BsHisComb,"Random Background from PV","l");
		legGen[i]->AddEntry(BsHisDecay,"Inclusive B_{s}^{0} Decay","l");
			
		legGen[i]->Draw("SAME");

		OutName = Form("PlotsInclusive/Gen/%s.png",SaveName[i].Data());

		c->SaveAs(OutName.Data());

	}



	TLegend * legAllComp[NVar];
	TLegend * legSigComp[NVar];
	TLegend * legSelComp[NVar];
	TLegend * legCombComp[NVar];
	TLegend * legDecayComp[NVar];

	

	for(int i = 0; i < NVar; i++){


		cout << "Now Working on Variable:  " << SaveName[i].Data() << endl;

		AllName = Form("%sAll",SaveName[i].Data());
		SigName = Form("%sSig",SaveName[i].Data());
		SelName = Form("%sSel",SaveName[i].Data());
		CombName = Form("%sComb",SaveName[i].Data());
		DecayName = Form("%sDecay",SaveName[i].Data());

	
		TH1D * BsHisAllCompGen = (TH1D *)  finGen->Get(AllName.Data());
		TH1D * BsHisSigCompGen = (TH1D *)  finGen->Get(SigName.Data());
		TH1D * BsHisSelCompGen = (TH1D *)  finGen->Get(SelName.Data());
		TH1D * BsHisCombCompGen = (TH1D *)  finGen->Get(CombName.Data());
		TH1D * BsHisDecayCompGen = (TH1D *)  finGen->Get(DecayName.Data());

	
		TH1D * BsHisAllComp = (TH1D *)  fin->Get(AllName.Data());
		TH1D * BsHisSigComp = (TH1D *)  fin->Get(SigName.Data());
		TH1D * BsHisSelComp = (TH1D *)  fin->Get(SelName.Data());
		TH1D * BsHisCombComp = (TH1D *)  fin->Get(CombName.Data());
		TH1D * BsHisDecayComp = (TH1D *)  fin->Get(DecayName.Data());


		BsHisAllCompGen->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisAllCompGen->GetYaxis()->SetTitle("Normalized Counts");
		BsHisAllCompGen->GetXaxis()->CenterTitle();
		BsHisAllCompGen->GetYaxis()->CenterTitle();
		BsHisAllCompGen->GetXaxis()->SetTitleOffset(1.2);
		BsHisAllCompGen->GetYaxis()->SetTitleOffset(1.2);

		BsHisAllCompGen->Scale(1.0/BsHisAllCompGen->Integral());
		BsHisAllCompGen->SetLineWidth(2);
		BsHisAllCompGen->SetLineColor(kBlue);
		


		BsHisAllComp->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisAllComp->GetYaxis()->SetTitle("Normalized Counts");
		BsHisAllComp->GetXaxis()->CenterTitle();
		BsHisAllComp->GetYaxis()->CenterTitle();
		BsHisAllComp->GetXaxis()->SetTitleOffset(1.2);
		BsHisAllComp->GetYaxis()->SetTitleOffset(1.2);

		BsHisAllComp->Scale(1.0/BsHisAllComp->Integral());
		BsHisAllComp->SetLineWidth(2);
		BsHisAllComp->SetLineColor(kRed);




		legAllComp[i] = new TLegend(0.36,0.65,0.75,0.85,NULL,"brNDC");
		legAllComp[i]->SetBorderSize(0);
		legAllComp[i]->SetTextSize(0.040);
		legAllComp[i]->SetTextFont(42);
		legAllComp[i]->SetFillStyle(0);
		legAllComp[i]->SetLineWidth(3);

		legAllComp[i]->AddEntry(BsHisAllCompGen,"GEN: Sig + Bkgd","l");
		legAllComp[i]->AddEntry(BsHisAllComp,"RECO: Sig + Bkgd","l");

		BsHisAllCompGen->Draw("hist");
		BsHisAllComp->Draw("histSAME");

		legAllComp[i]->Draw("SAME");

		OutName = Form("PlotsComparison/All/%s.png",SaveName[i].Data());

		c->SaveAs(OutName.Data());


		//Sig//

		BsHisSigCompGen->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisSigCompGen->GetYaxis()->SetTitle("Normalized Counts");
		BsHisSigCompGen->GetXaxis()->CenterTitle();
		BsHisSigCompGen->GetYaxis()->CenterTitle();
		BsHisSigCompGen->GetXaxis()->SetTitleOffset(1.2);
		BsHisSigCompGen->GetYaxis()->SetTitleOffset(1.2);

		BsHisSigCompGen->Scale(1.0/BsHisSigCompGen->Integral());
		BsHisSigCompGen->SetLineWidth(2);
		BsHisSigCompGen->SetLineColor(kBlue);
		


		BsHisSigComp->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisSigComp->GetYaxis()->SetTitle("Normalized Counts");
		BsHisSigComp->GetXaxis()->CenterTitle();
		BsHisSigComp->GetYaxis()->CenterTitle();
		BsHisSigComp->GetXaxis()->SetTitleOffset(1.2);
		BsHisSigComp->GetYaxis()->SetTitleOffset(1.2);

		BsHisSigComp->Scale(1.0/BsHisSigComp->Integral());
		BsHisSigComp->SetLineWidth(2);
		BsHisSigComp->SetLineColor(kRed);



		legSigComp[i] = new TLegend(0.36,0.65,0.75,0.85,NULL,"brNDC");
		legSigComp[i]->SetBorderSize(0);
		legSigComp[i]->SetTextSize(0.040);
		legSigComp[i]->SetTextFont(42);
		legSigComp[i]->SetFillStyle(0);
		legSigComp[i]->SetLineWidth(3);

		legSigComp[i]->AddEntry(BsHisSigCompGen,"GEN: Sig","l");
		legSigComp[i]->AddEntry(BsHisSigComp,"RECO: Sig","l");
	
		BsHisSigCompGen->Draw("hist");
		BsHisSigComp->Draw("histSAME");

		legSigComp[i]->Draw("SAME");

		OutName = Form("PlotsComparison/Sig/%s.png",SaveName[i].Data());

		c->SaveAs(OutName.Data());


		//Sel//

		BsHisSelCompGen->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisSelCompGen->GetYaxis()->SetTitle("Normalized Counts");
		BsHisSelCompGen->GetXaxis()->CenterTitle();
		BsHisSelCompGen->GetYaxis()->CenterTitle();
		BsHisSelCompGen->GetXaxis()->SetTitleOffset(1.2);
		BsHisSelCompGen->GetYaxis()->SetTitleOffset(1.2);

		BsHisSelCompGen->Scale(1.0/BsHisSelCompGen->Integral());
		BsHisSelCompGen->SetLineWidth(2);
		BsHisSelCompGen->SetLineColor(kBlue);
		


		BsHisSelComp->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisSelComp->GetYaxis()->SetTitle("Normalized Counts");
		BsHisSelComp->GetXaxis()->CenterTitle();
		BsHisSelComp->GetYaxis()->CenterTitle();
		BsHisSelComp->GetXaxis()->SetTitleOffset(1.2);
		BsHisSelComp->GetYaxis()->SetTitleOffset(1.2);

		BsHisSelComp->Scale(1.0/BsHisSelComp->Integral());
		BsHisSelComp->SetLineWidth(2);
		BsHisSelComp->SetLineColor(kRed);



		legSelComp[i] = new TLegend(0.36,0.65,0.75,0.85,NULL,"brNDC");
		legSelComp[i]->SetBorderSize(0);
		legSelComp[i]->SetTextSize(0.040);
		legSelComp[i]->SetTextFont(42);
		legSelComp[i]->SetFillStyle(0);
		legSelComp[i]->SetLineWidth(3);

		legSelComp[i]->AddEntry(BsHisSelCompGen,"GEN: Sel","l");
		legSelComp[i]->AddEntry(BsHisSelComp,"RECO: Sel","l");
	
	
		BsHisSelCompGen->Draw("hist");
		BsHisSelComp->Draw("histSAME");


		legSelComp[i]->Draw("SAME");

		OutName = Form("PlotsComparison/Sel/%s.png",SaveName[i].Data());

		c->SaveAs(OutName.Data());



		//Comb//

		BsHisCombCompGen->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisCombCompGen->GetYaxis()->SetTitle("Normalized Counts");
		BsHisCombCompGen->GetXaxis()->CenterTitle();
		BsHisCombCompGen->GetYaxis()->CenterTitle();
		BsHisCombCompGen->GetXaxis()->SetTitleOffset(1.2);
		BsHisCombCompGen->GetYaxis()->SetTitleOffset(1.2);

		BsHisCombCompGen->Scale(1.0/BsHisCombCompGen->Integral());
		BsHisCombCompGen->SetLineWidth(2);
		BsHisCombCompGen->SetLineColor(kBlue);
		


		BsHisCombComp->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisCombComp->GetYaxis()->SetTitle("Normalized Counts");
		BsHisCombComp->GetXaxis()->CenterTitle();
		BsHisCombComp->GetYaxis()->CenterTitle();
		BsHisCombComp->GetXaxis()->SetTitleOffset(1.2);
		BsHisCombComp->GetYaxis()->SetTitleOffset(1.2);

		BsHisCombComp->Scale(1.0/BsHisCombComp->Integral());
		BsHisCombComp->SetLineWidth(2);
		BsHisCombComp->SetLineColor(kRed);



		legCombComp[i] = new TLegend(0.36,0.65,0.75,0.85,NULL,"brNDC");
		legCombComp[i]->SetBorderSize(0);
		legCombComp[i]->SetTextSize(0.040);
		legCombComp[i]->SetTextFont(42);
		legCombComp[i]->SetFillStyle(0);
		legCombComp[i]->SetLineWidth(3);

		legCombComp[i]->AddEntry(BsHisCombCompGen,"GEN: Comb","l");
		legCombComp[i]->AddEntry(BsHisCombComp,"RECO: Comb","l");
	

		BsHisCombCompGen->Draw("hist");
		BsHisCombComp->Draw("histSAME");


		legCombComp[i]->Draw("SAME");

		OutName = Form("PlotsComparison/Comb/%s.png",SaveName[i].Data());

		c->SaveAs(OutName.Data());




		//Decay//

		BsHisDecayCompGen->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisDecayCompGen->GetYaxis()->SetTitle("Normalized Counts");
		BsHisDecayCompGen->GetXaxis()->CenterTitle();
		BsHisDecayCompGen->GetYaxis()->CenterTitle();
		BsHisDecayCompGen->GetXaxis()->SetTitleOffset(1.2);
		BsHisDecayCompGen->GetYaxis()->SetTitleOffset(1.2);

		BsHisDecayCompGen->Scale(1.0/BsHisDecayCompGen->Integral());
		BsHisDecayCompGen->SetLineWidth(2);
		BsHisDecayCompGen->SetLineColor(kBlue);
		


		BsHisDecayComp->GetXaxis()->SetTitle(XAxisName[i].Data());
		BsHisDecayComp->GetYaxis()->SetTitle("Normalized Counts");
		BsHisDecayComp->GetXaxis()->CenterTitle();
		BsHisDecayComp->GetYaxis()->CenterTitle();
		BsHisDecayComp->GetXaxis()->SetTitleOffset(1.2);
		BsHisDecayComp->GetYaxis()->SetTitleOffset(1.2);

		BsHisDecayComp->Scale(1.0/BsHisDecayComp->Integral());
		BsHisDecayComp->SetLineWidth(2);
		BsHisDecayComp->SetLineColor(kRed);



		legDecayComp[i] = new TLegend(0.36,0.65,0.75,0.85,NULL,"brNDC");
		legDecayComp[i]->SetBorderSize(0);
		legDecayComp[i]->SetTextSize(0.040);
		legDecayComp[i]->SetTextFont(42);
		legDecayComp[i]->SetFillStyle(0);
		legDecayComp[i]->SetLineWidth(3);

		legDecayComp[i]->AddEntry(BsHisDecayCompGen,"GEN: Decay","l");
		legDecayComp[i]->AddEntry(BsHisDecayComp,"RECO: Decay","l");

		BsHisDecayCompGen->Draw("hist");
		BsHisDecayComp->Draw("histSAME");

		legDecayComp[i]->Draw("SAME");

		OutName = Form("PlotsComparison/Decay/%s.png",SaveName[i].Data());

		c->SaveAs(OutName.Data());



	}




}
