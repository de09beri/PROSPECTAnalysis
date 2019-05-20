#include "PROSPECT_Style.cc"
#include "TROOT.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TPave.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TColor.h"
#include "TExec.h"
#include <sstream>
#include "TLatex.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "Header.C"


int PlotCutStudyResults(){

	setup_PROSPECT_style();
  	gROOT->ForceStyle();

	TFile *f = new TFile(Form("%s/CutStudiesResults.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));
	TGraph *grPromptE = (TGraph*)f->Get("PromptE");
	TGraph *grDelayE = (TGraph*)f->Get("DelayE");
	TGraph *grPromptPSD = (TGraph*)f->Get("PromptPSD");
	TGraph *grDelayPSD = (TGraph*)f->Get("DelayPSD");

	grPromptE->SetMarkerStyle(24);
	grDelayE->SetMarkerStyle(24);
	grPromptPSD->SetMarkerStyle(25);
	grDelayPSD->SetMarkerStyle(25);
	grPromptE->SetMarkerColor(1);
	grDelayE->SetMarkerColor(2);
	grPromptPSD->SetMarkerColor(4);
	grDelayPSD->SetMarkerColor(8);
	grPromptE->SetMarkerSize(1);
	grDelayE->SetMarkerSize(1);
	grPromptPSD->SetMarkerSize(1);
	grDelayPSD->SetMarkerSize(1);
	grPromptE->SetLineColor(1);
	grDelayE->SetLineColor(2);
	grPromptPSD->SetLineColor(4);
	grDelayPSD->SetLineColor(8);

	TLine *l = new TLine(2.35,0,4.15,0);
	l->SetLineStyle(2);
	l->SetLineColor(14);

	TCanvas *c = new TCanvas("c","c",1);
	grPromptE->Draw("APL");
	grPromptE->GetYaxis()->SetRangeUser(-40e-6,40e-6);
	grDelayE->Draw("PL");
	grPromptPSD->Draw("PL");
	grDelayPSD->Draw("PL");
	l->Draw("same");
	TLegend *leg = new TLegend(0.8,0.8,1,1);
	leg->AddEntry(grPromptE,"Prompt Energy","p");
	leg->AddEntry(grDelayE,"Delay Energy","p");
	leg->AddEntry(grPromptPSD,"Prompt PSD","p");
	leg->AddEntry(grDelayPSD,"Delay PSD","p");
	leg->Draw();	
	c->SaveAs(Form("%s/CutStudyResults.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	return 0;
}
