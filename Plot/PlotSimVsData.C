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
#include "TRandom.h"

int PlotSimVsData(){

	setup_PROSPECT_style();
  	gROOT->ForceStyle();

	int timeBin = 0;
	//-------------------------------------------------------------------------------------------------------
	TFile *fData = new TFile(Form("%s/Ac227_HistsPerTime_SimComparison.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));
	if(!fData){
		printf("File not found. Exiting. \n");
		return -1;
	}
cout<<"Loading hists"<<endl;

	TH1F *hRnPoDt_Data = (TH1F*)fData->Get(Form("hRnPoDt_%i",timeBin));
	TH1F *hRnPSD_Data  = (TH1F*)fData->Get(Form("hRnPSD_%i",timeBin));
	TH1F *hPoPSD_Data  = (TH1F*)fData->Get(Form("hPoPSD_%i",timeBin));
	TH1F *hRnEn_Data   = (TH1F*)fData->Get(Form("hRnEn_%i",timeBin));
	TH1F *hPoEn_Data   = (TH1F*)fData->Get(Form("hPoEn_%i",timeBin));
	TH1F *hRnPoDz_Data = (TH1F*)fData->Get(Form("hRnPoDz_%i",timeBin));

	//-------------------------------------------------------------------------------------------------------
	TFile *fSim = new TFile("/g/g20/berish1/Simulation/FinalRun_Results/Ac227_HistsPerTime.root");
	if(!fSim){
		printf("File not found. Exiting. \n");
		return -1;
	}
cout<<"Loading hists"<<endl;

	TH1F *hRnPoDt_Sim = (TH1F*)fSim->Get(Form("hRnPoDt_%i",timeBin));
	TH1F *hRnPSD_Sim  = (TH1F*)fSim->Get(Form("hRnPSD_%i",timeBin));
	TH1F *hPoPSD_Sim  = (TH1F*)fSim->Get(Form("hPoPSD_%i",timeBin));
	TH1F *hRnEn_Sim   = (TH1F*)fSim->Get(Form("hRnEn_%i",timeBin));
	TH1F *hPoEn_Sim   = (TH1F*)fSim->Get(Form("hPoEn_%i",timeBin));
	TH1F *hRnPoDz_Sim = (TH1F*)fSim->Get(Form("hRnPoDz_%i",timeBin));

	//-------------------------------------------------------------------------------------------------------
	int numBins = hPoEn_Sim->GetNbinsX();	

	//-------------------------------------------------------------------------------------------------------
	hRnEn_Data->Scale(1/hRnEn_Data->Integral());
	hRnEn_Sim->Scale(1/hRnEn_Sim->Integral());
	hPoEn_Data->Scale(1/hPoEn_Data->Integral());
	hPoEn_Sim->Scale(1/hPoEn_Sim->Integral());
	
	TH1F *hRnEn_Sub = (TH1F*)hRnEn_Data->Clone();
	hRnEn_Sub->Sumw2();
	hRnEn_Sub->Add(hRnEn_Sim,-1);

	TH1F *hPoEn_Sub = (TH1F*)hPoEn_Data->Clone();
	hPoEn_Sub->Sumw2();
	hPoEn_Sub->Add(hPoEn_Sim,-1);

	hRnEn_Data->SetLineColor(kRed);
	hRnEn_Data->SetMarkerColor(kRed);
	hRnEn_Sub->SetLineColor(kBlack);
	hRnEn_Sub->SetMarkerColor(kBlack);
	hPoEn_Data->SetLineColor(kRed);
	hPoEn_Data->SetMarkerColor(kRed);

	//-------------------------------------------------------------------------------------------------------
	TCanvas *cRnEn = new TCanvas("cRnEn","RnEn",700,700);
	TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
        TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
        pad1->Draw();
        pad2->Draw();

	pad1->cd();
	hRnEn_Data->GetYaxis()->SetRangeUser(0,0.05);
	hRnEn_Data->GetYaxis()->SetTitle("Normalized Counts");
	hRnEn_Data->Draw();
	hRnEn_Sim->Draw("same");
	TLegend *leg = new TLegend(0.82,0.82,0.97,0.95);
	leg->AddEntry(hRnEn_Data,"Data");
	leg->AddEntry(hRnEn_Sim,"Simulation");
	leg->Draw();
	pad2->cd();
	gPad->SetGrid();
	hRnEn_Sub->GetYaxis()->SetRangeUser(-0.005,0.005);
	hRnEn_Sub->GetYaxis()->SetTitle("Data-Sim");
	hRnEn_Sub->GetYaxis()->SetTitleSize(0.07);
	hRnEn_Sub->GetYaxis()->SetLabelSize(0.07);
	hRnEn_Sub->GetYaxis()->SetTitleOffset(0.4);
	hRnEn_Sub->GetXaxis()->SetTitleSize(0.07);
	hRnEn_Sub->GetXaxis()->SetLabelSize(0.07);
	hRnEn_Sub->GetXaxis()->SetTitleOffset(0.5);
	hRnEn_Sub->Draw();
	hRnEn_Sub->Fit("pol0");
	TPaveText *pt = new TPaveText(0.75,0.75,0.99,0.99,"NDCNB");
        pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hRnEn_Sub->GetFunction("pol0")->GetChisquare(),hRnEn_Sub->GetFunction("pol0")->GetNDF()));
        pt->AddText(Form("Prob   %.2f",hRnEn_Sub->GetFunction("pol0")->GetProb()));
        pt->AddText(Form("p0   %f #pm %f",hRnEn_Sub->GetFunction("pol0")->GetParameter(0),hRnEn_Sub->GetFunction("pol0")->GetParError(0)));
        pt->Draw();
	cRnEn->SaveAs(Form("%s/DataVsSim_RnEn.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoEn = new TCanvas("cPoEn","PoEn",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
        pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
        pad1->Draw();
        pad2->Draw();

	pad1->cd();
	hPoEn_Data->GetYaxis()->SetRangeUser(0,0.055);
	hPoEn_Data->GetYaxis()->SetTitle("Normalized Counts");
	hPoEn_Data->Draw();
	hPoEn_Sim->Draw("same");
	leg = new TLegend(0.82,0.82,0.97,0.95);
	leg->AddEntry(hPoEn_Data,"Data");
	leg->AddEntry(hPoEn_Sim,"Simulation");
	leg->Draw();
	pad2->cd();
	gPad->SetGrid();
	hPoEn_Sub->GetYaxis()->SetRangeUser(-0.0055,0.0055);
	hPoEn_Sub->GetYaxis()->SetTitle("Data-Sim");
	hPoEn_Sub->GetYaxis()->SetTitleSize(0.07);
	hPoEn_Sub->GetYaxis()->SetLabelSize(0.07);
	hPoEn_Sub->GetYaxis()->SetTitleOffset(0.4);
	hPoEn_Sub->GetXaxis()->SetTitleSize(0.07);
	hPoEn_Sub->GetXaxis()->SetLabelSize(0.07);
	hPoEn_Sub->GetXaxis()->SetTitleOffset(0.5);
	hPoEn_Sub->Draw();
	hPoEn_Sub->Fit("pol0");
	pt = new TPaveText(0.75,0.75,0.99,0.99,"NDCNB");
        pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hPoEn_Sub->GetFunction("pol0")->GetChisquare(),hPoEn_Sub->GetFunction("pol0")->GetNDF()));
        pt->AddText(Form("Prob   %.2f",hPoEn_Sub->GetFunction("pol0")->GetProb()));
        pt->AddText(Form("p0   %f #pm %f",hPoEn_Sub->GetFunction("pol0")->GetParameter(0),hPoEn_Sub->GetFunction("pol0")->GetParError(0)));
        pt->Draw();
	cPoEn->SaveAs(Form("%s/DataVsSim_PoEn.pdf",gSystem->Getenv("AD_AC227_PLOTS")));


	return 0;
} 	//end PlotDistributionsVsTime
