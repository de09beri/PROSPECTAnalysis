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

TH1F *makeResidHist(TH1F *h, TF1 *fit){
	TH1F *hResid = (TH1F*)h->Clone();
	int numPt = h->GetSize() - 2;

	double hx, hy, hyErr;
	double fity;
	double resid, residErr;

	for(int i=1;i<numPt+1;i++){
		hx = h->GetXaxis()->GetBinCenter(i);
		hy = h->GetBinContent(i);	
		hyErr = h->GetBinError(i);

		fity = fit->Eval(hx);
		
		resid = (hy-fity);
		residErr = hyErr;

		hResid->SetBinContent(i,resid);
		hResid->SetBinError(i,residErr);
	}

	hResid->GetYaxis()->SetTitle("Data - fit");
	hResid->GetYaxis()->SetLabelSize(0.07);
	hResid->GetYaxis()->SetTitleSize(0.09);
	hResid->GetYaxis()->SetTitleOffset(0.4);
	hResid->GetXaxis()->SetLabelSize(0.07);
	hResid->GetXaxis()->SetTitleSize(0.06);
	hResid->GetXaxis()->SetTitleOffset(0.7);


	return hResid;
}

int PlotDistributionsVsTime(int timeBin){

	int grIDX = timeBin;

	setup_PROSPECT_style();
  	gROOT->ForceStyle();

	TFile *fg = new TFile(Form("%s/Ac227_GraphsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"))); 
	if(!fg){
		printf("Graph file not found. Exiting. \n");
		return -1;
	}
	
	TGraphErrors *grRate      = (TGraphErrors*)fg->Get("grRate");
	TGraphErrors *grRnEnEff   = (TGraphErrors*)fg->Get("grPromptEnEff");
	TGraphErrors *grPoEnEff   = (TGraphErrors*)fg->Get("grDelayEnEff");
	TGraphErrors *grRnPSDEff  = (TGraphErrors*)fg->Get("grPromptPSDEff");
	TGraphErrors *grPoPSDEff  = (TGraphErrors*)fg->Get("grDelayPSDEff");
	TGraphErrors *grRnPoDzEff = (TGraphErrors*)fg->Get("grDzEff");

	fg->Close();

	//-------------------------------------------------------------------------------------------------------
	TFile *f = new TFile(Form("%s/Ac227_HistsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));
	if(!f){
		printf("File not found. Exiting. \n");
		return -1;
	}
cout<<"Loading hists"<<endl;
	TH1F *hSelectDt 	= (TH1F*)f->Get(Form("hSelectDt_%i",timeBin));
	TH1F *hBGDt 		= (TH1F*)f->Get(Form("hBGDt_%i",timeBin));
	TH1F *hRnPoDt 	 	= (TH1F*)f->Get(Form("hRnPoDt_%i",timeBin));
	TF1  *fRnPoDtExp 	= hRnPoDt->GetFunction("fRnPoDtExp");
	fRnPoDtExp->SetRange(0,2.57*5.0);

	TH1F *hSelectPromptPSD  = (TH1F*)f->Get(Form("hSelectPromptPSD_%i",timeBin));
	TH1F *hBGPromptPSD  	= (TH1F*)f->Get(Form("hBGPromptPSD_%i",timeBin));
	TH1F *hRnPSD 	 	= (TH1F*)f->Get(Form("hRnPSD_%i",timeBin));
	TF1  *fRnPSDGaus 	= hRnPSD->GetFunction("fRnPSDGaus");

	TH1F *hSelectDelayPSD   = (TH1F*)f->Get(Form("hSelectDelayPSD_%i",timeBin));
	TH1F *hBGDelayPSD   	= (TH1F*)f->Get(Form("hBGDelayPSD_%i",timeBin));
	TH1F *hPoPSD 	 	= (TH1F*)f->Get(Form("hPoPSD_%i",timeBin));
	TF1  *fPoPSDGaus 	= hPoPSD->GetFunction("fPoPSDGaus");

	TH1F *hSelectPromptEn 	= (TH1F*)f->Get(Form("hSelectPromptEn_%i",timeBin));
	TH1F *hBGPromptEn 	= (TH1F*)f->Get(Form("hBGPromptEn_%i",timeBin));
	TH1F *hRnEn   		= (TH1F*)f->Get(Form("hRnEn_%i",timeBin));
	TF1  *fRnEnGaus 	= hRnEn->GetFunction("fRnEnGaus");

	TH1F *hSelectDelayEn 	= (TH1F*)f->Get(Form("hSelectDelayEn_%i",timeBin));
	TH1F *hBGDelayEn 	= (TH1F*)f->Get(Form("hBGDelayEn_%i",timeBin));
	TH1F *hPoEn 		= (TH1F*)f->Get(Form("hPoEn_%i",timeBin));
	TF1  *fPoEnGaus 	= hPoEn->GetFunction("fPoEnGaus");

	TH1F *hSelectPromptPos	= (TH1F*)f->Get(Form("hSelectPromptPos_%i",timeBin));
	TH1F *hBGPromptPos	= (TH1F*)f->Get(Form("hBGPromptPos_%i",timeBin));
	TH1F *hRnPos 		= (TH1F*)f->Get(Form("hRnPos_%i",timeBin));

	TH1F *hSelectDelayPos	= (TH1F*)f->Get(Form("hSelectDelayPos_%i",timeBin));
	TH1F *hBGDelayPos	= (TH1F*)f->Get(Form("hBGDelayPos_%i",timeBin));
	TH1F *hPoPos 		= (TH1F*)f->Get(Form("hPoPos_%i",timeBin));

	TH1F *hSelectDz		= (TH1F*)f->Get(Form("hSelectDz_%i",timeBin));
	TH1F *hBGDz		= (TH1F*)f->Get(Form("hBGDz_%i",timeBin));
	TH1F *hRnPoDz 	  	= (TH1F*)f->Get(Form("hRnPoDz_%i",timeBin));
	TF1  *fRnPoDzGaus 	= hRnPoDz->GetFunction("fRnPoDzGaus");

	TH2F *hRnPoPSDvsEn 	= (TH2F*)f->Get(Form("hRnPoPSDvsEn_%i",timeBin));

	TH2F *hPoEnvsRnEn 	= (TH2F*)f->Get(Form("hPoEnvsRnEn_%i",timeBin));

	//-------------------------------------------------------------------------------------------------------
	TPaveText *pt;
	int p_col = 1, d_col = 4;
	int s_col = 8, bg_col = 1;
	
	TCanvas *cRnPoDt = new TCanvas("cRnPoDt","RnPo Dt",700,700);
	TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();
	gPad->SetLogy();
	hSelectDt->GetYaxis()->SetTitleOffset(1.1);
	hSelectDt->SetLineColor(8);
	hSelectDt->SetLineWidth(1);
	hSelectDt->SetMinimum(1e2);
	hSelectDt->Draw("HIST");
	hBGDt->SetLineColor(kBlack);
	hBGDt->SetLineWidth(1);
	hBGDt->Draw("same&HIST");
	hRnPoDt->SetLineWidth(1);
	hRnPoDt->Draw("same");
	hRnPoDt->GetXaxis()->SetTitle("t_{#alphaPo} - t_{#alphaRn} [ms]");
	fRnPoDtExp->SetLineStyle(2);
	fRnPoDtExp->SetRange(0,12.85);
	fRnPoDtExp->Draw("same");
	pt = new TPaveText(0.72,0.75,0.99,0.99,"NDCNB");
	pt->AddText(Form("Entries    %.0f",hRnPoDt->GetEntries()));
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPoDtExp->GetChisquare(),fRnPoDtExp->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnPoDtExp->GetProb()));
	pt->AddText(Form("N    %.2f #pm %.2f",fRnPoDtExp->GetParameter(0),fRnPoDtExp->GetParError(0)));
	pt->AddText(Form("#tau_{Po^{215}}    %.3f #pm %.3f ms",fRnPoDtExp->GetParameter(1),fRnPoDtExp->GetParError(1)));
	pt->Draw();

	pad2->cd();
	gPad->SetGrid();
	TH1F *hRnPoDtResid = makeResidHist(hRnPoDt,fRnPoDtExp);
	hRnPoDtResid->GetXaxis()->SetTitle("dt [ms]");
	hRnPoDtResid->GetYaxis()->SetRangeUser(-700,700);
	hRnPoDtResid->Draw();
	hRnPoDtResid->Fit("pol0","","",0.5,12.85);
	hRnPoDtResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hRnPoDtResid->GetFunction("pol0")->GetChisquare(),hRnPoDtResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hRnPoDtResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hRnPoDtResid->GetFunction("pol0")->GetParameter(0),hRnPoDtResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	

	cRnPoDt->SaveAs(Form("%s/RnPoDt_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));

	//----------------------------------------------------------
	TCanvas *cRnPSD = new TCanvas("cRnPSD","Rn PSD",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();
	hSelectPromptPSD->GetYaxis()->SetTitleOffset(1.1);
	hSelectPromptPSD->SetLineColor(s_col);
	hSelectPromptPSD->Draw("HIST");
	hBGPromptPSD->SetLineColor(bg_col);
	hBGPromptPSD->Draw("HIST&SAME");
	hRnPSD->Draw("same");
	fRnPSDGaus->SetLineStyle(2);
	fRnPSDGaus->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
	TText *tRn = pt->AddText("Rn^{219}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPSDGaus->GetChisquare(),fRnPSDGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnPSDGaus->GetProb()));
	pt->AddText(Form("Eff    %.3f #pm %.3f",grRnPSDEff->GetY()[grIDX],grRnPSDEff->GetEY()[grIDX]));
	pt->Draw();

	pad2->cd();
	gPad->SetGrid();
	TH1F *hRnPSDResid = makeResidHist(hRnPSD,fRnPSDGaus);
//	hRnPSDResid->GetYaxis()->SetRangeUser(-70,70);
	hRnPSDResid->Draw();
	hRnPSDResid->Fit("pol0");
	hRnPSDResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hRnPSDResid->GetFunction("pol0")->GetChisquare(),hRnPSDResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hRnPSDResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hRnPSDResid->GetFunction("pol0")->GetParameter(0),hRnPSDResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	

	cRnPSD->SaveAs(Form("%s/RnPSD_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));

	//----------------------------------------------------------
	TCanvas *cPoPSD = new TCanvas("cPoPSD","Po PSD",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();
	hSelectDelayPSD->GetYaxis()->SetTitleOffset(1.1);
	hSelectDelayPSD->SetLineColor(s_col);
	hSelectDelayPSD->Draw("HIST");
	hBGDelayPSD->SetLineColor(bg_col);
	hBGDelayPSD->Draw("HIST&SAME");
	hPoPSD->Draw("same");
	fPoPSDGaus->SetLineStyle(2);
	fPoPSDGaus->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
	TText *tPo = pt->AddText("Po^{215}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fPoPSDGaus->GetChisquare(),fPoPSDGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fPoPSDGaus->GetProb()));
	pt->AddText(Form("Eff    %.3f #pm %.3f",grPoPSDEff->GetY()[grIDX],grPoPSDEff->GetEY()[grIDX]));
	pt->Draw();

	pad2->cd();
	gPad->SetGrid();
	TH1F *hPoPSDResid = makeResidHist(hPoPSD,fPoPSDGaus);
//	hPoPSDResid->GetYaxis()->SetRangeUser(-70,70);
	hPoPSDResid->Draw();
	hPoPSDResid->Fit("pol0");
	hPoPSDResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hPoPSDResid->GetFunction("pol0")->GetChisquare(),hPoPSDResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hPoPSDResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hPoPSDResid->GetFunction("pol0")->GetParameter(0),hPoPSDResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	

	cPoPSD->SaveAs(Form("%s/PoPSD_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));

	//----------------------------------------------------------
	TCanvas *cRnEn = new TCanvas("cRnEn","Rn En",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();

	hSelectPromptEn->GetYaxis()->SetTitleOffset(1.1);
	hSelectPromptEn->SetLineColor(s_col);
	hSelectPromptEn->Draw("HIST");
	hBGPromptEn->SetLineColor(bg_col);
	hBGPromptEn->Draw("HIST&SAME");
	hRnEn->Draw("same");
	fRnEnGaus->SetLineStyle(2);	
	fRnEnGaus->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
	tRn = pt->AddText("Rn^{219}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnEnGaus->GetChisquare(),fRnEnGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnEnGaus->GetProb()));
	pt->AddText(Form("Eff    %.3f #pm %.3f",grRnEnEff->GetY()[grIDX],grRnEnEff->GetEY()[grIDX]));
	pt->Draw();

	pad2->cd();
	gPad->SetGrid();
	TH1F *hRnEnResid = makeResidHist(hRnEn,fRnEnGaus);
//	hRnEnResid->GetYaxis()->SetRangeUser(-70,70);
	hRnEnResid->Draw();
	hRnEnResid->Fit("pol0");
	hRnEnResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hRnEnResid->GetFunction("pol0")->GetChisquare(),hRnEnResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hRnEnResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hRnEnResid->GetFunction("pol0")->GetParameter(0),hRnEnResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	

	cRnEn->SaveAs(Form("%s/RnEn_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));

	//----------------------------------------------------------
	TCanvas *cPoEn = new TCanvas("cPoEn","Po En",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();
	hSelectDelayEn->GetYaxis()->SetTitleOffset(1.1);
	hSelectDelayEn->SetLineColor(s_col);
	hSelectDelayEn->Draw("HIST");
	hBGDelayEn->SetLineColor(bg_col);
	hBGDelayEn->Draw("HIST&SAME");
	hPoEn->Draw("same");
	fPoEnGaus->SetLineStyle(2);	
	fPoEnGaus->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
	tPo = pt->AddText("Po^{215}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fPoEnGaus->GetChisquare(),fPoEnGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fPoEnGaus->GetProb()));
	pt->AddText(Form("Eff    %.3f #pm %.3f",grPoEnEff->GetY()[grIDX],grPoEnEff->GetEY()[grIDX]));
	pt->Draw();

	pad2->cd();
	gPad->SetGrid();
	TH1F *hPoEnResid = makeResidHist(hPoEn,fPoEnGaus);
//	hPoEnResid->GetYaxis()->SetRangeUser(-70,70);
	hPoEnResid->Draw();
	hPoEnResid->Fit("pol0");
	hPoEnResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hPoEnResid->GetFunction("pol0")->GetChisquare(),hPoEnResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hPoEnResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hPoEnResid->GetFunction("pol0")->GetParameter(0),hPoEnResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	

	cPoEn->SaveAs(Form("%s/PoEn_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));

	//----------------------------------------------------------
	TCanvas *cRnPoDz = new TCanvas("cRnPoDz","RnPo Dz",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();
	hSelectDz->GetYaxis()->SetTitleOffset(1.1);
	hSelectDz->SetLineColor(s_col);
	hSelectDz->Draw("HIST");
	hBGDz->SetLineColor(bg_col);
	hBGDz->Draw("HIST&SAME");
	hRnPoDz->GetXaxis()->SetTitle("z_{#alphaPo} - z_{#alphaRn} [mm]");
	hRnPoDz->Draw("same");
	fRnPoDzGaus->SetLineStyle(2);
	fRnPoDzGaus->Draw("same");	
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPoDzGaus->GetChisquare(),fRnPoDzGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnPoDzGaus->GetProb()));
	pt->AddText(Form("Eff    %.2f #pm %.2f",grRnPoDzEff->GetY()[grIDX],grRnPoDzEff->GetEY()[grIDX]));
	pt->Draw();
	
	pad2->cd();
	gPad->SetGrid();
	TH1F *hRnPoDzResid = makeResidHist(hRnPoDz,fRnPoDzGaus);
//	hRnPoDzResid->GetYaxis()->SetRangeUser(-70,70);
	hRnPoDzResid->Draw();
	hRnPoDzResid->Fit("pol0");
	hRnPoDzResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hRnPoDzResid->GetFunction("pol0")->GetChisquare(),hRnPoDzResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hRnPoDzResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hRnPoDzResid->GetFunction("pol0")->GetParameter(0),hRnPoDzResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	
	cRnPoDz->SaveAs(Form("%s/RnPoDz_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));

	//----------------------------------------------------------
	TCanvas *cRnPos = new TCanvas("cRnPos","Rn Position",1);
	gPad->SetGrid();
	hSelectPromptPos->SetLineColor(s_col);
	hSelectPromptPos->Draw("HIST");
	hBGPromptPos->SetLineColor(bg_col);
	hBGPromptPos->Draw("HIST&SAME");
	hPoPos->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
    	tRn = pt->AddText("Rn^{219}");
  	pt->AddText(Form("#mu    %.2f #pm %.2f",hRnPos->GetMean(),hRnPos->GetMeanError()));
   	pt->AddText(Form("RMS    %.2f #pm %.2f",hRnPos->GetRMS(),hRnPos->GetRMSError()));
	pt->Draw();
	cRnPos->SaveAs(Form("%s/RnPos_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));

	TCanvas *cPoPos = new TCanvas("cPoPos","Po Position",1);
	gPad->SetGrid();
	hSelectDelayPos->SetLineColor(s_col);
	hSelectDelayPos->Draw("HIST");
	hBGDelayPos->SetLineColor(bg_col);
	hBGDelayPos->Draw("HIST&SAME");
	hPoPos->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
    	tPo = pt->AddText("Po^{215}");
  	pt->AddText(Form("#mu    %.2f #pm %.2f",hPoPos->GetMean(),hPoPos->GetMeanError()));
   	pt->AddText(Form("RMS    %.2f #pm %.2f",hPoPos->GetRMS(),hPoPos->GetRMSError()));
	pt->Draw();
	cPoPos->SaveAs(Form("%s/PoPos_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
	
	TCanvas *cRnPoPSDvsEn = new TCanvas("cRnPoPSDvsEn","RnPo PSD vs En",1);
	gPad->SetRightMargin(0.1);
	gPad->SetLeftMargin(0.12);
	hRnPoPSDvsEn->GetYaxis()->SetTitleOffset(1.1);
	hRnPoPSDvsEn->Draw("colz");
	cRnPoPSDvsEn->SaveAs(Form("%s/RnPoPSDvsEn_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
	
	TCanvas *cPoEnvsRnEn = new TCanvas("cPoEnvsRnEn","Po En vs Rn En",650,600);
	gPad->SetRightMargin(0.12);
	gPad->SetLeftMargin(0.12);
	hPoEnvsRnEn->GetYaxis()->SetTitleOffset(1.1);
	hPoEnvsRnEn->Draw("colz");
	cPoEnvsRnEn->SaveAs(Form("%s/PoEnvsRnEn_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));


	return 0;
} 	//end PlotDistributionsVsTime
