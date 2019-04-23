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

TGraphErrors *makeETGr(TGraphErrors *gr){
	const int numET = 46;
	int ET[numET] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,27,28,41,42,55,56,69,70,83,84,97,98,111,112,125,126,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153}; 	

	int numPt = gr->GetN();

	vector<double> vGrx, vGry, vGryErr;
	double grx, gry, gryErr;
	bool boolET;
	
	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);	

		boolET = find(begin(ET), end(ET), grx) != end(ET);
                if(boolET){
			vGrx.push_back(grx);
			vGry.push_back(gry);
			vGryErr.push_back(gryErr);
		}
	}

	int size = vGrx.size();
	TGraphErrors *grET = new TGraphErrors(size); 

	for(int i=0;i<size;i++){
		grET->SetPoint(i,vGrx[i],vGry[i]);
		grET->SetPointError(i,0,vGryErr[i]);
	}
	
	return grET;
}


TGraphErrors *makeRelGr(TGraphErrors *gr){
	gr->Fit("pol0","Q0");
        double mean = gr->GetFunction("pol0")->GetParameter(0);
        double meanErr = gr->GetFunction("pol0")->GetParError(0);

	TGraphErrors *grRel = (TGraphErrors*)gr->Clone();
	int numPt = gr->GetN();

	double rel, relErr; 
	double grx, gry, gryErr;

	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);

		rel = gry/mean;
		relErr = rel * sqrt(pow(gryErr/gry,2) + pow(meanErr/mean,2));	

		grRel->SetPoint(i,grx,rel);
		grRel->SetPointError(i,0,relErr);
	}

	return grRel;
}


void PlotRnPoVsCell(){

	setup_PROSPECT_style();
    	gROOT->ForceStyle();

	TFile *f = new TFile(Form("%s/Ac227_GraphsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));

	TGraphErrors *grRate 		= (TGraphErrors*)f->Get("grRate");
	TGraphErrors *grTotEff 		= (TGraphErrors*)f->Get("grTotEff");
	TGraphErrors *grLifetime 	= (TGraphErrors*)f->Get("grLifetime");

	TGraphErrors *grPoPSDMean 	= (TGraphErrors*)f->Get("grPoPSDMean");
	TGraphErrors *grPoPSDSigma 	= (TGraphErrors*)f->Get("grPoPSDSigma");
	TGraphErrors *grPoEnMean 	= (TGraphErrors*)f->Get("grPoEnMean");
	TGraphErrors *grPoEnSigma 	= (TGraphErrors*)f->Get("grPoEnSigma");
	TGraphErrors *grPoEnSmearMean 	= (TGraphErrors*)f->Get("grPoEnSmearMean");
	TGraphErrors *grPoEnSmearSigma 	= (TGraphErrors*)f->Get("grPoEnSmearSigma");
	TGraphErrors *grPoPosMean 	= (TGraphErrors*)f->Get("grPoPosMean");
	TGraphErrors *grPoPosSigma 	= (TGraphErrors*)f->Get("grPoPosSigma");
	TGraphErrors *grRnPoDzMean 	= (TGraphErrors*)f->Get("grRnPoDzMean");
	TGraphErrors *grRnPoDzSigma 	= (TGraphErrors*)f->Get("grRnPoDzSigma");

	TGraphErrors *grRnEnEff 	= (TGraphErrors*)f->Get("grPromptEnEff");
	TGraphErrors *grPoEnEff 	= (TGraphErrors*)f->Get("grDelayEnEff");
	TGraphErrors *grRnPSDEff 	= (TGraphErrors*)f->Get("grPromptPSDEff");
	TGraphErrors *grPoPSDEff 	= (TGraphErrors*)f->Get("grDelayPSDEff");
	TGraphErrors *grRnPoDzEff 	= (TGraphErrors*)f->Get("grDzEff");

	TGraph *grRnPSDChiSq = (TGraph*)f->Get("grRnPSDChiSq");
	TGraph *grPoPSDChiSq = (TGraph*)f->Get("grPoPSDChiSq");
	TGraph *grRnEnChiSq  = (TGraph*)f->Get("grRnEnChiSq");
	TGraph *grPoEnChiSq  = (TGraph*)f->Get("grPoEnChiSq");
	TGraph *grDzChiSq    = (TGraph*)f->Get("grDzChiSq");
	TGraph *grDtChiSq    = (TGraph*)f->Get("grDtChiSq");

	TGraph *grBGRate = (TGraph*)f->Get("grBGRate");

	TH1F *hRate = (TH1F*)f->Get("hRate");

//	f->Close();
	
	//-------------------------------------------------------------------------------------------------------
	double hRateMin = grRate->GetHistogram()->GetMinimum()-0.02;
	double hRateMax = grRate->GetHistogram()->GetMaximum()+0.02;

	TH2F *h2DRate = new TH2F("h2DRate","Rate per cell",14,0,14,11,0,11);

	int numPt = grRate->GetN();
	double grx,gry, gryErr;
	int binx,biny;
	for(int i=0;i<numPt;i++){
		grRate->GetPoint(i,grx,gry);
		gryErr = grRate->GetErrorY(i);
		binx = (int)grx%14 + 1;
		biny = ((int)grx/14) + 1;

		h2DRate->SetBinContent(binx,biny,gry);
	}

	//-------------------------------------------------------------------------------------------------------
	TGraphErrors *grRateET = makeETGr(grRate); 
	grRateET->SetMarkerColor(kRed);
	grRateET->SetLineColor(kRed);

	TGraphErrors *grPoPSDMeanET = makeETGr(grPoPSDMean);
	grPoPSDMeanET->SetMarkerColor(kRed);
	grPoPSDMeanET->SetLineColor(kRed);
	
	TGraphErrors *grPoPSDSigmaET = makeETGr(grPoPSDSigma);
	grPoPSDSigmaET->SetMarkerColor(kRed);
	grPoPSDSigmaET->SetLineColor(kRed);

	TGraphErrors *grPoEnMeanET = makeETGr(grPoEnMean);
	grPoEnMeanET->SetMarkerColor(kRed);
	grPoEnMeanET->SetLineColor(kRed);
	
	TGraphErrors *grPoEnSigmaET = makeETGr(grPoEnSigma);
	grPoEnSigmaET->SetMarkerColor(kRed);
	grPoEnSigmaET->SetLineColor(kRed);

	TGraphErrors *grPoPosMeanET = makeETGr(grPoPosMean);
	grPoPosMeanET->SetMarkerColor(kRed);
	grPoPosMeanET->SetLineColor(kRed);
	
	TGraphErrors *grPoPosSigmaET = makeETGr(grPoPosSigma);
	grPoPosSigmaET->SetMarkerColor(kRed);
	grPoPosSigmaET->SetLineColor(kRed);

	TGraphErrors *grRnPoDzMeanET = makeETGr(grRnPoDzMean);
	grRnPoDzMeanET->SetMarkerColor(kRed);
	grRnPoDzMeanET->SetLineColor(kRed);

	TGraphErrors *grRnPoDzSigmaET = makeETGr(grRnPoDzSigma);
	grRnPoDzSigmaET->SetMarkerColor(kRed);
	grRnPoDzSigmaET->SetLineColor(kRed);

	//-------------------------------------------------------------------------------------------------------
	const char *xLabel = "Cell";

	TCanvas *cRate = new TCanvas("cRate","Rate per cell",1000,400);
	grRate->GetXaxis()->SetTitle(xLabel);
	grRate->GetYaxis()->SetTitle("R_{RnPo} [mHz]");
	grRate->Draw("AP");
	grRateET->Draw("P");
	cRate->SaveAs(Form("%s/RatePerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *chRate = new TCanvas("chRate","Rate per cell",1);
	hRate->GetXaxis()->SetTitle("R_{RnPo} [mHz]");
	hRate->GetYaxis()->SetTitle("Counts");
	hRate->GetXaxis()->SetTitleOffset(1.35);
	hRate->GetXaxis()->SetTitleSize(0.04);
	hRate->Draw("HIST");
	chRate->SaveAs(Form("%s/HistRatePerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	gStyle->SetPaintTextFormat("2.2f");
	TCanvas *cRate2D = new TCanvas("cRate2D","Rate per cell",1);
	gPad->SetRightMargin(0.16);
	h2DRate->SetMarkerSize(1.25);
	h2DRate->SetMinimum(grRate->GetHistogram()->GetMinimum());
	h2DRate->Draw("colz && text");	
	cRate2D->SaveAs(Form("%s/2DRatePerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cBGRate = new TCanvas("cBGRate","BGRate per cell",1000,400);
	gPad->SetLogy();
	gPad->SetGrid();
	grBGRate->GetXaxis()->SetTitle(xLabel);
	grBGRate->GetYaxis()->SetTitle("BG Rate [mHz]");
	grBGRate->Draw("AP");
	cBGRate->SaveAs(Form("%s/BGRatePerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-------------------------
	TCanvas *cPoPSDMean = new TCanvas("cPoPSDMean","Po PSD Mean",1000,400);
	grPoPSDMean->GetXaxis()->SetTitle(xLabel);
	grPoPSDMean->GetYaxis()->SetTitle("PSD");
	grPoPSDMean->Draw("AP");
	grPoPSDMeanET->Draw("P");
	cPoPSDMean->SaveAs(Form("%s/PoPSDMeanPerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoPSDSigma = new TCanvas("cPoPSDSigma"," Po PSD Sigma",1000,400);
	grPoPSDSigma->GetXaxis()->SetTitle(xLabel);
	grPoPSDSigma->GetYaxis()->SetTitle("#sigma_{PSD}");
	grPoPSDSigma->Draw("AP");
	grPoPSDSigmaET->Draw("P");
	cPoPSDSigma->SaveAs(Form("%s/PoPSDSigmaPerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-------------------------
	TCanvas *cPoEnMean = new TCanvas("cPoEnMean"," Po Energy Mean",1000,400);
	grPoEnMean->GetXaxis()->SetTitle(xLabel);
	grPoEnMean->GetYaxis()->SetTitle("Energy [MeV]");
	grPoEnMean->Draw("AP");
	grPoEnMeanET->Draw("P");
	cPoEnMean->SaveAs(Form("%s/PoEnMeanPerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoEnSigma = new TCanvas("cPoEnSigma"," Po Energy Sigma",1000,400);
	grPoEnSigma->GetXaxis()->SetTitle(xLabel);
	grPoEnSigma->GetYaxis()->SetTitle("#sigma_{E} [MeV]");
	grPoEnSigma->Draw("AP");
	grPoEnSigmaET->Draw("P");
	cPoEnSigma->SaveAs(Form("%s/PoEnSigmaPerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-------------------------
/*	TCanvas *cPoPosMean = new TCanvas("cPoPosMean"," Po Position Mean",1000,400);
	grPoPosMean->GetXaxis()->SetTitle(xLabel);
	grPoPosMean->GetYaxis()->SetTitle("z [mm]");
	grPoPosMean->Draw("AP");
	grPoPosMeanET->Draw("P");
	cPoPosMean->SaveAs(Form("%s/PoPosMeanPerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoPosSigma = new TCanvas("cPoPosSigma"," Po Position Sigma",1000,400);
	grPoPosSigma->GetXaxis()->SetTitle(xLabel);
	grPoPosSigma->GetYaxis()->SetTitle("RMS_{z} [mm]");
	grPoPosSigma->Draw("AP");
	grPoPosSigmaET->Draw("P");
	cPoPosSigma->SaveAs(Form("%s/PoPosSigmaPerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
*/
	//-------------------------
	TCanvas *cRnPoDzMean = new TCanvas("cRnPoDzMean"," RnPo Dz Mean",1000,400);
	grRnPoDzMean->GetXaxis()->SetTitle(xLabel);
	grRnPoDzMean->GetYaxis()->SetTitle("#Deltaz [mm]");
	grRnPoDzMean->Draw("AP");
	grRnPoDzMeanET->Draw("P");
	cRnPoDzMean->SaveAs(Form("%s/RnPoDzMeanPerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRnPoDzSigma = new TCanvas("cRnPoDzSigma"," RnPo Dz Sigma",1000,400);
	grRnPoDzSigma->GetXaxis()->SetTitle(xLabel);
	grRnPoDzSigma->GetYaxis()->SetTitle("#sigma_{#Deltaz} [mm]");
	grRnPoDzSigma->Draw("AP");
	grRnPoDzSigmaET->Draw("P");
	cRnPoDzSigma->SaveAs(Form("%s/RnPoDzSigmaPerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-------------------------
	TCanvas *cEff = new TCanvas("cEff","Efficiency",1000,400);
	gPad->SetGrid();
	grTotEff->GetXaxis()->SetTitle("Segment");
	grTotEff->GetYaxis()->SetTitle("Efficiency");
	grTotEff->GetXaxis()->SetTitleOffset(0.89);
	grTotEff->GetYaxis()->SetRangeUser(0.998,1.0002);

	grRnPSDEff->SetMarkerStyle(26);
	grRnPSDEff->SetMarkerColor(kBlack);
	grRnPSDEff->SetLineColor(kBlack);
	grPoPSDEff->SetMarkerStyle(26);
	grPoPSDEff->SetMarkerColor(kRed);
	grPoPSDEff->SetLineColor(kRed);
	grRnEnEff->SetMarkerStyle(25);
	grRnEnEff->SetMarkerColor(kBlack);
	grRnEnEff->SetLineColor(kBlack);
	grPoEnEff->SetMarkerStyle(25);
	grPoEnEff->SetMarkerColor(kRed);
	grPoEnEff->SetLineColor(kRed);
	grRnPoDzEff->SetMarkerStyle(27);
	grRnPoDzEff->SetMarkerColor(8);
	grRnPoDzEff->SetLineColor(8);

	grTotEff->Draw("AP");
	grRnPSDEff->Draw("P");
	grPoPSDEff->Draw("P");
	grRnEnEff->Draw("P");
	grPoEnEff->Draw("P");
	grRnPoDzEff->Draw("P");

	TLegend *leg = new TLegend(0.90,0.67,0.99,0.99);
	leg->AddEntry(grRnPSDEff,"Rn PSD","p");
	leg->AddEntry(grPoPSDEff,"Po PSD","p");
	leg->AddEntry(grRnEnEff,"Rn E","p");
	leg->AddEntry(grPoEnEff,"Po E","p");
	leg->AddEntry(grRnPoDzEff,"Dz","p");
	leg->AddEntry(grTotEff,"Total","p");
	leg->Draw();	

	cEff->SaveAs(Form("%s/EfficiencyPerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//--------------------------------------------------
	double grLifetimexStart, grLifetimexEnd, grLifetimeyStart, grLifetimeyEnd;
	grLifetime->GetPoint(0,grLifetimexStart,grLifetimeyStart);
	grLifetime->GetPoint(grLifetime->GetN()-1,grLifetimexEnd,grLifetimeyEnd);
	grLifetimexStart = grLifetimexStart - 3;
	grLifetimexEnd = grLifetimexEnd + 3;

	TLine *lLifetime = new TLine(grLifetimexStart,2.569,grLifetimexEnd,2.569);
	lLifetime->SetLineStyle(2);
	lLifetime->SetLineColor(30);

	TGraph *grlLifetime = new TGraph(4);
	grlLifetime->SetPoint(0,grLifetimexStart,2.562);
	grlLifetime->SetPoint(1,grLifetimexEnd,2.562);
	grlLifetime->SetPoint(2,grLifetimexEnd,2.576);
	grlLifetime->SetPoint(3,grLifetimexStart,2.576);
	grlLifetime->SetFillStyle(3001);
	grlLifetime->SetFillColor(30);
	
	TCanvas *cLifetime = new TCanvas("cLifetime","Lifetime",1000,400);
	gPad->SetGrid();
	grLifetime->GetXaxis()->SetTitle(xLabel);
	grLifetime->GetYaxis()->SetTitle("#tau_{Po} [ms]");	
	grLifetime->Draw("AP");
	lLifetime->Draw("same");
	grlLifetime->Draw("f");
	grLifetime->Draw("P");
	grLifetime->Fit("pol0");
	grLifetime->GetFunction("pol0")->SetLineStyle(2);
	TPaveText *pv = new TPaveText(0.85,0.8,0.99,0.99,"NDC");
	pv->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grLifetime->GetFunction("pol0")->GetChisquare(),grLifetime->GetFunction("pol0")->GetNDF()));	
	pv->AddText(Form("Prob   %f",grLifetime->GetFunction("pol0")->GetProb()));
	pv->AddText(Form("p0   %.3f #pm %.3f",grLifetime->GetFunction("pol0")->GetParameter(0),grLifetime->GetFunction("pol0")->GetParError(0)));
	pv->Draw();
	cLifetime->SaveAs(Form("%s/LifetimePerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//---------------------------------
	grRnPSDChiSq->SetMarkerStyle(22);
	grRnPSDChiSq->SetMarkerSize(1.1);
	grRnPSDChiSq->SetMarkerColor(kBlack);
	grRnPSDChiSq->SetLineColor(kBlack);
	grPoPSDChiSq->SetMarkerStyle(22);
	grPoPSDChiSq->SetMarkerSize(1.1);
	grPoPSDChiSq->SetMarkerColor(kRed);
	grPoPSDChiSq->SetLineColor(kRed);

	grRnEnChiSq->SetMarkerStyle(34);
	grRnEnChiSq->SetMarkerSize(1.1);
	grRnEnChiSq->SetMarkerColor(kBlack);
	grRnEnChiSq->SetLineColor(kBlack);
	grPoEnChiSq->SetMarkerStyle(34);
	grPoEnChiSq->SetMarkerSize(1.1);
	grPoEnChiSq->SetMarkerColor(kRed);
	grPoEnChiSq->SetLineColor(kRed);

	grDzChiSq->SetMarkerStyle(20);
	grDzChiSq->SetMarkerSize(1.1);
	grDzChiSq->SetMarkerColor(8);
	grDzChiSq->SetLineColor(8);

	grDtChiSq->SetMarkerStyle(21);
	grDtChiSq->SetMarkerSize(1.1);
	
	TCanvas *cFitChiSq = new TCanvas("cFitChiSq","Fit Chisquared",1000,400);
	gPad->SetGrid();
	grRnPSDChiSq->GetXaxis()->SetTitle(xLabel);
	grRnPSDChiSq->GetYaxis()->SetTitle("#Chi^{2}/NDF");
	grRnPSDChiSq->GetYaxis()->SetRangeUser(0,4);
	grRnPSDChiSq->Draw("AP");
	grPoPSDChiSq->Draw("P");
	grRnEnChiSq->Draw("P");
	grPoEnChiSq->Draw("P");
	grDzChiSq->Draw("P");
	grDtChiSq->Draw("P");	
	leg = new TLegend(0.90,0.67,0.99,0.99);
	leg->AddEntry(grRnPSDChiSq,"Prompt PSD","p");
	leg->AddEntry(grPoPSDChiSq,"Delay PSD","p");
	leg->AddEntry(grRnEnChiSq,"Prompt E","p");
	leg->AddEntry(grPoEnChiSq,"Delay E","p");
	leg->AddEntry(grDzChiSq,"Dz","p");
	leg->AddEntry(grDtChiSq,"Dt","p");
	leg->Draw();	
	cFitChiSq->SaveAs(Form("%s/FitChiSqPerCell.pdf",gSystem->Getenv("AD_AC227_PLOTS")));



} 	//end PlotRnPoVsCell
