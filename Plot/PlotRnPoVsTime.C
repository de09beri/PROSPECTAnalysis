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

//Function for creating a relative to mean graph
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

//Function for creating Rate Graph corrected for the Ac227 lifetime
TGraphErrors *makeCorrGr(TGraphErrors *gr, TF1 *fExp){
	TGraphErrors *grCorr = (TGraphErrors*)gr->Clone();

	int numPt = gr->GetN();

	double t0, N0;
	gr->GetPoint(0,t0,N0);

	double grx,gry,gryErr;
	double grCorry,grCorryErr;

	double lambda = log(2)/(21.772*365.0*24.0*60.0*60.0);

	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);	

		grCorry = gry/(N0*exp(-((grx-t0)*lambda)));
		grCorryErr = gryErr/(N0*exp(-((grx-t0)*lambda)));
	
		grCorr->SetPoint(i,grx,grCorry);
		grCorr->SetPointError(i,0,grCorryErr);	
	}

	return grCorr;
}


void PlotRnPoVsTime(){

	setup_PROSPECT_style();
    	gROOT->ForceStyle();

	TFile *f = new TFile(Form("%s/Ac227_GraphsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));

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

	TGraph *grLivetime    = (TGraph*)f->Get("grLivetime");
	TGraph *grTotLivetime = (TGraph*)f->Get("grTotLivetime");
	TGraph *grPileupVeto  = (TGraph*)f->Get("grPileupVeto");
	TGraph *grPileupVetoFrac  = (TGraph*)f->Get("grPileupVetoFrac");
	
	TGraphErrors *grPromptEnEff  = (TGraphErrors*)f->Get("grPromptEnEff");
	TGraphErrors *grDelayEnEff   = (TGraphErrors*)f->Get("grDelayEnEff");
	TGraphErrors *grPromptPSDEff = (TGraphErrors*)f->Get("grPromptPSDEff");
	TGraphErrors *grDelayPSDEff  = (TGraphErrors*)f->Get("grDelayPSDEff");
	TGraphErrors *grDzEff  	     = (TGraphErrors*)f->Get("grDzEff");

	TGraph *grRnPSDChiSq = (TGraph*)f->Get("grRnPSDChiSq");
	TGraph *grPoPSDChiSq = (TGraph*)f->Get("grPoPSDChiSq");
	TGraph *grRnEnChiSq  = (TGraph*)f->Get("grRnEnChiSq");
	TGraph *grPoEnChiSq  = (TGraph*)f->Get("grPoEnChiSq");
	TGraph *grDzChiSq    = (TGraph*)f->Get("grDzChiSq");
	TGraph *grDtChiSq    = (TGraph*)f->Get("grDtChiSq");

	TGraph *grBGRate = (TGraph*)f->Get("grBGRate");	

	f->Close();

/*
grRate->SetPoint(87,1556726400,0.275617);
grRate->SetPointError(87,0,0.0015);
grRate->SetPoint(88,1556899200,0.275508);
grRate->SetPointError(88,0,0.0015);
grRate->SetPoint(89,1557072000,0.275400);
grRate->SetPointError(89,0,0.00155);
*/


	//-------------------------------------------------------------------------------------------------------
	const char *xLabel = "Date";
	TLegend *leg;

	int numPt = grRate->GetN();
	double grxStart, grxEnd, gryStart, gryEnd;
	grRate->GetPoint(0,grxStart,gryStart);
	grRate->GetPoint(numPt-1,grxEnd,gryEnd);

	double Feb20_2018 = 1519102800;		//Rx On
	double March16_2018 = 1521201600;	//Rx Off 
	double May01_2018  = 1525168800;	//Rx On
	double May25_2018  = 1527285600;	//Rx Off
	double June12_2018 = 1528797600;	//Rx On	
	double July06_2018 = 1530914400;	//Rx Off
	double July24_2018 = 1532426400;	//Rx On
	double Aug17_2018  = 1534534200;	//Rx Off
	double Sept04_2018 = 1536055200;	//Rx On
	double Sept28_2018 = 1538172000; 	//Rx Off

	TGraph *grMarchOn_2018 = new TGraph(4);
	TGraph *grMayOn_2018  = new TGraph(4);
	TGraph *grJuneOn_2018 = new TGraph(4);
	TGraph *grJulyOn_2018 = new TGraph(4);
	TGraph *grSeptOn_2018 = new TGraph(4);

	grMarchOn_2018->SetPoint(0,Feb20_2018,-400);
	grMarchOn_2018->SetPoint(1,Feb20_2018,400e9);
	grMarchOn_2018->SetPoint(2,March16_2018,400e9);
	grMarchOn_2018->SetPoint(3,March16_2018,-400);

	grMayOn_2018->SetPoint(0,May01_2018,-400);
	grMayOn_2018->SetPoint(1,May01_2018,400e9);
	grMayOn_2018->SetPoint(2,May25_2018,400e9);
	grMayOn_2018->SetPoint(3,May25_2018,-400);
	
	grJuneOn_2018->SetPoint(0,June12_2018,-400);
	grJuneOn_2018->SetPoint(1,June12_2018,400e9);
	grJuneOn_2018->SetPoint(2,July06_2018,400e9);
	grJuneOn_2018->SetPoint(3,July06_2018,-400);

	grJulyOn_2018->SetPoint(0,July24_2018,-400);
	grJulyOn_2018->SetPoint(1,July24_2018,400e9);
	grJulyOn_2018->SetPoint(2,Aug17_2018,400e9);
	grJulyOn_2018->SetPoint(3,Aug17_2018,-400);

	grSeptOn_2018->SetPoint(0,Sept04_2018,-400);
	grSeptOn_2018->SetPoint(1,Sept04_2018,400e9);
	grSeptOn_2018->SetPoint(2,Sept28_2018,400e9);
	grSeptOn_2018->SetPoint(3,Sept28_2018,-400);

	grMarchOn_2018->SetFillStyle(3002);
	grMarchOn_2018->SetFillColor(16);
	grMayOn_2018->SetFillStyle(3002);
	grMayOn_2018->SetFillColor(16);
	grJuneOn_2018->SetFillStyle(3002);
	grJuneOn_2018->SetFillColor(16);
	grJulyOn_2018->SetFillStyle(3002);
	grJulyOn_2018->SetFillColor(16);
	grSeptOn_2018->SetFillStyle(3002);
	grSeptOn_2018->SetFillColor(16);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	double convertT = 365.0*24.0*60.0*60.0;

	double useX,useY;
	grRate->GetPoint(2,useX,useY);

	TF1 *fAcExpDraw = new TF1("fAcExpDraw",Form("[0]*exp(-((x-%f)*log(2))/([1]*%f))",useX,convertT),grxStart,grxEnd);
	fAcExpDraw->FixParameter(0,useY);
	fAcExpDraw->FixParameter(1,21.772);
	fAcExpDraw->SetLineColor(kViolet);

	TCanvas *cRate = new TCanvas("cRate","Rate vs time",1000,400);
	grRate->GetXaxis()->SetTitle(xLabel);
	grRate->GetYaxis()->SetTitle("R_{RnPo} [Hz]");
	grRate->GetXaxis()->SetTimeDisplay(1);
	grRate->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grRate->GetXaxis()->SetTimeOffset(0,"est");
	grRate->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRate->Draw("P");
	//fAcExpDraw->Draw("same");
	cRate->SaveAs(Form("%s/RateVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-------------------
	//Fit rate 

	TF1 *fAcExp = new TF1("fAcExp",Form("[0]*exp(-((x-%f)*log(2))/([1]*%f))",grxStart,convertT),grxStart,grxEnd);
	fAcExp->SetParameters(0.32,21.772);
	fAcExp->SetParLimits(1,0,100);

	TF1 *fAcExpFixed = new TF1("fAcExpFixed",Form("[0]*exp(-((x-%f)*log(2))/([1]*%f))",grxStart,convertT),grxStart,grxEnd);
	fAcExp->SetParameters(0.32,21.772);
	fAcExpFixed->FixParameter(1,21.772);

	TCanvas *cRateFit = new TCanvas("cRate","Rate vs time Fit",1000,400);
	grRate->GetXaxis()->SetTitle();
	grRate->GetYaxis()->SetTitle("R_{RnPo} [Hz]");
	grRate->GetXaxis()->SetTimeDisplay(1);
	grRate->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grRate->GetXaxis()->SetTimeOffset(0,"est");
	grRate->Draw("AP");
	grRate->Fit(fAcExp,"0RB");
	grRate->Fit(fAcExpFixed,"0R");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRate->Draw("P");
	fAcExp->Draw("same");
	fAcExpFixed->SetLineColor(kBlack);
	fAcExpFixed->Draw("same");
	TPaveText *pv = new TPaveText(0.4,0.8,0.6,0.98,"NDC");
	pv->SetShadowColor(kRed);
	pv->AddText(Form("#Chi^{2}/NDF  %.1f/%d",fAcExp->GetChisquare(),fAcExp->GetNDF()));
	pv->AddText(Form("Prob  %.3f",fAcExp->GetProb()));
	pv->AddText(Form("R_{0}   %.2f #pm %.2f mHz",fAcExp->GetParameter(0)*1000.0,fAcExp->GetParError(0)*1000.0));
	pv->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExp->GetParameter(1),fAcExp->GetParError(1)));
	pv->Draw();
	TPaveText *pvFix = new TPaveText(0.65,0.8,0.85,0.98,"NDC");
	pvFix->SetShadowColor(kBlack);
	pvFix->AddText(Form("#Chi^{2}/NDF  %.1f/%d",fAcExpFixed->GetChisquare(),fAcExpFixed->GetNDF()));
	pvFix->AddText(Form("Prob  %.3f",fAcExpFixed->GetProb()));
	pvFix->AddText(Form("R_{0}   %.2f #pm %.2f mHz",fAcExpFixed->GetParameter(0)*1000.0,fAcExpFixed->GetParError(0)*1000.0));
	pvFix->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExpFixed->GetParameter(1),fAcExpFixed->GetParError(1)));
	pvFix->Draw();
	cRateFit->SaveAs(Form("%s/RateVsTime_Fit.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-------------------
	//Fit rate corrected for Ac227 lifetime
	TLine *lRateCorr = new TLine(grxStart,1.0,grxEnd,1.0);
	lRateCorr->SetLineStyle(2);
	lRateCorr->SetLineColor(kBlack);

	TGraphErrors *grRateCorr = makeCorrGr(grRate,fAcExpFixed);
	TCanvas *cRateCorr = new TCanvas("cRateCorr","Corrected rate",1000,400);
	gPad->SetGrid();
	grRateCorr->GetXaxis()->SetTitle(xLabel);	
	grRateCorr->GetYaxis()->SetTitle("R_{RnPo}/R_{Ac}");
	grRateCorr->GetXaxis()->SetTimeDisplay(1);
	grRateCorr->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grRateCorr->GetXaxis()->SetTimeOffset(0,"est");
	grRateCorr->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRateCorr->Draw("P");
	lRateCorr->Draw("same");
	grRateCorr->Fit("pol0");
	TPaveText *pvPol = new TPaveText(0.8,0.8,0.99,0.99,"NDC");
	pvPol->AddText(Form("#Chi^{2}/NDF  %.1f/%d",grRateCorr->GetFunction("pol0")->GetChisquare(),grRateCorr->GetFunction("pol0")->GetNDF()));
	pvPol->AddText(Form("Prob  %f",grRateCorr->GetFunction("pol0")->GetProb()));
	pvPol->AddText(Form("p0  %.3f #pm %.3f",grRateCorr->GetFunction("pol0")->GetParameter(0),grRateCorr->GetFunction("pol0")->GetParError(0)));	
	pvPol->Draw();
	cRateCorr->SaveAs(Form("%s/RateVsTime_Corrected.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	

	//======================================================================================
	TCanvas *cPoPSDMean = new TCanvas("cPoPSDMean","Po PSD Mean",1000,400);
	grPoPSDMean->GetYaxis()->SetTitle("PSD");
	grPoPSDMean->GetXaxis()->SetTitle(xLabel);
	grPoPSDMean->GetXaxis()->SetTimeDisplay(1);
	grPoPSDMean->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grPoPSDMean->GetXaxis()->SetTimeOffset(0,"est");
	grPoPSDMean->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoPSDMean->Draw("P");
	cPoPSDMean->SaveAs(Form("%s/PoPSDMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoPSDSigma = new TCanvas("cPoPSDSigma","Po PSD Sigma",1000,400);
	grPoPSDSigma->GetYaxis()->SetTitle("#sigma_{PSD}");
	grPoPSDSigma->GetXaxis()->SetTitle(xLabel);
	grPoPSDSigma->GetXaxis()->SetTimeDisplay(1);
	grPoPSDSigma->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grPoPSDSigma->GetXaxis()->SetTimeOffset(0,"est");
	grPoPSDSigma->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoPSDSigma->Draw("P");
	cPoPSDSigma->SaveAs(Form("%s/PoPSDSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-------------------------------
	TCanvas *cPoEnMean = new TCanvas("cPoEnMean","Po Energy Mean",1000,400);
	grPoEnMean->GetYaxis()->SetTitle("Energy [MeV]");
	grPoEnMean->GetXaxis()->SetTitle(xLabel);
	grPoEnMean->GetXaxis()->SetTimeDisplay(1);
	grPoEnMean->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grPoEnMean->GetXaxis()->SetTimeOffset(0,"est");
	grPoEnMean->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoEnMean->Draw("P");
	cPoEnMean->SaveAs(Form("%s/PoEnMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoEnSigma = new TCanvas("cPoEnSigma","Po Energy Sigma",1000,400);
	grPoEnSigma->GetYaxis()->SetTitle("#sigma_{E} [MeV]");
	grPoEnSigma->GetXaxis()->SetTitle(xLabel);
	grPoEnSigma->GetXaxis()->SetTimeDisplay(1);
	grPoEnSigma->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grPoEnSigma->GetXaxis()->SetTimeOffset(0,"est");
	grPoEnSigma->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoEnSigma->Draw("P");
	cPoEnSigma->SaveAs(Form("%s/PoEnSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoEnSmearMean = new TCanvas("cPoEnSmearMean","Po ESmear Mean",1000,400);
	grPoEnSmearMean->GetYaxis()->SetTitle("ESmear [MeV]");
	grPoEnSmearMean->GetXaxis()->SetTitle(xLabel);
	grPoEnSmearMean->GetXaxis()->SetTimeDisplay(1);
	grPoEnSmearMean->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grPoEnSmearMean->GetXaxis()->SetTimeOffset(0,"est");
	grPoEnSmearMean->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoEnSmearMean->Draw("P");
	cPoEnSmearMean->SaveAs(Form("%s/PoEnSmearMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoEnSmearSigma = new TCanvas("cPoEnSmearSigma","Po ESmear Sigma",1000,400);
	grPoEnSmearSigma->GetYaxis()->SetTitle("#sigma_{ESmear} [MeV]");
	grPoEnSmearSigma->GetXaxis()->SetTitle(xLabel);
	grPoEnSmearSigma->GetXaxis()->SetTimeDisplay(1);
	grPoEnSmearSigma->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grPoEnSmearSigma->GetXaxis()->SetTimeOffset(0,"est");
	grPoEnSmearSigma->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoEnSmearSigma->Draw("P");
	cPoEnSmearSigma->SaveAs(Form("%s/PoEnSmearSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoEnMeanCompare = new TCanvas("cPoEnMeanCompare","Po En Mean Compare",1000,400);
	grPoEnMean->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoEnMean->Draw("P");
	grPoEnSmearMean->SetMarkerColor(kRed);
	grPoEnSmearMean->SetLineColor(kRed);
	grPoEnSmearMean->Draw("P");	
	leg = new TLegend(0.90,0.8,0.99,0.95);
	leg->AddEntry(grPoEnMean,"E","p");
	leg->AddEntry(grPoEnSmearMean,"ESmear","p");
	leg->Draw();	
	cPoEnMeanCompare->SaveAs(Form("%s/PoEnMeanCompareVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoEnSigmaCompare = new TCanvas("cPoEnSigmaCompare","Po En Sigma Compare",1000,400);
	grPoEnSigma->GetYaxis()->SetRangeUser(38e-3,53e-3);
	grPoEnSigma->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoEnSigma->Draw("P");
	grPoEnSmearSigma->SetMarkerColor(kRed);
	grPoEnSmearSigma->SetLineColor(kRed);
	grPoEnSmearSigma->Draw("P");	
	leg = new TLegend(0.90,0.8,0.99,0.95);
	leg->AddEntry(grPoEnSigma,"E","p");
	leg->AddEntry(grPoEnSmearSigma,"ESmear","p");
	leg->Draw();	
	cPoEnSigmaCompare->SaveAs(Form("%s/PoEnSigmaCompareVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoEnMeanRelativeCompare = new TCanvas("cPoEnMeanRelativeCompare","Po En Mean Compare",1000,400);
	TGraphErrors *grPoEnMeanRelative = makeRelGr(grPoEnMean);
	grPoEnMeanRelative->GetYaxis()->SetTitle("E/#LTE#GT");
	grPoEnMeanRelative->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoEnMeanRelative->Draw("P");
	TGraphErrors *grPoEnSmearMeanRelative = makeRelGr(grPoEnSmearMean);
	grPoEnSmearMeanRelative->SetMarkerColor(kRed);
	grPoEnSmearMeanRelative->SetLineColor(kRed);
	grPoEnSmearMeanRelative->Draw("P");	
	leg = new TLegend(0.90,0.8,0.99,0.95);
	leg->AddEntry(grPoEnMeanRelative,"E","p");
	leg->AddEntry(grPoEnSmearMeanRelative,"ESmear","p");
	leg->Draw();	
	cPoEnMeanRelativeCompare->SaveAs(Form("%s/PoEnMeanRelativeCompareVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoEnSigmaRelativeCompare = new TCanvas("cPoEnSigmaRelativeCompare","Po En Sigma Compare",1000,400);
	TGraphErrors *grPoEnSigmaRelative = makeRelGr(grPoEnSigma);
	grPoEnSigmaRelative->GetYaxis()->SetTitle("#sigma_{E}/#LT#sigma_{E}#GT");
	grPoEnSigmaRelative->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoEnSigmaRelative->Draw("P");
	TGraphErrors *grPoEnSmearSigmaRelative = makeRelGr(grPoEnSmearSigma);
	grPoEnSmearSigmaRelative->SetMarkerColor(kRed);
	grPoEnSmearSigmaRelative->SetLineColor(kRed);
	grPoEnSmearSigmaRelative->Draw("P");	
	leg = new TLegend(0.90,0.8,0.99,0.95);
	leg->AddEntry(grPoEnSigmaRelative,"E","p");
	leg->AddEntry(grPoEnSmearSigmaRelative,"ESmear","p");
	leg->Draw();	
	cPoEnSigmaRelativeCompare->SaveAs(Form("%s/PoEnSigmaRelativeCompareVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	
	//-------------------------------
/*
	TCanvas *cPoPosMean = new TCanvas("cPoPosMean","Po Position Mean",1000,400);
	grPoPosMean->GetXaxis()->SetTitle(xLabel);
	grPoPosMean->GetXaxis()->SetTimeDisplay(1);
	grPoPosMean->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grPoPosMean->GetXaxis()->SetTimeOffset(0,"est");
	grPoPosMean->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoPosMean->Draw("P");
	cPoPosMean->SaveAs(Form("%s/PoPosMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoPosSigma = new TCanvas("cPoPosSigma","Po Position Sigma",1000,400);
	grPoPosSigma->GetXaxis()->SetTitle(xLabel);
	grPoPosSigma->GetXaxis()->SetTimeDisplay(1);
	grPoPosSigma->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grPoPosSigma->GetXaxis()->SetTimeOffset(0,"est");
	grPoPosSigma->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPoPosSigma->Draw("P");
	cPoPosSigma->SaveAs(Form("%s/PoPosSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
*/
	//-------------------------------
	TCanvas *cRnPoDzMean = new TCanvas("cRnPoDzMean","RnPo Dz Mean",1000,400);
	grRnPoDzMean->GetYaxis()->SetTitle("#Deltaz [mm]");
	grRnPoDzMean->GetXaxis()->SetTitle(xLabel);
	grRnPoDzMean->GetXaxis()->SetTimeDisplay(1);
	grRnPoDzMean->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grRnPoDzMean->GetXaxis()->SetTimeOffset(0,"est");
	grRnPoDzMean->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRnPoDzMean->Draw("P");
	cRnPoDzMean->SaveAs(Form("%s/RnPoDzMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRnPoDzSigma = new TCanvas("cRnPoDzSigma","RnPo Dz Sigma",1000,400);
	grRnPoDzSigma->GetYaxis()->SetTitle("#sigma_{#Deltaz} [mm]");
	grRnPoDzSigma->GetXaxis()->SetTitle(xLabel);
	grRnPoDzSigma->GetXaxis()->SetTimeDisplay(1);
	grRnPoDzSigma->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grRnPoDzSigma->GetXaxis()->SetTimeOffset(0,"est");
	grRnPoDzSigma->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRnPoDzSigma->Draw("P");
	cRnPoDzSigma->SaveAs(Form("%s/RnPoDzSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-------------------------------
	TCanvas *cEff = new TCanvas("cEff","Efficiency",1000,400);
	grTotEff->GetXaxis()->SetTitle(xLabel);
	grTotEff->GetYaxis()->SetTitle("Efficiency");
	grTotEff->GetXaxis()->SetTimeDisplay(1);
	grTotEff->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grTotEff->GetXaxis()->SetTimeOffset(0,"est");
	grTotEff->GetYaxis()->SetRangeUser(0.9992,1.0002);
	grTotEff->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grTotEff->Draw("P");

	grPromptPSDEff->SetMarkerStyle(26);
	grPromptPSDEff->SetMarkerColor(kBlack);
	grPromptPSDEff->SetLineColor(kBlack);
	grDelayPSDEff->SetMarkerStyle(26);
	grDelayPSDEff->SetMarkerColor(kRed);
	grDelayPSDEff->SetLineColor(kRed);
	grPromptEnEff->SetMarkerStyle(25);
	grPromptEnEff->SetMarkerColor(kBlack);
	grPromptEnEff->SetLineColor(kBlack);
	grDelayEnEff->SetMarkerStyle(25);
	grDelayEnEff->SetMarkerColor(kRed);
	grDelayEnEff->SetLineColor(kRed);
	grDzEff->SetMarkerStyle(27);
	grDzEff->SetMarkerColor(8);
	grDzEff->SetLineColor(8);

	grPromptPSDEff->Draw("P");
	grDelayPSDEff->Draw("P");
	grPromptEnEff->Draw("P");
	grDelayEnEff->Draw("P");
	grDzEff->Draw("P");

	leg = new TLegend(0.90,0.67,0.99,0.99);
	leg->AddEntry(grPromptPSDEff,"Prompt PSD","p");
	leg->AddEntry(grDelayPSDEff,"Delay PSD","p");
	leg->AddEntry(grPromptEnEff,"Prompt E","p");
	leg->AddEntry(grDelayEnEff,"Delay E","p");
	leg->AddEntry(grDzEff,"Dz","p");
	leg->AddEntry(grTotEff,"Total","p");
	leg->Draw();	
	cEff->SaveAs(Form("%s/EfficiencyVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-------------------------------
	TCanvas *cVetoTime = new TCanvas("cVetoTime","VetoTime",1000,400);
	gPad->SetGrid();
	grPileupVetoFrac->GetXaxis()->SetTitle(xLabel);
	grPileupVetoFrac->GetYaxis()->SetTitle("t_{veto}/t_{live}");
	grPileupVetoFrac->GetXaxis()->SetTimeDisplay(1);
	grPileupVetoFrac->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grPileupVetoFrac->GetXaxis()->SetTimeOffset(0,"est");
	grPileupVetoFrac->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPileupVetoFrac->Draw("P");
	cVetoTime->SaveAs(Form("%s/VetoTimeVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	
	//-----------------------
	//Fit Chisq
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
	grRnPSDChiSq->GetXaxis()->SetTimeDisplay(1);
	grRnPSDChiSq->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grRnPSDChiSq->GetXaxis()->SetTimeOffset(0,"est");
	grRnPSDChiSq->GetYaxis()->SetRangeUser(0,3);
	grRnPSDChiSq->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRnPSDChiSq->Draw("P");
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
	cFitChiSq->SaveAs(Form("%s/FitChiSqVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	
	//-------------------------------
	TCanvas *cBGRate = new TCanvas("cBGRate","BG Rate",1000,400);
	gPad->SetGrid();
	grBGRate->GetXaxis()->SetTitle(xLabel);
	grBGRate->GetYaxis()->SetTitle("BG Rate [Hz]");
	grBGRate->GetXaxis()->SetTimeDisplay(1);
	grBGRate->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grBGRate->GetXaxis()->SetTimeOffset(0,"est");
	grBGRate->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grBGRate->Draw("P");
	cBGRate->SaveAs(Form("%s/BGRateVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//----------------------------------------------------
	double grLifetimexStart, grLifetimexEnd, grLifetimeyStart, grLifetimeyEnd;
	grLifetime->GetPoint(0,grLifetimexStart,grLifetimeyStart);
	grLifetime->GetPoint(grLifetime->GetN()-1,grLifetimexEnd,grLifetimeyEnd);

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
	grLifetime->GetXaxis()->SetTimeDisplay(1);
	grLifetime->GetXaxis()->SetTimeOffset(0,"est");
	grLifetime->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grLifetime->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	lLifetime->Draw("same");
	grlLifetime->Draw("f");
	grLifetime->Draw("P");
	grLifetime->Fit("pol0");
	grLifetime->GetFunction("pol0")->SetLineStyle(2);
	pv = new TPaveText(0.85,0.8,0.99,0.99,"NDC");
	pv->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grLifetime->GetFunction("pol0")->GetChisquare(),grLifetime->GetFunction("pol0")->GetNDF()));	
	pv->AddText(Form("Prob   %.3f",grLifetime->GetFunction("pol0")->GetProb()));
	pv->AddText(Form("p0   %.3f #pm %.3f",grLifetime->GetFunction("pol0")->GetParameter(0),grLifetime->GetFunction("pol0")->GetParError(0)));
	pv->Draw();
	cLifetime->SaveAs(Form("%s/LifetimeVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-------------------------------
	TCanvas *cLivetime = new TCanvas("cLivetime","livetime",1000,400);
	gPad->SetGrid();
	grLivetime->GetXaxis()->SetTitle(xLabel);
	grLivetime->GetYaxis()->SetTitle("Livetime [ms]");
	grLivetime->GetXaxis()->SetTimeDisplay(1);
	grLivetime->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grLivetime->GetXaxis()->SetTimeOffset(0,"est");
	grLivetime->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grLivetime->Draw("P");
	cLivetime->SaveAs(Form("%s/LivetimeVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));


} 	//end PlotRnPoVsTime
