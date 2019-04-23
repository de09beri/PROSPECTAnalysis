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

TGraphErrors *makeRelGr(TGraphErrors *gr, double mean, double meanErr){
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

TGraphErrors *makeCorrGr(TGraphErrors *gr, TF1 *fExp){
	gr->RemovePoint(0);
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
	
	//-------------------------------------------------------------------------------------------------------
cout<<"\n RATE"<<endl;
	grRate->Fit("pol0","0");
	TGraphErrors *grRelRate = makeRelGr(grRate, grRate->GetFunction("pol0")->GetParameter(0), grRate->GetFunction("pol0")->GetParError(0));

cout<<"\n PSD MEAN"<<endl;
	grPoPSDMean->Fit("pol0","0");
	TGraphErrors *grRelPoPSDMean = makeRelGr(grPoPSDMean, grPoPSDMean->GetFunction("pol0")->GetParameter(0), grPoPSDMean->GetFunction("pol0")->GetParError(0));

cout<<"\n PSD SIGMA"<<endl;
	grPoPSDSigma->Fit("pol0","0");
	TGraphErrors *grRelPoPSDSigma = makeRelGr(grPoPSDSigma, grPoPSDSigma->GetFunction("pol0")->GetParameter(0), grPoPSDSigma->GetFunction("pol0")->GetParError(0));

cout<<"\n ENERGY MEAN"<<endl;
	grPoEnMean->Fit("pol0","0");
	TGraphErrors *grRelPoEnMean = makeRelGr(grPoEnMean, grPoEnMean->GetFunction("pol0")->GetParameter(0), grPoEnMean->GetFunction("pol0")->GetParError(0));

cout<<"\n ENERGY SIGMA"<<endl;
	grPoEnSigma->Fit("pol0","0");
	TGraphErrors *grRelPoEnSigma = makeRelGr(grPoEnSigma, grPoEnSigma->GetFunction("pol0")->GetParameter(0), grPoEnSigma->GetFunction("pol0")->GetParError(0));

cout<<"\n POSITION MEAN"<<endl;
	grPoPosMean->Fit("pol0","0");
	TGraphErrors *grRelPoPosMean = makeRelGr(grPoPosMean, grPoPosMean->GetFunction("pol0")->GetParameter(0), grPoPosMean->GetFunction("pol0")->GetParError(0));

cout<<"\n POSITION RMS"<<endl;
	grPoPosSigma->Fit("pol0","0");
	TGraphErrors *grRelPoPosSigma = makeRelGr(grPoPosSigma, grPoPosSigma->GetFunction("pol0")->GetParameter(0), grPoPosSigma->GetFunction("pol0")->GetParError(0));

cout<<"\n DZ MEAN"<<endl;
	grRnPoDzMean->Fit("pol0","0");
	TGraphErrors *grRelRnPoDzMean = makeRelGr(grRnPoDzMean, grRnPoDzMean->GetFunction("pol0")->GetParameter(0), grRnPoDzMean->GetFunction("pol0")->GetParError(0));

cout<<"\n DZ SIGMA"<<endl;
	grRnPoDzSigma->Fit("pol0","0");
	TGraphErrors *grRelRnPoDzSigma = makeRelGr(grRnPoDzSigma, grRnPoDzSigma->GetFunction("pol0")->GetParameter(0), grRnPoDzSigma->GetFunction("pol0")->GetParError(0));
	
	//-------------------------------------------------------------------------------------------------------
	const char *xLabel = "Time";

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
	TCanvas *cRate = new TCanvas("cRate","Rate vs time",1000,400);
	grRate->GetXaxis()->SetTitle(xLabel);
	grRate->GetYaxis()->SetTitle("R_{RnPo} [Hz]");
	grRate->GetXaxis()->SetTimeDisplay(1);
	grRate->GetXaxis()->SetTimeFormat("%m/%d");
	grRate->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRate->Draw("P");
	cRate->SaveAs(Form("%s/RateVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRate->SaveAs(Form("%s/RateVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRate = new TCanvas("cRelRate","Relative rate vs time",1000,400);
	grRelRate->GetXaxis()->SetTitle(xLabel);
	grRelRate->GetYaxis()->SetTitle("R_{RnPo}/#LTR_{RnPo}#GT");  
	grRelRate->GetXaxis()->SetTimeDisplay(1);
	grRelRate->GetXaxis()->SetTimeFormat("%m/%d");
	grRelRate->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRelRate->Draw("P");
	TPaveText *pvR = new TPaveText(0.80,0.87,0.99,0.99,"NDCNB");
	pvR->SetFillColor(19);
	pvR->AddText(Form("#LTR_{RnPo}#GT  %.3f #pm %.3f Hz",grRate->GetFunction("pol0")->GetParameter(0),grRate->GetFunction("pol0")->GetParError(0)));
	pvR->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grRate->GetFunction("pol0")->GetChisquare(),grRate->GetFunction("pol0")->GetNDF()));	
	pvR->Draw();
	cRelRate->SaveAs(Form("%s/RelativeRateVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRate->SaveAs(Form("%s/RelativeRateVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TF1 *fAcExp = new TF1("fAcExp","[0]*exp(-((x-[2])*log(2)/([1]*365.0*24.0*60.0*60.0)))",grxStart,grxEnd);
	fAcExp->SetParameters(0.32,21.772);
	fAcExp->FixParameter(2,grxStart);

	TF1 *fAcExpFixed = new TF1("fAcExp","[0]*exp(-((x-[2])*log(2)/([1]*365.0*24.0*60.0*60.0)))",grxStart,grxEnd);
	fAcExpFixed->FixParameter(1,21.772);
	fAcExpFixed->FixParameter(2,grxStart);

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
	pv->AddText(Form("Prob  %.2f",fAcExp->GetProb()));
	pv->AddText(Form("R_{0}   %.2f #pm %.2f mHz",fAcExp->GetParameter(0)*1000.0,fAcExp->GetParError(0)*1000.0));
	pv->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExp->GetParameter(1),fAcExp->GetParError(1)));
	pv->Draw();
	TPaveText *pvFix = new TPaveText(0.65,0.8,0.85,0.98,"NDC");
	pvFix->SetShadowColor(kBlack);
	pvFix->AddText(Form("#Chi^{2}/NDF  %.1f/%d",fAcExpFixed->GetChisquare(),fAcExpFixed->GetNDF()));
	pvFix->AddText(Form("Prob  %.2f",fAcExpFixed->GetProb()));
	pvFix->AddText(Form("R_{0}   %.2f #pm %.2f mHz",fAcExpFixed->GetParameter(0)*1000.0,fAcExpFixed->GetParError(0)*1000.0));
	pvFix->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExpFixed->GetParameter(1),fAcExpFixed->GetParError(1)));
	pvFix->Draw();
	cRateFit->SaveAs(Form("%s/RateVsTime_Fit.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRateFit->SaveAs(Form("%s/RateVsTime_Fit.png",gSystem->Getenv("AD_AC227_PLOTS")));

/*
	TF1 *fAcExpLine = new TF1("fAcExpLine","[0]*( exp(-((x-[2])*log(2)/([1]*365.0*24.0*60.0*60.0))) + [3]*(x-[2]) )",grxStart,grxEnd);
	fAcExpLine->SetParameters(0.32,21.772,grxStart,-1e-11);
	fAcExpLine->FixParameter(2,grxStart);

	TF1 *fAcExpLineFixed = new TF1("fAcExpLine","[0]*( exp(-((x-[2])*log(2)/([1]*365.0*24.0*60.0*60.0))) + [3]*(x-[2]) )",grxStart,grxEnd);
	fAcExpLineFixed->FixParameter(1,21.772);
	fAcExpLineFixed->FixParameter(2,grxStart);
	fAcExpLineFixed->SetParameter(3,-1e-11);

	TCanvas *cRateNewFit = new TCanvas("cRate","Rate vs time Fit",1000,400);
	grRate->GetXaxis()->SetTitle(xLabel);
	grRate->GetYaxis()->SetTitle("R_{RnPo} [Hz]");
	grRate->GetXaxis()->SetTimeDisplay(1);
	grRate->GetXaxis()->SetTimeFormat("%m/%d");
	grRate->Draw("AP");
	grRate->Fit(fAcExpLine,"0RB");
	grRate->Fit(fAcExpLineFixed,"0R");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRate->Draw("P");
	fAcExpLine->Draw("same");
	fAcExpLineFixed->SetLineColor(8);
	fAcExpLineFixed->Draw("same");
	pv = new TPaveText(0.4,0.8,0.6,0.98,"NDC");
	pv->SetShadowColor(kRed);
	pv->AddText(Form("#Chi^{2}/NDF  %.1f/%d",fAcExpLine->GetChisquare(),fAcExpLine->GetNDF()));
	pv->AddText(Form("Prob  %f",fAcExpLine->GetProb()));
	pv->AddText(Form("R_{0}   %.2f #pm %.2f mHz",fAcExpLine->GetParameter(0)*1000.0,fAcExpLine->GetParError(0)*1000.0));
	pv->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExpLine->GetParameter(1),fAcExpLine->GetParError(1)));
	pv->AddText(Form("R_{1}   %.4f #pm %.4f yr^{-1}",fAcExpLine->GetParameter(3)*60.0*60.0*24.0*365.0,fAcExpLine->GetParError(3)*60.0*60.0*24.0*365.0));
	pv->Draw();
	pvFix = new TPaveText(0.65,0.8,0.85,0.98,"NDC");
	pvFix->SetShadowColor(8);
	pvFix->AddText(Form("#Chi^{2}/NDF  %.1f/%d",fAcExpLineFixed->GetChisquare(),fAcExpLineFixed->GetNDF()));
	pvFix->AddText(Form("Prob  %f",fAcExpLineFixed->GetProb()));
	pvFix->AddText(Form("R_{0}   %.2f #pm %.2f mHz",fAcExpLineFixed->GetParameter(0)*1000.0,fAcExpLineFixed->GetParError(0)*1000.0));
	pvFix->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExpLineFixed->GetParameter(1),fAcExpLineFixed->GetParError(1)));
	pvFix->AddText(Form("R_{1}   %.4f #pm %.4f yr^{-1}",fAcExpLineFixed->GetParameter(3)*60.0*60.0*24.0*365.0,fAcExpLineFixed->GetParError(3)*60.0*60.0*24.0*365.0));
	pvFix->Draw();
	cRateNewFit->SaveAs(Form("%s/RateVsTime_Fit_New.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRateNewFit->SaveAs(Form("%s/RateVsTime_Fit_New.png",gSystem->Getenv("AD_AC227_PLOTS")));
*/
	TLine *lRateCorr = new TLine(grxStart,1.0,grxEnd,1.0);
	lRateCorr->SetLineStyle(2);
	lRateCorr->SetLineColor(kBlack);

	TGraphErrors *grRateCorr = makeCorrGr(grRate,fAcExpFixed);
	TCanvas *cRateCorr = new TCanvas("cRateCorr","Corrected rate",1000,400);
	gPad->SetGrid();
	grRateCorr->GetXaxis()->SetTitle(xLabel);	
	grRateCorr->GetYaxis()->SetTitle("R_{RnPo}/R_{Ac}");
	grRateCorr->GetXaxis()->SetTimeDisplay(1);
	grRateCorr->GetXaxis()->SetTimeFormat("%m/%d");
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
	cRateCorr->SaveAs(Form("%s/RateVsTime_Corrected.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRateCorr->SaveAs(Form("%s/RateVsTime_Corrected.png",gSystem->Getenv("AD_AC227_PLOTS")));
	

	//======================================================================================

	TCanvas *cRelPoPSDMean = new TCanvas("cRelPoPSDMean","Relative Po PSD Mean",1000,400);
	grRelPoPSDMean->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDMean->GetYaxis()->SetTitle("PSD_{Po}/#LTPSD_{Po}#GT");
	grRelPoPSDMean->GetXaxis()->SetTimeDisplay(1);
	grRelPoPSDMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPSDMean->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRelPoPSDMean->Draw("P");
	pvR = new TPaveText(0.80,0.87,0.99,0.99,"NDCNB");
	pvR->SetFillColor(19);
	pvR->AddText(Form("#LTPSD_{Po}#GT  %.3f #pm %.3f",grPoPSDMean->GetFunction("pol0")->GetParameter(0),grPoPSDMean->GetFunction("pol0")->GetParError(0)));
	pvR->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grPoPSDMean->GetFunction("pol0")->GetChisquare(),grPoPSDMean->GetFunction("pol0")->GetNDF()));	
	pvR->Draw();
	cRelPoPSDMean->SaveAs(Form("%s/RelativePoPSDMeanVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDMean->SaveAs(Form("%s/RelativePoPSDMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPSDSigma = new TCanvas("cRelPoPSDSigma","Relative Po PSD Sigma",1000,400);
	grRelPoPSDSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDSigma->GetYaxis()->SetTitle("#sigma_{PSD}/#LT#sigma_{PSD}#GT");
	grRelPoPSDSigma->GetXaxis()->SetTimeDisplay(1);
	grRelPoPSDSigma->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPSDSigma->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRelPoPSDSigma->Draw("P");
	pvR = new TPaveText(0.80,0.87,0.99,0.99,"NDCNB");
	pvR->SetFillColor(19);
	pvR->AddText(Form("#LT#sigma_{PSD}#GT  %.3f #pm %.3f",grPoPSDSigma->GetFunction("pol0")->GetParameter(0),grPoPSDSigma->GetFunction("pol0")->GetParError(0)));
	pvR->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grPoPSDSigma->GetFunction("pol0")->GetChisquare(),grPoPSDSigma->GetFunction("pol0")->GetNDF()));	
	pvR->Draw();
	cRelPoPSDSigma->SaveAs(Form("%s/RelativePoPSDSigmaVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDSigma->SaveAs(Form("%s/RelativePoPSDSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoEnMean = new TCanvas("cRelPoEnMean","Relative Po Energy Mean",1000,400);
	grRelPoEnMean->GetXaxis()->SetTitle(xLabel);
	grRelPoEnMean->GetYaxis()->SetTitle("E_{Po}/#LTE_{Po}#GT");
	grRelPoEnMean->GetXaxis()->SetTimeDisplay(1);
	grRelPoEnMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoEnMean->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRelPoEnMean->Draw("P");
	pvR = new TPaveText(0.80,0.87,0.99,0.99,"NDCNB");
	pvR->SetFillColor(19);
	pvR->AddText(Form("#LTE_{Po}#GT  %.3f #pm %.3f MeV",grPoEnMean->GetFunction("pol0")->GetParameter(0),grPoEnMean->GetFunction("pol0")->GetParError(0)));
	pvR->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grPoEnMean->GetFunction("pol0")->GetChisquare(),grPoEnMean->GetFunction("pol0")->GetNDF()));	
	pvR->Draw();
	cRelPoEnMean->SaveAs(Form("%s/RelativePoEnMeanVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnMean->SaveAs(Form("%s/RelativePoEnMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoEnSigma = new TCanvas("cRelPoEnSigma","Relative Po Energy Sigma",1000,400);
	//grRelPoEnSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoEnSigma->GetYaxis()->SetTitle("#sigma_{E}/#LT#sigma_{E}#GT");
	grRelPoEnSigma->GetXaxis()->SetTimeDisplay(1);
	grRelPoEnSigma->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grRelPoEnSigma->GetXaxis()->SetTimeOffset(0,"est");
	grRelPoEnSigma->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRelPoEnSigma->Draw("P");
	pvR = new TPaveText(0.80,0.87,0.99,0.99,"NDCNB");
	pvR->SetFillColor(19);
	pvR->AddText(Form("#LT#sigma_{E}#GT  %.3f #pm %.3f MeV",grPoEnSigma->GetFunction("pol0")->GetParameter(0),grPoEnSigma->GetFunction("pol0")->GetParError(0)));
	pvR->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grPoEnSigma->GetFunction("pol0")->GetChisquare(),grPoEnSigma->GetFunction("pol0")->GetNDF()));	
	pvR->Draw();
	cRelPoEnSigma->SaveAs(Form("%s/RelativePoEnSigmaVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnSigma->SaveAs(Form("%s/RelativePoEnSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPosMean = new TCanvas("cRelPoPosMean","Relative Po Position Mean",1000,400);
	grRelPoPosMean->GetXaxis()->SetTitle(xLabel);
	grRelPoPosMean->GetYaxis()->SetTitle("z_{Po}/#LTz_{Po}#GT");
	grRelPoPosMean->GetXaxis()->SetTimeDisplay(1);
	grRelPoPosMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPosMean->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRelPoPosMean->Draw("P");
	pvR = new TPaveText(0.80,0.87,0.99,0.99,"NDCNB");
	pvR->SetFillColor(19);
	pvR->AddText(Form("#LTz_{Po}#GT  %.3f #pm %.3f mm",grPoPosMean->GetFunction("pol0")->GetParameter(0),grPoPosMean->GetFunction("pol0")->GetParError(0)));
	pvR->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grPoPosMean->GetFunction("pol0")->GetChisquare(),grPoPosMean->GetFunction("pol0")->GetNDF()));	
	pvR->Draw();
	cRelPoPosMean->SaveAs(Form("%s/RelativePoPosMeanVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosMean->SaveAs(Form("%s/RelativePoPosMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPosSigma = new TCanvas("cRelPoPosSigma","Relative Po Position Sigma",1000,400);
	grRelPoPosSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoPosSigma->GetYaxis()->SetTitle("RMS_{z}/#LTRMS_{z}#GT");
	grRelPoPosSigma->GetXaxis()->SetTimeDisplay(1);
	grRelPoPosSigma->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPosSigma->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRelPoPosSigma->Draw("P");
	pvR = new TPaveText(0.80,0.87,0.99,0.99,"NDCNB");
	pvR->SetFillColor(19);
	pvR->AddText(Form("#LTRMS_{z}#GT  %.3f #pm %.3f mm",grPoPosSigma->GetFunction("pol0")->GetParameter(0),grPoPosSigma->GetFunction("pol0")->GetParError(0)));
	pvR->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grPoPosSigma->GetFunction("pol0")->GetChisquare(),grPoPosSigma->GetFunction("pol0")->GetNDF()));	
	pvR->Draw();
	cRelPoPosSigma->SaveAs(Form("%s/RelativePoPosSigmaVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosSigma->SaveAs(Form("%s/RelativePoPosSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRnPoDzMean = new TCanvas("cRelRnPoDzMean","Relative RnPo Dz Mean",1000,400);
	grRelRnPoDzMean->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzMean->GetYaxis()->SetTitle("dz_{RnPo}/#LTdz_{RnPo}#GT");
	grRelRnPoDzMean->GetXaxis()->SetTimeDisplay(1);
	grRelRnPoDzMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelRnPoDzMean->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRelRnPoDzMean->Draw("P");
	pvR = new TPaveText(0.80,0.87,0.99,0.99,"NDCNB");
	pvR->SetFillColor(19);
	pvR->AddText(Form("#LTdz_{RnPo}#GT  %.3f #pm %.3f mm",grRnPoDzMean->GetFunction("pol0")->GetParameter(0),grRnPoDzMean->GetFunction("pol0")->GetParError(0)));
	pvR->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grRnPoDzMean->GetFunction("pol0")->GetChisquare(),grRnPoDzMean->GetFunction("pol0")->GetNDF()));	
	pvR->Draw();
	cRelRnPoDzMean->SaveAs(Form("%s/RelativeRnPoDzMeanVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzMean->SaveAs(Form("%s/RelativeRnPoDzMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRnPoDzSigma = new TCanvas("cRelRnPoDzSigma","Relative RnPo Dz Sigma",1000,400);
	//grRelRnPoDzSigma->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzSigma->GetYaxis()->SetTitle("#sigma_{dz}/#LT#sigma_{dz}#GT");
	grRelRnPoDzSigma->GetXaxis()->SetTimeDisplay(1);
	grRelRnPoDzSigma->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grRelRnPoDzSigma->GetXaxis()->SetTimeOffset(0,"est");
	grRelRnPoDzSigma->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grRelRnPoDzSigma->Draw("P");
	pvR = new TPaveText(0.80,0.87,0.99,0.99,"NDCNB");
	pvR->SetFillColor(19);
	pvR->AddText(Form("#LT#sigma_{dz}#GT  %.3f #pm %.3f mm",grRnPoDzSigma->GetFunction("pol0")->GetParameter(0),grRnPoDzSigma->GetFunction("pol0")->GetParError(0)));
	pvR->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grRnPoDzSigma->GetFunction("pol0")->GetChisquare(),grRnPoDzSigma->GetFunction("pol0")->GetNDF()));	
	pvR->Draw();
	cRelRnPoDzSigma->SaveAs(Form("%s/RelativeRnPoDzSigmaVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzSigma->SaveAs(Form("%s/RelativeRnPoDzSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cEff = new TCanvas("cEff","Efficiency",1000,400);
//	grTotEff->GetXaxis()->SetTitle(xLabel);
	grTotEff->GetYaxis()->SetTitle("Efficiency");
	grTotEff->GetXaxis()->SetTimeDisplay(1);
	grTotEff->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grTotEff->GetXaxis()->SetTimeOffset(0,"est");
	grTotEff->GetYaxis()->SetRangeUser(0.94,1.0002);
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

	TLegend *leg = new TLegend(0.90,0.67,0.99,0.99);
	leg->AddEntry(grPromptPSDEff,"Prompt PSD","p");
	leg->AddEntry(grDelayPSDEff,"Delay PSD","p");
	leg->AddEntry(grPromptEnEff,"Prompt E","p");
	leg->AddEntry(grDelayEnEff,"Delay E","p");
	leg->AddEntry(grDzEff,"Dz","p");
	leg->AddEntry(grTotEff,"Total","p");
	leg->Draw();	
	cEff->SaveAs(Form("%s/EfficiencyVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cEff->SaveAs(Form("%s/EfficiencyVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

/*
	TCanvas *cEnEff = new TCanvas("cEnEff","En Eff",1000,400);
	gPad->SetGrid();
	grDelayEnEff->GetYaxis()->SetTitle("Efficiency");
	grDelayEnEff->GetXaxis()->SetTimeDisplay(1);
	grDelayEnEff->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grDelayEnEff->GetXaxis()->SetTimeOffset(0,"est");
	grDelayEnEff->GetYaxis()->SetRangeUser(0.965,1.0002);
	grDelayEnEff->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grDelayEnEff->Draw("P");
	grPromptEnEff->Draw("P");
	leg = new TLegend(0.90,0.80,0.99,0.99);
	leg->AddEntry(grPromptEnEff,"Rn E","p");
	leg->AddEntry(grDelayEnEff,"Po E","p");
	leg->Draw();	
	cEnEff->SaveAs(Form("%s/EnergyEfficiencyVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cEnEff->SaveAs(Form("%s/EnergyEfficiencyVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cDzEff = new TCanvas("cDzEff","dz eff",1000,400);
	gPad->SetGrid();
	grDzEff->GetYaxis()->SetTitle("Efficiency");
	grDzEff->GetXaxis()->SetTimeDisplay(1);
	grDzEff->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grDzEff->GetXaxis()->SetTimeOffset(0,"est");
	grDzEff->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grDzEff->Draw("P");
	cDzEff->SaveAs(Form("%s/DzEfficiencyVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cDzEff->SaveAs(Form("%s/DzEfficiencyVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));
*/


	TCanvas *cVetoTime = new TCanvas("cVetoTime","VetoTime",1000,400);
	gPad->SetGrid();
	grPileupVetoFrac->GetXaxis()->SetTitle(xLabel);
	grPileupVetoFrac->GetYaxis()->SetTitle("t_{veto}/t_{live}");
	grPileupVetoFrac->GetXaxis()->SetTimeDisplay(1);
	grPileupVetoFrac->GetXaxis()->SetTimeFormat("%m/%d");
	grPileupVetoFrac->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grPileupVetoFrac->Draw("P");
	cVetoTime->SaveAs(Form("%s/VetoTimeVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cVetoTime->SaveAs(Form("%s/VetoTimeVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));
	
	
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
	grRnPSDChiSq->GetXaxis()->SetTimeFormat("%m/%d");
	grRnPSDChiSq->GetYaxis()->SetRangeUser(0,4);
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
	cFitChiSq->SaveAs(Form("%s/FitChiSqVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cFitChiSq->SaveAs(Form("%s/FitChiSqVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));
	
	TCanvas *cBGRate = new TCanvas("cBGRate","BG Rate",1000,400);
	gPad->SetGrid();
	//grBGRate->GetXaxis()->SetTitle(xLabel);
	grBGRate->GetYaxis()->SetTitle("BG Rate [Hz]");
	grBGRate->GetXaxis()->SetTimeDisplay(1);
	grBGRate->GetXaxis()->SetTimeFormat("%m/%d/%Y");
	grBGRate->GetXaxis()->SetTimeOffset(0,"est");
	//grBGRate->GetYaxis()->SetRangeUser(0,62e-3);
	grBGRate->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grBGRate->Draw("P");
	cBGRate->SaveAs(Form("%s/BGRateVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cBGRate->SaveAs(Form("%s/BGRateVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

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
	
	TCanvas *cLifetime = new TCanvas("cLifetime","BG Rate",1000,400);
	gPad->SetGrid();
	grLifetime->GetXaxis()->SetTitle(xLabel);
	grLifetime->GetYaxis()->SetTitle("#tau_{Po} [ms]");
	grLifetime->GetXaxis()->SetTimeDisplay(1);
	grLifetime->GetXaxis()->SetTimeFormat("%m/%d");
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
	pv->AddText(Form("Prob   %f",grLifetime->GetFunction("pol0")->GetProb()));
	pv->AddText(Form("p0   %.3f #pm %.3f",grLifetime->GetFunction("pol0")->GetParameter(0),grLifetime->GetFunction("pol0")->GetParError(0)));
	pv->Draw();
	cLifetime->SaveAs(Form("%s/LifetimeVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cLifetime->SaveAs(Form("%s/LifetimeVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cLivetime = new TCanvas("cLivetime","livetime",1000,400);
	gPad->SetGrid();
	grLivetime->GetXaxis()->SetTitle(xLabel);
	grLivetime->GetYaxis()->SetTitle("Livetime [ms]");
	grLivetime->GetXaxis()->SetTimeDisplay(1);
	grLivetime->GetXaxis()->SetTimeFormat("%m/%d");
	grLivetime->Draw("AP");
	grMarchOn_2018->Draw("f");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grSeptOn_2018->Draw("f");
	grLivetime->Draw("P");
	cLivetime->SaveAs(Form("%s/LivetimeVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cLivetime->SaveAs(Form("%s/LivetimeVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));


} 	//end PlotRnPoVsTime
