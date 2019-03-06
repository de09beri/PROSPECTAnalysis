//This macro will create histograms for Ac227 coincidences
//according to cell number

#include "RNPO.C"

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TMath.h"
#include "TLatex.h"
#include "TVectorD.h"

#include "Header.C"

void RnPoVsTime_Calc(){
	//---------------------------------------------------------------------------------
	TFile *histFile = new TFile(Form("%s/Ac227_HistsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));

	vector<double> *vT;
	histFile->GetObject("vTimestamp",vT);

	vector<double> *vL;
	histFile->GetObject("vLivetime",vL);

	//---------------------------------------------------------------------------------
	TH1F *hRnPoDt;
	TH1F *hRnPSD,	*hPoPSD;
	TH1F *hRnEn, 	*hPoEn,	 *hRnTotEn, *hPoEnSmear;
	TH1F *hRnPos,	*hPoPos;
	TH1F *hRnPoDz;

	TH2F *hRnPoPSDvsEn;
	TH2F *hPoEnvsRnEn;

	//---------------------------------------------------------------------------------
	TF1 *fRnPoDtExp;
	TF1 *fRnPSDGaus, *fPoPSDGaus;
	TF1 *fRnEnGaus,  *fPoEnGaus;
	TF1 *fRnPoDzGaus;

	//---------------------------------------------------------------------------------
	vector<double> vRate,   vRateErr;
	vector<double> vTotEff, vTotEffErr;
	vector<double> vLifetime,   vLifetimeErr;
	vector<double> vPoPSDMean,  vPoPSDMeanErr,  vPoPSDSigma,  vPoPSDSigmaErr;
	vector<double> vPoEnMean,   vPoEnMeanErr,   vPoEnSigma,   vPoEnSigmaErr;
	vector<double> vPoPosMean,  vPoPosMeanErr,  vPoPosSigma,  vPoPosSigmaErr;
	vector<double> vRnPoDzMean, vRnPoDzMeanErr, vRnPoDzSigma, vRnPoDzSigmaErr;

	vector<double> vTotLivetime,  vPileupVetoT;
	vector<double> vPileupVetoFrac;
	
	vector<double> vPromptEnEff,  vPromptEnEffErr,  vDelayEnEff,  vDelayEnEffErr;
	vector<double> vPromptPSDEff, vPromptPSDEffErr, vDelayPSDEff, vDelayPSDEffErr; 
	vector<double> vDzEff, 	      vDzEffErr;

	vector<double> vRnPSDChiSq, vPoPSDChiSq, vRnEnChiSq, vPoEnChiSq, vDzChiSq, vDtChiSq;

	vector<double> vBGRate;
	
	//---------------------------------------------------------------------------------
	//Calculate results
	printf("=============== Calculating Results =============== \n"); 

	double dtBinWidth = (dtMax - dtMin)/(double)numDtBins;

	double promptPSDEff, delayPSDEff, promptPSDEffErr, delayPSDEffErr;
	double promptEnEff,  delayEnEff,  promptEnEffErr,  delayEnEffErr;
	double dzEff, dzEffErr;
	double totEff, totEffErr;

	double NAlpha, NAlphaErr, lifetime, lifetimeErr;
	double rate, rateErr;		

	//int numTimeBins = vLivetime.size();
	int numTimeBins = vL->size();

	for(int i=0;i<numTimeBins;i++){

		hRnPoDt = (TH1F*)histFile->Get(Form("hRnPoDt_%i",i));
		hRnPSD  = (TH1F*)histFile->Get(Form("hRnPSD_%i",i));
		hPoPSD  = (TH1F*)histFile->Get(Form("hPoPSD_%i",i));
		hRnEn   = (TH1F*)histFile->Get(Form("hRnEn_%i",i));
		hPoEn   = (TH1F*)histFile->Get(Form("hPoEn_%i",i));
		hRnPoDz = (TH1F*)histFile->Get(Form("hRnPoDz_%i",i));

		hPoPos   = (TH1F*)histFile->Get(Form("hPoPos_%i",i));

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Fit distributions
/*
		TF1 *fPoExp = new TF1("fPoExp",Form("[0]*exp(-x/[1])*(%f/[1])",dtBinWidth),2,4);
		fPoExp->SetParameter(1,POLIFETIME);
		hRnPoDt->Fit(fPoExp,"R0");

		fRnPoDtExp = new TF1("fRnPoDtExp",Form("[0]*exp(-x/[1])*(%f/[1]) + [2]*exp(-x/[3])",dtBinWidth),0,8);
		fRnPoDtExp->SetParameters(5e4,POLIFETIME,3.5e3,0.06);
		fRnPoDtExp->FixParameter(1,fPoExp->GetParameter(1));
		fRnPoDtExp->FixParameter(0,fPoExp->GetParameter(0));
		hRnPoDt->Fit(fRnPoDtExp,"R");
*/

		fRnPoDtExp = new TF1("fRnPoDtExp","[0]*exp(-[1]*x) + [2]*exp(-[3]*x)",0,6);
		fRnPoDtExp->SetParameters(1450,0.385,1950,21);
		hRnPoDt->Fit(fRnPoDtExp,"R");


/*		double fitPSDMin = promptLowPSDCut;
		double maxValue = hRnPSD->GetMaximum(), maxBin = hRnPSD->GetBinCenter(hRnPSD->GetMaximumBin());	
		fRnPSDGaus = new TF1("fRnPSDGaus","gaus(0)+gaus(3)",fitPSDMin,PSDMax);
		fRnPSDGaus->SetParameters(0.25*maxValue,maxBin+0.01,0.015,0.75*maxValue,maxBin-0.001,0.02);
		fRnPSDGaus->SetParLimits(0,0,maxValue);
		fRnPSDGaus->SetParLimits(3,0,maxValue);
		hRnPSD->Fit(fRnPSDGaus,"RQ0");
		fRnPSDGaus->SetRange(PSDMin,PSDMax);
*/
		//fitPSDMin = delayLowPSDCut;
		fPoPSDGaus = new TF1("fPoPSDGaus","gaus",PSDMin,PSDMax);
		hPoPSD->Fit(fPoPSDGaus,"RQ0");
		fPoPSDGaus->SetRange(PSDMin,PSDMax);
	
/*		double fitRnEnMin = EnMin;	
		fRnEnGaus = new TF1("fRnEnGaus","gaus(0)+gaus(3)+gaus(6)",fitRnEnMin,EnMax);
		maxValue = hRnEn->GetMaximum();
		maxBin = hRnEn->GetBinCenter(hRnEn->GetMaximumBin());	
		fRnEnGaus->SetParameters(0.8*maxValue,maxBin,0.035,0.05*maxValue,maxBin+0.14,0.08,0.15*maxValue,maxBin-0.03,0.05);
		hRnEn->Fit(fRnEnGaus,"RQ0");
		fRnEnGaus->SetRange(EnMin,EnMax);
*/
		//double fitPoEnMin = delayLowEnCut;
		fPoEnGaus = new TF1("fPoEnGaus","gaus",EnMin,EnMax);
		hPoEn->Fit(fPoEnGaus,"RQ0");
		fPoEnGaus->SetRange(EnMin,EnMax);

		fRnPoDzGaus = new TF1("fRnPoDzGaus","gaus",dzMin,dzMax);
		hRnPoDz->Fit(fRnPoDzGaus,"RQ0");
		fRnPoDzGaus->SetRange(dzMin,dzMax);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Calculate efficiencies
/*
		promptPSDEff = fRnPSDGaus->Integral(promptLowPSDCut,promptHighPSDCut)/fRnPSDGaus->Integral(PSDMin,PSDMax);
		promptPSDEffErr = sqrt((promptPSDEff*(1-promptPSDEff))/hRnPSD->GetEntries()); 
			
		delayPSDEff = fPoPSDGaus->Integral(delayLowPSDCut,delayHighPSDCut)/fPoPSDGaus->Integral(PSDMin,PSDMax);
		delayPSDEffErr = sqrt((delayPSDEff*(1-delayPSDEff))/hPoPSD->GetEntries()); 
		
		promptEnEff = fRnEnGaus->Integral(promptLowEnCut,promptHighEnCut)/fRnEnGaus->Integral(EnMin,EnMax);
		promptEnEffErr = sqrt((promptEnEff*(1-promptEnEff))/hRnEn->GetEntries());

		delayEnEff = fPoEnGaus->Integral(delayLowEnCut,delayHighEnCut)/fPoEnGaus->Integral(EnMin,EnMax);
		delayEnEffErr = sqrt((delayEnEff*(1-delayEnEff))/hPoEn->GetEntries()); 

		dzEff = fRnPoDzGaus->Integral(dzCutLow,dzCutHigh)/fRnPoDzGaus->Integral(dzMin,dzMax);
		dzEffErr = sqrt((dzEff*(1-dzEff))/hRnPoDz->GetEntries());
*/
		promptPSDEff = 1.0;
		promptPSDEffErr = 0.0;
		delayPSDEff = 1.0;
		delayPSDEffErr = 0.0;
		promptEnEff = 1.0;
		promptEnEffErr = 0.0;
		delayEnEff = 1.0;
		delayEnEffErr = 0.0;
		dzEff = 1.0;
		dzEffErr = 0.0;	
	
//		totEff = promptPSDEff * delayPSDEff * promptEnEff * delayEnEff * dzEff;	
//		totEffErr = totEff * sqrt( pow(promptPSDEffErr/promptPSDEff,2) + pow(delayPSDEffErr/delayPSDEff,2) + pow(promptEnEffErr/promptEnEff,2) + pow(delayEnEffErr/delayEnEff,2) + pow(dzEffErr/dzEff,2) );

		totEff = 1.0;
		totEffErr = 0.0;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Calculate rate
		double livetime = vL->at(i);

		NAlpha = (fRnPoDtExp->GetParameter(0)*fRnPoDtExp->GetParameter(1))/dtBinWidth;
		NAlphaErr = NAlpha*(1/dtBinWidth)*(pow(fRnPoDtExp->GetParError(0)/fRnPoDtExp->GetParameter(0),2) + pow(fRnPoDtExp->GetParError(1)/fRnPoDtExp->GetParameter(1),2));
		lifetime = log(2)/(double)fRnPoDtExp->GetParameter(1);
		lifetimeErr = log(2)/(double)fRnPoDtExp->GetParError(1);
	
		rate = (NAlpha/(livetime*totEff))*(1e3);		//Hz
		rateErr = rate * sqrt(pow(NAlphaErr/NAlpha,2) + pow(totEffErr/totEff,2));

		printf("Rate: %.4f +/- %.4f \n",rate,rateErr);

		//---------------------------------------------------------------------------------
		//Populate vectors
		vRate.push_back(rate);
		vRateErr.push_back(rateErr);
		vTotEff.push_back(totEff);
		vTotEffErr.push_back(totEffErr);
		vLifetime.push_back(lifetime);
		vLifetimeErr.push_back(lifetimeErr);
		
		vPoPSDMean.push_back(fPoPSDGaus->GetParameter(1));	
		vPoPSDMeanErr.push_back(fPoPSDGaus->GetParError(1));
		vPoPSDSigma.push_back(fPoPSDGaus->GetParameter(2));
		vPoPSDSigmaErr.push_back(fPoPSDGaus->GetParError(2));

		vPoEnMean.push_back(fPoEnGaus->GetParameter(1));
		vPoEnMeanErr.push_back(fPoEnGaus->GetParError(1));
		vPoEnSigma.push_back(fPoEnGaus->GetParameter(2));
		vPoEnSigmaErr.push_back(fPoEnGaus->GetParError(2));

		vPoPosMean.push_back(hPoPos->GetMean());
		vPoPosMeanErr.push_back(hPoPos->GetMeanError());
		vPoPosSigma.push_back(hPoPos->GetRMS());
		vPoPosSigmaErr.push_back(hPoPos->GetRMSError());

		vRnPoDzMean.push_back(fRnPoDzGaus->GetParameter(1));
		vRnPoDzMeanErr.push_back(fRnPoDzGaus->GetParError(1));
		vRnPoDzSigma.push_back(fRnPoDzGaus->GetParameter(2));
		vRnPoDzSigmaErr.push_back(fRnPoDzGaus->GetParError(2));

		vPromptEnEff.push_back(promptEnEff);
		vPromptEnEffErr.push_back(promptEnEffErr);
		vDelayEnEff.push_back(delayEnEff);
		vDelayEnEffErr.push_back(delayEnEffErr);

		vPromptPSDEff.push_back(promptPSDEff);
		vPromptPSDEffErr.push_back(promptPSDEffErr);
		vDelayPSDEff.push_back(delayPSDEff);
		vDelayPSDEffErr.push_back(delayPSDEffErr);

		vDzEff.push_back(dzEff);
		vDzEffErr.push_back(dzEffErr);			

		vDtChiSq.push_back(fRnPoDtExp->GetChisquare()/(double)fRnPoDtExp->GetNDF());

/*		vRnPSDChiSq.push_back(fRnPSDGaus->GetChisquare()/(double)fRnPSDGaus->GetNDF());
		vPoPSDChiSq.push_back(fPoPSDGaus->GetChisquare()/(double)fPoPSDGaus->GetNDF());
	
		vRnEnChiSq.push_back(fRnEnGaus->GetChisquare()/(double)fRnEnGaus->GetNDF());
		vPoEnChiSq.push_back(fPoEnGaus->GetChisquare()/(double)fPoEnGaus->GetNDF());
		
		vDzChiSq.push_back(fRnPoDzGaus->GetChisquare()/(double)fRnPoDzGaus->GetNDF());	
	
		vDtChiSq.push_back(fRnPoDtExp->GetChisquare()/(double)fRnPoDtExp->GetNDF());
*/
	}	//end while loop 

	histFile->Close();

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Initialize TGraphs for results

	int numPt = vRate.size();
	double x[numPt], xErr[numPt];
	double y[1], yErr[1];

	TGraphErrors *grRate 		= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grTotEff 		= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grLifetime 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPSDMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPSDSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoEnMean   	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoEnSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPosMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPosSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);

//	TGraph *grLivetime 	 = new TGraph(numPt,x,y);
//	TGraph *grTotLivetime 	 = new TGraph(numPt,x,y);
//	TGraph *grPileupVeto	 = new TGraph(numPt,x,y);
//	TGraph *grPileupVetoFrac = new TGraph(numPt,x,y);

	TGraphErrors *grPromptEnEff  = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDelayEnEff   = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPromptPSDEff = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDelayPSDEff  = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDzEff 	     = new TGraphErrors(numPt,x,y,xErr,yErr);

	TGraph *grDtChiSq    = new TGraph(numPt,x,y);

	TH1F *hDtChiSq = new TH1F("hDtChiSq","dt chisq",20,0.5,1.5);

/*	TGraph *grRnPSDChiSq = new TGraph(numPt,x,y);
	TGraph *grPoPSDChiSq = new TGraph(numPt,x,y);
	TGraph *grRnEnChiSq  = new TGraph(numPt,x,y);
	TGraph *grPoEnChiSq  = new TGraph(numPt,x,y);
	TGraph *grDzChiSq    = new TGraph(numPt,x,y);
	TGraph *grDtChiSq    = new TGraph(numPt,x,y);

	TGraph *grBGRate = new TGraph(numPt,x,y);
*/
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Fill TGraphs
	for(int i=0;i<numPt;i++){
		//double time = vTimestamp[i];
		double time = vT->at(i);

		grRate->SetPoint(i,time,vRate[i]);
		grRate->SetPointError(i,0,vRateErr[i]);

		grTotEff->SetPoint(i,time,vTotEff[i]);
		grTotEff->SetPointError(i,0,vTotEffErr[i]);

		grLifetime->SetPoint(i,time,vLifetime[i]);
		grLifetime->SetPointError(i,0,vLifetimeErr[i]);

		grPoPSDMean->SetPoint(i,time,vPoPSDMean[i]);
		grPoPSDMean->SetPointError(i,0,vPoPSDMeanErr[i]);

		grPoPSDSigma->SetPoint(i,time,vPoPSDSigma[i]);
		grPoPSDSigma->SetPointError(i,0,vPoPSDSigmaErr[i]);

		grPoEnMean->SetPoint(i,time,vPoEnMean[i]);
		grPoEnMean->SetPointError(i,0,vPoEnMeanErr[i]);
	
		grPoEnSigma->SetPoint(i,time,vPoEnSigma[i]);
		grPoEnSigma->SetPointError(i,0,vPoEnSigmaErr[i]);

		grPoPosMean->SetPoint(i,time,vPoPosMean[i]);
		grPoPosMean->SetPointError(i,0,vPoPosMeanErr[i]);

		grPoPosSigma->SetPoint(i,time,vPoPosSigma[i]);
		grPoPosSigma->SetPointError(i,0,vPoPosSigmaErr[i]);		

		grRnPoDzMean->SetPoint(i,time,vRnPoDzMean[i]);
		grRnPoDzMean->SetPointError(i,0,vRnPoDzMeanErr[i]);

		grRnPoDzSigma->SetPoint(i,time,vRnPoDzSigma[i]);
		grRnPoDzSigma->SetPointError(i,0,vRnPoDzSigmaErr[i]);


//		grLivetime->SetPoint(i,time,vLivetime[i]);
//		grTotLivetime->SetPoint(i,time,vTotLivetime[i]);
//		grPileupVeto->SetPoint(i,time,vPileupVetoT[i]);
//		grPileupVetoFrac->SetPoint(i,time,vPileupVetoFrac[i]);

		grPromptEnEff->SetPoint(i,time,vPromptEnEff[i]);
		grPromptEnEff->SetPointError(i,0,vPromptEnEffErr[i]);

		grDelayEnEff->SetPoint(i,time,vDelayEnEff[i]);
		grDelayEnEff->SetPointError(i,0,vDelayEnEffErr[i]);

		grPromptPSDEff->SetPoint(i,time,vPromptPSDEff[i]);
		grPromptPSDEff->SetPointError(i,0,vPromptPSDEffErr[i]);

		grDelayPSDEff->SetPoint(i,time,vDelayPSDEff[i]);
		grDelayPSDEff->SetPointError(i,0,vDelayPSDEffErr[i]);

		grDzEff->SetPoint(i,time,vDzEff[i]);
		grDzEff->SetPointError(i,0,vDzEffErr[i]);	

		grDtChiSq->SetPoint(i,time,vDtChiSq[i]);
		hDtChiSq->Fill(vDtChiSq[i]);

/*		
		grRnPSDChiSq->SetPoint(i,time,vRnPSDChiSq[i]);
		grPoPSDChiSq->SetPoint(i,time,vPoPSDChiSq[i]);
		grRnEnChiSq->SetPoint(i,time,vRnEnChiSq[i]);
		grPoEnChiSq->SetPoint(i,time,vPoEnChiSq[i]);
		grDzChiSq->SetPoint(i,time,vDzChiSq[i]);
		grDtChiSq->SetPoint(i,time,vDtChiSq[i]);

		grBGRate->SetPoint(i,time,vBGRate[i]);
*/
	}	//end for loop to populate TGraphs

	//---------------------------------------------------------------------------------
	//Write TGraphs to file
	TFile *graphFile = new TFile(Form("%s/Ac227_GraphsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");
	
	grRate->Write("grRate");
	grTotEff->Write("grTotEff");	
	grLifetime->Write("grLifetime");
	grPoPSDMean->Write("grPoPSDMean");
	grPoPSDSigma->Write("grPoPSDSigma");
	grPoEnMean->Write("grPoEnMean");
	grPoEnSigma->Write("grPoEnSigma");
	grPoPosMean->Write("grPoPosMean");
	grPoPosSigma->Write("grPoPosSigma");
	grRnPoDzMean->Write("grRnPoDzMean");
	grRnPoDzSigma->Write("grRnPoDzSigma");	
//	grLivetime->Write("grLivetime");
//	grTotLivetime->Write("grTotLivetime");
//	grPileupVeto->Write("grPileupVeto");
//	grPileupVetoFrac->Write("grPileupVetoFrac");
	grPromptEnEff->Write("grPromptEnEff");
	grDelayEnEff->Write("grDelayEnEff");
	grPromptPSDEff->Write("grPromptPSDEff");
	grDelayPSDEff->Write("grDelayPSDEff");
	grDzEff->Write("grDzEff");
	grDtChiSq->Write("grDtChiSq");
	hDtChiSq->Write();
/*
	grRnPSDChiSq->Write("grRnPSDChiSq");
	grPoPSDChiSq->Write("grPoPSDChiSq");
	grRnEnChiSq->Write("grRnEnChiSq");
	grPoEnChiSq->Write("grPoEnChiSq");
	grDzChiSq->Write("grDzChiSq");
	grDtChiSq->Write("grDtChiSq");
	grBGRate->Write("grBGRate");
*/
	graphFile->Close();

}	//end void RnPoVsTime

