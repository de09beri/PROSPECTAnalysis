//This macro will create histograms for Ac227 coincidences
//according to cell number

#include "TFile.h"
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
#include "TFitResult.h"
#include "TMatrixD.h"

#include "Header.C"

void RnPoVsTime_Calc(){
	//---------------------------------------------------------------------------------
	TFile *histFile = new TFile(Form("%s/Ac227_HistsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));

	vector<double> *vT;
	histFile->GetObject("vTimestamp",vT);

	vector<double> *vL;
	histFile->GetObject("vLivetime",vL);

	vector<double> *vTotLivetime;
	histFile->GetObject("vTotLivetime",vTotLivetime);	//minutes
	vector<double> *vPileupVetoTime;
	histFile->GetObject("vPileupVetoTime",vPileupVetoTime);	//minutes
	vector<double> *vPileupVetoFrac;
	histFile->GetObject("vPileupVetoFrac",vPileupVetoFrac);

	vector<double> *vRnPSDCutLow;
	histFile->GetObject("vRnPSDCutLow",vRnPSDCutLow);
	vector<double> *vRnPSDCutHigh;
	histFile->GetObject("vRnPSDCutHigh",vRnPSDCutHigh);

	vector<double> *vPoPSDCutLow;
	histFile->GetObject("vPoPSDCutLow",vPoPSDCutLow);
	vector<double> *vPoPSDCutHigh;
	histFile->GetObject("vPoPSDCutHigh",vPoPSDCutHigh);

	vector<double> *vRnEnCutLow;
	histFile->GetObject("vRnEnCutLow",vRnEnCutLow);
	vector<double> *vRnEnCutHigh;
	histFile->GetObject("vRnEnCutHigh",vRnEnCutHigh);

	vector<double> *vPoEnCutLow;
	histFile->GetObject("vPoEnCutLow",vPoEnCutLow);
	vector<double> *vPoEnCutHigh;
	histFile->GetObject("vPoEnCutHigh",vPoEnCutHigh);
	
	vector<double> *vRnPoDzCutLow;
	histFile->GetObject("vRnPoDzCutLow",vRnPoDzCutLow);
	vector<double> *vRnPoDzCutHigh;
	histFile->GetObject("vRnPoDzCutHigh",vRnPoDzCutHigh);

	//---------------------------------------------------------------------------------
	TH1F *hBGDt;
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

	vector<double> vPromptEnEff,  vPromptEnEffErr,  vDelayEnEff,  vDelayEnEffErr;
	vector<double> vPromptPSDEff, vPromptPSDEffErr, vDelayPSDEff, vDelayPSDEffErr; 
	vector<double> vDzEff, 	      vDzEffErr;

	vector<double> vRnPSDChiSq, vPoPSDChiSq, vRnEnChiSq, vPoEnChiSq, vDzChiSq, vDtChiSq;

	vector<double> vBGRate;
	
	//---------------------------------------------------------------------------------
	//Calculate results
	printf("=============== Calculating Results =============== \n"); 

	double dtBinWidth = (dtMax - dtMin)/(double)numDtBins;

	double promptLowPSDCut, promptHighPSDCut, promptLowEnCut, promptHighEnCut;
	double delayLowPSDCut, delayHighPSDCut, delayLowEnCut, delayHighEnCut;
	double dzCutLow, dzCutHigh;

	double promptPSDEff, delayPSDEff, promptPSDEffErr, delayPSDEffErr;
	double promptEnEff,  delayEnEff,  promptEnEffErr,  delayEnEffErr;
	double dzEff, dzEffErr;
	double totEff, totEffErr;

	double N0, N0Err, lifetime, lifetimeErr;
	double rate, rateErr;		
	double BGRate; 

	//int numTimeBins = vLivetime.size();
	int numTimeBins = vL->size();

	for(int i=0;i<numTimeBins;i++){

		hBGDt   = (TH1F*)histFile->Get(Form("hBGDt_%i",i));
		hRnPoDt = (TH1F*)histFile->Get(Form("hRnPoDt_%i",i));
		hRnPSD  = (TH1F*)histFile->Get(Form("hRnPSD_%i",i));
		hPoPSD  = (TH1F*)histFile->Get(Form("hPoPSD_%i",i));
		hRnEn   = (TH1F*)histFile->Get(Form("hRnEn_%i",i));
		hPoEn   = (TH1F*)histFile->Get(Form("hPoEn_%i",i));
		hRnPoDz = (TH1F*)histFile->Get(Form("hRnPoDz_%i",i));

		hPoPos   = (TH1F*)histFile->Get(Form("hPoPos_%i",i));

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Fit distributions
		printf("========= Time Bin %d ========= \n",i);
		printf("Fitting distributions \n");

		fRnPoDtExp = new TF1("fRnPoDtExp","[0]*exp(-x/[1])",0.5,8);
                fRnPoDtExp->SetParameter(1,POLIFETIME);
                TFitResultPtr fitr = hRnPoDt->Fit(fRnPoDtExp,"RS0Q");
                TMatrixD DtCov = fitr->GetCovarianceMatrix();

		//-----------------
		double fitSigma = 3.0;
                double fitMin, fitMax;
		
		//-----------------
		fitMin = hRnPSD->GetMean() - fitSigma*hRnPSD->GetStdDev();
                fitMax = hRnPSD->GetMean() + fitSigma*hRnPSD->GetStdDev();

                fRnPSDGaus = new TF1("fRnPSDGaus","gaus",fitMin,fitMax);
                hRnPSD->Fit(fRnPSDGaus,"R0LQ");
                fRnPSDGaus->SetRange(PSDMin,PSDMax);

		//-----------------
		fitMin = hPoPSD->GetMean() - fitSigma*hPoPSD->GetStdDev();
                fitMax = hPoPSD->GetMean() + fitSigma*hPoPSD->GetStdDev();

                fPoPSDGaus = new TF1("fPoPSDGaus","gaus",fitMin,fitMax);
                hPoPSD->Fit(fPoPSDGaus,"R0LQ");
                fPoPSDGaus->SetRange(PSDMin,PSDMax);
	
		//-----------------
		double maxValue = hRnEn->GetMaximum();
                double maxBinX = hRnEn->GetBinCenter(hRnEn->GetMaximumBin());

                TF1 *fRnMainGaus = new TF1("fRnMainGaus","gaus",maxBinX-0.05,maxBinX+0.05);
                hRnEn->Fit(fRnMainGaus,"R0LQ");

                TF1 *fRnGammaGaus = new TF1("fRnGammaGaus","gaus",maxBinX+0.15,maxBinX+0.25);
                hRnEn->Fit(fRnGammaGaus,"R0LQ");

                double RnGausPar[6];
                fRnMainGaus->GetParameters(&RnGausPar[0]);
                fRnGammaGaus->GetParameters(&RnGausPar[3]);

                fitMin = fRnMainGaus->GetParameter(1) - 2.5*fRnMainGaus->GetParameter(2);
                fitMax = fRnGammaGaus->GetParameter(1) + 2.5*fRnGammaGaus->GetParameter(2);

                fRnEnGaus = new TF1("fRnEnGaus","gaus(0)+gaus(3)",fitMin,fitMax);
                fRnEnGaus->SetParameters(RnGausPar);
                hRnEn->Fit(fRnEnGaus,"R0LQ");
                fRnEnGaus->SetRange(EnMin,EnMax);

		//-----------------
		fitMin = hPoEn->GetMean() - fitSigma*hPoEn->GetStdDev();
                fitMax = hPoEn->GetMean() + fitSigma*hPoEn->GetStdDev();

                fPoEnGaus = new TF1("fPoEnGaus","gaus",fitMin,fitMax);
                hPoEn->Fit(fPoEnGaus,"R0LQ");
                fPoEnGaus->SetRange(EnMin,EnMax);

		//-----------------
		fitMin = hRnPoDz->GetMean() - fitSigma*hRnPoDz->GetStdDev();
                fitMax = hRnPoDz->GetMean() + fitSigma*hRnPoDz->GetStdDev();

                fRnPoDzGaus = new TF1("fRnPoDzGaus","gaus",fitMin,fitMax);
                hRnPoDz->Fit(fRnPoDzGaus,"R0LQ");
                fRnPoDzGaus->SetRange(dzMin,dzMax);


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Calculate efficiencies
		promptLowPSDCut  = vRnPSDCutLow->at(i);
		promptHighPSDCut = vRnPSDCutHigh->at(i);	
		delayLowPSDCut   = vPoPSDCutLow->at(i);
		delayHighPSDCut  = vPoPSDCutHigh->at(i);	
		promptLowEnCut  = vRnEnCutLow->at(i);
		promptHighEnCut = vRnEnCutHigh->at(i);	
		delayLowEnCut   = vPoEnCutLow->at(i);
		delayHighEnCut  = vPoEnCutHigh->at(i);	
		dzCutLow  = vRnPoDzCutLow->at(i);
		dzCutHigh = vRnPoDzCutHigh->at(i);

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
	
		totEff = promptPSDEff * delayPSDEff * promptEnEff * delayEnEff * dzEff;	
		totEffErr = totEff * sqrt( pow(promptPSDEffErr/promptPSDEff,2) + pow(delayPSDEffErr/delayPSDEff,2) + pow(promptEnEffErr/promptEnEff,2) + pow(delayEnEffErr/delayEnEff,2) + pow(dzEffErr/dzEff,2) );


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Calculate rate
		double livetime = vL->at(i);

		N0 = fRnPoDtExp->GetParameter(0);
                N0Err = fRnPoDtExp->GetParError(0);
                lifetime = fRnPoDtExp->GetParameter(1);
                lifetimeErr = fRnPoDtExp->GetParError(1);

		rate = ((N0*lifetime)/(dtBinWidth*livetime*totEff))*(1e3);      //Hz
                rateErr = rate * sqrt( pow(N0Err/N0,2) + pow(lifetimeErr/lifetime,2) + 2*DtCov[0][1]/(N0*lifetime) + pow(totEffErr/totEff,2));
		BGRate = (hBGDt->GetEntries()/livetime)*(1e3);   //Hz

                printf("Rate: %.4f +/- %.4f \n",rate,rateErr);


		//---------------------------------------------------------------------------------
		//Populate vectors
		vRate.push_back(rate);
		vRateErr.push_back(rateErr);
		vTotEff.push_back(totEff);
		vTotEffErr.push_back(totEffErr);
		vLifetime.push_back(lifetime);
		vLifetimeErr.push_back(lifetimeErr);
                vBGRate.push_back(BGRate);
	
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

		vRnPSDChiSq.push_back(fRnPSDGaus->GetChisquare()/(double)fRnPSDGaus->GetNDF());
		vPoPSDChiSq.push_back(fPoPSDGaus->GetChisquare()/(double)fPoPSDGaus->GetNDF());
	
		vRnEnChiSq.push_back(fRnEnGaus->GetChisquare()/(double)fRnEnGaus->GetNDF());
		vPoEnChiSq.push_back(fPoEnGaus->GetChisquare()/(double)fPoEnGaus->GetNDF());
		
		vDzChiSq.push_back(fRnPoDzGaus->GetChisquare()/(double)fRnPoDzGaus->GetNDF());	
	}	//end while loop 

	histFile->Close();

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Initialize TGraphs for results
	printf("Creating TGraphs \n");

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

	TGraph *grLivetime 	 = new TGraph(numPt,x,y);
	TGraph *grTotLivetime 	 = new TGraph(numPt,x,y);
	TGraph *grPileupVeto	 = new TGraph(numPt,x,y);
	TGraph *grPileupVetoFrac = new TGraph(numPt,x,y);

	TGraphErrors *grPromptEnEff  = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDelayEnEff   = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPromptPSDEff = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDelayPSDEff  = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDzEff 	     = new TGraphErrors(numPt,x,y,xErr,yErr);

	TGraph *grDtChiSq    = new TGraph(numPt,x,y);
	TH1F *hDtChiSq 	     = new TH1F("hDtChiSq","dt chisq",20,0.5,1.5);

	TGraph *grRnPSDChiSq = new TGraph(numPt,x,y);
	TGraph *grPoPSDChiSq = new TGraph(numPt,x,y);
	TGraph *grRnEnChiSq  = new TGraph(numPt,x,y);
	TGraph *grPoEnChiSq  = new TGraph(numPt,x,y);
	TGraph *grDzChiSq    = new TGraph(numPt,x,y);

	TGraph *grBGRate = new TGraph(numPt,x,y);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Fill TGraphs
	for(int i=0;i<numPt;i++){
		double time = vT->at(i);
		printf("Time bin: %d ; Time: %f \n",i,time);

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

		grLivetime->SetPoint(i,time,vL->at(i));
		grTotLivetime->SetPoint(i,time,vTotLivetime->at(i));
		grPileupVeto->SetPoint(i,time,vPileupVetoTime->at(i));
		grPileupVetoFrac->SetPoint(i,time,vPileupVetoFrac->at(i));

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

		grRnPSDChiSq->SetPoint(i,time,vRnPSDChiSq[i]);
		grPoPSDChiSq->SetPoint(i,time,vPoPSDChiSq[i]);
		grRnEnChiSq->SetPoint(i,time,vRnEnChiSq[i]);
		grPoEnChiSq->SetPoint(i,time,vPoEnChiSq[i]);
		grDzChiSq->SetPoint(i,time,vDzChiSq[i]);
		grDtChiSq->SetPoint(i,time,vDtChiSq[i]);

		grBGRate->SetPoint(i,time,vBGRate[i]);
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
	grLivetime->Write("grLivetime");
	grTotLivetime->Write("grTotLivetime");
	grPileupVeto->Write("grPileupVeto");
	grPileupVetoFrac->Write("grPileupVetoFrac");
	grPromptEnEff->Write("grPromptEnEff");
	grDelayEnEff->Write("grDelayEnEff");
	grPromptPSDEff->Write("grPromptPSDEff");
	grDelayPSDEff->Write("grDelayPSDEff");
	grDzEff->Write("grDzEff");
	grRnPSDChiSq->Write("grRnPSDChiSq");
	grPoPSDChiSq->Write("grPoPSDChiSq");
	grRnEnChiSq->Write("grRnEnChiSq");
	grPoEnChiSq->Write("grPoEnChiSq");
	grDzChiSq->Write("grDzChiSq");
	grDtChiSq->Write("grDtChiSq");
	hDtChiSq->Write();
	grBGRate->Write("grBGRate");
	graphFile->Close();

}	//end void RnPoVsTime

