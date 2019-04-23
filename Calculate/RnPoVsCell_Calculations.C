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

void RnPoVsCell_Calc(string outputName){
	//---------------------------------------------------------------------------------
	TFile *histFile = new TFile(Form("%s/Ac227_HistsPerCell_%s.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"),outputName.c_str()),"UPDATE");

	vector<double> *vL;
	histFile->GetObject("vLivetime",vL);

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
	TF1 *fRnEnGaus,  *fPoEnGaus, *fPoEnSmearGaus;
	TF1 *fRnPoDzGaus;

	//---------------------------------------------------------------------------------
	vector<double> vSeg;
	vector<double> vRate,   vRateErr;
	vector<double> vTotEff, vTotEffErr;
	vector<double> vLifetime,   vLifetimeErr;
	vector<double> vPoPSDMean,  vPoPSDMeanErr,  vPoPSDSigma,  vPoPSDSigmaErr;
	vector<double> vPoEnMean,   vPoEnMeanErr,   vPoEnSigma,   vPoEnSigmaErr;
	vector<double> vPoEnSmearMean,   vPoEnSmearMeanErr,   vPoEnSmearSigma,   vPoEnSmearSigmaErr;
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

	bool exclude;

	for(int i=0;i<NUMCELLS;i++){
		if(i%10==0) printf("Cell: %d \n",i);
		
		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
                if(exclude) continue;

		hBGDt   = (TH1F*)histFile->Get(Form("hBGDt_%i",i));
		hRnPoDt = (TH1F*)histFile->Get(Form("hRnPoDt_%i",i));
		hRnPSD  = (TH1F*)histFile->Get(Form("hRnPSD_%i",i));
		hPoPSD  = (TH1F*)histFile->Get(Form("hPoPSD_%i",i));
		hRnEn   = (TH1F*)histFile->Get(Form("hRnEn_%i",i));
		hPoEn   = (TH1F*)histFile->Get(Form("hPoEn_%i",i));
		hPoEnSmear   = (TH1F*)histFile->Get(Form("hPoEnSmear_%i",i));
		hRnPoDz = (TH1F*)histFile->Get(Form("hRnPoDz_%i",i));

		hPoPos   = (TH1F*)histFile->Get(Form("hPoPos_%i",i));

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Fit distributions
		printf("========= Cell %d ========= \n",i);
		printf("Fitting distributions \n");

		fRnPoDtExp = new TF1("fRnPoDtExp","[0]*exp(-x/[1])",0.5,8);
                fRnPoDtExp->SetParameter(1,POLIFETIME);
                TFitResultPtr fitr = hRnPoDt->Fit(fRnPoDtExp,"RS0Q");
                TMatrixD DtCov = fitr->GetCovarianceMatrix();

		hRnPoDt->Write(Form("hRnPoDt_%i",i),TObject::kOverwrite);

		//-----------------
		double fitSigma = 2.0;
                double fitMin, fitMax;
		
		//-----------------
		fitMin = hRnPSD->GetMean() - fitSigma*hRnPSD->GetStdDev();
                fitMax = hRnPSD->GetMean() + fitSigma*hRnPSD->GetStdDev();

                fRnPSDGaus = new TF1("fRnPSDGaus","gaus",fitMin,fitMax);
                hRnPSD->Fit(fRnPSDGaus,"R0LQ");
                fRnPSDGaus->SetRange(PSDMin,PSDMax);

		hRnPSD->Write(Form("hRnPSD_%i",i),TObject::kOverwrite);

		//-----------------
		fitMin = hPoPSD->GetMean() - fitSigma*hPoPSD->GetStdDev();
                fitMax = hPoPSD->GetMean() + fitSigma*hPoPSD->GetStdDev();

                fPoPSDGaus = new TF1("fPoPSDGaus","gaus",fitMin,fitMax);
                hPoPSD->Fit(fPoPSDGaus,"R0LQ");
                fPoPSDGaus->SetRange(PSDMin,PSDMax);

		hPoPSD->Write(Form("hPoPSD_%i",i),TObject::kOverwrite);
	
		//-----------------
		fitMin = hRnEn->GetMean() - 1.3*hRnEn->GetStdDev();
                fitMax = hRnEn->GetMean() + 0.6*hRnEn->GetStdDev();

                fRnEnGaus = new TF1("fRnEnGaus","gaus",fitMin,fitMax);
	        hRnEn->Fit(fRnEnGaus,"R0LQ");
                fRnEnGaus->SetRange(EnMin,EnMax);

		hRnEn->Write(Form("hRnEn_%i",i),TObject::kOverwrite);

		//-----------------
		fitMin = hPoEn->GetMean() - fitSigma*hPoEn->GetStdDev();
                fitMax = hPoEn->GetMean() + fitSigma*hPoEn->GetStdDev();

                fPoEnGaus = new TF1("fPoEnGaus","gaus",fitMin,fitMax);
                hPoEn->Fit(fPoEnGaus,"R0LQ");
                fPoEnGaus->SetRange(EnMin,EnMax);

		hPoEn->Write(Form("hPoEn_%i",i),TObject::kOverwrite);

		//-----------------
		fitMin = hPoEnSmear->GetMean() - fitSigma*hPoEnSmear->GetStdDev();
                fitMax = hPoEnSmear->GetMean() + fitSigma*hPoEnSmear->GetStdDev();

                fPoEnSmearGaus = new TF1("fPoEnSmearGaus","gaus",fitMin,fitMax);
                hPoEnSmear->Fit(fPoEnSmearGaus,"R0LQ");
                fPoEnSmearGaus->SetRange(EnMin,EnMax);

		hPoEnSmear->Write(Form("hPoEnSmear_%i",i),TObject::kOverwrite);

		//-----------------
		fitMin = hRnPoDz->GetMean() - fitSigma*hRnPoDz->GetStdDev();
                fitMax = hRnPoDz->GetMean() + fitSigma*hRnPoDz->GetStdDev();

                fRnPoDzGaus = new TF1("fRnPoDzGaus","gaus",fitMin,fitMax);
                hRnPoDz->Fit(fRnPoDzGaus,"R0LQ");
                fRnPoDzGaus->SetRange(dzMin,dzMax);

		hRnPoDz->Write(Form("hRnPoDz_%i",i),TObject::kOverwrite);

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

		printf("%f < Rn PSD < %f \n",promptLowPSDCut,promptHighPSDCut);
		printf("%f < Po PSD < %f \n",delayLowPSDCut,delayHighPSDCut);
		printf("%f < Rn En < %f \n",promptLowEnCut,promptHighEnCut);
		printf("%f < Po En < %f \n",delayLowEnCut,delayHighEnCut);
		printf("%f < RnPo Dz < %f \n",dzCutLow,dzCutHigh);	

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

		printf("Efficiency: %f +/- %f \n",totEff,totEffErr);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Calculate rate
		double livetime = vL->at(0);

		N0 = fRnPoDtExp->GetParameter(0);
                N0Err = fRnPoDtExp->GetParError(0);
                lifetime = fRnPoDtExp->GetParameter(1);
                lifetimeErr = fRnPoDtExp->GetParError(1);

		rate = ((N0*lifetime)/(dtBinWidth*livetime*totEff))*(1e6);      //mHz
                rateErr = rate * sqrt( pow(N0Err/N0,2) + pow(lifetimeErr/lifetime,2) + 2*DtCov[0][1]/(N0*lifetime) + pow(totEffErr/totEff,2));
		BGRate = (hBGDt->GetEntries()/livetime)*(1e6);   //Hz

                printf("Rate: %.4f +/- %.4f \n",rate,rateErr);

		//---------------------------------------------------------------------------------
		//Populate vectors
		vSeg.push_back(i);
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

		vPoEnSmearMean.push_back(fPoEnSmearGaus->GetParameter(1));
		vPoEnSmearMeanErr.push_back(fPoEnSmearGaus->GetParError(1));
		vPoEnSmearSigma.push_back(fPoEnSmearGaus->GetParameter(2));
		vPoEnSmearSigmaErr.push_back(fPoEnSmearGaus->GetParError(2));

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

	int numPt = NUMCELLS - NUMEXCLUDECELLS;
        double x[numPt], xErr[numPt];
        double y[1], yErr[1];

	TGraphErrors *grRate 		= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grTotEff 		= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grLifetime 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPSDMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPSDSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoEnMean   	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoEnSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoEnSmearMean   = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoEnSmearSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPosMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPosSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);

	TGraphErrors *grPromptEnEff  = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDelayEnEff   = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPromptPSDEff = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDelayPSDEff  = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDzEff 	     = new TGraphErrors(numPt,x,y,xErr,yErr);

	TGraph *grDtChiSq    = new TGraph(numPt,x,y);

	TGraph *grRnPSDChiSq = new TGraph(numPt,x,y);
	TGraph *grPoPSDChiSq = new TGraph(numPt,x,y);
	TGraph *grRnEnChiSq  = new TGraph(numPt,x,y);
	TGraph *grPoEnChiSq  = new TGraph(numPt,x,y);
	TGraph *grDzChiSq    = new TGraph(numPt,x,y);

	TGraph *grBGRate = new TGraph(numPt,x,y);

	//-------------------
	TH1F *hRate = new TH1F("hRate","Rate",20,3.14,3.34);
	TH1F *hDtChiSq = new TH1F("hDtChiSq","dt chisq",20,0.5,1.5);


	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Fill TGraphs
	for(int i=0;i<numPt;i++){
		double seg = vSeg[i];

		grRate->SetPoint(i,seg,vRate[i]);
		grRate->SetPointError(i,0,vRateErr[i]);

		hRate->Fill(vRate[i]);

		grTotEff->SetPoint(i,seg,vTotEff[i]);
		grTotEff->SetPointError(i,0,vTotEffErr[i]);

		grLifetime->SetPoint(i,seg,vLifetime[i]);
		grLifetime->SetPointError(i,0,vLifetimeErr[i]);

		grPoPSDMean->SetPoint(i,seg,vPoPSDMean[i]);
		grPoPSDMean->SetPointError(i,0,vPoPSDMeanErr[i]);

		grPoPSDSigma->SetPoint(i,seg,vPoPSDSigma[i]);
		grPoPSDSigma->SetPointError(i,0,vPoPSDSigmaErr[i]);

		grPoEnMean->SetPoint(i,seg,vPoEnMean[i]);
		grPoEnMean->SetPointError(i,0,vPoEnMeanErr[i]);
	
		grPoEnSigma->SetPoint(i,seg,vPoEnSigma[i]);
		grPoEnSigma->SetPointError(i,0,vPoEnSigmaErr[i]);

		grPoEnSmearMean->SetPoint(i,seg,vPoEnSmearMean[i]);
		grPoEnSmearMean->SetPointError(i,0,vPoEnSmearMeanErr[i]);
	
		grPoEnSmearSigma->SetPoint(i,seg,vPoEnSmearSigma[i]);
		grPoEnSmearSigma->SetPointError(i,0,vPoEnSmearSigmaErr[i]);

		grPoPosMean->SetPoint(i,seg,vPoPosMean[i]);
		grPoPosMean->SetPointError(i,0,vPoPosMeanErr[i]);

		grPoPosSigma->SetPoint(i,seg,vPoPosSigma[i]);
		grPoPosSigma->SetPointError(i,0,vPoPosSigmaErr[i]);		

		grRnPoDzMean->SetPoint(i,seg,vRnPoDzMean[i]);
		grRnPoDzMean->SetPointError(i,0,vRnPoDzMeanErr[i]);

		grRnPoDzSigma->SetPoint(i,seg,vRnPoDzSigma[i]);
		grRnPoDzSigma->SetPointError(i,0,vRnPoDzSigmaErr[i]);

		grPromptEnEff->SetPoint(i,seg,vPromptEnEff[i]);
		grPromptEnEff->SetPointError(i,0,vPromptEnEffErr[i]);

		grDelayEnEff->SetPoint(i,seg,vDelayEnEff[i]);
		grDelayEnEff->SetPointError(i,0,vDelayEnEffErr[i]);

		grPromptPSDEff->SetPoint(i,seg,vPromptPSDEff[i]);
		grPromptPSDEff->SetPointError(i,0,vPromptPSDEffErr[i]);

		grDelayPSDEff->SetPoint(i,seg,vDelayPSDEff[i]);
		grDelayPSDEff->SetPointError(i,0,vDelayPSDEffErr[i]);

		grDzEff->SetPoint(i,seg,vDzEff[i]);
		grDzEff->SetPointError(i,0,vDzEffErr[i]);	

		grDtChiSq->SetPoint(i,seg,vDtChiSq[i]);
		hDtChiSq->Fill(vDtChiSq[i]);

		grRnPSDChiSq->SetPoint(i,seg,vRnPSDChiSq[i]);
		grPoPSDChiSq->SetPoint(i,seg,vPoPSDChiSq[i]);
		grRnEnChiSq->SetPoint(i,seg,vRnEnChiSq[i]);
		grPoEnChiSq->SetPoint(i,seg,vPoEnChiSq[i]);
		grDzChiSq->SetPoint(i,seg,vDzChiSq[i]);
		grDtChiSq->SetPoint(i,seg,vDtChiSq[i]);

		grBGRate->SetPoint(i,seg,vBGRate[i]);
	}	//end for loop to populate TGraphs

	//---------------------------------------------------------------------------------
	//Write TGraphs to file
	TFile *graphFile = new TFile(Form("%s/Ac227_GraphsPerCell_%s.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"),outputName.c_str()),"RECREATE");
	
	grRate->Write("grRate");
	grTotEff->Write("grTotEff");	
	grLifetime->Write("grLifetime");
	grPoPSDMean->Write("grPoPSDMean");
	grPoPSDSigma->Write("grPoPSDSigma");
	grPoEnMean->Write("grPoEnMean");
	grPoEnSigma->Write("grPoEnSigma");
	grPoEnSmearMean->Write("grPoEnSmearMean");
	grPoEnSmearSigma->Write("grPoEnSmearSigma");
	grPoPosMean->Write("grPoPosMean");
	grPoPosSigma->Write("grPoPosSigma");
	grRnPoDzMean->Write("grRnPoDzMean");
	grRnPoDzSigma->Write("grRnPoDzSigma");	
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
	grBGRate->Write("grBGRate");
	hDtChiSq->Write();
	hRate->Write();
	graphFile->Close();

}	//end void RnPoVsTime

