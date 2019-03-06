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


void CutsPerCell(double dtCut){

	ofstream cutFile("/g/g20/berish1/AD_Ac227Analysis/PROSPECTAD_Ac227/Calculate/CutParameterPerCell.txt");

	//---------------------------------------------------------------------------------
	//Initialize histograms for dt, PSD, E, dz, and position
	int N = NUMCELLS;

	TH1F *hSelectPromptPSD[N], *hBGPromptPSD[N], *hRnPSD[N];
	TH1F *hSelectDelayPSD[N],  *hBGDelayPSD[N],  *hPoPSD[N];
	TH1F *hSelectPromptEn[N],  *hBGPromptEn[N],  *hRnEn[N];
	TH1F *hSelectDelayEn[N],   *hBGDelayEn[N],   *hPoEn[N];
	TH1F *hSelectDz[N], 	   *hBGDz[N], 		 *hRnPoDz[N];

	for(int i=0;i<NUMCELLS;i++){
		hSelectPromptPSD[i] = new TH1F(Form("hSelectPromptPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
		hBGPromptPSD[i] 	= new TH1F(Form("hBGPromptPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
	
		hSelectDelayPSD[i] 	= new TH1F(Form("hSelectDelayPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
		hBGDelayPSD[i] 		= new TH1F(Form("hBGDelayPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
			
		hSelectPromptEn[i] 	= new TH1F(Form("hSelectPromptEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		hBGPromptEn[i] 		= new TH1F(Form("hBGPromptEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
	
		hSelectDelayEn[i] 	= new TH1F(Form("hSelectDelayEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		hBGDelayEn[i] 		= new TH1F(Form("hBGDelayEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);

		hSelectDz[i]		= new TH1F(Form("hSelectDz_%i",i),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);
		hBGDz[i]			= new TH1F(Form("hBGDz_%i",i),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);
	}	//end for loop creating histograms

	//---------------------------------------------------------------------------------
	//Fill histograms
	printf("=============== Filling Histograms =============== \n"); 

	RNPO *rnpo = new RNPO();
	
	bool exclude;
	Long64_t numEntries = Long64_t(rnpo->fChain->GetEntries());
	printf("Number of Ac-227 candidates: %lld \n",numEntries);

	double tstamp;
	double promptLowPSDCut, promptHighPSDCut, promptLowEnCut, promptHighEnCut;
	double delayLowPSDCut,  delayHighPSDCut,  delayLowEnCut,  delayHighEnCut;
	double dzCut;

	// Get cut values
	rnpo->GetEntry(0);
	tstamp 		 = rnpo->tstamp;
	//promptLowPSDCut  = rnpo->p_PSDCut[0]; 
	promptLowPSDCut  = 0.19; 
	promptHighPSDCut = rnpo->p_PSDCut[1];
	//delayLowPSDCut   = rnpo->d_PSDCut[0];
	delayLowPSDCut   = 0.19;
	delayHighPSDCut  = rnpo->d_PSDCut[1];
	promptLowEnCut   = rnpo->p_ECut[0];
	promptHighEnCut  = rnpo->p_ECut[1];
	delayLowEnCut    = rnpo->d_ECut[0];
	delayHighEnCut   = rnpo->d_ECut[1];
	//dzCut            = rnpo->dzCut;
	dzCut            = 200;


	//Initialize variables
	int seg;
	double dt, dz;

	for(Long64_t i=0;i<numEntries;i++){
		if(i%1000000==0) printf("Event: %lld \n",i);
		rnpo->GetEntry(i);

		if(rnpo->d_t*(1e-6) > ((double)((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1()*1000.0 - (TIMEWINDOW+TIMEOFFSET))) continue;

		seg = rnpo->d_seg;

		double rnpo_p_E = rnpo->p_E;	
		double rnpo_d_E = rnpo->d_E;
		double rnpo_f_E = rnpo->f_E;

		if(rnpo->d_PSD < delayLowPSDCut || rnpo->d_E < delayLowEnCut) continue;

		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), seg) != end(ExcludeCellArr);
		if(exclude) continue;


		//if prompt-delay pair
		dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
		dz = rnpo->d_z - rnpo->p_z;
		if(rnpo->p_seg > -1 && rnpo->p_PSD>promptLowPSDCut && rnpo->p_E>promptLowEnCut && dt>dtCut && rnpo->p_z>-500 && rnpo->p_z<500 && abs(dz)<dzCut){
			dz = rnpo->d_z - rnpo->p_z;
			hSelectPromptPSD[seg]->Fill(rnpo->p_PSD);
			hSelectDelayPSD[seg]->Fill(rnpo->d_PSD);
			hSelectPromptEn[seg]->Fill(rnpo_p_E);
			hSelectDelayEn[seg]->Fill(rnpo_d_E);
			hSelectDz[seg]->Fill(dz);
		}

		//if prompt-delay BG pair
		dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;
		dz = rnpo->d_z - rnpo->f_z;
		if(rnpo->f_seg > -1 && rnpo->f_PSD>promptLowPSDCut && rnpo->f_E>promptLowEnCut && dt>dtCut && rnpo->f_z>-500 && rnpo->f_z<500 && abs(dz)<dzCut){
			dz = rnpo->d_z - rnpo->f_z;
			hBGPromptPSD[seg]->Fill(rnpo->f_PSD);
			hBGDelayPSD[seg]->Fill(rnpo->d_PSD);
			hBGPromptEn[seg]->Fill(rnpo_f_E);
			hBGDelayEn[seg]->Fill(rnpo_d_E);
			hBGDz[seg]->Fill(dz);
		}
	}	//end for loop over TChain


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Subtract histograms
	printf("=============== Subtracting Histograms =============== \n"); 
	TF1 *fRnPSDGaus, *fPoPSDGaus;
	TF1 *fRnEnGaus,  *fPoEnGaus;
	TF1 *fRnPoDzGaus;

	for(int i=0;i<NUMCELLS;i++){
		double RnPSD = 0, RnPSDSigma = 0, PoPSD = 0, PoPSDSigma = 0;
		double RnEn = 0,  RnEnSigma = 0,  PoEn = 0,  PoEnSigma = 0;
		double RnPoDz = 0, RnPoDzSigma = 0;

		if(i%10==0) printf("Cell: %d \n",i);

		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
		if(exclude){ 
			cutFile<<std::fixed<<i<<" "<<std::setprecision(4)<<RnPSD<<" "<<std::setprecision(4)<<RnPSDSigma<<" "<<std::setprecision(4)<<PoPSD<<" "<<std::setprecision(4)<<PoPSDSigma<<" "<<std::setprecision(4)<<RnEn<<" "<<std::setprecision(4)<<RnEnSigma<<" "<<std::setprecision(4)<<PoEn<<" "<<std::setprecision(4)<<PoEnSigma<<" "<<std::setprecision(4)<<RnPoDz<<" "<<std::setprecision(4)<<RnPoDzSigma<<"\n";	
			continue;
		}

		hRnPSD[i] = (TH1F*)hSelectPromptPSD[i]->Clone();
		hRnPSD[i]->SetName(Form("hRnPSD_%i",i));
		hRnPSD[i]->Sumw2();
		hRnPSD[i]->Add(hBGPromptPSD[i],-1);

		hPoPSD[i] = (TH1F*)hSelectDelayPSD[i]->Clone();
		hPoPSD[i]->SetName(Form("hPoPSD_%i",i));
		hPoPSD[i]->Sumw2();
		hPoPSD[i]->Add(hBGDelayPSD[i],-1);

		hRnEn[i] = (TH1F*)hSelectPromptEn[i]->Clone();
		hRnEn[i]->SetName(Form("hRnEn_%i",i));
		hRnEn[i]->Sumw2();
		hRnEn[i]->Add(hBGPromptEn[i],-1);

		hPoEn[i] = (TH1F*)hSelectDelayEn[i]->Clone();
		hPoEn[i]->SetName(Form("hPoEn_%i",i));
		hPoEn[i]->Sumw2();
		hPoEn[i]->Add(hBGDelayEn[i],-1);
	
		hRnPoDz[i] = (TH1F*)hSelectDz[i]->Clone();
		hRnPoDz[i]->SetName(Form("hRnPoDz_%i",i));
		hRnPoDz[i]->Sumw2();
		hRnPoDz[i]->Add(hBGDz[i],-1);

		double RnPSDMaxBin = hRnPSD[i]->GetBinCenter(hRnPSD[i]->GetMaximumBin());
		fRnPSDGaus = new TF1("fRnPSDGaus","gaus",RnPSDMaxBin-0.02,RnPSDMaxBin+0.02);
		hRnPSD[i]->Fit(fRnPSDGaus,"RQ");
		RnPSD = fRnPSDGaus->GetParameter(1);
		RnPSDSigma = fRnPSDGaus->GetParameter(2);	

		double PoPSDMaxBin = hPoPSD[i]->GetBinCenter(hPoPSD[i]->GetMaximumBin());
		fPoPSDGaus = new TF1("fPoPSDGaus","gaus",PoPSDMaxBin-0.02,PoPSDMaxBin+0.02);
		hPoPSD[i]->Fit(fPoPSDGaus,"RQ");
		PoPSD = fPoPSDGaus->GetParameter(1);
		PoPSDSigma = fPoPSDGaus->GetParameter(2);	

		double RnEnMaxBin = hRnEn[i]->GetBinCenter(hRnEn[i]->GetMaximumBin());
		fRnEnGaus = new TF1("fRnEnGaus","gaus",RnEnMaxBin-0.04,RnEnMaxBin+0.04);	
		hRnEn[i]->Fit(fRnEnGaus,"RQ");
		RnEn = fRnEnGaus->GetParameter(1);
		RnEnSigma = fRnEnGaus->GetParameter(2);

		double PoEnMaxBin = hPoEn[i]->GetBinCenter(hPoEn[i]->GetMaximumBin());
		fPoEnGaus = new TF1("fPoEnGaus","gaus",PoEnMaxBin-0.04,PoEnMaxBin+0.04);	
		hPoEn[i]->Fit(fPoEnGaus,"RQ");
		PoEn = fPoEnGaus->GetParameter(1);
		PoEnSigma = fPoEnGaus->GetParameter(2);

		double RnPoDzMaxBin = hRnPoDz[i]->GetBinCenter(hRnPoDz[i]->GetMaximumBin());
		fRnPoDzGaus = new TF1("fRnPoDzGaus","gaus",RnPoDzMaxBin-50,RnPoDzMaxBin+50);
		hRnPoDz[i]->Fit(fRnPoDzGaus,"RQ");
		RnPoDz = fRnPoDzGaus->GetParameter(1);
		RnPoDzSigma = fRnPoDzGaus->GetParameter(2);	

		cutFile<<std::fixed<<i<<" "<<std::setprecision(4)<<RnPSD<<" "<<std::setprecision(4)<<RnPSDSigma<<" "<<std::setprecision(4)<<PoPSD<<" "<<std::setprecision(4)<<PoPSDSigma<<" "<<std::setprecision(4)<<RnEn<<" "<<std::setprecision(4)<<RnEnSigma<<" "<<std::setprecision(4)<<PoEn<<" "<<std::setprecision(4)<<PoEnSigma<<" "<<std::setprecision(4)<<RnPoDz<<" "<<std::setprecision(4)<<RnPoDzSigma<<"\n";	

	}	//end for loop to subtract hists

	cutFile.close();
}	//end void RnPoVsCell

