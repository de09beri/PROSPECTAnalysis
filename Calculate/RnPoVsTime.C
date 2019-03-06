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


void RnPoVsTime(double promptPSDStdDev, double delayPSDStdDev, double promptEnStdDev, double delayEnStdDev, double dzStdDev, double zLow, double zHigh, double timeBin, double dtCut, bool boolESmear){
	//---------------------------------------------------------------------------------
        //Read in text file with PSD and energy information
        //Creat bin by bin cuts for prompt and delay PSD and energy

	vector<int> vfileTimestamp;
        vector<double> vPromptPSDCutLow, vDelayPSDCutLow, vPromptEnCutLow, vDelayEnCutLow, vDzCutLow, vDzCutHigh;

        ifstream cutFile;
        cutFile.open("CutParameterVsTime.txt",ifstream::in);

        int fileTimebin, fileTimestamp;
        double RnPSD, RnPSDSigma, PoPSD, PoPSDSigma;
        double RnEn,  RnEnSigma,  PoEn,  PoEnSigma;
        double RnPoDz, RnPoDzSigma;

        double promptPSDCutLow, delayPSDCutLow, promptEnCutLow, delayEnCutLow, dzCutLow, dzCutHigh;

        string line;
        while(cutFile.good() & !cutFile.eof()){
                getline(cutFile,line);

                stringstream s(line);
                s>>fileTimebin>>fileTimestamp>>RnPSD>>RnPSDSigma>>PoPSD>>PoPSDSigma>>RnEn>>RnEnSigma>>PoEn>>PoEnSigma>>RnPoDz>>RnPoDzSigma;

                promptPSDCutLow = RnPSD - (promptPSDStdDev*RnPSDSigma);
                delayPSDCutLow  = PoPSD - (delayPSDStdDev*PoPSDSigma);
                promptEnCutLow  = RnEn  - (promptEnStdDev*RnEnSigma);
                delayEnCutLow   = PoEn  - (delayEnStdDev*PoEnSigma);
                dzCutLow        = RnPoDz - (dzStdDev*RnPoDzSigma);
                dzCutHigh       = RnPoDz + (dzStdDev*RnPoDzSigma);

		vfileTimestamp.push_back(fileTimestamp);
                vPromptPSDCutLow.push_back(promptPSDCutLow);
                vDelayPSDCutLow.push_back(delayPSDCutLow);
                vPromptEnCutLow.push_back(promptEnCutLow);
                vDelayEnCutLow.push_back(delayEnCutLow);
                vDzCutLow.push_back(dzCutLow);
                vDzCutHigh.push_back(dzCutHigh);
        }

        cutFile.close();


	//---------------------------------------------------------------------------------
	const double TIMEBREAK = timeBin*(3.6e6);	//[ms]

	TFile *histFile = new TFile(Form("%s/Ac227_HistsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");

	//---------------------------------------------------------------------------------
	TH1F *hSelectDt, 	*hBGDt,        *hRnPoDt;
	TH1F *hSelectPromptPSD, *hBGPromptPSD, *hRnPSD,	 *hSelectDelayPSD, *hBGDelayPSD, *hPoPSD;
	TH1F *hSelectPromptEn,  *hBGPromptEn,  *hRnEn,   *hSelectDelayEn,  *hBGDelayEn,  *hPoEn;
	TH1F *hSelectDelayEnSmear,  *hBGDelayEnSmear,  *hPoEnSmear;
	TH1F *hSelectPromptPos, *hBGPromptPos, *hRnPos,  *hSelectDelayPos, *hBGDelayPos, *hPoPos;
	TH1F *hSelectDz,        *hBGDz,        *hRnPoDz;

	TH1F *hSelectPromptTotEn, *hBGPromptTotEn, *hRnTotEn;

	TH2F *hSelectPSDvsEn, 			*hBGPSDvsEn, 		   *hRnPoPSDvsEn;
	TH2F *hSelectDelayEnvsPromptEn, *hBGDelayEnvsPromptEn, *hPoEnvsRnEn;

	//---------------------------------------------------------------------------------
	TF1 *fRnPoDtExp;
	TF1 *fRnPSDGaus, *fPoPSDGaus;
	TF1 *fRnEnGaus,  *fPoEnGaus;
	TF1 *fRnPoDzGaus;

	//---------------------------------------------------------------------------------
	vector<double> vLivetime,   vTimestamp;
	vector<double> vTotLivetime,  vPileupVetoT;
	vector<double> vPileupVetoFrac;

	//---------------------------------------------------------------------------------
	RNPO *rnpo = new RNPO();
	
	bool exclude;
	Long64_t numEntries = Long64_t(rnpo->fChain->GetEntries());
	printf("Number of Ac-227 candidates: %lld \n",numEntries);

	double setPromptLowPSDCut, setDelayLowPSDCut;
	double setPromptLowEnCut,  setDelayLowEnCut;
	double setDzCutLow, setDzCutHigh;

	double promptLowPSDCut, promptHighPSDCut, promptLowEnCut, promptHighEnCut;
	double delayLowPSDCut,  delayHighPSDCut,  delayLowEnCut,  delayHighEnCut;
	double dzCutLow, dzCutHigh;

	// Get cut values
	rnpo->GetEntry(0);
	setPromptLowPSDCut  = rnpo->p_PSDCut[0];
	promptHighPSDCut    = rnpo->p_PSDCut[1];
	setDelayLowPSDCut   = rnpo->d_PSDCut[0];
	delayHighPSDCut     = rnpo->d_PSDCut[1];
	setPromptLowEnCut   = rnpo->p_ECut[0];
    	promptHighEnCut     = rnpo->p_ECut[1];
    	setDelayLowEnCut    = rnpo->d_ECut[0];
    	delayHighEnCut      = rnpo->d_ECut[1];
	setDzCutLow  = -1.0*dzCut;
	setDzCutHigh = dzCut;

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Loop over all events, splitting into time periods of 24 hours

	Long64_t IDX = 0;
	int numTimeBin = 0;
	while(IDX<numEntries){

		promptLowPSDCut = (setPromptLowPSDCut > vPromptPSDCutLow[numTimeBin]) ? setPromptLowPSDCut : vPromptPSDCutLow[numTimeBin];
                delayLowPSDCut  = (setDelayLowPSDCut  > vDelayPSDCutLow[numTimeBin])  ? setDelayLowPSDCut  : vDelayPSDCutLow[numTimeBin];
                promptLowEnCut  = (setPromptLowEnCut  > vPromptEnCutLow[numTimeBin])  ? setPromptLowEnCut  : vPromptEnCutLow[numTimeBin];
                delayLowEnCut   = (setDelayLowEnCut   > vDelayEnCutLow[numTimeBin])   ? setDelayLowEnCut   : vDelayEnCutLow[numTimeBin];

                dzCutLow  = (setDzCutLow > vDzCutLow[numTimeBin])   ? setDzCutLow   : vDzCutLow[numTimeBin];
                dzCutHigh = (setDzCutHigh < vDzCutHigh[numTimeBin]) ? setDzHighzCut : vDzCutHigh[numTimeBin];

		//---------------------------------------------------------------------------------
		//Initialize histograms
		hSelectDt 		= new TH1F(Form("hSelectDt_%i",numTimeBin),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);
		hBGDt 			= new TH1F(Form("hBGDt_%i",numTimeBin),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);
	
		hSelectPromptPSD 	= new TH1F(Form("hSelectPromptPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	
		hBGPromptPSD 		= new TH1F(Form("hBGPromptPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	

		hSelectDelayPSD 	= new TH1F(Form("hSelectDelayPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	
		hBGDelayPSD 		= new TH1F(Form("hBGDelayPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	

		hSelectPromptEn 	= new TH1F(Form("hSelectPromptEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGPromptEn 		= new TH1F(Form("hBGPromptEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
	
		hSelectDelayEn 		= new TH1F(Form("hSelectDelayEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGDelayEn 		= new TH1F(Form("hBGDelayEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	

		hSelectPromptPos 	= new TH1F(Form("hSelectPromptPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);
		hBGPromptPos 		= new TH1F(Form("hBGPromptPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);

		hSelectDelayPos		= new TH1F(Form("hSelectDelayPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);
		hBGDelayPos 		= new TH1F(Form("hBGDelayPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);

		hSelectDz 		= new TH1F(Form("hSelectDz_%i",numTimeBin),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);
		hBGDz 			= new TH1F(Form("hBGDz_%i",numTimeBin),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);

		hSelectPromptTotEn 	= new TH1F(Form("hSelectPromptTotEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGPromptTotEn 		= new TH1F(Form("hBGPromptTotEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	

		hSelectDelayEnSmear 	= new TH1F(Form("hSelectDelayEnSmear_%i",numTimeBin),";ESmear [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGDelayEnSmear 	= new TH1F(Form("hBGDelayEnSmear_%i",numTimeBin),";ESmear [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	

		hSelectPSDvsEn 		= new TH2F(Form("hSelectPSDvsEn_%i",numTimeBin),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);
		hBGPSDvsEn 		= new TH2F(Form("hBGPSDvsEn_%i",numTimeBin),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);

		hSelectDelayEnvsPromptEn 	= new TH2F(Form("hSelectDelayEnvsPromptEn_%i",numTimeBin),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);
		hBGDelayEnvsPromptEn 		= new TH2F(Form("hBGDelayEnvsPromptEn_%i",numTimeBin),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);

		//---------------------------------------------------------------------------------
		//Fill histograms
		printf("=============== Filling Histograms =============== \n"); 
	
		int seg;
		double dt, dz;
		
		double livetime = 0.0, tstamp;
		double lastTime = 0.0, lastRunTime = 0.0, lastTimestamp = 0.0;	
		double sumWeightedTimestamp = 0.0, sumRunTime = 0.0;

		double lastNumClusts = 0.0, numClusts = 0.0;

		for(Long64_t i=IDX;i<numEntries;i++){
			if(i%100000==0) printf("Event: %lld \n",i);
			rnpo->GetEntry(i);

			if(rnpo->d_t*(1e-6) > ((double)((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1()*1000.0 - (TIMEWINDOW+TIMEOFFSET)) && i!=(numEntries-1)) continue;
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Calculate livetime and weighted timestamp
			if(rnpo->d_t < lastTime){ 
				livetime += lastTime*(1e-6);		//livetime in ms	

				sumWeightedTimestamp += lastRunTime * ((lastRunTime/2.0)+lastTimestamp);
				sumRunTime += lastRunTime;

				numClusts += lastNumClusts;
			}

			if(livetime>TIMEBREAK){

				tstamp = sumWeightedTimestamp/sumRunTime;
				vTimestamp.push_back(tstamp);	

				IDX = i;
				break;
			}

			lastTime = rnpo->d_t;
			lastRunTime = ((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1();		//seconds
			lastTimestamp = rnpo->tstamp;			//epoch seconds	
			lastNumClusts = rnpo->numClust;

			//if we are at the last entry	
			if(i == (numEntries-1)){ 
				livetime += lastTime*(1e-6);

				//if livetime is less than 12 hours 
				if(livetime*(2.778e-7) < 12) goto endloop;

				sumWeightedTimestamp += lastRunTime * ((lastRunTime/2.0)+lastTimestamp);
				sumRunTime += lastRunTime;

				tstamp = sumWeightedTimestamp/sumRunTime;
				vTimestamp.push_back(tstamp);	
	
				numClusts += lastNumClusts;

				IDX = i+1;
			}

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Fill histograms

			seg = rnpo->d_seg;

			if(rnpo->d_PSD < delayLowPSDCut || rnpo->d_E < delayLowEnCut) continue;	
			if(rnpo->d_z < zLow || rnpo->d_z > zHigh) continue;

			dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
			dz = rnpo->d_z - rnpo->p_z;
			if(rnpo->p_seg > -1 && rnpo->p_PSD>promptLowPSDCut && rnpo->p_E>promptLowEnCut && rnpo->p_z>zLow && rnpo->p_z<zHigh && dt>dtCut && dz>dzCutLow && dz<dzCutHigh){	
				dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
				dz = rnpo->d_z - rnpo->p_z;

				hSelectDt->Fill(dt);
				hSelectPromptPSD->Fill(rnpo->p_PSD);
				hSelectDelayPSD->Fill(rnpo->d_PSD);
				hSelectPromptEn->Fill(rnpo->p_E);
				hSelectDelayEn->Fill(rnpo->d_E);
				hSelectDelayEnSmear->Fill(rnpo->d_ESmear);
				hSelectPromptTotEn->Fill(rnpo->p_Etot);	
				hSelectPromptPos->Fill(rnpo->p_z);
				hSelectDelayPos->Fill(rnpo->d_z);
				hSelectDz->Fill(dz);

				hSelectPSDvsEn->Fill(rnpo->p_E,rnpo->p_PSD);
				hSelectPSDvsEn->Fill(rnpo->d_E,rnpo->d_PSD);
				hSelectDelayEnvsPromptEn->Fill(rnpo->p_E,rnpo->d_E);		
			}
			dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;	
			dz = rnpo->d_z - rnpo->f_z;
			if(rnpo->f_seg > -1 && rnpo->f_PSD>promptLowPSDCut && rnpo->f_E>promptLowEnCut && rnpo->f_z>zLow && rnpo->f_z<zHigh && dt>dtCut && dz>dzCutLow && dz<dzCutHigh){	
				dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;	
				dz = rnpo->d_z - rnpo->f_z;

				hBGDt->Fill(dt);
				hBGPromptPSD->Fill(rnpo->f_PSD);
				hBGDelayPSD->Fill(rnpo->d_PSD);
				hBGPromptEn->Fill(rnpo->f_E);
				hBGDelayEn->Fill(rnpo->d_E);
				hBGDelayEnSmear->Fill(rnpo->d_ESmear);
				hBGPromptTotEn->Fill(rnpo->f_Etot);	
				hBGPromptPos->Fill(rnpo->f_z);
				hBGDelayPos->Fill(rnpo->d_z);
				hBGDz->Fill(dz);

				hBGPSDvsEn->Fill(rnpo->f_E,rnpo->f_PSD);
				hBGPSDvsEn->Fill(rnpo->d_E,rnpo->d_PSD);
				hBGDelayEnvsPromptEn->Fill(rnpo->f_E,rnpo->d_E);		
			}

		}	//end for loop over events

	//===========================================================================================================

		double pileupVetoTime = numClusts*pileupVetoT;	//[ms]

		printf("Time bin: %i  |  Livetime: %f hrs \n",numTimeBin,livetime*(2.778e-7));
		printf("Pileup veto time: %f \n ms",pileupVetoTime);

		vTotLivetime.push_back(livetime/(1000.0*60.0));		//minutes
		vPileupVetoT.push_back(pileupVetoTime/(1000.0*60.0));	//minutes

		vPileupVetoFrac.push_back(pileupVetoTime/livetime);

		livetime = livetime - 2.0*pileupVetoTime;
		vLivetime.push_back(livetime);

		printf("Corrected Livetime: %f hours \n",livetime*(2.778e-7));

		//---------------------------------------------------------------------------------
		//Subtract histograms
		printf("=============== Subtracting Histograms =============== \n"); 

		hSelectDt->Sumw2();
		hRnPoDt = (TH1F*)hSelectDt->Clone();
		hRnPoDt->SetName(Form("hRnPoDt_%i",numTimeBin));
		hBGDt->Sumw2();
		hRnPoDt->Add(hBGDt,-1);

		hRnPSD = (TH1F*)hSelectPromptPSD->Clone();
		hRnPSD->SetName(Form("hRnPSD_%i",numTimeBin));
		hRnPSD->Sumw2();
		hRnPSD->Add(hBGPromptPSD,-1);

		hPoPSD = (TH1F*)hSelectDelayPSD->Clone();
		hPoPSD->SetName(Form("hPoPSD_%i",numTimeBin));
		hPoPSD->Sumw2();
		hPoPSD->Add(hBGDelayPSD,-1);

		hRnEn = (TH1F*)hSelectPromptEn->Clone();
		hRnEn->SetName(Form("hRnEn_%i",numTimeBin));
		hRnEn->Sumw2();
		hRnEn->Add(hBGPromptEn,-1);

		hRnTotEn = (TH1F*)hSelectPromptTotEn->Clone();
		hRnTotEn->SetName(Form("hRnTotEn_%i",numTimeBin));
		hRnTotEn->Sumw2();
		hRnTotEn->Add(hBGPromptTotEn,-1);

		hPoEn = (TH1F*)hSelectDelayEn->Clone();
		hPoEn->SetName(Form("hPoEn_%i",numTimeBin));
		hPoEn->Sumw2();
		hPoEn->Add(hBGDelayEn,-1);
	
		hPoEnSmear = (TH1F*)hSelectDelayEnSmear->Clone();
		hPoEnSmear->SetName(Form("hPoEnSmear_%i",numTimeBin));
		hPoEnSmear->Sumw2();
		hPoEnSmear->Add(hBGDelayEnSmear,-1);

		hRnPos = (TH1F*)hSelectPromptPos->Clone();
		hRnPos->SetName(Form("hRnPos_%i",numTimeBin));
		hRnPos->Sumw2();
		hRnPos->Add(hBGPromptPos,-1);

		hPoPos = (TH1F*)hSelectDelayPos->Clone();
		hPoPos->SetName(Form("hPoPos_%i",numTimeBin));
		hPoPos->Sumw2();
		hPoPos->Add(hBGDelayPos,-1);

		hRnPoDz = (TH1F*)hSelectDz->Clone();
		hRnPoDz->SetName(Form("hRnPoDz_%i",numTimeBin));
		hRnPoDz->Sumw2();
		hRnPoDz->Add(hBGDz,-1);

		hRnPoPSDvsEn = (TH2F*)hSelectPSDvsEn->Clone();
		hRnPoPSDvsEn->SetName(Form("hRnPoPSDvsEn_%i",numTimeBin));
		hRnPoPSDvsEn->Add(hBGPSDvsEn,-1);	

		hPoEnvsRnEn = (TH2F*)hSelectDelayEnvsPromptEn->Clone();
		hPoEnvsRnEn->SetName(Form("hPoEnvsRnEn_%i",numTimeBin));
		hPoEnvsRnEn->Add(hBGDelayEnvsPromptEn,-1);

		numTimeBin++;
	}	//end while loop IDX < numEntries

	endloop:
	histFile->Write();
	histFile->WriteObject(&vTimestamp,"vTimestamp");
	histFile->WriteObject(&vLivetime,"vLivetime");
	histFile->Close();

}	//end void RnPoVsTime

