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

void RnPoVsCell(double promptPSDStdDev, double delayPSDStdDev, double promptEnStdDev, double delayEnStdDev, double dzStdDev, double zLow, double zHigh, double dtCut, bool boolESmear){
	RNPO *rnpo = new RNPO();

	double setPromptLowPSDCut, setDelayLowPSDCut;
        double setPromptLowEnCut,  setDelayLowEnCut;
        double setDzCutLow, setDzCutHigh;

	double calcPromptLowPSDCut, calcDelayLowPSDCut;
        double calcPromptLowEnCut,  calcDelayLowEnCut;
        double calcDzCutLow, calcDzCutHigh;

	double promptLowPSDCut, promptHighPSDCut, promptLowEnCut, promptHighEnCut;
	double delayLowPSDCut,  delayHighPSDCut,  delayLowEnCut,  delayHighEnCut;
	double dzCut, dzCutLow, dzCutHigh;

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
	dzCut               = rnpo->dzCut;
	setDzCutLow  = -1.0*dzCut;
	setDzCutHigh = dzCut;

        //---------------------------------------------------------------------------------
        //Read in text file with PSD and energy information
        //Creat bin by bin cuts for prompt and delay PSD and energy
	vector<double> vPromptPSDCutLow, vDelayPSDCutLow, vPromptEnCutLow, vDelayEnCutLow, vDzCutLow, vDzCutHigh;
	vector<double> vPromptPSDCutHigh, vDelayPSDCutHigh, vPromptEnCutHigh, vDelayEnCutHigh;

        ifstream cutFile;
	cutFile.open("/g/g20/berish1/PROSPECTAnalysis/Calculate/CutParameterPerCell.txt",ifstream::in);

        int cellNum;
        double RnPSD, RnPSDSigma, PoPSD, PoPSDSigma;
        double RnEn,  RnEnSigma,  PoEn,  PoEnSigma;
        double RnPoDz, RnPoDzSigma;

        string line;
        while(cutFile.good() & !cutFile.eof()){
                getline(cutFile,line);

                stringstream s(line);
                s>>cellNum>>RnPSD>>RnPSDSigma>>PoPSD>>PoPSDSigma>>RnEn>>RnEnSigma>>PoEn>>PoEnSigma>>RnPoDz>>RnPoDzSigma;

                calcPromptLowPSDCut = RnPSD - (promptPSDStdDev*RnPSDSigma);
                calcDelayLowPSDCut  = PoPSD - (delayPSDStdDev*PoPSDSigma);
                calcPromptLowEnCut  = RnEn  - (promptEnStdDev*RnEnSigma);
                calcDelayLowEnCut   = PoEn  - (delayEnStdDev*PoEnSigma);
                calcDzCutLow        = RnPoDz - (dzStdDev*RnPoDzSigma);
                calcDzCutHigh       = RnPoDz + (dzStdDev*RnPoDzSigma);


		promptLowPSDCut = (setPromptLowPSDCut > calcPromptLowPSDCut) ? setPromptLowPSDCut : calcPromptLowPSDCut;
                delayLowPSDCut  = (setDelayLowPSDCut  > calcDelayLowPSDCut)  ? setDelayLowPSDCut  : calcDelayLowPSDCut;
                promptLowEnCut  = (setPromptLowEnCut  > calcPromptLowEnCut)  ? setPromptLowEnCut  : calcPromptLowEnCut;
                delayLowEnCut   = (setDelayLowEnCut   > calcDelayLowEnCut)   ? setDelayLowEnCut   : calcDelayLowEnCut;

                dzCutLow  = (setDzCutLow > calcDzCutLow)   ? setDzCutLow  : calcDzCutLow;
                dzCutHigh = (setDzCutHigh < calcDzCutHigh) ? setDzCutHigh : calcDzCutHigh;


                vPromptPSDCutLow.push_back(promptLowPSDCut);
                vDelayPSDCutLow.push_back(delayLowPSDCut);
                vPromptEnCutLow.push_back(promptLowEnCut);
                vDelayEnCutLow.push_back(delayLowEnCut);
                vDzCutLow.push_back(dzCutLow);
                vDzCutHigh.push_back(dzCutHigh);

		vPromptPSDCutHigh.push_back(promptHighPSDCut);
		vDelayPSDCutHigh.push_back(delayHighPSDCut);
		vPromptEnCutHigh.push_back(promptHighEnCut);
		vDelayEnCutHigh.push_back(delayHighEnCut);
        }

        cutFile.close();
        //---------------------------------------------------------------------------------

	TFile *histFile = new TFile(Form("%s/Ac227_HistsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");

	//---------------------------------------------------------------------------------
	//Initialize histograms for dt, PSD, E, dz, and position
	int N = NUMCELLS;

	TH1F *hSelectDt[N],	   *hBGDt[N],	     *hRnPoDt[N];
	TH1F *hSelectPromptPSD[N], *hBGPromptPSD[N], *hRnPSD[N];
	TH1F *hSelectDelayPSD[N],  *hBGDelayPSD[N],  *hPoPSD[N];
	TH1F *hSelectPromptEn[N],  *hBGPromptEn[N],  *hRnEn[N];
	TH1F *hSelectDelayEn[N],   *hBGDelayEn[N],   *hPoEn[N];
	TH1F *hSelectPromptPos[N], *hBGPromptPos[N], *hRnPos[N];
	TH1F *hSelectDelayPos[N],  *hBGDelayPos[N],  *hPoPos[N];
	TH1F *hSelectDz[N], 	   *hBGDz[N], 	     *hRnPoDz[N];

	TH1F *hSelectPromptTotEn[N],  *hBGPromptTotEn[N],  *hRnTotEn[N];
	TH1F *hSelectDelayEnSmear[N], *hBGDelayEnSmear[N], *hPoEnSmear[N];

	TH2F *hSelectPSDvsEn[N],   *hBGPSDvsEn[N], 	 *hRnPoPSDvsEn[N];
	TH2F *hSelectPSDvsPos[N],  *hBGPSDvsPos[N],      *hRnPoPSDvsPos[N];
	TH2F *hSelectEnvsPos[N],   *hBGEnvsPos[N],	 *hRnPoEnvsPos[N];
	TH2F *hSelectDelayEnvsPromptEn[N], *hBGDelayEnvsPromptEn[N], *hPoEnvsRnEn[N];

	for(int i=0;i<NUMCELLS;i++){
		hSelectDt[i] 		= new TH1F(Form("hSelectDt_%i",i),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);	
		hBGDt[i] 	 	= new TH1F(Form("hBGDt_%i",i),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);	
		
		hSelectPromptPSD[i] 	= new TH1F(Form("hSelectPromptPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
		hBGPromptPSD[i] 	= new TH1F(Form("hBGPromptPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
	
		hSelectDelayPSD[i] 	= new TH1F(Form("hSelectDelayPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
		hBGDelayPSD[i] 		= new TH1F(Form("hBGDelayPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
			
		hSelectPromptEn[i] 	= new TH1F(Form("hSelectPromptEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		hBGPromptEn[i] 		= new TH1F(Form("hBGPromptEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
	
		hSelectDelayEn[i] 	= new TH1F(Form("hSelectDelayEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		hBGDelayEn[i] 		= new TH1F(Form("hBGDelayEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);

		hSelectPromptPos[i] 	= new TH1F(Form("hSelectPromptPos_%i",i),";z [mm];Counts/cm",numPosBins,posMin,posMax);
		hBGPromptPos[i] 	= new TH1F(Form("hBGPromptPos_%i",i),";z [mm];Counts/cm",numPosBins,posMin,posMax);

		hSelectDelayPos[i] 	= new TH1F(Form("hSelectDelayPos_%i",i),";z [mm];Counts/cm",numPosBins,posMin,posMax);
		hBGDelayPos[i] 		= new TH1F(Form("hBGDelayPos_%i",i),";z [mm];Counts/cm",numPosBins,posMin,posMax);

		hSelectDz[i]		= new TH1F(Form("hSelectDz_%i",i),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);
		hBGDz[i]		= new TH1F(Form("hBGDz_%i",i),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);

		hSelectPromptTotEn[i] 	= new TH1F(Form("hSelectPromptTotEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		hBGPromptTotEn[i] 	= new TH1F(Form("hBGPromptTotEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		
		hSelectDelayEnSmear[i] 	= new TH1F(Form("hSelectDelayEnSmear_%i",i),";ESmear [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		hBGDelayEnSmear[i] 	= new TH1F(Form("hBGDelayEnSmear_%i",i),";ESmear [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);

		hSelectPSDvsEn[i] 	= new TH2F(Form("hSelectPSDvsEn_%i",i),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);
		hBGPSDvsEn[i] 		= new TH2F(Form("hBGPSDvsEn_%i",i),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);

		hSelectPSDvsPos[i] 	= new TH2F(Form("hSelectPSDvsPos_%i",i),";z [mm];PSD [arb]",numPosBins,posMin,posMax,numPSDBins,PSDMin,PSDMax);
		hBGPSDvsPos[i] 		= new TH2F(Form("hBGPSDvsPos_%i",i),";z [mm];PSD [arb]",numPosBins,posMin,posMax,numPSDBins,PSDMin,PSDMax);

		hSelectEnvsPos[i] 	= new TH2F(Form("hSelectEnvsPos_%i",i),";z [mm];En [MeVee]",numPosBins,posMin,posMax,numEnBins,EnMin,EnMax);
		hBGEnvsPos[i] 		= new TH2F(Form("hBGEnvsPos_%i",i),";z [mm];En [MeVee]",numPosBins,posMin,posMax,numEnBins,EnMin,EnMax);

		hSelectDelayEnvsPromptEn[i] = new TH2F(Form("hSelectDelayEnvsPromptEn_%i",i),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);
		hBGDelayEnvsPromptEn[i]     = new TH2F(Form("hBGDelayEnvsPromptEn_%i",i),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);


	}	//end for loop creating histograms

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Fill histograms
	printf("=============== Filling Histograms =============== \n"); 

//	RNPO *rnpo = new RNPO();
	
	bool exclude;
	Long64_t numEntries = Long64_t(rnpo->fChain->GetEntries());
	printf("Number of Ac-227 candidates: %lld \n",numEntries);

	double tstamp;
	rnpo->GetEntry(0);
	tstamp 		    = rnpo->tstamp;

/*	double setPromptLowPSDCut, setDelayLowPSDCut;
        double setPromptLowEnCut,  setDelayLowEnCut;
        double setDzCutLow, setDzCutHigh;

	double promptLowPSDCut, promptHighPSDCut, promptLowEnCut, promptHighEnCut;
	double delayLowPSDCut,  delayHighPSDCut,  delayLowEnCut,  delayHighEnCut;
	double dzCut, dzLowCut, dzHighCut;

	// Get cut values
	rnpo->GetEntry(0);
	tstamp 		    = rnpo->tstamp;
	setPromptLowPSDCut  = rnpo->p_PSDCut[0]; 
	promptHighPSDCut    = rnpo->p_PSDCut[1];
	setDelayLowPSDCut   = rnpo->d_PSDCut[0];
	delayHighPSDCut     = rnpo->d_PSDCut[1];
	setPromptLowEnCut   = rnpo->p_ECut[0];
	promptHighEnCut     = rnpo->p_ECut[1];
	setDelayLowEnCut    = rnpo->d_ECut[0];
	delayHighEnCut      = rnpo->d_ECut[1];
	dzCut               = rnpo->dzCut;
	setDzLowCut  = -1.0*dzCut;
	setDzHighCut = dzCut;
*/

	//Initialize variables
	int seg;
	double dt, dz;

	double livetime = 0.0;
	double lastTime = 0.0;

	double lastNumClusts = 0.0, numClusts = 0.0;
	double lastRuntime = 0.0, totRuntime = 0.0;

	//===============
	//Variables for keeping track of multiplicity
	int delay_mult = 0, BGdelay_mult = 0;
	int prompt_mult = 0, far_mult = 0;
	int count_mult = 0;
	int larger_mult = 0;

	double multDelayEn = 0 , multDelayPSD = 0;
	double multDelayEnSmear = 0;

	for(Long64_t i=0;i<numEntries;i++){
		if(i%1000000==0) printf("Event: %lld  Events to go: %lld \n",i,numEntries-i);
		rnpo->GetEntry(i);

		if(rnpo->d_t*(1e-6) < (TIMEWINDOW + TIMEOFFSET)) continue; 

		if(rnpo->d_t < lastTime){ 
			livetime += lastTime*(1e-6);		//livetime in ms	
			totRuntime += lastRuntime;
			numClusts += lastNumClusts;
		}
		lastTime = rnpo->d_t;
		lastNumClusts = rnpo->numClust;
		lastRuntime = ((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1();	//[s]	

		//------------------------
		if(count_mult == larger_mult && i>0){
                	if(delay_mult>0){
                        	hSelectDelayPSD[seg]->Fill(multDelayPSD,delay_mult);
                        	hSelectDelayEn[seg]->Fill(multDelayEn,delay_mult);
                        	hSelectDelayEnSmear[seg]->Fill(multDelayEnSmear,delay_mult);
                	}
                	if(BGdelay_mult>0){
                        	hBGDelayPSD[seg]->Fill(multDelayPSD,BGdelay_mult);
                        	hBGDelayEn[seg]->Fill(multDelayEn,BGdelay_mult);
                        	hBGDelayEnSmear[seg]->Fill(multDelayEnSmear,BGdelay_mult);
                	}

                	count_mult=0;
                	delay_mult = 0;
                	BGdelay_mult = 0;
        	}

        	prompt_mult = rnpo->p_mult;
        	far_mult = rnpo->f_mult;
        	larger_mult = (prompt_mult > far_mult) ? prompt_mult : far_mult;

        	count_mult++;

		//------------------------
		seg = rnpo->d_seg;

		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), seg) != end(ExcludeCellArr);
                if(exclude) continue;

		promptLowPSDCut = vPromptPSDCutLow[seg];
		delayLowPSDCut = vDelayPSDCutLow[seg];
		promptLowEnCut = vPromptEnCutLow[seg];
		delayLowEnCut = vDelayEnCutLow[seg];
		dzCutLow = vDzCutLow[seg];
		dzCutHigh = vDzCutHigh[seg];

		//------------------------
		if(rnpo->d_PSD < delayLowPSDCut || rnpo->d_E < delayLowEnCut) continue;
		if(rnpo->d_z < zLow || rnpo->d_z > zHigh) continue;

		//if prompt-delay pair
		dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
		dz = rnpo->d_z - rnpo->p_z;
		if(rnpo->p_seg > -1 && rnpo->p_PSD>promptLowPSDCut && rnpo->p_E>promptLowEnCut && rnpo->p_z>zLow && rnpo->p_z<zHigh && dt>dtCut && dz>dzCutLow && dz<dzCutHigh){
			delay_mult++;
                        multDelayEn = rnpo->d_E;
                        multDelayEnSmear = rnpo->d_ESmear;
                        multDelayPSD = rnpo->d_PSD;

			dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
			dz = rnpo->d_z - rnpo->p_z;

			hSelectDt[seg]->Fill(dt);
			hSelectPromptPSD[seg]->Fill(rnpo->p_PSD);
//			hSelectDelayPSD[seg]->Fill(rnpo->d_PSD);
			hSelectPromptEn[seg]->Fill(rnpo->p_E);
//			hSelectDelayEn[seg]->Fill(rnpo->d_E);
//			hSelectDelayEnSmear[seg]->Fill(rnpo->d_ESmear);
			hSelectPromptTotEn[seg]->Fill(rnpo->p_Etot);	
			hSelectPromptPos[seg]->Fill(rnpo->p_z);
			hSelectDelayPos[seg]->Fill(rnpo->d_z);
			hSelectDz[seg]->Fill(dz);

			hSelectPSDvsEn[seg]->Fill(rnpo->p_E,rnpo->p_PSD);
			hSelectPSDvsEn[seg]->Fill(rnpo->d_E,rnpo->d_PSD);
			hSelectPSDvsPos[seg]->Fill(rnpo->p_z,rnpo->p_PSD);
			hSelectPSDvsPos[seg]->Fill(rnpo->d_z,rnpo->d_PSD);
			hSelectEnvsPos[seg]->Fill(rnpo->p_z,rnpo->p_E);
			hSelectEnvsPos[seg]->Fill(rnpo->d_z,rnpo->d_E);
			hSelectDelayEnvsPromptEn[seg]->Fill(rnpo->p_E,rnpo->d_E);		
		}

		//if prompt-delay BG pair
		dt = (rnpo->d_t - rnpo->f_t)*(1e-6) - TIMEOFFSET;
		dz = rnpo->d_z - rnpo->f_z;
		if(rnpo->f_seg > -1 && rnpo->f_PSD>promptLowPSDCut && rnpo->f_E>promptLowEnCut && rnpo->f_z>zLow && rnpo->f_z<zHigh && dt>dtCut && dz>dzCutLow && dz<dzCutHigh){
			BGdelay_mult++;
                        multDelayEn = rnpo->d_E;
                        multDelayEnSmear = rnpo->d_ESmear;
                        multDelayPSD = rnpo->d_PSD;

			dt = (rnpo->d_t - rnpo->f_t)*(1e-6) - TIMEOFFSET;
			dz = rnpo->d_z - rnpo->f_z;

			hBGDt[seg]->Fill(dt);
			hBGPromptPSD[seg]->Fill(rnpo->f_PSD);
//			hBGDelayPSD[seg]->Fill(rnpo->d_PSD);
			hBGPromptEn[seg]->Fill(rnpo->f_E);
//			hBGDelayEn[seg]->Fill(rnpo->d_E);
//			hBGDelayEnSmear[seg]->Fill(rnpo->d_ESmear);
			hBGPromptTotEn[seg]->Fill(rnpo->f_Etot);
			hBGPromptPos[seg]->Fill(rnpo->f_z);
			hBGDelayPos[seg]->Fill(rnpo->d_z);
			hBGDz[seg]->Fill(dz);

			hBGPSDvsEn[seg]->Fill(rnpo->f_E,rnpo->f_PSD);
			hBGPSDvsEn[seg]->Fill(rnpo->d_E,rnpo->d_PSD);
			hBGPSDvsPos[seg]->Fill(rnpo->f_z,rnpo->f_PSD);
			hBGPSDvsPos[seg]->Fill(rnpo->d_z,rnpo->d_PSD);
			hBGEnvsPos[seg]->Fill(rnpo->f_z,rnpo->f_E);
			hBGEnvsPos[seg]->Fill(rnpo->d_z,rnpo->d_E);
			hBGDelayEnvsPromptEn[seg]->Fill(rnpo->f_E,rnpo->d_E);	
		}
	}	//end for loop over TChain


	livetime += lastTime*(1e-6);	//add time from last tree
	totRuntime += lastRuntime;
	numClusts += lastNumClusts;

	double pileupVetoTime = numClusts*pileupVetoT;	//[ms]

	printf("Total runtime: %f hours \n",totRuntime/(60.0*60.0));
	printf("Livetime: %f hours \n",livetime*(2.778e-7));
	printf("Pileup veto time: %f hours \n",pileupVetoTime*(2.778e-7));

	livetime = livetime - (2.0*pileupVetoTime);
	printf("Corrected Livetime: %f hours \n",livetime*(2.778e-7));
	
	TVectorD *vPileupTime = new TVectorD(1);
	TVectorD *vLivetime = new TVectorD(1);

	vPileupTime[0] = pileupVetoTime;
	vLivetime[0] = livetime;


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Subtract histograms
	printf("=============== Subtracting Histograms =============== \n"); 

	for(int i=0;i<NUMCELLS;i++){
		if(i%10==0) printf("Cell: %d \n",i);

		hSelectDt[i]->Sumw2();
		hRnPoDt[i] = (TH1F*)hSelectDt[i]->Clone();
		hRnPoDt[i]->SetName(Form("hRnPoDt_%i",i));
		hBGDt[i]->Sumw2();
		hRnPoDt[i]->Add(hBGDt[i],-1);

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

		hRnTotEn[i] = (TH1F*)hSelectPromptTotEn[i]->Clone();
		hRnTotEn[i]->SetName(Form("hRnTotEn_%i",i));
		hRnTotEn[i]->Sumw2();
		hRnTotEn[i]->Add(hBGPromptTotEn[i],-1);

		hPoEn[i] = (TH1F*)hSelectDelayEn[i]->Clone();
		hPoEn[i]->SetName(Form("hPoEn_%i",i));
		hPoEn[i]->Sumw2();
		hPoEn[i]->Add(hBGDelayEn[i],-1);
	
		hPoEnSmear[i] = (TH1F*)hSelectDelayEnSmear[i]->Clone();
		hPoEnSmear[i]->SetName(Form("hPoEnSmear_%i",i));
		hPoEnSmear[i]->Sumw2();
		hPoEnSmear[i]->Add(hBGDelayEnSmear[i],-1);

		hRnPos[i] = (TH1F*)hSelectPromptPos[i]->Clone();
		hRnPos[i]->SetName(Form("hRnPos_%i",i));
		hRnPos[i]->Sumw2();
		hRnPos[i]->Add(hBGPromptPos[i],-1);

		hPoPos[i] = (TH1F*)hSelectDelayPos[i]->Clone();
		hPoPos[i]->SetName(Form("hPoPos_%i",i));
		hPoPos[i]->Sumw2();
		hPoPos[i]->Add(hBGDelayPos[i],-1);
			
		hRnPoDz[i] = (TH1F*)hSelectDz[i]->Clone();
		hRnPoDz[i]->SetName(Form("hRnPoDz_%i",i));
		hRnPoDz[i]->Sumw2();
		hRnPoDz[i]->Add(hBGDz[i],-1);

		hRnPoPSDvsEn[i] = (TH2F*)hSelectPSDvsEn[i]->Clone();
		hRnPoPSDvsEn[i]->SetName(Form("hRnPoPSDvsEn_%i",i));
		hRnPoPSDvsEn[i]->Add(hBGPSDvsEn[i],-1);

		hRnPoPSDvsPos[i] = (TH2F*)hSelectPSDvsPos[i]->Clone();
		hRnPoPSDvsPos[i]->SetName(Form("hRnPoPSDvsPos_%i",i));
		hRnPoPSDvsPos[i]->Add(hBGPSDvsPos[i],-1);

		hRnPoEnvsPos[i] = (TH2F*)hSelectEnvsPos[i]->Clone();
		hRnPoEnvsPos[i]->SetName(Form("hRnPoEnvsPos_%i",i));
		hRnPoEnvsPos[i]->Add(hBGEnvsPos[i],-1);

		hPoEnvsRnEn[i] = (TH2F*)hSelectDelayEnvsPromptEn[i]->Clone();
		hPoEnvsRnEn[i]->SetName(Form("hPoEnvsRnEn_%i",i));
		hPoEnvsRnEn[i]->Add(hBGDelayEnvsPromptEn[i],-1);
	}	//end for loop to subtract hists

	histFile->Write();
	histFile->WriteObject(&vLivetime,"vLivetime");
	histFile->WriteObject(&vPileupTime,"vPileupTime");
	histFile->WriteObject(&vPromptPSDCutLow,"vRnPSDCutLow");
	histFile->WriteObject(&vPromptPSDCutHigh,"vRnPSDCutHigh");
	histFile->WriteObject(&vDelayPSDCutLow,"vPoPSDCutLow");
	histFile->WriteObject(&vDelayPSDCutHigh,"vPoPSDCutHigh");
	histFile->WriteObject(&vPromptEnCutLow,"vRnEnCutLow");
	histFile->WriteObject(&vPromptEnCutHigh,"vRnEnCutHigh");
	histFile->WriteObject(&vDelayEnCutLow,"vPoEnCutLow");
	histFile->WriteObject(&vDelayEnCutHigh,"vPoEnCutHigh");
	histFile->WriteObject(&vDzCutLow,"vRnPoDzCutLow");
	histFile->WriteObject(&vDzCutHigh,"vRnPoDzCutHigh");
	histFile->Close();



}	//end void RnPoVsCell

