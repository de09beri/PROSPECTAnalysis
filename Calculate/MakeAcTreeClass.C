#include "TChain.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include "TVectorD.h"
#include "TFile.h"

using namespace std;

int MakeAcTreeClass(){
	clock_t tStart = clock();	
	
	ifstream file;
	file.open("/g/g20/berish1/PROSPECTAnalysis/Calculate/2018C_GoodRuns_RxStatus.txt", ifstream::in);
	
	if(!(file.is_open() && file.good())){
		printf("Good runs file not found. Exiting. \n");
		return -1;
	}


	TChain *chain = new TChain("Ac227TreePlugin/TAc");

	string line;
	string filename;
	int RxFlag;

	int i=0;	
	int numRxOn = 0, numRxOff = 0;
	double RxOnTime = 0, RxOffTime = 0;

	while(file.good() & !file.eof()){
		getline(file, line);

		stringstream s(line);
		s>>filename>>RxFlag;

		TString str = Form("%s/%s/AD1_Extra_Phys.root",gSystem->Getenv("P2X_ANALYZED"),filename.data());

		if(TFile *f = TFile::Open((const char*)str.Data())){
		if(RxFlag==0){ 
			numRxOff++;
			RxOffTime += ((TVectorD*)f->Get("runtime"))->Norm1();
		}
		if(RxFlag==1){ 
			numRxOn++;
			RxOnTime += ((TVectorD*)f->Get("runtime"))->Norm1();
		}
		f->Close();
		}

		chain->Add(str.Data());

		if(i%100==0) printf("At file %i: %s \n",i,(const char*)str);

		file.peek();

		i++;
	}

	int nadded = chain->GetListOfFiles()->GetEntries();
	chain->Lookup(1);
	int nmissing = nadded - chain->GetListOfFiles()->GetEntries();

	printf("===========================================================================\n");
	printf("Attempted to add %d files from good run list to Ac227 TChain \n",nadded);
	printf("Unable to add %d files to TChain \n",nmissing);
	printf("Rx On files: %i, Time: %f s \n",numRxOn,RxOnTime);
	printf("Rx Off files: %i, Time: %f s \n",numRxOff,RxOffTime);
	printf("===========================================================================\n");

	printf("Making class \n");
	chain->MakeClass("RNPO");

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	return 0;
}
