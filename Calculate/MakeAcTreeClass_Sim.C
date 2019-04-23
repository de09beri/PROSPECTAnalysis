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
	
	TChain *chain = new TChain("Ac227TreePlugin/TAc");

	chain->Add("/g/g20/berish1/Simulation/AD_Ac227_FinalRun/AD1_Extra_Phys.root");

	printf("Making class \n");
	chain->MakeClass("RNPO");

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	return 0;
}
