//Macro containing variables used in RnPoVsCell.C and RnPoVsTime.C


using namespace std;
const int NUMCELLS = 154;
const double POLIFETIME = 2.569;   //[ms]
const double TIMEWINDOW = 5.0*POLIFETIME, TIMEOFFSET = 10.0*POLIFETIME;    //[ms]


//const int NUMEXCLUDECELLS = 64;
//int ExcludeCellArr[NUMEXCLUDECELLS] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 21, 23, 24, 27, 28, 29, 31, 32, 34, 36, 40, 41, 42, 43, 44, 46, 47, 48, 50, 52, 55, 56, 60, 63, 68, 69, 70, 73, 79, 83, 86, 87, 94, 97, 102, 107, 111, 115, 121, 122, 126, 127, 128, 130, 133, 136, 139, 141};

//No ET's
const int NUMEXCLUDECELLS = 82;
int ExcludeCellArr[NUMEXCLUDECELLS] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 21, 23, 24, 27, 28, 29, 31, 32, 34, 36, 40, 41, 42, 43, 44, 46, 47, 48, 50, 52, 55, 56, 60, 63, 68, 69, 70, 73, 79, 83, 84, 86, 87, 94, 97, 98, 102, 107, 111, 112, 115, 121, 122, 125, 126, 127, 128, 130, 133, 136, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153};


//Set up bins and ranges for histograms
double dtMin = 0.0, dtMax = 12.9;	//[ms]
//int numDtBins = (dtMax - dtMin)/0.1;
int numDtBins = (dtMax - dtMin)/0.05;
	
double PSDMin = 0.15, PSDMax = 0.37;
int numPSDBins = 100; 

double EnMin = 0.45, EnMax = 1.19;		//[MeVee]
int numEnBins = (EnMax - EnMin)/0.005;	// 5 keV bins
	
double dzMin = -250, dzMax = 250;	//[mm]
int numDzBins = (dzMax - dzMin)/2.5;	//0.25 cm bins

double posMin = -800, posMax = 800;	//[mm]
int numPosBins = (posMax - posMin)/10;	//1 cm bins 

double pileupVetoT = 800*(1e-6);	//[ms]



