//Macro containing variables used in RnPoVsCell.C and RnPoVsTime.C


using namespace std;
const int NUMCELLS = 154;
const double POLIFETIME = 2.569;   //[ms]
const double TIMEWINDOW = 5.0*POLIFETIME, TIMEOFFSET = 10.0*POLIFETIME;    //[ms]

//2018C exclude list
//const int NUMEXCLUDECELLS = 64;
//int ExcludeCellArr[NUMEXCLUDECELLS] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 21, 23, 24, 27, 28, 29, 31, 32, 34, 36, 40, 41, 42, 43, 44, 46, 47, 48, 50, 52, 55, 56, 60, 63, 68, 69, 70, 73, 79, 83, 86, 87, 94, 97, 102, 107, 111, 115, 121, 122, 126, 127, 128, 130, 133, 136, 139, 141};

//exclude top of detector
//const int NUMEXCLUDECELLS = 70 + 46;
//int ExcludeCellArr[NUMEXCLUDECELLS] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 21, 23, 24, 27, 28, 29, 31, 32, 34, 36, 40, 41, 42, 43, 44, 46, 47, 48, 50, 52, 55, 56, 60, 63, 68, 69, 70, 73, 79, 83 ,84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153}; 

//exclude bottom of detector
const int NUMEXCLUDECELLS = 84 + 18;
int ExcludeCellArr[NUMEXCLUDECELLS] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 86, 87, 94, 97, 102, 107, 111, 115, 121, 122, 126, 127, 128, 130, 133, 136, 139, 141};


//No ET's
//const int NUMEXCLUDECELLS = 82;
//int ExcludeCellArr[NUMEXCLUDECELLS] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 21, 23, 24, 27, 28, 29, 31, 32, 34, 36, 40, 41, 42, 43, 44, 46, 47, 48, 50, 52, 55, 56, 60, 63, 68, 69, 70, 73, 79, 83, 84, 86, 87, 94, 97, 98, 102, 107, 111, 112, 115, 121, 122, 125, 126, 127, 128, 130, 133, 136, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153};


//Set up bins and ranges for histograms
double dtMin = 0.0, dtMax = 12.9;	//[ms]
//int numDtBins = (dtMax - dtMin)/0.05;
int numDtBins = (dtMax - dtMin)/0.1;
	
double PSDMin = 0.15, PSDMax = 0.37;
int numPSDBins = 100; 

double EnMin = 0.45, EnMax = 1.19;		//[MeVee]
int numEnBins = (EnMax - EnMin)/0.005;	// 5 keV bins
	
double dzMin = -250, dzMax = 250;	//[mm]
int numDzBins = (dzMax - dzMin)/2.5;	//0.25 cm bins

double posMin = -1000, posMax = 1000;	//[mm]
int numPosBins = (posMax - posMin)/10;	//1 cm bins 

double pileupVetoT = 800*(1e-6);	//[ms]



