#include <vector>
#include <array>

// PbPb published binning
/* const std::vector<double> ptbinsvec = { */
/* 	/\* 5, 7, 10, 15, 20, 50, 60 *\/ */
/* 	7, 10, 15, 20, 50 */
/* }; */

// const unsigned nptBins = 4;
// const std::array<double, nptBins + 1> ptbinsvec = {
// 	7, 10, 15, 20, 50
// };
//


//const unsigned nptBins = 1;
//const std::array<double, nptBins + 1> ptbinsvec = {7, 10};

const unsigned nptBins = 4;
const std::array<double, nptBins + 1> ptbinsvec = { 7, 10, 15, 20, 50};


const unsigned nptBinsBP = 7;
const std::array<double, nptBinsBP + 1> ptbinsvecBP = { 5, 7, 10, 15, 20, 30, 50, 60};

const unsigned nptBinsBP_test  = 1;
const std::array<double, nptBinsBP_test  + 1> ptbinsvecBP_test  = {7,10 };

const unsigned nyBins_both = 8;
const std::array<double, nyBins_both + 1> ybinsvec = {-2.4,-1.5,-1.0,-0.5,0.0 ,0.5, 1.0, 1.5, 2.4};

//const unsigned nyBins_both = 2;
//const std::array<double, nyBins_both + 1> ybinsvec = {0.5, 1.0,1.5};

const unsigned nyBins_both_full = 1;
const std::array<double, nyBins_both_full + 1> ybinsvec_full = {-2.4, 2.4};

const unsigned nmBins_both = 7;
const std::array<double, nmBins_both + 1> nmbinsvec = { 0,20,30,40,50,60,70,100};

//const unsigned nmBins_both = 1;
//const std::array<double, nmBins_both + 1> nmbinsvec = {20,30};

const unsigned nmBins_both_full = 1;
const std::array<double, nmBins_both_full + 1> nmbinsvec_full = {0,100};


// 7 bins
/* float BPXsecPbPbY[NBins] = {4.82132e+06/11.1,311668,270167,64384.4,208537/11.1,28700.6/11.1,7000.73/11.1}; */
/* float BPXsecPbPbX[NBins] = {6,8.73,12.4,17.2,25,40,55}; */
/* float BPXSecPbPbXErrUp[NBins] = {1,1.27,2.6,2.8,5,10,5}; */
/* float BPXSecPbPbXErrDown[NBins] = {1,1.23,2.4,2.2,5,10,5}; */
/* float BPXSecPbPbYErrUpRatio[NBins] = {0.278198,0.159,0.041,0.0654,0.0690334,0.104543,0.24575}; */
/* float BPXSecPbPbYErrDownRatio[NBins] = {0.278198,0.145,0.0795,0.065,0.0690334,0.104543,0.24575}; */

const unsigned nptBinsPbPb = 4;
// const unsigned nptBinsPbPb = nptBins;
float BPXsecPbPbY[nptBinsPbPb] = {311668, 270167, 64384.4, 7704};
float BPXsecPbPbX[nptBinsPbPb] = {8.73,12.4,17.2,27.3};
float BPXSecPbPbXErrUp[nptBinsPbPb] = {1.27,2.6,2.8,22.7};
float BPXSecPbPbXErrDown[nptBinsPbPb] = {1.23,2.4,2.2,7.3};

float BPXSecPbPbYErrUpRatio[nptBinsPbPb] =   {0.159, 0.041, 0.0654, 0.0691};
float BPXSecPbPbYErrDownRatio[nptBinsPbPb] = {0.145, 0.0795, 0.065, 0.0526};
float BPXSecPbPbYSystUpRatio[nptBinsPbPb]   = {0.1404, 0.1714, 0.0775, 0.0715};
float BPXSecPbPbYSystDownRatio[nptBinsPbPb] = {0.1359, 0.1705, 0.0761, 0.0698};









    //TO BE USED IN Bmesons_Comparison.C (THESE ARE LEADLEAD INFO)

	 vector<float> vect_BPXsecPbPbY = {4.82132e+06/11.1,311668,270167,64384.4,208537/11.1,28700.6/11.1,7000.73/11.1};
	 vector<float> vect_BPXsecPbPbX = {6,8.73,12.4,17.2,25,40,55};
	 vector<float> vect_BPXSecPbPbXErrUp = {1,1.27,2.6,2.8,5,10,5};
	 vector<float> vect_BPXSecPbPbXErrDown = {1,1.23,2.4,2.2,5,10,5};
	 vector<float> vect_BPXSecPbPbYErrUpRatio = {0.278198,0.159,0.041,0.0654,0.0690334,0.104543,0.24575};
	 vector<float> vect_BPXSecPbPbYErrDownRatio = {0.278198,0.145,0.0795,0.065,0.0690334,0.104543,0.24575};
	 vector<float> vect_BPXSecPbPbYSystUpRatio = {0.3577,0.1404,0.1714,0.0775,0.0858,0.0715,0.1253};
	 vector<float> vect_BPXSecPbPbYSystDownRatio = {0.3210,0.1359,0.1705,0.0761,0.0843,0.0699,0.1220};
