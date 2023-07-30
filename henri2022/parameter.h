#include <vector>
#include <array>

const int nBins_full=1;
double ptBins_full[nBins_full+1] = {7,50}; 
double ptBins_fullBP[nBins_full+1] = {5,60};

const int nptBins = 4;
const std::array<double, nptBins + 1> ptbinsvec = { 7, 10, 15, 20, 50};

const int nptBinsBP = 7;
const std::array<double, nptBinsBP + 1> ptbinsvecBP = { 5, 7, 10, 15, 20, 30, 50, 60};

const int nptBinsBP_test  = 1;
const std::array<double, nptBinsBP_test  + 1> ptbinsvecBP_test  = {7,10};

//const unsigned nyBins_both = 8;
//const std::array<double, nyBins_both + 1> ybinsvec = {-2.4, -1.5, -1.0, -0.5, 0.0 ,0.5, 1.0, 1.5, 2.4};

const int nyBins_both = 5;
const std::array< double, nyBins_both + 1> ybinsvec = {0.0 ,0.5, 1.0, 1.5, 2, 2.4};

//const int nyBins_both_full = 1;
//const std::array<double, nyBins_both_full + 1> ybinsvec_full = {-2.4, 2.4};

const int nmBins_both = 7;
const std::array<double, nmBins_both + 1> nmbinsvec = { 0,20,30,40,50,60,70,100};



    //TO BE USED IN Bmesons_Comparison.C 
	//PbPb INFO
	// Bins
	double  BINS_PbPb[4]={7, 10, 15, 20, 50};	
	double  XsecPbPb_X_BP[4]={8.73, 12.4, 17.2, 27.3};
	double  XsecPbPb_X_Bs[4]={8.75, 12.6, 17.4, 27.3};
	double  XsecPbPb_XL_BP[4]={0,0,0,0};
	double  XsecPbPb_XR_BP[4]={0,0,0,0};
	double  XsecPbPb_XL_Bs[4]={0,0,0,0};
	double  XsecPbPb_XR_Bs[4]={0,0,0,0};
	// Bins

	//BP_PbPb
	double  XsecPbPb_Y_BP[4]={311668, 270167, 64384.4, 7704};
	double  XSecPbPb_BP_Y_StatUpRatio[4]={0.159, 0.041, 0.0654, 0.0691};
	double  XSecPbPb_BP_Y_StatDownRatio[4]={0.145, 0.0795, 0.065, 0.0526};
	double  XSecPbPb_BP_Y_SystUpRatio[4]={0.1404, 0.1714, 0.0775, 0.0715};
	double  XSecPbPb_BP_Y_SystDownRatio[4]={0.1359, 0.1705, 0.0761, 0.0698};
	//BP_PbPb
	//Bs_PbPb
	double  XsecPbPb_Y_Bs[4]={160432,75523.7,25354.5,2272.18};	
	double  XSecPbPb_Bs_Y_StatUpRatio[4]={0.513,0.224,0.216,0.216};
	double  XSecPbPb_Bs_Y_StatDownRatio[4]={0.483,0.256,0.207,0.163};
	double  XSecPbPb_Bs_Y_SystUpRatio[4]={0.4564,0.1482,0.1218,0.1647};
	double  XSecPbPb_Bs_Y_SystDownRatio[4]={0.4564,0.1454,0.1210,0.1640};
	//Bs_PbPb
	//PbPb INFO

 	//2015 INFO
	vector<float> vect_ptBins2015bp{7, 10, 15, 20, 30, 50};
	vector<float> vect_ptBins2015bs{7, 15, 20, 50};

	//Bs_2015_pp
	vector<float> vect_BsXsecPPX2015{11,17.5,35.0};
	vector<float> vect_BsXSecPPXErrDown2015{4,2.5,15};
	vector<float> vect_BsXSecPPXErrUp2015{4,2.5,15};
	vector<float> vect_BsXsecPPY2015{316000,34100,3830};
	vector<float> vect_BsXSecPPYErrDown2015{37000,6300,670};
	vector<float> vect_BsXSecPPYErrUp2015{37000,6300,670};
	vector<float> vect_BsXSecPPYSystDown2015{62000,3200,360};
	vector<float> vect_BsXSecPPYSystUp2015{62000,3200,360};
	//Bs_2015_pp
	//BP_2015_pp	
	vector<float> vect_BPXsecPPX2015{8.5,12.5,17.5,25,40};
	vector<float> vect_BPXSecPPXErrDown2015{1.5,2.5,2.5,5,10};
	vector<float> vect_BPXSecPPXErrUp2015{1.5,2.5,2.5,5,10};
	vector<float> vect_BPXsecPPY2015{2610000,744000,197000,46500,5300};
	vector<float> vect_BPXSecPPYErrDown2015{170000,29000,9000,2400,500};
	vector<float> vect_BPXSecPPYErrUp2015{170000,29000,9000,2400,500};
	vector<float> vect_BPXSecPPYSystDown2015{230000,59000,15000,3500,400};
	vector<float> vect_BPXSecPPYSystUp2015{230000,59000,15000,3500,400};
	//BP_2015_pp
 	//2015 INFO

	//TO BE USED IN Bmesons_Comparison.C 


















///////////// FOR MCEFF.C

enum Tracking{
  loose = 0,
  standard = 1,
  tight = 2,
};

std::map<Tracking, double> ptErr{
  {Tracking::standard, 0.1},
  {Tracking::tight, 0.05},
  {Tracking::loose, 0.15}
};

std::map<Tracking, double> chi2Nlayer{
  {Tracking::standard, 0.18},
  {Tracking::tight, 0.15},
  {Tracking::loose, 0.18}
};

///////////// FOR MCEFF.C