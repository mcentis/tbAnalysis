#include "iostream"
#include "fstream"
#include "stdlib.h"
#include "string"
#include "vector"

#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TProfile.h"
#include "TSpectrum.h"
#include "TSpline.h"
#include "TLegend.h"
#include "TMultiGraph.h"

#include "langauFit.hh"

#include "ConfigFileReader.hh"
#include "trackDef.h"

int findRunNumber(std::string in)
{
  int slashPos = in.find_last_of("/");
  std::string name = in.substr(slashPos + 1);

  char num[7];
  for(int i = 0; i < 6; ++i)
    num[i] = name[i];

  return atoi(num);
}

TF1* fitResiduals(TH1* inHist, double range1, double range2)
{
  double fr1 = inHist->GetMean() - inHist->GetRMS();
  double fr2 = inHist->GetMean() + inHist->GetRMS();
  TF1* gausFit = new TF1("gausFit", "gaus", fr1, fr2); // fit with gauss
  inHist->Fit(gausFit, "RQ");

  TF1* fitFunc = new TF1("fitFunc", "gaus + [3]", range1, range2);
  fitFunc->SetParNames("Const", "Mean", "Sigma", "Offset");
  fitFunc->SetLineColor(kRed);
  fitFunc->SetParameter(0, gausFit->GetParameter(0)); // use gauss parameter as start values
  fitFunc->SetParameter(1, gausFit->GetParameter(1));
  fitFunc->SetParameter(2, gausFit->GetParameter(2));
  fitFunc->SetParameter(3, 0); // no offset

  fitFunc->SetParLimits(0, 0, inHist->GetEntries());
  fitFunc->SetParLimits(1, inHist->GetMean() - inHist->GetRMS(), inHist->GetMean() + inHist->GetRMS());
  fitFunc->SetParLimits(2, 0, 2 * inHist->GetRMS());
  fitFunc->SetParLimits(3, 0, inHist->GetEntries());

  inHist->Fit(fitFunc, "RQ");

  return fitFunc;
}

std::vector<int> readGoodChFile(const char* chFileName)
{
  std::vector<int> goodChannels;
  const int nChannels = 256; // channels of the readout system

  std::ifstream chStr;
  chStr.open(chFileName, std::ifstream::in);

  if(chStr.is_open() == false)
    {
      std::cout << "Impossible to open the file " << chFileName << "\nAll channels will be used" << std::endl;
      for(int i = 0; i < nChannels; i++) goodChannels.push_back(i);
      return goodChannels;
    }

  std::string line;
  std::string::iterator newEnd; // end of the string after remove
  unsigned int sharpPos; // position of the #

  while(!chStr.eof())
    {
      line.clear();

      std::getline(chStr, line);

      if(line.length() == 0) continue; // empty line

      sharpPos = line.find('#');
      if(sharpPos == 0) continue; // full comment line, skip

      if(sharpPos < line.size()) // sharp found
      line.resize(sharpPos); // ignore what comes after the #

      // removes all the spaces from the string, remove does not change the length of the string, it is necessary to resize the string
      newEnd = std::remove(line.begin(), line.end(), ' ');
      if(newEnd == line.begin()) continue; // string of spaces and comments
      line.resize(newEnd - line.begin()); // resize the string to its new size

      // same treatment for the \t
      newEnd = std::remove(line.begin(), line.end(), '\t');
      if(newEnd == line.begin()) continue; // string of spaces, tabs and comments
      line.resize(newEnd - line.begin()); // resize the string to its new size

      goodChannels.push_back(atoi(line.c_str()));
    }

  chStr.close();

  return goodChannels;
}

TF1* fitPeak(TH1* hist, double negSig, double posSig)
{
  TF1* fitFunc = new TF1("fitFunc", "gaus");
  if(hist == 0)
    {
      fitFunc->SetParameters(0, 0, 0);
      return fitFunc;
    }

  TCanvas* fitCan = new TCanvas("fitCa"); // it is important that the canvases have different names, strange root...

  TSpectrum* spec = new TSpectrum(10, 5); // find the peaks (max 4) in the histo  
  int nPeaks = spec->Search(hist, 3, "nobackground new"); // sigma peak = 15, no backgrund considered

  Double_t* pos = spec->GetPositionX();
  Double_t* heigth = spec->GetPositionY();

  double maxPos = -1;
  double startHeigth = -1;
  for(int i = 0; i < nPeaks; ++i)
    {
      if(heigth[i] < 10) continue;
      if(pos[i] > maxPos)
	{
	  maxPos = pos[i];
	  startHeigth = heigth[i];
	}
    }

  fitFunc->SetParameter(0, startHeigth);
  fitFunc->SetParameter(1, maxPos);

  fitFunc->SetParameter(2, 15);
  double start = maxPos - 20;
  double stop = maxPos + 20;
  fitFunc->SetRange(start, stop);
  hist->Fit(fitFunc, "RQ");

  start = fitFunc->GetParameter(1) - negSig * fitFunc->GetParameter(2);
  stop = fitFunc->GetParameter(1) + posSig * fitFunc->GetParameter(2);
  fitFunc->SetRange(start, stop);
  hist->Fit(fitFunc, "RQ");

  delete fitCan;

  return fitFunc;
}

// global variables to be used in the findfracpos
TSpline3* sp;
TF1* level;

double toMin(double* x, double* par)
{
  return fabs(level->Eval(*x) - sp->Eval(*x)); 
}

struct findFracPos // different from the one in mergeresults
{
  TF1* func;
  double fracPos;
  double errHi;
  double errLo;
  double deltaEntries;

  findFracPos(TH1D* inHist, double frac)
  {
    TH1D* inte = new TH1D(*inHist); // CDF histogram
    int binStart = 1;
    int binStop = inHist->GetXaxis()->GetNbins();
    double totArea = inHist->Integral(binStart, binStop);
    double intBin = 0;

    for(int i = binStart; i < binStop + 1; ++i)
      {
	intBin = inHist->Integral(binStart, i);
	inte->SetBinContent(i, intBin);
      }

    sp = new TSpline3(inte); // sp is a global variable
    level = new TF1("level", "[0]");
    double xmin = inte->GetXaxis()->GetXmin();
    double xmax = inte->GetXaxis()->GetXmax();
    level->SetRange(xmin, xmax);
    level->SetParameter(0, frac * totArea);
    func = new TF1("func", toMin, xmin, xmax, 0);
    fracPos = func->GetMinimumX();
    deltaEntries = sqrt(totArea * frac * (1 - frac)); // std dev binomial
    level->SetParameter(0, totArea * frac + deltaEntries);
    errHi = func->GetMinimumX() - fracPos;
    level->SetParameter(0, totArea * frac - deltaEntries);
    errLo = fracPos - func->GetMinimumX();
    delete sp;
    delete inte;
  }

  ~findFracPos()
  {
    delete func;
  }
};

double Median(const TH1D * h1) {
  int n = h1->GetXaxis()->GetNbins(); 
  std::vector<double>  x(n);
  h1->GetXaxis()->GetCenter( &x[0] );
  const double * y = h1->GetArray();
  // exclude underflow/overflows from bin content array y
  return TMath::Median(n, &x[0], &y[1]);
}

int main(int argc, char* argv[])
{
  if(argc != 3)
    {
      std::cout << "Usage: chargeDistr inFile confFile" << std::endl;
      return 1;
    }

  gSystem->Load("libTree"); // necessary to use the program, root magic

  // big axis labels
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");

  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");

  gStyle->SetTitleOffset(0.95, "x");
  gStyle->SetTitleOffset(0.95, "y");

  TFile* inFile = TFile::Open(argv[1]);
  if (!inFile)
    {
      std::cout << "Impossible to open " << argv[1] << std::endl;
      return -1;
    }

  const int nChannels = 256; // channels of the readout system
  const int nPlanes = 7; // planes (6 telescope + 1 dut)
  const double pitch = 0.080; // mm
  char name[50]; // to be used in varios occasions
  char title[200];

  TDirectory* dir = (TDirectory*) inFile->Get("Ntuple");
  if(!dir)
    {
      std::cout << "TDirectory not found" << std::endl;
      return -1;
    }

  TTree* trkTree = (TTree*) dir->Get("EUFit");
  if(!dir)
    {
      std::cout << "TTree not found" << std::endl;
      return -1;
    }

  ConfigFileReader* conf = new ConfigFileReader(argv[2]);
  //conf->DumpConfMap();

  TString outFileName = conf->GetValue("outFilePath");
  outFileName += findRunNumber(argv[1]);
  outFileName += ".root";

  const int maxDist = atoi(conf->GetValue("hitMaxDist").c_str()); // maximum distance form the strip hit by a track in looking at the charge, in number of strips
  const float timeCut1 = atof(conf->GetValue("timeCut1").c_str()); // time cut for the events
  const float timeCut2 = atof(conf->GetValue("timeCut2").c_str()); // time cut for the events
  const float xCut1 = atof(conf->GetValue("xCut1").c_str()); // cut in the x direction, applied on the telescope ref frame
  const float xCut2 = atof(conf->GetValue("xCut2").c_str()); // cut in the x direction, applied on the telescope ref frame
  int readenPolarity = atoi(conf->GetValue("polarity").c_str()); // expected signal polarity
  const int polarity = readenPolarity / abs(readenPolarity);
  const float negSigmaFit = atof(conf->GetValue("negSigmaFit").c_str()); // number of sigma to which extend the landau gaussian fit in the negative direction
  const float posSigmaFit = atof(conf->GetValue("posSigmaFit").c_str()); // number of sigma to which extend the landau gaussian fit in the positive direction
  int maxEntryNum = atoi(conf->GetValue("maxEntryNum").c_str()); // max entry number to be processed
  int highestNeighbor = atoi(conf->GetValue("highestNeighbor").c_str()); // if not 0 only the highest neighbor is added to the hit charge

  double scaleFactor; // used in the track loop to apply temperature correction
  // if(atof(conf->GetValue("scaleFactor").c_str()) != 0)
  //   {
  //     scaleFactor = atof(conf->GetValue("scaleFactor").c_str());
  //     std::cout << "================================================================>>> WARNING You are going to apply a scaling factor to the signal, the single channel noise will not be scaled" << std::endl;
  //   }
  // else
  //   scaleFactor = 1;

  std::vector<int> goodChVec = readGoodChFile(conf->GetValue("goodChFile").c_str()); // good channels

  double ADCtoe = atof(conf->GetValue("ADCtoe").c_str()); // parameter to apply conversion of ADC to e-
  double ADCtoeErr = atof(conf->GetValue("ADCtoeErr").c_str()); // parameter to apply conversion of ADC to e-

  double tCorr_p0 = atof(conf->GetValue("tCorr_p0").c_str()); // parameters to apply temperature correction
  double tCorr_p1 = atof(conf->GetValue("tCorr_p1").c_str());
  double tCorr_p0Err = atof(conf->GetValue("tCorr_p0Err").c_str());
  double tCorr_p1Err = atof(conf->GetValue("tCorr_p1Err").c_str());
  double tCorr_p0p1Cov = atof(conf->GetValue("tCorr_p0p1Cov").c_str());
  double tCorr_p2 = atof(conf->GetValue("tCorr_p2").c_str());
  double tempErr = atof(conf->GetValue("tempErr").c_str());
  double targetChipTemp = atof(conf->GetValue("targetChipTemp").c_str());
  double targetGain;
  double targetGainErr;
  bool applyTcorr = false; // flag used in the loop on the tracks

  if(tCorr_p0 || tCorr_p1 || tCorr_p2)
    {
      targetGain = tCorr_p0 + tCorr_p1 * targetChipTemp + tCorr_p2 * targetChipTemp * targetChipTemp;
      targetGainErr = sqrt(pow(tCorr_p0Err, 2) + pow(tCorr_p1Err * targetChipTemp, 2) + 2 * tCorr_p0p1Cov * targetChipTemp);

      applyTcorr = true;

      std::cout << "=============================>>> WARNING a temperature correction is going to be applied to the signal, the single channel noise will not be corrected\n";
      std::cout << "Target chip temperature " << targetChipTemp << std::endl;
      std::cout << "Target gain " << targetGain << std::endl;
    }

  TFile* outFile = new TFile(outFileName, "RECREATE");

  //Declaration of leaves types
  Int_t           Event;
  Int_t           RunNr;
  Int_t           EvtNr;
  Int_t           Ndf;
  Float_t         Chi2;
  Double_t        measX[nPlanes];
  Double_t        measY[nPlanes];
  Double_t        measZ[nPlanes];
  Double_t        measQ[nPlanes];
  Double_t        fitX[nPlanes];
  Double_t        fitY[nPlanes];

  Double_t        dutTrackX_global; // extrapolation of the track on the DUT, global ref frame
  Double_t        dutTrackY_global;
  Double_t        dutTrackX_local;  // extrapolation of the track on the DUT, local ref frame
  Double_t        dutTrackY_local;
  Double_t        dutTrackX_pixel; // extrapolation of the track on the DUT, in pixels / strips
  Double_t        dutTrackY_pixel;
  Double_t        dutHitX_global; // measured DUT hit in global ref frame
  Double_t        dutHitY_global;
  Double_t        dutHitX_local;
  Double_t        dutHitY_local;
  Double_t        dutHitX_pixel;
  Double_t        dutHitY_pixel;
  Double_t        dutHitR; // distance between extrapolated and measured hit on DUT in global ref frame, mm
  Double_t        dutHitQ; // charge of the cluster hit on DUT, units of ADC * 100
  Float_t         alibava_TDC;
  Float_t         alibava_temp;
  Double_t        alibavaPH[nChannels] = {0.};

  // Set branch addresses.
  trkTree->SetBranchAddress("Event",&Event);
  trkTree->SetBranchAddress("RunNr",&RunNr);
  trkTree->SetBranchAddress("EvtNr",&EvtNr);
  trkTree->SetBranchAddress("Ndf",&Ndf);
  trkTree->SetBranchAddress("Chi2",&Chi2);

  for(int i = 0; i < nPlanes; ++i)
    {
      trkTree->SetBranchAddress(TString::Format("measX_%i", i), &measX[i]);
      trkTree->SetBranchAddress(TString::Format("measY_%i", i), &measY[i]);
      trkTree->SetBranchAddress(TString::Format("measZ_%i", i), &measZ[i]);
      trkTree->SetBranchAddress(TString::Format("measQ_%i", i), &measQ[i]);
      trkTree->SetBranchAddress(TString::Format("fitX_%i", i), &fitX[i]);
      trkTree->SetBranchAddress(TString::Format("fitY_%i", i), &fitY[i]);
    }

  trkTree->SetBranchAddress("dutTrackX_global",&dutTrackX_global); // position extrapolated from the track fit on the DUT, global ref frame
  trkTree->SetBranchAddress("dutTrackY_global",&dutTrackY_global);
  trkTree->SetBranchAddress("dutTrackX_local",&dutTrackX_local); // position extrapolated from the track fit on the DUT, local ref frame
  trkTree->SetBranchAddress("dutTrackY_local",&dutTrackY_local);
  trkTree->SetBranchAddress("dutTrackX_pixel",&dutTrackX_pixel); // position extrapolated on the DUT, dut ref frame, in pixel / strip number
  trkTree->SetBranchAddress("dutTrackY_pixel",&dutTrackY_pixel);
  trkTree->SetBranchAddress("dutHitX_global",&dutHitX_global); // cluster matched on the alibava
  trkTree->SetBranchAddress("dutHitY_global",&dutHitY_global);
  trkTree->SetBranchAddress("dutHitX_local",&dutHitX_local); // cluster matched on the alibava
  trkTree->SetBranchAddress("dutHitY_local",&dutHitY_local);
  trkTree->SetBranchAddress("dutHitX_pixel",&dutHitX_pixel); // cluster matched on the alibava
  trkTree->SetBranchAddress("dutHitY_pixel",&dutHitY_pixel);
  trkTree->SetBranchAddress("dutHitR",&dutHitR);
  trkTree->SetBranchAddress("dutHitQ",&dutHitQ);
  trkTree->SetBranchAddress("alibava_tdc",&alibava_TDC);
  trkTree->SetBranchAddress("alibava_temp",&alibava_temp);

  for(std::vector<int>::iterator ch = goodChVec.begin(); ch != goodChVec.end(); ch++) // associate just the good channels
    {
      sprintf(name, "alibava_reco_ch_%i", *ch);
      trkTree->SetBranchAddress(name,&alibavaPH[*ch]);
    }

  // activate the interesting branches
  trkTree->SetBranchStatus("*", 0);
  trkTree->SetBranchStatus("EvtNr", 1);
  trkTree->SetBranchStatus("alibava*", 1); // all the alibava info
  trkTree->SetBranchStatus("dutTrackX_*", 1);
  trkTree->SetBranchStatus("dutTrackY_*", 1);
  trkTree->SetBranchStatus("dutHitX_global", 1);
  trkTree->SetBranchStatus("dutHitY_global", 1);

  // tracks
  TH1I* trkEvt = new TH1I("traksEvt", "Number of tracks per event;Number of tracks;Entries", 21, -0.5, 20.5);
  TGraph* trkVsEvt = new TGraph();
  trkVsEvt->SetName("trkVsEvt");
  trkVsEvt->SetTitle("Number of tracks vs event");
  trkVsEvt->SetMarkerStyle(7);
  TH1I* trkEvtSelected = new TH1I("traksEvtSelected", "Number of tracks per event in events that pass the event selection;Number of tracks;Entries", 11, -0.5, 10.5);

  // check of time in the reconstruction
  TH1D* timeTrack = new TH1D("timeTrack", "Time of the reconstructed tracks;Time [ns];Entries", 60, 0, 120);
  TH2D* timeTrackPosY = new TH2D("timeTrackPosY", "Time and position in Y of the reconstructed tracks;Time [ns];y [mm]", 60, 0, 120, 100, -10, 10);
  TH2D* adcTime = new TH2D("adcTime", "ADC PH vs time for all the channels;Time [ns];PH [ADC]", 60, 0, 120, 1024, -511.5, 511.5);

  // residuals
  TH1D* residualsX = new TH1D("residualsX", "Difference between matched cluster and extrapolated position along x;x_{matched} - x_{extrapolated} [mm];Entries", 501, -1.5, 1.5);
  TH1D* residualsY = new TH1D("residualsY", "Difference between matched cluster and extrapolated position along y;y_{matched} - y_{extrapolated} [mm];Entries", 501, -1.5, 1.5);
  TH2D* residuals = new TH2D("residuals", "Difference between matched cluster and extrapolated position;x_{matched} - x_{extrapolated} [mm];y_{matched} - y_{extrapolated} [mm];Entries", 501, -1.5, 1.5, 501, -1.5, 1.5);
  TH2D* residualsXvsX = new TH2D("residualsXvsX", "Difference between matched cluster and extrapolated position along x as function of x;x_{extrapolated} [mm];x_{matched} - x_{extrapolated} [mm];Entries", 501, -20, 20, 501, -1.5, 1.5);
  TH2D* residualsYvsY = new TH2D("residualsYvsY", "Difference between matched cluster and extrapolated position along y as function of y;y_{extrapolated} [mm];y_{matched} - y_{extrapolated} [mm];Entries", 501, -10, 10, 501, -1.5, 1.5);
  TH2D* residualsYvsX = new TH2D("residualsYvsX", "Difference between matched cluster and extrapolated position along y as function of x;x_{extrapolated} [mm];y_{matched} - y_{extrapolated} [mm];Entries", 501, -20, 20, 501, -1.5, 1.5);
  TH2D* residualsXvsY = new TH2D("residualsXvsY", "Difference between matched cluster and extrapolated position along x as function of y;y_{extrapolated} [mm];x_{matched} - x_{extrapolated} [mm];Entries", 501, -10, 10, 501, -1.5, 1.5);
  TH2D* residualsYvsEvt = new TH2D("residualsYvsEvt", "Difference between matched cluster and extrapolated position along y as function of event number;Event number;y_{matched} - y_{extrapolated} [mm];Entries", 500, 0, 5e5, 501, -1.5, 1.5);
  TH2D* residualsYvsEntry = new TH2D("residualsYvsEntry", "Difference between matched cluster and extrapolated position along y as function of entry number;Entry number;y_{matched} - y_{extrapolated} [mm];Entries", 5000, 0, 5e6, 501, -1.5, 1.5);
  TH1D* residualsYselected = new TH1D("residualsYselected", "Difference between matched cluster and extrapolated position along y in events that pass the event selection;y_{matched} - y_{extrapolated} [mm];Entries", 501, -1.5, 1.5);

  // hitmaps
  TH2D* hitMapDUTtele = new TH2D("hitMapDUTtele", "Extrapolated position of the tracks on the strip sensor;x [mm];y [mm]", 200, -20, 20, 100, -10, 10);
  TH2D* hitMapMatched = new TH2D("hitMapMatched", "Matched hits on the strip sensor;x [mm];y [mm]", 200, -20, 20, 100, -10, 10);
  TH2D* hitMapDUTgoodCh = new TH2D("hitMapDUTgoodCh", "Extrapolated position of the tracks on the strip sensor, passing a good channel;x [mm];y [mm]", 200, -20, 20, 100, -10, 10);
  TH1I* extraChDistr = new TH1I("extraChDistr", "Distribution of the extrapolated position in channels;Channel;Entries", 513, -255.5, 255.5);
  TH1I* extraChDistrGoodCh = new TH1I("extraChDistrGoodCh", "Distribution of the extrapolated position in channels (only good channels shown);Channel;Entries", 256, -0.5, 255.5);

  // low ph entries investigation
  TH2D* hitMapLowPH = new TH2D("hitMapLowPH", "Position of tracks with a signal of less than 15 in the time cut;x [mm];y [mm]", 200, -20, 20, 100, -10, 10);
  TH1I* trkEvtLowPH = new TH1I("traksEvtLowPH", "Number of tracks per event in events that have a signal of less than 15 in the time cut;Number of tracks;Entries", 11, -0.5, 10.5);
  TH1D* stripHPHDiffExtraLowPH = new TH1D("stripHPHDiffExtraLowPH", "Difference in strip number between extracted and highest PH strip in the time cut, events with PH < 15;ExtraStr - HiPHSt  [Strip];Entries", 21, -10.5, 10.5);
  TH2D* diffExtraStripHPHvsPH = new TH2D("diffExtraStripHPHvsPH", "Difference extrapolated strip and strip HPH vs hit PH;Hit PH [ADC];ExtraStr - HiPHSt  [Strip];Entries",  562, -50.5, 511.5, 21, -10.5, 10.5);

  // signal and noise
  TH2D* signalTime = new TH2D("signalTime", "Hit signal vs time;Time [ns];Hit signal [ADC]", 60, 0, 120, 1024, -511.5, 511.5);
  TH2D* positivizedSignalTime = new TH2D("positivizedSignalTime", "Hit signal vs time (positivized);Time [ns];Hit signal [ADC]", 60, 0, 120, 151, -50.5, 511.5);
  TH1D* signalDistr = new TH1D("signalDistr", "Hit signal distribution (positivized);Hit signal[ADC];Entries", 562, -50.5, 511.5);
  TH1D* noiseDistr = new TH1D("noiseDistr", "Signal distribution (positivized) not associated with a hit;Signal [ADC];Entries", 201, -100.5, 100.5);
  sprintf(title, "Signal distribution (positivized) not associated with a hit, summed over %i channels ;Signal [ADC];Entries", maxDist * 2 + 1);
  TH1D* noiseDistrGroup = new TH1D("noiseDistrGroup", title, 201, -100.5, 100.5);
  sprintf(title, "Signal distribution (positivized) not associated with a hit, summed over %i channels ;Signal [e^{-}];Entries", maxDist * 2 + 1);
  TH1D* noiseDistrGroup_electrons = new TH1D("noiseDistrGroup_electrons", title, 201, -100.5 * ADCtoe, 100.5 * ADCtoe);
  TH1D* noiseDistrPair = new TH1D("noiseDistrPair", "Signal distribution (positivized) not associated with a hit, summed over 2 channels ;Signal [ADC];Entries", 201, -100.5, 100.5);

  TH1D* signalDistrTimeCut = new TH1D("signalDistrTimeCut", "Hit signal distribution (positivized) in the time cut;Hit signal[ADC];Entries", 500, -50.5, 600.5);
  TH1D* signalDistrTimeCutDistCutEta01 = new TH1D("signalDistrTimeCutDistCutEta01", "Hit signal distribution (positivized), time cut, dist cut, #eta > 1 or #eta < 0;Hit signal[ADC];Entries", 500, -50.5, 600.5);
  TH1D* signalDistrTimeCutDistCut = new TH1D("signalDistrTimeCutDistCut", "Hit signal distribution (positivized) in the time cut, highest PH strip neighboring the extrapolated one;Hit signal [ADC];Entries", 200, -50.5, 600.5);
  TH1D* signalDistrTimeCutDistCut_electrons = new TH1D("signalDistrTimeCutDistCut_electrons", "Hit signal distribution (positivized) in the time cut, highest PH strip neighboring the extrapolated one;Hit signal [e^{-}];Entries", 200, -50.5 * ADCtoe, 600.5 * ADCtoe);
  TH1D* signalDistrTimeCutDistCut_noisePeakSub = new TH1D(*signalDistrTimeCutDistCut); // to preserve binning
  signalDistrTimeCutDistCut_noisePeakSub->SetName("signalDistrTimeCutDistCut_noisePeakSub");
  signalDistrTimeCutDistCut_noisePeakSub->SetTitle("Hit signal distribution (positivized) in the time cut, highest PH strip neighboring the extrapolated one, noise peak subtracted;Hit signal [ADC];Entries");
  TH1D* signalDistrTimeCutDistCut_noisePeakSub_electrons = new TH1D(*signalDistrTimeCutDistCut_electrons); // to preserve binning
  signalDistrTimeCutDistCut_noisePeakSub_electrons->SetName("signalDistrTimeCutDistCut_noisePeakSub_electrons");
  signalDistrTimeCutDistCut_noisePeakSub_electrons->SetTitle("Hit signal distribution (positivized) in the time cut, highest PH strip neighboring the extrapolated one, noise peak subtracted;Hit signal [e^{-}];Entries");
  TH1D* noiseDistrTimeCut = new TH1D("noiseDistrTimeCut", "Signal distribution (positivized) not associated with a hit in the time cut;Signal [ADC];Entries", 201, -100.5, 100.5);
  TH1D* signalDistrTimeDistHPHcut = new TH1D("signalDistrTimeDistHPHcut", "Hit signal distribution (positivized), time cut, distance cut, strip HPH cut;Hit signal[ADC];Entries", 500, -50.5, 511.5);

  // signal on the strip with highest ph
  TH1D* stripHPHDistrTimeCut = new TH1D("stripHPHDistrTimeCut", "Hit signal distribution (positivized) in the time cut, for the strip with highest charge;Hit signal[ADC];Entries", 500, -50.5, 511.5);
  TH1D* stripHPHDistrTimeCutDistCut = new TH1D("stripHPHDistrTimeCutDistCut", "Hit signal distribution (positivized) in the time cut, for the strip with highest charge, highest PH strip neighboring the extrapolated one;Hit signal[ADC];Entries", 500, -50.5, 511.5);
  TH1D* backgroundDistrHPH = new TH1D("backgrounDistrHPH", "Backgorund distribution of the strip with highest PH;Signal [ADC];Entries", 500, -50.5, 511.5);
  TH1D* stripHPHDistrTimeCutDistCut_electrons = new TH1D("stripHPHDistrTimeCutDistCut_electrons", "Hit signal distribution (positivized) in the time cut, for the strip with highest charge, highest PH strip neighboring the extrapolated one;Hit signal[e^{-}];Entries", 500, -50.5 * ADCtoe, 511.5 * ADCtoe);
  TH1D* backgroundDistrHPH_electrons = new TH1D("backgrounDistrHPH_electrons", "Backgorund distribution of the strip with highest PH;Signal [e^{-}];Entries", 500, -50.5 * ADCtoe, 511.5 * ADCtoe);
  TH1D* stripHPH_plusNeigh_DistrTimeCutDistCut = new TH1D("stripHPH_plusNeigh_DistrTimeCutDistCut", "Hit signal distribution (positivized) in the time cut, for the strip with highest charge plus its highest neighbor, highest PH strip neighboring the extrapolated one;Hit signal[ADC];Entries", 500, -50.5, 511.5);
  TH1D* stripHPHDiffExtra = new TH1D("stripHPHDiffExtra", "Difference in strip number between extracted and highest PH strip in the time cut;ExtraStr - HiPHSt  [Strip];Entries", 21, -10.5, 10.5);
  TH2D* phAroundHPHstripTimeCut = new TH2D("phAroundHPHstripTimeCut", "PH of the hit centered on the strip with highest PH, in the time cut;Strip;PH [ADC]", 21, -10.5, 10.5, 562, -50.5, 511.5);
  TH2D* phAroundExtraStripTimeCut = new TH2D("phAroundExtraStripTimeCut", "PH of the hit centered on the extrapolated strip, in the time cut;Strip;PH [ADC]", 21, -10.5, 10.5, 562, -50.5, 511.5);
  TH2D* stripHPHSignalTime = new TH2D("stripHPHSignalTime", "Signal vs time (positivized), for the strip with highest charge;Time [ns];Hit signal [ADC]", 60, 0, 120, 562, -150.5, 411.5);
  TH2D* leftStripHPHSignalTime = new TH2D("leftStripHPHSignalTime", "Signal vs time (positivized), for the left strip neighboring the one with highest charge;Time [ns];Hit signal [ADC]", 60, 0, 120, 562, -150.5, 411.5);
  TH2D* rightStripHPHSignalTime = new TH2D("rightStripHPHSignalTime", "Signal vs time (positivized), for the right strip neighboring the one  with highest charge;Time [ns];Hit signal [ADC]", 60, 0, 120, 562, -150.5, 411.5);
  TH2D* correlationPHstripHPHhit = new TH2D("correlationPHstripHPHhit", "Correlation between PH on the strip Hi PH and the hit PH, time cut and dist cut;PH tot hit [ADC];PH strip Hi PH [ADC]", 562, -50.5, 511.5, 562, -50.5, 511.5);
 
  // chip temperature
  TGraph* tempEvt = new TGraph();
  tempEvt->SetName("tempEvt");
  tempEvt->SetTitle("Tempetrature of the beetle chip vs event number");

  TH1D* tempDistr = new TH1D("tempDistr", "Temperature of the beetle chip;Temperature [C];Events", 600, 0, 30);

  // temperature correction
  TH1D* scaleFactorDistr = new TH1D("scaleFactorDistr", "Factor for temperature correction;Factor;Entries", 4000, 0, 2);

  // three histos to get a map of charge collection over 2 strips in the time cut, distance cut
  double minX = -20000;
  double maxX = 20000;
  int binX = 100;
  double minY = -0.5;
  double maxY = 160.5;
  int binY = 31;
  TH2D* chargeMapMod160Normalized = new TH2D("chargeMapMod160Normalized", "Charge map in the time cut, distance cut, divided by the number of tracks;x [#mum];y mod 160 [#mum];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);
  TH2D* chargeMapMod160 = new TH2D("chargeMapMod160", "Charge map in the time cut, distance cut;x [#mum];y mod 160 [#mum];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);
  TH2D* hitMapMod160 = new TH2D("hitMapMod160", "Hit map in the time cut, distance cut;x [#mum];y mod 160 [#mum];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);

  // charge collection map over the sensor in the time cut, distance cut
  minX = -20;
  maxX = 20;
  binX = 200;
  minY = -10;
  maxY = 10;
  binY = 100;
  TH2D* chargeMapNormalized = new TH2D("chargeMapNormalized", "Charge map in the time cut, distance cut, divided by the number of tracks;x [mm];y [mm];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);
  TH2D* chargeMap = new TH2D("chargeMap", "Charge map in the time cut, distance cut;x [mm];y [mm];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);
  TH2D* hitMap = new TH2D("hitMap", "Hit map in the time cut, distance cut;x [mm];y [mm];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);

  // 2d histo to study signal in different parts of the strip
  minX = 0;
  maxX = 80;
  binX = 4; // normally 5
  minY = -50.5;
  maxY = 400.5;
  //  binY = 350;
  binY = 150;
  TH2D* signalStrip = new TH2D("signalStrip", "Signal in various strip parts (2 strips surrounding the hit position), time cut, dist cut;Position in the strip [#mum];Signal [ADC counts]", binX, minX, maxX, binY, minY, maxY);
  TGraphErrors* mpvStrip = new TGraphErrors(); // graph of the landau mpv for slices of the signalStrip
  mpvStrip->SetName("mpvStrip");
  mpvStrip->SetTitle("Landau MPV various strip parts, time cut, distance cut"); 

  TGraphErrors* sigmaStrip = new TGraphErrors(); // graph of the gaus sigma for slices of the signalStrip
  sigmaStrip->SetName("sigmaStrip");
  sigmaStrip->SetTitle("Gauss sigma various strip parts, time cut, distance cut");

  TGraphErrors* mpvStrip_norm = new TGraphErrors(); // graph of the landau mpv for slices of the signalStrip
  mpvStrip_norm->SetName("mpvStrip_norm");
  mpvStrip_norm->SetTitle("Landau MPV various strip parts, time cut, distance cut, normalized to the 5 strips MPV");

  TGraphErrors* maxDistrStrip_norm = new TGraphErrors(); // graph of the maximum of the charge distr for slices of the signalStrip
  maxDistrStrip_norm->SetName("maxDistrStrip_norm");
  maxDistrStrip_norm->SetTitle("Maximum of the charge distr various strip parts, time cut, distance cut, normalized to the one for 5 strips");

  TGraphErrors* maxDistrStrip = new TGraphErrors(); // graph of the maximum of the charge distr for slices of the signalStrip
  maxDistrStrip->SetName("maxDistrStrip");
  //maxDistrStrip->SetTitle("Maximum of the charge distr various strip parts, time cut, distance cut");
  maxDistrStrip->SetTitle("Sum signal L+R");

  TGraphAsymmErrors* medianDistrStrip = new TGraphAsymmErrors(); // graph of the median of the charge distr for slices of the signalStrip
  medianDistrStrip->SetName("medianDistrStrip");
  medianDistrStrip->SetTitle("Sum signal L+R");

  // single contributions of the 4 strips surrounding the hit point
  TH2D* signalStripL = new TH2D("signalStripL", "Signal on Left strip, dist cut;Position in the strip [#mum];Signal [ADC counts]", binX, minX, maxX, binY, minY, maxY);
  TH2D* signalStripR = new TH2D("signalStripR", "Signal on Right strip, dist cut;Position in the strip [#mum];Signal [ADC counts]", binX, minX, maxX, binY, minY, maxY);
  TH2D* signalStripLp1 = new TH2D("signalStripLp1", "Signal on Left+1 strip, dist cut;Position in the strip [#mum];Signal [ADC counts]", binX, minX, maxX, binY, minY, maxY);
  TH2D* signalStripRp1 = new TH2D("signalStripRp1", "Signal on Right+1 strip, dist cut;Position in the strip [#mum];Signal [ADC counts]", binX, minX, maxX, binY, minY, maxY);

  TGraphErrors* maxDistrStripL = new TGraphErrors();
  maxDistrStripL->SetName("maxDistrStripL");
  //maxDistrStripL->SetTitle("Maximum of the charge distr Left strip, time cut, distance cut");
  maxDistrStripL->SetTitle("Signal L");

  TGraphErrors* maxDistrStripR = new TGraphErrors();
  maxDistrStripR->SetName("maxDistrStripR");
  //maxDistrStripR->SetTitle("Maximum of the charge distr Right strip, time cut, distance cut");
  maxDistrStripR->SetTitle("Signal R");

  TGraphErrors* maxDistrStripLp1 = new TGraphErrors();
  maxDistrStripLp1->SetName("maxDistrStripLp1");
  //maxDistrStripLp1->SetTitle("Maximum of the charge distr Left+1 strip, time cut, distance cut");
  maxDistrStripLp1->SetTitle("Signal L+1");

  TGraphErrors* maxDistrStripRp1 = new TGraphErrors();
  maxDistrStripRp1->SetName("maxDistrStripRp1");
  //maxDistrStripRp1->SetTitle("Maximum of the charge distr Right+1 strip, time cut, distance cut");
  maxDistrStripRp1->SetTitle("Signal R+1");

  // eta distribution
  minX = -0.5;
  maxX = 1.5;
  binX = 100;
  TH1D* etaDistrTimeCutDistCut = new TH1D("etaDistrTimeCutDistCut", "#eta distribution in the time cut, dist cut;#eta;Entries", binX, minX, maxX);
  TH1D* CDFetaDistrTimeCutDistCut = new TH1D("CDFetaDistrTimeCutDistCut", "CDF #eta distribution in the time cut, dist cut;#eta;#int #eta", binX, minX, maxX);
  TH1D* etaDistrTrackTimeCut = new TH1D("etaDistrTrackTimeCut", "Track based #eta distribution in the time cut;#eta;Entries", binX, minX, maxX);
  TH1D* CDFetaDistrTrackTimeCut = new TH1D("CDFetaDistrTrackTimeCut", "CDF track based #eta distribution in the time cut;#eta;#int #eta", binX, minX, maxX);
  TH1D* etaDistrTrackTimeCutLowPH = new TH1D("etaDistrTrackTimeCutLowPH", "Track based #eta distribution in the time cut, PH < 15 ADC;#eta;Entries", binX, minX, maxX);

  // scatter plot eta
  TH2D* etaClustVsPos = new TH2D("etaClustVsPos", "Cluster #eta vs reduced track position;Position in the strip [AU];#eta cluster", 50, 0, 1, binX, minX, maxX);
  TH2D* etaTrackVsPos = new TH2D("etaTrackVsPos", "Track #eta vs reduced track position;Position in the strip [AU];#eta track", 50, 0, 1, binX, minX, maxX);

  TH1D* noiseHistCh[nChannels]; // calculation of the noise
  for(int i = 0; i < nChannels; ++i)
    {
      sprintf(name, "noiseDistrChannel_%i", i);
      sprintf(title, "Signal distribution (positivized) not associated with a hit, channel %i;Signal [ADC];Entries", i);
      noiseHistCh[i] = new TH1D(name, title, 111, -40.5, 40.5);
    }
  TH1D* fittedNoiseDistr = new TH1D("fittedNoiseDistr", "Distribution of the fitted noise;Noise [ADC];Entries", 61, -0.5, 20.5);

  TH1D* noiseHistCh_electrons[nChannels]; // calculation of the noise
  for(int i = 0; i < nChannels; ++i)
    {
      sprintf(name, "noiseDistrCh_electrons_%i", i);
      sprintf(title, "Signal distribution (positivized) not associated with a hit, channel %i;Signal [e^{-}];Entries", i);
      noiseHistCh_electrons[i] = new TH1D(name, title, 111, -40.5 * ADCtoe, 40.5 * ADCtoe);
    }
  TH1D* fittedNoiseDistr_electrons = new TH1D("fittedNoiseDistr_electrons", "Distribution of the fitted noise;Noise [e^{-}];Entries", 61, -0.5 * ADCtoe, 20.5 * ADCtoe);

  int nTrks = 0; // number of tracks in one event
  long int evtMrk = -1; // event marker

  std::vector<track> trkVec; // tracks in one event
  double evtAliPH[nChannels] = {0}; // alibava pulse height in on e event
  float evtAliTime = -1;
  float evtAliTemp = -275;
  track* trk;

  int extraCh = -1; // extrapolated channel number
  double hitCharge = 0; // charge on the hit
  int hiChargeCh = -1; // extrapolated ch num of the hit with highest charge in the event
  double highestCharge = 0; // highest hit charge in the event (believed to be the particle that passes the detector in time)
  int trackPos = 0; // position of the track in the track vector
  int iBin = 0; // bin to be filled
  double posX = 0; // store positions
  double posY = 0; // this one will be the mod of the position
  double oldContent = 0;

  double* intPart = new double(0); // used in the signalStrip histo

  bool goodNoise[nChannels]; // determine wether a channel was hit or not, for the noise analysis
  double noiseSum; // used for the noise over multiple channels
  int summed; // number of summed ch for the noise

  std::vector<double> hitSignal; // vectors to store the hit and strip hph signal to be used for a cut using strip hph
  std::vector<double> stripHPHsignal;

  int highestPHstrip = -1; // strip with the highest ph in the hit
  double phHighestStrip = 0; // charge on the strip with the highest charge

  double noiseHPH = 0; // to calculate the background for the signal distribution of the strip with highest ph

  double phR = -1; // variables for the eta distribution
  double phL = -1;
  double phLp1 = -1;
  double phRp1 = -1;

  bool analyzeEvent = false;

  long int nEntries = trkTree->GetEntries();

  if(maxEntryNum == 0 || maxEntryNum > nEntries) // if there is a maximum number of entries
    maxEntryNum = nEntries;

  for(int i = 0; i < maxEntryNum; ++i)
    {
      trkTree->GetEntry(i);

      // put the track info into a structure
      trk = new track();
      trk->extraPosDUT_global[0] = dutTrackX_global; // global reference frame
      trk->extraPosDUT_global[1] = dutTrackY_global;
      trk->extraPosDUT_local[0] = dutTrackX_local; // local reference frame
      trk->extraPosDUT_local[1] = dutTrackY_local;
      trk->extraPosDUT_pixel[0] = dutTrackX_pixel; // dut ref frame, in pixel / strip number
      trk->extraPosDUT_pixel[1] = dutTrackY_pixel;
      trk->measPosDUT_global[0] = dutHitX_global;
      trk->measPosDUT_global[1] = dutHitY_global;
      trk->measPosDUT_local[0] = dutHitX_local;
      trk->measPosDUT_local[1] = dutHitY_local;
      trk->measPosDUT_pixel[0] = dutHitX_pixel;
      trk->measPosDUT_pixel[1] = dutHitY_pixel;
      trk->entryNum = i; 

      hitMapDUTtele->Fill(dutTrackX_global, dutTrackY_global);
      extraChDistr->Fill(dutTrackY_pixel);

      timeTrack->Fill(alibava_TDC);
      timeTrackPosY->Fill(alibava_TDC, dutTrackY_global);

      if(evtMrk == EvtNr)
	{ // add track info to some container
	  nTrks++;

	  trkVec.push_back(*trk);
	}
      else // new event
	{ // analyze the old event, if there are tracks
	  if(trkVec.size() != 0)
	    {
	      trkEvt->Fill(nTrks);
	      trkVsEvt->SetPoint(trkVsEvt->GetN(), evtMrk, nTrks);

	      for(int iCh = 0; iCh < nChannels; ++iCh) goodNoise[iCh] = true; // reset the flag for the noise analysis

	      for(unsigned int iTrk = 0; iTrk < trkVec.size(); ++iTrk)
		{
		  if(trkVec.at(iTrk).measPosDUT_global[0] > -900 && trkVec.at(iTrk).measPosDUT_global[1] > -900) // if there is a matched hit
		    {
		      residualsYvsEntry->Fill(trkVec.at(iTrk).entryNum, trkVec.at(iTrk).measPosDUT_global[1] - trkVec.at(iTrk).extraPosDUT_global[1]);
		      hitMapMatched->Fill(trkVec.at(iTrk).measPosDUT_global[0], trkVec.at(iTrk).measPosDUT_global[1]);
		      residualsX->Fill(trkVec.at(iTrk).measPosDUT_global[0] - trkVec.at(iTrk).extraPosDUT_global[0]);
		      residualsY->Fill(trkVec.at(iTrk).measPosDUT_global[1] - trkVec.at(iTrk).extraPosDUT_global[1]);
		      residuals->Fill(trkVec.at(iTrk).measPosDUT_global[0] - trkVec.at(iTrk).extraPosDUT_global[0], trkVec.at(iTrk).measPosDUT_global[1] - trkVec.at(iTrk).extraPosDUT_global[1]);
		      residualsXvsX->Fill(trkVec.at(iTrk).extraPosDUT_global[0], trkVec.at(iTrk).measPosDUT_global[0] - trkVec.at(iTrk).extraPosDUT_global[0]);
		      residualsYvsY->Fill(trkVec.at(iTrk).extraPosDUT_global[1], trkVec.at(iTrk).measPosDUT_global[1] - trkVec.at(iTrk).extraPosDUT_global[1]);
		      residualsXvsY->Fill(trkVec.at(iTrk).extraPosDUT_global[1], trkVec.at(iTrk).measPosDUT_global[0] - trkVec.at(iTrk).extraPosDUT_global[0]);
		      residualsYvsX->Fill(trkVec.at(iTrk).extraPosDUT_global[0], trkVec.at(iTrk).measPosDUT_global[1] - trkVec.at(iTrk).extraPosDUT_global[1]);
		      residualsYvsEvt->Fill(evtMrk, trkVec.at(iTrk).measPosDUT_global[1] - trkVec.at(iTrk).extraPosDUT_global[1]);
		    }

		  // the center of the channel is at 0
		  if(modf(trkVec.at(iTrk).extraPosDUT_pixel[1], intPart) < 0.5)
		    extraCh = trkVec.at(iTrk).extraPosDUT_pixel[1];
		  else
		    extraCh = trkVec.at(iTrk).extraPosDUT_pixel[1] + 1;

		  for(int iCh = extraCh - maxDist; iCh <= extraCh + maxDist; ++iCh) // exclude channels that may have charge deposit from the noise analysis
		    if(iCh >= 0 && iCh < nChannels) // protect limits
		      goodNoise[iCh] = false;
		} // loop on the tracks

	      for(int iCh = 0; iCh < nChannels; ++iCh) // fill noise histos
	      	if(goodNoise[iCh] && evtAliPH[iCh] != 0) // no ph == 0 and no ch belonging toany hit
	      	  {
	      	    noiseDistr->Fill(evtAliPH[iCh] * polarity);
	      	    noiseHistCh[iCh]->Fill(evtAliPH[iCh] * polarity / scaleFactor);
	      	    noiseHistCh_electrons[iCh]->Fill(evtAliPH[iCh] * polarity / scaleFactor * ADCtoe);
	      	  }

	      noiseSum = 0;
	      noiseHPH = -1e3;
	      for(int iCh = 0; iCh < nChannels; ++iCh) // fill the histo of noise for grouped channels
	      	{
	      	  if(evtAliPH[iCh] == 0 || goodNoise[iCh] == false) // bad channels or channels belonging to a hit
		    {
		      summed = 0;
		      noiseSum = 0;
		      noiseHPH = -1e3;
		    }
		  else
		    {
		      noiseSum += evtAliPH[iCh];
		      summed++;

		      if(evtAliPH[iCh] * polarity > noiseHPH)
			noiseHPH = evtAliPH[iCh] * polarity;
		    }

		  if(summed == maxDist * 2 + 1)
		    {
		      noiseDistrGroup->Fill(noiseSum * polarity);
		      noiseDistrGroup_electrons->Fill(noiseSum * polarity * ADCtoe);
		      backgroundDistrHPH->Fill(noiseHPH);
		      backgroundDistrHPH_electrons->Fill(noiseHPH * ADCtoe);
		      summed = 0;
		      noiseSum = 0;
		      noiseHPH = -1e3;
		    }
	      	} // noise over group of channels

	      noiseSum = 0;
	      for(int iCh = 0; iCh < nChannels; ++iCh) // fill the histo of noise for pairs of channels
	      	{
	      	  if(evtAliPH[iCh] == 0 || goodNoise[iCh] == false) // bad channels or channels belonging to a hit
		    {
		      summed = 0;
		      noiseSum = 0;
		    }
		  else
		    {
		      noiseSum += evtAliPH[iCh];
		      summed++;
		    }

		  if(summed == 2)
		    {
		      noiseDistrPair->Fill(noiseSum * polarity);
		      summed = 0;
		      noiseSum = 0;
		    }
	      	} // noise for pairs of channels
	    }

	  analyzeEvent = true; // variable that determines wether an event will be analyzed for the charge
	  for(unsigned int iTrk = 0; iTrk < trkVec.size(); ++iTrk) // check that all the tracks are in the geom cut (in Y and X)
	    {
	      // the center of the channel is at 0
	      if(modf(trkVec.at(iTrk).extraPosDUT_pixel[1], intPart) < 0.5)
		extraCh = trkVec.at(iTrk).extraPosDUT_pixel[1];
	      else
		extraCh = trkVec.at(iTrk).extraPosDUT_pixel[1] + 1;

	      if(trkVec.at(iTrk).extraPosDUT_global[0] <= xCut1 || trkVec.at(iTrk).extraPosDUT_global[0] >= xCut2) // geom cut in X
		analyzeEvent = false;

	      for(int iCh = extraCh - maxDist; iCh <= extraCh + maxDist; ++iCh) // the track must be on the strips considered for the analysis
		if(iCh >= 0 && iCh < nChannels)
		  {
		    if(evtAliPH[iCh] == 0) analyzeEvent = false;
		  }
		else
		  analyzeEvent = false;

	      if(analyzeEvent == false) break; // if one track does not fullfill the cuts
	    } // geometry cuts to be fullfilled by all the tracks

	  if(analyzeEvent && trkVec.size() != 0) analyzeEvent = true; // check that there is at least a track 
	  else analyzeEvent = false;

	  if(analyzeEvent) // all the event tracks in a sensitive part of the sensor and at least one track
	    {
	      trkEvtSelected->Fill(nTrks);

	      for(unsigned int iTrk = 0; iTrk < trkVec.size(); ++iTrk) // residuals of selected events
		if(trkVec.at(iTrk).measPosDUT_global[0] > -900 && trkVec.at(iTrk).measPosDUT_global[1] > -900) // if there is a matched hit
		  residualsYselected->Fill(trkVec.at(iTrk).measPosDUT_global[1] - trkVec.at(iTrk).extraPosDUT_global[1]);

	      for(unsigned int iTrk = 0; iTrk < trkVec.size(); ++iTrk) // loop to select the right track
	      	{
		  // the center of the channel is at 0
		  if(modf(trkVec.at(iTrk).extraPosDUT_pixel[1], intPart) < 0.5)
		    extraCh = trkVec.at(iTrk).extraPosDUT_pixel[1];
		  else
		    extraCh = trkVec.at(iTrk).extraPosDUT_pixel[1] + 1;

		  hitMapDUTgoodCh->Fill(trkVec.at(iTrk).extraPosDUT_global[0], trkVec.at(iTrk).extraPosDUT_global[1]);
		  extraChDistrGoodCh->Fill(extraCh);

		  hitCharge = 0;

		  if(highestNeighbor) // if the variable is not 0, just the neighbor with the highest charge is added to the hit
		    {
		      if(fabs(evtAliPH[extraCh + 1]) > fabs(evtAliPH[extraCh - 1]))
			hitCharge = evtAliPH[extraCh + 1] + evtAliPH[extraCh];
		      else
			hitCharge = evtAliPH[extraCh - 1] + evtAliPH[extraCh];
		    }
		  else
		    {
		      for(int iCh = extraCh - maxDist; iCh <= extraCh + maxDist; ++iCh) // sum the charge of the strips around the extrapolated value
			if(iCh >=0 && iCh < nChannels) // protect array margins
			  hitCharge += evtAliPH[iCh];
		    }

		  if(fabs(hitCharge) > fabs(highestCharge))  // select the highest hit charge
		    {
		      highestCharge = hitCharge;
		      hiChargeCh = extraCh;
		      trackPos = iTrk;
		    }
		}// loop to select the right track

	      // fill histos old event
	      signalTime->Fill(evtAliTime, highestCharge);
	      positivizedSignalTime->Fill(evtAliTime, highestCharge * polarity);
	      signalDistr->Fill(highestCharge * polarity);

	      // strip with highest ph in the hit
	      phHighestStrip = -1e3;
	      highestPHstrip = -1;
	      for(int iCh = hiChargeCh - maxDist; iCh <= hiChargeCh + maxDist; ++iCh)// find the strip with the highest ph in the hit
		if(iCh >=0 && iCh < nChannels) // protect array margins
		  if(evtAliPH[iCh] * polarity > phHighestStrip)
		    {
		      phHighestStrip = evtAliPH[iCh] * polarity;
		      highestPHstrip = iCh;
		    }
	      stripHPHSignalTime->Fill(evtAliTime, phHighestStrip);
	      leftStripHPHSignalTime->Fill(evtAliTime, evtAliPH[highestPHstrip - 1] * polarity);
	      rightStripHPHSignalTime->Fill(evtAliTime, evtAliPH[highestPHstrip + 1] * polarity);

	      if(evtAliTime >= timeCut1 && evtAliTime <= timeCut2) // apply time cut
		{
		  signalDistrTimeCut->Fill(highestCharge * polarity);

		  if(abs(hiChargeCh - highestPHstrip) <= 1) // the strip with the highest PH is neighboring the hit one
		  //if(hiChargeCh == highestPHstrip) // the strip with the highest PH is the hit one
		    {
		      signalDistrTimeCutDistCut->Fill(highestCharge * polarity);
		      signalDistrTimeCutDistCut_electrons->Fill(highestCharge * polarity * ADCtoe);
		      stripHPHDistrTimeCutDistCut->Fill(phHighestStrip);
		      stripHPHDistrTimeCutDistCut_electrons->Fill(phHighestStrip * ADCtoe);
		      correlationPHstripHPHhit->Fill(highestCharge * polarity, phHighestStrip);

		      hitSignal.push_back(highestCharge * polarity);
		      stripHPHsignal.push_back(phHighestStrip);

		      // eta distribution
		      if(evtAliPH[highestPHstrip + 1] * polarity > evtAliPH[highestPHstrip - 1] * polarity)
			{
			  phL = phHighestStrip;
			  phR = evtAliPH[highestPHstrip + 1] * polarity;
			}
		      else
			{
			  phL = evtAliPH[highestPHstrip - 1] * polarity;
			  phR = phHighestStrip;
			}
		      etaDistrTimeCutDistCut->Fill(phR / (phR + phL));

		      etaClustVsPos->Fill(modf(trkVec.at(trackPos).extraPosDUT_pixel[1], intPart), phR / (phR + phL));

		      stripHPH_plusNeigh_DistrTimeCutDistCut->Fill(phR + phL);

		      // hitmap and charge map mod 160 (on 2 strips)
		      posX = 1000 * trkVec.at(trackPos).extraPosDUT_local[0]; // assign the positions in x and y
		      posY = abs((int)(1000 * trkVec.at(trackPos).extraPosDUT_local[1]) % (int)(2 * pitch * 1000));
		      //posY = abs((int)(1000 * (trkVec.at(trackPos).extraPosDUT_pixel[1] * pitch + 0.5 * pitch)) % (int)(2 * pitch * 1000));
		      iBin = chargeMapMod160->FindBin(posX, posY); // find the bin
		      oldContent = chargeMapMod160->GetBinContent(iBin);
		      chargeMapMod160->SetBinContent(iBin, oldContent + highestCharge * polarity);
		      
		      hitMapMod160->Fill(posX, posY);
		      //std::cout << hitMapMod160->GetEntries() << std::endl;
		      //hitmap and charge map over the sensor
		      posX = trkVec.at(trackPos).extraPosDUT_local[0];
		      posY = trkVec.at(trackPos).extraPosDUT_local[1];
		      iBin = chargeMap->FindBin(posX, posY); // find the bin
		      oldContent = chargeMap->GetBinContent(iBin);
		      chargeMap->SetBinContent(iBin, oldContent + highestCharge * polarity);
		      
		      hitMap->Fill(posX, posY);
		    } // distance cut

		  // totally track based eta distr
		  //the center of the channel is at 0 (hiChargeCh alone does not have enough information for this)
		  modf(trkVec.at(trackPos).extraPosDUT_pixel[1], intPart);
		  phL = evtAliPH[(int) *intPart];
		  phR = evtAliPH[(int) *intPart + 1];
		  phLp1 = evtAliPH[(int) *intPart - 1];
		  phRp1 = evtAliPH[(int) *intPart + 2];

		  etaDistrTrackTimeCut->Fill(phR / (phR + phL));

		  if(phR / (phR + phL) > 1 || phR / (phR + phL) < 0) // events outside expected eta
		    if(abs(hiChargeCh - highestPHstrip) <= 1) // dist cut
		      signalDistrTimeCutDistCutEta01->Fill(highestCharge * polarity);

		  etaTrackVsPos->Fill(modf(trkVec.at(trackPos).extraPosDUT_pixel[1], intPart), phR / (phR + phL));

		  if(abs(hiChargeCh - highestPHstrip) <= 1) // the strip with the highest PH is neighboring the hit one
		    {
		      posX = modf(trkVec.at(trackPos).extraPosDUT_pixel[1], intPart) * pitch * 1e3; // in um
		      signalStrip->Fill(posX, (phL + phR) * polarity); //distance cut for this guy
		      signalStripL->Fill(posX, phL * polarity);
		      signalStripR->Fill(posX, phR * polarity);
		      signalStripLp1->Fill(posX, phLp1 * polarity);
		      signalStripRp1->Fill(posX, phRp1 * polarity);
		    }

		  // output for test
		  // modf(trkVec.at(trackPos).extraPosDUT_pixel[1], intPart);
		  // if(hiChargeCh != *intPart) std::cout << "difference!!  " << trkVec.at(trackPos).extraPosDUT_pixel[1] << "   " << hiChargeCh << "    " << *intPart << "   " << hiChargeCh - *intPart << std::endl;

		  for(int iCh = 0; iCh < nChannels; ++iCh)
		    if(evtAliPH[iCh] != 0 && !(iCh >= hiChargeCh - maxDist && iCh <= hiChargeCh - maxDist)) // no ph == 0 and no ch belonging to the hit
		      noiseDistrTimeCut->Fill(evtAliPH[iCh] * polarity);

		  // strip with highest ph in the hit
		  stripHPHDistrTimeCut->Fill(phHighestStrip);
		  stripHPHDiffExtra->Fill(hiChargeCh - highestPHstrip);


		  if(highestCharge * polarity < 15) // investigation of the events with low ph
		    {
		      hitMapLowPH->Fill(trkVec.at(trackPos).extraPosDUT_global[0], trkVec.at(trackPos).extraPosDUT_global[1]);
		      trkEvtLowPH->Fill(trkVec.size());
		      stripHPHDiffExtraLowPH->Fill(hiChargeCh - highestPHstrip);
		      etaDistrTrackTimeCutLowPH->Fill(phR / (phR + phL));
		    }
		  diffExtraStripHPHvsPH->Fill(highestCharge * polarity, hiChargeCh - highestPHstrip);

		  for(int iCh = hiChargeCh - maxDist; iCh <= hiChargeCh + maxDist; ++iCh)// ph around strip hph and extrapolated strip
		    if(iCh >=0 && iCh < nChannels) // protect array margins
		      {
			phAroundHPHstripTimeCut->Fill(iCh - highestPHstrip, evtAliPH[iCh] * polarity);
			phAroundExtraStripTimeCut->Fill(iCh - hiChargeCh, evtAliPH[iCh] * polarity);
		      }		 
		} // time cut
	      tempEvt->SetPoint(tempEvt->GetN(), evtMrk, evtAliTemp);
	      tempDistr->Fill(evtAliTemp);
	      scaleFactorDistr->Fill(scaleFactor);
	    } // the analysis of the event should be contained in this scope

	  // set the counters
	  evtMrk = EvtNr;
	  trkVec.clear();
	  highestCharge = 0;
	  // store the first track of the event
	  trkVec.push_back(*trk);
	  nTrks = 1;
	  // store event conditions
	  evtAliTime = alibava_TDC;
	  evtAliTemp = alibava_temp;

	  if(applyTcorr) // determination of the temperature correction
	    scaleFactor = targetGain / (tCorr_p0 + tCorr_p1 * evtAliTemp + tCorr_p2 * evtAliTemp * evtAliTemp);
	  else 
	    scaleFactor = 1;

	  // store the PH from alibava
	  for(int iCh = 0; iCh < nChannels; ++iCh) evtAliPH[iCh] = alibavaPH[iCh] * scaleFactor; // a scale factor is applied on the need

	  for(int iCh = 0; iCh < nChannels; ++iCh) adcTime->Fill(evtAliTime, evtAliPH[iCh]); // fill the histo
	}

      delete trk;
    } // loop on the tracks (entries in the tree)

  double threshold = backgroundDistrHPH->GetMean() + backgroundDistrHPH->GetRMS();
  for(unsigned int i = 0; i < hitSignal.size(); ++i) // implementation of the cut using the strip with highest ph
    if(stripHPHsignal.at(i) > threshold)
      signalDistrTimeDistHPHcut->Fill(hitSignal.at(i));

  TProfile* positivizedSignalTimeProfile = positivizedSignalTime->ProfileX("positivizedSignalTimeProfile");
  positivizedSignalTimeProfile->SetTitle("Positivized signal time profile");

  TProfile* stripHPHSignalTimeProfile = stripHPHSignalTime->ProfileX("stripHPHSignalTimeProfile");
  stripHPHSignalTimeProfile->SetTitle("Time profile of the signal of the strip with highest PH;Time [ns];Signal [ADC]");

  TProfile* leftStripHPHSignalTimeProfile = leftStripHPHSignalTime->ProfileX("leftStripHPHSignalTimeProfile");
  leftStripHPHSignalTimeProfile->SetTitle("Time profile of the signal of the left strip neighboring the one with highest PH;Time [ns];Signal [ADC]");

  TProfile* rightStripHPHSignalTimeProfile = rightStripHPHSignalTime->ProfileX("rightStripHPHSignalTimeProfile");
  rightStripHPHSignalTimeProfile->SetTitle("Time profile of the signal of the right strip neighboring the one with highest PH;Time [ns];Signal [ADC]");

  TProfile* profileResYvsY = residualsYvsY->ProfileX("profileResYvsY");
  profileResYvsY->SetTitle("Profile histo of the residuals along y vs y");
  TProfile* profileResXvsX = residualsXvsX->ProfileX("profileResXvsX");
  profileResXvsX->SetTitle("Profile histo of the residuals along x vs x");
  TProfile* profileResYvsX = residualsYvsX->ProfileX("profileResYvsX");
  profileResYvsX->SetTitle("Profile histo of the residuals along y vs x");
  TProfile* profileResXvsY = residualsXvsY->ProfileX("profileResXvsY");
  profileResXvsY->SetTitle("Profile histo of the residuals along x vs y");

  TDirectory* timeSlicesDir = outFile->mkdir("timeSlices");
  timeSlicesDir->cd();   

  const int nPars = 4; // parameters of the langaus fit
  TGraphErrors* lanGausParVsTime[nPars];
  const char* names[nPars] = {"lanWidthTime", "mpvTime", "areaTime", "gausSigTime"};
  const char* titles[nPars] = {"Fitted Landau width vs time", "Fitted Landau MPV vs time", "Fitted area vs time", "Fitted Gaussian sigma vs time"};
  for(int iPar = 0; iPar < nPars; ++iPar)
    {
      lanGausParVsTime[iPar] = new TGraphErrors();
      lanGausParVsTime[iPar]->SetName(names[iPar]);
      lanGausParVsTime[iPar]->SetTitle(titles[iPar]);
      lanGausParVsTime[iPar]->SetMarkerStyle(8);
    }

  TGraph* maxLanGauFit = new TGraph();
  maxLanGauFit->SetName("maxLanGauFit");
  maxLanGauFit->SetTitle("Maximum of the fitted Landau Gaussian convolution");
  maxLanGauFit->SetMarkerStyle(8);

  TGraph* chi2SliceFit = new TGraph();
  chi2SliceFit->SetName("chi2SliceFit");
  chi2SliceFit->SetTitle("Reduced #chi^{2} of the landau gaussian fits");
  chi2SliceFit->SetMarkerStyle(8);

  int nBins = positivizedSignalTime->GetXaxis()->GetNbins();
  double time = 0;
  double binW = 0;
  TH1D* slice = NULL;
  TF1* fit = NULL;

  TCanvas* fitCan = new TCanvas("fitCan");

  fitResiduals(residualsX, -0.5, 0.5);
  fitResiduals(residualsY, -0.25, 0.25);
  fitResiduals(residualsYselected, -0.15, 0.15);

  // fit noise distr of the noise in groups of channels
  fit = new TF1("gausFit", "gaus", noiseDistrGroup->GetMean() - 2 * noiseDistrGroup->GetRMS(), noiseDistrGroup->GetMean() + 1 * noiseDistrGroup->GetRMS());
  noiseDistrGroup->Fit(fit, "RQ");

  // TF1* lanGausFitFunc = gausLanGausFitFixGaus(signalDistrTimeCutDistCut, negSigmaFit, posSigmaFit,
  // 					      fit->GetParameter(1), fit->GetParameter(2)); // gaus mean and all the sigma determined from the noise distr and landau gauss convolution fitted simultaneously

  TF1* lanGausFitFunc = gausLanGausFitFixGausNoise(signalDistrTimeCutDistCut, negSigmaFit, posSigmaFit,
   						   fit->GetParameter(1), fit->GetParameter(2)); // gaus mean and sigma determined from the noise distr and landau gauss convolution fitted simultaneously

  // fit noise distr of the noise in groups of channels
  fit = new TF1("gausFit", "gaus", noiseDistrGroup_electrons->GetMean() - 2 * noiseDistrGroup_electrons->GetRMS(), noiseDistrGroup_electrons->GetMean() + 1 * noiseDistrGroup_electrons->GetRMS());
  noiseDistrGroup_electrons->Fit(fit, "RQ");

  TF1* lanGausFitFunc_electrons = gausLanGausFitFixGausNoise(signalDistrTimeCutDistCut_electrons, negSigmaFit, posSigmaFit,
							     fit->GetParameter(1), fit->GetParameter(2)); // gaus mean and sigma determined from the noise distr and landau gauss convolution fitted simultaneously

  //lanGausFit(signalDistrTimeCut, negSigmaFit, posSigmaFit);
  //lanGausFit(signalDistrTimeCutDistCut, negSigmaFit, posSigmaFit);
  // lanGausFit(stripHPHDistrTimeCutDistCut, negSigmaFit, posSigmaFit);
  lanGausFit(signalDistrTimeDistHPHcut, 10, 10); // fixed range !!!!!!
  lanGausFit(stripHPH_plusNeigh_DistrTimeCutDistCut, negSigmaFit, posSigmaFit);

  //gausLanGausFit(signalDistrTimeCutDistCut, negSigmaFit, posSigmaFit); // peack at 0 and landau gauss convolution fitted simultaneously
  // gausLanGausFit(signalDistrTimeCut, negSigmaFit, posSigmaFit);
  // gausLanGausFit(stripHPHDistrTimeCutDistCut, negSigmaFit, posSigmaFit); 
  // gausLanGausFit(stripHPH_plusNeigh_DistrTimeCutDistCut, negSigmaFit, posSigmaFit);

  for(int iBin = 1; iBin <= nBins; ++iBin) // fit the signal distributions of the various times
    {
      time = positivizedSignalTime->GetXaxis()->GetBinCenter(iBin);
      binW = positivizedSignalTime->GetXaxis()->GetBinWidth(iBin);

      sprintf(name, "timeSlice_%f", time);
      sprintf(title, "Hit charge distribution (positivized) at time %f;Hit charge [ADC];Entries", time);

      slice = positivizedSignalTime->ProjectionY(name, iBin, iBin);
      slice->SetTitle(title);

      fit = lanGausFit(slice, negSigmaFit, posSigmaFit);

      for(int iPar = 0; iPar < nPars; ++iPar)
	{
	  if(fit->GetParameter(iPar) > -200) // exclude nonsense
	    {
	      //std::cout << iPar << "   " <<  fit->GetParameter(iPar) << "  +-  " << fit->GetParError(iPar) << std::endl;
	      lanGausParVsTime[iPar]->SetPoint(lanGausParVsTime[iPar]->GetN(), time, fit->GetParameter(iPar));
	      lanGausParVsTime[iPar]->SetPointError(lanGausParVsTime[iPar]->GetN() - 1, binW / 2, fit->GetParError(iPar));
	    }
	}

      if(fit->GetMaximumX() > -200) // exclude nonsense
	maxLanGauFit->SetPoint(maxLanGauFit->GetN(), time, fit->GetMaximumX());

      if(fit->GetChisquare() > 0 && fit->GetChisquare() < 1e6 && fit->GetNDF() > 0) // avoid to put nonsense in the graph
	chi2SliceFit->SetPoint(chi2SliceFit->GetN(), time, fit->GetChisquare() / fit->GetNDF());

      slice->Write();
    }

  // maximum of the charge distr in adc for the normalized signal maximum over the strips
  double rangeConstNeg = 1.8;
  double rangeConstPos = 2.5;
  //TF1* fitMax = new TF1("fitMax", "[0] + [1] * TMath::Power([2] - x, 2)");
  TF1* fitMax = new TF1("fitMax", "gaus");
  fitMax->SetRange(lanGausFitFunc->GetParameter(4) - rangeConstNeg * lanGausFitFunc->GetParameter(6), lanGausFitFunc->GetParameter(4) + rangeConstPos * lanGausFitFunc->GetParameter(6));
  fitMax->SetParameter(0, lanGausFitFunc->GetMaximum());
  fitMax->SetParameter(2, sqrt(pow(lanGausFitFunc->GetParameter(6), 2) + pow(lanGausFitFunc->GetParameter(3), 2)));
  fitMax->SetParameter(1, lanGausFitFunc->GetParameter(4));  
  signalDistrTimeCutDistCut->Fit(fitMax, "RQN");
  double max5strip = fitMax->GetParameter(1);
  double max5stripErr = fitMax->GetParError(1);

  TDirectory* posSlicesDir = outFile->mkdir("positionSlices");
  posSlicesDir->cd();   

  nBins = signalStrip->GetXaxis()->GetNbins();
  double pos;
  for(int iBin = 1; iBin <= nBins; ++iBin) // fit the signal distributions in various strip positions
    {
      pos = signalStrip->GetXaxis()->GetBinCenter(iBin);
      binW = signalStrip->GetXaxis()->GetBinWidth(iBin);

      sprintf(name, "posSlice_%f", pos);
      sprintf(title, "Hit charge distribution (positivized) at position %f;Hit charge [ADC];Entries", pos);

      slice = signalStrip->ProjectionY(name, iBin, iBin);
      // slice->Sumw2();
      // slice->Add(signalStrip->ProjectionY(name, nBins - iBin + 1, nBins - iBin + 1)); // add symmetric bin to increase statistics
      // slice->Sumw2();
      slice->SetTitle(title);

      // use noise from pair strip, since the charge is now summed from 2 strips
      fit = gausLanGausFitFixGausNoise(slice, negSigmaFit, posSigmaFit,
				       noiseDistrPair->GetMean(), noiseDistrPair->GetRMS()); // gaus mean and sigma determined from the noise distr and landau gauss convolution fitted simultaneously

      double gainMeas;
      double gainMeasErr;
      double tempTotErr = sqrt(pow(tempErr, 2) + pow(tempDistr->GetMeanError(), 2));
      if(tCorr_p1 || tCorr_p0)
	{
	  gainMeas = tCorr_p0 + tCorr_p1 * tempDistr->GetMean();
	  gainMeasErr = sqrt(pow(tCorr_p0Err, 2) + pow(tCorr_p1Err * tempDistr->GetMean(), 2) + pow(tCorr_p1 * tempTotErr, 2) + 2 * tCorr_p0p1Cov * tempDistr->GetMean());
	}

      double error = fit->GetParameter(4) * sqrt(pow(targetGainErr / targetGain, 2) + pow(gainMeasErr / gainMeas, 2) + pow(fit->GetParError(4) / fit->GetParameter(4), 2));

      mpvStrip->SetPoint(mpvStrip->GetN(), pos, fit->GetParameter(4));
      mpvStrip->SetPointError(mpvStrip->GetN() - 1, binW / 2, error);

      error = fit->GetParameter(6) * sqrt(pow(targetGainErr / targetGain, 2) + pow(gainMeasErr / gainMeas, 2) + pow(fit->GetParError(6) / fit->GetParameter(6), 2));
      sigmaStrip->SetPoint(sigmaStrip->GetN(), pos, fit->GetParameter(6));
      sigmaStrip->SetPointError(sigmaStrip->GetN() - 1, binW / 2, error);

      // normalization to the 5 strip MPV
      error = fit->GetParameter(4) / lanGausFitFunc->GetParameter(4) * sqrt(pow(fit->GetParError(4) / fit->GetParameter(4), 2) + pow(lanGausFitFunc->GetParError(4) / lanGausFitFunc->GetParameter(4), 2));

      mpvStrip_norm->SetPoint(mpvStrip_norm->GetN(), pos, fit->GetParameter(4) / lanGausFitFunc->GetParameter(4));
      mpvStrip_norm->SetPointError(mpvStrip_norm->GetN() - 1, binW / 2, error);

      // fitMax->SetRange(fit->GetParameter(4) - rangeConstNeg * fit->GetParameter(6), fit->GetParameter(4) + rangeConstPos * fit->GetParameter(6));
      //double maxPos = slice->GetXaxis()->GetBinCenter(slice->GetMaximumBin());
      // fitMax->SetRange(maxPos - rangeConstNeg * slice->GetRMS(), maxPos + rangeConstPos * slice->GetRMS());

      // fitMax->SetParameter(0, fit->GetMaximum());
      // fitMax->SetParameter(2, sqrt(pow(fit->GetParameter(6), 2) + pow(fit->GetParameter(3), 2)));
      // fitMax->SetParameter(1, fit->GetParameter(4));
      // slice->Fit(fitMax, "RQ");
      fitMax = fitPeak(slice, rangeConstNeg, rangeConstPos);
      double maxSlice = fitMax->GetParameter(1);
      double maxSliceErr = fitMax->GetParError(1);

      error = maxSlice / max5strip * sqrt(pow(maxSliceErr / maxSlice, 2) + pow(max5stripErr / max5strip, 2));

      maxDistrStrip_norm->SetPoint(maxDistrStrip_norm->GetN(), pos, maxSlice / max5strip);
      maxDistrStrip_norm->SetPointError(maxDistrStrip_norm->GetN() - 1, binW / 2, error);

      maxDistrStrip->SetPoint(maxDistrStrip->GetN(), pos, maxSlice);
      maxDistrStrip->SetPointError(maxDistrStrip->GetN() - 1, binW / 2, maxSliceErr);
      //maxDistrStrip->SetPointError(maxDistrStrip->GetN() - 1, 0, maxSliceErr);

      findFracPos med(slice, 0.5);
      medianDistrStrip->SetPoint(medianDistrStrip->GetN(), pos, med.fracPos);
      medianDistrStrip->SetPointError(medianDistrStrip->GetN() - 1, binW / 2, binW / 2, med.errLo, med.errHi);

      slice->Write();
    }

  TDirectory* posSlicesDirL = outFile->mkdir("positionSlicesL");
  TDirectory* posSlicesDirR = outFile->mkdir("positionSlicesR");
  TDirectory* posSlicesDirLp1 = outFile->mkdir("positionSlicesLp1");
  TDirectory* posSlicesDirRp1 = outFile->mkdir("positionSlicesRp1");

  for(int iBin = 1; iBin <= nBins; ++iBin) // fit the signal distributions in various strip positions
    {
      pos = signalStrip->GetXaxis()->GetBinCenter(iBin);
      binW = signalStrip->GetXaxis()->GetBinWidth(iBin);

      sprintf(name, "posSlice_%f", pos);
      sprintf(title, "Hit charge distribution (positivized) at position %f;Hit charge [ADC];Entries", pos);

      slice = signalStripL->ProjectionY(name, iBin, iBin);
      slice->SetTitle(title);
      fitMax = fitPeak(slice, rangeConstNeg, rangeConstPos);
      maxDistrStripL->SetPoint(maxDistrStripL->GetN(), pos, fitMax->GetParameter(1));
      maxDistrStripL->SetPointError(maxDistrStripL->GetN() - 1, binW / 2, fitMax->GetParError(1));
      posSlicesDirL->cd();   
      slice->Write();

      slice = signalStripR->ProjectionY(name, iBin, iBin);
      slice->SetTitle(title);
      fitMax = fitPeak(slice, rangeConstNeg, rangeConstPos);
      maxDistrStripR->SetPoint(maxDistrStripR->GetN(), pos, fitMax->GetParameter(1));
      maxDistrStripR->SetPointError(maxDistrStripR->GetN() - 1, binW / 2, fitMax->GetParError(1));
      posSlicesDirR->cd();   
      slice->Write();

      slice = signalStripLp1->ProjectionY(name, iBin, iBin);
      slice->SetTitle(title);
      fitMax = fitPeak(slice, rangeConstNeg, rangeConstPos);
      maxDistrStripLp1->SetPoint(maxDistrStripLp1->GetN(), pos, fitMax->GetParameter(1));
      maxDistrStripLp1->SetPointError(maxDistrStripLp1->GetN() - 1, binW / 2, fitMax->GetParError(1));
      posSlicesDirLp1->cd();   
      slice->Write();

      slice = signalStripRp1->ProjectionY(name, iBin, iBin);
      slice->SetTitle(title);
      fitMax = fitPeak(slice, rangeConstNeg, rangeConstPos);
      maxDistrStripRp1->SetPoint(maxDistrStripRp1->GetN(), pos, fitMax->GetParameter(1));
      maxDistrStripRp1->SetPointError(maxDistrStripRp1->GetN() - 1, binW / 2, fitMax->GetParError(1));
      posSlicesDirRp1->cd();   
      slice->Write();
    }

  maxDistrStrip->SetFillColor(kWhite);
  maxDistrStripL->SetLineColor(kRed);
  maxDistrStripL->SetFillColor(kWhite);
  maxDistrStripR->SetLineColor(kGreen);
  maxDistrStripR->SetFillColor(kWhite);
  maxDistrStripLp1->SetLineColor(kBlue);
  maxDistrStripLp1->SetFillColor(kWhite);
  maxDistrStripRp1->SetLineColor(kViolet);
  maxDistrStripRp1->SetFillColor(kWhite);

  TMultiGraph* sigLR = new TMultiGraph();
  sigLR->SetName("signalLeftRightSum");
  sigLR->SetTitle("Max signal distr. Errors from fit (no T corr errors)");
  sigLR->Add(maxDistrStrip);
  sigLR->Add(maxDistrStripL);
  sigLR->Add(maxDistrStripR);
  sigLR->Add(maxDistrStripLp1);
  sigLR->Add(maxDistrStripRp1);

  TH1D* entriesStrip = signalStrip->ProjectionX();
  entriesStrip->SetName("entriesStrip");
  entriesStrip->SetTitle("Entries vs position between strips;Position [A.U.];Entries");

  TGraphErrors* noiseMeanCh = new TGraphErrors();
  noiseMeanCh->SetName("noiseMeanCh");
  noiseMeanCh->SetTitle("Mean of the noise fit vs channel");
  noiseMeanCh->SetMarkerStyle(8);

  TGraphErrors* noiseCh = new TGraphErrors();
  noiseCh->SetName("noiseCh");
  noiseCh->SetTitle("Sigma of the noise fit vs channel");
  noiseCh->SetMarkerStyle(8);

  TGraph* noiseChiCh = new TGraph();
  noiseChiCh->SetName("noiseChiCh");
  noiseChiCh->SetTitle("#chi^{2} / ndf of the noise fit vs channel");
  noiseChiCh->SetMarkerStyle(8);

  TGraphErrors* noiseMeanCh_electrons = new TGraphErrors();
  noiseMeanCh_electrons->SetName("noiseMeanCh_electrons");
  noiseMeanCh_electrons->SetTitle("Mean of the noise fit vs channel");
  noiseMeanCh_electrons->SetMarkerStyle(8);

  TGraphErrors* noiseCh_electrons = new TGraphErrors();
  noiseCh_electrons->SetName("noiseCh_electrons");
  noiseCh_electrons->SetTitle("Sigma of the noise fit vs channel");
  noiseCh_electrons->SetMarkerStyle(8);

  for(int i = 0; i < nChannels; ++i) // fit of the noise distribution for all the channels
    if(noiseHistCh[i]->GetEntries() != 0)
      {
	fit = new TF1("gausFit", "gaus", noiseHistCh[i]->GetMean() - 3 * noiseHistCh[i]->GetRMS(), noiseHistCh[i]->GetMean() + 2 * noiseHistCh[i]->GetRMS());
	noiseHistCh[i]->Fit(fit, "RQ");

	noiseMeanCh->SetPoint(noiseMeanCh->GetN(), i, fit->GetParameter(1));
	noiseMeanCh->SetPointError(noiseMeanCh->GetN() - 1, 0, fit->GetParError(1));

	noiseCh->SetPoint(noiseCh->GetN(), i, fit->GetParameter(2));
	noiseCh->SetPointError(noiseCh->GetN() - 1, 0, fit->GetParError(2));

	noiseChiCh->SetPoint(noiseChiCh->GetN(), i, fit->GetChisquare() / fit->GetNDF());

	fittedNoiseDistr->Fill(fit->GetParameter(2));

	// noise in electrons
	fit = new TF1("gausFit", "gaus", noiseHistCh_electrons[i]->GetMean() - 3 * noiseHistCh_electrons[i]->GetRMS(), noiseHistCh_electrons[i]->GetMean() + 2 * noiseHistCh_electrons[i]->GetRMS());
	noiseHistCh_electrons[i]->Fit(fit, "RQ");

	double error = fit->GetParameter(1) * sqrt(pow(ADCtoeErr / ADCtoe, 2) + pow(fit->GetParError(1) / fit->GetParameter(1), 2));
	noiseMeanCh_electrons->SetPoint(noiseMeanCh_electrons->GetN(), i, fit->GetParameter(1));
	noiseMeanCh_electrons->SetPointError(noiseMeanCh_electrons->GetN() - 1, 0, error);

	error = fit->GetParameter(2) * sqrt(pow(ADCtoeErr / ADCtoe, 2) + pow(fit->GetParError(2) / fit->GetParameter(2), 2));
	noiseCh_electrons->SetPoint(noiseCh_electrons->GetN(), i, fit->GetParameter(2));
	noiseCh_electrons->SetPointError(noiseCh_electrons->GetN() - 1, 0, fit->GetParError(2));

	fittedNoiseDistr_electrons->Fill(fit->GetParameter(2));
      }

  // fit noise distr for a noise on pairs of channels
  fit = new TF1("gausFit", "gaus", noiseDistrPair->GetMean() - 2 * noiseDistrPair->GetRMS(), noiseDistrPair->GetMean() + 1 * noiseDistrPair->GetRMS());
  noiseDistrPair->Fit(fit, "RQ");

  // fit of the background strip hph distribution, for derivation see logbook entry 16.07.2014
  TF1* fit_bg = new TF1("fit_bg", "[2] * [1] * TMath::Gaus(x, 0, [0], 1) * TMath::Power(0.5 * (1 + TMath::Erf(x / ([0] * TMath::Sqrt(2)))), [1] - 1)");
  fit_bg->SetParameter(0, fittedNoiseDistr->GetMean());
  fit_bg->FixParameter(1, 2 * maxDist + 1); // fixed parameter
  fit_bg->SetParameter(2, backgroundDistrHPH->GetMaximum());
  fit_bg->SetRange(backgroundDistrHPH->GetMean() - 2 * backgroundDistrHPH->GetRMS(), backgroundDistrHPH->GetMean() + 1 * backgroundDistrHPH->GetRMS());
  // backgroundDistrHPH->Fit(fit_bg, "RQ");
  // double fitSig = fit_bg->GetParameter(0);
  // fit_bg->FixParameter(0, fitSig);
  stripHPHDistrTimeCutDistCut->Fit(fit_bg, "RQ");

  // noise peak subtraction for the histo of charge
  TF1* noiseGaus = new TF1("noiseGaus", "gaus", signalDistrTimeCutDistCut->GetXaxis()->GetXmin(), signalDistrTimeCutDistCut->GetXaxis()->GetXmax());
  noiseGaus->SetParameter(0, lanGausFitFunc->GetParameter(0));
  noiseGaus->SetParameter(1, lanGausFitFunc->GetParameter(1));
  noiseGaus->SetParameter(2, lanGausFitFunc->GetParameter(2));
  double cont;
  for(int iBin = 1; iBin < signalDistrTimeCutDistCut->GetNbinsX(); iBin++) // under flow and over flow bins not considered
    {
      cont = signalDistrTimeCutDistCut->GetBinContent(iBin) - noiseGaus->Eval(signalDistrTimeCutDistCut->GetBinCenter(iBin));
      if(cont < 0) cont = 0; // avoid to go negative
      signalDistrTimeCutDistCut_noisePeakSub->SetBinContent(iBin, cont);
    }

  //lanGausFit(signalDistrTimeCutDistCut_noisePeakSub, negSigmaFit * 4, posSigmaFit * 5);

  noiseGaus->SetParameter(0, lanGausFitFunc_electrons->GetParameter(0));
  noiseGaus->SetParameter(1, lanGausFitFunc_electrons->GetParameter(1));
  noiseGaus->SetParameter(2, lanGausFitFunc_electrons->GetParameter(2));

  for(int iBin = 1; iBin < signalDistrTimeCutDistCut_electrons->GetNbinsX(); iBin++) // under flow and over flow bins not considered
    {
      cont = signalDistrTimeCutDistCut_electrons->GetBinContent(iBin) - noiseGaus->Eval(signalDistrTimeCutDistCut_electrons->GetBinCenter(iBin));
      if(cont < 0) cont = 0; // avoid to go negative
      signalDistrTimeCutDistCut_noisePeakSub_electrons->SetBinContent(iBin, cont);
    }

  delete fitCan;

  double chargeSum;
  int nTracks;
  for(int iBinX = 1; iBinX < chargeMapMod160->GetNbinsX(); ++iBinX) // normalize charge map in the time cut
    for(int iBinY = 1; iBinY < chargeMapMod160->GetNbinsY(); ++iBinY)
      {
	chargeSum = chargeMapMod160->GetBinContent(iBinX, iBinY);
	nTracks = hitMapMod160->GetBinContent(iBinX, iBinY);

	if(nTracks != 0)
	  chargeMapMod160Normalized->SetBinContent(iBinX, iBinY, chargeSum / nTracks);
      }

  TH1D* chargeMapMod160NormProjX = chargeMapMod160Normalized->ProjectionX("chargeMapMod160NormProjX"); // projections in x and y
  chargeMapMod160NormProjX->SetTitle("X projection of the normalized charge distribution");
  TH1D* chargeMapMod160NormProjY = chargeMapMod160Normalized->ProjectionY("chargeMapMod160NormProjY");
  chargeMapMod160NormProjY->SetTitle("Y projection of the normalized charge distribution");

  for(int iBinX = 1; iBinX < chargeMap->GetNbinsX(); ++iBinX) // normalize charge map in the time cut
    for(int iBinY = 1; iBinY < chargeMap->GetNbinsY(); ++iBinY)
      {
	chargeSum = chargeMap->GetBinContent(iBinX, iBinY);
	nTracks = hitMap->GetBinContent(iBinX, iBinY);

	if(nTracks != 0)
	  chargeMapNormalized->SetBinContent(iBinX, iBinY, chargeSum / nTracks);

      }

  TH1D* chargeMapNormProjX = chargeMapNormalized->ProjectionX("chargeMapNormProjX"); // projections in x and y
  chargeMapNormProjX->SetTitle("X projection of the normalized charge distribution");
  TH1D* chargeMapNormProjY = chargeMapNormalized->ProjectionY("chargeMapNormProjY");
  chargeMapNormProjY->SetTitle("Y projection of the normalized charge distribution");

  // background subtraction for the strip with highest PH
  double bgMean = backgroundDistrHPH->GetMean();
  double bgRMS = backgroundDistrHPH->GetRMS();
  int binStart = backgroundDistrHPH->GetXaxis()->FindBin(bgMean - 2 * bgRMS);
  int binStop = backgroundDistrHPH->GetXaxis()->FindBin(bgMean + 0.7 * bgRMS);
  double bgInt = backgroundDistrHPH->Integral(binStart, binStop);
  double sigPlusBgInt = stripHPHDistrTimeCutDistCut->Integral(binStart, binStop);
  TH1D* stripHPHDistrTimeCutDistCut_BGsub = new TH1D(*stripHPHDistrTimeCutDistCut);
  stripHPHDistrTimeCutDistCut_BGsub->SetName("stripHPHDistrTimeCutDistCut_BGsub");
  stripHPHDistrTimeCutDistCut_BGsub->SetTitle("Signal distr strip Hi PH time cut dist cut, bg subtracted;Signal [ADC];Entries");

  backgroundDistrHPH->Sumw2();
  stripHPHDistrTimeCutDistCut_BGsub->Sumw2();
  stripHPHDistrTimeCutDistCut_BGsub->Add(backgroundDistrHPH, -1 * sigPlusBgInt / bgInt);

  for(int iBin = 1; iBin < stripHPHDistrTimeCutDistCut_BGsub->GetNbinsX(); iBin++) // avoid negatives
    if(stripHPHDistrTimeCutDistCut_BGsub->GetBinContent(iBin) < 0)
      stripHPHDistrTimeCutDistCut_BGsub->SetBinContent(iBin, 0);

  // background subtraction for the strip with highest PH in electrons
  bgMean = backgroundDistrHPH_electrons->GetMean();
  bgRMS = backgroundDistrHPH_electrons->GetRMS();
  binStart = backgroundDistrHPH_electrons->GetXaxis()->FindBin(bgMean - 2 * bgRMS);
  binStop = backgroundDistrHPH_electrons->GetXaxis()->FindBin(bgMean + 0.7 * bgRMS);
  bgInt = backgroundDistrHPH_electrons->Integral(binStart, binStop);
  sigPlusBgInt = stripHPHDistrTimeCutDistCut_electrons->Integral(binStart, binStop);
  TH1D* stripHPHDistrTimeCutDistCut_BGsub_electrons = new TH1D(*stripHPHDistrTimeCutDistCut_electrons);
  stripHPHDistrTimeCutDistCut_BGsub_electrons->SetName("stripHPHDistrTimeCutDistCut_BGsub_electrons");
  stripHPHDistrTimeCutDistCut_BGsub_electrons->SetTitle("Signal distr strip Hi PH time cut dist cut, bg subtracted;Signal [e^{-}];Entries");

  backgroundDistrHPH_electrons->Sumw2();
  stripHPHDistrTimeCutDistCut_BGsub_electrons->Sumw2();
  stripHPHDistrTimeCutDistCut_BGsub_electrons->Add(backgroundDistrHPH_electrons, -1 * sigPlusBgInt / bgInt);

  for(int iBin = 1; iBin < stripHPHDistrTimeCutDistCut_BGsub_electrons->GetNbinsX(); iBin++) // avoid negatives
    if(stripHPHDistrTimeCutDistCut_BGsub_electrons->GetBinContent(iBin) < 0)
      stripHPHDistrTimeCutDistCut_BGsub_electrons->SetBinContent(iBin, 0);

  // background subtraction using fitted function
  TH1D* stripHPHDistrTimeCutDistCut_BGsub_fromFit = new TH1D(*stripHPHDistrTimeCutDistCut);
  stripHPHDistrTimeCutDistCut_BGsub_fromFit->SetName("stripHPHDistrTimeCutDistCut_BGsub_fromFit");
  stripHPHDistrTimeCutDistCut_BGsub_fromFit->SetTitle("Signal distr strip Hi PH time cut dist cut, bg subtracted usinf the fitted bg;Signal [ADC];Entries");
  stripHPHDistrTimeCutDistCut_BGsub_fromFit->Sumw2();
  fit_bg->SetRange(-1000, 1000);
  stripHPHDistrTimeCutDistCut_BGsub_fromFit->Add(fit_bg, -1);

  // normalized integral from the stripHPH with cuts and bg subtraction
  TH1D* stripHPH_BGsub_integral = new TH1D(*stripHPHDistrTimeCutDistCut_BGsub);
  stripHPH_BGsub_integral->SetName("stripHPH_BGsub_integral");
  stripHPH_BGsub_integral->SetTitle("Normalized integral of the strip with Hi PH, with cuts and BG subtracted;Signal [ADC];Integral");
  binStart = 1;
  binStop = stripHPHDistrTimeCutDistCut_BGsub->GetXaxis()->GetNbins();
  double totArea = stripHPHDistrTimeCutDistCut_BGsub->Integral(binStart, binStop);
  double intBin = 0;

  for(int i = binStart; i < binStop + 1; ++i)
    {
      intBin = stripHPHDistrTimeCutDistCut_BGsub->Integral(binStart, i);
      stripHPH_BGsub_integral->SetBinContent(i, intBin / totArea);
    }

  // normalized integral from the stripHPH with cuts and bg subtraction in electrons
  TH1D* stripHPH_BGsub_integral_electrons = new TH1D(*stripHPHDistrTimeCutDistCut_BGsub_electrons);
  stripHPH_BGsub_integral_electrons->SetName("stripHPH_BGsub_integral_electrons");
  stripHPH_BGsub_integral_electrons->SetTitle("Normalized integral of the strip with Hi PH, with cuts and BG subtracted;Signal [e^{-}];Integral");
  binStart = 1;
  binStop = stripHPHDistrTimeCutDistCut_BGsub_electrons->GetXaxis()->GetNbins();
  totArea = stripHPHDistrTimeCutDistCut_BGsub_electrons->Integral(binStart, binStop);
  intBin = 0;

  for(int i = binStart; i < binStop + 1; ++i)
    {
      intBin = stripHPHDistrTimeCutDistCut_BGsub_electrons->Integral(binStart, i);
      stripHPH_BGsub_integral_electrons->SetBinContent(i, intBin / totArea);
    }

  // normalized histo of the noise distribution
  TH1D* noiseDistr_integral = new TH1D(*noiseDistr);
  noiseDistr_integral->SetName("noiseDistr_integral");
  noiseDistr_integral->SetTitle("Normalized integral of the single strip noise (all the strips of the sensor);Noise [ADC];Integral");
  binStart = 1;
  binStop = noiseDistr->GetXaxis()->GetNbins();
  totArea = noiseDistr->Integral(binStart, binStop);

  for(int i = binStart; i < binStop + 1; ++i)
    {
      intBin = noiseDistr->Integral(binStart, i);
      noiseDistr_integral->SetBinContent(i, intBin / totArea);
    }

  // nomalized CDF of the eta distros
  binStart = 1;
  binStop = etaDistrTimeCutDistCut->GetXaxis()->GetNbins();
  totArea = etaDistrTimeCutDistCut->Integral(binStart, binStop);

  for(int i = binStart; i < binStop + 1; ++i)
    {
      intBin = etaDistrTimeCutDistCut->Integral(binStart, i);
      CDFetaDistrTimeCutDistCut->SetBinContent(i, intBin / totArea);
    }

  binStart = 1;
  binStop = etaDistrTrackTimeCut->GetXaxis()->GetNbins();
  totArea = etaDistrTrackTimeCut->Integral(binStart, binStop);

  for(int i = binStart; i < binStop + 1; ++i)
    {
      intBin = etaDistrTrackTimeCut->Integral(binStart, i);
      CDFetaDistrTrackTimeCut->SetBinContent(i, intBin / totArea);
    }

  // profiles of the eta distr
 TProfile* profileEtaTrack = etaTrackVsPos->ProfileX("profileEtaTrack");
 profileEtaTrack->SetTitle("Profile of #eta track vs poisition");
 TProfile* profileEtaClust = etaClustVsPos->ProfileX("profileEtaClust");
 profileEtaClust->SetTitle("Profile of #eta cluster vs poisition");

  // draw graphs to name the axis
  TCanvas* servCan = new TCanvas("servCan");
  servCan->cd();

  trkVsEvt->Draw("AP");
  trkVsEvt->GetXaxis()->SetTitle("Event number");
  trkVsEvt->GetYaxis()->SetTitle("Number of tracks");

  tempEvt->Draw("AP");
  tempEvt->GetXaxis()->SetTitle("Event number");
  tempEvt->GetYaxis()->SetTitle("Temperature [C]");

  const char* axTitles[nPars] = {"Width [ADC]", "MPV [ADC]", "Area [Entries]", "Gaus sigma [ADC]"};
  for(int iPar = 0; iPar < nPars; ++iPar)
    {
      lanGausParVsTime[iPar]->Draw("AP");
      lanGausParVsTime[iPar]->GetXaxis()->SetTitle("Time [ns]");
      lanGausParVsTime[iPar]->GetYaxis()->SetTitle(axTitles[iPar]);
    }

  maxLanGauFit->Draw("AP");
  maxLanGauFit->GetXaxis()->SetTitle("Time [ns]");
  maxLanGauFit->GetYaxis()->SetTitle("Function max [ADC]");

  chi2SliceFit->Draw("AP");
  chi2SliceFit->GetXaxis()->SetTitle("Time [ns]");
  chi2SliceFit->GetYaxis()->SetTitle("#chi^{2} / ndf");

  noiseMeanCh->Draw("AP");
  noiseMeanCh->GetXaxis()->SetTitle("Channel number");
  noiseMeanCh->GetYaxis()->SetTitle("Mean [ADC]");

  noiseCh->Draw("AP");
  noiseCh->GetXaxis()->SetTitle("Channel number");
  noiseCh->GetYaxis()->SetTitle("Sigma [ADC]");

  noiseChiCh->Draw("AP");
  noiseChiCh->GetXaxis()->SetTitle("Channel number");
  noiseChiCh->GetYaxis()->SetTitle("#chi^{2} / ndf [ADC]");

  noiseMeanCh_electrons->Draw("AP");
  noiseMeanCh_electrons->GetXaxis()->SetTitle("Channel number");
  noiseMeanCh_electrons->GetYaxis()->SetTitle("Mean [e^{-}]");

  noiseCh_electrons->Draw("AP");
  noiseCh_electrons->GetXaxis()->SetTitle("Channel number");
  noiseCh_electrons->GetYaxis()->SetTitle("Sigma [e^{-}]");

  mpvStrip->Draw("AP");
  mpvStrip->GetXaxis()->SetTitle("Position [AU]");
  mpvStrip->GetYaxis()->SetTitle("Landau MPV [ADC]");

  sigmaStrip->Draw("AP");
  sigmaStrip->GetXaxis()->SetTitle("Position [AU]");
  sigmaStrip->GetYaxis()->SetTitle("Gaus sigma [ADC]");

  mpvStrip_norm->Draw("AP");
  mpvStrip_norm->GetXaxis()->SetTitle("Position [AU]");
  mpvStrip_norm->GetYaxis()->SetTitle("MPV slice / MPV 5 strips");

  maxDistrStrip_norm->Draw("AP");
  maxDistrStrip_norm->GetXaxis()->SetTitle("Position [AU]");
  maxDistrStrip_norm->GetYaxis()->SetTitle("Max distr slice / Max distr 5 strips");

  maxDistrStrip->Draw("AP");
  maxDistrStrip->GetXaxis()->SetTitle("Position [#mum]");
  maxDistrStrip->GetYaxis()->SetTitle("Max distr slice [ADC counts]");

  medianDistrStrip->Draw("AP");
  medianDistrStrip->GetXaxis()->SetTitle("Position [#mum]");
  medianDistrStrip->GetYaxis()->SetTitle("Median distr slice [ADC counts]");

  maxDistrStripL->Draw("AP");
  maxDistrStripL->GetXaxis()->SetTitle("Position [#mum]");
  maxDistrStripL->GetYaxis()->SetTitle("Max distr slice [ADC counts]");

  maxDistrStripR->Draw("AP");
  maxDistrStripR->GetXaxis()->SetTitle("Position [#mum]");
  maxDistrStripR->GetYaxis()->SetTitle("Max distr slice [ADC counts]");

  maxDistrStripLp1->Draw("AP");
  maxDistrStripLp1->GetXaxis()->SetTitle("Position [#mum]");
  maxDistrStripLp1->GetYaxis()->SetTitle("Max distr slice [ADC counts]");

  maxDistrStripRp1->Draw("AP");
  maxDistrStripRp1->GetXaxis()->SetTitle("Position [#mum]");
  maxDistrStripRp1->GetYaxis()->SetTitle("Max distr slice [ADC counts]");

  delete servCan;

  TCanvas* sigLRcan = new TCanvas("signalLeftRightSumCan");
  sigLRcan->SetGridx();
  sigLRcan->SetGridy();
  sigLR->Draw("APL");
  sigLR->GetXaxis()->SetTitle("Position [#mum]");
  sigLR->GetYaxis()->SetTitle("Max distr slice [ADC counts]");
  sigLRcan->BuildLegend()->SetLineColor(kWhite);
  sigLRcan->Modified();
  sigLRcan->Update();

  outFile->cd();

  trkEvt->Write();
  trkVsEvt->Write();
  trkEvtSelected->Write();
  timeTrack->Write();
  timeTrackPosY->Write();
  adcTime->Write();

  TDirectory* resDir = outFile->mkdir("Residuals");
  resDir->cd();
  residualsX->Write();
  residualsY->Write();
  residuals->Write();
  residualsXvsX->Write();
  profileResXvsX->Write();
  residualsYvsY->Write();
  profileResYvsY->Write();
  residualsXvsY->Write();
  profileResXvsY->Write();
  residualsYvsX->Write();
  profileResYvsX->Write();
  residualsYvsEvt->Write();
  residualsYvsEntry->Write();
  residualsYselected->Write();

  outFile->cd();
  hitMapDUTtele->Write();
  hitMapMatched->Write();
  hitMapDUTgoodCh->Write();
  hitMapLowPH->Write();
  trkEvtLowPH->Write();
  stripHPHDiffExtraLowPH->Write();
  diffExtraStripHPHvsPH->Write();
  extraChDistr->Write();
  extraChDistrGoodCh->Write();
  signalTime->Write();
  positivizedSignalTime->Write();
  positivizedSignalTimeProfile->Write();

  for(int iPar = 0; iPar < nPars; ++iPar)
    lanGausParVsTime[iPar]->Write();
   
  maxLanGauFit->Write();
  chi2SliceFit->Write();
  signalDistr->Write();
  noiseDistr->Write();
  noiseDistr_integral->Write();
  noiseDistrGroup->Write();
  noiseDistrGroup_electrons->Write();
  noiseDistrPair->Write();
  signalDistrTimeCut->Write();
  signalDistrTimeCutDistCut->Write();
  signalDistrTimeCutDistCut_electrons->Write();
  signalDistrTimeCutDistCutEta01->Write();
  signalDistrTimeCutDistCut_noisePeakSub->Write();
  signalDistrTimeCutDistCut_noisePeakSub_electrons->Write();
  signalDistrTimeDistHPHcut->Write();
  noiseDistrTimeCut->Write();
  backgroundDistrHPH->Write();
  stripHPHDistrTimeCut->Write();
  stripHPHDistrTimeCutDistCut->Write();
  stripHPHDistrTimeCutDistCut_BGsub->Write();
  stripHPH_BGsub_integral->Write();
  stripHPHDistrTimeCutDistCut_electrons->Write();
  backgroundDistrHPH_electrons->Write();
  stripHPHDistrTimeCutDistCut_BGsub_electrons->Write();
  stripHPH_BGsub_integral_electrons->Write();
  stripHPHDistrTimeCutDistCut_BGsub_fromFit->Write();
  stripHPH_plusNeigh_DistrTimeCutDistCut->Write();
  stripHPHDiffExtra->Write();
  stripHPHSignalTime->Write();
  leftStripHPHSignalTime->Write();
  rightStripHPHSignalTime->Write();
  stripHPHSignalTimeProfile->Write();
  leftStripHPHSignalTimeProfile->Write();
  rightStripHPHSignalTimeProfile->Write();
  phAroundHPHstripTimeCut->Write();
  phAroundExtraStripTimeCut->Write();
  correlationPHstripHPHhit->Write();
  tempEvt->Write();
  tempDistr->Write();
  scaleFactorDistr->Write();
  hitMapMod160->Write();
  chargeMapMod160->Write();
  chargeMapMod160Normalized->Write();
  chargeMapMod160NormProjX->Write();
  chargeMapMod160NormProjY->Write();
  hitMap->Write();
  chargeMap->Write();
  chargeMapNormalized->Write();
  chargeMapNormProjX->Write();
  chargeMapNormProjY->Write();
  signalStrip->Write();
  signalStripL->Write();
  signalStripR->Write();
  signalStripLp1->Write();
  signalStripRp1->Write();
  entriesStrip->Write();
  mpvStrip->Write();
  sigmaStrip->Write();
  mpvStrip_norm->Write();
  maxDistrStrip_norm->Write();
  maxDistrStrip->Write();
  medianDistrStrip->Write();
  maxDistrStripL->Write();
  maxDistrStripR->Write();
  maxDistrStripLp1->Write();
  maxDistrStripRp1->Write();
  sigLR->Write();
  sigLRcan->Write();
  etaDistrTimeCutDistCut->Write();
  etaDistrTrackTimeCut->Write();
  etaDistrTrackTimeCutLowPH->Write();
  CDFetaDistrTimeCutDistCut->Write();
  CDFetaDistrTrackTimeCut->Write();
  etaClustVsPos->Write();
  etaTrackVsPos->Write();
  profileEtaClust->Write();
  profileEtaTrack->Write();

  TDirectory* noiseDir = outFile->mkdir("noiseChannels");
  noiseDir->cd();
  for(int i = 0; i < nChannels; ++i)
    noiseHistCh[i]->Write();

  TDirectory* noiseDir_electrons = outFile->mkdir("noiseChannels_electrons");
  noiseDir_electrons->cd();
  for(int i = 0; i < nChannels; ++i)
    noiseHistCh_electrons[i]->Write();

  outFile->cd();
  noiseMeanCh->Write();
  noiseCh->Write();
  noiseChiCh->Write();
  fittedNoiseDistr->Write();
  noiseMeanCh_electrons->Write();
  noiseCh_electrons->Write();
  fittedNoiseDistr_electrons->Write();

  outFile->Close();
  /*
  sprintf(name, "chargeBetweenStripsRun_%i.dat", findRunNumber(argv[1]));
  std::ofstream outGr(name, std::ofstream::out);
  int np = maxDistrStrip->GetN();
  double* x = maxDistrStrip->GetX();
  double* y = maxDistrStrip->GetY();
  double* Ex = maxDistrStrip->GetEX();
  double* Ey = maxDistrStrip->GetEY();

  for(int i = 0; i < np; ++i)
    outGr << x[i] << ' ' << y[i] << ' ' << Ex[i] << ' ' << Ey[i] << '\n';

  outGr.close();
  */
  return 0;
}
