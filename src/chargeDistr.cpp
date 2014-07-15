#include "iostream"
#include "stdlib.h"
#include "string"
#include "vector"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TProfile.h"

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
  TF1* fitFunc = new TF1("fitFunc", "gaus + [3]", range1, range2);
  fitFunc->SetParNames("Const", "Mean", "Sigma", "Offset");
  fitFunc->SetLineColor(kRed);
  fitFunc->SetParameter(0, inHist->GetMaximum());
  fitFunc->SetParameter(1, inHist->GetMean()); // the distr is expected to be symmetric
  fitFunc->SetParameter(2, inHist->GetRMS() * 0.5);
  fitFunc->SetParameter(3, 0); // no offset

  fitFunc->SetParLimits(0, 0, inHist->GetEntries());
  fitFunc->SetParLimits(1, inHist->GetMean() - inHist->GetRMS(), inHist->GetMean() + inHist->GetRMS());
  fitFunc->SetParLimits(2, 0, 2 * inHist->GetRMS());
  fitFunc->SetParLimits(3, 0, inHist->GetEntries());

  inHist->Fit(fitFunc, "RQ");

  return fitFunc;
}

int main(int argc, char* argv[])
{
  gSystem->Load("libTree"); // necessary to use the program, root magic

  if(argc != 3)
    {
      std::cout << "Usage: chargeDistr inFile confFile" << std::endl;
      return 1;
    }

  TFile* inFile = TFile::Open(argv[1]);
  if (!inFile)
    {
      std::cout << "Impossible to open " << argv[1] << std::endl;
      return -1;
    }

  const int nChannels = 256;
  const double pitch = 0.080; // mm
  char name[50]; // to be used in varios occasions
  char title[200];

  TDirectory* dir = (TDirectory*) inFile->Get("WriteTracksToNTuple");
  TTree* trkTree = (TTree*) dir->Get("EUFit");

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

  TFile* outFile = new TFile(outFileName, "RECREATE");

  //Declaration of leaves types
  Int_t           Event;
  Int_t           RunNr;
  Int_t           EvtNr;
  Int_t           Ndf;
  Float_t         Chi2;
  Double_t        measX_0;
  Double_t        measY_0;
  Double_t        measZ_0;
  Double_t        measQ_0;
  Double_t        fitX_0;
  Double_t        fitY_0;
  Double_t        measX_1;
  Double_t        measY_1;
  Double_t        measZ_1;
  Double_t        measQ_1;
  Double_t        fitX_1;
  Double_t        fitY_1;
  Double_t        measX_2;
  Double_t        measY_2;
  Double_t        measZ_2;
  Double_t        measQ_2;
  Double_t        fitX_2;
  Double_t        fitY_2;
  Double_t        measX_3;
  Double_t        measY_3;
  Double_t        measZ_3;
  Double_t        measQ_3;
  Double_t        fitX_3;
  Double_t        fitY_3;
  Double_t        measX_4;
  Double_t        measY_4;
  Double_t        measZ_4;
  Double_t        measQ_4;
  Double_t        fitX_4;
  Double_t        fitY_4;
  Double_t        measX_5;
  Double_t        measY_5;
  Double_t        measZ_5;
  Double_t        measQ_5;
  Double_t        fitX_5;
  Double_t        fitY_5;
  Double_t        measX_6;
  Double_t        measY_6;
  Double_t        measZ_6;
  Double_t        measQ_6;
  Double_t        fitX_6;
  Double_t        fitY_6;
  Double_t        dutTrackX;
  Double_t        dutTrackY;
  Double_t        dutPixelX;
  Double_t        dutPixelY;
  Double_t        dutHitX;
  Double_t        dutHitY;
  Double_t        dutHitR;
  Double_t        dutHitQ;
  Float_t         alibava_TDC;
  Float_t         alibava_temp;
  Double_t        alibavaPH[nChannels];

  // Set branch addresses.
  trkTree->SetBranchAddress("Event",&Event);
  trkTree->SetBranchAddress("RunNr",&RunNr);
  trkTree->SetBranchAddress("EvtNr",&EvtNr);
  trkTree->SetBranchAddress("Ndf",&Ndf);
  trkTree->SetBranchAddress("Chi2",&Chi2);
  trkTree->SetBranchAddress("measX_0",&measX_0);
  trkTree->SetBranchAddress("measY_0",&measY_0);
  trkTree->SetBranchAddress("measZ_0",&measZ_0);
  trkTree->SetBranchAddress("measQ_0",&measQ_0);
  trkTree->SetBranchAddress("fitX_0",&fitX_0);
  trkTree->SetBranchAddress("fitY_0",&fitY_0);
  trkTree->SetBranchAddress("measX_1",&measX_1);
  trkTree->SetBranchAddress("measY_1",&measY_1);
  trkTree->SetBranchAddress("measZ_1",&measZ_1);
  trkTree->SetBranchAddress("measQ_1",&measQ_1);
  trkTree->SetBranchAddress("fitX_1",&fitX_1);
  trkTree->SetBranchAddress("fitY_1",&fitY_1);
  trkTree->SetBranchAddress("measX_2",&measX_2);
  trkTree->SetBranchAddress("measY_2",&measY_2);
  trkTree->SetBranchAddress("measZ_2",&measZ_2);
  trkTree->SetBranchAddress("measQ_2",&measQ_2);
  trkTree->SetBranchAddress("fitX_2",&fitX_2);
  trkTree->SetBranchAddress("fitY_2",&fitY_2);
  trkTree->SetBranchAddress("measX_3",&measX_3);
  trkTree->SetBranchAddress("measY_3",&measY_3);
  trkTree->SetBranchAddress("measZ_3",&measZ_3);
  trkTree->SetBranchAddress("measQ_3",&measQ_3);
  trkTree->SetBranchAddress("fitX_3",&fitX_3);
  trkTree->SetBranchAddress("fitY_3",&fitY_3);
  trkTree->SetBranchAddress("measX_4",&measX_4);
  trkTree->SetBranchAddress("measY_4",&measY_4);
  trkTree->SetBranchAddress("measZ_4",&measZ_4);
  trkTree->SetBranchAddress("measQ_4",&measQ_4);
  trkTree->SetBranchAddress("fitX_4",&fitX_4);
  trkTree->SetBranchAddress("fitY_4",&fitY_4);
  trkTree->SetBranchAddress("measX_5",&measX_5);
  trkTree->SetBranchAddress("measY_5",&measY_5);
  trkTree->SetBranchAddress("measZ_5",&measZ_5);
  trkTree->SetBranchAddress("measQ_5",&measQ_5);
  trkTree->SetBranchAddress("fitX_5",&fitX_5);
  trkTree->SetBranchAddress("fitY_5",&fitY_5);
  trkTree->SetBranchAddress("measX_6",&measX_6);
  trkTree->SetBranchAddress("measY_6",&measY_6);
  trkTree->SetBranchAddress("measZ_6",&measZ_6);
  trkTree->SetBranchAddress("measQ_6",&measQ_6);
  trkTree->SetBranchAddress("fitX_6",&fitX_6);
  trkTree->SetBranchAddress("fitY_6",&fitY_6);
  trkTree->SetBranchAddress("dutTrackX",&dutTrackX); // position extrapolated from the track fit on the DUT, global ref frame
  trkTree->SetBranchAddress("dutTrackY",&dutTrackY);
  trkTree->SetBranchAddress("dutPixelX",&dutPixelX); // position extrapolated on the DUT, dut ref frame, in pixel / strip number
  trkTree->SetBranchAddress("dutPixelY",&dutPixelY);
  trkTree->SetBranchAddress("dutHitX",&dutHitX); // cluster matched on the alibava
  trkTree->SetBranchAddress("dutHitY",&dutHitY);
  trkTree->SetBranchAddress("dutHitR",&dutHitR);
  trkTree->SetBranchAddress("dutHitQ",&dutHitQ);
  trkTree->SetBranchAddress("alibava_tdc",&alibava_TDC);
  trkTree->SetBranchAddress("alibava_temp",&alibava_temp);

  for(int i = 0; i < nChannels; ++i)
    {
      sprintf(name, "alibava_reco_ch_%i", i);
      trkTree->SetBranchAddress(name,&alibavaPH[i]);
    }

  // activate the interesting branches
  trkTree->SetBranchStatus("*", 0);
  trkTree->SetBranchStatus("EvtNr", 1);
  trkTree->SetBranchStatus("alibava*", 1); // all the alibava info
  trkTree->SetBranchStatus("dutTrackX", 1);
  trkTree->SetBranchStatus("dutTrackY", 1);
  trkTree->SetBranchStatus("dutPixelX", 1);
  trkTree->SetBranchStatus("dutPixelY", 1);
  trkTree->SetBranchStatus("dutHitX", 1);
  trkTree->SetBranchStatus("dutHitY", 1);

  // tracks
  TH1I* trkEvt = new TH1I("traksEvt", "Number of tracks per event;Number of tracks;Entries", 21, -0.5, 20.5);
  TGraph* trkVsEvt = new TGraph();
  trkVsEvt->SetName("trkVsEvt");
  trkVsEvt->SetTitle("Number of tracks vs event");
  trkVsEvt->SetMarkerStyle(7);
  TH1I* trkEvtSelected = new TH1I("traksEvtSelected", "Number of tracks per event in events that pass the event selection;Number of tracks;Entries", 11, -0.5, 10.5);

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
  TH1D* noiseDistrPair = new TH1D("noiseDistrPair", "Signal distribution (positivized) not associated with a hit, summed over 2 channels ;Signal [ADC];Entries", 201, -100.5, 100.5);
  TH1D* signalDistrTimeCut = new TH1D("signalDistrTimeCut", "Hit signal distribution (positivized) in the time cut;Hit signal[ADC];Entries", 562, -50.5, 511.5);
  TH1D* signalDistrTimeCutDistCut = new TH1D("signalDistrTimeCutDistCut", "Hit signal distribution (positivized) in the time cut, highest PH strip neighboring the extrapolated one;Hit signal[ADC];Entries", 562, -50.5, 511.5);
  TH1D* noiseDistrTimeCut = new TH1D("noiseDistrTimeCut", "Signal distribution (positivized) not associated with a hit in the time cut;Signal [ADC];Entries", 201, -100.5, 100.5);

  // signal on the strip with highest ph
  TH1D* stripHPHDistrTimeCut = new TH1D("stripHPHDistrTimeCut", "Hit signal distribution (positivized) in the time cut, for the strip with highest charge;Hit signal[ADC];Entries", 562, -50.5, 511.5);
  TH1D* stripHPHDistrTimeCutDistCut = new TH1D("stripHPHDistrTimeCutDistCut", "Hit signal distribution (positivized) in the time cut, for the strip with highest charge, highest PH strip neighboring the extrapolated one;Hit signal[ADC];Entries", 562, -50.5, 511.5);
  TH1D* backgroundDistrHPH = new TH1D("backgrounDistrHPH", "Backgorund distribution of the strip with highest PH;Signal [ADC];Entries", 562, -50.5, 511.5);
  TH1D* stripHPH_plusNeigh_DistrTimeCutDistCut = new TH1D("stripHPH_plusNeigh_DistrTimeCutDistCut", "Hit signal distribution (positivized) in the time cut, for the strip with highest charge plus its highest neighbor, highest PH strip neighboring the extrapolated one;Hit signal[ADC];Entries", 562, -50.5, 511.5);
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

  // three histos to get a map of charge collection over 2 strips in the time cut
  double minX = -20000;
  double maxX = 20000;
  int binX = 100;
  double minY = -0.5;
  double maxY = 160.5;
  int binY = 31;
  TH2D* chargeMapMod160Normalized = new TH2D("chargeMapMod160Normalized", "Charge map in the time cut divided by the number of tracks;x [#mum];y mod 160 [#mum];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);
  TH2D* chargeMapMod160 = new TH2D("chargeMapMod160", "Charge map in the time cut;x [#mum];y mod 160 [#mum];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);
  TH2D* hitMapMod160 = new TH2D("hitMapMod160", "Hit map in the time cut;x [#mum];y mod 160 [#mum];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);

  TH1D* etaDistrTimeCut = new TH1D("etaDistrTimeCut", "#eta distribution in the time cut;#eta;Entries", 200, -0.5, 1.5);

  // charge collection map over the sensor in the time cut
  minX = -20;
  maxX = 20;
  binX = 200;
  minY = -10;
  maxY = 10;
  binY = 100;
  TH2D* chargeMapNormalized = new TH2D("chargeMapNormalized", "Charge map in the time cut divided by the number of tracks;x [mm];y [mm];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);
  TH2D* chargeMap = new TH2D("chargeMap", "Charge map in the time cut;x [mm];y [mm];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);
  TH2D* hitMap = new TH2D("hitMap", "Hit map in the time cut;x [mm];y [mm];Charge [ADC]", binX, minX, maxX, binY, minY, maxY);

  TH1D* noiseHistCh[nChannels]; // calculation of the noise
  for(int i = 0; i < nChannels; ++i)
    {
      sprintf(name, "noiseDistrChannel_%i", i);
      sprintf(title, "Signal distribution (positivized) not associated with a hit, channel %i;Signal [ADC];Entries", i);
      noiseHistCh[i] = new TH1D(name, title, 111, -40.5, 40.5);
    }
  TH1D* fittedNoiseDistr = new TH1D("fittedNoiseDistr", "Distribution of the fitted noise;Noise [ADC];Entries", 61, -0.5, 20.5);

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

  bool goodNoise[nChannels]; // determine wether a channel was hit or not, for the noise analysis
  double noiseSum; // used for the noise over multiple channels
  int summed; // number of summed ch for the noise

  int highestPHstrip = -1; // strip with the highest ph in the hit
  double phHighestStrip = 0; // charge on the strip with the highest charge

  double noiseHPH = 0; // to calculate the background for the signal distribution of the strip with highest ph

  double phR = -1; // variables for the eta distribution
  double phL = -1;

  bool analyzeEvent = false;

  long int nEntries = trkTree->GetEntries();

  if(maxEntryNum == 0 || maxEntryNum > nEntries) // if there is a maximum number of entries
    maxEntryNum = nEntries;

  for(int i = 0; i < maxEntryNum; ++i)
    {
      trkTree->GetEntry(i);

      // put the track info into a structure
      trk = new track();
      trk->extraPosDUT[0] = dutTrackX; // global reference frame
      trk->extraPosDUT[1] = dutTrackY;
      trk->extraPosDUTpix[0] = dutPixelX; // dut ref frame, in pixel / strip number
      trk->extraPosDUTpix[1] = dutPixelY;
      trk->measPosDUT[0] = dutHitX;
      trk->measPosDUT[1] = dutHitY;
      trk->entryNum = i; 

      hitMapDUTtele->Fill(dutTrackX, dutTrackY);
      extraChDistr->Fill(dutPixelY);

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
		  if(trkVec.at(iTrk).measPosDUT[0] > -900 && trkVec.at(iTrk).measPosDUT[1] > -900) // if there is a matched hit
		    {
		      residualsYvsEntry->Fill(trkVec.at(iTrk).entryNum, trkVec.at(iTrk).measPosDUT[1] - trkVec.at(iTrk).extraPosDUT[1]);
		      hitMapMatched->Fill(trkVec.at(iTrk).measPosDUT[0], trkVec.at(iTrk).measPosDUT[1]);
		      residualsX->Fill(trkVec.at(iTrk).measPosDUT[0] - trkVec.at(iTrk).extraPosDUT[0]);
		      residualsY->Fill(trkVec.at(iTrk).measPosDUT[1] - trkVec.at(iTrk).extraPosDUT[1]);
		      residuals->Fill(trkVec.at(iTrk).measPosDUT[0] - trkVec.at(iTrk).extraPosDUT[0], trkVec.at(iTrk).measPosDUT[1] - trkVec.at(iTrk).extraPosDUT[1]);
		      residualsXvsX->Fill(trkVec.at(iTrk).extraPosDUT[0], trkVec.at(iTrk).measPosDUT[0] - trkVec.at(iTrk).extraPosDUT[0]);
		      residualsYvsY->Fill(trkVec.at(iTrk).extraPosDUT[1], trkVec.at(iTrk).measPosDUT[1] - trkVec.at(iTrk).extraPosDUT[1]);
		      residualsXvsY->Fill(trkVec.at(iTrk).extraPosDUT[1], trkVec.at(iTrk).measPosDUT[0] - trkVec.at(iTrk).extraPosDUT[0]);
		      residualsYvsX->Fill(trkVec.at(iTrk).extraPosDUT[0], trkVec.at(iTrk).measPosDUT[1] - trkVec.at(iTrk).extraPosDUT[1]);
		      residualsYvsEvt->Fill(evtMrk, trkVec.at(iTrk).measPosDUT[1] - trkVec.at(iTrk).extraPosDUT[1]);
		    }

		  extraCh = trkVec.at(iTrk).extraPosDUTpix[1];
		  for(int iCh = extraCh - maxDist; iCh <= extraCh + maxDist; ++iCh) // exclude channels that may have charge deposit from the noise analysis
		    if(iCh >= 0 && iCh < nChannels) // protect limits
		      goodNoise[iCh] = false;
		} // loop on the tracks

	      for(int iCh = 0; iCh < nChannels; ++iCh) // fill noise histos
	      	if(goodNoise[iCh] && evtAliPH[iCh] != 0) // no ph == 0 and no ch belonging toany hit
	      	  {
	      	    noiseDistr->Fill(evtAliPH[iCh] * polarity);
	      	    noiseHistCh[iCh]->Fill(evtAliPH[iCh] * polarity);
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
		      backgroundDistrHPH->Fill(noiseHPH);
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
	      extraCh = trkVec.at(iTrk).extraPosDUTpix[1];

	      if(trkVec.at(iTrk).extraPosDUT[0] <= xCut1 || trkVec.at(iTrk).extraPosDUT[0] >= xCut2) // geom cut in X
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
	      for(unsigned int iTrk = 0; iTrk < trkVec.size(); ++iTrk) // loop to select the right track
	      	{
		  extraCh = trkVec.at(iTrk).extraPosDUTpix[1]; // this should be right thing

		  hitMapDUTgoodCh->Fill(trkVec.at(iTrk).extraPosDUT[0], trkVec.at(iTrk).extraPosDUT[1]);
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
		    {
		      signalDistrTimeCutDistCut->Fill(highestCharge * polarity);
		      stripHPHDistrTimeCutDistCut->Fill(phHighestStrip);
		      correlationPHstripHPHhit->Fill(highestCharge * polarity, phHighestStrip);

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
		      etaDistrTimeCut->Fill(phR / (phR + phL));

		      stripHPH_plusNeigh_DistrTimeCutDistCut->Fill(phR + phL);
		    }

		  for(int iCh = 0; iCh < nChannels; ++iCh)
		    if(evtAliPH[iCh] != 0 && !(iCh >= hiChargeCh - maxDist && iCh <= hiChargeCh - maxDist)) // no ph == 0 and no ch belonging to the hit
		      noiseDistrTimeCut->Fill(evtAliPH[iCh] * polarity);

		  // strip with highest ph in the hit
		  stripHPHDistrTimeCut->Fill(phHighestStrip);
		  stripHPHDiffExtra->Fill(hiChargeCh - highestPHstrip);


		  if(highestCharge * polarity < 15) // investigation of the events with low ph
		    {
		      hitMapLowPH->Fill(trkVec.at(trackPos).extraPosDUT[0], trkVec.at(trackPos).extraPosDUT[1]);
		      trkEvtLowPH->Fill(trkVec.size());
		      stripHPHDiffExtraLowPH->Fill(hiChargeCh - highestPHstrip);
		    }
		  diffExtraStripHPHvsPH->Fill(highestCharge * polarity, hiChargeCh - highestPHstrip);

		  for(int iCh = hiChargeCh - maxDist; iCh <= hiChargeCh + maxDist; ++iCh)// ph around strip hph and extrapolated strip
		    if(iCh >=0 && iCh < nChannels) // protect array margins
		      {
			phAroundHPHstripTimeCut->Fill(iCh - highestPHstrip, evtAliPH[iCh] * polarity);
			phAroundExtraStripTimeCut->Fill(iCh - hiChargeCh, evtAliPH[iCh] * polarity);
		      }

		  // hitmap and charge map mod 160 (on 2 strips)
		  posX = 1000 * trkVec.at(trackPos).extraPosDUT[0]; // assign the positions in x and y
		  posY = abs((int)(1000 * trkVec.at(trackPos).extraPosDUT[1]) % (int)(2 * pitch * 1000));
		  iBin = chargeMapMod160->FindBin(posX, posY); // find the bin
		  oldContent = chargeMapMod160->GetBinContent(iBin);
		  chargeMapMod160->SetBinContent(iBin, oldContent + highestCharge * polarity);
		       
		  hitMapMod160->Fill(posX, posY);
		       
		  //hitmap and charge map over the sensor
		  posX = trkVec.at(trackPos).extraPosDUT[0];
		  posY = trkVec.at(trackPos).extraPosDUT[1];
		  iBin = chargeMap->FindBin(posX, posY); // find the bin
		  oldContent = chargeMap->GetBinContent(iBin);
		  chargeMap->SetBinContent(iBin, oldContent + highestCharge * polarity);
		       
		  hitMap->Fill(posX, posY);
		       
		}
	      tempEvt->SetPoint(tempEvt->GetN(), evtMrk, evtAliTemp);
	    } // the analysis of the event should be contained in this scope

	  // set the counters
	  evtMrk = EvtNr;
	  trkVec.clear();
	  highestCharge = 0;
	  // store the first track of the event
	  trkVec.push_back(*trk);
	  nTrks = 1;
	  // store the alibava info of the event
	  for(int iCh = 0; iCh < nChannels; ++iCh) evtAliPH[iCh] = alibavaPH[iCh];
	  evtAliTime = alibava_TDC;
	  evtAliTemp = alibava_temp;
	}

      delete trk;
    } // loop on the tracks (entries in the tree)

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

  // gausLanGausFitFixGaus(signalDistrTimeCutDistCut, negSigmaFit, posSigmaFit,
  //           		noiseDistrGroup->GetMean(), noiseDistrGroup->GetRMS()); // gaus mean and sigma determined from the noise distr and landau gauss convolution fitted simultaneously

  //lanGausFit(signalDistrTimeCut, negSigmaFit, posSigmaFit);
  //lanGausFit(signalDistrTimeCutDistCut, negSigmaFit, posSigmaFit);
  lanGausFit(stripHPHDistrTimeCutDistCut, negSigmaFit, posSigmaFit);
  lanGausFit(stripHPH_plusNeigh_DistrTimeCutDistCut, negSigmaFit, posSigmaFit);

  gausLanGausFit(signalDistrTimeCutDistCut, negSigmaFit, posSigmaFit); // peack at 0 and landau gauss convolution fitted simultaneously
  gausLanGausFit(signalDistrTimeCut, negSigmaFit, posSigmaFit);
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
      }

  // fit noise distr of the noise in groups of channels
  fit = new TF1("gausFit", "gaus", noiseDistrGroup->GetMean() - 2 * noiseDistrGroup->GetRMS(), noiseDistrGroup->GetMean() + 1 * noiseDistrGroup->GetRMS());
  noiseDistrGroup->Fit(fit, "RQ");

  // fit noise distr for a noise on pairs of channels
  fit = new TF1("gausFit", "gaus", noiseDistrPair->GetMean() - 2 * noiseDistrPair->GetRMS(), noiseDistrPair->GetMean() + 1 * noiseDistrPair->GetRMS());
  noiseDistrPair->Fit(fit, "RQ");

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
  int binStop = backgroundDistrHPH->GetXaxis()->FindBin(bgMean + bgRMS);
  double bgInt = backgroundDistrHPH->Integral(binStart, binStop);
  double sigPlusBgInt = stripHPHDistrTimeCutDistCut->Integral(binStart, binStop);
  TH1D* stripHPHDistrTimeCutDistCut_BGsub = new TH1D(*stripHPHDistrTimeCutDistCut);
  stripHPHDistrTimeCutDistCut_BGsub->SetName("stripHPHDistrTimeCutDistCut_BGsub");
  stripHPHDistrTimeCutDistCut_BGsub->SetTitle("Signal distr strip Hi PH time cut dist cut, bg subtracted;Signal [ADC];Entries");

  backgroundDistrHPH->Sumw2();
  stripHPHDistrTimeCutDistCut_BGsub->Sumw2();
  stripHPHDistrTimeCutDistCut_BGsub->Add(backgroundDistrHPH, -1 * sigPlusBgInt / bgInt);

  // normalized histo from the stripHPH with cuts and bg subtraction
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

  delete servCan;

  outFile->cd();

  trkEvt->Write();
  trkVsEvt->Write();
  trkEvtSelected->Write();

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
  noiseDistrGroup->Write();
  noiseDistrPair->Write();
  signalDistrTimeCut->Write();
  signalDistrTimeCutDistCut->Write();
  noiseDistrTimeCut->Write();
  backgroundDistrHPH->Write();
  stripHPHDistrTimeCut->Write();
  stripHPHDistrTimeCutDistCut->Write();
  stripHPHDistrTimeCutDistCut_BGsub->Write();
  stripHPH_BGsub_integral->Write();
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
  etaDistrTimeCut->Write();

  TDirectory* noiseDir = outFile->mkdir("noiseChannels");
  noiseDir->cd();
  for(int i = 0; i < nChannels; ++i)
    noiseHistCh[i]->Write();

  outFile->cd();
  noiseMeanCh->Write();
  noiseCh->Write();
  noiseChiCh->Write();
  fittedNoiseDistr->Write();

  outFile->Close();

  return 0;
}
