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

   char name[50];
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

   long int nEntries = trkTree->GetEntries();

   TH1I* trkEvt = new TH1I("traksEvt", "Number of tracks per event;Number of tracks;Entries", 11, -0.5, 10.5);
   TH2D* hitMapDUTtele = new TH2D("hitMapDUTtele", "Extrapolated position of the tracks on the strip sensor;x [mm];y [mm]", 200, -20, 20, 100, -10, 10);
   TH2D* hitMapMatched = new TH2D("hitMapMatched", "Matched hits on the strip sensor;x [mm];y [mm]", 200, -20, 20, 100, -10, 10);
   TH2D* hitMapDUTgoodCh = new TH2D("hitMapDUTgoodCh", "Extrapolated position of the tracks on the strip sensor, passing a good channel;x [mm];y [mm]", 200, -20, 20, 100, -10, 10);
   TH1I* extraChDistr = new TH1I("extraChDistr", "Distribution of the extrapolated position in channels;Channel;Entries", 513, -255.5, 255.5);
   TH1I* extraChDistrGoodCh = new TH1I("extraChDistrGoodCh", "Distribution of the extrapolated position in channels (only good channels shown);Channel;Entries", 256, -0.5, 255.5);
   TH2D* signalTime = new TH2D("signalTime", "Hit signal vs time;Time [ns];Hit signal [ADC]", 60, 0, 120, 1024, -511.5, 511.5);
   TH2D* positivizedSignalTime = new TH2D("positivizedSignalTime", "Hit signal vs time (positivized);Time [ns];Hit signal [ADC]", 60, 0, 120, 151, -50.5, 511.5);
   TH1D* signalDistr = new TH1D("signalDistr", "Hit signal distribution (positivized);Hit signal[ADC];Entries", 562, -50.5, 511.5);
   TH1D* signalDistrTimeCut = new TH1D("signalDistrTimeCut", "Hit signal distribution (positivized) in the time cut;Hit signal[ADC];Entries", 151, -50.5, 511.5);
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


   int nTrks = 0; // number of tracks in one event
   long int evtMrk = -1; // event marker

   std::vector<track> trkVec; // tracks in one event
   double evtAliPH[nChannels] = {0}; // alibava pulse height in on e event
   float evtAliTime = -1;
   float evtAliTemp = -275;
   track* trk;

   int extraCh = -1; // extrapolated channel number
   double hitCharge = 0; // charge on the hit
   double highestCharge = 0; // highest hit charge in the event (believed to be the particle that passes the detector in time)
   int trackPos = 0; // position of the track in the track vector
   int iBin = 0; // bin to be filled
   double posX = 0; // store positions
   double posY = 0; // this one will be the mod of the position
   double oldContent = 0;

   bool analyzeEvent = false;

   for(int i = 0; i < nEntries; ++i)
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

       hitMapDUTtele->Fill(dutTrackX, dutTrackY);
       if(dutHitX > -900 && dutHitY > -900) hitMapMatched->Fill(dutHitX, dutHitY);
       extraChDistr->Fill(dutPixelY);

       if(evtMrk == EvtNr)
	 { // add track info to some container
	   nTrks++;

	   trkVec.push_back(*trk);
	 }
       else // new event
	 { // analyze the old event, if there are tracks
	   analyzeEvent = true;
	   for(unsigned int iTrk = 0; iTrk < trkVec.size(); ++iTrk) // check that all the tracks are in the geom cut (in Y)
	     {
	       extraCh = trkVec.at(iTrk).extraPosDUTpix[1];
	       if(extraCh < 0 || extraCh >= nChannels) // protect array margins (no seg violation hopefully)
		 {
		   analyzeEvent = false;
		   break;
		 }

	       if(evtAliPH[extraCh] != 0 && evtAliPH[extraCh + 1] != 0 && evtAliPH[extraCh - 1] != 0) continue; // the strip traversed is not a border strip or a not bonded one
	       else
		 {
		   analyzeEvent = false;
		   break;
		 }
	     }

	   if(analyzeEvent && trkVec.size() != 0) analyzeEvent = true; // check that there is at least a track 
	   else analyzeEvent = false;

	   if(analyzeEvent) // all the event tracks in a sensitive part of the sensor and at least one track
	     {
	       trkEvt->Fill(nTrks);
	       for(unsigned int iTrk = 0; iTrk < trkVec.size(); ++iTrk)
		 {
		   extraCh = trkVec.at(iTrk).extraPosDUTpix[1]; // this should be right thing

		   hitMapDUTgoodCh->Fill(trkVec.at(iTrk).extraPosDUT[0], trkVec.at(iTrk).extraPosDUT[1]);
		   extraChDistrGoodCh->Fill(extraCh);

		   hitCharge = 0;

		   for(int dist = 0; dist <= maxDist; ++dist) // sum the charge of the strips around the extrapolated value
		     {
		       if(dist == 0) hitCharge += evtAliPH[extraCh];
		       else
			 {
			   hitCharge += evtAliPH[extraCh + dist];
			   hitCharge += evtAliPH[extraCh - dist];
			 }
		     }

		   if(fabs(hitCharge) > fabs(highestCharge))  // select the highest hit charge
		     {
		       highestCharge = hitCharge;
		       trackPos = iTrk;
		     }
		 }
	       // fill histos old event
	       if(trkVec.at(trackPos).extraPosDUT[0] > xCut1 && trkVec.at(trackPos).extraPosDUT[0] < xCut2)
		 {
		   signalTime->Fill(evtAliTime, highestCharge);
		   positivizedSignalTime->Fill(evtAliTime, highestCharge * polarity);
		   signalDistr->Fill(highestCharge * polarity);
		   if(evtAliTime >= timeCut1 && evtAliTime <= timeCut2) // apply time cut
		     {
		       signalDistrTimeCut->Fill(highestCharge * polarity);
		       
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
     }

   TProfile* positivizedSignalTimeProfile = positivizedSignalTime->ProfileX("positivizedSignalTimeProfile");
   positivizedSignalTimeProfile->SetTitle("Positivized signal time profile");

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
   char title[200];

   TCanvas* fitCan = new TCanvas("fitCan");

   //lanGausFit(signalDistrTimeCut, negSigmaFit, posSigmaFit);
   gausLanGausFit(signalDistrTimeCut, negSigmaFit, posSigmaFit); // peack at 0 and landau gauss convolution fitted simultaneously

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
	   lanGausParVsTime[iPar]->SetPoint(iBin - 1, time, fit->GetParameter(iPar));
	   lanGausParVsTime[iPar]->SetPointError(iBin - 1, binW / 2, fit->GetParError(iPar));
	 }

       if(fit->GetMaximumX() > -200) // exclude nonsense
	 maxLanGauFit->SetPoint(maxLanGauFit->GetN(), time, fit->GetMaximumX());

       if(fit->GetChisquare() > 0 && fit->GetNDF() > 0) // avoid to put nonsense in the graph
	 chi2SliceFit->SetPoint(chi2SliceFit->GetN(), time, fit->GetChisquare() / fit->GetNDF());

       slice->Write();
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

   // draw graphs to name the axis
   TCanvas* servCan = new TCanvas("servCan");
   servCan->cd();

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

   delete servCan;

   outFile->cd();

   trkEvt->Write();
   hitMapDUTtele->Write();
   hitMapMatched->Write();
   hitMapDUTgoodCh->Write();
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
   signalDistrTimeCut->Write();
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

   outFile->Close();

  return 0;
}
