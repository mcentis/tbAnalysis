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
  //const double pitch = 0.080; // mm

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
   TH2D* positivizedSignalTime = new TH2D("positivizedSignalTime", "Hit signal vs time (positivized);Time [ns];Hit signal [ADC]", 60, 0, 120, 251, -50.5, 511.5);
   TH1D* signalDistr = new TH1D("signalDistr", "Hit signal distribution (positivized);Hit signal[ADC];Entries", 562, -50.5, 511.5);
   TH1D* signalDistrTimeCut = new TH1D("signalDistrTimeCut", "Hit signal distribution (positivized) in the time cut;Hit signal[ADC];Entries", 251, -50.5, 511.5);
   TGraph* tempEvt = new TGraph();
   tempEvt->SetName("tempEvt");
   tempEvt->SetTitle("Tempetrature of the beetle chip vs event number");

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

		   if(fabs(hitCharge) > fabs(highestCharge)) highestCharge = hitCharge; // select the highest hit charge
		 }
	       // fill histos old event
	       signalTime->Fill(evtAliTime, highestCharge);
	       positivizedSignalTime->Fill(evtAliTime, highestCharge * polarity);
	       signalDistr->Fill(highestCharge * polarity);
	       if(evtAliTime >= timeCut1 && evtAliTime <= timeCut2) signalDistrTimeCut->Fill(highestCharge * polarity);
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

   TDirectory* timeSlicesDir = outFile->mkdir("timeSlices");
   timeSlicesDir->cd();   

   TGraphErrors* mpvTime = new TGraphErrors();
   mpvTime->SetName("mpvTime");
   mpvTime->SetTitle("Fitted Landau MPV vs time");

   TGraph* chi2SliceFit = new TGraph();
   chi2SliceFit->SetName("chi2SliceFit");
   chi2SliceFit->SetTitle("Reduced #chi^{2} of the landau gaussian fits");

   int nBins = positivizedSignalTime->GetXaxis()->GetNbins();
   double time = 0;
   double binW = 0;
   TH1D* slice = NULL;
   TF1* fit = NULL;
   char title[200];

   TCanvas* fitCan = new TCanvas("fitCan");

   lanGausFit(signalDistrTimeCut, negSigmaFit, posSigmaFit);

   for(int iBin = 1; iBin <= nBins; ++iBin) // fit the signal distributions of the various times
     {
       time = positivizedSignalTime->GetXaxis()->GetBinCenter(iBin);
       binW = positivizedSignalTime->GetXaxis()->GetBinWidth(iBin);

       sprintf(name, "timeSlice_%f", time);
       sprintf(title, "Hit charge distribution (positivized) at time %f;Hit charge [ADC];Entries", time);

       slice = positivizedSignalTime->ProjectionY(name, iBin, iBin);
       slice->SetTitle(title);

       fit = lanGausFit(slice, negSigmaFit, posSigmaFit);
       mpvTime->SetPoint(iBin - 1, time, fit->GetParameter(1));
       mpvTime->SetPointError(iBin - 1, binW / 2, fit->GetParError(1));

       if(fit->GetChisquare() > 0 && fit->GetNDF() > 0) // avoid to put nonsense in the graph
	 chi2SliceFit->SetPoint(chi2SliceFit->GetN(), time, fit->GetChisquare() / fit->GetNDF());

       slice->Write();
     }

   delete fitCan;

   // draw graphs to name the axis
   TCanvas* servCan = new TCanvas("servCan");
   servCan->cd();

   tempEvt->Draw("AP");
   tempEvt->GetXaxis()->SetTitle("Event number");
   tempEvt->GetYaxis()->SetTitle("Temperature [C]");

   mpvTime->Draw("AP");
   mpvTime->GetXaxis()->SetTitle("Time [ns]");
   mpvTime->GetYaxis()->SetTitle("MPV [ADC]");

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
   mpvTime->Write();
   chi2SliceFit->Write();
   signalDistr->Write();
   signalDistrTimeCut->Write();
   tempEvt->Write();

   outFile->Close();

  return 0;
}
