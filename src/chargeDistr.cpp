#include "iostream"
#include "stdlib.h"
#include "string"
#include "vector"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH2D.h"
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
  const double pitch = 0.080; // mm

  TDirectory* dir = (TDirectory*) inFile->Get("WriteTracksToNTuple");
  TTree* trkTree = (TTree*) dir->Get("EUFit");

  ConfigFileReader* conf = new ConfigFileReader(argv[2]);
  //conf->DumpConfMap();

  TString outFileName = conf->GetValue("outFilePath");
  outFileName += findRunNumber(argv[1]);
  outFileName += ".root";

  int maxDist = atoi(conf->GetValue("hitMaxDist").c_str()); // maximum distance form the strip hit by a track in looking at the charge, in number of strips

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
   trkTree->SetBranchAddress("dutTrackX",&dutTrackX); // position extrapolated from the track fit
   trkTree->SetBranchAddress("dutTrackY",&dutTrackY);
   trkTree->SetBranchAddress("dutHitX",&dutHitX); // cluster matched on the alibava
   trkTree->SetBranchAddress("dutHitY",&dutHitY);
   trkTree->SetBranchAddress("dutHitR",&dutHitR);
   trkTree->SetBranchAddress("dutHitQ",&dutHitQ);
   trkTree->SetBranchAddress("alibava_TDC",&alibava_TDC);
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
   trkTree->SetBranchStatus("dutHitX", 1);
   trkTree->SetBranchStatus("dutHitY", 1);

   long int nEntries = trkTree->GetEntries();

   TH1I* trkEvt = new TH1I("traksEvt", "Number of tracks per event;Number of tracks;Entries", 11, -0.5, 10.5);
   TH2D* hitMapDUTtele = new TH2D("hitMapDUTtele", "Extrapolated position of the tracks on the strip sensor;x [mm];y [mm]", 200, -20, 20, 100, -10, 10);
   TH2D* hitMapMatched = new TH2D("hitMapMatched", "Matched hits on the strip sensor;x [mm];y [mm]", 200, -20, 20, 100, -10, 10);
   TH2D* hitMapDUTgoodCh = new TH2D("hitMapDUTgoodCh", "Extrapolated position of the tracks on the strip sensor, passing a good channel;x [mm];y [mm]", 200, -20, 20, 100, -10, 10);
   TH1I* extraChDistr = new TH1I("extraChDistr", "Distribution of the extrapolated position in channels;Channel;Entries", 256, -0.5, 255.5);
   TH1I* extraChDistrGoodCh = new TH1I("extraChDistrGoodCh", "Distribution of the extrapolated position in channels (only good channels shown);Channel;Entries", 256, -0.5, 255.5);

   int nTrks = 0; // number of tracks in one event
   long int evtMrk = -1; // event marker

   std::vector<track> trkVec; // tracks in one event
   double evtAliPH[nChannels] = {0}; // alibava pulse height in on e event
   float evtAliTime = -1;
   float evtAliTemp = -275;
   track* trk;

   int extraCh = -1; // extrapolated channel number

   for(int i = 0; i < nEntries; ++i)
     {
       trkTree->GetEntry(i);

       // put the track info into a structure
       trk = new track();
       trk->extraPosDUT[0] = dutTrackX;
       trk->extraPosDUT[1] = dutTrackY;
       trk->measPosDUT[0] = dutHitX;
       trk->measPosDUT[1] = dutHitY;

       if(evtMrk == EvtNr)
	 { // add track info to some container
	   nTrks++;
	   hitMapDUTtele->Fill(dutTrackX, dutTrackY);
	   if(dutHitX > -900 && dutHitY > -900) hitMapMatched->Fill(dutHitX, dutHitY);

	   trkVec.push_back(*trk);
	 }
       else // new event
	 { // analyze the old event
	   trkEvt->Fill(nTrks);
	   for(unsigned int iTrk = 0; iTrk < trkVec.size(); ++iTrk)
	     {
	       // extraCh = (right variable for this) / pitch; // check this!!!
	       extraCh = trkVec[iTrk].extraPosDUT[1] / pitch; // to be changed

	       extraChDistr->Fill(extraCh);

	       if(evtAliPH[extraCh] != 0)
		 {
		   hitMapDUTgoodCh->Fill(trkVec[iTrk].extraPosDUT[0], trkVec[iTrk].extraPosDUT[1]);
		   extraChDistrGoodCh->Fill(extraCh);
		 }


	     }

	   // set the counters
	   evtMrk = EvtNr;
	   trkVec.clear();
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

   trkEvt->Write();
   hitMapDUTtele->Write();
   hitMapMatched->Write();
   hitMapDUTgoodCh->Write();
   extraChDistr->Write();
   extraChDistrGoodCh->Write();

   outFile->Close();

  return 0;
}
