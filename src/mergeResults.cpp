#include "iostream"
#include "fstream"
#include "vector"
#include "math.h"

#include "TStyle.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

#include "ConfigFileReader.hh"

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
  if(argc != 4)
    {
      std::cout << "Usage: mergeResults listOfLists dataDir confFile" << std::endl;
      return 1;
    }

  ConfigFileReader* conf = new ConfigFileReader(argv[3]);
  conf->DumpConfMap();

  double ADCtoe = atof(conf->GetValue("ADCtoe").c_str()); // parameter to apply conversion of ADC to e-
  double ADCtoeErr = atof(conf->GetValue("ADCtoeErr").c_str());

  double tCorr_p0 = atof(conf->GetValue("tCorr_p0").c_str()); // parameters to apply temperature correction
  double tCorr_p1 = atof(conf->GetValue("tCorr_p1").c_str());
  double tCorr_p0Err = atof(conf->GetValue("tCorr_p0Err").c_str());
  double tCorr_p1Err = atof(conf->GetValue("tCorr_p1Err").c_str());

  double targetChipTemp = atof(conf->GetValue("targetChipTemp").c_str());
  double tempErr = atof(conf->GetValue("tempErr").c_str());
  double targetGain;
  double targetGainErr;

  if(tCorr_p1 || tCorr_p0)
    {
      targetGain = tCorr_p0 + tCorr_p1 * targetChipTemp;
      targetGainErr = sqrt(pow(tCorr_p0Err, 2) + pow(tCorr_p1Err * targetChipTemp, 2));
    }
  else
    {
      targetGain = 1;
      targetGainErr = 1;
      std::cout << "================================================> no gain error calculated!!!" << std::endl;
    }

  // big axis labels
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");

  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");

  gStyle->SetTitleOffset(0.95, "x");
  gStyle->SetTitleOffset(0.95, "y");

  gStyle->SetMarkerSize(2);
  // gStyle->SetLineWidth(1);

  char fileName[200];
  sprintf(fileName, "%s/mergedData.root", argv[2]);
  TFile* outFile = new TFile(fileName, "RECREATE");

  // sensors that has to be merged
  std::vector<std::string> sensorType; // es epi100P
  std::vector<double> fluences; // in neq / cm^2
  std::vector<std::string> runList; // file with the run list for each sensor

  std::string type;
  double flu;
  std::string list;

  // read the list of lists
  std::ifstream listStr(argv[1], std::ifstream::in);
  if(listStr.is_open() == false)
    {
      std::cout << "Could not open " << argv[1] << std::endl;
      return -1;
    }

  do
    {
      listStr >> type >> flu >> list;
      sensorType.push_back(type);
      fluences.push_back(flu);
      runList.push_back(list);
    } while(listStr.good());
    
  listStr.close();
  
  std::cout << "Lists for sensors and fluences\n";
  for(unsigned int i = 0; i < sensorType.size(); ++i)
    std::cout << sensorType.at(i) << '\t' << fluences.at(i) << '\t' << runList.at(i) << std::endl;
  std::cout << "If too many sensors appear, have a look for empty lines in the lists" << std::endl;

  std::vector<std::string> sensorLabel; // labels for captions
  bool longLabel = false;
  for(unsigned int i = 0; i < sensorType.size(); ++i)
    {
      if(sensorType.at(i)[6] == 'P' || sensorType.at(i)[6] == 'p')
	{
	  if(longLabel) sensorLabel.push_back("p-bulk, p-stop");
	  else sensorLabel.push_back("P");
	}
      else if(sensorType.at(i)[6] == 'Y' || sensorType.at(i)[6] == 'y')
	{
	  if(longLabel) sensorLabel.push_back("p-bulk, p-spray");
	  else sensorLabel.push_back("Y");
	}
      else if(sensorType.at(i)[6] == 'N' || sensorType.at(i)[6] == 'n')
	{
	  if(longLabel) sensorLabel.push_back("n-bulk         ");
	  else sensorLabel.push_back("N");
	}
      else
	sensorLabel.push_back("   no label    ");
    }

  std::vector<std::string> sensorThickness; // determine thickness and material
  std::vector<std::string> sensorMaterial;
  std::string str;
  for(unsigned int i = 0; i < sensorType.size(); ++i)
    {
      str = std::string(sensorType.at(i), 0, 3);
      sensorMaterial.push_back(str);

      str = std::string(sensorType.at(i), 3, 3);
      sensorThickness.push_back(str);
    }

  std::vector<TGraphErrors*> medianBiasVec;
  std::vector<TGraphErrors*> meanBiasVec;
  std::vector<TGraphErrors*> maxFitBiasVec;
  std::vector<TGraphErrors*> mpvBiasVec;
  std::vector<TGraphErrors*> mpvBiasVec_electrons;
  std::vector<TGraphErrors*> lanWBiasVec;
  std::vector<TGraphErrors*> gSigBiasVec;
  std::vector<TGraphErrors*> noiseBiasVec;
  std::vector<TGraphErrors*> noiseBiasVec_electrons;
  std::vector<TGraphErrors*> noiseGroupBiasVec;
  std::vector<TGraphErrors*> noisePairBiasVec;
  std::vector<TGraphErrors*> eff95BiasVec;
  std::vector<TGraphErrors*> chargeSharingBiasVec;
  std::vector<TGraphErrors*> resYBiasVec;
  std::vector<TGraphErrors*> chipTempBiasVec;

  std::vector<TCanvas*> histSupVec; // canvases with the superimposition of the used histos
  std::vector<TCanvas*> etaSupVec; // canvases with the superimposition of the used eta distr

  std::vector<TMultiGraph*> chargePosSupVec; // suiperimposition of the charge vs pos

  std::vector<double> bias;
  std::vector<double> angle;
  std::vector<int> run;
  // variables used for reading
  double bi;
  double an;
  int ru;

  TFile* inFile;
  TH1* chDistr;
  TF1* func;
  TH1* noiseDistr;
  TH1* noiseGroupDistr;
  TH1* noisePairDistr;
  TH1* intSeedDistr;
  TH1* etaDistr;
  TH1* resYDistr;
  TH1* tempDistr;
  TDirectory* resDir;
  TGraphErrors* medianGr;
  TGraphErrors* meanGr;
  TGraphErrors* maxFitGr;
  TGraphErrors* mpvGr;
  TGraphErrors* mpvGr_electrons;
  TGraphErrors* lanWGr;
  TGraphErrors* gSigGr;
  TGraphErrors* noiseGr;
  TGraphErrors* noiseGr_electrons;
  TGraphErrors* noiseGroupGr;
  TGraphErrors* noisePairGr;
  TGraphErrors* eff95Gr;
  TGraphErrors* chargeSharingGr;
  TGraphErrors* resYGr;
  TGraphErrors* chipTempGr;
  TGraphErrors* chargePos;
  TCanvas* histSupCan; // each sensor gets a canvas (charge distribution)
  TCanvas* etaSupCan; // each sensor gets a canvas (eta distribution)
  TMultiGraph* chargePosSup; // each sensor gets a graph

  char name[200];
  char title[500];
  int linStyle = 1;
  int linWidth = 1;
  int iColor = 0; // color of the graphs
  int mrkStyle = 0; // marker style
  int iColSuper = 0; // color for superimpositions

  int startBin;
  int endBin;
  double chargeSharing;
  double EchargeSharing;

  for(unsigned int i = 0; i < sensorType.size(); ++i) // loop on the sensors
    {
      bias.clear();
      angle.clear();
      run.clear();

      std::ifstream runListStr(runList.at(i).c_str(), std::ifstream::in);
      if(runListStr.is_open() == false)
	{
	  std::cout << "Could not open " << runList.at(i) << std::endl;
	  continue;
	}

      do
	{
	  runListStr >> bi >> an >> ru;
	  bias.push_back(bi);
	  angle.push_back(an);
	  run.push_back(ru);
	} while(runListStr.good());

      runListStr.close();

      std::cout << "=======================================================\n";
      std::cout << "Lists of runs for sensor " << sensorType.at(i) << "    " << fluences.at(i) << std::endl;
      for(unsigned int j = 0; j < bias.size(); ++j)
	std::cout << bias.at(j) << '\t' << angle.at(j) << '\t' << run.at(j) << std::endl;
      std::cout << "If too many runs appear, have a look for empty lines in the lists" << std::endl;

      // iColor = i % 6 + 1; // use "good" colors
      // if(iColor >= 5) ++iColor; // skip yellow
      // if(iColor - 1 == 0) linStyle++;

      linStyle = 1;
      // the color will represent the irradiation
      if(fluences.at(i) == 0)
	{
	  iColor = kRed;
	  linStyle = 2; // dotted line for 0 fluence
	}
      else if(fluences.at(i) == 1.e15)
	iColor = kGreen;// in case it is too bright use 8
      else if(fluences.at(i) == 1.5e15)
	iColor = kBlack;
      else if(fluences.at(i) == 3e15)
	iColor = kBlue;
      else if(fluences.at(i) == 1.3e16)
	iColor = kRed;
      else
	{
	  iColor = 6;
	  std::cout << "-------------------------------------------------------------------- Fluence not recognized for " << sensorType.at(i) << std::endl;
	}

      // the symbol represents the material
      if(sensorMaterial.at(i) == "Epi" || sensorMaterial.at(i) == "epi")
	{
	  mrkStyle = 20;
	  linWidth = 1;
	}
      else if(sensorMaterial.at(i) == "Mcz" || sensorMaterial.at(i) == "mcz")
	{
	  mrkStyle = 22;
	  linWidth = 2;
	}
      else if(sensorMaterial.at(i) == "Fth" || sensorMaterial.at(i) == "fth")
	{
	  mrkStyle = 21;
	  linWidth = 2;
	}
      else
	{
	  mrkStyle = 3;
	  std::cout << "-------------------------------------------------------------------- Material not recognized for " << sensorType.at(i) << std::endl;
	}

      // open symbols for p-stop, full for p-spray, full for n type
      if(sensorType.at(i)[6] == 'P' || sensorType.at(i)[6] == 'p')
	mrkStyle += 4;
      else if(sensorType.at(i)[6] != 'Y' && sensorType.at(i)[6] != 'y' && sensorType.at(i)[6] != 'N' && sensorType.at(i)[6] != 'n')
	{
	  mrkStyle = 27;
	  std::cout << "-------------------------------------------------------------------- Sensor type not recognized for " << sensorType.at(i) << std::endl;
	}

      // if(sensorType.at(i)[6] == 'P' || sensorType.at(i)[6] == 'p')
      // 	; // do nothing
      // else if(sensorType.at(i)[6] == 'Y' || sensorType.at(i)[6] == 'y')
      // 	mrkStyle += 4;
      // else if(sensorType.at(i)[6] == 'N' || sensorType.at(i)[6] == 'n')
      // 	mrkStyle = 29;
      // else
      // 	{
      // 	  mrkStyle = 27;
      // 	  std::cout << "-------------------------------------------------------------------- Sensor type not recognized for " << sensorType.at(i) << std::endl;
      // 	}

      sprintf(name, "mpv_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      mpvGr = new TGraphErrors();
      mpvGr->SetName(name);
      mpvGr->SetTitle(title);
      mpvGr->SetMarkerStyle(mrkStyle);
      mpvGr->SetFillColor(kWhite);
      mpvGr->SetLineColor(iColor); // set line color and style
      mpvGr->SetLineWidth(linWidth);
      mpvGr->SetMarkerColor(iColor);
      mpvGr->SetLineStyle(linStyle);

      sprintf(name, "mpv_electrons_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      mpvGr_electrons = new TGraphErrors();
      mpvGr_electrons->SetName(name);
      mpvGr_electrons->SetTitle(title);
      mpvGr_electrons->SetMarkerStyle(mrkStyle);
      mpvGr_electrons->SetFillColor(kWhite);
      mpvGr_electrons->SetLineColor(iColor); // set line color and style
      mpvGr_electrons->SetLineWidth(linWidth);
      mpvGr_electrons->SetMarkerColor(iColor);
      mpvGr_electrons->SetLineStyle(linStyle);

      sprintf(name, "mean_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      meanGr = new TGraphErrors();
      meanGr->SetName(name);
      meanGr->SetTitle(title);
      meanGr->SetMarkerStyle(mrkStyle);
      meanGr->SetFillColor(kWhite);
      meanGr->SetLineColor(iColor); // set line color and style
      meanGr->SetLineWidth(linWidth);
      meanGr->SetMarkerColor(iColor);
      meanGr->SetLineStyle(linStyle);

      sprintf(name, "median_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      medianGr = new TGraphErrors();
      medianGr->SetName(name);
      medianGr->SetTitle(title);
      medianGr->SetMarkerStyle(mrkStyle);
      medianGr->SetFillColor(kWhite);
      medianGr->SetLineColor(iColor); // set line color and style
      medianGr->SetLineWidth(linWidth);
      medianGr->SetMarkerColor(iColor);
      medianGr->SetLineStyle(linStyle);

      sprintf(name, "maxFit_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      maxFitGr = new TGraphErrors();
      maxFitGr->SetName(name);
      maxFitGr->SetTitle(title);
      maxFitGr->SetMarkerStyle(mrkStyle);
      maxFitGr->SetFillColor(kWhite);
      maxFitGr->SetLineColor(iColor); // set line color and style
      maxFitGr->SetLineWidth(linWidth);
      maxFitGr->SetMarkerColor(iColor);
      maxFitGr->SetLineStyle(linStyle);

      sprintf(name, "lanW_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      lanWGr = new TGraphErrors();
      lanWGr->SetName(name);
      lanWGr->SetTitle(title);
      lanWGr->SetMarkerStyle(mrkStyle);
      lanWGr->SetFillColor(kWhite);
      lanWGr->SetLineColor(iColor); // set line color and style
      lanWGr->SetLineWidth(linWidth);
      lanWGr->SetMarkerColor(iColor);
      lanWGr->SetLineStyle(linStyle);

      sprintf(name, "gSig_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      gSigGr = new TGraphErrors();
      gSigGr->SetName(name);
      gSigGr->SetTitle(title);
      gSigGr->SetMarkerStyle(mrkStyle);
      gSigGr->SetFillColor(kWhite);
      gSigGr->SetLineColor(iColor); // set line color and style
      gSigGr->SetLineWidth(linWidth);
      gSigGr->SetMarkerColor(iColor);
      gSigGr->SetLineStyle(linStyle);

      sprintf(name, "noise_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      noiseGr = new TGraphErrors();
      noiseGr->SetName(name);
      noiseGr->SetTitle(title);
      noiseGr->SetMarkerStyle(mrkStyle);
      noiseGr->SetFillColor(kWhite);
      noiseGr->SetLineColor(iColor); // set line color and style
      noiseGr->SetLineWidth(linWidth);
      noiseGr->SetMarkerColor(iColor);
      noiseGr->SetLineStyle(linStyle);

      sprintf(name, "noise_electrons_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      noiseGr_electrons = new TGraphErrors();
      noiseGr_electrons->SetName(name);
      noiseGr_electrons->SetTitle(title);
      noiseGr_electrons->SetMarkerStyle(mrkStyle);
      noiseGr_electrons->SetFillColor(kWhite);
      noiseGr_electrons->SetLineColor(iColor); // set line color and style
      noiseGr_electrons->SetLineWidth(linWidth);
      noiseGr_electrons->SetMarkerColor(iColor);
      noiseGr_electrons->SetLineStyle(linStyle);

      sprintf(name, "noiseGroup_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      noiseGroupGr = new TGraphErrors();
      noiseGroupGr->SetName(name);
      noiseGroupGr->SetTitle(title);
      noiseGroupGr->SetMarkerStyle(mrkStyle);
      noiseGroupGr->SetFillColor(kWhite);
      noiseGroupGr->SetLineColor(iColor); // set line color and style
      noiseGroupGr->SetLineWidth(linWidth);
      noiseGroupGr->SetMarkerColor(iColor);
      noiseGroupGr->SetLineStyle(linStyle);

      sprintf(name, "noisePair_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      noisePairGr = new TGraphErrors();
      noisePairGr->SetName(name);
      noisePairGr->SetTitle(title);
      noisePairGr->SetMarkerStyle(mrkStyle);
      noisePairGr->SetFillColor(kWhite);
      noisePairGr->SetLineColor(iColor); // set line color and style
      noisePairGr->SetLineWidth(linWidth);
      noisePairGr->SetMarkerColor(iColor);
      noisePairGr->SetLineStyle(linStyle);

      sprintf(name, "eff95_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      eff95Gr = new TGraphErrors();
      eff95Gr->SetName(name);
      eff95Gr->SetTitle(title);
      eff95Gr->SetMarkerStyle(mrkStyle);
      eff95Gr->SetFillColor(kWhite);
      eff95Gr->SetLineColor(iColor); // set line color and style
      eff95Gr->SetLineWidth(linWidth);
      eff95Gr->SetMarkerColor(iColor);
      eff95Gr->SetLineStyle(linStyle);

      sprintf(name, "chargeSharing_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      chargeSharingGr = new TGraphErrors();
      chargeSharingGr->SetName(name);
      chargeSharingGr->SetTitle(title);
      chargeSharingGr->SetMarkerStyle(mrkStyle);
      chargeSharingGr->SetFillColor(kWhite);
      chargeSharingGr->SetLineColor(iColor); // set line color and style
      chargeSharingGr->SetLineWidth(linWidth);
      chargeSharingGr->SetMarkerColor(iColor);
      chargeSharingGr->SetLineStyle(linStyle);

      sprintf(name, "resY_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      resYGr = new TGraphErrors();
      resYGr->SetName(name);
      resYGr->SetTitle(title);
      resYGr->SetMarkerStyle(mrkStyle);
      resYGr->SetFillColor(kWhite);
      resYGr->SetLineColor(iColor); // set line color and style
      resYGr->SetLineWidth(linWidth);
      resYGr->SetMarkerColor(iColor);
      resYGr->SetLineStyle(linStyle);

      sprintf(name, "chipTemp_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      chipTempGr = new TGraphErrors();
      chipTempGr->SetName(name);
      chipTempGr->SetTitle(title);
      chipTempGr->SetMarkerStyle(mrkStyle);
      chipTempGr->SetFillColor(kWhite);
      chipTempGr->SetLineColor(iColor); // set line color and style
      chipTempGr->SetLineWidth(linWidth);
      chipTempGr->SetMarkerColor(iColor);
      chipTempGr->SetLineStyle(linStyle);

      sprintf(name, "histSup_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "Distributions at different biases %s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      histSupCan = new TCanvas(name, title);
      histSupCan->SetGridx();
      histSupCan->SetGridy();

      sprintf(name, "etaSup_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "#eta distributions at different biases %s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      etaSupCan = new TCanvas(name, title);
      etaSupCan->SetGridx();
      etaSupCan->SetGridy();

      sprintf(name, "chargePosSup_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "Charge vs position at different biases %s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      chargePosSup = new TMultiGraph(name, title);

      for(unsigned int iRun = 0; iRun < bias.size(); ++iRun) // loop on the runs for a sensor
	{
	  sprintf(name, "%s/%d.root", argv[2], run.at(iRun));
	  inFile = TFile::Open(name);

	  chDistr = (TH1*) inFile->Get("signalDistrTimeCutDistCut"); // full hit 
	  //chDistr = (TH1*) inFile->Get("stripHPHDistrTimeCutDistCut"); // strip highest ph
	  //chDistr = (TH1*) inFile->Get("stripHPH_plusNeigh_DistrTimeCutDistCut"); // strip highest ph and highest neighbor

	  // func = chDistr->GetFunction("lanGausFit");

	  // mpvGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(1));
	  // mpvGr->SetPointError(iRun, 0, func->GetParError(1));

	  // lanWGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(0));
	  // lanWGr->SetPointError(iRun, 0, func->GetParError(0));

	  // gSigGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(3));
	  // gSigGr->SetPointError(iRun, 0, func->GetParError(3));

	  func = chDistr->GetFunction("gausLang");

	  mpvGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(4));
	  mpvGr->SetPointError(iRun, 0, func->GetParError(4));

	  maxFitGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetMaximumX());

	  lanWGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(3));
	  lanWGr->SetPointError(iRun, 0, func->GetParError(3));

	  gSigGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(6));
	  gSigGr->SetPointError(iRun, 0, func->GetParError(6));

	  noiseDistr = (TH1*) inFile->Get("fittedNoiseDistr");

	  noiseGr->SetPoint(iRun, fabs(bias.at(iRun)), noiseDistr->GetMean());
	  noiseGr->SetPointError(iRun, 0, noiseDistr->GetRMS() / sqrt(noiseDistr->GetEntries()));

	  noiseGroupDistr = (TH1*) inFile->Get("noiseDistrGroup");

	  noiseGroupGr->SetPoint(iRun, fabs(bias.at(iRun)), noiseGroupDistr->GetRMS());
	  noiseGroupGr->SetPointError(iRun, 0, noiseGroupDistr->GetRMSError());

	  noisePairDistr = (TH1*) inFile->Get("noiseDistrPair");

	  noisePairGr->SetPoint(iRun, fabs(bias.at(iRun)), noisePairDistr->GetRMS());
	  noisePairGr->SetPointError(iRun, 0, noisePairDistr->GetRMSError());

	  intSeedDistr = (TH1*) inFile->Get("stripHPH_BGsub_integral"); // get the threshold to have 95 % efficiency
	  for(int iBin = 1; iBin < intSeedDistr->GetNbinsX(); ++iBin)
	    if(intSeedDistr->GetBinContent(iBin) > 0.05)
	      {
		eff95Gr->SetPoint(iRun, fabs(bias.at(iRun)), intSeedDistr->GetXaxis()->GetBinCenter(iBin));
		eff95Gr->SetPointError(iRun, 0, intSeedDistr->GetXaxis()->GetBinWidth(iBin));

		break;
	      }

	  //etaDistr = (TH1*) inFile->Get("etaDistrTimeCutDistCut"); // from the "clusters"
	  etaDistr = (TH1*) inFile->Get("etaDistrTrackTimeCut"); // track based
	  etaDistr->Sumw2();
	  etaDistr->Scale(1 / etaDistr->Integral());

	  startBin = etaDistr->GetXaxis()->FindBin(0.2);
	  endBin = etaDistr->GetXaxis()->FindBin(0.8);
	  chargeSharing = etaDistr->IntegralAndError(startBin, endBin, EchargeSharing) / 0.6;
	  EchargeSharing /= 0.6;

	  chargeSharingGr->SetPoint(iRun, fabs(bias.at(iRun)), chargeSharing);
	  chargeSharingGr->SetPointError(iRun, 0, EchargeSharing);

	  iColSuper = iRun % 9 + 1;
	  if(iColSuper >= 5) ++iColSuper; // skip yellow
	  etaDistr->SetLineColor(iColSuper);
	  sprintf(title, "%.00f", bias.at(iRun));
	  etaDistr->SetTitle(title);
	  //etaDistr->SetLineWidth(2);
	  etaSupCan->cd();
	  if(iRun == 0) etaDistr->Draw("E");
	  else etaDistr->Draw("Esame");
	  // if(iRun == 0) etaDistr->Draw("hist");
	  // else etaDistr->Draw("histsame");

	  chargePos = (TGraphErrors*) inFile->Get("mpvStrip");
	  chargePos->SetLineColor(iColSuper);
	  chargePos->SetFillColor(kWhite);
	  sprintf(title, "%.00f", bias.at(iRun));
	  chargePos->SetTitle(title);
	  chargePosSup->Add(chargePos);

	  resDir = (TDirectory*) inFile->Get("Residuals");
	  resYDistr = (TH1*) resDir->Get("residualsYselected");
	  func = resYDistr->GetFunction("fitFunc");

	  resYGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(2));
	  resYGr->SetPointError(iRun, 0, func->GetParError(2));

	  tempDistr = (TH1*) inFile->Get("tempDistr");
	  chipTempGr->SetPoint(iRun, fabs(bias.at(iRun)), tempDistr->GetMean());
	  chipTempGr->SetPointError(iRun, 0, tempDistr->GetRMS());

	  chDistr->SetLineColor(iColSuper);
	  sprintf(title, "%.00f", bias.at(iRun));
	  chDistr->SetTitle(title);
	  //chDistr->SetLineWidth(2);
	  histSupCan->cd();
	  if(iRun == 0) chDistr->Draw("E");
	  else chDistr->Draw("Esame");

	  // mean and median
	  TH1D* chDistr_nSub = (TH1D*) inFile->Get("signalDistrTimeCutDistCut_noisePeakSub"); // noise peak subtracted
	  meanGr->SetPoint(iRun, fabs(bias.at(iRun)), chDistr_nSub->GetMean());
	  meanGr->SetPointError(iRun, 0, chDistr_nSub->GetMeanError());

	  medianGr->SetPoint(iRun, fabs(bias.at(iRun)), Median(chDistr_nSub));

	  // charge in electrons
	  chDistr = (TH1*) inFile->Get("signalDistrTimeCutDistCut_electrons"); // full hit 
	  func = chDistr->GetFunction("gausLang");

	  double gainMeas;
	  double gainMeasErr;
	  double tempTotErr = sqrt(pow(tempErr, 2) + pow(tempDistr->GetMeanError(), 2));
	  if(tCorr_p1 || tCorr_p0)
	    {
	      gainMeas = tCorr_p0 + tCorr_p1 * tempDistr->GetMean();
	      gainMeasErr = sqrt(pow(tCorr_p0Err, 2) + pow(tCorr_p1Err * tempDistr->GetMean(), 2) + pow(tCorr_p1 * tempTotErr, 2));
	    }

	  double error = func->GetParameter(4) * sqrt(pow(ADCtoeErr / ADCtoe, 2) + pow(targetGainErr / targetGain, 2) + pow(gainMeasErr / gainMeas, 2) + pow(func->GetParError(4) / func->GetParameter(4), 2));

	  mpvGr_electrons->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(4));
	  mpvGr_electrons->SetPointError(iRun, 0, error);

	  noiseDistr = (TH1*) inFile->Get("fittedNoiseDistr_electrons");

	  error =  noiseDistr->GetMean() * sqrt(pow(ADCtoeErr / ADCtoe, 2) + pow(noiseDistr->GetMeanError() / noiseDistr->GetMean(), 2));
	  noiseGr_electrons->SetPoint(iRun, fabs(bias.at(iRun)), noiseDistr->GetMean());
	  noiseGr_electrons->SetPointError(iRun, 0, error);

	  //inFile->Close(); // do not close the files, otherwise the canvas does not find the histos to draw and save
	} // loop on the runs for a sensor


      mpvBiasVec.push_back(mpvGr);
      mpvBiasVec_electrons.push_back(mpvGr_electrons);
      meanBiasVec.push_back(meanGr);
      medianBiasVec.push_back(medianGr);
      maxFitBiasVec.push_back(maxFitGr);
      lanWBiasVec.push_back(lanWGr);
      gSigBiasVec.push_back(gSigGr);
      noiseBiasVec.push_back(noiseGr);
      noiseBiasVec_electrons.push_back(noiseGr_electrons);
      noiseGroupBiasVec.push_back(noiseGroupGr);
      noisePairBiasVec.push_back(noisePairGr);
      eff95BiasVec.push_back(eff95Gr);
      chargeSharingBiasVec.push_back(chargeSharingGr);
      resYBiasVec.push_back(resYGr);
      chipTempBiasVec.push_back(chipTempGr);
      histSupVec.push_back(histSupCan);
      etaSupVec.push_back(etaSupCan);
      chargePosSupVec.push_back(chargePosSup);

    } // loop on the sensors

  std::vector<TGraphErrors*> snrBiasVec;
  TGraphErrors* snrGr;

  double* bia;
  double* mpv;
  double* errMpv;
  double* noise;
  double* errNoise;

  double snr;
  double errSnr;

  for(unsigned int i = 0; i < sensorType.size(); ++i) // loop on the sensors
    {
      mpvGr = mpvBiasVec.at(i);
      noiseGr = noiseBiasVec.at(i); // mean noise on single channel
      //noiseGr = noiseGroupBiasVec.at(i); // noise on a group of strips
      //noiseGr = noisePairBiasVec.at(i); // noise on a pair of strips

      sprintf(name, "snr_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      snrGr = new TGraphErrors();
      snrGr->SetName(name);
      snrGr->SetTitle(title);
      snrGr->SetMarkerStyle(mpvGr->GetMarkerStyle());
      snrGr->SetFillColor(kWhite);
      snrGr->SetLineColor(mpvGr->GetLineColor()); // set line color and style
      snrGr->SetLineWidth(linWidth);
      snrGr->SetMarkerColor(mpvGr->GetMarkerColor());
      snrGr->SetLineStyle(mpvGr->GetLineStyle());
      snrGr->SetLineWidth(2);

      bia = mpvGr->GetX();
      mpv = mpvGr->GetY();
      errMpv = mpvGr->GetEY();
      noise = noiseGr->GetY();
      errNoise = noiseGr->GetEY();

      for(int iPoint = 0; iPoint < mpvGr->GetN(); ++iPoint)
      	{
      	  snr = mpv[iPoint] / noise[iPoint];
      	  errSnr = snr * sqrt(pow(errMpv[iPoint] / mpv[iPoint], 2) + pow(errNoise[iPoint] / noise[iPoint], 2));

      	  snrGr->SetPoint(iPoint, bia[iPoint], snr);
      	  snrGr->SetPointError(iPoint, 0, errSnr);
      	}

      snrBiasVec.push_back(snrGr);
    } // loop on the sensors

  // study of correlations
  std::vector<TGraphErrors*> corrMPVmaxFitVec;
  std::vector<TGraphErrors*> corrWidthGsigFitVec;
  std::vector<TGraphErrors*> corrNoiseGrGsigFitVec;

  TMultiGraph* corrMPVmaxFit = new TMultiGraph();
  corrMPVmaxFit->SetName("corrMPVmaxFit");
  corrMPVmaxFit->SetTitle("Max fit vs Landau MPV");

  TMultiGraph* corrWidthGsigFit = new TMultiGraph();
  corrWidthGsigFit->SetName("corrWidthGsigFit");
  corrWidthGsigFit->SetTitle("Gaus sigma fit vs Landau width");

  TMultiGraph* corrNoiseGrGsigFit = new TMultiGraph();
  corrNoiseGrGsigFit->SetName("corrNoiseGrGsigFit");
  corrNoiseGrGsigFit->SetTitle("Gaus sigma fit vs noise group (RMS)");

  double* y1;
  double* y2;
  double* eY1;
  double* eY2;
  TGraphErrors* corrGr;
  TGraphErrors* gr1;
  TGraphErrors* gr2;

  for(unsigned int i = 0; i < sensorType.size(); ++i) // loop on the sensors, for correlations
    {
      //------------------------------------------- landau mpv vs max of the fit
      gr1 = mpvBiasVec.at(i);
      gr2 = maxFitBiasVec.at(i);

      sprintf(name, "corrMPVmaxFit_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      corrGr = new TGraphErrors();
      corrGr->SetName(name);
      corrGr->SetTitle(title);
      corrGr->SetMarkerStyle(gr1->GetMarkerStyle());
      corrGr->SetFillColor(kWhite);
      corrGr->SetLineColor(gr1->GetLineColor()); // set line color and style
      corrGr->SetLineWidth(gr1->GetLineWidth());
      corrGr->SetMarkerColor(gr1->GetMarkerColor());
      corrGr->SetLineStyle(gr1->GetLineStyle());
      //corrGr->SetLineWidth(2);

      y1 = gr1->GetY();
      eY1 = gr1->GetEY();
      y2 = gr2->GetY();
      eY2 = gr2->GetEY();

      for(int iPoint = 0; iPoint < gr1->GetN(); ++iPoint)
	{
	  corrGr->SetPoint(iPoint, y1[iPoint], y2[iPoint]);
	  corrGr->SetPointError(iPoint, eY1[iPoint], eY2[iPoint]);
	}

      corrMPVmaxFitVec.push_back(corrGr);
      corrMPVmaxFit->Add(corrGr);

      // -------------------------------------------------- landau width gaus sigma of the fit
      gr1 = lanWBiasVec.at(i);
      gr2 = gSigBiasVec.at(i);

      sprintf(name, "corrLanWgausSigFit_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      corrGr = new TGraphErrors();
      corrGr->SetName(name);
      corrGr->SetTitle(title);
      corrGr->SetMarkerStyle(gr1->GetMarkerStyle());
      corrGr->SetFillColor(kWhite);
      corrGr->SetLineColor(gr1->GetLineColor()); // set line color and style
      corrGr->SetLineWidth(gr1->GetLineWidth());
      corrGr->SetMarkerColor(gr1->GetMarkerColor());
      corrGr->SetLineStyle(gr1->GetLineStyle());
      //corrGr->SetLineWidth(2);

      y1 = gr1->GetY();
      eY1 = gr1->GetEY();
      y2 = gr2->GetY();
      eY2 = gr2->GetEY();

      for(int iPoint = 0; iPoint < gr1->GetN(); ++iPoint)
	{
	  corrGr->SetPoint(iPoint, y1[iPoint], y2[iPoint]);
	  corrGr->SetPointError(iPoint, eY1[iPoint], eY2[iPoint]);
	}

      corrWidthGsigFitVec.push_back(corrGr);
      corrWidthGsigFit->Add(corrGr);

      // -------------------------------------------- noise of group strips gaus sigma fit
      gr1 = gSigBiasVec.at(i);
      gr2 = noiseGroupBiasVec.at(i);

      sprintf(name, "corrgSigNoiseGroup_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      corrGr = new TGraphErrors();
      corrGr->SetName(name);
      corrGr->SetTitle(title);
      corrGr->SetMarkerStyle(gr1->GetMarkerStyle());
      corrGr->SetFillColor(kWhite);
      corrGr->SetLineColor(gr1->GetLineColor()); // set line color and style
      corrGr->SetLineWidth(gr1->GetLineWidth());
      corrGr->SetMarkerColor(gr1->GetMarkerColor());
      corrGr->SetLineStyle(gr1->GetLineStyle());
      //corrGr->SetLineWidth(2);

      y1 = gr1->GetY();
      eY1 = gr1->GetEY();
      y2 = gr2->GetY();
      eY2 = gr2->GetEY();

      for(int iPoint = 0; iPoint < gr1->GetN(); ++iPoint)
	{
	  corrGr->SetPoint(iPoint, y1[iPoint], y2[iPoint]);
	  corrGr->SetPointError(iPoint, eY1[iPoint], eY2[iPoint]);
	}

      corrNoiseGrGsigFitVec.push_back(corrGr);
      corrNoiseGrGsigFit->Add(corrGr);
    }

  // normalization graph
  TGraphErrors* normGraph = new TGraphErrors();
  normGraph->SetName("normGraph");
  normGraph->SetTitle("Normalization graph, using non irradiated sensors");
  TF1* normFit = new TF1("normFit", "pol0", 150, 450); // function to fit the normalization constant
  int point;
  for(unsigned int i = 0; i < sensorType.size(); ++i) // loop on the sensors, for correlations
    if(fluences.at(i) == 0)
      {
	bia = mpvBiasVec.at(i)->GetX();
	mpv = mpvBiasVec.at(i)->GetY();
	errMpv = mpvBiasVec.at(i)->GetEY();

	for(int iPoint = 0; iPoint < mpvBiasVec.at(i)->GetN(); ++iPoint)
	  {
	    point = normGraph->GetN();
	    normGraph->SetPoint(point, bia[iPoint], mpv[iPoint]);
	    normGraph->SetPointError(point, 0, errMpv[iPoint]);
	  }
      }
  normGraph->Fit(normFit, "R");

  // normalized mpv graph
  TMultiGraph* normMpvGraph = new TMultiGraph("normMpvGraph", "MPV normalizet to non irradiated devices between 200 and 400 V");
  std::vector<TGraphErrors*> normMpvGrVec;
  TGraphErrors* normMpvGr;
  double norm = normFit->GetParameter(0);
  double normErr = normFit->GetParError(0);
  double normMpv;
  double normMpvErr;
  for(unsigned int i = 0; i < sensorType.size(); ++i) // loop on the sensors, for correlations
    {
      sprintf(name, "normMpvGr_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      normMpvGr = new TGraphErrors();
      normMpvGr->SetName(name);
      normMpvGr->SetTitle(title);
      normMpvGr->SetMarkerStyle(mpvBiasVec.at(i)->GetMarkerStyle());
      normMpvGr->SetFillColor(kWhite);
      normMpvGr->SetLineColor(mpvBiasVec.at(i)->GetLineColor()); // set line color and style
      normMpvGr->SetLineWidth(mpvBiasVec.at(i)->GetLineWidth());
      normMpvGr->SetMarkerColor(mpvBiasVec.at(i)->GetMarkerColor());
      normMpvGr->SetLineStyle(mpvBiasVec.at(i)->GetLineStyle());

      bia = mpvBiasVec.at(i)->GetX();
      mpv = mpvBiasVec.at(i)->GetY();
      errMpv = mpvBiasVec.at(i)->GetEY();

      for(int iPoint = 0; iPoint < mpvBiasVec.at(i)->GetN(); ++iPoint)
	{
	  normMpv = mpv[iPoint] / norm;
	  normMpvErr = sqrt(pow(errMpv[iPoint] / mpv[iPoint], 2) + pow(normErr / norm, 2)) * normMpv;
	  normMpvGr->SetPoint(iPoint, bia[iPoint], normMpv);
	  normMpvGr->SetPointError(iPoint, 0, normMpvErr);
	}

      normMpvGrVec.push_back(normMpvGr);
      normMpvGraph->Add(normMpvGr);
    }

  TMultiGraph* mpvAllSensors = new TMultiGraph();
  mpvAllSensors->SetName("mpvAllSensors");
  mpvAllSensors->SetTitle("MPV vs bias");

  TMultiGraph* mpvAll_electrons = new TMultiGraph();
  mpvAll_electrons->SetName("mpvAll_electrons");
  mpvAll_electrons->SetTitle("MPV vs bias");

  TMultiGraph* meanAllSensors = new TMultiGraph();
  meanAllSensors->SetName("meanAllSensors");
  meanAllSensors->SetTitle("Mean of the charge distribution (noise sub)");

  TMultiGraph* medianAllSensors = new TMultiGraph();
  medianAllSensors->SetName("medianAllSensors");
  medianAllSensors->SetTitle("Median of the charge distribution (noise sub)");

  TMultiGraph* maxFitAllSensors = new TMultiGraph();
  maxFitAllSensors->SetName("maxFitAllSensors");
  maxFitAllSensors->SetTitle("Maximum of the fit function");

  TMultiGraph* lanWAllSensors = new TMultiGraph();
  lanWAllSensors->SetName("lanWAllSensors");
  lanWAllSensors->SetTitle("Landau width vs bias");

  TMultiGraph* gSigAllSensors = new TMultiGraph();
  gSigAllSensors->SetName("gSigAllSensors");
  gSigAllSensors->SetTitle("Gaus #sigma vs bias");

  TMultiGraph* noiseAllSensors = new TMultiGraph();
  noiseAllSensors->SetName("noiseAllSensors");
  noiseAllSensors->SetTitle("Noise vs bias");

  TMultiGraph* noiseAll_electrons = new TMultiGraph();
  noiseAll_electrons->SetName("noiseAll_electrons");
  noiseAll_electrons->SetTitle("Noise vs bias");

  TMultiGraph* noiseGroupAllSensors = new TMultiGraph();
  noiseGroupAllSensors->SetName("noiseGroupAllSensors");
  noiseGroupAllSensors->SetTitle("Noise over multiple channels vs bias");

  TMultiGraph* noisePairAllSensors = new TMultiGraph();
  noisePairAllSensors->SetName("noisePairAllSensors");
  noisePairAllSensors->SetTitle("Noise over multiple channels vs bias");

  TMultiGraph* snrAllSensors = new TMultiGraph();
  snrAllSensors->SetName("snrAllSensors");
  snrAllSensors->SetTitle("SNR vs bias");

  TMultiGraph* eff95AllSensors = new TMultiGraph();
  eff95AllSensors->SetName("eff95AllSensors");
  eff95AllSensors->SetTitle("Threshold to get 95% efficiency vs bias");

  TMultiGraph* chargeSharingAllSensors = new TMultiGraph();
  chargeSharingAllSensors->SetName("chargeSharingAllSensors");
  chargeSharingAllSensors->SetTitle("Charge sharing vs bias");

  TMultiGraph* resYAllSensors = new TMultiGraph();
  resYAllSensors->SetName("resYAllSensors");
  resYAllSensors->SetTitle("#sigma resduals Y vs bias");

  TMultiGraph* chipTempAllSensors = new TMultiGraph();
  chipTempAllSensors->SetName("chipTempAllSensors");
  chipTempAllSensors->SetTitle("Chip temperature");

  TCanvas* servCan = new TCanvas();
  servCan->SetName("servCan");

  // limits of the x axis for the graphs with all sensors
  double xmin = 0;
  double xmax = 1100;

  // test of legend
  TLegend* legMpv = new TLegend(0.5, 0.5, 0.7, 0.7);
  legMpv->SetFillColor(kWhite);
  double prevFlu = -1;

  for(unsigned int i = 0; i < mpvBiasVec.size(); ++i) // loop on the graphs
    {
      mpvBiasVec.at(i)->Draw("AP");
      mpvBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      mpvBiasVec.at(i)->GetYaxis()->SetTitle("Landau MPV [ADC counts]");

      mpvAllSensors->Add(mpvBiasVec.at(i));

      if(fluences.at(i) != prevFlu)
	{
	  // TH1F* labFlu = new TH1F("labFlu","", 10, 0, 1);
	  // labFlu->SetLineColor(mpvBiasVec.at(i)->GetLineColor());
	  // labFlu->SetFillColor(mpvBiasVec.at(i)->GetLineColor());
	  // if(fluences.at(i) == 0) labFlu->SetFillColor(kWhite);
	  sprintf(title, "%.01e cm^{-2}", fluences.at(i));
	  // legMpv->AddEntry(labFlu, title);
	  legMpv->AddEntry((TObject*)0, title, "");
	  prevFlu = fluences.at(i);
	}

      sprintf(title, "%s %s %s", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str());
      legMpv->AddEntry(mpvBiasVec.at(i), title);
    }

  mpvAllSensors->Draw("AP");
  mpvAllSensors->GetXaxis()->SetTitle("Bias [V]");
  mpvAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  mpvAllSensors->GetYaxis()->SetTitle("Landau MPV [ADC counts]");

  for(unsigned int i = 0; i < mpvBiasVec_electrons.size(); ++i) // loop on the graphs
      mpvAll_electrons->Add(mpvBiasVec_electrons.at(i));

  mpvAll_electrons->Draw("AP");
  mpvAll_electrons->GetXaxis()->SetTitle("Bias [V]");
  mpvAll_electrons->GetXaxis()->SetLimits(xmin, xmax);
  mpvAll_electrons->GetYaxis()->SetTitle("Landau MPV [e^{-}]");

  for(unsigned int i = 0; i < meanBiasVec.size(); ++i) // loop on the graphs
    meanAllSensors->Add(meanBiasVec.at(i));

  meanAllSensors->Draw("AP");
  meanAllSensors->GetXaxis()->SetTitle("Bias [V]");
  meanAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  meanAllSensors->GetYaxis()->SetTitle("Mean [ADC counts]");

  for(unsigned int i = 0; i < medianBiasVec.size(); ++i) // loop on the graphs
    medianAllSensors->Add(medianBiasVec.at(i));

  medianAllSensors->Draw("AP");
  medianAllSensors->GetXaxis()->SetTitle("Bias [V]");
  medianAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  medianAllSensors->GetYaxis()->SetTitle("Median [ADC counts]");

  for(unsigned int i = 0; i < maxFitBiasVec.size(); ++i) // loop on the graphs
    {
      maxFitBiasVec.at(i)->Draw("AP");
      maxFitBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      maxFitBiasVec.at(i)->GetYaxis()->SetTitle("Max of the fit [ADC counts]");

      maxFitAllSensors->Add(maxFitBiasVec.at(i));
    }

  maxFitAllSensors->Draw("AP");
  maxFitAllSensors->GetXaxis()->SetTitle("Bias [V]");
  maxFitAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  maxFitAllSensors->GetYaxis()->SetTitle("Max of the fit [ADC counts]");

  for(unsigned int i = 0; i < lanWBiasVec.size(); ++i) // loop on the graphs
    {
      lanWBiasVec.at(i)->Draw("AP");
      lanWBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      lanWBiasVec.at(i)->GetYaxis()->SetTitle("Landau width [ADC counts]");

      lanWAllSensors->Add(lanWBiasVec.at(i));
    }

  lanWAllSensors->Draw("AP");
  lanWAllSensors->GetXaxis()->SetTitle("Bias [V]");
  lanWAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  lanWAllSensors->GetYaxis()->SetTitle("Landau width [ADC counts]");

  for(unsigned int i = 0; i < gSigBiasVec.size(); ++i) // loop on the graphs
    {
      gSigBiasVec.at(i)->Draw("AP");
      gSigBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      gSigBiasVec.at(i)->GetYaxis()->SetTitle("Gaus #sigma [ADC counts]");

      gSigAllSensors->Add(gSigBiasVec.at(i));
    }

  gSigAllSensors->Draw("AP");
  gSigAllSensors->GetXaxis()->SetTitle("Bias [V]");
  gSigAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  gSigAllSensors->GetYaxis()->SetTitle("Gaus #sigma [ADC counts]");

  for(unsigned int i = 0; i < noiseBiasVec.size(); ++i) // loop on the graphs
    {
      noiseBiasVec.at(i)->Draw("AP");
      noiseBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      noiseBiasVec.at(i)->GetYaxis()->SetTitle("Noise [ADC counts]");

      noiseAllSensors->Add(noiseBiasVec.at(i));
    }

  noiseAllSensors->Draw("AP");
  noiseAllSensors->GetXaxis()->SetTitle("Bias [V]");
  noiseAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  noiseAllSensors->GetYaxis()->SetTitle("Noise [ADC counts]");

  for(unsigned int i = 0; i < noiseBiasVec.size(); ++i) // loop on the graphs
    noiseAll_electrons->Add(noiseBiasVec_electrons.at(i));

  noiseAll_electrons->Draw("AP");
  noiseAll_electrons->GetXaxis()->SetTitle("Bias [V]");
  noiseAll_electrons->GetXaxis()->SetLimits(xmin, xmax);
  noiseAll_electrons->GetYaxis()->SetTitle("Noise [e^{-}]");

  for(unsigned int i = 0; i < noiseGroupBiasVec.size(); ++i) // loop on the graphs
    {
      noiseGroupBiasVec.at(i)->Draw("AP");
      noiseGroupBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      noiseGroupBiasVec.at(i)->GetYaxis()->SetTitle("Noise (RMS) [ADC counts]");

      noiseGroupAllSensors->Add(noiseGroupBiasVec.at(i));
    }

  noiseGroupAllSensors->Draw("AP");
  noiseGroupAllSensors->GetXaxis()->SetTitle("Bias [V]");
  noiseGroupAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  noiseGroupAllSensors->GetYaxis()->SetTitle("Noise (RMS) [ADC counts]");

  for(unsigned int i = 0; i < noisePairBiasVec.size(); ++i) // loop on the graphs
    {
      noisePairBiasVec.at(i)->Draw("AP");
      noisePairBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      noisePairBiasVec.at(i)->GetYaxis()->SetTitle("Noise (RMS) [ADC counts]");

      noisePairAllSensors->Add(noisePairBiasVec.at(i));
    }

  noisePairAllSensors->Draw("AP");
  noisePairAllSensors->GetXaxis()->SetTitle("Bias [V]");
  noisePairAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  noisePairAllSensors->GetYaxis()->SetTitle("Noise (RMS) [ADC counts]");

  for(unsigned int i = 0; i < snrBiasVec.size(); ++i) // loop on the graphs
    {
      snrBiasVec.at(i)->Draw("AP");
      snrBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      snrBiasVec.at(i)->GetYaxis()->SetTitle("SNR");

      snrAllSensors->Add(snrBiasVec.at(i));
    }

  snrAllSensors->Draw("AP");
  snrAllSensors->GetXaxis()->SetTitle("Bias [V]");
  snrAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  snrAllSensors->GetYaxis()->SetTitle("SNR");

  for(unsigned int i = 0; i < chargePosSupVec.size(); ++i) // loop on the graphs
    {
      chargePosSupVec.at(i)->Draw("AP");
      chargePosSupVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      chargePosSupVec.at(i)->GetYaxis()->SetTitle("Landau MPV [ADC counts]");
    }

  for(unsigned int i = 0; i < eff95BiasVec.size(); ++i) // loop on the graphs
    {
      eff95BiasVec.at(i)->Draw("AP");
      eff95BiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      eff95BiasVec.at(i)->GetYaxis()->SetTitle("Threshold 95% eff [ADC counts]");

      eff95AllSensors->Add(eff95BiasVec.at(i));
    }

  eff95AllSensors->Draw("AP");
  eff95AllSensors->GetXaxis()->SetTitle("Bias [V]");
  eff95AllSensors->GetXaxis()->SetLimits(xmin, xmax);
  eff95AllSensors->GetYaxis()->SetTitle("Threshold 95% eff [ADC counts]");

  for(unsigned int i = 0; i < chargeSharingBiasVec.size(); ++i) // loop on the graphs
    {
      chargeSharingBiasVec.at(i)->Draw("AP");
      chargeSharingBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      chargeSharingBiasVec.at(i)->GetYaxis()->SetTitle("Charge sharing");

      chargeSharingAllSensors->Add(chargeSharingBiasVec.at(i));
    }

  chargeSharingAllSensors->Draw("AP");
  chargeSharingAllSensors->GetXaxis()->SetTitle("Bias [V]");
  chargeSharingAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  chargeSharingAllSensors->GetYaxis()->SetTitle("Charge sharing");

  for(unsigned int i = 0; i < resYBiasVec.size(); ++i) // loop on the graphs
    {
      resYBiasVec.at(i)->Draw("AP");
      resYBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      resYBiasVec.at(i)->GetYaxis()->SetTitle("#sigma res Y [mm]");

      resYAllSensors->Add(resYBiasVec.at(i));
    }

  resYAllSensors->Draw("AP");
  resYAllSensors->GetXaxis()->SetTitle("Bias [V]");
  resYAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  resYAllSensors->GetYaxis()->SetTitle("#sigma res Y [mm]");

  for(unsigned int i = 0; i < chipTempBiasVec.size(); ++i) // loop on the graphs
    {
      chipTempBiasVec.at(i)->Draw("AP");
      chipTempBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      chipTempBiasVec.at(i)->GetYaxis()->SetTitle("Temperature [C]");

      chipTempAllSensors->Add(chipTempBiasVec.at(i));
    }

  chipTempAllSensors->Draw("AP");
  chipTempAllSensors->GetXaxis()->SetTitle("Bias [V]");
  chipTempAllSensors->GetXaxis()->SetLimits(xmin, xmax);
  chipTempAllSensors->GetYaxis()->SetTitle("Temperature [C]");

  // correlations graphs
  corrMPVmaxFit->Draw("AP");
  corrMPVmaxFit->GetXaxis()->SetTitle("Landau MPV [ADC counts]");
  corrMPVmaxFit->GetYaxis()->SetTitle("Max fit function [ADC counts]");

  corrWidthGsigFit->Draw("AP");
  corrWidthGsigFit->GetXaxis()->SetTitle("Landau width [ADC counts]");
  corrWidthGsigFit->GetYaxis()->SetTitle("G sigma fit function [ADC counts]");

  corrNoiseGrGsigFit->Draw("AP");
  corrNoiseGrGsigFit->GetXaxis()->SetTitle("G sigma fit function [ADC counts]");
  corrNoiseGrGsigFit->GetYaxis()->SetTitle("Noise group (RMS) [ADC counts]");

  normGraph->Draw("AP");
  normGraph->GetXaxis()->SetTitle("Bias [V]");
  normGraph->GetYaxis()->SetTitle("MPV [ADC counts]");

  normMpvGraph->Draw("AP");
  normMpvGraph->GetXaxis()->SetTitle("Bias [V]");
  normMpvGraph->GetYaxis()->SetTitle("CCE");

  delete servCan;

  TLegend* leg;
  outFile->cd();
  // for(unsigned int i = 0; i < mpvBiasVec.size(); ++i) // loop on the graphs
  //   mpvBiasVec.at(i)->Write();
  mpvAllSensors->Write();
/*
  // build test legend
  TLegend* testLeg = new TLegend(0.7, 0.15, 0.9, 0.6);
  testLeg->SetLineColor(kWhite);
  testLeg->SetFillColor(kWhite);
  TH1I* h1 = new TH1I(); h1->SetLineColor(kBlue); h1->SetFillColor(kBlue);
  testLeg->AddEntry(h1, "Epi 100 #mum");
  TH1I* h2 = new TH1I(); h2->SetLineColor(kRed); h2->SetFillColor(kRed);
  testLeg->AddEntry(h2, "MCZ 200 #mum");
  TH1I* h3 = new TH1I(); h3->SetLineColor(kBlack); h3->SetFillColor(kBlack);
  testLeg->AddEntry(h3, "FTH 200 #mum");
  testLeg->AddEntry((TObject*) NULL, "P-stop open symbols", "");
  testLeg->AddEntry((TObject*) NULL, "P-spray full symbols", "");

  // build legend for the fluences
  TLegend* lineLeg = new TLegend(0.3, 0.15, 0.5, 0.6);
  lineLeg->SetFillColor(kWhite);
  lineLeg->SetLineColor(kWhite);
  TGraph* a = new TGraph(); a->SetLineStyle(0);
  TGraph* b = new TGraph(); b->SetLineStyle(9);
  TGraph* c = new TGraph(); c->SetLineStyle(2);
  lineLeg->AddEntry(a, "0 cm^{-2}", "L");
  lineLeg->AddEntry(b, "1e15 cm^{-2}", "L");
  lineLeg->AddEntry(c, "1.5e15 cm^{-2}", "L");
  lineLeg->AddEntry(a, "3e15 cm^{-2}", "L");
  lineLeg->AddEntry(b, "1.3e16 cm^{-2}", "L");
*/
  // next legend attempt
  TLegend* legend = new TLegend(0.77, 0.21, 1, 0.8);
  legend->SetLineColor(kWhite);
  legend->SetFillColor(kWhite);
  TGraph* epiLeg = new TGraph(); epiLeg->SetMarkerStyle(20); epiLeg->SetMarkerSize(2); epiLeg->SetLineWidth(1);
  legend->AddEntry(epiLeg, "Epi 100 #mum", "PL");
  TGraph* mczLeg = new TGraph(); mczLeg->SetMarkerStyle(22); mczLeg->SetMarkerSize(2); mczLeg->SetLineWidth(2);
  legend->AddEntry(mczLeg, "MCz 200 #mum", "PL");
  TGraph* fthLeg = new TGraph(); fthLeg->SetMarkerStyle(21); fthLeg->SetMarkerSize(2); fthLeg->SetLineWidth(2);
  legend->AddEntry(fthLeg, "Fth 200 #mum", "PL");
  legend->AddEntry((TObject*) NULL, "P-stop open symbols", "");
  legend->AddEntry((TObject*) NULL, "P-spray full symbols", "");
  TGraph* f0 = new TGraph(); f0->SetLineColor(kRed); f0->SetLineStyle(2);
  legend->AddEntry(f0, "0 cm^{-2}", "L");
  TGraph* f1e15 = new TGraph(); f1e15->SetLineColor(kGreen);
  legend->AddEntry(f1e15, "1 10^{15} cm^{-2}", "L");
  TGraph* f15e14 = new TGraph(); f15e14->SetLineColor(kBlack);
  legend->AddEntry(f15e14, "1.5 10^{15} cm^{-2}", "L");
  TGraph* f3e15 = new TGraph(); f3e15->SetLineColor(kBlue);
  legend->AddEntry(f3e15, "3 10^{15} cm^{-2}", "L");
  TGraph* f13e15 = new TGraph(); f13e15->SetLineColor(kRed);
  legend->AddEntry(f13e15, "1.3 10^{16} cm^{-2}", "L");

  TCanvas* mpvAllSenCan = new TCanvas("mpvAllSenCan");
  mpvAllSensors->Draw("APL");
  // leg = mpvAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  // legMpv->Draw(); // different test legend
  // testLeg->Draw(); // different test legend
  // lineLeg->Draw(); // different test legend
  legend->Draw();
  mpvAllSenCan->SetGridx();
  mpvAllSenCan->SetGridy();
  mpvAllSenCan->Modified();
  mpvAllSenCan->Update();
  mpvAllSenCan->Write();

  meanAllSensors->Write();

  TCanvas* meanAllSenCan = new TCanvas("meanAllSenCan");
  meanAllSensors->Draw("APL");
  legend->Draw();
  meanAllSenCan->SetGridx();
  meanAllSenCan->SetGridy();
  meanAllSenCan->Modified();
  meanAllSenCan->Update();
  meanAllSenCan->Write();

  medianAllSensors->Write();

  TCanvas* medianAllSenCan = new TCanvas("medianAllSenCan");
  medianAllSensors->Draw("APL");
  legend->Draw();
  medianAllSenCan->SetGridx();
  medianAllSenCan->SetGridy();
  medianAllSenCan->Modified();
  medianAllSenCan->Update();
  medianAllSenCan->Write();

  // for(unsigned int i = 0; i < maxFitBiasVec.size(); ++i) // loop on the graphs
  //   maxFitBiasVec.at(i)->Write();
  maxFitAllSensors->Write();

  TCanvas* maxFitAllSenCan = new TCanvas("maxFitAllSenCan");
  maxFitAllSensors->Draw("APL");
  // leg = maxFitAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  maxFitAllSenCan->SetGridx();
  maxFitAllSenCan->SetGridy();
  maxFitAllSenCan->Modified();
  maxFitAllSenCan->Update();
  maxFitAllSenCan->Write();

  // for(unsigned int i = 0; i < lanWBiasVec.size(); ++i) // loop on the graphs
  //   lanWBiasVec.at(i)->Write();
  lanWAllSensors->Write();

  TCanvas* lanWAllSenCan = new TCanvas("lanWAllSenCan");
  lanWAllSensors->Draw("APL");
  // leg = lanWAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  lanWAllSenCan->SetGridx();
  lanWAllSenCan->SetGridy();
  lanWAllSenCan->Modified();
  lanWAllSenCan->Update();
  lanWAllSenCan->Write();

  // for(unsigned int i = 0; i < gSigBiasVec.size(); ++i) // loop on the graphs
  //   gSigBiasVec.at(i)->Write();
  gSigAllSensors->Write();

  TCanvas* gSigAllSenCan = new TCanvas("gSigAllSenCan");
  gSigAllSensors->Draw("APL");
  // leg = gSigAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  gSigAllSenCan->SetGridx();
  gSigAllSenCan->SetGridy();
  gSigAllSenCan->Modified();
  gSigAllSenCan->Update();
  gSigAllSenCan->Write();

  // for(unsigned int i = 0; i < noiseBiasVec.size(); ++i) // loop on the graphs
  //   noiseBiasVec.at(i)->Write();
  noiseAllSensors->Write();

  TCanvas* noiseAllSenCan = new TCanvas("noiseAllSenCan");
  noiseAllSensors->Draw("APL");
  // leg = noiseAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  noiseAllSenCan->SetGridx();
  noiseAllSenCan->SetGridy();
  noiseAllSenCan->Modified();
  noiseAllSenCan->Update();
  noiseAllSenCan->Write();

  // for(unsigned int i = 0; i < noiseGroupBiasVec.size(); ++i) // loop on the graphs
  //   noiseGroupBiasVec.at(i)->Write();
  noiseGroupAllSensors->Write();

  TCanvas* noiseGroupAllSenCan = new TCanvas("noiseGroupAllSenCan");
  noiseGroupAllSensors->Draw("APL");
  // leg = noiseGroupAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  noiseGroupAllSenCan->SetGridx();
  noiseGroupAllSenCan->SetGridy();
  noiseGroupAllSenCan->Modified();
  noiseGroupAllSenCan->Update();
  noiseGroupAllSenCan->Write();

  // for(unsigned int i = 0; i < noisePairBiasVec.size(); ++i) // loop on the graphs
  //   noisePairBiasVec.at(i)->Write();
  noisePairAllSensors->Write();

  TCanvas* noisePairAllSenCan = new TCanvas("noisePairAllSenCan");
  noisePairAllSensors->Draw("APL");
  // leg = noisePairAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  noisePairAllSenCan->SetGridx();
  noisePairAllSenCan->SetGridy();
  noisePairAllSenCan->Modified();
  noisePairAllSenCan->Update();
  noisePairAllSenCan->Write();

  // for(unsigned int i = 0; i < snrBiasVec.size(); ++i) // loop on the graphs
  //   snrBiasVec.at(i)->Write();
  snrAllSensors->Write();

  TCanvas* snrAllSenCan = new TCanvas("snrAllSenCan");
  snrAllSensors->Draw("APL");
  // leg = snrAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  snrAllSenCan->SetGridx();
  snrAllSenCan->SetGridy();
  snrAllSenCan->Modified();
  snrAllSenCan->Update();
  snrAllSenCan->Write();

  // for(unsigned int i = 0; i < eff95BiasVec.size(); ++i) // loop on the graphs
  //   eff95BiasVec.at(i)->Write();
  eff95AllSensors->Write();

  TCanvas* eff95AllSenCan = new TCanvas("eff95AllSenCan");
  eff95AllSensors->Draw("APL");
  // leg = eff95AllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  eff95AllSenCan->SetGridx();
  eff95AllSenCan->SetGridy();
  eff95AllSenCan->Modified();
  eff95AllSenCan->Update();
  eff95AllSenCan->Write();

  // for(unsigned int i = 0; i < chargeSharingBiasVec.size(); ++i) // loop on the graphs
  //   chargeSharingBiasVec.at(i)->Write();
  chargeSharingAllSensors->Write();

  TCanvas* chargeSharingAllSenCan = new TCanvas("chargeSharingAllSenCan");
  chargeSharingAllSensors->Draw("APL");
  // leg = chargeSharingAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  chargeSharingAllSenCan->SetGridx();
  chargeSharingAllSenCan->SetGridy();
  chargeSharingAllSenCan->Modified();
  chargeSharingAllSenCan->Update();
  chargeSharingAllSenCan->Write();

  // for(unsigned int i = 0; i < resYBiasVec.size(); ++i) // loop on the graphs
  //   resYBiasVec.at(i)->Write();
  resYAllSensors->Write();

  TCanvas* resYAllSenCan = new TCanvas("resYAllSenCan");
  resYAllSensors->Draw("APL");
  // leg = resYAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  resYAllSenCan->SetGridx();
  resYAllSenCan->SetGridy();
  resYAllSenCan->Modified();
  resYAllSenCan->Update();
  resYAllSenCan->Write();

  // for(unsigned int i = 0; i < chipTempBiasVec.size(); ++i) // loop on the graphs
  //   chipTempBiasVec.at(i)->Write();
  chipTempAllSensors->Write();

  TCanvas* chipTempAllSenCan = new TCanvas("chipTempAllSenCan");
  chipTempAllSensors->Draw("APL");
  // leg = chipTempAllSenCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  chipTempAllSenCan->SetGridx();
  chipTempAllSenCan->SetGridy();
  chipTempAllSenCan->Modified();
  chipTempAllSenCan->Update();
  chipTempAllSenCan->Write();

  for(unsigned int i = 0; i < histSupVec.size(); ++i)
    {
      leg = histSupVec.at(i)->BuildLegend();
      leg->SetFillColor(kWhite);
      histSupVec.at(i)->Write();
    }

  for(unsigned int i = 0; i < etaSupVec.size(); ++i)
    {
      leg = etaSupVec.at(i)->BuildLegend();
      leg->SetFillColor(kWhite);
      etaSupVec.at(i)->Write();
    }

  TCanvas* chargePosCan;
  for(unsigned int i = 0; i < chargePosSupVec.size(); ++i)
    {
      sprintf(name, "%s_can", chargePosSupVec.at(i)->GetName());
      chargePosCan = new TCanvas(name);
      chargePosCan->SetGridx();
      chargePosCan->SetGridy();
      chargePosSupVec.at(i)->Draw("APL");
      leg = chargePosCan->BuildLegend();
      leg->SetFillColor(kWhite);
      chargePosCan->Write();
      //chargePosSupVec.at(i)->Write();
    }

  TCanvas* corrMPVmaxFitAllSenCan = new TCanvas("corrMPVmaxFitAllSenCan");
  corrMPVmaxFit->Draw("AP");
  leg = corrMPVmaxFitAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  corrMPVmaxFitAllSenCan->SetGridx();
  corrMPVmaxFitAllSenCan->SetGridy();
  corrMPVmaxFitAllSenCan->Modified();
  corrMPVmaxFitAllSenCan->Update();
  corrMPVmaxFitAllSenCan->Write();

  TCanvas* corrWidthGsigFitAllSenCan = new TCanvas("corrWidthGsigFitAllSenCan");
  corrWidthGsigFit->Draw("AP");
  leg = corrWidthGsigFitAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  corrWidthGsigFitAllSenCan->SetGridx();
  corrWidthGsigFitAllSenCan->SetGridy();
  corrWidthGsigFitAllSenCan->Modified();
  corrWidthGsigFitAllSenCan->Update();
  corrWidthGsigFitAllSenCan->Write();

  TCanvas* corrNoiseGrGsigFitAllSenCan = new TCanvas("corrNoiseGrGsigFitAllSenCan");
  corrNoiseGrGsigFit->Draw("AP");
  leg = corrNoiseGrGsigFitAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  corrNoiseGrGsigFitAllSenCan->SetGridx();
  corrNoiseGrGsigFitAllSenCan->SetGridy();
  corrNoiseGrGsigFitAllSenCan->Modified();
  corrNoiseGrGsigFitAllSenCan->Update();
  corrNoiseGrGsigFitAllSenCan->Write();

  normGraph->Write();
  normMpvGraph->Write();

  TCanvas* normMpvGraphCan = new TCanvas("normMpvGraphCan");
  normMpvGraph->Draw("APL");
  // leg = normMpvGraphCan->BuildLegend();
  // leg->SetFillColor(kWhite);
  legend->Draw();
  normMpvGraphCan->SetGridx();
  normMpvGraphCan->SetGridy();
  normMpvGraphCan->Modified();
  normMpvGraphCan->Update();
  normMpvGraphCan->Write();

  mpvAll_electrons->Write();

  TCanvas* mpvAllSenCan_electrons = new TCanvas("mpvAllSenCan_electrons");
  mpvAll_electrons->Draw("APL");
  legend->Draw();
  mpvAllSenCan_electrons->SetGridx();
  mpvAllSenCan_electrons->SetGridy();
  mpvAllSenCan_electrons->Modified();
  mpvAllSenCan_electrons->Update();
  mpvAllSenCan_electrons->Write();

  noiseAll_electrons->Write();

  TCanvas* noiseAllSenCan_electrons = new TCanvas("noiseAllSenCan_electrons");
  noiseAll_electrons->Draw("APL");
  legend->Draw();
  noiseAllSenCan_electrons->SetGridx();
  noiseAllSenCan_electrons->SetGridy();
  noiseAllSenCan_electrons->Modified();
  noiseAllSenCan_electrons->Update();
  noiseAllSenCan_electrons->Write();

  outFile->Close();

  return 0;
}
