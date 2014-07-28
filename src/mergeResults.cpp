#include "iostream"
#include "fstream"
#include "vector"
#include "math.h"

#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

int main(int argc, char* argv[])
{
  if(argc != 3)
    {
      std::cout << "Usage: mergeResults listOfLists dataDir" << std::endl;
      return 1;
    }

  // big axis labels
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");

  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");

  gStyle->SetTitleOffset(0.95, "x");
  gStyle->SetTitleOffset(0.95, "y");

  gStyle->SetMarkerSize(2);
  gStyle->SetLineWidth(1);

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

  std::vector<TGraphErrors*> maxFitBiasVec;
  std::vector<TGraphErrors*> mpvBiasVec;
  std::vector<TGraphErrors*> lanWBiasVec;
  std::vector<TGraphErrors*> gSigBiasVec;
  std::vector<TGraphErrors*> noiseBiasVec;
  std::vector<TGraphErrors*> noiseGroupBiasVec;
  std::vector<TGraphErrors*> noisePairBiasVec;
  std::vector<TGraphErrors*> eff95BiasVec;
  std::vector<TGraphErrors*> resYBiasVec;
  std::vector<TGraphErrors*> chipTempBiasVec;

  std::vector<TCanvas*> histSupVec; // canvases with the superimposition of the used histos

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
  TH1* resYDistr;
  TH1* tempDistr;
  TDirectory* resDir;
  TGraphErrors* maxFitGr;
  TGraphErrors* mpvGr;
  TGraphErrors* lanWGr;
  TGraphErrors* gSigGr;
  TGraphErrors* noiseGr;
  TGraphErrors* noiseGroupGr;
  TGraphErrors* noisePairGr;
  TGraphErrors* eff95Gr;
  TGraphErrors* resYGr;
  TGraphErrors* chipTempGr;
  TCanvas* histSupCan; // each sensor gets a canvas

  char name[200];
  char title[500];
  int linStyle = 0;
  int iColor = 0; // color of the graphs
  int mrkStyle = 0; // marker style
  int iColHist = 0; // color for the histograms

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

      // the symbol represents the sensor type: n, y, p
      if(sensorType.at(i)[6] == 'P' || sensorType.at(i)[6] == 'p')
	mrkStyle = 20;
      else if(sensorType.at(i)[6] == 'Y' || sensorType.at(i)[6] == 'y')
	mrkStyle = 21;
      else if(sensorType.at(i)[6] == 'N' || sensorType.at(i)[6] == 'n')
	mrkStyle = 22;
      else
	{
	  mrkStyle = 27;
	  std::cout << "-------------------------------------------------------------------- Sensor type not recognized for " << sensorType.at(i) << std::endl;
	}

      // the color will represent the irradiation
      if(fluences.at(i) == 0)
	{
	  iColor = kRed;
	  mrkStyle += 4; // open symbols for the non irradiated devices (works if NPY was recognized...)
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

      // the line represents the thickness, if needed
      linStyle = 9;

      sprintf(name, "mpv_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      mpvGr = new TGraphErrors();
      mpvGr->SetName(name);
      mpvGr->SetTitle(title);
      mpvGr->SetMarkerStyle(mrkStyle);
      mpvGr->SetFillColor(kWhite);
      mpvGr->SetLineColor(iColor); // set line color and style
      mpvGr->SetMarkerColor(iColor);
      mpvGr->SetLineStyle(linStyle);

      sprintf(name, "maxFit_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      maxFitGr = new TGraphErrors();
      maxFitGr->SetName(name);
      maxFitGr->SetTitle(title);
      maxFitGr->SetMarkerStyle(mrkStyle);
      maxFitGr->SetFillColor(kWhite);
      maxFitGr->SetLineColor(iColor); // set line color and style
      maxFitGr->SetMarkerColor(iColor);
      maxFitGr->SetLineStyle(linStyle);

      sprintf(name, "lanW_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      lanWGr = new TGraphErrors();
      lanWGr->SetName(name);
      lanWGr->SetTitle(title);
      lanWGr->SetMarkerStyle(mrkStyle);
      lanWGr->SetFillColor(kWhite);
      lanWGr->SetLineColor(iColor); // set line color and style
      lanWGr->SetMarkerColor(iColor);
      lanWGr->SetLineStyle(linStyle);

      sprintf(name, "gSig_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      gSigGr = new TGraphErrors();
      gSigGr->SetName(name);
      gSigGr->SetTitle(title);
      gSigGr->SetMarkerStyle(mrkStyle);
      gSigGr->SetFillColor(kWhite);
      gSigGr->SetLineColor(iColor); // set line color and style
      gSigGr->SetMarkerColor(iColor);
      gSigGr->SetLineStyle(linStyle);

      sprintf(name, "noise_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      noiseGr = new TGraphErrors();
      noiseGr->SetName(name);
      noiseGr->SetTitle(title);
      noiseGr->SetMarkerStyle(mrkStyle);
      noiseGr->SetFillColor(kWhite);
      noiseGr->SetLineColor(iColor); // set line color and style
      noiseGr->SetMarkerColor(iColor);
      noiseGr->SetLineStyle(linStyle);

      sprintf(name, "noiseGroup_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      noiseGroupGr = new TGraphErrors();
      noiseGroupGr->SetName(name);
      noiseGroupGr->SetTitle(title);
      noiseGroupGr->SetMarkerStyle(mrkStyle);
      noiseGroupGr->SetFillColor(kWhite);
      noiseGroupGr->SetLineColor(iColor); // set line color and style
      noiseGroupGr->SetMarkerColor(iColor);
      noiseGroupGr->SetLineStyle(linStyle);

      sprintf(name, "noisePair_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      noisePairGr = new TGraphErrors();
      noisePairGr->SetName(name);
      noisePairGr->SetTitle(title);
      noisePairGr->SetMarkerStyle(mrkStyle);
      noisePairGr->SetFillColor(kWhite);
      noisePairGr->SetLineColor(iColor); // set line color and style
      noisePairGr->SetMarkerColor(iColor);
      noisePairGr->SetLineStyle(linStyle);

      sprintf(name, "eff95_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      eff95Gr = new TGraphErrors();
      eff95Gr->SetName(name);
      eff95Gr->SetTitle(title);
      eff95Gr->SetMarkerStyle(mrkStyle);
      eff95Gr->SetFillColor(kWhite);
      eff95Gr->SetLineColor(iColor); // set line color and style
      eff95Gr->SetMarkerColor(iColor);
      eff95Gr->SetLineStyle(linStyle);

      sprintf(name, "resY_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      resYGr = new TGraphErrors();
      resYGr->SetName(name);
      resYGr->SetTitle(title);
      resYGr->SetMarkerStyle(mrkStyle);
      resYGr->SetFillColor(kWhite);
      resYGr->SetLineColor(iColor); // set line color and style
      resYGr->SetMarkerColor(iColor);
      resYGr->SetLineStyle(linStyle);

      sprintf(name, "chipTemp_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      chipTempGr = new TGraphErrors();
      chipTempGr->SetName(name);
      chipTempGr->SetTitle(title);
      chipTempGr->SetMarkerStyle(mrkStyle);
      chipTempGr->SetFillColor(kWhite);
      chipTempGr->SetLineColor(iColor); // set line color and style
      chipTempGr->SetMarkerColor(iColor);
      chipTempGr->SetLineStyle(linStyle);

      sprintf(name, "histSup_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "Distributions at different biases %s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      histSupCan = new TCanvas(name, title);
      histSupCan->SetGridx();
      histSupCan->SetGridy();

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

	  mpvGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(3));
	  mpvGr->SetPointError(iRun, 0, func->GetParError(3));

	  maxFitGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetMaximumX());

	  lanWGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(2));
	  lanWGr->SetPointError(iRun, 0, func->GetParError(2));

	  gSigGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(5));
	  gSigGr->SetPointError(iRun, 0, func->GetParError(5));

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

	  resDir = (TDirectory*) inFile->Get("Residuals");
	  resYDistr = (TH1*) resDir->Get("residualsY");
	  func = resYDistr->GetFunction("fitFunc");

	  resYGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(2));
	  resYGr->SetPointError(iRun, 0, func->GetParError(2));

	  tempDistr = (TH1*) inFile->Get("tempDistr");
	  chipTempGr->SetPoint(iRun, fabs(bias.at(iRun)), tempDistr->GetMean());
	  chipTempGr->SetPointError(iRun, 0, tempDistr->GetRMS());


	  iColHist = iRun % 9 + 1;
	  if(iColHist == 5) ++iColHist; // skip yellow
	  chDistr->SetLineColor(iColHist);
	  sprintf(title, "%.00f", bias.at(iRun));
	  chDistr->SetTitle(title);
	  //chDistr->SetLineWidth(2);
	  histSupCan->cd();
	  if(iRun == 0) chDistr->Draw("E");
	  else chDistr->Draw("Esame");

	  //inFile->Close(); // do not close the files, otherwise the canvas does not find the histos to draw and save
	} // loop on the runs for a sensor


      mpvBiasVec.push_back(mpvGr);
      maxFitBiasVec.push_back(maxFitGr);
      lanWBiasVec.push_back(lanWGr);
      gSigBiasVec.push_back(gSigGr);
      noiseBiasVec.push_back(noiseGr);
      noiseGroupBiasVec.push_back(noiseGroupGr);
      noisePairBiasVec.push_back(noisePairGr);
      eff95BiasVec.push_back(eff95Gr);
      resYBiasVec.push_back(resYGr);
      chipTempBiasVec.push_back(chipTempGr);
      histSupVec.push_back(histSupCan);

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
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      snrGr = new TGraphErrors();
      snrGr->SetName(name);
      snrGr->SetTitle(title);
      snrGr->SetMarkerStyle(mpvGr->GetMarkerStyle());
      snrGr->SetFillColor(kWhite);
      snrGr->SetLineColor(mpvGr->GetLineColor()); // set line color and style
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

  TMultiGraph* mpvAllSensors = new TMultiGraph();
  mpvAllSensors->SetName("mpvAllSensors");
  mpvAllSensors->SetTitle("MPV vs bias");

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

  TMultiGraph* resYAllSensors = new TMultiGraph();
  resYAllSensors->SetName("resYAllSensors");
  resYAllSensors->SetTitle("#sigma resduals Y vs bias");

  TMultiGraph* chipTempAllSensors = new TMultiGraph();
  chipTempAllSensors->SetName("chipTempAllSensors");
  chipTempAllSensors->SetTitle("Chip temperature");

  TCanvas* servCan = new TCanvas();
  servCan->SetName("servCan");

  for(unsigned int i = 0; i < mpvBiasVec.size(); ++i) // loop on the graphs
    {
      mpvBiasVec.at(i)->Draw("AP");
      mpvBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      mpvBiasVec.at(i)->GetYaxis()->SetTitle("Landau MPV [ADC]");

      mpvAllSensors->Add(mpvBiasVec.at(i));
    }

  mpvAllSensors->Draw("AP");
  mpvAllSensors->GetXaxis()->SetTitle("Bias [V]");
  mpvAllSensors->GetYaxis()->SetTitle("Landau MPV [ADC]");

  for(unsigned int i = 0; i < maxFitBiasVec.size(); ++i) // loop on the graphs
    {
      maxFitBiasVec.at(i)->Draw("AP");
      maxFitBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      maxFitBiasVec.at(i)->GetYaxis()->SetTitle("Max of the fit [ADC]");

      maxFitAllSensors->Add(maxFitBiasVec.at(i));
    }

  maxFitAllSensors->Draw("AP");
  maxFitAllSensors->GetXaxis()->SetTitle("Bias [V]");
  maxFitAllSensors->GetYaxis()->SetTitle("Max of the fit [ADC]");

  for(unsigned int i = 0; i < lanWBiasVec.size(); ++i) // loop on the graphs
    {
      lanWBiasVec.at(i)->Draw("AP");
      lanWBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      lanWBiasVec.at(i)->GetYaxis()->SetTitle("Landau width [ADC]");

      lanWAllSensors->Add(lanWBiasVec.at(i));
    }

  lanWAllSensors->Draw("AP");
  lanWAllSensors->GetXaxis()->SetTitle("Bias [V]");
  lanWAllSensors->GetYaxis()->SetTitle("Landau width [ADC]");

  for(unsigned int i = 0; i < gSigBiasVec.size(); ++i) // loop on the graphs
    {
      gSigBiasVec.at(i)->Draw("AP");
      gSigBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      gSigBiasVec.at(i)->GetYaxis()->SetTitle("Gaus #sigma [ADC]");

      gSigAllSensors->Add(gSigBiasVec.at(i));
    }

  gSigAllSensors->Draw("AP");
  gSigAllSensors->GetXaxis()->SetTitle("Bias [V]");
  gSigAllSensors->GetYaxis()->SetTitle("Gaus #sigma [ADC]");

  for(unsigned int i = 0; i < noiseBiasVec.size(); ++i) // loop on the graphs
    {
      noiseBiasVec.at(i)->Draw("AP");
      noiseBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      noiseBiasVec.at(i)->GetYaxis()->SetTitle("Noise [ADC]");

      noiseAllSensors->Add(noiseBiasVec.at(i));
    }

  noiseAllSensors->Draw("AP");
  noiseAllSensors->GetXaxis()->SetTitle("Bias [V]");
  noiseAllSensors->GetYaxis()->SetTitle("Noise [ADC]");

  for(unsigned int i = 0; i < noiseGroupBiasVec.size(); ++i) // loop on the graphs
    {
      noiseGroupBiasVec.at(i)->Draw("AP");
      noiseGroupBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      noiseGroupBiasVec.at(i)->GetYaxis()->SetTitle("Noise (RMS) [ADC]");

      noiseGroupAllSensors->Add(noiseGroupBiasVec.at(i));
    }

  noiseGroupAllSensors->Draw("AP");
  noiseGroupAllSensors->GetXaxis()->SetTitle("Bias [V]");
  noiseGroupAllSensors->GetYaxis()->SetTitle("Noise (RMS) [ADC]");

  for(unsigned int i = 0; i < noisePairBiasVec.size(); ++i) // loop on the graphs
    {
      noisePairBiasVec.at(i)->Draw("AP");
      noisePairBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      noisePairBiasVec.at(i)->GetYaxis()->SetTitle("Noise (RMS) [ADC]");

      noisePairAllSensors->Add(noisePairBiasVec.at(i));
    }

  noisePairAllSensors->Draw("AP");
  noisePairAllSensors->GetXaxis()->SetTitle("Bias [V]");
  noisePairAllSensors->GetYaxis()->SetTitle("Noise (RMS) [ADC]");

  for(unsigned int i = 0; i < snrBiasVec.size(); ++i) // loop on the graphs
    {
      snrBiasVec.at(i)->Draw("AP");
      snrBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      snrBiasVec.at(i)->GetYaxis()->SetTitle("SNR");

      snrAllSensors->Add(snrBiasVec.at(i));
    }

  snrAllSensors->Draw("AP");
  snrAllSensors->GetXaxis()->SetTitle("Bias [V]");
  snrAllSensors->GetYaxis()->SetTitle("SNR");

  for(unsigned int i = 0; i < eff95BiasVec.size(); ++i) // loop on the graphs
    {
      eff95BiasVec.at(i)->Draw("AP");
      eff95BiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      eff95BiasVec.at(i)->GetYaxis()->SetTitle("Threshold 95% eff [ADC]");

      eff95AllSensors->Add(eff95BiasVec.at(i));
    }

  eff95AllSensors->Draw("AP");
  eff95AllSensors->GetXaxis()->SetTitle("Bias [V]");
  eff95AllSensors->GetYaxis()->SetTitle("Threshold 95% eff [ADC]");

  for(unsigned int i = 0; i < resYBiasVec.size(); ++i) // loop on the graphs
    {
      resYBiasVec.at(i)->Draw("AP");
      resYBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      resYBiasVec.at(i)->GetYaxis()->SetTitle("#sigma res Y [mm]");

      resYAllSensors->Add(resYBiasVec.at(i));
    }

  resYAllSensors->Draw("AP");
  resYAllSensors->GetXaxis()->SetTitle("Bias [V]");
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
  chipTempAllSensors->GetYaxis()->SetTitle("Temperature [C]");

  delete servCan;

  outFile->cd();
  for(unsigned int i = 0; i < mpvBiasVec.size(); ++i) // loop on the graphs
    mpvBiasVec.at(i)->Write();
  mpvAllSensors->Write();

  TCanvas* mpvAllSenCan = new TCanvas("mpvAllSenCan");
  mpvAllSensors->Draw("APL");
  TLegend* leg = mpvAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  mpvAllSenCan->SetGridx();
  mpvAllSenCan->SetGridy();
  mpvAllSenCan->Modified();
  mpvAllSenCan->Update();
  mpvAllSenCan->Write();

  for(unsigned int i = 0; i < maxFitBiasVec.size(); ++i) // loop on the graphs
    maxFitBiasVec.at(i)->Write();
  maxFitAllSensors->Write();

  TCanvas* maxFitAllSenCan = new TCanvas("maxFitAllSenCan");
  maxFitAllSensors->Draw("APL");
  leg = maxFitAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  maxFitAllSenCan->SetGridx();
  maxFitAllSenCan->SetGridy();
  maxFitAllSenCan->Modified();
  maxFitAllSenCan->Update();
  maxFitAllSenCan->Write();

  for(unsigned int i = 0; i < lanWBiasVec.size(); ++i) // loop on the graphs
    lanWBiasVec.at(i)->Write();
  lanWAllSensors->Write();

  TCanvas* lanWAllSenCan = new TCanvas("lanWAllSenCan");
  lanWAllSensors->Draw("APL");
  leg = lanWAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  lanWAllSenCan->SetGridx();
  lanWAllSenCan->SetGridy();
  lanWAllSenCan->Modified();
  lanWAllSenCan->Update();
  lanWAllSenCan->Write();

  for(unsigned int i = 0; i < gSigBiasVec.size(); ++i) // loop on the graphs
    gSigBiasVec.at(i)->Write();
  gSigAllSensors->Write();

  TCanvas* gSigAllSenCan = new TCanvas("gSigAllSenCan");
  gSigAllSensors->Draw("APL");
  leg = gSigAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  gSigAllSenCan->SetGridx();
  gSigAllSenCan->SetGridy();
  gSigAllSenCan->Modified();
  gSigAllSenCan->Update();
  gSigAllSenCan->Write();

  for(unsigned int i = 0; i < noiseBiasVec.size(); ++i) // loop on the graphs
    noiseBiasVec.at(i)->Write();
  noiseAllSensors->Write();

  TCanvas* noiseAllSenCan = new TCanvas("noiseAllSenCan");
  noiseAllSensors->Draw("APL");
  leg = noiseAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  noiseAllSenCan->SetGridx();
  noiseAllSenCan->SetGridy();
  noiseAllSenCan->Modified();
  noiseAllSenCan->Update();
  noiseAllSenCan->Write();

  for(unsigned int i = 0; i < noiseGroupBiasVec.size(); ++i) // loop on the graphs
    noiseGroupBiasVec.at(i)->Write();
  noiseGroupAllSensors->Write();

  TCanvas* noiseGroupAllSenCan = new TCanvas("noiseGroupAllSenCan");
  noiseGroupAllSensors->Draw("APL");
  leg = noiseGroupAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  noiseGroupAllSenCan->SetGridx();
  noiseGroupAllSenCan->SetGridy();
  noiseGroupAllSenCan->Modified();
  noiseGroupAllSenCan->Update();
  noiseGroupAllSenCan->Write();

  for(unsigned int i = 0; i < noisePairBiasVec.size(); ++i) // loop on the graphs
    noisePairBiasVec.at(i)->Write();
  noisePairAllSensors->Write();

  TCanvas* noisePairAllSenCan = new TCanvas("noisePairAllSenCan");
  noisePairAllSensors->Draw("APL");
  leg = noisePairAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  noisePairAllSenCan->SetGridx();
  noisePairAllSenCan->SetGridy();
  noisePairAllSenCan->Modified();
  noisePairAllSenCan->Update();
  noisePairAllSenCan->Write();

  for(unsigned int i = 0; i < snrBiasVec.size(); ++i) // loop on the graphs
    snrBiasVec.at(i)->Write();
  snrAllSensors->Write();

  TCanvas* snrAllSenCan = new TCanvas("snrAllSenCan");
  snrAllSensors->Draw("APL");
  leg = snrAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  snrAllSenCan->SetGridx();
  snrAllSenCan->SetGridy();
  snrAllSenCan->Modified();
  snrAllSenCan->Update();
  snrAllSenCan->Write();

  for(unsigned int i = 0; i < eff95BiasVec.size(); ++i) // loop on the graphs
    eff95BiasVec.at(i)->Write();
  eff95AllSensors->Write();

  TCanvas* eff95AllSenCan = new TCanvas("eff95AllSenCan");
  eff95AllSensors->Draw("APL");
  leg = eff95AllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  eff95AllSenCan->SetGridx();
  eff95AllSenCan->SetGridy();
  eff95AllSenCan->Modified();
  eff95AllSenCan->Update();
  eff95AllSenCan->Write();

  for(unsigned int i = 0; i < resYBiasVec.size(); ++i) // loop on the graphs
    resYBiasVec.at(i)->Write();
  resYAllSensors->Write();

  TCanvas* resYAllSenCan = new TCanvas("resYAllSenCan");
  resYAllSensors->Draw("APL");
  leg = resYAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  resYAllSenCan->SetGridx();
  resYAllSenCan->SetGridy();
  resYAllSenCan->Modified();
  resYAllSenCan->Update();
  resYAllSenCan->Write();

  for(unsigned int i = 0; i < chipTempBiasVec.size(); ++i) // loop on the graphs
    chipTempBiasVec.at(i)->Write();
  chipTempAllSensors->Write();

  TCanvas* chipTempAllSenCan = new TCanvas("chipTempAllSenCan");
  chipTempAllSensors->Draw("APL");
  leg = chipTempAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
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

  outFile->Close();

  return 0;
}
