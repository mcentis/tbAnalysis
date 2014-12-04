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

  std::vector<std::string> sensorLabel; // labels for captions
  for(unsigned int i = 0; i < sensorType.size(); ++i)
    {
      if(sensorType.at(i)[6] == 'P' || sensorType.at(i)[6] == 'p')
	sensorLabel.push_back("p-bulk, p-stop");
      else if(sensorType.at(i)[6] == 'Y' || sensorType.at(i)[6] == 'y')
	sensorLabel.push_back("p-bulk, p-spray");
      else if(sensorType.at(i)[6] == 'N' || sensorType.at(i)[6] == 'n')
	sensorLabel.push_back("n-bulk         ");
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

  std::vector<TGraphErrors*> maxFitBiasVec;
  std::vector<TGraphErrors*> mpvBiasVec;
  std::vector<TGraphErrors*> lanWBiasVec;
  std::vector<TGraphErrors*> gSigBiasVec;
  std::vector<TGraphErrors*> noiseBiasVec;
  std::vector<TGraphErrors*> noiseGroupBiasVec;
  std::vector<TGraphErrors*> noisePairBiasVec;
  std::vector<TGraphErrors*> eff95BiasVec;
  std::vector<TGraphErrors*> chargeSharingBiasVec;
  std::vector<TGraphErrors*> resYBiasVec;
  std::vector<TGraphErrors*> chipTempBiasVec;

  std::vector<TCanvas*> histSupVec; // canvases with the superimposition of the used histos
  std::vector<TCanvas*> etaSupVec; // canvases with the superimposition of the used eta distr

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
  TGraphErrors* maxFitGr;
  TGraphErrors* mpvGr;
  TGraphErrors* lanWGr;
  TGraphErrors* gSigGr;
  TGraphErrors* noiseGr;
  TGraphErrors* noiseGroupGr;
  TGraphErrors* noisePairGr;
  TGraphErrors* eff95Gr;
  TGraphErrors* chargeSharingGr;
  TGraphErrors* resYGr;
  TGraphErrors* chipTempGr;
  TCanvas* histSupCan; // each sensor gets a canvas (charge distribution)
  TCanvas* etaSupCan; // each sensor gets a canvas (eta distribution)

  char name[200];
  char title[500];
  int linStyle = 0;
  int iColor = 0; // color of the graphs
  int mrkStyle = 0; // marker style
  int iColHist = 0; // color for the histograms
  int iColEta = 0; // color for the histograms

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

      // the line represents the material
      if(sensorMaterial.at(i) == "Epi" || sensorMaterial.at(i) == "epi")
	linStyle = 9;
      else if(sensorMaterial.at(i) == "Mcz" || sensorMaterial.at(i) == "mcz")
	linStyle = 2;
      else if(sensorMaterial.at(i) == "Fth" || sensorMaterial.at(i) == "fth")
	linStyle = 10;
      else
	{
	  linStyle = 1;
	  std::cout << "-------------------------------------------------------------------- Material not recognized for " << sensorType.at(i) << std::endl;
	}

      sprintf(name, "mpv_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      mpvGr = new TGraphErrors();
      mpvGr->SetName(name);
      mpvGr->SetTitle(title);
      mpvGr->SetMarkerStyle(mrkStyle);
      mpvGr->SetFillColor(kWhite);
      mpvGr->SetLineColor(iColor); // set line color and style
      mpvGr->SetMarkerColor(iColor);
      mpvGr->SetLineStyle(linStyle);

      sprintf(name, "maxFit_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      maxFitGr = new TGraphErrors();
      maxFitGr->SetName(name);
      maxFitGr->SetTitle(title);
      maxFitGr->SetMarkerStyle(mrkStyle);
      maxFitGr->SetFillColor(kWhite);
      maxFitGr->SetLineColor(iColor); // set line color and style
      maxFitGr->SetMarkerColor(iColor);
      maxFitGr->SetLineStyle(linStyle);

      sprintf(name, "lanW_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      lanWGr = new TGraphErrors();
      lanWGr->SetName(name);
      lanWGr->SetTitle(title);
      lanWGr->SetMarkerStyle(mrkStyle);
      lanWGr->SetFillColor(kWhite);
      lanWGr->SetLineColor(iColor); // set line color and style
      lanWGr->SetMarkerColor(iColor);
      lanWGr->SetLineStyle(linStyle);

      sprintf(name, "gSig_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      gSigGr = new TGraphErrors();
      gSigGr->SetName(name);
      gSigGr->SetTitle(title);
      gSigGr->SetMarkerStyle(mrkStyle);
      gSigGr->SetFillColor(kWhite);
      gSigGr->SetLineColor(iColor); // set line color and style
      gSigGr->SetMarkerColor(iColor);
      gSigGr->SetLineStyle(linStyle);

      sprintf(name, "noise_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      noiseGr = new TGraphErrors();
      noiseGr->SetName(name);
      noiseGr->SetTitle(title);
      noiseGr->SetMarkerStyle(mrkStyle);
      noiseGr->SetFillColor(kWhite);
      noiseGr->SetLineColor(iColor); // set line color and style
      noiseGr->SetMarkerColor(iColor);
      noiseGr->SetLineStyle(linStyle);

      sprintf(name, "noiseGroup_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      noiseGroupGr = new TGraphErrors();
      noiseGroupGr->SetName(name);
      noiseGroupGr->SetTitle(title);
      noiseGroupGr->SetMarkerStyle(mrkStyle);
      noiseGroupGr->SetFillColor(kWhite);
      noiseGroupGr->SetLineColor(iColor); // set line color and style
      noiseGroupGr->SetMarkerColor(iColor);
      noiseGroupGr->SetLineStyle(linStyle);

      sprintf(name, "noisePair_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      noisePairGr = new TGraphErrors();
      noisePairGr->SetName(name);
      noisePairGr->SetTitle(title);
      noisePairGr->SetMarkerStyle(mrkStyle);
      noisePairGr->SetFillColor(kWhite);
      noisePairGr->SetLineColor(iColor); // set line color and style
      noisePairGr->SetMarkerColor(iColor);
      noisePairGr->SetLineStyle(linStyle);

      sprintf(name, "eff95_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      eff95Gr = new TGraphErrors();
      eff95Gr->SetName(name);
      eff95Gr->SetTitle(title);
      eff95Gr->SetMarkerStyle(mrkStyle);
      eff95Gr->SetFillColor(kWhite);
      eff95Gr->SetLineColor(iColor); // set line color and style
      eff95Gr->SetMarkerColor(iColor);
      eff95Gr->SetLineStyle(linStyle);

      sprintf(name, "chargeSharing_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      chargeSharingGr = new TGraphErrors();
      chargeSharingGr->SetName(name);
      chargeSharingGr->SetTitle(title);
      chargeSharingGr->SetMarkerStyle(mrkStyle);
      chargeSharingGr->SetFillColor(kWhite);
      chargeSharingGr->SetLineColor(iColor); // set line color and style
      chargeSharingGr->SetMarkerColor(iColor);
      chargeSharingGr->SetLineStyle(linStyle);

      sprintf(name, "resY_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      resYGr = new TGraphErrors();
      resYGr->SetName(name);
      resYGr->SetTitle(title);
      resYGr->SetMarkerStyle(mrkStyle);
      resYGr->SetFillColor(kWhite);
      resYGr->SetLineColor(iColor); // set line color and style
      resYGr->SetMarkerColor(iColor);
      resYGr->SetLineStyle(linStyle);

      sprintf(name, "chipTemp_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      chipTempGr = new TGraphErrors();
      chipTempGr->SetName(name);
      chipTempGr->SetTitle(title);
      chipTempGr->SetMarkerStyle(mrkStyle);
      chipTempGr->SetFillColor(kWhite);
      chipTempGr->SetLineColor(iColor); // set line color and style
      chipTempGr->SetMarkerColor(iColor);
      chipTempGr->SetLineStyle(linStyle);

      sprintf(name, "histSup_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "Distributions at different biases %s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      histSupCan = new TCanvas(name, title);
      histSupCan->SetGridx();
      histSupCan->SetGridy();

      sprintf(name, "etaSup_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "#eta distributions at different biases %s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      etaSupCan = new TCanvas(name, title);
      etaSupCan->SetGridx();
      etaSupCan->SetGridy();

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

	  etaDistr = (TH1*) inFile->Get("etaDistrTimeCutDistCut"); // from the "clusters"
	  //etaDistr = (TH1*) inFile->Get("etaDistrTrackTimeCut"); // track based
	  etaDistr->Sumw2();
	  etaDistr->Scale(1 / etaDistr->Integral());

	  startBin = etaDistr->GetXaxis()->FindBin(0.2);
	  endBin = etaDistr->GetXaxis()->FindBin(0.8);
	  chargeSharing = etaDistr->IntegralAndError(startBin, endBin, EchargeSharing) / 0.6;
	  EchargeSharing /= 0.6;

	  chargeSharingGr->SetPoint(iRun, fabs(bias.at(iRun)), chargeSharing);
	  chargeSharingGr->SetPointError(iRun, 0, EchargeSharing);

	  iColEta = iRun % 9 + 1;
	  if(iColEta == 5) ++iColEta; // skip yellow
	  etaDistr->SetLineColor(iColEta);
	  sprintf(title, "%.00f", bias.at(iRun));
	  etaDistr->SetTitle(title);
	  //etaDistr->SetLineWidth(2);
	  etaSupCan->cd();
	  if(iRun == 0) etaDistr->Draw("E");
	  else etaDistr->Draw("Esame");
	  // if(iRun == 0) etaDistr->Draw("hist");
	  // else etaDistr->Draw("histsame");

	  resDir = (TDirectory*) inFile->Get("Residuals");
	  resYDistr = (TH1*) resDir->Get("residualsYselected");
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
      chargeSharingBiasVec.push_back(chargeSharingGr);
      resYBiasVec.push_back(resYGr);
      chipTempBiasVec.push_back(chipTempGr);
      histSupVec.push_back(histSupCan);
      etaSupVec.push_back(etaSupCan);

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
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
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
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      corrGr = new TGraphErrors();
      corrGr->SetName(name);
      corrGr->SetTitle(title);
      corrGr->SetMarkerStyle(gr1->GetMarkerStyle());
      corrGr->SetFillColor(kWhite);
      corrGr->SetLineColor(gr1->GetLineColor()); // set line color and style
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
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      corrGr = new TGraphErrors();
      corrGr->SetName(name);
      corrGr->SetTitle(title);
      corrGr->SetMarkerStyle(gr1->GetMarkerStyle());
      corrGr->SetFillColor(kWhite);
      corrGr->SetLineColor(gr1->GetLineColor()); // set line color and style
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
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      corrGr = new TGraphErrors();
      corrGr->SetName(name);
      corrGr->SetTitle(title);
      corrGr->SetMarkerStyle(gr1->GetMarkerStyle());
      corrGr->SetFillColor(kWhite);
      corrGr->SetLineColor(gr1->GetLineColor()); // set line color and style
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
      sprintf(title, "%s %s %s %.01e n_{eq} cm^{-2}", sensorMaterial.at(i).c_str(), sensorThickness.at(i).c_str(), sensorLabel.at(i).c_str(), fluences.at(i));
      normMpvGr = new TGraphErrors();
      normMpvGr->SetName(name);
      normMpvGr->SetTitle(title);
      normMpvGr->SetMarkerStyle(mpvBiasVec.at(i)->GetMarkerStyle());
      normMpvGr->SetFillColor(kWhite);
      normMpvGr->SetLineColor(mpvBiasVec.at(i)->GetLineColor()); // set line color and style
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

  for(unsigned int i = 0; i < chargeSharingBiasVec.size(); ++i) // loop on the graphs
    {
      chargeSharingBiasVec.at(i)->Draw("AP");
      chargeSharingBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      chargeSharingBiasVec.at(i)->GetYaxis()->SetTitle("Charge sharing");

      chargeSharingAllSensors->Add(chargeSharingBiasVec.at(i));
    }

  chargeSharingAllSensors->Draw("AP");
  chargeSharingAllSensors->GetXaxis()->SetTitle("Bias [V]");
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

  // correlations graphs
  corrMPVmaxFit->Draw("AP");
  corrMPVmaxFit->GetXaxis()->SetTitle("Landau MPV [ADC]");
  corrMPVmaxFit->GetYaxis()->SetTitle("Max fit function [ADC]");

  corrWidthGsigFit->Draw("AP");
  corrWidthGsigFit->GetXaxis()->SetTitle("Landau width [ADC]");
  corrWidthGsigFit->GetYaxis()->SetTitle("G sigma fit function [ADC]");

  corrNoiseGrGsigFit->Draw("AP");
  corrNoiseGrGsigFit->GetXaxis()->SetTitle("G sigma fit function [ADC]");
  corrNoiseGrGsigFit->GetYaxis()->SetTitle("Noise group (RMS) [ADC]");

  normGraph->Draw("AP");
  normGraph->GetXaxis()->SetTitle("Bias [V]");
  normGraph->GetYaxis()->SetTitle("MPV [ADC]");

  normMpvGraph->Draw("AP");
  normMpvGraph->GetXaxis()->SetTitle("Bias [V]");
  normMpvGraph->GetYaxis()->SetTitle("CCE");

  delete servCan;

  outFile->cd();
  // for(unsigned int i = 0; i < mpvBiasVec.size(); ++i) // loop on the graphs
  //   mpvBiasVec.at(i)->Write();
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

  // for(unsigned int i = 0; i < maxFitBiasVec.size(); ++i) // loop on the graphs
  //   maxFitBiasVec.at(i)->Write();
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

  // for(unsigned int i = 0; i < lanWBiasVec.size(); ++i) // loop on the graphs
  //   lanWBiasVec.at(i)->Write();
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

  // for(unsigned int i = 0; i < gSigBiasVec.size(); ++i) // loop on the graphs
  //   gSigBiasVec.at(i)->Write();
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

  // for(unsigned int i = 0; i < noiseBiasVec.size(); ++i) // loop on the graphs
  //   noiseBiasVec.at(i)->Write();
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

  // for(unsigned int i = 0; i < noiseGroupBiasVec.size(); ++i) // loop on the graphs
  //   noiseGroupBiasVec.at(i)->Write();
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

  // for(unsigned int i = 0; i < noisePairBiasVec.size(); ++i) // loop on the graphs
  //   noisePairBiasVec.at(i)->Write();
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

  // for(unsigned int i = 0; i < snrBiasVec.size(); ++i) // loop on the graphs
  //   snrBiasVec.at(i)->Write();
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

  // for(unsigned int i = 0; i < eff95BiasVec.size(); ++i) // loop on the graphs
  //   eff95BiasVec.at(i)->Write();
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

  // for(unsigned int i = 0; i < chargeSharingBiasVec.size(); ++i) // loop on the graphs
  //   chargeSharingBiasVec.at(i)->Write();
  chargeSharingAllSensors->Write();

  TCanvas* chargeSharingAllSenCan = new TCanvas("chargeSharingAllSenCan");
  chargeSharingAllSensors->Draw("APL");
  leg = chargeSharingAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
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
  leg = resYAllSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
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

  for(unsigned int i = 0; i < etaSupVec.size(); ++i)
    {
      leg = etaSupVec.at(i)->BuildLegend();
      leg->SetFillColor(kWhite);
      etaSupVec.at(i)->Write();
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
  leg = normMpvGraphCan->BuildLegend();
  leg->SetFillColor(kWhite);
  normMpvGraphCan->SetGridx();
  normMpvGraphCan->SetGridy();
  normMpvGraphCan->Modified();
  normMpvGraphCan->Update();
  normMpvGraphCan->Write();

  outFile->Close();

  return 0;
}
