#include "iostream"
#include "fstream"
#include "vector"
#include "math.h"

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

  std::vector<TGraphErrors*> mpvBiasVec;
  std::vector<TGraphErrors*> noiseBiasVec;
  std::vector<TGraphErrors*> resYBiasVec;
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
  TH1* resYDistr;
  TDirectory* resDir;
  TGraphErrors* mpvGr;
  TGraphErrors* noiseGr;
  TGraphErrors* resYGr;

  char name[200];
  char title[500];
  int linStyle = 0;
  int iColor = 0; // color of the graphs

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

      sprintf(name, "mpv_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      mpvGr = new TGraphErrors();
      mpvGr->SetName(name);
      mpvGr->SetTitle(title);
      mpvGr->SetMarkerStyle(8);
      mpvGr->SetFillColor(kWhite);
      iColor = i % 9 + 1;
      if(iColor == 5) ++iColor; // skip yellow
      mpvGr->SetLineColor(iColor); // set line color and style
      mpvGr->SetMarkerColor(iColor);
      if(iColor - 1 == 0) linStyle++;
      mpvGr->SetLineStyle(linStyle);
      mpvGr->SetLineWidth(2);

      sprintf(name, "noise_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      noiseGr = new TGraphErrors();
      noiseGr->SetName(name);
      noiseGr->SetTitle(title);
      noiseGr->SetMarkerStyle(8);
      noiseGr->SetFillColor(kWhite);
      noiseGr->SetLineColor(iColor); // set line color and style
      noiseGr->SetMarkerColor(iColor);
      noiseGr->SetLineStyle(linStyle);
      noiseGr->SetLineWidth(2);

      sprintf(name, "resY_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      resYGr = new TGraphErrors();
      resYGr->SetName(name);
      resYGr->SetTitle(title);
      resYGr->SetMarkerStyle(8);
      resYGr->SetFillColor(kWhite);
      resYGr->SetLineColor(iColor); // set line color and style
      resYGr->SetMarkerColor(iColor);
      resYGr->SetLineStyle(linStyle);
      resYGr->SetLineWidth(2);

      for(unsigned int iRun = 0; iRun < bias.size(); ++iRun) // loop on the runs for a sensor
	{
	  sprintf(name, "%s/%d.root", argv[2], run.at(iRun));
	  inFile = TFile::Open(name);

	  chDistr = (TH1*) inFile->Get("signalDistrTimeCut");
	  func = chDistr->GetFunction("lanGausFit_0");

	  mpvGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(1));
	  mpvGr->SetPointError(iRun, 0, func->GetParError(1));

	  noiseDistr =  (TH1*) inFile->Get("fittedNoiseDistr");

	  noiseGr->SetPoint(iRun, fabs(bias.at(iRun)), noiseDistr->GetMean());
	  noiseGr->SetPointError(iRun, 0, noiseDistr->GetRMS());

	  resDir = (TDirectory*) inFile->Get("Residuals");
	  resYDistr = (TH1*) resDir->Get("residualsY");
	  func = resYDistr->GetFunction("fitFunc");

	  resYGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(2));
	  resYGr->SetPointError(iRun, 0, func->GetParError(2));

	  inFile->Close();
	} // loop on the runs for a sensor


      mpvBiasVec.push_back(mpvGr);
      noiseBiasVec.push_back(noiseGr);
      resYBiasVec.push_back(resYGr);

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
      noiseGr = noiseBiasVec.at(i);

      sprintf(name, "snr_%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      snrGr = new TGraphErrors();
      snrGr->SetName(name);
      snrGr->SetTitle(title);
      snrGr->SetMarkerStyle(8);
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

  TMultiGraph* noiseAllSensors = new TMultiGraph();
  noiseAllSensors->SetName("noiseAllSensors");
  noiseAllSensors->SetTitle("Noise vs bias");

  TMultiGraph* snrAllSensors = new TMultiGraph();
  snrAllSensors->SetName("snrAllSensors");
  snrAllSensors->SetTitle("SNR vs bias");

  TMultiGraph* resYAllSensors = new TMultiGraph();
  resYAllSensors->SetName("resYAllSensors");
  resYAllSensors->SetTitle("#sigma resduals Y vs bias");

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

  outFile->Close();

  return 0;
}
