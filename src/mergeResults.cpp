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
  TGraphErrors* mpvGr;

  char name[200];
  char title[500];
  int linStyle = 0;

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

      sprintf(name, "%s_%.01e", sensorType.at(i).c_str(), fluences.at(i));
      sprintf(title, "%s %.01e n_{eq} cm^{-2}", sensorType.at(i).c_str(), fluences.at(i));
      mpvGr = new TGraphErrors();
      mpvGr->SetName(name);
      mpvGr->SetTitle(title);
      mpvGr->SetMarkerStyle(8);
      mpvGr->SetFillColor(kWhite);
      mpvGr->SetLineColor(i % 9 + 1); // set line color and style
      mpvGr->SetMarkerColor(i % 9 + 1);
      if(i % 9 == 0) linStyle++;
      mpvGr->SetLineStyle(linStyle);
      mpvGr->SetLineWidth(2);

      for(unsigned int iRun = 0; iRun < bias.size(); ++iRun) // loop on the runs for a sensor
	{
	  sprintf(name, "%s/%d.root", argv[2], run.at(iRun));
	  inFile = TFile::Open(name);

	  chDistr = (TH1*) inFile->Get("signalDistrTimeCut");
	  func = chDistr->GetFunction("lanGausFit_0");

	  mpvGr->SetPoint(iRun, fabs(bias.at(iRun)), func->GetParameter(1));
	  mpvGr->SetPointError(iRun, 0, func->GetParError(1));

	  inFile->Close();
	} // loop on the runs for a sensor


      mpvBiasVec.push_back(mpvGr);

    } // loop on the sensors

  TMultiGraph* allSensors = new TMultiGraph();
  allSensors->SetName("allSensors");
  allSensors->SetTitle("MPV vs bias");

  TCanvas* servCan = new TCanvas();
  servCan->SetName("servCan");

  for(unsigned int i = 0; i < mpvBiasVec.size(); ++i) // loop on the graphs
    {
      mpvBiasVec.at(i)->Draw("AP");
      mpvBiasVec.at(i)->GetXaxis()->SetTitle("Bias [V]");
      mpvBiasVec.at(i)->GetYaxis()->SetTitle("Landau MPV [ADC]");

      allSensors->Add(mpvBiasVec.at(i));
    }

  allSensors->Draw("AP");
  allSensors->GetXaxis()->SetTitle("Bias [V]");
  allSensors->GetYaxis()->SetTitle("Landau MPV [ADC]");

  delete servCan;

  outFile->cd();
  for(unsigned int i = 0; i < mpvBiasVec.size(); ++i) // loop on the graphs
    mpvBiasVec.at(i)->Write();
  allSensors->Write();

  TCanvas* allSenCan = new TCanvas("allSenCan");
  allSensors->Draw("APL");
  TLegend* leg = allSenCan->BuildLegend();
  leg->SetFillColor(kWhite);
  allSenCan->Modified();
  allSenCan->Update();
  allSenCan->Write();

  outFile->Close();

  return 0;
}
