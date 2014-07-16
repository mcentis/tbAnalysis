void backgroundHPH(double sigma, int nChannels)
{
  // range and binning
  const double loR = -20.5;
  const double hiR = 20.5;
  const double binW = 0.1;

  double noise;
  double highestNoise = 0;

  TH1D* backgroundDistrHPH = new TH1D("backgrounDistrHPH", "Backgorund distribution of the strip with highest PH;Signal [ADC];Entries", (hiR - loR) / binW, loR, hiR);

  for(int iEvt = 0; iEvt < 100000; iEvt++)
    {
      highestNoise = -1e3;
      for(int i = 0; i < nChannels; ++i)
	{
	  noise = gRandom->Gaus(0, sigma);
	  if(noise > highestNoise)
	    highestNoise = noise;
	}
      backgroundDistrHPH->Fill(highestNoise);
    }

  TH1D* bgNotNorm = new TH1D(*backgroundDistrHPH); // copy to test the fit function
  bgNotNorm->SetName("bgNotNorm");
  bgNotNorm->SetTitle("Backgorund distribution of the strip with highest PH;Signal [ADC];Entries");

  backgroundDistrHPH->Sumw2();
  backgroundDistrHPH->Scale(1 / backgroundDistrHPH->GetEntries()); // normalize the thing

  TCanvas* numCan = new TCanvas("numCan", "numerical generation");
  backgroundDistrHPH->Draw("hist");

  TF1* backgroundProbHPH = new TF1("backgroundProbHPH", "[1] * TMath::Gaus(x, 0, [0], 1) * TMath::Power(0.5 * (1 + TMath::Erf(x / ([0] * TMath::Sqrt(2)))), [1] - 1)", loR, hiR);
  backgroundProbHPH->SetParameter(0, sigma);
  backgroundProbHPH->SetParameter(1, nChannels);
  backgroundProbHPH->SetNpx(1e5);

  TCanvas* anaCan = new TCanvas("anaCan", "analitical derivation");
  backgroundProbHPH->Draw();

  std::cout << "Integral of the analitical form: " << backgroundProbHPH->Integral(loR, hiR) << std::endl;

  TF1* fit_bg = new TF1("fit_bg", "[2] * [1] * TMath::Gaus(x, 0, [0], 1) * TMath::Power(0.5 * (1 + TMath::Erf(x / ([0] * TMath::Sqrt(2)))), [1] - 1)", loR, hiR);
  fit_bg->SetParameter(0, 1);
  fit_bg->SetParameter(1, nChannels);
  fit_bg->SetParameter(2, bgNotNorm->GetMaximum());

  TCanvas* fitCan = new TCanvas("fitCan", "test fit");
  bgNotNorm->Draw();
  bgNotNorm->Fit(fit_bg);

  return;
}
