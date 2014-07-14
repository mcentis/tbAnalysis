void backgroundHPH(double sigma, int nChannels)
{
  double noise;
  double highestNoise = 0;

  TH1D* backgroundDistrHPH = new TH1D("backgrounDistrHPH", "Backgorund distribution of the strip with highest PH;Signal [ADC];Entries", 201, -100.5, 100.5);

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

  backgroundDistrHPH->Draw();

  return;
}
