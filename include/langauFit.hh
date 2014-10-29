#ifndef LANGAUFIT_HH
#define LANGAUFIT_HH

//-----------------------------------------------------------------------
//
//	Convoluted Landau and Gaussian Fitting Function
//         (using ROOT's Landau and Gauss functions)
//
//  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
//  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
//   Markus Friedl (Markus.Friedl@cern.ch)
//
//  to execute this example, do:
//  root > .x langaus.C
// or
//  root > .x langaus.C++
//
//-----------------------------------------------------------------------

#include "TH1.h"
#include "TF1.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

Double_t langaufun(Double_t* x, Double_t* par) // the landau gaussian convolution
{

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

// the function performs a gaussian fit around the highest bin content of the histo and uses the gaus sigma and center to set the 
// fit range and starting parameters of the landau gaus fit
// the parameters are: histo to be fitted, number of sigma in the negative direction for the fit range and the same in positive direction
TF1* lanGausFit(TH1* inHist, double negSigmaFit, double posSigmaFit) // function to be used
{
  int histMaxBin = inHist->GetMaximumBin();
  double histMax = inHist->GetXaxis()->GetBinCenter(histMaxBin);

  const double minMaxPosition = 0.5; // minimum position accepted for the max
  if(histMax < minMaxPosition) // if the max is too low it probably it is in the noise
    {
      int startBin = 1 + (minMaxPosition - inHist->GetXaxis()->GetXmin()) / inHist->GetXaxis()->GetBinWidth(5);
      int endBin = inHist->GetXaxis()->GetNbins() - 2;
      double yMax = 0;
      for(int iBin = startBin; iBin < endBin; ++iBin)
	if(inHist->GetBinContent(iBin) >= yMax)
	  {
	    yMax = inHist->GetBinContent(iBin);
	    histMax = inHist->GetXaxis()->GetBinCenter(iBin);
	  }
    }

  double halfRange = 20; // guess for the range of the gaus fit

  TF1* gausFit = new TF1("gausFit", "gaus", histMax - halfRange, histMax + halfRange);
  gausFit->SetLineColor(kBlue);

  inHist->Fit(gausFit, "RQL");

  double* Gpar = gausFit->GetParameters();

  double fitR1 = Gpar[1] - Gpar[2] * negSigmaFit; // fit range from the gaus mean and sigma
  double fitR2 = Gpar[1] + Gpar[2] * posSigmaFit; //...
  double gausSig = Gpar[2]; // initial guess of the gaus sigma

  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4];
  fr[0]=fitR1;
  fr[1]=fitR2;

  // find mpv and integral start value
  int binMin; // max and min bin number (corresponding to the range)
  int binMax;
  double intStart = 0; // start value of the integral
  double mpvStart = Gpar[1]; // start value of the mpv

  double binW = inHist->GetXaxis()->GetBinWidth(5); // bin width from a random bin
  double xMin = inHist->GetXaxis()->GetXmin();
  binMin = 1 + (fitR1 - xMin) / binW;
  binMax = 1 + (fitR2 - xMin) / binW;

  double binCont;
  // double yMax = 0; // variable used to look for the maximum (mpv start)
  for(int iBn = binMin; iBn < binMax; ++iBn) // integral in the fit range to get the start value
    {
      binCont = inHist->GetBinContent(iBn);
      intStart += binCont;
      // if(binCont > yMax) 
      // 	 {
      // 	   yMax = binCont;
      // 	   mpvStart = inHist->GetXaxis()->GetBinCenter(iBn);
      // 	 }
    }

  // starting parameters
  sv[0] = 5;//landau width
  sv[1] = mpvStart; // mpv landau
  sv[2] = intStart; // integral
  if(gausSig > sv[0])
    sv[3] = sqrt(gausSig * gausSig - sv[0] * sv[0]); // gaussian width
  else
    sv[3] = gausSig;
    
  // std::cout << "Fitting histogram " << inHist->GetName() << std::endl;
  // std::cout << "Starting parameters" << std::endl;
  // std::cout << "Landau width " << sv[0] << std::endl;
  // std::cout << "MPV          " << sv[1] << std::endl;
  // std::cout << "Area         " << sv[2] << std::endl;
  // std::cout << "Gaus sigma   " << sv[3] << std::endl;

  // parameter limits
  pllo[0]=0.01; pllo[1]=-15.0; pllo[2]=1.0; pllo[3]=gausSig * 0.1;
  plhi[0]=20.0; plhi[1]=200.0; plhi[2]=10000000.0; plhi[3]=gausSig;

  TF1* ffit = new TF1("lanGausFit", langaufun, fr[0], fr[1], 4);
  //  ffit->SetNpx(1e4);
  ffit->SetParameters(sv);
  ffit->SetParNames("Width","MPV","Area","GSigma");
  ffit->SetLineColor(kRed);
   
  for (int i = 0; i < 4; i++)
    {
      ffit->SetParLimits(i, pllo[i], plhi[i]);
    }

  inHist->Fit(ffit,"RQL"); // fit within specified range

  return ffit;
}

Double_t gausLangaufun(Double_t* x, Double_t* par) // a peak at 0 and a landau gaussian convolution
{
  Double_t gausPart = par[0] * TMath::Gaus(*x, par[1], par[5]);
  Double_t langauPart = langaufun(x, &par[2]); // par 0 to 2 belong to the gauss part

  return langauPart + gausPart;
}

// the sigma of the convolution and the one of the noise are constrained to be the same
TF1* gausLanGausFit(TH1* inHist, double negSigmaFit, double posSigmaFit)
{
  TF1* gausFunc = new TF1("gausFunc", "gaus", -30, 5); // referred as g0 in the next comments
  inHist->Fit(gausFunc, "RQL");

  TF1* langauFunc = lanGausFit(inHist, negSigmaFit, posSigmaFit);

  const int nPars = 6;
  double par[nPars] = {0};

  for(int i = 0; i < 2; ++i) par[i] = gausFunc->GetParameter(i); // get the gaus fit parameters
  for(int i = 2; i < nPars; ++i) par[i] = langauFunc->GetParameter(i - 2); // get parameters form langaus fit

  double parLimHi[nPars] = {0};
  double parLimLo[nPars] = {0};

  parLimLo[0] = -5; // g0 const
  parLimHi[0] = 1000000;
  parLimLo[1] = par[1] - 0.5 * gausFunc->GetParameter(2); // g0 mean
  parLimHi[1] = par[1] + 0.5 * gausFunc->GetParameter(2);

  for(int i = 2; i < nPars; ++i) // allow a 50% variation on the already fitted parameters
    {
      parLimLo[i] = par[i] - 0.5 * fabs(par[i]);
      parLimHi[i] = par[i] + 0.5 * fabs(par[i]);
    }

  const char* parNames[nPars] = {"ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma"};
  std::cout << "Start parameters and limits\n";
  for(int i = 0; i < nPars; ++i)
    std::cout << parNames[i] << "\t\t" << par[i] << "    " << parLimLo[i] << "   " << parLimHi[i] << " \n";
  std::cout << std::endl;

  double fitR1 = inHist->GetXaxis()->GetXmin();
  double fitR2 = inHist->GetXaxis()->GetXmax();

  TF1* gausLang = new TF1("gausLang", gausLangaufun, fitR1, fitR2, nPars);
  gausLang->SetNpx(1e4);
  gausLang->SetParameters(par);
  gausLang->SetParNames("ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma");
  for(int i = 0; i < nPars; ++i)
    gausLang->SetParLimits(i, parLimLo[i], parLimHi[i]);

  inHist->Fit(gausLang, "RL");

  return gausLang;
}

// the sigma of the convolution and the one of the noise are constrained to be the same, mean and sigma of the noise are fixed
TF1* gausLanGausFitFixGaus(TH1* inHist, double negSigmaFit, double posSigmaFit, double mean, double sigma) // gauss parameters (mean and sigma) from another histo
{
  TF1* langauFunc = lanGausFit(inHist, negSigmaFit, posSigmaFit);

  const int nPars = 6;
  double par[nPars] = {0};

  // set the gaus fit parameters
  par[0] = inHist->GetBinContent(inHist->FindBin(0)); // constant gets the value of the bin at 0
  par[1] = mean; // these 2 remain fixed
  par[5] = sigma;
  for(int i = 2; i < nPars - 1; ++i) par[i] = langauFunc->GetParameter(i - 2); // get parameters form langaus fit, except gaus sigma

  double parLimHi[nPars] = {0};
  double parLimLo[nPars] = {0};

  parLimLo[0] = -5; // g0 const
  parLimHi[0] = 1000000;
  parLimLo[1] = mean; // g0 mean
  parLimHi[1] = mean;
  parLimLo[5] = sigma; // g0 sigma
  parLimHi[5] = sigma;

  for(int i = 2; i < nPars - 1; ++i) // allow a 50% variation on the already fitted parameters, fix gaus sigma
    {
      parLimLo[i] = par[i] - 0.5 * fabs(par[i]);
      parLimHi[i] = par[i] + 0.5 * fabs(par[i]);
    }

  const char* parNames[nPars] = {"ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma"};
  std::cout << "Start parameters and limits\n";
  for(int i = 0; i < nPars; ++i)
    std::cout << parNames[i] << "\t\t" << par[i] << "    " << parLimLo[i] << "   " << parLimHi[i] << " \n";
  std::cout << std::endl;

  double fitR1 = inHist->GetXaxis()->GetXmin();
  double fitR2 = inHist->GetXaxis()->GetXmax();

  TF1* gausLang = new TF1("gausLang", gausLangaufun, fitR1, fitR2, nPars);
  gausLang->SetNpx(1e4);
  gausLang->SetParameters(par);
  gausLang->SetParNames("ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma");
  for(int i = 0; i < nPars; ++i)
    gausLang->SetParLimits(i, parLimLo[i], parLimHi[i]);

  inHist->Fit(gausLang, "RL");

  return gausLang;
}

Double_t gausNoiseLangaufun(Double_t* x, Double_t* par) // a peak at 0 and a landau gaussian convolution
{
  Double_t gausPart = par[0] * TMath::Gaus(*x, par[1], par[2]);
  Double_t langauPart = langaufun(x, &par[3]); // par 0 to 2 belong to the gauss part

  return langauPart + gausPart;
}

// just fix the gaus parameters of the gaussian close to 0, the rest is free
TF1* gausLanGausFitFixGausNoise(TH1* inHist, double negSigmaFit, double posSigmaFit, double mean, double sigma) // gauss parameters (mean and sigma) from another histo
{
  // subtract the noise from the starting histo to determine the landau gauss parameters
  // TF1* noise = new TF1("noise", "gaus", inHist->GetXaxis()->GetXmin(), inHist->GetXaxis()->GetXmax());
  // noise->SetParameter(0, inHist->GetBinContent(inHist->FindBin(0)));
  // noise->SetParameter(1, mean);
  // noise->SetParameter(2, sigma);
  // TH1D* subHist = new TH1D(*inHist);
  // subHist->Add(noise , -1);

  // TF1* langauFunc = lanGausFit(subHist, negSigmaFit, posSigmaFit);
  TF1* langauFunc = lanGausFit(inHist, negSigmaFit, posSigmaFit);

  const int nPars = 7;
  double par[nPars] = {0};

  // set the gaus fit parameters
  par[0] = inHist->GetBinContent(inHist->FindBin(0)); // constant gets the value of the bin at 0
  par[1] = mean; // these 2 remain fixed
  par[2] = sigma;
  for(int i = 3; i < nPars; ++i) par[i] = langauFunc->GetParameter(i - 3); // get parameters form langaus fit, except gaus sigma

  double parLimHi[nPars] = {0};
  double parLimLo[nPars] = {0};

  parLimLo[0] = -5; // g0 const
  parLimHi[0] = 1000000;
  parLimLo[1] = mean; // g0 mean
  parLimHi[1] = mean;
  parLimLo[2] = sigma; // g0 sigma
  parLimHi[2] = sigma;

  for(int i = 3; i < nPars; ++i) // allow a 50% variation on the already fitted parameters, fix gaus sigma
    {
      parLimLo[i] = par[i] - 0.5 * fabs(par[i]);
      parLimHi[i] = par[i] + 0.5 * fabs(par[i]);
    }

  const char* parNames[nPars] = {"ConstG0", "MeanG0", "SigmaG0", "Width", "MPV", "Area", "GSigma"};
  std::cout << "Start parameters and limits\n";
  for(int i = 0; i < nPars; ++i)
    std::cout << parNames[i] << "\t\t" << par[i] << "    " << parLimLo[i] << "   " << parLimHi[i] << " \n";
  std::cout << std::endl;

  double fitR1 = inHist->GetXaxis()->GetXmin();
  double fitR2 = inHist->GetXaxis()->GetXmax();

  TF1* gausLang = new TF1("gausLang", gausNoiseLangaufun, fitR1, fitR2, nPars);
  gausLang->SetNpx(1e4);
  gausLang->SetParameters(par);
  gausLang->SetParNames("ConstG0", "MeanG0", "SigmaG0", "Width", "MPV", "Area", "GSigma");
  for(int i = 0; i < nPars; ++i)
    gausLang->SetParLimits(i, parLimLo[i], parLimHi[i]);

  inHist->Fit(gausLang, "RL");

  return gausLang;
}

#endif // #ifndef LANGAUFIT_HH
