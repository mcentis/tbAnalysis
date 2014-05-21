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

TF1* lanGausFit(TH1* inHist, double fitR1, double fitR2, double gausSig = 1) // function to be used
{
   // Setting fit range and start values
   Double_t fr[2];
   Double_t sv[4], pllo[4], plhi[4];
   fr[0]=fitR1;
   fr[1]=fitR2;

   // find mpv and integral start value
   int binMin; // max and min bin number (corresponding to the range)
   int binMax;
   double intStart = 0; // start value of the integral
   double mpvStart = 0; // start value of the mpv

   double binW = inHist->GetXaxis()->GetBinWidth(5); // bin width from a random bin
   double xMin = inHist->GetXaxis()->GetXmin();
   binMin = 1 + (fitR1 - xMin) / binW;
   binMax = 1 + (fitR2 - xMin) / binW;

   double binCont;
   double yMax = 0; // variable used to look for the maximum (mpv start)
   for(int iBn = binMin; iBn < binMax; ++iBn)
     {
       binCont = inHist->GetBinContent(iBn);
       intStart += binCont;
       if(binCont > yMax) 
	 {
	   yMax = binCont;
	   mpvStart = inHist->GetXaxis()->GetBinCenter(iBn);
	 }
     }

   // starting parameters
   sv[0] = 5;//landau width
   sv[1] = mpvStart; // mpv landau
   sv[2] = intStart; // integral
   sv[3] = gausSig + gausSig * 0.5; // gaussian width

   // std::cout << "Fitting histogram " << inHist->GetName() << std::endl;
   // std::cout << "Starting parameters" << std::endl;
   // std::cout << "Landau width " << sv[0] << std::endl;
   // std::cout << "MPV          " << sv[1] << std::endl;
   // std::cout << "Area         " << sv[2] << std::endl;
   // std::cout << "Gaus sigma   " << sv[0] << std::endl;

   // parameter limits
   pllo[0]=0.01; pllo[1]=0.5; pllo[2]=1.0; pllo[3]=gausSig;
   plhi[0]=20.0; plhi[1]=200.0; plhi[2]=10000000.0; plhi[3]=20.0;

   static int funcCall = 0;
   char name[50];
   sprintf(name, "lanGausFit_%i", funcCall); // different names in order not to overwrite oder fit functions
   funcCall++;

   TF1* ffit = new TF1(name, langaufun, fr[0], fr[1], 4);
   ffit->SetParameters(sv);
   ffit->SetParNames("Width","MPV","Area","GSigma");
   
   for (int i = 0; i < 4; i++)
     {
       ffit->SetParLimits(i, pllo[i], plhi[i]);
     }

   inHist->Fit(ffit,"RBQ");   // fit within specified range, use ParLimits

   return ffit;
}


#endif // #ifndef LANGAUFIT_HH
