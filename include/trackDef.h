#ifndef TRACKDEF_H
#define TRACKDEF_H

struct track
{
  static const int nTelPlanes = 6; // total number of the telescope planes
  static const int nParameters = 5; // parameters of the track model

  double trackParameters[nParameters]; // track parameters
  double measPosTelescope[nTelPlanes][3]; // measured hits on the telescope
  double extraPosTelescope[nTelPlanes][3]; // extrapolated hits on the telescope
  double measPosDUT[3]; // measured hit in the DUT
  double extraPosDUT[3]; // extrapolated hit on the DUT
  bool recHit[nTelPlanes + 1]; // true if a hit was found in the plane, +1 for DUT
  int ndf; // degrees of freedom of the track fit
  double chi2; // chi2 of the track fit

  track() // default constructor
  {
    for(int i = 0; i < nParameters; ++i) trackParameters[i] = 0;

    for(int i = 0; i < nTelPlanes; ++i)
      for(int j = 0; j < 3; ++j)
	{
	  measPosTelescope[i][j] = -10000;
	  extraPosTelescope[i][j] = -10000;
	}

    for(int i = 0; i < 3; ++i)
      {
	measPosDUT[i] = -10000;
	extraPosDUT[i] = -10000;
      }

    for(int i = 0; i < nTelPlanes + 1; ++i) recHit[i] = false;

    ndf = 0;
    chi2 = 0;
  };

};

#endif // #ifndef TRACKDEF_H

