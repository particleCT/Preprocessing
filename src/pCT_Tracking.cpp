// Tracker pattern recognition routines
// R.P. Johnson 5/22/2016

#include "pCT_Tracking.h"

// This is the pattern recognition, which works only in 2D, separately for the V-U and T-U views.
std::vector<Tkr2D> pCT_Tracking::Tracking2D(int Idx, TkrHits &pCThits, pCTgeo *Geometry) {

  /* Idx = 0 for V-U view
     Idx = 1 for T-U view
    Algorithm summary:
    Find all reasonable vectors in the front tracker
    Find all reasonable vectors in the rear tracker
    If there are no vectors in just one T tracker
       extrapolate or interpolate the position
       if the position is close to a crack, place the coordinate in the middle of the crack
    If there still are no vectors in just the front tracker
       Use the one hit to make a vector pointing back toward the vertex
    If there still are no vectors in just the rear tracker
       Extrapolate or interpolate to infer the missing hit
    Cut on the displacement at u=0 between front and back vectors
  */

  const double mxSlopeFront[2] = { 0.03, 0.09 }; // Cut on the slope of the front tracker vector, separately for V and T
  const double mxSlopeBack[2] = { 0.1, 0.15 };   // Cut on the slope of the rear tracker vector, separately for V and T
  const double deltaMx = 6.0;                    // Cut on how far the two vectors miss each other at u=0, in mm

  // Create lists of all vectors in the front and back trackers, within the specified slope cuts
  std::vector<Tkr2D> tmp;
  std::vector<vct> front, back;
  for (int i = 0; i < pCThits.Lyr[0].N[Idx]; i++) {
    for (int j = 0; j < pCThits.Lyr[1].N[Idx]; j++) {
      double slope = (pCThits.Lyr[1].Y[Idx].at(j) - pCThits.Lyr[0].Y[Idx].at(i)) /
                     (pCThits.Lyr[1].U[Idx].at(j) - pCThits.Lyr[0].U[Idx].at(i));
      if (abs(slope) < mxSlopeFront[Idx]) {
        vct vctmp;
        vctmp.slope = slope;
        vctmp.Y[0] = pCThits.Lyr[0].Y[Idx].at(i);
        vctmp.U[0] = pCThits.Lyr[0].U[Idx].at(i);
        vctmp.I[0] = i;
        vctmp.Y[1] = pCThits.Lyr[1].Y[Idx].at(j);
        vctmp.U[1] = pCThits.Lyr[1].U[Idx].at(j);
        vctmp.I[1] = j;
        vctmp.intercept = vctmp.Y[0] - slope * vctmp.U[0];
        front.push_back(vctmp);
      }
    }
  }
  for (int i = 0; i < pCThits.Lyr[2].N[Idx]; i++) {
    for (int j = 0; j < pCThits.Lyr[3].N[Idx]; j++) {
      double slope = (pCThits.Lyr[3].Y[Idx].at(j) - pCThits.Lyr[2].Y[Idx].at(i)) /
                     (pCThits.Lyr[3].U[Idx].at(j) - pCThits.Lyr[2].U[Idx].at(i));
      if (abs(slope) < mxSlopeBack[Idx]) {
        vct vctmp;
        vctmp.slope = slope;
        vctmp.Y[0] = pCThits.Lyr[2].Y[Idx].at(i);
        vctmp.U[0] = pCThits.Lyr[2].U[Idx].at(i);
        vctmp.I[0] = i;
        vctmp.Y[1] = pCThits.Lyr[3].Y[Idx].at(j);
        vctmp.U[1] = pCThits.Lyr[3].U[Idx].at(j);
        vctmp.I[1] = j;
        vctmp.intercept = vctmp.Y[0] - slope * vctmp.U[0];
        back.push_back(vctmp);
      }
    }
  }
  if (front.size() == 0 && back.size() == 0)
    return tmp; // Hopeless case; no recovery, so give up.

  // Look at all combinations of front and rear vectors for tracks that match in the center
  for (int i = 0; i < front.size(); i++) {
    for (int j = 0; j < back.size(); j++) {
      double miss = front[i].intercept - back[j].intercept;
      if (abs(miss) < deltaMx) {
        // Check whether any of the hits have already been used. Do not allow hit sharing.
        if (tmp.size() > 0) { // Don't waste time if there are no tracks already (the usual situation)
          bool reject = false;
          std::vector<int> shareList;
          for (int k = 0; k < 2; k++) {
            int tk1 = pCThits.Lyr[k].F[Idx].at(front[i].I[k]);
            if (tk1 >= 0) {
              if (tmp[tk1].Good) {
                if (abs(miss) > abs(tmp[tk1].Miss)) {
                  reject = true;
                  break; // This track is inferior, so don't check further.
                }
                shareList.push_back(tk1);
              }
            }
            tk1 = pCThits.Lyr[k + 2].F[Idx].at(back[j].I[k]);
            if (tk1 >= 0) {
              if (tmp[tk1].Good) {
                if (abs(miss) > abs(tmp[tk1].Miss)) {
                  reject = true;
                  break;
                }
                shareList.push_back(tk1);
              }
            }
          }
          if (reject)
            break; // This track is worse than some tracks already found that shares hits, so skip this one.
          for (int shr = 0; shr < shareList.size(); shr++) {
            tmp[shareList[shr]].Good = false; // Get rid of earlier inferior tracks sharing hits with this one
          }
        }
        Tkr2D tmpTk;
        tmpTk.X[0] = front[i].Y[0];
        tmpTk.U[0] = front[i].U[0];
        tmpTk.Qh[0] = 3;
        pCThits.Lyr[0].F[Idx].at(front[i].I[0]) = tmp.size(); // Flags hit as used and makes a pointer from hit to track
        tmpTk.X[1] = front[i].Y[1];
        tmpTk.U[1] = front[i].U[1];
        tmpTk.Qh[1] = 3;
        pCThits.Lyr[1].F[Idx].at(front[i].I[1]) = tmp.size();
        tmpTk.X[2] = back[j].Y[0];
        tmpTk.U[2] = back[j].U[0];
        tmpTk.Qh[2] = 3;
        pCThits.Lyr[2].F[Idx].at(back[j].I[0]) = tmp.size();
        tmpTk.X[3] = back[j].Y[1];
        tmpTk.U[3] = back[j].U[1];
        tmpTk.Qh[3] = 3;
        pCThits.Lyr[3].F[Idx].at(back[j].I[1]) = tmp.size();
        tmpTk.Miss = miss;
        tmpTk.Q = 3;
        tmpTk.Good = true;
        tmp.push_back(tmpTk);
      }
    }
  }

  // Now consider tracks passing through gaps in the T layers.
  // With unused vectors in the back, try to make vectors in the front with an unused hit and a gap, and look for a
  // match.
  // Then repeat, using gaps in the back and unused vectors in the front.
  if (Idx == 1) { // T layers only, where the cracks are parallel to the strips. This won't work for V layers.
    for (int i = 0; i < back.size(); i++) {
      if (pCThits.Lyr[2].F[Idx].at(back[i].I[0]) < 0 && pCThits.Lyr[3].F[Idx].at(back[i].I[1]) < 0) {
        int otherLyr[2] = { 1, 0 };
        for (int lyr = 0; lyr < 2; lyr++) {
          for (int j = 0; j < pCThits.Lyr[lyr].N[Idx]; j++) {
            int tk = pCThits.Lyr[lyr].F[Idx].at(j);
            if (tk >= 0) {
              if (tmp[tk].Good)
                break;
            }
            double t1 = pCThits.Lyr[lyr].Y[Idx].at(j);
            double u1 = pCThits.Lyr[lyr].U[Idx].at(j);
            for (int gap = 0; gap < 3; gap++) {
              double t2 = Geometry->tBoardGap(gap, otherLyr[lyr]);
              double u2 = Geometry->uT(otherLyr[lyr]);
              double slope = (t2 - t1) / (u2 - u1);
              if (abs(slope) < mxSlopeFront[Idx]) {
                double intercept = t1 - slope * u1;
                double miss = intercept - back[i].intercept;
                if (abs(miss) < deltaMx) {
                  Tkr2D tmpTk;
                  tmpTk.X[lyr] = t1;
                  tmpTk.U[lyr] = u1;
                  tmpTk.Qh[lyr] = 3;
                  pCThits.Lyr[lyr].F[Idx].at(j) = tmp.size(); // Flags hit as used and makes a pointer from hit to track
                  tmpTk.X[otherLyr[lyr]] = t2;
                  tmpTk.U[otherLyr[lyr]] = u2;
                  tmpTk.Qh[otherLyr[lyr]] = 2;
                  tmpTk.X[2] = back[i].Y[0];
                  tmpTk.U[2] = back[i].U[0];
                  tmpTk.Qh[2] = 3;
                  pCThits.Lyr[2].F[Idx].at(back[i].I[0]) = tmp.size();
                  tmpTk.X[3] = back[i].Y[1];
                  tmpTk.U[3] = back[i].U[1];
                  tmpTk.Qh[3] = 3;
                  pCThits.Lyr[3].F[Idx].at(back[i].I[1]) = tmp.size();
                  tmpTk.Miss = miss;
                  tmpTk.Q = 2;
                  tmpTk.Good = true;
                  tmp.push_back(tmpTk);
                  goto nextBackVector;
                }
              }
            }
          }
        }
      }
    nextBackVector:
      ;
    }

    for (int i = 0; i < front.size(); i++) {
      if (pCThits.Lyr[0].F[Idx].at(front[i].I[0]) < 0 && pCThits.Lyr[1].F[Idx].at(front[i].I[1]) < 0) {
        int otherLyr[4] = { 0, 0, 3, 2 };
        for (int lyr = 2; lyr < 4; lyr++) {
          for (int j = 0; j < pCThits.Lyr[lyr].N[Idx]; j++) {
            int tk = pCThits.Lyr[lyr].F[Idx].at(j);
            if (tk >= 0) {
              if (tmp[tk].Good)
                break;
            }
            double t1 = pCThits.Lyr[lyr].Y[Idx].at(j);
            double u1 = pCThits.Lyr[lyr].U[Idx].at(j);
            for (int gap = 0; gap < 3; gap++) {
              double t2 = Geometry->tBoardGap(gap, otherLyr[lyr]);
              double u2 = Geometry->uT(otherLyr[lyr]);
              double slope = (t2 - t1) / (u2 - u1);
              if (abs(slope) < mxSlopeBack[Idx]) {
                double intercept = t1 - slope * u1;
                double miss = intercept - front[i].intercept;
                if (abs(miss) < deltaMx) {
                  Tkr2D tmpTk;
                  tmpTk.X[0] = front[i].Y[0];
                  tmpTk.U[0] = front[i].U[0];
                  tmpTk.Qh[0] = 3;
                  pCThits.Lyr[0].F[Idx].at(front[i].I[0]) = tmp.size();
                  tmpTk.X[1] = front[i].Y[1];
                  tmpTk.U[1] = front[i].U[1];
                  tmpTk.Qh[1] = 3;
                  pCThits.Lyr[1].F[Idx].at(front[i].I[1]) = tmp.size();
                  tmpTk.X[lyr] = t1;
                  tmpTk.U[lyr] = u1;
                  tmpTk.Qh[lyr] = 3;
                  pCThits.Lyr[lyr].F[Idx].at(j) = tmp.size(); // Flags hit as used and makes a pointer from hit to track
                  tmpTk.X[otherLyr[lyr]] = t2;
                  tmpTk.U[otherLyr[lyr]] = u2;
                  tmpTk.Qh[otherLyr[lyr]] = 2;
                  tmpTk.Miss = miss;
                  tmpTk.Q = 2;
                  tmpTk.Good = true;
                  tmp.push_back(tmpTk);
                  goto nextFrontVector;
                }
              }
            }
          }
        }
      }
    nextFrontVector:
      ;
    }
  }

  // Now consider track candidates missing one hit in the front tracker.  Make a vector from the hit that
  // is present that points back to the putative beam origin (this works a bit better at LLUMC, using the
  // lead foil, than at the Chicago Proton Center).
  for (int i = 0; i < back.size(); i++) {
    if (pCThits.Lyr[2].F[Idx].at(back[i].I[0]) < 0 && pCThits.Lyr[3].F[Idx].at(back[i].I[1]) < 0) {
      int otherLyr[2] = { 1, 0 };
      for (int lyr = 0; lyr < 2; lyr++) {
        for (int j = 0; j < pCThits.Lyr[lyr].N[Idx]; j++) {
          int tk = pCThits.Lyr[lyr].F[Idx].at(j);
          if (tk >= 0) {
            if (tmp[tk].Good)
              break;
          }
          double t1 = pCThits.Lyr[lyr].Y[Idx].at(j);
          double u1 = pCThits.Lyr[lyr].U[Idx].at(j);
          double u2;
          if (Idx == 1)
            u2 = Geometry->uT(otherLyr[lyr]);
          else
            u2 = Geometry->uV(otherLyr[lyr]);
          double slope = (t1 - Geometry->BeamVertex(Idx)[Idx]) / (u1 - Geometry->BeamVertex(Idx)[2]);
          double intercept = Geometry->BeamVertex(Idx)[Idx] - slope * Geometry->BeamVertex(Idx)[2];
          double miss = intercept - back[i].intercept;
          if (abs(miss) < deltaMx) {
            double t2 = intercept + slope * u2;
            Tkr2D tmpTk;
            tmpTk.X[lyr] = t1;
            tmpTk.U[lyr] = u1;
            tmpTk.Qh[lyr] = 3;
            pCThits.Lyr[lyr].F[Idx].at(j) = tmp.size(); // Flags hit as used and makes a pointer from hit to track
            tmpTk.X[otherLyr[lyr]] = t2;
            tmpTk.U[otherLyr[lyr]] = u2;
            tmpTk.Qh[otherLyr[lyr]] = 1;
            tmpTk.X[2] = back[i].Y[0];
            tmpTk.U[2] = back[i].U[0];
            tmpTk.Qh[2] = 3;
            pCThits.Lyr[2].F[Idx].at(back[i].I[0]) = tmp.size();
            tmpTk.X[3] = back[i].Y[1];
            tmpTk.U[3] = back[i].U[1];
            tmpTk.Qh[3] = 3;
            pCThits.Lyr[3].F[Idx].at(back[i].I[1]) = tmp.size();
            tmpTk.Q = 1;
            tmpTk.Miss = miss;
            tmpTk.Good = true;
            tmp.push_back(tmpTk);
            goto nextBkVctr;
          }
        }
      }
    }
  nextBkVctr:
    ;
  }

  // Finally, if there is an unused front vector but only one unused hit in the back, try to make a track by
  // extrapolation or interpolation
  // These will be the lowest quality tracks.
  for (int i = 0; i < front.size(); i++) {
    if (pCThits.Lyr[0].F[Idx].at(front[i].I[0]) < 0 && pCThits.Lyr[1].F[Idx].at(front[i].I[1]) < 0) {
      int otherLyr[4] = { 0, 0, 3, 2 };
      for (int lyr = 2; lyr < 4; lyr++) {
        for (int j = 0; j < pCThits.Lyr[lyr].N[Idx]; j++) {
          int tk = pCThits.Lyr[lyr].F[Idx].at(j);
          if (tk >= 0) {
            if (tmp[tk].Good && tmp[tk].Q > 0)
              break; // Skip hits already used on better tracks
          }
          double t1 = front[i].Y[0];
          double u1 = front[i].U[0];
          double t2 = front[i].Y[1];
          double u2 = front[i].U[1];
          double t3 = pCThits.Lyr[lyr].Y[Idx].at(j);
          double u3 = pCThits.Lyr[lyr].U[Idx].at(j);
          double u4;
          if (Idx == 1)
            u4 = Geometry->uT(otherLyr[lyr]);
          else
            u4 = Geometry->uV(otherLyr[lyr]);
          double t4 = quadExtrap(u1, u2, u3, u4, t1, t2, t3);
          double slope = (t4 - t3) / (u4 - u3);
          if (abs(slope) < mxSlopeBack[Idx]) {
            double intercept = t3 - slope * u3;
            double miss = intercept - front[i].intercept;
            if (abs(miss) < deltaMx) {
              if (tk >= 0) {
                if (tmp[tk].Good) {
                  if (abs(tmp[tk].Miss) > abs(miss)) {
                    tmp[tk].Good = false; // This track is better than an existing track using the same hit
                  } else {
                    continue; // An existing track using this hit is better
                  }
                }
              }
              Tkr2D tmpTk;
              tmpTk.X[0] = t1;
              tmpTk.U[0] = u1;
              tmpTk.Qh[0] = 3;
              pCThits.Lyr[0].F[Idx].at(front[i].I[0]) = tmp.size();
              tmpTk.X[1] = t2;
              tmpTk.U[1] = u2;
              tmpTk.Qh[1] = 3;
              pCThits.Lyr[1].F[Idx].at(front[i].I[1]) = tmp.size();
              tmpTk.X[lyr] = t3;
              tmpTk.U[lyr] = u3;
              tmpTk.Qh[lyr] = 3;
              pCThits.Lyr[lyr].F[Idx].at(j) = tmp.size(); // Flags hit as used and makes a pointer from hit to track
              tmpTk.X[otherLyr[lyr]] = t4;
              tmpTk.U[otherLyr[lyr]] = u4;
              tmpTk.Qh[otherLyr[lyr]] = 0;
              tmpTk.Q = 0;
              tmpTk.Miss = miss;
              tmpTk.Good = true;
              tmp.push_back(tmpTk);
            }
          }
        }
      }
    }
  }

  return tmp;
} // End of Tracking2D

pCT_Tracking::pCT_Tracking(TkrHits &pCThits, pCTgeo* Geometry) { // This class constructor finds all the tracks

  nTracks = 0;

  int nVht = 0;
  int nTht = 0;
  for (int lyr = 0; lyr < 4; lyr++) {
    if (pCThits.Lyr[lyr].N[0] > 0)
      nVht++;
    if (pCThits.Lyr[lyr].N[1] > 0)
      nTht++;
  }
  if (nVht < 3 || nTht < 3)
    return; // Abort if more than one layer is missing a hit in either view

  VTracks = Tracking2D(0, pCThits, Geometry);

  if (VTracks.size() == 0)
    return; // Don't waste time with the T view if there is no track in V

  TTracks = Tracking2D(1, pCThits, Geometry);

  if (TTracks.size() == 0)
    return;

  // Count up just the good tracks in V and T.  Bad tracks were rejected because they overlapped (shared hits with)
  // superior tracks.
  int nVtkrs = 0;
  int nTtkrs = 0;
  itkV = -1;
  itkT = -1;
  for (int i = 0; i < VTracks.size(); i++) {
    if (VTracks[i].Good) {
      nVtkrs++;
      if (itkV < 0)
        itkV = i; // The first good V track
    }
  }
  for (int i = 0; i < TTracks.size(); i++) {
    if (TTracks[i].Good) {
      nTtkrs++;
      if (itkT < 0)
        itkT = i; // The first good T track
    }
  }
  if (nVtkrs == 0 || nTtkrs == 0)
    nTracks = 0;
  else
    nTracks = (nVtkrs > nTtkrs) ? nVtkrs
                                : nTtkrs; // Generally we will reject events if this number is not exactly equal to 1.
};

// Method to print out the list of tracks
void pCT_Tracking::dumpTracks(int eventNumber) {
  std::cout << "pCT_Tracking::dumpTracks: dump of the list of detected tracks for event " << eventNumber
            << ".  Number of tracks=" << nTracks << std::endl;
  std::cout << "    The first good V track is " << itkV << " and the first good T track is " << itkT << std::endl;
  std::cout << "    There are " << VTracks.size() << " tracks in the V view:" << std::endl;
  for (int i = 0; i < VTracks.size(); i++) {
    if (true) {
      //            if (VTracks[i].Good) {
      std::cout << "    Track number " << i << " Quality=" << VTracks[i].Q << "  Mismatch at u=0 is " << VTracks[i].Miss
                << " mm "
                << " good=" << VTracks[i].Good << std::endl;
      for (int lyr = 0; lyr < 4; lyr++) {
        std::cout << "          Hit on layer " << lyr << "   U=" << VTracks[i].U[lyr] << "   V=" << VTracks[i].X[lyr]
                  << "  Q=" << VTracks[i].Qh[lyr] << std::endl;
      }
    }
  }
  std::cout << "    There are " << TTracks.size() << " tracks in the T view:" << std::endl;
  for (int i = 0; i < TTracks.size(); i++) {
    if (true) {
      //            if (TTracks[i].Good) {
      std::cout << "    Track number " << i << " Quality=" << TTracks[i].Q << "  Mismatch at u=0 is " << TTracks[i].Miss
                << " mm "
                << " good=" << TTracks[i].Good << std::endl;
      for (int lyr = 0; lyr < 4; lyr++) {
        std::cout << "          Hit on layer " << lyr << "   U=" << TTracks[i].U[lyr] << "   T=" << TTracks[i].X[lyr]
                  << "  Q=" << TTracks[i].Qh[lyr] << std::endl;
      }
    }
  }
};

// Method to display all good tracks in an event, using Gnuplot.
void pCT_Tracking::displayEvent(int eventNumber, TkrHits &pCThits, std::string outputDir) {
  char fn[80] = "";
  char Q[5] = "ebg ";
  sprintf(fn, "%s/Event%d.gp", outputDir.c_str(), eventNumber);
  std::cout << "**** pCT_Tracking::displayEvent: Writing a single event display for event number" << eventNumber
            << " to file " << fn << std::endl;
  FILE *oFile = fopen(fn, "w");
  if (oFile != NULL) {
    fprintf(oFile, "#*** This file is intended to be displayed by gnuplot.\n");
    fprintf(oFile, "#*** Either double click on the file (works in Windows at least),\n");
    fprintf(oFile, "#*** or else start up gnuplot and use the load command to display the plot.\n");
    fprintf(
        oFile,
        "set multiplot title 'pCT Tracking Single Event Display for Event %d' layout 2,1 columnsfirst scale 1.0,1.0\n",
        eventNumber);
    fprintf(oFile, "set title 'Tracker V View'\n");
    fprintf(oFile, "set xlabel 'U (mm)'\n");
    fprintf(oFile, "set ylabel 'V (mm)'\n");
    fprintf(oFile, "set xrange [-250.:250.]\n");
    fprintf(oFile, "set yrange [-50.:50.]\n");
    fprintf(oFile, "set nokey\n");
    fprintf(oFile, "plot '-' with labels\n");
    for (int tkr = 0; tkr < VTracks.size(); tkr++) {
      if (VTracks[tkr].Good) {
        for (int lyr = 0; lyr < 4; lyr++) {
          fprintf(oFile, "%7.2f %7.2f %d%c\n", VTracks[tkr].U[lyr], VTracks[tkr].X[lyr], tkr,
                  Q[VTracks[tkr].Qh[lyr]]); // plot hits on tracks
        }
      }
    }
    for (int lyr = 0; lyr < 4; lyr++) {
      for (int hit = 0; hit < pCThits.Lyr[lyr].N[0]; hit++) {
        int tk = pCThits.Lyr[lyr].F[0].at(hit);
        if (tk >= 0) {
          if (!VTracks[tk].Good)
            tk = -1;
        }
        if (tk < 0) { // plot unused hits
          fprintf(oFile, "%7.2f %7.2f x\n", pCThits.Lyr[lyr].U[0].at(hit), pCThits.Lyr[lyr].Y[0].at(hit));
        }
      }
    }
    fprintf(oFile, "e\n");

    fprintf(oFile, "set title 'Tracker T View'\n");
    fprintf(oFile, "set xlabel 'U (mm)'\n");
    fprintf(oFile, "set ylabel 'T (mm)'\n");
    fprintf(oFile, "set xrange [-250.:250.]\n");
    fprintf(oFile, "set yrange [-150.:150.]\n");
    fprintf(oFile, "set nokey\n");
    fprintf(oFile, "plot '-' with labels\n");
    for (int tkr = 0; tkr < TTracks.size(); tkr++) {
      if (TTracks[tkr].Good) {
        for (int lyr = 0; lyr < 4; lyr++) {
          fprintf(oFile, "%7.2f %7.2f %d%c\n", TTracks[tkr].U[lyr], TTracks[tkr].X[lyr], tkr, Q[TTracks[tkr].Qh[lyr]]);
        }
      }
    }
    for (int lyr = 0; lyr < 4; lyr++) {
      for (int hit = 0; hit < pCThits.Lyr[lyr].N[1]; hit++) {
        int tk = pCThits.Lyr[lyr].F[1].at(hit);
        if (tk >= 0) {
          if (!TTracks[tk].Good)
            tk = -1;
        }
        if (tk < 0) { // plot unused hits
          fprintf(oFile, "%7.2f %7.2f x\n", pCThits.Lyr[lyr].U[1].at(hit), pCThits.Lyr[lyr].Y[1].at(hit));
        }
      }
    }
    fprintf(oFile, "e\n");
    fprintf(oFile, "unset multiplot\n");
    fprintf(oFile, "show label");
    fclose(oFile);
  } else
    std::cout << "pCT_Tracking::displayEvent: unable to open file for display of event " << eventNumber << std::endl;
}
