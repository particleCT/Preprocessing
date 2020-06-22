// Tracker pattern recognition routines
// R.P. Johnson 5/22/2016
// ****************************************************
// This class first creates so called tracks from all 
// Hits recorded for an event on all 4 tracker layers
// If multiple Hits are found, the track with the 
// smallest difference in backprojected front and rear
// tracker vectors at isocenter is "good". 
// If there is missing front/back vectors: try to
// use the existing back vector with either gaps or
// the vertex position.
// 
// Lenny June 2020 
// 
//
//
//
//*****************************************************
#include "pCT_Tracking.h"
#include "pCTconfig.h" 
#include "pCTcut.h"


pCT_Tracking::pCT_Tracking(TkrHits &pCThits, pCTgeo* Geometry) { // This class constructor finds all the tracks

  theCuts = pCTcut::GetInstance();  
  theConfig= pCTconfig::GetInstance();

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



//Reprint from the TkrHits class for undertsanding
/*
    int N[2];                       // Number of hits in each of the V and T views for a given layer
    std::vector<double> Y[2], U[2]; // 0=V and 1=T
    std::vector<int> F[2];          // Track number; -1 if not use
*/
//End For Undertsanding only 


  // Create lists of all vectors in the front and back trackers, within the specified slope cuts
  std::vector<Tkr2D> tmp;
  std::vector<vct> front, back;
  bool GoodTk = false; // See if there is any useable tracks, if there is any, the code terminates as soon as it is found (to avoid loosing good tracks by also looking for hit+gap tracks)


  Y0_candidate,U0;
  Y1_candidate,U1;  
  Y2_candidate,U2;  
  Y3_candidate,U3;  
 
  
   
//-----------------------------------------------------------------------------
// First list all possible front/rear vectors from the recorded hit clusters 
//----------------------------------------------------------------------------
  for (int i = 0; i < pCThits.Lyr[0].N[Idx]; i++) {
    for (int j = 0; j < pCThits.Lyr[1].N[Idx]; j++) {

      U0 = pCThits.Lyr[0].U[Idx].at(i); // U position of the layers always the same for all hits
      U1 = pCThits.Lyr[1].U[Idx].at(j); // U position of the layers always the same for all hits

      Y0_candidate = pCThits.Lyr[0].Y[Idx].at(i);
      Y1_candidate = pCThits.Lyr[1].Y[Idx].at(j);

 
      double slope = (Y1_candidate - Y0_candidate) /
                     (U1 - U0);

      if(theCuts->cutHitSlope(0,Idx,slope)) continue; // Check if the slope of the hits matches the expected maximum slope (front=0/back=1,Idx,slope)

      //save this entrance vector candidate
      vct vctmp;
      vctmp.slope = slope;
      vctmp.Y[0] = Y0_candidate;
      vctmp.U[0] = U0; 
      vctmp.I[0] = i; // Hit id 
      vctmp.Y[1] = Y1_candidate; 
      vctmp.U[1] = U1; 
      vctmp.I[1] = j; // Hit id
      vctmp.intercept = vctmp.Y[0] - slope * vctmp.U[0];
      front.push_back(vctmp);

      
    }
  }
  for (int i = 0; i < pCThits.Lyr[2].N[Idx]; i++) {
    for (int j = 0; j < pCThits.Lyr[3].N[Idx]; j++) {

      //U view (always the same but reset here, to accomodate the case U[Idx].size()==0)
      U2 = pCThits.Lyr[2].U[Idx].at(i); 
      U3 = pCThits.Lyr[3].U[Idx].at(j); 

      //T or V view (Idx=0 = T)
      Y2_candidate = pCThits.Lyr[2].Y[Idx].at(i);
      Y3_candidate = pCThits.Lyr[3].Y[Idx].at(j);



      double slope = (Y3_candidate - Y2_candidate) /
                     (U3 - U2);


      if(theCuts->cutHitSlope(1,Idx,slope)) continue;
    
      //save this exit vector candidate
      vct vctmp;
      vctmp.slope = slope;
      vctmp.Y[0] = Y2_candidate; 
      vctmp.U[0] = U2; 
      vctmp.I[0] = i;
      vctmp.Y[1] = Y3_candidate; 
      vctmp.U[1] = U3; 
      vctmp.I[1] = j;
      vctmp.intercept = vctmp.Y[0] - slope * vctmp.U[0];
      back.push_back(vctmp);

    }
  }
  //-----------------------------------------------------------------------------
  // End grouping of hits for track assembling. Results are stored in back and front vectors. 
  //-----------------------------------------------------------------------------


  if (front.size() == 0 && back.size() == 0) return tmp; // No useful tracks front and back, skip the event!

  if(theConfig->item_int["MultiTrackReject"]){ //Added, Lenny June 2020
    if (front.size()>1 || back.size()>1) return tmp; //Test // If there is more than one reasonable track possible for either front or rear, skip this event.
  }

  //-----------------------------------------------------------------------------
  // Find any front/back vector combination from the hits
  //----------------------------------------------------------------------------
    
  // Look at all combinations of front and rear vectors for tracks that match in the center
  for (int i = 0; i < front.size(); i++) {
    for (int j = 0; j < back.size(); j++) {


      miss = front[i].intercept - back[j].intercept; // distance of forward/backprojected tracks at isocenter

      if(theCuts->cutTrackIsocenterIntercept(miss)) continue;

        // Check whether any of the hits have already been used. Do not allow hit sharing.
        if (tmp.size() > 0) { // Don't waste time if there are no tracks already (the usual situation)

          bool reject = false;
          std::vector<int> shareList;

          for (int k = 0; k < 2; k++) { // k goes over the layers
            int tk1 = pCThits.Lyr[k].F[Idx].at(front[i].I[k]); //Check for any track that used a hit of teh proposed track already

            if (tk1 >= 0) { // F is set down below if a track is found, else it is -1
              if (tmp[tk1].Good) { 
                if (abs(miss) > abs(tmp[tk1].Miss)) { // Use the track with the best match
                  reject = true;
                  break; // This track is inferior
                }
                shareList.push_back(tk1); // store the inferior track to reject it later on
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

          if (reject) continue; // This track is worse than some tracks already found that shares hits, so skip this one. // FORMERLY break -> Pot. Bug

          for (int shr = 0; shr < shareList.size(); shr++) {
            tmp[shareList[shr]].Good = false; // Get rid of earlier inferior tracks sharing hits with this one
          }
        }

        Tkr2D tmpTk;
        tmpTk.X[0] = front[i].Y[0];
        tmpTk.U[0] = front[i].U[0]; 

        tmpTk.Qh[0] = 3; // Hit quality! 
        pCThits.Lyr[0].F[Idx].at(front[i].I[0]) = tmp.size(); // Flags hit as used and enables to search for the track corresponding to hit // not used anymore! 
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

	// Temporary track storage
        tmp.push_back(tmpTk);

	GoodTk = true; // There is at least 1 good track 
      
    }
  }
  //-----------------------------------------------------------------------------
  //Done with track matching
  //-----------------------------------------------------------------------------

  //Check if there was any good track already. 
  //If there was, there is no need for dubious gap+hit or vertex plus hit tracks
  //Added by Lenny
  if(GoodTk) return tmp;


  //-----------------------------------------------------------------------------
  // if No good tracks could be made, make the tracks from gaps or vertex extrapolation
  //-----------------------------------------------------------------------------

  //No possible front tracks, but rear tracks
  
  if(front.size()==0 && back.size()>0){

  //---------------------------------------------------------------------------
  // Try to make vectors in the front with an unused hit and a gap
  //---------------------------------------------------------------------------

  if (Idx == 1) { // T layers only, where the cracks are parallel to the strips. This won't work for V layers.
    for (int i = 0; i < back.size(); i++) {
      if (pCThits.Lyr[2].F[Idx].at(back[i].I[0]) >= 0 || pCThits.Lyr[3].F[Idx].at(back[i].I[1]) >= 0) continue; // Unused back vectors only (F = -1)

        int otherLyr[2] = { 1, 0 };
	
        for (int lyr = 0; lyr < 2; lyr++) { // Loop over both front layers

          for (int j = 0; j < pCThits.Lyr[lyr].N[Idx]; j++) { // number of Hits

/*            int tk = pCThits.Lyr[lyr].F[Idx].at(j); // is there a track with any of these hits already?

            if (tk >= 0) {
              if (tmp[tk].Good)
                continue; // if a possible track with a hit in this layer is good, don't look for further hit+"gap". //formerly break
            }
*///Can be skipped with condition introduced above
	    
            double t1 = pCThits.Lyr[lyr].Y[Idx].at(j);
            double u1 = pCThits.Lyr[lyr].U[Idx].at(j);

            for (int gap = 0; gap < 3; gap++) {

              double t2 = Geometry->tBoardGap(gap, otherLyr[lyr]);
              double u2 = Geometry->uT(otherLyr[lyr]);

              double slope = (t2 - t1) / (u2 - u1);

              if (theCuts->cutHitSlope(0,Idx,slope)) continue; // front = 0 in the cuts function

              double intercept = t1 - slope * u1;
              miss = intercept - back[i].intercept;

              if (theCuts->cutTrackIsocenterIntercept(miss)) continue;
	      

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
 	      GoodTk = true; 
	      goto nextBackVector; //found a potential solution for this one already so skip any further possibilities. //FIXME get rid of goto's
	    } 
          }
        }
    nextBackVector:
      ;
    }
  }
  //--------------------------------------------------------------------------
  // End find possible front+gap + back vector tracks
  //--------------------------------------------------------------------------

  //For now I assume if a hit+gap track is possible, that it is more likely than a particle passing close to the gap and having a missing hit
  //This can also be updated by using either one of the hit+gap or hit+vertex solutions that give the smallest distance at isocenter. 
  if(GoodTk) return tmp;  

  //-----------------------------------------------------------------------------
  // Consider all possibile tracks using missing hits and the beam vertex
  //-----------------------------------------------------------------------------

  // Now consider track candidates missing one hit in the front tracker.  Make a vector from the hit that
  // is present that points back to the putative beam origin (this works a bit better at LLUMC, using the
  // lead foil, than at the Chicago Proton Center).
  for (int i = 0; i < back.size(); i++) {

    if (pCThits.Lyr[2].F[Idx].at(back[i].I[0]) >= 0 || pCThits.Lyr[3].F[Idx].at(back[i].I[1]) >= 0) continue; //skip if there was a track with that back vector already

      int otherLyr[2] = { 1, 0 };

      for (int lyr = 0; lyr < 2; lyr++) {
        for (int j = 0; j < pCThits.Lyr[lyr].N[Idx]; j++) {
/*          int tk = pCThits.Lyr[lyr].F[Idx].at(j);
          if (tk >= 0) {
            if (tmp[tk].Good)
              continue; // there was a track with that Hit already in combination //formerly break
          }
*/

          double t1 = pCThits.Lyr[lyr].Y[Idx].at(j);
          double u1 = pCThits.Lyr[lyr].U[Idx].at(j);
          double u2;

          if (Idx == 1) u2 = Geometry->uT(otherLyr[lyr]);
          else u2 = Geometry->uV(otherLyr[lyr]);

          double slope = (t1 - Geometry->BeamVertex(Idx)[Idx]) / (u1 - Geometry->BeamVertex(Idx)[2]);
          double intercept = Geometry->BeamVertex(Idx)[Idx] - slope * Geometry->BeamVertex(Idx)[2];

          miss = intercept - back[i].intercept;
          if (theCuts->cutTrackIsocenterIntercept(miss)) continue;

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
    nextBkVctr:
    ;
    }
  }


  //Existing front vector but no back tracker vector
  if( back.size()==0 && front.size()>0){

  //----------------------------------------------------------------------------
  // Match up unused front tracker vectors with hits from the back+gap
  //----------------------------------------------------------------------------
    for (int i = 0; i < front.size(); i++) {
      if (pCThits.Lyr[0].F[Idx].at(front[i].I[0]) >=0 || pCThits.Lyr[1].F[Idx].at(front[i].I[1]) >= 0) continue; // this front vector was connnetced to a good track already

        int otherLyr[4] = { 0, 0, 3, 2 };

        for (int lyr = 2; lyr < 4; lyr++) {
          for (int j = 0; j < pCThits.Lyr[lyr].N[Idx]; j++) {

/*            int tk = pCThits.Lyr[lyr].F[Idx].at(j); // Was there any track with the hit already? 
            if (tk >= 0) {
              if (tmp[tk].Good)
                break; //continue; //If there was any good hit in this layer, no need to look into gaps // FORMERLY loop break!
            }
*///Can be skipped with above condition

            double t1 = pCThits.Lyr[lyr].Y[Idx].at(j);
            double u1 = pCThits.Lyr[lyr].U[Idx].at(j);

            for (int gap = 0; gap < 3; gap++) { //Check for all combinations of gaps with the unused hit

              double t2 = Geometry->tBoardGap(gap, otherLyr[lyr]);
              double u2 = Geometry->uT(otherLyr[lyr]);
              double slope = (t2 - t1) / (u2 - u1);

              if (theCuts->cutHitSlope(1,Idx,slope)) continue;

              double intercept = t1 - slope * u1;
              double miss = intercept - front[i].intercept;
              if (theCuts->cutTrackIsocenterIntercept(miss)) continue;

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
	      GoodTk = true; 
 	      
	      goto nextFrontVector;	
          }
        }
      }
    nextFrontVector:
      ;
    }
  
    if(GoodTk) return tmp; //Don't bother with anything more if a good track was found //same as above, this might be improved by forcing the track with the smallest isocenter distance.

    // Finally, if there is an unused front vector but only one unused hit in the back, try to make a track by
    // extrapolation or interpolation
    // These will be the lowest quality tracks.
    for (int i = 0; i < front.size(); i++) {
      if (pCThits.Lyr[0].F[Idx].at(front[i].I[0]) < 0 && pCThits.Lyr[1].F[Idx].at(front[i].I[1]) < 0) {
        int otherLyr[4] = { 0, 0, 3, 2 };
        for (int lyr = 2; lyr < 4; lyr++) {
          for (int j = 0; j < pCThits.Lyr[lyr].N[Idx]; j++) {
/*            int tk = pCThits.Lyr[lyr].F[Idx].at(j);
            if (tk >= 0) {
              if (tmp[tk].Good && tmp[tk].Q > 0)
                break; // Skip hits already used on better tracks
            }*/
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
            if (theCuts->cutHitSlope(1,Idx,slope)) continue;
            double intercept = t3 - slope * u3;
            miss = intercept - front[i].intercept;
            if (theCuts->cutTrackIsocenterIntercept(miss)) continue;

/*                if (tk >= 0) {
                  if (tmp[tk].Good) {
                    if (abs(tmp[tk].Miss) > abs(miss)) {
                      tmp[tk].Good = false; // This track is better than an existing track using the same hit
                    } else {
                      continue; // An existing track using this hit is better
                    }
                  }
                }
*/  
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
  return tmp;
} // End of Tracking2D

// Method to print out the list of tracks
void pCT_Tracking::dumpTracks(int eventNumber) {
  std::cout << "pCT_Tracking::dumpTracks: dump of the list of detected tracks for event " << eventNumber
            << ".  Number of tracks=" << nTracks << std::endl;
  std::cout << "    The first good V track is " << itkV << " and the first good T track is " << itkT << std::endl;
  std::cout << "    There are " << VTracks.size() << " tracks in the V view:" << std::endl;
  for (int i = 0; i < VTracks.size(); i++) {
  //            if (VTracks[i].Good) {
    std::cout << "    Track number " << i << " Quality=" << VTracks[i].Q << "  Mismatch at u=0 is " << VTracks[i].Miss
              << " mm "
              << " good=" << VTracks[i].Good << std::endl;
    for (int lyr = 0; lyr < 4; lyr++) {
      std::cout << "          Hit on layer " << lyr << "   U=" << VTracks[i].U[lyr] << "   V=" << VTracks[i].X[lyr]
                << "  Q=" << VTracks[i].Qh[lyr] << std::endl;
    }
  }
  std::cout << "    There are " << TTracks.size() << " tracks in the T view:" << std::endl;
  for (int i = 0; i < TTracks.size(); i++) {
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
