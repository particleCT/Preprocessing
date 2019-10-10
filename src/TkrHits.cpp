// A class to find and store all of the tracker hits for an event, starting from
// raw clusters
#include "TkrHits.h"
//  Convert the raw tracker strip information into coordinates
TkrHits::TkrHits(pCTraw &pCTEvent, const pCTgeo &Geometry, bool print) {
  // Some arrays for temporary storage of merged clusters
  int cl1[110];
  int cln[110];
  int clcp[110];

  for (int lyr = 0; lyr < 4; lyr++) {
    Lyr[lyr].N[0] = 0;
    Lyr[lyr].N[1] = 0;
  }

  // Loop over the raw data, sorted by FPGA, merge clusters when broken between
  // chips, and calculate coordinates per layer
  for (int iFPGA = 0; iFPGA < 12; iFPGA++) {
    pCTraw::TrackerFPGA fpga = pCTEvent.tkr_fpga[iFPGA];
    if (fpga.num_chips != 0) {
      if (print) std::cout << "FPGA= " << iFPGA << "\n";
      int ncl = 0;
      int idx0new = -1;
      pCTraw::TrackerFPGA::TrackerChip prevChip = fpga.chip[0];
      for (int ichip = 0; ichip < fpga.num_chips; ichip++) {
        int idx0 = idx0new;
        idx0new = -1;
        pCTraw::TrackerFPGA::TrackerChip chip = fpga.chip[ichip];
        int chip_ID = chip.address;
        if (chip_ID < 0 || chip_ID > 11) {
          std::cout << "TkrHits **********Data Error********* Chip ID " << chip_ID << " is out of range!!\n";
          break;
        }
        if (print)
          std::cout << "  Chip " << chip_ID << " number of clusters=" << chip.num_clusts << "\n";
        if (chip.num_clusts > 10) {
          std::cout << "TkrHits ************Data Error********* Number of clusters = " << chip.num_clusts << "!!\n";
          break;
        }
        for (int icluster = chip.num_clusts - 1; icluster >= 0; icluster--) {
          pCTraw::TrackerFPGA::TrackerChip::TrackerCluster cluster = chip.cluster[icluster];
          // Check if the cluster is a continuation from the
          // previous chip, as chip boundaries can break up valid clusters
          bool merge = false;
          if (ichip != 0) {
            if (print)
              std::cout << "      Cluster " << icluster << "  First strip=" << cluster.first
                        << "  Number strips=" << cluster.length << "\n";
            if (chip_ID != 0 && chip_ID != 6) {
              pCTraw::TrackerFPGA::TrackerChip::TrackerCluster prevClus = prevChip.cluster[0];
              int lastStrip = cluster.first + cluster.length;
              int firstStrip = prevClus.first;
              int prevChip_ID = prevChip.address;
              if (print)
                std::cout << "        Last strip=" << lastStrip << " first strip previous chip=" << firstStrip
                          << " Chip=" << chip_ID << " Previous chip=" << prevChip_ID << "\n";
              if (lastStrip == 63 && prevChip_ID == chip_ID - 1 && firstStrip == 0 && cluster.length != 63) {
                if (idx0 < 0 || cl1[idx0] != 0 || chip_ID != clcp[idx0] + 1) {
                  std::cout << "**TkrHits: OOPS!  Logic screwed up, idx0=" << idx0 << " cl1=" << cl1[idx0]
                            << " ID=" << chip_ID << " clcp=" << clcp[idx0] << "\n";
                  std::cout << "**FPGA=" << iFPGA << " chip_ID=" << chip_ID << " ichip=" << ichip
                            << " icluster=" << icluster << "\n";
                  pCTEvent.dumpEvt();
                } else {
                  cl1[idx0] = cluster.first;
                  cln[idx0] = cluster.length + prevClus.length + 2;
                  clcp[idx0] = chip_ID;
                  merge = true;
                  if (print)
                    std::cout << "          Merging FPGA " << iFPGA << " chip " << chip_ID << " cluster "
                              << cluster.first << " " << cluster.length << " idx0=" << idx0 << "\n";
                }
              }
            }
          }
          if (!merge) {
            cl1[ncl] = cluster.first;
            cln[ncl] = cluster.length + 1;
            clcp[ncl] = chip_ID;
            if (cluster.first == 0)
              idx0new = ncl; // To facilitate overwriting a cluster when two are
                             // merged
            ncl++;
          }
        }
        prevChip = chip;
      }

      // Now calculate coordinates from the (merged) clusters for this FPGA
      for (int icl = 0; icl < ncl; icl++) {
        if (print)
          std::cout << "          Merged cluster " << icl << " 1st strip " << cl1[icl] << " length " << cln[icl]
                    << "  chip " << clcp[icl] << "\n";
        double Strip = float(cl1[icl]) + (float(cln[icl]) / 2.0) - 0.5;
        if (iFPGA < 4) {
          double Vcoord = Geometry.strip_position_V(iFPGA, clcp[icl], Strip);
          Lyr[Geometry.Lyr(iFPGA)].Y[0].push_back(Vcoord);
          Lyr[Geometry.Lyr(iFPGA)].U[0].push_back(Geometry.uV(Geometry.Lyr(iFPGA)));
          Lyr[Geometry.Lyr(iFPGA)].N[0]++;
          Lyr[Geometry.Lyr(iFPGA)].F[0].push_back(-1);
        } else {
          double Tcoord = Geometry.strip_position_T(iFPGA, clcp[icl], Strip);
          Lyr[Geometry.Lyr(iFPGA)].Y[1].push_back(Tcoord);
          Lyr[Geometry.Lyr(iFPGA)].U[1].push_back(Geometry.uT(Geometry.Lyr(iFPGA)));
          Lyr[Geometry.Lyr(iFPGA)].N[1]++;
          Lyr[Geometry.Lyr(iFPGA)].F[1].push_back(-1);
        }
      }
    }
  }
};
// Method to print the list of hits
void TkrHits::dumpHits(int eventNumber) {
  std::cout << "List of the coordinates in V and T layers for event " << eventNumber << ":" << std::endl;
  for (int lyr = 0; lyr < 4; lyr++) {
    std::cout << "    Layer number " << lyr << " has " << Lyr[lyr].N[0] << " V hits and " << Lyr[lyr].N[1]
              << " T hits\n";
    if (Lyr[lyr].N[0] <= 0)
      std::cout << "        No V coordinates found." << std::endl;
    else {
      std::cout << "        V coordinates: U=" << Lyr[lyr].U[0].at(0) << " V=";
      for (int i = 0; i < Lyr[lyr].N[0]; i++) {
        std::cout << " " << Lyr[lyr].Y[0].at(i) << " (" << Lyr[lyr].F[0].at(i) << ")";
      }
      std::cout << std::endl;
    }
    if (Lyr[lyr].N[1] <= 0)
      std::cout << "        No T coordinates found." << std::endl;
    else {
      std::cout << "        T coordinates: U=" << Lyr[lyr].U[1].at(0) << " T=";
      for (int i = 0; i < Lyr[lyr].N[1]; i++) {
        std::cout << " " << Lyr[lyr].Y[1].at(i) << " (" << Lyr[lyr].F[1].at(i) << ")";
      }
      std::cout << std::endl;
    }
  }
};
