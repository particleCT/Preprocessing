#include <iostream>
#include <fstream>
#include <math.h>       /* sqrt */
#include <cmath>        // std::abs
#include <stdlib.h>  
#include "TVector3.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMatrixD.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TBrowser.h"
#include "TRandom.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;
struct Particle{
  Float_t x0,y0,z0,px0,py0,pz0;
  Float_t x1,y1,z1,px1,py1,pz1;
  Float_t E[5];
  Float_t E_tot;
  Float_t wepl, WET_prob;
  Int_t MaxEnergyTransFilter;
  Int_t ThresholdFilter;
  Int_t dEEFilter;
  Int_t TS =0;
  Int_t PR =0;
  Float_t t[4];
  Float_t v[4];
};

int NStep = 250;
double det_Xmin = -16.43; // cm
double det_Xmax =  16.43; // cm
double det_Zmin = -8.76/2; // cm
double det_Zmax =  8.76/2; // cm
double det_Ymin = -35.2/2; // cm
double det_Ymax =  35.2/2; // cm

double reco_Xmin = -16.00/2.; // cm
double reco_Xmax =  16.00/2.; // cm
double reco_Ymin = -16.00/2.; // cm
double reco_Ymax =  16.00/2.; // cm
double reco_Zmin = -3.00/2.; // cm
double reco_Zmax =  3.00/2.; // cm

double midX = (det_Xmax+det_Xmin)/2.; // cm
double midY = (det_Ymax+det_Ymin)/2.; // cm
double midZ = (det_Zmax+det_Zmin)/2.; // cm   
double findWET(double, double);
bool getLineIntersection(double, double, double, double, double, double, double, double, double&, double&);
double EstimateLength(Particle Point, TH3D* Prior);
double extrap2D(double X[2], double Y[2], double Xnew);
bool GetLength(Particle Point, int Nbricks, double &Length,TProfile*);
int main(int, char** argv){
  Particle Point;
  TGraph* calWEPL;
  char*  priorfilename = argv[2];
  TFile* priorfile = new TFile(priorfilename,"update");
  TH3D*  Prior = (TH3D*)priorfile->Get("test");

  char * datafilename = argv[1];
  TFile* datafile = new TFile(datafilename,"update");
  
  // Phasespace
  TTree* t = (TTree*)datafile->Get("phase");
  t->SetBranchAddress("x0",&Point.x0);
  t->SetBranchAddress("y0",&Point.y0);
  t->SetBranchAddress("z0",&Point.z0);

  t->SetBranchAddress("x1",&Point.x1);
  t->SetBranchAddress("y1",&Point.y1);
  t->SetBranchAddress("z1",&Point.z1);
  
  t->SetBranchAddress("px0",&Point.px0);
  t->SetBranchAddress("py0",&Point.py0);
  t->SetBranchAddress("pz0",&Point.pz0);

  t->SetBranchAddress("px1",&Point.px1);
  t->SetBranchAddress("py1",&Point.py1);
  t->SetBranchAddress("pz1",&Point.pz1);

  t->SetBranchAddress("t",&Point.t);
  t->SetBranchAddress("v",&Point.v);
  t->SetBranchAddress("E_tot",&Point.E_tot);
  
  t->SetBranchAddress("wepl",&Point.wepl);
  t->SetBranchAddress("MaxEnergyTransFilter",&Point.MaxEnergyTransFilter);
  t->SetBranchAddress("ThresholdFilter",&Point.ThresholdFilter);
  t->SetBranchAddress("dEEFilter",&Point.dEEFilter);

  TBranch* bpt = t->Branch("PriorFilter",&Point.PR,"Prior/I");  

  calWEPL = (TGraph*)datafile->Get("calWEPL/RangeVsEnergy");
  TProfile* hProf = new TProfile("Length","",500, -150, 150);
  TH2D* REhist_test_Normal = new TH2D("Test_Normal","Test",250,0,260, 1000, 0,260*4);
  TH2D* REhist_test_Prior = new TH2D("Test_Prior","Test",250,0,260, 1000, 0,260*4);
  TH2D* REhist_test_TS = new TH2D("Test_TS","Test",250,0,260, 1000, 0,260*4);
  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    Point.PR = false;

    // Base Filter
    if(i%100000==0) cout<<i<<endl;
    //if(Point.dEEFilter==1 &&
    //Point.ThresholdFilter==1 &&
    //Point.MaxEnergyTransFilter==1 &&
    //  Point.wepl>0 && Point.wepl < 260){
    //for(int j =0; j<5; j++)  Point.E_tot+= Point.E[j];
    double Length = 0;

    //For the wedge
    //if(GetLength(Point, atoi(argv[2]), Length, hProf)){

    //For the cylinder
    double WET = EstimateLength(Point,Prior);
    
    //double WET = Length*1.030;	
    double Epred = calWEPL->Eval(WET);
    Point.PR = 0;
    //if(Point.TS==1) REhist_test_TS->Fill(Length,Point.E_tot);
    
    if( abs(Point.E_tot - Epred) <30){
	Point.PR = true;
	REhist_test_Prior->Fill(WET,Point.E_tot);
	//if(abs(Point.z0)<10) 
	
    }
    if(abs(Point.z0)<10) REhist_test_Normal->Fill(WET,Point.E_tot);
	
	
	//} // For the wedge
  //}
    bpt->Fill();
  }
  REhist_test_Normal->Write("",TObject::kOverwrite);
  REhist_test_Prior->Write("",TObject::kOverwrite);
  REhist_test_TS->Write("",TObject::kOverwrite);
  hProf->Write("",TObject::kOverwrite);
  t->Write("",TObject::kOverwrite);
  datafile->Close();
}
//////////////////////////////////////////////////////////////////////
// 2D equation to find line intersection
//////////////////////////////////////////////////////////////////////
bool getLineIntersection(double p0u, double p0t, double p1u, double p1t,  // Two points on the first line
			 double p2u, double p2t, double p3u, double p3t,  // Two points on the second line
			 double &ipu, double &ipt) { // Return the intersection point
  if (p0u == p1u || p2u == p3u) return false;

  double slope1 = (p1t - p0t) / (p1u - p0u); // Slope of the first line
  double slope2 = (p3t - p2t) / (p3u - p2u); // Slope of the second line
  if (slope1 == slope2) return false;                   // There is no intersection
  double int1 = p0t - slope1 * p0u; // Intercept of the first line
  double int2 = p2t - slope2 * p2u; // Intercept of the second line

  ipu = (int1 - int2) / (slope2 - slope1); // intersection point U
  ipt = (slope2 * int1 - slope1 * int2) / (slope2 - slope1); // Intersection Point T
  return true;
}
//////////////////////////////////////////////////////////////////////
// 2D equation to project forward
//////////////////////////////////////////////////////////////////////
double extrap2D(double X[2], double Y[2], double Xnew) {
    double dX = X[1] - X[0];
    if (dX == 0.) {
      cout << "pCTgeo::extrap2D, division by zero; x values must be different." << endl;
      dX = 1.0e-23;
    }
    double slope = (Y[1] - Y[0]) / (dX);
    return Y[1] + (Xnew - X[1]) * slope;
}
//////////////////////////////////////////////////////////////////////
// 2D equation to project forward
//////////////////////////////////////////////////////////////////////
bool GetLength(Particle Point, int Nbricks, double &Length,TProfile* hProf){
  /// CALCULATE THE LENGTH
  //double Ut[2] , Uv[2] , T[2] , V[2];
  double Uft[2], Ufv[2], Tf[2], Vf[2];
  double Ust[2];
  double Tw1, Tw2, Tw3, Tw4;
  double tWedgeOffset = -50.8;
  Tw1 = -104.50 + tWedgeOffset; // Start of the wedge slope
  Tw2 = -4.75   + tWedgeOffset; // End of the slope, start of the flat
  Tw3 =  4.75   + tWedgeOffset; // End of the flat, start of the opposite slope
  Tw4 =  104.50 + tWedgeOffset; // End of the opposite slope.

  double BrickThickness =  50.8; // Brick thickness; also the wedge thickness
  Ust[0] = -2.5 * BrickThickness; // centered around the middle of the bricks, 4 bricks + 1 wedge
  Ust[1] = Ust[0] + BrickThickness;
  double Wbricks = BrickThickness * Nbricks;
  double Uout = Ust[1] + Wbricks; // U coordinate of the last brick, downstream side
  double tBrickEnd = Tw4 + 25.0;
  
  Uft[0] = -211.4; Uft[1] = -161.4;// Front tracker coordinates
  Ufv[0] = -217.3; Ufv[1] = -167.2;
  //Ut[0]  =  161.4; Ut[1]  =  211.4;// Rear tracker coordinates  
  //Uv[0]  =  167.2; Uv[1]  =  217.3;
  
  Tf[0] = Point.t[0] -1; Tf[1] = Point.t[1] - 1;// Front tracker coordinates
  //Tf[0] = Point.t[0]; Tf[1] = Point.t[1] ;// Front tracker coordinates  
  Vf[0] = Point.v[0]; Vf[1] = Point.v[1];
  //T[0]  = Point.t[2];  T[1] = Point.t[3];// Rear tracker coordinates
  //V[0]  = Point.v[2];  V[1] = Point.v[3];
  // First estimate of the u coordinate at entry into the phantom wedge
  double Tin  = extrap2D(Uft, Tf, Ust[0]);
  double Vin  = extrap2D(Ufv, Vf, Ust[0]);
  
  // First brick past the wedge
  double TinB = extrap2D(Uft, Tf, Ust[1]);
  double VinB = extrap2D(Ufv, Vf, Ust[1]);  
  //Back of bricks     
  double Tout = extrap2D(Uft, Tf, Uout);
  double Vout = extrap2D(Ufv, Vf, Uout);
  //double Tout = extrap2D(Ut, T, Uout);
  //double Vout = extrap2D(Uv, V, Uout);  
  double Uin  = Ust[0];
  // Before the negative slope
  if(TinB < Tw1 || Tin < Tw1) return false;
  
  // Negative Slope
  else if (TinB >= Tw1 &&TinB<=Tw2 && Tin>=Tw1 && Tin<=Tw2) {
    //return false;
    if(!getLineIntersection(-211.4, Point.t[0], -161.4, Point.t[1], Uin + BrickThickness, Tw1, Uin, Tw2, Uin, Tin)) return false;
    
  }
  // Flat part of the wedge // No need to change Vin
  else if(Tin > Tw2 && Tin < Tw3 && TinB > Tw2 && TinB < Tw3){
    //return false;
    Uin = Ust[0];

  }
  
  // Positive Slope
  else if (Tin<= Tw4 && Tin >= Tw3 && TinB <= Tw4 && TinB >=Tw3) {  
    if (!getLineIntersection( -211.4, Point.t[0], -161.4, Point.t[1], Uin, Tw3 ,Uin + BrickThickness, Tw4, Uin, Tin)) return false;
    //return false;
  }
  // Past the positive slope
  else if(TinB > Tw4 && TinB < tBrickEnd){
    //return false;
    Uin = Ust[1];
    Tin = TinB;
    Vin = VinB;
  }
  
  // Past the end of the brick
  else if(TinB > tBrickEnd){
    //return false;
    Tin = Tout;
    Vin = Vout;
    Uin = Uout;
  }
  else return false;
  Vin  = extrap2D(Ufv, Vf, Uin);    
  Length = sqrt( pow((Tin - Tout),2) + pow((Vin - Vout),2) + pow((Uin - Uout),2) );
  hProf->Fill(TinB,Length*1.03);
  return true;
}


////////////////////////////////////////////
// Compute Spline and Estimate Estop
////////////////////////////////////////////
double EstimateLength(Particle Point, TH3D* Prior){
  TVector3 p0(Point.px0,   Point.py0,   Point.pz0);
  TVector3 p1(Point.px1,   Point.py1,   Point.pz1);
  TVector3 m0(Point.x0, Point.y0, Point.z0);
  TVector3 m1(Point.x1, Point.y1, Point.z1);

  TVector3 m, m_old; // position of the path iterated
  TVector3 m_entry, m_exit; // position of the path at the entry and exit of the object
  TVector3 m_entry_proj, m_exit_proj; // Projection from the entry/exit position to the opposite tracker
  double t, t_entry = 0., t_exit = 1.; // fraction of the path at which point the spline enter the object
  m_entry = m0;
  m_exit  = m1;

  m_entry_proj.SetX( m1.x());
  m_entry_proj.SetY( m0.y() + p0.y()*(Point.x1-Point.x0) );
  m_entry_proj.SetZ( m0.z() + p0.z()*(Point.x1-Point.x0) );

  m_exit_proj.SetX( m0.x());
  m_exit_proj.SetY( m1.y() - p1.y()*(Point.x1-Point.x0) );
  m_exit_proj.SetZ( m1.z() - p1.z()*(Point.x1-Point.x0) );

  // Propagate from the entrance to the radius
  for(int k =1; k<NStep-1; k++){
    t_entry = double(k)/NStep;
    m_entry = m0 + t_entry*(m_entry_proj - m0);
    double RSP = Prior->GetBinContent(Prior->FindBin(m_entry.x(), m_entry.y(), m_entry.z()));
    if(RSP>1) break;
  }
  // Retro Propagate from the exit to the Hull
  for(int k =NStep; k>0; k--){
    t_exit = double(k)/NStep;
    m_exit = m_exit_proj + t_exit*(m1 - m_exit_proj);
    double RSP = Prior->GetBinContent(Prior->FindBin(m_exit.x(), m_exit.y(), m_exit.z()));
    if(RSP>1) break;    
  }

  if(t_entry >= t_exit) { // No Hull we are in air
    return 0;
    }
  double WER        = 25.693404; // 200 MeV -- cm
  double wepl       = Point.wepl/10; // mm -> cm
  double alpha1     = 1.01+0.43*pow(wepl/WER,2);
  double alpha2     = 0.99-0.46*pow(wepl/WER,2);
  double HullLength = TVector3(m_entry-m_exit).Mag();
  
  p0.SetMag(alpha1*HullLength);
  p1.SetMag(alpha2*HullLength);
  TVector3 A,B,C,D;
  A       =    2*m_entry - 2*m_exit + p0+p1;
  B       =   -3*m_entry + 3*m_exit - 2*p0-p1;
  C       =    p0;
  D       =    m_entry;
  
  m_old   =    m_entry;
  double WET = 0;
  double RSP;
  for(int i=0;i<NStep;i++){
    t = double(i)/NStep;
    m = D+t*(C+t*(B+t*A));
    RSP = Prior->GetBinContent(Prior->FindBin(m.x(), m.y(), m.z()));

    float L    = TVector3(m-m_old).Mag();
    WET       += L*RSP;
    /*cout<<"T:"<<t<<" t_entry:"<<t_entry<<" t_exit:"<<t_exit<<" WET:"<<WET<<" L:"<<L<<" RSP:"<<RSP<<" radius:"<<radius<<endl;
    cout<<"Vector X-Y-Z"<<m.x()<<" "<<m.y()<<" "<<m.z()<<endl;
    cout<<"Vector Entry X-Y-Z"<<m_entry.x()<<" "<<m_entry.y()<<" "<<m_entry.z()<<endl;
    cout<<"Vector Exit X-Y-Z"<<m_exit.x()<<" "<<m_exit.y()<<" "<<m_exit.z()<<endl;
    cout<<"\n"<<endl;*/
    m_old      = m;
  }
  return WET;

}
