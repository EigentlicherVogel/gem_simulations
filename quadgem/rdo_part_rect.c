#define UNUSED(expr)                                                           \
  do {                                                                         \
    (void)(expr);                                                              \
  } while (0)
#include "AvalancheMC.hh"
#include "AvalancheMicroscopic.hh"
#include "ComponentAnsys123.hh"
#include "MediumMagboltz.hh"
#include "Plotting.hh"
#include "Random.hh"
#include "Sensor.hh"
#include "TObjArray.h"
#include "ViewField.hh"
#include "ViewGeometry.hh"
#include "ViewSignal.hh"
#include "DriftLineRKF.hh"
#include <TApplication.h>
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TROOT.h>
#include <TVirtualViewer3D.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

#ifdef __MAKECINT__
#pragma link C++ class vector < float> + ;
#endif

using namespace Garfield;
using namespace std;
using namespace std::chrono;

auto start = high_resolution_clock::now();
/* Tree to store variables produced by Garfield++
 */
int main(int argc, char *argv[]) {
  UNUSED(argc);

  //TApplication app("app", &argc, argv);
  //plottingEngine.SetDefaultStyle();

  // ANSYS geometry variables
  Float_t kapton = 0.0050; // Thickness of the kapton layer, in cm
  Float_t metal = 0.0005;  // Thickness of the meta layers, in cm
  Float_t halfind = 0.05; // inputs?


  // Define variables.

  int nevents = 0; // Number of events.

  if (nevents < 1)
    nevents = 1000;

  vector<float> *ex2 = 0;
  vector<float> *ey2 = 0;
  vector<float> *ez2 = 0;
  vector<float> *et2 = 0;
  vector<float> *ee2 = 0;

  Double_t xi, yi, zi, ei, ti;
  Int_t rnele, rnelend;
  TString etreeName, itreeName;
  Double_t rx0, ry0, rz0, re0, rt0;
  std::vector<float> rx1, ry1, rz1, re1, rt1;
  std::vector<float> rx2, ry2, rz2, re2, rt2;

  std::vector<float> rixf, riyf, rizf, ritf;
  std::vector<float> rixi, riyi, rizi, riti;
  std::vector<float> rstat;
  std::vector<float> es, is;

  Int_t rentries, revent;
  Double_t rgain;
  Int_t rntotalend = 0;

  TString filename;
  filename = Form("OUTFILE/crosschev_single_%s.root", argv[1]);

  TFile *treefile = new TFile(filename, "RECREATE");
  treefile->cd();

  etreeName = "retree";
  TTree *retree = new TTree(etreeName, "RDO e- variables");

  // electron tree
  // retree->Branch("rentries",&rentries,"rentries/I");
  retree->Branch("revent", &revent, "revent/I");
  retree->Branch("rgain", &rgain, "rgain/D");
  retree->Branch("rnele", &rnele, "rnele/I");
  retree->Branch("rnelend", &rnelend, "rnelend/I");
  retree->Branch("rntotalend", &rntotalend, "rntotalend/I");
  retree->Branch("rez0", &rz0, "rez0/D");
  retree->Branch("rex0", &rx0, "rex0/D");
  retree->Branch("rey0", &ry0, "rey0/D");
  retree->Branch("ree0", &re0, "ree0/D");
  retree->Branch("ret0", &rt0, "ret0/D");

  retree->Branch("rstatus", &rstat);
  retree->Branch("rex1", &rx1);
  retree->Branch("rey1", &ry1);
  retree->Branch("ree1", &re1);
  retree->Branch("ret1", &rt1);
  retree->Branch("rez1", &rz1);
  retree->Branch("rex2", &rx2);
  retree->Branch("rey2", &ry2);
  retree->Branch("ree2", &re2);
  retree->Branch("ret2", &rt2);
  retree->Branch("rez2", &rz2);

  // Ion tree
  itreeName = "ritree";
  TTree *ritree = new TTree(itreeName, "rdo_board ions variables");

  // ritree->Branch("rentries",&rentries,"rentries/I");
  ritree->Branch("revent", &revent, "revent/I");
  ritree->Branch("rixi", &rixi);
  ritree->Branch("riyi", &riyi);
  ritree->Branch("riti", &riti);
  ritree->Branch("rizi", &rizi);
  ritree->Branch("rixf", &rixf);
  ritree->Branch("riyf", &riyf);
  ritree->Branch("ritf", &ritf);
  ritree->Branch("rizf", &rizf);

    //! Inputs
  cout << "Reading input" << endl;
  
  
  TString fileinp;
  fileinp = Form("OUTFILE/quadgem_%s_out.root", argv[1]);
  TFile *fin = new TFile(fileinp, "READ");
  fin->cd();

  TTree *etree = (TTree *)fin->Get("etree");

  etree->SetBranchAddress("ex2", &ex2);
  etree->SetBranchAddress("ey2", &ey2);
  etree->SetBranchAddress("ez2", &ez2);
  etree->SetBranchAddress("ee2", &ee2);
  etree->SetBranchAddress("et2", &et2);
  cout << "Read inputs" << endl;
  nevents = etree->GetEntries();

  // Gain histogram
  TH1D *h_gain = new TH1D("h_gain", "", 50, -0.1, 49.9);

  // Loading the Ansys files (weighting afterwards)
  ComponentAnsys123 *rdo_board = new ComponentAnsys123();

  rdo_board->Initialise("ANSYS_geo/ELIST_ro_pvoltf_kap.lis",
                  "ANSYS_geo/NLIST_ro_pvoltf_kap.lis",
                  "ANSYS_geo/MPLIST_ro_pvoltf_kap.lis",
                  "ANSYS_geo/PRNSOL_ro_pvoltf_kap.lis", "mm");
  rdo_board->EnableMirrorPeriodicityX();
  rdo_board->EnableMirrorPeriodicityY();
  //rdo_board->EnableMirrorPeriodicityY();
  rdo_board->PrintRange();
  
    //=========================================================
  // Create histograms for aspect ratio and element size by retrieving the dimesions of each mesh elemennt. This will enable us to know about the quality of the mesh
  //=======================================================
  TH1F* hAspectRatio = new TH1F("hAspectRatioE", "Aspect RatioE",100, 0., 50.);
  TH1F* hSize = new TH1F("hSizeE", "Element SizeE",100, 0., 30.);
  const int nel = rdo_board->GetNumberOfElements();
  double volume; // volume of the element 
  double dmin, dmax; // Min and max distance between nodal points
  for (int i = 0; i< nel; i++) {
    rdo_board->GetElement(i, volume, dmin, dmax);
    if (dmin > 0.) hAspectRatio->Fill(dmax / dmin);
    hSize->Fill(volume * 1.e9);
  }
  cout << "Mesh element analysis created" << endl;


  //=========================================================
  // Create histograms for aspect ratio and element size by retrieving the
  // dimesions of each mesh elemennt. This will enable us to know about the
  // quality of the mesh
  //=======================================================

  

 // rdo_board->SetMagneticField(0,by,bz);
  TCanvas* cp1 = new TCanvas("cp1", "", 600, 600);
      TCanvas* cf1 = new TCanvas("cf1", "", 600, 600);
      TCanvas* cf1_0 = new TCanvas("cf1_0", "", 600, 600);
      TCanvas* cf1_1 = new TCanvas("cf1_1", "", 600, 600);
      TCanvas* cp2 = new TCanvas("cp2", "", 600, 600);
      TCanvas* cf2 = new TCanvas("cf2", "", 600, 600);
      TCanvas* cf2_1 = new TCanvas("cf2_1", "", 600, 600);
      const double pitch=0.0200; // in cm
      const bool plotField =true;
      if (plotField) {
        ViewField* fieldView = new ViewField();
        fieldView->SetComponent(rdo_board);
        // Plot the potential along the hole axis.
        //TCanvas* cp = new TCanvas("cp", "", 600, 600);
        cp1->SetLeftMargin(0.16);
        fieldView->SetCanvas(cp1);
        fieldView->PlotProfile(0., 0., 0.02, 0., 0., -0.02);
        
       
        fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
        //fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);
        fieldView->SetArea(-2.*pitch , 0., 2.*pitch , 0.08);
	fieldView->SetElectricFieldRange(0., 2000); // gem example
        //TCanvas* cf = new TCanvas("cf", "", 600, 600);
        cf1->SetLeftMargin(0.16);
        fieldView->SetCanvas(cf1);
        fieldView->PlotContour("e");

	ViewField* fieldView0= new ViewField();
        fieldView0->SetComponent(rdo_board);
        // Plot the potential along the hole axis.                                                                       
	//TCanvas* cp = new TCanvas("cp", "", 600, 600);                                                            
        fieldView0->SetPlane(-1., 0., 0., 0., 0., 0.);
        //fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);
	fieldView0->SetArea(-2.*pitch , 0., 2.*pitch , 0.08);
        fieldView0->SetElectricFieldRange(0., 2000); // gem example
                                                         
        //TCanvas* cf = new TCanvas("cf", "", 600, 600);        
	cf1_0->SetLeftMargin(0.16);
        fieldView0->SetCanvas(cf1_0);
        fieldView0->PlotContour("e");


	
	ViewField* fieldView1 = new ViewField();
        fieldView1->SetComponent(rdo_board);
	// Set the viewing plane (normal vector) and ranges                                                         
	fieldView1->SetPlane(0., 0., -1., 0., 0., 0.005);
        //fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);                                                   
	fieldView1->SetArea(-45.*pitch , -10.*pitch, 45.*pitch , 10.0*pitch);
        fieldView1->SetElectricFieldRange(0., 2000.);                                                          
	//TCanvas* cf = new TCanvas("cf", "", 600, 600);    
	cf1_1->SetLeftMargin(0.16);
        fieldView1->SetCanvas(cf1_1);
        fieldView1->PlotContour("e");



        ViewField* fieldView2 = new ViewField();
        fieldView2->SetComponent(rdo_board);
        // Plot the potential along the hole axis.
        cp2->SetLeftMargin(0.16);
        fieldView2->SetCanvas(cp2);
        fieldView2->PlotProfile(0., 0., 0.15, 0., 0., -0.15);
        
        // Set the viewing plane (normal vector) and ranges
        fieldView2->SetPlane(0., -1., 0., 0., 0., 0.);
        //fieldView2->SetArea(-pitch / 2., -0.11, pitch / 2., 0.1);
        fieldView2->SetArea(-2.*pitch, 0., 2.*pitch , 0.08);
        fieldView2->SetVoltageRange(0., -200);
        //TCanvas* cf = new TCanvas("cf", "", 600, 600);
        cf2->SetLeftMargin(0.16);
        fieldView2->SetCanvas(cf2);
        fieldView2->PlotContour();

	ViewField* fieldView2_1 = new ViewField();
        fieldView2_1->SetComponent(rdo_board);
	fieldView2_1->SetPlane(0., 0, -1., 0., 0., 0.);
        //fieldView2->SetArea(-pitch / 2., -0.11, pitch / 2., 0.1);                                                       
        fieldView2_1->SetArea(-45.*pitch , -10.*pitch, 45.*pitch , 10.0*pitch);
        fieldView2_1->SetVoltageRange(0., -200);
        //TCanvas* cf = new TCanvas("cf", "", 600, 600);                                                                  
        cf2_1->SetLeftMargin(0.16);
        fieldView2_1->SetCanvas(cf2_1);
        fieldView2_1->PlotContour();

      }
  cout << "Mesh quality analysis created" << endl;

  // Gas setup
  MediumMagboltz *gas = new MediumMagboltz();
  gas->SetComposition("ar", 70., "co2", 30.);
  // gas->SetComposition("ar",45.,"co2",15.,"cf4",40.0);
  // gas->SetComposition("ar",80.,"cf4",20.);
  // gas->SetComposition("cf4",100.);
  // gas->SetComposition("ar",95.,"ch4",5.)
  gas->SetTemperature(293.5);
  gas->SetPressure(760.0);
  // gas->SetMaxElectronEnergy(200.);
  // gas->EnablePenningTransfer(0.51, 0.0, "Ne");
  gas->EnablePenningTransfer(0.46, 0.0, "ar");
  gas->Initialise(true);
  // gas->Initialise();
  int numofmat = rdo_board->GetNumberOfMaterials();
  cout << "Number of materials: " << numofmat << endl;

  //const std::string path = getenv("GARFIELD_HOME");
  gas->LoadIonMobility("Ion_Data/IonMobility_Ar+_Ar.txt");

  for (int i = 0; i < numofmat; i++) {
    const double eps = rdo_board->GetPermittivity(i);
    cout << "The permittivity of material " << i << " is " << eps << endl;
    if (fabs(eps - 1) < 1.e-13) {
      rdo_board->SetMedium(i, gas);
    }
    // if (eps == 1.) rdo_board->SetMedium(i, gas);
  }

  Sensor *sensor = new Sensor();
  sensor->AddComponent(rdo_board);

  /*
  sensor->SetArea(-5 * pitch, -5 * pitch, -0.02,
  5 * pitch,  5 * pitch,   0.02  );
  */
    
  AvalancheMicroscopic *aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);
  //drift_e->SetMaximumStepSize(0.001);

  // Set the ion transport
  AvalancheMC *drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);

      const bool plotDrift = false;
      ViewDrift* driftView = new ViewDrift();
      if (plotDrift) {
	driftView->SetArea(-1.0* pitch, - 1.0*pitch, -0.02,
			   1 * pitch,  1 * pitch,  0.02);
	aval->EnablePlotting(driftView);
	drift->EnablePlotting(driftView);
      }
  
  


  // aval->EnableDebugging();
  // aval->EnableMagneticField();
  // aval->EnableAvalancheSizeLimit(1000);
  rentries = 0;
  for (int i = 0; i < nevents; i++) {

    etree->GetEntry(i);
    revent = i;

    for (int j = 0; j < ex2->size(); j++) {

      if ((ez2->at(j) > -0.4))
        continue;
      xi = ex2->at(j);
      yi = ey2->at(j);
      // shift up half ind offset + half ind remaining gap
      zi = ez2->at(j) - ez2->at(j) + 0.049; 
      ei = ee2->at(j);
      ti = et2->at(j);
      cout << " z3land :" << ez2->at(j) << " input z gem3 :" << zi
           << " xi :" << ex2->at(j) << " yi :" << ey2->at(j) << endl;
      double dx0 = 0.0;
      double dy0 = 0.0;
      double dz0 = 0.0;

      float timeused = (float)clock() / CLOCKS_PER_SEC;
      if (timeused > 25000000.0) {
        cout << "Time limit reached, " << timeused << " sec" << endl;
        break;
      }
      // Start e- avalanche
      //int ne, ni;
      // Get # of e-s and ions produced in avalanche
      //aval->GetAvalancheSize(ne, ni);
      //rnele = ne;

      // Get the history of each avalanche e- (starting and end point of
      // avalanche e-)

    aval->AvalancheElectron(xi, yi, zi, ti, ei, dx0, dy0, dz0);
    int ne, ni;
    
    aval->GetAvalancheSize(ne, ni);
    rnele = ne;
      double xa0, ya0, za0, ea0, ta0;
      double xa1, ya1, za1, ea1, ta1;
      double xi1, yi1, zi1, ti1;
      double xi2, yi2, zi2, ti2;
      int istat;
    rnelend = aval->GetNumberOfElectronEndpoints();
    cout << "nelend :" << rnelend << endl;

      // entries ++;
      Int_t e_anode = 0;
      for (int k = 0; k < rnelend; k++) {

      aval->GetElectronEndpoint(k, xa0, ya0, za0, ta0, ea0, xa1, ya1, za1, ta1,
                                ea1, istat);

        rx1.push_back(xa0);
        ry1.push_back(ya0);
        rz1.push_back(za0);
        rt1.push_back(ta0);
        re1.push_back(ea0);

        // Store end point of avalanche e-s and primaries
        rx2.push_back(xa1);
        ry2.push_back(ya1);
        rz2.push_back(za1);
        rt2.push_back(ta1);
        re2.push_back(ea1);
        // status of avalanche e-
        rstat.push_back(istat);

        drift->DriftIon(xa0, ya0, za0, ta0);
        drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, istat);
        rixi.push_back(xi1);
        riyi.push_back(yi1);
        rizi.push_back(zi1);
        riti.push_back(ti1);
        rixf.push_back(xi2);
        riyf.push_back(yi2);
        rizf.push_back(zi2);
        ritf.push_back(ti2);

        rntotalend++;
        // if(za1 < -1.0* (metal+ kapton/2.+0.04))e_anode++ ; //rdo_board bottom?

        if (za1 < -0.03)
          e_anode++; // rdo only

      } // e- loop

      if (e_anode != 0) {
        h_gain->Fill(e_anode);
      }
      rgain = e_anode;
      cout << " e-  in anode " << e_anode << endl;

      retree->Fill();
      ritree->Fill();

      rx1.clear();
      rx2.clear();
      ry1.clear();
      ry2.clear();
      rz1.clear();
      rz2.clear();
      re1.clear();
      re2.clear();
      rt1.clear();
      rt2.clear();
      rixi.clear();
      riyi.clear();
      rizi.clear();
      rixf.clear();
      riyf.clear();
      rizf.clear();
      riti.clear();
      ritf.clear();

      //entries++;
    }
  } // event loop end

  TCanvas *cd = new TCanvas();
  if (plotDrift) {
    driftView->SetCanvas(cd);
    driftView->Plot("e");
  }

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<minutes>(stop - start);
  cout << "Time taken: " << duration.count() << " minutes" << '\n';

  // tree->Fill();
  cout << rntotalend << endl;

  // Write out the Tree File
  treefile->cd();
  retree->Write();
  ritree->Write();
  // sigtree->Write();
  cf1->Write();
  cp1->Write();
  cf1_0->Write();
  cp2->Write();
  cf2->Write();
  cf2_1->Write();
  cd->Write();
  hAspectRatio->Write();
  hSize->Write();
  h_gain->Write();

  treefile->Close();

  delete gas;
  delete sensor;
  delete aval;
  delete rdo_board;

  return 0;

  // app.Run(kTRUE);
}
