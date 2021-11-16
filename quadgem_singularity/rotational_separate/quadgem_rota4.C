#define UNUSED(expr) do { (void)(expr); } while (0)
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <stdbool.h>
#include <chrono>
#include <stdio.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <algorithm>
#include <TH2F.h>
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualViewer3D.h>
#include <TFile.h>
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/Random.hh"
#include <TChain.h>
#include "TObjArray.h"
#include <time.h>
/*
  #ifdef __MAKECINT__
  #pragma link C++ class vector<float>+;
  #endif
*/
using namespace Garfield;
using namespace std;
using namespace std::chrono;

auto start = high_resolution_clock::now();
/* Tree to store variables produced by Garfield++
 */
int main(int argc, char * argv[]) {
  UNUSED(argc);

  //TApplication app("app", &argc, argv);
  //plottingEngine.SetDefaultStyle();
      
  //ANSYS geometry variables
  Float_t kapton = 0.0050;   //Thickness of the kapton layer, in cm
  Float_t metal =  0.0005;   // Thickness of the meta layers, in cm
  Float_t indgap =   0.1;    //Induction gap
  Float_t tr_gap = 0.2; 

  int gnevents = 1 ; // Number of events.


  if (gnevents<1) gnevents=1000;

  //Previous stage variables
  Double_t ex0,ey0,ez0,ee0;
  vector <float> *ex2 = 0 ;
  vector <float> *ey2 = 0 ;
  vector <float> *ez2 = 0 ;
  vector <float> *et2 = 0 ;
  vector <float> *ee2 = 0 ;
 
  //Current stage variables
  Double_t xi,yi,zi,ei,ti;
  Int_t rnele,rnelend;
  TString etreeName, itreeName;
  // Double_t rx0,ry0,rz0,re0,rt0;
  
  std::vector <float> rx1,ry1,rz1, re1, rt1;
  std::vector <float> rx2,ry2,rz2, re2, rt2;
  std::vector <float> rx0,ry0,rz0, re0, rt0;
  
  
  //Double_t rx2,ry2,rz2,re2,rt2;
  
  std::vector <float> rixf,riyf,rizf,ritf;
  std::vector <float> rixi,riyi,rizi,riti;
  std::vector <float> rstat;
  std::vector <float> es, is;
      
  Int_t entries, revent;
  Double_t gain;
  Int_t rntotalend = 0;

  //Get the history of each avalanche e- (starting and end point of avalanche e-)
        
  double xa0, ya0, za0, ea0, ta0;
  double xa1, ya1, za1, ea1, ta1;
  double xi1,yi1,zi1,ti1;
  double xi2,yi2,zi2,ti2;
  int istat;

  TString filename;
  //filename = Form("./OUTFILE/xyrdo_rdo45ite_quadgem350_6mmtotal_30um_%s.root", argv[1]);
  
//   filename = Form("./OUTFILE/chevrdo_rdoxite_quadgem350_expcld1_5_40um300c_%s.root", argv[1]);
   filename = Form("OUTFILE/quadgem_gemout_73ArCo2_350dv_itlay_4_%s.root", argv[1]);

  TFile *treefile = new TFile(filename, "RECREATE");
  treefile->cd();
  
    etreeName = "retree";
  
  TTree *retree = new TTree(etreeName,"e- variables");
  cout << "Will write to File output ROOT file " << endl;
      
  // electron tree
  retree->Branch("entries",&entries,"entries/I");
  retree->Branch("revent",&revent,"revent/I");
  retree->Branch("gain",&gain,"gain/D");
  retree->Branch("rnele",&rnele,"rnele/I");
  retree->Branch("rnelend",&rnelend,"rnelend/I");
  retree->Branch("rntotalend",&rntotalend, "rntotalend/I");
  /*  
  retree->Branch("rez0",&rz0, "rez0/D");
  retree->Branch("rex0",&rx0, "rex0/D");
  retree->Branch("rey0",&ry0, "rey0/D");
  retree->Branch("ree0",&re0, "ree0/D");
  retree->Branch("ret0",&rt0, "ret0/D");
  */    
      
  retree->Branch("rstatus",&rstat);
  /*   
       retree->Branch("rex1",&rx1);
       retree->Branch("rey1",&ry1);
       retree->Branch("ree1",&re1);
       retree->Branch("ret1",&rt1);
       retree->Branch("rez1",&rz1);
  */
  retree->Branch("rex0",&rx0);
  retree->Branch("rey0",&ry0);
  retree->Branch("ree0",&re0);
  retree->Branch("ret0",&rt0);
  retree->Branch("rez0",&rz0);
  
  retree->Branch("rex2",&rx2);
  retree->Branch("rey2",&ry2);
  retree->Branch("ree2",&re2);
  retree->Branch("ret2",&rt2);
  retree->Branch("rez2",&rz2);
   
  //Gain histogram
  TH1D *h_gain = new TH1D("h_gain","",50, -0.1, 49.9);
      
  // Loading the Ansys files (weighting afterwards)
      ComponentAnsys123* gem = new ComponentAnsys123();
     
        TString lsrs = "quadgem_bot"; 
	TString elistf, nlistf, mplistf, prnsolf;
      elistf = Form("./ANSYS_geo/ELIST_%s.lis", lsrs.Data());
      nlistf = Form("./ANSYS_geo/NLIST_%s.lis", lsrs.Data());
      mplistf = Form("./ANSYS_geo/MPLIST_%s.lis", lsrs.Data());
      prnsolf = Form("./ANSYS_geo/PRNSOL_%s.lis", lsrs.Data());
        cout << "Field file formats: " << elistf << endl;
      std::string elistfile, nlistfile, mplistfile, prnsolfile;
      elistfile = elistf;
      nlistfile = nlistf;
      mplistfile = mplistf;
      prnsolfile = prnsolf;
      
      
      
  gem->Initialise(elistfile,		  nlistfile,		  mplistfile,                  prnsolfile, "mm");
  gem->EnableMirrorPeriodicityY();
  //gem->EnablePeriodicityX();
  gem->EnableMirrorPeriodicityX();
  gem->PrintRange();
      
  //=========================================================
  // Create histograms for aspect ratio and element size by retrieving the dimesions of each mesh elemennt. This will enable us to know about the quality of the mesh
  //=======================================================
      double ang = 0.0; // Angle between E&B fields
      double bfield = 0.0; // Magnetic field
      cout<<bfield<<endl;
      double bx= bfield;
      double by=(bfield)*sin(ang*(22.0/7.0)/180.0);
      double bz=(bfield)*cos(ang*(22.0/7.0)/180.0);
      cout<<bx<<" "<<by<<" "<<bz<<endl;
      gem->SetMagneticField(bx,0,0);
      
      
  TCanvas* cp1 = new TCanvas("cp1", "", 600, 600);
  TCanvas* cf1 = new TCanvas("cf1", "", 600, 600);
  TCanvas* cf1_0 = new TCanvas("cf1_0", "", 600, 600);
  TCanvas* cf1_1 = new TCanvas("cf1_1", "", 600, 600);
  TCanvas* cp2 = new TCanvas("cp2", "", 600, 600);
  TCanvas* cf2 = new TCanvas("cf2", "", 600, 600);
  TCanvas* cf2_1 = new TCanvas("cf2_1", "", 600, 600);
  const double pitch=0.0200; // in cm
  const bool plotField = false;
  
  if (plotField) {
    ViewField* fieldView = new ViewField();
    fieldView->SetComponent(gem);
    // Plot the potential along the hole axis.
    //TCanvas* cp = new TCanvas("cp", "", 600, 600);
    cp1->SetLeftMargin(0.16);
    fieldView->SetCanvas(cp1);
    fieldView->PlotProfile(0., 0., 0.02, 0., 0., -0.02);
        
       
    //fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
    fieldView->SetPlaneXZ();
    //fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);
    fieldView->SetArea(-2.*pitch , 0., 2.*pitch , 0.08);
    fieldView->SetElectricFieldRange(0., 2000); // gem example
    //TCanvas* cf = new TCanvas("cf", "", 600, 600);
    cf1->SetLeftMargin(0.16);
    fieldView->SetCanvas(cf1);
    fieldView->PlotContour("e");

    ViewField* fieldView0= new ViewField();
    fieldView0->SetComponent(gem);
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
    fieldView1->SetComponent(gem);
    // Set the viewing plane (normal vector) and ranges                                                         
    //fieldView1->SetPlane(0., 0., -1., 0., 0., 0.005);
    fieldView->SetPlaneXY();
    //fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);                                                   
    fieldView1->SetArea(-5.*pitch , -5.*pitch, 5.*pitch , 5.0*pitch);
    //fieldView1->SetElectricFieldRange(0., 2000.);                                                          
    //TCanvas* cf = new TCanvas("cf", "", 600, 600);    
    cf1_1->SetLeftMargin(0.16);
    fieldView1->SetCanvas(cf1_1);
    fieldView1->PlotContour("e");



    ViewField* fieldView2 = new ViewField();
    fieldView2->SetComponent(gem);
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
    fieldView2_1->SetComponent(gem);
    fieldView2_1->SetPlane(0., 0, -1., 0., 0., 0.);
    //fieldView2->SetArea(-pitch / 2., -0.11, pitch / 2., 0.1);                                                       
    fieldView2_1->SetArea(-4.*pitch , -4.*pitch, 4.*pitch , 4.0*pitch);
    //fieldView2_1->SetVoltageRange(0., -200);
    //TCanvas* cf = new TCanvas("cf", "", 600, 600);                                                                  
    cf2_1->SetLeftMargin(0.16);
    fieldView2_1->SetCanvas(cf2_1);
    fieldView2_1->PlotContour();

  }
      
  //Gas setup
  MediumMagboltz* gas=new MediumMagboltz();
  gas->SetComposition("ar",70.,"co2",30.);
  // gas->SetComposition("ar",45.,"co2",15.,"cf4",40.0);
  // gas->SetComposition("ar",80.,"cf4",20.);
  // gas->SetComposition("cf4",100.);
  // gas->SetComposition("ar",95.,"ch4",5.)
  gas->SetTemperature(293.5);
  gas->SetPressure(760.0);
  // gas->SetMaxElectronEnergy(200.);
  // gas->EnablePenningTransfer(0.51, 0.0, "Ne");
  gas->EnablePenningTransfer(0.55, 0.0, "ar");
  gas->Initialise(true);
  // gas->Initialise();
  int numofmat=gem->GetNumberOfMaterials();
  //cout << "Number of materials: " << numofmat << endl;
      
  //const std::string path = getenv("GARFIELD_HOME");
  //gas->LoadIonMobility("../Ion_Data/IonMobility_Ar+_Ar.txt");
  const std::string path = std::getenv("GARFIELD_INSTALL");
  gas->LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
  
  for(int i=0;i<numofmat;i++){
    const double eps = gem->GetPermittivity(i);
    cout<<"The permittivity of material "<<i<<" is "<< eps<<endl;
    if(fabs(eps-1)<1.e-13){gem->SetMedium(i,gas);}
    //if (eps == 1.) gem->SetMedium(i, gas);
  }
      
      
      
  Sensor* sensor = new Sensor();
  sensor->AddComponent(gem);
      
  /*
    sensor->SetArea(-5 * pitch, -5 * pitch, -0.02,
    5 * pitch,  5 * pitch,   0.02  );
  */
      
  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);
      
  //Set the ion transport
  AvalancheMC* drift = new AvalancheMC();
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
  aval->EnableMagneticField();
  // aval->EnableAvalancheSizeLimit(1000);
  // Loop over primary electrons
  //const double smear=pitch/2;

  //sprintf(filein, "../singlegem/OUTFILE/singlegem_4RO_dvgem400_dr1_ind2kvpcm_%s.root", argv[1]);
  //sprintf(filein, "../singlegem/OUTFILE/singlegem_4RO_einitmorespread_v2_dvgem400_dr1_ind2kvpcm_%s.root", argv[1]);
  TString fileinp;
  //fileinp = Form("./OUTFILE/quadgem_gemout_73ArCo2_350dv_xbfield15_IonF_%s.root", argv[1]);
  fileinp = Form("OUTFILE/quadgem_gemout_73ArCo2_350dv_itlay_3_%s.root", argv[1]);
  //
  //fileinp = "./OUTFILE/single_halfind.root";
  TFile *fin = new TFile(fileinp, "READ");
  fin->cd();

  //TTree *etree = (TTree*)fin->Get("ietree");
  TTree* etree;
  etree=(TTree*)fin->Get("etree");
    
  //From top gem 
  etree->SetBranchAddress("ex0", &ex0);
  etree->SetBranchAddress("ey0", &ey0);
  etree->SetBranchAddress("ez0", &ez0);
  etree->SetBranchAddress("ee0", &ee0);
    
  etree->SetBranchAddress("ex2", &ex2);
  etree->SetBranchAddress("ey2", &ey2);
  etree->SetBranchAddress("ez2", &ez2);
  etree->SetBranchAddress("ee2", &ee2);
  etree->SetBranchAddress("et2", &et2);
   

  gnevents = etree->GetEntries();
  treefile->cd();
  


  
  
  entries = 0;
  int it = 1;
  for(int l = 1; l <= it; l++){
     cout << "Starting iteration " << l << endl;

    
    for(int i=0; i < gnevents; i++) {
      
      etree->GetEntry(i);
      revent = i;
	
      Double_t e_drift;
      Double_t multip_cloudwidth = 1.0;
      Double_t iteration_movement = 0.0; //30u 45d
      Double_t const_y_offset = 0.0; //cm, to mid-strip
      //e_drift =  -1.0*(kapton/2. + metal + 0.005) ; // Single GEM
      //e_drift = -1.0*(kapton/2. + metal + tr_gap + kapton + 2.*metal + 0.005  ); // Triple GEM
     e_drift = -1.0*(tr_gap/2. + tr_gap + 2.*kapton + 4.*metal + 0.005  ); // 4GEM
      cout << " e- drift init z :" << e_drift << endl;
 
          for(int j = 0 ; j < ex2->size(); j++) {
               
            	//if((ez2->at(j) > -0.09)) continue; // For single GEM
            	if((ez2->at(j) > e_drift)) continue; // For triple GEM
            // 	xi = (ex2->at(j)*multip_cloudwidth)+iteration_movement*l;
            // 	yi = (ey2->at(j)*multip_cloudwidth)+const_y_offset;
            	xi = (ex2->at(j)*multip_cloudwidth)+iteration_movement*l-0.5*pitch;
            	yi = (ey2->at(j)*multip_cloudwidth)+iteration_movement*l;
            	//xi = ex2->at(j);
            	//yi = ey2->at(j)+0.1*l; // displace along Y
            	
            	// zi = (ez2->at(j))- ind+0.2; // This is 0.3 cm drift (0.1 mm previous stage + 0.2 mm shifting in z here
            	zi = (ez2->at(j)) - (ez2->at(j)) + 0.0950; // This is 0.1 cm drift (+0.1 cm to move the previous stage electron to R0 @z=0 and +0.1 cm to move them further up to drift region of RO geometry
            	ei = ee2->at(j);
            	ti = et2->at(j);
            	//	  cout << " z2land :" << ez2->at(j) << " input z ro :" << zi << " xi :" << ex2->at(j) << " yi :" << ey2->at(j) << endl;
            	double dx0 = 0.0; double dy0 = 0.0; double dz0 = 0.0;
            
                    
            	float timeused = (float)clock()/CLOCKS_PER_SEC;
            	if (timeused > 25000000.0) {
            	  cout << "Time limit reached, " << timeused << " sec" << endl;
            	  break;
            	}
            	//Start e- avalanche
            	aval->AvalancheElectron(xi,yi,zi,ti,ei,dx0,dy0,dz0);
            	int ne,ni;
            	//Get # of e-s and ions produced in avalanche
            	aval->GetAvalancheSize(ne,ni);
            	rnele = ne;
            	
                    
            	rnelend = aval->GetNumberOfElectronEndpoints();
            	//        cout << "nelend :" << rnelend << endl;
                    
            	//entries ++;
            	Int_t e_anode = 0 ;
            	for (int j= 0 ; j < rnelend; j++) {
                    
            	  aval->GetElectronEndpoint(j,xa0,ya0,za0,ta0,ea0,xa1,ya1,za1,ta1,ea1,istat);
                      
            	  rntotalend++       ;
            	  if(za1 < 0.05)e_anode++ ; // 0.04 cm added to make sure that we are not counting electrons stuck at bottom of GEM
            	  
                     
            	}// e- loop
            	if(e_anode != 0) {
            	  h_gain->Fill(e_anode);
                      
            	}
            	gain = e_anode;
            	//        cout << " e-  in anode " << e_anode << endl;
            
            	rx0.push_back(xi);
            	ry0.push_back(yi);
            	rz0.push_back(zi);
            	rt0.push_back(ti);
            	re0.push_back(ei);
            	
            	rx2.push_back(xa1);
            	ry2.push_back(ya1);
            	rz2.push_back(za1);
            	rt2.push_back(ta1);
            	re2.push_back(ea1);
            	//status of avalanche e-
            	rstat.push_back(istat);
                    
            
          }
      //entries++;
    }//event loop end
    
    retree->Fill();

    rx0.clear(); 
    ry0.clear();
    rz0.clear();
    re0.clear();
    rt0.clear();
    
    rx2.clear(); 
    ry2.clear();
    rz2.clear();
    re2.clear();
    rt2.clear();
    rstat.clear();
    entries++;
  } // iteration loop  end
  
  TCanvas* cd = new TCanvas();
  if (plotDrift) {
    driftView->SetCanvas(cd);
    driftView->Plot();
  }
  cd->Write();
      
      
      
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<minutes>(stop-start);
  cout << "Time taken: " << duration.count() << " minutes" << '\n';
      
  //tree->Fill();
  cout << rntotalend << endl;
      
  // Write out the Tree File

  treefile->cd();
  retree->Write();
  /*
  for(int m = 1; m <= it; m++){
    retree[m]->Write();
  }*/
  cp1->Write();
  cf1->Write();
  cf1_0->Write();
  cf1_1->Write();
  cp2->Write();
  cf2->Write();      
  cf2_1->Write();
  cd->Write();
  h_gain->Write();
      
     
  treefile->Close();
      
  delete gas;
  delete sensor;
  delete aval;
  delete gem;
      
      
  return 0;
      
  //app.Run(kTRUE);
      
}
