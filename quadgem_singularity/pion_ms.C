#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <stdlib.h>
#include <stdbool.h>

#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/Random.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  //TFile *fout = new TFile("../OUTFILE/gem_dep_piplus_womagfield_DCX1_drfitlineplot_v2_singlepion.root","RECREATE");
  //TFile *fout = new TFile("../OUTFILE/gem_dep_piplus_womagfield_DCX1_drfitlineplot_xzplane_v2_singlepion.root","RECREATE");
  //TFile *fout = new TFile("../OUTFILE/gem_dep_piplus_womagfield_DCX1_drfitlineplot_xzplane_v2_triplepion.root","RECREATE");
  //TFile *fout = new TFile("../OUTFILE/gem_dep_piplus_womagfield_DC_allzero_drfitlineplot_xzplane_v2_triplepion.root","RECREATE");
  //TFile *fout = new TFile("../OUTFILE/gem_dep_piplus_womagfield_DCX_o5_drfitlineplot_xzplane_v2_triplepion.root","RECREATE");
  // TFile *fout = new TFile("../OUTFILE/gem_dep_piminus_womagfield_DCZ_1_drfitlineplot_xzplane_v2_triplepion.root","RECREATE");
  //TFile *fout = new TFile("../OUTFILE/gem_dep_piminus_womagfield_DCZ_o5_drfitlineplot_xzplane_v2_triplepion.root","RECREATE");
  //TFile *fout = new TFile("../OUTFILE/gem_dep_300e6eVpiminus_womagfield_DCZ_neg_1_500pions.root","RECREATE");
  //TFile *fout = new TFile("../OUTFILE/gem_dep_300e6eVpiminus_wmagfield_DCZ_neg_1_500pions.root","RECREATE");
  //TFile *fout = new TFile("../OUTFILE/gem_dep_300e6eVpiminus_wmagfield_enableetrkfield_DCZ_neg_1_500pions.root","RECREATE");
  TFile *fout = new TFile("pion3benchmark.root","RECREATE");
  
  TCanvas* cf = new TCanvas("cf", "", 600, 600);
  TCanvas* cf1 = new TCanvas("cf1", "", 600, 600);
  TCanvas* cxy = new TCanvas("cxy", "", 600, 600);
  TCanvas* cxy1 = new TCanvas("cxy1", "", 600, 600);
  TCanvas* cd = new TCanvas("cd", "", 600, 600);

  // Histograms
  TH1::StatOverflows(true); 
  TH1F hElectrons("hElectrons", "Number of electrons", 200, 0, 200);
  TH1F hEdep("hEdep", "Energy Loss", 100, 0., 2.);

    
  // Load the field map.
  ComponentAnsys123 fm;
  //fm.Initialise("ELIST.lis", "NLIST.lis", "MPLIST.lis", "PRNSOL.lis", "mm");
    TString lsrs = "4gemf_350dv";
    TString elistf, nlistf, mplistf, prnsolf;
    elistf = Form("./ANSYS_geo/ELIST_%s.lis", lsrs.Data());
    nlistf = Form("./ANSYS_geo/NLIST_%s.lis", lsrs.Data());
    mplistf = Form("./ANSYS_geo/MPLIST_%s.lis", lsrs.Data());
    prnsolf = Form("./ANSYS_geo/PRNSOL_%s.lis", lsrs.Data());
    std::cout << "Field file formats: " << elistf << std::endl;
    std::string elistfile, nlistfile, mplistfile, prnsolfile;
    elistfile = elistf;
    nlistfile = nlistf;
    mplistfile = mplistf;
    prnsolfile = prnsolf;

    fm.Initialise(elistfile, nlistfile, mplistfile, prnsolfile, "mm");
  fm.EnableMirrorPeriodicityX();
  fm.EnableMirrorPeriodicityY();
  fm.PrintRange();

  // Dimensions of the GEM [cm]
  constexpr double pitch = 0.014;

  ViewField fieldView;
  ViewField fieldView1;
  ViewField fieldView2;
  ViewField fieldView3;

  
  constexpr bool plotField = true;
  if (plotField) {
    fieldView.SetComponent(&fm);
    // Set the normal vector of the viewing plane (xz plane).
    //fieldView.SetPlane(0, -1, 0, 0, 0, 0);
    fieldView.SetPlaneXZ();
    // Set the plot limits in the current viewing plane.
    //fieldView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02);
    fieldView.SetArea(-2. * pitch, -0.05, 2.0 * pitch, 0.05);
    fieldView.SetVoltageRange(-100., -900.);
    //TCanvas* cf = new TCanvas("cf", "", 600, 600);
    cf->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf);
    fieldView.PlotContour();

    fieldView1.SetComponent(&fm);
    //fieldView1.SetPlaneXZ();
    //fieldView1.SetArea(-2. * pitch, -0.05, 2.0 * pitch, 0.05);
    fieldView1.SetPlaneYZ();
    fieldView1.SetElectricFieldRange(0., 2000);
    cf1->SetLeftMargin(0.16);
    fieldView1.SetCanvas(cf1);
    fieldView1.PlotContour("e");

    fieldView2.SetComponent(&fm);
    fieldView2.SetPlaneXY();
    fieldView2.SetArea(-2.*pitch, 2.*pitch, -2.0*pitch, 2.*pitch);
    //fieldView1.SetVoltageRange(-160., 160.);
    cxy->SetLeftMargin(0.16);
    fieldView2.SetCanvas(cxy);
    fieldView2.PlotContour();

    fieldView3.SetComponent(&fm);
    fieldView3.SetPlaneXY();
    fieldView3.SetArea(-2.*pitch, 2.*pitch, -2.0*pitch, 2.*pitch);
    //fieldView1.SetElectricFieldRange(0., 20000);
    cxy1->SetLeftMargin(0.16);
    fieldView3.SetCanvas(cxy1);
    fieldView3.PlotContour("e");
    
    
  }

  // Setup the gas.
  MediumMagboltz gas;
  gas.SetComposition("ar", 80., "co2", 20.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.Initialise(true);  
  // Set the Penning transfer efficiency.
  constexpr double rPenning = 0.51;
  constexpr double lambdaPenning = 0.;
  gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  const std::string path = std::getenv("GARFIELD_INSTALL");
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
  // Associate the gas with the corresponding field map material. 
  const unsigned int nMaterials = fm.GetNumberOfMaterials();
  for (unsigned int i = 0; i < nMaterials; ++i) {
    const double eps = fm.GetPermittivity(i);
    if (eps == 1.) fm.SetMedium(i, &gas);
  }
  fm.PrintMaterials();
  /*
  // Set the magnetic field [Tesla].
  constexpr double bfield = 3.; 
  fm.SetMagneticField(bfield,  0, 0);
  */
  
  // Create the sensor.
  Sensor sensor;
  sensor.AddComponent(&fm);
  sensor.SetArea(-5 * pitch, -5 * pitch, -0.2,
                  5 * pitch,  5 * pitch,  0.3);

  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);
  aval.EnableMagneticField();
  aval.EnableAvalancheSizeLimit(1000);
  
  AvalancheMC drift;
  drift.SetSensor(&sensor);
  drift.SetDistanceSteps(2.e-4);

  // Simulate an ionizing particle (negative pion) using Heed.
  TrackHeed track;
  track.SetParticle("pi-");
  constexpr double momentum = 300.e6; // [eV / c]
  track.SetMomentum(momentum);
  track.SetSensor(&sensor);
  //track.EnableMagneticField();
  //track.EnableElectricField();
  
  ViewDrift driftView;
  constexpr bool plotDrift = true;
  if (plotDrift) {
    //aval.EnablePlotting(&driftView);
    //drift.EnablePlotting(&driftView);

    track.EnablePlotting(&driftView);
    aval.EnablePlotting(&driftView);
    aval.EnableExcitationMarkers(false);
    aval.EnableIonisationMarkers(false);
    
  }
  /*
  constexpr unsigned int nEvents = 10;
  for (unsigned int i = 0; i < nEvents; ++i) { 
    std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position. 
    const double x0 = -0.5 * pitch + RndmUniform() * pitch;
    const double y0 = -0.5 * pitch + RndmUniform() * pitch;
    const double z0 = 0.02; 
    const double t0 = 0.;
    const double e0 = 0.1;
    aval.AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    int ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);
    const unsigned int np = aval.GetNumberOfElectronEndpoints();
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1;
    double xi2, yi2, zi2, ti2;
    int status;
    for (unsigned int j = 0; j < np; ++j) {
      aval.GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
                                  xe2, ye2, ze2, te2, e2, status);
      drift.DriftIon(xe1, ye1, ze1, te1);
      drift.GetIonEndpoint(0, xi1, yi1, zi1, ti1, 
                              xi2, yi2, zi2, ti2, status);
    }
  }
  */

  const int nEvents = 2;
  track.EnableDebugging();

  for (int i = 0; i < nEvents; ++i) {
    if (i == 1) track.DisableDebugging();
    if (i % 1 == 0) std::cout << i << "/" << nEvents << "\n";

  // Set the starting point of the incoming ionising track.
    //For single pion to visualize
    //double x0 =  0.; // [cm]
    //double y0 =  0.;
    //For higher statistics randon distribution of multiple pions
  double x0 =  -0.5 * pitch + RndmUniform() * pitch; // [cm]
  double y0 =  -0.5 * pitch + RndmUniform() * pitch;
  double z0 = 0.29;
  double t0 = 0.0;
  double dx0 = 0.0;
  double dy0 = 0.0;
  double dz0 = -1.;

  // Now simulate a track.
  track.NewTrack(x0, y0, z0, t0, dx0, dy0 ,dz0 );
  // Loop over the clusters.
  double xc =0., yc =0, zc =0., tc=0., ec=0., extra;
  int nc =0; // # of electrons produced in collision
  //Total energy loss along the track
  double esum = 0 ;
  // Total number of electrons produced along the track
    int nsum = 0;
    //Loop over clusters
  while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    esum += ec;
    nsum += nc;
    //Loop over total electrons 
    for (int j = 0; j < nc; ++j) { 
      double xe, ye, ze, te, ee, dxe, dye, dze;
      track.GetElectron(j, xe, ye, ze, te, ee, dxe, dye, dze);
      // Simulate the drift/avalanche of this electron.
      aval.AvalancheElectron(xe, ye, ze, te, 0.1, dxe, dye, dze);
    }
  }
  hElectrons.Fill(nsum);
  hEdep.Fill(esum * 1.e-3);
  }

  if (plotDrift) {
    //TCanvas* cd = new TCanvas();
    constexpr bool plotMesh = false;
    if (plotMesh) {
      ViewFEMesh* meshView = new ViewFEMesh();
      
      meshView->SetCanvas(cd);
      meshView->SetComponent(&fm);
      // x-z projection.
      //meshView->SetPlane(0, -1, 0, 0, 0, 0);
      //meshView->SetArea(-2 * pitch, -2 * pitch, -0.2, 2 * pitch,  2 * pitch, 0.3);
      //y-z
      meshView->SetPlane(1, 0, 0, 0, 0, 0);
      meshView->SetArea(-5 * pitch, -5 * pitch, -0.2, 5 * pitch,  5 * pitch, 0.3);
      meshView->SetFillMesh(true);
      meshView->SetColor(0, kGray);
      // Set the color of the kapton.
      meshView->SetColor(2, kYellow + 3);
      meshView->EnableAxes();
      meshView->SetViewDrift(&driftView);
      meshView->Plot();
    } else {
      /*
      //For xz plane
      
      driftView.SetPlane(0, -1, 0, 0, 0, 0); //when dc z0 = -1
      //driftView.SetPlane(0, 1, 0, 0, 0, 0);   // when dc z0 = +1
      driftView.SetArea(-2 * pitch, -0.2, 2 * pitch, 0.3);
      */
      //For yz plane
      //driftView.SetPlane(1, 0, 0, 0, 0, 0);
      driftView.SetPlaneYZ();
      //driftView.SetArea(-0.2, -5 * pitch ,0.3,  5 * pitch);
      
      driftView.SetCanvas(cd);
      constexpr bool twod = true;
      driftView.Plot(twod);
    }
  }
  
  fout->cd();
  cf->Write();
  cf1->Write();
  cxy->Write();
  cxy1->Write();
  hElectrons.Write();
  hEdep.Write();
  
  cd->Write();
  fout->Close();
  /*
  //cd->SaveAs("driftplot_gem_dep_piplus_womagfield_DCX1_drfitlineplot_xzplane_v2_triplepion.pdf");
  //cd->SaveAs("driftplot_gem_dep_piplus_womagfield_DCX_o5_drfitlineplot_xzplane_v2_triplepion.pdf");
  //cd->SaveAs("driftplot_gem_dep_piplus_womagfield_DCX_o5_drfitlineplot_yzplane_v2_triplepion.pdf");
  //cd->SaveAs("driftplot_meshview_gem_dep_piplus_womagfield_DCX_o5_drfitlineplot_yzplane_v2_triplepion.pdf");
  //cd->SaveAs("driftplot_gem_dep_piminus_womagfield_DCZ_1_drfitlineplot_xzplane_v2_triplepion.pdf");
  //cd->SaveAs("gem_dep_piminus_womagfield_DCZ_o5_drfitlineplot_xzplane_v2_triplepion.pdf");
  //cd->SaveAs("driftplot_meshview_gem_dep_piminus_womagfield_DCZ_o5_drfitlineplot_xzplane_v2_triplepion.pdf");
  //cd->SaveAs("driftplot_gem_dep_piminus_womagfield_DCZ_o5_drfitlineplot_yzplane_v2_triplepion.pdf");
  //cd->SaveAs("driftplot_gem_dep_piminus_womagfield_DCZ_+1_drfitlineplot_yzplane_v2_triplepion.pdf");
  //cd->SaveAs("driftplot_gem_dep_piminus_womagfield_DCZ_+1_drfitlineplot_xzplane_dc_-1_triplepion_it2.pdf");
  //cd->SaveAs("driftplot_gem_dep_piminus_womagfield_DCZ_-1_drfitlineplot_xzplane_dc_-1_triplepion.pdf");
  //cd->SaveAs("driftplot_gem_dep_piminus_womagfield_DCZ_+1_drfitlineplot_xzplane_dc_+1_triplepion.pdf");
  //cd->SaveAs("../OUTFILE/driftplot_gem_dep_piminus_w3Tmagfield_DCZ_+1_drfitlineplot_xzplane_dc_+1_triplepion.pdf");
  //cd->SaveAs("../OUTFILE/driftplot_gem_dep_piminus_w3Tmagfield_trkefieldenable_DCZ_+1_drfitlineplot_xzplane_dc_+1_triplepion.pdf");
  fout->cd();
  cf->Write();
  cf1->Write();
  cxy->Write();
  cxy1->Write();
  hElectrons.Write();
  hEdep.Write();
  
  //cd->Write();
  //fout->Close();
  */
  app.Run(kTRUE);
  return 0;
  //app.Run(kTRUE);

}
