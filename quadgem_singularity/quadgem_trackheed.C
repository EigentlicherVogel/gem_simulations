#define UNUSED(expr)  \
    do                \
    {                 \
        (void)(expr); \
    } while (0)
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
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/Random.hh"
#include <TChain.h>
#include "TObjArray.h"
#include <time.h>

#ifdef __MAKECINT__
#pragma link C++ class vector < float> + ;
#endif

using namespace Garfield;
using namespace std;
using namespace std::chrono;

auto start = high_resolution_clock::now();
/* Tree to store variables produced by Garfield++
     */
int main(int argc, char *argv[])
{
    UNUSED(argc);

    //TApplication app("app", &argc, argv);
    //plottingEngine.SetDefaultStyle();

    //ANSYS geometry variables
    Float_t kapton = 0.0050; //Thickness of the kapton layer, in cm
    Float_t metal = 0.0005;  // Thickness of the meta layers, in cm
    Float_t drift_gap = 0.3;
    Float_t transfer_gap = 0.2;
    Float_t induction_gap = 0.1;
    Float_t mesh = 1.5 * transfer_gap + 2 * (kapton + 2. * metal) + drift_gap;             //Coordinate of drift electrode/cathode along Z  in cm
    Float_t pad = -1.0 * (1.5 * transfer_gap + 2 * (kapton + 2. * metal) + induction_gap); // Coordinate of readout pad/ anode along Z in cm

    //Float_t drift_gap = drift-metal-kapton/2 // drift gap in cm
    //Float_t induction_gap =-1.0* (induct + metal + kapton/2) //induction gap in cm

    //Field applied in different regions (V/cm) and potential difference (V) across each surface of a gem
    Float_t induction_field = 2000.;
    Float_t transfer_field = 2000.;
    Float_t dv_gem = 350;

    // Define variables.
    Double_t xi, yi, zi, ei, ti;
    Int_t nele, nelend;

    int nevents = 50; // Number of events.
    //if (nevents<1) nevents=1000;

    TString filename;
    filename = Form("OUTFILE/quadgem_gemout_73ArCo2_350dv_50rxy01pi-_24f15bf_%s.root", argv[1]);
    TFile *treefile = new TFile(filename, "RECREATE");
    treefile->cd();

    TString etreeName, itreeName, ptreeName;
    Double_t x0, y0, z0, e0, t0;
    std::vector<float> x1, y1, z1, e1, t1;
    std::vector<float> x2, y2, z2, e2, t2;

    std::vector<float> ixf, iyf, izf, itf;
    std::vector<float> ixi, iyi, izi, iti;
    std::vector<float> stat;
    std::vector<float> es, is;

    //std::vector<float> pxi, pyi, pzi, pti;
    Int_t nclusters;
    Int_t neleprod;
    Double_t eloss;
    std::vector<int> neleprodc;
    std::vector<float> elossc;
    std::vector<float> pxc, pyc, pzc, ptc;

    Int_t entries, event;
    // Int_t entriespi, eventpi;
    Double_t gain;
    Int_t ntotalend = 0;

    etreeName = "etree";
    TTree *etree = new TTree(etreeName, "e- variables");
    cout << "Will write to File output ROOT file " << endl;

    // electron tree
    etree->Branch("entries", &entries, "entries/I");
    etree->Branch("event", &event, "event/I");
    etree->Branch("gain", &gain, "gain/D");
    etree->Branch("nele", &nele, "nele/I");
    etree->Branch("nelend", &nelend, "nelend/I");
    etree->Branch("ntotalend", &ntotalend, "ntotalend/I");

    //   etree->Branch("ez0",&z0, "ez0/D");
    //   etree->Branch("ex0",&x0, "ex0/D");
    //   etree->Branch("ey0",&y0, "ey0/D");
    //   etree->Branch("ee0",&e0, "ee0/D");
    //   etree->Branch("et0",&t0, "et0/D");

    etree->Branch("status", &stat);

    etree->Branch("ex1", &x1);
    etree->Branch("ey1", &y1);
    etree->Branch("ee1", &e1);
    etree->Branch("et1", &t1);
    etree->Branch("ez1", &z1);

    etree->Branch("ex2", &x2);
    etree->Branch("ey2", &y2);
    etree->Branch("ee2", &e2);
    etree->Branch("et2", &t2);
    etree->Branch("ez2", &z2);

    //Ion tree
    itreeName = "itree";
    TTree *itree = new TTree(itreeName, "ions variables");

    itree->Branch("entries", &entries, "entries/I");
    itree->Branch("event", &event, "event/I");
    itree->Branch("ixi", &ixi);
    itree->Branch("iyi", &iyi);
    itree->Branch("iti", &iti);
    itree->Branch("izi", &izi);
    itree->Branch("ixf", &ixf);
    itree->Branch("iyf", &iyf);
    itree->Branch("itf", &itf);
    itree->Branch("izf", &izf);

    ptreeName = "piontree";
    TTree *ptree = new TTree(ptreeName, "pions variables");

    // ptree->Branch("entriespi", &entriespi, "entriespi/I");
    // ptree->Branch("eventpi", &eventpi, "eventpi/I");

    ptree->Branch("z0", &z0, "z0/D");
    ptree->Branch("x0", &x0, "x0/D");
    ptree->Branch("y0", &y0, "y0/D");
    ptree->Branch("e0", &e0, "e0/D");
    ptree->Branch("t0", &t0, "t0/D");
    ptree->Branch("x0", &x0, "x0/D");
    ptree->Branch("nclusters", &nclusters, "nclusters/I");
    ptree->Branch("neleprod", &neleprod, "neleprod/I");

    // ptree->Branch("pxi", &pxi);
    // ptree->Branch("pyi", &pyi);
    // ptree->Branch("pti", &pti);
    // ptree->Branch("pzi", &pzi);
    ptree->Branch("pxc", &pxc);
    ptree->Branch("pyc", &pyc);
    ptree->Branch("ptc", &ptc);
    ptree->Branch("pzc", &pzc);
    ptree->Branch("eloss", &eloss, "eloss/D");
    ptree->Branch("elossc", &elossc);
    ptree->Branch("neleprodc", &neleprodc);

    //Gain histogram
    TH1D *h_gain = new TH1D("h_gain", "", 50, -0.1, 49.9);

    // Loading the Ansys files (weighting afterwards)
    ComponentAnsys123 *gem = new ComponentAnsys123();

    TString lsrs = "4gemf_350dv";
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

    gem->Initialise(elistfile, nlistfile, mplistfile, prnsolfile, "mm");
    gem->EnableMirrorPeriodicityX();
    gem->EnableMirrorPeriodicityY();
    gem->PrintRange();

    double ang = 0.0;    // Angle between E&B fields
    double bfield = 1.5; // Magnetic field
    cout << bfield << endl;
    double bx = bfield;
    double by = (bfield)*sin(ang * (22.0 / 7.0) / 180.0);
    double bz = (bfield)*cos(ang * (22.0 / 7.0) / 180.0);
    cout << bx << " " << by << " " << bz << endl;
    gem->SetMagneticField(bx, 0, 0);

    //=========================================================
    // Create histograms for aspect ratio and element size by retrieving the dimesions of each mesh elemennt. This will enable us to know about the quality of the mesh
    //=======================================================

    TCanvas *cp1 = new TCanvas("cp1", "", 600, 600);
    TCanvas *cp2 = new TCanvas("cp2", "", 600, 600);
    TCanvas *cf1 = new TCanvas("cf1", "", 600, 600);
    TCanvas *cf2 = new TCanvas("cf2", "", 600, 600);
    TCanvas *cf3 = new TCanvas("cf3", "", 600, 600);
    TCanvas *cfb = new TCanvas("cfb", "", 600, 600);
    const double pitch = 0.014; // in cm
    const bool plotField = false;
    if (plotField)
    {
        ViewField *fieldView = new ViewField();
        fieldView->SetComponent(gem);
        // Plot the potential along the hole axis.
        //TCanvas* cp = new TCanvas("cp", "", 600, 600);
        cp1->SetLeftMargin(0.16);
        fieldView->SetCanvas(cp1);
        fieldView->PlotProfile(0., 0., 0.02, 0., 0., -0.02);

        Float_t mkm = metal + kapton + metal;

        //gem4 coordinate + voltage
        Float_t g4_t = -1.0 * (mkm + 1.5 * transfer_gap - 0.02);        // 3rd gem top surface z coordinate offset by -0.3 mm
        Float_t g4_b = -1.0 * (mkm + 1.5 * transfer_gap + mkm + 0.02);  // 3rd gem bottom surface z coordinate offset by 0.3 mm
        Float_t g4_b_v = -1.0 * (induction_field * induction_gap - 50); // 3rd  gem bottom surface voltage offset by 50 V
        Float_t g4_t_v = g4_b_v - dv_gem - 50;
        cout << " g4b :" << g4_b << " " << g4_b_v << endl;
        cout << " g4t :" << g4_t << " " << g4_t_v << endl;

        //gem3 coordinate + voltage
        Float_t g3_t = -1.0 * (0.5 * transfer_gap - 0.02);                // 3rd gem top surface z coordinate offset by -0.3 mm
        Float_t g3_b = -1.0 * (0.5 * transfer_gap + mkm + 0.02);          // 3rd gem bottom surface z coordinate offset by 0.3 mm
        Float_t g3_b_v = -1.0 * (transfer_field * transfer_gap) + g4_t_v; // 3rd  gem bottom surface voltage offset by 50 V
        Float_t g3_t_v = g3_b_v - dv_gem - 50;
        cout << " g3b :" << g3_b << " " << g3_b_v << endl;
        cout << " g3t :" << g3_t << " " << g3_t_v << endl;

        //gem2 coordinate + voltage
        Float_t g2_t = (0.5 * transfer_gap + 0.02);                       // 2nd  gem top surface z coordinate offset by 0.3 mm
        Float_t g2_b = (0.5 * transfer_gap + mkm - 0.02);                 // 2nd gem bottom surface z coordinate offset by 0.3 mm
        Float_t g2_b_v = -1.0 * (transfer_field * transfer_gap) + g3_t_v; // 2nd  gem bottom surface voltage offset by 50 V
        Float_t g2_t_v = g2_b_v - dv_gem - 50;
        cout << " g2b :" << g2_b << " " << g2_b_v << endl;
        cout << " g2t :" << g2_t << " " << g2_t_v << endl;

        //Calculate top gem spatial coordinate in Z (offset by few mm) and voltages applied
        Float_t g1_b = mkm + 1.5 * transfer_gap - 0.02;       // Top gem bottom surface z coordinate offset by -0.3 mm
        Float_t g1_t = mkm + 1.5 * transfer_gap + mkm + 0.02; // Top gem top offset by +0.3 mm
        Float_t g1_b_v = -1.0 * (transfer_field * transfer_gap) + g2_t_v;
        Float_t g1_t_v = g1_b_v - dv_gem - 50;
        cout << " g1b :" << g1_b << " " << g1_b_v << endl;
        cout << " g1t :" << g1_t << " " << g1_t_v << endl;

        // Set the viewing plane (normal vector) and ranges
        fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
        //fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);
        fieldView->SetArea(-pitch / 2., g1_b, pitch / 2., g1_t);
        fieldView->SetVoltageRange(g1_t_v, g1_b_v); // gem example
        //TCanvas* cf = new TCanvas("cf", "", 600, 600);
        cf1->SetLeftMargin(0.16);
        fieldView->SetCanvas(cf1);
        fieldView->PlotContour();

        ViewField *fieldView2 = new ViewField();
        fieldView2->SetComponent(gem);
        // Plot the potential along the hole axis.
        cp2->SetLeftMargin(0.16);
        fieldView2->SetCanvas(cp2);
        fieldView2->PlotProfile(0., 0., 0.15, 0., 0., -0.15);

        // Set the viewing plane (normal vector) and ranges
        fieldView2->SetPlane(0., -1., 0., 0., 0., 0.);
        fieldView2->SetArea(-pitch / 2., g2_b, pitch / 2., g2_t);
        //fieldView2->SetArea(-5.*pitch, -0.11, 5*pitch , 0.1);
        fieldView2->SetVoltageRange(g2_t_v, g2_b_v);
        //TCanvas* cf = new TCanvas("cf", "", 600, 600);
        cf2->SetLeftMargin(0.16);
        fieldView2->SetCanvas(cf2);
        fieldView2->PlotContour();

        ViewField *fieldView3 = new ViewField();
        fieldView3->SetComponent(gem);
        fieldView3->SetPlane(0., -1., 0., 0., 0., 0.);
        fieldView3->SetArea(-pitch / 2., g3_b, pitch / 2., g3_t);
        //fieldView2->SetArea(-5.*pitch, -0.11, 5*pitch , 0.1);
        fieldView3->SetVoltageRange(g3_t_v, g3_b_v);
        //TCanvas* cf = new TCanvas("cf", "", 600, 600);
        cf3->SetLeftMargin(0.16);
        fieldView3->SetCanvas(cf3);
        fieldView3->PlotContour();

        ViewField *fieldView4 = new ViewField();
        fieldView4->SetComponent(gem);
        fieldView4->SetPlane(0., -1., 0., 0., 0., 0.);
        fieldView4->SetArea(-pitch / 2., g4_b, pitch / 2., g4_t);
        //fieldView2->SetArea(-5.*pitch, -0.11, 5*pitch , 0.1);
        fieldView4->SetVoltageRange(g4_t_v, g4_b_v);
        //TCanvas* cf = new TCanvas("cf", "", 600, 600);
        cf3->SetLeftMargin(0.16);
        fieldView4->SetCanvas(cf3);
        fieldView4->PlotContour();

        ViewField *fieldViewB = new ViewField();
        fieldViewB->SetComponent(gem);
        fieldViewB->SetPlane(0., -1., 0., 0., 0., 0.);
        //fieldView2->SetArea(-pitch / 2., -0.11, pitch / 2., 0.1);
        fieldViewB->SetArea(-4. * pitch, -4. * pitch, 4. * pitch, 4.0 * pitch);
        //fieldView2_1->SetVoltageRange(0., -200);
        //TCanvas* cf = new TCanvas("cf", "", 600, 600);
        cfb->SetLeftMargin(0.16);
        fieldViewB->SetCanvas(cfb);
        fieldViewB->PlotContour();
    }

    TH1F hElectrons("hElectrons", "Number of electrons", 200, 0, 200);
    TH1F hEdep("hEdep", "Energy Loss", 100, 0., 10.);

    //Gas setup
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
      const std::string path = std::getenv("GARFIELD_INSTALL");
      gas->LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
    gas->EnablePenningTransfer(0.55, 0.0, "ar");
    gas->Initialise(true);
    // gas->Initialise();
    int numofmat = gem->GetNumberOfMaterials();
    //cout << "Number of materials: " << numofmat << endl;

    //const std::string path = getenv("GARFIELD_HOME");

    for (int i = 0; i < numofmat; i++)
    {
        const double eps = gem->GetPermittivity(i);
        cout << "The permittivity of material " << i << " is " << eps << endl;
        if (fabs(eps - 1) < 1.e-13)
        {
            gem->SetMedium(i, gas);
        }
        //if (eps == 1.) gem->SetMedium(i, gas);
    }
    gem->PrintMaterials();

    Sensor *sensor = new Sensor();
    sensor->AddComponent(gem);
      sensor->SetArea(-24 * pitch, -48 * pitch, -0.6,
                  24 * pitch,  24 * pitch,  0.7);
    /*
      sensor->SetArea(-5 * pitch, -5 * pitch, -0.02,
              5 * pitch,  5 * pitch,   0.02  );
      */
      
          //pi-
    TrackHeed *track = new TrackHeed();
    track->SetSensor(sensor);
    track->SetParticle("pi-");
    constexpr double momentum = 300.e6;
    track->SetMomentum(momentum);
    track->EnableMagneticField();
    track->EnableElectricField();

    //e-
    AvalancheMicroscopic *aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);
    //aval->EnableDebugging();
    aval->EnableMagneticField();
    //aval->EnableAvalancheSizeLimit(1000);
    //aval->SetCollisionSteps(1000);
    //aval->EnableAvalancheSizeLimit(500);

    //ion
    AvalancheMC *drift = new AvalancheMC();
    drift->SetSensor(sensor);
    drift->SetDistanceSteps(2.e-4);

    // Loop over primary electrons
    //const double smear=pitch/2;



    const bool plotDrift = true;
    ViewDrift *driftView = new ViewDrift();
    if (plotDrift)
    {
        driftView->SetArea(-1.0 * pitch, -1.0 * pitch, -0.3,
                           1 * pitch, 1 * pitch, 0.3);
        track->EnablePlotting(driftView);
        aval->EnablePlotting(driftView);
    }

    entries = 0;
    for (int i = 0; i < nevents; i++)
    {
       cout << "pion event loop " << i << endl;
       track->EnableDebugging();
        event = i;
        //pxi.clear();
        pxc.clear();
        //pyi.clear();
        pyc.clear();
        //pzi.clear();
        pzc.clear();
        //pti.clear();
        ptc.clear();

        x0 = -0.5 * pitch + RndmUniform() * pitch; // [cm]
        y0 = -0.5 * pitch + RndmUniform() * pitch;
        z0 = 0.6 ;
        t0 = 0.0;
        double dx0 = 0.1 * (RndmUniform() - 0.5);
        double dy0 = 0.1 * (RndmUniform() - 0.5);
        double dz0 = -0.5;
        cout << x0 << " : " << y0 << ":" << z0 << endl;

        // const double smear = pitch/2;
        // xi = -smear + RndmUniform() * pitch;
        // yi = -smear + RndmUniform() * pitch;
        // zi = mesh - 0.05 ;
        // cout << "zi :" << zi << endl;
        // ei = 1.0;
        // ti = 0.0;
        // double dx0 = 0.0; double dy0 = 0.0; double dz0 = 0.0;

        xi = x0;
        yi = y0;
        zi = z0;
        ei = e0;
        ti = t0;
        track->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
        // Loop over the clusters.
        double xc = 0., yc = 0, zc = 0., tc = 0., ec = 0., extra;
        int nc = 0; // # of electrons produced in collision

        //Total energy loss along the track
        double esum = 0.;
        // Total number of electrons produced along the track
        int nsum = 0;

        // pxi.push_back(x0);
        // pyi.push_back(y0);
        // pzi.push_back(z0);
        // pti.push_back(t0);
        int clcnt = 0;
        while (track->GetCluster(xc, yc, zc, tc, nc, ec, extra))
        {
            clcnt++;
            cout << "pion cluster loop" << clcnt << " of pion " << i << endl;
            cout << "coord: " << xc << " : " << yc << " : " << zc << endl;
            pxc.push_back(xc);
            pyc.push_back(yc);
            pzc.push_back(zc);
            ptc.push_back(tc);
            elossc.push_back(ec);
            neleprodc.push_back(nc);

            esum += ec;
            nsum += nc;


            x1.clear();
            x2.clear();
            y1.clear();
            y2.clear();
            z1.clear();
            z2.clear();
            e1.clear();
            e2.clear();
            t1.clear();
            t2.clear();

            ixi.clear();
            iyi.clear();
            izi.clear();
            ixf.clear();
            iyf.clear();
            izf.clear();
            iti.clear();
            itf.clear();

            //Loop over total electrons
            for (int j = 0; j < nc; ++j)
            {

            cout << "electron loop " << j << endl;
                double xe, ye, ze, te, ee, dxe, dye, dze;
                track->GetElectron(j, xe, ye, ze, te, ee, dxe, dye, dze);

                aval->AvalancheElectron(xe, ye, ze, te, 0.1, dxe, dye, dze);
                int ne, ni;
                //Get # of e-s and ions produced in avalanche
                aval->GetAvalancheSize(ne, ni);
                nele = ne;

                //Get the history of each avalanche e- (starting and end point of avalanche e-)

                double xa0, ya0, za0, ea0, ta0;
                double xa1, ya1, za1, ea1, ta1;
                double xi1, yi1, zi1, ti1;
                double xi2, yi2, zi2, ti2;
                int istat;
                nelend = aval->GetNumberOfElectronEndpoints();
                cout << "nelend :" << nelend << endl;

                Int_t e_anode = 0;

                for (int j = 0; j < nelend; j++)
                {

                    aval->GetElectronEndpoint(j, xa0, ya0, za0, ta0, ea0, xa1, ya1, za1, ta1, ea1, istat);
                    //Store starting point of avalanche e-s and primaries
                    x1.push_back(xa0);
                    y1.push_back(ya0);
                    z1.push_back(za0);
                    t1.push_back(ta0);
                    e1.push_back(ea0);

                    //Store end point of avalanche e-s and primaries
                    x2.push_back(xa1);
                    y2.push_back(ya1);
                    z2.push_back(za1);
                    t2.push_back(ta1);
                    e2.push_back(ea1);
                    //status of avalanche e-
                    stat.push_back(istat);

                    drift->DriftIon(xa0, ya0, za0, ta0);
                    drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, istat);
                    ixi.push_back(xi1);
                    iyi.push_back(yi1);
                    izi.push_back(zi1);
                    iti.push_back(ti1);
                    ixf.push_back(xi2);
                    iyf.push_back(yi2);
                    izf.push_back(zi2);
                    itf.push_back(ti2);

                    ntotalend++;
                    if (za1 < -1.0 * (kapton / 2 + metal + transfer_gap + kapton + 2. * metal + 0.04))
                        e_anode++; // 0.04 cm added to make sure that we are not counting electrons stuck at bottom of GEM


                } // e- loop
                gain = e_anode;
                cout << " e-  in anode " << e_anode << endl;
            }

            etree->Fill();
            itree->Fill();
        }

        eloss = esum;
        nclusters = nc;
        neleprod = nsum;
        ptree->Fill();
        hElectrons.Fill(nsum);
        hEdep.Fill(esum * 1.e-3);
    }

//     float timeused = (float)clock() / CLOCKS_PER_SEC;
//     if (timeused > 25000000.0)
//     {
//         cout << "Time limit reached, " << timeused << " sec" << endl;
//         break;
//     }
    entries++;
// } //event loop end

TCanvas *cd = new TCanvas();
if (plotDrift)
{
    driftView->SetCanvas(cd);
    driftView->Plot();
}

//sensor->ClearSignal();

auto stop = high_resolution_clock::now();
auto duration = duration_cast<minutes>(stop - start);
cout << "Time taken: " << duration.count() << " minutes" << '\n';

//tree->Fill();
cout << ntotalend << endl;

// Write out the Tree File
treefile->cd();
etree->Write();
itree->Write();
ptree->Write();
//sigtree->Write();
cp1->Write();
cf1->Write();
cp2->Write();
cf2->Write();
cf3->Write();
cfb->Write();
cd->Write();
h_gain->Write();

treefile->Close();

delete track;
delete gas;
delete sensor;
delete aval;
delete gem;

return 0;

//app.Run(kTRUE);
}
