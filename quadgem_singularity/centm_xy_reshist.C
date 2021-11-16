void centm_xy_reshist()
{


  gROOT->SetBatch(kTRUE);
  TH1* finalhist_x = new TH1S("txrx", "truth x - recon x", 160, -0.02, 0.02);
  TH1* finalhist_y = new TH1S("tyry", "truth y - recon y", 160, -0.02, 0.02);

  Float_t kapton = 0.0050; // Thickness of the kapton layer, in cm
  Float_t metal = 0.0005;  // Thickness of the meta layers, in cm
  Int_t binct = 300;
  Int_t layer = 0;

  Float_t pitch = 0.0200;
  Float_t xstrip = 0.0150;
  Float_t ystrip = 0.0050;

  Float_t z_init = ((kapton + metal) * 1) + kapton + 0.0001;
  Float_t z_fin = z_init + metal;
  char cuttingX[100];
  sprintf(cuttingX, "rez2 > %f && rez2 < %f", z_init, z_fin);

  Float_t z_inity = kapton;
  Float_t z_finy = z_inity + metal;
  char cuttingY[100];
  sprintf(cuttingY, "rez2 > %f && rez2 < %f", z_inity, z_finy);
  //puts(cutting);

  cout << "histogram cut string: " << cuttingX << endl;
  TCut hist_cutX("metal layer", cuttingY);
  TCut hist_cutY("metal layer Y", cuttingX);


  bool iter_dir_x = false;
  bool iter_dir_45 = true;
  int file_ct = 6;

  TChain* tch_fln = new TChain("retree","chain of e- trees");
  tch_fln->AddFile("xyrdo_73ArCo2_350dv_50pi-_xy80i3u_xbfield15_1.root");
  int max_iter = tch_fln->GetEntries();

    tch_fln->Draw("rex0>>initer_elecsX", hist_cutX, "", 1, 0);
    tch_fln->Draw("rey0>>initer_elecsY", hist_cutY, "", 1, 0);
    
    TH1F *htempx = (TH1F*)gDirectory->FindObject("initer_elecsX");
    TH1F *htempy = (TH1F*)gDirectory->FindObject("initer_elecsY");

    
    double offsetXmax = htempx->GetXaxis()->GetBinCenter(htempx->GetMaximumBin());
    double offsetXmin = htempx->GetXaxis()->GetBinCenter(htempx->GetMinimumBin());
    
    double offsetYmax = htempy->GetXaxis()->GetBinCenter(htempy->GetMaximumBin());
    double offsetYmin = htempy->GetXaxis()->GetBinCenter(htempy->GetMinimumBin());
    
    double Yoffrange = (offsetYmax - offsetYmin) / 2;
    double Xoffrange = (offsetXmax - offsetXmin) / 2;
    cout << "X init will be from " << -Xoffrange << " to " << Xoffrange << endl; 
    cout << "Y init will be from " << -Yoffrange << " to " << Yoffrange << endl; 
    cout << "Rgx " << offsetXmin << " - " << offsetXmax << endl; 
    cout << "Rgy " << offsetYmin << " - " << offsetYmax << endl; 
    delete htempx;
    delete htempy;

    
    

  TFile *store_file = new TFile("xy_iterativexy_pi-offs_aggr.root","RECREATE");



  Double_t truth_cent;
  Double_t recon_cent;
  for(int fl = 1;fl <= file_ct;fl++){
  
    char filename_iter[50];
    
    // sprintf(filename_iter, "xyrdo_xite_quadgem350_samev_ncont_%d.root", file_ct);
    //sprintf(filename_iter, "xyrdo_xite_quadgem350_6mmtotal_20um_c%d.root", fl);
    sprintf(filename_iter, "xyrdo_73ArCo2_350dv_50pi-_xy80i3u_xbfield15_%d.root", fl);

    TFile *fin =  new TFile(filename_iter);
    fin->cd();

    TTree *out_etree = (TTree*)fin->Get("retree");
    //TTree etree;
    //etree=(TTree*)fin->Get("etree");
    //out_etree->Print();

    TGraph* recons_strx = new TGraph(max_iter);
    TGraph* recons_stry = new TGraph(max_iter);

    //Double_t iteration_distance = 0.003;
    Double_t iteration_distance = 0.00212;
    
    for(int cur_iter = 1;cur_iter <= max_iter;cur_iter++){
        int treeiter;
        Int_t iterationx;
        Int_t iterationy;
        if(iter_dir_45){
            
          iterationx = cur_iter;
          iterationy = cur_iter;
          treeiter = iterationx - 1;
          
        }else if (iter_dir_x)
        { 
          iterationx = cur_iter;
          iterationy = 0;
          treeiter = iterationx - 1;
        }
        else
        {
          iterationy = cur_iter;
          iterationx = 0;
          treeiter = iterationy - 1;
        }


        //middle of electron cloud = zero point
        // Double_t e_centerX_r = iteration_distance * iterationx;
        // Double_t e_centerY_r = iteration_distance * iterationy;
        // Double_t x_it_min = -0.12 + iteration_distance * iterationx;
        // Double_t x_it_max = 0.12 + iteration_distance * iterationx;
        // Double_t y_it_min = -0.12 + iteration_distance * iterationy;
        // Double_t y_it_max = 0.12 + iteration_distance * iterationy;

        Double_t x_it_min = offsetXmin + iteration_distance * iterationx;
        Double_t x_it_max = offsetXmax + iteration_distance * iterationx;
        Double_t y_it_min = offsetYmin + iteration_distance * iterationy;
        Double_t y_it_max = offsetYmax + iteration_distance * iterationy;
        
        //auto rexcanva = new TCanvas("rex", "rex");
        //rexcanva->cd();
        TH1 *h_allX = new TH1D("all_elecsX", "electrons X", binct, x_it_min, x_it_max);
        out_etree->Draw("rex0>>all_elecsX", hist_cutX, "", 1, treeiter);
        TF1 *fit_hX = new TF1("fit_h_x","gaus",x_it_min,x_it_max);
        fit_hX->SetName("fit_h_x");
        h_allX->Fit("fit_h_x","R");
        
        cout << "Iter: " << cur_iter << " - Fit Center X " << fit_hX->GetParameter(1) << endl;
        
        cout <<  x_it_min << " X Range " << x_it_max << endl;
        Double_t e_centerX = fit_hX->GetParameter(1);
        

        //auto reycanva = new TCanvas("rey", "rey");
        //reycanva->cd();
        TH1 *h_allY = new TH1D("all_elecsY", "electrons Y", binct, y_it_min, y_it_max);
        out_etree->Draw("rey0>>all_elecsY", hist_cutY, "", 1, treeiter);
        TF1 *fit_hY = new TF1("fit_h_y","gaus",y_it_min,y_it_max);
        fit_hY->SetName("fit_h_y");
        h_allY->Fit("fit_h_y","R");
        
        
        cout << "Iter: " << cur_iter << " - Fit Center Y " << fit_hY->GetParameter(1) << endl;
       
        cout <<  y_it_min << " Y Range " << y_it_max << endl;
       	Double_t e_centerY = fit_hY->GetParameter(1);
        
        
        //cout << "read histograms" << endl;

        Double_t initial_offset_x = pitch / 2 - xstrip / 2 - pitch * 2;
        Double_t xs_iter = initial_offset_x;
        Double_t xs_iter_con = 0.0;
        std::list<Double_t> contentlistX;
        std::list<Double_t> coordlistX;
        while (xs_iter < e_centerX || xs_iter_con > 0.0)
        {
          Double_t xs_iter_up = xs_iter + xstrip + 0.001;
          char cuttingStX[100];
          sprintf(cuttingStX, "rex2 > %f && rex2 < %f", xs_iter, xs_iter_up);
          //cout << cuttingStX << endl;

          TCut hist_cut_itX("X iterative cutter", cuttingStX);

          TH1 *x_temphist = new TH1S("xtemphist", "iterator histogram", 50, xs_iter, xs_iter_up);
          out_etree->Draw("rex2>>xtemphist", hist_cutX && hist_cut_itX, "colz", 1, treeiter);


          xs_iter_con = x_temphist->Integral();
          Double_t stripmean = xs_iter + (xstrip / 2);

          //cout << "strip at: " << stripmean << endl;
          //cout << "strip content: " << xs_iter_con << endl;
          if (xs_iter_con > 2.0)
          {
            contentlistX.push_back(xs_iter_con);
            coordlistX.push_back(stripmean);
          }
          xs_iter += pitch;
          delete x_temphist;
        }

        Double_t initial_offset_y = pitch / 2 - ystrip / 2 - pitch * 2;
        Double_t ys_iter = initial_offset_y;
        Double_t ys_iter_con = 0.0;
        std::list<Double_t> contentlistY;
        std::list<Double_t> coordlistY;
        while (ys_iter < e_centerY || ys_iter_con > 0.0)
        {
          Double_t ys_iter_up = ys_iter + ystrip + 0.001;
          char cuttingStY[100];
          sprintf(cuttingStY, "rey2 > %f && rey2 < %f", ys_iter, ys_iter_up);
          //cout << cuttingStY << endl;

          TCut hist_cut_itY("Y iterative cutter", cuttingStY);

          TH1 *y_temphist = new TH1S("ytemphist", "iterator histogram", 50, ys_iter, ys_iter_up);
          out_etree->Draw("rey2>>ytemphist", hist_cutY && hist_cut_itY, "colz", 1, treeiter);

          ys_iter_con = y_temphist->Integral();
          Double_t stripmean = ys_iter + (ystrip / 2);

          //cout << "strip at: " << stripmean << endl;
          //cout << "strip content: " << ys_iter_con << endl;
          if (ys_iter_con > 2.0)
          {
            contentlistY.push_back(ys_iter_con);
            coordlistY.push_back(stripmean);
          }
          ys_iter += pitch;
          delete y_temphist;
        }
        std::list<Double_t>::iterator ictx = contentlistX.begin();
        std::list<Double_t>::iterator icox = coordlistX.begin();

        Double_t xsum = 0.0;
        Double_t x_t_coord = 0.0;

        while (ictx != contentlistX.end())
        {
          xsum += *ictx;
          x_t_coord += *ictx * *icox;
          ictx++;
          icox++;
        }
        Double_t xcentroid = x_t_coord / xsum;

        std::list<Double_t>::iterator icty = contentlistY.begin();
        std::list<Double_t>::iterator icoy = coordlistY.begin();

        Double_t ysum = 0.0;
        Double_t y_t_coord = 0.0;

        while (icty != contentlistY.end())
        {
          ysum += *icty;
          y_t_coord += *icty * *icoy;
          icty++;
          icoy++;
        }
        Double_t ycentroid = y_t_coord / ysum;

        cout << "centroid x det " << xcentroid << endl;
        cout << "centroid y det " << ycentroid << endl;
        cout << "cloud center x " << e_centerX << endl;
        cout << "cloud center y " << e_centerY << endl;
        cout << xcentroid - e_centerX << " // " << ycentroid - e_centerY << endl;

        
        recons_strx->SetPoint(cur_iter-1,e_centerX,xcentroid);
        recons_stry->SetPoint(cur_iter-1,e_centerY,ycentroid);

        Double_t xvx = xcentroid - e_centerX;
        finalhist_x->Fill(xvx);
        Double_t yvy = ycentroid - e_centerY;
        finalhist_y->Fill(yvy);
        
        delete h_allX;
        delete h_allY;
        delete fit_hX;
        delete fit_hY;
      }


      char graphxname[20];
      char graphyname[20];
      sprintf(graphxname, "xre_xt_%d", fl);
      sprintf(graphyname, "yre_yt_%d", fl);

      store_file->cd();
      recons_strx->SetName(graphxname);
      recons_stry->SetName(graphyname);
      recons_strx->SetTitle("cloud reconstructed X vs actual X");
      recons_stry->SetTitle("cloud reconstructed Y vs actual Y");
      recons_strx->SetMarkerStyle(3);
      recons_stry->SetMarkerStyle(3);
      recons_strx->Write();
      recons_stry->Write();
      delete recons_strx;
      delete recons_stry;

    }
  gROOT->SetBatch(kFALSE);
  //TCanvas* c10 = new TCanvas("c10");
  //c10->cd();
  finalhist_x->SetMarkerStyle(20);
  finalhist_x->SetMarkerSize(1);
  finalhist_y->SetMarkerStyle(20);
  finalhist_y->SetMarkerSize(1);
  //finalhist_x->Draw("P");
  //c10->SaveAs("histby.png");
  //tch_fln->Draw("rey2:rex2","","colz",1,1);

  store_file->cd();
  finalhist_x->Write();
  finalhist_y->Write();
  store_file->Close();


  // cout << "contents X:" << endl;
  // for (Double_t a : contentlistX)
  // {
  //   cout << a << ",";
  // }
  // cout << endl;
  // cout << "coordinates X" << endl;
  // for (Double_t b : coordlistX)
  // {
  //   cout << b << ",";
  // }
  // cout << endl;

  // cout << "content Y:" << endl;
  // for (Double_t a : contentlistY)
  // {
  //   cout << a << ",";
  // }
  // cout << endl;
  // cout << "coordinates Y:" << endl;
  // for (Double_t b : coordlistY)
  // {
  //   cout << b << ",";
  // }

  //delete h_allX;
  //delete h_allY;
}
