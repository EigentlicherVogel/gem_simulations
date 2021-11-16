void centm_chev_reshist()
{

  gROOT->SetBatch(kTRUE);
  TH1 *finalhist_x = new TH1S("txrx", "truth x - recon x", 400, -0.5, 0.8);
  TH1 *finalhist_y = new TH1S("tyry", "truth y - recon y", 300, -0.5, 0.5);

  Double_t kapton = 0.0050; // Thickness of the kapton layer, in cm
  Double_t metal = 0.0005;  // Thickness of the meta layers, in cm
  Int_t binct = 160;
  Int_t layer = 0;

  //Parameters for geometry, in cm
  Double_t geocel_half_y = 0.0250;
  Double_t geocel_x = 0.3252;
  Double_t geoperpdis = 0.0107;
  Double_t geoperpwid = 0.0249;
  Double_t geoapxwid = 0.2275;
  Double_t geoapxgap = 0.0977;
  Double_t geoapxpitch = geoapxwid + geoapxgap;
  Double_t charangle = 6.285; //degrees
  Double_t charanglerad = charangle * (M_PI / 180.0);

  Double_t yxslope = tan(charanglerad);
  cout << charanglerad << " < angle , tan > " << yxslope << endl;

  Double_t geoperpgap = geoapxgap * yxslope;
  Double_t geoperpwidth = geoapxwid * yxslope;
  Double_t geoperppitch = geoapxpitch * yxslope;
  Double_t geoperprestw = geocel_half_y - geoperpgap;

  Double_t yxcos = cos(charanglerad);
  Double_t geoxgap = geoperpwid * yxcos;
  Double_t dist_xopposcorner = geocel_x + geoxgap;

  Double_t z_init = kapton + 0.0001;
  Double_t z_fin = z_init + metal;
  char cuttingX[100];
  sprintf(cuttingX, "rez2 > %f && rez2 < %f", z_init, z_fin);
  
  //cout << "histogram cut string: " << cuttingX << endl;
  TCut hist_cutX("metal layer",cuttingX);

  bool iter_dir_x = true;

  int file_ct = 7;

  TChain *tch_fln = new TChain("retree", "chain of e- trees");
  tch_fln->AddFile("chevrdo_xite_quadgem350_6mmtotal_20um_1.root");
  int max_iter = tch_fln->GetEntries();

  TFile *store_file = new TFile("chev_iterative_aggr_ynr.root", "RECREATE");

  Double_t truth_cent;
  Double_t recon_cent;
  for (int fl = 1; fl <= file_ct; fl++)
  {

    char filename_iter[50];
      cout << " file starting: " << filename_iter << endl;

    // sprintf(filename_iter, "xyrdo_xite_quadgem350_samev_ncont_%d.root", file_ct);
    //sprintf(filename_iter, "chevrdo_xite_quadgem350_expcld1_5_40um300c_m%d.root", fl);
    
    sprintf(filename_iter, "chevrdo_xite_quadgem350_6mmtotal_20um_%d.root", fl);

    TFile *fin = new TFile(filename_iter);
    fin->cd();

    TTree *out_etree = (TTree *)fin->Get("retree");
    //TTree etree;
    //etree=(TTree*)fin->Get("etree");
    //out_etree->Print();

    TGraph *recons_strx = new TGraph(max_iter);
    TGraph* recons_stry = new TGraph(max_iter);

    Double_t iteration_distance = 0.002;

    //Double_t iteration_distance = 0.004;
    for (int cur_iter = 1; cur_iter <= max_iter; cur_iter++)
    {
      int treeiter;
      Int_t iterationx;
      Int_t iterationy;
      if (iter_dir_x)
      {
        iterationx = cur_iter;
        iterationy = 1;
        treeiter = iterationx - 1;
      }
      else
      {
        iterationy = cur_iter;
        iterationx = 1;
        treeiter = iterationy - 1;
      }
      cout << " iteration x: " << iterationx << endl;
      cout << " iteration y: " << iterationy << endl;

      //middle of electron cloud = zero point
        Double_t e_centerX_r = iteration_distance * iterationy;
        Double_t e_centerY_r = iteration_distance * iterationx;

        Double_t x_it_min = -0.12 + iteration_distance * iterationx;
        Double_t x_it_max = 0.12 + iteration_distance * iterationx;
        Double_t y_it_min = -0.12 + iteration_distance * iterationy;
        Double_t y_it_max = 0.12 + iteration_distance * iterationy;
        
        // auto rexcanva = new TCanvas("rex", "rex");
        // rexcanva->cd();
        TH1 *h_allX = new TH1D("all_elecsX", "electrons X", binct, x_it_min, x_it_max);
        out_etree->Draw("rex0>>all_elecsX", hist_cutX, "", 1, treeiter);
        TF1 *fit_hX = new TF1("fit_h_x","gaus",x_it_min,x_it_max);
        fit_hX->SetName("fit_h_x");
        h_allX->Fit("fit_h_x","R");
        
        
        Double_t e_centerY = fit_hX->GetParameter(1);
        

        // auto reycanva = new TCanvas("rey", "rey");
        // reycanva->cd();
        TH1 *h_allY = new TH1D("all_elecsY", "electrons Y", binct, y_it_min, y_it_max);
        out_etree->Draw("rey0>>all_elecsY", hist_cutX, "", 1, treeiter);
        TF1 *fit_hY = new TF1("fit_h_y","gaus",y_it_min,y_it_max);
        fit_hY->SetName("fit_h_y");
        h_allY->Fit("fit_h_y","R");
        
        
        Double_t e_centerX = fit_hY->GetParameter(1);
      // auto rexcanva = new TCanvas("rex", "rex");
      // rexcanva->cd();
      // TH2 *h_allX = new TH2D("all_elecsX", "electrons X", binct, x_it_min, x_it_max, binct, y_it_min, y_it_max);
      // out_etree->Draw("rex2:rey2>>all_elecsX", hist_cutX, "colz", 1, treeiter);

      // auto reycanva = new TCanvas("rey", "rey");
      // reycanva->cd();
      // TH2 *h_allY = new TH2D("all_elecsY", "electrons Y", binct, x_it_min, x_it_max, binct, y_it_min, y_it_max);
      // out_etree->Draw("rex2:rey2>>all_elecsY", hist_cutY, "colz", 1, treeiter);
      //cout << "read histograms" << endl;

      bool direction_ltor = true;
      int y_dir_width = 3;
      int x_dir_width = 3;
      bool geo_up = false;
      if (y_dir_width % 2 == 0)
      {
        geo_up = true;
      }
      // std::list<Double_t> contentlist;
      // std::list<Double_t> coordlist;

      std::map<Double_t, Double_t> con_coord_map;
      Double_t ys_iter_init = -1 * y_dir_width * geocel_half_y;
      Double_t ys_iter = ys_iter_init;
      Double_t ys_iter_con = 0.0;
      Int_t initgeo = -y_dir_width;

      int n_graphs = 0;

      // TF1 *cutfnu;
      // TF1 *cutfnd;

      while (ys_iter < e_centerY || ys_iter_con > 0.0)
      {
        Double_t ys_iter_up = ys_iter + geocel_half_y;
        char cuttingStY[100];

        sprintf(cuttingStY, "rey2 > %f && rey2 < %f", ys_iter, ys_iter_up);
        // cout << "Y geometry number: " << initgeo << endl;
        // cout << cuttingStY << endl;

        TCut hist_cut_itY("Y iterative cutter", cuttingStY);
        TH1 *y_temphist = new TH1S("ytemphist", "iterator histogram", 50, ys_iter, ys_iter_up);
       out_etree->Draw("rey2>>ytemphist", hist_cutX && hist_cut_itY, "colz", 1, treeiter);
        ys_iter_con = y_temphist->Integral();

        //cout << "row sum sh: " << ys_iter_con << endl;

        Double_t xlimiter_iter = -1 * x_dir_width * geocel_x;
        Double_t xs_iter_con = 0.0;
        Double_t e_stopX = (x_dir_width * geocel_x) + e_centerX;

        while (xlimiter_iter < e_stopX || xs_iter_con > 0.0)
        {

          Double_t xlimiter_u = xlimiter_iter + geoapxwid + 0.0001;
          Double_t xmid = xlimiter_iter + (dist_xopposcorner / 2);

          Double_t y_const_d;
          Double_t y_const_u;

          TH2 *h_cuttr = new TH2D("cut_line_hist", "cutting histogram", 50, xlimiter_iter, xlimiter_u, 50, ys_iter, ys_iter_up);
          //reycanva->cd();
          char cuttingLN1[200];
          if (geo_up)
          {
            y_const_d = ys_iter - (xlimiter_iter * yxslope); //y >/< mx + [b]
            y_const_u = ys_iter - (xlimiter_u * yxslope);    //[b] = y - mx
            sprintf(cuttingLN1, "rez2 > %f && rez2 < %f "
                                " && rey2 > %f && rey2 < %f "
                                " && rey2 < %f * rex2 + %f "
                                " && rey2 > %f * rex2 + %f ",
                    z_init, z_fin, ys_iter, ys_iter_up, yxslope, y_const_d, yxslope, y_const_u);

            // cutfnu = new TF1("upper_cut", "[0]*x+[1]", xlimiter_iter, xlimiter_u);
            // cutfnd = new TF1("lower_cut", "[0]*x+[1]", xlimiter_iter, xlimiter_u);

            // cutfnu->SetLineColor(4);
            // cutfnd->SetLineColor(4);
            // cutfnu->SetLineStyle(2);

            // cutfnu->SetParameters(yxslope, y_const_u);
            // cutfnd->SetParameters(yxslope, y_const_d);
            // cutfnu->SetRange(x_it_min, ys_iter, x_it_max, ys_iter_up);
            // cutfnd->SetRange(x_it_min, ys_iter, x_it_max, ys_iter_up);
          }
          else
          {
            y_const_u = ys_iter + geocel_half_y - (xlimiter_iter * -yxslope);
            y_const_d = ys_iter + geocel_half_y - (xlimiter_u * -yxslope);
            sprintf(cuttingLN1, "rez2 > %f && rez2 < %f "
                                " && rey2 > %f && rey2 < %f "
                                " && rey2 < %f * rex2 + %f "
                                " && rey2 > %f * rex2 + %f ",
                    z_init, z_fin, ys_iter, ys_iter_up, -yxslope, y_const_d, -yxslope, y_const_u);
            // cutfnu = new TF1("upper_cut", "[0]*x+[1]", xlimiter_iter, xlimiter_u);
            // cutfnd = new TF1("lower_cut", "[0]*x+[1]", xlimiter_iter, xlimiter_u);

            // cutfnu->SetLineColor(2);
            // cutfnd->SetLineColor(2);
            // cutfnu->SetLineStyle(2);

            // cutfnu->SetParameters(-yxslope, y_const_d);
            // cutfnd->SetParameters(-yxslope, y_const_u);
            // cutfnu->SetRange(x_it_min, ys_iter, x_it_max, ys_iter_up);
            // cutfnd->SetRange(x_it_min, ys_iter, x_it_max, ys_iter_up);
          }
          TCut hist_cutL("limiter1", cuttingLN1);
         out_etree->Draw("rey2:rex2>>cut_line_hist", hist_cutL, "colz", 1, treeiter);

          Double_t integz = h_cuttr->ProjectionY()->Integral();
          xs_iter_con = integz;
          // TObject *clonefnu = cutfnu->Clone();
          // TObject *clonefnd = cutfnd->Clone();

          if (integz > 10.0)
          {
            //cout << "cut content at: " << cuttingLN1 << endl;
            //cout << integz << endl;
            // fncanva->cd();
            // clonefnu->Draw("same");
            // clonefnd->Draw("same");

            // if (geo_up)
            // {
            //   cout << "up" << endl;
            // }
            // else
            // {
            //   cout << "down" << endl;
            // }
            // // reycanva->cd();
            // cout << "x between: " << xlimiter_iter << " and  " << xlimiter_u << endl;

            // contentlist.push_back(integz);
            // coordlist.push_back(xmid);
            if (con_coord_map.count(xmid) == 0)
            {
              con_coord_map.insert({xmid, integz});
            }
            else
            {
              con_coord_map[xmid] += integz;
            }

            //  n_graphs++;
            //  char graphname[30];
            //  sprintf(graphname, "debugraph%d.png" , n_graphs);
            //  reycanva->SaveAs(graphname);
          }
          xlimiter_iter += geocel_x;
          delete h_cuttr;
        }
        //ys_iter_con = xs_iter_con;

        // TH1* y_temphist = new TH1S("ytemphist","iterator histogram",50,ys_iter, ys_iter_up);
        // out_etree->Draw("rex2>>ytemphist",hist_cutX&&hist_cut_itY,"colz");
        // ys_iter_con = y_temphist->Integral();
        // cout << y_temphist->Integral() << endl;

        ys_iter += geocel_half_y;
        geo_up = !geo_up;
        initgeo++;
        delete y_temphist;
      }
      Double_t xsum = 0.0;
      Double_t x_t_coord = 0.0;
      cout << " coordinates : content " << endl;
      for (auto const &x : con_coord_map)
      {
        cout << x.first << ':' << x.second << endl;
        xsum += x.second;
        x_t_coord += x.first * x.second;
      }
      Double_t xcentroid = x_t_coord / xsum;
      cout << "centroid x det " << xcentroid << endl;
      cout << "cloud center x " << e_centerX << endl;
      cout << xcentroid - e_centerX << " // " << endl;

      recons_strx->SetPoint(cur_iter - 1, e_centerY, xcentroid);
      recons_stry->SetPoint(cur_iter-1,e_centerX, xcentroid);
      Double_t xvx = xcentroid - e_centerX;
      Double_t yvy = e_centerY - xcentroid;
      finalhist_x->Fill(yvy);
      finalhist_y->Fill(xvx);
      
        delete h_allX;
        delete h_allY;
        delete fit_hX;
        delete fit_hY;
      //delete y_temphist;
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

 

    // }
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


}
