void centm_udxy_reshist()
{


  gROOT->SetBatch(kTRUE);
  
  //storage of final results
  TH1* finalhist_u = new TH1S("tur", "truth u - recon u", 160, -0.01, 0.01);
  TH1* finalhist_l = new TH1S("tlr", "truth l - recon l", 160, -0.1, 0.1);
    
  Double_t truth_cent;
  Double_t recon_cent;
  

  Float_t kapton = 0.0050; // Thickness of the kapton layer, in cm
  Float_t metal = 0.0005;  // Thickness of the meta layers, in cm
  Int_t binct = 320;
  Int_t layer = 0;

  Float_t pitch = 0.0200;
  Float_t dstrip = 0.0150;
  Float_t ustrip = 0.0050;

  //Cut for upper strips
  Float_t z_init = ((kapton + metal) * 1) + kapton + 0.0001; 
  Float_t z_fin = z_init + metal;
  char cuttingU[100]; //Convert to string
  sprintf(cuttingU, "rez2 > %f && rez2 < %f", z_init, z_fin);

  //Cut for lower strips
  Float_t z_init1 = kapton + 0.0001;
  Float_t z_fin1 = z_init1 + metal;
  char cuttingL[100]; //Convert to string
  sprintf(cuttingL, "rez2 > %f && rez2 < %f", z_init1, z_fin1);

  cout << "histogram cut string u: " << cuttingU << endl;
  cout << "histogram cut string l: " << cuttingL << endl;
  
  //Convert strings to TCut
  TCut hist_cutU("metal layer upper", cuttingU);
  TCut hist_cutL("metal layer lower", cuttingL);

  //Configuration
  bool iter_dir_x = false;
  bool iter_dir_45 = true;
  int file_ct = 6;        
  bool u_gauss = true;
  bool l_gauss = true;


   //Chain for obtaining the initial cloud position and size
   TChain* tch_fln = new TChain("retree","chain of e- trees");

   tch_fln->Add("xyrdo_73ArCo2_350dv_xy150i3u_xbfield15_1.root");
    
    //tch_fln->Add("xyrdo_xite_quadgem350_6mmtotal_20um_c1.root");
    //tch_fln->Add("xyrdo_73ArCo2_350dv_5pi-_300mev_xy50i30u_xbfield15_1.root");
    //tch_fln->Add("xyrdo_73ArCo2_350dv_60i30ux_4iterot_propstag_1.root");
    // tch_fln->Add("xyrdo_73ArCo2_350dv_1pi-_5mev_x50i50u_xbfield15_1.root");
    //tch_fln->Add(" xyrdo_73ArCo2_350dv_80i30ux_4iterot_bf15_1.root");
    //tch_fln->Add("xyrdo_73ArCo2_350dv_x100i3u_xbfield15_1.root");
    // tch_fln->Add("xyrdo_xite_quadgem350_xat0_30um250c_c1.root");
    //tch_fln->Add("xyrdo_73ArCo2_350dv_60i30u_4iterot_4.root");
  
  
    //int max_iter = tch_fln->GetEntries() - 2;
    //int max_iter = 10;

    //Draw cloud to histograms of electron x and y single-axis histograms, and find them
    tch_fln->Draw("rex0>>initer_elecs0", hist_cutL, "", 1, 0);
    tch_fln->Draw("rey0>>initer_elecs1", hist_cutU, "", 1, 0);
    TH1F *htemp0 = (TH1F*)gDirectory->FindObject("initer_elecs0");
    TH1F *htemp1 = (TH1F*)gDirectory->FindObject("initer_elecs1");


    //Obtain the maximum and minimum on the axes, calculate range
    double offsetLmax = htemp0->GetXaxis()->GetXmax();
    double offsetLmin = htemp0->GetXaxis()->GetXmin();
    
    double offsetUmax = htemp1->GetXaxis()->GetXmax();
    double offsetUmin = htemp1->GetXaxis()->GetXmin();
    
    double Uoffrange = (offsetUmax - offsetUmin) / 2;
    double Xoffrange = (offsetLmax - offsetLmin) / 2;

    cout << "lower strips init will be from " << offsetLmin << " - " << offsetLmax << endl; 
    cout << "upper strips init will be from " << offsetUmin << " - " << offsetUmax << endl; 
    
    delete htemp0;
    delete htemp1;

    
    
    //Make the storage file
  TFile *store_file = new TFile("iterative_e-stgb_fit_aggr.root","RECREATE");
  

  //Loop 1: File
  for(int fl = 1;fl <= file_ct;fl++){
  
    char filename_iter[50];
    
    // sprintf(filename_iter, "xyrdo_xite_quadgem350_samev_ncont_%d.root", file_ct);
    //sprintf(filename_iter, "xyrdo_xite_quadgem350_6mmtotal_20um_c%d.root", fl);
    // sprintf(filename_iter, "xyrdo_xite_quadgem350_xat0_30um250c_c%d.root", fl);
    //  sprintf(filename_iter, "xyrdo_73ArCo2_350dv_1pi-_5mev_x50i50u_xbfield15_%d.root", fl);
    //sprintf(filename_iter, "xyrdo_73ArCo2_350dv_50pi-xyot_xy5i30u_xbfield15_%d.root", fl);
    sprintf(filename_iter, "xyrdo_73ArCo2_350dv_xy150i3u_xbfield15_%d.root", fl);
    //tch_fln->Add("xyrdo_73ArCo2_350dv_x100i3u_xbfield15_1.root");
    //sprintf(filename_iter, "xyrdo_73ArCo2_350dv_60i30ux_4iterot_propstag_%d.root", fl);
    //sprintf(filename_iter, "xyrdo_73ArCo2_350dv_60i30u_4iterot_%d.root", fl);
    

    TFile *fin =  new TFile(filename_iter);
    fin->cd();

    TTree *out_etree = (TTree*)fin->Get("retree");
    //TTree etree;
    //etree=(TTree*)fin->Get("etree");
    //out_etree->Print();
    
    //Get number of iterations
    int max_iter = out_etree->GetEntries() - 1;
    // Double_t iteration_distance = 0.005;
    // Double_t iteration_distance = 0.003;
    Double_t iteration_distance = 0.0212;


    //Declaration of misc result histograms
    TGraph* recons_strx = new TGraph(max_iter);
    TGraph* recons_stry = new TGraph(max_iter);
    TH2* cutelecsl = new TH2D("cutelecsL", "L strips alle", binct, offsetLmin, offsetLmax + max_iter*iteration_distance, binct, offsetUmin, offsetUmax + max_iter*iteration_distance);
    TH2* cutelecsu = new TH2D("cutelecsU", "U strips alle", binct, offsetLmin, offsetLmax + max_iter*iteration_distance, binct, offsetUmin, offsetUmax + max_iter*iteration_distance);
    TH2* cutelecsla = new TH2D("cutelecsLA", "L strips singular cloud", binct, offsetLmin, offsetLmax, binct, offsetUmin, offsetUmax);
    TH2* cutelecsua = new TH2D("cutelecsUA", "U strips singular cloud", binct, offsetLmin, offsetLmax, binct, offsetUmin, offsetUmax);
    TH2* cutelecsa = new TH2D("cutelecsA", "singular cloud", binct, offsetLmin, offsetLmax, binct, offsetUmin, offsetUmax);
    TH2* cutelecsta = new TH2D("cutelecsS", "singular cloud strips", binct, offsetLmin, offsetLmax, binct, offsetUmin, offsetUmax);

    //Loop 2: Iteration
    for(int cur_iter = 1;cur_iter <= max_iter;cur_iter++){
        int treeiter;
        Int_t iterationx;
        Int_t iterationy;
        
        //Iteration direction from configurations earlier
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
        //Calculate the electron cloud's upper and lower limits for the current iteration
        Double_t l_it_min = offsetLmin - 0.1 + iteration_distance * iterationx;
        Double_t l_it_max = offsetLmax + 0.1 + iteration_distance * iterationx;
        Double_t u_it_min = offsetUmin - 0.1 + iteration_distance * iterationy;
        Double_t u_it_max = offsetUmax + 0.1 + iteration_distance * iterationy;
        
	
        Double_t e_centerL;
        //Draw electrons on RO board layer to histogram
        TH1 *h_allL = new TH1D("all_elecsL", "electrons L", binct, l_it_min, l_it_max);
        out_etree->Draw("rex2>>all_elecsL", hist_cutL, "", 1, treeiter);
        //If electron cloud is gaussian
        if(l_gauss){
            //Create gaussian fit, get center of the fit
            TF1 *fit_hL = new TF1("fit_h_l","gaus",l_it_min,l_it_max);
            fit_hL->SetName("fit_h_l");
            h_allL->Fit("fit_h_l","R");
            e_centerL = fit_hL->GetParameter(1);
            delete fit_hL;
            
        }else{
            //Else just find histogram mean as cloud center
           	e_centerL = h_allL->GetMean();
        }
        cout << "Iter: " << cur_iter << " - Fit Center L " << e_centerL << endl;
        cout <<  l_it_min << " Lower strip Range " << l_it_max << endl;  

        
        //Same process for Y
        Double_t e_centerU;
        TH1 *h_allU = new TH1D("all_elecsU", "electrons U", binct, u_it_min, u_it_max);
        out_etree->Draw("rey2>>all_elecsU", hist_cutU, "", 1, treeiter);
        if(u_gauss){
           	//Double_t e_centerU = h_allU->GetMean();
            TF1 *fit_hU = new TF1("fit_h_u","gaus",u_it_min,u_it_max);
            fit_hU->SetName("fit_h_u");
            h_allU->Fit("fit_h_u","R");
           	e_centerU = fit_hU->GetParameter(1);
           	delete fit_hU;
        }else{
           	e_centerU = h_allU->GetMean();
        }


       	    cout << "Iter: " << cur_iter << " - Fit Center U " << e_centerU << endl;
            cout <<  u_it_min << " Upper strip Range " << u_it_max << endl;
       	
        
        //If empty histograms are detected, alert, throw away result and continue
        if(h_allU->Integral() == 0.0){
       	    cout << "empty histograms detected on upper strips" << endl;
       	}
       	if(h_allL->Integral() == 0.0){
       	    cout << "empty histograms detected on lower strips" << endl;
       	}       	
       	if(h_allU->Integral() == 0.0 || h_allL->Integral() == 0.0){
           	delete h_allU;
            delete h_allL;
       	    continue;
       	}
        

        //Use calculated lower bound of the cloud to find which multiplier of (strip+pitch) should the iteration start at
        Double_t xpitch_mult = floor(abs((l_it_min + (pitch / 2 - dstrip / 2)) / pitch)); 
        //Calculate the initial offset from the multiplier
        Double_t initial_offset_x = pitch / 2 - dstrip / 2 - pitch * xpitch_mult;
        cout << "lower starts at " << xpitch_mult << " units wide and at " << initial_offset_x  << endl;
        
        
        //Iterations loop starts at the initial offset
        Double_t xs_iter = initial_offset_x;
        //Blank double for the contents of one strip
        Double_t xs_iter_con = 0.0;
        //Declare lists
        std::list<Double_t> contentlistL;
        std::list<Double_t> coordlistL;
        
        //Loop 3(d): Strips
        while (xs_iter < e_centerL || xs_iter_con > 0.0)
        {
          //Make cut for "between two strips" using current iterator
          Double_t xs_iter_up = xs_iter + dstrip + 0.001;
          char cuttingStL[100];
          sprintf(cuttingStL, "rex2 > %f && rex2 < %f", xs_iter, xs_iter_up);
          //cout << cuttingStL << endl;
          TCut hist_cut_itL("L iterative cutter", cuttingStL);

            
          //Combine the between-strip cut and the layer-only cut to get the current strip's 1D electron distribution
          TH1 *x_temphist = new TH1S("xtemphist", "iterator histogram", 50, xs_iter, xs_iter_up);
          out_etree->Draw("rex2>>xtemphist", hist_cutL && hist_cut_itL, "colz", 1, treeiter);

          //Obtain content of strip
          xs_iter_con = x_temphist->Integral();
          //Obtain center of strip
          Double_t stripmean = xs_iter + (dstrip / 2);

          cout << "l strip at: " << stripmean << endl;
          cout << "l strip content: " << xs_iter_con << endl;
          cout << "l cut used: " << cuttingStL << " && " << hist_cutL << endl;
          
          //Put the content and center into the respective lists ONLY if there is significant amount of electrons
          if (xs_iter_con > 10.0)
          {
            contentlistL.push_back(xs_iter_con);
            coordlistL.push_back(stripmean);
          }
          
          //Advance the iterator by pitch
          xs_iter += pitch;
          delete x_temphist;
        }
        //Loop 3(d) Ends

        Double_t ypitch_mult = floor(abs((u_it_min + (pitch / 2 - ustrip / 2)) / pitch)); 
        Double_t initial_offset_y = pitch / 2 - ustrip / 2 - pitch * ypitch_mult;
        
        cout << "upper starts at " << ypitch_mult << " units wide and at " << initial_offset_y <<  endl;
        Double_t ys_iter = initial_offset_y;
        Double_t ys_iter_con = 0.0;
        std::list<Double_t> contentlistU;
        std::list<Double_t> coordlistU;
        //Loop 3(u) starts
        while (ys_iter < e_centerU || ys_iter_con > 0.0)
        {
          Double_t ys_iter_up = ys_iter + ustrip + 0.001;
          char cuttingStU[100];
          sprintf(cuttingStU, "rey2 > %f && rey2 < %f", ys_iter, ys_iter_up);
          //cout << cuttingStU << endl;

          TCut hist_cut_itU("U iterative cutter", cuttingStU);

          TH1 *y_temphist = new TH1S("ytemphist", "iterator histogram", 50, ys_iter, ys_iter_up);
          out_etree->Draw("rey2>>ytemphist", hist_cutU && hist_cut_itU, "colz", 1, treeiter);

          ys_iter_con = y_temphist->Integral();
          Double_t stripmean = ys_iter + (ustrip / 2);

          cout << "u strip at: " << stripmean << endl;
          cout << "u strip content: " << ys_iter_con << endl;
          cout << "u cut used: " << cuttingStU << " && " << hist_cutU << endl;
          if (ys_iter_con > 10.0)
          {
            contentlistU.push_back(ys_iter_con);
            coordlistU.push_back(stripmean);
          }
          ys_iter += pitch;
          delete y_temphist;
        }
        //Loop 3(u) ends
        
        //Iterate through the filled content/coordinate lists again
        std::list<Double_t>::iterator ictx = contentlistL.begin();
        std::list<Double_t>::iterator icox = coordlistL.begin();
        //Storage for sum of electrons and sum of (e*coordinates)
        Double_t xsum = 0.0;
        Double_t x_t_coord = 0.0;

        //Loop 4: centroid
        while (ictx != contentlistL.end())
        {
          xsum += *ictx;
          x_t_coord += *ictx * *icox;
          //cout << "lowere cd: " << xsum << " " << x_t_coord << " " << *ictx * *icox << endl;
          ictx++;
          icox++;
        }
        //Loop 4 ends
        Double_t xcentroid = x_t_coord / xsum;



        std::list<Double_t>::iterator icty = contentlistU.begin();
        std::list<Double_t>::iterator icoy = coordlistU.begin();

        Double_t ysum = 0.0;
        Double_t y_t_coord = 0.0;

        while (icty != contentlistU.end())
        {
          ysum += *icty;
          y_t_coord += *icty * *icoy;
          //cout << "upper cd: " << ysum << " " << y_t_coord << " " << *icty * *icoy << endl;
          icty++;
          icoy++;
        }
        Double_t ycentroid = y_t_coord / ysum;


        
        cout << "centroid lower det " << xcentroid << endl;
        cout << "centroid upper det " << ycentroid << endl;
        cout << "cloud center lower " << e_centerL << endl;
        cout << "cloud center upper " << e_centerU << endl;
        cout << xcentroid - e_centerL << " // " << ycentroid - e_centerU << endl;


        //If both the cloud center and the centroid shows up to 0.0 (very small probability unless something goes wrong) then omit cloud from final histogram
        if(e_centerL == 0.0 || xcentroid == 0.0){
            cout << "l omission" << endl;
        //Else put into histogram
        }else{
            recons_strx->SetPoint(cur_iter-1,e_centerL,xcentroid);
            Double_t xvx = xcentroid - e_centerL;  
            finalhist_l->Fill(xvx);
        }
        
        
        if(e_centerU == 0.0 || ycentroid == 0.0){
            cout << "u omission" << endl;

        }else{
            recons_stry->SetPoint(cur_iter-1,e_centerU,ycentroid);
            Double_t yvy = ycentroid - e_centerU;
            finalhist_u->Fill(yvy);
            
        }


        
        delete h_allL;
        delete h_allU;
        //delete fit_hL;
        //delete fit_hU;
      }
      //Loop 2 ends


      char graphxname[20];
      char graphyname[20];
      sprintf(graphxname, "xre_xt_%d", fl);
      sprintf(graphyname, "yre_yt_%d", fl);
      
      char h2xname[20];
      char h2yname[20];
      sprintf(h2xname, "lstrips_%d", fl);
      sprintf(h2yname, "ustrips_%d", fl);
      
      char h2xnamea[20];
      char h2ynamea[20];
      sprintf(h2xnamea, "lstrips_a_%d", fl);
      sprintf(h2ynamea, "ustrips_a_%d", fl);
      
      
      out_etree->Draw("rey2:rex2>>cutelecsL", hist_cutL, "COLZ2");
      out_etree->Draw("rey2:rex2>>cutelecsU", hist_cutU, "COLZ2");
      cutelecsl->SetName(h2xname);
      cutelecsu->SetName(h2yname);
      
      out_etree->Draw("rey2:rex2>>cutelecsLA", hist_cutL, "COLZ2", 1, 0);
      out_etree->Draw("rey2:rex2>>cutelecsUA", hist_cutU, "COLZ2", 1, 0);
      cutelecsla->SetName(h2xnamea);
      cutelecsua->SetName(h2ynamea);
      
      
      char h2aname[10];
      char h2tname[15];
      sprintf(h2aname, "st_a_%d", fl);
      sprintf(h2tname, "strs_a_%d", fl);
      out_etree->Draw("rey2:rex2>>cutelecsA", "" , "COLZ2", 1, 0);
      cutelecsa->SetName(h2aname);
      out_etree->Draw("rey2:rex2>>cutelecsS", hist_cutL || hist_cutU, "COLZ2", 1, 0);
      cutelecsta->SetName(h2tname);

      store_file->cd();
      recons_strx->SetName(graphxname);
      recons_stry->SetName(graphyname);
      recons_strx->SetTitle("cloud reconstructed L vs actual L");
      recons_stry->SetTitle("cloud reconstructed U vs actual U");
      recons_strx->SetMarkerStyle(3);
      recons_stry->SetMarkerStyle(3);
      recons_strx->Write();
      recons_stry->Write();
      cutelecsl->Write();
      cutelecsu->Write();
      cutelecsla->Write();
      cutelecsua->Write();
      cutelecsa->Write();
      delete cutelecsa;
      
      cutelecsta->Write();
      delete cutelecsta;
      delete cutelecsl;
      delete cutelecsu;      
      delete cutelecsla;
      delete cutelecsua;
      delete recons_strx;
      delete recons_stry;

    }
    //Loop 1 ends
  gROOT->SetBatch(kFALSE);
  //TCanvas* c10 = new TCanvas("c10");
  //c10->cd();
  finalhist_u->SetMarkerStyle(20);
  finalhist_u->SetMarkerSize(1);
  finalhist_l->SetMarkerStyle(20);
  finalhist_l->SetMarkerSize(1);
  //finalhist_x->Draw("P");
  //c10->SaveAs("histby.png");
  //tch_fln->Draw("rey2:rex2","","colz",1,1);

  store_file->cd();
  finalhist_u->Write();
  finalhist_l->Write();
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

  // cout << "content U:" << endl;
  // for (Double_t a : contentlistU)
  // {
  //   cout << a << ",";
  // }
  // cout << endl;
  // cout << "coordinates U:" << endl;
  // for (Double_t b : coordlistU)
  // {
  //   cout << b << ",";
  // }

  //delete h_allX;
  //delete h_allU;
}
