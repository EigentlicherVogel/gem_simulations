THStack* single_use_hist()
{
THStack *hst = new THStack("hst","");
TChain* retree_iter_1 = new TChain("retree_iter_1","tchainiter1");
TChain* retree_iter_2 = new TChain("retree_iter_2","tchainiter2");
TChain* retree_iter_3 = new TChain("retree_iter_3","tchainiter3");
TChain* retree_iter_4 = new TChain("retree_iter_4","tchainiter4");
TChain* retree_iter_5 = new TChain("retree_iter_5","tchainiter5");
TChain* retree_iter_6 = new TChain("retree_iter_6","tchainiter6");
retree_iter_1->Add("chev_quadgem_73p350_r4_*.root");
retree_iter_2->Add("chev_quadgem_73p350_r4_*.root");
retree_iter_3->Add("chev_quadgem_73p350_r4_*.root");
retree_iter_4->Add("chev_quadgem_73p350_r4_*.root");
retree_iter_5->Add("chev_quadgem_73p350_r4_*.root");
retree_iter_6->Add("chev_quadgem_73p350_r4_*.root");

TH2* rtt1 = new TH2D("ro_xy_iter1", "Iterative RDO", 160, -0.1, 0.3, 80, -0.05, 0.15);
TH2* rtt2 = new TH2D("ro_xy_iter2", "Iterative RDO", 160, -0.05, 0.35, 80, -0.05, 0.15);
TH2* rtt3 = new TH2D("ro_xy_iter3", "Iterative RDO", 160, 0.0, 0.4, 80, -0.05, 0.15);
TH2* rtt4 = new TH2D("ro_xy_iter4", "Iterative RDO", 160, 0.05, 0.45, 80, -0.05, 0.15);
TH2* rtt5 = new TH2D("ro_xy_iter5", "Iterative RDO", 160, 0.1, 0.5, 80, -0.05, 0.15);
TH2* rtt6 = new TH2D("ro_xy_iter6", "Iterative RDO", 160, 0.15, 0.55, 80, -0.05, 0.15);

retree_iter_1->Draw("rey2:rex2>>ro_xy_iter1","rez2<0.01","colz");
retree_iter_2->Draw("rey2:rex2>>ro_xy_iter2","rez2<0.01","colz");
retree_iter_3->Draw("rey2:rex2>>ro_xy_iter3","rez2<0.01","colz");
retree_iter_4->Draw("rey2:rex2>>ro_xy_iter4","rez2<0.01","colz");
retree_iter_5->Draw("rey2:rex2>>ro_xy_iter5","rez2<0.01","colz");
retree_iter_6->Draw("rey2:rex2>>ro_xy_iter6","rez2<0.01","colz");

hst->Add(ro_xy_iter1);
hst->Add(ro_xy_iter2);
hst->Add(ro_xy_iter3);
hst->Add(ro_xy_iter4);
hst->Add(ro_xy_iter5);
hst->Add(ro_xy_iter6);

return hst;
}
