// fraction of SS-filter-pass events whose leading qualifying SS pair (filter defn) exceeds thresholds
void ss_pair_pt_threshold_fractions(const char* g){ TChain ch("HeavyIonD3PD"); ch.Add(g);
  std::vector<float> *pt1=0,*eta1=0,*pt2=0,*eta2=0,*m=0,*pt=0; std::vector<int> *c1=0,*c2=0;
  ch.SetBranchAddress("truth_mupair_pt1",&pt1); ch.SetBranchAddress("truth_mupair_eta1",&eta1); ch.SetBranchAddress("truth_mupair_ch1",&c1);
  ch.SetBranchAddress("truth_mupair_pt2",&pt2); ch.SetBranchAddress("truth_mupair_eta2",&eta2); ch.SetBranchAddress("truth_mupair_ch2",&c2);
  ch.SetBranchAddress("truth_mupair_m",&m); ch.SetBranchAddress("truth_mupair_pt",&pt);
  double T[4]={40,45,50,60}; int n[4]={0}; long N=ch.GetEntries();
  for(long i=0;i<N;++i){ch.GetEntry(i); double lead=-1;
    for(size_t k=0;k<m->size();++k){ if(fabs(pt1->at(k))*1e-3<=3.7||fabs(pt2->at(k))*1e-3<=3.7||fabs(eta1->at(k))>=2.5||fabs(eta2->at(k))>=2.5) continue;
      if(c1->at(k)*c2->at(k)<=0) continue; double mm=m->at(k)*1e-3; if(mm>0&&mm<10) lead=std::max(lead,(double)pt->at(k)*1e-3);}
    for(int j=0;j<4;++j) if(lead>T[j]) n[j]++;}
  for(int j=0;j<4;++j) printf("  pT^mumu > %2.0f: %4d / %ld = %.3f\n",T[j],n[j],N,double(n[j])/N);}
