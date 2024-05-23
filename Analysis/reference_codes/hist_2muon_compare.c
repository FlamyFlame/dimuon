{
  const int nPair = 3;
  const int nPrtGp = 5;
  
  std::ifstream fin;
  fin.open("muon_parent_2muon.txt");

  if (!fin){
    std::cout << "Unable to open file muon_parent_2muon.txt";
    exit(1);
  }

  std::string line, word;
  std::vector<std::string> row;
  int row_size = 12;
  
  std::vector<int> eventn_list;
  std::vector<int> muon_parent;
  std::vector<int> muon_id;
  std::vector<double> muon_pt;
  std::vector<double> muon_eta;
  std::vector<double> muon_phi;
  std::vector<bool> tau_tag;
  std::vector<bool> c_tag;
  std::vector<double> weight_list;
  std::vector<double> xgamma_list;
  std::vector<double> ygamma_list;
  std::vector<double> xa_list;
  
  while(std::getline(fin,line)){
    row.clear();
    if (line[0] == '#')
      continue;

    std::stringstream s(line);
    while(std::getline(s,word,'\t')){
      row.push_back(word);
    }
    assert (row.size() == row_size);
    
    eventn_list.push_back(std::stoi(row[0]));
    muon_parent.push_back(std::stoi(row[1]));
    muon_id.push_back(std::stoi(row[2]));
    muon_pt.push_back(std::stod(row[3]));
    muon_eta.push_back(std::stod(row[4]));
    muon_phi.push_back(std::stod(row[5]));
    tau_tag.push_back(std::stoi(row[6]));
    c_tag.push_back(std::stoi(row[7]));
    weight_list.push_back(std::stod(row[8]));
    xgamma_list.push_back(std::stod(row[9]));
    ygamma_list.push_back(std::stod(row[10]));
    xa_list.push_back(std::stod(row[11]));
  }

  fin.close();

  TCanvas *c1 = new TCanvas();
  c1->SetWindowSize(1100,720);
  c1->Divide(2,2);
  auto pi = TMath::Pi();
   
  double leadpTBins[101];
  double subleadpTBins[101];
  double deltaEtaBins[101];
  for(int i = 0; i <= 100; i++) // Index should span 5 orders of magnitude. 10e-5 -> 10e0             
  {
    leadpTBins[i] = 53 * pow(10.0, ((float)(i-100))*0.0113);
    subleadpTBins[i] = 27 * pow(10.0, ((float)(i-100))*0.00837);
    deltaEtaBins[i] = 4.38 * pow(10.0, ((float)(i-100))*0.04);
  }

  //TH2F *h1 = new TH2F("h1","P_{T}_{#mu,lead} versus P_{T}_{#mu,sublead}",100,subleadpTBins,100,leadpTBins);
  //TH2F *h2 = new TH2F("h2","#eta_{#mu,lead} versus #eta_{#mu,sublead}",28,-2.7,2.7,28,-2.7,2.7);
  //TH2F *h3 = new TH2F("h3","#bar{#eta} versus |#Delta #eta|",100,deltaEtaBins,28,-2.7,2.7);
  TH2F* h [nPair];
  h[0] = new TH2F("h1","",100,subleadpTBins,100,leadpTBins);
  h[1] = new TH2F("h2","",28,-2.7,2.7,28,-2.7,2.7);
  h[2] = new TH2F("h3","",100,deltaEtaBins,28,-2.7,2.7);
  
  int cur_eventn = eventn_list[0];
  int cur_event_ind = 0; // the index of the first muon in the current event
  double cur_weight = weight_list[0];
  std::vector<double> cur_muon_pt;
  std::vector<double> cur_muon_eta;
  std::vector<std::pair<double,int>> pt_ind_arr;
  for (int i = 0; i < eventn_list.size(); i++){
    if (eventn_list[i] == cur_eventn){
      // push back the muon kinematics for the current event
      cur_muon_pt.push_back(muon_pt[i]);
      cur_muon_eta.push_back(muon_eta[i]);
      pt_ind_arr.push_back(std::make_pair(muon_pt[i],i-cur_event_ind));
    }
    else{
      // add pairs to histograms
      unsigned int n = cur_muon_pt.size();
      std::sort(cur_muon_pt.begin(), cur_muon_pt.end()); //sorting of cur_muon_pt doesn't affect the vector pt_ind_arr
      std::sort(pt_ind_arr.begin(),pt_ind_arr.end()); //sort by the first row of pt_ind_arr (muon pt)
      //the second row now contains the ordering of muon pt in the current event (from low to high)
      std::vector<double> eta_sort(n); //muon eta in the current event sorted by pt (from low pt to high pt)
      for (unsigned int j = 0; j < n; j++)
         eta_sort[j] = cur_muon_eta.at(pt_ind_arr[j].second);
      h[0]->Fill(cur_muon_pt.at(n-2), cur_muon_pt.at(n-1), cur_weight); //(pt_sublead, pt_lead, weight)
      h[1]->Fill(eta_sort[n-2], eta_sort[n-1], cur_weight); //(eta_sublead, eta_lead, weight)
      double sum = std::accumulate(cur_muon_eta.begin(),cur_muon_eta.end(),0.0); //0.0: init val
      h[2]->Fill(abs(eta_sort[n-2] - eta_sort[n-1]), sum/n, cur_weight); //(|Delta eta|, eta-avg, weight)

      // clear the kinematic arrays & update cur_eventn
      cur_muon_pt.clear();
      cur_muon_eta.clear();
      pt_ind_arr.clear();
      cur_eventn = eventn_list[i];
      cur_event_ind = i;
      cur_weight = weight_list[i];
      
      //push back kinematics of the current muon
      cur_muon_pt.push_back(muon_pt[i]);
      cur_muon_eta.push_back(muon_eta[i]);
      pt_ind_arr.push_back(std::make_pair(muon_pt[i],i-cur_event_ind));

    }
  }
  

  const char *xtitle_arr[nPair] = {"p_{T #mu,sublead}","#eta_{#mu,sublead}","|#Delta #eta|"};
  const char *ytitle_arr[nPair] = {"p_{T #mu,lead}","#eta_{#mu,lead}","#bar{#eta}"};
   //float ytitle_offset_arr[nPair] = {1.0,1.0,1.0};
  bool logx[nPair] = {true,false,true};
  bool logy[nPair] = {true,false,false};
  for (int i = 0; i < nPair; i++){
    c1->cd(i+1);
    h[i]->SetStats(0);
    //h[i]->SetTitle(title_arr[j]);
    h[i]->GetXaxis()->SetTitle(xtitle_arr[i]);
    h[i]->GetYaxis()->SetTitle(ytitle_arr[i]);
    h[i]->GetXaxis()->SetTitleSize(h[i]->GetXaxis()->GetTitleSize()*2);
    h[i]->GetXaxis()->SetLabelSize(h[i]->GetXaxis()->GetLabelSize()*1.35);
    h[i]->GetYaxis()->SetTitleSize(h[i]->GetYaxis()->GetTitleSize()*2);
    h[i]->GetYaxis()->SetLabelSize(h[i]->GetYaxis()->GetLabelSize()*1.35);
    h[i]->GetZaxis()->SetLabelSize(h[i]->GetZaxis()->GetLabelSize()*1.35);
    h[i]->GetYaxis()->SetTitleOffset(0.45);
    h[i]->GetXaxis()->SetTitleOffset(0.6);
    gPad->SetLogx(logx[i]);
    gPad->SetLogy(logy[i]);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.135);
    gPad->SetRightMargin(0.135);
    h[i]->Scale(1.0/h[i]->Integral("width"));
    h[i]->Draw("COLZ");
  }
  c1->Print("2muon_kin_compare_4gev.pdf");
}
