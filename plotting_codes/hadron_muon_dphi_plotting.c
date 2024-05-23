ROOT::VecOps::RVec<float> HadronMuonDphiCalc(float m1_phi, float m2_phi, ROOT::VecOps::RVec<float> hadron1_pt_eta_phi_m, ROOT::VecOps::RVec<float> hadron2_pt_eta_phi_m){
	// function to be applied to either both-b case or both-c case
	ROOT::VecOps::RVec<float> dphi_vec = {};
	if (hadron1_pt_eta_phi_m.size() == 0 || hadron2_pt_eta_phi_m.size() == 0){
		return dphi_vec; // return empty vector if either of the two muons aren't originally from b
	}
	float dphi1 = m1_phi - hadron1_pt_eta_phi_m[2];
	dphi1 = atan2(sin(dphi1),cos(dphi1));
	float dphi2 = m2_phi - hadron2_pt_eta_phi_m[2];
	dphi2 = atan2(sin(dphi2),cos(dphi2));
	
	dphi_vec.push_back(dphi1);
	dphi_vec.push_back(dphi2);
	return dphi_vec;
}

ROOT::VecOps::RVec<float> CHadronMuonDphiCalc(bool both_from_c, float m1_phi, float m2_phi, ROOT::VecOps::RVec<float> hadron1_pt_eta_phi_m, ROOT::VecOps::RVec<float> hadron2_pt_eta_phi_m){
	ROOT::VecOps::RVec<float> c_dphi_vec = {};
	if (!both_from_c){
		return c_dphi_vec; // return empty if not both from c
	}
	return HadronMuonDphiCalc(m1_phi, m2_phi, hadron1_pt_eta_phi_m, hadron2_pt_eta_phi_m);
}

void hadron_muon_dphi_plotting(){
	ROOT::RDataFrame df_powheg_bb_ss("muon_pair_tree_sign1", "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/muon_pairs_mc_truth_bb_1-5.root");
	ROOT::RDataFrame df_powheg_bb_op("muon_pair_tree_sign2", "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/bb_full_sample/muon_pairs_mc_truth_bb_1-5.root");
	ROOT::RDataFrame df_powheg_cc_ss("muon_pair_tree_sign1", "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/muon_pairs_mc_truth_cc_1-5.root");
	ROOT::RDataFrame df_powheg_cc_op("muon_pair_tree_sign2", "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/cc_full_sample/muon_pairs_mc_truth_cc_1-5.root");
	ROOT::RDataFrame df_pythia_ss("muon_pair_tree_sign1", "/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_after0322.root");
	ROOT::RDataFrame df_pythia_op("muon_pair_tree_sign2", "/usatlas/u/yuhanguo/usatlasdata/pythia/muon_pairs_pythia_after0322.root");

	ROOT::RDF::RNode df_powheg_bb_ss_updated = df_powheg_bb_ss.Filter("m1.pt + m2.pt > 10").Filter("abs(deta) > 0.8").Define("ptavg",[](float pt1, float pt2){return static_cast<float>((pt1 + pt2)/2.);}, {"m1.pt", "m2.pt"}).Define("b_hadron_muon_dphi", HadronMuonDphiCalc, {"m1.phi", "m2.phi", "m1_last_b_hadron_prt_pt_eta_phi_m", "m2_last_b_hadron_prt_pt_eta_phi_m"}).Redefine("weight", [](double w){return static_cast<double>(w / 1000000.);},{"weight"});
	ROOT::RDF::RNode df_powheg_bb_op_updated = df_powheg_bb_op.Filter("m1.pt + m2.pt > 10").Filter("abs(deta) > 0.8").Define("ptavg",[](float pt1, float pt2){return static_cast<float>((pt1 + pt2)/2.);}, {"m1.pt", "m2.pt"}).Define("b_hadron_muon_dphi", HadronMuonDphiCalc, {"m1.phi", "m2.phi", "m1_last_b_hadron_prt_pt_eta_phi_m", "m2_last_b_hadron_prt_pt_eta_phi_m"}).Redefine("weight", [](double w){return static_cast<double>(w / 1000000.);},{"weight"});
	ROOT::RDF::RNode df_powheg_cc_ss_updated = df_powheg_cc_ss.Filter("m1.pt + m2.pt > 10").Filter("abs(deta) > 0.8").Define("ptavg",[](float pt1, float pt2){return static_cast<float>((pt1 + pt2)/2.);}, {"m1.pt", "m2.pt"}).Define("c_hadron_muon_dphi", CHadronMuonDphiCalc, {"both_from_c", "m1.phi", "m2.phi", "m1_last_hf_hadron_prt_pt_eta_phi_m", "m2_last_hf_hadron_prt_pt_eta_phi_m"}).Redefine("weight", [](double w){return static_cast<double>(w / 1000000.);},{"weight"});
	ROOT::RDF::RNode df_powheg_cc_op_updated = df_powheg_cc_op.Filter("m1.pt + m2.pt > 10").Filter("abs(deta) > 0.8").Define("ptavg",[](float pt1, float pt2){return static_cast<float>((pt1 + pt2)/2.);}, {"m1.pt", "m2.pt"}).Define("c_hadron_muon_dphi", CHadronMuonDphiCalc, {"both_from_c", "m1.phi", "m2.phi", "m1_last_hf_hadron_prt_pt_eta_phi_m", "m2_last_hf_hadron_prt_pt_eta_phi_m"}).Redefine("weight", [](double w){return static_cast<double>(w / 1000000.);},{"weight"});
	ROOT::RDF::RNode df_pythia_ss_updated = df_pythia_ss.Filter("m1.pt + m2.pt > 10").Filter("abs(deta) > 0.8").Define("ptavg",[](float pt1, float pt2){return static_cast<float>((pt1 + pt2)/2.);}, {"m1.pt", "m2.pt"}).Define("b_hadron_muon_dphi", HadronMuonDphiCalc, {"m1.phi", "m2.phi", "m1_last_b_hadron_prt_pt_eta_phi_m", "m2_last_b_hadron_prt_pt_eta_phi_m"}).Define("c_hadron_muon_dphi", CHadronMuonDphiCalc, {"both_from_c", "m1.phi", "m2.phi", "m1_last_hf_hadron_prt_pt_eta_phi_m", "m2_last_hf_hadron_prt_pt_eta_phi_m"}).Redefine("weight", [](double w){return static_cast<double>(w * 1000000.);},{"weight"});
	ROOT::RDF::RNode df_pythia_op_updated = df_pythia_op.Filter("m1.pt + m2.pt > 10").Filter("abs(deta) > 0.8").Define("ptavg",[](float pt1, float pt2){return static_cast<float>((pt1 + pt2)/2.);}, {"m1.pt", "m2.pt"}).Define("b_hadron_muon_dphi", HadronMuonDphiCalc, {"m1.phi", "m2.phi", "m1_last_b_hadron_prt_pt_eta_phi_m", "m2_last_b_hadron_prt_pt_eta_phi_m"}).Define("c_hadron_muon_dphi", CHadronMuonDphiCalc, {"both_from_c", "m1.phi", "m2.phi", "m1_last_hf_hadron_prt_pt_eta_phi_m", "m2_last_hf_hadron_prt_pt_eta_phi_m"}).Redefine("weight", [](double w){return static_cast<double>(w * 1000000.);},{"weight"});
	auto h_b_hadron_powheg_bb_ss = df_powheg_bb_ss_updated.Histo1D<ROOT::VecOps::RVec<float>,double>({"h_b_hadron_muon_dphi_powheg_bb_ss","h_b_hadron_muon_dphi_powheg_bb_ss",32, -1,1},"b_hadron_muon_dphi","weight");
	auto h_b_hadron_powheg_bb_op = df_powheg_bb_op_updated.Histo1D<ROOT::VecOps::RVec<float>,double>({"h_b_hadron_muon_dphi_powheg_bb_op","h_b_hadron_muon_dphi_powheg_bb_op",32, -1,1},"b_hadron_muon_dphi","weight");
	auto h_c_hadron_powheg_cc_ss = df_powheg_cc_ss_updated.Histo1D<ROOT::VecOps::RVec<float>,double>({"h_c_hadron_muon_dphi_powheg_cc_ss","h_c_hadron_muon_dphi_powheg_cc_ss",32, -1,1},"c_hadron_muon_dphi","weight");
	auto h_c_hadron_powheg_cc_op = df_powheg_cc_op_updated.Histo1D<ROOT::VecOps::RVec<float>,double>({"h_c_hadron_muon_dphi_powheg_cc_op","h_c_hadron_muon_dphi_powheg_cc_op",32, -1,1},"c_hadron_muon_dphi","weight");
	auto h_b_hadron_pythia_ss 	= df_pythia_ss_updated.Histo1D<ROOT::VecOps::RVec<float>,double>({"h_b_hadron_muon_dphi_pythia_ss","h_b_hadron_muon_dphi_pythia_ss",32, -1,1},"b_hadron_muon_dphi","weight");
	auto h_b_hadron_pythia_op 	= df_pythia_op_updated.Histo1D<ROOT::VecOps::RVec<float>,double>({"h_b_hadron_muon_dphi_pythia_op","h_b_hadron_muon_dphi_pythia_op",32, -1,1},"b_hadron_muon_dphi","weight");
	auto h_c_hadron_pythia_ss 	= df_pythia_ss_updated.Histo1D<ROOT::VecOps::RVec<float>,double>({"h_c_hadron_muon_dphi_pythia_ss","h_c_hadron_muon_dphi_pythia_ss",32, -1,1},"c_hadron_muon_dphi","weight");
	auto h_c_hadron_pythia_op 	= df_pythia_op_updated.Histo1D<ROOT::VecOps::RVec<float>,double>({"h_c_hadron_muon_dphi_pythia_op","h_c_hadron_muon_dphi_pythia_op",32, -1,1},"c_hadron_muon_dphi","weight");
	auto h_ptavg_powheg_bb_ss = df_powheg_bb_ss_updated.Histo1D({"h_ptavg_powheg_bb_ss","h_ptavg_powheg_bb_ss",50,0,20},"ptavg","weight");
	auto h_ptavg_powheg_bb_op = df_powheg_bb_op_updated.Histo1D({"h_ptavg_powheg_bb_op","h_ptavg_powheg_bb_op",50,0,20},"ptavg","weight");
	auto h_ptavg_powheg_cc_ss = df_powheg_cc_ss_updated.Histo1D({"h_ptavg_powheg_cc_ss","h_ptavg_powheg_cc_ss",50,0,20},"ptavg","weight");
	auto h_ptavg_powheg_cc_op = df_powheg_cc_op_updated.Histo1D({"h_ptavg_powheg_cc_op","h_ptavg_powheg_cc_op",50,0,20},"ptavg","weight");
	auto h_ptavg_pythia_ss = df_pythia_ss_updated.Histo1D({"h_ptavg_pythia_ss","h_ptavg_pythia_ss",50,0,20},"ptavg","weight");
	auto h_ptavg_pythia_op = df_pythia_op_updated.Histo1D({"h_ptavg_pythia_op","h_ptavg_pythia_op",50,0,20},"ptavg","weight");
	
	TFile * m_outfile=new TFile("/usatlas/u/yuhanguo/workarea/dimuon_codes/plots/hadron_muon_dphi_plots.root", "recreate");
	h_b_hadron_powheg_bb_ss->Write();
	h_b_hadron_powheg_bb_op->Write();
	h_c_hadron_powheg_cc_ss->Write();
	h_c_hadron_powheg_cc_op->Write();
	h_b_hadron_pythia_ss->Write();
	h_b_hadron_pythia_op->Write();
	h_c_hadron_pythia_ss->Write();
	h_c_hadron_pythia_op->Write();
	h_ptavg_powheg_bb_ss->Write();
	h_ptavg_powheg_bb_op->Write();
	h_ptavg_powheg_cc_ss->Write();
	h_ptavg_powheg_cc_op->Write();
	h_ptavg_pythia_ss->Write();
	h_ptavg_pythia_op->Write();
	delete m_outfile;
}

