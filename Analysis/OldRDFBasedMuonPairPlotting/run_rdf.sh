root -b -l << EOF
	.L RDFBasedMuonPairPlotting.c
	// cout << "Have loaded the code." << endl;

	// RDFBasedMuonPairPlottingPbPb rdf_pbpb
	// rdf_pbpb.Run()

	// RDFBasedMuonPairPlottingPP rdf_pp
	// rdf_pp.Run()

	RDFBasedMuonPairPlottingPowheg rdf_powheg_bb
	rdf_powheg_bb.powheg_mode = RDFBasedMuonPairPlottingPowheg::PowhegMode::bb
	rdf_powheg_bb.Run()

	// RDFBasedMuonPairPlottingPowheg rdf_powheg_cc
	// rdf_powheg_cc.powheg_mode = RDFBasedMuonPairPlottingPowheg::PowhegMode::cc
	// rdf_powheg_cc.Run()

	// RDFBasedMuonPairPlottingPythia rdf_pythia
	// rdf_pythia.Run()
	.q

EOF
