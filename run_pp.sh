root -b -l << EOF
	.L MuonNTupleFirstPassPP.c

	MuonNTupleFirstPassPP mmedium;
	mmedium.isTight = false;
	mmedium.Run();

	// MuonNTupleFirstPassPP mtight;
	// mtight.isTight = true;
	// mtight.Run();
	.q;

EOF

