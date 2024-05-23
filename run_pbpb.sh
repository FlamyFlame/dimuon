root -b -l << EOF
	.L MuonNTupleFirstPass.c

	MuonNTupleFirstPass mmedium;
	mmedium.mode = 4;
	mmedium.Run();

	// MuonNTupleFirstPass mtight;
	// mtight.isTight = true;
	// mtight.Run();
	
	.q;

EOF

