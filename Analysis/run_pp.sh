root -b -l << EOF
	.L MuonNTupleFirstPassPP.c

	MuonNTupleFirstPassPP mmedium;
	mmedium.isTight = false;
	mmedium.mode = 4;
	mmedium.Run();

	MuonNTupleFirstPassPP mtight;
	mtight.isTight = true;
	mtight.mode = 4;
	mtight.Run();
	
	.q;

EOF

