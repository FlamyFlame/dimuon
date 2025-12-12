root -b -l << EOF
	.L MuonPairPlottingPP.c

	// MuonPairPlottingPP pp1
	// pp1.Run()

	MuonPairPlottingPP pp2;
	pp2.isScram = true;
	pp2.Run();

	// MuonPairPlottingPP pp3;
	// pp3.isTight = true;
	// pp3.Run();

	.q;
EOF

