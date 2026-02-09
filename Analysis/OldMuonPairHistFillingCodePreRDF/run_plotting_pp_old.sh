root -b -l << EOF
	.L MuonPairPlottingOLDPP.c

	MuonPairPlottingOLDPP pp1
	pp1.Run()

	// MuonPairPlottingOLDPP pp2;
	// pp2.isScram = true;
	// pp2.Run();

	// MuonPairPlottingOLDPP pp3;
	// pp3.isTight = true;
	// pp3.Run();

	.q;
EOF

