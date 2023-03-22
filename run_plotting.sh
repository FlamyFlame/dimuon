root -b -l << EOF
	.L MuonPairPlotting.c
	MuonPairPlotting pp1
	pp1.isScram = false;
	pp1.Run()

	// MuonPairPlotting pp2;
	// pp2.isScram = true;
	// pp2.Run();
	.q;
EOF

