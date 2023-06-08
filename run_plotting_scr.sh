root -b -l << EOF
	.L MuonPairPlotting.c
	// MuonPairPlotting pbpb1
	// pbpb1.isScram = false;
	// pbpb1.Run()

	MuonPairPlotting pbpb2;
	pbpb2.isScram = true;
	pbpb2.Run();
	.q;
EOF

