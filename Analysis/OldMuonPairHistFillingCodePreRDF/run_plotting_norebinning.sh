root -b -l << EOF
	.L MuonPairPlottingNoRebinning.c
	MuonPairPlottingNoRebinning pp1;
	pp1.isScram = false;
	pp1.Run();
	MuonPairPlottingNoRebinning pp2;
	pp2.isScram = true;
	pp2.Run();
	.q;
EOF

