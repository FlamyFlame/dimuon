root -b -l << EOF
	.L MuonPairPlottingPP.c

	MuonPairPlottingPP pp1
	pp1.isScram = false;
    pp1.isMCTruthBB = false;
    pp1.isMCTruthCC = true;
	pp1.Run()

	MuonPairPlottingPP pp2;
	pp2.isScram = false;
    pp2.isMCTruthBB = true;
    pp2.isMCTruthCC = false;
	pp2.Run();
	
	.q;
EOF

