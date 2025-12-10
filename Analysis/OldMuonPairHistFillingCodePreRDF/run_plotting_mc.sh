root -b -l << EOF
	.L MuonPairPlottingPowheg.c
	
	MuonPairPlottingPowheg pp1
    pp1.isMCTruthBB = false;
    pp1.isMCTruthCC = true;
	pp1.Run()

	MuonPairPlottingPowheg pp2;
    pp2.isMCTruthBB = true;
    pp2.isMCTruthCC = false;
	pp2.Run();
	
	.q;
EOF

