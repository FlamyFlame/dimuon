root -b -l << EOF
	.L MuonPairPlottingOLDPP.c

	MuonPairPlottingOLDPP pp1
	pp1.isScram = false;
    pp1.isMCTruthBB = false;
    pp1.isMCTruthCC = true;
	pp1.Run()

	MuonPairPlottingOLDPP pp2;
	pp2.isScram = false;
    pp2.isMCTruthBB = true;
    pp2.isMCTruthCC = false;
	pp2.Run();
	
	.q;
EOF

