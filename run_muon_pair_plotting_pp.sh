root -b -l << EOF
	.L MuonPairPlottingPP.c
	MuonPairPlottingPP pp1
	pp1.isScram = false;
	pp1.Run()
	MuonPairPlottingPP pp2;
	pp2.isScram = true;
	pp2.Run();
	.q;
EOF

