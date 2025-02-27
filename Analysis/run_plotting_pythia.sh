root -b -l << EOF
	.L MuonPairPlottingPythia.c

	MuonPairPlottingPythia py
	py.Run()
	
	.q;
EOF

