root -b -l << EOF
	.L SimpleOutput.c

	SimpleOutput cc
	cc.mc_mode="mc_truth_cc"
	cc.Run()

	// SimpleOutput bb;
	// bb.mc_mode="mc_truth_bb"
	// bb.Run();
	.q;
EOF


