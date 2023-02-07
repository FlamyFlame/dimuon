root -b -l << EOF
	.L MCNTupleFirstPass.c

	MCNTupleFirstPass cc
	cc.mc_mode="mc_truth_cc"
	cc.Run()

	MCNTupleFirstPass bb;
	bb.mc_mode="mc_truth_bb"
	bb.Run();
	.q;
EOF

