root -b -l << EOF
	.L MCNTupleFirstPass.c

	MCNTupleFirstPass* cc = new MCNTupleFirstPass();
	cc->mc_mode = "mc_truth_cc"
	cc->Run()
	delete cc

	MCNTupleFirstPass* bb = new MCNTupleFirstPass();
	bb->mc_mode = "mc_truth_bb"
	bb->Run();
	delete bb
	.q;
EOF

