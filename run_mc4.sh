root -b -l << EOF
	.L MCNTupleFirstPass.c

	MCNTupleFirstPass* cc = new MCNTupleFirstPass();
	cc->is_full_sample = true;
	cc->full_sample_batch_num = 4;
	mc_mode = "mc_truth_cc"
	cc->Run()
	delete cc

	MCNTupleFirstPass* bb = new MCNTupleFirstPass();
	bb->is_full_sample = true;
	bb->full_sample_batch_num = 4;
	bb->mc_mode = "mc_truth_bb"
	bb->Run();
	delete bb
		
	.q;
EOF

