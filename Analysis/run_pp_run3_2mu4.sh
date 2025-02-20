root -b -l << EOF
	.L MuonNTupleFirstPassPP.c

	MuonNTupleFirstPassPP p1;
	p1.run3_use_mu4_mu4_noL1 = false;
	p1.run3_file_batch = 1;
	p1.Run()
	
	.q;

EOF

root -b -l << EOF
	.L MuonNTupleFirstPassPP.c

	MuonNTupleFirstPassPP p2;
	p2.run3_use_mu4_mu4_noL1 = false;
	p2.run3_file_batch = 2;
	p2.Run()

	.q;

EOF

root -b -l << EOF
	.L MuonNTupleFirstPassPP.c

	MuonNTupleFirstPassPP p3;
	p3.run3_use_mu4_mu4_noL1 = false;
	p3.run3_file_batch = 3;
	p3.Run()
	
	.q;

EOF

root -b -l << EOF
	.L MuonNTupleFirstPassPP.c

	MuonNTupleFirstPassPP p4;
	p4.run3_use_mu4_mu4_noL1 = false;
	p4.run3_file_batch = 4;
	p4.Run()

	.q;

EOF