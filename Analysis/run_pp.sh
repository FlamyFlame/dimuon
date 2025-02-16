# root -b -l << EOF
# 	.L MuonNTupleFirstPassPP.c

# 	MuonNTupleFirstPassPP mmedium;
# 	mmedium.isTight = false;
# 	mmedium.mode = 4;
# 	mmedium.Run();

# 	MuonNTupleFirstPassPP mtight;
# 	mtight.isTight = true;
# 	mtight.mode = 4;
# 	mtight.Run();
	
	# .q;

# EOF


# root -b -l << EOF
# 	.L MuonNTupleFirstPassPP.c

# 	MuonNTupleFirstPassPP p1;
# 	p1.run3_file_batch = 1;
# 	p1.Run()
	
# 	MuonNTupleFirstPassPP p2;
# 	p2.run3_file_batch = 2;
# 	p2.Run()
	
# 	MuonNTupleFirstPassPP p3;
# 	p3.run3_file_batch = 3;
# 	p3.Run()
	
# 	MuonNTupleFirstPassPP p4;
# 	p4.run3_file_batch = 4;
# 	p4.Run()
#	.q;

# EOF


root -b -l << EOF
	.L MuonNTupleFirstPassPP.c

	MuonNTupleFirstPassPP p1;
	p1.run3_file_batch = 1;
	p1.Run()
	
	.q;

EOF

root -b -l << EOF
	.L MuonNTupleFirstPassPP.c

	MuonNTupleFirstPassPP p2;
	p2.run3_file_batch = 2;
	p2.Run()

	.q;

EOF

root -b -l << EOF
	.L MuonNTupleFirstPassPP.c

	MuonNTupleFirstPassPP p3;
	p3.run3_file_batch = 3;
	p3.Run()
	
	.q;

EOF

root -b -l << EOF
	.L MuonNTupleFirstPassPP.c

	MuonNTupleFirstPassPP p4;
	p4.run3_file_batch = 4;
	p4.Run()

	.q;

EOF