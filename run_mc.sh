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

# root -b -l << EOF
# 	.L MuonPairPlottingPP.c
	
# 	MuonPairPlottingPP cc
# 	cc.isScram = false
# 	cc.isMCTruthBB = false
# 	cc.isMCTruthCC = true
# 	cc.Run()

# 	MuonPairPlottingPP bb;
# 	bb.isScram = false
# 	bb.isMCTruthBB = true
# 	bb.isMCTruthCC = false
# 	bb.Run();
# 	.q;
# EOF

