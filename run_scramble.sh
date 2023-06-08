root -b -l << EOF
	.L ScrambGen.c

	ScrambSampleGen gg;
	gg.keep_spill = true;
	gg.Run();

	// ScrambSampleGen gg_nospill;
	// gg_nospill.keep_spill = false;
	// gg_nospill.Run();
	
	.q;

EOF

