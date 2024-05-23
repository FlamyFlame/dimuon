root -b -l << EOF
	.L ScrambGenPP.c
	ScrambSampleGenPP pp;
	pp.Run();

	.q;

EOF