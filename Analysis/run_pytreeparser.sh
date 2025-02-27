#!/bin/bash
root -b -l << EOF
	.L PythiaNTupleFirstPass.c

	PythiaNTupleFirstPass* py = new PythiaNTupleFirstPass();
	py->Run()
	delete py

	.q;
EOF

