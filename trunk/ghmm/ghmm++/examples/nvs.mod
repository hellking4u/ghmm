SHMM = {
	M = 1;
	N = 2;
	cos = 1;
	density = 0;
	fix_state = vector {
	0,0;
	};
	Pi = vector {
	0.2, 0.8;
	};
#	A = matrix {random} ;
	A = matrix {
	0.2, 0.8;
	0.2, 0.8;
	};

	# Gewichtsmatrix der Verteilungen, NxM
	C = matrix {
		1;
		1;
	};
	# Mittelwerte mue, NxM
	Mue = matrix {
		0;
		0.2;
	};
	# Varianz u, NxM
	U = matrix {
		2;
		2;
	};			

};