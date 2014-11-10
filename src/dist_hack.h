
#ifdef FAST
#define NAME distMatrix
#define P_OUTER _Pragma("omp parallel for num_threads( THREADS)")
#define P_INNER
#else
#undef NAME
#undef P_OUTER
#undef P_INNER
#define NAME distMatrixLM
#define P_OUTER 
#define P_INNER _Pragma("omp parallel for num_threads( THREADS)")
#endif

double *NAME( seq_t* sequences, size_t n, data_t *M){
	errno = 0;
	double *D = (double*) malloc( n * n * sizeof(double));

	if( !D){
		err(1, "Could not allocate enough memory for the comparison matrix. Try using --join or --low-memory.");
	}

	size_t i;

	//#pragma
	P_OUTER
	for(i=0;i<n;i++){
		esa_t E;

		if( esa_init( &E, &(sequences[i])) != 0){
			warnx("Failed to create index for %s.", sequences[i].name);

			for( size_t j=0; j< n; j++){
				D(i,j) = (i==j) ? 0.0 : NAN;
				if( M) M(i,j).coverage = 0.0;
			}

			continue;
		}

		// now compare every other sequence to i
		size_t j;

		P_INNER
		for(j=0; j<n; j++){
			if( j == i) {
				D(i,j) = 0.0;
				if( M) M(i,j).coverage = 0.0;
				continue;
			}

			// TODO: Provide a nicer progress indicator.
			if( FLAGS & F_EXTRA_VERBOSE ){
				#pragma omp critical
				{
					fprintf( stderr, "comparing %lu and %lu\n", i, j);
				}
			}

			size_t ql = sequences[j].len;

			data_t datum = dist_anchor( &E, sequences[j].S, ql, sequences[i].gc);

			if( !(FLAGS & F_RAW)){
				datum.distance = -0.75 * log(1.0- (4.0 / 3.0) * datum.distance ); // jukes cantor
			}
			// fix negative zero
			if( datum.distance <= 0.0 ){
				datum.distance = 0.0;
			}

			D(i,j) = datum.distance;
			if( M){
				M(i,j) = datum;
			}
		}

		esa_free(&E);
	}

	return D;
}