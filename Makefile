sequential:
	mpicxx -O2 main_seq.cpp iter.cpp mesh.cpp -o sequential
parallel:
	mpicxx -O2 main_par.cpp iter.cpp mesh.cpp -o parallel
parallel_omp:
	mpixlcxx_r -O2 main_par_omp.cpp iter.cpp mesh.cpp -qsmp=omp -o parallel_omp
