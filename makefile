TEXT = "Robert Robert Robert merge la magazin pentru ca Robert nu are altceva de facut decat sa mearga la magazinut"
PATTERN = "Robert"

build:
	gcc grep.c -o grep -g
	gcc grep_nr.c -o grep -g
	gcc grep_nr_omp.c -o grep -g -fopenmp
	gcc grep_nr_pthread.c -o grep -g -lpthread
	mpicc grep_nr_mpi.c -o grep -g
	#mpicc grep_nr_mpi_omp.c -o grep -g -fopenmp
	#mpicc grep_nr_mpi_pThread.c -o grep -g -lpthread
debug:
	gdb -tui --args ./grep $(TEXT) $(PATTERN)
run_omp:
	./grep $(TEXT) $(PATTERN)
run_mpi:
	mpirun -np 2 ./grep $(TEXT) $(PATTERN)
