TEXT = "Robert Robert Robert merge la magazin pentru ca Robert nu are altceva de facut decat sa mearga la magazinut"
PATTERN = "Robert"

build:
	gcc grep.c -o grep -g
	gcc grep_nr.c -o grep -g
	gcc grep_nr_omp.c -o grep -g -fopenmp
debug:
	gdb -tui --args ./grep $(TEXT) $(PATTERN)
run:
	./grep $(TEXT) $(PATTERN)
