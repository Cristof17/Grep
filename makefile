TEXT = CCTTTTGCGCTTCTGCTACCTTTTGCGCGCGCGCGGAACCTTTTGCCCTTTTGC
PATTERN = CCTTTTGC

build:
	gcc grep.c -o grep -g
run:
	./grep $(TEXT) $(PATTERN) 
debug:
	gdb -tui --args ./grep $(TEXT) $(PATTERN)
