TEXT = CCTTTTGCGCTTCTGCTACCTTTTGCGCGCGCGCGGAACCTTTTGCCCTTTTGC
PATTERN = CCTTTTGC

build:
	gcc grep.c -o grep -g
	gcc grep_nr.c -o grep -g
debug_p:
	gdb -tui --args ./grep $(TEXT) $(PATTERN)
