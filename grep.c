#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define TRUE 1
#define FALSe 0

void print(char *v){
	puts(v);
}

/*
 * Function how_many_positions
 * ---------------------
 * Returns the distance from pointer p to the character c in the preffix of v
 * The preffix of v is the followind: 
 * Say v is an array of size N and at any give position k in the v
 * the preffix is v[0..k] and the suffix is v[k+1..N]
 */
int how_many_positions(char *p, char searched, int pattern_size){
	int nr_pos=0;
	char *curr = p;
	int i = 0;
	while(*curr != searched && i < pattern_size){
		nr_pos++;
		curr--;
		i++;
	}
	/*
	 * If the letter searched does not fit in our pattern, 
	 * unlike the algorithm where you move the pattern to the right,
	 * we move the text to the left, changing the missalligned letter
 	 */
	if (nr_pos == pattern_size)
		return 1;
	return nr_pos;
}

int main(int argc, char **argv){
	if (argc < 2){
		printf("grep <pattern> <text>\n");
		return 1;
	}

	char *p, *t;
	int found = 0;

	p = argv[2];
	t = argv[1];

	long processed=0;
	long total_processed=0;
	t += strlen(p)-1;
	p += strlen(p)-1;
	long text_len = strlen(argv[1]);
	
	clock_t start = clock();
	while(total_processed != text_len){
		if (strlen(p) == strlen(argv[2])){
		   found = TRUE; //we found a pattern
		   //go find the next one
		   printf("Found one pattern\n");
		   t += strlen(argv[2]);
		}
		//if we found the pattern and the pattern is at the begining
		//skip the text, don't stop, search for more
		if (found && (strlen(t) == strlen(argv[1]))){
			t += strlen(argv[2]);
		}
		if (*t == *p){
			t--;
			p--;
			processed++;
		}
		else {
			if (strlen(t) < strlen(argv[2]))
				break;
			//shift t to the left				
			int pattern_length = strlen(argv[2]);
			int num = how_many_positions(p, *t, pattern_length);
			t+= processed;
			t+= num;
			total_processed += processed;
			total_processed += num;
			processed = 0;
			p = argv[2];
			p += strlen(p)-1;
		}
	}
	clock_t stop = clock();
	printf("Executed in %f\n", ((float)stop-start)/CLOCKS_PER_SEC);
	return 0;
}
