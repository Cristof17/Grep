#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

	/**
 	 * Starting index part
 	 *
 	 */
	int start_p = 0;
	int start_t = 0;
	int stop_p = strlen(argv[2]);
	int stop_t = strlen(argv[1]);

	start_p += strlen(argv[2])-1;
	start_t += strlen(argv[2])-1;
	
	while(start_t < stop_t){
		if (start_p == 0){
		   found = TRUE; //we found a pattern
		   //go find the next one
		   printf("Found one pattern\n");
		   start_t += processed;
		   start_p += processed;
		   start_t += stop_p;
		   processed = 0;
		}

		if (t[start_t] == p[start_p]){
			start_t --;
			start_p --;
			processed++;
		}
		else {
			//if we don't do this, make sure we
			//don't skip the end of the t by
			//adding the lengt of p
			if (stop_t-start_t < stop_p)
				break;
			//shift t to the left				
			int num = how_many_positions(p + (stop_p - processed -1), t[start_t], stop_p);
			start_t+= processed;
			start_t+= num;
			total_processed += processed;
			total_processed += num;
			processed = 0;
			start_p = stop_p-1;
		}
	}
	return 0;
}
