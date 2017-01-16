#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <pthread.h>
#include <time.h>

#define TRUE 1
#define FALSe 0

#define THREADS 32

void print(char *v){
	puts(v);
}

typedef struct params_t{
	char *t;
	char *p;
	int start_t;
	int stop_t;
	int start_p;
	int stop_p;
	int thread_id;
	int padding;
	short *found_pos;
	pthread_mutex_t *mutex;
}params;


/*
 * Function how_many_positions
 * ---------------------------
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

/*
 * Function process_text
 * ---------------------
 * Discovers matching pattern p of size stop_p in the t text of size stop_t
 * and puts TRUE in the start_t position of the found_pos array when a 
 * matching pattern was found (this is necessary for multithreading because
 * there may be multiple which find the same pattern at the same time, and 
 * if they found it there is no need to signal a match because a previous
 * thread did it before
 */

void* process_text(void *par){
	
	params *pmtrs = (params*)par;
	int thread_id = pmtrs->thread_id;
	char *t = pmtrs->t;
	char *p = pmtrs->p;
	int start_t = pmtrs->start_t;
	int stop_t = pmtrs-> stop_t;
	int start_p = pmtrs-> start_p;
	int stop_p = pmtrs-> stop_p;
	short *found_pos = pmtrs-> found_pos;
	pthread_mutex_t *mutex = pmtrs->mutex;
	int found = 0;
	long processed=0;
	long total_processed=0;
	start_p += stop_p-1;
	start_t += stop_p-1;

	while(start_t <= stop_t){
		if (start_p == 0){
			found = TRUE; //we found a pattern and then go find the next one
			//acquire the mutex
			pthread_mutex_lock(mutex);
			if (found_pos[start_t] != TRUE){
				printf("%d Found one pattern\n", thread_id);
				found_pos[start_t] = TRUE;
			}
			//release the mutex
			pthread_mutex_unlock(mutex);
			//reset positions
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
}

int main(int argc, char **argv){
	if (argc < 2){
		printf("grep <pattern> <text>\n");
		return 1;
	}
	
	/*
	 * pThread flavour
	 *
	 */
	int id = 0;
	pthread_t threads[THREADS];
	pthread_t threads2[THREADS];
	params prms[THREADS];
	pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

	/*
	 * P = pattern searched
	 * T = text in which to find the pattern
	 */
	char *p, *t;

	p = argv[2];
	t = argv[1];
	int start_p = 0; //current character in the pattern examined
	int start_t = 0; //current character in the text examined
	int stop_p = strlen(argv[2]);
	int stop_t = strlen(argv[1]);

	int NR_THREADS = THREADS;
	int chunk = stop_t/NR_THREADS; //how much each thread will be processing
	int remainder = stop_t %NR_THREADS; //the last thread gets the uneven part
	//found contains TRUE on each start position of pattern p in text t
	short *found = (short*)calloc(stop_t, sizeof(short));

	clock_t start = clock();
	for (id = 0; id < NR_THREADS; ++id)
	{
		start_p = 0; //omp does not use the previous value, except firstprivate
		start_t = id * chunk;
		stop_t = (id + 1) * chunk -1;
		if (id == NR_THREADS-1)
			stop_t += remainder;//give the last chunk to the last array
		//go find the pattern and flag the beginning of the patterns found
		//in the found array
		prms[id].t = t;
		prms[id].p = p;
		prms[id].start_t = start_t;
		prms[id].stop_t = stop_t;
		prms[id].start_p = start_p;
		prms[id].stop_p = stop_p;
		prms[id].found_pos = found;
		prms[id].thread_id = id;
		prms[id].mutex = &mutex;
		int rc = pthread_create(&threads[id], NULL, process_text, (void *)&prms[id]);
	}

	for(id = 0; id < NR_THREADS; ++id){
		int rc = pthread_join(threads[id], NULL);
	}
	
	pthread_mutex_destroy(&mutex);

	pthread_mutex_init(&mutex, NULL);
	/*
	 * If the splitting of the chunks was in the middle of the pattern p
	 * in a text, it is needed to check the (-pattern_size-1;+pattern_size-1)
	 * region of each point of split to see if we find any more patterns
	 */
	 
	
	for (id = 0; id < NR_THREADS; ++id)
	{
		//when searching at the border between two thread chunks, the 
		//last thread has no neighbor so there is no other chunk to check
		//for fragmentation except the chunk the second to last thread is
		//processing
		int last_thread_id = NR_THREADS-1;
		if (id != NR_THREADS-1){
			start_p = 0;
			start_t = (id + 1) * chunk;
			stop_t = (id + 1) * chunk;
			//look around the delimitation between two chunks
			start_t -= stop_p;
			stop_t += stop_p-1;
			//go look for patterns
			prms[id].t = t;
			prms[id].p = p;
			prms[id].start_t = start_t;
			prms[id].stop_t = stop_t;
			prms[id].start_p = start_p;
			prms[id].stop_p = stop_p;
			prms[id].found_pos = found;
			prms[id].thread_id = id;
			prms[id].mutex = &mutex;
			int rc = pthread_create(&threads[id], NULL, process_text, (void *)&prms[id]);
		}
	}

	for(id = 0; id < NR_THREADS; ++id){
		//NR_THREADS-1 thread did not start so there is no need to stop it
		//it will get blocked
		if (id == NR_THREADS-1)
			continue;
		int rc = pthread_join(threads[id], NULL);
	}
	
	pthread_mutex_destroy(&mutex);
	
	clock_t stop = clock();
	printf("Executed in %f\n",((float)stop-start)/CLOCKS_PER_SEC);
			return 0;
}
