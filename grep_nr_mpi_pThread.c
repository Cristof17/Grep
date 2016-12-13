#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <pthread.h>

#define TRUE 1
#define FALSE 0

#define THREADS 4

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
			if (found_pos[start_t] == -1){
				//printf("%d Found one pattern\n", thread_id);
				found_pos[start_t] = thread_id;
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

void print(short *found, int size){
	int i =0;
	for (i = 0; i < size; ++i){
		printf("%d", found[i]);
	}
	printf("\n");
}

void init(short *v, int size){
	int i= 0;
	for(i = 0; i < size; ++i){
		v[i] = -1;
	}
}

int main(int argc, char **argv){
	if (argc < 2){
		printf("grep <pattern> <text>\n");
		return 1;
	}
	/*
	 * openMPI flavour
	 *
	 */
	int rank;
	int numtasks;
	int rc;
	int id;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	/*
	 * pThread flavour
	 *
	 */
	int idp = 0;
	pthread_t threads[THREADS];
	pthread_t threads2[THREADS];
	params prms[THREADS];
	pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

	if (rank == 0){
		if (numtasks < 2 ) {
			printf("Need at least two MPI tasks. Quitting...\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(1);
		}
	}
		
	/*
	 * P = pattern searched
	 * T = text in which to find the pattern
	 */
	char *p, *t;
	
	p = argv[2];
	t = argv[1];
	int start_p; //current character in the pattern examined
	int start_t; //current character in the text examined
	int stop_p;
	int stop_t;
	int chunk;
	int remainder;
	short *found;
	int len_t; //send the length of t
	
	if (rank == 0){
		/*
		 * Init-tot 
		 * 
		 */
		start_t = 0;
		start_p = 0;
		stop_p = strlen(argv[2]);
		stop_t = strlen(argv[1]);
		len_t = strlen(argv[1]);
		chunk = stop_t/numtasks; //how much each thread will be processing
		remainder = stop_t %numtasks; //the last thread gets the uneven part
		//found contains TRUE on each start position of pattern p in text t
		found = (short*)calloc(stop_t, sizeof(short));
		init(found, stop_t);
		/*
		 * Send chunks
		 */
		for (id = 1; id < numtasks; ++id)
		{
			start_p = 0; //omp does not use the previous value, except firstprivate
			start_t = id * chunk;
			stop_t = (id + 1) * chunk -1;
			if (id == numtasks-1)
				stop_t += remainder;//give the last chunk to the last array
			/*
			 * Send integers
			 */
			MPI_Send(&start_t, 1, MPI_INT, id, 0, MPI_COMM_WORLD);
			MPI_Send(&stop_t, 1, MPI_INT, id, 0, MPI_COMM_WORLD);
			MPI_Send(&start_p, 1, MPI_INT, id, 0, MPI_COMM_WORLD);
			MPI_Send(&stop_p, 1, MPI_INT, id, 0, MPI_COMM_WORLD);
			MPI_Send(&len_t, 1, MPI_INT, id, 0, MPI_COMM_WORLD);
			MPI_Send(&p[0], stop_p-start_p +1, MPI_CHAR, id, 0, MPI_COMM_WORLD);
			MPI_Send(&t[start_t], stop_t-start_t + 1, MPI_CHAR, id, 0, MPI_COMM_WORLD);
		}
		/*
 		 * Reset values
 		 */
		start_t = 0;
		start_p = 0;
		stop_p = strlen(argv[2]);
		stop_t = chunk - 1;
		/*
		 * Debug
		 */
		//printf("Task %d received: start_t:%d, stop_t:%d, start_p:%d, stop_p:%d\n", rank, start_t, stop_t, start_p, stop_p);
		//puts(t);

		/*
		 * Process
		 */
	int NR_THREADS = THREADS;
	int chunk_size = stop_t/NR_THREADS; //how much each thread will be processing
	int rest = stop_t %NR_THREADS; //the last thread gets the uneven part
	//found contains TRUE on each start position of pattern p in text t

	for (idp = 0; idp < NR_THREADS; ++idp)
	{
		start_p = 0;
		start_t += idp * chunk_size;
		stop_t += ((idp + 1) * chunk_size) -1;
		if (idp == NR_THREADS-1)
			stop_t += rest;//give the last chunk to the last array
		//go find the pattern and flag the beginning of the patterns found
		//in the found array
		prms[idp].t = t;
		prms[idp].p = p;
		prms[idp].start_t = start_t;
		prms[idp].stop_t = stop_t;
		prms[idp].start_p = start_p;
		prms[idp].stop_p = stop_p;
		prms[idp].found_pos = found;
		prms[idp].thread_id = rank;
		prms[idp].mutex = &mutex;
		int rc = pthread_create(&threads[idp], NULL, process_text, (void *)&prms[idp]);
	}

	for(idp = 0; idp < NR_THREADS; ++idp){
		int rc = pthread_join(threads[idp], NULL);
	}
	
	pthread_mutex_destroy(&mutex);

	pthread_mutex_init(&mutex, NULL);
	/*
	 * If the splitting of the chunks was in the middle of the pattern p
	 * in a text, it is needed to check the (-pattern_size-1;+pattern_size-1)
	 * region of each point of split to see if we find any more patterns
	 */
	 
	
	for (idp = 0; idp < NR_THREADS; ++idp)
	{
		//when searching at the border between two thread chunks, the 
		//last thread has no neighbor so there is no other chunk to check
		//for fragmentation except the chunk the second to last thread is
		//processing
		int last_thread_id = NR_THREADS-1;
		if (idp != NR_THREADS-1){
			start_p = 0;
			start_t = (idp+ 1) * chunk_size;
			stop_t = (idp + 1) * chunk_size;
			//look around the delimitation between two chunks
			start_t -= stop_p;
			stop_t += stop_p-1;
			//go look for patterns
			prms[idp].t = t;
			prms[idp].p = p;
			prms[idp].start_t = start_t;
			prms[idp].stop_t = stop_t;
			prms[idp].start_p = start_p;
			prms[idp].stop_p = stop_p;
			prms[idp].found_pos = found;
			prms[idp].thread_id = rank;
			prms[idp].mutex = &mutex;
			int rc = pthread_create(&threads[idp], NULL, process_text, (void *)&prms[idp]);
		}
	}

	for(idp = 0; idp < NR_THREADS; ++idp){
		//NR_THREADS-1 thread did not start so there is no need to stop it
		//it will get blocked
		if (idp == NR_THREADS-1)
			continue;
		int rc = pthread_join(threads[idp], NULL);
	}
	
	pthread_mutex_destroy(&mutex);
		
		/*
		 * Receive
		 */
		for (id = 1; id < numtasks; ++id){
			start_t = id * chunk;
			stop_t = (id + 1) * chunk;
			if (id == numtasks-1)
				stop_t += remainder;
			MPI_Recv(&found[start_t], stop_t-start_t +1, MPI_SHORT, id, 0, MPI_COMM_WORLD, NULL); 
		}
	//	print(found, stop_t);
	}
	else{
		//receive all the parameters
		/*
		 *  Receive
		 */
		MPI_Recv(&start_t, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
		MPI_Recv(&stop_t, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
		MPI_Recv(&start_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
		MPI_Recv(&stop_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
		MPI_Recv(&len_t, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
		/*
		 * Init array
		 */
		found = (short*)calloc(stop_t, sizeof(short));
		init(found, stop_t);
		p =(char*)malloc(((stop_p-start_p + 1) * sizeof(char)) +1);
		t = (char*)malloc((len_t * sizeof(char))) ;
		MPI_Recv(&p[0], stop_p-start_p +1, MPI_CHAR, 0, 0,MPI_COMM_WORLD, NULL); 
		MPI_Recv(&t[start_t], stop_t-start_t+1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, NULL);
		/*
 		 * Debug stuff
 		 */
		//printf("Task %d received: start_t:%d, stop_t:%d, start_p:%d, stop_p:%d\n", rank, start_t, stop_t, start_p, stop_p);
		//puts(t);
		/*
		 * Process
		 */
	int NR_THREADS = THREADS;
	int chunk_size = stop_t/NR_THREADS; //how much each thread will be processing
	int rest = stop_t %NR_THREADS; //the last thread gets the uneven part
	//found contains TRUE on each start position of pattern p in text t

	for (idp = 0; idp < NR_THREADS; ++idp)
	{
		start_p = 0;
		start_t += idp * chunk_size;
		stop_t += ((idp + 1) * chunk_size) -1;
		if (idp == NR_THREADS-1)
			stop_t += rest;//give the last chunk to the last array
		//go find the pattern and flag the beginning of the patterns found
		//in the found array
		prms[idp].t = t;
		prms[idp].p = p;
		prms[idp].start_t = start_t;
		prms[idp].stop_t = stop_t;
		prms[idp].start_p = start_p;
		prms[idp].stop_p = stop_p;
		prms[idp].found_pos = found;
		prms[idp].thread_id = rank;
		prms[idp].mutex = &mutex;
		int rc = pthread_create(&threads[idp], NULL, process_text, (void *)&prms[idp]);
	}

	for(idp = 0; idp < NR_THREADS; ++idp){
		int rc = pthread_join(threads[idp], NULL);
	}
	
	pthread_mutex_destroy(&mutex);

	pthread_mutex_init(&mutex, NULL);
	/*
	 * If the splitting of the chunks was in the middle of the pattern p
	 * in a text, it is needed to check the (-pattern_size-1;+pattern_size-1)
	 * region of each point of split to see if we find any more patterns
	 */
	 
	
	for (idp = 0; idp < NR_THREADS; ++idp)
	{
		//when searching at the border between two thread chunks, the 
		//last thread has no neighbor so there is no other chunk to check
		//for fragmentation except the chunk the second to last thread is
		//processing
		int last_thread_id = NR_THREADS-1;
		if (idp != NR_THREADS-1){
			start_p = 0;
			start_t = (idp+ 1) * chunk_size;
			stop_t = (idp + 1) * chunk_size;
			//look around the delimitation between two chunks
			start_t -= stop_p;
			stop_t += stop_p-1;
			//go look for patterns
			prms[idp].t = t;
			prms[idp].p = p;
			prms[idp].start_t = start_t;
			prms[idp].stop_t = stop_t;
			prms[idp].start_p = start_p;
			prms[idp].stop_p = stop_p;
			prms[idp].found_pos = found;
			prms[idp].thread_id = rank;
			prms[idp].mutex = &mutex;
			int rc = pthread_create(&threads[idp], NULL, process_text, (void *)&prms[idp]);
		}
	}

	for(idp = 0; idp < NR_THREADS; ++idp){
		//NR_THREADS-1 thread did not start so there is no need to stop it
		//it will get blocked
		if (idp == NR_THREADS-1)
			continue;
		int rc = pthread_join(threads[idp], NULL);
	}
	
	pthread_mutex_destroy(&mutex);		
		/*
		 * Send
		 */
		MPI_Send(found, stop_t - start_t +1, MPI_SHORT, 0, 0, MPI_COMM_WORLD);
		free(found);
	}

	MPI_Barrier(MPI_COMM_WORLD);
		
	/*
	 * If the splitting of the chunks was in the middle of the pattern p
	 * in a text, it is needed to check the (-pattern_size-1;+pattern_size-1)
	 * region of each point of split to see if we find any more patterns
	 */

	if (rank == 0){
 		//Send 
		for (id = 1; id < numtasks; ++id){
			if (id != numtasks-1){
				start_t = ((id + 1) * chunk)-1;
				stop_t = ((id + 1) * chunk);
				start_t -= (stop_p);
				stop_t += (stop_p);
				MPI_Send(&start_t, 1, MPI_INT, id, 0, MPI_COMM_WORLD);
				MPI_Send(&stop_t, 1, MPI_INT, id, 0, MPI_COMM_WORLD);
				MPI_Send(&t[start_t], stop_t-start_t+1, MPI_CHAR, id, 0, MPI_COMM_WORLD);
			}
		}
 		//Init
		start_t = chunk-1;
		stop_t = chunk;
		start_t -= (stop_p);
		stop_t += (stop_p); 
 		//Process
	int NR_THREADS = THREADS;
	int chunk_size = stop_t/NR_THREADS; //how much each thread will be processing
	int rest = stop_t %NR_THREADS; //the last thread gets the uneven part
	//found contains TRUE on each start position of pattern p in text t

	for (idp = 0; idp < NR_THREADS; ++idp)
	{
		start_p = 0;
		start_t += idp * chunk_size;
		stop_t += ((idp + 1) * chunk_size) -1;
		if (idp == NR_THREADS-1)
			stop_t += rest;//give the last chunk to the last array
		//go find the pattern and flag the beginning of the patterns found
		//in the found array
		prms[idp].t = t;
		prms[idp].p = p;
		prms[idp].start_t = start_t;
		prms[idp].stop_t = stop_t;
		prms[idp].start_p = start_p;
		prms[idp].stop_p = stop_p;
		prms[idp].found_pos = found;
		prms[idp].thread_id = rank;
		prms[idp].mutex = &mutex;
		int rc = pthread_create(&threads[idp], NULL, process_text, (void *)&prms[idp]);
	}

	for(idp = 0; idp < NR_THREADS; ++idp){
		int rc = pthread_join(threads[idp], NULL);
	}
	
	pthread_mutex_destroy(&mutex);

	pthread_mutex_init(&mutex, NULL);
	/*
	 * If the splitting of the chunks was in the middle of the pattern p
	 * in a text, it is needed to check the (-pattern_size-1;+pattern_size-1)
	 * region of each point of split to see if we find any more patterns
	 */
	 
	
	for (idp = 0; idp < NR_THREADS; ++idp)
	{
		//when searching at the border between two thread chunks, the 
		//last thread has no neighbor so there is no other chunk to check
		//for fragmentation except the chunk the second to last thread is
		//processing
		int last_thread_id = NR_THREADS-1;
		if (idp != NR_THREADS-1){
			start_p = 0;
			start_t = (idp+ 1) * chunk_size;
			stop_t = (idp + 1) * chunk_size;
			//look around the delimitation between two chunks
			start_t -= stop_p;
			stop_t += stop_p-1;
			//go look for patterns
			prms[idp].t = t;
			prms[idp].p = p;
			prms[idp].start_t = start_t;
			prms[idp].stop_t = stop_t;
			prms[idp].start_p = start_p;
			prms[idp].stop_p = stop_p;
			prms[idp].found_pos = found;
			prms[idp].thread_id = rank;
			prms[idp].mutex = &mutex;
			int rc = pthread_create(&threads[idp], NULL, process_text, (void *)&prms[idp]);
		}
	}

	for(idp = 0; idp < NR_THREADS; ++idp){
		//NR_THREADS-1 thread did not start so there is no need to stop it
		//it will get blocked
		if (idp == NR_THREADS-1)
			continue;
		int rc = pthread_join(threads[idp], NULL);
	}
	
	pthread_mutex_destroy(&mutex);		
		//Receive
		for (id = 1; id < numtasks; ++id){
			if (id != numtasks-1){
				start_t = ((id + 1)* chunk)-1;
				stop_t = ((id + 1)* chunk);
				start_t -= (stop_p);
				stop_t += (stop_p);
				MPI_Recv(&found[start_t -1], stop_t-start_t+1, MPI_SHORT, id, 0, MPI_COMM_WORLD, NULL);
			}
		}
	}
	else if (rank != numtasks-1 && rank != 0){
		//Debug
		//puts(t);
 		//Receive
		MPI_Recv(&start_t, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
		MPI_Recv(&stop_t, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
		MPI_Recv(&t[start_t], stop_t-start_t+1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, NULL);
 		//Init
		found = (short*)calloc(stop_t, sizeof(short));
		init(found, stop_t);
 		//Process
	int NR_THREADS = THREADS;
	int chunk_size = stop_t/NR_THREADS; //how much each thread will be processing
	int rest = stop_t %NR_THREADS; //the last thread gets the uneven part
	//found contains TRUE on each start position of pattern p in text t

	for (idp = 0; idp < NR_THREADS; ++idp)
	{
		start_p = 0;
		start_t += idp * chunk_size;
		stop_t += ((idp + 1) * chunk_size) -1;
		if (idp == NR_THREADS-1)
			stop_t += rest;//give the last chunk to the last array
		//go find the pattern and flag the beginning of the patterns found
		//in the found array
		prms[idp].t = t;
		prms[idp].p = p;
		prms[idp].start_t = start_t;
		prms[idp].stop_t = stop_t;
		prms[idp].start_p = start_p;
		prms[idp].stop_p = stop_p;
		prms[idp].found_pos = found;
		prms[idp].thread_id = rank;
		prms[idp].mutex = &mutex;
		int rc = pthread_create(&threads[idp], NULL, process_text, (void *)&prms[idp]);
	}

	for(idp = 0; idp < NR_THREADS; ++idp){
		int rc = pthread_join(threads[idp], NULL);
	}
	
	pthread_mutex_destroy(&mutex);

	pthread_mutex_init(&mutex, NULL);
	/*
	 * If the splitting of the chunks was in the middle of the pattern p
	 * in a text, it is needed to check the (-pattern_size-1;+pattern_size-1)
	 * region of each point of split to see if we find any more patterns
	 */
	 
	
	for (idp = 0; idp < NR_THREADS; ++idp)
	{
		//when searching at the border between two thread chunks, the 
		//last thread has no neighbor so there is no other chunk to check
		//for fragmentation except the chunk the second to last thread is
		//processing
		int last_thread_id = NR_THREADS-1;
		if (idp != NR_THREADS-1){
			start_p = 0;
			start_t = (idp+ 1) * chunk_size;
			stop_t = (idp + 1) * chunk_size;
			//look around the delimitation between two chunks
			start_t -= stop_p;
			stop_t += stop_p-1;
			//go look for patterns
			prms[idp].t = t;
			prms[idp].p = p;
			prms[idp].start_t = start_t;
			prms[idp].stop_t = stop_t;
			prms[idp].start_p = start_p;
			prms[idp].stop_p = stop_p;
			prms[idp].found_pos = found;
			prms[idp].thread_id = rank;
			prms[idp].mutex = &mutex;
			int rc = pthread_create(&threads[idp], NULL, process_text, (void *)&prms[idp]);
		}
	}

	for(idp = 0; idp < NR_THREADS; ++idp){
		//NR_THREADS-1 thread did not start so there is no need to stop it
		//it will get blocked
		if (idp == NR_THREADS-1)
			continue;
		int rc = pthread_join(threads[idp], NULL);
	}
	
	pthread_mutex_destroy(&mutex);		
 		//Send
		MPI_Send(&found[start_t-1], stop_t-start_t+1, MPI_SHORT, 0, 0, MPI_COMM_WORLD);
	}

 	//Debug
 	
	if (rank == 0){
	
		int i = 0;
		for (i = 0; i < len_t; ++i){
			if (found[i] != -1)
				printf("%d found pattern\n", found[i]);
		}
		printf("\n");

		/*
		for (i = 0; i < len_t; ++i){
			printf("%d", found[i]);
		}
		printf("\n");
		*/
	}

	MPI_Finalize();

	return 0;
}
