#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include "structs.h"

void terror(char *s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

Queue * generate_queue(Head * queue_head, uint64_t t_reads, uint64_t n_threads, uint64_t levels){
    uint64_t i, j;
    uint64_t reads_per_thread;
    if(levels > t_reads) levels = t_reads;
    uint64_t pieces = t_reads/levels;
    uint64_t from, to, t_queues = 0, current_queue = 0;
    for(i=0;i<levels;i++) t_queues += ((i+1)*n_threads);
    Queue * queues = (Queue *) malloc(t_queues * sizeof(Queue));
    if(queues == NULL) terror("Could not allocate queue tasks");
    queue_head->head = &queues[0];
    
    for(i=0;i<levels;i++){

        //reads_per_thread = (uint64_t) (floorl((long double) pieces / (long double) ((i+1)*n_threads)));
        reads_per_thread = (uint64_t) (ceill((long double) pieces / (long double) ((i+1)*n_threads)));
        

        for(j=0;j<(i+1)*n_threads;j++){
            from = j * reads_per_thread + (pieces*i);
            to = (j + 1) * reads_per_thread + (pieces*i);
            
            if(j == (i+1)*n_threads - 1) to = pieces*(i+1);


            if(i == levels - 1 && j == (i+1)*n_threads - 1){
                //If its the last 
                queues[current_queue].next = NULL;
            }else{
                //Else add it to the queue
                queues[current_queue].next = &queues[current_queue+1];
            }
            

            queues[current_queue].r1 = from;
            queues[current_queue].r2 = to;
            current_queue++;
            //printf("current_piece: %"PRIu64"-%"PRIu64" diff: %"PRIu64"\n", from, to, to - from);

        }

    }
    //printf("TREADS was %"PRIu64"\n", t_reads);    
    return &queues[0];
}

void print_queue(Queue * q){
    fprintf(stdout, "Task: %"PRIu64"-%"PRIu64"\n", q->r1, q->r2);
}

Queue * get_task_from_queue(Head * queue_head, pthread_mutex_t * lock){
    pthread_mutex_lock(lock);

    Queue * ptr = queue_head->head;
    if(queue_head->head != NULL) queue_head->head = queue_head->head->next;
    //if(ptr != NULL){ printf("Taking "); /*print_queue(ptr);*/ }

    pthread_mutex_unlock(lock);


    return ptr;
}

uint64_t quick_pow4(uint64_t n){
    static uint64_t pow4[33]={1L, 4L, 16L, 64L, 256L, 1024L, 4096L, 16384L, 65536L, 
    262144L, 1048576L, 4194304L, 16777216L, 67108864L, 268435456L, 1073741824L, 4294967296L, 
    17179869184L, 68719476736L, 274877906944L, 1099511627776L, 4398046511104L, 17592186044416L, 
    70368744177664L, 281474976710656L, 1125899906842624L, 4503599627370496L, 18014398509481984L, 
    72057594037927936L, 288230376151711744L, 1152921504606846976L, 4611686018427387904L};
    return pow4[n];
}

uint64_t quick_pow4byLetter(uint64_t n, const char c){
    static uint64_t pow4_G[33]={2*1L, 2*4L, 2*16L, 2*64L, 2*256L, 2*1024L, 2*4096L, 2*16384L, 2*65536L, 
    (uint64_t)2*262144L, (uint64_t)2*1048576L,(uint64_t)2*4194304L, (uint64_t)2*16777216L, (uint64_t)2*67108864L, (uint64_t)2*268435456L, (uint64_t)2*1073741824L, (uint64_t)2*4294967296L, 
    (uint64_t)2*17179869184L, (uint64_t)2*68719476736L, (uint64_t)2*274877906944L, (uint64_t)2*1099511627776L, (uint64_t)2*4398046511104L, (uint64_t)2*17592186044416L, 
    (uint64_t)2*70368744177664L, (uint64_t)2*281474976710656L, (uint64_t)2*1125899906842624L, (uint64_t)2*4503599627370496L, (uint64_t)2*18014398509481984L, 
    (uint64_t)2*72057594037927936L, (uint64_t) 2*288230376151711744L, (uint64_t) 2*1152921504606846976L, (uint64_t) 2*4611686018427387904L};
    
    static uint64_t pow4_T[33]={3*1L, 3*4L, 3*16L, 3*64L, 3*256L, 3*1024L, 3*4096L, 3*16384L, 3*65536L, 
    (uint64_t)3*262144L, (uint64_t) 3*1048576L, (uint64_t)3*4194304L, (uint64_t)3*16777216L, (uint64_t)3*67108864L, (uint64_t)3*268435456L, (uint64_t)3*1073741824L, (uint64_t)3*4294967296L, 
    (uint64_t)3*17179869184L, (uint64_t)3*68719476736L, (uint64_t)3*274877906944L, (uint64_t)3*1099511627776L, (uint64_t)3*4398046511104L, (uint64_t)3*17592186044416L, 
    (uint64_t)3*70368744177664L, (uint64_t)3*281474976710656L, (uint64_t)3*1125899906842624L, (uint64_t)3*4503599627370496L, (uint64_t)3*18014398509481984L, 
    (uint64_t)3*72057594037927936L, (uint64_t) 3*288230376151711744L, (uint64_t) 3*1152921504606846976L, (uint64_t) 3*4611686018427387904L};
    
    if(c == 'A') return 0;
    if(c == 'C') return quick_pow4(n);
    if(c == 'G') return pow4_G[n];
    if(c == 'T') return pow4_T[n];
    return 0;
}

uint64_t hashOfWord(const unsigned char * word, uint32_t k, uint64_t offset){
    
    uint64_t value = 0, jIdx;
    for(jIdx=0;jIdx<k-offset;jIdx++){
        value += quick_pow4byLetter(k-(jIdx+1), (char) word[jIdx]);
    }
    return value;
}

void perfect_hash_to_word(unsigned char * word, uint64_t hash, uint32_t k){
    uint64_t jIdx = (uint64_t) (k-1), upIdx = 0;
    uint64_t v;
    while(jIdx >= 0){
        v = (uint64_t) floor(hash / (pow(4, jIdx)));
        if(v == 0){ word[upIdx++] = (unsigned char) 'A'; hash -= ((uint64_t) pow(4, jIdx) * 0); }
        if(v == 1){ word[upIdx++] = (unsigned char) 'C'; hash -= ((uint64_t) pow(4, jIdx) * 1); }
        if(v == 2){ word[upIdx++] = (unsigned char) 'G'; hash -= ((uint64_t) pow(4, jIdx) * 2); }
        if(v == 3){ word[upIdx++] = (unsigned char) 'T'; hash -= ((uint64_t) pow(4, jIdx) * 3); }
        
        if(jIdx == 0) break;
        --jIdx;
    }
}

uint64_t collisioned_hash(const unsigned char * word, uint32_t k){
    uint64_t jIdx, value = 0;
    for(jIdx=0;jIdx<k;jIdx+=BYTES_IN_MER){
        value += quick_pow4byLetter(k - (jIdx+1), (char) word[jIdx]);
    }
    return value;
}

void decomposed_hash_of_word(const unsigned char * word, unsigned char * vector, uint32_t k){
    uint64_t jIdx;
    for(jIdx=0;jIdx<k/BYTES_IN_MER;jIdx++){
        vector[jIdx] = (unsigned char) quick_pow4byLetter(0, (char) word[jIdx*BYTES_IN_MER]);
        vector[jIdx] += (unsigned char) quick_pow4byLetter(1, (char) word[jIdx*BYTES_IN_MER+1]);
        vector[jIdx] += (unsigned char) quick_pow4byLetter(2, (char) word[jIdx*BYTES_IN_MER+2]);
        vector[jIdx] += (unsigned char) quick_pow4byLetter(3, (char) word[jIdx*BYTES_IN_MER+3]);
    }
}

uint64_t xor_decomposed_hash(unsigned char * vector1, unsigned char * vector2, uint32_t k){
    uint64_t jIdx;
    uint64_t diffs = 0;
    for(jIdx=0;jIdx<k/BYTES_IN_MER;jIdx++){
        diffs += (uint64_t) ((vector1[jIdx] ^ vector2[jIdx]) & 0x1);
    }
    return diffs;
}

uint64_t asciiToUint64(const char *text){
	uint64_t number=0;

	for(;*text;text++){
		char digit=*text-'0';           
		number=(number*10)+digit;
	}
	return number;
}

unsigned char complement(unsigned char c){
    
    switch(c){
    case ((unsigned)'A'): return (unsigned)'T';
    break;
    case ((unsigned)'C'): return (unsigned)'G';
    break;
    case ((unsigned)'G'): return (unsigned)'C';
    break;
    case ((unsigned)'T'): return (unsigned)'A';
    break;
    case ((unsigned)'N'): return (unsigned)'N';
    break;
    case ((unsigned)'-'): return (unsigned)'-'; // Gaps
    break;
    case ((unsigned)'\0'): return (unsigned)'\0';
    
    }
    
    printf("T B: (%c) ____\n", c);
    terror("Unrecognized nucleotide in hit");
    return 0; 
}

void inplace_reverse_and_complement(unsigned char *d, uint64_t l){
    uint64_t i;
    char c;
    for(i=0;i<l/2;i++){
        c = d[i];
        d[i] = d[l-i-1];
        d[l-i-1] = c;
        //printf("%c%c", d[i], d[l-i-1]);
    }
    
    //printf("\n");
    for(i=0;i<l;i++){
        //printf("eeee: %d", (int)d[i]);
        d[i] = complement(d[i]);
    }
}

