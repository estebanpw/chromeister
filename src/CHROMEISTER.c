/*********

File        CHROMEISTER.c
Author      EPW <estebanpw@uma.es>
Description Computes hits and generates a dotplot

USAGE       Usage is described by calling ./CHROMEISTER --help



**********/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))
#define STARTING_SEQS 1000
#define PIECE_OF_DB_REALLOC 3200000 //half a gigabyte if divided by 8 bytes

uint64_t custom_kmer = 12; // Defined as external in structs.h

void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, uint64_t * custom_kmer, uint64_t * dimension);

int VERBOSE_ACTIVE = 0;

int main(int argc, char ** av){
    

    clock_t begin, end;

    uint64_t i, j;

    //query to read kmers from, database to find seeds
    FILE * query = NULL, * database = NULL, * out_database = NULL;
    uint64_t dimension = 1000; // Default 1000 * 1000
    
    
    init_args(argc, av, &query, &database, &out_database, &custom_kmer, &dimension);

    unsigned char char_converter[91];
    char_converter[(unsigned char)'A'] = 0;
    char_converter[(unsigned char)'C'] = 1;
    char_converter[(unsigned char)'G'] = 2;
    char_converter[(unsigned char)'T'] = 3;

    
    //Variables to account for positions
    //Print info
    fprintf(stdout, "[INFO] Loading database\n");
    //Variables to read kmers
    char c = 'N'; //Char to read character
    //Current length of array and variables for the buffer
    uint64_t idx = 0, r = 0;
    
    //Vector to read in batches
    char * temp_seq_buffer = NULL;
    if ((temp_seq_buffer = calloc(READBUF, sizeof(char))) == NULL) {
        terror("Could not allocate memory for read buffer");
    }

    //Dimensional matrix
    unsigned char ** representation = (unsigned char **) calloc(dimension+1, sizeof(unsigned char *));
    if(representation == NULL) terror("Could not allocate representation");
    for(i=0; i<dimension+1; i++){
        representation[i] = (unsigned char *) calloc(dimension+1, sizeof(unsigned char));
        if(representation[i] == NULL) terror("Could not allocate second loop representation");
    }

    unsigned char curr_kmer[custom_kmer], reverse_kmer[custom_kmer];
    curr_kmer[0] = reverse_kmer[0] = '\0';
    uint64_t word_size = 0, word_size_rev = 0, pos_in_database = 0;

    //To hold all information related to database
    SeqInfo data_database;
    SeqInfo data_query;
    data_database.total_len = 0;
    data_database.n_seqs = 0;
    
    //To force reading from the buffer
    idx = READBUF + 1;

    //Store positions of kmers
    uint64_t n_pools_used = 0;
    //Mempool_l * mp = (Mempool_l *) malloc(MAX_MEM_POOLS*sizeof(Mempool_l));
    //if(mp == NULL) terror("Could not allocate vectors for memory pools");
    Mempool_l mp[MAX_MEM_POOLS];
    init_mem_pool_llpos(&mp[n_pools_used]);
    llpos * aux, * pointer;

    unsigned char aux_kmer[custom_kmer+1];
    
    //Vector to store query seq
    unsigned char * seq_vector_query = (unsigned char *) malloc(READBUF*sizeof(unsigned char));
    if(seq_vector_query == NULL) terror("Could not allocate memory for query vector");
    uint64_t pos_in_query = 0;


    Container * ct = (Container *) calloc(1, sizeof(Container));
    if(ct == NULL) terror("Could not allocate container");
    
    

    begin = clock();

    c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
    while((!feof(database) || (feof(database) && idx < r))){

        if(c == '>'){
            data_database.n_seqs++;
            // data_database.start_pos[data_database.n_seqs++] = pos_in_database;
            
            // if(pos_in_database == READBUF*n_realloc_database){ 
            //     n_realloc_database++; data_database.sequences = (unsigned char *) realloc(data_database.sequences, READBUF*n_realloc_database*sizeof(unsigned char));
            //     if(data_database.sequences == NULL) terror("Could not reallocate temporary database");
            // }

            // if(data_database.n_seqs == INITSEQS*n_seqs_database_realloc){
            //     n_seqs_database_realloc++; data_database.start_pos =  (uint64_t *) realloc(data_database.start_pos, INITSEQS*n_seqs_database_realloc*sizeof(uint64_t));
            // }


            while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);  //Skip ID
                

            while(c != '>' && (!feof(database) || (feof(database) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    curr_kmer[word_size] = (unsigned char) c;
                    if(word_size < custom_kmer) word_size++;
                    pos_in_database++;
                    //data_database.sequences[pos_in_database++] = (unsigned char) c;
            
                    // if(pos_in_database == READBUF*n_realloc_database){ 
                    //     n_realloc_database++; data_database.sequences = (unsigned char *) realloc(data_database.sequences, READBUF*n_realloc_database*sizeof(unsigned char));
                    //     if(data_database.sequences == NULL) terror("Could not reallocate temporary database");
                    // }


                }else{ //It can be anything (including N, Y, X ...)

                    if(c != '\n' && c != '\r' && c != '>'){
                        word_size = 0;
                        // data_database.sequences[pos_in_database++] = (unsigned char) 'N'; //Convert to N
                        pos_in_database++;
                        // if(pos_in_database == READBUF*n_realloc_database){ 
                        //     n_realloc_database++; data_database.sequences = (unsigned char *) realloc(data_database.sequences, READBUF*n_realloc_database*sizeof(unsigned char));
                        // if(data_database.sequences == NULL) terror("Could not reallocate temporary database");
                        // }
                    } 
                }
                if(word_size == custom_kmer){
                    //write to hash table
                    
		
                    pointer = ct->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                    [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                    [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                    [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];

                    

                    if(pointer == NULL){

                        pointer = getNewLocationllpos(mp, &n_pools_used);
                        

                        pointer->pos = pos_in_database;

                        pointer->extended_hash = hashOfWord(&curr_kmer[FIXED_K], custom_kmer - FIXED_K);

                        pointer->next = NULL;

                    

                    }else{

                        
                        aux = ct->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                        [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                        [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                        [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];

                        pointer = getNewLocationllpos(mp, &n_pools_used);

                        pointer->pos = pos_in_database;
                        pointer->extended_hash = hashOfWord(&curr_kmer[FIXED_K], custom_kmer - FIXED_K);
                        pointer->next = aux;

                    }

                    ct->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                    [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                    [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                    [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]] = pointer;
		
		            memcpy(aux_kmer, &curr_kmer[1], custom_kmer-1);
                    memcpy(curr_kmer, aux_kmer, custom_kmer-1);
                    word_size--;
                }
            }
            word_size = 0;
            
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);    
        }
        
    }

    end = clock();

    data_database.total_len = pos_in_database;

    fprintf(stdout, "[INFO] Database loaded and of length %"PRIu64". Hash table building took %e seconds\n", data_database.total_len, (double)(end-begin)/CLOCKS_PER_SEC);
    //close database
    fclose(database);

    
    
    begin = clock();

    

    //To force reading from the buffer
    idx = READBUF + 1;

    
    //data_query.sequences = seq_vector_query;
    data_query.total_len = 0;
    data_query.n_seqs = 0;
    
    
    //Print info
    fprintf(stdout, "[INFO] Generating hits\n");   

    // Get file length
    fseek(query, 0, SEEK_END);
    uint64_t aprox_len_query = ftell(query);
    rewind(query);
    long double pixel_size_query = (long double) dimension / (long double) aprox_len_query;
    long double pixel_size_db = (long double) dimension / (long double) data_database.total_len;
    
    long double ratio_query = (long double) aprox_len_query / dimension;
    long double ratio_db = (long double) data_database.total_len / dimension;

    uint64_t hash_forward, hash_reverse;
    
    fprintf(stdout, "[INFO] Ratios: Q [%Le] D [%Le]. Lenghts: Q [%"PRIu64"] D [%"PRIu64"]\n", ratio_query, ratio_db, aprox_len_query, data_database.total_len);

    c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);    
    
    while((!feof(query) || (feof(query) && idx < r))){

        if(c == '>'){
            data_query.n_seqs++;
            word_size = 0;
            word_size_rev = custom_kmer-1;
            



            while(c != '\n'){ c = buffered_fgetc(temp_seq_buffer, &idx, &r, query); } //Skip ID
                

            while(c != '>' && (!feof(query) || (feof(query) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    
                    pos_in_query++;
                    if(pos_in_query % (aprox_len_query/10) == 0) printf("%"PRIu64"%%..", 1+100*pos_in_query/aprox_len_query);
                    curr_kmer[word_size++] = (unsigned char) c;

                    switch(c){
                        case ('A'): reverse_kmer[word_size_rev] = (unsigned)'T';
                        break;
                        case ('C'): reverse_kmer[word_size_rev] = (unsigned)'G';
                        break;
                        case ('G'): reverse_kmer[word_size_rev] = (unsigned)'C';
                        break;
                        case ('T'): reverse_kmer[word_size_rev] = (unsigned)'A';
                        break;
                    }
                    if(word_size_rev != 0) word_size_rev--;

                    
                    if(word_size == custom_kmer){

                        aux = ct->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                        [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                        [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                        [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];


                        hash_forward = hashOfWord(&curr_kmer[FIXED_K], custom_kmer - FIXED_K);
                        hash_reverse = hashOfWord(&reverse_kmer[FIXED_K], custom_kmer - FIXED_K);

                        while(aux != NULL){

                            if(aux->extended_hash == hash_forward){
                                

                                // Convert scale to representation
                                uint64_t redir_db = (uint64_t) (aux->pos / (ratio_db));
                                uint64_t redir_query = (uint64_t) (pos_in_query / (ratio_query));

                                long double i_r = MAX(1.0, custom_kmer * pixel_size_query);
                                long double j_r = MAX(1.0, custom_kmer * pixel_size_db);
                                while((uint64_t) i_r >= 1 && (uint64_t) j_r >= 1){
                                    // printf("I have %Le %Le which is %"PRIu64" %"PRIu64"\n", i_r, j_r, (uint64_t) i_r, (uint64_t) j_r);
                                    // printf("Hit is at %"PRIu64", %"PRIu64"\n", pos_in_query, aux->pos);
                                    // printf("redir : %"PRIu64" %"PRIu64"\n", redir_query, redir_db);
                                    // getchar();
                                     if((int64_t) redir_query - (int64_t) i_r > 0 && (int64_t) redir_db - (int64_t) j_r > 0){
                                         representation[(int64_t) redir_query - (int64_t) i_r][(int64_t) redir_db - (int64_t) j_r] = 1;
                                    }else{
                                        representation[redir_query][redir_db] = 1;
                                        break;
                                    }
                                    if(pixel_size_query == 0 || pixel_size_db == 0) break;
                                    
                                    i_r -= pixel_size_query;
                                    j_r -= pixel_size_db;

                                }
                            }

                            
                            
                            if(aux->extended_hash == hash_reverse){
                                //printf("enter\n");
                    

                                // Convert scale to representation
                                uint64_t redir_db = (uint64_t) ( (aux->pos + custom_kmer) / (ratio_db));
                                uint64_t redir_query = (uint64_t) ((pos_in_query ) / (ratio_query));

                                long double i_r = MAX(1.0, custom_kmer * pixel_size_query);
                                long double j_r = MAX(1.0, custom_kmer * pixel_size_db);
                                while((uint64_t) i_r >= 1 && (uint64_t) j_r >= 1){
                                    // printf("I have %Le %Le which is %"PRIu64" %"PRIu64"\n", i_r, j_r, (uint64_t) i_r, (uint64_t) j_r);
                                    // printf("Hit is at %"PRIu64", %"PRIu64"\n", pos_in_query, aux->pos);
                                    // printf("redir : %"PRIu64" %"PRIu64"\n", redir_query, redir_db);
                                    // printf("drwaing: %"PRId64", %"PRId64"\n", (int64_t) redir_query - (int64_t) i_r, (int64_t) redir_db - (int64_t) j_r);
                                    // getchar();
                                    if((int64_t) redir_query + (int64_t) i_r < dimension && (int64_t) redir_db - (int64_t) j_r > 0){
                                         representation[(int64_t) redir_query + (int64_t) i_r][(int64_t) redir_db - (int64_t) j_r] = 1;
                                    }else{
                                        representation[redir_query][redir_db] = 1;
                                        break;
                                    }
                                    if(pixel_size_query == 0 || pixel_size_db == 0) break;
                                    i_r -= pixel_size_query;
                                    j_r -= pixel_size_db;

                                }
                            }
                            
                            aux = aux->next;
                        }

                        memmove(&curr_kmer[0], &curr_kmer[1], custom_kmer-1);
                        memmove(&reverse_kmer[1], &reverse_kmer[0], custom_kmer-1);
                        word_size--;
                    }


            
                    // if(pos_in_query == READBUF*n_realloc_database){ 
                    //     n_realloc_database++; data_query.sequences = (unsigned char *) realloc(data_query.sequences, READBUF*n_realloc_database*sizeof(unsigned char));
                    //     if(data_query.sequences == NULL) terror("Could not reallocate temporary query");
                    // }
                }else{
                    if(c != '\n' && c != '\r' && c != '>'){
                        word_size = 0;
                        word_size_rev = custom_kmer-1;
                        pos_in_query++;
                        // data_query.sequences[pos_in_query++] = (unsigned char) 'N'; //Convert to N
                        // if(pos_in_query == READBUF*n_realloc_database){ 
                        //     n_realloc_database++; data_query.sequences = (unsigned char *) realloc(data_query.sequences, READBUF*n_realloc_database*sizeof(unsigned char));
                        //     if(data_query.sequences == NULL) terror("Could not reallocate temporary query");
                        // }
                    }
                }
            }
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);    
        }
        
    }

    end = clock();

    
    
    data_query.total_len = pos_in_query;

    fprintf(stdout, "\n[INFO] Query length %"PRIu64". Hits completed. Took %e seconds\n", data_query.total_len, (double)(end-begin)/CLOCKS_PER_SEC);

    begin = clock();


    /*
    Queue * traverse = queue_head.head;
    while(traverse != NULL){
        printf("current_piece: %"PRIu64"-%"PRIu64"\n", traverse->r1, traverse->r2);
        traverse = traverse->next;
    }
    */


    //reads_per_thread = (uint64_t) (floorl((long double) data_query.n_seqs / (long double) n_threads));
    
    fprintf(stdout, "[INFO] Writing matrix\n");

    // representation[0][0] = 1;
    // representation[1][0] = 1;
    // representation[2][0] = 1;
    // representation[3][0] = 1;
    // representation[4][0] = 1;
    // representation[5][0] = 1;


    for(i=0; i<dimension+1; i++){
        for(j=0; j<dimension+1; j++){
            fprintf(out_database, "%d ", (int) representation[i][j]);
        }
        fprintf(out_database, "\n");
    }

    /*
    uint64_t z;
    for(z=0; z<POOL_SIZE; z++){
        aux = mp[0].base + z;
        fprintf(stdout, "%p\n", aux);
        fflush(stdout);
    }
    */
    
    free(ct->table);

    for(i=0;i<=n_pools_used;i++){
        free(mp[i].base);
    }
    for(i=0;i<dimension;i++){
        free(representation[i]);
    }
    free(representation);
    if(out_database != NULL) fclose(out_database);
    
    return 0;
}

void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, uint64_t * custom_kmer, uint64_t * dimension){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--verbose") == 0) VERBOSE_ACTIVE = 1;
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           CHROMEISTER -query [query] -db [database] -out [outfile]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -kmer       [Integer:   k>1 (default 12)]\n");
            fprintf(stdout, "           -dimension  Size of the output [Integer:   d>0 (default 1000)]\n");
            fprintf(stdout, "           -out        [File path]\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            fprintf(stdout, "PLEASE NOTICE: The reverse complementary is calculated for the QUERY.");
            exit(1);
        }
        if(strcmp(av[pNum], "-query") == 0){
            *query = fopen64(av[pNum+1], "rt");
            if(query==NULL) terror("Could not open query file");
        }
        if(strcmp(av[pNum], "-db") == 0){
            *database = fopen64(av[pNum+1], "rt");
            if(database==NULL) terror("Could not open database file");
        }
        if(strcmp(av[pNum], "-out") == 0){
            *out_database = fopen64(av[pNum+1], "wt");
            if(out_database==NULL) terror("Could not open output database file");
        }
        if(strcmp(av[pNum], "-kmer") == 0){
            *custom_kmer = (uint64_t) atoi(av[pNum+1]);
            if(*custom_kmer < 2) terror("K-mer size must be larger than 1");
        }
        if(strcmp(av[pNum], "-dimension") == 0){
            *dimension = (uint64_t) atoi(av[pNum+1]);
            if(*dimension < 1) terror("Dimension must be a positive integer");
        }
        pNum++;
    }
    
    if(*query==NULL || *database==NULL || *out_database==NULL) terror("A query, database and output file is required");
}

