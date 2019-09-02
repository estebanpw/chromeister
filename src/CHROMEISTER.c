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
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"
#define STARTING_SEQS 1000
#define PIECE_OF_DB_REALLOC 3200000 //half a gigabyte if divided by 8 bytes
#define RANGE 2

uint64_t custom_kmer = 32; // Defined as external in structs.h
uint64_t diffuse_z = 4; // Defined as external in structs.h

uint64_t get_seq_len(FILE * f);


void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, uint64_t * custom_kmer, uint64_t * dimension, uint64_t * diffuse_z);

int main(int argc, char ** av){



    
    //Store positions of kmers
    uint64_t n_pools_used = 0;
    //Mempool_l * mp = (Mempool_l *) malloc(MAX_MEM_POOLS*sizeof(Mempool_l));
    //if(mp == NULL) terror("Could not allocate vectors for memory pools");
    Mempool_l mp[MAX_MEM_POOLS];
    init_mem_pool_llpos(&mp[n_pools_used]);
    //llpos * aux;

    uint64_t n_pools_used_AVL = 0;
    Mempool_AVL mp_AVL[MAX_MEM_POOLS];
    init_mem_pool_AVL(&mp_AVL[n_pools_used_AVL]);
    
    //Tuple_hits * thit;
    
    AVLTree * root = NULL;

    /*
    AVLTree * root = NULL;
    root = insert_AVLTree(root, 10, mp_AVL, &n_pools_used_AVL, 0, mp, &n_pools_used);
    
    llpos * some = find_AVLTree(root, 25);
    while(some != NULL){
        printf("#%"PRIu64", ", some->pos); some = some->next;
    }
    */

    uint64_t i, j;

    //query to read kmers from, database to find seeds
    FILE * query = NULL, * database = NULL, * out_database = NULL;
    uint64_t dimension = 1000; // Default 1000 * 1000
    
    
    init_args(argc, av, &query, &database, &out_database, &custom_kmer, &dimension, &diffuse_z);



    unsigned char char_converter[91];
    char_converter[(unsigned char)'A'] = 0;
    char_converter[(unsigned char)'C'] = 1;
    char_converter[(unsigned char)'G'] = 2;
    char_converter[(unsigned char)'T'] = 3;

    
    // Variables to account for positions
    // Print info
    fprintf(stdout, "[INFO] Loading database\n");
    // Variables to read kmers
    char c = 'N'; //Char to read character
    // Current length of array and variables for the buffer
    uint64_t idx = 0, r = 0;
    
    //Vector to read in batches
    char * temp_seq_buffer = NULL;
    if ((temp_seq_buffer = calloc(READBUF, sizeof(char))) == NULL) {
        terror("Could not allocate memory for read buffer");
    }

    //Dimensional matrix
    uint64_t ** representation = (uint64_t **) calloc(dimension+1, sizeof(uint64_t *));
    if(representation == NULL) terror("Could not allocate representation");
    for(i=0; i<dimension+1; i++){
        representation[i] = (uint64_t *) calloc(dimension+1, sizeof(uint64_t));
        if(representation[i] == NULL) terror("Could not allocate second loop representation");
    }

    
    uint64_t aprox_len_query = get_seq_len(database);
    uint64_t aprox_len_db = aprox_len_query;

    uint64_t a_hundreth = (aprox_len_query/100);

    unsigned char curr_kmer[custom_kmer], reverse_kmer[custom_kmer];
    curr_kmer[0] = reverse_kmer[0] = '\0';
    uint64_t word_size = 0, word_size_rev = 0;

    //To hold all information related to database
    uint64_t current_len = 0;
    
    //To force reading from the buffer
    idx = READBUF + 1;

    //unsigned char aux_kmer[custom_kmer+1];
    
    //Vector to store query seq
    unsigned char * seq_vector_query = (unsigned char *) malloc(READBUF*sizeof(unsigned char));
    if(seq_vector_query == NULL) terror("Could not allocate memory for query vector");

    


    Index * ctidx = (Index *) calloc(1, sizeof(Index));
    if(ctidx == NULL) terror("Could not allocate container");
    

    c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
    while((!feof(database) || (feof(database) && idx < r))){

        if(c == '>'){
            


            while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);  //Skip ID
                

            while(c != '>' && (!feof(database) || (feof(database) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    curr_kmer[word_size] = (unsigned char) c;
                    if(word_size < custom_kmer) ++word_size;
                    ++current_len;
                    if(current_len % a_hundreth == 0){ 
                        fprintf(stdout, "\r%"PRIu64"%%...", 1+100*current_len/aprox_len_query); 
                        fflush(stdout);
                    }
                    


                }else{ //It can be anything (including N, Y, X ...)

                    if(c != '\n' && c != '>'){
                        word_size = 0;
                        ++current_len;

                    } 
                }
                if(word_size == custom_kmer){
                    //write to hash table
                    
                    /*
                    thit = &ctidx->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                    [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                    [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                    [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];
                    
                    

                    if(thit->repetition == UNSET){
                        // Then we can insert
                        thit->hit_count = 0;
                        thit->key = collisioned_hash(&curr_kmer[0], custom_kmer);
                        thit->pos = current_len;
                        thit->repetition = FALSE;
                    }else{
                        // Otherwise we break it
                        thit->repetition = TRUE;
                    }

                    */

                    root = insert_AVLTree_x(root, hashOfWord(&curr_kmer[0], custom_kmer, 0), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used);
                    //root = insert_AVLTree(root, collisioned_hash(&curr_kmer[0], custom_kmer), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used);
                    
                    

                    // Non overlapping
                    //word_size = 0;
                    
		
		            // Overlapping
                    memmove(&curr_kmer[0], &curr_kmer[1], custom_kmer-1);
                    --word_size;
                }
            }
            word_size = 0;
            
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);    
        }
        
    }

    fprintf(stdout, "[INFO] Database loaded and of length %"PRIu64".\n", current_len);
    //close database
    fclose(database);

    
    
    //begin = clock();
    

    
    
    
    // Get file length
    
    aprox_len_query = get_seq_len(query);

    //uint64_t reallocs_hash_holder_table = 1;
    //uint64_t n_items_hash_holder_table = aprox_len_query / 5;

    //Hash_holder * hash_holder_table = (Hash_holder *) calloc(n_items_hash_holder_table, sizeof(Hash_holder));
    //if(hash_holder_table == NULL) terror("Could not allocate hash holding table");

    a_hundreth = (aprox_len_query/100);


    fprintf(stdout, "[INFO] Computing absolute hit numbers.\n");


    current_len = 0;

    //llpos * the_original_hit;

    //To force reading from the buffer
    idx = READBUF + 1;
    c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);

    //uint64_t c_hash_holder = 0;
    
    while((!feof(query) || (feof(query) && idx < r))){

        if(c == '>'){
            word_size = 0;
            word_size_rev = custom_kmer-1;
            



            while(c != '\n'){ c = buffered_fgetc(temp_seq_buffer, &idx, &r, query); } //Skip ID
                

            while(c != '>' && (!feof(query) || (feof(query) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    
                    ++current_len;
                    if(current_len % a_hundreth == 0){ 
                        fprintf(stdout, "\r%"PRIu64"%%...", 1+100*current_len/aprox_len_query); 
                        fflush(stdout);
                    }
                    curr_kmer[word_size] = (unsigned char) c;
                    ++word_size;

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
                    if(word_size_rev != 0) --word_size_rev;



                    
                    if(word_size == custom_kmer){


                        //hash_forward = hashOfWord(&curr_kmer[0], custom_kmer, FIXED_K);
                        //hash_reverse = hashOfWord(&reverse_kmer[0], custom_kmer, FIXED_K);
                        uint64_t hash_forward, hash_reverse;
                        //hash_forward = collisioned_hash(&curr_kmer[0], custom_kmer);
                        //hash_reverse = collisioned_hash(&reverse_kmer[0], custom_kmer);
                        
                        /*
                        thit = &ctidx->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                        [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                        [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                        [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];

                        */

                        root = insert_AVLTree_y(root, hashOfWord(&curr_kmer[0], custom_kmer, 0), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used);
                        root = insert_AVLTree_y(root, hashOfWord(&reverse_kmer[0], custom_kmer, 0), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used);

                        // Overlapping
                        
                        memmove(&curr_kmer[0], &curr_kmer[1], custom_kmer-1);
                        memmove(&reverse_kmer[1], &reverse_kmer[0], custom_kmer-1);
                        --word_size;

                        // Non overlapping
                        //word_size = 0;
                        //word_size_rev = custom_kmer-1;
                    }
                }else{
                    if(c != '\n' && c != '>'){
                        word_size = 0;
                        word_size_rev = custom_kmer-1;
                        ++current_len;
                    }
                }
            }
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);    
        }
        
    }






    return 0;
}

void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, uint64_t * custom_kmer, uint64_t * dimension, uint64_t * diffuse_z){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           CHROMEISTER -query [query] -db [database] -out [outfile]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -kmer       [Integer:   k>1 (default 32)]\n");
	    fprintf(stdout, "		-diffuse    [Integer:   z>0 (default 4)]\n");
            fprintf(stdout, "           -dimension  Size of the output [Integer:   d>0 (default 1000)]\n");
            fprintf(stdout, "           -out        [File path]\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            fprintf(stdout, "PLEASE NOTICE: The reverse complementary is calculated for the QUERY.\n");
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
            *out_database = fopen(av[pNum+1], "wt");
            if(out_database==NULL) terror("Could not open output database file");
        }
        if(strcmp(av[pNum], "-kmer") == 0){
            *custom_kmer = (uint64_t) atoi(av[pNum+1]);
            if(*custom_kmer < BYTES_IN_MER) terror("K-mer size must be larger than 4");
            if(*custom_kmer % BYTES_IN_MER != 0) terror("K-mer size must be a multiple of 4");

        }
        if(strcmp(av[pNum], "-diffuse") == 0){
            *diffuse_z = (uint64_t) atoi(av[pNum+1]);
            if(*diffuse_z == 0 || *diffuse_z > 32) terror("Z must satisfy 0<z<=32");

        }

        if(strcmp(av[pNum], "-dimension") == 0){
            *dimension = (uint64_t) atoi(av[pNum+1]);
            if(*dimension < 1) terror("Dimension must be a positive integer");
        }
        
        pNum++;
    }
    
    if(*query==NULL || *database==NULL || *out_database==NULL) terror("A query, database and output file is required");
}

uint64_t get_seq_len(FILE * f) {
    char c = '\0';
    uint64_t l = 0;

    while(!feof(f)){
        c = getc(f);

        if(c == '>'){
            while(c != '\n') c = getc(f);
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z'){
            ++l;
        }
    }


    rewind(f);
    return l;
}
