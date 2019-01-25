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



    /*
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
    */
    Tuple_hits * thit;
    
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

    /*
    fseek(database, 0, SEEK_END);
    uint64_t aprox_len_query = ftell(database);
    uint64_t aprox_len_db = aprox_len_query;
    rewind(database);
    */
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

    /*
    Container * ct = (Container *) calloc(1, sizeof(Container));
    if(ct == NULL) terror("Could not allocate container");    
    */



    Index * ctidx = (Index *) calloc(1, sizeof(Index));
    if(ctidx == NULL) terror("Could not allocate container");
    

    //begin = clock();
    

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
                        //printf("%"PRIu64"%%..wasted: (%e) (%e)", 1+100*pos_in_query/aprox_len_query, (double)(wasted_cycles_forward)/CLOCKS_PER_SEC, (double)(wasted_cycles_reverse)/CLOCKS_PER_SEC); 
                        fflush(stdout);
                    }
                    


                }else{ //It can be anything (including N, Y, X ...)

                    if(c != '\n' && c != '>'){
                        word_size = 0;
                        // data_database.sequences[pos_in_database++] = (unsigned char) 'N'; //Convert to N
                        ++current_len;

                    } 
                }
                //if(current_len % 1000000 == 0) printf(" curr len %" PRIu64"\n", current_len);
                if(word_size == custom_kmer){
                    //write to hash table
                    

                    thit = &ctidx->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                    [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                    [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                    [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];
                    
                    /*
                    typedef struct tuple_hits{
                        int repetition;
                        int hit_count;
                        uint64_t key;
                        uint64_t pos;
                    } Tuple_hits;
                    */

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

                    //thit->root = insert_AVLTree(thit->root, hashOfWord(&curr_kmer[0], custom_kmer, FIXED_K), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used);
                    //thit->root = insert_AVLTree(thit->root, collisioned_hash(&curr_kmer[0], custom_kmer), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used);
                    
                    

                    // Non overlapping
                    word_size = 0;
                    
		
		    // Overlapping
                    //memmove(&curr_kmer[0], &curr_kmer[1], custom_kmer-1);
                    //--word_size;
                }
            }
            word_size = 0;
            
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);    
        }
        
    }


    //end = clock();

    // data_database.total_len = pos_in_database;

    //fprintf(stdout, "[INFO] Database loaded and of length %"PRIu64". Hash table building took %e seconds\n", data_database.total_len, (double)(end-begin)/CLOCKS_PER_SEC);
    fprintf(stdout, "[INFO] Database loaded and of length %"PRIu64".\n", current_len);
    //close database
    fclose(database);

    
    
    //begin = clock();
    

    
    
    
    double pixel_size_db = (double) dimension / (double) current_len;
    double ratio_db = (double) current_len / dimension;


    // Get file length
    
    //fseek(query, 0, SEEK_END);
    //aprox_len_query = ftell(query);
    //rewind(query);
    aprox_len_query = get_seq_len(query);

    //uint64_t reallocs_hash_holder_table = 1;
    //uint64_t n_items_hash_holder_table = aprox_len_query / 5;

    //Hash_holder * hash_holder_table = (Hash_holder *) calloc(n_items_hash_holder_table, sizeof(Hash_holder));
    //if(hash_holder_table == NULL) terror("Could not allocate hash holding table");

    a_hundreth = (aprox_len_query/100);
    double pixel_size_query = (double) dimension / (double) aprox_len_query;
    double ratio_query = (double) aprox_len_query / dimension;
    

    double i_r_fix = MAX(1.0, custom_kmer * pixel_size_query);
    double j_r_fix = MAX(1.0, custom_kmer * pixel_size_db);

    
    
    fprintf(stdout, "[INFO] Ratios: Q [%e] D [%e]. Lenghts: Q [%"PRIu64"] D [%"PRIu64"]\n", ratio_query, ratio_db, aprox_len_query, current_len);
    fprintf(stdout, "[INFO] Pixel size: Q [%e] D [%e].\n", pixel_size_query, pixel_size_db);


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
                        
			
                        thit = &ctidx->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                        [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                        [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                        [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];
			

                        //AVLTree * search = find_AVLTree(thit->root, hash_forward);

                        if(thit->repetition == FALSE){

                        	hash_forward = collisioned_hash(&curr_kmer[0], custom_kmer);

				            if( hash_forward == thit->key){
                                // Attention ::::: you were not removing the ones with count==1 earlier 
                                thit->pos_in_y = current_len;
                                thit->hit_count++;
				            }
                        }

                        thit = &ctidx->table[char_converter[reverse_kmer[0]]][char_converter[reverse_kmer[1]]][char_converter[reverse_kmer[2]]]
                        [char_converter[reverse_kmer[3]]][char_converter[reverse_kmer[4]]][char_converter[reverse_kmer[5]]]
                        [char_converter[reverse_kmer[6]]][char_converter[reverse_kmer[7]]][char_converter[reverse_kmer[8]]]
                        [char_converter[reverse_kmer[9]]][char_converter[reverse_kmer[10]]][char_converter[reverse_kmer[11]]];

                        if(thit->repetition == FALSE){ 


                        	hash_reverse = collisioned_hash(&reverse_kmer[0], custom_kmer);

                            if(hash_reverse == thit->key){
        	                    // Attention ::::: you were not removing the ones with count==1 earlier 
	                            thit->pos_in_y = current_len;
	                            thit->hit_count++;
				            }
                        }

                        /*
                        if(search != NULL && search->count == 1){ //If count is two, then it is a rep
                            thit->hit_count += search->count;
                            
                            hash_holder_table[c_hash_holder].pos = current_len;
                            hash_holder_table[c_hash_holder].node = search;
                            hash_holder_table[c_hash_holder].th = thit;
                            ++c_hash_holder;
                            if(c_hash_holder == n_items_hash_holder_table*reallocs_hash_holder_table){
                                ++reallocs_hash_holder_table;
                                hash_holder_table = (Hash_holder *) realloc(hash_holder_table, n_items_hash_holder_table*reallocs_hash_holder_table*sizeof(Hash_holder));
                                if(hash_holder_table == NULL) terror("Could not realloc hash holder table");
                            }
                        }
                        */

                        

                        

                        //search = find_AVLTree(thit->root, hash_reverse);
                        /*
                        if(search != NULL && search->count == 1){ //If count is two, then it is a rep
                            
                            thit->hit_count += search->count;
                            hash_holder_table[c_hash_holder].pos = current_len;
                            hash_holder_table[c_hash_holder].node = search;
                            hash_holder_table[c_hash_holder].th = thit;
                            ++c_hash_holder;
                            if(c_hash_holder == n_items_hash_holder_table*reallocs_hash_holder_table){
                                ++reallocs_hash_holder_table;
                                hash_holder_table = (Hash_holder *) realloc(hash_holder_table, n_items_hash_holder_table*reallocs_hash_holder_table*sizeof(Hash_holder));
                                if(hash_holder_table == NULL) terror("Could not realloc hash holder table");
                            }
                        }
                        */

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


    /// Out

    fprintf(stdout, "Scanning hits table.\n");

    a_hundreth = MAX(1, TOTAL_ENTRIES/100);
    uint64_t t_computed = 0;
    uint64_t w0,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11;
    for(w0=0;w0<4;w0++){
        for(w1=0;w1<4;w1++){
            for(w2=0;w2<4;w2++){
                for(w3=0;w3<4;w3++){
                    for(w4=0;w4<4;w4++){
                        for(w5=0;w5<4;w5++){
                            for(w6=0;w6<4;w6++){
                                for(w7=0;w7<4;w7++){
                                    for(w8=0;w8<4;w8++){
                                        for(w9=0;w9<4;w9++){
                                            for(w10=0;w10<4;w10++){
                                                for(w11=0;w11<4;w11++){

                                                    if(t_computed % a_hundreth == 0){ 
                                                        fprintf(stdout, "\r%"PRIu64"%%...", 1+100*t_computed/TOTAL_ENTRIES); 
                                                        fflush(stdout);
                                                    }
                                                    ++t_computed;
                                                    Tuple_hits * taux = &ctidx->table[w0][w1][w2][w3][w4][w5][w6][w7][w8][w9][w10][w11];
                                                    if(taux->hit_count == 1){
                                                        // We plot it   
                                                        // Convert scale to representation
                                                        uint64_t redir_db = (uint64_t) (taux->pos / (ratio_db));
                                                        uint64_t redir_query = (uint64_t) (taux->pos_in_y / (ratio_query));
                                                        double i_r = i_r_fix; double j_r = j_r_fix;
                                                        while((uint64_t) i_r >= 1 && (uint64_t) j_r >= 1){
                                                            if((int64_t) redir_query - (int64_t) i_r > 0 && (int64_t) redir_db - (int64_t) j_r > 0){
                                                                representation[(int64_t) redir_query - (int64_t) i_r][(int64_t) redir_db - (int64_t) j_r]++;
                                                            }else{
                                                                representation[redir_query][redir_db]++;
                                                                break;
                                                            }
                                                            i_r -= MIN(1.0, pixel_size_query);
                                                            j_r -= MIN(1.0, pixel_size_db);
                                                        }                                                     
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    

    //double average_hit = ((double) total_hits / (double) table_size);
    //average_hit = 2.2;

    /*
    //fprintf(stdout, "[INFO] Total hit count is %"PRIu64" on a size of %"PRIu64" Avg = %e.\n", total_hits, table_size, average_hit);
    fprintf(stdout, "[INFO] Total hit count is %"PRIu64" on a size of %"PRIu64" E' = %e.\n", total_hits, table_size, Eprime);

    

    

    a_hundreth = MAX(1,c_hash_holder/100);
    
    for(current_len = 0; current_len < c_hash_holder; current_len++){

        if(current_len % a_hundreth == 0){ 
            fprintf(stdout, "\r%"PRIu64"%%...", 1+100*current_len/c_hash_holder); 
            fflush(stdout);
        }

        aux = hash_holder_table[current_len].node->next;

        //if(hash_holder_table[current_len].th->hit_count < (uint64_t) average_hit){
        if(hash_holder_table[current_len].th->hit_count < (uint64_t) Eprime){
            while(aux != NULL){
                // Convert scale to representation
                uint64_t redir_db = (uint64_t) (aux->pos / (ratio_db));
                uint64_t redir_query = (uint64_t) (hash_holder_table[current_len].pos / (ratio_query));
                double i_r = i_r_fix; double j_r = j_r_fix;
                while((uint64_t) i_r >= 1 && (uint64_t) j_r >= 1){
                    if((int64_t) redir_query - (int64_t) i_r > 0 && (int64_t) redir_db - (int64_t) j_r > 0){
                        representation[(int64_t) redir_query - (int64_t) i_r][(int64_t) redir_db - (int64_t) j_r]++;
                    }else{
                        representation[redir_query][redir_db]++;
                        break;
                    }
                    i_r -= MIN(1.0, pixel_size_query);
                    j_r -= MIN(1.0, pixel_size_db);
                }
                aux = aux->next;
            }
        }
    }
    */


    //end = clock();

    
    

    //fprintf(stdout, "\n[INFO] Query length %"PRIu64". Hits completed. Took %e seconds\n", data_query.total_len, (double)(end-begin)/CLOCKS_PER_SEC);
    fprintf(stdout, "\n[INFO] Query length %"PRIu64".\n", current_len);

    //begin = clock();

    //reads_per_thread = (uint64_t) (floorl((long double) data_query.n_seqs / (long double) n_threads));
    
    fprintf(stdout, "[INFO] Writing matrix.\n");


    uint64_t unique_diffuse = 0;
	fprintf(out_database, "%"PRIu64"\n", aprox_len_query);
    fprintf(out_database, "%"PRIu64"\n", aprox_len_db);
     // And replace 2's with 1's 
	
    for(i=0; i<dimension+1; i++){
        for(j=0; j<dimension; j++){
            fprintf(out_database, "%"PRIu64" ", representation[i][j]);
	    unique_diffuse += representation[i][j];
        }
        fprintf(out_database, "%"PRIu64"\n",  representation[i][dimension]);
	unique_diffuse += representation[i][dimension];
    }

    fprintf(stdout, "[INFO] Found %"PRIu64" unique hits for z = %"PRIu64".\n", unique_diffuse, diffuse_z);
    

    
    //free(ct->table);
    //free(hash_holder_table);
    /*
    for(i=0;i<=n_pools_used_AVL;i++){
        free(mp_AVL[i].base);
    }

    for(i=0;i<=n_pools_used;i++){
        free(mp[i].base);
    }
    */
    for(i=0;i<dimension;i++){
        free(representation[i]);
    }
    free(representation);
    if(out_database != NULL) fclose(out_database);
    
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
