/*********

File        lowcompfilter.c
Author      EPW <estebanpw@uma.es>
Description blablalbal

USAGE       ./lowcompfilter --help



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
#define ENTROPY_WINDOW 12

uint64_t custom_kmer = 32; // Defined as external in structs.h
uint64_t diffuse_z = 4; // Defined as external in structs.h

typedef struct initial_D{
    unsigned char r;
    unsigned char p;
    unsigned char first;
} Initial_D;

void init_args(int argc, char ** av, FILE ** database, unsigned char * cutoff);

int main(int argc, char ** av){
    

    

    uint64_t i, j;

    //query to read kmers from, database to find seeds
    FILE * database = NULL, * out_database = NULL;

    unsigned char cutoff = 5;
    unsigned char current_cut = 0;
    
    
    init_args(argc, av, &database, &cutoff);
    out_database = fopen("seq-filtered.fasta", "wt");
    if(out_database == NULL) terror("Could not open output");



    unsigned char char_converter[91];
    char_converter[(unsigned char)'A'] = 0;
    char_converter[(unsigned char)'C'] = 1;
    char_converter[(unsigned char)'G'] = 2;
    char_converter[(unsigned char)'T'] = 3;

    double chances[91];
    chances[(unsigned char)'A'] = 0.29;
    chances[(unsigned char)'C'] = 0.21;
    chances[(unsigned char)'G'] = 0.21;
    chances[(unsigned char)'T'] = 0.29;

    double log_chances[91];
    log_chances[(unsigned char)'A'] = log(chances[(unsigned char)'A']) / log(2);
    log_chances[(unsigned char)'C'] = log(chances[(unsigned char)'C']) / log(2);
    log_chances[(unsigned char)'G'] = log(chances[(unsigned char)'G']) / log(2);
    log_chances[(unsigned char)'T'] = log(chances[(unsigned char)'T']) / log(2);

    double first_entropies[ENTROPY_WINDOW];
    double current_entropy = 0;
    int added_entropy;
    
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

    custom_kmer = 12;

    unsigned char curr_kmer[custom_kmer];
    curr_kmer[0] = '\0';
    uint64_t word_size = 0;

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

    //begin = clock();
    char * huge_seq = (char *) malloc(MAX_SEQ);
    uint64_t write_pos = 0;
    if(huge_seq == NULL) terror("Could not allocate huge sequence");
    // Tables

    uint64_t ts1 = (uint64_t) pow(4,1);
    uint64_t ts2 = (uint64_t) pow(4,2);
    uint64_t ts3 = (uint64_t) pow(4,3);
    uint64_t ts4 = (uint64_t) pow(4,4);
    uint64_t ts5 = (uint64_t) pow(4,5);
    uint64_t ts6 = (uint64_t) pow(4,6);

    Initial_D * t1 = (Initial_D *) calloc(ts1, sizeof(Initial_D));
    Initial_D * t2 = (Initial_D *) calloc(ts2, sizeof(Initial_D));
    Initial_D * t3 = (Initial_D *) calloc(ts3, sizeof(Initial_D));
    Initial_D * t4 = (Initial_D *) calloc(ts4, sizeof(Initial_D));
    Initial_D * t5 = (Initial_D *) calloc(ts5, sizeof(Initial_D));
    Initial_D * t6 = (Initial_D *) calloc(ts6, sizeof(Initial_D));

    // Add entropy filter
    // AAATAAAATAAAT

    c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
    huge_seq[write_pos++] = c;
    while((!feof(database) || (feof(database) && idx < r))){

        if(c == '>'){
            current_entropy = 0;
            added_entropy = 0;
            while(c != '\n'){ c = buffered_fgetc(temp_seq_buffer, &idx, &r, database); huge_seq[write_pos++] = c; }
            while(c != '>' && (!feof(database) || (feof(database) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
                c = toupper(c);
                huge_seq[write_pos++] = c;
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    curr_kmer[word_size] = (unsigned char) c;
                    first_entropies[added_entropy] = chances[(unsigned char) c] * log_chances[(unsigned char)c];
                    current_entropy += first_entropies[added_entropy];
                    ++added_entropy;
                    if(current_cut == 0) current_cut++;
                    else{
                        if(c == curr_kmer[word_size-1]) current_cut++;
                        else current_cut = 1;
                    }
                    if(current_cut == cutoff){
                        i=1; j=0;
                        while(i-j<=(uint64_t) cutoff){
                            
                            if(huge_seq[write_pos-(i)] != '\n') huge_seq[write_pos-(i)] = 'N';
                            if(huge_seq[write_pos-(i)] == '\n') j++;
                            i++;

                        }
                        current_cut--;
                        
                    }
                    if(word_size < custom_kmer) ++word_size;
                    ++current_len;
                // Wasted kmers 
                }else{ if(c != '\n' && c != '>'){ word_size = 0; ++current_len; current_cut = 0; current_entropy = 0; added_entropy = 0; } }
                

                if(added_entropy == ENTROPY_WINDOW){
                    current_entropy = 0;
                    for(i=0;i<ENTROPY_WINDOW;i++){
                        current_entropy += first_entropies[i];
                    }
                    printf("Entropy of %.12s -> %e\n", curr_kmer, -current_entropy);

                    current_entropy -= first_entropies[0];
                    memmove(&first_entropies[0], &first_entropies[1], ENTROPY_WINDOW-1);
                    --added_entropy;
                    getchar();
                }
                
                
                if(word_size == 12){
                    curr_kmer[12] = '\0';
                    memset(t1, 0x0, ts1*sizeof(Initial_D));
                    memset(t2, 0x0, ts2*sizeof(Initial_D));
                    memset(t3, 0x0, ts3*sizeof(Initial_D));
                    memset(t4, 0x0, ts4*sizeof(Initial_D));
                    memset(t5, 0x0, ts5*sizeof(Initial_D));
                    memset(t6, 0x0, ts6*sizeof(Initial_D));
                    
                    unsigned char n_rep;
                    uint64_t hash;
                    int debug = 0, death = 0;
                    if(current_len % 5000000 == 0) printf("Computed 5 millions\n");
                    //if(strcmp("ACACACACACAC", (char *) curr_kmer) == 0){ printf("FOUND %s\n", curr_kmer); debug = 1; }
                    for(i=1; i<7; i++){
                        if(death != 0) break;

                        for(j=0; j<12; j++){
                            
                            
                            if(j + i <= 10 && death == 0){
                                switch(i){
                                    case 1: {
                                        hash = hashOfWord(&curr_kmer[j], 1, 0);
                                        //printf("(T1)hash: %"PRIu64" ->attempting insert %.1s at %u prior pos is (%u) result pos is = %u; curr_rep = %u\n", hash, &curr_kmer[j], j, t5[hash].p, t5[hash].p + 4, t5[hash].r);
                                        if(t1[hash].first == 0) n_rep = ++(t1[hash].r); else if(t1[hash].p <= (unsigned char) j) ++(t1[hash].r); 
                                        t1[hash].first = 1;
                                        n_rep = (t1[hash].r);
                                        t1[hash].p = (unsigned char) (j+1);
                                        //if(n_rep == 10){ printf("%s\n", curr_kmer); printf("Got longest repeat of size 1"); getchar(); }
                                        if(n_rep == 10) death = 1;
                                    } break;
                                    case 2: { 
                                        hash = hashOfWord(&curr_kmer[j], 2, 0);
                                        //printf("(T2)hash: %"PRIu64" ->attempting insert %.2s at %u prior pos is (%u) result pos is = %u; curr_rep = %u\n", hash, &curr_kmer[j], j, t2[hash].p, t2[hash].p + 2, t2[hash].r); getchar();
                                        if(t2[hash].first == 0) { n_rep = ++(t2[hash].r); if(debug==1) printf("Done\n"); }else if(t2[hash].p <= (unsigned char) j){  ++(t2[hash].r); if(debug==1) printf("DONE\n") ;} 
                                        t2[hash].first = 1;
                                        n_rep = (t2[hash].r);
                                        t2[hash].p = (unsigned char) (j+2);
                                        //if(n_rep == 5){ printf("%s\n", curr_kmer); printf("Got longest repeat of size 2"); getchar(); }
                                        if(n_rep == 5) death = 1;
                                    } break;    
                                    case 5: { 
                                        
                                        hash = hashOfWord(&curr_kmer[j], 5, 0);
                                        //printf("(T5)hash: %"PRIu64" ->attempting insert %.5s at %u prior pos is (%u) result pos is = %u; curr_rep = %u\n", hash, &curr_kmer[j], j, t5[hash].p, t5[hash].p + 4, t5[hash].r);
                                        if(t5[hash].first == 0) n_rep = ++(t5[hash].r); else if(t5[hash].p <= (unsigned char) j) ++(t5[hash].r); 
                                        t5[hash].first = 1;
                                        n_rep = (t5[hash].r);
                                        t5[hash].p = (unsigned char) (j+5);
                                        //if(n_rep == 2){ printf("%s\n", curr_kmer); printf("Got longest repeat of size 5"); getchar(); }
                                        if(n_rep == 2) death = 1;
                                    } break;    
                                }
                            }

                            
                            if(j + i < 12 && death == 0){
                                switch(i){
                                    case 3: { 
                                        hash = hashOfWord(&curr_kmer[j], 3, 0);
                                        //printf("(T3)hash: %"PRIu64" ->attempting insert %.3s at %u prior pos is (%u) result pos is = %u; curr_rep = %u\n", hash, &curr_kmer[j], j, t5[hash].p, t5[hash].p + 4, t5[hash].r);
                                        if(t3[hash].first == 0) n_rep = ++(t3[hash].r); else if(t3[hash].p <= (unsigned char) j) ++(t3[hash].r); 
                                        t3[hash].first = 1;
                                        n_rep = (t3[hash].r);
                                        t3[hash].p = (unsigned char) (j+3);
                                        //if(n_rep == 4){ printf("%s\n", curr_kmer); printf("Got longest repeat of size 3"); getchar(); }
                                        if(n_rep == 4) death = 2;
                                    } break;
                                    case 4: { 
                                        hash = hashOfWord(&curr_kmer[j], 4, 0);
                                        //printf("(T4)hash: %"PRIu64" ->attempting insert %.4s at %u prior pos is (%u) result pos is = %u; curr_rep = %u\n", hash, &curr_kmer[j], j, t5[hash].p, t5[hash].p + 4, t5[hash].r);
                                        if(t4[hash].first == 0) n_rep = ++(t4[hash].r); else if(t4[hash].p <= (unsigned char) j) ++(t4[hash].r); 
                                        t4[hash].first = 1;
                                        n_rep = (t4[hash].r);
                                        t4[hash].p = (unsigned char) (j+4);
                                        //if(n_rep == 3){ printf("%s\n", curr_kmer); printf("Got longest repeat of size 4"); getchar(); }
                                        if(n_rep == 3) death = 2;
                                    } break;
                                    case 6: { 
                                        hash = hashOfWord(&curr_kmer[j], 6, 0);
                                        //printf("(T6)hash: %"PRIu64" ->attempting insert %.6s at %u prior pos is (%u) result pos is = %u; curr_rep = %u\n", hash, &curr_kmer[j], j, t5[hash].p, t5[hash].p + 4, t5[hash].r);
                                        if(t6[hash].first == 0) n_rep = ++(t6[hash].r); else if(t6[hash].p <= (unsigned char) j) ++(t6[hash].r); 
                                        t6[hash].first = 1;
                                        n_rep = (t6[hash].r);
                                        t6[hash].p = (unsigned char) (j+6);
                                        //if(n_rep == 2){ printf("%s\n", curr_kmer); printf("Got longest repeat of size 6"); getchar(); }
                                        if(n_rep == 2) death = 2;
                                    } break;
                                }
                            }
                            

                            if(death == 1){
                                i=0; j=0;
                                while(j<10){
                                    if(huge_seq[write_pos-i-2] != '\n') huge_seq[write_pos-i-2] = 'N';
                                    if(huge_seq[write_pos-i-2] != '\n') ++j;
                                    i++;
                                }
                            } 
                            if(death == 2){
                                i=0; j=0;
                                while(j<12){ 
                                    if(huge_seq[write_pos-i] != '\n') huge_seq[write_pos-i] = 'N';
                                    if(huge_seq[write_pos-i] != '\n') ++j;
                                    i++;
                                }
                            } 
                            if(j + i >= 12) break;
                            
                        }
                    }
                    
                    /*
                    printf("For seq: %s\n", curr_kmer);
                    for(i=0; i<4; i++){
                        printf("%u\n", t1[i]);
                    }
                    printf("--------\n");
                    for(i=0; i<16; i++){
                        printf("%u\n", t2[i]);
                    }
                    */
		
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

    fwrite(huge_seq, 1, write_pos, out_database);


    fprintf(stdout, "[INFO] Database loaded and of length %"PRIu64".\n", current_len);
    fclose(database);

    
    if(out_database != NULL) fclose(out_database);
    
    return 0;
}

void init_args(int argc, char ** av, FILE ** database, unsigned char * cutoff){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           lowcompfilter -db [database]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -cut       [Integer:   c>1 (default 5)]\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }
        
        if(strcmp(av[pNum], "-db") == 0){
            *database = fopen64(av[pNum+1], "rt");
            if(database==NULL) terror("Could not open database file");
        }
        if(strcmp(av[pNum], "-cut") == 0){
            *cutoff = (unsigned char) atoi(av[pNum+1]);
            if(database==NULL) terror("Could not open database file");
        }

        pNum++;
    }

}

