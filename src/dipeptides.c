/*********

File        dipeptides.c
Author      EPW <estebanpw@uma.es>
Description Calculates frequencies of dipeptides

USAGE       Usage is described by calling ./dipeptides --help



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

uint64_t custom_kmer = 2; // Defined as external in structs.h
uint64_t diffuse_z = 4; // Defined as external in structs.h

void init_args(int argc, char ** av, FILE ** database);

int main(int argc, char ** av){
    

    uint64_t dipeptides[16];

    uint64_t i;
    for(i=0;i<16;i++){
        dipeptides[i] = 0;
    }

    //query to read kmers from, database to find seeds
    FILE * database = NULL;
    
    
    init_args(argc, av, &database);



    uint64_t char_converter[91];
    char_converter[(unsigned char)'A'] = 0;
    char_converter[(unsigned char)'C'] = 1;
    char_converter[(unsigned char)'G'] = 2;
    char_converter[(unsigned char)'T'] = 3;

    
    //Variables to account for positions
    //Print info
    //fprintf(stdout, "[INFO] Loading sequence\n");
    //Variables to read kmers
    char c = 'N'; //Char to read character
    //Current length of array and variables for the buffer
    uint64_t idx = 0, r = 0;
    
    //Vector to read in batches
    char * temp_seq_buffer = NULL;
    if ((temp_seq_buffer = calloc(READBUF, sizeof(char))) == NULL) {
        terror("Could not allocate memory for read buffer");
    }


    fseek(database, 0, SEEK_END);
    uint64_t aprox_len_query = ftell(database);
    rewind(database);

    uint64_t a_hundreth = (aprox_len_query/100);

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

    long double total = 0;

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
                        //fprintf(stdout, "\r%"PRIu64"%%...", 1+100*current_len/aprox_len_query); 
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
                if(word_size == 2){
                    //write to hash table
                    

                    dipeptides[char_converter[curr_kmer[1]] * 4 + char_converter[curr_kmer[0]]]++;
                    total += 1.0;
                    
                    
		
		            // Overlapping
                    memmove(&curr_kmer[0], &curr_kmer[1], 1);
                    --word_size;
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
    //fprintf(stdout, "[INFO] Database loaded and of length %"PRIu64".\n", current_len);
    //close database
    fclose(database);

    //fprintf(stdout, "[INFO] Frequency profile generated.\n");

    for(i=0;i<16;i++){
        fprintf(stdout, "$\t%Le\t->\t%"PRIu64"\n", (long double)dipeptides[i]/(long double)total, dipeptides[i]);
    }
    
    //begin = clock();
    

    
    
    
    
    return 0;
}

void init_args(int argc, char ** av, FILE ** database){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           dipeptides -db [sequence]\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }
        
        if(strcmp(av[pNum], "-db") == 0){
            *database = fopen64(av[pNum+1], "rt");
            if(database==NULL) terror("Could not open database file");
        }
        
        
        pNum++;
    }
    
    if(*database==NULL) terror("A sequence is required");
}

