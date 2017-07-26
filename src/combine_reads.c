#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "structs.h"
#include "commonFunctions.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))
#define STARTING_SEQS 1000
#define PIECE_OF_DB_REALLOC 3200000 //half a gigabyte if divided by 8 bytes
#define MATVAL 101

int main(int argc, char ** av){
    
    if(argc != 2) terror("USE:: combine_reads <file>");

    FILE * results, * data;
    data = fopen(av[1], "rt");
    if(data == NULL) terror("Could not open input file");

    uint64_t * mat = (uint64_t *) calloc(MATVAL, sizeof(uint64_t));
    if(mat == NULL) terror("Could not allocate matrix array");

    char buffer[MAXLID];
    if ((results = fopen("accu.log", "r")) == NULL){
        results = fopen("accu.log", "wt");
    }else{
        // Load the matrix
        uint64_t i;
        
        for(i=0;i<100;i++){
            if(0 == fgets(buffer, MAXLID, results)) terror("Missing number on load");
            
            //fprintf(stdout, "Have %s\n", buffer);
            buffer[strlen(buffer)-1] = '\0';
            mat[i] = asciiToUint64(buffer);
            //fprintf(stdout, "%"PRIu64"\n", mat[i]);
            //getchar();
        }
        fclose(results);
        results = fopen("accu.log", "wt"); // Re open
    }

    // Read file 
    uint64_t read_id_1, read_id_2, coverage, identity, length, current, currmax, j;
    while(!feof(data)){
        if(0 == fgets(buffer, MAXLID, data) && !feof(data)) terror("Missing values");
        // 2 77277 89 64 213
        sscanf(buffer, "%"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64, &read_id_1, &read_id_2, &coverage, &identity, &length);
        //fprintf(stdout, "%"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64"\n", read_id_1, read_id_2, coverage, identity, length);
        currmax = MIN(identity, coverage);
        //fprintf(stdout, "%"PRIu64"\n", currmax);
        current = read_id_1;
        /*
        for(j=currmax; j > 1; j--){
            if(current != lasts[j]){
                mat[j]++;
                lasts[j] = current;
            }
        }
        */
        mat[currmax] += 1;

    }

    for(j=99; j>0; j--){
        mat[j] += mat[j+1];
    }
    mat[0] = mat[1];

    for(j=0; j<100; j++){
        fprintf(stdout, "%"PRIu64"\n", mat[j]);
        fprintf(results, "%"PRIu64"\n", mat[j]);
    }


    fclose(results);
    fclose(data);
    free(mat);
    return 0;
}