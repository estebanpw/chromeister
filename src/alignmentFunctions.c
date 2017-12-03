#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))


int64_t compare_letters(unsigned char a, unsigned char b){
    if(a != (unsigned char) 'N' && a != (unsigned char) '>') return (a == b) ? POINT : -POINT;
    return -POINT;
}

void reset_llpos(Mempool_l * mp, uint64_t * n_pools_used, uint64_t n_reset_llpos){
    if(mp[*n_pools_used].current >= n_reset_llpos){
        mp[*n_pools_used].current -= n_reset_llpos;
    }else{
        if(*n_pools_used == 0){ mp[*n_pools_used].current = 0; return; }

        n_reset_llpos -= mp[*n_pools_used].current;
        --(*n_pools_used);
        mp[*n_pools_used].current -= n_reset_llpos;
    }
}

llpos * getNewLocationllpos(Mempool_l * mp, uint64_t * n_pools_used){

    if(mp[*n_pools_used].current == POOL_SIZE){
        *n_pools_used += 1;
        if(*n_pools_used == MAX_MEM_POOLS) terror("Reached max pools");
        init_mem_pool_llpos(&mp[*n_pools_used]);
        
    }

    llpos * new_pos = mp[*n_pools_used].base + mp[*n_pools_used].current;
    mp[*n_pools_used].current++;

    
    return new_pos;
}

void init_mem_pool_llpos(Mempool_l * mp){
    mp->base = (llpos *) calloc(POOL_SIZE, sizeof(llpos));
    if(mp->base == NULL) terror("Could not request memory pool");
    mp->current = 0;
}

/*
typedef struct linked_list_kmer{
    unsigned char kmer[32];
    linked_list_kmer * next;
} llpos_kmer;

typedef struct static_table{
    llpos_kmer * table;
} Static_table;


*/

llpos_kmer * getNewLocationllpos_kmer(Mempool_l_kmer * mp, uint64_t * n_pools_used){

    if(mp[*n_pools_used].current == POOL_SIZE){
        *n_pools_used += 1;
        if(*n_pools_used == MAX_MEM_POOLS) terror("Reached max pools");
        init_mem_pool_llpos_kmer(&mp[*n_pools_used]);
        
    }

    llpos_kmer * new_pos = mp[*n_pools_used].base + mp[*n_pools_used].current;
    mp[*n_pools_used].current++;

    
    return new_pos;
}

void init_mem_pool_llpos_kmer(Mempool_l_kmer * mp){
    mp->base = (llpos_kmer *) calloc(POOL_SIZE, sizeof(llpos_kmer));
    if(mp->base == NULL) terror("Could not request memory pool");
    mp->current = 0;
}