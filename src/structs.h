#ifndef STRUCTS_H
#define STRUCTS_H

#include <inttypes.h>

//Structs required for the dotplot workflow
#define MAXLID 200
//#define READBUF 2000000000 //2 GB
#define READBUF 500000000 //50MB
#define INITSEQS 3000 //Number of starting sequences (in database)
#define POINT 4

#define FIXED_K 12

#define MAXLID 200
#define ALIGN_LEN 60 //For NW alignment
#define MAX_READ_SIZE 20000 //Maximum length of read to have a portion of the table allocated
#define MAX_WINDOW_SIZE 500 //Maximum window length to explore NW table
//#define POOL_SIZE 2500000000 // by 16 bytes it is 40 GB
#define POOL_SIZE 12500000 // 1 GB if 16 bytes
#define MAX_MEM_POOLS 256 

#define NHASH 6400

#define BYTES_IN_MER 4
#define MAX_DECOMP_HASH 10
#define FALSE 0
#define TRUE 1

extern uint64_t custom_kmer;

extern uint64_t hash_functions[NHASH];

typedef struct linked_list_kmer{
    unsigned char kmer[32];
    uint64_t hash;
    uint64_t pos;
    struct linked_list_kmer * next;
} llpos_kmer;

//Struct for linked list of positions
typedef struct linked_list_pos{
    uint64_t pos;
    struct linked_list_pos * next;
} llpos;

typedef struct def_pos{
    uint64_t full_hash;
    llpos * next;
} Def_pos;

typedef struct static_table{
    llpos_kmer * table;
    uint64_t * count;
    uint64_t * picked_hash;
    Def_pos ** second_level;

} Static_table;


//Struct for memory pool por lists
typedef struct mempool_l{
    llpos * base;
    uint64_t current;
} Mempool_l;

// Mempool for lists of kmers
typedef struct mempool_l_kmer{
    llpos_kmer * base;
    uint64_t current;
} Mempool_l_kmer;

//Struct for a whole sequence(s) data
typedef struct seqinfo{
    unsigned char * sequences;
    uint64_t * start_pos;
    uint64_t total_len;
    uint64_t n_seqs;
} SeqInfo;

//Struct for the alignment of a quick hit
typedef struct quickfrag{
    uint64_t x_start;
    uint64_t y_start;
    uint64_t t_len;
    long double coverage;
    long double e_value;
} Quickfrag;

typedef struct point{
    uint64_t x;
    uint64_t y;
} Point;



struct cell{
    int64_t score;
    uint64_t xfrom;
    uint64_t yfrom;
};

struct positioned_cell{
    int64_t score;
    uint64_t xpos;
    uint64_t ypos;
};

struct best_cell{
    struct positioned_cell c;
    uint64_t j_prime;
};

typedef struct{
    uint64_t identities;
    uint64_t length;
    uint64_t igaps;
    uint64_t egaps;
} BasicAlignment;

typedef struct queue{
    uint64_t r1; //reads region
    uint64_t r2;
    struct queue * next;
} Queue;

typedef struct{
    Queue * head;
} Head;


#endif
