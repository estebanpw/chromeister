#define QF_LAMBDA 0.275
#define QF_KARLIN 0.333

typedef struct container{
    llpos * table[4][4][4][4][4][4][4][4][4][4][4][4];
} Container;



typedef struct {
    uint64_t id;        //The thread id
    SeqInfo * database; //Database sequence and lengths
    SeqInfo * query;    //Query sequence and lengths
    uint64_t from;      //Starting READ to compute alignments from
    uint64_t to;        //End READ to compute alignments from
    Container * container; //Container to hold the multidimensional array
    uint64_t accepted_query_reads; //Number of reads that have a fragment with evalue less than specified
    long double min_e_value;    //Minimum evalue to accept read
    long double min_coverage;    //Minimum coverage percentage to accept read
    long double min_identity;    //Minimum identity percentage to accept read
    long double window;         //Percentage of window that will be explored (+-)
    FILE * out; //File to write alignments out
    int igap;
    int egap;
    uint64_t * hits;        // To work in hits mode only
    struct positioned_cell * mc;
    struct cell ** table;
    char * reconstruct_X;
    char * reconstruct_Y;
    char * writing_buffer_alignment;
    unsigned char * my_x;
    unsigned char * my_y;
    Head * queue_head;  //To tell where the queue starts after modifications
    pthread_mutex_t * lock;
    unsigned char full_comp; // Tells whether read reporting should stop at first match or keep reporting
    unsigned char * markers; // To tell which sequences were already used
} HashTableArgs;


/*
    Nucleotides matching function
*/
int64_t compare_letters(unsigned char a, unsigned char b);


/**

Reset a number of llpos from a pool

*/
void reset_llpos(Mempool_l * mp, uint64_t * n_pools_used, uint64_t n_reset_llpos);
/**
 * Initialize the memory pool to later retrieve individual memory addresses for llpos
 * 
 */
void init_mem_pool_llpos(Mempool_l * mp);

/**
 * Get a new memory address from the pool mp for a type llpos
 * 
 */
llpos * getNewLocationllpos(Mempool_l * mp, uint64_t * n_pools_used);

/*
llpos for kmers
*/

llpos_kmer * getNewLocationllpos_kmer(Mempool_l_kmer * mp, uint64_t * n_pools_used);

void init_mem_pool_llpos_kmer(Mempool_l_kmer * mp);