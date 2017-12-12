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



int64_t compare_letters(unsigned char a, unsigned char b){
    if(a != (unsigned char) 'N' && a != (unsigned char) '>') return (a == b) ? POINT : -POINT;
    return -POINT;
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

AVLTree * getNewLocationAVLTree(Mempool_AVL * mp, uint64_t * n_pools_used, uint64_t key){

    if(mp[*n_pools_used].current == POOL_SIZE){
        *n_pools_used += 1;
        if(*n_pools_used == MAX_MEM_POOLS) terror("Reached max pools");
        init_mem_pool_AVL(&mp[*n_pools_used]);
        
    }

    AVLTree * new_pos = mp[*n_pools_used].base + mp[*n_pools_used].current;
    mp[*n_pools_used].current++;

    new_pos->key = key;
    new_pos->height = 1;

    return new_pos;
}

void init_mem_pool_AVL(Mempool_AVL * mp){
    mp->base = (AVLTree *) calloc(POOL_SIZE, sizeof(AVLTree));
    if(mp->base == NULL) terror("Could not request memory pool");
    mp->current = 0;
}



/*
// An AVL tree node
typedef struct AVL_Node{
    uint64_t key;
    struct AVL_Node * left;
    struct AVL_Node * right;
    uint64_t height;
    llpos * next;
} AVLTree;
*/
 
// A utility function to get height of the tree
/*
uint64_t height(AVLTree * N){
    if (N == NULL)
        return 0;
    return N->height;
}
*/
/* Substituted by (x == NULL) ? (0) : (x->height) */
 
/* Helper function that allocates a new node with the given key and
    NULL left and right pointers. */

/* This one is substituted by AVLTree * getNewLocationAVLTree(Mempool_AVL * mp, uint64_t * n_pools_used, uint64_t key) */
 
// A utility function to right rotate subtree rooted with y
// See the diagram given above.
AVLTree * right_rotate(AVLTree * y){
    AVLTree * x = y->left;
    AVLTree * T2 = x->right;
 
    // Perform rotation
    x->right = y;
    y->left = T2;
 
    // Update heights
    x->height = MAX((x == NULL) ? (0) : (x->left->height), (x == NULL) ? (0) : (x->right->height))+1;
    y->height = MAX((y == NULL) ? (0) : (y->left->height), (y == NULL) ? (0) : (y->right->height))+1;
 
    // Return new root
    return x;
}
 
// A utility function to left rotate subtree rooted with x
// See the diagram given above.
AVLTree * left_rotate(AVLTree * x){
    AVLTree * y = x->right;
    AVLTree * T2 = y->left;
 
    // Perform rotation
    y->left = x;
    x->right = T2;
 
    //  Update heights
    x->height = MAX((x == NULL) ? (0) : (x->left->height), (x == NULL) ? (0) : (x->right->height))+1;
    y->height = MAX((y == NULL) ? (0) : (y->left->height), (y == NULL) ? (0) : (y->right->height))+1;
 
    // Return new root
    return y;
}
 
// Get Balance factor of node N
/*
int64_t get_balance(AVLTree * N){
    if (N == NULL)
        return 0;
    return height(N->left) - height(N->right);
}
*/
/* Substituted by (node == NULL) ? (0) : ((int64_t) node->left->height - (int64_t) node->right->height) */
 
// Recursive function to insert key in subtree rooted
// with node and returns new root of subtree.
AVLTree * insert_AVLTree(AVLTree * node, uint64_t key, Mempool_AVL * mp, uint64_t * n_pools_used){
    /* 1.  Perform the normal BST insertion */
    if (node == NULL)
        return(getNewLocationAVLTree(mp, n_pools_used, key));
 
    if (key < node->key)
        node->left  = insert_AVLTree(node->left, key, mp, n_pools_used);
    else if (key > node->key)
        node->right = insert_AVLTree(node->right, key, mp, n_pools_used);
    else // Equal keys are not allowed in BST
        return node;
 
    /* 2. Update height of this ancestor node */
    node->height = 1 + MAX((node->left == NULL) ? (0) : (node->left->height), (node->right == NULL) ? (0) : (node->right->height));
 
    /* 3. Get the balance factor of this ancestor
          node to check whether this node became
          unbalanced */
    int64_t balance = (node->left == NULL || node->right == NULL) ? (0) : ((int64_t) node->left->height - (int64_t) node->right->height);
 
    // If this node becomes unbalanced, then
    // there are 4 cases
 
    // Left Left Case
    if (balance > 1 && key < node->left->key)
        return right_rotate(node);
 
    // Right Right Case
    if (balance < -1 && key > node->right->key)
        return left_rotate(node);
 
    // Left Right Case
    if (balance > 1 && key > node->left->key)
    {
        node->left =  left_rotate(node->left);
        return right_rotate(node);
    }
 
    // Right Left Case
    if (balance < -1 && key < node->right->key)
    {
        node->right = right_rotate(node->right);
        return left_rotate(node);
    }
 
    /* return the (unchanged) node pointer */
    return node;
}
 
// A utility function to print preorder traversal
// of the tree.
// The function also prints height of every node

void pre_order(AVLTree * root){
    if(root != NULL){
        printf("%"PRIu64" ", root->key);
        pre_order(root->left);
        pre_order(root->right);
    }
}
