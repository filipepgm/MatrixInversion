#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <pthreads.h>
#include <omp.h>
#include <mpi.h>

#define TRUE 1
#define FALSE 0

#define MAX_LINE 1024
#define MAX_RANDOM 2147483647

#define CONF_FILE "conf_inv_mat_MC"
#define DISTRIB_FILE "distrib.conf"
//#define MAXSIZE 100000000
#define MAXSIZE (1024*1024*1024)
#define DEFAULT_N_STEPS 30
#define MAX_STEPS 100
#define DEFAULT_N_PLAYS 1000
#define MAX_PLAYS 1000000
#define RAND_SEED 0                     // time(NULL)

#define fltype float
#define TYPE_FLOAT MPI_FLOAT

#define M_MAT 111
#define M_NROWS 112
#define M_NCOLS 113
#define M_ROWMIN 114
#define M_ROWMAX 115
#define M_COLMIN 116
#define M_COLMAX 117
#define TOKEN 120
#define ALLDONE 121
#define SAMPLE 122
#define PRINT_IN_TURN 123
#define M_SIZES 130
#define M_COLNUMS 131
#define M_VALUES 132


#define BLOCK_LOW(id,p,n)  ((unsigned long long)(id)*(unsigned long long)(n)/(unsigned long long)(p))
//#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
//#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
//#define BLOCK_OWNER(index,p,n) ((((unsigned long long) p)*((unsigned long long)( (unsigned long long) index)+1)-1)/( (unsigned long long)n))

#define OWNED_ROW(index) ( ( (index) >= taskrowmin ) && ( (index) <= taskrowmax ) )

typedef struct _playInfo {
    int startRow;
    int curRow;
    int stepsLeft;
    fltype value;
    fltype sum;
} playInfo;
