#include "matInvMCBinMPI.h"
//#define DEBUG
#define OUTBUFFERSIZE 2048
#define INQUEUESIZE 1024

int nsteps = DEFAULT_N_STEPS;
int nplays = DEFAULT_N_PLAYS;
int rowMin, rowMax, colMin, colMax;
int pid, ntasks, nrows, ncolumns, ntaskrows, taskrowmin, taskrowmax;

int bufferSize = 0;
int *sizeList, *colNums;
float *valueList;

volatile int nthreads_idle;
volatile char done = 0, all_samples_given = 0, token_holder = 0;
long long doneMsg = 0, totalDones;
long long **doneMsgCounter = NULL;
char msg[MAX_LINE];
fltype *V, **count;
int nthreads = 0;

typedef struct thread_public_info {
    volatile int to_fill;
    volatile int to_process;
    volatile int waiting;
    playInfo * in_sample;
    MPI_Status *status;
    MPI_Request req_sample;
} thread_public_info;

thread_public_info **thread_in;

typedef struct out_struct {
    playInfo **outArray;
    int *outFilled;
    char *outBit;
    char *outWaiting;
    MPI_Request *outRequests;
} out_struct; 

out_struct global_out;

typedef struct thread_info {
    unsigned int rand_val;
    int gen_r, gen_play;
    out_struct out;
    char thread_gen_done;
    int cache_row; float cache_value;
    long long thread_messages_done;
} thread_info, *lp_thread_info;

long long log_full_rcv_buffer = 0;

void printSample(playInfo *sample)
{
    printf("start=%d, cur=%d, stepsLeft=%d, value=%f",
	   sample->startRow, sample->curRow, sample->stepsLeft,
	   sample->value);
}

void cleanExit(char *msg, int code)
{
    fprintf(stderr, "%s", msg);
    MPI_Abort(MPI_COMM_WORLD, code);
    exit(code);
}


//===========================================================
// Uneven distribution
int *row_distribution;

void parse_distrib_file(FILE* fpc) {
    int t;
    row_distribution = malloc( (ntasks+1)*sizeof(int));
    for( t = 0; t <= ntasks;  t++) {
        if(1 != fscanf(fpc, "%d", &(row_distribution[t])) ) {
            sprintf(msg, "Error in distrib file: %d\n", t);
            cleanExit(msg, -3);
        }
    }
    fclose(fpc);
}

void default_distrib() {
    int t;
    row_distribution = malloc((ntasks+1)*sizeof(int) );
    for( t = 0; t <= ntasks;  t++) {
        row_distribution[t] = BLOCK_LOW(t, ntasks, nrows);
    }
}

int block_owner(int row) {
    int t =0;
    for (t = 0; t < ntasks; t++) {
        if (row_distribution[t+1] > row)
            return t;
    }
    return -1;
}

//===========================================================

void parse_conf_file(FILE *fpc)
{
    char line[MAX_LINE];
    
    // ignore 1st line
    fgets(line, MAX_LINE, fpc);
    
    // number of plays per step
    fscanf(fpc, "%d", &nplays);
    if((nplays < 1) || (nplays >= MAX_PLAYS)){
        sprintf(msg, "Error in number of plays in configuration file, %d.\nMust be a natural number below %d.\n", nplays, MAX_PLAYS);
        cleanExit(msg, -2);
    }
    
    // number of steps
    fscanf(fpc, "%d", &nsteps);
    if((nsteps < 1) || (nsteps >= MAX_STEPS)){
        sprintf(msg, "Error in number of steps in configuration file, %d.\nMust be a natural number below %d.\n", nsteps, MAX_STEPS);
        cleanExit(msg, -3);
    }
    
    fclose(fpc);
}

// Allocates and returns the array where to store the results
fltype **countMatMalloc(int nthreads)
{
    fltype **m1 = calloc(nthreads, sizeof(fltype*));
    if((m1 == NULL)){
        sprintf(msg, "Error allocating count matrix memory space.\n");
        cleanExit(msg, -4);
    }
    
    /*for(t=0; t<nthreads; t++) {
        m1[t] = calloc(nrows, sizeof(fltype));
        if((m1[t] == NULL)){
            sprintf(msg, "Error allocating count matrix memory space.\n");
            cleanExit(msg, -4);
        }
    }*/
    return m1;
}

void mpi_receive_and_check(void* buffer, int size, MPI_Datatype type, int origin, int tag){
    MPI_Status status;
    int count;

    MPI_Recv(buffer, size, type, origin, tag, MPI_COMM_WORLD, &status);

    MPI_Get_count(&status, type, &count);
    if (count != size)
    {
        sprintf(msg, "Ooops, task %d expected %d elements and received %d.\n",
                pid, size, (int)(count));
        cleanExit(msg, -5);
    }
}

FILE * parse_arguments(int argc, char *argv[]){
    FILE *fpm;
	int scanres;
	
    if((argc < 2) || (argc > 4)){
        sprintf(msg, "Usage: %s \"matrix file\" <line to compute> <column to compute>\n", argv[0]);
        cleanExit(msg, -1);
    }
    
    // 1st argument: matrix file
    fpm = fopen(argv[1], "r");
    if(fpm == NULL){
        sprintf(msg, "Error reading file %s.\n", argv[1]);
        cleanExit(msg, -2);
    }
    
    // two integers, #rows #columns
    scanres = fscanf(fpm, "%d %d", &nrows, &ncolumns);
    if(scanres != 2 || (nrows < 1) || (nrows >= MAXSIZE)){
        sprintf(msg, "Error in number of rows, %d. Must be a natural number below %d.\n", nrows, MAXSIZE);
        cleanExit(msg, -2);
    }
    if((ncolumns < 1) || (ncolumns >= MAXSIZE)){
        sprintf(msg, "Error in number of columns, %d. Must be a natural number below %d.\n", ncolumns, MAXSIZE);
        cleanExit(msg, -2);
    }
    
    if(argc > 2){
        rowMin = rowMax = atoi(argv[2]);
        if((rowMin < 0) || (rowMin > nrows)){
            sprintf(msg, "Error in row to compute, %d. Must be a natural number below %d.\n", rowMin, nrows);
            cleanExit(msg, -2);
        }
    }
    else{
        rowMin = 0;
        rowMax = nrows - 1;
    }
    
    if(argc == 4){
        colMin = colMax = atoi(argv[3]);
        if((colMin < 0) || (colMin > ncolumns)){
            sprintf(msg, "Error in column to compute, %d. Must be a natural number below %d.\n", colMin, ncolumns);
            cleanExit(msg, -2);
        }
    }
    else{
        colMin = 0;
        colMax = ncolumns - 1;
    }

    FILE *fpc;
    // check if distrib file exists
    fpc = fopen(DISTRIB_FILE, "r");
    if(fpc!= NULL)
        parse_distrib_file(fpc);
    else
        default_distrib();

    taskrowmin = row_distribution[pid];//BLOCK_LOW(pid, ntasks, nrows);
    taskrowmax = row_distribution[pid+1] - 1 ;//BLOCK_HIGH(pid, ntasks, nrows);
    ntaskrows = taskrowmax - taskrowmin + 1;//BLOCK_SIZE(pid,ntasks,nrows);

    return fpm;
}


//============================================================
//===================   READ INPUT   =========================
//============================================================
int addValue(int row, int col, float value, int bufferFilled) {
    if(bufferSize == bufferFilled){
        bufferSize *= 2;
        colNums = realloc(colNums, sizeof(int)*(bufferSize));
        valueList = realloc(valueList, sizeof(fltype)*(bufferSize));
        if( !colNums || !valueList) {
            sprintf(msg, "Error reallocating memory space for matrix\n");
            cleanExit(msg, -13);
        }
    }

    colNums[bufferFilled] = col;
    valueList[bufferFilled] = value;
    bufferFilled++;
    sizeList[row-taskrowmin+1] = bufferFilled;
    return bufferFilled;
}

int load_from_partial_file(FILE * pfpm, int bufferFilled, int start, int end) {
    int row=-100, col=-100;
    float value;
    int lastRow = -100;

    int scanres = fscanf(pfpm, "%d %d %f\n", &row, &col, &value);
    row--; col--;
    if (start < taskrowmin)
        start = taskrowmin;

    while(scanres == 3 && row <= taskrowmax){
        if(value != 0.0 && row >=taskrowmin) {
            //================================================
            // Independent of input file format
            if(row != lastRow){
                if(start == lastRow) // lastRow was read
                    start++;
                while( start < row){ // Add a 0 to every non read line
                    sizeList[start-taskrowmin] = bufferFilled;
                    bufferFilled = addValue(start, start, 0.0, bufferFilled);
                    start++;
                }
                sizeList[row-taskrowmin] = bufferFilled;
			}
            bufferFilled = addValue(row, col, value, bufferFilled);
			lastRow = row;
            //================================================
        }
		scanres = fscanf(pfpm, "%d %d %f\n", &row, &col, &value);
        row--; col--;
    }
    //sizeList[lastRow-taskrowmin+1] = bufferFilled;
    if(start == lastRow) // lastRow was read
        start++;
    while (start <= taskrowmax && start < end) { //Add a 0 for every non read line
        sizeList[start-taskrowmin] = bufferFilled;
        bufferFilled = addValue(start, start, 0.0, bufferFilled);
        start++;
    }
    return bufferFilled;
}

// Input files start at 1, internal representation starts at 0!!!!!
int load_from_index(FILE * fpm) {
    int start, end, declaredSize, scanres;
    char filename[1024];
    int bufferFilled = 0;

    #ifdef DEBUG
    printf("P%d.main: reading from file, rows %d to %d\n", pid, taskrowmin, taskrowmax);
    fflush(stdout);
    #endif

    while( 3 == fscanf(fpm, "%d %d %1023s", &start, &end, filename ) ) {
        start--; end--;
        if(taskrowmax < start) {
            break; // No more interesting files
        }
        if(start <= taskrowmax && end > taskrowmin) {
            // Read this file
            FILE * pfpm;
            pfpm = fopen(filename, "r");
            scanres = fscanf(pfpm, "%d %d", &declaredSize, &declaredSize);
            if(scanres != 2 || declaredSize != nrows)
                fprintf(stderr, "ERROR: file declares wrong size! %s", filename);
            bufferFilled  = load_from_partial_file(pfpm, bufferFilled, start, end );
            fclose(pfpm);
        }
    }
    return bufferFilled;
}

void read_matrix(int argc, char *argv[]) {
    FILE *fpm;

    fpm = parse_arguments(argc, argv);

    bufferSize = ntaskrows;
    sizeList = malloc(sizeof(int)*(ntaskrows+2));
    colNums = malloc(sizeof(int)*(bufferSize));
    valueList = malloc(sizeof(fltype)*(bufferSize));

    if( !sizeList || !colNums || !valueList) {
        sprintf(msg, "Error allocating memory space for matrix\n");
        cleanExit(msg, -12);
    }

    load_from_index(fpm);
}


fltype *normalize_accumulate()
{
    int row, col;
	fltype fValue;
    fltype *V = malloc(ntaskrows * sizeof(fltype));
    if(V == NULL){
        sprintf(msg, "Error allocating memory space for normalization vector V.\n");
        cleanExit(msg, -4);
    }

    for(row = 0; row < ntaskrows; row++){

#ifdef DEBUG
	printf("P%d: 1, Processing row: %d\n", pid, row);
	fflush(stdout);
#endif

        V[row] = 0.0;
		for(col= sizeList[row]; col < sizeList[row+1]; col++)
			V[row] += valueList[col];

#ifdef DEBUG
	printf("P%d: V[%d] = %f\n", pid, row, V[row]);
	fflush(stdout);
#endif
        fValue = 0.0;
		for(col= sizeList[row]; col < sizeList[row+1]; col++){
            if(V[row] != 0.0)
                fValue = valueList[col] / V[row] + fValue;
            valueList[col]  = fValue;
        }
		valueList[col-1] = 1.0; // Avoid float imprecision
    }
    return V;
}

// Change to sparse
int selectColumn(int curRow, int rndVal, int ncolumns) {
	int nLow = sizeList[curRow], nUp = sizeList[curRow+1], nMid;
    int select = rndVal % (nUp-nLow);
    /*
	while(valueList[nLow] < rndVal)
	{
		nMid = (nUp + nLow) >>1;
		if (valueList[nMid] < rndVal)
			nLow = nMid+1;
		else
			nUp = nMid;
	}*/
	return colNums[nLow + select];
}


//============================================================
//===================   PROCESSOR    =========================
//============================================================
int gen_row=0;
playInfo *get_next_sample(lp_thread_info data, playInfo *sample)
{
#ifdef DEBUG
    printf("P%d.main: GET sample", pid);
    fflush(stdout);
#endif

    if (data->gen_play >= nplays) {
        // Get new row
        data->gen_r += nthreads-1;
        if(data->gen_r > taskrowmax)
            data->thread_gen_done = TRUE;
        
        /*#pragma omp critical (gen_row)
        {
            if(all_samples_given) {
                data->thread_gen_done = TRUE;
            } else {
                data->gen_r = gen_row;
                gen_row++;
                if(gen_row > taskrowmax)
                    all_samples_given = TRUE;
            }
        }*/
        data->gen_play = 0;
    }

    if (!data->thread_gen_done) {
        sample->startRow = sample->curRow = data->gen_r;
        sample->stepsLeft = nsteps;
        sample->value = 1.0; sample->sum = 0.0;
        data->gen_play++;
        return sample;
    } else {
        return NULL;
    }
}

inline void send_to(out_struct* out, int owner) {
    MPI_Status status;
    int outbufferIndex = owner*2+out->outBit[owner];
    int pos = out->outFilled[owner];

    if(out->outWaiting[owner])
        MPI_Wait(&(out->outRequests[owner]), &status); // Wait before placing the next request

    MPI_Isend(out->outArray[outbufferIndex], pos*sizeof(playInfo), MPI_BYTE,
            owner, SAMPLE, MPI_COMM_WORLD, &(out->outRequests[owner]));
    //MPI_Send(out->outArray[outbufferIndex], pos*sizeof(playInfo), MPI_BYTE,
            //owner, SAMPLE, MPI_COMM_WORLD);
    out->outFilled[owner] = 0;
    out->outBit[owner] ^= 1; // Flip the bit 1 <-> 0
    out->outWaiting[owner] = 1;
}

void push_out_struct(out_struct* out, playInfo * sample, int owner) {
    int pos = out->outFilled[owner];
    int outbufferIndex = owner*2+out->outBit[owner];
    
    out->outArray[outbufferIndex][pos].startRow = sample->startRow;
    out->outArray[outbufferIndex][pos].curRow = sample->curRow;
    out->outArray[outbufferIndex][pos].stepsLeft= sample->stepsLeft;
    out->outArray[outbufferIndex][pos].value= sample->value;
    out->outArray[outbufferIndex][pos].sum= sample->sum;
    out->outFilled[owner]++;
}

void send_sample(out_struct* out, playInfo * sample, int owner){
    {
        // Add new entry
        push_out_struct(out, sample, owner);

        // Send if needed
        if( ! (out->outFilled[owner] < OUTBUFFERSIZE)){
            send_to(out, owner);
        }
    }
}

inline void flush_thread_cache(int threadId, lp_thread_info data) {
    if (data->cache_row != -1) {
        //#pragma omp atomic
        count[threadId][data->cache_row - taskrowmin] += data->cache_value;
    }
}

void compute_sample(int threadId, lp_thread_info data, playInfo *sample)
{
    out_struct* out = & (data->out);
    int rndVal;
    fltype partialSum = 0.0;
#ifdef DEBUG
    printf("P%d.main: handling sample: ", pid);
    printSample(sample);
    fflush(stdout);
#endif
    while (OWNED_ROW(sample->curRow)  && (sample->stepsLeft > 0))
    {
        #ifdef DEBUG
        //printf("P%d.main: accounting sample: ", pid);
        //printSample(sample);
        //fflush(stdout);
        #endif
        sample->sum += sample->value;

        sample->value *= V[sample->curRow - taskrowmin];
        rndVal = rand_r(& (data->rand_val));//random(); //rndVal = ((fltype) random()) / MAX_RANDOM;
        sample->curRow = selectColumn(sample->curRow - taskrowmin,
                                      rndVal, ncolumns);
        sample->stepsLeft--;
#ifdef DEBUG
        printf("P%d.main: new row: %d\n", pid, sample->curRow);
        fflush(stdout);
#endif
    }

    if ( sample->stepsLeft == 0 && OWNED_ROW(sample->startRow) )
    {

        if (data->cache_row != sample->startRow) {
            flush_thread_cache(threadId, data);
            
            data->cache_row = sample->startRow;
            data->cache_value = 0;
        }
        data->cache_value += sample->sum;
        
        //#pragma omp atomic
        //count[sample->startRow - taskrowmin] += sample->sum;
        //#pragma omp atomic
        data->thread_messages_done++;
    }
    else {
        if ( sample->stepsLeft == 0)
            sample->curRow = sample->startRow;
#ifdef DEBUG
        printf("P%d.main: sending sample to %d: ", pid, block_owner(sample->curRow));
        printSample(sample);
        fflush(stdout);
#endif
        send_sample(out, sample, block_owner(sample->curRow));
    }
}

void prepare_idle(out_struct* out) {
    int owner, thread, i;
    long long token_value;
    // Send incomplete buffers to global buffer
    #pragma omp critical (idle)
    for(owner = 0; owner< ntasks; owner++){
        if(out->outFilled[owner] > 0){
            int outbufferIndex = owner*2+out->outBit[owner];
            for(i=0; i < out->outFilled[owner]; i++) {
                send_sample(&global_out, &(out->outArray[outbufferIndex][i]), owner);
            }
            out->outFilled[owner] = 0;
        }
    }

    #pragma omp atomic
    nthreads_idle++;

    //FIXME repeat if
    #pragma omp critical (token)
    {if(token_holder ){
        if(pid != 0)
            token_value = doneMsg;
        else
            token_value = 0;
        for(thread = 1; thread < nthreads; thread++) {
            token_value += *(doneMsgCounter[thread]);
        }
        #ifdef DEBUG
            printf("P%d.main: sending token: %d\n", pid, token_value); fflush(stdout);
        #endif
		MPI_Send(&token_value, 1, MPI_LONG_LONG_INT, (pid+1) % ntasks, TOKEN, MPI_COMM_WORLD);
		token_holder = 0;
	}}
}

//============================================================
//===================   MESSENGER    =========================
//============================================================

void initialize_out_struct(out_struct* out) {
    int i;
    
    out->outFilled = calloc(sizeof(int), ntasks); 
    out->outBit = calloc(sizeof(char), ntasks); // Select the buffer to use
    out->outWaiting = calloc(sizeof(char), ntasks); // Select the buffer to use
    out->outRequests = calloc(sizeof(MPI_Request), ntasks);
    out->outArray = calloc(sizeof(playInfo *), ntasks*2); // Double out buffers
    for(i = 0; i < ntasks*2; i++) {
        out->outArray[i] = calloc(sizeof(playInfo), OUTBUFFERSIZE);
    }
}

void process_sample_batch(int threadId, lp_thread_info data, playInfo *sample, int count)
{
	int nSamples;
for(nSamples = count; nSamples > 0; nSamples--, sample++){
        compute_sample(threadId, data, sample);
    }
}

void finish_computation()
{
    if(pid != ntasks - 1)
		MPI_Send(&done, 1, MPI_CHAR, pid+1, ALLDONE,
		 MPI_COMM_WORLD);
    done = 1;
    //omp_unset_lock(&fifo_empty);   // release other thread
}

void handle_token(playInfo *sample)
{
    doneMsg = *((long long*)sample);
    #ifdef DEBUG
    printf("P%d.helper: received token, value=%d\n",pid, doneMsg);
    fflush(stdout);
    #endif
    if(pid == 0 && doneMsg == totalDones ){
		finish_computation();
    } else {
        #pragma omp critical (token)
		token_holder = 1;
	}
}

void handle_recv_msg(int threadId, lp_thread_info data, MPI_Status *status, playInfo *sample)
{
	int count;
#ifdef DEBUG
    printf("P%d.helper: received msg", pid);
    fflush(stdout);
    printf("P%d.helper: received msg from pid=%d, tag=%d\n",
	   pid, status->MPI_SOURCE, status->MPI_TAG);
    fflush(stdout);
#endif
	MPI_Get_count(status, MPI_BYTE, &count);
    switch(status->MPI_TAG){
      case SAMPLE:
	process_sample_batch(threadId, data, sample, count/sizeof(playInfo));
	break;
      case TOKEN:
	handle_token(sample);
	break;
      case ALLDONE:
	finish_computation();
	break;
      default:
    	sprintf(msg, "Unknown tag for sample, %d!\n", status->MPI_TAG);
    	cleanExit(msg, -3);
    }
}

void communication_handler(int threadId) {
    count[threadId] = calloc(ntaskrows, sizeof(fltype)); //For now..
    doneMsgCounter[threadId] = calloc(sizeof(long long), 1);

    initialize_out_struct(&global_out);

    #pragma omp barrier
    ;

    int thread;
    while(!done) {
        for (thread = 1; thread < nthreads; thread ++) {
            if( thread_in[thread]->waiting ) {
                int flag = FALSE;
                #ifdef DEBUG
                printf("P%d.%d.comm: test %d %d\n", pid, thread, thread_in[thread]->to_process, thread_in[thread]->to_fill); fflush(stdout);
                #endif
                MPI_Test(&(thread_in[thread]->req_sample), &flag, &thread_in[thread]->status[thread_in[thread]->to_fill]);
                if(flag) {
                    thread_in[thread]->waiting = FALSE;
                    thread_in[thread]->to_fill = ((thread_in[thread]->to_fill) +1) % INQUEUESIZE;
                }
            }

            if( (! (thread_in[thread]->waiting)) && ((thread_in[thread]->to_fill+1) % INQUEUESIZE != thread_in[thread]->to_process) ) {
                #ifdef DEBUG
                printf("P%d.%d.comm: irecv %d %d\n", pid, thread, thread_in[thread]->to_process, thread_in[thread]->to_fill); fflush(stdout);
                #endif
                MPI_Irecv(thread_in[thread]->in_sample + (thread_in[thread]->to_fill*OUTBUFFERSIZE), OUTBUFFERSIZE*sizeof(playInfo), MPI_BYTE, MPI_ANY_SOURCE,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &thread_in[thread]->req_sample);
                thread_in[thread]->waiting = TRUE;
            }
        }

        if(nthreads_idle == nthreads-1) { // flush global buffers
            int owner;
            #pragma omp critical (idle)
            for(owner = 0; owner< ntasks; owner++){
                if(global_out.outFilled[owner] > 0){
                    send_to(&global_out, owner);
                }
            }
        }

    }
    for (thread = 1; thread < nthreads; thread ++) {
        if( thread_in[thread]->waiting )
            MPI_Cancel(&(thread_in[thread]->req_sample));
        ;
    }
    #pragma omp barrier
    ;
}

void sample_handler(int threadId)
{
    thread_info data;
    thread_public_info in_data;
    out_struct* out = &(data.out);
    unsigned int rand = RAND_SEED + threadId;

    playInfo * in_sample = malloc(sizeof(playInfo) * OUTBUFFERSIZE*INQUEUESIZE);
    MPI_Status *status = malloc(sizeof(MPI_Status) * INQUEUESIZE);
    count[threadId] = calloc(ntaskrows, sizeof(fltype));

    //MPI_Request req_sample;
    playInfo sample;
    int i;

    // Thread data initialization
    thread_in[threadId] = &in_data;
    in_data.to_fill = 0;
    in_data.waiting = FALSE;
    in_data.to_process = 0;
    in_data.in_sample = in_sample;
    in_data.status = status;


    data.thread_gen_done = FALSE;
    data.rand_val = RAND_SEED;
    data.gen_play = nplays+2;
    data.gen_r = taskrowmin - nthreads+1 + threadId-1;
    data.cache_row = -1;
    data.cache_value = 0.0;
    data.thread_messages_done = 0;
    doneMsgCounter[threadId] = &data.thread_messages_done;

    initialize_out_struct(out);

    #pragma omp barrier
    ;

    while(!done){
        // Process arriving samples
        while(in_data.to_process != in_data.to_fill ) {
            handle_recv_msg(threadId, &data, &status[in_data.to_process], &in_sample[in_data.to_process*OUTBUFFERSIZE]);
            in_data.to_process = (in_data.to_process+1) % INQUEUESIZE;
        }

        // Process new generated samples
        while( (in_data.to_process == in_data.to_fill) && !data.thread_gen_done) {
            for(i = 0; i < OUTBUFFERSIZE; i ++) {
                if(NULL != get_next_sample(&data, &sample)) // if there is new sample
                    compute_sample(threadId, &data, &sample);
                else
                    break;
            }
        }

        if ( (in_data.to_process == in_data.to_fill) ) {
            // Wait for activity
            #ifdef DEBUG
            printf("P%d.main: prepare idle\n", pid); fflush(stdout);
            #endif
            prepare_idle(out);
            #ifdef DEBUG
            printf("P%d.main: active waiting\n", pid); fflush(stdout);
            #endif
            while (!done && (in_data.to_process == in_data.to_fill) ) {
                ;
            }
            #pragma omp atomic
            nthreads_idle--;
        }
    }
    #ifdef DEBUG
        printf("P%d.helper: 3. flushing\n", pid);fflush(stdout);
    #endif
    flush_thread_cache(threadId, &data);
    #pragma omp barrier
    ;
#ifdef DEBUG
    printf("P%d.helper: 3. done\n", pid);fflush(stdout);
#endif
}

//======================================================================
//======================================================================
//======================================================================
void print_output(int argc, char *argv[], fltype* res, double time_diff)
{
    int r;
    char outfn[1024];
    FILE * fp;
    sprintf(outfn, "%s.%d.%d.%d_%d_%d.out", argv[1], pid, ntasks, nthreads, nsteps, nplays);
    fp = fopen(outfn, "w");
    for(r = 0; r < ntaskrows; r++){
    /*if(BLOCK_OWNER(r, ntasks, nrows) != pid)
        continue;*/
    fprintf(fp, "%.3f\n", res[r]);
    }
    fflush(stdout);
    fclose(fp);

    if(pid == 0) {
        sprintf(outfn, "%s.%d.%d.time", argv[1], ntasks, nthreads);
        fp = fopen(outfn, "w");
        fprintf(fp, "%f\n", time_diff);
        fclose(fp);
    }
}

int main(int argc, char *argv[])
{
    int r, level;
    FILE *fpc;
    double time_start, time_end;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &level);
    if(level != MPI_THREAD_MULTIPLE){
	sprintf(msg, "Requested level %d, got %d.\n", MPI_THREAD_MULTIPLE, level);
	cleanExit(msg, -8);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    // check if configuration file exists
    fpc = fopen(CONF_FILE, "r");
    if(fpc != NULL)
    	parse_conf_file(fpc);

#ifdef DEBUG
    printf("P%d: 1, after parsing standard conf file\n", pid);
    fflush(stdout);
#endif

    read_matrix(argc, argv);
    totalDones = ((long long)nrows) * ((long long)nplays); //* nsteps

#ifdef DEBUG
    printf("P%d: 2, after reading matrix\n", pid);
    fflush(stdout);
#endif

    /*int send_buffer_count = 1024*1024*1024;//100*OUTBUFFERSIZE*nrows*sizeof(playInfo) ;
    void *buffer = malloc(send_buffer_count);
    if(NULL == buffer) {
        sprintf(msg, "Failed to alloc msg buffer of %d bytes.\n",send_buffer_count);
        cleanExit(msg, -8);
    }
    MPI_Buffer_attach( buffer, send_buffer_count); //Acepts an !!int!!´*/

    MPI_Barrier(MPI_COMM_WORLD);
    
    if (pid == 0)
        time_start = omp_get_wtime();

    V = normalize_accumulate();

    srandom(RAND_SEED);

    done = 0;                         // global flag
    token_holder = (pid == 0);        // pid 0 holds initially the token

    // Prepare generation: move to right row
    gen_row = taskrowmin;

    #pragma omp parallel
    {
        #pragma omp single
        {
            nthreads = omp_get_num_threads();
            doneMsgCounter = calloc(nthreads, sizeof(long long*));
            // count starts at all 0s
            count = countMatMalloc(nthreads);
            thread_in = malloc(sizeof(thread_public_info *) * nthreads);
        }
        //printf("%d %d threadnum = %d\n", nplays, nsteps, omp_get_thread_num());
        if (omp_get_thread_num() == 0)
            communication_handler(0);
        else
            sample_handler(omp_get_thread_num());
        #pragma omp barrier
        ;
        #pragma omp for private(r)
        for(r = 0; r < ntaskrows; r++) {
            int t;
            for(t=1; t < nthreads; t++) {
                count[0][r] += count[t][r];
            }
        }
    }

    if (pid == 0)
        time_end = omp_get_wtime();

    //printf("p%d was full %lld times\n", pid, log_full_rcv_buffer); fflush(stdout);

#ifdef DEBUG
	printf("P%d: 1. calculating results!\n", pid);
	fflush(stdout);
#endif

    //Gather results from all processes
    //fltype *res = calloc(nrows, sizeof(fltype));
    //MPI_Allreduce(count, res, nrows, TYPE_FLOAT, MPI_SUM, MPI_COMM_WORLD); // Can be optimized

    for(r = 0; r < ntaskrows; r++){
        /*if(BLOCK_OWNER(r + taskrowmin, ntasks, nrows) != pid)
            continue;*/
        count[0][r] /= (nplays);
    }

#ifdef DEBUG
    printf("P%d: 2. printing results!\n", pid);
    fflush(stdout);
#endif
    print_output(argc, argv, count[0], time_end - time_start);

#ifdef DEBUG
    printf("P%d: 3. results printed and process going to finish\n", pid);
    fflush(stdout);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
	return 0;
}
