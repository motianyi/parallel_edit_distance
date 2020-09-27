// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/ 
// with many fixes and changes for multiple sequence alignment and to include an MPI driver
#include <mpi.h>
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include "sha512.hh"


using namespace std;

std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap, int *penalties);
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans);
void do_MPI_task(int rank);

/*
Examples of sha512 which returns a std::string
sw::sha512::calculate("SHA512 of std::string") // hash of a string, or
sw::sha512::file(path) // hash of a file specified by its path, or
sw::sha512::calculate(&data, sizeof(data)) // hash of any block of data
*/

// Return current time, for performance measurement
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;

// Driver code
int main(int argc, char **argv){
	int rank;
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
	MPI_Comm_rank(comm, &rank);
	if(rank==root){
		int misMatchPenalty;
		int gapPenalty;
		int k;
		std::cin >> misMatchPenalty;
		std::cin >> gapPenalty;
		std::cin >> k;	
		std::string genes[k];
		for(int i=0;i<k;i++) std::cin >> genes[i];

		int numPairs= k*(k-1)/2;

		int penalties[numPairs];
		
		uint64_t start = GetTimeStamp ();

		// return all the penalties and the hash of all allignments
		std::string alignmentHash = getMinimumPenalties(genes,
			k,misMatchPenalty, gapPenalty,
			penalties);
		
		// print the time taken to do the computation
		printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
		
		// print the alginment hash
		std::cout<<alignmentHash<<std::endl;

		for(int i=0;i<numPairs;i++){
			std::cout<<penalties[i] << " ";
		}
		std::cout << std::endl;
	} else {
		// do stuff for MPI tasks that are not rank==root
		do_MPI_task(rank);
	}
	MPI_Finalize();
	return 0;
}

/******************************************************************************/
/* Do not change any lines above here.            */
/* All of your changes should be below this line. */
/******************************************************************************/
#include <queue>
#include "omp.h"
#include "math.h"

struct Problem { 
    int i; 
    int j;
	int probNum; 
	int size;
    bool operator<(const Problem& anotherProblem) const
    {
		// order by the size (Descending order)
        return size < anotherProblem.size;
    }
}; 

struct Result { 
	int probNum; 
	int penality;
	char problemHash[129];
    bool operator<(const Result& anotherResult) const
    {
		// order by the probNum (Ascending order)
        return probNum > anotherResult.probNum;
    }
}; 

void compute(string* genes, int i, int j, int pxy, int pgap, Result* result);


int min3(int a, int b, int c) {
	if (a <= b && a <= c) {
		return a;
	} else if (b <= a && b <= c) {
		return b;
	} else {
		return c;
	}
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
int **new2d (int width, int height)
{
	int **dp = new int *[width];
	size_t size = width;
	size *= height;
	int *dp0 = new int [size];
	if (!dp || !dp0)
	{
	    std::cerr << "getMinimumPenalty: new failed" << std::endl;
	    exit(1);
	}
	dp[0] = dp0;
	for (int i = 1; i < width; i++)
	    dp[i] = dp[i-1] + height;

	return dp;
}

// called by the root MPI task only
// this procedure should distribute work to other MPI tasks
// and put together results, etc.
std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap,
	int *penalties)
{
	// MPI_Comm shmcomm;
	// MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
    //                 MPI_INFO_NULL, &shmcomm);
	// int shmrank;
	// MPI_Comm_rank(shmcomm, &shmrank);
//	std::cout << "shmrank: "<< shmrank <<" \n";

	
	// int rank;
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	std::cout << "rank: "<< rank <<" \n";
	
	// send pxy and pgap
	MPI_Bcast(&pxy, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&pgap, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// send all strings to all process
	MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
	for(int i = 0; i < k; i++){
		int length = (genes[i]).length();
		MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
		char *buffer = new char[length + 1];
		memcpy(buffer, genes[i].c_str(), length);
		buffer[length] = 0;
		MPI_Bcast(buffer, length + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
	}

	// Problem priority
	std::priority_queue<Problem> problem_queue;
	int probNum=0;
	for(int i=1;i<k;i++){
		for(int j=0;j<i;j++){
			int size = (genes[i].length()/100) * (genes[j].length()/100);
			Problem p = {i, j ,probNum, size};
			problem_queue.push(p);
			probNum ++;
		}
	}

	// Result priority queue
	std::priority_queue<Result> result_queue;


	// while(probNum > 0){
	// 	Problem p = problem_queue.top();
	// 	std::cout << "i: "<< p.i <<"  j: "<< p.j <<"  ProbNum: "<< p.probNum<<"  Problem Size: "<< p.size  <<" \n" ;
	// 	problem_queue.pop();
	// 	probNum --;
	// }

	MPI_Datatype MPI_PROBLEM;
	MPI_Type_contiguous(4, MPI_INT, &MPI_PROBLEM);
	MPI_Type_commit(&MPI_PROBLEM);

	// MPI_Datatype MPI_HASH;
	// MPI_Type_contiguous(128+1, MPI_CHAR, &MPI_HASH);

	const int elements = 3;
	MPI_Aint offsets[elements] = {0, 4, 8};
	int blocklengths[elements] = {1, 1, 128 + 1};
	MPI_Datatype types[elements] = {MPI_INT, MPI_INT, MPI_CHAR};
	MPI_Datatype MPI_RESULT;
	MPI_Type_create_struct(elements, blocklengths, offsets, types, &MPI_RESULT);
	MPI_Type_commit(&MPI_RESULT);


	// while (!problem_queue.empty()){
	// 	Problem p = problem_queue.top();
	// 	std::cout << "i: "<< p.i <<"  j: "<< p.j <<"  ProbNum: "<< p.probNum<<"  Problem Size: "<< p.size  <<" \n" ;
	// 	problem_queue.pop();
	// 	MPI_Send(&p, 1, MPI_PROBLEM, 0, 1111, MPI_COMM_WORLD);
	// }
	 
	//send initial job to each process (EXCEPT root)
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int workers = min(world_size - 1, (k*(k-1))/2);
	for(int i = 0; i < workers; i++){
		Problem p = problem_queue.top();
		problem_queue.pop();
		// (i + 1) is the destination, tag is 0
		MPI_Send(&p, 1, MPI_PROBLEM, i + 1, 0, MPI_COMM_WORLD);
	}

	while(!problem_queue.empty()){
		Result result;
		MPI_Status status;
		// send to rank 0
		MPI_Recv(&result, 1, MPI_RESULT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		result_queue.push(result);
//		std::cout<<"RECEIVE result\n";

		// send new task to it
		Problem p = problem_queue.top();
		problem_queue.pop();
		MPI_Send(&p, 1, MPI_PROBLEM, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
//		std::cout<<"Send Another\n";
	}

	// receive rest of the results
	for(int i = 0; i < workers; i++){
		Result result;
		MPI_Status status;
		MPI_Recv(&result, 1, MPI_RESULT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		result_queue.push(result);
//		std::cout<<"RECEIVE result\n";
		Problem end_exec = {-1, -1, -1, -1};
		MPI_Send(&end_exec, 1, MPI_PROBLEM, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	}

	// combine and calculate the result
	std::string alignmentHash="";
	while(!result_queue.empty()){
		Result result = result_queue.top();
		result_queue.pop();
		// std::cout <<"Collected at Rank0: ProbNum "<< result.probNum << ", Penality"<< result.penality<< ",result hash: "<< result.problemHash <<" \n" ;
		alignmentHash=sw::sha512::calculate(alignmentHash.append(result.problemHash));
		penalties[result.probNum] = result.penality;
	}	
	return alignmentHash;
}

// called for all tasks with rank!=root
// do stuff for each MPI task based on rank
void do_MPI_task(int rank)
{	
	// MPI_Comm shmcomm;
	// MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
    //                 MPI_INFO_NULL, &shmcomm);
	// int shmrank;
	// MPI_Comm_rank(shmcomm, &shmrank);
//	std::cout << "shmrank: "<< shmrank <<" \n";

	int root = 0;
	// int rank2;
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank2);
//	std::cout << "rank: "<< rank2 <<" \n";
	int pxy, pgap;
	MPI_Bcast(&pxy, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&pgap, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//receive number of strings
	int k;
	MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);

	std::string genes[k];

	//receive strings
	for(int i = 0; i < k; i++){
		int length;
		MPI_Bcast(&length, 1, MPI_INT, root, MPI_COMM_WORLD);
		char *buffer = new char[length + 1];
		buffer[length] = 0;
		MPI_Bcast(buffer, length + 1, MPI_CHAR, root, MPI_COMM_WORLD);
		genes[i] = buffer;
	}

	// MPI Datatypes
	MPI_Datatype MPI_PROBLEM;
	MPI_Type_contiguous(4, MPI_INT, &MPI_PROBLEM);
	MPI_Type_commit(&MPI_PROBLEM);
	const int elements = 3;
	MPI_Aint offsets[elements] = {0, 4, 8};
	int blocklengths[elements] = {1, 1, 128 + 1};
	MPI_Datatype types[elements] = {MPI_INT, MPI_INT, MPI_CHAR};
	MPI_Datatype MPI_RESULT;
	MPI_Type_create_struct(elements, blocklengths, offsets, types, &MPI_RESULT);
	MPI_Type_commit(&MPI_RESULT);

	while(true){
		// receive task
		Problem p;
		MPI_Status status;
		MPI_Recv(&p, 1, MPI_PROBLEM, root, 0, MPI_COMM_WORLD, &status);
		if(p.probNum < 0){
			break;
		}
//		std::cout <<"Rank"<< rank<< ", i: "<< p.i <<"  j: "<< p.j <<"  ProbNum: "<< p.probNum<<"  Problem Size: "<< p.size  <<" \n" ;

		//compute
		Result result;
		result.probNum = p.probNum;
		compute(genes, p.i, p.j, pxy, pgap, &result);
//		std::cout <<"ProbNum "<< result.probNum << ", Penality"<< result.penality<< ",result hash: "<< result.problemHash <<" \n" ;

		// send result to rank 0
		MPI_Send(&result, 1, MPI_RESULT, root, 0, MPI_COMM_WORLD);
		// std::cout <<"Sent\n ";
	}
}



void compute(string* genes, int i, int j, int pxy, int pgap, Result* result)
{
	std::string gene1 = genes[i];
	std::string gene2 = genes[j];

	int m = gene1.length(); // length of gene1
	int n = gene2.length(); // length of gene2
	int l = m+n;
	// penalties[probNum]=
	int xans[l+1], yans[l+1];

	// Assign result
	result->penality = getMinimumPenalty(gene1,gene2,pxy,pgap,xans,yans);
	// Since we have assumed the answer to be n+m long,
	// we need to remove the extra gaps in the starting
	// id represents the index from which the arrays
	// xans, yans are useful
	int id = 1;
	int a;
	for (a = l; a >= 1; a--)
	{
		if ((char)yans[a] == '_' && (char)xans[a] == '_')
		{
			id = a + 1;
			break;
		}
	}
	std::string align1="";
	std::string align2="";
	for (a = id; a <= l; a++)
	{
		align1.append(1,(char)xans[a]);
	}
	for (a = id; a <= l; a++)
	{
		align2.append(1,(char)yans[a]);
	}
	std::string align1hash = sw::sha512::calculate(align1);
	std::string align2hash = sw::sha512::calculate(align2);
	std::string problemhash = sw::sha512::calculate(align1hash.append(align2hash));

	// Assign result
	strcpy(result -> problemHash, problemhash.c_str());
}


// // function to find out the minimum penalty
// // return the minimum penalty and put the aligned sequences in xans and yans
// int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans)
// {
	
// 	int i, j; // intialising variables

// 	int m = x.length(); // length of gene1
// 	int n = y.length(); // length of gene2
// 	std::cout <<"m: "<< m<<", "<<"n: "<< n<< "\n";
	
// 	// table for storing optimal substructure answers
// 	int **dp = new2d (m+1, n+1);
// 	size_t size = m + 1;
// 	size *= n + 1;
// 	memset (dp[0], 0, size);

// 	// intialising the table
// 	for (i = 0; i <= m; i++)
// 	{
// 		dp[i][0] = i * pgap;
// 	}
// 	for (i = 0; i <= n; i++)
// 	{
// 		dp[0][i] = i * pgap;
// 	}

// 	// calcuting the minimum penalty
// 	for (i = 1; i <= m; i++)
// 	{
// 		for (j = 1; j <= n; j++)
// 		{
// 			if (x[i - 1] == y[j - 1])
// 			{
// 				dp[i][j] = dp[i - 1][j - 1];
// 			}
// 			else
// 			{
// 				dp[i][j] = min3(dp[i - 1][j - 1] + pxy ,
// 						dp[i - 1][j] + pgap ,
// 						dp[i][j - 1] + pgap);
// 			}
// 		}
// 	}

// 	// Reconstructing the solution
// 	int l = n + m; // maximum possible length
	
// 	i = m; j = n;
	
// 	int xpos = l;
// 	int ypos = l;
	
// 	while ( !(i == 0 || j == 0))
// 	{
// 		if (x[i - 1] == y[j - 1])
// 		{
// 			xans[xpos--] = (int)x[i - 1];
// 			yans[ypos--] = (int)y[j - 1];
// 			i--; j--;
// 		}
// 		else if (dp[i - 1][j - 1] + pxy == dp[i][j])
// 		{
// 			xans[xpos--] = (int)x[i - 1];
// 			yans[ypos--] = (int)y[j - 1];
// 			i--; j--;
// 		}
// 		else if (dp[i - 1][j] + pgap == dp[i][j])
// 		{
// 			xans[xpos--] = (int)x[i - 1];
// 			yans[ypos--] = (int)'_';
// 			i--;
// 		}
// 		else if (dp[i][j - 1] + pgap == dp[i][j])
// 		{
// 			xans[xpos--] = (int)'_';
// 			yans[ypos--] = (int)y[j - 1];
// 			j--;
// 		}
// 	}
// 	while (xpos > 0)
// 	{
// 		if (i > 0) xans[xpos--] = (int)x[--i];
// 		else xans[xpos--] = (int)'_';
// 	}
// 	while (ypos > 0)
// 	{
// 		if (j > 0) yans[ypos--] = (int)y[--j];
// 		else yans[ypos--] = (int)'_';
// 	}

// 	int ret = dp[m][n];

// 	delete[] dp[0];
// 	delete[] dp;
	
// 	return ret;
// }


// function to find out the minimum penalty
// return the maximum penalty and put the aligned sequences in xans and yans
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap,
	int* xans, int* yans)
{
	int i, j; // intialising variables
	
	int m = x.length(); // length of gene1
	int n = y.length(); // length of gene2
	
	// table for storing optimal substructure answers
	int **dp = new2d (m+1, n+1);
//	size_t size = m + 1;
//	size *= n + 1;
//	memset (dp[0], 0, size);

	// intialising the table
	#pragma omp parallel for
	for (i = 0; i <= m; i++)
	{
		dp[i][0] = i * pgap;
	}
	#pragma omp parallel for
	for (i = 0; i <= n; i++)
	{
		dp[0][i] = i * pgap;
	}

	int number_threads = 64;
    omp_set_num_threads(16);
	//integer devision
	//calculate the deminsion of each block
	int block_width = m / number_threads;
	int block_height = n / number_threads;
	int width_remainder = m % number_threads;
	int height_remainder = n % number_threads;

	for(int k = 0; k < (2 * number_threads - 1); k++ ){

		int column_min = max(0, k - number_threads + 1);
        int column_max = min(k, number_threads - 1);

		#pragma omp parallel for
		for (int column = column_min; column <= column_max; column++){
			int row = k - column;

			//calculate the elements index in the block
			int block_min_column = column * block_height + 1;
			int block_max_column = (column + 1) * block_height;
			int block_min_row = row * block_width + 1;
			int block_max_row = (row + 1) * block_width;
			
			// handle row and column remainder
			if(column == number_threads-1){
				block_max_column += height_remainder;
			}
			if(row == number_threads-1){
				block_max_row += width_remainder;
			}
			
			// sequential calculation within each block
			if(block_max_column >= block_min_column && block_max_row >= block_min_row){
				
				for (int i = block_min_row; i <= block_max_row; i++){
					for (int j = block_min_column; j <= block_max_column; j++){
						// printf("i = %d\n", i);
						if (x[i - 1] == y[j - 1]){
							dp[i][j] = dp[i - 1][j - 1];
						}else{
							dp[i][j] = min(min(dp[i - 1][j - 1] + pxy ,
									dp[i - 1][j] + pgap) ,
									dp[i][j - 1] + pgap);
						}
					}
				}
				
			}
		}
	}



	// Reconstructing the solution
	int l = n + m; // maximum possible length
	
	i = m; j = n;
	
	int xpos = l;
	int ypos = l;
	
	while ( !(i == 0 || j == 0))
	{
		if (x[i - 1] == y[j - 1])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j - 1] + pxy == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)'_';
			i--;
		}
		else if (dp[i][j - 1] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)'_';
			yans[ypos--] = (int)y[j - 1];
			j--;
		}
	}
	omp_set_num_threads(omp_get_max_threads());
    int x_diff = xpos - i, y_diff = ypos-j;
#pragma omp parallel for
    for (int ii = i; ii>0; --ii){
        xans[ii+x_diff] = (int)x[ii-1];
    }
#pragma omp parallel for
    for (int x_pos2=xpos-i; x_pos2>0; --x_pos2){
        xans[x_pos2] = (int)'_';
    }

#pragma omp parallel for
    for (int jj = j; jj>0; --jj){
        yans[jj+y_diff] = (int)y[jj-1];
        if (jj==0){
        }
    }

#pragma omp parallel for
    for (int y_pos2=ypos-j; y_pos2>0; --y_pos2){
        yans[y_pos2] = (int)'_';
    }

	int ret = dp[m][n];

	delete[] dp[0];
	delete[] dp;
	return ret;
}
// mpicxx -std=c++14 -fopenmp -o tianyim-seqalkway tianyim-seqalkway.cpp -O3
