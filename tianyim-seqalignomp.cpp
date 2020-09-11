// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/ and
// fixed an error when initializing the dp array :-)
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>

using namespace std;

int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int* xans, int* yans);

// Return current time, for performance measurement
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}


// Driver code
int main(){
	int misMatchPenalty;
	int gapPenalty;
	std::string gene1;
	std::string gene2;
	std::cin >> misMatchPenalty;
	std::cin >> gapPenalty;
	std::cin >> gene1;
	std::cin >> gene2;
	std::cout << "misMatchPenalty=" << misMatchPenalty << std::endl;
	std::cout << "gapPenalty=" << gapPenalty << std::endl;

	int m = gene1.length(); // length of gene1
	int n = gene2.length(); // length of gene2
	int l = m+n;
	int xans[l+1], yans[l+1];

	uint64_t start = GetTimeStamp ();

	// calling the function to calculate the result
	int penalty = getMinimumPenalty(gene1, gene2,
		misMatchPenalty, gapPenalty,
		xans,yans);
	
	// print the time taken to do the computation
	printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
	
	// postprocessing of the answer, for printing results

	// Since we have assumed the answer to be n+m long,
	// we need to remove the extra gaps in the starting
	// id represents the index from which the arrays
	// xans, yans are useful
	int id = 1;
	int i;
	for (i = l; i >= 1; i--)
	{
		if ((char)yans[i] == '_' && (char)xans[i] == '_')
		{
			id = i + 1;
			break;
		}
	}
	
	// Printing the final answer
	std::cout << "Minimum Penalty in aligning the genes = ";
	std::cout << penalty << std::endl;
	std::cout << "The aligned genes are :" << std::endl;
	for (i = id; i <= l; i++)
	{
		std::cout<<(char)xans[i];
	}
	std::cout << "\n";
	for (i = id; i <= l; i++)
	{
		std::cout << (char)yans[i];
	}
	std::cout << "\n";

	return 0;
}

/******************************************************************************/
/* Do not change any lines above here.            */
/* All of your changes should be below this line. */
/******************************************************************************/

int inline min3(int a, int b, int c) {
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
	size_t size = m + 1;
	size *= n + 1;
	memset (dp[0], 0, size);

	// intialising the table
	for (i = 0; i <= m; i++)
	{
		dp[i][0] = i * pgap;
	}
	for (i = 0; i <= n; i++)
	{
		dp[0][i] = i * pgap;
	}

	// original sequential version
	// calcuting the minimum penalty
	// for (i = 1; i <= m; i++)
	// {
	// 	for (j = 1; j <= n; j++)
	// 	{
	// 		if (x[i - 1] == y[j - 1])
	// 		{
	// 			dp[i][j] = dp[i - 1][j - 1];
	// 		}
	// 		else
	// 		{
	// 			dp[i][j] = min3(dp[i - 1][j - 1] + pxy ,
	// 					dp[i - 1][j] + pgap ,
	// 					dp[i][j - 1] + pgap);
	// 		}
	// 	}
	// }

	// first parallel version
	// int k, i_min, i_max;
	// for (k = 2; k <= m + n; k++){
    //     int i_min = max(1, k - n);
    //     int i_max   = min(k - 1, m);
	// 	// printf("k value: %d\n", k);
	// 	#pragma omp parallel for shared(dp)
    //     for (int i = i_min; i <= i_max; i++){
    //         int j = k - i;
	// 		// printf("i value: %d, j value: %d, %d,%d,%d\n", i, j,dp[i-1][j-1] + (x[j-1] == y[i-1] ? 0 : pxy),dp[i-1][j] + pgap, dp[i][j-1] + pgap);
    //         dp[i][j] = min3(dp[i-1][j-1] + (x[i-1] == y[j-1] ? 0 : pxy), dp[i-1][j] + pgap, dp[i][j-1] + pgap);
			
		
    //     }
    // } 

	int k;

	int number_threads = 24;
	//integer devision
	int block_width = m / number_threads;
	int block_height = n / number_threads;
	int width_remainder = m % number_threads;
	int height_remainder = n % number_threads;

	for(k = 0; k < (2 * number_threads - 1); k++ ){
		// printf("KKKK");
		int column_min = max(0, k - number_threads + 1);
        int column_max = min(k, number_threads - 1);
		// printf("cloumn-min-max %d, %d \n", column_min, column_max);
		#pragma omp parallel for shared(dp, x,y)
		for (int column = column_min; column <= column_max; column++){
			int row = k - column;
			int block_min_column = column * block_height + 1;
			int block_max_column = (column + 1) * block_height;
			int block_min_row = row * block_width + 1;
			int block_max_row = (row + 1) * block_width;
			// printf("%d, %d \n", row, column);
			// handle last row
			if(column == number_threads-1){
				block_max_column += height_remainder;
			}
			if(row == number_threads-1){
				block_max_row += width_remainder;
			}
			// printf("%d, %d, %d, %d \n", block_min_column, block_max_column, block_min_row, block_max_row);
			// sequential block calculation
			for (int i = block_min_column; i <= block_max_column; i++){
				// printf("i = %d\n", i);
				for (int j = block_min_row; j <= block_max_row; j++){
					if (x[i - 1] == y[j - 1]){
						dp[i][j] = dp[i - 1][j - 1];
					}else{
						dp[i][j] = min3(dp[i - 1][j - 1] + pxy ,
								dp[i - 1][j] + pgap ,
								dp[i][j - 1] + pgap);
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
	while (xpos > 0)
	{
		if (i > 0) xans[xpos--] = (int)x[--i];
		else xans[xpos--] = (int)'_';
	}
	while (ypos > 0)
	{
		if (j > 0) yans[ypos--] = (int)y[--j];
		else yans[ypos--] = (int)'_';
	}

	int ret = dp[m][n];

	delete[] dp[0];
	delete[] dp;
	
	return ret;
}

// int getMinimumPenalty2(std::string x, std::string y, int pxy, int pgap,int* xans, int* yans){
// 	int i, j; // intialising variables
// 	int m = x.length(); // length of gene1
// 	int n = y.length(); // length of gene2
	

// 	// second parallel version
// 	int **dp_v = new2d(m + n + 1, min(n, m) +1);
// 	int width = min(n, m);
// 	size_t size = (m + n + 1);
// 	size *= (min(n, m) +1);
// 	memset (dp_v[0], 0, size);
// 	dp_v[0][0] = 0;
// 	dp_v[1][0] = pgap;
// 	dp_v[1][1] = pgap;
// 	int k;
// 	for (k = 2; k <= m + n; k++){
//         int i_min = max(1, k - n);
//         int i_max = min(k - 1, m);
// 		int *last_last_row = dp_v[k-2];
// 		int *last_row = dp_v[k-1];
// 		int *this_row = dp_v[k];
//         //Fill in Diagonal with base case (if the diagonal passes through a base case cell
//         if(k <= min(m, n)){
//         	this_row[k] = k * pgap;
//         	this_row[0] = k * pgap;
//         }

						
//        	#pragma omp parallel for firstprivate(i_min, i_max, this_row, last_row, last_last_row) shared(x, y)
//         for (int i = i_min; i <= i_max; i++){
// 			int j = k - i;
//             this_row[i] = min3(last_row[i] + pgap, last_row[i - 1] + pgap, last_last_row[i - 1] + (y[j - 1] == x[i - 1] ? 0 : pxy));
// 		}
//     }
   
// 	// Reconstructing the solution
// 	int l = n + m; // maximum possible length
	
// 	i = m; j = n;
	
// 	int xpos = l;
// 	int ypos = l;
	
// 	while ( !(i == 0 || j == 0))
// 	{
// 		k = i + j;
// 		if (x[i - 1] == y[j - 1])
// 		{
// 			xans[xpos--] = (int)x[i - 1];
// 			yans[ypos--] = (int)y[j - 1];
// 			i--; j--;
// 		}
// 		else if (dp_v[k - 2][i - 1] + pxy == dp_v[k][i])
// 		{
// 			xans[xpos--] = (int)x[i - 1];
// 			yans[ypos--] = (int)y[j - 1];
// 			i--; j--;
// 		}
// 		else if (dp_v[k-1][i - 1] + pgap == dp_v[k][i])
// 		{
// 			xans[xpos--] = (int)x[i - 1];
// 			yans[ypos--] = (int)'_';
// 			i--;
// 		}
// 		else if (dp_v[k - 1][i] + pgap == dp_v[k][i])
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

// 	// int ret = dp[m][n];
// 	int ret = dp_v[m + n][min(n, m)];

// 	delete[] dp_v[0];
// 	delete[] dp_v;
	
// 	return ret;
// }

// int inline min(int x, int y){
// 	if(x > y){
// 		return y;
// 	}
// 	return x;
// }

// int inline max(int x, int y){
// 	if(x < y){
// 		return y;
// 	}
// 	return x;
// }

//g++ -fopenmp tianyim-seqalignomp.cpp -o tianyim-seqalignomp -O3