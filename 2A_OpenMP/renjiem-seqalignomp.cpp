// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/ and
// fixed an error when initializing the dp array :-)
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>

using namespace std;

int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans);

// Return current time, for performance measurement
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t) 1000000 + tv.tv_usec;
}


// Driver code
int main() {
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
    int l = m + n;
    int xans[l + 1], yans[l + 1];

    uint64_t start = GetTimeStamp();

    // calling the function to calculate the result
    int penalty = getMinimumPenalty(gene1, gene2,
                                    misMatchPenalty, gapPenalty,
                                    xans, yans);

    // print the time taken to do the computation
    printf("Time: %ld us\n", (uint64_t)(GetTimeStamp() - start));

    // postprocessing of the answer, for printing results

    // Since we have assumed the answer to be n+m long,
    // we need to remove the extra gaps in the starting
    // id represents the index from which the arrays
    // xans, yans are useful
    int id = 1;
    int i;
    for (i = l; i >= 1; i--) {
        if ((char) yans[i] == '_' && (char) xans[i] == '_') {
            id = i + 1;
            break;
        }
    }

    // Printing the final answer
    std::cout << "Minimum Penalty in aligning the genes = ";
    std::cout << penalty << std::endl;
    std::cout << "The aligned genes are :" << std::endl;
    for (i = id; i <= l; i++) {
        std::cout << (char) xans[i];
    }
    std::cout << "\n";
    for (i = id; i <= l; i++) {
        std::cout << (char) yans[i];
    }
    std::cout << "\n";

    return 0;
}

/******************************************************************************/
/* Do not change any lines above here.            */
/* All of your changes should be below this line. */
/******************************************************************************/
#include "omp.h"
#include "math.h"

inline int min3(int a, int b, int c) {
    if (a <= b) {
        if (a <= c) {
            return a;
        }
        return c;
    } else if (b <= c) {
        return b;
    }
    return c;
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
inline int **new2d(int width, int height) {
    int **dp = new int *[width];
    size_t size = width;
    size *= height;
    int *dp0 = new int[size];
    if (!dp || !dp0) {
        std::cerr << "getMinimumPenalty: new failed" << std::endl;
        exit(1);
    }
    dp[0] = dp0;


    for (int i = 1; i < width; i++)
        dp[i] = dp[i - 1] + height;

    return dp;
}

// function to find out the minimum penalty
// return the maximum penalty and put the aligned sequences in xans and yans
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap,
                      int *xans, int *yans) {
    int i, j; // intialising variables

    int m = x.length(); // length of gene1
    int n = y.length(); // length of gene2

    // table for storing optimal substructure answers
    omp_set_num_threads(omp_get_max_threads());
    int **dp = new2d(m + 1, n + 1);


    // remove unnecessary memset
//    size_t size = m + 1;
//    size *= n + 1;
//    memset (dp[0], 0, size);

    // intialising the table
#pragma omp parallel for
    for (i = 0; i <= m; ++i) {
        dp[i][0] = i * pgap;
    }
#pragma omp parallel for
    for (i = 0; i <= n; ++i) {
        dp[0][i] = i * pgap;
    }

    int n_threads = 22;
    omp_set_num_threads(n_threads);
    // calculating the minimum penalty with the tiling technique in an anti-diagonal version
    int tile_row_size = (int) ceil((1.0 * m) / n_threads); // Number of dp elements in row of each tile
    int tile_col_size = (int) ceil((1.0 * n) / n_threads); // Number of dp elements in column of each tile

//    int tile_row_size = 256; // Number of dp elements in row of each tile
//    int tile_col_size = 256; // Number of dp elements in column of each tile
    int tile_m = (int) ceil((1.0 * m) / tile_row_size); // Number of tiles in row of the dp matrix
    int tile_n = (int) ceil((1.0 * n) / tile_col_size); // Number of tile in column of the dp matrix

    int total_diagonal = tile_m + tile_n - 1;
    int row_min, row_max, diagonal_index, k;
//    cout << "tile_row_size: " << tile_row_size << ", tile_col_size: " << tile_col_size << endl;
//    cout << "tile_m: " << tile_m << ", tile_n: " << tile_n << endl;
//    cout << "total_diagonal: " << total_diagonal << endl;
    for (diagonal_index = 1; diagonal_index <= total_diagonal; ++diagonal_index) {
        row_min = max(1, diagonal_index - tile_n + 1);
        row_max = min(diagonal_index, tile_m);
#pragma omp parallel for
        for (k = row_min; k <= row_max; ++k) {
            int tile_row_start = 1 + (k - 1) * tile_row_size; // index inclusive
            int tile_row_end = min(tile_row_start + tile_row_size, m + 1); // index exclusive
            int tile_col_start = 1 + (diagonal_index - k) * tile_col_size; // index inclusive
            int tile_col_end = min(tile_col_start + tile_col_size, n + 1); // index exclusive

//            cout << "(" << tile_row_start<< "," << tile_col_start << ")" << " | ";
//            cout << "-> (" << tile_row_end << "," << tile_col_end << ")" << '|';
            for (int ii = tile_row_start; ii < tile_row_end; ++ii) {
                for (int jj = tile_col_start; jj < tile_col_end; ++jj) {
                    if (x[ii - 1] == y[jj - 1]) {
                        dp[ii][jj] = dp[ii - 1][jj - 1];
                    } else {
                        dp[ii][jj] = min3(dp[ii - 1][jj - 1] + pxy,
                                          dp[ii - 1][jj] + pgap,
                                          dp[ii][jj - 1] + pgap);
                    }
                }
            }
        }
//        cout << "done" << endl;
    }


    // Reconstructing the solution
    int l = n + m; // maximum possible length

    i = m;
    j = n;

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
        else
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

//    delete[] dp[0];
//    delete[] dp;

    return ret;
}
//g++ -fopenmp -o renjiem-seqalignomp renjiem-seqalignomp.cpp -O3