#include<vector>
#include<limits>
#include<algorithm>
#include<cmath>
#include "mex.h"
#include "matrix.h"
#include <iostream>
using namespace std;

constexpr double INF{std::numeric_limits<double>::infinity()};

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band


double* compute_dtw(double* A, double* B, double *cb, int m, int r, double best_so_far) {
    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i, j, k;
    double x, y, z, min_cost;
    int ea = 0;
    double* res;
    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two arrays of size O(r).
    res = (double*)calloc(2, sizeof(double));
    cost = (double*) calloc(2 * r + 1, sizeof(double));
    cost_prev = (double*) calloc(2 * r + 1, sizeof(double));
    for (k = 0; k < 2 * r + 1; k++) {
        cost[k] = INF;
        cost_prev[k] = INF;
    }

    for (i = 0; i < m; i++) {
        k = std::max(0, r - i);
        min_cost = INF;

        for (j = std::max(0, i - r); j <= std::min(m - 1, i + r); j++, k++) {
            /// Initialize all row and column
            if ((i == 0) && (j == 0)) {
                double c = (A[0] - B[0]);
                cost[k] = c * c;
                min_cost = cost[k];
                continue;
            }

            if ((j - 1 < 0) || (k - 1 < 0)) {
                y = INF;
            } else {
                y = cost[k - 1];
            }
            if ((i < 1) || (k > 2 * r - 1)) {
                x = INF;
            } else {
                x = cost_prev[k + 1];
            }
            if ((i < 1) || (j < 1)) {
                z = INF;
            } else {
                z = cost_prev[k];
            }

            /// Classic DTW calculation
            double d = A[i] - B[j];
            cost[k] = std::min(std::min( x, y) , z) + d * d;

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost) {
                min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distance with lower bound together are larger than best_so_far
        if (double((i+r+1))/m<=0.5 && min_cost + cb[i + r + 1] >= best_so_far) {
            free(cost);
            free(cost_prev);
           // cout << "Early abondoned " << i+r+1 << endl;
           // if (((i+r+1)/m) <= 0.01)
                ea = ea + 1;//double((i+r+1))/m;
            
            res[0] = min_cost + cb[i + r + 1];
            res[1] = ea;
            return res;
        }

        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    res[0] = final_dtw;
    res[1] = ea;
    //return final_dtw;
    return res;
}

bool is_scalar(const mxArray* p){
    return (mxGetM(p) == 1) && (mxGetN(p) == 1);
}

bool is_real_doub_vec(const mxArray* p){
    return (std::min(mxGetM(p), mxGetN(p)) == 1) && mxIsDouble(p) && !mxIsComplex(p) && !mxIsSparse(p);
}
        
bool is_real_scalar(const mxArray* p){
    return (mxGetM(p) == 1) && (mxGetN(p) == 1);  
    // Matlab defaults to double. If we check for an integral type here
    // it would require passing int64(arg) from matlab, which is annoying.
}

bool is_real_scalar_doub(const mxArray* p){
    return (mxGetM(p) == 1) && (mxGetN(p) == 1) && mxIsDouble(p);
}

bool matching_len(const mxArray* p, const mxArray* q){
    return (mxGetM(p) == mxGetM(q)) && (mxGetN(p) == mxGetN(q));
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[] ){
    // This expects 2 sequences, a max warping parameter as a 64 bit integer, and 
    if(nrhs < 3 or nrhs > 4){
        mexErrMsgIdAndTxt( "MATLAB:invalidNumInputs", "(required signature is (seq_x, seq_y, warpmax, bsf).");
    }
//     else if (nlhs > 1){
//         mexErrMsgIdAndTxt( "MATLAB:invalidNumOutputs", "More than 1 output argument requested by the calling function.");
//     }
    else if(!is_real_doub_vec(prhs[0])){
        mexErrMsgIdAndTxt( "MATLAB:invalidArguments", "wrong data type, first argument");
    }
    else if(!is_real_doub_vec(prhs[1])){
        mexErrMsgIdAndTxt( "MATLAB:invalidArguments", "wrong data type, second argument");
    }
    else if(!matching_len(prhs[0], prhs[1])){
        mexErrMsgIdAndTxt("MATLAB:invalidArguments", "mismatched input dimensions");
    }
    else if(!is_real_scalar(prhs[2])){
        mexErrMsgIdAndTxt("MATLAB:invalidArguments", "invalid warping factor parameter, expected int64 type.");
    }
    
    // look up check for scalar int on third 
    if(!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0])){
        mexErrMsgIdAndTxt( "MATLAB:invalidArguments", "wrong data type");
    }
    long long warpmax = mxGetScalar(prhs[2]);
    long long seqlen = mxGetNumberOfElements(prhs[0]); 
    
    if(warpmax > seqlen){
        mexErrMsgIdAndTxt("MATLAB:invalid parameter", "max warping window size exceeds sequence length");
    }
    
    double* x = mxGetDoubles(prhs[0]);
    double* y = mxGetDoubles(prhs[1]);
    std::vector<double> cb(seqlen, 0);
    double max_threshold;
    if(nrhs == 4){
        if(!is_real_scalar_doub(prhs[3])){
            mexErrMsgIdAndTxt("MATLAB:invalidArguments", "mismatched input dimensions");
        }
        max_threshold = mxGetScalar(prhs[3]);
        max_threshold *= max_threshold;  // this uses squared distance internally
    }
    else{
        max_threshold = std::numeric_limits<double>::infinity(); //INF; 
    }
    double* res = compute_dtw(x, y, cb.data(), seqlen, warpmax, max_threshold);
    double dist = res[0];
    plhs[0] = mxCreateDoubleScalar(sqrt(dist));
    plhs[1] = mxCreateDoubleScalar(res[1]);
    return;
}

