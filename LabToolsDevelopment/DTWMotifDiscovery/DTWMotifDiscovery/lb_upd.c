#include "mex.h"
#include "matrix.h"


double LB_Keogh_ea(double* restrict x, double* restrict U, double* restrict L, long long seqlen, double bsf){
    double dist = 0.0;
    for(long long i = 0; i < seqlen; i++){
        if(x[i] > U[i]){
            dist += (x[i] - U[i]) * (x[i] - U[i]);
        }
        else if(x[i] < L[i]){
            dist += (x[i] - L[i]) * (x[i] - L[i]);
        }
        if(dist >= bsf){
            return dist;
        }
    }
    return dist;
}


double LB_Keogh(double* restrict x, double* restrict U, double* restrict L, long long seqlen){
    
    double dist = 0.0;
    for(long long i = 0; i < seqlen; i++){
        if(x[i] > U[i]){
            dist += (x[i] - U[i]) * (x[i] - U[i]);
        }
        else if(x[i] < L[i]){
            dist += (x[i] - L[i]) * (x[i] - L[i]);
        }
    }
    return dist;
}


bool is_scalar(const mxArray* p){
    return (mxGetM(p) == 1) && (mxGetN(p) == 1);
}


bool is_real_doub_vec(const mxArray* p){
    return (mxGetM(p) == 1 ||  mxGetN(p) == 1) && mxIsDouble(p) && !mxIsComplex(p) && !mxIsSparse(p);
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
    
    if(nrhs < 3 || nrhs > 4){
        mexErrMsgIdAndTxt( "MATLAB:invalidNumInputs", "(required signature is (seq_x, Uy, Ly, best_so_far).");
    }
    else if (nlhs > 1){
        mexErrMsgIdAndTxt( "MATLAB:invalidNumOutputs", "More than 1 output argument requested by the calling function.");
    }
    else if(!is_real_doub_vec(prhs[0])){
        mexErrMsgIdAndTxt( "MATLAB:invalidArguments", "wrong data type, first argument");
    }
    else if(!is_real_doub_vec(prhs[1])){
        mexErrMsgIdAndTxt( "MATLAB:invalidArguments", "wrong data type, second argument");
    }
    else if(!is_real_doub_vec(prhs[2])){
        mexErrMsgIdAndTxt( "MATLAB:invalidArguments", "wrong data type, second argument");
    }
    
    else if(!matching_len(prhs[0], prhs[1]) || !matching_len(prhs[1], prhs[2])){
        mexErrMsgIdAndTxt("MATLAB:invalidArguments", "mismatched input dimensions");
    }
    long long seqlen = mxGetNumberOfElements(prhs[0]);
    double dist;
    if(nrhs == 4){
        dist = LB_Keogh_ea(mxGetDoubles(prhs[0]), mxGetDoubles(prhs[1]), mxGetDoubles(prhs[2]), seqlen, mxGetScalar(prhs[3]));
    }
    else{
        dist = LB_Keogh(mxGetDoubles(prhs[0]), mxGetDoubles(prhs[1]), mxGetDoubles(prhs[2]), seqlen);
    }
    plhs[0] = mxCreateDoubleScalar(dist);
    return;
}

