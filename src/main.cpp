#define CONSTANT 1

#include "fwMultivariate.h"
#include "matrix.h"
#include "optim.h"
#include "numerical.h"


// TODO(jonas): check Xcode build settings with BLAS, maybe use LAPACK
// TODO(jonas): make arena for multivariate


FVar logNormalPDF(FVar mu, FVar sigma) {

    FVar fv = FVDiv( FVMul(mu, mu), FVMul(sigma, sigma) );

#if CONSTANT
    fv = FVAddD( fv, log( 2 * M_PI ) );
#endif

    fv = FVMulD( fv, -0.5 );

    fv = FVSub( fv, FVLog( sigma ) );

    return fv;
}


FVar target( FVar *input ) {
    return logNormalPDF(input[0], input[1]);
}


int main(int argc, const char * argv[]) {
    
//    u32 num = 2;
//
//    FVar mu    = FVMake(num, 5.0, 0);
//    FVar sigma = FVMake(num, 2.0, 1);
//
//    FVar l = logNormalPDF(mu, sigma);
//    FVPrint(l, "l");
//
//    FVar params[2];
//    params[0] = mu;
//    params[1] = sigma;
//
//    l = target( params );
//    FVPrint(l, "l - 2");
//
//
//    FVar df1 = FDiff( target, params, num, 1E-10 );
//    FVPrint(df1, "df1");
//
//    FVar df2 = CDiff( target, params, num, 1E-10 );
//    FVPrint(df2, "df2");
    
    
//    FVar a = FVMake(num, 5.0,  0);
//    FVar b = FVMake(num, -1.0, 1);
//
//
//    double c = 3.0;
//
//    FVPrint(a, "a");
//    FVPrint(b, "b");
//
//    FVar out = FVAdd( FVMulD(a, c), b );
//
//    FVPrint(out, "out");
    
    Mat A, B, C;
    u32 m = 10, k = 10, n = 5, i, j;
    
    A = MatMake(DefaultAllocator, m, k);
    B = MatMake(DefaultAllocator, k, n);

    printf (" Intializing matrix data \n\n");
    
    for (i=0; i<m; ++i) {
        for (j=0; j<k; ++j) {
            SetElement(A, i, j, (f64) i+j);
        }
    }
    MatPrint(A, "A");

    for (i=0; i<k; ++i) {
        for (j=0; j<n; ++j) {
            SetElement(B, i, j, -(f64)i-(f64)j);
        }
    }
    MatPrint(B, "B");
    
    C = MatMul(A, B);
    
    MatPrint(A, "A");
    MatPrint(B, "B");
    MatPrint(C, "C");
    
    return 0;
}
