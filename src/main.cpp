#define CONSTANT 1

#include "fwMultivariate.h"
#include "optim.h"
#include "numerical.h"
//#include "BLAS/cblas_f77.h"
//#include "BLAS/cblas.h"

#include <time.h>


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
    
    u32 num = 2;

    FVar mu    = FVMake(num, 5.0, 0);
    FVar sigma = FVMake(num, 2.0, 1);

    FVar l = logNormalPDF(mu, sigma);
    FVPrint(l, "l");

    FVar params[2];
    params[0] = mu;
    params[1] = sigma;

    l = target( params );
    FVPrint(l, "l - 2");


    FVar df1 = FDiff( target, params, num, 1E-10 );
    FVPrint(df1, "df1");

    FVar df2 = CDiff( target, params, num, 1E-10 );
    FVPrint(df2, "df2");
    
    
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

    return 0;
}
