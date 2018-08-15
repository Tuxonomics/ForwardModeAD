
#define CONSTANT 1


#include "utilities.h"
#include "grad.h"
//#include "optim.h"


// TODO(jonas): make arena for multivariate


FVar logNormalPDF(FVar mu, FVar sigma)
{

    FVar fv = FVDiv( FVMul(mu, mu), FVMul(sigma, sigma) );

#if CONSTANT
    fv = FVAddD( fv, log( 2 * M_PI ) );
#endif

    fv = FVMulD( fv, -0.5 );

    fv = FVSub( fv, FVLog( sigma ) );

    return fv;
}


Inline
FVar target( FVar *input )
{
    return logNormalPDF(input[0], input[1]);
}


Inline
FVar target2( FVarMat input )
{
    return logNormalPDF( input.data[0], input.data[1] );
}



FVar bimodal( FVar x, FVar y )
{
    FVar out;
    
    f64 mu = 2.5;
    f64 s  = 1;
    
    out = FVExp( FVDMul( -0.5*s, FVAdd( FVPow( FVDAdd( mu, y), 2), FVPow( FVDAdd( -mu, x), 2) ) ) );
    
    out = FVAdd( out, FVExp( FVDMul( -0.5*s, FVAdd( FVPow( FVDAdd( -mu, y), 2), FVPow( FVDAdd( mu, x), 2) ) ) ) );
    
    return FVLog( out );
}

Inline
FVar bmTarget( FVarMat input )
{
    return bimodal( input.data[0], input.data[1] );
}



#ifndef TEST
int main(int argc, const char * argv[]) {

    FVar mu    = FVMake(5.0, 0);
    FVar sigma = FVMake(2.0, 1);

    FVar l = logNormalPDF(mu, sigma);
    FVPRINT(l);
    
    
    FVar params[2];
    params[0] = mu;
    params[1] = sigma;

    l = target( params );
    FVPrint(l, "l - 2");

    FVarMat mat = FVarMatMake( DefaultAllocator, 2, 1 );
    
    FVarMatSetElement( mat, 0, 0, mu );
    FVarMatSetElement( mat, 1, 0, sigma );
    
    FVarMatPrint( mat, "mat" );
    
    l = target2( mat );
    FVPrint( l, "l - 3" );
    
    
    FVarMat grad = FVarMatMake( DefaultAllocator, 2, 1 );
    
    FVarGradient( DefaultAllocator, &target2, mat, grad );
    
    FVarMatPrint( grad, "grad" );
    
    
    FVarFDiff( DefaultAllocator, &target2, mat, grad, 0.000001 );
    
    FVarMatPrint( grad, "FDiff" );

    
    FVarCDiff( DefaultAllocator, &target2, mat, grad, 0.000001 );
    
    FVarMatPrint( grad, "CDiff" );
    
    
//    Xorshift1024 x = Xorshift1024Init( 37473 );
//
//    printf( "new number: %llu\n", rngXorshift1024Next( &x ) );
//    printf("new number: %.8f\n", rngXorshift1024NextFloat( &x ) );
//
//
//    Xorshift1024 xx = XORSHIFT_INIT;
//    Rng rng = RngInitXorshift1024( &xx );
//
//    printf( "new number: %llu\n", RngNext( rng ) );
//    printf("new number: %.8f\n", RngNextFloat( rng ) );
//
//    printf("c(\n");
//    for ( u32 i=0; i<1000; ++i )
//        printf( "%.8f,\n", RngNormal( rng ) );
//
//    printf( "%.8f\n", RngNormal( rng ) );
//    printf(")\n");
    
////    u32 num = 2;
////
////    FVar mu    = FVMake(num, 5.0, 0);
////    FVar sigma = FVMake(num, 2.0, 1);
////
////    FVar l = logNormalPDF(mu, sigma);
////    FVPRINT(l);
////
////    FVar params[2];
////    params[0] = mu;
////    params[1] = sigma;
////
////    l = target( params );
////    FVPrint(l, "l - 2");
////
////
////    FVar df1 = FDiff( target, params, num, 1E-10 );
////    FVPRINT(df1);
////
////    FVar df2 = CDiff( target, params, num, 1E-10 );
////    FVPRINT(df2);
//
//
////    FVar a = FVMake(num, 5.0,  0);
////    FVar b = FVMake(num, -1.0, 1);
////
////
////    double c = 3.0;
////
////    FVPRINT(a);
////    FVPRINT(b);
////
////    FVar out = FVAdd( FVMulD(a, c), b );
////
////    FVPRINT(out);

//    Mat A, B, C;
//    u32 n = 5, m = 2, p = 3, i, j;
//
//    A = MatMake(DefaultAllocator, n, m);
//    B = MatMake(DefaultAllocator, m, p);
//    C = MatMake(DefaultAllocator, n, p);
//
//    printf (" Intializing matrix data \n\n");
//
//    for (i=0; i<n; ++i) {
//        for (j=0; j<m; ++j) {
//            SetElement(A, i, j, (f64) i+j);
//        }
//    }
//    MATPRINT(A);
//
//    for (i=0; i<m; ++i) {
//        for (j=0; j<p; ++j) {
//            SetElement(B, i, j, -(f64)i-(f64)j);
//        }
//    }
//    MATPRINT(B);
//
//    MatMul(A, B, C);
//    MATPRINT(C);
//
//    MatMul_BLAS(A, B, C);
//    MATPRINT(C);

    return 0;
}
#endif
