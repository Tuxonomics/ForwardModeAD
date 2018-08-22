#import "fw_univariate.h"


/* Finite Difference */
void f64FVarFDiff( Allocator al, f64FVar f( f64FVarMat ), f64Mat input, f64Mat grad, f64 h )
{

    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );

    f64FVar tmp;
    f64FVar tmpF;

    u32 N = input.dim0;

    f64FVarMat xCpy  = f64FVarMatMake( al, N, 1 );
    f64FVarMat xCpy2 = f64FVarMatMake( al, N, 1 );

    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i]  = f64FVConst( input.data[i] );
        xCpy2.data[i] = f64FVConst( input.data[i] );
    }

    tmpF = f(xCpy2);

    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i].val += h;

        tmp = f64FVMulf64( f64FVSub( f(xCpy), tmpF ), 1/h );

        grad.data[i] = tmp.val;

        xCpy.data[i].val -= h;
    }

    f64FVarMatFree( al, &xCpy );
    f64FVarMatFree( al, &xCpy2 );
}


/* Central Difference */
void f64FVarCDiff( Allocator al, f64FVar f( f64FVarMat ), f64Mat input, f64Mat grad, f64 h )
{

    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );

    f64FVar tmp;
    u32 N = input.dim0;

    f64FVarMat xCpy  = f64FVarMatMake( al, N, 1 );
    f64FVarMat xCpy2 = f64FVarMatMake( al, N, 1 );

    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i]  = f64FVConst( input.data[i] );
        xCpy2.data[i] = f64FVConst( input.data[i] );
    }

    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i].val  += h;
        xCpy2.data[i].val -= h;

        tmp = f64FVMulf64( f64FVSub( f(xCpy), f(xCpy2) ), 1/(2*h) );

        grad.data[i]     = tmp.val;

        xCpy.data[i].val  -= h;
        xCpy2.data[i].val += h;
    }

    f64FVarMatFree( al, &xCpy  );
    f64FVarMatFree( al, &xCpy2 );
}


/* AD gradient */
void f64FVarGradient( Allocator al, f64FVar f( f64FVarMat ), f64Mat input, f64Mat grad )
{
    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );

    f64FVar tmp;
    u32 N = input.dim0;

    f64FVarMat xCpy = f64FVarMatMake( al, N, 1 );

    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i] = f64FVConst( input.data[i] );
    }

    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i].dot = 1.0;

        tmp     = f( xCpy );

        grad.data[i]     = tmp.dot;
        xCpy.data[i].dot = 0;
    }

    f64FVarMatFree( al, &xCpy );
}


/* test function, will be refactored once the test tool is updated */
f64FVar test_f( f64FVarMat input )
{
    return f64FVAdd( f64FVTanh( input.data[0] ), f64FVSin( input.data[1] ) );
}

#if TEST
void test_grad()
{
#define P 1E-8

    u32 N = 2;

    f64Mat input  = f64MatMake( DefaultAllocator, N, 1 );
    f64Mat gradC  = f64MatMake( DefaultAllocator, N, 1 );
    f64Mat gradF  = f64MatMake( DefaultAllocator, N, 1 );
    f64Mat gradAD = f64MatMake( DefaultAllocator, N, 1 );

    input.data[0] = 0.5;
    input.data[1] = 0.0;

    f64FVarFDiff( DefaultAllocator, test_f, input, gradF, P );

    f64FVarCDiff( DefaultAllocator, test_f, input, gradC, P );

    f64FVarGradient( DefaultAllocator, test_f, input, gradAD );


    TEST_ASSERT( f64MatEqual( gradF, gradC,  P ) );
    TEST_ASSERT( f64MatEqual( gradC, gradAD, P ) );

//    f64MatPrint( gradAD,  "gradAD" );

    f64MatFree( DefaultAllocator, &input );
    f64MatFree( DefaultAllocator, &gradC );
    f64MatFree( DefaultAllocator, &gradF );
    f64MatFree( DefaultAllocator, &gradAD );

#undef P
}
#endif


/* numerical hessian based on finite differences */
void f64FVarNumHess( Allocator al, f64FVarFVar f( f64FVarFVarMat ), f64Mat input, f64Mat hess, f64 h )
{
#define Hess(i,j) hess.data[i*hess.dim1 + j]

    ASSERT( input.dim1 == 1 );
    ASSERT( hess.dim0 == hess.dim1 && hess.dim0 == MAX(input.dim0, input.dim1) );

    f64FVarFVar tmpF;
    f64FVarFVar tmp0;
    f64FVarFVar tmp1;
    f64FVarFVar tmp01;
    f64FVarFVar tmp;

    u32 N = input.dim0;

    f64FVarFVarMat xCpy = f64FVarFVarMatMake( al, N, 1 );

    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i] = f64FVarFVMake( f64FVConst( input.data[i] ), f64FVConst( 0 ) );
    }

    tmpF = f( xCpy );

    for ( u32 i=0; i<N; ++i ) {


        for ( u32 j=i; j<N; ++j ) {
            xCpy.data[i].val.val += h;
            tmp0 = f( xCpy );

            xCpy.data[j].val.val += h;
            tmp01 = f( xCpy );

            xCpy.data[i].val.val -= h;
            tmp1 = f( xCpy );

            xCpy.data[j].val.val -= h;


            tmp = f64FVarFVSub( f64FVarFVAdd( tmpF, tmp01 ), f64FVarFVAdd( tmp0, tmp1 ) );

            tmp = f64FVarFVMulf64FVar( tmp, f64FVConst( 1/(h*h) ) );


            Hess(i, j) = tmp.val.val;
            Hess(j, i) = Hess(i, j);
        }

    }

    f64FVarFVarMatFree( al, &xCpy );

#undef Hess
}


/* hessian and gradient with second-order forward AD variables */
void f64FVarHessian( Allocator al, f64FVarFVar f( f64FVarFVarMat ), f64Mat input, f64Mat grad, f64Mat hess )
{
#define Hess(i,j) hess.data[i*hess.dim1 + j]

    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
    ASSERT( hess.dim0 == hess.dim1 && hess.dim0 == MAX(grad.dim0, grad.dim1) );

    f64FVarFVar tmp;
    u32 N = input.dim0;

    f64FVarFVarMat xCpy = f64FVarFVarMatMake( al, N, 1 );

    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i] = f64FVarFVMake( f64FVConst( input.data[i] ), f64FVConst( 0 ) );
    }

    for ( u32 i=0; i<N; ++i ) {
        for ( u32 j=i; j<N; ++j ) {

            for (u32 k=0; k<N; ++k) {
                xCpy.data[k].val.dot = (i32) (j == k);
                xCpy.data[k].dot.val = (i32) (i == k);
            }

            tmp = f( xCpy );

            if ( i == j ) {
                grad.data[i] = tmp.dot.val;
            }

            Hess(i, j) = tmp.dot.dot;
            Hess(j, i) = Hess(i, j);

        }
    }

    f64FVarFVarMatFree( al, &xCpy );

#undef Hess
}



/* test function, will be refactored once the test tool is updated */
f64FVarFVar test_f2( f64FVarFVarMat input )
{
    return f64FVarFVAdd(
        f64FVarFVTanh( input.data[0] ),
        f64FVarFVMul( input.data[0], f64FVarFVPow( input.data[1], 2 ) )
    );
}


#if TEST
void test_hessian()
{
#define P 1E-5

    u32 N = 2;

    f64Mat input  = f64MatMake( DefaultAllocator, N, 1 );
    f64Mat grad   = f64MatMake( DefaultAllocator, N, 1 );
    f64Mat hess   = f64MatMake( DefaultAllocator, N, N );
    f64Mat num    = f64MatMake( DefaultAllocator, N, N );

    input.data[0] = 0.5;
    input.data[1] = 2.0;


    f64FVarHessian( DefaultAllocator, test_f2, input, grad, hess );

    f64FVarNumHess( DefaultAllocator, test_f2, input, num, P );


    TEST_ASSERT( f64MatEqual( hess, num, P ) );


    f64MatFree( DefaultAllocator, &input );
    f64MatFree( DefaultAllocator, &grad );
    f64MatFree( DefaultAllocator, &hess );

#undef P
}
#endif





