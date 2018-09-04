// Author:  https://github.com/Tuxonomics
// Created: Aug, 2018
//

#include "fw_dod.h"

///* Finite Difference */
//void f64FVarFDiff( Allocator al, f64FVar f( f64FVarMat ), f64Mat input, f64Mat grad, f64 h )
//{
//
//    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
//
//    f64FVar tmp;
//    f64FVar tmpF;
//
//    u32 N = input.dim0;
//
//    f64FVarMat xCpy  = f64FVarMatMake( al, N, 1 );
//    f64FVarMat xCpy2 = f64FVarMatMake( al, N, 1 );
//
//    for ( u32 i=0; i<N; ++i ) {
//        xCpy.data[i]  = f64FVConst( input.data[i] );
//        xCpy2.data[i] = f64FVConst( input.data[i] );
//    }
//
//    tmpF = f(xCpy2);
//
//    for ( u32 i=0; i<N; ++i ) {
//        xCpy.data[i].val += h;
//
//        tmp = f64FVMulf64( f64FVSub( f(xCpy), tmpF ), 1/h );
//
//        grad.data[i] = tmp.val;
//
//        xCpy.data[i].val -= h;
//    }
//
//    f64FVarMatFree( al, &xCpy );
//    f64FVarMatFree( al, &xCpy2 );
//}
//
//
///* Central Difference */
//void f64FVarCDiff( Allocator al, f64FVar f( f64FVarMat ), f64Mat input, f64Mat grad, f64 h )
//{
//
//    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
//
//    f64FVar tmp;
//    u32 N = input.dim0;
//
//    f64FVarMat xCpy  = f64FVarMatMake( al, N, 1 );
//    f64FVarMat xCpy2 = f64FVarMatMake( al, N, 1 );
//
//    for ( u32 i=0; i<N; ++i ) {
//        xCpy.data[i]  = f64FVConst( input.data[i] );
//        xCpy2.data[i] = f64FVConst( input.data[i] );
//    }
//
//    for ( u32 i=0; i<N; ++i ) {
//        xCpy.data[i].val  += h;
//        xCpy2.data[i].val -= h;
//
//        tmp = f64FVMulf64( f64FVSub( f(xCpy), f(xCpy2) ), 1/(2*h) );
//
//        grad.data[i]     = tmp.val;
//
//        xCpy.data[i].val  -= h;
//        xCpy2.data[i].val += h;
//    }
//
//    f64FVarMatFree( al, &xCpy  );
//    f64FVarMatFree( al, &xCpy2 );
//}

/* AD gradient */
void f64FVGradient( Allocator al, f64FVar f( f64FVar ), f64Mat input, f64Mat grad )
{
    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
    
    f64FVar tmp;
    u32 N = input.dim0;
    
    f64FVar xCpy = f64FVMake( al, N, 1 );
    
    for ( u32 i=0; i<N; ++i ) {
        f64FVSetElement( xCpy, i, 0, input.data[i], 0.0 );
    }
    
    for ( u32 i=0; i<N; ++i ) {
        xCpy.dot.data[i] = 1.0;

        tmp     = f( xCpy );

        grad.data[i]     = tmp.dot.data[0];
        xCpy.dot.data[i] = 0.0;
    }
    
    f64FVFree( al, &xCpy );
}

//
///* test function, will be refactored once the test tool is updated */
//f64FVar test_f( f64FVarMat input )
//{
//    return f64FVAdd( f64FVTanh( input.data[0] ), f64FVSin( input.data[1] ) );
//}
//
//#if TEST
//void test_grad()
//{
//#define EPS 1E-8
//
//    u32 N = 2;
//
//    f64Mat input  = f64MatMake( DefaultAllocator, N, 1 );
//    f64Mat gradC  = f64MatMake( DefaultAllocator, N, 1 );
//    f64Mat gradF  = f64MatMake( DefaultAllocator, N, 1 );
//    f64Mat gradAD = f64MatMake( DefaultAllocator, N, 1 );
//
//    input.data[0] = 0.5;
//    input.data[1] = 0.0;
//
//    f64FVarFDiff( DefaultAllocator, test_f, input, gradF, EPS );
//
//    f64FVarCDiff( DefaultAllocator, test_f, input, gradC, EPS );
//
//    f64FVarGradient( DefaultAllocator, test_f, input, gradAD );
//
//
//    TEST_ASSERT( f64MatEqual( gradF, gradC,  EPS ) );
//    TEST_ASSERT( f64MatEqual( gradC, gradAD, EPS ) );
//
//
//    f64MatFree( DefaultAllocator, &input );
//    f64MatFree( DefaultAllocator, &gradC );
//    f64MatFree( DefaultAllocator, &gradF );
//    f64MatFree( DefaultAllocator, &gradAD );
//
//#undef EPS
//}
//#endif

