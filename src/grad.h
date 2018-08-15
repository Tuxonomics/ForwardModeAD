#import "fw_univariate.h"


/* Finite Difference */
void FVarFDiff( Allocator al, FVar f( FVarMat ), FVarMat input, FVarMat grad, f64 h ) {
    
    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
    
    FVar tmp;
    u32 N = input.dim0;
    
    FVarMat xCpy = FVarMatMake( al, N, 1 );
    
    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i]     = input.data[i];
        xCpy.data[i].dot = 0;
    }
    
    
    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i].val += h;
        
        tmp = FVMulD( FVSub( f(xCpy), f(input) ), 1/h );
        tmp.dot = 0;
        
        grad.data[i] = tmp;
        
        xCpy.data[i].val -= h;
    }
    
    FVarMatFree( al, &xCpy );
}


/* Central Difference */
void FVarCDiff( Allocator al, FVar f( FVarMat ), FVarMat input, FVarMat grad, f64 h ) {

    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
    
    FVar tmp;
    u32 N = input.dim0;
    
    FVarMat xCpy  = FVarMatMake( al, N, 1 );
    FVarMat xCpy2 = FVarMatMake( al, N, 1 );
    
    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i]      = input.data[i];
        xCpy.data[i].dot  = 0;
        xCpy2.data[i]     = input.data[i];
        xCpy2.data[i].dot = 0;
    }
    
    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i].val  += h;
        xCpy2.data[i].val -= h;
        
        tmp = FVMulD( FVSub( f(xCpy), f(xCpy2) ), 1/(2*h) );
        tmp.dot = 0;
        
        grad.data[i]     = tmp;
        
        xCpy.data[i].val  -= h;
        xCpy2.data[i].val += h;
    }
    
    FVarMatFree( al, &xCpy  );
    FVarMatFree( al, &xCpy2 );
}


/* AD gradient */
void FVarGradient( Allocator al, FVar f( FVarMat ), FVarMat input, FVarMat grad )
{
    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
    
    FVar tmp;
    u32 N = input.dim0;
    
    FVarMat xCpy = FVarMatMake( al, N, 1 );
    
    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i]     = input.data[i];
        xCpy.data[i].dot = 0;
    }
    
    for ( u32 i=0; i<N; ++i ) {
        xCpy.data[i].dot = 1.0;
        
        tmp     = f( xCpy );
        tmp.val = tmp.dot;
        tmp.dot = 0;
        
        grad.data[i]     = tmp;
        xCpy.data[i].dot = 0;
    }
    
    FVarMatFree( al, &xCpy );
}