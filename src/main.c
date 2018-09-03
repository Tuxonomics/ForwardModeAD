#include "fw_dod_grad.h"





#ifndef TEST
int main(int argc, const char * argv[]) {

    InitializeFV( 4, 'l' );
    
    u32 N = 4000;
    
    f64FVar fv1 = f64FVMake( DefaultAllocator, N, N );
    f64FVar fv2 = f64FVMake( DefaultAllocator, N, N );
    f64FVar fv3 = f64FVMake( DefaultAllocator, N, N );
    
    
    for ( u32 i=0; i<N; ++i ) {
        for ( u32 j=0; j<N; ++j ) {
            f64FVSetElement( fv1, i, j, (f64) (i+0.5) / (j+1), 0.0 );
            f64FVSetElement( fv2, i, j, (f64) (j+0.5) / (i+1), 0.0 );
        }
    }
    
    
    f64FVMatMul( fv1, fv2, fv3 );

    printf("%.4f\n", fv3.val.data[0]);
    
//    f64FVPrint( fv3, "fv3" );
    
//    f64FVFree( DefaultAllocator, &fv1 );
//    f64FVFree( DefaultAllocator, &fv2 );
//    f64FVFree( DefaultAllocator, &fv3 );

    return 0;
}
#endif
