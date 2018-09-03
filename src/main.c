#include "utilities.h"
#include "forward/grad.h"
//#include "optim.h"




#ifndef TEST
int main(int argc, const char * argv[]) {


//    f64FVarFVar a = f64FVarFVMake(
//        f64FVMake( 0.5, 1.0 ),
//        f64FVMake( 1.0, 0.0)
//    );
//
//    f64FVarFVPrint( a, "a" );
//
//    f64FVarFVar b = f64FVarFVSin( a );
//
//    f64FVarFVPrint( b, "b" );
//
//    f64FVarFVar c = f64FVarFVCos( a );
//
//    f64FVarFVPrint( c, "c" );


    f64FVarFVar a = f64FVarFVMake(
        f64FVMake( 0.5, 0.0 ),
        f64FVMake( 0.0, 0.0)
    );

    f64FVarFVar b = f64FVarFVMake(
        f64FVMake( 0.0, 1.0 ),
        f64FVMake( 1.0, 0.0)
    );

    f64FVarFVarMat m = f64FVarFVarMatMake( DefaultAllocator, 2, 1 );
    m.data[0] = a;
    m.data[1] = b;

    f64FVarFVarMatPrint( m, "m" );

    f64FVarFVar c = test_f2( m );

    f64FVarFVPrint( c, "c" );

    return 0;
}
#endif
