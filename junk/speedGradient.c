#include "../src/utilities.h"
#include "../src/grad.h"
#include <time.h>

int main(int argn, const char ** argv) {

#define P 1E-8

    u32 N = 50000;


    f64Mat input  = f64MatMake( DefaultAllocator, N, 1 );
    f64Mat gradC  = f64MatMake( DefaultAllocator, N, 1 );
    f64Mat gradF  = f64MatMake( DefaultAllocator, N, 1 );
    f64Mat gradAD = f64MatMake( DefaultAllocator, N, 1 );

    input.data[0] = 0.5;
    input.data[1] = 0.0;

    // one
    f64 dump = 0;
    clock_t begin0 = clock();

    for (u32 i=0; i<N; ++i) {
        f64FVarFDiff( DefaultAllocator, test_f, input, gradF, P );
        dump += gradF.data[0];
    }

    clock_t end0 = clock();
    double time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;

    printf("dump = %.4f\n", dump);
    printf("Time of f64FVarFDiff function call: %.4f sec.\n", time_spent0);


    // two
    dump = 0;
    begin0 = clock();

    for (u32 i=0; i<N; ++i) {
        f64FVarCDiff( DefaultAllocator, test_f, input, gradC, P );
        dump += gradC.data[0];
    }

    end0 = clock();
    time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;

    printf("dump = %.4f\n", dump);
    printf("Time of f64FVarCDiff function call: %.4f sec.\n", time_spent0);


    // three
    dump = 0;
    begin0 = clock();

    for (u32 i=0; i<N; ++i) {
        f64FVarGradient( DefaultAllocator, test_f, input, gradAD );
        dump += gradAD.data[0];
    }

    end0 = clock();
    time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;

    printf("dump = %.4f\n", dump);
    printf("Time of f64FVarGradient function call: %.4f sec.\n", time_spent0);


    f64MatFree( DefaultAllocator, &input );
    f64MatFree( DefaultAllocator, &gradC );
    f64MatFree( DefaultAllocator, &gradF );
    f64MatFree( DefaultAllocator, &gradAD );

#undef P

    return 0;
}
