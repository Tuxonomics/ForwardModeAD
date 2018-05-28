
// Function from https://www.youtube.com/watch?v=r2hhRSHiQwY at 10:10

#include "fwUnivariate.h"
#include <time.h>

int main(int argn, const char ** argv) {
    
    u32 N = 1410065408;
    
    
    FVar x = FVMake(0.0, 1.0);
    
    FVar out;
    f64 dump = 0;
    
    clock_t begin0 = clock();
    
    for (u32 i=0; i<N; ++i) {
        out = FVDiv( FVExp(x), FVSqrt( FVAdd( FVPow( FVCos(x), 3.0) , FVPow( FVSin(x), 3.0) ) ) );
        dump += out.val;
    }
    
    clock_t end0 = clock();
    double time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;
    
    printf("dump = %.4f\n", dump);
    printf("Time of AD function call: %.4f sec.\n", time_spent0);
    
    
    dump = 0;
    f64 y = 0.0;
    
    clock_t begin1 = clock();
    
    for (u32 i=0; i<N; ++i) {
        dump += exp(y) / sqrt( pow(cos(y), 3.0) + pow(sin(y), 3.0) );
    }
    
    clock_t end1 = clock();
    double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
    
    printf("dump = %.4f\n", dump);
    printf("Time of normal function call: %.4f sec.\n", time_spent1);
    
    return 0;
}
