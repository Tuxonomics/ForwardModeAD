#include "../src/utilities.h"
#include "../src/grad.h"
#include <time.h>

#define CLS 64

void fastMM( f64Mat a, f64Mat b, f64Mat c )
{
#define SM (CLS / sizeof (f64))
#define A(i, j) a.data[i*m + j]
#define B(i, j) b.data[i*p + j]
#define C(i, j) c.data[i*p + j]

    ASSERT(a.dim0 == c.dim0 && a.dim1 == b.dim0 && b.dim1 == c.dim1); \

    // A is n x m
    // B is m x p
    // C is n x p

    /* NOTE(jonas): for row-major ordering the algorithm needs to increment by
    row size */

    u32 n = a.dim0;
    u32 m = a.dim1;
    u32 p = b.dim1;

    i32 i, i2, j, j2, k, k2;

    u32 count = 0;

    for (i = 0; i < n; i += SM) {

        for (j = 0; j < p; j += SM) {

            for (k = 0; k < m; k += SM) {

                for ( i2=i; i2<MIN(i+SM, n); ++i2 ) {

                    for ( j2=j; j2<MIN(j+SM, p); ++j2 ) {
                        f64 sum = 0;

                        for ( k2=k; k2<MIN(k+SM, m); ++k2 ) {
                            sum += A(i2, k2) * B(k2, j2);
                        }
                        C(i2, j2) += sum;
                    }
                }
            }
        }
    }

#undef SM
#undef A
#undef B
#undef C
}


int main(int argn, const char ** argv) {

    u32 N = 1;

    u32 S = 4000;

    f64Mat a = f64MatMake( DefaultAllocator, S, S );
    f64Mat b = f64MatMake( DefaultAllocator, S, S );
    f64Mat c = f64MatMake( DefaultAllocator, S, S );
    f64Mat d = f64MatZeroMake( DefaultAllocator, S, S );


    Xorshift1024 x = Xorshift1024Init( 37473 );

    for ( u32 i=0; i<S*S; ++i ) {
        a.data[i] = rngXorshift1024NextFloat( &x );
        b.data[i] = rngXorshift1024NextFloat( &x ) * 2;
    }


//    f64MatMul( a, b, c );
//
//    f64MatPrint( c, "c" );
//
//
//    fastMM( a, b, d );
//
//    f64MatPrint( d, "d" );
//    memset( d.data, 0, d.dim0 * d.dim1 * sizeof(f64) );



    // one
    f64 dump = 0;
    clock_t begin0 = clock();

//    for (u32 i=0; i<N; ++i) {
//        f64MatMul( a, b, c );
//        dump += c.data[0];
//    }

    clock_t end0 = clock();
    double time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;

    printf("dump = %.4f\n", dump);
    printf("Time of f64MatMul function call: %.4f sec.\n", time_spent0);

//    f64MatPrint( c, "c" );


    // two
    dump = 0;
    begin0 = clock();

    for (u32 i=0; i<N; ++i) {
//        memset( d.data, 0, d.dim0 * d.dim1 * sizeof(f64) );
        fastMM( a, b, d );
        dump += d.data[0];
    }

    end0 = clock();
    time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;

    printf("dump = %.4f\n", dump);
    printf("Time of fastMM function call: %.4f sec.\n", time_spent0);

//    f64MatPrint( d, "d" );


    return 0;
}

