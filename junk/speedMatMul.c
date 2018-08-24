#include "../src/utilities.h"
#include "../src/grad.h"
#include <time.h>

//#include <immintrin.h>
//#include <emmintrin.h>

#if defined (__has_include) && (__has_include(<x86intrin.h>))
    #include <x86intrin.h>
#else
    #error "No Intel intrinsics file found..."
#endif

//#include <x86intrin.h>

#include "build-cblas/blas/CBLAS/include/cblas.h"



void print256_num(__m256 var)
{
    f64 *val = (f64*) &var;
    printf("%.4f %.4f %.4f %.4f \n",
           val[0], val[1], val[2], val[3]);
}


void print256_mul(__m256 var1, __m256 var2)
{
    f64 *val1 = (f64*) &var1;
    f64 *val2 = (f64*) &var2;
    printf("%.4f %.4f %.4f %.4f \n",
           val1[0]*val2[0], val1[1]*val2[1], val1[2]*val2[2], val1[3]*val2[3]);
}



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

    u32 i, i2, j, j2, k, k2;

    f64 *ra = a.data;
    f64 *rb = b.data;
    f64 *rc = c.data;
    
    for (i = 0; i < n; i += SM) {

        for (j = 0; j < p; j += SM) {

            for (k = 0; k < m; k += SM) {

                for ( i2=i; i2<MIN(i+SM, n); ++i2 ) {

                    rc = &C(i2, j);

                    for ( j2=j; j2<MIN(j+SM, p); ++j2 ) {
                        f64 sum = 0;

                        for ( k2=k; k2<MIN(k+SM, m); ++k2 ) {
                            sum += A(i2, k2) * B(k2, j2);
                        }
                        rc[j2] += sum;
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


void fastMM2( f64Mat a, f64Mat b, f64Mat c )
{
    ASSERT(a.dim0 == c.dim0 && a.dim1 == b.dim0 && b.dim1 == c.dim1);
    
    // A is n x m
    // B is m x p
    // C is n x p
    
    u32 n = a.dim0;
    u32 m = a.dim1;
    u32 p = b.dim1;
    
    i32 i, j, k;
    
    f64 tmpA;
    
    f64 *ra  = a.data;
    f64 *rb  = b.data;
    f64 *rc  = c.data;
//    f64 *rc2 = c.data;

//    for ( i = 0; i < n; ++i ) {
//        for ( k = 0; k < m; ++k ) {
//
//            tmpA = (ra++)[0];
//            // tmpA = A(i, k);
//
//            rc = rc2;
//
//            for ( j = 0; j < p; ++j ) {
//
//                rc[0] = f64Add( rc[0], f64Mul( tmpA, (rb++)[0] ) );
//
//                rc++;
//                // C(i, j) += tmpA * B(k, j);
//
//            }
//
//
//        }
//
//        rc2 = rc;
//        rb = b.data;
//    }

    for ( i = 0; i < n; ++i ) {
        for ( k = 0; k < m; ++k ) {

            //            tmpA = ra[k];


            for ( j = 0; j < p; ++j ) {

                rc[j] = f64Add( rc[j], f64Mul( ra[k], rb[j] ) );

//                rc[j] += ra[k] * rb[j];
                //                rc[j] += tmpA * rb[j];

            }

            rb += p;
        }

        ra += k;
        rc += p;
        rb -= m*p;
    }

}


void fastMM3( f64Mat a, f64Mat b, f64Mat c )
{
    ASSERT(a.dim0 == c.dim0 && a.dim1 == b.dim0 && b.dim1 == c.dim1);

    // A is n x m
    // B is m x p
    // C is n x p

    u32 n = a.dim0;
    u32 m = a.dim1;
    u32 p = b.dim1;

    i32 i, j, k;

    f64 tmpA;

    f64 *ra = a.data;
    f64 *rb = b.data;
    f64 *rc = c.data;

    for ( i = 0; i < n; ++i ) {
        for ( k = 0; k < m; ++k ) {

            tmpA = (ra++)[0];

            for ( j = 0; (j + 4) < p; j += 4) {

                rc[j]   += tmpA * (rb)[0];
                rc[j+1] += tmpA * (rb+1)[0];
                rc[j+2] += tmpA * (rb+2)[0];
                rc[j+3] += tmpA * (rb+3)[0];

                rb += 4;

            }

            for ( ; j < p; ++j ) {
                rc[j] += tmpA * (rb++)[0];
            }
        }

        rc = rc + p;
        rb = b.data;
    }
}


void fastMM4( f64Mat a, f64Mat b, f64Mat c )
{
    ASSERT(a.dim0 == c.dim0 && a.dim1 == b.dim0 && b.dim1 == c.dim1);
    
    // A is n x m
    // B is m x p
    // C is n x p
    
    u32 n = a.dim0;
    u32 m = a.dim1;
    u32 p = b.dim1;
    
    i32 i, j, k;
    
    f64 tmpA;
    
    f64 *ra  = a.data;
    f64 *rb  = b.data;
    f64 *rc  = c.data;
    f64 *rc2 = c.data;
    
    for ( i = 0; i < n; ++i ) {
        for ( k = 0; k < m; ++k ) {
            
//            tmpA = (ra++)[0];
            
            __m256 tmpSIMD = _mm256_set1_pd( (ra++)[0] );
            
            for ( j = 0; j < p; j += 4, rb += 4, rc += 4 ) {
                __m256 bSIMD = _mm256_load_pd( rb );
                __m256 cSIMD = _mm256_load_pd( rc );

                __m256 mul = _mm256_mul_pd( tmpSIMD, bSIMD );
                
                __m256 add = _mm256_add_pd( cSIMD, mul );
                
                _mm256_store_pd( rc, add );
            }
            
            rc = rc2;
        }
        
        
        rc = rc + p;
        rc2 = rc;
        rb = b.data;
    }
}


Inline
void microKernel( f64 * ra, f64 * rb, f64 * rc )

//Inline
//void microKernel( f64 *restrict ra, f64 *restrict rb, f64 *restrict rc )
{
    rc[0] += ra[0] * rb[0];
    rc[1] += ra[0] * rb[1];
    rc[2] += ra[0] * rb[2];
    rc[3] += ra[0] * rb[3];
}


Inline
void macroKernel( f64 * ra, f64 * rb, f64 * rc, u32 n, u32 m, u32 p )
//Inline
//void macroKernel( f64 *restrict ra, f64 *restrict rb, f64 *restrict rc, u32 n, u32 m, u32 p )
{

    u32 i, j, k;

//    f64 tmpA;

    for ( i = 0; i < n; ++i ) {
        for ( k = 0; k < m; ++k ) {

//            tmpA = ra[k];


            for ( j = 0; j < p; j += 4 ) {

                microKernel( ra+k, rb+j, rc+j );

//                rc[j] += ra[k] * rb[j];
//                rc[j] += tmpA * rb[j];

            }

            rb += p;
        }

        ra += k;
        rc += p;
        rb -= m*p;
    }

}


void fastMM5( f64Mat a, f64Mat b, f64Mat c )
{
    ASSERT(a.dim0 == c.dim0 && a.dim1 == b.dim0 && b.dim1 == c.dim1);

    // A is n x m
    // B is m x p
    // C is n x p

    macroKernel( a.data, b.data, c.data, a.dim0, a.dim1, b.dim1 );

}



void FVFastMM( f64FVarMat a, f64FVarMat b, f64FVarMat c )
{
    ASSERT(a.dim0 == c.dim0 && a.dim1 == b.dim0 && b.dim1 == c.dim1);

    // A is n x m
    // B is m x p
    // C is n x p

    u32 n = a.dim0;
    u32 m = a.dim1;
    u32 p = b.dim1;

    i32 i, j, k;

    f64FVar tmpA;

    f64FVar *ra = a.data;
    f64FVar *rb = b.data;
    f64FVar *rc = c.data;

    for ( i = 0; i < n; ++i ) {
        for ( k = 0; k < m; ++k ) {

            tmpA = (ra++)[0];
            // tmpA = A(i, k);

            for ( j = 0; j < p; ++j ) {

                rc[j] = f64FVAdd( rc[j], f64FVMul( tmpA, (rb++)[0] ) );

//                rc[j] += tmpA * (rb++)[0];
                // C(i, j) += tmpA * B(k, j);

            }
        }

        rc = rc + p;

//        __builtin_prefetch( &rc, 1, 1 ); // buffer overflow in the end?

        rb = b.data;
    }
}


int main(int argn, const char ** argv) {
    

//    u32 S = 4000;
//
//    f64FVarMat a = f64FVarMatMake( DefaultAllocator, S, S );
//    f64FVarMat b = f64FVarMatMake( DefaultAllocator, S, S );
//    f64FVarMat c = f64FVarMatZeroMake( DefaultAllocator, S, S );
//    f64FVarMat d = f64FVarMatZeroMake( DefaultAllocator, S, S );
//
//
//    Xorshift1024 x = Xorshift1024Init( 37473 );
//
//    for ( u32 i=0; i<S*S; ++i ) {
//        a.data[i] = f64FVConst( rngXorshift1024NextFloat( &x ) );
//        b.data[i] = f64FVConst( rngXorshift1024NextFloat( &x ) * 2  );
//    }
//
//
//    //    f64MatMul( a, b, c );
//    //
//    //    fastMM2( a, b, c );
//    //
//    //    f64MatPrint( c, "c" );
//    //
//    //
//    //    fastMM( a, b, d );
//    //
//    //    f64MatPrint( d, "d" );
//    //    memset( d.data, 0, d.dim0 * d.dim1 * sizeof(f64) );
//
//
//    f64 dump;
//    clock_t begin0;
//    clock_t end0;
//    double time_spent0;
//
//    // one
//    dump = 0;
//    begin0 = clock();
//
//    FVFastMM( a, b, c );
//
//    dump = c.data[0].val;
//
//    end0 = clock();
//    time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;
//
//    printf("dump = %.4f\n", dump);
//    printf("Time of fastMM2 function call: %.4f sec.\n", time_spent0);





    u32 S = 4000;
    u32 S2 = S;

    f64Mat a = f64MatMake( DefaultAllocator, S, S2 );
    f64Mat b = f64MatMake( DefaultAllocator, S2, S );
    f64Mat c = f64MatZeroMake( DefaultAllocator, S, S );
    f64Mat d = f64MatZeroMake( DefaultAllocator, S, S );


    Xorshift1024 x = Xorshift1024Init( 37473 );

    for ( u32 i=0; i<S*S2; ++i ) {
        a.data[i] = rngXorshift1024NextFloat( &x );
        b.data[i] = rngXorshift1024NextFloat( &x ) * 2;
    }

//    f64MatPrint( a, "a" );
//    f64MatPrint( b, "b" );

//    f64MatMul( a, b, c );
//
////    fastMM3( a, b, c );
////
//    f64MatPrint( c, "c" );
//
////
//    fastMM4( a, b, d );
////
//    f64MatPrint( d, "d" );
//    memset( d.data, 0, d.dim0 * d.dim1 * sizeof(f64) );
//
//
    f64 dump;
    clock_t begin0;
    clock_t end0;
    double time_spent0;

    // one
    dump = 0;
    begin0 = clock();

    fastMM2( a, b, c );

    dump = c.data[0];

    end0 = clock();
    time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;

    printf("dump = %.4f\n", dump);
    printf("Time of fastMM2 function call: %.4f sec.\n", time_spent0);

//    f64MatPrint( c, "c" );


    // two
    dump = 0;
    begin0 = clock();

    fastMM5( a, b, d );

    dump = d.data[0];

    end0 = clock();
    time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;

    printf("dump = %.4f\n", dump);
    printf("Time of fastMM5 function call: %.4f sec.\n", time_spent0);

//    f64MatPrint( d, "d" );

    
    printf("c == d : %d\n", f64MatEqual( c, d, 1E-10 ) );
    
    // three
    memset( c.data, 0, c.dim0 * c.dim1 * sizeof(f64) );

    dump = 0;
    begin0 = clock();

    cblas_dgemm(
        CblasRowMajor, CblasNoTrans, CblasNoTrans,
        a.dim0, b.dim1, a.dim1, 1.0,
        a.data, a.dim1, b.data, b.dim0,
        0.0, c.data, b.dim1
    );

    dump = c.data[0];

    end0 = clock();
    time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;

    printf("dump = %.4f\n", dump);
    printf("Time of cblas_dgemm function call: %.4f sec.\n", time_spent0);

    
//    // four
//    memset( c.data, 0, c.dim0 * c.dim1 * sizeof(f64) );
//
//    dump = 0;
//    begin0 = clock();
//
//    fastMM2( a, b, c );
//
//    dump = c.data[0];
//
//
//    end0 = clock();
//    time_spent0 = (double)(end0 - begin0) / CLOCKS_PER_SEC;
//
//    printf("dump = %.4f\n", dump);
//    printf("Time of fastMM2 function call: %.4f sec.\n", time_spent0);

    return 0;
}


