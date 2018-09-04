// Author:  https://github.com/Tuxonomics
// Created: Aug, 2018
//

#include "dependencies/utilities.h"

static Arena     FVScratchArena;
static Allocator FVScratchBuffer;


void InitializeFV( u32 numThreads, char threadScope )
{
    ArenaInit( &FVScratchArena, DefaultAllocator, KB(1) );
    FVScratchBuffer = ArenaAllocatorMake( &FVScratchArena );
    
    InitializeMatrices( numThreads, threadScope );
    
}


void TerminateFV( void )
{
    ArenaDestroy( &FVScratchArena );
    
    TerminateMatrices();
}


#define EPS 1E-10


#define FVAR_DECL(type) typedef struct type##FVar type##FVar; \
    struct type##FVar { \
        u32 dim0; \
        u32 dim1; \
        type##Mat val; \
        type##Mat dot; \
    };

#define FVAR_MAKE(type) type##FVar \
    type##FVMake(Allocator al, u32 dim0, u32 dim1) \
    { \
        type##FVar fv;\
        fv.dim0 = dim0; \
        fv.dim1 = dim1; \
        fv.val  = type##MatMake( al, dim0, dim1 ); \
        fv.dot  = type##MatMake( al, dim0, dim1 ); \
        return fv; \
    }

#define FVAR_FREE(type) void \
    type##FVFree( Allocator al, type##FVar *fv ) \
    { \
        type##MatFree( al, &fv->val );\
        type##MatFree( al, &fv->dot );\
    }

#define FVAR_CONST(type, baseFun) type##FVar \
    type##FVConst(Allocator al, u32 dim0, u32 dim1, type val) \
    { \
    type##FVar fv = type##FVMake( al, dim0, dim1 ); \
    for ( u32 i=0; i<(dim0*dim1); ++i ) { \
        fv.val.data[i] = val; \
        fv.dot.data[i] = baseFun( 0 ); \
    } \
    return fv; \
    }

#define FVAR_SETELEMENT(type) void \
    type##FVSetElement( type##FVar fv, u32 dim0, u32 dim1, type val, type dot ) \
    { \
        fv.val.data[ dim0 * fv.dim1 + dim1 ] = val; \
        fv.dot.data[ dim0 * fv.dim1 + dim1 ] = dot; \
    }

#define FVAR_GETELEMENT(type) type##FVar \
    type##FVGetElement(Allocator al, type##FVar fv, u32 dim0, u32 dim1) \
    { \
        type##FVar el = type##FVMake( al, dim0, dim1 );\
        el.val.data[0] = fv.val.data[ dim0 * fv.dim1 + dim1 ]; \
        el.dot.data[0] = fv.dot.data[ dim0 * fv.dim1 + dim1 ]; \
        return el; \
    }

#define FVAR_PRINT(type, baseFun) void \
    type##FVPrint( type##FVar fv, const char* name) \
    { \
        printf("FVar (%s): {\n", name); \
        printf("\t.dim0 = %u\n", fv.dim0); \
        printf("\t.dim1 = %u\n", fv.dim1); \
        printf("\t.data = {\n\n"); \
        u32 dim; \
        for (u32 i=0; i<fv.dim0; ++i) { \
            for (u32 j=0; j<fv.dim1; ++j) { \
                dim = i*fv.dim1 + j; \
                printf("[%d] ", dim); \
                baseFun(fv.val.data[dim]); \
                printf("\t"); \
                baseFun(fv.dot.data[dim]); \
                printf("  "); \
            } \
            printf("\n"); \
        } \
        printf("\n\t}\n"); \
        printf("}\n"); \
    }

#define FVAR_EQUAL(type, equalFun) b32 \
    type##FVEqual( type##FVar a, type##FVar b, f64 eps ) \
    { \
        if ( a.dim0 != b.dim0 ) \
            return 0; \
        if ( a.dim1 != b.dim1 ) \
            return 0; \
        \
        for ( u32 i=0; i<(a.dim0 * a.dim1); ++i ) { \
            if ( ! equalFun( a.val.data[i], b.val.data[i], eps ) ) \
                return 0; \
        } \
        for ( u32 i=0; i<(a.dim0 * a.dim1); ++i ) { \
            if ( ! equalFun( a.dot.data[i], b.dot.data[i], eps ) ) \
                return 0; \
        } \
        return 1; \
    }

#define FVAR_COPY(type, baseFun) void \
    type##FVCopy( type##FVar src, type##FVar dst ) \
    { \
        ASSERT( src.dim0 == dst.dim0 ); \
        ASSERT( src.dim1 == dst.dim1 ); \
        \
        u32 copySize = src.dim0 * src.dim1; \
        \
        for ( u32 i=0; i<copySize; ++i ) { \
            baseFun( src.val.data[i], dst.val.data+i ); \
        } \
        for ( u32 i=0; i<copySize; ++i ) { \
            baseFun( src.dot.data[i], dst.dot.data+i ); \
        } \
    }

#define FVAR_ADD(type, addFun) void \
    type##FVAdd( type##FVar src, type##FVar dst ) \
    { \
        ASSERT( src.dim0 == dst.dim0 ); \
        ASSERT( src.dim1 == dst.dim1 ); \
        \
        for ( u32 i=0; i<(src.dim0*src.dim1); ++i ) { \
            dst.val.data[i] = addFun(dst.val.data[i], src.val.data[i]); \
        } \
        for ( u32 i=0; i<(src.dim0*src.dim1); ++i ) { \
            dst.dot.data[i] = addFun(dst.dot.data[i], src.dot.data[i]); \
        } \
    }

#define FVAR_ADD_TYPE(type, addFun) void \
    type##FVAdd##type( type##FVar fv, type val ) \
    { \
        for ( u32 i=0; i<(fv.dim0*fv.dim1); ++i ) { \
            fv.val.data[i] = addFun(fv.val.data[i], val); \
        } \
    }

#define FVAR_SUB(type, subFun) void \
    type##FVSub( type##FVar src, type##FVar dst ) \
    { \
        ASSERT( src.dim0 == dst.dim0 ); \
        ASSERT( src.dim1 == dst.dim1 ); \
        \
        for ( u32 i=0; i<(src.dim0*src.dim1); ++i ) { \
            dst.val.data[i] = subFun(dst.val.data[i], src.val.data[i]); \
        } \
        for ( u32 i=0; i<(src.dim0*src.dim1); ++i ) { \
            dst.dot.data[i] = subFun(dst.dot.data[i], src.dot.data[i]); \
        } \
    }

#define FVAR_MUL(type, addFun, mulFun) void \
    type##FVMul( type##FVar src, type##FVar dst ) \
    { \
        ASSERT( src.dim0 == dst.dim0 ); \
        ASSERT( src.dim1 == dst.dim1 ); \
        \
        type tmp; \
        for ( u32 i=0; i<(src.dim0*src.dim1); ++i ) { \
            tmp             = dst.val.data[i]; \
            \
            dst.val.data[i] = mulFun(src.val.data[i], tmp); \
            dst.dot.data[i] = addFun( \
                mulFun( src.val.data[i], dst.dot.data[i] ), \
                mulFun( src.dot.data[i], tmp ) \
            ); \
        } \
    }

#define FVAR_MUL_TYPE(type, mulFun) void \
    type##FVMul##type( type##FVar fv, type val ) \
    { \
        for ( u32 i=0; i<(fv.dim0*fv.dim1); ++i ) { \
            fv.val.data[i] = mulFun(fv.val.data[i], val); \
        } \
        for ( u32 i=0; i<(fv.dim0*fv.dim1); ++i ) { \
            fv.dot.data[i] = mulFun(fv.dot.data[i], val); \
        } \
    }

#define FVAR_DIV(type, subFun, mulFun, divFun) void \
    type##FVDiv( type##FVar src, type##FVar dst ) \
    { \
        ASSERT( src.dim0 == dst.dim0 ); \
        ASSERT( src.dim1 == dst.dim1 ); \
        \
        type tmp; \
        for ( u32 i=0; i<(src.dim0*src.dim1); ++i ) { \
            tmp             = dst.val.data[i]; \
            \
            dst.val.data[i] = divFun(src.val.data[i], tmp); \
            dst.dot.data[i] = divFun( \
                subFun( \
                    mulFun( src.dot.data[i], tmp ), \
                    mulFun( src.val.data[i], dst.dot.data[i] ) \
                ), \
                mulFun( tmp, tmp ) \
            ); \
        } \
    }

#define FVAR_DIV_TYPE(type, divFun) void \
    type##FVDiv##type( type##FVar fv, type val ) \
    { \
        for ( u32 i=0; i<(fv.dim0*fv.dim1); ++i ) { \
            fv.val.data[i] = divFun(fv.val.data[i], val); \
        } \
        for ( u32 i=0; i<(fv.dim0*fv.dim1); ++i ) { \
            fv.dot.data[i] = divFun(fv.dot.data[i], val); \
        } \
    }

#define FVAR_NEG(type, negFun) void \
    type##FVNeg( type##FVar fv ) \
    { \
        for ( u32 i=0; i<(fv.dim0*fv.dim1); ++i ) { \
            fv.val.data[i] = negFun( fv.val.data[i] ); \
        } \
        for ( u32 i=0; i<(fv.dim0*fv.dim1); ++i ) { \
            fv.dot.data[i] = negFun( fv.dot.data[i] ); \
        } \
    }


/* elementary element-wise functions */

#define FVAR_EXP(type, mulFun, expFun) void \
    type##FVExp( type##FVar fv ) \
    { \
        type tmp; \
        for ( u32 i=0; i<(fv.dim0*fv.dim1); ++i ) { \
            tmp = expFun( fv.val.data[i] ); \
             \
            fv.val.data[i] = tmp; \
            fv.dot.data[i] = mulFun( fv.dot.data[i], tmp ); \
        } \
    }

#define FVAR_LOG(type, divFun, logFun) void \
    type##FVLog( type##FVar fv ) \
    { \
        type tmp; \
        for ( u32 i=0; i<(fv.dim0*fv.dim1); ++i ) { \
            tmp = fv.val.data[i]; \
            \
            fv.val.data[i] = logFun( tmp ); \
            fv.dot.data[i] = divFun( fv.dot.data[i], tmp ); \
        } \
    }


/* matrix-based functions */

#define FVAR_MATADD(type, addFun) void \
    type##FVMatAdd( type##FVar a, type##FVar b, type##FVar dst ) \
    { \
        ASSERT(a.dim0 == b.dim0 && a.dim1 == b.dim1 && a.dim0 == dst.dim0 && a.dim1 == dst.dim1); \
        \
        addFun( a.val, b.val, dst.val ); \
        addFun( a.dot, b.dot, dst.dot ); \
    }

#define FVAR_MATSUB(type, subFun) void \
    type##FVMatSub( type##FVar a, type##FVar b, type##FVar dst ) \
    { \
        ASSERT(a.dim0 == b.dim0 && a.dim1 == b.dim1 && a.dim0 == dst.dim0 && a.dim1 == dst.dim1); \
        \
        subFun( a.val, b.val, dst.val ); \
        subFun( a.dot, b.dot, dst.dot ); \
    }

/* in-place addition to dst.dot */
#define FVAR_MATMUL(type, mulFun, mulFunIP) void \
    type##FVMatMul( type##FVar a, type##FVar b, type##FVar dst )\
    { \
        ASSERT(a.dim0 == dst.dim0 && a.dim1 == b.dim0 && b.dim1 == dst.dim1); \
        \
        mulFun(   a.dot, b.val, dst.dot ); \
        mulFunIP( a.val, b.dot, dst.dot ); \
        \
        mulFun( a.val, b.val, dst.val ); \
    }


FVAR_DECL(f64);
FVAR_MAKE(f64);
FVAR_FREE(f64);
FVAR_CONST(f64, f64Const);
FVAR_EQUAL(f64, f64Equal);
FVAR_SETELEMENT(f64);
FVAR_GETELEMENT(f64);
FVAR_PRINT(f64, f64Print);
FVAR_COPY(f64, f64Copy);

void f64FVCopyFast( f64FVar src, f64FVar dst )
{
    ASSERT( src.dim0 == dst.dim0 );
    ASSERT( src.dim1 == dst.dim1 );
    
    u32 copySize = src.dim0 * src.dim1 * sizeof(f64);
    
    memcpy( src.val.data, dst.val.data, copySize );
    memcpy( src.dot.data, dst.dot.data, copySize );
}


#if TEST
void test_fvar_basics()
{
    u32 N = 3;
    
    f64FVar fv1 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv1, 0, 0, 3.0, 0.0 );
    f64FVSetElement( fv1, 1, 0, 2.5, 1.0 );
    f64FVSetElement( fv1, 2, 0, 0.5, 0.0 );
    
    
    f64FVar fv2 = f64FVConst( DefaultAllocator, N, 1, 4.0 );
    f64FVar fv3 = f64FVConst( DefaultAllocator, N, 1, 5.0 );
    
    TEST_ASSERT( f64FVEqual( fv2, fv2, EPS ) );
    TEST_ASSERT( ! f64FVEqual( fv2, fv3, EPS ) );
    
    f64FVCopy( fv1, fv2 );
    
    TEST_ASSERT( f64FVEqual( fv1, fv2, EPS ) );
    
    f64FVFree( DefaultAllocator, &fv1 );
    f64FVFree( DefaultAllocator, &fv2 );
    f64FVFree( DefaultAllocator, &fv3 );
}
#endif


FVAR_ADD(f64, f64Add);
FVAR_ADD_TYPE(f64, f64Add);
FVAR_SUB(f64, f64Sub);
FVAR_MUL(f64, f64Add, f64Mul);
FVAR_MUL_TYPE(f64, f64Mul);
FVAR_DIV(f64, f64Sub, f64Mul, f64Div);
FVAR_DIV_TYPE(f64, f64Div);
FVAR_NEG(f64, f64Neg);


#if TEST
void test_fvar_add_sub()
{
    u32 N = 3;
    
    f64FVar fv1 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv1, 0, 0, 3.0, 0.0 );
    f64FVSetElement( fv1, 1, 0, 2.5, 1.0 );
    f64FVSetElement( fv1, 2, 0, 0.5, 0.0 );
    
    f64FVar fv2 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv2, 0, 0, 6.0, 0.0 );
    f64FVSetElement( fv2, 1, 0, 5.0, 2.0 );
    f64FVSetElement( fv2, 2, 0, 1.0, 0.0 );
    
    f64FVAdd( fv1, fv1 );
    
    TEST_ASSERT( f64FVEqual( fv1, fv2, EPS ) );
    
    f64FVSetElement( fv2, 0, 0, 6.5, 0.0 );
    f64FVSetElement( fv2, 1, 0, 5.5, 2.0 );
    f64FVSetElement( fv2, 2, 0, 1.5, 0.0 );
    
    f64FVAddf64( fv1, 0.5 );
    
    TEST_ASSERT( f64FVEqual( fv1, fv2, EPS ) );
    
    f64FVar fv3 = f64FVConst( DefaultAllocator, N, 1, 0.0 );
    
    f64FVSub( fv1, fv2 );
    
    TEST_ASSERT( f64FVEqual( fv2, fv3, EPS ) );
    
    f64FVFree( DefaultAllocator, &fv1 );
    f64FVFree( DefaultAllocator, &fv2 );
    f64FVFree( DefaultAllocator, &fv3 );
}
#endif

#if TEST
void test_fvar_mul()
{
    u32 N = 3;
    
    f64FVar fv1 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv1, 0, 0, 3.0, 0.0 );
    f64FVSetElement( fv1, 1, 0, 2.5, 1.0 );
    f64FVSetElement( fv1, 2, 0, 0.5, 0.0 );
    
    f64FVar fv2 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv2, 0, 0, 2.0, 0.0 );
    f64FVSetElement( fv2, 1, 0, 2.0, 0.0 );
    f64FVSetElement( fv2, 2, 0, 2.0, 0.0 );

    f64FVMul( fv1, fv2 );
    
    f64FVMulf64( fv1, 2.0 );
    
    TEST_ASSERT( f64FVEqual( fv1, fv2, EPS ) );
    
    f64FVFree( DefaultAllocator, &fv1 );
    f64FVFree( DefaultAllocator, &fv2 );
}
#endif

#if TEST
void test_fvar_div()
{
    u32 N = 3;
    
    f64FVar fv1 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv1, 0, 0, 3.0, 0.0 );
    f64FVSetElement( fv1, 1, 0, 2.5, 1.0 );
    f64FVSetElement( fv1, 2, 0, 0.5, 0.0 );
    
    f64FVar fv2 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv2, 0, 0, 2.0, 0.0 );
    f64FVSetElement( fv2, 1, 0, 2.0, 0.0 );
    f64FVSetElement( fv2, 2, 0, 2.0, 0.0 );
    
    f64FVDiv( fv1, fv2 );
    
    f64FVDivf64( fv1, 2.0 );
    
    TEST_ASSERT( f64FVEqual( fv1, fv2, EPS ) );
    
    f64FVFree( DefaultAllocator, &fv1 );
    f64FVFree( DefaultAllocator, &fv2 );
}
#endif

#if TEST
void test_fvar_neg()
{
    u32 N = 3;
    
    f64FVar fv1 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv1, 0, 0, 3.0, 0.0 );
    f64FVSetElement( fv1, 1, 0, 2.5, 1.0 );
    f64FVSetElement( fv1, 2, 0, 0.5, 0.0 );
    
    f64FVar fv2 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv2, 0, 0, -3.0,  0.0 );
    f64FVSetElement( fv2, 1, 0, -2.5, -1.0 );
    f64FVSetElement( fv2, 2, 0, -0.5,  0.0 );
    
    f64FVNeg( fv1 );
    
    TEST_ASSERT( f64FVEqual( fv1, fv2, EPS ) );
    
    f64FVFree( DefaultAllocator, &fv1 );
    f64FVFree( DefaultAllocator, &fv2 );
}
#endif


FVAR_EXP(f64, f64Mul, exp);
FVAR_LOG(f64, f64Div, log);

#if TEST
void test_fvar_exp()
{
    u32 N = 3;
    
    f64FVar fv1 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv1, 0, 0, 3.0, 0.0 );
    f64FVSetElement( fv1, 1, 0, 2.5, 1.0 );
    f64FVSetElement( fv1, 2, 0, 0.5, 0.0 );
    
    f64FVar fv2 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv2, 0, 0, exp(3.0),  0.0 );
    f64FVSetElement( fv2, 1, 0, exp(2.5),  exp(2.5) );
    f64FVSetElement( fv2, 2, 0, exp(0.5),  0.0 );
    
    f64FVExp( fv1 );
    
    TEST_ASSERT( f64FVEqual( fv1, fv2, EPS ) );
    
    
    f64FVFree( DefaultAllocator, &fv1 );
    f64FVFree( DefaultAllocator, &fv2 );
}
#endif

#if TEST
void test_fvar_log()
{
    u32 N = 3;
    
    f64FVar fv1 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv1, 0, 0, 3.0, 0.0 );
    f64FVSetElement( fv1, 1, 0, 2.5, 1.0 );
    f64FVSetElement( fv1, 2, 0, 0.5, 0.0 );
    
    f64FVar fv2 = f64FVMake( DefaultAllocator, N, 1 );
    
    f64FVSetElement( fv2, 0, 0, log(3.0),  0.0 );
    f64FVSetElement( fv2, 1, 0, log(2.5),  1 / 2.5 );
    f64FVSetElement( fv2, 2, 0, log(0.5),  0.0 );
    
    f64FVLog( fv1 );
    
    TEST_ASSERT( f64FVEqual( fv1, fv2, EPS ) );
    
    
    f64FVFree( DefaultAllocator, &fv1 );
    f64FVFree( DefaultAllocator, &fv2 );
}
#endif


FVAR_MATADD(f64, f64MatAdd);
FVAR_MATSUB(f64, f64MatSub);
FVAR_MATMUL(f64, f64MatMul, f64MatMulIP);


#if TEST
void test_fvar_matadd_matsub()
{
    u32 N = 2;
    
    f64FVar fv1 = f64FVMake( DefaultAllocator, N, N );
    f64FVar fv2 = f64FVMake( DefaultAllocator, N, N );
    f64FVar fv3 = f64FVMake( DefaultAllocator, N, N );
    
    
    for ( u32 i=0; i<N; ++i ) {
        for ( u32 j=0; j<N; ++j ) {
            f64FVSetElement( fv1, i, j, (f64) (i+0.5) / (j+1), 0.0 );
            f64FVSetElement( fv2, i, j, (f64) (j+0.5) / (i+1), 0.0 );
        }
    }
    
    f64FVMatAdd( fv1, fv2, fv3 );
    
    f64FVSetElement( fv1, 0, 0, 1.00, 0.0 );
    f64FVSetElement( fv1, 0, 1, 1.75, 0.0 );
    f64FVSetElement( fv1, 1, 0, 1.75, 0.0 );
    f64FVSetElement( fv1, 1, 1, 1.50, 0.0 );
    
    TEST_ASSERT( f64FVEqual( fv1, fv3, EPS ) );
    
    
    f64FVMatSub( fv3, fv2, fv1 );
    
    f64FVSetElement( fv2, 0, 0, 0.5,  0.0 );
    f64FVSetElement( fv2, 0, 1, 0.25, 0.0 );
    f64FVSetElement( fv2, 1, 0, 1.5,  0.0 );
    f64FVSetElement( fv2, 1, 1, 0.75, 0.0 );
    
    TEST_ASSERT( f64FVEqual( fv1, fv2, EPS ) );
    
    f64FVFree( DefaultAllocator, &fv1 );
    f64FVFree( DefaultAllocator, &fv2 );
    f64FVFree( DefaultAllocator, &fv3 );
}
#endif


#if TEST
void test_fvar_matmul()
{
    u32 N = 2;
    
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
    
    f64FVSetElement( fv1, 0, 0, 0.3125, 0.0 );
    f64FVSetElement( fv1, 0, 1, 0.9375, 0.0 );
    f64FVSetElement( fv1, 1, 0, 0.9375, 0.0 );
    f64FVSetElement( fv1, 1, 1, 2.8125, 0.0 );
    
    
    TEST_ASSERT( f64FVEqual( fv1, fv3, EPS ) );
    
    
    f64FVFree( DefaultAllocator, &fv1 );
    f64FVFree( DefaultAllocator, &fv2 );
    f64FVFree( DefaultAllocator, &fv3 );
}
#endif


#undef EPS



