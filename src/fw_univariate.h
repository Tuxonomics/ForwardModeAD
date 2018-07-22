
typedef struct FVar FVar;
struct FVar {
    f64 val;
    f64 dot;
};

FVar FVMake(f64 val, f64 dot)  /* make AD number */
{
    FVar fv;
    fv.val = val;
    fv.dot = dot;
    return fv;
}

void FVPrint(FVar x, const char* name) /* print AD number */
{
    printf("FVar (%s): {\n", name);
    printf("\t.val = %.4f\n", x.val);
    printf("\t.dot = %.4f\n", x.dot);
    printf("}\n");
}

#define FVPRINT(x) FVPrint(x, #x)
#define print_FVar(x) FVPrint(x, NULL)

f64 FVReal(FVar x) /* get value part */
{
    return x.val;
}

f64 FVDot(FVar x) /* get derivative part */
{
    return x.dot;
}

FVar FVAdd(FVar x, FVar y)  /* add two AD numbers */
{
    return (FVar) {
        .val = x.val + y.val,
        .dot = x.dot + y.dot
    };
}

FVar FVAddD(FVar x, f64 a)  /* add AD with double */
{
    return (FVar) {
        .val = x.val + a,
        .dot = x.dot
    };
}

FVar FVDAdd(f64 a, FVar x)  /* add AD with double */
{
    return (FVar) {
        .val = x.val + a,
        .dot = x.dot
    };
}

FVar FVSub(FVar x, FVar y)  /* subtract y from x */
{
    return (FVar) {
        .val = x.val - y.val,
        .dot = x.dot - y.dot
    };
}

FVar FVMul(FVar x, FVar y)  /* multiply two AD numbers */
{
    return (FVar) {
        .val = x.val * y.val,
        .dot = x.val * y.dot + x.dot * y.val
    };
}

FVar FVMulD(FVar x, f64 a)  /* multiply AD number times double */
{
    return (FVar) {
        .val = x.val * a,
        .dot = x.dot * a
    };
}

FVar FVDMul(f64 a, FVar x)  /* multiply double times AD number */
{
    return (FVar) {
        .val = x.val * a,
        .dot = x.dot * a
    };
}

FVar FVDiv(FVar x, FVar y)  /* divide AD by AD */
{
    return (FVar) {
        .val = x.val / y.val,
        .dot = (x.dot * y.val - x.val * y.dot) / (y.val * y.val)
    };
}

FVar FDDivD(FVar x, double a)  /* divide AD by double */
{
    return (FVar) {
        .val = x.val / a,
        .dot = x.dot / a
    };
}

FVar FVNeg(FVar x) /* negate AD */
{
    return (FVar) {
        .val = -x.val,
        .dot = -x.dot
    };
}


/* elementary forward AD functions */

FVar FVSqrt(FVar x) /* square root of a AD number */
{
    f64 tmp = sqrt(x.val);
    return (FVar) {
        .val = tmp,
        .dot = 0.5 * x.dot / tmp
    };
}

FVar FVPow(FVar x, f64 a) /* x^a, power of AD number */
{
    f64 tmp = pow(x.val, a - 1.0);
    return (FVar) {
        .val = tmp * x.val,
        .dot = a * x.dot * tmp
    };
}

FVar FVSin(FVar x) /* sine of AD number */
{
    return (FVar) {
        .val = sin(x.val),
        .dot = x.dot * cos(x.val)
    };
}

FVar FVCos(FVar x) /* cosine of AD number */
{
    return (FVar) {
        .val = cos(x.val),
        .dot = -x.dot * sin(x.val)
    };
}

FVar FVTan(FVar x) /* tangent of AD number */
{
    f64 tmp = cos(x.val);
    return (FVar) {
        .val = tan(x.val),
        .dot = x.dot / ( tmp * tmp )
    };
}

FVar FVAtan(FVar x) /* arctangent of AD number */
{
    return (FVar) {
        .val = atan(x.val),
        .dot = x.dot / ( 1.0 + (x.val * x.val) )
    };
}

FVar FVExp(FVar x) /* e^a, exp of AD number */
{
    f64 tmp = exp(x.val);
    return (FVar) {
        .val = tmp,
        .dot = x.dot * tmp
    };
}

FVar FVLog(FVar x) /* log base e of AD number */
{
    return (FVar) {
        .val = log(x.val),
        .dot = x.dot / x.val
    };
}

FVar FVLogAbs(FVar x) /* log base e of absolute value of AD number */
{
    return (FVar) {
        .val = log(fabs(x.val)),
        .dot = x.dot / x.val
    };
}

FVar FVSinh(FVar x) /* hyperbolic sine of AD number */
{
    return (FVar) {
        .val = sinh(x.val),
        .dot = x.dot * cosh(x.val)
    };
}

FVar FVCosh(FVar x) /* hyperbolic cosine of AD number */
{
    return (FVar) {
        .val = cosh(x.val),
        .dot = x.dot * sinh(x.val)
    };
}

FVar FVTanh(FVar x) /* hyperbolic tangent of AD number */
{
    f64 tmp = tanh(x.val);
    return (FVar) {
        .val = tmp,
        .dot = x.dot * ( 1.0 - (tmp * tmp) )
    };
}

FVar FVAtanh(FVar x) /* hyperbolic arctangent of AD number */
{
    return (FVar) {
        .val = atanh(x.val),
        .dot = x.dot / ( 1.0 - (x.val * x.val) )
    };
}


#ifndef MAT_DECL

#define MAT_DECL(type) typedef struct Mat_##type Mat_##type; \
    struct Mat_##type { \
        u32 dim0; \
        u32 dim1; \
        type *data; \
    };

#define MAT_MAKE(type) Mat_##type \
    MatMake_##type(Allocator al, u32 dim0, u32 dim1) \
        { \
            Mat_##type m; \
            m.dim0  = dim0; \
            m.dim1  = dim1; \
            m.data  = (type *) Alloc(al, dim0 * dim1 * sizeof(type)); \
            return m; \
        }

#define MAT_FREE(type) void \
    MatFree_##type(Allocator al, Mat_##type *m) \
        { \
            ASSERT(m->data); \
            Free(al, m->data); \
        }

#define MAT_PRINT(type) void \
    MatPrint_##type(Mat_##type m, const char* name) \
        { \
            printf("Mat (%s): {\n", name); \
            printf("\t.dim0 = %u\n", m.dim0); \
            printf("\t.dim1 = %u\n", m.dim1); \
            printf("\t.data = {\n"); \
            u32 dim; \
            for (u32 i=0; i<m.dim0; ++i) { \
                for (u32 j=0; j<m.dim1; ++j) { \
                    dim = i*m.dim1 + j; \
                    printf("[%d]\n", dim); \
                    print_##type(m.data[dim]); \
                    printf("\n"); \
                } \
                printf("\n"); \
            } \
            printf("\n\t}\n\n"); \
            printf("}\n"); \
        }

#define MATPRINT(type, x) MatPrint_##type(x, #x)

#define MAT_ZERO(type) Mat_##type \
    MatZeroMake_##type(Allocator al, u32 dim0, u32 dim1) \
        { \
            Mat_##type m = MatMake_##type(al, dim0, dim1); \
            memset(m.data, 0, dim0 * dim1 * sizeof(type)); \
            return m; \
        } \

#define MAT_SETELEMENT(type) Inline void \
    MatSetElement_##type(Mat_##type m, u32 dim0, u32 dim1, type val) \
        { \
            ASSERT( m.data ); \
            ASSERT( dim0 <= m.dim0 ); \
            ASSERT( dim1 <= m.dim1 ); \
            u32 dim = dim0*m.dim1 + dim1; \
            m.data[dim] = val; \
        }

#define MAT_GETELEMENT(type) Inline type \
    MatGetElement_##type(Mat_##type m, u32 dim0, u32 dim1) \
        { \
            ASSERT( dim0 <= m.dim0 ); \
            ASSERT( dim1 <= m.dim1 ); \
            u32 dim = dim0 * m.dim1 + dim1; \
            return m.data[dim]; \
        }

#define MAT_ADD(type, fun) Inline void \
    MatAdd_##type(Mat_##type a, Mat_##type b, Mat_##type c) /* c = a + b */ \
    { \
        ASSERT(a.dim0 == b.dim0 && a.dim1 == b.dim1 && a.dim0 == c.dim0 && a.dim1 == c.dim1); \
     \
        for ( u32 i=0; i<(a.dim0*a.dim1); ++i ) { \
            c.data[i] = fun( a.data[i], b.data[i] ); \
        } \
    }

#define MAT_SUB(type, fun) Inline void \
    MatSub_##type(Mat_##type a, Mat_##type b, Mat_##type c) /* c = a - b */ \
        { \
            ASSERT(a.dim0 == b.dim0 && a.dim1 == b.dim1 && a.dim0 == c.dim0 && a.dim1 == c.dim1); \
         \
            for ( u32 i=0; i<(a.dim0*a.dim1); ++i ) { \
                c.data[i] = fun( a.data[i], b.data[i] ); \
            } \
        }


// NOTE(jonas): naive implementation, change when necessary
#define MAT_MUL(type, addFun, mulFun) Inline void \
MatMul_##type(Mat_##type a, Mat_##type b, Mat_##type c) \
{ \
    ASSERT(a.dim0 == c.dim0 && a.dim1 == b.dim0 && b.dim1 == c.dim1); \
 \
    /* type val; */ \
 \
    for ( u32 i = 0; i < c.dim0; ++i ) { \
        for ( u32 j = 0; j < c.dim1; ++j ) { \
            type val = {0}; \
            for ( u32 k = 0; k < a.dim1; ++k ) { \
                val = addFun( \
                    val, \
                    mulFun( MatGetElement_##type( a, i, k ), MatGetElement_##type( b, k, j ) \
                    ) \
                ); \
            } \
            MatSetElement_##type( c, i, j, val ); \
        } \
    } \
}

#endif



MAT_DECL(FVar);
MAT_MAKE(FVar);
MAT_FREE(FVar);
MAT_PRINT(FVar);
MAT_ZERO(FVar);
MAT_SETELEMENT(FVar);
MAT_GETELEMENT(FVar);
MAT_ADD(FVar, FVAdd);
MAT_SUB(FVar, FVSub);
MAT_MUL(FVar, FVAdd, FVMul);




void FVarGradient( FVar f( Mat_FVar ), Mat_FVar input, Mat_FVar grad )
{
    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
    
    FVar tmp;
    u32 N = input.dim0;
    
    Mat_FVar xCpy = MatMake_FVar( DefaultAllocator, N, 1 );
    
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
}





