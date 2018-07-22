
#include "BLAS/cblas_f77.h"
#include "BLAS/cblas.h"

#define ORDER  enum CBLAS_ORDER
#define ROWMAJ CblasRowMajor
#define COLMAJ CblasColMajor

#define DEFAULT_MAJ ROWMAJ


// NOTE(jonas): all row major for now

typedef struct Mat Mat;
struct Mat {
//    ORDER order;
    
    u32 dim0;
    u32 dim1;
    
    f64 *data;
};


// row-major
// u32 idx = i * nCols + j;
// now mat[idx] corresponds to m(i, j)


Mat MatMake(Allocator al, u32 dim0, u32 dim1)
{
    Mat m;
//    m.order = DEFAULT_MAJ;
    m.dim0  = dim0;
    m.dim1  = dim1;
    m.data  = (f64 *) Alloc(al, dim0 * dim1 * sizeof(f64));
    return m;
}

void MatFree(Allocator al, Mat *m)
{
    ASSERT(m->data);
    Free(al, m->data);
}

void MatPrint(Mat m, const char* name)
{
    printf("Mat (%s): {\n", name);
//    printf("\t.order = %s\n\n", (m.order == ROWMAJ) ? "RowMajor" : "ColMajor");
    printf("\t.dim0 = %u\n", m.dim0);
    printf("\t.dim1 = %u\n", m.dim1);
    printf("\t.data = {\n");
    u32 dim;
    for (u32 i=0; i<m.dim0; ++i) {
        for (u32 j=0; j<m.dim1; ++j) {
            dim = i*m.dim1 + j;
            printf("\t%.4f", m.data[dim]);
        }
        printf("\n");
    }
    printf("\n\t}\n\n");
    printf("}\n");
}

#define MATPRINT(x) MatPrint(x, #x)

Mat MatZeroMake(Allocator al, u32 dim0, u32 dim1)
{
    Mat m = MatMake(al, dim0, dim1);
    memset(m.data, 0, dim0 * dim1 * sizeof(f64));
    return m;
}

Inline
void MatSetElement(Mat m, u32 dim0, u32 dim1, f64 val)
{
    ASSERT( m.data );
    ASSERT( dim0 <= m.dim0 );
    ASSERT( dim1 <= m.dim1 );
    u32 dim = dim0*m.dim1 + dim1;
    m.data[dim] = val;
}

Inline
f64 MatGetElement(Mat m, u32 dim0, u32 dim1)
{
    ASSERT( dim0 <= m.dim0 );
    ASSERT( dim1 <= m.dim1 );
    u32 dim = dim0 * m.dim1 + dim1;
    return m.data[dim];
}


// TODO(jonas): make addition and subtraction cache aware
Inline
void MatAdd(Mat a, Mat b, Mat c) /* c = a + b */
{
    ASSERT(a.dim0 == b.dim0 && a.dim1 == b.dim1 && a.dim0 == c.dim0 && a.dim1 == c.dim1);
    
    for ( u32 i=0; i<(a.dim0*a.dim1); ++i ) {
        c.data[i] = a.data[i] + b.data[i];
    }
}

Inline
void MatSub(Mat a, Mat b, Mat c) /* c = a - b */
{
    ASSERT(a.dim0 == b.dim0 && a.dim1 == b.dim1 && a.dim0 == c.dim0 && a.dim1 == c.dim1);
    
    for ( u32 i=0; i<(a.dim0*a.dim1); ++i ) {
        c.data[i] = a.data[i] - b.data[i];
    }
}


// NOTE(jonas): naive implementation, change when necessary
Inline
void MatMul(Mat a, Mat b, Mat c)
{
    ASSERT(a.dim0 == c.dim0 && a.dim1 == b.dim0 && b.dim1 == c.dim1);
    
    f64 val;
    
    for ( u32 i = 0; i < c.dim0; ++i ) {
        for ( u32 j = 0; j < c.dim1; ++j ) {
            val = 0;
            for ( u32 k = 0; k < a.dim1; ++k ) {
                val += MatGetElement( a, i, k ) * MatGetElement( b, k, j );
            }
            MatSetElement( c, i, j, val );
        }
    }
}


//#define CACHE_LINE_SIZE 64
//#define TILE_SIZE ( CACHE_LINE_SIZE / sizeof(f64) )
//void MatMul1(Mat a, Mat b, Mat c)
//{
//    ASSERT(a.dim0 == c.dim0 && a.dim1 == b.dim0 && b.dim1 == c.dim1);
//    u32 N = a.dim0, M = a.dim1, P = b.dim1;
//    u32 i, j, k, i2, j2, k2;
//    f64 tmp;
//
//    // A(i,k), B(k, j), C(i, j)
//
//    f64 *restrict rC;
//    f64 *restrict rA;
//    f64 *restrict rB;
//
//    for ( i=0; i<N; i+=TILE_SIZE ) {
//        for ( j=0; j<P; j+=TILE_SIZE ) {
//            for ( k=0; k<P; k+= TILE_SIZE) {
//
//                for ( i2=i; i2<(i+TILE_SIZE) && i<N; ++i2 ) {
//                    for ( j2=j; j2<(j+TILE_SIZE) && j2<P; ++j2 ) {
//                        tmp = 0;
//                        for ( k2=k; k2<(k+TILE_SIZE) && k2<M; ++k2 ) {
//                            tmp += a.data[i2 * M + k2] * b.data[k2 * P + j2];
//                        }
//                        c.data[i2 * P + j2] = tmp;
//                    }
//                }
//
//            }
//        }
//    }
//}


void MatMul_BLAS(Mat a, Mat b, Mat c)
{
    ASSERT(a.dim0 == c.dim0 && a.dim1 == b.dim0 && b.dim1 == c.dim1);
    
    cblas_dgemm(
        DEFAULT_MAJ, CblasNoTrans, CblasNoTrans,
        a.dim0, b.dim1, a.dim1, 1.0,
        a.data, a.dim1, b.data, b.dim1,
        0.0, c.data, b.dim1);
}


