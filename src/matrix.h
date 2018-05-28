
#include "BLAS/cblas_f77.h"
#include "BLAS/cblas.h"

#define ORDER  enum CBLAS_ORDER
#define ROWMAJ CblasRowMajor
#define COLMAJ CblasColMajor

#define DEFAULT_MAJ ROWMAJ

typedef struct Mat Mat;
struct Mat {
    ORDER order;
    
    u32 dim0;
    u32 dim1;
    
    f64 *data;
};


Mat MatMake(Allocator al, u32 dim0, u32 dim1)
{
    Mat m;
    m.order = DEFAULT_MAJ;
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

void MatPrint(Mat m, const char* name = NULL)
{
    printf("Mat (%s): {\n", name);
    printf("\t.order = %s\n\n", (m.order == ROWMAJ) ? "RowMajor" : "ColMajor");
    printf("\t.dim0 = %u\n", m.dim0);
    printf("\t.dim1 = %u\n", m.dim1);
    printf("\t.data = {\n");
    u32 dim;
    for (u32 i=0; i<m.dim0; ++i) {
        for (u32 j=0; j<m.dim1; ++j) {
            dim = i*m.dim1 + j;
            printf("\t%.4f", dim, m.data[dim]);
        }
        printf("\n");
    }
    printf("\n\t}\n\n");
    printf("}\n");
}

Mat MatZeroMake(Allocator al, u32 dim0, u32 dim1)
{
    Mat m = MatMake(al, dim0, dim1);
    memset(m.data, 0, dim0 * dim1 * sizeof(f64));
    return m;
}

void SetElement(Mat m, u32 dim0, u32 dim1, f64 val)
{
    ASSERT(m.data);
    ASSERT(dim0 >= 0 && dim0 <= m.dim0 );
    ASSERT(dim1 >= 0 && dim1 <= m.dim1 );
    u32 dim = dim0*m.dim1 + dim1;
    m.data[dim] = val;
}

Mat MatMul(Mat a, Mat b)
{
    Mat c = MatMake(DefaultAllocator, a.dim0, b.dim1);
    
    cblas_dgemm(
        DEFAULT_MAJ, CblasNoTrans, CblasNoTrans,
        a.dim0, b.dim1, a.dim1, 1.0,
        a.data, a.dim1, b.data, b.dim0,
        0.0, c.data, b.dim1);

    return c;
}



