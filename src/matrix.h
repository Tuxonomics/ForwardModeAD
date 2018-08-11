
#ifndef MAT_DECL

#define MAT_DECL(type) typedef struct type##Mat type##Mat; \
    struct type##Mat { \
            u32 dim0; \
            u32 dim1; \
            type *data; \
        };

#define MAT_MAKE(type) type##Mat \
    type##MatMake(Allocator al, u32 dim0, u32 dim1) \
        { \
            type##Mat m; \
            m.dim0  = dim0; \
            m.dim1  = dim1; \
            m.data  = (type *) Alloc(al, dim0 * dim1 * sizeof(type)); \
            return m; \
        }

#define MAT_FREE(type) void \
    type##MatFree(Allocator al, type##Mat *m) \
        { \
            ASSERT(m->data); \
            Free(al, m->data); \
        }

#define MAT_PRINT(type) void \
    type##MatPrint(type##Mat m, const char* name) \
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

#define MATPRINT(type, x) type##MatPrint(x, #x)

#define MAT_ZERO(type) type##Mat \
    type##MatZeroMake(Allocator al, u32 dim0, u32 dim1) \
        { \
            type##Mat m = type##MatMake(al, dim0, dim1); \
            memset(m.data, 0, dim0 * dim1 * sizeof(type)); \
            return m; \
        } \

#define MAT_SETELEMENT(type) Inline void \
    type##MatSetElement(type##Mat m, u32 dim0, u32 dim1, type val) \
        { \
            ASSERT( m.data ); \
            ASSERT( dim0 <= m.dim0 ); \
            ASSERT( dim1 <= m.dim1 ); \
            u32 dim = dim0*m.dim1 + dim1; \
            m.data[dim] = val; \
        }

#define MAT_GETELEMENT(type) Inline type \
    type##MatGetElement(type##Mat m, u32 dim0, u32 dim1) \
        { \
            ASSERT( dim0 <= m.dim0 ); \
            ASSERT( dim1 <= m.dim1 ); \
            u32 dim = dim0 * m.dim1 + dim1; \
            return m.data[dim]; \
        }

#define MAT_SETCOL(type) void \
    type##MatSetCol(type##Mat m, u32 dim, type##Mat newCol) \
        { \
            ASSERT( m.data ); \
            ASSERT( newCol.dim0 * newCol.dim1 == m.dim0 ); \
             \
            u32 idx = dim*m.dim1; \
            memcpy( &m.data[idx], newCol.data, sizeof( type ) * m.dim0 ); \
        }

#define MAT_GETCOL(type) void \
    type##MatGetCol(type##Mat m, u32 dim, type##Mat newCol) \
        { \
            ASSERT( m.data ); \
            ASSERT( newCol.dim0 * newCol.dim1 == m.dim0 ); \
            \
            u32 idx = dim*m.dim1; \
            memcpy( newCol.data, &m.data[idx], sizeof( type ) * m.dim0 ); \
        }

#define MAT_ADD(type, fun) Inline void \
    type##MatAdd(type##Mat a, type##Mat b, type##Mat c) /* c = a + b */ \
        { \
            ASSERT(a.dim0 == b.dim0 && a.dim1 == b.dim1 && a.dim0 == c.dim0 && a.dim1 == c.dim1); \
            \
            for ( u32 i=0; i<(a.dim0*a.dim1); ++i ) { \
                c.data[i] = fun( a.data[i], b.data[i] ); \
            } \
        }

#define MAT_SUB(type, fun) Inline void \
    type##MatSub(type##Mat a, type##Mat b, type##Mat c) /* c = a - b */ \
        { \
            ASSERT(a.dim0 == b.dim0 && a.dim1 == b.dim1 && a.dim0 == c.dim0 && a.dim1 == c.dim1); \
            \
            for ( u32 i=0; i<(a.dim0*a.dim1); ++i ) { \
                c.data[i] = fun( a.data[i], b.data[i] ); \
            } \
        }


// NOTE(jonas): naive implementation, change when necessary
#define MAT_MUL(type, addFun, mulFun) Inline void \
    type##MatMul(type##Mat a, type##Mat b, type##Mat c) \
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
                            mulFun( type##MatGetElement( a, i, k ), type##MatGetElement( b, k, j ) ) \
                        ); \
                    } \
                    type##MatSetElement( c, i, j, val ); \
                } \
            } \
        }

#endif
