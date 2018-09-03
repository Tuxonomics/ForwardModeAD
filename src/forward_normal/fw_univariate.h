#include "fw_matrix.h"

#ifndef FVAR_DECL

#define FVAR_DECL(type) typedef struct type##FVar type##FVar; \
    struct type##FVar { \
        type val; \
        type dot; \
    }

#define FVAR_MAKE(type) type##FVar \
    type##FVMake(type val, type dot) \
    { \
        type##FVar fv;\
        fv.val = val; \
        fv.dot = dot; \
        return fv; \
    }

#define FVAR_CONST(type) type##FVar \
    type##FVConst(type val) \
    { \
        type##FVar fv; \
        fv.val = val; \
        memset(&fv.dot, 0, sizeof(type)); \
        return fv; \
    }

#define FVAR_PRINT(type, printFun) void \
    type##FVPrint(type##FVar x, const char* name) \
    { \
        printf("FVar (%s): {\n", name); \
        printf("\t.val = "); \
        printFun(x.val); \
        printf("\n"); \
        printf("\t.dot = "); \
        printFun(x.dot); \
        printf("\n"); \
        printf("}\n"); \
    }

/* numerical equality of two AD numbers */
#define FVAR_EQUAL(type, baseFun) b32 \
    type##FVEqual( type##FVar x, type##FVar y, f64 eps ) \
    { \
        return baseFun( x.val, y.val, eps ) && baseFun( x.dot, y.dot, eps ); \
    }

/* add two AD numbers */
#define FVAR_ADD(type, addFun) type##FVar \
    type##FVAdd(type##FVar x, type##FVar y) \
    { \
        return type##FVMake( \
            addFun( x.val, y.val ), \
            addFun( x.dot, y.dot ) \
        ); \
    }

/* add AD with double */
#define FVAR_ADD_TYPE(type, addFun) type##FVar \
    type##FVAdd##type(type##FVar x, type a) \
    { \
        return type##FVMake( addFun( x.val, a ), x.dot ); \
    } \

/* add AD with double */
#define FVAR_TYPE_ADD(type, addFun) type##FVar \
    type##FV##type##Add(type a, type##FVar x) \
    { \
        return type##FVMake( addFun( x.val, a ), x.dot ); \
    }

/* subtract y from x */
#define FVAR_SUB(type, subFun) type##FVar \
    type##FVSub(type##FVar x, type##FVar y) \
    { \
        return type##FVMake( \
            subFun( x.val, y.val), \
            subFun( x.dot, y.dot) \
        ); \
    }

/* multiply two AD numbers */
#define FVAR_MUL(type, mulFun, addFun) type##FVar \
    type##FVMul(type##FVar x, type##FVar y) \
    { \
        return type##FVMake( \
            mulFun( x.val, y.val ), \
            addFun( mulFun( x.val, y.dot), mulFun( x.dot, y.val) ) \
        ); \
    }

/* multiply AD number times double */
#define FVAR_MUL_TYPE(type, mulFun) type##FVar \
    type##FVMul##type(type##FVar x, type a) \
    { \
        return type##FVMake( \
            mulFun( x.val, a ), \
            mulFun( x.dot, a ) \
        ); \
    }

/* multiply double times AD number */
#define FVAR_TYPE_MUL(type, mulFun) type##FVar \
    type##FV##type##Mul(type a, type##FVar x) \
    { \
        return type##FVMake( \
            mulFun( x.val, a ), \
            mulFun( x.dot, a ) \
        ); \
    }

/* divide AD by AD */
#define FVAR_DIV(type, subFun, mulFun, divFun) type##FVar \
type##FVDiv(type##FVar x, type##FVar y) \
{ \
    return type##FVMake( \
        divFun(x.val, y.val), \
        divFun(subFun(mulFun(x.dot, y.val), mulFun(x.val, y.dot)), mulFun(y.val, y.val) ) \
    ); \
}

/* divide AD by f64 */
#define FVAR_DIV_TYPE(type, divFun) type##FVar \
    type##FVDiv##type(type##FVar x, type a) \
    { \
        return type##FVMake( \
            divFun( x.val, a ), \
            divFun( x.dot, a ) \
        ); \
    }

/* divide f64 by AD */
#define FVAR_TYPE_DIV(type, mulFun, divFun) type##FVar \
    type##FV##typeDiv(type a, type##FVar x) \
    { \
        return type##FVMake( \
            divFun( a, x.val ), \
            divFun( a, mulFun( x.val, x.val ) ) \
        ); \
    }

/* negate AD */
#define FVAR_NEG(type, negFun) type##FVar \
    type##FVNeg(type##FVar x) \
    { \
        return type##FVMake( \
            negFun( x.val ), \
            negFun( x.dot ) \
        ); \
    }


/* elementary forward AD functions */


/* square root of a AD number */
#define FVAR_SQRT(type, constFun, mulFun, divFun, sqrtFun) type##FVar \
    type##FVSqrt(type##FVar x) \
    { \
        type c = constFun( 0.5 ); \
        type tmp = sqrtFun(x.val); \
        return type##FVMake( \
            tmp, \
            divFun( mulFun( c, x.dot ), tmp ) \
        ); \
    }

/* x^a, power of AD number */
#define FVAR_POW(type, constFun, mulFun, powFun) type##FVar \
    type##FVPow(type##FVar x, f64 a) \
    { \
        type c = constFun( a ); \
        type tmp = powFun(x.val, a - 1.0); \
        return type##FVMake( \
            mulFun( tmp, x.val ), \
            mulFun( c, mulFun( x.dot, tmp ) ) \
        ); \
    }

/* sine of AD number */
#define FVAR_SIN(type, mulFun, sinFun, cosFun) type##FVar \
    type##FVSin(type##FVar x) \
    { \
        return type##FVMake( \
            sinFun(x.val), \
            mulFun( x.dot, cosFun(x.val) ) \
        ); \
    }

/* cosine of AD number */
#define FVAR_COS(type, negFun, mulFun, sinFun, cosFun) type##FVar \
    type##FVCos(type##FVar x) \
    { \
        return type##FVMake( \
            cosFun( x.val ), \
            mulFun( negFun(x.dot), sinFun(x.val) ) \
        ); \
    }

/* tangent of AD number */
#define FVAR_TAN(type, mulFun, divFun, cosFun, tanFun) type##FVar \
    type##FVTan(type##FVar x) \
    { \
        type tmp = cosFun(x.val); \
        return type##FVMake( \
            tanFun( x.val ), \
            divFun( x.dot, mulFun(tmp, tmp) ) \
        ); \
    }

/* arctangent of AD number */
#define FVAR_ATAN(type, constFun, addFun, mulFun, divFun, atanFun) type##FVar \
    type##FVAtan(type##FVar x) \
    { \
        type c = constFun( 1.0 ); \
        return type##FVMake( \
            atanFun( x.val ), \
            divFun( x.dot, addFun(c, mulFun(x.val, x.val)) ) \
        ); \
    }

/* e^a, exp of AD number */
#define FVAR_EXP(type, mulFun, expFun) type##FVar \
    type##FVExp(type##FVar x) \
    { \
        type tmp = expFun(x.val); \
        return type##FVMake( \
            tmp, \
            mulFun( x.dot, tmp ) \
        ); \
    }

/* log base e of AD number */
#define FVAR_LOG(type, divFun, logFun) type##FVar \
    type##FVLog(type##FVar x) \
    { \
        return type##FVMake( \
            logFun( x.val ), \
            divFun( x.dot, x.val ) \
        ); \
    }

/* log base e of absolute value of AD number */
#define FVAR_LOGABS(type, divFun, absFun, logFun) type##FVar \
    type##FVLogAbs(type##FVar x) \
    { \
        return type##FVMake( \
            logFun( absFun(x.val) ), \
            divFun( x.dot, x.val ) \
        ); \
    }

/* hyperbolic sine of AD number */
#define FVAR_SINH(type, mulFun, sinhFun, coshFun) type##FVar \
    type##FVSinh(type##FVar x) \
    { \
        return type##FVMake( \
            sinhFun( x.val ), \
            mulFun( x.dot, coshFun(x.val) ) \
        ); \
    }

/* hyperbolic cosine of AD number */
#define FVAR_COSH(type, mulFun, sinhFun, coshFun) type##FVar \
    type##FVCosh(type##FVar x) \
    { \
        return type##FVMake( \
            coshFun(x.val), \
            mulFun( x.dot, sinhFun(x.val) ) \
        ); \
    }

/* hyperbolic tangent of AD number */
#define FVAR_TANH(type, constFun, subFun, mulFun, tanhFun) type##FVar \
    type##FVTanh(type##FVar x) \
    { \
        type c = constFun( 1.0 ); \
        type tmp = tanhFun( x.val ); \
        return type##FVMake( \
            tmp, \
            mulFun( x.dot, subFun( c, mulFun(tmp, tmp) ) ) \
        ); \
    }

/* hyperbolic arctangent of AD number */
#define FVAR_ATANH(type, constFun, subFun, divFun, mulFun, atanhFun) type##FVar \
    type##FVAtanh(type##FVar x) \
    { \
        type c = constFun( 1.0 ); \
        return type##FVMake( \
            atanhFun(x.val), \
            divFun( c, subFun( c, mulFun(x.val, x.val) ) )\
        ); \
    }

#endif

/* first-order derivatives */

FVAR_DECL(f64);
FVAR_MAKE(f64);
FVAR_CONST(f64);
FVAR_PRINT(f64, f64Print);

#define f64FVPRINT0(x) f64FVPrint( x, NULL );

FVAR_EQUAL(f64, f64Equal);
FVAR_ADD(f64, f64Add);
FVAR_ADD_TYPE(f64, f64Add);
FVAR_TYPE_ADD(f64, f64Add);
FVAR_SUB(f64, f64Sub);
FVAR_MUL(f64, f64Mul, f64Add);
FVAR_MUL_TYPE(f64, f64Mul);
FVAR_TYPE_MUL(f64, f64Mul);
FVAR_DIV(f64, f64Sub, f64Mul, f64Div);
FVAR_DIV_TYPE(f64, f64Div);
FVAR_TYPE_DIV(f64, f64Mul, f64Div);
FVAR_NEG(f64, f64Neg);

#if TEST
void test_fvar_basic_functions()
{
#define EPS 1E-10
    
    f64FVar a = f64FVMake( 1.1, 1.0 );
    f64FVar b = f64FVConst( 2.0 );
    f64FVar c = a;

    TEST_ASSERT( f64FVEqual( a, c, EPS ) );
    TEST_ASSERT( ! f64FVEqual( a, b, EPS ) );

    TEST_ASSERT( f64FVEqual( f64FVAdd( a, b ), f64FVAdd( b, a ), EPS ) );

    TEST_ASSERT( f64FVEqual( f64FVAddf64( a, 5 ), f64FVf64Add( 5, a ), EPS ) );

    TEST_ASSERT( f64FVEqual( f64FVSub( a, b ), f64FVNeg( f64FVSub( b, a ) ), EPS ) );

    TEST_ASSERT( f64FVEqual( f64FVMul( a, b ), f64FVMul( b, a ), EPS ) );

    TEST_ASSERT( f64FVEqual( f64FVMulf64( a, 5 ), f64FVf64Mul( 5, a ), EPS ) );

    TEST_ASSERT( f64FVEqual( f64FVDiv( a, b ), f64FVMake( 0.55, 0.5 ), EPS ) );

    TEST_ASSERT( f64FVEqual( f64FVDivf64( a, 2 ), f64FVMake( 0.55, 0.5 ), EPS ) );

#undef EPS
}
#endif


FVAR_SQRT(f64, f64Const, f64Mul, f64Div, sqrt);
FVAR_POW(f64, f64Const, f64Mul, pow);
FVAR_SIN(f64, f64Mul, sin, cos);
FVAR_COS(f64, f64Neg, f64Mul, sin, cos);
FVAR_TAN(f64, f64Mul, f64Div, cos, tan);
FVAR_ATAN(f64, f64Const, f64Add, f64Mul, f64Div, atan);
FVAR_EXP(f64, f64Mul, exp);
FVAR_LOG(f64, f64Div, log);
FVAR_LOGABS(f64, f64Div, fabs, log);
FVAR_SINH(f64, f64Mul, sinh, cosh);
FVAR_COSH(f64, f64Mul, sinh, cosh);
FVAR_TANH(f64, f64Const, f64Sub, f64Mul, tanh);
FVAR_ATANH(f64, f64Const, f64Sub, f64Div, f64Mul, atanh);


#if TEST
void test_fv_elementary_functions()
{
#define EPS 1E-10
    
    f64FVar a = f64FVMake( 2.0, 1.0 );
    f64FVar b = f64FVMake( 0.5, 1.0 );
    
    
    TEST_ASSERT(
        f64FVEqual(
            f64FVSqrt(a),
            f64FVMake( sqrt(a.val), 1 / (2 * sqrt(a.val)) ),
            EPS
        )
    );
    
    TEST_ASSERT(
        f64FVEqual(
            f64FVPow( a, 2 ),
            f64FVMake( pow(a.val, 2.0), 2 * a.val ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVSin( a ),
            f64FVMake( sin(a.val), cos(a.val) ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVCos( a ),
            f64FVMake( cos(a.val), -sin(a.val) ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVTan( a ),
            f64FVDiv( f64FVSin(a), f64FVCos(a) ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVAtan( a ),
            f64FVMake( atan(a.val), 1 / ( 1 + (a.val * a.val) ) ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVExp( a ),
            f64FVMake( exp(a.val), exp(a.val) ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVLog( a ),
            f64FVMake( log(a.val), 1 / a.val ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVLogAbs( f64FVNeg(a) ),
            f64FVMake( log( a.val ), 1 / a.val ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVSinh( a ),
            f64FVMake( sinh( a.val ), cosh( a.val ) ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVCosh( a ),
            f64FVMake( cosh( a.val ), sinh( a.val ) ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVTanh( a ),
            f64FVMake( tanh( a.val ), 1 - tanh( a.val ) * tanh( a.val ) ),
            EPS
        )
    );

    TEST_ASSERT(
        f64FVEqual(
            f64FVAtanh( b ),
            f64FVMake( atanh( b.val ), 1.0 / ( 1.0 - (b.val * b.val) ) ),
            EPS
        )
    );

    
    
#undef EPS
}
#endif


MAT_DECL(f64FVar);
MAT_MAKE(f64FVar);
MAT_FREE(f64FVar);
MAT_PRINT(f64FVar, f64FVPRINT0);
MAT_EQUAL(f64FVar, f64FVEqual);
MAT_ZERO(f64FVar);
MAT_SETELEMENT(f64FVar);
MAT_GETELEMENT(f64FVar);
MAT_SETCOL(f64FVar);
MAT_GETCOL(f64FVar);
MAT_ADD(f64FVar, f64FVAdd);
MAT_SUB(f64FVar, f64FVSub);
MAT_MUL(f64FVar, f64FVAdd, f64FVMul);



/* second-order derivative */

FVAR_DECL(f64FVar);
FVAR_MAKE(f64FVar);
FVAR_CONST(f64FVar);
FVAR_PRINT(f64FVar, f64FVPRINT0);

#define f64FVarFVPRINT0(x) f64FVarFVPrint( x, NULL );

FVAR_EQUAL(f64FVar, f64FVEqual);
FVAR_ADD(f64FVar, f64FVAdd);
FVAR_ADD_TYPE(f64FVar, f64FVAdd);
FVAR_TYPE_ADD(f64FVar, f64FVAdd);
FVAR_SUB(f64FVar, f64FVSub);
FVAR_MUL(f64FVar, f64FVMul, f64FVAdd);
FVAR_MUL_TYPE(f64FVar, f64FVMul);
FVAR_TYPE_MUL(f64FVar, f64FVMul);
FVAR_DIV(f64FVar, f64FVSub, f64FVMul, f64FVDiv);
FVAR_DIV_TYPE(f64FVar, f64FVDiv);
FVAR_TYPE_DIV(f64FVar, f64FVMul, f64FVDiv);
FVAR_NEG(f64FVar, f64FVNeg);

FVAR_SQRT(f64FVar, f64FVConst, f64FVMul, f64FVDiv, f64FVSqrt);
FVAR_POW(f64FVar, f64FVConst, f64FVMul, f64FVPow);
FVAR_SIN(f64FVar, f64FVMul, f64FVSin, f64FVCos);
FVAR_COS(f64FVar, f64FVNeg, f64FVMul, f64FVSin, f64FVCos);
FVAR_TAN(f64FVar, f64FVMul, f64FVDiv, f64FVCos, f64FVTan);
FVAR_ATAN(f64FVar, f64FVConst, f64FVAdd, f64FVMul, f64FVDiv, f64FVAtan);
FVAR_EXP(f64FVar, f64FVMul, f64FVExp);
FVAR_LOG(f64FVar, f64FVDiv, f64FVLog);
//FVAR_LOGABS(f64, f64FVDiv, f64FVAbs, f64FVLog);
FVAR_SINH(f64FVar, f64FVMul, f64FVSinh, f64FVCosh);
FVAR_COSH(f64FVar, f64FVMul, f64FVSinh, f64FVCosh);
FVAR_TANH(f64FVar, f64FVConst, f64FVSub, f64FVMul, f64FVTanh);
FVAR_ATANH(f64FVar, f64FVConst, f64FVSub, f64FVDiv, f64FVMul, f64FVAtanh);


MAT_DECL(f64FVarFVar);
MAT_MAKE(f64FVarFVar);
MAT_FREE(f64FVarFVar);
MAT_PRINT(f64FVarFVar, f64FVarFVPRINT0);
MAT_EQUAL(f64FVarFVar, f64FVarFVEqual);
MAT_ZERO(f64FVarFVar);
MAT_SETELEMENT(f64FVarFVar);
MAT_GETELEMENT(f64FVarFVar);
MAT_SETCOL(f64FVarFVar);
MAT_GETCOL(f64FVarFVar);
MAT_ADD(f64FVarFVar, f64FVarFVAdd);
MAT_SUB(f64FVarFVar, f64FVarFVSub);
MAT_MUL(f64FVarFVar, f64FVarFVAdd, f64FVarFVMul);



/* third-order derivative */

//FVAR_DECL(f64FVarFVar);
//FVAR_MAKE(f64FVarFVar);
//FVAR_CONST(f64FVarFVar);
//FVAR_PRINT(f64FVarFVar, f64FVarFVPRINT0);
//
//#define f64FVarFVarFVPRINT0(x) f64FVarFVarFVPrint( x, NULL );
//
//FVAR_EQUAL(f64FVarFVar, f64FVarFVEqual);
//FVAR_ADD(f64FVarFVar, f64FVarFVAdd);
//FVAR_ADD_TYPE(f64FVarFVar, f64FVarFVAdd);
//FVAR_TYPE_ADD(f64FVarFVar, f64FVarFVAdd);
//FVAR_SUB(f64FVarFVar, f64FVarFVSub);
//FVAR_MUL(f64FVarFVar, f64FVarFVMul, f64FVarFVAdd);
//FVAR_MUL_TYPE(f64FVarFVar, f64FVarFVMul);
//FVAR_TYPE_MUL(f64FVarFVar, f64FVarFVMul);
//FVAR_DIV(f64FVarFVar, f64FVarFVSub, f64FVarFVMul, f64FVarFVDiv);
//FVAR_DIV_TYPE(f64FVarFVar, f64FVarFVDiv);
//FVAR_TYPE_DIV(f64FVarFVar, f64FVarFVMul, f64FVarFVDiv);
//FVAR_NEG(f64FVarFVar, f64FVarFVNeg);

//FVAR_SQRT(f64FVar, f64FVConst, f64FVMul, f64FVDiv, f64FVSqrt);
//FVAR_POW(f64FVar, f64FVConst, f64FVMul, f64FVPow);
//FVAR_SIN(f64FVar, f64FVMul, f64FVSin, f64FVCos);
//FVAR_COS(f64FVar, f64FVNeg, f64FVMul, f64FVSin, f64FVCos);
//FVAR_TAN(f64FVar, f64FVMul, f64FVDiv, f64FVCos, f64FVTan);
//FVAR_ATAN(f64FVar, f64FVConst, f64FVAdd, f64FVMul, f64FVDiv, f64FVAtan);
//FVAR_EXP(f64FVar, f64FVMul, f64FVExp);
//FVAR_LOG(f64FVar, f64FVDiv, f64FVLog);
////FVAR_LOGABS(f64, f64FVDiv, f64FVAbs, f64FVLog);
//FVAR_SINH(f64FVar, f64FVMul, f64FVSinh, f64FVCosh);
//FVAR_COSH(f64FVar, f64FVMul, f64FVSinh, f64FVCosh);
//FVAR_TANH(f64FVar, f64FVConst, f64FVSub, f64FVMul, f64FVTanh);
//FVAR_ATANH(f64FVar, f64FVConst, f64FVSub, f64FVDiv, f64FVMul, f64FVAtanh);
//
//
//MAT_DECL(f64FVarFVar);
//MAT_MAKE(f64FVarFVar);
//MAT_FREE(f64FVarFVar);
//MAT_PRINT(f64FVarFVar, f64FVarFVPRINT0);
//MAT_EQUAL(f64FVarFVar, f64FVarFVEqual);
//MAT_ZERO(f64FVarFVar);
//MAT_SETELEMENT(f64FVarFVar);
//MAT_GETELEMENT(f64FVarFVar);
//MAT_SETCOL(f64FVarFVar);
//MAT_GETCOL(f64FVarFVar);
//MAT_ADD(f64FVarFVar, f64FVarFVAdd);
//MAT_SUB(f64FVarFVar, f64FVarFVSub);
//MAT_MUL(f64FVarFVar, f64FVarFVAdd, f64FVarFVMul);



//typedef struct FVar FVar;
//struct FVar {
//    f64 val;
//    f64 dot;
//};
//
//FVar FVMake(f64 val, f64 dot)  /* make AD number */
//{
//    FVar fv;
//    fv.val = val;
//    fv.dot = dot;
//    return fv;
//}
//
//void FVPrint(FVar x, const char* name) /* print AD number */
//{
//    printf("FVar (%s): {\n", name);
//    printf("\t.val = %.4f\n", x.val);
//    printf("\t.dot = %.4f\n", x.dot);
//    printf("}\n");
//}
//
//#define FVPRINT(x) FVPrint(x, #x)
//#define print_FVar(x) FVPrint(x, NULL)
//
//b32 FVEqual( FVar x, FVar y, f64 eps ) /* numerical equality of two AD numbers */
//{
//    return f64Equal( x.val, y.val, eps ) && f64Equal( x.dot, y.dot, eps );
//}
//
//FVar FVAdd(FVar x, FVar y)  /* add two AD numbers */
//{
//    return (FVar) {
//        .val = x.val + y.val,
//        .dot = x.dot + y.dot
//    };
//}
//
//FVar FVAddD(FVar x, f64 a)  /* add AD with double */
//{
//    return (FVar) {
//        .val = x.val + a,
//        .dot = x.dot
//    };
//}
//
//FVar FVDAdd(f64 a, FVar x)  /* add AD with double */
//{
//    return (FVar) {
//        .val = x.val + a,
//        .dot = x.dot
//    };
//}
//
//FVar FVSub(FVar x, FVar y)  /* subtract y from x */
//{
//    return (FVar) {
//        .val = x.val - y.val,
//        .dot = x.dot - y.dot
//    };
//}
//
//FVar FVMul(FVar x, FVar y)  /* multiply two AD numbers */
//{
//    return (FVar) {
//        .val = x.val * y.val,
//        .dot = x.val * y.dot + x.dot * y.val
//    };
//}
//
//FVar FVMulD(FVar x, f64 a)  /* multiply AD number times double */
//{
//    return (FVar) {
//        .val = x.val * a,
//        .dot = x.dot * a
//    };
//}
//
//FVar FVDMul(f64 a, FVar x)  /* multiply double times AD number */
//{
//    return (FVar) {
//        .val = x.val * a,
//        .dot = x.dot * a
//    };
//}
//
//FVar FVDiv(FVar x, FVar y)  /* divide AD by AD */
//{
//    return (FVar) {
//        .val = x.val / y.val,
//        .dot = (x.dot * y.val - x.val * y.dot) / (y.val * y.val)
//    };
//}
//
//FVar FVDivD(FVar x, f64 a)  /* divide AD by f64 */
//{
//    return (FVar) {
//        .val = x.val / a,
//        .dot = x.dot / a
//    };
//}
//
//FVar FVDDiv(f64 a, FVar x)  /* divide f64 by AD */
//{
//    return (FVar) {
//        .val = a / x.val,
//        .dot = a / ( x.val * x.val )
//    };
//}
//
//FVar FVNeg(FVar x) /* negate AD */
//{
//    return (FVar) {
//        .val = -x.val,
//        .dot = -x.dot
//    };
//}
//
//
//#if TES
//void test_fv_basic_functions()
//{
//    FVar a = { .val = 1.1, .dot = 1.0 };
//    FVar b = { .val = 2.0, .dot = 0.0 };
//    FVar c = a;
//
//    TEST_ASSERT( FVEqual( a, c, 1E-10 ) );
//    TEST_ASSERT( ! FVEqual( a, b, 1E-10 ) );
//
//    TEST_ASSERT( FVEqual( FVAdd( a, b ), FVAdd( b, a ), 1E-10 ) );
//
//    TEST_ASSERT( FVEqual( FVAddD( a, 5 ), FVDAdd( 5, a ), 1E-10 ) );
//
//    TEST_ASSERT( FVEqual( FVSub( a, b ), FVNeg( FVSub( b, a ) ), 1E-10 ) );
//
//    TEST_ASSERT( FVEqual( FVMul( a, b ), FVMul( b, a ), 1E-10 ) );
//
//    TEST_ASSERT( FVEqual( FVMulD( a, 5 ), FVDMul( 5, a ), 1E-10 ) );
//
//    TEST_ASSERT( FVEqual( FVDiv( a, b ), (FVar) { .val = 0.55, .dot = 0.5 }, 1E-10 ) );
//
//    TEST_ASSERT( FVEqual( FVDivD( a, 2 ), (FVar) { .val = 0.55, .dot = 0.5 }, 1E-10 ) );
//}
//#endif


///* elementary forward AD functions */
//
//FVar FVSqrt(FVar x) /* square root of a AD number */
//{
//    f64 tmp = sqrt(x.val);
//    return (FVar) {
//        .val = tmp,
//        .dot = 0.5 * x.dot / tmp
//    };
//}
//
//FVar FVPow(FVar x, f64 a) /* x^a, power of AD number */
//{
//    f64 tmp = pow(x.val, a - 1.0);
//    return (FVar) {
//        .val = tmp * x.val,
//        .dot = a * x.dot * tmp
//    };
//}
//
//FVar FVSin(FVar x) /* sine of AD number */
//{
//    return (FVar) {
//        .val = sin(x.val),
//        .dot = x.dot * cos(x.val)
//    };
//}
//
//FVar FVCos(FVar x) /* cosine of AD number */
//{
//    return (FVar) {
//        .val = cos(x.val),
//        .dot = -x.dot * sin(x.val)
//    };
//}
//
//FVar FVTan(FVar x) /* tangent of AD number */
//{
//    f64 tmp = cos(x.val);
//    return (FVar) {
//        .val = tan(x.val),
//        .dot = x.dot / ( tmp * tmp )
//    };
//}
//
//FVar FVAtan(FVar x) /* arctangent of AD number */
//{
//    return (FVar) {
//        .val = atan(x.val),
//        .dot = x.dot / ( 1.0 + (x.val * x.val) )
//    };
//}
//
//FVar FVExp(FVar x) /* e^a, exp of AD number */
//{
//    f64 tmp = exp(x.val);
//    return (FVar) {
//        .val = tmp,
//        .dot = x.dot * tmp
//    };
//}
//
//FVar FVLog(FVar x) /* log base e of AD number */
//{
//    return (FVar) {
//        .val = log(x.val),
//        .dot = x.dot / x.val
//    };
//}
//
//FVar FVLogAbs(FVar x) /* log base e of absolute value of AD number */
//{
//    return (FVar) {
//        .val = log(fabs(x.val)),
//        .dot = x.dot / x.val
//    };
//}
//
//FVar FVSinh(FVar x) /* hyperbolic sine of AD number */
//{
//    return (FVar) {
//        .val = sinh(x.val),
//        .dot = x.dot * cosh(x.val)
//    };
//}
//
//FVar FVCosh(FVar x) /* hyperbolic cosine of AD number */
//{
//    return (FVar) {
//        .val = cosh(x.val),
//        .dot = x.dot * sinh(x.val)
//    };
//}
//
//FVar FVTanh(FVar x) /* hyperbolic tangent of AD number */
//{
//    f64 tmp = tanh(x.val);
//    return (FVar) {
//        .val = tmp,
//        .dot = x.dot * ( 1.0 - (tmp * tmp) )
//    };
//}
//
//FVar FVAtanh(FVar x) /* hyperbolic arctangent of AD number */
//{
//    return (FVar) {
//        .val = atanh(x.val),
//        .dot = 1.0 / ( 1.0 - (x.val * x.val) )
//    };
//}
//
//
//
//
//
//#if TES
//void test_fv_elementary_functions()
//{
//#define P 1E-10
//
//    FVar a = { .val = 2.0, .dot = 1.0 };
//    FVar b = { .val = 0.5, .dot = 1.0 };
//
//
//    TEST_ASSERT(
//        FVEqual(
//            FVSqrt(a),
//            (FVar) { .val = sqrt(a.val), .dot = 1 / (2 * sqrt(a.val)) },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVPow( a, 2 ),
//            (FVar) { .val = pow(a.val, 2.0), .dot = 2 * a.val },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVSin( a ),
//            (FVar) { .val = sin(a.val), .dot = cos(a.val) },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVCos( a ),
//            (FVar) { .val = cos(a.val), .dot = -sin(a.val) },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVTan( a ),
//            FVDiv( FVSin(a), FVCos(a) ),
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVAtan( a ),
//            (FVar) { .val = atan(a.val), .dot = 1 / ( 1 + (a.val * a.val) ) },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVExp( a ),
//            (FVar) { .val = exp(a.val), .dot = exp(a.val) },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVLog( a ),
//            (FVar) { .val = log(a.val), .dot = 1 / a.val },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVLogAbs( FVNeg(a) ),
//            (FVar) { .val = log( a.val ), .dot = 1 / a.val },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVSinh( a ),
//            (FVar) { .val = sinh( a.val ), .dot = cosh( a.val ) },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVCosh( a ),
//            (FVar) { .val = cosh( a.val ), .dot = sinh( a.val ) },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVTanh( a ),
//            (FVar) { .val = tanh( a.val ), .dot = 1 - tanh( a.val ) * tanh( a.val ) },
//            P
//        )
//    );
//
//    TEST_ASSERT(
//        FVEqual(
//            FVAtanh( b ),
//            (FVar) { .val = atanh( b.val ), .dot = 1.0 / ( 1.0 - (b.val * b.val) ) },
//            P
//        )
//    );
//
//
//
//#undef P
//}
//#endif
//
//
//MAT_DECL(FVar);
//MAT_MAKE(FVar);
//MAT_FREE(FVar);
//MAT_PRINT(FVar, FVPRINT);
//MAT_EQUAL(FVar, FVEqual);
//MAT_ZERO(FVar);
//MAT_SETELEMENT(FVar);
//MAT_GETELEMENT(FVar);
//MAT_SETCOL(FVar);
//MAT_GETCOL(FVar);
//MAT_ADD(FVar, FVAdd);
//MAT_SUB(FVar, FVSub);
//MAT_MUL(FVar, FVAdd, FVMul);





