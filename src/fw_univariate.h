
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

b32 FVEqual( FVar x, FVar y, f64 eps ) /* numerical equality of two AD numbers */
{
    return F64Equal( x.val, y.val, eps ) && F64Equal( x.dot, y.dot, eps );
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

FVar FVDivD(FVar x, double a)  /* divide AD by double */
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


#if TEST
void test_fv_basic_functions()
{
    FVar a = { .val = 1.1, .dot = 1.0 };
    FVar b = { .val = 2.0, .dot = 0.0 };
    FVar c = a;
    
    TEST_ASSERT( FVEqual( a, c, 1E-10 ) );
    TEST_ASSERT( ! FVEqual( a, b, 1E-10 ) );
    
    TEST_ASSERT( FVEqual( FVAdd( a, b ), FVAdd( b, a ), 1E-10 ) );
    
    TEST_ASSERT( FVEqual( FVAddD( a, 5 ), FVDAdd( 5, a ), 1E-10 ) );
    
    TEST_ASSERT( FVEqual( FVSub( a, b ), FVNeg( FVSub( b, a ) ), 1E-10 ) );
    
    TEST_ASSERT( FVEqual( FVMul( a, b ), FVMul( b, a ), 1E-10 ) );
    
    TEST_ASSERT( FVEqual( FVMulD( a, 5 ), FVDMul( 5, a ), 1E-10 ) );
    
    TEST_ASSERT( FVEqual( FVDiv( a, b ), (FVar) { .val = 0.55, .dot = 0.5 }, 1E-10 ) );
    
    TEST_ASSERT( FVEqual( FVDivD( a, 2 ), (FVar) { .val = 0.55, .dot = 0.5 }, 1E-10 ) );
}
#endif


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


#if TEST
void test_fv_elementary_functions()
{
    FVar a = { .val = 2.0, .dot = 1.0 };
    
    TEST_ASSERT( FVEqual( FVSqrt(a), (FVar) { .val = 1.414, .dot = 0.353 }, 1E-3 ) );
    
    TEST_ASSERT( FVEqual( FVPow( a, 2 ), (FVar) { .val = 4.0, .dot = 4.0 }, 1E-10 ) );
    
//    FVSin
    
//    FVCos
    
//    FVTan
    
//    FVAtan
    
//    FVExp
    
//    FVExp
    
//    FVLog
    
//    FVLogAbs
    
//    FVSinh
    
//    FVCosh
    
//    FVTanh
    
//    FVAtanh
}
#endif


MAT_DECL(FVar);
MAT_MAKE(FVar);
MAT_FREE(FVar);
MAT_PRINT(FVar);
MAT_ZERO(FVar);
MAT_SETELEMENT(FVar);
MAT_GETELEMENT(FVar);
MAT_SETCOL(FVar);
MAT_GETCOL(FVar);
MAT_ADD(FVar, FVAdd);
MAT_SUB(FVar, FVSub);
MAT_MUL(FVar, FVAdd, FVMul);


