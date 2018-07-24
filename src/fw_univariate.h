
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



void FVarGradient( FVar f( FVarMat ), FVarMat input, FVarMat grad )
{
    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
    
    FVar tmp;
    u32 N = input.dim0;
    
    FVarMat xCpy = FVarMatMake( DefaultAllocator, N, 1 );
    
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





