
typedef struct FVar FVar;
struct FVar {
    f64 val;
    u32 dim;
    f64 *dot;
};

FVar FVInit(u32 dim)  /* init AD number */
{
    FVar fv;
    fv.dim = dim;
    fv.dot = Calloc(DefaultAllocator, dim, sizeof(f64));
    return fv;
}

FVar FVMake(u32 dim, f64 val, u32 idx)  /* make AD number */
{
    FVar fv = FVInit(dim);
    fv.val = val;
    
    fv.dot[idx] = 1.0;

    return fv;
}

void FVFree(FVar x) {
    ASSERT(x.dot);
    Free(DefaultAllocator, x.dot);
}

void FVPrint(FVar x, const char* name) /* print AD number */
{
    printf("FVar (%s): {\n", name);
    printf("\t.val = %.4f\n", x.val);
    printf("\t.dim = %u\n", x.dim);
    printf("\t.dot = {\n");
    for (u32 i=0; i<x.dim; ++i) {
        printf("\t\t[%d] = %.4f\n", i, x.dot[i]);
    }
    printf("\t}\n");
    printf("}\n");
}

#define FVPRINT(x) FVPrint(x, #x)

FVar FVAdd(FVar x, FVar y)  /* add two AD numbers */
{
    u32 dim = x.dim;
    
    ASSERT(dim == y.dim);
    
    FVar fv = FVMake(dim, x.val + y.val, 0);
    
    for (u32 i=0; i<dim; ++i) {
        fv.dot[i] = x.dot[i] + y.dot[i];
    }
    
    return fv;
}

FVar FVAddD(FVar x, f64 a)  /* add AD with double */
{
    FVar fv = FVInit(x.dim);
    fv.val = x.val + a;
    
    memcpy(fv.dot, x.dot, x.dim * sizeof(f64));
    
    return fv;
}

FVar FVDAdd(f64 a, FVar x)  /* add AD with double */
{
    FVar fv = FVInit(x.dim);
    fv.val = x.val + a;
    
    memcpy(fv.dot, x.dot, x.dim * sizeof(f64));
    
    return fv;
}

FVar FVSub(FVar x, FVar y)  /* subtract y from x */
{
    u32 dim = x.dim;
    
    ASSERT(dim == y.dim);
    
    FVar fv = FVInit(dim);
    fv.val = x.val - y.val;
    
    for (u32 i=0; i<dim; ++i) {
        fv.dot[i] = x.dot[i] - y.dot[i];
    }
    
    return fv;
}

FVar FVMul(FVar x, FVar y)  /* multiply two AD numbers */
{
    u32 dim = x.dim;
    
    ASSERT(dim == y.dim);
    
    FVar fv = FVInit(dim);
    fv.val = x.val * y.val;
    
    for (u32 i=0; i<dim; ++i) {
        fv.dot[i] = x.val * y.dot[i] + x.dot[i] * y.val;
    }
    
    return fv;
}

FVar FVMulD(FVar x, f64 a)  /* multiply complex times double */
{
    FVar fv = FVInit(x.dim);
    fv.val = x.val * a;
    memcpy(fv.dot, x.dot, x.dim * sizeof(f64));
    
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = fv.dot[i] * a;
    }
    
    return fv;
}

FVar FVDiv(FVar x, FVar y)  /* divide AD by AD */
{
    u32 dim = x.dim;
    
    ASSERT(dim == y.dim);
    
    FVar fv = FVInit(dim);
    fv.val = x.val / y.val;

    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = (x.dot[i] * y.val - x.val * y.dot[i]) / (y.val * y.val);
    }

    return fv;
}

FVar FDDivD(FVar x, f64 a)  /* multiply complex times double */
{
    f64 tmp = 1 / a;
    
    FVar fv = FVInit(x.dim);
    fv.val = x.val * tmp;
    
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVNeg(FVar x) /* negate AD */
{
    FVar fv = FVInit(x.dim);
    fv.val = x.val;
    
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = -x.dot[i];
    }
    
    return fv;
}


/* elementary forward AD functions */

FVar FVSqrt(FVar x) /* square root of a AD number */
{
    f64 tmp = sqrt(x.val);
    
    FVar fv = FVInit(x.dim);
    fv.val = tmp;
    
    tmp = 1 / tmp;
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = 0.5 * x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVPow(FVar x, f64 a) /* x^a, power of AD number */
{
    f64 tmp = pow(x.val, a - 1.0);
    
    FVar fv = FVInit(x.dim);
    fv.val = tmp * x.val;
    
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = a * x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVSin(FVar x) /* sine of AD number */
{
    FVar fv = FVInit(x.dim);
    fv.val = sin(x.val);
    
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * cos(x.val);
    }
    
    return fv;
}

FVar FVCos(FVar x) /* cosine of AD number */
{
    FVar fv = FVInit(x.dim);
    fv.val = cos(x.val);
    
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = -x.dot[i] * sin(x.val);
    }
    
    return fv;
}

FVar FVTan(FVar x) /* tangent of AD number */
{
    f64 tmp = cos(x.val);
    tmp = 1 / ( tmp * tmp );
    
    FVar fv = FVInit(x.dim);
    fv.val = tan(x.val);
    
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVAtan(FVar x) /* arctangent of AD number */
{
    FVar fv = FVInit(x.dim);
    fv.val = atan(x.val);
    
    f64 tmp = 1.0 / ( 1.0 + (x.val * x.val) );
    
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVExp(FVar x) /* e^a, exp of AD number */
{
    f64 tmp = exp(x.val);
    
    FVar fv = FVInit(x.dim);
    fv.val = tmp;
    
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVLog(FVar x) /* log base e of AD number */
{
    FVar fv = FVInit(x.dim);
    fv.val = log(x.val);
    
    f64 tmp = 1 / x.val;
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVLogAbs(FVar x) /* log base e of absolute value of AD number */
{
    FVar fv = FVInit(x.dim);
    fv.val = log(fabs(x.val));
    
    f64 tmp = 1 / x.val;
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVSinh(FVar x) /* hyperbolic sine of AD number */
{
    FVar fv = FVInit(x.dim);
    fv.val = sinh(x.val);
    
    f64 tmp = cosh(x.val);
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVCosh(FVar x) /* hyperbolic cosine of AD number */
{
    FVar fv = FVInit(x.dim);
    fv.val = cosh(x.val);
    
    f64 tmp = sinh(x.val);
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVTanh(FVar x) /* hyperbolic tangent of AD number */
{
    f64 tmp = tanh(x.val);
    
    FVar fv = FVInit(x.dim);
    fv.val = tmp;
    
    tmp = 1.0 - (tmp * tmp);
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * tmp;
    }
    
    return fv;
}

FVar FVAtanh(FVar x) /* hyperbolic arctangent of AD number */
{
    FVar fv = FVInit(x.dim);
    fv.val = atanh(x.val);
    
    f64 tmp = 1.0 / ( 1.0 - (x.val * x.val) );
    for (u32 i=0; i<x.dim; ++i) {
        fv.dot[i] = x.dot[i] * tmp;
    }
    
    return fv;
}



/* Finite Difference */
FVar FDiff( Allocator al, FVar target(FVar *), FVar *params, u32 dim, f64 h ) {
    
    FVar *params2 = Calloc( al, dim, sizeof(FVar) );
    for ( u32 i=0; i<dim; ++i ) {
        params2[i] = FVMake( dim, params[i].val, i );
    }
    
    FVar t2;
    FVar t = target( params );
    
    FVar res = FVInit( dim );
    res.val = t.val;
    
    for ( u32 i=0; i<dim; ++i ) {
        params2[i].val += h;
        
        t2 = FVMulD( FVSub( target(params2), target(params) ), 1/h );
        
        res.dot[i] = t2.val;
        
        params2[i].val -= h;
    }
    
    return res;
}


/* Central Difference */
FVar CDiff( Allocator al, FVar target(FVar *), FVar *params, u32 dim, f64 h ) {
    
    FVar *params2 = Calloc( al, dim, sizeof(FVar) );
    FVar *params3 = Calloc( al, dim, sizeof(FVar) );
    for ( u32 i=0; i<dim; ++i ) {
        params2[i] = FVMake( dim, params[i].val, i );
        params3[i] = FVMake( dim, params[i].val, i );
    }
    
    FVar t2;
    FVar t = target( params );
    
    FVar res = FVInit( dim );
    res.val = t.val;
    
    for ( u32 i=0; i<dim; ++i ) {
        params2[i].val += h;
        params3[i].val -= h;
        
        t2 = FVMulD( FVSub( target(params2), target(params3) ), 1/(2*h) );
        
        res.dot[i] = t2.val;
        
        params2[i].val -= h;
        params3[i].val += h;
    }
    
    return res;
}


