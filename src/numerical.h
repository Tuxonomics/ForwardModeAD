#ifndef NUMERICAL_H
#define NUMERICAL_H

#include "utilities.h"

/* Finite Difference */
FVar FDiff( FVar target(FVar *), FVar *params, u32 dim, f64 h ) {
    
    FVar *params2 = (FVar *) Calloc(DefaultAllocator, dim, sizeof(FVar));
    for (u32 i=0; i<dim; ++i) {
        params2[i] = FVMake(dim, params[i].val, i);
    }
    
    FVar t2;
    FVar t = target( params );
    
    FVar res = FVInit(dim);
    res.val = t.val;
    
    for (u32 i=0; i<dim; ++i) {
        params2[i].val += h;
        
        t2 = FVMulD( FVSub( target(params2), target(params) ), 1/h );
        
        res.dot[i] = t2.val;
        
        params2[i].val -= h;
    }
    
    return res;
}


/* Central Difference */
FVar CDiff( FVar target(FVar *), FVar *params, u32 dim, f64 h ) {

    FVar *params2 = (FVar *) Calloc(DefaultAllocator, dim, sizeof(FVar));
    FVar *params3 = (FVar *) Calloc(DefaultAllocator, dim, sizeof(FVar));
    for (u32 i=0; i<dim; ++i) {
        params2[i] = FVMake(dim, params[i].val, i);
        params3[i] = FVMake(dim, params[i].val, i);
    }
    
    FVar t2;
    FVar t = target( params );
    
    FVar res = FVInit(dim);
    res.val = t.val;
    
    for (u32 i=0; i<dim; ++i) {
        params2[i].val += h;
        params3[i].val -= h;
        
        t2 = FVMulD( FVSub( target(params2), target(params3) ), 1/(2*h) );
        
        res.dot[i] = t2.val;
        
        params2[i].val -= h;
        params3[i].val += h;
    }
    
    return res;
}


#endif /* NUMERICAL_H */