#import "fw_univariate.h"


///* Finite Difference */
//void FVFDiff( Allocator al, FVar f( FVarMat ), FVarMat input, FVarMat grad, f64 h ) {
//    
//    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
//    
//    FVar tmp;
//    u32 N = input.dim0;
//    
//    FVarMat xCpy = FVarMatMake( al, N, 1 );
//    
//    for ( u32 i=0; i<N; ++i ) {
//        xCpy.data[i]     = input.data[i];
//        xCpy.data[i].dot = 0;
//    }
//    
//    
//    for ( u32 i=0; i<N; ++i ) {
//        xCpy.data[i].val += h;
//        
//        tmp = FVMulD( FVSub( f(xCpy), f(input) ), 1/h );
//        tmp.dot = 0;
//        
//        grad.data[i] = tmp;
//        
//        xCpy.data[i].val -= h;
//    }
//    
//    FVarMatFree( al, &xCpy );
//}
//
//
///* Central Difference */
//void FVCDiff( Allocator al, FVar f( FVarMat ), FVarMat input, FVarMat grad, f64 h ) {
//
//    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
//    
//    FVar tmp;
//    u32 N = input.dim0;
//    
//    FVarMat xCpy  = FVarMatMake( al, N, 1 );
//    FVarMat xCpy2 = FVarMatMake( al, N, 1 );
//    
//    for ( u32 i=0; i<N; ++i ) {
//        xCpy.data[i]      = input.data[i];
//        xCpy.data[i].dot  = 0;
//        xCpy2.data[i]     = input.data[i];
//        xCpy2.data[i].dot = 0;
//    }
//    
//    for ( u32 i=0; i<N; ++i ) {
//        xCpy.data[i].val  += h;
//        xCpy2.data[i].val -= h;
//        
//        tmp = FVMulD( FVSub( f(xCpy), f(xCpy2) ), 1/(2*h) );
//        tmp.dot = 0;
//        
//        grad.data[i]     = tmp;
//        
//        xCpy.data[i].val  -= h;
//        xCpy2.data[i].val += h;
//    }
//    
//    FVarMatFree( al, &xCpy  );
//    FVarMatFree( al, &xCpy2 );
//}
//
//
///* AD gradient */
//void FVGradient( Allocator al, FVar f( FVarMat ), FVarMat input, FVarMat grad )
//{
//    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
//    
//    FVar tmp;
//    u32 N = input.dim0;
//    
//    FVarMat xCpy = FVarMatMake( al, N, 1 );
//    
//    for ( u32 i=0; i<N; ++i ) {
//        xCpy.data[i]     = input.data[i];
//        xCpy.data[i].dot = 0;
//    }
//    
//    for ( u32 i=0; i<N; ++i ) {
//        xCpy.data[i].dot = 1.0;
//        
//        tmp     = f( xCpy );
//        tmp.val = tmp.dot;
//        tmp.dot = 0;
//        
//        grad.data[i]     = tmp;
//        xCpy.data[i].dot = 0;
//    }
//    
//    FVarMatFree( al, &xCpy );
//}
//
//
//FVar test_test_f( FVarMat input )
//{
//    return FVTanh( input.data[0] );
//}
//
//#if TEST
//void test_grad()
//{
//#define P 1E-8
//    
//    FVarMat input  = FVarMatMake( DefaultAllocator, 1, 1 );
//    FVarMat gradC  = FVarMatMake( DefaultAllocator, 1, 1 );
//    FVarMat gradF  = FVarMatMake( DefaultAllocator, 1, 1 );
//    FVarMat gradAD = FVarMatMake( DefaultAllocator, 1, 1 );
//    
//    input.data[0] = (FVar) { .val = 0.5, .dot = 1.0 };
//    
//    FVFDiff( DefaultAllocator, test_test_f, input, gradF, P );
//    
//    FVCDiff( DefaultAllocator, test_test_f, input, gradC, P );
//    
//    FVGradient( DefaultAllocator, test_test_f, input, gradAD );
//    
//    TEST_ASSERT( FVarMatEqual( gradF, gradC,  P ) );
//    TEST_ASSERT( FVarMatEqual( gradC, gradAD, P ) );
//    
//#undef P
//}
//#endif
//
//
////void hessian( Allocator al, FVar f( FVarMat ), FVarMat input, FVarMat grad, FVarMat hess )
////{
////    ASSERT( input.dim0 == grad.dim0 && input.dim1 == grad.dim1 && input.dim1 == 1 );
////    ASSERT( input.dim0 == hess.dim0 && input.dim0 == hess.dim1 );
////
////
////
////}
//
////template <typename T, typename F>
////void hessian(const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& x, T& fx,
////             Eigen::Matrix<T, Eigen::Dynamic, 1>& grad,
////             Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H) {
////    H.resize(x.size(), x.size());
////    grad.resize(x.size());
////    // size 0 separate because nothing to loop over in main body
////    if (x.size() == 0) {
////        fx = f(x);
////        return;
////    }
////    Eigen::Matrix<fvar<fvar<T> >, Eigen::Dynamic, 1> x_fvar(x.size());
////    for (int i = 0; i < x.size(); ++i) {
////        for (int j = i; j < x.size(); ++j) {
////            for (int k = 0; k < x.size(); ++k)
////                x_fvar(k) = fvar<fvar<T> >(fvar<T>(x(k), j == k), fvar<T>(i == k, 0));
////            fvar<fvar<T> > fx_fvar = f(x_fvar);
////            if (j == 0)
////                fx = fx_fvar.val_.val_;
////            if (i == j)
////                grad(i) = fx_fvar.d_.val_;
////            H(i, j) = fx_fvar.d_.d_;
////            H(j, i) = H(i, j);
////        }
////    }
////}
////
////}  // namespace math
////}  // namespace stan
//
