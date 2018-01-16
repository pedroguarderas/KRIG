# Test 2 -------------------------------------------------------------------------------------------
context( "Sensitivity analysis" )

test_that( "Sum var verification", {
  options( stringsAsFactors = FALSE )
  
  kernel_1<-function( x, y ) exp( -0.5*(x-y)^2)
  kernel_2<-function( x, y ) exp( -0.7*(x-y)^2)
  kernel_3<-function( x, y ) exp( -0.1*(x-y)^2)
  
  Kernels<-data.frame( kernel = c( 'kernel_1', 'kernel_2', 'kernel_3' ), 
                       min = c( 0, -1, -5 ), 
                       max = c( 1, 1, 5 ),
                       n = c( 500, 500, 500 ) )
  
  n<-20
  X<-matrix( c( seq( 0, 1, length.out = n ), 
                seq( -1, 1, length.out = n ),
                seq( -5, 5, length.out = n ) ), n, 3 )
  
  KI<-vector_integrate_kernel( Kernels, X )
  
  # Returns Gamma array and Kanova matrix
  GK<-Kanova( Kernels, KI, X )
  
  f<-function( x ) abs( x[1] + 30 * x[2] + 60 * x[3] )
  Func<-apply( X, 1, FUN = f )
  
  KF<-solve( GK$Kanova + diag( 1e-8, n, n ), Func )
  
  SbI<-NULL
  for ( j in 1:3 ) {
    CB<-combn( 1:3, j )  
    for ( l in 1:ncol( CB ) ) {
      SbI<-c( SbI, sens_idx( KF, CB[,l], X, GK$Gamma ) )
      names(SbI)[length(SbI)]<-paste( 'C.', paste( CB[,l], collapse='.' ), sep = '' )
    }
  }
  
  Var<-sens_var( KF, GK$Gamma )
  
  SVar<-sum( SbI / Var )

  
  expect_less_than( abs( SVar - 1 ), 1e-6 )
  
} )

