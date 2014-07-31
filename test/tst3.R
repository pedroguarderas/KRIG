#___________________________________________________________________________________________________
library(RKHSENS)
options( stringsAsFactors = FALSE )

kernel<-function( x, y ) min(x,y)
Xi<-data.frame( x1 = c( 1.0, 0.5 ) , x2 = c( 1.0, 2.0 ) )
gamma_integrate( Xi, kernel, 0, 1 )
