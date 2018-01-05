# Reading data -------------------------------------------------------------------------------------

# Loading libraries --------------------------------------------------------------------------------
library( data.table )

path<-'/home/aju/Documents/Lectures/spatial_statistics/data/'

# Simulated example
data_file<-paste0( path, 'samples.dat' )

SG<-read.csv( file = data_file, header = FALSE, skip = 3, sep = ' ', dec = '.', 
              colClasses = c( 'numeric', 'numeric', 'numeric' ) )
SG<-as.data.table( SG )
SG[ , V1 := NULL ]
setnames( SG, c( 'x1', 'x2', 'Z' ) )                       

# Brenda mines
data_file<-paste0( path, 'BrendaMines.dat' )

BM<-read.csv( file = data_file, header = FALSE, skip = 4, sep = ',', dec = '.' )
BM<-as.data.table( BM )
BM[ , V1 := NULL ]
setnames( BM, c( 'a', 's', 'x1', 'x2', 'x3', 'Z1', 'Z2', 'C1', 'C2', 'C3' ) )  

save( SG, BM, file = 'RData/mining_data.RData' )

# Copper mine
data_file<-paste0( path, 'Copper.dat' )

Copper<-read.csv( file = data_file, header = FALSE, skip = 2, sep = '\t', dec = '.',
                  colClasses = c( 'numeric', 'character', 'numeric', 'numeric', 'numeric', 
                                  'numeric', 'numeric', 'numeric' ) )
Copper<-as.data.table( Copper )
Copper[ , V1 := NULL ]
setnames( Copper, c( 'a', 's', 'x1', 'x2', 'x3', 'Z', 'C' ) )  

save( SG, BM, Copper, file = 'RData/mining_data.RData' )
rm( list =  ls() )
gc()
