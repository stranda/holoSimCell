#Not necessary to include this in the package, just for the MSU cluster
#These lines need to be run for each user before trying to run our simulations

#Setting up the libraries and packages...
module load GCC/6.4.0-2.28  
module load OpenMPI/2.1.2
module load GSL/2.4
module load R/3.5.0-X11-20180131
module load fastsimcoal 

#Start R
R

#Now in R, install packages
install.packages("cubature")  	
install.packages("Rcpp")  	
install.packages("RcppGSL")	
install.packages("matrixStats")	
install.packages("tripack")		
install.packages("vegan")		
install.packages("geoR")		
install.packages("LDcorSV")		
install.packages("sna")			
install.packages("igraph")		
install.packages("strataG")		

#Then test...
library("cubature")  	
library("Rcpp")  		
library("RcppGSL")
library("matrixStats")	
library("tripack")		
library("vegan")		
library("geoR")			
library("LDcorSV")		
library("sna")			
library("igraph")		
library("strataG")		
