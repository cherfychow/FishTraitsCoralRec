######################################################################################
# Function to compute trait richness index TOP (Fontana et al. 2015)                 
# Partially adapted from Villeger et al. 2008
#
#                                                                                     
# It requires R libraries 'geozoo', 'flexmix' and 'geometry'                                                              #                                                                              
# input:                                                                             
# - data frame 'traits': every row represents an individual, every column a trait
# - number of individuals must be higher than number of traits    
# - NA are not allowed                                                                         
# - ADVICE: it is generally meaningful to standardise trait values (mean=0 and sd=1)                                      
######################################################################################


##################
##  TOP index   ##
##################

library(geometry)


TOP.index <- function(traitdat){
  
  # TOP
  
  dim1 <- ncol(traitdat)
  
  #definitions: index i, area as empty vector
  
  i=0
  area<-matrix(ncol=2,nrow=nrow(traitdat))
  
  while(nrow(traitdat)>dim1){
    i=i+1
    
    # use of convhulln function
    
    # area
    area[i,2] <- convhulln(traitdat,"FA")$area
    
    # identity of vertices
    vert0<-convhulln(traitdat,"Fx TO 'vert.txt'")
    vert1<-scan("vert.txt",quiet=T)
    vert2<-vert1+1
    
    vertices <- vert2[-1]
    
    traitdat <- traitdat[-vertices,]
    
    area[i,1] <- length(vertices)
    
  }
  
  
  area<-na.omit(area)
  
  # Output (2 numbers): Number of points touched by areas; Sum of the areas (TOP index)
  colSums(area)
  
}
