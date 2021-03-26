
###########################################
# Convex hull function (for plotting)
# This function calculates a two dimensional convex hull for plotting
# returns a data frame of the hull vertices in clockwise order

# Dependencies: requires the geometry package convhull function
###########################################

convhull.vert <- function(x) {
  
  # check for dependencies
  require(geometry)
  
  # check x
  if (is.data.frame(x)==FALSE) {stop('x must be a data frame or tibble object')}
  if (ncol(x)!=2) {
    stop('x contains an incorrect number of dimensions- it must contain 2. Use convhulln() if you require hull calculation for > 2 dimensions')
  }
  
  # where x is the dataframe of multivariate analysis generated point coordinates
  # one column each for axis a and b
  
  hull <- convhulln(x, 'FA')
  
  # use the convex hull results to create x,y vertices coordinates
  # because the hull list is by row number, we use the hull object to index from points dataframe
  data.frame(
    do.call(
      rbind,
      lapply(1:nrow(hull$hull), function(i) {
        rbind(x[hull$hull[i,1],], x[hull$hull[i,2],])
      })
    )
  ) -> hull.v
  
  # reorder clockwise
  hull.v <- hull.v[order(-1 * atan2(hull.v[,2] - mean(range(hull.v[,2])), hull.v[,1] - mean(range(hull.v[,1])))),]
  hull.v <- rbind(hull.v, hull.v[1,])
  return(hull.v)
}