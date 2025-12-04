# This function finds the datanugget mean of a given datanugget object
# It returns the dataframe of the datanugget means.

# Function Inputs:
# x is the original large dataset
# DN is the given datanugget


getDN.means <- function(x, DN){
  
  
  
  ## -------------------------------
  ## Argument checks
  ## -------------------------------
  
  # check x
  if (!any(class(x) %in% c("matrix", "data.frame", "data.table"))){
    stop(" x must be of class \"matrix\", \"data.frame\", or \"data.table\" ")
  }
  
  
  # check DN
  if (!inherits(DN,"datanugget")){
    stop('DN must be of class "datanugget"')
  }
  
  
  ## ----------------------------------
  ## Pre processing and Initialization
  ## ----------------------------------
  
  
  # convert the original dataset to matrix
  x_mat <- as.matrix(x)
  
  
  # number of data nuggets
  n.dns <- nrow(DN$"Data Nuggets")
  
  
  # number of features(excluding datanugget number, scale and weight)
  n.features <- ncol(DN$"Data Nuggets") - 3
  
  
  # datanugget centers
  DN.centers <- DN$"Data Nuggets"[, 2:(n.features + 1)]
  
  
  # datanugget assignments
  DN.assignments <- DN$"Data Nugget Assignments"
  
  
  ## ----------------------
  ## Working functions
  ## ----------------------
  
  
  out <- t(sapply(1:n.dns, function(i){
    
    
    index_vector <- which(DN.assignments == i)
    
    
    obs.dn <- x_mat[index_vector, , drop = FALSE]
    
    
    center <- colMeans(obs.dn)
  }))
  
  colnames(out) <- colnames(DN.centers)
  
  out <- as.data.frame(out)
  
  return(out)
  
}



