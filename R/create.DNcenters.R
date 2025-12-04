# This function creates the centers of data nuggets(DN) from a random sample.
# It returns the dataframe with DN.num of DN centers.

# Function inputs
# RS: random sample of observations.
# delete.percent: proportion of data points to be deleted at each iteration.
# DN.num: final number of DNs to retain.
# dist.metric: pairwise distance measure (e.g., "euclidean" or "manhattan").
# make.pbs: logical; whether to show a progress bar while the function runs.


create.DNcenters = function(RS, 
                            delete.percent = .1, 
                            DN.num,
                            dist.metric = "euclidean", 
                            make.pbs = FALSE){ 
  
  
  ## ------------------------------- 
  ## Argument checks 
  ## ------------------------------- 
  
  # check RS 
  if (!any(class(RS) %in% c("matrix", "data.frame", "data.table"))) { 
    stop("RS must be of class \"matrix\", \"data.frame\", or \"data.table\"") 
  } 
  
  
  # check delete.percent 
  if (!is.numeric(delete.percent)) { 
    stop("delete.percent must be of class \"numeric\"") 
  } 
  
  
  # make sure delete.percent is between 0 and 1 
  if (delete.percent <= 0 | delete.percent >= 1){ 
    stop("delete.percent must be within (0,1)") 
  } 
  
  
  # check DN.num 
  if (!(class(DN.num) %in% c("numeric", "integer"))) { 
    stop("DN.num must be of class \"numeric\" or \"integer\"") 
  } 
  
  
  # check dist.metric 
  if (dist.metric != "euclidean" & dist.metric != "manhattan") { 
    stop("dist.metric must be \"euclidean\" or \"manhattan\"") 
  }
  
  
  # check make.pbs 
  if (!is.logical(make.pbs)) { 
    stop("make.pbs must be TRUE OR FALSE") 
  } 
  
  
  
  ## ---------------------------------- 
  ## Pre processing and Initialization 
  ## ---------------------------------- 
  
  
  # Convert RS to dataframe 
  RS.data <- as.data.frame(RS) 
  
  
  # reset the rownames to numbers 
  rownames(RS.data) <- 1:nrow(RS.data) 
  
  
  # storing the generated DN centers
  DN.data <- RS.data 
  
  
  # storing the remaining no of rows after removal until completion 
  tmp.num <- nrow(DN.data) 
  
  
  
  ## ------------------------- 
  ## Pairwise Distance Matrix 
  ## ------------------------- 
  
  
  # check that nrow(RS) is not smaller than DN.num
  if (tmp.num <= DN.num){
    stop("DN.num cannot be more than the number of observations in the random sample")
  }
  
  
  # Compute the pairwise distance matrix
  DN.dist.matrix <- as.matrix(Rfast::Dist(DN.data, method = dist.metric)) 
  
  
  
  # eliminate the diagonals from becoming a closest neighbor choice 
  diag(DN.dist.matrix) <- max(DN.dist.matrix) + 1 
  
  
  
  # check if user wants a progress bar 
  if (make.pbs){ 
    
    pbs.max <- tmp.num - DN.num 
    
    # initialize progress bar 
    pb <- txtProgressBar(min = 0, max = pbs.max) 
    
    # initialize the value for the progress bar 
    for.prog <- 0 
    
  } 
  
  
  
  ## ------------------------- 
  ## Working Function 
  ## ------------------------- 
  
  
  # eliminate the desired number data points that are closest together 
  while (tmp.num > DN.num){ 
    
    
    # retrieve the number of sample points to delete 
    delete.num <- floor(tmp.num * delete.percent) 
    
    
    # check if the current delete.num would delete too many entries 
    if((tmp.num - DN.num) <= delete.num) {
      
      delete.num = tmp.num - DN.num 
    
    }
    
    # break if delete.num equals 0
    if (delete.num <= 0) break 
    
    
    # find the smallest pairwise distance for each sample point 
    row_min <- Rfast::rowMins(DN.dist.matrix) 
    
    
    # selected row indices for removal 
    delete.entries <- order(row_min)[1:delete.num] 
    
    
    # delete the selected rows 
    DN.data <- DN.data[-delete.entries, ,drop = FALSE] 
    
    
    # delete the selected rows from the distance matrix 
    DN.dist.matrix <- DN.dist.matrix[-delete.entries, -delete.entries, 
                                     drop = FALSE] 
    
    
    # remaining no of datapoints 
    tmp.num = nrow(DN.data) 
    
    
    # check if user wants a progress bar 
    if (make.pbs == TRUE){ 
      
      for.prog <- for.prog + length(delete.entries) 
      
      # update the progress bar 
      setTxtProgressBar(pb, for.prog) 
    } 
  } 
  
  if (make.pbs) close(pb) 
  
  return(DN.data) 
}




