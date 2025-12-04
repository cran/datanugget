# This function creates the datanuggets.
# It returns a list of two objects:
# 1. Data Nuggets is a dataframe with cols DN no., DN centers, Scale and Weight.
# 2. Data Nugget Assignments is the vector containing the data nugget assignment of each observation in x.


# Function inputs
# x: original dataset.
# center.method: the method used for choosing data nugget centers.
# R: the number of observations to sample from the data matrix when creating the initial data nugget centers.
# delete.percent: proportion of data points to be deleted at each iteration.
# DN.num1: the number of initial data nugget centers to create.
# DN.num2: the number of final data nuggets to create.
# dist.metric: pairwise distance measure (e.g., "euclidean" or "manhattan").
# seed: random seed for replication.
# no.cores: number of cores used for parallel processing.
# make.pbs: logical; whether to show a progress bar while the function runs.


create.DN = function(x,
                     center.method = "original",
                     R = 5000,
                     delete.percent = .1,
                     DN.num1 = 10^4,
                     DN.num2 = 2000,
                     dist.metric = "euclidean", 
                     seed = 291102,
                     no.cores = (detectCores() - 1),
                     make.pbs = FALSE){
  
  
  ## -------------------------------
  ## Argument checks
  ## -------------------------------
  
  
  # check x
  if (!any(class(x) %in% c("matrix", "data.frame", "data.table"))){
    stop("x must be of class \"matrix\", \"data.frame\", or \"data.table\"")
  }
  
  
  
  # check center.method
  if (!(center.method %in% c("mean", "random", "original"))){
    stop("center.method must be \"mean\" or \"random\" or \"original\"")
  }
  
  
  
  # check R
  if (!(class(R) %in% c("numeric", "integer"))){
    stop("R must be of class \"numeric\" or \"integer\"")
  }
  
  
  
  # make sure R is between 100 and 10000
  if (R < 100 | R > 10000){
    stop("R must be within [100,10000]")
  }
  
  
  
  # check delete.percent 
  if (!is.numeric(delete.percent)){
    stop("delete.percent must be of class \"numeric\"")
  }
  
  
  
  # make sure delete.percent is between 0 and 1
  if (delete.percent <= 0 | delete.percent >= 1){
    stop("delete.percent must be within (0,1)")
  }
  
  
  
  # check DN.num1
  if (!(class(DN.num1) %in% c("numeric", "integer"))){
    stop("DN.num1 must be of class \"numeric\" or \"integer\"")
  }
  
  
  
  # check DN.num2 
  if (!(class(DN.num2) %in% c("numeric", "integer"))){
    stop("DN.num2 must be of class \"numeric\" or \"integer\"")
  }
  
  
  
  # make sure DN.num2 <= DN.num1
  if (DN.num2 > DN.num1){
    stop("DN.num2 must be less than DN.num1")
  }
  
  
  
  # check dist.metric 
  if (dist.metric != "euclidean" & dist.metric != "manhattan"){
    stop("dist.metric must be \"euclidean\" or \"manhattan\"")
  }
  
  
  
  # check seed 
  if (!is.numeric(seed)){
    stop("seed must be of class \"numeric\"")
  }
  
  
  
  # check no.cores 
  if (!(class(no.cores) %in% c("numeric", "integer"))){
    stop("no.cores must be of class \"numeric\" or \"integer\"")
  }
  
  
  
  # check make.pbs 
  if (!is.logical(make.pbs)){
    stop("make.pbs must be TRUE OR FALSE")
  }
  
  
  
  ## ----------------------------------
  ## Pre processing and Initialization
  ## ----------------------------------
  
  
  # set the seed
  set.seed(seed)
  
  
  
  # number of observations in the data
  obs.num <- nrow(x)
  RS.num <- obs.num
  x_mat <- as.matrix(x)
  
  
  
  # create random sample ordering
  RS.loc <- sample(obs.num, RS.num, replace = FALSE)
  
  
  
  # batch size
  set.num <- R 
  
  
  # number of splits needed to divide RS.num into batches of size set.num
  RS.splits <- ceiling(RS.num/set.num)
  
  
  
  # check if user wants to use parallel processing
  if (no.cores > 0){
    
    
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Package 'parallel' is required for parallel processing")
    }
    
    
    
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Package 'foreach' is required for parallel processing")
    }
    
    
    
    if (!requireNamespace("doSNOW", quietly = TRUE)) {
      stop("Package 'doSNOW' is required for parallel processing with progress bars")
    }
    
    
    # create the cluster for parallel processing
    cl <- parallel::makeCluster(no.cores)
    
    # engage the cluster for parallel processing
    doSNOW::registerDoSNOW(cl)
    
  }
  
  
  
  ## ------------------------------------------
  ## Helper functions
  ## -------------------------------------------
  
  
  # Function to call create.DNcenters and get the DN centers for each RS.splits
  get_centers_for_split <- function(tmp.RS, DN_num_for_split){
    
    
    out <- create.DNcenters(RS = tmp.RS, 
                            delete.percent = delete.percent, 
                            DN.num = DN_num_for_split, 
                            dist.metric = dist.metric, 
                            make.pbs = FALSE)
    
    
    return(as.data.frame(out))
  }
  
  
  
  # Function to assign each observation to the nearest DN center
  assign_each <- function(obs, DN_centers_data){
    
    if(dist.metric == "euclidean"){
      
      diffs <- sweep(DN_centers_data, 2, obs)
      dists <- rowSums(diffs^2)
      which.min(dists)
      
      }else{
        
        diffs <- sweep(DN_centers_data, 2, obs)
        dists <- rowSums(abs(diffs))
        which.min(dists)
      }
    
  }
  
  
  # Function for recalculating center and calculating scale for each data nugget
  calculate_center_and_scale <- function(i, DN_centers_data){
    
    index_vector <- which(DN.assignments == i)
    
    obs.dn <- x_mat[index_vector, , drop = FALSE]
    
    if(nrow(obs.dn) == 1){
      
      return(as.vector(c(obs.dn, 0)))
    }
    
    
    center <- 
      if (center.method == "mean") {
        
        colMeans(obs.dn)
      
      }else if (center.method == "random"){
          
        obs.dn[sample(nrow(obs.dn), 1), ]
        
      }else{
        
        DN_centers_data[i, ]
      }
    
    
    scale <- 
     if (ncol(x_mat) == 1) {
       
       as.numeric(var(obs.dn))
       
       }else{
         
         sum(diag(cov(obs.dn)))/ncol(x_mat)
       
       }
    
    return(as.vector(c(center, scale)))
  }
  
  
  
  ## ------------------------------------------
  ## Create initial DN centers of size DN.num1
  ## -------------------------------------------
  
  
  # create initial set of data nugget centers
  message("Creating initial set of data nugget centers... ") 
  
  
  # Storage
  DN.centers1 <- NULL

  
  # check if user wants to use parallel processing
  if (no.cores > 0){
    
    # check if user wants a progress bar
    if (make.pbs){
      
      # initialize progress bar
      pb <- txtProgressBar(min = 0, max = RS.splits)
      
      # update the progress bar
      progress <- function(n){setTxtProgressBar(pb, n)}
      opts <- list(progress = progress)
      
    }else{
        
      opts <- NULL
    }
    
    DN.centers1 <- foreach::foreach(
      i = 1:RS.splits,
      .combine = rbind,
      .export = c("create.DNcenters"),
      .packages = c("Rfast", "parallel"),
      .options.snow = opts
      ) %dopar% {
        
        
        # set the beginning and endpoint indices for this batch
        RS.loc.start <- ((i-1) * set.num) + 1
        # take the min to avoid going past the end if the last split is smaller
        RS.loc.end <- min((set.num * i), length(RS.loc)) 
          
          
          
        # selected random sample for this batch
        tmp.RS <- x_mat[RS.loc[RS.loc.start:RS.loc.end], ,drop = FALSE]
          
          
          
        # size of the DN for this batch
        DN_num_for_split = min(ceiling(DN.num1/RS.splits), 
                               RS.loc.end - RS.loc.start + 1)
          
          
          
        # reduce the batch to size DN.num1/RS.splits
        get_centers_for_split(tmp.RS = tmp.RS, 
                              DN_num_for_split = DN_num_for_split)
        }
      
      # close the progress bar
      if (make.pbs) close(pb)
    
    }else{
      
      # check if user wants a progress bar
      if (make.pbs){
        
        # initialize progress bar
        pb <- txtProgressBar(min = 0, max = RS.splits)
      }
      
      
      
      # loop through each batch
      for (i in 1:RS.splits){
        
        # set the beginning and endpoint indices for this batch
        RS.loc.start <- ((i-1) * set.num) + 1
        # take the min to avoid going past the end if last split is smaller
        RS.loc.end <- min((set.num * i), length(RS.loc)) 
        
        
        
        # selected random sample for this batch
        tmp.RS <- x_mat[RS.loc[RS.loc.start:RS.loc.end], ,drop = FALSE]
        
        
        # size of the DN for this batch
        DN_num_for_split = min(ceiling(DN.num1/RS.splits), 
                               RS.loc.end - RS.loc.start + 1)
        
        
        # reduce the batch to size DN.num1/RS.splits
        out <- get_centers_for_split(tmp.RS = tmp.RS, 
                                     DN_num_for_split = DN_num_for_split)
        
        
        DN.centers1 <- rbind.data.frame(DN.centers1, out)
          
        if (make.pbs) setTxtProgressBar(pb,i)
        
      }
      
      # close the progress bar
      if (make.pbs) close(pb)
      
    }
  
  message("completed!")
  
  
  ## ---------------------------------------------------------
  ## Create final DN centers (reduce from DN.num1 to DN.num2)
  ## --------------------------------------------------------
  
  
  
  # create final set of data nugget centers(DN.num1 to DN.num2)
  message("Creating final set of data nugget centers...")
  
  
  DN.data <- create.DNcenters(RS = DN.centers1,
                              delete.percent = delete.percent,
                              DN.num = DN.num2,
                              dist.metric = dist.metric,
                              make.pbs = make.pbs)
  
  
  
  DN.data <- as.data.frame(DN.data)
  
  
  
  message("completed!")
  
  
  
  ## -----------------------------------------
  ## Assign observations to the nearest DN
  ## -----------------------------------------
  
  message("Assigning observations to the nearest data nuggets...")
  
  
  DN.data_mat <- as.matrix(DN.data)
  ncols_dn <- ncol(DN.data_mat)
  
  
  if (no.cores > 0){
    
    parallel::clusterExport(
      cl = cl, 
      varlist = c("dist.metric", "DN.data_mat", "assign_each"),
      envir = environment())
    
    
    
    DN.assignments <- parallel::parApply(
      cl = cl,
      X = x_mat,
      MARGIN = 1,
      FUN = function(row){
        assign_each(obs = row,
                    DN_centers_data = DN.data_mat)})
    
    }else{
      
        DN.assignments <- apply(
          X = x_mat,
          MARGIN = 1,
          FUN = function(row){
            assign_each(obs = row,
                        DN_centers_data = DN.data_mat)})
    }
  



  message("completed!")
  
  
  
  ## ---------------------------------------------------------------------------------------------------
  ## Recalculating the Data nugget centers based on the given method and calculating Data nugget scales
  ## ---------------------------------------------------------------------------------------------------
  
  
  message("Recalculating the data nugget centers and calculating the data nugget scales...")
  
  
  nugget_params <- t(sapply(1:DN.num2, function(i) {
    
    calculate_center_and_scale(i = i, DN_centers_data = DN.data_mat)
    }))
  
  if (no.cores > 0) parallel::stopCluster(cl)
  
  message("completed!")
  
  
  ## -----------------------------------------
  ## Create Data nugget information 
  ## -----------------------------------------
  
  
  # Create the data nugget information data frame
  DN.information <- data.frame("Data Nugget" = 1:DN.num2, 
                               nugget_params)

  
  # give column names for the centers and scale
  colnames(DN.information)[1:(ncols_dn + 2)] <- 
    c("Data Nugget", 
      paste0("Center", 1:ncols_dn), 
      "Scale")
  
  
  # weight for each data nugget
  DN.information$Weight = as.vector(table(DN.assignments))
  
  # Reorder
  DN.information <- DN.information[, c(
    setdiff(names(DN.information), c("Scale", "Weight")),
    "Weight",
    "Scale"
  )]
  
  
  # change the rownames
  rownames(DN.information) <- 1:DN.num2
  
  
  # replace the scales of those datanuggets having weight 1 with the (min scale)/100
  DN.information$Scale[DN.information$Weight == 1] <- min(DN.information$Scale)/100
  
  
  # create output dataframe
  output <- list("Data Nuggets" = DN.information,
                "Data Nugget Assignments" = DN.assignments)
  
  
  # assign the data nugget class to the output
  class(output) <- "datanugget"
  
  # return the data nugget
  return(output)
  
}
