create.DNcenters = function(RS,
                            
                            delete.percent,
                            
                            DN.num,
                            
                            dist.metric,
                            
                            make.pb = FALSE){
  
  
  
  # retrieve the random sample data
  
  RS.data = as.data.frame(RS)
  
  
  
  # change the rownames
  
  rownames(RS.data) = 1:nrow(RS.data)
  
  
  
  # initialize the first iteration of data nugget data
  
  DN.data = RS.data
  
  if (is.null(nrow(DN.data))){
    
    DN.data = as.data.frame(DN.data)
    
  }
  
  
  
  # initialize the total number of observations
  
  original.num = tmp.num = nrow(DN.data)
  
  
  
  # check that the DN.data is not already below the desired number of observations
  
  if (tmp.num > DN.num){
    
    # create the distance matrix for the random sample data
    
    # convert the distance matrix to a matrix
    
    # computing distance matrix could be time-consuming
    
    # using Rfast::Dist could make this faster
    
    DN.dist.matrix = as.matrix(Rfast::Dist(RS.data, method = dist.metric))
    
    
    
    # eliminate the diagonal as a choice for the minimum
    
    diag(DN.dist.matrix) = max(DN.dist.matrix) + 1
    
    
    
    # check if user wants a progress bar
    
    if (make.pb == TRUE){
      
      
      
      pb.max = nrow(RS.data) - DN.num
      
      
      
      # initialize progress bar
      
      pb = txtProgressBar(min = 0,
                          
                          max = pb.max)
      
    }
    
    
    
    # initialize the value for the progress bar
    
    for.prog = 0
    
    
    
    # initial delete indexes
    
    delete.ind <- NULL
    
    
    
    # eliminate the desired number data points that are closest together
    
    while (tmp.num > DN.num){
      
      # retrieve the number of entries to delete (i.e. 10% of the current number of entries)
      
      delete.num = tmp.num%/%(1/delete.percent)
      
      
      
      # check if the current delete.num would cause too many entries to be deleted
      
      if((tmp.num - DN.num) <= delete.num){
        
        # switch the delete num to delete just enough values
        
        delete.num = tmp.num - DN.num
        
      }
      
      
      
      # retrieve the entries for deletion
      
      if(delete.num!=0){
        
        if(delete.num%%2==0){
          
          delete_threshold = nth(DN.dist.matrix, k=delete.num, descending = FALSE)
          
          delete.entries = which(rowMins(DN.dist.matrix, value=TRUE) <= delete_threshold)
          
          
          
        } else{
          
          delete_threshold = nth(DN.dist.matrix, k=delete.num, descending = FALSE)
          
          DN.dist.matrix.rowmin = rowMins(DN.dist.matrix, value=TRUE)
          
          delete.entries = which(DN.dist.matrix.rowmin < delete_threshold)
          
          
          
          delete_entries_extract = which(DN.dist.matrix.rowmin == delete_threshold)
          
          if (length(delete_entries_extract) != 0) {
            
            delete_entries_extract = max(delete_entries_extract)
            
          }
          
          delete.entries = c(delete.entries, delete_entries_extract)
          
        }
        
      }else{
        
        delete.entries = which(rowMins(DN.dist.matrix, value=TRUE) <= min(DN.dist.matrix))
        
      }
      
      
      
      delete.entries[delete.entries == 0] = tmp.num
      
      #delete the data and dist matrix row and col related to the data
      
      #print(rownames(DN.data))
      
      DN.data = DN.data[-delete.entries, ]
      
      DN.dist.matrix = DN.dist.matrix[-delete.entries, -delete.entries]
      
      
      
      DN.data=as.data.frame(DN.data)
      
      # update the value for the progress bar
      
      for.prog = for.prog + (tmp.num - nrow(DN.data))
      
      # check if user wants a progress bar
      
      if (make.pb == TRUE){
        
        # update the progress bar
        
        setTxtProgressBar(pb, for.prog)
        
      }
      
      # store the total number of observations
      
      DN.data=as.data.frame(DN.data)
      
      tmp.num = nrow(DN.data)
      
    }
    
    
    
    if (make.pb == TRUE){ 
      
      close(pb)
      
      
      
    }
    
    
    
  }
  
  return(DN.data)
  
}

