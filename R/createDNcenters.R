create.DNcenters = function(RS,
                            DN.num,
                            dist.metric,
                            make.pb = FALSE){

  # retrieve the random sample data
  RS.data = as.data.frame(RS)

  # change the rownames
  rownames(RS.data) = 1:nrow(RS.data)

  # initialize the first iteration of data nugget data
  DN.data = RS.data

  # initialize the total number of observations
  tmp.num = nrow(DN.data)

  # check that the DN.data is not already below the desired number of observations
  if (tmp.num > DN.num){

    # create the distance matrix for the random sample data
    DN.dist.matrix  = dist(RS.data,
                        method = dist.metric)

    # convert the distance matrix to a matrix
    DN.dist.matrix = as.matrix(DN.dist.matrix)

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

    # eliminate the desired number data points that are closest together
    while (tmp.num > DN.num){

      # retrieve the number of entries to delete (i.e. 10% of the current number of entries)
      delete.num = tmp.num%/%10

      # check if the current delete.num would cause too many entries to be deleted
      if((tmp.num - DN.num) <= delete.num){

        # switch the delete num to delete just enough values
        delete.num = tmp.num - DN.num

      }

      # retrieve the entries for deletion
      delete.entries = (sort.list(DN.dist.matrix)[1:delete.num])%%tmp.num

      # convert any 0 entries to the proper number
      delete.entries[delete.entries == 0] = tmp.num

      # delete the entries marked for deletion from the dataset and the distance matrix
      DN.data = DN.data[rownames(DN.data)[-delete.entries], ]
      DN.dist.matrix = DN.dist.matrix[-delete.entries,-delete.entries]

      if (is.null(nrow(DN.data))){

        DN.data = as.data.frame(DN.data)

      }

      # update the value for the progress bar
      for.prog = for.prog + (tmp.num - nrow(DN.data))

      # check if user wants a progress bar
      if (make.pb == TRUE){

        # update the progress bar
        setTxtProgressBar(pb, for.prog)

      }

      # store the total number of observations
      tmp.num = nrow(DN.data)

    }

    # check if user wants a progress bar
    if (make.pb == TRUE){

      # close the progress bar
      close(pb)

    }

  }

  return(DN.data)

}
