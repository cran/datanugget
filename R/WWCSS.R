WWCSS = function(x,
                 k,
                 P){

  # initialize the output
  output = rep(0,k)

  # cycle through the different clusters
  for (i in 1:k){

    # retrieve the cluster observations
    clus.obs = x[x[, "Cluster"] == i, 1:P]

    # check if the cluster is non-empty
    if (!is.null(nrow(clus.obs))){

      # find the cluster center
      cl.center = colMeans(clus.obs)

      # retrieve the euclidean distances
      eucl.dist = rowSums(t(apply(X = clus.obs,
                                  MARGIN = 1,
                                  FUN = function(input){

                                    (cl.center - input)^2

                                  })))

      # find the weigehted within cluster sum of squares for the current cluster
      output[i] = t(x[x[, "Cluster"] == i, "Weight"])%*%eucl.dist

    }else{

      output[i] = 0

    }

  }

  # give the sums of squares names
  names(output) = 1:k

  # find the combined weighted within cluster sum of squares
  sum.output = sum(output)

  # return the WWCSS information
  return(list(output,
              sum.output))

}

