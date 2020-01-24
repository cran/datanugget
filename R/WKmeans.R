WKmeans = function(dataset,
                   k,
                   cl.centers = NULL,
                   obs.weights,
                   num.init = 1,
                   max.iterations = 10,
                   print.progress = TRUE,
                   seed = 291102,
                   reassign.prop = .25){

  # Argument checking/fixing ####

  # make sure k is of class numeric or integer
  if (class(k) != "numeric" &
      class(k) != "integer"){

    stop('k must be of class "numeric" or "integer"')

  }

  # make sure cl.centers is of proper dimension if it is not NULL
  if (!is.null(cl.centers)){

    if (nrow(cl.centers) != k |
        ncol(cl.centers) != ncol(dataset)){

      stop('cl.centers must be a k by ncol(dataset) numeric matrix')

    }

  }

  # make sure obs.weights is of class numeric or table
  if (class(obs.weights) != "numeric" &
      class(obs.weights) != "table"){

    stop('obs.weights must be of class "numeric" or "table"')

  }

  # make sure obs.weights is length nrow(dataset)
  if (length(obs.weights) != nrow(dataset)){

    stop(paste('obs.weights must be length ',
               nrow(dataset),
               sep = ""))

  }

  # make sure num.init is of class numeric or integer
  if (class(num.init) != "numeric" &
      class(num.init) != "integer"){

    stop('num.init must be of class "numeric" or "integer"')

  }

  # make sure max.iterations is of class numeric or integer
  if (class(max.iterations) != "numeric" &
      class(max.iterations) != "integer"){

    stop('max.iterations must be of class "numeric" or "integer"')

  }

  # make sure print.progress is TRUE or FALSE
  if (class(print.progress) != "logical"){

    stop('print.progress must be TRUE OR FALSE')

  }

  # make sure seed is numeric
  if (class(seed) != "numeric"){

    stop('seed must be of class "numeric"')

  }

  # make sure reassign.prop is numeric
  if (class(reassign.prop) != "numeric"){

    stop('reassign.prop must be of class "numeric"')

  }

  # make sure reassign.prop is within (0,1]
  if (reassign.prop > 1 |
      reassign.prop <= 0){

    stop('reassign.prop must be within (0,1]')

  }

  # Function ####

  # set the random seed
  set.seed(seed)

  # convert k, num.init, and max.itertations to integers
  k = floor(k)
  num.init = floor(num.init)
  max.iterations = floor(max.iterations)

  # convert the data set to a data frame
  dataset = as.data.frame(dataset)

  # find the number of observations and the
  # number of variables in the dataset
  N = nrow(dataset)
  P = ncol(dataset)

  # initalize the weights variable
  dataset[, "Weight"] = obs.weights

  # initialize the cluster variable
  dataset[, "Cluster"] = 0

  # convert the dataset to a matrix
  dataset = as.matrix(dataset)

  # initialize list for all possible centers
  all.centers = list()

  # check if the the cluster centers have been initialized yet
  if (is.null(cl.centers)){

    # create multiple different initial centers through
    # the different possible starting points
    for (i in 1:num.init){

      # choose random cluster centers to start with
      cl.centers = dataset[sample(1:nrow(dataset),k), 1:P]

      # name the cluster centers
      rownames(cl.centers) = 1:k

      # update the list of all centers
      all.centers[[i]] = cl.centers

    }

  }else{

    all.centers[[1]] = cl.centers

    num.init = 1

  }

  # initialize the best set of clusters
  best.clusters = NULL

  # cycle through the initial cluster centers
  for (U in 1:num.init){

    # retrieve the current cluster center
    cl.centers = all.centers[[U]]

    # find the closest cluster center for the current observation
    dataset[, "Cluster"] = apply(X = dataset[, 1:P],
                                 MARGIN = 1,
                                 FUN = function(input){

                                   return(as.numeric(which.min(colSums((t(cl.centers) - input)^2))))

                                 })

    # recompute the cluster center
    for (i in 1:k){

      # check if there is more than 1 observation in the current cluster
      if (sum(dataset[, "Cluster"] == i) > 1){

        # if so, compute the means
        cl.centers[i, ] = colMeans(dataset[dataset[, "Cluster"] == i, 1:P])

      }else{

        # otherwise, the cluster center is the cluster
        cl.centers[i, ] = dataset[dataset[, "Cluster"] == i, 1:P]

      }

    }

    # initialize the number of iterations and whether a switch occurred
    A = 1
    switch.occur = TRUE

    # while the maximum number of iterations has not been reached and a switch occurred on the previous iteration
    while (A <= max.iterations &
           switch.occur == TRUE){

      # change the switch occur to false
      switch.occur = FALSE

      # cycle through the observations of the dataset
      for (i in sample(nrow(dataset),
                       ceiling(reassign.prop*nrow(dataset)),
                       replace = FALSE)){

        # retrieve the cluster of the current observation
        curr.cluster = dataset[i, "Cluster"]

        # initialize the WWCSS for each cluster
        WWCSS.vec = rep(0,k)

        # cycle through the clusters
        for (j in 1:k){

          # switch the cluster for the current observation to the current cluster
          dataset[i, "Cluster"] = j

          # calculate the WWCSS if current observation is in the curren cluster
          WWCSS.vec[j] = WWCSS(x = dataset,
                               k = k,
                               P = P)[[2]]

        }

        # check if placing the observation in the current cluster yields a greater WWCSS than
        # placing the observation in a different cluster
        if (curr.cluster != which.min(WWCSS.vec)){

          # place the observation in the better cluster
          dataset[i, "Cluster"] = which.min(WWCSS.vec)

          # note that a switch occurred
          switch.occur = TRUE

          # recompute the cluster center
          for (l in 1:k){

            # check if there is more than 1 observation in the current cluster
            if (sum(dataset[, "Cluster"] == l) > 1){

              # if so, compute the means
              cl.centers[l, ] = colMeans(dataset[dataset[, "Cluster"] == l, 1:P])

            }else{

              # otherwise, the cluster center is the cluster
              cl.centers[l, ] = dataset[dataset[, "Cluster"] == l, 1:P]

            }

          }

        }else{

          dataset[i, "Cluster"] = curr.cluster

        }

      }

      # print the progress of the algorithm if the user desires
      if (print.progress == TRUE){

        message(paste("Iteration ",
                  A,
                  " of Initialization ",
                  U,
                  ":",
                  sep = ""))

        print(table(dataset[, "Cluster"]))

      }

      # increase the iteration counter
      A = A+1

    }

    # retrieve the current WWCSS results
    WWCSS.results = WWCSS(x = dataset,
                          k = k,
                          P = P)

    # check if the best results exist yet
    if (is.null(best.clusters) == TRUE){

      # if so, create them
      best.clusters = dataset[, "Cluster"]
      best.cl.centers = cl.centers
      best.WWCSS.results = WWCSS.results

      # check if current results are better
    }else if(best.WWCSS.results[[2]] > WWCSS.results[[2]]){

      # if so, replace them
      best.clusters = dataset[, "Cluster"]
      best.cl.centers = cl.centers
      best.WWCSS.results = WWCSS.results

    }

  }

  # Create output ####

  # retrieve the final WWCSS results and clusters
  output = list(best.clusters,
                best.cl.centers,
                best.WWCSS.results)

  # give names to the output
  names(output) = c("Cluster Assignments",
                    "Cluster Centers",
                    "Weighted WCSS")

  # return the output
  return(output)

}
