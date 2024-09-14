refine.DN = function(x,
                     DN,
                     EV.tol = .9,
                     max.splits = 5,
                     min.nugget.size = 2,
                     delta = 2,
                     seed = 291102,
                     no.cores = (detectCores() - 1),
                     make.pbs = TRUE){

  # Argument checking/fixing ####

  # bind global variable 'i'
  i = NULL

  # make sure DN is of class datanugget
  if (!inherits(DN,"datanugget")){

    stop('DN must be of class "datanugget"')

  }

  # make sure EV.tol is of class numeric
  if (!is.numeric(EV.tol)){

    stop('EV.tol must be of class "numeric"')

  }

  # make sure EV.tol is within (0,1)
  if (EV.tol >= 1 |
      EV.tol <= 0){

    stop('EV.tol must be within (0,1)')

  }

  # make sure min.nugget.size is of class numeric
  if (!is.numeric(min.nugget.size)){

    stop('min.nugget.size must be of class "numeric"')

  }

  # convert min.nugget size to an integer
  min.nugget.size = floor(min.nugget.size)

  # make sure min.nugget.size is greater than 1 is of class datanugget
  if (min.nugget.size < 2){

    stop('min.nugget.size must be a number greater than 2')

  }

  # make sure max.splits is a number
  if (!(class(max.splits) %in% c("numeric",
                               "integer"))){

    stop('max.splits must be of class "numeric" or "integer"')

  }

  # convert max.splits to an integer
  max.splits = floor(max.splits)

  # make sure seed is numeric
  if (!is.numeric(seed)){

    stop('seed must be of class "numeric"')

  }

  # make sure no.cores is numeric
  if (!is.numeric(no.cores)){

    stop('no.cores must be of class "numeric"')

  }

  # make sure make.pbs is TRUE or FALSE
  if (!is.logical(make.pbs)){

    stop('make.pbs must be TRUE OR FALSE')

  }

  # Function ####

  # set the random seed
  set.seed(seed)

  # convert no.cores to integer
  no.cores = floor(no.cores)

  # check if user wants to use parallel processing
  if (no.cores != 0){

    # create the cluster for parallel processing
    cl = makeCluster(no.cores)

    # engage the cluster for parallel processing
    registerDoSNOW(cl)

  }

  # retrieve the data nuggets and data nugget assignments
  DN.information = DN$`Data Nuggets`
  DN.assignments = DN$`Data Nugget Assignments`

  # check if the max EV of the data nugget is relevant
  if (!is.null(ncol(x))){

    # check if user wants to use parallel processing
    if (no.cores != 0){

      # find the data nuggets with the largest max EVs
      tmp.max.EVs = foreach(i = 1:nrow(DN.information),
                                 .combine = c)  %dopar%

        {

          # retrieve locations of observations assigned to this data nugget
          assign.locs = which(DN.assignments == i)

          # check that there is more than one observation in this data nugget
          if (length(assign.locs) > 1 &
              DN.information[DN.information[, "Data Nugget"] ==
                             DN.information[i, "Data Nugget"], "Scale"] > 0){

            # retrieve observations
            assign.obs = x[assign.locs, ]

            # retrieve the eigenvalues of the covariance matrix
            eigen.values = eigen(cov(assign.obs))$values

            # retrieve the top eigenvalue
            max.EV = eigen.values[1]

          }else{

            max.EV = 0

          }

          # return the eigenvalue ratio
          return(max.EV)

        }

    }else{

      # initialize the vector of data nugget max EVs
      tmp.max.EVs = NULL

      # cycle through the data nuggets
      for (i in 1:nrow(DN.information)){

        # retrieve locations of observations assigned to this data nugget
        assign.locs = which(DN.assignments == i)

        # check that there is more than one observation in this data nugget
        if (length(assign.locs) > 1 &
            DN.information[DN.information[, "Data Nugget"] ==
                           DN.information[i, "Data Nugget"], "Scale"] > 0){

          # retrieve observations
          assign.obs = x[assign.locs, ]

          # retrieve the eigenvalues of the covariance matrix
          eigen.values = eigen(cov(assign.obs))$values

          # retrieve the top eigenvalue
          max.EV = eigen.values[1]

        }else{

          max.EV = 0

        }

        # return the eigenvalue ratio
        tmp.max.EVs = c(tmp.max.EVs,
                             max.EV)

      }

    }

    # retrieve the original threshold for shape
    orig.threshold = as.numeric(quantile(tmp.max.EVs[tmp.max.EVs!= 0],
                                         probs = EV.tol))

    # retrieve the data nuggets with max EVs larger than the max EV tolerance
    large.max.EVs = which(tmp.max.EVs > orig.threshold)

    # retrieve the amount of data nuggets with large shapes
    cur.LS.num = length(large.max.EVs)

    # initialize new amount of data nuggets with large shapes
    new.LS.num = cur.LS.num + 1

    # initialize the number of split attempts
    split.attempts = 0

    # initialize progress bar
    message("splitting data nuggets according to max EV...")

    # perform the following operations while the number of data nuggets with large max EVs is
    # both > 0 and different from the previous number
    while(cur.LS.num > 0 &
          new.LS.num > 0 &
          cur.LS.num != new.LS.num &
          split.attempts < max.splits){

      # increase the number of split attempts
      split.attempts = split.attempts + 1

      # print split attempt
      message(paste("split #",
                    split.attempts,
                    sep = ""))

      # reset the current ammount of data nuggets with large shapes
      cur.LS.num = new.LS.num

      # initialize shift for new nuggets
      shift = max(DN.information[, "Data Nugget"])

      # check if user wants to use parallel processing
      if (no.cores != 0){

        # check if user wants a progress bar
        if (make.pbs == TRUE){

          # initialize progress bar
          pb = txtProgressBar(min = 0,
                              max = length(large.max.EVs))

          # update the progress bar
          progress = function(n){setTxtProgressBar(pb, n)}
          opts = list(progress = progress)

          # split data nuggets with large shapes
          tmp.new.DN.information = foreach(i = large.max.EVs,
                                           .combine = rbind,
                                           .options.snow = opts)  %dopar%

            {

              # retrieve the temporary assignments
              tmp.DN.assignments = DN.assignments

              # retrieve locations of observations assigned to this data nugget
              assign.locs = which(tmp.DN.assignments == i)

              # retrieve observations
              assign.obs = x[assign.locs, ]

              # check if the current data nugget can be split
              if (nrow(assign.obs) >= 2*min.nugget.size){

                # check if there are more than 2 observations assigned to this nugget
                if (nrow(assign.obs) > 2){

                  # if so, split the nugget into two nuggets with k means
                  new.nugget.info = kmeans(x = assign.obs,
                                           centers = 2,
                                           iter.max = 1000,
                                           nstart = 25)

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = shift + 2*(which(large.max.EVs == i)-1) + new.nugget.info$cluster

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.max.EVs == i)-1) + 1:2

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2:(ncol(x)+1)] = new.nugget.info$centers[1, ]
                  new.DN.information[2, 2:(ncol(x)+1)] = new.nugget.info$centers[2, ]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                        sum(new.nugget.info$cluster == 2))

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # retrieve the shape for the new data nuggets
                  if (length(new.assign.locs1) > 1){

                    if (ncol(x) == 1){

                      new.DN.information[1, "Scale"] = var(x[new.assign.locs1, ])

                    }else{

                      new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[1, "Scale"] = 0

                  }

                  if (length(new.assign.locs2) > 1){

                    if (ncol(x) == 1){

                      new.DN.information[2, "Scale"] = var(x[new.assign.locs2, ])

                    }else{

                      new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[2, "Scale"] = 0

                  }

                  # check if either of the new nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }else{

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2:(ncol(x)+1)] = assign.obs[1, ]
                  new.DN.information[2, 2:(ncol(x)+1)] = assign.obs[2, ]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = 1

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # check if the new data nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }

                return(new.DN.information)

              }

            }

          # close the progress bar
          close(pb)

        }else{

          # split data nuggets with large shapes
          tmp.new.DN.information = foreach(i = large.max.EVs,
                                           .combine = rbind)  %dopar%

            {

              # retrieve the temporary assignments
              tmp.DN.assignments = DN.assignments

              # retrieve locations of observations assigned to this data nugget
              assign.locs = which(tmp.DN.assignments == i)

              # retrieve observations
              assign.obs = x[assign.locs, ]

              # check if the current data nugget can be split
              if (nrow(assign.obs) >= 2*min.nugget.size){

                # check if there are more than 2 observations assigned to this nugget
                if (nrow(assign.obs) > 2){

                  # if so, split the nugget into two nuggets with k means
                  new.nugget.info = kmeans(x = assign.obs,
                                           centers = 2,
                                           iter.max = 1000,
                                           nstart = 25)

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = shift + 2*(which(large.max.EVs == i)-1) + new.nugget.info$cluster

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.max.EVs == i)-1) + 1:2

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2:(ncol(x)+1)] = new.nugget.info$centers[1, ]
                  new.DN.information[2, 2:(ncol(x)+1)] = new.nugget.info$centers[2, ]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                        sum(new.nugget.info$cluster == 2))

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # retrieve the shape for the new data nuggets
                  if (length(new.assign.locs1) > 1){

                    if (ncol(x) == 1){

                      new.DN.information[1, "Scale"] = var(x[new.assign.locs1, ])

                    }else{

                      new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[1, "Scale"] = 0

                  }

                  if (length(new.assign.locs2) > 1){

                    if (ncol(x) == 1){

                      new.DN.information[2, "Scale"] = var(x[new.assign.locs2, ])

                    }else{

                      new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[2, "Scale"] = 0

                  }

                  # check if either of the new nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }else{

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2:(ncol(x)+1)] = assign.obs[1, ]
                  new.DN.information[2, 2:(ncol(x)+1)] = assign.obs[2, ]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = 1

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # check if the new data nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }

                return(new.DN.information)

              }

            }

        }

      }else{

        # check if user wants a progress bar
        if (make.pbs == TRUE){

          # initialize progress bar
          pb = txtProgressBar(min = 0,
                              max = length(large.max.EVs))

          # initialize data for split data nuggets
          tmp.new.DN.information = NULL

          # cycle through the data nuggets with large max EVs
          for (i in large.max.EVs){



            # retrieve the temporary assignments
            tmp.DN.assignments = DN.assignments

            # retrieve locations of observations assigned to this data nugget
            assign.locs = which(tmp.DN.assignments == i)

            # retrieve observations
            assign.obs = x[assign.locs, ]

            # check if the current data nugget can be split
            if (nrow(assign.obs) >= 2*min.nugget.size){

              # check if there are more than 2 observations assigned to this nugget
              if (nrow(assign.obs) > 2){

                # if so, split the nugget into two nuggets with k means
                new.nugget.info = kmeans(x = assign.obs,
                                         centers = 2,
                                         iter.max = 1000,
                                         nstart = 25)

                # reassign the observations to their new nuggets
                tmp.DN.assignments[assign.locs] = shift + 2*(which(large.max.EVs == i)-1) + new.nugget.info$cluster

                # initialize data frame for new data nuggets
                new.DN.information = DN.information[1:2, ]

                # create new data nugget names
                new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.max.EVs == i)-1) + 1:2

                # retrieve the new data nugget centers
                new.DN.information[1, 2:(ncol(x)+1)] = new.nugget.info$centers[1, ]
                new.DN.information[2, 2:(ncol(x)+1)] = new.nugget.info$centers[2, ]

                # retrieve the weights for these data nuggets
                new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                      sum(new.nugget.info$cluster == 2))

                # retrieve the assignments for the new data nuggets
                new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                # retrieve the shape for the new data nuggets
                if (length(new.assign.locs1) > 1){

                  if (ncol(x) == 1){

                    new.DN.information[1, "Scale"] = var(x[new.assign.locs1, ])

                  }else{

                    new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                  }

                }else{

                  new.DN.information[1, "Scale"] = 0

                }

                if (length(new.assign.locs2) > 1){

                  if (ncol(x) == 1){

                    new.DN.information[2, "Scale"] = var(x[new.assign.locs2, ])

                  }else{

                    new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                  }

                }else{

                  new.DN.information[2, "Scale"] = 0

                }

                # check if either of the new nuggets meet the minimum nugget size criteria
                if (length(new.assign.locs1) >= min.nugget.size &
                    length(new.assign.locs2) >= min.nugget.size){

                  new.DN.information[, "Delete"] =
                    DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                  new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                            collapse = ",")

                  new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                            collapse = ",")

                }else{

                  new.DN.information = NULL

                }

              }else{

                # reassign the observations to their new nuggets
                tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                # initialize data frame for new data nuggets
                new.DN.information = DN.information[1:2, ]

                # create new data nugget names
                new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                # retrieve the new data nugget centers
                new.DN.information[1, 2:(ncol(x)+1)] = assign.obs[1, ]
                new.DN.information[2, 2:(ncol(x)+1)] = assign.obs[2, ]

                # retrieve the weights for these data nuggets
                new.DN.information[1:2, "Weight"] = 1

                # retrieve the assignments for the new data nuggets
                new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                # check if the new data nuggets meet the minimum nugget size criteria
                if (length(new.assign.locs1) >= min.nugget.size &
                    length(new.assign.locs2) >= min.nugget.size){

                  new.DN.information[, "Delete"] =
                    DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                  new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                            collapse = ",")

                  new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                            collapse = ",")

                }else{

                  new.DN.information = NULL

                }

              }

              rbind(tmp.new.DN.information,
                    new.DN.information)

            }

            setTxtProgressBar(pb, i)

          }

          # close the progress bar
          close(pb)

        }else{

          # initialize data for split data nuggets
          tmp.new.DN.information = NULL

          # cycle through the data nuggets with large max EVs
          for (i in large.max.EVs){



            # retrieve the temporary assignments
            tmp.DN.assignments = DN.assignments

            # retrieve locations of observations assigned to this data nugget
            assign.locs = which(tmp.DN.assignments == i)

            # retrieve observations
            assign.obs = x[assign.locs, ]

            # check if the current data nugget can be split
            if (nrow(assign.obs) >= 2*min.nugget.size){

              # check if there are more than 2 observations assigned to this nugget
              if (nrow(assign.obs) > 2){

                # if so, split the nugget into two nuggets with k means
                new.nugget.info = kmeans(x = assign.obs,
                                         centers = 2,
                                         iter.max = 1000,
                                         nstart = 25)

                # reassign the observations to their new nuggets
                tmp.DN.assignments[assign.locs] = shift + 2*(which(large.max.EVs == i)-1) + new.nugget.info$cluster

                # initialize data frame for new data nuggets
                new.DN.information = DN.information[1:2, ]

                # create new data nugget names
                new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.max.EVs == i)-1) + 1:2

                # retrieve the new data nugget centers
                new.DN.information[1, 2:(ncol(x)+1)] = new.nugget.info$centers[1, ]
                new.DN.information[2, 2:(ncol(x)+1)] = new.nugget.info$centers[2, ]

                # retrieve the weights for these data nuggets
                new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                      sum(new.nugget.info$cluster == 2))

                # retrieve the assignments for the new data nuggets
                new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                # retrieve the shape for the new data nuggets
                if (length(new.assign.locs1) > 1){

                  if (ncol(x) == 1){

                    new.DN.information[1, "Scale"] = var(x[new.assign.locs1, ])

                  }else{

                    new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                  }

                }else{

                  new.DN.information[1, "Scale"] = 0

                }

                if (length(new.assign.locs2) > 1){

                  if (ncol(x) == 1){

                    new.DN.information[2, "Scale"] = var(x[new.assign.locs2, ])

                  }else{

                    new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                  }

                }else{

                  new.DN.information[2, "Scale"] = 0

                }

                # check if either of the new nuggets meet the minimum nugget size criteria
                if (length(new.assign.locs1) >= min.nugget.size &
                    length(new.assign.locs2) >= min.nugget.size){

                  new.DN.information[, "Delete"] =
                    DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                  new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                            collapse = ",")

                  new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                            collapse = ",")

                }else{

                  new.DN.information = NULL

                }

              }else{

                # reassign the observations to their new nuggets
                tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                # initialize data frame for new data nuggets
                new.DN.information = DN.information[1:2, ]

                # create new data nugget names
                new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                # retrieve the new data nugget centers
                new.DN.information[1, 2:(ncol(x)+1)] = assign.obs[1, ]
                new.DN.information[2, 2:(ncol(x)+1)] = assign.obs[2, ]

                # retrieve the weights for these data nuggets
                new.DN.information[1:2, "Weight"] = 1

                # retrieve the assignments for the new data nuggets
                new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                # check if the new data nuggets meet the minimum nugget size criteria
                if (length(new.assign.locs1) >= min.nugget.size &
                    length(new.assign.locs2) >= min.nugget.size){

                  new.DN.information[, "Delete"] =
                    DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                  new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                            collapse = ",")

                  new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                            collapse = ",")

                }else{

                  new.DN.information = NULL

                }

              }

              rbind(tmp.new.DN.information,
                    new.DN.information)

            }

          }

        }

      }

      # check if any data nuggets were split
      if (!is.null(tmp.new.DN.information)){

        # bind the new data nuggets with the old data nuggets
        DN.information = rbind.data.frame(DN.information,
                                          tmp.new.DN.information[, -c(ncol(tmp.new.DN.information)-1,
                                                                      ncol(tmp.new.DN.information))])

        # delete the data nuggets that were split
        DN.information = DN.information[!(DN.information[, "Data Nugget"] %in% unique(tmp.new.DN.information[, "Delete"])), ]

        # cycle through the new data nuggets
        for (j in 1:nrow(tmp.new.DN.information)){

          # reassign the data nuggets
          DN.assignments[as.numeric(unlist(strsplit(tmp.new.DN.information[j, "New Assignment Locations"],
                                                    split = ",")))] = tmp.new.DN.information[j, "Data Nugget"]

        }

      }

      # check if user wants to use parallel processing
      if (no.cores != 0){

        # find the data nuggets with the largest max EVs
        tmp.max.EVs = foreach(i = 1:nrow(DN.information),
                                   .combine = c)  %dopar%

          {

            # retrieve locations of observations assigned to this data nugget
            assign.locs = which(DN.assignments == i)

            # check that there is more than one observation in this data nugget
            if (length(assign.locs) > 1 &
                DN.information[DN.information[, "Data Nugget"] ==
                               DN.information[i, "Data Nugget"], "Scale"] > 0){

              # retrieve observations
              assign.obs = x[assign.locs, ]

              # retrieve the eigenvalues of the covariance matrix
              eigen.values = eigen(cov(assign.obs))$values

              # retrieve the top eigenvalue
              max.EV = eigen.values[1]

            }else{

              max.EV = 0

            }

            # return the eigenvalue ratio
            return(max.EV)

          }

      }else{

        # initialize the vector of data nugget max EVs
        tmp.max.EVs = NULL

        # cycle through the data nuggets
        for (i in 1:nrow(DN.information)){

          # retrieve locations of observations assigned to this data nugget
          assign.locs = which(DN.assignments == i)

          # check that there is more than one observation in this data nugget
          if (length(assign.locs) > 1 &
              DN.information[DN.information[, "Data Nugget"] ==
                             DN.information[i, "Data Nugget"], "Scale"] > 0){

            # retrieve observations
            assign.obs = x[assign.locs, ]

            # retrieve the eigenvalues of the covariance matrix
            eigen.values = eigen(cov(assign.obs))$values

            # retrieve the top eigenvalue
            max.EV = eigen.values[1]

          }else{

            max.EV = 0

          }

          # return the eigenvalue ratio
          tmp.max.EVs = c(tmp.max.EVs,
                               max.EV)

        }

      }

      # retrieve the data nuggets with max EVs larger than the max EV tolerance
      large.max.EVs = which(tmp.max.EVs > orig.threshold)

      # retrieve the amount of data nuggets with large shapes
      new.LS.num = length(large.max.EVs)

    }

    message("complete!")

    # rename the rows of DN information
    rownames(DN.information) = 1:nrow(DN.information)

    # rename the data nuggets
    DN.information[, "Data Nugget"] = 1:nrow(DN.information)

    # retrieve list of data nuggets
    DN.list = unique(DN.assignments)[order(unique(DN.assignments))]

    # cycle through the list of data nuggets
    for (j in DN.list){

      # change the data nugget number
      DN.assignments[DN.assignments == j] = which(DN.list[order(DN.list)] == j)

    }

  }else{

    # retrieve the original threshold for shape
    orig.threshold = as.numeric(quantile(DN.information[DN.information[, "Scale"] != 0, "Scale"],
                                         probs = EV.tol))

    # retrieve the data nuggets with shapes larger than the shape tolerance
    large.shapes = DN.information[which(DN.information[, "Scale"] > orig.threshold), "Data Nugget"]

    # retrieve the amount of data nuggets with large shapes
    cur.LS.num = length(large.shapes)

    # initialize new amount of data nuggets with large shapes
    new.LS.num = cur.LS.num + 1

    # initialize the number of split attempts
    split.attempts = 0

    # initialize progress bar
    message("splitting data nuggets according to shape...")

    # perform the following operations while the number of data nuggets with large shapes is
    # both > 0 and different from the previous number
    while(cur.LS.num > 0 &
          new.LS.num > 0 &
          cur.LS.num != new.LS.num &
          split.attempts < max.splits){

      # increase the number of split attempts
      split.attempts = split.attempts + 1

      # print split attempt
      message(paste("split #",
                    split.attempts,
                    sep = ""))

      # reset the current ammount of data nuggets with large shapes
      cur.LS.num = new.LS.num

      # initialize shift for new nuggets
      shift = max(DN.information[, "Data Nugget"])

      # check if user wants to use parallel processing
      if (no.cores != 0){

        # check if user wants a progress bar
        if (make.pbs == TRUE){

          # initialize progress bar
          pb = txtProgressBar(min = 0,
                              max = length(large.shapes))

          # update the progress bar
          progress = function(n){setTxtProgressBar(pb, n)}
          opts = list(progress = progress)

          # split data nuggets with large shapes
          tmp.new.DN.information = foreach(i = large.shapes,
                                           .combine = rbind,
                                           .options.snow = opts)  %dopar%

            {

              # retrieve the temporary assignments
              tmp.DN.assignments = DN.assignments

              # retrieve locations of observations assigned to this data nugget
              assign.locs = which(tmp.DN.assignments == i)

              if (!is.null(ncol(x))){

                # retrieve observations
                assign.obs = x[assign.locs, ]

                # check if the current data nugget can be split
                if (nrow(assign.obs) >= 2*min.nugget.size){

                  # check if there are more than 2 observations assigned to this nugget
                  if (nrow(assign.obs) > 2){

                    # if so, split the nugget into two nuggets with k means
                    new.nugget.info = kmeans(x = assign.obs,
                                             centers = 2,
                                             iter.max = 1000,
                                             nstart = 25)

                    # reassign the observations to their new nuggets
                    tmp.DN.assignments[assign.locs] = shift + 2*(which(large.shapes == i)-1) + new.nugget.info$cluster

                    # initialize data frame for new data nuggets
                    new.DN.information = DN.information[1:2, ]

                    # create new data nugget names
                    new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.shapes == i)-1) + 1:2

                    # retrieve the new data nugget centers
                    new.DN.information[1, 2:(ncol(x)+1)] = new.nugget.info$centers[1, ]
                    new.DN.information[2, 2:(ncol(x)+1)] = new.nugget.info$centers[2, ]

                    # retrieve the weights for these data nuggets
                    new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                          sum(new.nugget.info$cluster == 2))

                    # retrieve the assignments for the new data nuggets
                    new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                    new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                    # retrieve the shape for the new data nuggets
                    if (length(new.assign.locs1) > 1){

                      if (ncol(x) == 1){

                        new.DN.information[1, "Scale"] = var(x[new.assign.locs1, ])

                      }else{

                        new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                      }

                    }else{

                      new.DN.information[1, "Scale"] = 0

                    }

                    if (length(new.assign.locs2) > 1){

                      if (ncol(x) == 1){

                        new.DN.information[2, "Scale"] = var(x[new.assign.locs2, ])

                      }else{

                        new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                      }

                    }else{

                      new.DN.information[2, "Scale"] = 0

                    }

                    # check if either of the new nuggets meet the minimum nugget size criteria
                    if (length(new.assign.locs1) >= min.nugget.size &
                        length(new.assign.locs2) >= min.nugget.size){

                      new.DN.information[, "Delete"] =
                        DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                      new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                                collapse = ",")

                      new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                                collapse = ",")

                    }else{

                      new.DN.information = NULL

                    }

                  }else{

                    # reassign the observations to their new nuggets
                    tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                    # initialize data frame for new data nuggets
                    new.DN.information = DN.information[1:2, ]

                    # create new data nugget names
                    new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                    # retrieve the new data nugget centers
                    new.DN.information[1, 2:(ncol(x)+1)] = assign.obs[1, ]
                    new.DN.information[2, 2:(ncol(x)+1)] = assign.obs[2, ]

                    # retrieve the weights for these data nuggets
                    new.DN.information[1:2, "Weight"] = 1

                    # retrieve the assignments for the new data nuggets
                    new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                    new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                    # check if the new data nuggets meet the minimum nugget size criteria
                    if (length(new.assign.locs1) >= min.nugget.size &
                        length(new.assign.locs2) >= min.nugget.size){

                      new.DN.information[, "Delete"] =
                        DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                      new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                                collapse = ",")

                      new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                                collapse = ",")

                    }else{

                      new.DN.information = NULL

                    }

                  }

                  return(new.DN.information)

                }

              }else{

                # retrieve observations
                assign.obs = x[assign.locs]

                # check if the current data nugget can be split
                if (length(assign.obs) >= 2*min.nugget.size){

                  # check if there are more than 2 observations assigned to this nugget
                  if (length(assign.obs) > 2){

                    # if so, split the nugget into two nuggets with k means
                    new.nugget.info = kmeans(x = assign.obs,
                                             centers = 2,
                                             iter.max = 1000,
                                             nstart = 25)

                    # reassign the observations to their new nuggets
                    tmp.DN.assignments[assign.locs] = shift + 2*(which(large.shapes == i)-1) + new.nugget.info$cluster

                    # initialize data frame for new data nuggets
                    new.DN.information = DN.information[1:2, ]

                    # create new data nugget names
                    new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.shapes == i)-1) + 1:2

                    # retrieve the new data nugget centers
                    new.DN.information[1, 2] = new.nugget.info$centers[1]
                    new.DN.information[2, 2] = new.nugget.info$centers[2]

                    # retrieve the weights for these data nuggets
                    new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                          sum(new.nugget.info$cluster == 2))

                    # retrieve the assignments for the new data nuggets
                    new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                    new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                    # retrieve the shape for the new data nuggets
                    if (length(new.assign.locs1) > 1){

                      if (is.null(ncol(x))){

                        new.DN.information[1, "Scale"] = var(x[new.assign.locs1])

                      }else{

                        new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                      }

                    }else{

                      new.DN.information[1, "Scale"] = 0

                    }

                    if (length(new.assign.locs2) > 1){

                      if (is.null(ncol(x))){

                        new.DN.information[2, "Scale"] = var(x[new.assign.locs2])

                      }else{

                        new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                      }

                    }else{

                      new.DN.information[2, "Scale"] = 0

                    }

                    # check if either of the new nuggets meet the minimum nugget size criteria
                    if (length(new.assign.locs1) >= min.nugget.size &
                        length(new.assign.locs2) >= min.nugget.size){

                      new.DN.information[, "Delete"] =
                        DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                      new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                                collapse = ",")

                      new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                                collapse = ",")

                    }else{

                      new.DN.information = NULL

                    }

                  }else{

                    # reassign the observations to their new nuggets
                    tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                    # initialize data frame for new data nuggets
                    new.DN.information = DN.information[1:2, ]

                    # create new data nugget names
                    new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                    # retrieve the new data nugget centers
                    new.DN.information[1, 2] = assign.obs[1]
                    new.DN.information[2, 2] = assign.obs[2]

                    # retrieve the weights for these data nuggets
                    new.DN.information[1:2, "Weight"] = 1

                    # retrieve the assignments for the new data nuggets
                    new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                    new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                    # check if the new data nuggets meet the minimum nugget size criteria
                    if (length(new.assign.locs1) >= min.nugget.size &
                        length(new.assign.locs2) >= min.nugget.size){

                      new.DN.information[, "Delete"] =
                        DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                      new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                                collapse = ",")

                      new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                                collapse = ",")

                    }else{

                      new.DN.information = NULL

                    }

                  }

                  return(new.DN.information)

                }

              }

            }

          # close the progress bar
          close(pb)

        }else{

          # split data nuggets with large shapes
          tmp.new.DN.information = foreach(i = large.shapes,
                                           .combine = rbind)  %dopar%

            {

              # retrieve the temporary assignments
              tmp.DN.assignments = DN.assignments

              # retrieve locations of observations assigned to this data nugget
              assign.locs = which(tmp.DN.assignments == i)

              if (!is.null(ncol(x))){

                # retrieve observations
                assign.obs = x[assign.locs, ]

                # check if the current data nugget can be split
                if (nrow(assign.obs) >= 2*min.nugget.size){

                  # check if there are more than 2 observations assigned to this nugget
                  if (nrow(assign.obs) > 2){

                    # if so, split the nugget into two nuggets with k means
                    new.nugget.info = kmeans(x = assign.obs,
                                             centers = 2,
                                             iter.max = 1000,
                                             nstart = 25)

                    # reassign the observations to their new nuggets
                    tmp.DN.assignments[assign.locs] = shift + 2*(which(large.shapes == i)-1) + new.nugget.info$cluster

                    # initialize data frame for new data nuggets
                    new.DN.information = DN.information[1:2, ]

                    # create new data nugget names
                    new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.shapes == i)-1) + 1:2

                    # retrieve the new data nugget centers
                    new.DN.information[1, 2:(ncol(x)+1)] = new.nugget.info$centers[1, ]
                    new.DN.information[2, 2:(ncol(x)+1)] = new.nugget.info$centers[2, ]

                    # retrieve the weights for these data nuggets
                    new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                          sum(new.nugget.info$cluster == 2))

                    # retrieve the assignments for the new data nuggets
                    new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                    new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                    # retrieve the shape for the new data nuggets
                    if (length(new.assign.locs1) > 1){

                      if (ncol(x) == 1){

                        new.DN.information[1, "Scale"] = var(x[new.assign.locs1, ])

                      }else{

                        new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                      }

                    }else{

                      new.DN.information[1, "Scale"] = 0

                    }

                    if (length(new.assign.locs2) > 1){

                      if (ncol(x) == 1){

                        new.DN.information[2, "Scale"] = var(x[new.assign.locs2, ])

                      }else{

                        new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                      }

                    }else{

                      new.DN.information[2, "Scale"] = 0

                    }

                    # check if either of the new nuggets meet the minimum nugget size criteria
                    if (length(new.assign.locs1) >= min.nugget.size &
                        length(new.assign.locs2) >= min.nugget.size){

                      new.DN.information[, "Delete"] =
                        DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                      new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                                collapse = ",")

                      new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                                collapse = ",")

                    }else{

                      new.DN.information = NULL

                    }

                  }else{

                    # reassign the observations to their new nuggets
                    tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                    # initialize data frame for new data nuggets
                    new.DN.information = DN.information[1:2, ]

                    # create new data nugget names
                    new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                    # retrieve the new data nugget centers
                    new.DN.information[1, 2:(ncol(x)+1)] = assign.obs[1, ]
                    new.DN.information[2, 2:(ncol(x)+1)] = assign.obs[2, ]

                    # retrieve the weights for these data nuggets
                    new.DN.information[1:2, "Weight"] = 1

                    # retrieve the assignments for the new data nuggets
                    new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                    new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                    # check if the new data nuggets meet the minimum nugget size criteria
                    if (length(new.assign.locs1) >= min.nugget.size &
                        length(new.assign.locs2) >= min.nugget.size){

                      new.DN.information[, "Delete"] =
                        DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                      new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                                collapse = ",")

                      new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                                collapse = ",")

                    }else{

                      new.DN.information = NULL

                    }

                  }

                  return(new.DN.information)

                }

              }else{

                # retrieve observations
                assign.obs = x[assign.locs]

                # check if the current data nugget can be split
                if (length(assign.obs) >= 2*min.nugget.size){

                  # check if there are more than 2 observations assigned to this nugget
                  if (length(assign.obs) > 2){

                    # if so, split the nugget into two nuggets with k means
                    new.nugget.info = kmeans(x = assign.obs,
                                             centers = 2,
                                             iter.max = 1000,
                                             nstart = 25)

                    # reassign the observations to their new nuggets
                    tmp.DN.assignments[assign.locs] = shift + 2*(which(large.shapes == i)-1) + new.nugget.info$cluster

                    # initialize data frame for new data nuggets
                    new.DN.information = DN.information[1:2, ]

                    # create new data nugget names
                    new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.shapes == i)-1) + 1:2

                    # retrieve the new data nugget centers
                    new.DN.information[1, 2] = new.nugget.info$centers[1]
                    new.DN.information[2, 2] = new.nugget.info$centers[2]

                    # retrieve the weights for these data nuggets
                    new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                          sum(new.nugget.info$cluster == 2))

                    # retrieve the assignments for the new data nuggets
                    new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                    new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                    # retrieve the shape for the new data nuggets
                    if (length(new.assign.locs1) > 1){

                      if (is.null(ncol(x))){

                        new.DN.information[1, "Scale"] = var(x[new.assign.locs1])

                      }else{

                        new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                      }

                    }else{

                      new.DN.information[1, "Scale"] = 0

                    }

                    if (length(new.assign.locs2) > 1){

                      if (is.null(ncol(x))){

                        new.DN.information[2, "Scale"] = var(x[new.assign.locs2])

                      }else{

                        new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                      }

                    }else{

                      new.DN.information[2, "Scale"] = 0

                    }

                    # check if either of the new nuggets meet the minimum nugget size criteria
                    if (length(new.assign.locs1) >= min.nugget.size &
                        length(new.assign.locs2) >= min.nugget.size){

                      new.DN.information[, "Delete"] =
                        DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                      new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                                collapse = ",")

                      new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                                collapse = ",")

                    }else{

                      new.DN.information = NULL

                    }

                  }else{

                    # reassign the observations to their new nuggets
                    tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                    # initialize data frame for new data nuggets
                    new.DN.information = DN.information[1:2, ]

                    # create new data nugget names
                    new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                    # retrieve the new data nugget centers
                    new.DN.information[1, 2] = assign.obs[1]
                    new.DN.information[2, 2] = assign.obs[2]

                    # retrieve the weights for these data nuggets
                    new.DN.information[1:2, "Weight"] = 1

                    # retrieve the assignments for the new data nuggets
                    new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                    new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                    # check if the new data nuggets meet the minimum nugget size criteria
                    if (length(new.assign.locs1) >= min.nugget.size &
                        length(new.assign.locs2) >= min.nugget.size){

                      new.DN.information[, "Delete"] =
                        DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                      new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                                collapse = ",")

                      new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                                collapse = ",")

                    }else{

                      new.DN.information = NULL

                    }

                  }

                  return(new.DN.information)

                }

              }

            }

        }

      }else{

        # check if user wants a progress bar
        if (make.pbs == TRUE){

          # initialize progress bar
          pb = txtProgressBar(min = 0,
                              max = length(large.shapes))

          # initialize the information for the split data nugget
          tmp.new.DN.information = NULL

          # cycle through the data nuggets with large shapes
          for (i in large.shapes){

            # retrieve the temporary assignments
            tmp.DN.assignments = DN.assignments

            # retrieve locations of observations assigned to this data nugget
            assign.locs = which(tmp.DN.assignments == i)

            if (!is.null(ncol(x))){

              # retrieve observations
              assign.obs = x[assign.locs, ]

              # check if the current data nugget can be split
              if (nrow(assign.obs) >= 2*min.nugget.size){

                # check if there are more than 2 observations assigned to this nugget
                if (nrow(assign.obs) > 2){

                  # if so, split the nugget into two nuggets with k means
                  new.nugget.info = kmeans(x = assign.obs,
                                           centers = 2,
                                           iter.max = 1000,
                                           nstart = 25)

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = shift + 2*(which(large.shapes == i)-1) + new.nugget.info$cluster

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.shapes == i)-1) + 1:2

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2:(ncol(x)+1)] = new.nugget.info$centers[1, ]
                  new.DN.information[2, 2:(ncol(x)+1)] = new.nugget.info$centers[2, ]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                        sum(new.nugget.info$cluster == 2))

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # retrieve the shape for the new data nuggets
                  if (length(new.assign.locs1) > 1){

                    if (ncol(x) == 1){

                      new.DN.information[1, "Scale"] = var(x[new.assign.locs1, ])

                    }else{

                      new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[1, "Scale"] = 0

                  }

                  if (length(new.assign.locs2) > 1){

                    if (ncol(x) == 1){

                      new.DN.information[2, "Scale"] = var(x[new.assign.locs2, ])

                    }else{

                      new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[2, "Scale"] = 0

                  }

                  # check if either of the new nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }else{

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2:(ncol(x)+1)] = assign.obs[1, ]
                  new.DN.information[2, 2:(ncol(x)+1)] = assign.obs[2, ]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = 1

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # check if the new data nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }

                rbind(tmp.new.DN.information,
                      new.DN.information)

              }

            }else{

              # retrieve observations
              assign.obs = x[assign.locs]

              # check if the current data nugget can be split
              if (length(assign.obs) >= 2*min.nugget.size){

                # check if there are more than 2 observations assigned to this nugget
                if (length(assign.obs) > 2){

                  # if so, split the nugget into two nuggets with k means
                  new.nugget.info = kmeans(x = assign.obs,
                                           centers = 2,
                                           iter.max = 1000,
                                           nstart = 25)

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = shift + 2*(which(large.shapes == i)-1) + new.nugget.info$cluster

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.shapes == i)-1) + 1:2

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2] = new.nugget.info$centers[1]
                  new.DN.information[2, 2] = new.nugget.info$centers[2]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                        sum(new.nugget.info$cluster == 2))

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # retrieve the shape for the new data nuggets
                  if (length(new.assign.locs1) > 1){

                    if (is.null(ncol(x))){

                      new.DN.information[1, "Scale"] = var(x[new.assign.locs1])

                    }else{

                      new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[1, "Scale"] = 0

                  }

                  if (length(new.assign.locs2) > 1){

                    if (is.null(ncol(x))){

                      new.DN.information[2, "Scale"] = var(x[new.assign.locs2])

                    }else{

                      new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[2, "Scale"] = 0

                  }

                  # check if either of the new nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }else{

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2] = assign.obs[1]
                  new.DN.information[2, 2] = assign.obs[2]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = 1

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # check if the new data nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }

                rbind(tmp.new.DN.information,
                      new.DN.information)

              }

            }

            setTxtProgressBar(pb, i)

          }

          # close the progress bar
          close(pb)

        }else{

          # initialize the information for the split data nugget
          tmp.new.DN.information = NULL

          # cycle through the data nuggets with large shapes
          for (i in large.shapes){

            # retrieve the temporary assignments
            tmp.DN.assignments = DN.assignments

            # retrieve locations of observations assigned to this data nugget
            assign.locs = which(tmp.DN.assignments == i)

            if (!is.null(ncol(x))){

              # retrieve observations
              assign.obs = x[assign.locs, ]

              # check if the current data nugget can be split
              if (nrow(assign.obs) >= 2*min.nugget.size){

                # check if there are more than 2 observations assigned to this nugget
                if (nrow(assign.obs) > 2){

                  # if so, split the nugget into two nuggets with k means
                  new.nugget.info = kmeans(x = assign.obs,
                                           centers = 2,
                                           iter.max = 1000,
                                           nstart = 25)

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = shift + 2*(which(large.shapes == i)-1) + new.nugget.info$cluster

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.shapes == i)-1) + 1:2

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2:(ncol(x)+1)] = new.nugget.info$centers[1, ]
                  new.DN.information[2, 2:(ncol(x)+1)] = new.nugget.info$centers[2, ]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                        sum(new.nugget.info$cluster == 2))

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # retrieve the shape for the new data nuggets
                  if (length(new.assign.locs1) > 1){

                    if (ncol(x) == 1){

                      new.DN.information[1, "Scale"] = var(x[new.assign.locs1, ])

                    }else{

                      new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[1, "Scale"] = 0

                  }

                  if (length(new.assign.locs2) > 1){

                    if (ncol(x) == 1){

                      new.DN.information[2, "Scale"] = var(x[new.assign.locs2, ])

                    }else{

                      new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[2, "Scale"] = 0

                  }

                  # check if either of the new nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }else{

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2:(ncol(x)+1)] = assign.obs[1, ]
                  new.DN.information[2, 2:(ncol(x)+1)] = assign.obs[2, ]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = 1

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # check if the new data nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }

                rbind(tmp.new.DN.information,
                      new.DN.information)

              }

            }else{

              # retrieve observations
              assign.obs = x[assign.locs]

              # check if the current data nugget can be split
              if (length(assign.obs) >= 2*min.nugget.size){

                # check if there are more than 2 observations assigned to this nugget
                if (length(assign.obs) > 2){

                  # if so, split the nugget into two nuggets with k means
                  new.nugget.info = kmeans(x = assign.obs,
                                           centers = 2,
                                           iter.max = 1000,
                                           nstart = 25)

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = shift + 2*(which(large.shapes == i)-1) + new.nugget.info$cluster

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = shift + 2*(which(large.shapes == i)-1) + 1:2

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2] = new.nugget.info$centers[1]
                  new.DN.information[2, 2] = new.nugget.info$centers[2]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = c(sum(new.nugget.info$cluster == 1),
                                                        sum(new.nugget.info$cluster == 2))

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # retrieve the shape for the new data nuggets
                  if (length(new.assign.locs1) > 1){

                    if (is.null(ncol(x))){

                      new.DN.information[1, "Scale"] = var(x[new.assign.locs1])

                    }else{

                      new.DN.information[1, "Scale"] = sum(diag(cov(x[new.assign.locs1, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[1, "Scale"] = 0

                  }

                  if (length(new.assign.locs2) > 1){

                    if (is.null(ncol(x))){

                      new.DN.information[2, "Scale"] = var(x[new.assign.locs2])

                    }else{

                      new.DN.information[2, "Scale"] = sum(diag(cov(x[new.assign.locs2, ])))/ncol(x)

                    }

                  }else{

                    new.DN.information[2, "Scale"] = 0

                  }

                  # check if either of the new nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }else{

                  # reassign the observations to their new nuggets
                  tmp.DN.assignments[assign.locs] = max(DN.information[, "Data Nugget"]) + 1:2

                  # initialize data frame for new data nuggets
                  new.DN.information = DN.information[1:2, ]

                  # create new data nugget names
                  new.DN.information[1:2, "Data Nugget"] = c(max(DN.information[, "Data Nugget"]) + 1:2)

                  # retrieve the new data nugget centers
                  new.DN.information[1, 2] = assign.obs[1]
                  new.DN.information[2, 2] = assign.obs[2]

                  # retrieve the weights for these data nuggets
                  new.DN.information[1:2, "Weight"] = 1

                  # retrieve the assignments for the new data nuggets
                  new.assign.locs1 = which(tmp.DN.assignments == min(tmp.DN.assignments[assign.locs]))
                  new.assign.locs2 = which(tmp.DN.assignments == max(tmp.DN.assignments[assign.locs]))

                  # check if the new data nuggets meet the minimum nugget size criteria
                  if (length(new.assign.locs1) >= min.nugget.size &
                      length(new.assign.locs2) >= min.nugget.size){

                    new.DN.information[, "Delete"] =
                      DN.information[DN.information[, "Data Nugget"] == i, "Data Nugget"]

                    new.DN.information[1, "New Assignment Locations"] = paste(new.assign.locs1,
                                                                              collapse = ",")

                    new.DN.information[2, "New Assignment Locations"] = paste(new.assign.locs2,
                                                                              collapse = ",")

                  }else{

                    new.DN.information = NULL

                  }

                }

                rbind(tmp.new.DN.information,
                      new.DN.information)

              }

            }

          }

        }

      }

      # check if any data nuggets were split
      if (!is.null(tmp.new.DN.information)){

        # bind the new data nuggets with the old data nuggets
        DN.information = rbind.data.frame(DN.information,
                                          tmp.new.DN.information[, -c(ncol(tmp.new.DN.information)-1,
                                                                      ncol(tmp.new.DN.information))])

        # delete the data nuggets that were split
        DN.information = DN.information[!(DN.information[, "Data Nugget"] %in% unique(tmp.new.DN.information[, "Delete"])), ]

        # cycle through the new data nuggets
        for (j in 1:nrow(tmp.new.DN.information)){

          # reassign the data nuggets
          DN.assignments[as.numeric(unlist(strsplit(tmp.new.DN.information[j, "New Assignment Locations"],
                                                    split = ",")))] = tmp.new.DN.information[j, "Data Nugget"]

        }

      }

      # retrieve the data nuggets with shapes larger than the shape tolerance
      large.shapes = DN.information[which(DN.information[, "Scale"] >
                                            orig.threshold), "Data Nugget"]

      # retrieve the amount of data nuggets with large shapes
      new.LS.num = length(large.shapes)

    }

    message("complete!")

    # rename the rows of DN information
    rownames(DN.information) = 1:nrow(DN.information)

    # rename the data nuggets
    DN.information[, "Data Nugget"] = 1:nrow(DN.information)

    # retrieve list of data nuggets
    DN.list = unique(DN.assignments)[order(unique(DN.assignments))]

    # cycle through the list of data nuggets
    for (j in DN.list){

      # change the data nugget number
      DN.assignments[DN.assignments == j] = which(DN.list[order(DN.list)] == j)

    }

  }

  # check if user wants to use parallel processing
  if (no.cores != 0){

    # stop the cluster
    stopCluster(cl)

  }

  # Create output ####
  output = list(DN.information,
                DN.assignments)

  # assign names to the output
  names(output) = c("Data Nuggets",
                    "Data Nugget Assignments")

  # assign the data nugget class to the output
  class(output) = "datanugget"

  # return(output)
  refineeigen.DN=function(xx,DN,delta=2,minsize=10){

    DN1mm=DN
    nc = ncol(xx)
    zz=as.matrix(DN[[1]][,2:(nc+1)])
    nn= as.vector(DN[[1]][,nc+2])
    sn= as.vector(DN[[1]][,nc+3])
    v = as.vector(DN[[2]])
    sthreshold= median(sn)+2*mad(sn)
    n = nrow(zz)

    ii  = (1:n)[nn>=minsize]
    rr= cbind(0*nn,0)
    for(i in ii ) {
      zz1=xx[v==i,,drop=F]
      rr=eigen(var(zz1))
      if(nrow(zz1)>2)
        if(rr$val[1]/rr$val[2] > delta | rr$val[1]> sthreshold){
          iii=kmeans(c(as.matrix(zz1)%*%rr$vectors[,1,drop=F]),2)$clus
          if (min(table(iii)>1)) {
            rr=eigen(var(zz1[iii==1,]))$values[1:2]
            rr1=eigen(var(zz1[iii==2,]))$values[1:2]
            if( 0.5*(rr[1]/rr[2]+rr1[1]/rr1[2])<delta) {
              # DN1mm[[1]][i,2:(nc+3)] = c(apply(zz1[iii==1,],2,mean),sum(iii==1), sqrt(rr[1]))
              DN1mm[[1]][i,2:(nc+3)] = c(apply(zz1[iii==1,],2,mean),sum(iii==1), sum(diag(cov(zz1[iii==1,])))/nc)
              n= nrow(DN1mm[[1]])
              # DN1mm[[1]] =rbind(DN1mm[[1]], c(n+1,apply(zz1[iii==2,],2,mean),sum(iii==2),sqrt(rr1[1])))
              DN1mm[[1]] =rbind(DN1mm[[1]], c(n+1,apply(zz1[iii==2,],2,mean),sum(iii==2),sum(diag(cov(zz1[iii==2,])))/nc))
              DN1mm[[2]][v==i][iii==2] <- n+1
              # print(i)
            }}}}
    DN1mm[[1]][,nc+2] = as.integer(DN1mm[[1]][,nc+2])
    DN1mm
  }

  # return the data nugget dataset
  return(refineeigen.DN(x, output,delta = delta))

}


