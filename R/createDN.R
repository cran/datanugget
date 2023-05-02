create.DN = function(x,
                     center.method = "mean",
                     R = 5000,
                     delete.percent = .1,
                     DN.num1 = 10^4,
                     DN.num2 = 2000,
                     dist.metric = "euclidean",
                     seed = 291102,
                     no.cores = (detectCores() - 1),
                     make.pbs = TRUE){

  # Argument checking/fixing ####

  # check if x is a vector
  if (!(class(x) %in% c("matrix",
                        "data.frame",
                        "data.table"))){

    stop('x must be of class "matrix", "data.frame", or "data.table"')

  }

  # make sure R is a number
  if (!(class(R) %in% c("numeric",
                        "integer"))){

    stop('R must be of class "numeric" or "integer"')

  }

  # make sure R is between 100 and 10000
  if (R < 100 |
      R > 10000){

    stop('R must be within [100,10000]')

  }

  # make sure center.method is mean or random
  if (!(center.method %in% c("mean",
                           "random"))){

    stop('center.method must be "mean" or "random"')

  }

  # make sure delete.percent is a number
  if (!is.numeric(delete.percent)){

    stop('delete.percent must be of class "numeric"')

  }

  # make sure delete.percent is between 0 and 1
  if (delete.percent <= 0 |
      delete.percent >= 1){

    stop('delete.percent must be within (0,1)')

  }

  # make sure DN.num1 is a number
  if (!(class(DN.num1) %in% c("numeric",
                            "integer"))){

    stop('DN.num1 must be of class "numeric" or "integer"')

  }

  # make sure DN.num2 is a number
  if (!(class(DN.num2) %in% c("numeric",
                            "integer"))){

    stop('DN.num2 must be of class "numeric" or "integer"')

  }

  if (DN.num2 > DN.num1){

    stop('DN.num2 must be less than DN.num1')

  }

  # convert RS.num and DN.num to integers
  RS.num = floor(nrow(x))
  R = floor(R)
  DN.num1 = floor(DN.num1)
  DN.num2 = floor(DN.num2)

  # make sure given distance metric can be used
  if (dist.metric != "euclidean" &
      dist.metric != "manhattan" &
      dist.metric != "minkowski"){

    stop('dist.metric must be "euclidean" or "manhattan"')


  }

  # # make sure minkowksi.p is numeric
  # if (class(minkowski.p) != "numeric"){
  #
  #   stop('minkowski.p must be of class "numeric"')
  #
  # }

  # make sure seed is numeric
  if (!is.numeric(seed)){

    stop('seed must be of class "numeric"')

  }

  # make sure no.cores is numeric or integer
  if (!(class(no.cores) %in% c("numeric",
                              "integer"))){

    stop('no.cores must be of class "numeric" or "integer"')

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

  # retrieve the number of observations
  obs.num = nrow(x)

  # create set.num
  set.num = R

  # count the number of sets of set.num the random sample will be split into
  RS.splits = RS.num%/%set.num + (RS.num%%set.num > 0)

  # find the locations of the observations that will compose the random sample
  RS.loc = sample(obs.num,
                  RS.num,
                  replace = FALSE)

  # create initial set of data nugget centers
  message("creating initial set of data nugget centers...")

  # check if user wants to use parallel processing
  if (no.cores != 0){

    # check if user wants a progress bar
    if (make.pbs == TRUE){

      # initialize progress bar
      pb = txtProgressBar(min = 0,
                          max = RS.splits)

      # update the progress bar
      progress = function(n){setTxtProgressBar(pb, n)}
      opts = list(progress = progress)

      DN.centers1 = foreach(i = 1:RS.splits,
                            .combine = rbind,
                            .export = "create.DNcenters",
                            .options.snow = opts) %dopar%

        {

          # set the beginning and ends of the random sample selection
          RS.loc.start = ((i-1)*set.num)+1
          RS.loc.end = ifelse(set.num*i > length(RS.loc),
                              length(RS.loc),
                              set.num*i)

          tmp.RS = x[RS.loc[RS.loc.start:RS.loc.end], ]

          # retrieve the data nugget centers
          output = create.DNcenters(RS = tmp.RS,
                                    delete.percent = delete.percent,
                                    DN.num = min(c(ceiling(DN.num1/RS.splits),
                                                   length(RS.loc.start:RS.loc.end))),
                                    dist.metric = dist.metric,
                                    make.pb = FALSE)

          output = as.data.frame(output)

          # return the data nugget centers
          return(output)

        }

      # close the progress bar
      close(pb)

    }else{

      DN.centers1 = foreach(i = 1:RS.splits,
                            .combine = rbind) %dopar%

        {

          # set the beginning and ends of the random sample selection
          RS.loc.start = ((i-1)*set.num)+1
          RS.loc.end = ifelse(set.num*i > length(RS.loc),
                              length(RS.loc),
                              set.num*i)

          tmp.RS = x[RS.loc[RS.loc.start:RS.loc.end], ]

          # retrieve the data nugget centers
          output = create.DNcenters(RS = tmp.RS,
                                    delete.percent = delete.percent,
                                    DN.num = min(c(ceiling(DN.num1/RS.splits),
                                                   length(RS.loc.start:RS.loc.end))),
                                    dist.metric = dist.metric,
                                    make.pb = FALSE)

          output = as.data.frame(output)

          # return the data nugget centers
          return(output)

        }

    }

  }else{

    # check if user wants a progress bar
    if (make.pbs == TRUE){

      # initialize progress bar
      pb = txtProgressBar(min = 0,
                          max = RS.splits)

      # initialize the initial centers of the data nuggets
      DN.centers1 = NULL

      # cycle through the random sample splits
      for (i in 1:RS.splits){

        # set the beginning and ends of the random sample selection
        RS.loc.start = ((i-1)*set.num)+1
        RS.loc.end = ifelse(set.num*i > length(RS.loc),
                            length(RS.loc),
                            set.num*i)

        tmp.RS = x[RS.loc[RS.loc.start:RS.loc.end], ]

        # retrieve the data nugget centers
        output = create.DNcenters(RS = tmp.RS,
                                  delete.percent = delete.percent,
                                  DN.num = min(c(ceiling(DN.num1/RS.splits),
                                                 length(RS.loc.start:RS.loc.end))),
                                  dist.metric = dist.metric,
                                  make.pb = FALSE)

        DN.centers1 = rbind.data.frame(DN.centers1,
                                       output)

        setTxtProgressBar(pb,i)

      }

      # close the progress bar
      close(pb)

    }else{

      # initialize the initial centers of the data nuggets
      DN.centers1 = NULL

      # cycle through the random sample splits
      for (i in 1:RS.splits){

        # set the beginning and ends of the random sample selection
        RS.loc.start = ((i-1)*set.num)+1
        RS.loc.end = ifelse(set.num*i > length(RS.loc),
                            length(RS.loc),
                            set.num*i)

        tmp.RS = x[RS.loc[RS.loc.start:RS.loc.end], ]

        # retrieve the data nugget centers
        output = create.DNcenters(RS = tmp.RS,
                                  delete.percent = delete.percent,
                                  DN.num = min(c(ceiling(DN.num1/RS.splits),
                                                 length(RS.loc.start:RS.loc.end))),
                                  dist.metric = dist.metric,
                                  make.pb = FALSE)

        DN.centers1 = rbind.data.frame(DN.centers1,
                                       output)

      }

    }

  }

  message("completed!")

  # create final set of data nugget centers
  message("creating final set of data nugget centers...")

  DN.data = create.DNcenters(RS = DN.centers1,
                             delete.percent = delete.percent,
                             DN.num = DN.num2,
                             dist.metric = dist.metric,
                             make.pb = make.pbs)

  DN.data = as.data.frame(DN.data)

  message("completed!")

  # assign the observations to data nuggets
  message("assigning observations to data nuggets...")

  if (no.cores != 0){

    parallel::clusterExport(cl = cl,
                  varlist = c("dist.metric",
                              "DN.data"),
                  envir = environment())

    if (ncol(x) == 1){

      DN.assignments = parallel::parApply(cl = cl,
                                X = x,
                                MARGIN = 1,
                                FUN = function(input){

                                  # apply the given distance metric for assignment
                                  if (dist.metric == "euclidean"){

                                    return(as.numeric(which.min((t(DN.data) - input)^2)))

                                  }else if (dist.metric == "manhattan"){

                                    return(as.numeric(which.min(abs((t(DN.data) - input)))))

                                    # }else if (dist.metric == "minkowski"){
                                    #
                                    #   return(as.numeric(which.min((colSums(abs((t(DN.data) - input))^minkowski.p))^(1/minkowski.p))))

                                  }

                                })

    }else{

      DN.assignments = parallel::parApply(cl = cl,
                                X = x,
                                MARGIN = 1,
                                FUN = function(input){

                                  # apply the given distance metric for assignment
                                  if (dist.metric == "euclidean"){

                                    return(as.numeric(which.min(colSums((t(DN.data) - input)^2))))

                                  }else if (dist.metric == "manhattan"){

                                    return(as.numeric(which.min(colSums(abs((t(DN.data) - input))))))

                                    # }else if (dist.metric == "minkowski"){
                                    #
                                    #   return(as.numeric(which.min((colSums(abs((t(DN.data) - input))^minkowski.p))^(1/minkowski.p))))

                                  }

                                })

    }

  }else{

    if (ncol(x) == 1){

      DN.assignments = apply(X = x,
                             MARGIN = 1,
                             FUN = function(input){

                               # apply the given distance metric for assignment
                               if (dist.metric == "euclidean"){

                                 return(as.numeric(which.min((t(DN.data) - input)^2)))

                               }else if (dist.metric == "manhattan"){

                                 return(as.numeric(which.min(abs((t(DN.data) - input)))))

                                 # }else if (dist.metric == "minkowski"){
                                 #
                                 #   return(as.numeric(which.min((colSums(abs((t(DN.data) - input))^minkowski.p))^(1/minkowski.p))))

                               }

                             })

    }else{

      DN.assignments = apply(X = x,
                             MARGIN = 1,
                             FUN = function(input){

                               # apply the given distance metric for assignment
                               if (dist.metric == "euclidean"){

                                 return(as.numeric(which.min(colSums((t(DN.data) - input)^2))))

                               }else if (dist.metric == "manhattan"){

                                 return(as.numeric(which.min(colSums(abs((t(DN.data) - input))))))

                                 # }else if (dist.metric == "minkowski"){
                                 #
                                 #   return(as.numeric(which.min((colSums(abs((t(DN.data) - input))^minkowski.p))^(1/minkowski.p))))

                               }

                             })

    }

  }

  message("completed!")

  # initialize the dataset holding th edata nugget information
  DN.information = data.frame(Data.Nugget = 1:DN.num2)

  # bind the data nugget centers the information dataset
  DN.information = cbind.data.frame(DN.information, DN.data)

  # give column names for the centers
  if (ncol(DN.data) != 1){

    colnames(DN.information)[1:(ncol(DN.data)+1)] = c("Data Nugget",
                                                      paste("Center",
                                                            1:ncol(DN.data),
                                                            sep = ""))

  }else{

    colnames(DN.information)[1:2] = c("Data Nugget",
                                      "Center")

  }

  # find the weight for each data nugget
  DN.information[, "Weight"] = table(DN.assignments)

  message("recalculating data nugget centers and calculating data nugget scales...")

  # check if user wants to use parallel processing
  if (no.cores != 0){

    # check if user wants a progress bar
    if (make.pbs == TRUE){

      # initialize progress bar
      pb = txtProgressBar(min = 0,
                          max = DN.num2)

      # update the progress bar
      progress = function(n){setTxtProgressBar(pb, n)}
      opts = list(progress = progress)

      # recalculate the center and calculate the shape parameter for each data nugget
      for.params = foreach(i = 1:DN.num2,
                           .combine = rbind,
                           .export = "create.DNcenters",
                           .options.snow = opts)  %dopar%

        {

          # extract the row locations of the data points assigned to the current data nugget
          tmp.extract = which(DN.assignments == i)

          if (ncol(x) == 1){

            # check if there is only one data point assigned to this data nugget
            if (length(tmp.extract) == 1){

              # if so, make the center equal to the lone observation
              output1 = x[tmp.extract, ]

              # if so, make the shape parameter missing
              output2 = 0

            }else{

              # otherwise, find the new center
              if (center.method == "mean"){

                output1 = mean(x[tmp.extract, ])

              }else if (center.method == "random"){

                output1 = x[sample(tmp.extract,1), ]

              }

              # otherwise, find the scale parameter
              output2 = var(x[tmp.extract, ])

            }

          }else{

            # check if there is only one data point assigned to this data nugget
            if (length(tmp.extract) == 1){

              # if so, make the center equal to the lone observation
              output1 = x[tmp.extract, ]

              # if so, make the shape parameter missing
              output2 = 0

            }else{

              # otherwise, find the new center
              if (center.method == "mean"){

                output1 = colMeans(x[tmp.extract, ])

              }else if (center.method == "random"){

                output1 = x[sample(tmp.extract,1), ]

              }

              # otherwise, find the scale parameter
              output2 = sum(diag(cov(x[tmp.extract, ])))/ncol(x)

            }

          }


          # return the center and shape parameters
          return(c(as.vector(output1),
                   output2))

        }

      # close the progress bar
      close(pb)

    }else{

      # recalculate the center and calculate the shape parameter for each data nugget
      for.params = foreach(i = 1:DN.num2,
                           .combine = rbind)  %dopar%

        {

          # extract the row locations of the data points assigned to the current data nugget
          tmp.extract = which(DN.assignments == i)

          if (ncol(x) == 1){

            # check if there is only one data point assigned to this data nugget
            if (length(tmp.extract) == 1){

              # if so, make the center equal to the lone observation
              output1 = x[tmp.extract, ]

              # if so, make the shape parameter missing
              output2 = 0

            }else{

              # otherwise, find the new center
              if (center.method == "mean"){

                output1 = mean(x[tmp.extract, ])

              }else if (center.method == "random"){

                output1 = x[sample(tmp.extract,1), ]

              }

              # otherwise, find the scale parameter
              output2 = var(x[tmp.extract, ])

            }

          }else{

            # check if there is only one data point assigned to this data nugget
            if (length(tmp.extract) == 1){

              # if so, make the center equal to the lone observation
              output1 = x[tmp.extract, ]

              # if so, make the shape parameter missing
              output2 = 0

            }else{

              # otherwise, find the new center
              if (center.method == "mean"){

                output1 = colMeans(x[tmp.extract, ])

              }else if (center.method == "random"){

                output1 = x[sample(tmp.extract,1), ]

              }

              # otherwise, find the scale parameter
              output2 = sum(diag(cov(x[tmp.extract, ])))/ncol(x)

            }

          }


          # return the center and shape parameters
          return(c(as.vector(output1),
                   output2))

        }

    }

  }else{

    # check if user wants a progress bar
    if (make.pbs == TRUE){

      # initialize progress bar
      pb = txtProgressBar(min = 0,
                          max = DN.num2)

      # initialize matrix holding new centers and scales
      for.params = NULL

      # cycle through the final data nuggets
      for (i in 1:DN.num2){

        # extract the row locations of the data points assigned to the current data nugget
        tmp.extract = which(DN.assignments == i)

        if (ncol(x) == 1){

          # check if there is only one data point assigned to this data nugget
          if (length(tmp.extract) == 1){

            # if so, make the center equal to the lone observation
            output1 = x[tmp.extract, ]

            # if so, make the shape parameter missing
            output2 = 0

          }else{

            # otherwise, find the new center
            output1 = mean(x[tmp.extract, ])

            # otherwise, find the scale parameter
            output2 = var(x[tmp.extract, ])

          }

        }else{

          # check if there is only one data point assigned to this data nugget
          if (length(tmp.extract) == 1){

            # if so, make the center equal to the lone observation
            output1 = x[tmp.extract, ]

            # if so, make the shape parameter missing
            output2 = 0

          }else{

            # otherwise, find the new center
            output1 = colMeans(x[tmp.extract, ])

            # otherwise, find the scale parameter
            output2 = sum(diag(cov(x[tmp.extract, ])))/ncol(x)

          }

        }


        # return the center and shape parameters
        for.params = rbind(for.params,
                           c(as.vector(output1),
                             output2))

        setTxtProgressBar(pb, i)

      }

      # close the progress bar
      close(pb)

    }else{

      # initialize matrix holding new centers and scales
      for.params = NULL

      # cycle through the final data nuggets
      for (i in 1:DN.num2){

        # extract the row locations of the data points assigned to the current data nugget
        tmp.extract = which(DN.assignments == i)

        if (ncol(x) == 1){

          # check if there is only one data point assigned to this data nugget
          if (length(tmp.extract) == 1){

            # if so, make the center equal to the lone observation
            output1 = x[tmp.extract, ]

            # if so, make the shape parameter missing
            output2 = 0

          }else{

            # otherwise, find the new center
            output1 = mean(x[tmp.extract, ])

            # otherwise, find the scale parameter
            output2 = var(x[tmp.extract, ])

          }

        }else{

          # check if there is only one data point assigned to this data nugget
          if (length(tmp.extract) == 1){

            # if so, make the center equal to the lone observation
            output1 = x[tmp.extract, ]

            # if so, make the shape parameter missing
            output2 = 0

          }else{

            # otherwise, find the new center
            output1 = colMeans(x[tmp.extract, ])

            # otherwise, find the scale parameter
            output2 = sum(diag(cov(x[tmp.extract, ])))/ncol(x)

          }

        }


        # return the center and shape parameters
        for.params = rbind(for.params,
                           c(as.vector(output1),
                             output2))

      }

    }

  }

  message("completed!")

  # check if user wants to use parallel processing
  if (no.cores != 0){

    # stop the cluster
    stopCluster(cl)

  }

  # place the center parameter info in the data nugget dataset
  for (i in 2:(ncol(DN.information)-1)){

    DN.information[, i] = unlist(as.data.frame(for.params)[, i-1])

  }

  # place the shape parameter info in the data nugget dataset
  DN.information[, "Scale"] = unlist(as.data.frame(for.params)[, ncol(for.params)])

  # change the rownames
  rownames(DN.information) = 1:DN.num2

  # Create output ####
  output = list(DN.information,
                DN.assignments)

  # assign names to the output
  names(output) = c("Data Nuggets",
                    "Data Nugget Assignments")

  # assign the data nugget class to the output
  class(output) = "datanugget"

  # return the data nugget dataset
  return(output)

}
