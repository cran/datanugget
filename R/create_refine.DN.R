# This function creates and refines a dataset using create.DN() and refine.DN() functions
# It returns the final data nugget object

# Function inputs
# x: original dataset.
# center.method: this method used for choosing data nugget centers.
# R: the number of observations to sample from the data matrix when creating the initial data nugget centers.
# delete.percent: proportion of data points to be deleted at each iteration.
# DN.num1: the number of initial data nugget centers to create.
# DN.num2: the number of final data nuggets to create.
# dist.metric: pairwise distance measure (e.g., "euclidean" or "manhattan").
# eps: a very small number taken to be the scale for datanuggets with only one observation.
# seed: random seed for replication.
# no.cores: number of cores used for parallel processing.
# make.pbs: logical; whether to show a progress bar while the function runs.
# EV.tol: a value designating the percentile for finding the corresponding quantile that will designate how large the largest eigenvalue of the covariance matrix of a data nugget can be before it must be split. 
# max.splits:	a value designating the maximum amount of attempts that will be made to split data nuggets according to their largest eigenvalue before the algorithm breaks
# min.nugget.size: a value designating the minimum amount of observations a data nugget created from a split must contain. 
# delta: ratio between the first and second eigenvalues of the covariance matrix of a data nugget to force its split.


create_refine.DN = function(x,
                     center.method = "original",
                     R = 5000,
                     delete.percent = .1,
                     DN.num1 = 10^4,
                     DN.num2 = 2000,
                     dist.metric = "euclidean",
                     seed = 291102,
                     no.cores = (detectCores() - 1),
                     make.pbs = TRUE,
                     EV.tol = .9,
                     max.splits = 5,
                     min.nugget.size = 2,
                     delta = 2){

  my.DN = create.DN(x = x,
                    center.method = center.method,
                    R = R,
                    delete.percent = delete.percent,
                    DN.num1 = DN.num1,
                    DN.num2 = DN.num2,
                    dist.metric = dist.metric,
                    seed = seed,
                    no.cores = no.cores,
                    make.pbs = make.pbs)

  my.DN2 = refine.DN(x = x,
                     DN = my.DN,
                     EV.tol = EV.tol,
                     max.splits = max.splits,
                     min.nugget.size = min.nugget.size,
                     delta = delta,
                     seed = seed,
                     no.cores = no.cores,
                     make.pbs = make.pbs)

  return(my.DN2)
}
