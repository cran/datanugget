create_refine.DN = function(x,
                     center.method = "mean",
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
                     delta = 2
                     ){

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
