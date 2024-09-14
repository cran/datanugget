AC = function(x,
              R,
              delete.percent,
              DN.num1,
              DN.num2){

  # check if x is a vector
  if (!(class(x) %in% c("matrix",
                        "data.frame",
                        "data.table"))){

    stop('x must be of class "matrix", "data.frame", or "data.table"')

  }

  # make sure R is a number
  if (!is.numeric(R)){

    stop('R must be of class "numeric"')

  }

  # make sure R is between 100 and 10000
  if (R < 100 |
      R > 10000){

    stop('R must be within [100,10000]')

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

  # make sure DN.num1 is numeric
  if (!is.numeric(DN.num1)){

    stop('DN.num1 must be of class "numeric"')

  }

  # make sure DN.num2 is numeric
  if (!is.numeric(DN.num2)){

    stop('DN.num2 must be of class "numeric"')

  }

  if (DN.num2 > DN.num1){

    stop('DN.num2 must be less than DN.num1')

  }

  R = floor(R)
  DN.num1 = floor(DN.num1)
  DN.num2 = floor(DN.num2)

  N = nrow(x)
  P = ncol(x)
  G = ceiling(N/R)

  A = 1

  while (R != DN.num1/G){

    if (R - floor(delete.percent*R) < DN.num1/G){

      R = R - (R-DN.num1/G)

    }else if (floor(delete.percent*R) == 0){

      R = R - 1

    }else{

      R = R - floor(delete.percent*R)

    }

    A = A+1

  }

  B = 1

  while (DN.num1 != DN.num2){

    if (DN.num1 - floor(delete.percent*DN.num1) < DN.num2){

      DN.num1 = DN.num1 - (DN.num1-DN.num2)

    }else if (floor(delete.percent*DN.num1) == 0){

      DN.num1 = DN.num1 - 1

    }else{

      DN.num1 = DN.num1 - floor(delete.percent*DN.num1)

    }

    B = B+1

  }

  my.AC = log10(P*(G*((R*(R-1)/2)*(A)) +
                     (DN.num1*(DN.num1-1)/2)*(B) +
                     (N+1)*DN.num2))

  return(my.AC)

}
