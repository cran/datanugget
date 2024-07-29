\name{refine.DN}
\alias{refine.DN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Refine Data Nuggets
}
\description{
This function refines the data nuggets found in an object of class datanugget created using the create.DN function.
}
\usage{
refine.DN(x,
          DN,
          EV.tol = .9,
          max.splits = 5,
          min.nugget.size = 2,
          delta = 2,
          seed = 291102,
          no.cores = (detectCores() - 1),
          make.pbs = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{
A data matrix (data frame, data table, matrix, etc.) containing only entries of class numeric.
}
  \item{DN}{
An object of class data nugget created using the create.DN function.
}

  \item{EV.tol}{
A value designating the percentile for finding the corresponding quantile that will designate how large the largest eigenvalue of the covariance matrix of a data nugget can be before it must be split. Must be of class numeric and within (0,1).
}

  \item{max.splits}{
A value designating the maximum amount of attempts that will be made to split data nuggets according to their largest eigenvalue before the algorithm breaks. Must be of class numeric.
}


  \item{min.nugget.size}{
A value designating the minimum amount of observations a data nugget created from a split must contain. Must be of class numeric and greater than 1.
}

 \item{delta}{
Ratio between the first and second eigenvalues of the covariance matrix of a data nugget to force its split. Default is 2.
}

  \item{seed}{
Random seed for replication. Must be of class numeric.
}

  \item{no.cores}{
Number of cores used for parallel processing. If '0' then parallel processing is not used. Must be of class numeric.
}

  \item{make.pbs}{
Print progress bars? Must be TRUE or FALSE.
}


}
\details{
Data nuggets can be refined by attempting to make all of the data nugget shapes as spherical as possible. This is achieved by designating an eigenvalue tolerance (EV.tol) which is used to give a lower threshold for a data nugget's deviation from sphericity, respectively.

If the largest eigenvalue of a data nugget's covariance matrix has a ratio greater than the quantile associated with the percentile given by EV.tol, this data nugget is split into two smaller data nuggets using K-means clustering.

However, if either of the two data nuggets created by this split have less than the designated minimum data nugget size (min.nugget.size), then the split is cancelled and the data nugget remains as is. This function refines data nuggets using Algorithm 2 provided in the reference.

Updated: When data nuggets are not spherical, with the ratio between the first and second eigenvalues of the covariance matrix of the data nugget is greater than \code{delta} (its default value is 2), the data nugget is split.
}
\value{
An object of class datanugget:

\item{Data Nuggets}{DN.num by (ncol(x)+3) data frame containing the information for the data nuggets created (index, center, weight, scale).}

\item{Data Nugget Assignments}{Vector of length nrow(x) containing the data nugget assignment of each observation in x.}

}
\references{
Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, 1-21.

Cherasia, K. E., Cabrera, J., Fernholz, L. T., & Fernholz, R. (2022). Data Nuggets in Supervised Learning. \emph{In Robust and Multivariate Statistical Methods: Festschrift in Honor of David E. Tyler} (pp. 429-449). Cham: Springer International Publishing.
}
\author{
Traymon Beavers, Javier Cabrera, Mariusz Lubomirski
}
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }

\examples{

    ## small example
    X = cbind.data.frame(rnorm(10^3),
                         rnorm(10^3),
                         rnorm(10^3))

    suppressMessages({

      my.DN = create.DN(x = X,
                        R = 500,
                        delete.percent = .1,
                        DN.num1 = 500,
                        DN.num2 = 250,
                        no.cores = 0,
                        make.pbs = FALSE)

      my.DN2 = refine.DN(x = X,
                         DN = my.DN,
                         EV.tol = .9,
                         min.nugget.size = 2,
                         max.splits = 5,
                         no.cores = 0,
                         make.pbs = FALSE)

    })

    my.DN2$`Data Nuggets`
    my.DN2$`Data Nugget Assignments`

    \donttest{

      ## large example
      X = cbind.data.frame(rnorm(5*10^4),
                           rnorm(5*10^4),
                           rnorm(5*10^4),
                           rnorm(5*10^4),
                           rnorm(5*10^4))

      my.DN = create.DN(x = X,
                        R = 5000,
                        delete.percent = .9,
                        DN.num1 = 10^4,
                        DN.num2 = 2000,
                        no.cores = 2)

      my.DN2 = refine.DN(x = X,
                         DN = my.DN,
                         EV.tol = .9,
                         min.nugget.size = 2,
                         max.splits = 5,
                         no.cores = 2)

      my.DN2$`Data Nuggets`
      my.DN2$`Data Nugget Assignments`

    }

}

%% ## The function is currently defined as
%% function (x)
%% {
%% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
