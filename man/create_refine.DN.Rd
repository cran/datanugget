\name{create_refine.DN}
\alias{create_refine.DN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create and Refine Data Nuggets in one function
}
\description{
This function combines creating and refining data nuggets in one function. It's a wrapper function for \code{create.DN} and \code{refine.DN}.
}
\usage{
create_refine.DN(x,
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
                 delta = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{
A data matrix (of class matrix, data.frame, or data.table) containing only entries of class numeric.
}

  \item{center.method}{
The method used for choosing data nugget centers. Must be 'mean' or 'random'. 'mean' chooses the data nugget center to be the mean of all observations within that data nugget, while 'random' chooses the data nugget center to be some random observation within that data nugget.
}

  \item{R}{
The number of observations to sample from the data matrix when creating the initial data nugget centers. Must be of class numeric within [100,10000].
}

  \item{delete.percent}{
The proportion of observations to remove from the data matrix at each iteration when finding data nugget centers. Must be of class numeric and within (0,1).
}

  \item{DN.num1}{
The number of initial data nugget centers to create. Must be of class numeric.
}

  \item{DN.num2}{
The number of data nuggets to create. Must be of class numeric.
}

  \item{dist.metric}{
The distance metric used to create the initial centers of data nuggets. Must be 'euclidean' or 'manhattan'.
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

}
\details{
Data nuggets are a representative sample meant to summarize Big Data by reducing a large dataset to a much smaller dataset by eliminating redundant points while also preserving the peripheries of the dataset. Each data nugget is defined by a center (location), weight (importance), and scale (internal variability). This function combines creating and refining data nuggets in one function. It's a wrapper function for \code{create.DN} and \code{refine.DN}.
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

        my.DN = create_refine.DN(x = X,
                          R = 500,
                          delete.percent = .1,
                          DN.num1 = 500,
                          DN.num2 = 250,
                          no.cores = 0,
                          make.pbs = FALSE,
                          EV.tol = .9,
                          min.nugget.size = 2,
                          max.splits = 5,
                          delta = 2)

      })

      my.DN$`Data Nuggets`
      my.DN$`Data Nugget Assignments`

    \donttest{

      ## large example
      X = cbind.data.frame(rnorm(5*10^4),
                           rnorm(5*10^4),
                           rnorm(5*10^4),
                           rnorm(5*10^4),
                           rnorm(5*10^4))

      my.DN = create_refine.DN(x = X,
                          R = 5000,
                          delete.percent = .9,
                          DN.num1 = 10^4,
                          DN.num2 = 2000,
                          no.cores = 2,
                          EV.tol = .9,
                          min.nugget.size = 2,
                          max.splits = 5,
                          delta = 2)

      my.DN$`Data Nuggets`
      my.DN$`Data Nugget Assignments`

    }

}
