\name{WKmeans}
\alias{WKmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted K-means Clustering of Data Nuggets
}
\description{
This function clusters data nuggets using a form of weighted K-means clustering.
}
\usage{
WKmeans(dataset,
        k,
        cl.centers = NULL,
        obs.weights,
        num.init = 1,
        max.iterations = 10,
        print.progress = TRUE,
        seed = 291102,
        reassign.prop = .25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{dataset}{
A data matrix (data frame, data table, matrix, etc) containing only entries of class numeric (i.e. matrix of data nugget centers).
}

  \item{k}{
Number of desired clusters. Must be of class numeric.
}

  \item{cl.centers}{
Chosen cluster centers. If not NULL, must be a k by ncol(dataset) matrix containing only entries of class numeric.
}

  \item{obs.weights}{
Vector of length nrow(dataset) of weights for each observation in the dataset. Must be of class numeric.
}

  \item{num.init}{
Number of initial clusters to attempt. Ignored if cl.centers is not NULL. Must be of class numeric.
}

  \item{max.iterations}{
Maximum number of iterations attempted for convergence before quitting. Must be of class numeric.
}

  \item{print.progress}{
Print progress of algorithm? Must be TRUE or FALSE.
}

  \item{seed}{
Random seed for replication. Must be of class numeric.
}

  \item{reassign.prop}{
Proportion of data to attempt to reassign during each iteration. Must be of class numeric and within (0,1].
}

}
\details{
Weighted K-means clustering can be used as an unsupervised learning technique to cluster observations contained in datasets that also have a measure of importance (e.g. weight) associated with them. In the case of data nuggets, this is the weight parameter associated with the data nuggets, so the centers of data nuggets are clustered using their weight parameters. The objective of the algorithm which performs this method of clustering is to minimize the weighted within cluster sum of squares (WWCSS). This function clusters data nuggets using Algorithm 3 provided in the reference.

Note that although this method was designed for use with data nuggets, there is no obvious reason to suggest that it cannot be used to perform clustering for other datasets which also have some weighting scheme.
}
\value{

\item{Cluster Assignments}{Vector of length nrow(dataset) containing the cluster assignment for each observation.}

\item{Cluster Centers}{k by ncol(dataset) matrix containing the cluster centers for each cluster.}

\item{Weighted WCSS}{List containing the individual WWCSS for each cluster and the combined sum of all individual WWCSS's.}

}
\references{
Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure (Submitted for Publication, 2019)
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
    X = cbind.data.frame(rnorm(10^4),
                         rnorm(10^4),
                         rnorm(10^4))

    suppressMessages({

      my.DN = create.DN(x = X,
                        RS.num = 10^3,
                        DN.num1 = 500,
                        DN.num2 = 250,
                        no.cores = 0,
                        make.pbs = FALSE)

      my.DN2 = refine.DN(x = X,
                         DN = my.DN,
                         scale.tol = .9,
                         shape.tol = .9,
                         min.nugget.size = 2,
                         max.nuggets = 1000,
                         scale.max.splits = 5,
                         shape.max.splits = 5,
                         no.cores = 0,
                         make.pbs = FALSE)


      DN.clus = WKmeans(dataset = my.DN2$`Data Nuggets`[, c("Center1",
                                                            "Center2",
                                                            "Center3")],
                        k = 3,
                        obs.weights = my.DN2$`Data Nuggets`[, "Weight"],
                        num.init = 1,
                        max.iterations = 3,
                        reassign.prop = .33,
                        print.progress = FALSE)

    })

    DN.clus$`Cluster Assignments`
    DN.clus$`Cluster Centers`
    DN.clus$`Weighted WCSS`

    \donttest{

      ## large example
      X = cbind.data.frame(rnorm(10^6),
                           rnorm(10^6),
                           rnorm(10^6),
                           rnorm(10^6),
                           rnorm(10^6))

      my.DN = create.DN(x = X,
                        RS.num = 10^5,
                        DN.num1 = 10^4,
                        DN.num2 = 2000)

      my.DN$`Data Nuggets`
      my.DN$`Data Nugget Assignments`

      my.DN2 = refine.DN(x = X,
                         DN = my.DN,
                         scale.tol = .9,
                         shape.tol = .9,
                         min.nugget.size = 2,
                         max.nuggets = 10000,
                         scale.max.splits = 5,
                         shape.max.splits = 5)

      my.DN2$`Data Nuggets`
      my.DN2$`Data Nugget Assignments`

      DN.clus = WKmeans(dataset = my.DN2$`Data Nuggets`[, c("Center1",
                                                            "Center2",
                                                            "Center3")],
                        k = 3,
                        obs.weights = my.DN2$`Data Nuggets`[, "Weight"],
                        num.init = 1,
                        max.iterations = 3,
                        reassign.prop = .33)

      DN.clus$`Cluster Assignments`
      DN.clus$`Cluster Centers`
      DN.clus$`Weighted WCSS`

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
