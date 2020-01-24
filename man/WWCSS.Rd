\name{WWCSS}
\alias{WWCSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Within Cluster Sum of Squares
}
\description{
This function computes the weighted within cluster sum of squares (WWCSS) for a set of cluster assignments provided to a dataset with some weighting scheme.
}
\usage{
WWCSS(x,
      k,
      P)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix (data frame, data table, matrix, etc) containing only entries of class numeric (i.e. matrix of data nugget centers). Must contain the location, weight, and cluster assignment for each observation.
}
  \item{k}{
The number of possible clusters. Must be of class numeric.
}
  \item{P}{
The number of columns from the original dataset before clustering and without weight variable. Must be of class numeric.
}
}
\details{
The WWCSS is used for optimizing the cluster assignments in the WKmeans function. NOTE THAT THIS FUNCTION IS NOT DESIGNED FOR USE OUTSIDE OF THE WKmeans FUNCTION.
}
\value{

\item{output}{Vector of individual WWCSS's for each cluster}

\item{sum.output}{Combined sum of all individual WWCSS's.}
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
%% \examples{
%% ##---- Should be DIRECTLY executable !! ----
%% ##-- ==>  Define data, use random,
%% ##--	or do  help(data=index)  for the standard data sets.

%% ## The function is currently defined as
%% function (x)
%% {
%% }
%% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
