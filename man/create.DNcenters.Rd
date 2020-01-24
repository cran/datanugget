\name{create.DNcenters}
\alias{create.DNcenters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Data Nugget Centers
}
\description{
This function creates the centers of data nuggets from a random sample.
}
\usage{
create.DNcenters(RS,
                 DN.num,
                 dist.metric,
                 make.pb = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{RS}{
A data matrix (data frame, data table, matrix, etc) containing only entries of class numeric.
}

  \item{DN.num}{
The number of data nuggets to create. Must be of class numeric.
}

  \item{dist.metric}{
The distance metric used to create the initial centers of data nuggets. Must be 'euclidean' or 'manhattan'.
}

  \item{make.pb}{
Print progress bar? Must be TRUE or FALSE.
}


}
\details{
This function is used for reducing a random sample to data nugget centers in the create.DN function. NOTE THAT THIS FUNCTION IS NOT DESIGNED FOR USE OUTSIDE OF THE create.DN FUNCTION.
}
\value{

\item{DN.data}{DN.num by (ncol(RS)) data frame containing the data nugget centers.}

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
