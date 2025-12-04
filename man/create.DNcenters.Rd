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
                 delete.percent = 0.1,
                 DN.num,
                 dist.metric = "euclidean",
                 make.pbs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{RS}{
A data matrix (data frame, data table, matrix, etc) containing only entries of class numeric.
}

  \item{delete.percent}{
The proportion of observations to remove from the data matrix at each iteration when finding data nugget centers. Must be of class numeric and within (0,1). Default value is 0.1.
}


  \item{DN.num}{
The number of data nuggets to create. Must be of class numeric.
}

  \item{dist.metric}{
The distance metric used to create the initial centers of data nuggets. Must be 'euclidean' or 'manhattan'. Default is "euclidean".
}

  \item{make.pbs}{
Logical; whether to show a progress bar while the function runs. Default is FALSE.
}


}
\details{
This function is used for reducing a random sample to data nugget centers in the create.DN function. NOTE THAT THIS FUNCTION IS NOT DESIGNED FOR USE OUTSIDE OF THE create.DN FUNCTION.
}
\value{

\item{DN.data}{DN.num by ncol(RS) data frame containing the data nugget centers.}

}
\references{
Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, 1-21.


Cherasia, K. E., Cabrera, J., Fernholz, L. T., & Fernholz, R. (2022). Data Nuggets in Supervised Learning. \emph{In Robust and Multivariate Statistical Methods: Festschrift in Honor of David E. Tyler} (pp. 429-449). Cham: Springer International Publishing.
}
\author{
Rituparna Dey, Traymon Beavers, Javier Cabrera, Mariusz Lubomirski
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
