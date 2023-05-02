\name{AC}
\alias{AC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Arithmetic Complexicity of the Algorithm That Creates Data Nuggets
}
\description{
This function creates the centers of data nuggets from a random sample.
}
\usage{
AC(x,
   R,
   delete.percent,
   DN.num1,
   DN.num2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{
A data matrix (of class matrix, data.frame, or data.table) containing only entries of class numeric.
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

}
\details{
This function is used for calculating the arithmetic complexicity of the algorithm behind the create.DN function for the given parameter choices.
}
\value{

\item{my.AC}{The arithmetic complexicity of the algorithm behind the create.DN function for the given parameter choices on a log10 scale.}

}
\references{
Cherasia, K. E., Cabrera, J., Fernholz, L. T., & Fernholz, R. (2022). Data Nuggets in Supervised Learning. \emph{In Robust and Multivariate Statistical Methods: Festschrift in Honor of David E. Tyler} (pp. 429-449). Cham: Springer International Publishing.

Beavers, T., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., Teigler, J. (2023). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure (Submitted for Publication)
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

      X = cbind.data.frame(rnorm(10^6),
                           rnorm(10^6),
                           rnorm(10^6),
                           rnorm(10^6),
                           rnorm(10^6))

      my.AC = AC(x = X,
                 R = 5000,
                 delete.percent = .1,
                 DN.num1 = 10^4,
                 DN.num2 = 2000)

}

%% ## The function is currently defined as
%% function (x)
%% {
%% }
%% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
