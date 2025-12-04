\name{getDN.means}
\alias{getDN.means}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the means of a given Data Nugget
}
\description{
This function calculates the datanugget means of a given datanugget object.
}
\usage{
getDN.means(x, 
            DN)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix (of class matrix, data.frame, or data.table) containing only entries of class numeric.
}


  \item{DN}{
A datanugget object created using the create.DN or refined using the refine.DN function. 
  }
}
\details{
create.DN function generates Data nuggets with centers based on different centering methods viz., mean, random and original. For any datanugget object created by create.DN or refined by refine.DN, this function computes the means of the datanuggets.
}
\value{

\item{out}{nrow(DN$'Data Nuggets') by ncol(x) data frame containing the means of the data nuggets.}
}
\references{
Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, 1-21.

Cherasia, K. E., Cabrera, J., Fernholz, L. T., & Fernholz, R. (2022). Data Nuggets in Supervised Learning. \emph{In Robust and Multivariate Statistical Methods: Festschrift in Honor of David E. Tyler} (pp. 429-449). Cham: Springer International Publishing.
}
\author{
Rituparna Dey, Traymon Beavers, Javier Cabrera, Mariusz Lubomirski
}

\examples{
      
      set.seed(55)
      X = cbind.data.frame(rnorm(10^3),
                           rnorm(10^3),
                           rnorm(10^3))
      
      
      # with center.method = "original"
      my.DN = create.DN(x = X,
                        R = 500,
                        delete.percent = .1,
                        DN.num1 = 500,
                        DN.num2 = 250,
                        no.cores = 0,
                        make.pbs = FALSE)

      

      my.DN.means = getDN.means(x = X, 
                                DN = my.DN)
      my.DN.means
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
