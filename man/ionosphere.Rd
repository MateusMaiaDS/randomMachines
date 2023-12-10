\name{ionosphere}
\title{Ionosphere Dataset}
\alias{ionosphere}
\description{
  The 'ionosphere' dataset contains radar data for the classification of radar returns as either 'good' or 'bad'.
}
\usage{
  data(ionosphere)
}
\format{
  A data frame with 351 rows and 35 columns.
}
\details{
  This dataset includes the following columns:
  \describe{
    \item{X1-X34}{Features extracted from radar signals.}
    \item{y}{Class label indicating whether the radar return is 'g' (good) or 'b' (bad).}
  }
}
\source{
  The 'ionosphere' dataset is sourced from the UCI Machine Learning Repository:
  \url{https://archive.ics.uci.edu/ml/datasets/ionosphere}
}
\examples{
    data(ionosphere)
    head(ionosphere)
}
\keyword{datasets}
