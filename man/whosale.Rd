\name{whosale}
\title{Wholesale Dataset}
\alias{whosale}
\description{
  The 'whosale' dataset contains information about wholesale customers' annual spending on various product categories.
}
\usage{
  data(whosale)
}
\format{
  A data frame with 440 rows and 8 columns.
}
\details{
  This dataset includes the following columns:
  \describe{
    \item{y}{Type of customer, either 'Horeca' (Hotel/Restaurant/Cafe), coded as \code{1} or 'Retail' coded as \code{2}.}
    \item{Region}{Geographic region of the customer, either 'Lisbon', 'Oporto', or 'Other'. Coded as \code{{1,2,3}}, respectively.}
    \item{Fresh}{Annual spending (in monetary units) on fresh products.}
    \item{Milk}{Annual spending on milk products.}
    \item{Grocery}{Annual spending on grocery products.}
    \item{Frozen}{Annual spending on frozen products.}
    \item{Detergents Paper}{Annual spending on detergents and paper products.}
    \item{Delicassen}{Annual spending on delicatessen products.}
  }
}
\source{
  The 'whosale' dataset is sourced from the UCI Machine Learning Repository:
  \url{https://archive.ics.uci.edu/ml/datasets/wholesale+customers}
}
\examples{
    data(whosale)
    head(whosale)
}
\keyword{datasets}
