\name{bolsafam}
\title{Bolsa Família Dataset}
\alias{bolsafam}
\description{
  The 'bolsafam' dataset contains information about the utilization rate of the Bolsa Família program in Brazilian municipalities. The utilization rate \eqn{y_{i}} is defined as the number of people benefiting from the assistance divided by the total population of the city.
}
\usage{
  data(bolsafam)
}
\format{
  A data frame with 5564 rows and 11 columns.
}
\details{
  This dataset includes the following columns:
  \describe{
    \item{y}{Rate of use of the social assistance program by municipality.}
    \item{COD_UF}{Code to identify the Brazilian state to which the city belongs.}
    \item{T_DENS}{Percentage of the population living in households with a density greater than 2 people per bedroom.}
    \item{TRABSC}{Percentage of employed persons aged 18 or over who are employed without a formal contract.}
    \item{PPOB}{Proportion of people vulnerable to poverty.}
    \item{T_NESTUDA_NTRAB_MMEIO}{Percentage of people aged 15 to 24 who do not study or work and are vulnerable to poverty.}
    \item{T_FUND15A17}{Percentage of the population aged 15 to 17 with complete primary education.}
    \item{RAZDEP}{Dependency ratio.}
    \item{T_ATRASO_0_BASICO}{Percentage of the population aged 6 to 17 years attending basic education that does not have an age-grade delay.}
    \item{T_AGUA}{Percentage of the population living in households with running water.}
    \item{REGIAO}{Aggregation of states according to the regions defined by IBGE.}
  }
}
\source{
  The 'bolsafam' dataset is sourced from the Brazilian organizational site called \emph{Transparency Portal}.
}
\references{

Mateus Maia & Anderson Ara (2023). rmachines: Random Machines: a package for a support vector ensemble based on random kernel space. R package version 0.1.0.

}
\examples{
    data(bolsafam)
    head(bolsafam)
}
\keyword{datasets}
