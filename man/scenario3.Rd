% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Scenarios.R
\name{scenario3}
\alias{scenario3}
\title{Scenario-3:  Large location shift, and small covariance difference from an additive logistics normal distribtion}
\usage{
scenario3(n1 = 30, n2 = 30, dms_ = 75, seed)
}
\arguments{
\item{n1}{number of sample in class 1}

\item{n2}{number of samples in class 2}

\item{dms_}{number of features}

\item{seed}{random seed}
}
\value{
A data frame with simulated data from two classes
}
\description{
Simulate data for two classes using scenario-2(Hinton (2021)). Data are sampled from the additive logistics normal distribution (see Aitichinson (1986))
}
\examples{
\dontrun{
scenario3(n1=40,n2=40,dms_ = 75)
}


}
\references{
Aitchinson, J 1986
}
