% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mst_fun.R
\name{mst_fun}
\alias{mst_fun}
\title{A function to update the model parameters.}
\usage{
mst_fun(m, V, G, K, z)
}
\arguments{
\item{m}{variational parameter m from last iteration.}

\item{V}{variational parameter Vfrom last iteration.}

\item{G}{number of component.}

\item{K}{number of taxa.}

\item{z}{class membership matrix at the current iteration.}
}
\description{
This function to update the model parameters for the latent Gaussian variables.
}
\examples{
mst_fun()

}
\keyword{mst}
