% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/custom_pi_prior.R
\name{custom_pi_prior}
\alias{custom_pi_prior}
\title{Create a vector of pi Dirichlet priors with specified values for one or more collections}
\usage{
custom_pi_prior(P, C)
}
\arguments{
\item{P}{A data frame of one or more desired pi prior parameters. One column, "collection",
is a character vector, with valid values including the names of any reference collections,
or the special value "DEFAULT_PI". The second column, "pi_param" is the prior value to be
used for each collection.}

\item{C}{a tibble with a column "collection" collection names}
}
\description{
This handles a case in which the user provides a data frame for \code{pi_prior}. The
data frame lists desired Dirichlet parameter priors for at least one reference collection,
and/or a default value for all unspecified collections.
}
\details{
Input checking is currently done in the early stages of \code{infer_mixture} in order to
throw errors before long processing times, and avoid re-checking during \code{bootstrap_rho}.
}
\keyword{internal}
