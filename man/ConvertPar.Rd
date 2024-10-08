% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConvertPar.R
\name{ConvertPar}
\alias{ConvertPar}
\alias{BiCopEta2Par}
\alias{BiCopPar2Eta}
\alias{BiCopEta2Tau}
\alias{BiCopTau2Eta}
\title{Conversions between various bivariate copula parametrizations.}
\usage{
BiCopEta2Par(family, eta, eta2 = 0)

BiCopPar2Eta(family, par, par2 = 0)

BiCopEta2Tau(family, eta, eta2 = 0)

BiCopTau2Eta(family, tau)
}
\arguments{
\item{family}{An integer defining the bivariate copula family to use.  See \strong{Details}.}

\item{eta, eta2}{Vector of parameters on the \code{eta} scale.  See \strong{Details}.}

\item{par, par2}{Vector of parameters on the \code{par} scale.}

\item{tau}{Vector of parameters on the \code{tau} scale.}
}
\value{
Vector of converted parameters.
}
\description{
Conversions between various bivariate copula parametrizations.
}
\details{
The copula family integer codes are identical to those of the \pkg{VineCopula} package. Currently, the following families are implemented:
\describe{
\item{1}{Gaussian copula.}
\item{2}{Student-t copula.}
\item{3}{Clayton copula.}
\item{4}{Gumbel copula.}
\item{5}{Frank copula.}
\item{13}{Clayton copula -- rotated 180 degrees.}
\item{14}{Gumbel copula -- rotated 180 degrees.}
\item{23}{Clayton copula -- rotated 90 degrees.}
\item{24}{Gumbel copula -- rotated 90 degrees.}
\item{33}{Clayton copula -- rotated 270 degrees.}
\item{34}{Gumbel copula -- rotated 270 degrees.}
}
}
