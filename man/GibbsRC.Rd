% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GibbsRC.R
\name{GibbsRC}
\alias{GibbsRC}
\title{Gibbs Sampler for the Random Cluster model}
\usage{
GibbsRC(t, p, q, n.nodes, burnin = 50)
}
\arguments{
\item{t}{The number of iterations to run the algorithm}

\item{p}{The wiring probability}

\item{q}{The clustering coefficient of the Random Cluster model}

\item{n.nodes}{The number of nodes in the configuration}

\item{burnin}{Burnin period for the markov chain}
}
\value{
The sampled configurations; minus the burnin
}
\description{
Function to sample configurations from the Random cluster model using a Gibbs sampler.
}
