% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATEN.R
\name{evolTree}
\alias{evolTree}
\title{Compute the output of the tree if the data matrix contain less columns}
\usage{
evolTree(tree, input, nodes)
}
\arguments{
\item{tree}{An And/Or tree}

\item{input}{a data matrix, each row is a time point, each column corresponds to a node/gene of the tree}

\item{nodes}{the number of nodes in the Boolean network}
}
\value{
a vector of the output of the tree
}
\description{
A helper function for computing the output of the tree according to the Boolean functions (And/Or tree) and the data matrix
}
