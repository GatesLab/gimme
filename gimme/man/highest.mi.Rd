% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highest.mi.R
\name{highest.mi}
\alias{highest.mi}
\title{Identifies highest MI from list of MIs.}
\usage{
highest.mi(
  mi_list,
  indices,
  elig_paths,
  prop_cutoff,
  n_subj,
  chisq_cutoff,
  allow.mult,
  ms_tol,
  hybrid,
  dir_prop_cutoff
)
}
\arguments{
\item{mi_list}{A list of MIs across individuals}

\item{indices}{A list of fit indices. Only relevant at the individual-level.}

\item{elig_paths}{A character vector containing eligible paths that
gimme is allowed to add to a model (e.g., no nonsense paths).}

\item{prop_cutoff}{The proportion of individuals for whom a path must be
significant in order for it to be added to the models. NULL if used 
at the individual-level.}

\item{n_subj}{The number of subjects in a given stage of the search. If
in the group stage, n_subj equals the number of subjects. If in the subgroup
stage, n_subj equals the number of individuals in a given subgroup. At the 
individual stage, n_subj = 1.}

\item{chisq_cutoff}{Cutoff used in order for MI to be considered significant.
Value varies depending on stage of search (e.g., group, subgroup, 
individual).}
}
\value{
Returns name of parameter associated with highest MI. If no MI meets 
the criteria, returns NA.
}
\description{
Identifies highest MI from list of MIs.
}
\keyword{internal}
