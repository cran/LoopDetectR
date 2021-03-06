% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{li08_solution}
\alias{li08_solution}
\title{Solution for the cell cycle model related to func_li08}
\format{
A dataframe with 634 rows and 19 columns
\describe{
  \item{time}{time variable}
  \item{y1}{first variable value}
  \item{y2}{second variable value, }
  \item{etc.}{etc. }
  \item{y18}{18th variable value}
}
}
\source{
The Caulobacter cell cycle model was proposed in 
Li S, Brazhnik P, Sobral B, Tyson JJ. A Quantitative Study of the 
Division Cycle of Caulobacter crescentus Stalked Cells. Plos Comput Biol. 
2008;4(1):e9. The solutions were generated with MATLAB using the
functions accompanying the above reference on 
\url{http://mpf.biol.vt.edu/research/caulobacter/SWST/pp/}.
}
\usage{
li08_solution
}
\description{
The file contains the solution over time (3 oscillatory cycles) for the 
ordinary differential equation model as given in func_li08. In addition, 
events as constructed in the 
original publication [Li et al., 2008] are considered.
}
\keyword{datasets}
