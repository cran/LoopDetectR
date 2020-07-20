
#' Solution for the cell cycle model related to func_li08
#' 
#' The file contains the solution over time (3 oscillatory cycles) for the 
#' ordinary differential equation model as given in func_li08. In addition, 
#' events as constructed in the 
#' original publication [Li et al., 2008] are considered. 
#' 
#' @format A dataframe with 634 rows and 19 columns
#' \describe{
#'   \item{time}{time variable}
#'   \item{y1}{first variable value}
#'   \item{y2}{second variable value, }
#'   \item{etc.}{etc. }
#'   \item{y18}{18th variable value}
#' }
#' 
#' @source The Caulobacter cell cycle model was proposed in 
#' Li S, Brazhnik P, Sobral B, Tyson JJ. A Quantitative Study of the 
#' Division Cycle of Caulobacter crescentus Stalked Cells. Plos Comput Biol. 
#' 2008;4(1):e9. The solutions were generated with MATLAB using the
#' functions accompanying the above reference on 
#' \url{http://mpf.biol.vt.edu/research/caulobacter/SWST/pp/}.
#' 
"li08_solution"



#' Example ODE function: chain model with positive regulation.
#'
#' The file contains the function definition an ordinary differential equation
#' model of a chain model of 4 variables with positive feedback. 
#'
#' @format R file with definition of function func_POSm4 that takes as
#' input arguments time t (dimension 1), variable values x (dimension 4), and 
#' kinetic parameter values klin (dimension 8) and knonlin (dimension 2).
#' 
#' @source 
#' The chain model was used in 
#' Baum K, Politi AZ, Kofahl B, Steuer R, Wolf J. Feedback, Mass 
#'  Conservation and Reaction Kinetics Impact the Robustness of Cellular 
#'  Oscillations. PLoS Comput Biol. 2016;12(12):e1005298.
#' 
"func_POSm4"



#' Example ODE function: bacterial cell cycle.
#' 
#' The file contains the function definition an ordinary differential equation
#' model of Caulobacter crescentus cell cycle as proposed by Li et al., 2008. 
#' It has 18 variables. 
#'
#' @format R file with definition of function func_li08 that takes as
#' input arguments time t (dimension 1), and variable values y (dimension 18).
#' the kinetic parameters are defined within the function.
#' 
#' @source The Caulobacter cell cycle model was proposed in 
#' Li S, Brazhnik P, Sobral B, Tyson JJ. A Quantitative Study of the 
#' Division Cycle of Caulobacter crescentus Stalked Cells. Plos Comput Biol. 
#' 2008;4(1):e9. The function corresponds to the MATLAB function modelwtin(t,y)
#' as given on \url{http://mpf.biol.vt.edu/research/caulobacter/SWST/pp/}.
#' 
#' @details
#' The Caulobacter cell cycle model function will only give the solution as 
#' shown in the publication [Li et al., 2008] if the change in variables at 
#' defined events are taken into account. Please refer to the original 
#' reference for details.
#' 
"func_li08"








