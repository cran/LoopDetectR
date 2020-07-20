## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install_fromCRAN, eval=FALSE---------------------------------------------
#  # Download and install
#  install.packages("LoopDetectR")

## ----install_fromgit, message=FALSE, results='hide', eval=FALSE---------------
#  # Install the package including the vignette and the manual
#  remotes::install_gitlab("kabaum/LoopDetectR", build_manual=TRUE,
#                          build_vignettes = TRUE)
#  

## ----load_package-------------------------------------------------------------
# Load package
library("LoopDetectR")

## ----quick_example------------------------------------------------------------

# Load example ODE system with function func_POSm4, 4 variables
data("func_POSm4")
# Example variable values
s_star <- rep(1,4)
# Further arguments of func_POSm4, in addition: time t as argument
klin <- rep(1,8)
knonlin <- c(2.5,3)

# compute loops
res_tab <- find_loops_vset(func_POSm4,vset=list(s_star),t=1,klin=klin,
                           knonlin=knonlin,max_num_loops=10)
# The loop list is reported
res_tab$loop_rep[[1]]
# This is the sixth loop of the list. It is a positive feedback loop (sign in 
# the loop list equals +1) of length 3 in that variable 3 regulates variable 4, 
# variable 4 regulates variable 2, and variable 2 regulates variable 3. 
res_tab$loop_rep[[1]][6,]
# The corresponding signed Jacobian matrix
res_tab$jac_rep[[1]]

## ----solving_ODE--------------------------------------------------------------
# Load example ODE system with function func_POSm4, 
# Positive feedback chain model from [Baum et al., 2016], 4 variables
data(func_POSm4)
# The function func_POSm4 returns a vector, but deSolve needs the vector within
# a list as output. Therefore, we define a function that simply puts the output 
# of func_POSm4 into a list:
func_POSm4_list <- function(t,x,klin,knonlin){list(func_POSm4(t,x,klin,knonlin))}
# Kinetic parameters of the model, supplied as arguments to func_POSm4
klin <- c(165,0.044,0.27,550,5000,78,4.4,5.1)
knonlin <- c(0.3,2)
# Solve the system using deSolve
sol <-  deSolve::ode(y = rep(1,4), times = seq(0,15,0.1), func = func_POSm4_list, 
                     parms=klin, knonlin=knonlin)
# The solution of the 4-variable system is oscillatory, showing only the first 
# variable here
plot(sol[,1],sol[,2],type='l',xlab='time',ylab ='variable 1')
# Set the last point of the numeric solution as point of interest, omit the 
# first column (it contains the time)
s_star <- sol[dim(sol)[1],2:dim(sol)[2]]

## ----compute_jacobian---------------------------------------------------------
klin <- c(165,0.044,0.27,550,5000,78,4.4,5.1)
knonlin <- c(0.3,2)
j_matrix <- numDeriv::jacobian(func_POSm4,s_star,method="complex",
                               t=1,klin=klin,knonlin=knonlin,)
j_matrix
signed_jacobian <- sign(j_matrix)

## ----find_loops---------------------------------------------------------------
# Determine the loop_list from Jacobian matrix j_matrix
loop_list <- find_loops(j_matrix)
loop_list
# The signed Jacobian matrix can be supplied instead, delivering the same results
loop_list <- find_loops(signed_jacobian)
loop_list

## ----retrieve_loops-----------------------------------------------------------
# Retrieve fifth loop
loop_list[5,1] 
# In order to obtain the vector of the order of variables of a loop, you have to 
# call the list element
loop_list[5,1][[1]]
# Retrieve all loops of length 4
loop_list[loop_list$length==4,] 

## ----loop_summary-------------------------------------------------------------
loop_summary(loop_list)

## ----loops_with_specific_node-------------------------------------------------
# Index of node of interest
noi <- 2
# Return all loops from loop_list containing node 2
loop_list[vapply(loop_list$loop,function(x){noi %in% x},logical(1)),]

## ----find_edge----------------------------------------------------------------
# Obtain the indices of the loops with edge '2 regulates 3' in loop_list 
loop_edge_ind <- find_edge(loop_list,source_node=2,target_node=3);
loops_with_edge_2_to_3 <- loop_list[loop_edge_ind,]
loops_with_edge_2_to_3

## ----save_files---------------------------------------------------------------
# This writes a file 'looplist_func_POSm4.RData'
save(loop_list,file=file.path(tempdir(),'looplist_func_POSm4.RData')) 
# This loads the object 'loop_list' into the workspace
load(file.path(tempdir(),'looplist_func_POSm4.RData'))

# This saves the loop list into a table with separation by tab
write.table(loop_list,file=file.path(tempdir(),'looplist_func_POSm4.txt'),
            sep='\t',quote=F,row.names=F)
# This reads the loop list back into the R session, as object ll
ll <- read.table(file.path(tempdir(),'looplist_func_POSm4.txt'),
                 header=T,sep='\t')
ll
# This transforms the "c(x,y,z)"-entries into the format used by LoopDetectR
ll$loop <- lapply(ll$loop,function(x){eval(parse(text=x))})
ll

## ----li_example---------------------------------------------------------------
# Load in the example ODE model function func_li08, the kinetic parameters are
# defined within the function, and the function returns a vector of the time
# derivatives (in the same order as the modelled variables in the arguments)
data("func_li08")
# Load sets of variable values (solution to the ODE over time)
data("li08_solution") #loads the data.frame li08_solution, columns: variables
# Cast the solution into the correct list format and remove time (first column)
li08_sol_list <- as.list(as.data.frame(t(li08_solution[,-1])))
# Compute all different loop lists along the solution
res_tab <- find_loops_vset(func_li08,vset=li08_sol_list,t=1,
                           compute_full_list=FALSE)

## ----example_loop_lists_li----------------------------------------------------
loop_list_2 <- res_tab$loop_rep[[2]][res_tab$loop_rep[[2]]$length>1,]
loop_list_2
loop_list_7 <- res_tab$loop_rep[[7]][res_tab$loop_rep[[7]]$length>1,]
loop_list_7

## ----example_li_index---------------------------------------------------------
# Determine the loop list belonging to the 10th up to 20th input state. 
# It is loop_list_2 for all of these input states.
res_tab$loop_rep_index[10:20]

## ----compare_loops------------------------------------------------------------
# Load the ODE function of the positive feedback chain model
data(func_POSm4)
# Set kinetic parameters
klin <- c(165,0.044,0.27,550,5000,78,4.4,5.1)
knonlin <- c(0.3,2)
# Compute the Jacobian matrix of the system at the state [1,1,1,1]
j_matrix <- numDeriv::jacobian(func_POSm4,rep(1,4),method="complex",
                               t=1, klin=klin, knonlin=knonlin)
# Compute all loops for this Jacobian
loop_list_pos <- find_loops(j_matrix)

# Function with negative regulation. The altered regulation affects 
# two entries of the Jacobian matrix. Parameter values and the set of 
# variable values remain identical.
j_matrix_neg <- sign(j_matrix)
j_matrix_neg[1:2,4] <- -sign(j_matrix[1:2,4])
# Compute the loop list for this Jacobian with altered regulation
loop_list_neg <- find_loops(j_matrix_neg);

# Compute comparison
comp_inds <- compare_loop_list(loop_list_pos,loop_list_neg)
# Only the four self-loops remain identical in both systems (saved in ind_a_id).
loop_list_pos[comp_inds$ind_a_id,]
# ind_b_id saves the indices of the corresponding loops in the negatively 
# regulated system.
loop_list_neg[comp_inds$ind_b_id,]
# Two loops are the same in both systems but they have switched their signs. 
# Their indices in the first list are saved in ind_a_switch
loop_list_pos[comp_inds$ind_a_switch,]
# Their indices in the second list are saved in ind_b_switch 
loop_list_neg[comp_inds$ind_b_switch,]
# All loops in the positively regulated system do also occur in the negatively
# regulated system, i.e. ind_a_notin is empty.
comp_inds$ind_a_notin

