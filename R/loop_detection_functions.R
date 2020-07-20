
find_loops_noscc <- function(jacobian, max_num_loops=100000){
    #' @title Loop detection in a matrix
    #'
    #' @description Given the Jacobian matrix of an ODE system or the adjacency
    #'  matrix of a graph, this function determines all loops in the system up
    #'  to the maximal number supplied. No
    #'  decomposition into strongly connected components is performed.
    #'
    #' @details The input matrix delivers the directed interactions in the ODE
    #'  system; if entry \code{(i,j)} is non-zero it means that variable (or
    #'  node) \code{i} is regulated by variable (node) \code{j}. Johnson's
    #'  algorithm for path detections
    #'  from the igraph package (function:
    #'  \code{\link[igraph]{all_simple_paths}}) is used.
    #'  No decomposition into strongly connected components is employed which
    #'  could be beneficial for smaller systems (compared to
    #'  \code{\link{find_loops}}).
    #'  The queried graph is increased stepwise leading to the output of loops
    #'  in a certain order determined by the order of occurrence in the Jacobian
    #'  matrix:
    #'  \itemize{
    #'  \item first the self-loops,
    #'  \item then feedback loops incorporating only the first and second species of
    #'  the Jacobian,
    #'  \item then feedback loops incorporating the third and at most also the first
    #'  and second species of the jacobian, etc.}
    #'  If the maximal number of loops, max_num_loops, is reached, no warning is
    #'  issued. It is very probable that not all feedback loops of the system
    #'  have been found.
    #'  Up to which species this function searched before stopping due to
    #'  reaching the maximal allowed loop number can be inferred from the last
    #'  exported feedback loop.
    #'  Running the function multiple times with re-ordered jacobian as input
    #'  can enable detection of alternative feedback loops while limiting the
    #'  runtime and output size of single runs.
    #'  If columns of the Jacobian are named,
    #'  the identification is given by the attribute \code{node_ids},
    #'  \code{attr(result,"node_ids")}.
    #'
    #' @param jacobian Square Jacobian matrix of an ODE system or the adjacency
    #'  matrix of a graph; captures interactions such that entry \code{(i,j)} is
    #'  negative (positive) if variable \code{j} regulates variable \code{i}
    #'  negatively (positively).
    #'
    #' @param max_num_loops Positive numeric value indicating the maximal number
    #' of loops that are reported. Default: \eqn{10^5}.
    #'
    #' @return A data.frame with three columns: \code{loop}, \code{length}, \code{sign}
    #'  containing up to \code{max_num_loops} loops of the systems defined
    #'  by matrix \code{jacobian}. Each entry in the loop column is a list of
    #'  identifiers that correspond to the indices of the variable in the
    #'  Jacobian matrix and denote
    #'  in which order the variables form the loop.
    #'
    #' @examples
    #' #sample Jacobian matrix of a system with 4 variables
    #' jac_matrix <- rbind(c(-1,0,0,-1),c(1,-1,0,1),c(0,1,-1,0),c(0,0,1,-1))
    #' #find the first 5 feedback loops of the system
    #' loop_list <- find_loops_noscc(jac_matrix,5)
    #'
    #' @seealso  \code{\link{find_loops}},  \code{\link{find_loops_vset}}
    #'
    #' @export
    #'

  #initialize variables
  all_loops <- vector('list',max_num_loops) #list to store all loops
  loop_signs <- numeric(max_num_loops) #vector to store if the loops are positive or negative
  loop_length <- numeric(max_num_loops) #vector to store the length of all loops
  loopcount <- 0

  #determine the sign of the entries only once, this is important to avoid overflow for large entries
  sig_jacobian <- sign(jacobian)

  #determine self-loops
  for (i in which(diag(sig_jacobian)!=0)){ #iterate over all nonzero entries of the diagonal, these are self-loops
    loopcount <- loopcount+1
    loop_signs[loopcount] <- sig_jacobian[i,i] #sign of the self-loop
    loop_length[loopcount] <-  1
    all_loops[loopcount] <- list(c(i,i))
  }

  #only include these nodes into the examination for which we have at least one outgoing and one ingoing edge
  #this makes computation a lot faster for sparse matrices. Might increase runtime for smaller matrices.
  node_inds <- which(sapply(1:dim(sig_jacobian)[1],function(x){any(sig_jacobian[x,][-x]!=0) & any(sig_jacobian[,x][-x]!=0)}))
  node_count <- 1

  #determine loops with length > 1
  for (i in node_inds[-1]){ #initialize a loop from the first row/column of the jacobian matrix to the last
    node_count <- node_count + 1 #keeps track of the node count
    if (loopcount < max_num_loops){
      sign_tmp_jacobian <- sig_jacobian[node_inds[1:node_count],node_inds[1:node_count]] #store the current sub-jacobian in a temporary matrix
      if (i==node_inds[2]) { # for the first graph with more than one node
        #build the directed graph from the current jacobian matrix, assign the id in the jacobian to the nodes
        jac_graph <- igraph::graph_from_adjacency_matrix(t(abs(sign_tmp_jacobian)),mode="directed",diag=FALSE, add.colnames = NA)
        #find the nodes which affect the most recently added node (not the node itself)
        secondlastnodes <- which(sign_tmp_jacobian[node_count,1:(node_count-1)] != 0)
      } else {# if it is not the first node
        #just add the one node (this will be the node with index node_count)
        jac_graph <- igraph::add_vertices(jac_graph,nv=1)
        #species acting on node i, omit self-edges:
        spec_act_on_i <- which(sign_tmp_jacobian[node_count,1:(node_count-1)] != 0)
        #odes which affect the most recently added node are used further
        secondlastnodes <- spec_act_on_i
        #add the edges to the graph
        if (length(spec_act_on_i)>0){
          jac_graph <- igraph::add_edges(jac_graph, c(rbind(spec_act_on_i, node_count)))
        }
        #species on which node i acts:
        spec_i_acts_on <- which(sign_tmp_jacobian[1:(node_count-1),node_count]!=0)
        #add the edges to the graph
        if (length(spec_i_acts_on)>0){
          jac_graph <- igraph::add_edges(jac_graph, c(rbind(node_count, spec_i_acts_on)))
        }

      }

      #search all simple paths from the current node to the nodes affecting node i
      if (length(secondlastnodes)!=0){ #only iterate over the nodes if the list is not empty
        #loop over all nodes which have an edge to node i
        for (k in secondlastnodes) {
          #ensure that the limit of numbers of loops is not yet exceeded
          if (loopcount < max_num_loops){
            #find and store all paths from the current node i to node k
            tmp_path <- igraph::all_simple_paths(jac_graph,from=node_count,to=k)
            if (length(tmp_path)!=0) { # if there are simple paths found
              #iterate over all paths in the list
              for (j in 1:length(tmp_path)){
                #ensure that the limit of numbers of loops is not yet exceeded
                if (loopcount < max_num_loops){
                  #store the current loop in the variable
                  #(e.g.  k = 1, tmp_path[[j]] = 4 3 2 1 --> loop = 1 4 3 2 1)
                  loop <- c(k, igraph::as_ids(tmp_path[[j]]))
                  #count how many loops have been found
                  loopcount <- loopcount+1
                  #store the loop in the list of all paths (and reformulate to original node ids)
                  all_loops[loopcount] <- list(node_inds[loop])
                  #store the loop length (it is loop length-1, because the edges should be counted and not the nodes in the loop)
                  loop_length[loopcount] <- length(loop)-1
                  #edge indices of the loop path
                  edge_vec <- cbind(loop[-1],loop[-length(loop)])
                  #check the product of all weights along the path to determine the sign of the loop
                  loop_signs[loopcount] <- sign(prod(sign_tmp_jacobian[edge_vec]))
                } else { #stop execution if we reach the upper limit to the number of nodes
                  break
                }
              }
            }
          } else { #stop execution if we reach the upper limit to the number of nodes
            break
          }
        }
      }
    } else { #stop execution if we reach the upper limit to the number of nodes
      break
    }
  }

  #save results
  result <- data.frame(loop=I(all_loops[1:loopcount]), length=loop_length[1:loopcount],sign=loop_signs[1:loopcount]) #store the results in a dataframe
  #add node identifiers if column names are provided in the Jacobian
  if (!is.null(colnames(jacobian))){
    name_vec <- (1:dim(jacobian)[2])
    names(name_vec) <-  colnames(jacobian)
    attr(result,'node_ids') <- name_vec
  }
  return(result) #return the data frame
}



find_loops <- function(jacobian, max_num_loops=100000){
    #' @title Loop detection in a matrix
    #'
    #' @description Given the Jacobian matrix of an ODE system or the adjacency
    #'  matrix of a graph, this function determines all loops in the system up
    #'  to the maximal number supplied.
    #'
    #' @details The input matrix delivers the directed interactions in the ODE
    #'  system; if entry \code{(i,j)} is non-zero it means that variable (or node) \code{i} is
    #'  regulated by variable (node) \code{j}. Johnson's algorithm for path detection
    #'  as well as Tarjan's algorithm for detecting strongly connected
    #'  components are used as implemented in the igraph package (functions:
    #'  \code{\link[igraph]{all_simple_paths}}, \code{\link[igraph]{components}}) .
    #'  If the maximal number of loops, max_num_loops, is reached, no warning is
    #'  issued. It is very probable that not all feedback loops of the system
    #'  have been found.
    #'  Running the function multiple times with re-ordered jacobian as input
    #'  can enable detection of alternative feedback loops while limiting the
    #'  runtime and output size of single runs.
    #'  If columns of the Jacobian are named,
    #'  the identification is given by the attribute \code{node_ids},
    #'  \code{attr(result,"node_ids")}.
    #'
    #'
    #' @param jacobian Square Jacobian matrix of an ODE system or the adjacency
    #'  matrix of a graph; captures interactions such that entry \code{(i,j)} is
    #'  negative (positive) if variable \code{j} regulates variable \code{i} negatively
    #'  (positively).
    #'
    #' @param max_num_loops Positive numeric value indicating the maximal number
    #' of loops that are reported. Default: \eqn{10^5}.
    #'
    #' @return A data.frame with three columns: \code{loop}, \code{length}, \code{sign}
    #'  containing up to \code{max_num_loops} loops of the systems defined
    #'  by matrix \code{jacobian}. Each entry in the loop column is a list of
    #'  identifiers that correspond to the indices of the variable in the
    #'  Jacobian matrix and denote
    #'  in which order the variables form the loop.
    #'
    #' @examples
    #' #sample Jacobian matrix of a system with 4 variables
    #' jac_matrix <- rbind(c(-1,0,0,-1),c(1,-1,0,1),c(0,1,-1,0),c(0,0,1,-1))
    #' #find the first 5 feedback loops of the system
    #' loop_list <- find_loops(jac_matrix,5)
    #'
    #' @seealso  \code{\link{find_loops_noscc}},  \code{\link{find_loops_vset}}
    #'
    #' @export


  #initialize variables
  all_loops <- vector('list',max_num_loops)#list to store all loops
  loop_signs <- numeric(max_num_loops) #vector to store if the loops are positive or negative
  loop_length <- numeric(max_num_loops) #vector to store the length of all loops
  #number of detected loops
  loopcount <- 0

  #determine the sign of the entries only once, this is important to avoid overflow for large entries
  sig_jacobian <- sign(jacobian)

  #determine self-loops
  for (i in which(diag(jacobian)!=0)){ #iterate over all nonzero entries of the diagonal, these are self-loops
    loopcount <- loopcount+1
    loop_signs[loopcount] <- sig_jacobian[i,i] #sign of the self-loop
    loop_length[loopcount] <-  1
    all_loops[loopcount] <- list(c(i,i))
  }

  #only include these nodes into the examination for which we have at least one outgoing and one ingoing edge
  #this makes computation a lot faster for sparse matrices. Might increase runtime for smaller matrices.
  glob_node_inds <- which(sapply(1:dim(sig_jacobian)[1],function(x){any(sig_jacobian[x,][-x]!=0) & any(sig_jacobian[,x][-x]!=0)}))

  #determine strongly connected components in the associated graph
  jac_graph <- igraph::graph_from_adjacency_matrix(t(abs(sig_jacobian[glob_node_inds,glob_node_inds])),mode="directed",diag=FALSE, add.colnames = NA)
  str_comp <- igraph::components(jac_graph, mode="strong")

  #for all clusters of size>1
  for (clus_id in which(str_comp$csize>1)) {
    #determine node_inds
    node_inds <- glob_node_inds[str_comp$membership==clus_id]
    #
    node_count <- 1
    #determine loops with length > 1
    for (i in node_inds[-1]){ #initialize a loop from the first row/column of the jacobian matrix to the last
      node_count <- node_count + 1 #keeps track of the node count
      if (loopcount < max_num_loops){
        sign_tmp_jacobian <- sig_jacobian[node_inds[1:node_count],node_inds[1:node_count]] #store the current sub-jacobian in a temporary matrix
        if (i==node_inds[2]) { # for the first graph with more than one node
          #build the directed graph from the current jacobian matrix, assign the id in the jacobian to the nodes
          jac_graph <- igraph::graph_from_adjacency_matrix(t(abs(sign_tmp_jacobian)),mode="directed",diag=FALSE, add.colnames = NA)
          #find the nodes which affect the most recently added node (not the node itself)
          secondlastnodes <- which(sign_tmp_jacobian[node_count,1:(node_count-1)] != 0)
        } else {# if it is not the first node
          #just add the one node
          jac_graph <- igraph::add_vertices(jac_graph,nv=1)
          #species acting on node i, omit self-edges:
          spec_act_on_i <- which(sign_tmp_jacobian[node_count,1:(node_count-1)]!=0)
          #the nodes which affect the most recently added node are used further
          secondlastnodes <- spec_act_on_i
          #add the edges to the graph
          if (length(spec_act_on_i)>0){
            jac_graph <- igraph::add_edges(jac_graph, c(rbind(spec_act_on_i, node_count)))
          }
          #species on which node i acts:
          spec_i_acts_on <- which(sign_tmp_jacobian[1:(node_count-1),node_count]!=0)
          #add the edges to the graph
          if (length(spec_i_acts_on)>0){
            jac_graph <- igraph::add_edges(jac_graph, c(rbind(node_count, spec_i_acts_on)))
          }

        }

        #search all simple paths from the current node to the nodes affecting the ith node
        if (length(secondlastnodes)!=0){ #only iterate over the nodes if the list is not empty
          #loop over all nodes which have an edge to node i
          for (k in secondlastnodes) {
            #ensure that the limit of numbers of loops is not yet exceeded
            if (loopcount < max_num_loops){
              #find and store all paths from the current node node_count to node k
              tmp_path <- igraph::all_simple_paths(jac_graph,from=node_count,to=k)
              if (length(tmp_path)!=0) { # if there are simple paths found
                #iterate over all paths in the list
                for (j in 1:length(tmp_path)){
                  #ensure that the limit of numbers of loops is not yet exceeded
                  if (loopcount < max_num_loops){
                    #store the current loop in the variable
                    #(e.g.  k = 1, tmp_path[[j]] = 4 3 2 1 --> loop = 1 4 3 2 1)
                    loop <- c(k, igraph::as_ids(tmp_path[[j]]))
                    #count how many loops have been found
                    loopcount <- loopcount+1
                    #store the loop in the list of all paths (and reformulate to original node ids)
                    all_loops[loopcount] <- list(node_inds[loop])
                    #store the loop length (it is loop length-1, because the edges should be counted and not the nodes in the loop)
                    loop_length[loopcount] <- length(loop)-1
                    #edge indices of the loop path
                    edge_vec <- cbind(loop[-1],loop[-length(loop)])
                    #check the product of all weights along the path to determine the sign of the loop
                    loop_signs[loopcount] <- sign(prod(sign_tmp_jacobian[edge_vec]))
                  } else { #stop execution if we reach the upper limit to the number of nodes
                    break
                  }
                }
              }
            } else { #stop execution if we reach the upper limit to the number of nodes
              break
            }
          }
        }
      } else { #stop execution if we reach the upper limit to the number of nodes
        break
      }
    }
  }
  #save results
  result <- data.frame(loop=I(all_loops[1:loopcount]), length=loop_length[1:loopcount],sign=loop_signs[1:loopcount]) #store the results in a dataframe
  #add node identifiers if column names are provided in the Jacobian
  if (!is.null(colnames(jacobian))){
    name_vec <- (1:dim(jacobian)[2])
    names(name_vec) <-  colnames(jacobian)
    attr(result,'node_ids') <- name_vec
  }
  return(result) #return the data frame
}



find_loops_vset <- function(fun, vset, ..., max_num_loops=100000, compute_full_list=TRUE){
    #' @title Loop detection for an ODE model at multiple sets of variables
    #'
    #' @description Determines loop lists for an ODE system given by a function
    #' and at multiple sets of variables. Loop lists are reported if signs of
    #' Jacobian matrix have changed.
    #'
    #' @param fun Function defining the ODE system, returns the vector \eqn{dx/dt}.
    #' May depend on further parameters in \code{...}.
    #' @param vset List of variable values at which the loops are determined.
    #' @param ... Further parameters except variable values to the function \code{fun},
    #' none called \code{x}.
    #' @param max_num_loops Positive numeric value indicating the maximal number
    #' of loops that are reported in a loop list. Default: \eqn{10^5}.
    #' @param compute_full_list Logical value indicating whether for each
    #' Jacobian matrix with any different sign the loop list is computed (\code{TRUE},
    #' default), or whether further checks are performed to ensure that loops
    #' may be altered.
    #'
    #' @details
    #' The supplied function can take more arguments, but only the variables
    #' are allowed to be named \code{x} (they can also be named differently).
    #' The Jacobian matrices are computed for each of the variable
    #' values defined in vset using the \code{\link[numDeriv]{jacobian}} function
    #' from the \code{NumDeriv} package with option \code{method = 'complex'}, i.e. using a
    #' complex-step approach.
    #' If \code{compute_full_list = TRUE} (default), loop lists are not re-computed
    #' for Jacobians that clearly do not allow for altered loop lists. This is
    #' the case if no new regulation appear and only signs of regulations are
    #' altered that are not member of any loop. Loop lists can still be
    #' identical for different Jacobians, e.g. if two sign switches occur that
    #' are both affecting the same loops.
    #'
    #' @return A list with four entries:
    #' \itemize{
    #' \item \code{loop_rep} List of loop lists.
    #' \item \code{loop_rep_index} Vector of integer numbers returning the index of the
    #' loop list in loop_rep belonging to each entry in \code{vset}.
    #' \item \code{jac_rep} List of signed Jacobian matrices.
    #' \item \code{jac_rep_index} Vector of integer numbers returning the index of the
    #' Jacobian matrix in jac_rep belonging to each entry in \code{vset}.
    #' }
    #'
    #' @details
    #' If there is only one class of Jacobian matrix (i.e. the signs of the
    #' Jacobian matrix are the same for all entries in \code{vset}), \code{loop_rep} and
    #' \code{jac_rep} will have only one entry each. The number of entries for
    #' \code{loop_rep_index} and \code{jac_rep_index} corresponds to the length of \code{vset}.
    #' Only if \code{compute_full_list} is set to \code{FALSE}, \code{loop_rep} can contain
    #' fewer elements than \code{jac_rep}, otherwise both have the same number of
    #' elements.
    #'
    #' @examples
    #' #default call to determine loops from an ODE model given by a function
    #' #read in example functions
    #' data("func_POSm4")
    #' #the loaded function func_POSm4 takes arguments t, x, klin, knonlin
    #' res_tab <- find_loops_vset(func_POSm4,vset=list(c(1,1,1,1)),t=1,
    #' klin=c(1,2,0.5,1,2,0.1,3,2,3),knonlin=c(1,2))
    #' #computed loop list:
    #' res_tab$loop_rep[[1]] #or res_tab[[1]][[1]]
    #'
    #' #determine loops from an ODE model over the course of a solution
    #' #read in the example function defining the bacterial cell cycle
    #' data("func_li08")
    #' #kinetic parameter values are defined within the function
    #' #read in a set of variable values (the solution of func_li08 with events)
    #' data("li08_solution")
    #' #transform the solution (columns: variables) to the correct list format
    #' #and remove the time (first column)
    #' li08_sol_list <- as.list(as.data.frame(t(li08_solution[,-1])))
    #' res_tab <- find_loops_vset(func_li08,vset=li08_sol_list,t=1,
    #' compute_full_list=FALSE)
    #'
    #' @export

    #compute all jacobians
    jacob_list <- lapply(vset,function(y){
        sign(numDeriv::jacobian(fun,x=y,method="complex",...))})
    #find different signed Jacobians
    jac_rep <- unique(jacob_list)

    #identify the indices back
    jac_rep_index <- match(jacob_list,jac_rep)

    if (compute_full_list) { #if we compute loops for all Jacobians
        loop_rep <- lapply(jac_rep,function(y){
            find_loops(y,max_num_loops)
        })
        loop_rep_index <- jac_rep_index
    } else { #if we search for different patterns before
        #loop list for the first Jacobian
        loop_rep <- list(find_loops(jac_rep[[1]],max_num_loops))
        loop_rep_index <- rep(0,length(jac_rep_index))
        loop_rep_index[jac_rep_index==1] <- 1
        #we go through every Jacobian separately
        if (length(jac_rep)>1) { #if we have more than one Jacobian
        for (i in 2:length(jac_rep)) { #check each Jacobian (if there are more
            #than 1)
            J_temp <- jac_rep[[i]]
            #determine if there is a switch from zero to nonzero
            switch_to_nonzero <- vapply(jac_rep[1:i-1],function(x)
                {any(abs(J_temp[x==0])>0)},logical(1))
            existing_loop_changed <- rep(FALSE,i-1)
            for (j in which(switch_to_nonzero==0)) { #for the Jacobians in which
                #no additional regulation came up: determine if the sign
                #switches of an edge that is contained in a loop
                #target-source node pairs of altered Jacobian entries
                tn_sn_ind <- which(!((J_temp-jac_rep[[j]])==0),arr.ind=T)
                #find loops that contain any of the edges
                for (k in 1:length(tn_sn_ind)) {
                    if (length(find_edge(loop_rep[[j]],
                                          tn_sn_ind[k,2],tn_sn_ind[k,1]))>0){
                        #in case we found an edge that is affected
                        existing_loop_changed[j] <- TRUE
                        break
                    }
                }
            }
            #if there is a relevant difference for all existing Jacobians
            if (all(switch_to_nonzero|existing_loop_changed)) {
                loop_rep[[length(loop_rep)+1]] <-
                    find_loops(J_temp,max_num_loops)
                loop_rep_index[jac_rep_index==i] <- max(loop_rep_index)+1
            } else { #if there is one Jacobian for which no interaction came up
                # AND no edge of a loop changed - assign index of the first
                # similar Jacobian (in principle, there only should be one)
                loop_rep_index[jac_rep_index==i] <-
                    loop_rep_index[jac_rep_index==which(!(switch_to_nonzero|existing_loop_changed))[1]][1]

            }
        }
        }
    }
    return(list('loop_rep'=loop_rep,'loop_rep_index'=loop_rep_index,
                'jac_rep'=jac_rep,'jac_rep_index'=jac_rep_index))
}
