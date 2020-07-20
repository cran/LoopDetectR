
sort_loop_index <- function(loop_list){
    #' Sort loop indices
    #'
    #' Changes the loop representation such that every loop starts with the
    #' smallest node index. Returns a loop list of the same dimensions, only
    #' column \code{loop} will be altered.
    #'
    #' @param loop_list Dataframe with a column \code{loop} that contains the lists
    #' of loops, e.g. obtained from \code{find_loops()}.
    #'
    #' @examples
    #' #sample Jacobian matrix of a system with 4 variables
    #' jac_matrix <- rbind(c(-1,0,0,-1),c(1,-1,0,1),c(0,1,-1,0),c(0,0,1,-1))
    #' #find the feedback loops of the system
    #' loop_list <- find_loops(jac_matrix,10)
    #' #sort the loop indices to start with the smallest
    #' sorted_loop_list <- sort_loop_index(loop_list)
    #'
    #' @seealso \code{\link{compare_loop_list}}
    #'
    #' @export

    sorted_loop_list <- loop_list
    sorted_loops <- lapply(sorted_loop_list$loop,function(x){
        ind_min <- which.min(x)
        if (ind_min>1){
            return(x[c(ind_min:length(x),(2:ind_min))])
        } else {
            return(x)
        }
    })
    sorted_loop_list$loop <- sorted_loops
    return(sorted_loop_list)
}



find_edge <- function(loop_list,source_node,target_node){
    #' Detecting loops with a certain edge
    #'
    #' @description Finds those loops in a loop list that contain a regulation
    #' from a certain variable (source node) to a certain variable (target node).
    #'
    #' @param loop_list Dataframe with a column \code{loop} that contains the lists
    #' of loops, e.g. obtained from \code{\link{find_loops}}.
    #' @param source_node Index of the variable that is the source of the
    #' queried interaction, i.e. that regulates the target node.
    #' @param target_node Index of the variable that is the target of the
    #' queried interaction, i.e. that is regulated by the source node.
    #'
    #' @return A vector that gives the indices in the loop list of those loops
    #' that contain the indicated edge.
    #'
    #' @examples
    #' #sample Jacobian matrix of a system with 4 variables
    #' jac_matrix <- rbind(c(-1,0,0,-1),c(1,-1,0,1),c(0,1,-1,0),c(0,0,1,-1))
    #' #find the feedback loops of the system
    #' loop_list <- find_loops(jac_matrix,10)
    #' #find the loops containing the regulation from variable 3 to variable 4
    #' inds_3_to_4 <- find_edge(loop_list,3,4)
    #'
    #' @export

    vec_is_edge <- vapply(loop_list$loop,function(x){
        source_ind <- which(x==source_node) #find index of source node
        if (length(source_ind)==0){ #if the source node was not found
            return(FALSE)
        } else {
            if (target_node==x[source_ind[1]+1]){ #check whether the subsequent
                #node of the source node is the target node
                return(TRUE)
            } else {
                return(FALSE)
            }
        }},logical(1))

    return(which(vec_is_edge))
}



loop_summary <- function(loop_list,column_val='length'){
    #' Summary of a loop list
    #'
    #' Summarizes the loops in a loop list by their length and sign, returns an
    #' overview table of the numbers of all, negative and positive loops divided
    #' by their lengths.
    #'
    #' @details Lengths are abbreviated by \code{len_1}, \code{len_2},
    #' \code{len_3} etc., signs are
    #' abbreviated by \code{pos} for positive, \code{neg} for negative loops. The table
    #' contains entries for each loop length from 1 to the maximal loop length
    #' encountered in the table, and zeros are filled in if no loops of a
    #' certain length exist in the table.
    #'
    #' @param loop_list List of loops as dataframe with columns \code{length},
    #'  \code{sign}.
    #' @param column_val String indicating the orientation of the summary table.
    #' By default, rows of the results table are the sign of the loops, columns
    #' are loop lengths. If \code{column_val} is set to \code{"sign"}, columns
    #' and rows are exchanged.
    #'
    #' @examples
    #' #sample Jacobian matrix of a system with 4 variables
    #' jac_matrix <- rbind(c(-1,0,0,-1),c(1,-1,0,1),c(0,1,-1,0),c(0,0,1,-1))
    #' #find the feedback loops of the system
    #' loop_list <- find_loops(jac_matrix,10)
    #' #loop summary table
    #' loop_sum_tab <- loop_summary(loop_list)
    #'
    #' @export

    res_mat <- table(loop_list$length)
    loop_tab <- data.frame(
        matrix(0,ncol=as.numeric(names(res_mat))[length(names(res_mat))],nrow=3),
        row.names=c('all','pos','neg'))
    loop_tab[1,as.numeric(row.names(res_mat))] <- res_mat
    res_mat_temp <- table(loop_list$length[loop_list$sign==1])
    loop_tab[2,as.numeric(row.names(res_mat_temp))] <- res_mat_temp
    res_mat_temp <- table(loop_list$length[loop_list$sign==-1])
    loop_tab[3,as.numeric(row.names(res_mat_temp))] <- res_mat_temp
    colnames(loop_tab) <- paste0('len_',1:length(loop_tab))

    if (column_val == 'sign') {
        loop_tab <- t(loop_tab[2:3,])
    }

    return(loop_tab)
}



compare_loop_list <- function(loop_list_a, loop_list_b){
    #' Compare two loop lists
    #'
    #' Compared two loop lists and returns the indices of those loops that are
    #' identical in both lists, that switch only the sign or that do not occur
    #' in both lists
    #'
    #' @details Indices of loops are given with respect to the order of the
    #' loops in the first supplied loop list as well as for the second loop
    #' list. The loops are sorted to represent their loops starting from the
    #' smallest variable index (using the function \code{\link{sort_loop_index}}).
    #'
    #' @param loop_list_a,loop_list_b Loop lists with columns \code{loop} and
    #' \code{sign}, for example generated from \code{\link{find_loops}}.
    #'
    #' @return A list with 5 (possible empty) vectors as entries.
    #' \itemize{
    #' \item \code{ind_a_id} - indices of the loops in the first loop list that occur
    #' identically in the second loop list
    #' \item \code{ind_a_switch} - indices of the loops in the first loop list that occur
    #' in the second loop list with a different sign
    #' \item \code{ind_a_notin} - indices of the loops in the first loop list that do not
    #' occur in the second loop list
    #' \item \code{ind_b_id} - indices of loops in the second loop list corresponding to
    #' the loops reported in \code{ind_a_id}
    #' \item \code{ind_b_switch} - indices of loops in the second loop list corresponding
    #' to loops reported in \code{ind_a_switch}.
    #' }
    #'
    #' @examples
    #' #sample Jacobian matrix of a system with 4 variables
    #' jac_matrix <- rbind(c(-1,0,0,-1),c(1,-1,0,1),c(0,1,-1,0),c(0,0,1,-1))
    #' #find the feedback loops of the system
    #' loop_list <- find_loops(jac_matrix,10)
    #' #a slightly different Jacobian matrix of the system with 4 variables
    #' jac_matrix_alt <- rbind(c(-1,0,0,1),c(1,-1,0,-1),c(0,1,-1,0),c(0,0,1,-1))
    #' #find the feedback loops of the system
    #' loop_list_alt <- find_loops(jac_matrix_alt,10)
    #' #compare the loop lists
    #' comp_loop_list <- compare_loop_list(loop_list,loop_list_alt)
    #' #loops that switch sign
    #' comp_loop_list[['ind_a_switch']]
    #'
    #' @export

    sorted_loop_list_a <- sort_loop_index(loop_list_a)
    sorted_loop_list_b <- sort_loop_index(loop_list_b)

    inds_a_in_b <- match(sorted_loop_list_a$loop,sorted_loop_list_b$loop)
    ind_a_notin <- which(is.na(inds_a_in_b))
    ind_a_in <- which(!is.na(inds_a_in_b))
    #check whether signs are the same in the matching loops
    logvec_a_id <- sorted_loop_list_a$sign[ind_a_in] == sorted_loop_list_b$sign[inds_a_in_b[ind_a_in]]

    return(list('ind_a_id' = ind_a_in[logvec_a_id],
                    'ind_a_switch' = ind_a_in[!logvec_a_id],
                    'ind_a_notin' = ind_a_notin,
                    'ind_b_id' = inds_a_in_b[ind_a_in[logvec_a_id]],
                    'ind_b_switch' = inds_a_in_b[ind_a_in[!logvec_a_id]]))

}
