

func_POSm4 <- function(t,x,klin,knonlin){
    # example function: chain with positive regulation [Baum et al., 2016]
    #
    # Baum K, Politi AZ, Kofahl B, Steuer R, Wolf J. Feedback, Mass
    # Conservation and Reaction Kinetics Impact the Robustness of Cellular
    # Oscillations. PLoS Comput Biol. 2016;12(12):e1005298
    #
    dx <- rep(0,4)
    dx[1] <- klin[1]-(klin[2]*(1 + (x[4]/knonlin[1])^knonlin[2]) +
                          klin[3])*x[1]
    dx[2] <- klin[2]*(1 + (x[4]/knonlin[1])^knonlin[2])*x[1] -
        (klin[4] + klin[5])*x[2]
    dx[3] <- klin[4]*x[2] - (klin[6] + klin[7])*x[3];
    dx[4] <- klin[6]*x[3] - klin[8]*x[4];
    return(dx)
}
