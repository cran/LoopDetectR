


func_li08 <- function(t,y) {
    # example function: bacterial cell cycle [modelwtin(t,y), Li et al. 2008]
    # Note: To obtain the solution as published, also events have to be 
    # considered, i.e. certain conditions lead to a change in certain variable 
    # values; see Li et al., 2008 for details.
    # 
    # Li S, Brazhnik P, Sobral B, Tyson JJ. A Quantitative Study of the 
    # Division Cycle of Caulobacter crescentus Stalked Cells. Plos Comput Biol. 
    # 2008;4(1):e9.
    #
    
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
    #Parameters values for the equations

ksCtrAP1  <- 0.0083 
JiCtrACtrA  <- 0.4 	
niCtrACtrA  <- 2 
ksCtrAP2  <- 0.073 	
JaCtrACtrA  <- 0.45 	
naCtrACtrA  <- 2 
kdCtrA1  <- 0.002 
kdCtrA2  <- 0.15 	
ndCtrA2  <- 2 
JdCtrADivKP  <- 0.55 
ksGcrA  <- 0.045 
JiGcrACtrA  <- 0.2 	
niGcrACtrA  <- 2 
kdGcrA  <- 0.022 
ksFts  <- 0.063 	
kdFts  <- 0.035 
kzringopen  <- 0.8 	
Jaopen  <- 0.01 
kzringclosed1  <- 0.0001 	
Jaclosed1  <- 0.1 
kzringclosed2  <- 0.6 	
nzringclosed2  <- 4 
JZringFts  <- 0.78 
ksDivK  <- 0.0054 	
ktransDivKP  <- 0.0295 	
ktransDivK  <- 0.5 	
kdDivK  <- 0.002 
ksI  <- 0.08 	
kdI  <- 0.04 
ksCcrM  <- 0.072 	
kdCcrM  <- 0.07 
kaDnaA  <- 0.0165 	
JiDnaAGcrA  <- 0.5 	
niDnaAGcrA  <- 2 
kdDnaA  <- 0.007 
kaIni  <- 0.01 	
JaIni  <- 1 		
naIni  <- 4 
thetaCtrA  <- 0.2 	
nthetaCtrA  <- 4 
thetaDnaA  <- 0.6 	
nthetaDnaA  <- 4 
thetaGcrA  <- 0.45 	
nthetaGcrA  <- 4 
thetaCori  <- 0.0002 	
nthetaCori  <- 1 
kmcori  <- 0.4 	
Jmcori  <- 0.95 	
nmcori  <- 4 
kmccrM  <- 0.4 	
JmccrM  <- 0.95 	
nmccrM  <- 4 
kmctrA  <- 0.4 	
JmctrA  <- 0.95 	
nmctrA  <- 4 
kmfts  <- 0.4 	
Jmfts  <- 0.95 	
nmfts  <- 4 
kelong  <- 0.95/160 
nelong  <- 4 

#$end of parameters
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
    ##Differential equations of the model
##variable map
#y[1]  <- [CtrA]
#y[2]  <- [GcrA]
#y[3]  <- [DnaA]
#y[4]  <- [Fts]
#y[5]  <- [Zring]
#y[6]  <- [DivK]
#y[7]  <- [DivK~P]
#y[8]  <- [Total DivK]
#y[9]  <- [I] (Intermediate)]
#y[10]  <- [CcrM]
#y[11]  <- [hCori] (hemimethylated Cori)
#y[12]  <- [hccrM] (hemimethylated ccrM)
#y[13]  <- [hctrA] (hemimethylated ctrA)
#y[14]  <- [hfts](hemimethylated fts)
#y[15]  <- [Ini] (Initiation)
#y[16]  <- [Elong](Elongation)
#y[17]  <- [DNA] (Total DNA)
#y[18]  <- Count (# of Chromosome) 
    
    dydt  <- rep(0,18)   
    dydt[1]  <- (ksCtrAP1*JiCtrACtrA^niCtrACtrA/(JiCtrACtrA^niCtrACtrA+y[1]^niCtrACtrA)*y[2]+ksCtrAP2*y[1]^naCtrACtrA/(JaCtrACtrA^naCtrACtrA+y[1]^naCtrACtrA))*y[13]-(kdCtrA1+kdCtrA2*y[7]^ndCtrA2/(JdCtrADivKP^ndCtrA2+y[7]^ndCtrA2))*y[1] 
    dydt[2]  <- (ksGcrA*JiGcrACtrA^niGcrACtrA/(JiGcrACtrA^niGcrACtrA+y[1]^niGcrACtrA)*y[3]-kdGcrA*y[2]) 
    dydt[3]  <- kaDnaA*JiDnaAGcrA^niDnaAGcrA/(JiDnaAGcrA^niDnaAGcrA+y[2]^niDnaAGcrA)*y[1]*(2-y[11])-kdDnaA*y[3] 
    dydt[4]  <- ksFts*y[1]*y[14]-kdFts*y[4] 
    dydt[5]  <- (kzringopen*(1-y[5])/(0.01+(1-y[5]))-(kzringclosed1+kzringclosed2*(y[4]/JZringFts)^nzringclosed2)*y[5]/(0.05+y[5])) 
    dydt[6]  <- (ksDivK*y[1]+ktransDivKP*y[7]-ktransDivK*(1-y[5])*y[6]-kdDivK*y[6]) 
    dydt[7]  <- (-ktransDivKP*y[7]+ktransDivK*(1-y[5])*y[6]-kdDivK*y[7]) 
    dydt[8]  <- (ksDivK*y[1]-kdDivK*y[8]) 
    dydt[9]  <- ksI*y[12]*y[1]-kdI*y[9] 
    dydt[10]  <- ksCcrM*y[9]-kdCcrM*y[10] 
    dydt[11]  <- -kmcori*y[10]^nmcori/(Jmcori^nmcori+y[10]^nmcori)*y[11] 
    dydt[12]  <- -kmccrM*y[10]^nmccrM/(JmccrM^nmccrM+y[10]^nmccrM)*y[12] 
    dydt[13]  <- -kmctrA*y[10]^nmctrA/(JmctrA^nmctrA+y[10]^nmctrA)*y[13] 
    dydt[14]  <- -kmfts*y[10]^nmfts/(Jmfts^nmfts+y[10]^nmfts)*y[14] 
    dydt[15]  <- kaIni*(y[3]/thetaDnaA)^nthetaDnaA*(y[2]/thetaGcrA)^4/(JaIni^naIni+(y[1]/thetaCtrA)^nthetaCtrA+(y[3]/thetaDnaA)^nthetaDnaA+(y[2]/thetaGcrA)^nthetaGcrA+(y[11]/thetaCori)^nthetaCori) 
    dydt[16]  <- kelong*y[16]^nelong/(y[16]^nelong+0.05^nelong)*y[18] 
    dydt[17]  <- kelong*y[16]^nelong/(y[16]^nelong+0.05^nelong)*y[18] 
    dydt[18]  <- 0  #count (of chromosome - is only altered at an event
    
    return(dydt)
##end of equations
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
}