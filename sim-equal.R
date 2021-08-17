#!/usr/bin/env Rscript
options(scipen=999)
args <- commandArgs(TRUE)
if(length(args)!=9) {
  cat("usage: <#points-to-simulate> <v=alt-freq> <q-concentration> <p-concentration> <#DNA-reps> <#RNA-reps> <effect size> <total-DNA-reads> <total-RNA-reads>\n");
  q(status=1)
}
numPoints <- as.numeric(args[1])
v <- as.numeric(args[2])
Q_CONCENTRATION <- as.numeric(args[3])
P_CONCENTRATION <- as.numeric(args[4])
NUM_DNA_REPS <- as.integer(args[5])
NUM_RNA_REPS <- as.integer(args[6])
THETA <- as.numeric(args[7])
TOTAL_DNA <- as.integer(args[8])
TOTAL_RNA <- as.integer(args[9])
    
#TOTAL_DNA <- as.integer(TOTAL_DNA/NUM_DNA_REPS)
#TOTAL_RNA <- as.integer(TOTAL_RNA/NUM_RNA_REPS)

betaModeConc <- function(mu,conc) {
    alpha <- mu*(conc-2)+1
    beta <- (1-mu)*(conc-2)+1
    return(rbeta(1,alpha,beta))
}

DNAtotals <- rep(TOTAL_DNA/NUM_DNA_REPS,NUM_DNA_REPS)
RNAtotals <- rep(TOTAL_RNA/NUM_RNA_REPS,NUM_RNA_REPS)

# Outputs:  v p q theta qj DNAreps (DNAref,DNAalt) RNAreps (RNAref,RNAalt)
for(i in 1:numPoints) {
    # Generate unobserved parameters
    p <- 0
    while(p==0 || p==1) { p <- betaModeConc(v,P_CONCENTRATION) }
    theta <- THETA
    q <- theta*p/(1-p+theta*p)
    cat(paste(v,round(p,6),round(q,6),round(theta,6),sep="\t"))

    # Generate DNA & RNA counts for all replicates
    cat(paste("\t",NUM_DNA_REPS))
    for(i in 1:NUM_DNA_REPS) {
        #DNA_N <- TOTAL_DNA
        DNA_N <- DNAtotals[i]
        DNAalt <- rbinom(1,DNA_N,p)
        DNAref <- DNA_N-DNAalt
        cat(paste("\t",DNAref,"\t",DNAalt))
    }
    cat(paste("\t",NUM_RNA_REPS))
    #debug=vector()
    for(i in 1:NUM_RNA_REPS) {
        qj <- betaModeConc(q,Q_CONCENTRATION)
        #RNA_N <- TOTAL_RNA
        RNA_N <- RNAtotals[i]
        RNAalt <- rbinom(1,RNA_N,qj)
        RNAref <- RNA_N-RNAalt
        cat(paste("\t",RNAref,"\t",RNAalt))
        #debug=c(debug,RNAref,RNAalt)
    }
    #print(debug)
    #print(sum(debug))
    cat("\n")
}

