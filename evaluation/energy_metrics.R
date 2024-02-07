# NES = mean(abs(samples - true))**beta - 0.5 * mean(abs(samples - samples[perm]))**beta
# If sample is vector need euclid distance instead of -

score <- function(sample, value, perm, beta){
    if (is.vector(sample)){
        permu <- sample[perm]
    } else {
        permu <- sample[perm,]
    }
    # Find the ES for a single value
    within <- mean(rowSums(t(t(sample) - permu)**2)**0.5)**beta
    between <- mean(rowSums(t(t(sample)-value)**2)**0.5)**beta
    return (between - 0.5*within)
}

energy_score <- function(Pdist, Qdist, beta=1){
    # ES is a strictly proper scoring rule, that is ES(P, Q) = ES(Q, Q) iff 
    # P = Q
    if (is.vector(Pdist)){
        n <- length(Pdist)
        m <- length(Qdist)
    } else {
        n <- dim(Pdist)[1]
        m <- dim(Qdist)[1]
    }
    # Generate permutation vector
    perm <- sample(c(1:n))
    escore <- 0
    for (i in c(1:m)){
        if (is.vector(Pdist)){
            # Single param
            escore = escore + score(Pdist, Qdist[i], perm, beta)
        } else {
            # Multiple params
            escore = escore + score(Pdist, Qdist[i,], perm, beta)
        }
    }
    return(escore/m)
}

energy_stat <- function(Pdist, Qdist){
    # Similar to the score above, however this is used for a null-hypothesis:
    # H_0 = P and Q are from the same distribution iff energy_stat == 0
    # Setup holder for three terms
    between <- 0
    withP <- 0
    withQ <- 0
    if (is.vector(Pdist)){
        n <- length(Pdist)
        m <- length(Qdist)
        for (i in c(1:max(n,m))){
            if (i <= n){
                withP <- withP + rowSums(((t(t(Pdist) - Pdist[i])**2)**0.5))
            }
            if (i <= m){
                withQ <- withQ + rowSums(((t(t(Qdist) - Qdist[i])**2)**0.5))
                between <- between + rowSums(((t(t(Pdist)-Qdist[i])**2)**0.5))
            }
        }
    } else {
        n <- dim(Pdist)[1]
        m <- dim(Qdist)[1]
        for (i in c(1:max(n,m))){
            if (i <= n){
                withP <- withP + rowSums(((t(t(Pdist) - Pdist[i,])**2)**0.5))
            }
            if (i <= m){
                withQ <- withQ + rowSums(((t(t(Qdist) - Qdist[i,])**2)**0.5))
                between <- between + rowSums(((t(t(Pdist)-Qdist[i,])**2)**0.5))
            }
        }
    }
    return (2*sum(between)/(n*m) - sum(withP)/(n*n) - sum(withQ)/(m*m))
}
