library(rms)
model <- function(f, x, do.print=TRUE,dropfirst=TRUE) {
#    print(f)
    if(!("formula" %in% class(f)))
        f <- as.formula(f)
    mod <- ols(f, x=TRUE, y=TRUE, data=x)
    mod.r <- robcov(mod, cluster=x$baitID)
    if(do.print)
        return(modelprint(mod.r,dropfirst=dropfirst))
    return(mod.r)
}
modelpred <- function(f, x) {
#    print(f)
    mod <- ols(f, x=TRUE, y=TRUE, data=x)
    mod.r <- robcov(mod, cluster=x$baitID)
}
modelprint <- function(mod.r,dropfirst=TRUE) {
    cf <- coefficients(mod.r)
    ## A <- anova(mod.r)
    ## use <- 1:(nrow(A)-2)
    ## P <- split(A[use,"P"], rownames(A[use,]))
    ## DF <- split(A[use,"d.f."], rownames(A[use,])) %>% unlist()
    ## if(any(DF>1)) {
    ##     for(i in which(DF>1)) {
    ##         P[[i]] <- c(P[[i]],rep(NA,DF[i]-1))
    ##     }
    ## }
    ## p <- p[ 1:(length(p)-2) ]
    se <- ( vcov(mod.r) %>% diag() %>% sqrt() )
    z <- cf/se
    p <- pnorm(abs(z),lower.tail=FALSE) * 2
    ret <- cbind(Coef=cf, Lower.95=cf-1.96*se, Upper.95=cf+1.96*se, P=p)
    if(dropfirst)
        return(ret[-1,])
    return(ret)
}
modelcmp <- function(L, x, do.print=TRUE) {
    L <- lapply(L, as.formula)
    mods <- lapply(L,function(l) { #message(l);
        tryCatch(ols(l,x=TRUE, y=TRUE, data=x))})
    mods.r <- tryCatch(lapply(mods,robcov, cluster=x$baitID))
    if(do.print) {
        lapply(mods.r,modelprint) %>% print()
        message("BIC:")
    }
    #do.call("anova",mods.r) %>% print()
    sapply(mods.r,BIC)
}
lmodel <- function(f, x, do.print=TRUE) {
#    print(f)
    mod <- lrm(f, x=TRUE, y=TRUE, data=x)
    mod.r <- robcov(mod, cluster=x$baitID)
    if(do.print)
        return(modelprint(mod.r))
    return(mod.r)
}
