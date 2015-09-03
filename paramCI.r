gpdret = function (z, p)
{
    thres <- z$thres
    sig <- z$mle[1]
    xi <- z$mle[2]
    rate <- z$rate
    npy <- z$npy
    thres + (sig/xi) * ((p * npy * rate)^xi - 1)
}

paramCI = function (z, m, conf = 0.95, nint = 100, rl.xup = NULL, rl.xlow = NULL,
    xi.xup = NULL, xi.xlow = NULL, rl.only = FALSE, xi.only = FALSE,
    make.plot = FALSE)
{
#Function to calculate confidence intervals for given return
#period, using the Profile Likelihood method.
#This function is called by Function 'cnfInt.r'.
#Function developed by AS from the orginal presented in
#'gpd.parameterCI' (package 'ismev')  11/04/08:3:41 PM
    library("ismev")
    retp_ID <- paste("RP (yrs) = ",m)        #Given RP
    lmts <- list()
    xdat <- z$data
    ma <- -z$nllh
    u <- z$threshold
    la <- z$rate
    npy <- z$npy
    v <- numeric(nint)
    eps <- 1e-06
    if (make.plot) {
        if (!rl.only & !xi.only)
            par(mfrow = c(2, 1))
        else par(mfrow = c(1, 1))
    }
    if (!xi.only) {     #Calc. CI for given return period:
        est.rl.xup <- is.null(rl.xup)
        est.rl.xlow <- is.null(rl.xlow)
        mstar <- m
        m <- m * npy
        rl.mle <- gpdret(z, mstar)
        q <- gpdq2(z$mle, u, la, m)
        d1 <- rep(0, length(q))
        d2 <- (gpdq2(c(z$mle[1] + eps, z$mle[2]), u, la, m) -
            q)/eps
        d3 <- (gpdq2(c(z$mle[1], z$mle[2]), u, la, m) - q)/eps
        d <- cbind(d1, d2, d3)
        mat <- z$cov
        mat <- matrix(c((la * (1 - la))/z$n, 0, 0, 0, mat[1,
            1], mat[1, 2], 0, mat[2, 1], mat[2, 2]), nc = 3)
        vv <- apply(d, 1, q.form, m = mat)
        if (vv == 0){
          l <- 0
          l
        }
        if (est.rl.xlow)
            rl.xlow <- rl.mle - 1.5 * qnorm((1 - conf)/2, lower.tail = FALSE) *
                sqrt(vv)
        if (est.rl.xup)
            rl.xup <- rl.mle + 1.5 * qnorm((1 - conf)/2, lower.tail = FALSE) *
                sqrt(vv)
        x <- seq(rl.xlow, rl.xup, length = nint)
        sol <- z$mle[2]
        gpd.plik <- function(a) {
            if (m != Inf && m != 1/la ){
              sc <- (a * (xp - u))/((m * la)^a - 1)
            } else if (m == 1/la) {
              sc <- 0
            }
            else {
              sc <- (u - xp)/a
            }
            if (abs(a) < 10^(-4))
              l <- length(xdat) * log(sc) + sum(xdat - u)/sc
            else {
              y <- (xdat - u)/sc
              y <- 1 + a * y
              if (any(y <= 0) || sc <= 0)
                l <- 10^6
              else l <- length(xdat) * log(sc) + sum(log(y)) *
                (1/a + 1)
            }
            l
        }
        for (i in 1:nint) {
            xp <- x[i]
            opt <- optim(sol, gpd.plik, method = "BFGS")
            sol <- opt$par
            v[i] <- opt$value
        }
        lmts$upcross.level <- ma - 0.5 * qchisq(conf, 1)
        lmts$rl$mle <- rl.mle
        sfun <- splinefun(x, -v)
        lmts$rl$sfun <- sfun
        x1 <- order(sfun(c(rl.mle, rl.xup)))
        if (x1[1] == 1)
            lmts$rl$up <- bisearch(rl.mle, rl.xup, f = sfun,
                upcross.level = ma - 0.5 * qchisq(conf, 1))$x
        else lmts$rl$up <- bisearch(rl.xup, rl.mle, f = sfun,
            upcross.level = ma - 0.5 * qchisq(conf, 1))$x
        x2 <- order(sfun(c(rl.mle, rl.xlow)))
        if (x2[1] == 1)
            lmts$rl$dn <- bisearch(rl.mle, rl.xlow, f = sfun,
                upcross.level = ma - 0.5 * qchisq(conf, 1))$x
        else lmts$rl$dn <- bisearch(rl.xlow, rl.mle, f = sfun,
            upcross.level = ma - 0.5 * qchisq(conf, 1))$x
        if (make.plot) {    #Plot Profile Likelihood Curves
            plot(x, -v, type = "l", xlab = "Return Level", ylab = "Profile Log-likelihood",sub=retp_ID)
            abline(h = ma, col = "green")
            abline(h = ma - 0.5 * qchisq(conf, 1), col = "green")
            abline(v = c(lmts$rl$dn, lmts$rl$up), lty = 2)
            readline("Press a key to continue ...")
        }
    }
    if (!rl.only) {    #Calc. CI for shape parameter
        est.xi.xup <- is.null(xi.xup)
        est.xi.xlow <- is.null(xi.xlow)
        xdat <- z$data
        xi.mle <- z$mle[2]
        if (est.xi.xlow)
            xi.xlow <- xi.mle - 1.5 * qnorm((1 - conf)/2, lower.tail = FALSE) *
                z$se[2]
        if (est.xi.xup)
            xi.xup <- xi.mle + 1.5 * qnorm((1 - conf)/2, lower.tail = FALSE) *
                z$se[2]
        u <- z$threshold
        v <- numeric(nint)
        x <- seq(xi.xup + qnorm((1 - conf)/2, lower.tail = FALSE) *
            z$se[2], xi.xlow - qnorm((1 - conf)/2, lower.tail = FALSE) *
            z$se[2], length = nint)
        sol <- z$mle[1]
        gpd.plikxi <- function(a) {
            if (abs(xi) < 10^(-4))
                l <- length(xdat) * log(a) + sum(xdat - u)/a
            else {
                y <- (xdat - u)/a
                y <- 1 + xi * y
                if (any(y <= 0) || a <= 0)
                  l <- 10^6
                else l <- length(xdat) * log(a) + sum(log(y)) *
                  (1/xi + 1)
            }
            l
        }
        for (i in 1:nint) {
            xi <- x[i]
            opt <- optim(sol, gpd.plikxi, method = "BFGS")
            sol <- opt$par
            v[i] <- opt$value
        }
        sfun <- splinefun(x, -v)
        lmts$xi$sfun <- sfun
        x1 <- order(sfun(c(xi.mle, xi.xup)))
        x2 <- order(sfun(c(xi.mle, xi.xlow)))
        if (x1[1] == 1)
            lmts$xi$up <- bisearch(xi.mle, xi.xup, f = sfun,
                upcross.level = ma - 0.5 * qchisq(conf, 1))$x
        else lmts$xi$up <- bisearch(xi.xup, xi.mle, f = sfun,
            upcross.level = ma - 0.5 * qchisq(conf, 1))$x
        if (x2[1] == 1)
            lmts$xi$dn <- bisearch(xi.mle, xi.xlow, f = sfun,
                upcross.level = ma - 0.5 * qchisq(conf, 1))$x
        else lmts$xi$dn <- bisearch(xi.xlow, xi.mle, f = sfun,
            upcross.level = ma - 0.5 * qchisq(conf, 1))$x
        if (make.plot) {
            plot(x, -v, type = "l", xlab = "Shape Parameter",
                ylab = "Profile Log-likelihood")
            abline(h = ma, lty = 1, col = "green")
            abline(h = ma - 0.5 * qchisq(conf, 1), col = "green")
            abline(v = c(lmts$xi$dn, lmts$xi$up), lty = 2)
        }
    }
    class(lmts) <- "gpd.parameterCI.obj"
    invisible(lmts)
}
