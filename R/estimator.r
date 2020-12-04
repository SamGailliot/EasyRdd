estimator <- function (formula, data, subset = NULL, cutpoint = NULL, bw = NULL,
          kernel = "triangular", se.type = "HC1", cluster = NULL, verbose = FALSE,
          model = FALSE, frame = FALSE)
{
  call <- match.call()
  if (missing(data))
    data <- environment(formula)
  formula <- as.Formula(formula)
  X <- model.frame(formula, rhs = 1, lhs = 0, data = data,
                   na.action = na.pass)[[1]]
  Y <- model.frame(formula, rhs = 0, lhs = NULL, data = data,
                   na.action = na.pass)[[1]]
  if (!is.null(subset)) {
    X <- X[subset]
    Y <- Y[subset]
    if (!is.null(cluster))
      cluster <- cluster[subset]
  }
  if (!is.null(cluster)) {
    cluster <- as.character(cluster)
    robust.se <- function(model, cluster) {
      M <- length(unique(cluster))
      N <- length(cluster)
      K <- model$rank
      dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
      uj <- apply(estfun(model), 2, function(x) tapply(x,
                                                       cluster, sum))
      rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
      rcse.se <- coeftest(model, rcse.cov)
      return(rcse.se[2, 2])
    }
  }
  na.ok <- complete.cases(X) & complete.cases(Y)
  if (length(all.vars(formula(formula, rhs = 1, lhs = F))) >
      1) {
    type <- "fuzzy"
    Z <- model.frame(formula, rhs = 1, lhs = 0, data = data,
                     na.action = na.pass)[[2]]
    if (!is.null(subset))
      Z <- Z[subset]
    na.ok <- na.ok & complete.cases(Z)
    if (length(all.vars(formula(formula, rhs = 1, lhs = F))) >
        2)
      stop("Invalid formula. Read ?RDestimate for proper syntax")
  }
  else {
    type = "sharp"
  }
  covs <- NULL
  if (length(formula)[2] > 1) {
    covs <- model.frame(formula, rhs = 2, lhs = 0, data = data,
                        na.action = na.pass)
    if (!is.null(subset))
      covs <- subset(covs, subset)
    na.ok <- na.ok & complete.cases(covs)
    covs <- subset(covs, na.ok)
  }
  X <- X[na.ok]
  Y <- Y[na.ok]
  if (type == "fuzzy")
    Z <- as.double(Z[na.ok])
  if (is.null(cutpoint)) {
    cutpoint <- 0
    if (verbose)
      cat("No cutpoint provided. Using default cutpoint of zero.\n")
  }
  if (frame) {
    if (type == "sharp") {
      if (!is.null(covs))
        dat.out <- data.frame(X, Y, covs)
      else dat.out <- data.frame(X, Y)
    }
    else {
      if (!is.null(covs))
        dat.out <- data.frame(X, Y, Z, covs)
      else dat.out <- data.frame(X, Y, Z)
    }
  }
  if (is.null(bw)) {
    bw <- IKbandwidth(X = X, Y = Y, cutpoint = cutpoint,
                      kernel = kernel, verbose = verbose)
    bws <- c(bw, 0.5 * bw, 2 * bw)
    names(bws) <- c("LATE", "Half-BW", "Double-BW")
  }
  else if (length(bw) == 1) {
    bws <- c(bw, 0.5 * bw, 2 * bw)
    names(bws) <- c("LATE", "Half-BW", "Double-BW")
  }
  else {
    bws <- bw
  }
  o <- list()
  o$type <- type
  o$call <- call
  o$est <- vector(length = length(bws), mode = "numeric")
  names(o$est) <- names(bws)
  o$bw <- as.vector(bws)
  o$se <- vector(mode = "numeric")
  o$z <- vector(mode = "numeric")
  o$p <- vector(mode = "numeric")
  o$obs <- vector(mode = "numeric")
  o$ci <- matrix(NA, nrow = length(bws), ncol = 2)
  o$model <- list()
  if (type == "fuzzy") {
    o$model$firststage <- list()
    o$model$iv <- list()
  }
  o$frame <- list()
  o$na.action <- which(na.ok == FALSE)
  class(o) <- "RD"
  X <- X - cutpoint
  Xl <- (X < 0) * X
  Xr <- (X >= 0) * X
  Tr <- as.integer(X >= 0)
  for (bw in bws) {
    ibw <- which(bw == bws)
    sub <- X >= (-bw) & X <= (+bw)
    if (kernel == "gaussian")
      sub <- TRUE
    w <- kernelwts(X, 0, bw, kernel = kernel)
    o$obs[ibw] <- sum(w > 0)
    if (type == "sharp") {
      if (verbose) {
        cat("Running Sharp RD\n")
        cat("Running variable:", all.vars(formula(formula,
                                                  rhs = 1, lhs = F))[1], "\n")
        cat("Outcome variable:", all.vars(formula(formula,
                                                  rhs = F, lhs = 1))[1], "\n")
        if (!is.null(covs))
          cat("Covariates:", paste(names(covs), collapse = ", "),
              "\n")
      }
      if (!is.null(covs)) {
        data <- data.frame(Y, Tr, Xl, Xr, covs, w)
        form <- as.formula(paste("Y~Tr+Xl+Xr+", paste(names(covs),
                                                      collapse = "+", sep = ""), sep = ""))
      }
      else {
        data <- data.frame(Y, Tr, Xl, Xr, w)
        form <- as.formula(Y ~ Tr + Xl + Xr)
      }
      mod <- lm(form, weights = w, data = subset(data,
                                                 w > 0))
      if (verbose == TRUE) {
        cat("Model:\n")
        print(summary(mod))
      }
      o$est[ibw] <- coef(mod)["Tr"]
      if (is.null(cluster)) {
        o$se[ibw] <- coeftest(mod, vcovHC(mod, type = se.type))[2,
                                                                2]
      }
      else {
        o$se[ibw] <- robust.se(mod, cluster[na.ok][w >
                                                     0])
      }
      o$z[ibw] <- o$est[ibw]/o$se[ibw]
      o$p[ibw] <- 2 * pnorm(abs(o$z[ibw]), lower.tail = F)
      o$ci[ibw, ] <- c(o$est[ibw] - qnorm(0.975) * o$se[ibw],
                       o$est[ibw] + qnorm(0.975) * o$se[ibw])
      if (model)
        o$model[[ibw]] = mod
      if (frame)
        o$frame[[ibw]] = dat.out
    }
    else {
      if (verbose) {
        cat("Running Fuzzy RD\n")
        cat("Running variable:", all.vars(formula(formula,
                                                  rhs = 1, lhs = F))[1], "\n")
        cat("Outcome variable:", all.vars(formula(formula,
                                                  rhs = F, lhs = 1))[1], "\n")
        cat("Treatment variable:", all.vars(formula(formula,
                                                    rhs = 1, lhs = F))[2], "\n")
        if (!is.null(covs))
          cat("Covariates:", paste(names(covs), collapse = ", "),
              "\n")
      }
      if (!is.null(covs)) {
        data <- data.frame(Y, Tr, Xl, Xr, Z, covs, w)
        form <- as.Formula(paste("Y~Z+Xl+Xr+", paste(names(covs),
                                                     collapse = "+"), "|Tr+Xl+Xr+", paste(names(covs),
                                                                                          collapse = "+"), sep = ""))
        form1 <- as.Formula(paste("Z~Tr+Xl+Xr+", paste(names(covs),
                                                       collapse = "+", sep = "")))
      }
      else {
        data <- data.frame(Y, Tr, Xl, Xr, Z, w)
        form <- as.Formula(Y ~ Z + Xl + Xr | Tr + Xl +
                             Xr)
        form1 <- as.formula(Z ~ Tr + Xl + Xr)
      }
      mod1 <- lm(form1, weights = w, data = subset(data,
                                                   w > 0))
      mod <- ivreg(form, weights = w, data = subset(data,
                                                    w > 0))
      if (verbose == TRUE) {
        cat("First stage:\n")
        print(summary(mod1))
        cat("IV-RD:\n")
        print(summary(mod))
      }
      o$est[ibw] <- coef(mod)["Z"]
      if (is.null(cluster)) {
        o$se[ibw] <- coeftest(mod, vcovHC(mod, type = se.type))[2,
                                                                2]
      }
      else {
        o$se[ibw] <- robust.se(mod, cluster[na.ok][w >
                                                     0])
      }
      o$z[ibw] <- o$est[ibw]/o$se[ibw]
      o$p[ibw] <- 2 * pnorm(abs(o$z[ibw]), lower.tail = F)
      o$ci[ibw, ] <- c(o$est[ibw] - qnorm(0.975) * o$se[ibw],
                       o$est[ibw] + qnorm(0.975) * o$se[ibw])
      if (model) {
        o$model$firststage[[ibw]] <- mod1
        o$model$iv[[ibw]] = mod
      }
      if (frame)
        o$frame = dat.out
    }
  }
  return(o)
}

estimator()
