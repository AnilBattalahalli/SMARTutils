#' @importFrom stats rbinom
#' @importFrom stats var
#' @importFrom stats terms
#' @importFrom stats pnorm
#' @importFrom stats model.frame
#' @importFrom dplyr %>%
#' @importFrom utils combn
NULL

Ii <- function(A1, R, A2, a1, a2) {
  ifelse(A1 == a1, 1, 0) * (R + (ifelse(A2 == a2, 1, 0) * (1 - R)))
}

Wi <- function(R) {
  2 * (R + (2 * (1 - R)))
}

IWi <- function(A1, R, A2, a1, a2) {
  2 * ifelse(A1 == a1, 1, 0) * (R + (2 * ifelse(A2 == a2, 1, 0) * (1 - R)))
}

padding_str <- function(x, n) {
  x <- as.character(x)
  charlen <- nchar(x)
  trails <- n - charlen

  return(paste(c(x, rep(" ", trails)), collapse = ""))
}

pval <- function(x) {
  2 * pnorm(q = abs(x), lower.tail = FALSE)
}

perm_mul_sum <- function(l) {
  l <- rbind(t(combn(l, 2)), t(combn(l, 2)))
  l <- data.frame(l)
  return(sum(l$X1 * l$X2))
}

D <- function(data_subset, a1, a2, formula) {
  data_subset$a1 <- rep(a1, nrow(data_subset))
  data_subset$a2 <- rep(a2, nrow(data_subset))
  df <- model.frame(formula, data_subset)
  df <- within(df, rm("Y"))
  I <- rep(1, nrow(df))
  df <- as.matrix(cbind(I, df)) %>% unname()
  return(df)
}

mmm <- function(designmat, betas) {
  mu <- designmat %*% matrix(betas)
  return(mu)
}

get_residuals <- function(betas, data, formula) {
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <-
    list(
      "1,1" = c(1, 1),
      "1,-1" = c(1, -1),
      "-1,1" = c(-1, 1),
      "-1,-1" = c(-1, -1)
    )
  N <- data$i %>%
    unique() %>%
    length()
  epsilon <- list()
  for (d in dtrs) {
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    temp <- list()
    for (ci in unique(data$i)) {
      data_sub <- data[data$i == ci, ]
      temp[ci] <-
        list(as.numeric(data_sub$Y - mmm(D(
          data_sub, a1, a2, formula
        ), betas)))
    }
    epsilon[d] <- list(temp)
  }
  return(epsilon)
}

get_varrho_star <- function(data, epsilon, weights) {
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <-
    list(
      "1,1" = c(1, 1),
      "1,-1" = c(1, -1),
      "-1,1" = c(-1, 1),
      "-1,-1" = c(-1, -1)
    )
  var_star <- list()
  rho_star <- list()
  N <- data$i %>%
    unique() %>%
    length()

  for (d in dtrs) {
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    s <- 0
    den <- 0
    for (ci in unique(data$i)) {
      data_sub <- data[data$i == ci, ]
      ni <- nrow(data_sub)
      A1 <- data_sub$A1[1]
      R <- data_sub$R[1]
      A2 <- data_sub$A2[1]
      eps_sq <- epsilon[[d]][[ci]]^2
      if (is.null(weights) == F) {
        W <- weights[weights$i == ci, ]$W[1] * Ii(A1, R, A2, a1, a2)
      } else {
        W <- IWi(A1, R, A2, a1, a2)
      }
      if (W != 0) {
        s <- s + (W * sum(eps_sq))
        den <- den + (W * ni)
      }
    }
    var_star[d] <- s / den
  }

  for (d in dtrs) {
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    s <- 0
    den <- 0
    for (ci in unique(data$i)) {
      data_sub <- data[data$i == ci, ]
      ni <- nrow(data_sub)
      if (ni > 1) {
        A1 <- data_sub$A1[1]
        R <- data_sub$R[1]
        A2 <- data_sub$A2[1]
        eps_mul_sum <- perm_mul_sum(epsilon[[d]][[ci]])
        if (is.null(weights) == F) {
          W <- weights[weights$i == ci, ]$W[1] * Ii(A1, R, A2, a1, a2)
        } else {
          W <- IWi(A1, R, A2, a1, a2)
        }
        if (W != 0) {
          s <- s + (W * eps_mul_sum)
          den <- den + (W * ni * (ni - 1))
        }
      }
    }
    den <- den * var_star[[d]]
    rho_star[d] <- s / den
  }
  rhovar_stars <- list(rho_star = rho_star, var_star = var_star)
  return(rhovar_stars)
}

get_beta <- function(var, rho, data, formula, weights) {
  params <- (formula %>% terms() %>% labels() %>% length()) + 1
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <-
    list(
      "1,1" = c(1, 1),
      "1,-1" = c(1, -1),
      "-1,1" = c(-1, 1),
      "-1,-1" = c(-1, -1)
    )
  f <- matrix(0, params, params)
  N <- data$i %>%
    unique() %>%
    length()
  for (ci in unique(data$i)) {
    data_sub <- data[data$i == ci, ]
    A1 <- data_sub$A1[1]
    R <- data_sub$R[1]
    A2 <- data_sub$A2[1]
    ni <- nrow(data_sub)
    for (d in dtrs) {
      a1 <- dtrenc[[d]][1]
      a2 <- dtrenc[[d]][2]
      if (is.null(weights) == F) {
        W <- weights[weights$i == ci, ]$W[1] * Ii(A1, R, A2, a1, a2)
      } else {
        W <- IWi(A1, R, A2, a1, a2)
      }
      if (W != 0) {
        V <- var[[d]] * (diag(ni) + (rho[[d]] * (matrix(1, ni, ni) - diag(ni))))
        f <- f + (W * t(D(data_sub, a1, a2, formula)) %*% solve(V) %*% D(data_sub, a1, a2, formula))
      }
    }
  }

  s <- matrix(0, params, 1)
  for (ci in unique(data$i)) {
    data_sub <- data[data$i == ci, ]
    A1 <- data_sub$A1[1]
    R <- data_sub$R[1]
    A2 <- data_sub$A2[1]
    ni <- nrow(data_sub)
    Y <- matrix(data_sub$Y)
    for (d in dtrs) {
      a1 <- dtrenc[[d]][1]
      a2 <- dtrenc[[d]][2]
      if (is.null(weights) == F) {
        W <- weights[weights$i == ci, ]$W[1] * Ii(A1, R, A2, a1, a2)
      } else {
        W <- IWi(A1, R, A2, a1, a2)
      }
      if (W != 0) {
        V <- var[[d]] * (diag(ni) + (rho[[d]] * (matrix(1, ni, ni) - diag(ni))))
        s <-
          s + (IWi(A1, R, A2, a1, a2) * t(D(data_sub, a1, a2, formula)) %*% solve(V) %*% Y)
      }
    }
  }
  beta_hat <- solve(f) %*% s
  return(beta_hat)
}

clustered.estimate_var_betas <-
  function(data, var, rho, betas, formula, weights) {
    params <- (formula %>% terms() %>% labels() %>% length()) + 1
    dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
    dtrenc <-
      list(
        "1,1" = c(1, 1),
        "1,-1" = c(1, -1),
        "-1,1" = c(-1, 1),
        "-1,-1" = c(-1, -1)
      )
    J_hat <- matrix(0, params, params)
    N <- data$i %>%
      unique() %>%
      length()
    for (d in dtrs) {
      a1 <- dtrenc[[d]][1]
      a2 <- dtrenc[[d]][2]
      temp <- matrix(0, params, params)
      for (ci in unique(data$i)) {
        data_sub <- data[data$i == ci, ]
        ni <- nrow(data_sub)
        A1 <- data_sub$A1[1]
        R <- data_sub$R[1]
        A2 <- data_sub$A2[1]
        if (is.null(weights) == F) {
          W <- weights[weights$i == ci, ]$W[1] * Ii(A1, R, A2, a1, a2)
        } else {
          W <- IWi(A1, R, A2, a1, a2)
        }
        if (W != 0) {
          V <- var[[d]] * (diag(ni) + (rho[[d]] * (matrix(1, ni, ni) - diag(ni))))
          temp <-
            temp + (W * t(D(data_sub, a1, a2, formula)) %*% solve(V) %*% D(data_sub, a1, a2, formula))
        }
      }
      J_hat <- J_hat + temp
    }
    J_hat <- J_hat / N

    A_hat <- matrix(0, params, params)
    for (ci in unique(data$i)) {
      for (d in dtrs) {
        a1 <- dtrenc[[d]][1]
        a2 <- dtrenc[[d]][2]

        data_sub <- data[data$i == ci, ]
        A1 <- data_sub$A1[1]
        R <- data_sub$R[1]
        A2 <- data_sub$A2[1]
        ni <- nrow(data_sub)
        Y <- matrix(data_sub$Y)
        diff <- Y - mmm(D(data_sub, a1, a2, formula), betas)
        V <- var[[d]] * (diag(ni) + (rho[[d]] * (matrix(1, ni, ni) - diag(ni))))
        if (is.null(weights) == F) {
          W <- weights[weights$i == ci, ]$W[1] * Ii(A1, R, A2, a1, a2)
        } else {
          W <- IWi(A1, R, A2, a1, a2)
        }
        if (W != 0) {
          Ui <- W * t(D(data_sub, a1, a2, formula)) %*% solve(V) %*% diff
          A_hat <- A_hat + (Ui %*% t(Ui))
        }
      }
    }
    A_hat <- A_hat / N

    var_hat_beta_hat <-
      (1 / N) * (solve(J_hat) %*% A_hat %*% solve(J_hat))
    return(var_hat_beta_hat)
  }

clustered.prepare_data <-
  function(data,
           i = "i",
           A1 = "A1",
           R = "R",
           A2 = "A2",
           Y = "Y") {
    names(data)[names(data) == i] <- "i"
    names(data)[names(data) == A1] <- "A1"
    names(data)[names(data) == R] <- "R"
    names(data)[names(data) == A2] <- "A2"
    names(data)[names(data) == Y] <- "Y"
    return(data)
  }

#' @title Function to estimate the marginal mean model parameters
#' @param formula formula for the estimation
#' @param data a data.frame object
#' @param i cluster ID column name
#' @param A1 A1 column name
#' @param R R column name
#' @param A2 A2 column name
#' @param Y outcome variable column name
#' @param covstr structure of the working covariance matrix 'IND' or 'EXCH' or 'EXCH|AI'
#' @param verbose T/F for the model fit reports
#' @param weights a data frame containing cluster level weights with columns `i` and `W` for cluster ID and weights respectively
#' @return estimates object containing all the estimated parameters.
#' @export cSMART.mm
cSMART.mm <-
  function(formula,
           data,
           i = "i",
           A1 = "A1",
           R = "R",
           A2 = "A2",
           Y = "Y",
           covstr = "EXCH",
           verbose = T,
           weights = NULL) {
    data <- clustered.prepare_data(data, i, A1, R, A2, Y)
    data$A2[is.na(data$A2)] <- 0

    if (covstr == "IND") {
      var <- list(
        "1,1" = 1,
        "1,-1" = 1,
        "-1,1" = 1,
        "-1,-1" = 1
      )
      rho <- list(
        "1,1" = 0,
        "1,-1" = 0,
        "-1,1" = 0,
        "-1,-1" = 0
      )

      betas <- get_beta(var, rho, data, formula, weights)
      epsilon <- get_residuals(betas, data, formula)
      varrho_stars <- get_varrho_star(data, epsilon, weights)
      var <- varrho_stars$var_star

      cov_hat <-
        clustered.estimate_var_betas(data, var, rho, betas, formula, weights)
    }

    if (covstr == "EXCH") {
      var <- list(
        "1,1" = 1,
        "1,-1" = 1,
        "-1,1" = 1,
        "-1,-1" = 1
      )
      rho <- list(
        "1,1" = 0,
        "1,-1" = 0,
        "-1,1" = 0,
        "-1,-1" = 0
      )

      betas <- get_beta(var, rho, data, formula, weights)
      epsilon <- get_residuals(betas, data, formula)
      varrho_stars <- get_varrho_star(data, epsilon, weights)
      var <- varrho_stars$var_star
      rho <- varrho_stars$rho_star
      var_mean <- mean(unlist(var))
      rho_mean <- mean(unlist(rho))
      var$`1,1` <- var_mean
      var$`1,-1` <- var_mean
      var$`-1,1` <- var_mean
      var$`-1,-1` <- var_mean
      rho$`1,1` <- rho_mean
      rho$`1,-1` <- rho_mean
      rho$`-1,1` <- rho_mean
      rho$`-1,-1` <- rho_mean

      betas <- get_beta(var, rho, data, formula, weights)
      epsilon <- get_residuals(betas, data, formula)
      varrho_stars <- get_varrho_star(data, epsilon, weights)
      var <- varrho_stars$var_star
      rho <- varrho_stars$rho_star
      var_mean <- mean(unlist(var))
      rho_mean <- mean(unlist(rho))
      var$`1,1` <- var_mean
      var$`1,-1` <- var_mean
      var$`-1,1` <- var_mean
      var$`-1,-1` <- var_mean
      rho$`1,1` <- rho_mean
      rho$`1,-1` <- rho_mean
      rho$`-1,1` <- rho_mean
      rho$`-1,-1` <- rho_mean

      betas <- get_beta(var, rho, data, formula, weights)
      epsilon <- get_residuals(betas, data, formula)
      varrho_stars <- get_varrho_star(data, epsilon, weights)
      var <- varrho_stars$var_star
      rho <- varrho_stars$rho_star
      var_mean <- mean(unlist(var))
      rho_mean <- mean(unlist(rho))
      var$`1,1` <- var_mean
      var$`1,-1` <- var_mean
      var$`-1,1` <- var_mean
      var$`-1,-1` <- var_mean
      rho$`1,1` <- rho_mean
      rho$`1,-1` <- rho_mean
      rho$`-1,1` <- rho_mean
      rho$`-1,-1` <- rho_mean


      cov_hat <-
        clustered.estimate_var_betas(data, var, rho, betas, formula, weights)
    }

    if (covstr == "EXCH|AI") {
      var <- list(
        "1,1" = 1,
        "1,-1" = 1,
        "-1,1" = 1,
        "-1,-1" = 1
      )
      rho <- list(
        "1,1" = 0,
        "1,-1" = 0,
        "-1,1" = 0,
        "-1,-1" = 0
      )

      betas <- get_beta(var, rho, data, formula, weights)
      epsilon <- get_residuals(betas, data, formula)
      varrho_stars <- get_varrho_star(data, epsilon, weights)
      var <- varrho_stars$var_star
      rho <- varrho_stars$rho_star

      betas <- get_beta(var, rho, data, formula, weights)
      epsilon <- get_residuals(betas, data, formula)
      varrho_stars <- get_varrho_star(data, epsilon, weights)
      var <- varrho_stars$var_star
      rho <- varrho_stars$rho_star

      betas <- get_beta(var, rho, data, formula, weights)
      epsilon <- get_residuals(betas, data, formula)
      varrho_stars <- get_varrho_star(data, epsilon, weights)
      var <- varrho_stars$var_star
      rho <- varrho_stars$rho_star


      cov_hat <-
        clustered.estimate_var_betas(data, var, rho, betas, formula, weights)
    }

    report <- list(
      beta_hat = betas,
      var_hat = var,
      rho_hat = rho,
      cov_hat_beta_hat = cov_hat,
      formula = formula
    )
    formula <- report$formula %>% terms()
    terms <- c("(Intercept)", attr(terms(formula), "term.labels"))
    rownames(report$cov_hat_beta_hat) <- terms
    colnames(report$cov_hat_beta_hat) <- terms

    if (verbose) {
      n <- lapply(terms, nchar) %>%
        unlist() %>%
        max()
      terms_padded <- lapply(terms, padding_str, n = n) %>% unlist()

      cat(
        padding_str("Parameter", n),
        "\t",
        "Estimate",
        "\t",
        "Std.Err",
        "\t",
        "Z Score",
        "\t",
        "Pr(>|z|)",
        "\n"
      )
      cat(
        paste(rep("_", n), collapse = ""),
        "\t",
        "--------",
        "\t",
        "-------",
        "\t",
        "-------",
        "\t",
        "-------",
        "\n"
      )
      for (i in 1:length(terms)) {
        cat(
          terms_padded[i],
          "\t",
          sprintf("%.5f", report$beta_hat[i]),
          "\t",
          sprintf("%.5f", sqrt(report$cov_hat_beta_hat[i, i])),
          "\t",
          report$beta_hat[i] / sqrt(report$cov_hat_beta_hat[i, i]),
          "\t",
          pval(report$beta_hat[i] / sqrt(report$cov_hat_beta_hat[i, i])),
          "\n"
        )
      }
      cat(
        paste(rep("_", n), collapse = ""),
        "\t",
        "--------",
        "\t",
        "-------",
        "\t",
        "-------",
        "\t",
        "-------",
        "\n"
      )

      cat("\n")

      cat(paste("Marginal Mean Model: ", format(formula)), "\n\n")

      if (covstr == "IND") {
        cat(
          "Working covariance structure: 'IND' ",
          "(Independent-Homogeneous covariance structure)",
          "\n"
        )
        cat("Variance", "\t", mean(unlist(report$var_hat)))
        cat("\n")
      }

      if (covstr == "EXCH") {
        cat(
          "Working covariance structure: 'EXCH' ",
          "(Homogeneous-Exchangeable covariance structure)"
        )
        cat("\n")
        cat("Variance", "\t", mean(unlist(report$var_hat)))
        cat("\n")
        cat("Correlation", "\t", mean(unlist(report$rho_hat)))
        cat("\n")
      }

      if (covstr == "EXCH|AI") {
        cat(
          "Working covariance structure: 'EXCH|AI' ",
          "(Exchangeable by adaptive intervention covariance structure)",
          "\n\n"
        )
        cat("  AI ", "\t", "Variance", "\n")
        cat("-----", "\t", "--------", "\n")
        cat(" 1, 1", "\t", report$var_hat$`1,1`, "\n")
        cat(" 1,-1", "\t", report$var_hat$`1,-1`, "\n")
        cat("-1, 1", "\t", report$var_hat$`-1,1`, "\n")
        cat("-1,-1", "\t", report$var_hat$`-1,-1`, "\n")
        cat("\n")
        cat("  AI ", "\t", "Correlation", "\n")
        cat("-----", "\t", "--------", "\n")
        cat(" 1, 1", "\t", report$rho_hat$`1,1`, "\n")
        cat(" 1,-1", "\t", report$rho_hat$`1,-1`, "\n")
        cat("-1, 1", "\t", report$rho_hat$`-1,1`, "\n")
        cat("-1,-1", "\t", report$rho_hat$`-1,-1`, "\n")
      }
      cat("\n")
      cat("Variance-Covariance matrix of the estimates")
      cat("\n")
      print(report$cov_hat_beta_hat)
    }
    return(report)
  }

nocovariates.true_effect_size <- function(recipe, d1, d2) {
  treat_mu <- recipe$treat_mu
  treat_var <- recipe$treat_var
  dtrnum <- list(
    "1,1" = 1,
    "1,-1" = 2,
    "-1,1" = 3,
    "-1,-1" = 4
  )
  numerator <- treat_mu[dtrnum[[d1]]] - treat_mu[dtrnum[[d2]]]
  denominator <-
    (0.5 * treat_var[dtrnum[[d1]]]) + (0.5 * treat_var[dtrnum[[d2]]])
  denominator <- sqrt(denominator)
  return(numerator / denominator)
}

nocovariates.treat <- function(recipe) {
  encoder <-
    list(
      "1,1,." = 1,
      "1,0,1" = 2,
      "1,0,-1" = 3,
      "-1,1,." = 4,
      "-1,0,1" = 5,
      "-1,0,-1" = 6
    )
  N <- recipe$N
  m_vec <- recipe$m
  cell_mu <- recipe$cell_mu
  cell_var <- recipe$cell_var
  cell_cor <- recipe$cell_cor
  cell_cov <- recipe$cell_cov
  p_1 <- recipe$p_1
  p_2 <- recipe$p_2
  p_A1 <- recipe$p_A1
  p_A2 <- recipe$p_A2

  data <- matrix(nrow = 0, ncol = 6)

  for (clusID in 1:N) {
    m <- m_vec[clusID]
    A1 <- (2 * rbinom(1, 1, p_A1)) - 1

    if (A1 == 1) {
      R <- rbinom(1, 1, p_1)
      if (R == 1) {
        treat_code <- "1,1,."
        A2 <- NA
      } else if (R == 0) {
        A2 <- (2 * rbinom(1, 1, p_A2)) - 1
        treat_code <- sprintf("1,0,%d", A2)
      }
    } else if (A1 == -1) {
      R <- rbinom(1, 1, p_2)
      if (R == 1) {
        treat_code <- "-1,1,."
        A2 <- NA
      } else if (R == 0) {
        A2 <- (2 * rbinom(1, 1, p_A2)) - 1
        treat_code <- sprintf("-1,0,%d", A2)
      }
    }

    cell_index <- as.numeric(encoder[treat_code])
    V <-
      cell_var[cell_index] * (diag(m) + (cell_cor[cell_index] * (matrix(1, m, m) -
        diag(m))))
    Y <-
      t(cell_mu[cell_index] + mvtnorm::rmvnorm(1, rep(0, m), V, method = "chol"))

    tempmat <-
      cbind(rep(clusID, m), seq(1, m), rep(A1, m), rep(R, m), rep(A2, m), Y)
    data <- rbind(data, tempmat)
  }
  full_data_frame <- as.data.frame(data)
  names(full_data_frame) <- c("i", "j", "A1", "R", "A2", "Y")
  return(full_data_frame)
}

#' @title Function to generate object from the conditional parameters
#' @description Returns the data.frame object from the conditional parameters
#' (mean, variance, correlation), number of clusters, number of participants in each cluster,
#' response probabilities of the different first level treatment options.
#' conditional mean, variance, and correlation must be numerical lists of length 6 with each
#' element corresponding to the cells (1, 1, .), (1, 0, 1), (1, 0, -1), (-1, 1, .), (-1, 0, 1), and (-1, 0, -1)
#' encoded in (A1, R, A2) form.
#' @param cell_mu a vector of conditional means.
#' @param cell_var a vector of conditional variances.
#' @param cell_cor a vector of conditional correlations.
#' @param N number of clusters to be generated.
#' @param m number of participants in a cluster.
#' This can be a numeric variable \code{m_1} in which case, it is assumed that the
#' number of participants in each cluster is \code{m_1}. This can also be a list of length N
#' where each element in the list corresponds to the number of participants in the cluster.
#' @param p_1 probability of response of the first level treatment A_1 = 1
#' @param p_2 probability of response of the first level treatment A_1 = -1
#' @param treatment.summary T/F to return the summary object
#' @return a data.frame object of the generated data
#' @export cSMART.dgen

cSMART.dgen <-
  function(cell_mu = c(10, 8, 2, 4, -2, 3),
           cell_var = c(100, 100, 100, 100, 100, 100),
           cell_cor = c(0.1, 0.2, 0.15, 0.12, 0.11, 0.07),
           N = 200,
           m = 5,
           p_1 = 0.3,
           p_2 = 0.2,
           treatment.summary = F) {
    if (length(m) != N) {
      m <- rep(m[1], N)
    }

    p_A1 <- 0.5
    p_A2 <- 0.5

    q_1 <- 1 - p_1
    q_2 <- 1 - p_2

    cell_cov <- cell_var * cell_cor

    treat_mu <-
      c(
        (cell_mu[1] * p_1) + (cell_mu[2] * q_1),
        (cell_mu[1] * p_1) + (cell_mu[3] * q_1),
        (cell_mu[4] * p_2) + (cell_mu[5] * q_2),
        (cell_mu[4] * p_2) + (cell_mu[6] * q_2)
      )

    treat_var <-
      c(
        (cell_var[1] * p_1 + cell_var[2] * q_1),
        (cell_var[1] * p_1 + cell_var[3] * q_1),
        (cell_var[4] * p_2 + cell_var[5] * q_2),
        (cell_var[4] * p_2 + cell_var[6] * q_2)
      )

    mean_diff_term <-
      c(
        ((cell_mu[1] - cell_mu[2])^2) * p_1 * q_1, ((cell_mu[1] - cell_mu[3])^
          2) * p_1 * q_1,
        ((cell_mu[4] - cell_mu[5])^2) * p_2 * q_2, ((cell_mu[4] -
          cell_mu[6])^2) * p_2 * q_2
      )

    treat_var <- treat_var + mean_diff_term

    treat_cov <-
      c(
        (cell_cov[1] * p_1) + (cell_cov[2] * q_1),
        (cell_cov[1] * p_1) + (cell_cov[3] * q_1),
        (cell_cov[4] * p_2) + (cell_cov[5] * q_2),
        (cell_cov[4] * p_2 + cell_cov[6] * q_2)
      )

    treat_cov <- treat_cov + mean_diff_term
    treat_cor <- treat_cov / treat_var
    A <-
      matrix(c(1, 1, 1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1), 4, 4)
    B <- matrix(treat_mu, 4, 1)
    betas <- solve(A, B)

    treat_summary <-
      data.frame(treat_mu, treat_var, treat_cor, betas)
    row.names(treat_summary) <-
      c("(1, 1)", "(1, -1)", "(-1, 1)", "(-1,-1)")

    recipe <-
      list(
        cell_mu = cell_mu,
        cell_var = cell_var,
        cell_cor = cell_cor,
        cell_cov = cell_cov,
        treat_mu = treat_mu,
        treat_var = treat_var,
        treat_cov = treat_cov,
        treat_cor = treat_cor,
        p_1 = p_1,
        p_2 = p_2,
        p_A1 = p_A1,
        p_A2 = p_A2,
        N = N,
        m = m,
        betas = betas,
        treat_summary = treat_summary
      )
    dtrs <- c(" 1, 1", " 1,-1", "-1, 1", "-1,-1")


    data <- nocovariates.treat(recipe)

    if (treatment.summary == T) {
      return(list(data = data, summary = treat_summary))
    } else {
      cat(
        " DTR",
        "\t",
        "Treatment Mean",
        "\t",
        "Treatment Variance",
        "\t",
        "Treatment Correlation",
        "\n"
      )
      cat(
        "______",
        "\t",
        "______________",
        "\t",
        "__________________",
        "\t",
        "_____________________",
        "\n"
      )
      for (i in 1:4) {
        cat(
          dtrs[i],
          "\t\t",
          treat_mu[i],
          "\t\t",
          " ",
          treat_var[i],
          "\t\t",
          " ",
          treat_cor[i],
          "\n"
        )
      }
      cat(
        "______",
        "\t",
        "______________",
        "\t",
        "__________________",
        "\t",
        "_____________________",
        "\n"
      )

      cat("\nTrue Beta:\n\n")
      for (i in 1:4) {
        cat(sprintf("Beta[%d]", i - 1), "-", betas[i], "\n")
      }

      return(data)
    }
  }