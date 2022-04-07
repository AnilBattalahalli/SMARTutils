#' @importFrom stats rbinom
#' @importFrom stats var
#' @importFrom utils combn
#' @importFrom pbmcapply pbmclapply
NULL

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

IWi <- function(A1, R, A2, a1, a2){
  2*ifelse(A1==a1,1,0)*(R+(2*ifelse(A2==a2,1,0)*(1-R)))
}

perm_mul_sum <- function(l){
  l <- rbind(t(combn(l, 2)), t(combn(l, 2)))
  l <- data.frame(l)
  return(sum(l$X1*l$X2))
}

D <- function(n, a1, a2){
  temp <- c(1, a1, a2, a1*a2)
  rep.row(temp, n)
}

mmm <- function(y, betas, a1, a2){
  l <- length(y)
  mu <- rep.row(betas[1] + (a1 * betas[2]) + (a2 * betas[3]) + (a1 * a2 * betas[4]), l)
  return(mu)
}

get_reiduals <- function(betas, data){
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <- list("1,1" = c(1,1), "1,-1"=c(1,-1), "-1,1"=c(-1,1), "-1,-1"=c(-1,-1))
  N <- max(data$i)
  epsilon <- list()
  for (d in dtrs){
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    temp <- list()
    for (ci in 1:N){
      data_sub <- data[data$i == ci, ]
      temp[ci] <- list(as.numeric(data_sub$Y - mmm(data_sub$Y, betas, a1, a2)))
    }
    epsilon[d] <- list(temp)
  }
  return(epsilon)
}

get_varrho_star <- function(data, epsilon){
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <- list("1,1" = c(1,1), "1,-1"=c(1,-1), "-1,1"=c(-1,1), "-1,-1"=c(-1,-1))
  var_star <- list()
  rho_star <- list()
  N <- max(data$i)

  for (d in dtrs){
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    s <- 0
    den <- 0
    for (ci in 1:N){
      data_sub <- data[data$i == ci, ]
      ni <- nrow(data_sub)
      A1 <- data_sub$A1[1]
      R <- data_sub$R[1]
      A2 <- data_sub$A2[1]
      eps_sq <- epsilon[[d]][[ci]]^2
      if (IWi(A1, R, A2, a1, a2) != 0){
        s <- s + (IWi(A1, R, A2, a1, a2) * sum(eps_sq))
        den <- den + (IWi(A1, R, A2, a1, a2) * ni)
      }
    }
    var_star[d] <- s/den
  }

  for (d in dtrs){
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    s <- 0
    den <- 0
    for (ci in 1:N){
      data_sub <- data[data$i == ci, ]
      ni <- nrow(data_sub)
      A1 <- data_sub$A1[1]
      R <- data_sub$R[1]
      A2 <- data_sub$A2[1]
      eps_mul_sum <- perm_mul_sum(epsilon[[d]][[ci]])
      if (IWi(A1, R, A2, a1, a2) != 0){
        s <- s + (IWi(A1, R, A2, a1, a2) * eps_mul_sum)
        den <- den + (IWi(A1, R, A2, a1, a2) * ni * (ni - 1))
      }
    }
    den <- den * var_star[[d]]
    rho_star[d] <- s/den
  }
  rhovar_stars <- list(rho_star=rho_star, var_star=var_star)
  return(rhovar_stars)
}

get_beta <- function(var, rho, data){
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <- list("1,1" = c(1,1), "1,-1"=c(1,-1), "-1,1"=c(-1,1), "-1,-1"=c(-1,-1))
  f <- matrix(0, 4, 4)
  N <- max(data$i)
  for (ci in 1:N){
    data_sub <- data[data$i == ci, ]
    A1 <- data_sub$A1[1]
    R <- data_sub$R[1]
    A2 <- data_sub$A2[1]
    ni <- nrow(data_sub)
    for (d in dtrs){
      a1 <- dtrenc[[d]][1]
      a2 <- dtrenc[[d]][2]
      if (IWi(A1, R, A2, a1, a2) != 0){
        V = var[[d]] * (diag(ni)+(rho[[d]]*(matrix(1,ni,ni)-diag(ni))))
        f = f + (IWi(A1, R, A2, a1, a2) * t(D(ni, a1, a2)) %*% solve(V) %*% D(ni, a1, a2))
      }
    }
  }

  s <- matrix(0, 4, 1)
  for (ci in 1:N){
    data_sub <- data[data$i == ci, ]
    A1 <- data_sub$A1[1]
    R <- data_sub$R[1]
    A2 <- data_sub$A2[1]
    ni <- nrow(data_sub)
    Y <- matrix(data_sub$Y)
    for (d in dtrs){
      a1 <- dtrenc[[d]][1]
      a2 <- dtrenc[[d]][2]
      if (IWi(A1, R, A2, a1, a2) != 0){
        V = var[[d]] * (diag(ni)+(rho[[d]]*(matrix(1,ni,ni)-diag(ni))))
        s <- s + (IWi(A1, R, A2, a1, a2) * t(D(ni, a1, a2)) %*% solve(V) %*% Y)
      }
    }
  }
  beta_hat <- solve(f) %*% s
  return(beta_hat)
}

nocovariates.estimate_var_betas <- function(data, var, rho, betas){
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <- list("1,1" = c(1,1), "1,-1"=c(1,-1), "-1,1"=c(-1,1), "-1,-1"=c(-1,-1))
  J_hat <- matrix(0, 4, 4)
  N <- max(data$i)
  for (d in dtrs){
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    temp <- matrix(0, 4, 4)
    for (ci in 1:N){
      data_sub <- data[data$i == ci, ]
      ni <- nrow(data_sub)
      A1 <- data_sub$A1[1]
      R <- data_sub$R[1]
      A2 <- data_sub$A2[1]
      if (IWi(A1, R, A2, a1, a2) != 0){
        V = var[[d]] * (diag(ni)+(rho[[d]]*(matrix(1,ni,ni)-diag(ni))))
        temp <- temp + (IWi(A1, R, A2, a1, a2) * t(D(ni, a1, a2)) %*% solve(V) %*% D(ni, a1, a2))
      }
    }
    J_hat <- J_hat+temp
  }
  J_hat <-  J_hat/N

  A_hat <- matrix(0, 4, 4)
  for (d in dtrs){
    temp <- matrix(0, 4, 4)
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    for (ci in 1:N){
      data_sub <- data[data$i == ci, ]
      A1 <- data_sub$A1[1]
      R <- data_sub$R[1]
      A2 <- data_sub$A2[1]
      ni <- nrow(data_sub)
      Y <- matrix(data_sub$Y)
      diff <- Y - mmm(Y, betas, a1, a2)
      V = var[[d]] * (diag(ni)+(rho[[d]]*(matrix(1,ni,ni)-diag(ni))))
      if (IWi(A1, R, A2, a1, a2) != 0){
        Ui <- IWi(A1, R, A2, a1, a2) * t(D(ni, a1, a2)) %*% solve(V) %*% diff
        temp <- temp + (Ui %*% t(Ui))
      }
    }
    A_hat <- A_hat + temp
  }
  A_hat <- A_hat/N

  var_hat_beta_hat <- (1/N) * (solve(J_hat) %*% A_hat %*% solve(J_hat))
  return(var_hat_beta_hat)
}

#' @title Function to get treatment recipe object from the conditional parameters
#' @description Returns the treatment recipe object from the conditional parameters
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
#' @return a treatment recipe object containing all the conditional and marginal parameters.
#' @export nocovairiates.from_conditional
nocovairiates.from_conditional <- function(cell_mu = c(10,8,2,4,-2,3), cell_var = c(100,100,100,100,100,100),
                                           cell_cor = c(0.1, 0.2, 0.15, 0.12, 0.11, 0.07), N=200, m=5, p_1=0.3, p_2=0.2){

  if (length(m) != N){
    m <- rep(m[1], N)
  }

  p_A1 <- 0.5
  p_A2 <- 0.5

  q_1 <-  1 - p_1
  q_2 <-  1 - p_2

  cell_cov <-  cell_var*cell_cor

  treat_mu <-  c((cell_mu[1]*p_1) + (cell_mu[2]*q_1), (cell_mu[1]*p_1) + (cell_mu[3]*q_1),
                 (cell_mu[4]*p_2) + (cell_mu[5]*q_2), (cell_mu[4]*p_2) + (cell_mu[6]*q_2))

  treat_var <-  c((cell_var[1]*p_1 + cell_var[2]*q_1), (cell_var[1]*p_1 + cell_var[3]*q_1),
                (cell_var[4]*p_2 + cell_var[5]*q_2), (cell_var[4]*p_2 + cell_var[6]*q_2))

  mean_diff_term <-  c(((cell_mu[1]-cell_mu[2])^2)*p_1*q_1, ((cell_mu[1]-cell_mu[3])^2)*p_1*q_1,
                     ((cell_mu[4]-cell_mu[5])^2)*p_2*q_2, ((cell_mu[4]-cell_mu[6])^2)*p_2*q_2)

  treat_var <-  treat_var + mean_diff_term

  treat_cov <-  c((cell_cov[1]*p_1) + (cell_cov[2]*q_1), (cell_cov[1]*p_1) + (cell_cov[3]*q_1),
                (cell_cov[4]*p_2) + (cell_cov[5]*q_2), (cell_cov[4]*p_2 + cell_cov[6]*q_2))

  treat_cov <-  treat_cov + mean_diff_term
  treat_cor <-  treat_cov/treat_var
  A <-  matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1), 4, 4)
  B <-  matrix(treat_mu, 4, 1)
  betas <- solve(A, B)

  treat_summary <-  data.frame(treat_mu, treat_var, treat_cor, betas)
  row.names(treat_summary) <-  c("(1, 1)", "(1, -1)", "(-1, 1)", "(-1,-1)" )

  recipe <- list(cell_mu=cell_mu, cell_var=cell_var, cell_cor=cell_cor, cell_cov=cell_cov,
                 treat_mu=treat_mu, treat_var=treat_var, treat_cov=treat_cov, treat_cor=treat_cor,
                 p_1=p_1, p_2=p_2, p_A1=p_A1, p_A2=p_A2, N=N, m=m, betas=betas, treat_summary=treat_summary)
  return(recipe)
}

#' @title Function to output the true effect size from the recipe
#' @param recipe the generated recipe from nocovairiates.from_conditional() function
#' @param d1 DTR 1
#' @param d2 DTR 2
#' @return the true effect size of two DTRs
#' @export nocovariates.true_effect_size
nocovariates.true_effect_size <- function(recipe, d1, d2){
  treat_mu <- recipe$treat_mu
  treat_var <- recipe$treat_var
  dtrnum <- list("1,1" = 1, "1,-1"=2, "-1,1"=3, "-1,-1"=4)
  numerator <- treat_mu[dtrnum[[d1]]] - treat_mu[dtrnum[[d2]]]
  denominator <- (0.5 * treat_var[dtrnum[[d1]]])+(0.5 * treat_var[dtrnum[[d2]]])
  denominator <- sqrt(denominator)
  return(numerator/denominator)
}


#' @title Function to generate data from the treatment recipe object
#' @description Generates data from the treatment recipe object. This function can be
#' reused with the same treatment recipe object to generate more samples with the same parameters.
#' @param recipe the generated recipe from nocovairiates.from_conditional() function
#' @return a data.frame object of the generated data with cluster ID, participant ID, A1, R, A2, and response
#' @export nocovariates.treat
nocovariates.treat <- function(recipe){
  encoder <- list("1,1,."=1, "1,0,1"=2, "1,0,-1"=3, "-1,1,."=4, "-1,0,1"=5, "-1,0,-1"=6)
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

  data <- matrix(nrow=0, ncol=6)

  for (clusID in 1:N){
    m <- m_vec[clusID]
    A1 <-  (2*rbinom(1, 1, p_A1)) - 1

    if (A1 == 1){
      R <-  rbinom(1,1,p_1)
      if (R == 1){
        treat_code <- "1,1,."
        A2 <- 1
      } else if (R == 0){
        A2 <-  (2*rbinom(1, 1, p_A2)) - 1
        treat_code <- sprintf("1,0,%d", A2)
      }
    } else if (A1 == -1){
      R <-  rbinom(1, 1 ,p_1)
      if (R == 1){
        treat_code <- "-1,1,."
        A2 <- 1
      }
      else if (R == 0){
        A2 <-  (2*rbinom(1, 1, p_A2)) - 1
        treat_code <- sprintf("-1,0,%d", A2)
      }
    }

    cell_index <- as.numeric(encoder[treat_code])
    V <-  cell_var[cell_index] * (diag(m)+(cell_cor[cell_index]*(matrix(1,m,m)-diag(m))))
    Y <- t(cell_mu[cell_index]+mvtnorm::rmvnorm(1, rep(0,m), V, method='chol'))

    tempmat <- cbind(rep(clusID, m), seq(1,m), rep(A1, m), rep(R, m), rep(A2, m), Y)
    data <- rbind(data, tempmat)
  }
  full_data_frame = as.data.frame(data)
  names(full_data_frame) = c("i", "j", "A1", "R", "A2", "Y")
  return(full_data_frame)
}

#' @title Function to estimate beta, variance of beta_hat, and ICC
#' @param data a data.frame object containing cluster ID, participant ID, A1, R, A2, Y with
#' column names being i, j, A1, R, A2, Y
#' @param mode if set to 'all', the function returns the estimate for the covariance matrix of beta_hat
#' @return estimates object containing all the estimated parameters.
#' @export nocovariates.estimate
nocovariates.estimate <- function(data, mode='all'){
  var <- list("1,1" = 1, "1,-1"=1, "-1,1"=1, "-1,-1"=1)
  rho <- list("1,1" = 0, "1,-1"=0, "-1,1"=0, "-1,-1"=0)
  betas <- get_beta(var, rho, data)
  epsilon <- get_reiduals(betas, data)
  varrho_stars <- get_varrho_star(data, epsilon)
  var <- varrho_stars$var_star
  rho <- varrho_stars$rho_star

  betas <- get_beta(var, rho, data)
  epsilon <- get_reiduals(betas, data)
  varrho_stars <- get_varrho_star(data, epsilon)
  var <- varrho_stars$var_star
  rho <- varrho_stars$rho_star

  betas <- get_beta(var, rho, data)
  epsilon <- get_reiduals(betas, data)
  varrho_stars <- get_varrho_star(data, epsilon)
  var <- varrho_stars$var_star
  rho <- varrho_stars$rho_star
  B <- matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1),4,4)
  cov_hat <- nocovariates.estimate_var_betas(data, var, rho, betas)
  mu_hat_dtr <- B %*% betas
  return(list(beta_hat=betas, var_hat=var, rho_hat=rho, cov_hat_beta_hat=cov_hat, mu_hat_dtr=mu_hat_dtr))
}

rename_clusters <- function(df, id){
  ni <- nrow(df)
  df$i <- rep(id, ni)
  df$j <- 1:ni
  return(df)
}

clusterBoot <- function(data, n=NULL){
  N <- unique(data$i)
  if (is.null(n)){
    n <- length(N)
  }
  bootdata <- data.frame()
  clusters <- sample(N, n, replace=TRUE)
  for (i in 1:n){
    c <- clusters[i]
    bootdata <- rbind(bootdata, rename_clusters(data[data$i == c, ], i))
  }
  rownames(bootdata) <- 1:nrow(bootdata)
  return(bootdata)
}

#' @title Function to prepare data for estimation functions
#' @description nocovariates.estimate requires data as a data.frame object containing cluster ID, participant ID, A1, R, A2, Y with
#' column names being i, j, A1, R, A2, Y with cluster ID (i) being a contiguous integer from 1 to N and participant ID (j) being a contiguous integer from 1 to m_i
#' where N is the number of clusters and m_i is the number of participants in cluster i.
#' This helper function prepares data for estimation.
#' @param data data frame to be modified
#' @param cluster_ID the column name corresponding to the cluster ID in the data
#' @param participant_ID the column name corresponding to the participant ID in the data
#' @param A1 the column name corresponding to A1 in the data
#' @param R the column name corresponding to R in the data
#' @param A2 the column name corresponding to A2 in the data
#' @param Y the column name corresponding to the outcome variable in the data
#' @return a data.frame object with renamed columns, cluster ID and participant ID.
#' @export nocovariates.data_prep
nocovariates.data_prep <- function(data, cluster_ID, participant_ID, A1, R, A2, Y){
  names(data)[names(data) == cluster_ID] <- "i"
  names(data)[names(data) == participant_ID] <- "j"
  names(data)[names(data) == A1] <- "A1"
  names(data)[names(data) == R] <- "R"
  names(data)[names(data) == A2] <- "A2"
  names(data)[names(data) == Y] <- "Y"
  N_array <- unique(data$i)
  N <- length(N_array)
  data_renamed <- data.frame()
  for (i in 1:N){
    data_renamed <- rbind(data_renamed, rename_clusters(data[data$i == N_array[i], ], i))
  }
  rownames(data_renamed) <- 1:nrow(data_renamed)
  return(data_renamed)
}

#' @title Function to estimate effect size for two DTRs
#' @param estimates estimates object returned by the nocovariates.estimate() function
#' @param d1 DTR 1
#' @param d2 DTR 2
#' @return estimated effect size of type numeric
#' @export nocovariates.estimate_effectsize
nocovariates.estimate_effectsize <- function(estimates, d1, d2){
  var_hat <- estimates$var_hat
  mu_hat_dtr <- estimates$mu_hat_dtr
  dtrnum <- list("1,1" = 1, "1,-1"=2, "-1,1"=3, "-1,-1"=4)
  numerator <- mu_hat_dtr[dtrnum[[d1]]] - mu_hat_dtr[dtrnum[[d2]]]
  denominator <- (0.5 * var_hat[[d1]])+(0.5 * var_hat[[d2]])
  denominator <- sqrt(denominator)
  return(numerator/denominator)
}

#' @title Function to estimate effect size for all DTR pairs
#' @param estimates estimates object returned by the nocovariates.estimate() function
#' @return a named list of estimated effect sizes
#' @export nocovariates.estimate_effectsizes
nocovariates.estimate_effectsizes <- function(estimates){
  all_pairs <- t(combn(c("1,1", "1,-1", "-1,1", "-1,-1"), 2))
  effect_sizes <- list()
  for (i in 1:6){
    d1 <- all_pairs[i,1]
    d2 <- all_pairs[i,2]
    es <- nocovariates.estimate_effectsize(estimates, d1, d2)
    effect_sizes[sprintf("(%s), (%s)", d1, d2)] <- es
  }
  return(effect_sizes)
}

boot_var_es <- function(i, data, d1, d2){
  data_boot <- clusterBoot(data)
  estimates <- nocovariates.estimate(data_boot)
  return(nocovariates.estimate_effectsize(estimates, d1, d2))
}

#' @title Function to estimate variance of the effect size for two DTRs
#' @param estimates estimates object returned by the nocovariates.estimate() function
#' @param d1 DTR 1
#' @param d2 DTR 2
#' @param data data to be passed as an argument for 'bootstrap' method
#' @param method can be either 'naive' or 'bootstrap'
#' @param numBoot the number of bootstrap simulations to be performed
#' @param numCores number of parallel cores to be passed for parallel processing
#' @return estimated variance of the effect size of type numeric
#' @export nocovariates.estimate_var_effectsize
nocovariates.estimate_var_effectsize <- function(estimates, d1, d2, data, method='naive', numBoot=1000, numCores=1){
  var_hat <- estimates$var_hat
  dtrnum <- list("1,1" = 1, "1,-1"=2, "-1,1"=3, "-1,-1"=4)
  i1 <- dtrnum[[d1]]
  i2 <- dtrnum[[d2]]
  if (method == 'naive'){
    B <- matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1),4,4)
    M <- estimates$cov_hat_beta_hat
    var_hat_mu_dtr <- B%*%M%*%B
    numerator <- var_hat_mu_dtr[i1, i1] + var_hat_mu_dtr[i2, i2] - (2* var_hat_mu_dtr[i1, i2])
    denominator <- (0.5 * var_hat[[d1]])+(0.5 * var_hat[[d2]])
    return(numerator/denominator)

  } else if (method == 'bootstrap'){
    results_1 <- pbmcapply::pbmclapply(1:numBoot, boot_var_es, data=data, d1='1,1', d2='1,-1', mc.cores = numCores)
    #results_1 <- lapply(1:1000, boot_var_es, data=data, d1='1,1', d2='1,-1')
    return(var(unlist(results_1)))
  }
}

