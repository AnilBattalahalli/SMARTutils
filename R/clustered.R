#' @importFrom stats rbinom
#' @importFrom stats var
#' @importFrom stats terms
#' @importFrom dplyr %>%
#' @importFrom utils combn
#' @importFrom pbmcapply pbmclapply
NULL

Ii <- function(A1, R, A2, a1, a2){
  ifelse(A1==a1,1,0)*(R+(ifelse(A2==a2,1,0)*(1-R)))
}

Wi <- function(R){
  2*(R+(2*(1-R)))
}

IWi <- function(A1, R, A2, a1, a2){
  2*ifelse(A1==a1,1,0)*(R+(2*ifelse(A2==a2,1,0)*(1-R)))
}

perm_mul_sum <- function(l){
  l <- rbind(t(combn(l, 2)), t(combn(l, 2)))
  l <- data.frame(l)
  return(sum(l$X1*l$X2))
}

D <- function(data_subset, a1, a2, formula){
  data_subset$a1 <- rep(a1, nrow(data_subset))
  data_subset$a2 <- rep(a2, nrow(data_subset))
  df <- model.frame(formula, data_subset)
  df <- within(df, rm("Y"))
  I <- rep(1, nrow(df))
  df <- as.matrix(cbind(I, df)) %>% unname()
  return(df)
}

mmm <- function(designmat, betas){
  mu <- designmat %*% matrix(betas)
  return(mu)
}

get_residuals <- function(betas, data, formula){
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <- list("1,1" = c(1,1), "1,-1"=c(1,-1), "-1,1"=c(-1,1), "-1,-1"=c(-1,-1))
  N <- data$i %>% unique %>% length
  epsilon <- list()
  for (d in dtrs){
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    temp <- list()
    for (ci in unique(data$i)){
      data_sub <- data[data$i == ci, ]
      temp[ci] <- list(as.numeric(data_sub$Y - mmm(D(data_sub, a1, a2, formula), betas)))
    }
    epsilon[d] <- list(temp)
  }
  return(epsilon)
}

get_varrho_star <- function(data, epsilon, weights){
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <- list("1,1" = c(1,1), "1,-1"=c(1,-1), "-1,1"=c(-1,1), "-1,-1"=c(-1,-1))
  var_star <- list()
  rho_star <- list()
  N <- data$i %>% unique %>% length

  for (d in dtrs){
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    s <- 0
    den <- 0
    for (ci in unique(data$i)){
      data_sub <- data[data$i == ci, ]
      ni <- nrow(data_sub)
      A1 <- data_sub$A1[1]
      R <- data_sub$R[1]
      A2 <- data_sub$A2[1]
      eps_sq <- epsilon[[d]][[ci]]^2
      if (is.null(weights) == F){
        W = weights[weights$i == ci,]$W[1] * Ii(A1, R, A2, a1, a2)
      } else{
        W = IWi(A1, R, A2, a1, a2)
      }
      if (W != 0){
        s <- s + (W * sum(eps_sq))
        den <- den + (W * ni)
      }
    }
    var_star[d] <- s/den
  }

  for (d in dtrs){
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    s <- 0
    den <- 0
    for (ci in unique(data$i)){
      data_sub <- data[data$i == ci, ]
      ni <- nrow(data_sub)
      if (ni > 1){
        A1 <- data_sub$A1[1]
        R <- data_sub$R[1]
        A2 <- data_sub$A2[1]
        eps_mul_sum <- perm_mul_sum(epsilon[[d]][[ci]])
        if (is.null(weights) == F){
          W = weights[weights$i == ci,]$W[1] * Ii(A1, R, A2, a1, a2)
        } else{
          W = IWi(A1, R, A2, a1, a2)
        }
        if (W != 0){
          s <- s + (W * eps_mul_sum)
          den <- den + (W * ni * (ni - 1))
        }
      }
    }
    den <- den * var_star[[d]]
    rho_star[d] <- s/den
  }
  rhovar_stars <- list(rho_star=rho_star, var_star=var_star)
  return(rhovar_stars)
}

get_beta <- function(var, rho, data, formula, weights){
  params <- (formula %>% terms %>% labels %>% length)+1
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <- list("1,1" = c(1,1), "1,-1"=c(1,-1), "-1,1"=c(-1,1), "-1,-1"=c(-1,-1))
  f <- matrix(0, params, params)
  N <- data$i %>% unique %>% length
  for (ci in unique(data$i)){
    data_sub <- data[data$i == ci, ]
    A1 <- data_sub$A1[1]
    R <- data_sub$R[1]
    A2 <- data_sub$A2[1]
    ni <- nrow(data_sub)
    for (d in dtrs){
      a1 <- dtrenc[[d]][1]
      a2 <- dtrenc[[d]][2]
      if (is.null(weights) == F){
        W = weights[weights$i == ci,]$W[1] * Ii(A1, R, A2, a1, a2)
      } else{
        W = IWi(A1, R, A2, a1, a2)
      }
      if (W != 0){
        V = var[[d]] * (diag(ni)+(rho[[d]]*(matrix(1,ni,ni)-diag(ni))))
        f = f + (W * t(D(data_sub, a1, a2, formula)) %*% solve(V) %*% D(data_sub, a1, a2, formula))
      }
    }
  }

  s <- matrix(0, params, 1)
  for (ci in unique(data$i)){
    data_sub <- data[data$i == ci, ]
    A1 <- data_sub$A1[1]
    R <- data_sub$R[1]
    A2 <- data_sub$A2[1]
    ni <- nrow(data_sub)
    Y <- matrix(data_sub$Y)
    for (d in dtrs){
      a1 <- dtrenc[[d]][1]
      a2 <- dtrenc[[d]][2]
      if (is.null(weights) == F){
        W = weights[weights$i == ci,]$W[1] * Ii(A1, R, A2, a1, a2)
      } else{
        W = IWi(A1, R, A2, a1, a2)
      }
      if (W != 0){
        V = var[[d]] * (diag(ni)+(rho[[d]]*(matrix(1,ni,ni)-diag(ni))))
        s <- s + (IWi(A1, R, A2, a1, a2) * t(D(data_sub, a1, a2, formula)) %*% solve(V) %*% Y)
      }
    }
  }
  beta_hat <- solve(f) %*% s
  return(beta_hat)
}

clustered.estimate_var_betas <- function(data, var, rho, betas, formula, weights){
  params <- (formula %>% terms %>% labels %>% length)+1
  dtrs <- c("1,1", "1,-1", "-1,1", "-1,-1")
  dtrenc <- list("1,1" = c(1,1), "1,-1"=c(1,-1), "-1,1"=c(-1,1), "-1,-1"=c(-1,-1))
  J_hat <- matrix(0, params, params)
  N <- data$i %>% unique %>% length
  for (d in dtrs){
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    temp <- matrix(0, params, params)
    for (ci in unique(data$i)){
      data_sub <- data[data$i == ci, ]
      ni <- nrow(data_sub)
      A1 <- data_sub$A1[1]
      R <- data_sub$R[1]
      A2 <- data_sub$A2[1]
      if (is.null(weights) == F){
        W = weights[weights$i == ci,]$W[1] * Ii(A1, R, A2, a1, a2)
      } else{
        W = IWi(A1, R, A2, a1, a2)
      }
      if (W != 0){
        V = var[[d]] * (diag(ni)+(rho[[d]]*(matrix(1,ni,ni)-diag(ni))))
        temp <- temp + (W * t(D(data_sub, a1, a2, formula)) %*% solve(V) %*% D(data_sub, a1, a2, formula))
      }
    }
    J_hat <- J_hat+temp
  }
  J_hat <-  J_hat/N

  A_hat <- matrix(0, params, params)
  for (d in dtrs){
    temp <- matrix(0, params, params)
    a1 <- dtrenc[[d]][1]
    a2 <- dtrenc[[d]][2]
    for (ci in unique(data$i)){
      data_sub <- data[data$i == ci, ]
      A1 <- data_sub$A1[1]
      R <- data_sub$R[1]
      A2 <- data_sub$A2[1]
      ni <- nrow(data_sub)
      Y <- matrix(data_sub$Y)
      diff <- Y - mmm(D(data_sub, a1, a2, formula), betas)
      V = var[[d]] * (diag(ni)+(rho[[d]]*(matrix(1,ni,ni)-diag(ni))))
      if (is.null(weights) == F){
        W = weights[weights$i == ci,]$W[1] * Ii(A1, R, A2, a1, a2)
      } else{
        W = IWi(A1, R, A2, a1, a2)
      }
      if (W != 0){
        Ui <- W * t(D(data_sub, a1, a2, formula)) %*% solve(V) %*% diff
        temp <- temp + (Ui %*% t(Ui))
      }
    }
    A_hat <- A_hat + temp
  }
  A_hat <- A_hat/N

  var_hat_beta_hat <- (1/N) * (solve(J_hat) %*% A_hat %*% solve(J_hat))
  return(var_hat_beta_hat)
}

#' @title Function to prepare data for estimation
#' @param data data frame to be prepared for estimation
#' @param i cluster ID column name
#' @param A1 A1 column name
#' @param R R column name
#' @param A2 A2 column name
#' @param Y outcome variable column name
#' @return a data frame ready to be plugging into the estimation function
#' @export clustered.prepare_data
clustered.prepare_data <- function(data, i='i', A1='A1', R='R', A2='A2', Y='Y'){
  names(data)[names(data) == i] <- 'i'
  names(data)[names(data) == A1] <- 'A1'
  names(data)[names(data) == R] <- 'R'
  names(data)[names(data) == A2] <- 'A2'
  names(data)[names(data) == Y] <- 'Y'
  return(data)
}

#' @title Function to estimate beta, variance of beta_hat, and ICC
#' @param formula formula for the estimation
#' @param data a data.frame object containing cluster ID, participant ID, A1, R, A2, Y with
#' column names being i, j, A1, R, A2, Y
#' @param est_var T/F to estimate the variance or not
#' @param weights a data frame containing cluster level weights with columns `i` and `W` for cluster ID and weights respectively
#' @return estimates object containing all the estimated parameters.
#' @export clustered.estimate
clustered.estimate <- function(formula, data, est_var = T, weights=NULL){
  var <- list("1,1" = 1, "1,-1"=1, "-1,1"=1, "-1,-1"=1)
  rho <- list("1,1" = 0, "1,-1"=0, "-1,1"=0, "-1,-1"=0)

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

  if (est_var){
    cov_hat <- clustered.estimate_var_betas(data, var, rho, betas, formula, weights)
    return(list(beta_hat=betas, var_hat=var, rho_hat=rho, cov_hat_beta_hat=cov_hat, formula=formula))
  } else{
    return(list(beta_hat=betas, var_hat=var, rho_hat=rho, formula=formula))
  }

}
