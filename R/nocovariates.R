#' @importFrom stats rbinom
NULL

#' @title Function to get treatment recipe from the conditionals
#'
#' @description Does what it says
#' @param cell_mu a vector of conditional means
#' @param cell_var a vector of conditional vars
#' @param cell_cor a vector of conditional correlation
#' @param N number of clusters
#' @param m number of participants in a cluster
#' @param p_1 probability of response of A_1 = 1
#' @param p_2 probability of response of A_1 = -1
#' @export nocovairiates.from_conditional
nocovairiates.from_conditional <- function(cell_mu = c(10,8,2,4,-2,3), cell_var = c(100,100,100,100,100,100),
                                           cell_cor = c(0.1,0.2,0.15,0.12,0.11,0.07), N=100, m=10, p_1=0.3, p_2=0.2){

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

#' @title Function to generate data from treatment
#' @param recipe the generated recipe from nocovairiates.from_conditional()
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
    Y <- t(cell_mu[cell_index]+mvtnorm::rmvnorm(1, rep(0,m), V*0.5, method='chol'))

    tempmat <- cbind(rep(clusID, m), seq(1,m), rep(A1, m), rep(R, m), rep(A2, m), Y)
    data <- rbind(data, tempmat)
  }
  full_data_frame = as.data.frame(data)
  names(full_data_frame) = c("i", "j", "A1", "R", "A2", "Y")
  return(full_data_frame)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

IWi <- function(A1, R, A2, a1, a2){
  2*ifelse(A1==a1,1,0)*(R+(2*ifelse(A2==a2,1,0)*(1-R)))
}

D <- function(n, a1, a2){
  temp <- c(1, a1, a2, a1*a2)
  rep.row(temp, n)
}

#' @title Function to estimate betas given the data, variance and correlation
#' @param var list of variances
#' @param rho list of correlations
#' @param data data with cluster ID, A1, R, A2, Y
#' @export nocovariates.estimate_betas
nocovariates.estimate_betas <- function(var, rho, data){
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
      V = var[[d]] * (diag(ni)+(rho[[d]]*(matrix(1,ni,ni)-diag(ni))))
      f = f + (IWi(A1, R, A2, a1, a2) * t(D(ni, a1, a2)) %*% solve(V) %*% D(ni, a1, a2))
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
      V = var[[d]] * (diag(ni)+(rho[[d]]*(matrix(1,ni,ni)-diag(ni))))
      s <- s + (IWi(A1, R, A2, a1, a2) * t(D(ni, a1, a2)) %*% solve(V) %*% Y)
    }
  }
  beta_hat <- solve(f) %*% s
  return(beta_hat)
}
