#' @importFrom stats rbinom
#' @importFrom stats var
#' @importFrom stats terms
#' @importFrom stats model.frame
#' @importFrom dplyr %>%
#' @importFrom utils combn
#' @importFrom pbmcapply pbmclapply
NULL

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
#' @export nocovariates.from_conditional
nocovariates.from_conditional <- function(cell_mu = c(10,8,2,4,-2,3), cell_var = c(100,100,100,100,100,100),
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
