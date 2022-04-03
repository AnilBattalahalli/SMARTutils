rename_clusters <- function(df, id){
  ni <- nrow(df)
  df$i <- rep(id, ni)
  return(df)
}

#' @title Function to draw cluster-level bootstrap samples
#' @param data data.frame object containing cluster ID, participant ID, A1, R, A2, Y with
#' column names being i, j, A1, R, A2, Y.
#' @param n number of bootstrap samples to be drawn
#' @return data.frame object with the bootstrapped samples
#' @export clusterBoot

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
