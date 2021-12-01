#' @title Residualization
#' @importFrom RcppEigen fastLm
#' @author Yiqing Xu, Eddie Yang
#' @examples
#' @export

fLM <- function(y, x){
	out <- RcppEigen::fastLm(X=x, y=y)
	return(out$residuals)
}


getResidual <- function(P, pos.list=NULL){
	n <- length(pos.list)-1
  for (i in 1:n){
    old <- sum(pos.list[1:i])
    new <- old + pos.list[i+1]
    residual_list <- apply(as.matrix(P[,(old+1):new]), 2, fLM, x=as.matrix(P[,1:old]))
    P[,(old+1):new] <- as.matrix(residual_list)
  }
	#P <- scale(P)
	P <- as.data.frame(P)
	return(P)
}


loss_func0 <- function(coefs, co_x, tr_total, Q, alpha){
  XS = - co_x %*% coefs
  norm = abs(coefs)
  maxXS = max(XS)
  loss <- maxXS + log(t(Q) %*% exp( XS - maxXS)) + t(tr_total) %*% coefs + sum(alpha*norm)
  return(c(loss))
}


entbal_wts <- function(Q, co_x, coefs){
  V <- - co_x %*% coefs
  maxV <- max(V)
  norm_c <- t(Q) %*% exp( V - maxV )
  return (Q * exp( V - maxV ) / c(norm_c))
}





