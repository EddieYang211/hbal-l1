#' @title Ridge Penalty Selection through Cross Validation
#' @aliases CrossValidate
#' @description Internal function called by \code{hbal} to select ridge penalties through cross-validation.
#' @param alpha                   alpha. Controls degree of regularization.
#' @param grouping                different groupings of the covariates.
#' @param folds                   number of folds to perform cross validation. 
#' @param treatment               covariate matrix for treatment group.
#' @param fold.co                 fold assignments for control units.
#' @param fold.tr                 fold assignments for treated units.
#' @param coefs                   starting coefficients (lambda).
#' @param control                 covariate matrix for control group.
#' @param constraint.tolerance    tolerance level for imbalance.
#' @param print.level             details of printed output.
#' @param base.weight             target weight distribution for the control units.
#' @param full.t                  (unresidualized) ovariate matrix for treatment group.
#' @param full.c                  (unresidualized) ovariate matrix for control group.
#' @param shuffle.treat           whether to create folds for the treated units
#' @return alpha, lambda
#' @importFrom stats na.omit coef optim
#' @importFrom glmnet cv.glmnet
#' @importFrom nloptr nloptr
#' @author Yiqing Xu, Eddie Yang

crossValidate <- function(
	alpha=NULL,
	grouping=NULL,
	folds=NULL,
	treatment = NULL,
	fold.co = NULL,
	fold.tr=NULL,
	coefs=NULL,
	vcontrol = NULL,
	constraint.tolerance = NULL,
	print.level = NULL,
	base.weight = NULL,
	full.t=NULL,
	full.c=NULL,
	shuffle.treat=NULL){

	if (any(!is.finite(alpha))){
		return(Inf)
	}
	
	res <- rep(NA, folds) #store cross validation results
	coe <- NULL
	opt_constraints <- c(-5, 5)

	# loop over each alpha value

	penalty <- rep(c(0, alpha), times=grouping)
	sub.coef <- Coefs <- rep(0, length(coefs))

	counter <- 0
	# loop over each fold for each alpha
	for (k in 1:folds){
		co.test.k <- which(fold.co==k)
		tr.test.k <- which(fold.tr==k)
		base.w <- base.weight[-co.test.k]
		train.control <- vcontrol[-co.test.k,]
		train.treat <- treatment[-tr.test.k,]
		train.total <- c(1, colMeans(train.treat))
		test.control <- vcontrol[co.test.k,]
		test.treat <- full.t[tr.test.k,]
		test.cc <- full.c[co.test.k,]
		if(!shuffle.treat){
			train.total <- c(1, colMeans(treatment))
			test.treat <- full.t
		}
		if(!is.null(sub.coef) && all(is.finite(sub.coef)) && max(abs(sub.coef))<=10) Coefs <- sub.coef

#		out <- try(optim(par=as.matrix(Coefs), 
#			loss_func0, 
#			co_x=train.control, 
#			tr_total=as.matrix(train.total), 
#			Q=base.w, 
#			alpha = as.matrix(penalty),
#			method='BFGS',
#			lower = opt_constraints[1], 
#			upper = opt_constraints[2],
#			control = list(lmm = 25, factr = 1e-10, maxit = 500, trace=1)
#			),
#		silent=TRUE
#		)

#		out <- try(nmkb(par=as.matrix(Coefs), 
#			loss_func0, 
#			control = list(tol=1e-4),
#			lower = opt_constraints[1], 
#			upper = opt_constraints[2],
#			co_x=train.control, 
#			tr_total=as.matrix(train.total), 
#			Q=base.w, 
#			alpha = as.matrix(penalty)
#			),
#		silent=TRUE
#		)		
		out <- nloptr(x0 = as.matrix(Coefs), 
                eval_f = loss_func0,
                lb = rep(opt_constraints[1], length(Coefs)),
                ub = rep(opt_constraints[2], length(Coefs)), 
                opts = list('algorithm'='NLOPT_LN_BOBYQA',
                            'maxeval' =500,
                            'print_level'=print.level),
                co_x=train.control, 
                tr_total=as.matrix(train.total), 
                Q=base.w, 
                alpha = as.matrix(penalty))

		if (class(out)!="try-error"){
			coe <- as.matrix(out$solution)
			weights <- c(entbal_wts(base.weight[co.test.k], test.control, coe))
			test.weight <- all(is.finite(weights)) 
			test.coef <- all(is.finite(coe))
			if(test.weight && test.coef){
				test.cr.mean <- c(weights %*% test.cc)
				tr.t <- c(1, colMeans(test.treat))
				dif <- mean(abs(tr.t-test.cr.mean))
				res[k] <- dif
				counter <- counter + 1
			}
		}
		sub.coef <- coe
	}#end of inner loop

	oo <- ifelse(is.finite(mean(res)), mean(res), Inf)
	return(oo)
}



