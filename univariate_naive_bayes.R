library(glmnet)
library(methods)
library(stats)
library(utils)
library(Biobase)
library(pROC)
library(ggplot2)

#' Posterior probabilities of FR given G and E.
#'
#' \code{getFuncRvPosteriors} computes posterior probabilities of functionality
#'         of regulatory variant (FR) given genomic features (G) and outlier
#'         status (E) with current estimate of beta (parameters between FR and G)
#'         and theta (parameters between FR and E).
#'
#' @param Out Binary values of outlier status (E).
#' @param probFuncRvFeat probabilities of FR given genomic features and estimated
#'         beta, P(FR | G, beta), from \code{getFuncRvFeat}.
#' @param theta Current estimate of theta.
#'
#' @return posterior probabilities of FR (P(FR | G, E, beta, theta)) and probable
#'         status of FR.
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' theta.init<-matrix(c(.99, .01, .3, .7), nrow=2)
#' costs<-c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha=0, nfolds=10)
#' probFuncRvFeat <- getFuncRvFeat(Feat, logisticAllCV$glmnet.fit,
#'         logisticAllCV$lambda.min)
#' posteriors <- getFuncRvPosteriors(Out, probFuncRvFeat, theta=theta.init)
#'
#' @export


observed_naive_bayes_likelihood <- function(x, theta) {
  return(theta[1,x[1]])
}


getFuncRvPosteriors <- function(Out, probFuncRvFeat, theta_cur) {

  probOut_FuncRv <- matrix(NA,dim(Out)[1], 2)
  probOut <- matrix(NA, dim(Out)[1], 1)

  probOut_FuncRvInlier <- apply(Out, 1, observed_naive_bayes_likelihood, theta_cur$inlier_component)
  probOut_FuncRvOutlier <- apply(Out, 1, observed_naive_bayes_likelihood, theta_cur$outlier_component)
  probOut_FuncRv <- cbind(probOut_FuncRvInlier, probOut_FuncRvOutlier)

  probOut <- rowSums(probOut_FuncRv*cbind(1.0-probFuncRvFeat,probFuncRvFeat))
  post <- probOut_FuncRv*c(1-probFuncRvFeat,probFuncRvFeat)/c(probOut, probOut)
  list(posterior=post, mle = max.col(post)-1)
}



mle_theta <- function(Out, FuncRv, pseudoc,num_bins) {
  dimension <- dim(Out)[2]
  theta_outlier <- matrix(1,dimension,num_bins)
  theta_inlier <- matrix(1,dimension,num_bins)
  for (bin_number in 1:num_bins) {
    theta_outlier[,bin_number] <- colSums(((Out==bin_number)*FuncRv[,2]),na.rm=TRUE)
    theta_inlier[,bin_number] <- colSums(((Out==bin_number)*FuncRv[,1]),na.rm=TRUE)
  }
  theta_outlier <- theta_outlier + pseudoc
  theta_inlier <- theta_inlier + pseudoc

  theta_outlier <- theta_outlier/rowSums(theta_outlier)
  theta_inlier <- theta_inlier/rowSums(theta_inlier)
  theta_cur <- list(inlier_component = theta_inlier, outlier_component = theta_outlier)
  return(theta_cur)
}

#' Maximum likelihoood estimate of beta.
#'
#' \code{mleBeta} computes maximum likelihoood estimate of beta (parameters
#'         between FR (functionality of regulatory variant) and G (genomic
#'         annotations); multivariate logistic regression).
#'
#' @param Feat Genomic features (G)
#' @param FuncRv Soft-assignments of FR from E-step
#' @param costs Candidate penalty parameter values for L2-regularization within
#'         logistic regression.
#'
#' @return MLE of beta
#'
#' @section Warning: To input a vector of candidate penalty values makes
#'         \code{glmnet} faster than to input a single penalty value
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' theta.init <- matrix(c(.99, .01, .3, .7), nrow=2)
#' costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha=0, nfolds=10)
#' probFuncRvFeat <- getFuncRvFeat(Feat, logisticAllCV$glmnet.fit, logisticAllCV$lambda.min)
#' posteriors <- getFuncRvPosteriors(Out, probFuncRvFeat, theta=theta.init)
#' logistic.cur <- mleBeta(Feat, FuncRv=posteriors$posterior, costs)
#'
#' @seealso \code{\link[glmnet]{glmnet}}
#'
#' @export

mleBeta <- function(Feat, FuncRv, costs) {
  glmnet(Feat, FuncRv, lambda=costs, family="binomial", alpha = 0)
}

#' Posterior probabilities of FR given G
#'
#' \code{getFuncRvFeat} computes posterior probabilities of FR (functionality of
#'         regulatory variant) given G (genomic features) and current estimate
#'         of beta (parameters between FR and G).
#'
#' @param Feat Genomic features (G)
#' @param logistic.model Logistic regression model with current estimate of beta
#' @param lambda Selected lambda
#'
#' @return probabilities of FR given genomic features, P(FR | G)
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha = 0, nfolds=10)
#' probFuncRvFeat <- getFuncRvFeat(Feat, logistic.model=logisticAllCV$glmnet.fit,
#'         lambda=logisticAllCV$lambda.min)
#'
#' @seealso \code{\link{predict}}
#'
#' @export

getFuncRvFeat <- function(Feat, logistic.model, lambda) {
  predict(logistic.model, Feat, s=lambda, type="response")
}

#' Test posterior probabilities of FR given G and E
#'
#' \code{testPosteriors} computes posterior probabilities of FR (functionality
#'         of regulatory variant) given G (genomic annotations) and E (outlier
#'         status) with estimate of beta (parameters between FR and G) and
#'         theta (parameters between FR and E).
#'
#' @param Feat Genomic features (G)
#' @param Out Binary values of outlier status (E).
#' @param emModel Estimated parameters including beta and theta via EM and
#'         selected lambdas
#'
#' @return test posterior probabilities of FR given new outlier status (E)
#'         and genomic features (G), P(FR | G, E, beta, theta), and probable
#'         status of FR.
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' theta.init <- matrix(c(.99, .01, .3, .7), nrow=2)
#' costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha = 0, nfolds=10)
#' emModelAll <- integratedEM(Feat, Out, logisticAllCV$lambda.min, logisticAllCV$glmnet.fit,
#'         pseudoc=50, theta.init, costs, verbose=FALSE)
#' trainedpost <- testPosteriors(Feat, Out, emModel=emModelAll)
#'
#' @seealso \code{\link{getFuncRvFeat}} and \code{\link{getFuncRvPosteriors}}
#'
#' @export

testPosteriors <- function(Feat, Out, emModel) {
  probFuncRvFeat <- getFuncRvFeat(Feat, emModel$logistic.model, emModel$lambda)
  getFuncRvPosteriors(Out, probFuncRvFeat, emModel$theta)
}




compute_log_probability <- function(Feat, Out, logistic, beta, lambda, theta_cur,pseudoc) {
    probFuncRvFeat <- getFuncRvFeat(Feat, logistic, lambda)

    probOut_FuncRvInlier <- apply(Out, 1, observed_naive_bayes_likelihood, theta_cur$inlier_component)
    probOut_FuncRvOutlier <- apply(Out, 1, observed_naive_bayes_likelihood, theta_cur$outlier_component)

    probOut_FuncRv <- cbind(probOut_FuncRvInlier, probOut_FuncRvOutlier)

    probOut <- log(rowSums(probOut_FuncRv*cbind(1.0-probFuncRvFeat,probFuncRvFeat)))
    regularizer_term <- 0
    for (j in 1:length(beta)) {
        regularizer_term <- regularizer_term + dnorm(beta[j],mean=0,sd=sqrt((1.0/lambda)),log=TRUE)
    }
    log_prob <- sum(probOut) + regularizer_term
    return(log_prob)
}


#' An iterative expectation-maximization algorithm for RIVER
#'
#' \code{integratedEM} iteratively executes e-step and m-step until it
#'         converges. This is a main function of RIVER.
#'
#' @param Feat Genomic features (G).
#' @param Out Binary values of outlier status (E).
#' @param lambda Selected lambda.
#' @param logistic.init Smart initialization of beta (parameters between
#'         FR and G) from estimate of beta with E via multivariate logistic
#'         regression.
#' @param pseudoc Pseudo count.
#' @param theta.init Initial theta (parameters between FR (functionality
#'         of regulatory variant) and E).
#' @param costs Candidate penalty parameter values for L2-regularization
#'         within logistic regression.
#' @param verbose Logical option for showing extra information on progress.
#'
#' @return Best estimate of beta and theta, final multivariate logistic
#'         regression model, and posterior probabilities of FR.
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @seealso \code{\link{getFuncRvFeat}}, \code{\link{getFuncRvPosteriors}},
#'         \code{\link{mleTheta}}, \code{\link{mleBeta}},
#'         \code{\link[glmnet]{cv.glmnet}},
#'         and \url{https://github.com/ipw012/RIVER}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' theta.init=matrix(c(.99, .01, .3, .7), nrow=2)
#' costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha = 0, nfolds=10)
#' emModelAll <- integratedEM(Feat, Out, lambda=logisticAllCV$lambda.min,
#'         logistic.init=logisticAllCV$glmnet.fit, pseudoc=50, theta=theta.init,
#'         costs, verbose=FALSE)
#'
#' @export

integratedEM <- function(Feat, Out, lambda, logistic.init,
                         theta_init, num_bins, costs,pseudoc,
                         verbose=FALSE){
  theta_cur <- theta_init
  beta.cur <- logistic.init$beta[,which(logistic.init$lambda == lambda)]
  logistic.cur <- logistic.init

  steps <- 1
  maxIter <- 1000  
  converged <- 0
  for (iter in 1:maxIter) {
    if (verbose) {
      cat(' *** STREAM: EM step ',steps,'\n',sep="")
    }

    ## E-step:
    ## Compute expected posterior probabilities
    ##           given current parameters and data
    probFuncRvFeat <- getFuncRvFeat(Feat, logistic.cur, lambda)
    posteriors <- getFuncRvPosteriors(Out, probFuncRvFeat, theta_cur)
    if (verbose) {
      cat('     E-step: Top 10 % Threshold of expected P(FR=1 | G, E): ',
          round(quantile(posteriors$posterior[,2], .9),4),'\n',sep='')
    }
    ## M-step:
    ## Update theta and beta
    theta_old <- theta_cur
    ##HEREEE
    # ML estimate of mu
    theta_cur <- mle_theta(Out, posteriors$posterior, pseudoc, num_bins)

    beta.old <- beta.cur
    # ML estimate of beta
    logistic.cur <- mleBeta(Feat, posteriors$posterior, costs)
    beta.cur <- logistic.cur$beta[,which(logistic.cur$lambda == lambda)]

    # Compute observed log probability
    log_prob <- compute_log_probability(Feat, Out, logistic.cur, beta.cur, lambda, theta_cur,pseudoc)
    cat('    Current log probability: ', log_prob,'\n',sep='')
    print(theta_cur)
    # Print convergence info
    if (verbose) {
      cat('     M-step: norm(theta_inlier difference) = ',
          round(norm(matrix(theta_cur$inlier_component)-matrix(theta_old$inlier_component)),4),
          ', norm(theta_outlier difference) = ',
          round(norm(matrix(theta_cur$outlier_component)-matrix(theta_old$outlier_component)),4),
          ', norm(beta difference) = ',
          round(norm(matrix(beta.cur)-matrix(beta.old)),4),
          " *** \n\n", sep="")
    }

    ## Check convergence
    if ((norm(matrix(beta.cur) - matrix(beta.old)) < 5e-2) &
        (norm(theta_cur$inlier_component - theta_old$inlier_component) < 5e-2) &
        (norm(theta_cur$outlier_component - theta_old$outlier_component) < 5e-2) 
        ) {
      converged <- 1
      break
    }
    steps <- steps + 1
  }

  if (converged == 1) {
    cat(" ::: EM iteration is terminated since it converges within a
        predefined tolerance (0.001) ::: \n\n\n",sep="")
  } else if ((converged == 0) && (iter == maxIter)) {
    cat(" ::: EM iteration is terminated since it reaches a
        predefined maximum value (1000) ::: \n\n\n",sep="")
  }

  list(logistic.model=logistic.cur, beta=beta.cur, theta = theta_cur,
      posteriors=posteriors, lambda=lambda)
}



load_data <- function(input_file, ZscoreThrd=1.5) {
    expData <- read.table(input_file, header=TRUE)
    Feat <- expData[,3:(ncol(expData)-46)] # genomic features
    # sample name as SubjectID:GeneName
    rownames(Feat) <- paste(expData[,"SubjectID"], ":",
    expData[,"GeneName"],sep="")
    Feat <- as.matrix(t(Feat)) # feature x sample
    # outlier status, N2 pairs
    pData <-
        data.frame(Outlier=factor(ifelse(abs(expData[,"neg_log10_median_pvalue"])>=ZscoreThrd,1,0),
        levels=c(0,1)),
        N2pair=factor(expData[,"N2pair"],
        levels=unique(expData[,"N2pair"])))
    rownames(pData) <-
        paste(expData[,"SubjectID"],":",expData[,"GeneName"],sep="")

    # descrition of outlier status and N2 pairs
    metadata <-
        data.frame(labelDescription=c("Outlier status based on Z-scores",
                                  "Pairs of samples having same rare variants"),
               row.names=c("Outlier","N2pair"))
        phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
    dataInput <- ExpressionSet(assayData=Feat, phenoData=phenoData)
    tbt_expression_data <- expData[,(ncol(expData)-45):(ncol(expData)-2)]
    median_expression_data <- expData[,"neg_log10_median_pvalue"]
    return(list(dataInput,median_expression_data, tbt_expression_data))
}

plot_roc <- function(evaROC, ZscoreThrd, figure_file_name) {
  pdf(figure_file_name)
  par(mar=c(6.1, 6.1, 4.1, 4.1))
  plot(NULL, xlim=c(0,1), ylim=c(0,1), 
     xlab="False positive rate", ylab="True positive rate", 
     cex.axis=1.3, cex.lab=1.6)
  abline(0, 1, col="gray")
  lines(1-evaROC$RIVER_spec, evaROC$RIVER_sens, 
      type="s", col='dodgerblue', lwd=2)
  lines(1-evaROC$GAM_spec, evaROC$GAM_sens, 
      type="s", col='mediumpurple', lwd=2)
  legend(0.7,0.2,c("RIVER","GAM"), lty=c(1,1), lwd=c(2,2),
       col=c("dodgerblue","mediumpurple"), cex=1.2, 
       pt.cex=1.2, bty="n")
  title(main=paste("Threshold = ", ZscoreThrd," / AUC: RIVER = ", round(evaROC$RIVER_auc,3), 
                 ", GAM = ", round(evaROC$GAM_auc,3),sep=""))

  dev.off()
}

scatter_plot_fill <- function(data_framer, emModel, figure_output_file) {



  scatter <- ggplot(data = data_framer, mapping = aes(x = neg_log10_median_pvalue, y = GAM_posterior)) + geom_point(aes(colour = posterior), shape = 18,size=3) +
                scale_colour_gradient2(low = "khaki1",midpoint=.5, mid="dodgerblue",high="mediumpurple", guide = "colourbar")
  scatter <- scatter + theme(axis.text.x = element_text(size=20,hjust=.5,vjust=.5,face="plain"),axis.text.y = element_text(size=20,hjust=.5,vjust=.5,face="plain"))
  scatter <- scatter + labs(x = "median(-log10(pvalue))", y = "GAM Posterior") + theme(axis.title=element_text(size=20,face="bold"))

  ggsave(scatter, file=figure_output_file,width = 30,height=15,units="cm")


}



full_data_visualization_driver <- function(input_file, ZscoreThrd, output_root, dimensions, theta_init,num_bins, costs, verbose,pseudoc) {
    ## Extract required data
    ## Extract required data
    all_data <- load_data(input_file, ZscoreThrd)
    dataInput <- all_data[[1]]
    E_all <- all_data[[2]]
    E_tbt_real_valued <- t(as.matrix(all_data[[2]]))

    E_tbt <- discritize_expression_data(E_tbt_real_valued,dim,num_bins)


    ## All genomic features (G)
    FeatAll <- t(exprs(dataInput))
    ## All median log 10 (pvalue) expression data
  
    OutAll <- t(E_tbt)
    OutAll_binary <- as.numeric(unlist(dataInput$Outlier))-1
 
    ## Search a best lambda with outlier status based on 10 cross-validation
    logisticAllCV <- cv.glmnet(FeatAll, as.vector(OutAll_binary), lambda=costs,
                        family="binomial", alpha=0, nfolds=10) # GAM
    if (verbose) {
      cat(' *** best lambda = ',logisticAllCV$lambda.min,' *** \n\n', sep='')
    }

    ## Compute P(FR=1 | G)
    postporbGAM <- predict(logisticAllCV, FeatAll, s="lambda.min", type="response")
    ## Train RIVER with all data for application
    emModelAll <- integratedEM(FeatAll, OutAll, logisticAllCV$lambda.min,
                logisticAllCV$glmnet.fit, theta_init,
                 num_bins,costs, pseudoc, verbose)


    ## Compute P(FR | G, E)
    postprobRIVER <- testPosteriors(FeatAll, OutAll, emModelAll)
    ## Output: postprobs
    postprobs<-list(
        indiv_name=unlist(strsplit(rownames(FeatAll),":"))
        [seq(1,2*nrow(FeatAll),by=2)],
        gene_name=unlist(strsplit(rownames(FeatAll),":"))
        [seq(2,2*nrow(FeatAll),by=2)],
        RIVER_posterior=postprobRIVER$posterior[,2],
        GAM_posterior=as.numeric(postporbGAM),
        fitRIVER=emModelAll)
    class(postprobs) <- "appl"

    data_framer = data.frame(neg_log10_median_pvalue = E_all, GAM_posterior = postprobs$GAM_posterior, posterior = postprobs$RIVER_posterior)
    scatter_plot_fill(data_framer, emModelAll, paste0(output_root,"scatter_fill_", ZscoreThrd,".png"))
}

max_ignore_na <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)


min_ignore_na <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)


discritize_expression_data <- function(E_tbt, dim, num_bins) {
  maxy <- max_ignore_na(E_tbt) + .00001  #add .00001 in order to create bins
  miny <- min_ignore_na(E_tbt)
  bin_vector <- seq(miny,maxy,length.out=num_bins+1)
  discretized_mat <- matrix(0,dim(E_tbt)[1],dim(E_tbt)[2])
  for (bin_number in 1:num_bins) {
    bin_start <- bin_vector[bin_number]
    bin_end <- bin_vector[bin_number+1]
    temp_matrix <- (E_tbt>=bin_start & E_tbt<bin_end)*bin_number
    discretized_mat <- discretized_mat + temp_matrix
  }
  return(discretized_mat) 
}


# Train on all non-N2 pairs
# Used trained model to do N2 pair prediction
roc_analysis_driver <- function(input_file, ZscoreThrd, output_root, dimensions, theta_init, num_bins, costs, verbose, pseudoc) {
  ## Extract required data
  all_data <- load_data(input_file, 3.3)

  # Load in expression data
  dataInput <- all_data[[1]]
  E_all <- all_data[[2]]
  E_tbt_real_valued <- t(as.matrix(all_data[[2]]))

  E_tbt <- discritize_expression_data(E_tbt_real_valued,dim,num_bins)

  # all genomic features (G)
  FeatAll <- t(exprs(dataInput))
  # all outlier status (E)
  OutAll <- as.numeric(unlist(dataInput$Outlier))-1

  # G for training models
  FeatTrng <- t(exprs(dataInput[,is.na(dataInput$N2pair)]))
  # E for training models
  OutTrng <- t(t(E_tbt[,is.na(dataInput$N2pair)]))
  OutTrng_binary <- as.numeric(unlist(dataInput$Outlier[is.na(dataInput$N2pair)]))-1

  # G for test
  FeatTest <-
    t(cbind(exprs(dataInput[,!is.na(dataInput$N2pair)])
    [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
    exprs(dataInput[,!is.na(dataInput$N2pair)])
    [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))

  # E for test (1st and then 2nd individuals from N2 pairs)



  OutTest1 <-
    as.numeric(unlist(
      c(E_tbt[!is.na(dataInput$N2pair)]
      [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
      E_tbt[!is.na(dataInput$N2pair)]
      [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)])))

  OutTest1 <- as.matrix(OutTest1)


  # E for test (2nd and then 1st individuals from N2 pairs)
  OutTest2 <-
    as.numeric(unlist(
      c(dataInput$Outlier[!is.na(dataInput$N2pair)]
      [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)],
      dataInput$Outlier[!is.na(dataInput$N2pair)]
      [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1

  ## Standardization
  meanFeat <- apply(FeatAll, 2, mean)
  sdFeat <- apply(FeatAll,2,sd)
  FeatAll <- scale(FeatAll, center=meanFeat, scale=sdFeat)
  FeatTrng <- scale(FeatTrng, center=meanFeat, scale=sdFeat)

  ## Search a best lambda from a multivariate logistic regression
  ##         with outlier status with 10 cross-validation
  ## GAM (genomeic annotation model)
  logisticCV <- cv.glmnet(FeatTrng, as.vector(OutTrng_binary), lambda=costs,
                          family="binomial", alpha=0, nfolds=10)
  if (verbose) {
    cat(' *** best lambda = ',logisticCV$lambda.min,' *** \n\n', sep='')
  }

  ## Compute a P(FR | G) for all data
  postprobTest <- predict(logisticCV, FeatTest, s="lambda.min", type="response")


  ## Train RIVER on training data
  emModel <- integratedEM(FeatTrng, OutTrng, logisticCV$lambda.min,
          logisticCV$glmnet.fit, theta_init,
          num_bins, costs,pseudoc, verbose)

  # ## Generate G data for test data (Revised)
  FeatTest <- scale(FeatTest, center=meanFeat, scale=sdFeat)

  ## Compute P(FR | G, E)
  dup.post <- testPosteriors(FeatTest, OutTest1, emModel)

  all_data <- load_data(input_file, ZscoreThrd)
  dataInput <- all_data[[1]]


  # E for test (2nd and then 1st individuals from N2 pairs)
  OutTest2 <-
    as.numeric(unlist(
      c(dataInput$Outlier[!is.na(dataInput$N2pair)]
      [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)],
      dataInput$Outlier[!is.na(dataInput$N2pair)]
      [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1




  ## Check performance of models with N2 pairs
  RIVER.roc <- roc(OutTest2, dup.post$posterior[,2]) # RIVER
  GAM.roc <- roc(OutTest2, as.numeric(postprobTest)) # GAM

  if (verbose) {
    cat('*** AUC (GAM - genomic annotation model): ',round(GAM.roc$auc,3),
      '\n    AUC (RIVER): ',round(RIVER.roc$auc,3),'\n     P-value: ',
      format.pval(roc.test(RIVER.roc, GAM.roc)$p.value,digits=2,eps=0.001),
      '***\n\n')
  }

  evaROC <-
    list(RIVER_sens=RIVER.roc$sensitivities,
         RIVER_spec=RIVER.roc$specificities,
         RIVER_auc=RIVER.roc$auc[1],
         GAM_sens=GAM.roc$sensitivities,
         GAM_spec=GAM.roc$specificities,
         GAM_auc=GAM.roc$auc[1],
         pvalue=roc.test(RIVER.roc, GAM.roc)$p.value)
  class(evaROC) <- "eval"

  plot_roc(evaROC, ZscoreThrd, paste0(output_root, "roc_curve_", ZscoreThrd,".pdf"))


}

# Initialize multinomial distribution
initialize_theta <- function(num_bins,dim) {
  theta_outlier <- matrix(1,dim,num_bins)
  theta_inlier <- matrix(1,dim,num_bins)
  theta_inlier[,1] = .4
  theta_inlier[,2] = .3
  theta_inlier[,3] = .1
  theta_inlier[,4] = .05
  theta_inlier[,5] = .05
  theta_inlier[,6] = .05
  theta_inlier[,7] = .05

  theta_outlier[,1] = .05
  theta_outlier[,2] = .05
  theta_outlier[,3] = .05
  theta_outlier[,4] = .05
  theta_outlier[,5] = .1
  theta_outlier[,6] = .3
  theta_outlier[,7] = .4
  theta_init <- list(inlier_component = theta_inlier, outlier_component = theta_outlier)

}


print("UNIVARIATE BINNED NAIVE BAYES")
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_root = args[2]


pseudoc=100  # hyperparameter for prior on multinomial distribution
dimensions=1 # number of tissues (just the median.. so only 1)

costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)  # Grid-space to search for lambda hyperparameters
verbose=TRUE
ZscoreThrd = 3.3 # Really a -log10pvalue threshold
num_bins <- 7  # Number of dimensions of multinomial distribution


# Initialize multinomial distribution
theta_init <- initialize_theta(num_bins,dimensions) 


# Train on all non-N2 pairs
# Used trained model to do N2 pair prediction
roc_analysis_driver(input_file, ZscoreThrd, output_root, dimensions, theta_init, num_bins, costs, verbose, pseudoc)





# Train on all data. Make nice visualization of posteriors
full_data_visualization_driver(input_file, ZscoreThrd, output_root, dimensions, theta_init,num_bins, costs, verbose,pseudoc)
