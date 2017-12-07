library(glmnet)
library(methods)
library(stats)
library(utils)
library(Biobase)
library(pROC)
library(ggplot2)
library(cowplot)
library(sigmoid)
library(Rcpp)
library(RColorBrewer)
library(ggthemes)
library(lbfgs)
library(optimx)
library(numDeriv)
sourceCpp("crf_vi_updates.cpp")





getFuncRvPosteriors <- function(Out,Feat, model_params) {
  observed_indices <- !is.na(Out)
  unobserved_indices <- is.na(Out)


  # Now compute posterior for every (tissue, sample) where we do have expression data

  probOut_FuncRvInlier <- matrix(0, model_params$N, model_params$T)
  probOut_FuncRvOutlier <- matrix(0, model_params$N, model_params$T)

  for (n in 1:model_params$N) {
    for (t in 1:model_params$T) {
      probOut_FuncRvInlier[n,t] <- model_params$phi$inlier_component[1,Out[n,t]]
      probOut_FuncRvOutlier[n,t] <- model_params$phi$outlier_component[1, Out[n,t]]
    }
  }


  probOut_FuncRvInlier[unobserved_indices] <- 1.0
  probOut_FuncRvOutlier[unobserved_indices] <- 1.0

  ll_plus <- log(probOut_FuncRvOutlier)
  ll_minus <- log(probOut_FuncRvInlier)

  vectorized_coefficients <- extract_crf_coefficient_vector(model_params$beta, model_params$theta_singleton, model_params$theta_pair, model_params$T)
  model_params$posterior <- variational_inference_posterior_updates(Feat, model_params$posterior, model_params$posterior_mu, vectorized_coefficients$singleton_vec, vectorized_coefficients$beta_vec, vectorized_coefficients$pair_vec, ll_plus, ll_minus, model_params$vi_damping_factor)
  model_params$posterior_mu <- 2.0*model_params$posterior - 1.0


  return(model_params)
}



mle_phi <- function(Out, model_params, pseudoc, num_bins) {
  dimension <- dim(Out)[2]
  phi_outlier <- matrix(1,dimension,num_bins)
  phi_inlier <- matrix(1,dimension,num_bins)
  for (bin_number in 1:num_bins) {
    phi_outlier[,bin_number] <- colSums(((Out==bin_number)*model_params$posterior),na.rm=TRUE)
    phi_inlier[,bin_number] <- colSums(((Out==bin_number)*(1-model_params$posterior)),na.rm=TRUE)
  }
  phi_outlier <- colSums(phi_outlier) + pseudoc
  phi_inlier <- colSums(phi_inlier) + pseudoc

  phi_outlier <- phi_outlier/sum(phi_outlier)
  phi_inlier <- phi_inlier/sum(phi_inlier)
  model_params$phi$inlier_component <- t(as.matrix(phi_inlier))
  model_params$phi$outlier_component <- t(as.matrix(phi_outlier))

  return(model_params)
}




compute_elbo <- function(mu, T, theta_singleton, theta_pair, G, beta) {
  elbo <- 0
  for (t in 1:T) {
     elbo <- elbo + theta_singleton[t]*mu[t] + sum(G*beta)*mu[t] - mu[t]*log(mu[t]) - (1-mu[t])*log(1-mu[t])
    for (t2 in (t+1):(T-1)) {
      if (t2 > T) {
        break
      }
      if (t2 > T) {
        break
      }
      elbo <- elbo + theta_pair[t,t2]*mu[t]*mu[t2]
    }
  }
  return(elbo)
}


mean_field_vi_update <- function(G, mu, T, theta_singleton, theta_pair, beta) {
  maxIter <- 1000  
  converged <- 0
  for (iter in 1:maxIter) {
    # print(compute_elbo(mu, T, theta_singleton, theta_pair, G, beta))
    mu_old <- mu
    for (t in 1:T) {
      mu[t] <- sigmoid( sum(mu * theta_pair[t,]) - mu[t]*theta_pair[t,t] + theta_singleton[t] + sum(G*beta) )
    }
    if (sum(abs(mu - mu_old)) < 5e-5) {
      converged <- 1
      break
    }

  }
  if (converged == 0) {
    print('ERROR!! LACK OF CONVERGENCE')
  }  
  return(mu)
}


getFuncRvFeat <- function(Feat, model_params) {
  vectorized_coefficients <- extract_crf_coefficient_vector(model_params$beta, model_params$theta_singleton, model_params$theta_pair, model_params$T)
  model_params$q <- variational_inference_mu_updates(Feat, model_params$q, model_params$mu, vectorized_coefficients$singleton_vec, vectorized_coefficients$beta_vec, vectorized_coefficients$pair_vec, model_params$vi_damping_factor)
  model_params$mu <- 2.0*model_params$q - 1.0
  return(model_params)
}

update_mu_sample <- function(G_sample, mu_sample, T, gradient_params) {
  N <- dim(G_sample)[1]
  theta_singleton <- as.matrix(gradient_params$singleton_vec)
  theta_pair <- matrix(1,T,T)
  counter <- 1
  for (t in 1:(T-1)) {
    for (t2 in (t+1):T) {
      theta_pair[t,t2] <- gradient_params$pair_vec[counter]
      theta_pair[t2,t] <- gradient_params$pair_vec[counter]
      counter <- counter + 1
    }
  }
  for (n in 1:N) {
    mu_sample[n,] <- mean_field_vi_update(G_sample[n,], mu_sample[n,], T, theta_singleton, theta_pair,gradient_params$beta_vec)
  }
  return(mu_sample)
}


testPosteriors <- function(Feat, Out, emModel) {
  model_params <- emModel$model_params
  model_params$N <- dim(Out)[1]
  model_params$mu <- matrix(0, model_params$N, model_params$T)
  model_params$q <- matrix(.5, model_params$N, model_params$T)
  model_params$posterior <- matrix(.5,model_params$N,model_params$T)
  model_params$posterior_mu <- matrix(0,model_params$N,model_params$T)
  model_params <- getFuncRvFeat(Feat, model_params)
  model_params <- getFuncRvPosteriors(Out, Feat, model_params)
  median_observed_posterior <- calculate_median_observed_posterior(Out, model_params$posterior)
  return(median_observed_posterior)
}

testPosteriors_tbt <- function(Feat, Out, emModel) {
  model_params <- emModel$model_params
  model_params$N <- dim(Out)[1]
  model_params$mu <- matrix(0, model_params$N, model_params$T)
  model_params$q <- matrix(.5, model_params$N, model_params$T)
  model_params$posterior <- matrix(.5,model_params$N,model_params$T)
  model_params$posterior_mu <- matrix(0,model_params$N,model_params$T)
  model_params <- getFuncRvFeat(Feat, model_params)
  model_params <- getFuncRvPosteriors(Out, Feat, model_params)
  return(model_params$posterior)
}

testPosteriors_exact <- function(Feat, Out, emModel) {
  model_params <- emModel$model_params
  model_params$N <- dim(Out)[1]
  model_params$mu <- matrix(.5, model_params$N, model_params$T)
  model_params$posterior <- matrix(.5,model_params$N,model_params$T)
  model_params <- get_marginal_exact_posteriors(Feat, Out, model_params)
  median_observed_posterior <- calculate_median_observed_posterior(Out, model_params$posterior)
  return(median_observed_posterior)
}




compute_log_probability <- function(Feat, Out, model_params) {
    probFuncRvFeat <- model_params$mu



    observed_indices <- !is.na(Out)
    unobserved_indices <- is.na(Out)

    probOut_FuncRvInlier <- matrix(0, model_params$N, model_params$T)
    probOut_FuncRvOutlier <- matrix(0, model_params$N, model_params$T)

    for (n in 1:model_params$N) {
      for (t in 1:model_params$T) {
        probOut_FuncRvInlier[n,t] <- model_params$phi$inlier_component[1,Out[n,t]]
        probOut_FuncRvOutlier[n,t] <- model_params$phi$outlier_component[1, Out[n,t]]
      }
    }


    inlier <- probOut_FuncRvInlier*(1-probFuncRvFeat)
    outlier <- probOut_FuncRvOutlier*probFuncRvFeat

    inlier[unobserved_indices] <- (1-probFuncRvFeat)[unobserved_indices]
    outlier[unobserved_indices] <- probFuncRvFeat[unobserved_indices]

    observed_data_log_likelihood <- sum(log(inlier + outlier))

    return(observed_data_log_likelihood)
}



temp_initialization <- function(Feat, Out, phi_init, beta_init, lambda,singleton_init, vi_damping_factor) {
   T <- dim(Out)[2]
   N <- dim(Out)[1]
   D <- dim(Feat)[2]

   model_params <- list(theta_singleton = matrix(singleton_init, T, 1), theta_pair = matrix(.01,T,T), q = matrix(.5,N,T),
                   mu = matrix(0, N, T), posterior=matrix(.5, N, T), posterior_mu=matrix(0, N, T), N = N, D = D, beta = beta_init, 
                   T = T, phi = phi_init, lambda = lambda,lambda_theta_pair=lambda, lambda_theta_singleton=0, vi_damping_factor=vi_damping_factor,
                   pseudo_q = matrix(0,N,T), pseudo_mu = matrix(0,N,T))
   return(model_params)
}


extract_crf_coefficient_vector <- function(beta, theta_singleton, theta_pair, T) {

  beta_num <- length(beta)
  beta_vec <- beta

  singleton_num <- T
  singleton_vec <- theta_singleton[,1]

  pair_num <- (T*T - T)/2
  pair_vec <- numeric(length=pair_num)

  counter <- 1
  for (t in 1:(T-1)) {
    for (t2 in (t+1):T) {
      pair_vec[counter] <- theta_pair[t,t2]
      counter <- counter + 1
    }
  }
  gradient_parameters <- list(beta_vec=beta_vec, singleton_vec=singleton_vec, pair_vec=pair_vec)
  return(gradient_parameters)
}

restore_crf_coefficient_vector <- function(model_params, gradient_params) {
  model_params$beta <- gradient_params$beta_vec

  model_params$theta_singleton[,1] = gradient_params$singleton_vec

  counter <- 1
  T <- model_params$T
  for (t in 1:(T-1)) {
    for (t2 in (t+1):T) {
      model_params$theta_pair[t,t2] <- gradient_params$pair_vec[counter]
      model_params$theta_pair[t2,t] <- gradient_params$pair_vec[counter]
      counter <- counter + 1
    }
  }

  return(model_params)
}

compute_singleton_gradient <- function(singleton_vec, posterior_sample, mu_sample, lambda, T, batch_size, N) {
  grad <- (colSums(posterior_sample) - colSums(mu_sample))*(1/batch_size) - 2*lambda*singleton_vec
  return(grad)
}

compute_exact_singleton_gradient <- function(singleton_vec, posterior_sample, mu_sample, lambda, T, batch_size, N) {
  grad <- (colSums(posterior_sample) - colSums(mu_sample))*(1/batch_size) - 2*lambda*singleton_vec
  return(grad)
}

compute_beta_gradient <- function(beta, posterior_sample, mu_sample, G_sample, lambda, T, batch_size, N) {
  grad <- numeric(length = length(beta))
  for (t in 1:T) {
    grad <- grad + colSums(posterior_sample[,t]*G_sample) - colSums(mu_sample[,t]*G_sample)
  }
  regularizer_term <- lambda*beta
  # regularizer_term[1] <- 0  # Intercept should not be penalized
  grad <- (grad)*(1/(batch_size)) - regularizer_term
  return(grad)
}

compute_pair_gradient <- function(pair_vec, posterior_sample, mu_sample, lambda, T, batch_size, N) {
  cross_posterior <- t(posterior_sample) %*% posterior_sample
  cross_mu <- t(mu_sample) %*% mu_sample

  grad <- numeric(length(pair_vec))
  num_samples <- dim(posterior_sample)[1]
  counter <- 1
  for (t in 1:(T-1)) {
    for (t2 in (t+1):T) {
      grad[counter] <- (cross_posterior[t,t2] - cross_mu[t,t2])
      counter <- counter + 1
    }
  }
  grad <- grad*(1/batch_size) - lambda*pair_vec
  return(grad)
}



compute_crf_log_likelihood <- function(singleton_vec, beta, pair_vec, lambda, posterior_mu, mu, T, G, N, batch_size, lambda_theta_singleton,lambda_theta_pair) {
  q <- (mu + 1.0)/2.0  

  singleton_term <- sum(posterior_mu %*% as.matrix(singleton_vec)) - sum(mu %*% as.matrix(singleton_vec))

  cross_posterior <- t(posterior_mu) %*% posterior_mu
  cross_mu <- t(mu) %*% mu 
  pair_term <- 0
  counter <- 1
  for (t in 1:(T-1)) {
    for (t2 in (t+1):T) {
      pair_term <- pair_term + pair_vec[counter]*cross_posterior[t,t2] - pair_vec[counter]*cross_mu[t,t2]
      counter <- counter + 1
    }
  }

  beta_G <- (G %*% as.matrix(beta))[,1]
  beta_term <- sum(posterior_mu*beta_G) - sum(mu*beta_G)

  non_differentiable_term_1 <- q*log(q)
  non_differentiable_term_1[is.na(non_differentiable_term_1)] <- 0

  non_differentiable_term_2 <- (1-q)*log(1-q)
  non_differentiable_term_2[is.na(non_differentiable_term_2)] <- 0


  non_differentiable_term <- sum(non_differentiable_term_1) + sum(non_differentiable_term_2)
  regularization_term <- (lambda/2)*sum(beta*beta) + (lambda_theta_pair/2)*sum(pair_vec*pair_vec) + (lambda_theta_singleton/2)*sum(singleton_vec*singleton_vec)
  log_likelihood <- singleton_term + pair_term + beta_term - regularization_term*batch_size + non_differentiable_term
  return(log_likelihood)
}



gradient_ascent <- function(gradient_params, model_params, Feat, eta) {
  log_likelihood <- -100000000000000000000
  converged <- 0
  G_sample <- Feat
  posterior_mu_sample <- model_params$posterior_mu
  mu_sample <- model_params$mu
  q_sample <- model_params$q
  batch_size <- dim(Feat)[1]
  counter <- 0
  while (converged == 0) {
    # Randomly sample 1 batch
    prev_log_likelihood <- log_likelihood
    # mu_sample <- update_mu_sample(G_sample, mu_sample, model_params$T, gradient_params)

    q_sample <- variational_inference_mu_updates(G_sample, q_sample, mu_sample, gradient_params$singleton_vec, gradient_params$beta_vec, gradient_params$pair_vec)
    mu_sample <- 2.0*q_sample - 1.0
    
    log_likelihood <- compute_crf_log_likelihood(gradient_params$singleton_vec, gradient_params$beta_vec, gradient_params$pair_vec, model_params$lambda, posterior_mu_sample, mu_sample, model_params$T, G_sample, model_params$N, batch_size,model_params$lambda_theta_singleton,model_params$lambda_theta_pair)



    #grad_singleton <- compute_singleton_gradient(gradient_params$singleton_vec, posterior_mu_sample, mu_sample, model_params$lambda_theta_singleton, model_params$T, batch_size, model_params$N)
    grad_singleton <- compute_singleton_gradient_cpp(gradient_params$singleton_vec, posterior_mu_sample, mu_sample)
    gradient_params$singleton_vec <- gradient_params$singleton_vec + eta*grad_singleton
    #grad_beta <- compute_beta_gradient(gradient_params$beta_vec, posterior_mu_sample, mu_sample, G_sample, model_params$lambda, model_params$T, batch_size, model_params$N)
    grad_beta <- compute_beta_gradient_cpp(gradient_params$beta_vec, posterior_mu_sample, mu_sample, G_sample, model_params$lambda)
    gradient_params$beta_vec <- gradient_params$beta_vec + eta*grad_beta
    # grad_pair_old <- compute_pair_gradient(gradient_params$pair_vec, posterior_mu_sample, mu_sample, model_params$lambda_theta_pair, model_params$T,batch_size, model_params$N)
    grad_pair <- compute_pair_gradient_cpp(gradient_params$pair_vec, posterior_mu_sample, mu_sample, model_params$lambda_theta_pair)
    gradient_params$pair_vec <- gradient_params$pair_vec + eta*grad_pair

    if (counter == 0) {
      print("Hello")
      print(log_likelihood)
      print(grad_singleton)
      print(grad_beta)
      print(grad_pair)
    }  


    if (counter %% 10 == 0) {
      log_likelihood <- compute_crf_log_likelihood(gradient_params$singleton_vec, gradient_params$beta_vec, gradient_params$pair_vec, model_params$lambda, posterior_mu_sample, mu_sample, model_params$T, G_sample, model_params$N, batch_size,model_params$lambda_theta_singleton,model_params$lambda_theta_pair)
      print(log_likelihood)

      if (abs(prev_log_likelihood - log_likelihood) < .1) {
        converged <- 1
        break
       }
        prev_log_likelihood <- log_likelihood
    }
    counter <- counter + 1
    # print(paste(iter,abs(sum(grad_singleton)),abs(sum(grad_beta)), abs(sum(grad_pair))))

  }
  if (converged == 0) {
    print('CONVERGENCE ERROR (gradient ascent)')
  }
  return(gradient_params)
}


mle_beta <- function(model_params, Feat) {
  eta <- .0008
  #eta <- .001
  batch_size <- 1000
  gradient_params <- extract_crf_coefficient_vector(model_params$beta, model_params$theta_singleton, model_params$theta_pair, model_params$T)
  #gradient_params <- mini_batch_gradient_ascent(gradient_params, model_params, Feat, eta, batch_size)
  gradient_params <- gradient_ascent(gradient_params, model_params, Feat, eta)
  model_params <- restore_crf_coefficient_vector(model_params, gradient_params)
  return(model_params)
}

gradient_lbfgs <- function(x, Feat, posterior, posterior_mu, q, mu, pair_length, num_feat, T, lambda, lambda_theta_pair, N, vi_damping_factor) {
  gradient <- numeric(length(x))
  intercept <- x[1:T]
  pair_vec <- x[(T+1):(T+pair_length)]
  beta <- x[(T+pair_length+1):length(x)]

  q <- pseudo_mu_updates(Feat, q, mu, intercept, beta, pair_vec,posterior_mu)
  mu <- 2.0*q - 1.0

  #gradient[1:T] <- compute_singleton_gradient_cpp_2(intercept, posterior_mu, mu)
  gradient[1:T] <- compute_singleton_gradient_pseudo_cpp(intercept, posterior_mu, mu)

  #gradient[(T+1):(T+pair_length)] <- compute_pair_gradient_cpp_2(pair_vec, posterior_mu, mu, lambda_theta_pair)
  gradient[(T+1):(T+pair_length)] <- compute_pair_gradient_pseudo_cpp(pair_vec, posterior_mu, mu, lambda_theta_pair)

  #gradient[(T+pair_length+1):length(x)] <- compute_beta_gradient_cpp_2(beta, posterior_mu, mu, Feat, lambda)
  gradient[(T+pair_length+1):length(x)] <- compute_beta_gradient_pseudo_cpp(beta, posterior_mu, mu, Feat, lambda)

  return(-gradient)
}


compute_crf_log_likelihood_lbfgs <- function(x, Feat, posterior, posterior_mu, q, mu, pair_length, num_feat, T, lambda, lambda_theta_pair, N, vi_damping_factor) {
  intercept <- x[1:T]
  pair_vec <- x[(T+1):(T+pair_length)]
  beta <- x[(T+pair_length+1):length(x)]

  q <- pseudo_mu_updates(Feat, q, mu, intercept, beta, pair_vec,posterior_mu)
  mu <- 2.0*q - 1.0

  log_likelihood <- crf_log_likelihood_helper_pseudo_cpp(intercept, pair_vec, beta, mu, q, posterior_mu,posterior, Feat, T, lambda, lambda_theta_pair)
  print(round(log_likelihood,digits=22))
  return(-log_likelihood)
}


mle_beta_lbfgs <- function(model_params, Feat) {
  gradient_params <- extract_crf_coefficient_vector(model_params$beta, model_params$theta_singleton, model_params$theta_pair, model_params$T)
  
  x <- c(gradient_params$singleton_vec, gradient_params$pair_vec, gradient_params$beta_vec)
  
  T <- model_params$T
  N <- dim(Feat)[1]
  pair_length <- (T*T - T)/2.0
  num_feat <- dim(Feat)[2]

  #grad <- gradient_lbfgs(x, Feat, model_params$posterior, model_params$posterior_mu, model_params$pseudo_q, model_params$pseudo_mu, pair_length, num_feat, T, model_params$lambda, model_params$lambda_theta_pair, N, model_params$vi_damping_factor)
  #ll <- compute_crf_log_likelihood_lbfgs(x, Feat, model_params$posterior, model_params$posterior_mu, model_params$pseudo_q, model_params$pseudo_mu, pair_length, num_feat, T, model_params$lambda, model_params$lambda_theta_pair, N, model_params$vi_damping_factor)
  #numer_grad <- grad(compute_crf_log_likelihood_lbfgs,x,Feat = Feat, posterior = model_params$posterior, posterior_mu = model_params$posterior_mu, q = model_params$pseudo_q, mu = model_params$pseudo_mu, pair_length = pair_length, num_feat = num_feat, T =T, lambda=model_params$lambda, lambda_theta_pair = model_params$lambda_theta_pair, N = N, vi_damping_factor = model_params$vi_damping_factor)
  #print(numer_grad)
  #print(grad)
  #print(ll)
  output <- lbfgs(compute_crf_log_likelihood_lbfgs, gradient_lbfgs, x, Feat=Feat, posterior=model_params$posterior,posterior_mu=model_params$posterior_mu,q=model_params$pseudo_q, mu=model_params$pseudo_mu, pair_length = pair_length, num_feat=num_feat, T = T, lambda=model_params$lambda, lambda_theta_pair=model_params$lambda_theta_pair, N = N, vi_damping_factor=model_params$vi_damping_factor, invisible=0,epsilon=5e-2,m=10)
  new_x <- output$par
  if (output$convergence != 0) {
    print("ERRROOROR IN LBFGS!!!")
  }



  gradient_params$singleton_vec <- new_x[1:T]
  gradient_params$pair_vec <- new_x[(T+1):(T+pair_length)]
  gradient_params$beta_vec <- new_x[(T+pair_length+1):length(new_x)]

  model_params <- restore_crf_coefficient_vector(model_params, gradient_params)
  return(model_params)
}

calculate_median_observed_posterior <- function(Out, posterior) {
  temp_Out <- Out
  temp_Out[!is.na(temp_Out)] = 1
  new_posterior <- temp_Out*posterior
  median_observed_posterior <- apply(new_posterior,1,median,na.rm=TRUE)
  return(median_observed_posterior)
}

compute_expected_complete_log_likelihood <- function(Feat, Out, model_params) {
  vectorized_coefficients <- extract_crf_coefficient_vector(model_params$beta, model_params$theta_singleton, model_params$theta_pair, model_params$T)
  model_params$q <- variational_inference_mu_updates(Feat, model_params$q, model_params$mu, vectorized_coefficients$singleton_vec, vectorized_coefficients$beta_vec, vectorized_coefficients$pair_vec, model_params$vi_damping_factor)

  model_params$mu <- 2.0*model_params$q - 1.0

  observed_indices <- !is.na(Out)
  unobserved_indices <- is.na(Out)

  probOut_FuncRvInlier <- matrix(0, model_params$N, model_params$T)
  probOut_FuncRvOutlier <- matrix(0, model_params$N, model_params$T)

  for (n in 1:model_params$N) {
    for (t in 1:model_params$T) {
      probOut_FuncRvInlier[n,t] <- model_params$phi$inlier_component[1,Out[n,t]]
      probOut_FuncRvOutlier[n,t] <- model_params$phi$outlier_component[1, Out[n,t]]
    }
  }

  probOut_FuncRvInlier[unobserved_indices] = 1
  probOut_FuncRvOutlier[unobserved_indices] = 1

  inlier_term <- (1.0-model_params$posterior)*log(probOut_FuncRvInlier*(1.0-model_params$q))
  outlier_term <- (model_params$posterior)*log(probOut_FuncRvOutlier*model_params$q)

  return(sum(inlier_term) + sum(outlier_term))

}


compute_observed_data_log_like <- function(Feat, Out, model_params_exact) {
  gradient_params <- extract_crf_coefficient_vector(model_params_exact$beta, model_params_exact$theta_singleton, model_params_exact$theta_pair, model_params_exact$T)
  phi_temp <- rbind(model_params_exact$phi$inlier_component,model_params_exact$phi$outlier_component)
  log_like <- 0
  for (n in 1:model_params_exact$N) {
    normalization_constant <- compute_mu_normalization_constant(Feat[n,], gradient_params$singleton_vec, gradient_params$beta_vec, gradient_params$pair_vec)
    summer <- 0
    for (z1 in 0:1) {
      for (z2 in 0:1) {
        for (z3 in 0:1) {
          crf_prob <- exact_mu_prob(z1, z2, z3,normalization_constant, Feat[n,], gradient_params$singleton_vec, gradient_params$beta_vec, gradient_params$pair_vec)
          joint_prob <- crf_prob*phi_temp[(z1+1),Out[n,1]]*phi_temp[(z2+1),Out[n,2]]*phi_temp[(z3+1),Out[n,3]]
          summer <- summer + joint_prob
        }
      }
    }
    log_like <- log_like + log(summer)
  }


  return(log_like)
}


compute_observed_data_log_like3 <- function(Feat, Out, model_params_exact) {
  vectorized_coefficients <- extract_crf_coefficient_vector(model_params_exact$beta, model_params_exact$theta_singleton, model_params_exact$theta_pair, model_params_exact$T)
  phi_temp <- rbind(model_params_exact$phi$inlier_component,model_params_exact$phi$outlier_component)
  model_params_exact$q <- variational_inference_mu_updates(Feat, model_params_exact$q, model_params_exact$mu, vectorized_coefficients$singleton_vec, vectorized_coefficients$beta_vec, vectorized_coefficients$pair_vec)
  model_params_exact$mu <- 2.0*model_params_exact$q - 1.0
  log_like <- 0
  for (n in 1:model_params_exact$N) {
    summer <- 0
    for (z1 in 0:1) {
      for (z2 in 0:1) {
        for (z3 in 0:1) {
          crf_prob <- (z1*model_params_exact$q[n,1] + (1-z1)*(1-model_params_exact$q[n,1]))*(z2*model_params_exact$q[n,2] + (1-z2)*(1-model_params_exact$q[n,2]))*(z3*model_params_exact$q[n,3] + (1-z3)*(1-model_params_exact$q[n,3]))
          if (is.na(Out[n,1]) == FALSE) {
            crf_prob <- crf_prob*phi_temp[(z1+1),Out[n,1]]
          }
          if (is.na(Out[n,2]) == FALSE) {
            crf_prob <- crf_prob*phi_temp[(z2+1),Out[n,2]]
          }
          if (is.na(Out[n,3]) == FALSE) {
            crf_prob <- crf_prob*phi_temp[(z3+1),Out[n,3]]
          }
          summer <- summer + crf_prob
        }
      }
    }
    log_like <- log_like + log(summer)
  }


  return(log_like)
}




integratedEM <- function(Feat, Out, lambda_ideal,
                         phi_init, beta_init, singleton_init, num_bins, costs,pseudoc,
                         verbose, OutTest1, OutTest2_tbt, FeatTest, Zscore_initialization_threshold, Zscore_n2_threshold, output_root, vi_damping_factor){


  model_params <- temp_initialization(Feat, Out, phi_init, beta_init, lambda_ideal, singleton_init, vi_damping_factor)
  prev_expected_log_like <- -1000000000000000000
  prev_observed_log_like <- -1000000000000000000
  steps <- 1
  maxIter <- 1000  
  converged <- 0
  for (iter in 1:maxIter) {
    if (verbose) {
      cat(' *** STREAM: EM step ',steps,'\n',sep="")
    }




    # model_params_exact <- get_marginal_exact_posteriors(Feat, Out, model_params_exact)

    ## E-step:
    ## Compute expected posterior probabilities
    ##           given current parameters and data

    model_params <- getFuncRvFeat(Feat, model_params)
    model_params <- getFuncRvPosteriors(Out, Feat, model_params)

    #expected_data_log_like <- compute_expected_complete_log_likelihood(Feat, Out, model_params)
    #cat('    Current expected data log probability after E step: ', expected_data_log_like,'\n',sep='')
    #observed_data_log_like <- compute_observed_data_log_like3(Feat, Out, model_params)
    #cat('    Current observed data log probability after E Step: ', observed_data_log_like,'\n',sep='')
    expected_data_log_like <- compute_expected_complete_log_likelihood(Feat, Out, model_params)
    cat('    Current expected data log probability after E Step: ', expected_data_log_like,'\n',sep='')


    if (verbose) {
      cat('E-step: complete', '\n',sep='')
    }
    ## M-step:

    #  Keep track of previous iteration's parameters in order to check for convergence
    phi_old <- model_params$phi
    beta_old <- model_params$beta
    theta_singleton_old <- model_params$theta_singleton
    theta_pair_old <- model_params$theta_pair

    # Maximum Likelihood Estimate (really Map...) of Phi
    model_params <- mle_phi(Out, model_params, pseudoc, num_bins)


    # Maximum Likelihood estimate (really MAP...) of Beta and Theta
    model_params <- mle_beta_lbfgs(model_params, Feat)

    #model_params <- mle_beta(model_params, Feat)
    #model_params <- mle_beta_lbfgs(model_params, Feat)


    #model_params_exact <- mle_phi(Out, model_params_exact,pseudoc,num_bins)

    #model_params_exact <- mle_beta_exact(model_params_exact, Feat)


    # Compute observed log probability


    print(head(model_params$theta_pair))
    print(model_params$theta_singleton)
    print(model_params$beta)
    print(model_params$phi)


    # Print convergence info
    if (verbose) {
      cat('     M-step: norm(phi_inlier difference) = ',
          round(norm(matrix(model_params$phi$inlier_component)-matrix(phi_old$inlier_component)),4),
          ', norm(phi_outlier difference) = ',
          round(norm(matrix(model_params$phi$outlier_component)-matrix(phi_old$outlier_component)),4),
          ', norm(beta difference) = ',
          round(norm(matrix(model_params$beta)-matrix(beta_old)),4),
          ', norm(theta singleton difference) = ',
          round(norm(matrix(model_params$theta_singleton)-matrix(theta_singleton_old)),4),
          ', norm(theta pairs difference) = ',
          round(norm(matrix(model_params$theta_pair)-matrix(theta_pair_old)),4),
          " *** \n\n", sep="")
    }

    #expected_data_log_like <- compute_expected_complete_log_likelihood(Feat, Out, model_params)
    #cat('    Current expected data log probability after M step: ', expected_data_log_like,'\n',sep='')
    #observed_data_log_like <- compute_observed_data_log_like3(Feat, Out, model_params)
    #cat('    Current observed data log probability after M Step: ', observed_data_log_like,'\n',sep='')
    expected_data_log_like <- compute_expected_complete_log_likelihood(Feat, Out, model_params)
    cat('    Current expected data log probability after M Step: ', expected_data_log_like,'\n',sep='')


    ## Check convergence
    if ((norm(matrix(model_params$beta) - matrix(beta_old)) < .05 ) &
        (norm(model_params$phi$inlier_component - phi_old$inlier_component) < .05) &
        (norm(model_params$phi$outlier_component - phi_old$outlier_component) < .05) &
        (norm(matrix(model_params$theta_singleton) - matrix(theta_singleton_old)) < .05) &
        (norm(matrix(model_params$theta_pair) - matrix(theta_pair_old)) < .5)
        ) {
      converged <- 1
      break
    }

    if (abs(prev_expected_log_like - expected_data_log_like) < .5) {
      converged <- 1
      break
    }
    prev_expected_log_like <- expected_data_log_like
    #if (abs(prev_observed_log_like - observed_data_log_like) < .5) {
    #  converged <- 1
    #  break
    #}
    #prev_observed_log_like <- observed_data_log_like
    dup.post <- testPosteriors_tbt(FeatTest, OutTest1, list(model_params=model_params))
    RIVER.roc_observed <- tbt_roc(dup.post, OutTest1, OutTest2_tbt, "observed_only")
    print(RIVER.roc_observed)

    RIVER.roc_inferred <- tbt_roc(dup.post, OutTest1, OutTest2_tbt, "inferred")

    print(RIVER.roc_inferred)
    write.table(dup.post,paste0(output_root, "test_posterior_",Zscore_initialization_threshold,"_",Zscore_n2_threshold, "_",pseudoc,"_",steps,".txt"),quote=FALSE,sep="\t")
    write.table(model_params$theta_pair,paste0(output_root, "theta_pair_",Zscore_initialization_threshold,"_",Zscore_n2_threshold, "_",pseudoc,"_",steps,".txt"),quote=FALSE,sep="\t")

     evaROC <-
      list(RIVER_sens=RIVER.roc_observed$sensitivities,
          RIVER_spec=RIVER.roc_observed$specificities,
          RIVER_auc=RIVER.roc_observed$auc[1],
          GAM_sens=RIVER.roc_inferred$sensitivities,
          GAM_spec=RIVER.roc_inferred$specificities,
          GAM_auc=RIVER.roc_inferred$auc[1],
          pvalue=roc.test(RIVER.roc_observed, RIVER.roc_inferred)$p.value)
    class(evaROC) <- "eval"
    plot_roc(evaROC, Zscore_n2_threshold, paste0(output_root, "roc_curve_",Zscore_initialization_threshold,"_",Zscore_n2_threshold, "_",pseudoc,"_",steps,".pdf"))



    steps <- steps + 1
  }

  if (converged == 1) {
    cat(" ::: EM iteration is terminated since it converges within a
        predefined tolerance (0.001) ::: \n\n\n",sep="")
  } else if ((converged == 0) && (iter == maxIter)) {
    cat(" ::: EM iteration is terminated since it reaches a
        predefined maximum value (1000) ::: \n\n\n",sep="")
  }

  median_observed_posterior <- calculate_median_observed_posterior(Out, model_params$posterior)
  median_posterior <- apply(model_params$posterior,1,median)
  list(model_params=model_params,
      posteriors=median_observed_posterior)
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
  legend(0.7,0.2,c("w-shed","GAM"), lty=c(1,1), lwd=c(2,2),
       col=c("dodgerblue","mediumpurple"), cex=1.2, 
       pt.cex=1.2, bty="n")
  title(main=paste("Threshold = ", ZscoreThrd," / AUC: Watershed = ", round(evaROC$RIVER_auc,3), 
                 ", GAM = ", round(evaROC$GAM_auc,3),sep=""))

  dev.off()
}


tbt_roc <- function(dup.post, OutTest1, OutTest2_tbt, version) {

  if (version == "observed_only") {
    valid_instances <- !is.na(OutTest1*OutTest2_tbt)
  }
  if (version == "inferred") {
    valid_instances <- !is.na(OutTest2_tbt)
  }

  posterior_vec <- dup.post[valid_instances]
  pseudo_gold_standard <- OutTest2_tbt[valid_instances]

  print(sum(pseudo_gold_standard))
  print(length(pseudo_gold_standard))


  return(roc(pseudo_gold_standard, posterior_vec))

}


scatter_plot_fill <- function(data_framer, emModel, figure_output_file) {



  scatter <- ggplot(data = data_framer, mapping = aes(x = neg_log10_median_pvalue, y = GAM_posterior)) + geom_point(aes(colour = posterior), shape = 18,size=2) +
                scale_colour_gradient2(low = "khaki1",midpoint=.5, mid="dodgerblue",high="mediumpurple", guide = "colourbar")
  scatter <- scatter + theme(axis.text.x = element_text(size=20,hjust=.5,vjust=.5,face="plain"),axis.text.y = element_text(size=20,hjust=.5,vjust=.5,face="plain"))
  scatter <- scatter + labs(x = "median(-log10(pvalue))", y = "GAM Posterior", colour="Median observed posterior") + theme(axis.title=element_text(size=20,face="bold"))

  ggsave(scatter, file=figure_output_file,width = 30,height=15,units="cm")


}



energy_function_visualization <- function(model_params, figure_output_file) {
  theta_pair <- model_params$theta_pair[1,2]
  theta1 <- model_params$theta_singleton[1,1]
  theta2 <- model_params$theta_singleton[2,1]

  x_range <- seq(0,1,.01)
  y_range <- seq(0,1,.01)

  num_entries <- length(x_range)*length(y_range)

  x <- numeric(num_entries)
  y <- numeric(num_entries)
  e <- numeric(num_entries)

  counter <- 1
  for (x_iter in 1:length(x_range)) {
    x_val <- 2.0*x_range[x_iter] - 1.0 # x_val <- 2.0*x_range[x_iter] - 1.0
    for (y_iter in 1:length(y_range)) {
      y_val <- 2.0*y_range[y_iter] - 1.0 # y_val <- 2.0*y_range[y_iter] - 1.0



      energy <- theta1*(x_val) + theta2*(y_val) + theta_pair*x_val*y_val
      x[counter] <- x_range[x_iter]
      y[counter] <- y_range[y_iter]
      e[counter] <- energy
      counter <- counter + 1
    }
  }
  df <- data.frame(x=x, y=y, energy=e)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

  energy_map <- ggplot(df, aes(x,y, fill = energy)) + geom_raster() + 
                coord_fixed(ratio = 1)  + 
                scale_fill_gradientn(colours = myPalette(4)) + 
                labs(x = "Z1") + 
                labs(y = "Z2") +
                labs(title = "Energy Function")

  ggsave(energy_map, file=figure_output_file)      
}

posterior_probability_historgram <- function(array, output_file) {
  df <- data.frame(prob = array)
  histo <-ggplot(df, aes(x=prob)) + 
          geom_histogram(color="black", fill="white")

  ggsave(histo,file=output_file)
}

full_data_visualization_driver <- function(input_file, ZscoreThrd, output_root, dimensions, phi_init,num_bins, costs, verbose,pseudoc, initialization) {
    ## Extract required data
    ## Extract required data
    all_data <- load_data(input_file, ZscoreThrd)
    dataInput <- all_data[[1]]
    E_all <- all_data[[2]]

    E_tbt_real_valued <- t(as.matrix(all_data[[3]]))

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

    # Under smart initializtion assumptions
    beta_init <- logisticAllCV$glmnet.fit$beta[,logisticAllCV$lambda == logisticAllCV$lambda.min]
    singleton_init <- logisticAllCV$glmnet.fit$a0[logisticAllCV$lambda == logisticAllCV$lambda.min]


    ## Compute P(FR=1 | G)
    postporbGAM <- predict(logisticAllCV, FeatAll, s="lambda.min", type="response")
    ## Train RIVER with all data for application

    # If not performing "smart initializtion", set all parameters to zeros
    if (initialization == "zeros") {
      beta_init <- numeric(length(beta_init))
      singleton_init <- 0
    }



    emModelAll <- integratedEM(FeatAll, OutAll, logisticAllCV$lambda.min,
              phi_init, beta_init, singleton_init, num_bins, costs, pseudoc, verbose)


    ## Compute P(FR | G, E)
    postprobRIVER <- testPosteriors(FeatAll, OutAll, emModelAll)

    ## Output: postprobs
    postprobs<-list(
        indiv_name=unlist(strsplit(rownames(FeatAll),":"))
        [seq(1,2*nrow(FeatAll),by=2)],
        gene_name=unlist(strsplit(rownames(FeatAll),":"))
        [seq(2,2*nrow(FeatAll),by=2)],
        RIVER_posterior=postprobRIVER,
        GAM_posterior=as.numeric(postporbGAM),
        fitRIVER=emModelAll)
    class(postprobs) <- "appl"


    # Compare relation between posterior probability and number of observed tissues
    num_observed_tissues <- apply(OutAll, 1, function(x) length(which(!is.na(x))))
    temp_matrix <- cbind(postprobs$RIVER_posterior, postprobs$GAM_posterior, E_all, num_observed_tissues)
    write.table(temp_matrix, file=paste0(output_root,"summarize_results.txt"), quote=FALSE, sep="\t", col.names=FALSE,row.names=FALSE)

    write.table(emModelAll$model_params$theta_pair,file=paste0(output_root,"theta_pair_table.txt"), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

    posterior_probability_historgram(as.vector(emModelAll$model_params$posterior), paste0(output_root,"posterior_histogram_", ZscoreThrd,"_",pseudoc,".png"))

    posterior_probability_historgram(as.vector(emModelAll$model_params$q), paste0(output_root,"crf_posterior_histogram_", ZscoreThrd,"_",pseudoc,".png"))

    energy_function_visualization(emModelAll$model_params, paste0(output_root,"energy_viz_", ZscoreThrd,"_",pseudoc,".png"))

    data_framer = data.frame(neg_log10_median_pvalue = E_all, GAM_posterior = postprobs$GAM_posterior, posterior = postprobs$RIVER_posterior)
    scatter_plot_fill(data_framer, emModelAll, paste0(output_root,"scatter_fill_", ZscoreThrd,"_",pseudoc,".png"))
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



roc_analysis_driver <- function(input_file, Zscore_initialization_threshold, Zscore_n2_threshold, output_root, dimensions, phi_init, num_bins, costs, verbose, pseudoc,initialization, vi_damping_factor) {
  ## Extract required data
  all_data <- load_data(input_file, Zscore_initialization_threshold)
  dataInput <- all_data[[1]]
  E_all <- all_data[[2]]
  E_tbt_real_valued <- t(as.matrix(all_data[[3]]))

  E_tbt <- discritize_expression_data(E_tbt_real_valued,dim,num_bins)

  # all genomic features (G)
  FeatAll <- t(exprs(dataInput))
  # all outlier status (E)
  OutAll <- as.numeric(unlist(dataInput$Outlier))-1

  # G for training models
  FeatTrng <- t(exprs(dataInput[,is.na(dataInput$N2pair)]))
  # E for training models
  OutTrng <- t(E_tbt[,is.na(dataInput$N2pair)])
  OutTrng_binary <- as.numeric(unlist(dataInput$Outlier[is.na(dataInput$N2pair)]))-1

  # G for test
  FeatTest <-
    t(cbind(exprs(dataInput[,!is.na(dataInput$N2pair)])
    [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
    exprs(dataInput[,!is.na(dataInput$N2pair)])
    [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))

  # E for test (1st and then 2nd individuals from N2 pairs)
  OutTest1 <-
      t(cbind(E_tbt[,!is.na(dataInput$N2pair)]
      [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
      E_tbt[,!is.na(dataInput$N2pair)]
      [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))

  # E for test (2nd and then 1st individuals from N2 pairs)
  all_data_n2 <- load_data(input_file, Zscore_n2_threshold)
  dataInput_n2 <- all_data_n2[[1]]
  OutTest2 <-
    as.numeric(unlist(
      c(dataInput_n2$Outlier[!is.na(dataInput_n2$N2pair)]
      [seq(from=2,to=sum(!is.na(dataInput_n2$N2pair)),by=2)],
      dataInput_n2$Outlier[!is.na(dataInput_n2$N2pair)]
      [seq(from=1,to=sum(!is.na(dataInput_n2$N2pair)),by=2)])))-1


  OutTest2_tbt_real_valued <- 
      t(cbind(E_tbt_real_valued[,!is.na(dataInput$N2pair)]
      [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)],
      E_tbt_real_valued[,!is.na(dataInput$N2pair)]
      [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)]))

  OutTest2_tbt <- (OutTest2_tbt_real_valued > Zscore_n2_threshold)*1.0
  print(head(OutTest2_tbt))
  OutTest2_tbt_real_valued_cp <- OutTest2_tbt_real_valued
  OutTest2_tbt_real_valued_cp[!is.nan(OutTest2_tbt_real_valued_cp)] = 1
  OutTest2_tbt <- OutTest2_tbt*OutTest2_tbt_real_valued_cp
  print(head(OutTest2_tbt))



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
  beta_init <- logisticCV$glmnet.fit$beta[,logisticCV$lambda == logisticCV$lambda.min]
  singleton_init <- logisticCV$glmnet.fit$a0[logisticCV$lambda == logisticCV$lambda.min]
  if (verbose) {
    cat(' *** best lambda = ',logisticCV$lambda.min,' *** \n\n', sep='')
  }

  ## Compute a P(FR | G) for all data
  postprobTest <- predict(logisticCV, FeatTest, s="lambda.min", type="response")

      # If not performing "smart initializtion", set all parameters to zeros
  if (initialization == "zeros") {
      beta_init <- numeric(length(beta_init))
      singleton_init <- 0
  }


  ## Train RIVER on training data
  emModel <- integratedEM(FeatTrng, OutTrng, logisticCV$lambda.min,
              phi_init, beta_init, singleton_init, num_bins, costs, pseudoc, verbose, OutTest1, OutTest2_tbt, FeatTest, Zscore_initialization_threshold, Zscore_n2_threshold, output_root, vi_damping_factor)





  # ## Generate G data for test data (Revised)
  FeatTest <- scale(FeatTest, center=meanFeat, scale=sdFeat)

  ## Compute P(FR | G, E)
  dup.post <- testPosteriors_tbt(FeatTest, OutTest1, emModel)


  RIVER.roc_observed <- tbt_roc(dup.post, OutTest1, OutTest2_tbt, "observed_only")
  RIVER.roc_inferred <- tbt_roc(dup.post, OutTest1, OutTest2_tbt, "inferred")

  print(RIVER.roc_observed)
  print(RIVER.roc_inferred)

  ## Check performance of models with N2 pairs
  # RIVER.roc <- roc(OutTest2, dup.post) # RIVER
  #GAM.roc <- roc(OutTest2, as.numeric(postprobTest)) # GAM

  #if (verbose) {
  #  cat('*** AUC (GAM - genomic annotation model): ',round(GAM.roc$auc,3),
  #    '\n    AUC (RIVER): ',round(RIVER.roc$auc,3),'\n     P-value: ',
  #    format.pval(roc.test(RIVER.roc, GAM.roc)$p.value,digits=2,eps=0.001),
 #     '***\n\n')
  #}

  evaROC <-
    list(RIVER_sens=RIVER.roc_observed$sensitivities,
         RIVER_spec=RIVER.roc_observed$specificities,
         RIVER_auc=RIVER.roc_observed$auc[1],
         GAM_sens=RIVER.roc_inferred$sensitivities,
         GAM_spec=RIVER.roc_inferred$specificities,
         GAM_auc=RIVER.roc_inferred$auc[1],
         pvalue=roc.test(RIVER.roc_observed, RIVER.roc_inferred)$p.value)
  class(evaROC) <- "eval"
  plot_roc(evaROC, Zscore_n2_threshold, paste0(output_root, "roc_curve_",Zscore_initialization_threshold,"_",Zscore_n2_threshold, "_",pseudoc,".pdf"))

}


initialize_phi<- function(num_bins,dim) {
  phi_outlier <- matrix(1,dim,num_bins)
  phi_inlier <- matrix(1,dim,num_bins)
  phi_inlier[,1] = .5
  phi_inlier[,2] = .2
  phi_inlier[,3] = .1
  phi_inlier[,4] = .1
  phi_inlier[,5] = .05
  phi_inlier[,6] = .05


  phi_outlier[,1] = .05
  phi_outlier[,2] = .05
  phi_outlier[,3] = .1
  phi_outlier[,4] = .1
  phi_outlier[,5] = .2
  phi_outlier[,6] = .5


  phi_init <- list(inlier_component = phi_inlier, outlier_component = phi_outlier)

}


print("Watershed Ising Model :)")
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_root = args[2]
version = args[3]  # Run "full_data_visualization" or "roc_analysis"
initialization = args[4]  # How to initialize parameters ("smart_initialization" or "zeros")
Zscore_n2_threshold = as.numeric(args[5])
vi_damping_factor = as.numeric(args[6])

print(version)
print(initialization)
print(Zscore_n2_threshold)
print(vi_damping_factor)

pseudoc=100
dimensions=1 # number of tissues

costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)
verbose=TRUE
ZscoreThrd = 3.3
num_bins <- 6

phi_init <- initialize_phi(num_bins,dimensions) 

if (version == "full_data_visualization") {
  full_data_visualization_driver(input_file, ZscoreThrd, output_root, dimensions, phi_init,num_bins, costs, verbose, pseudoc, initialization)
}

Zscore_initialization_threshold=3.3
if (version == "roc_analysis") {
  roc_analysis_driver(input_file, Zscore_initialization_threshold, Zscore_n2_threshold, output_root, dimensions, phi_init, num_bins, costs, verbose, pseudoc, initialization, vi_damping_factor)
}

