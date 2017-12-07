#include <vector>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <limits>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;

double sigmoid_fxn(double input) {
    double output = 0;
    output = 1.0 / (1.0 + exp(-input));
    return output;
}

double abs_difference(NumericMatrix mu_mat, int n, NumericVector prev_mu, int T) {
    double difference = 0;
    for (int t = 0; t < T; t++) {
        difference += fabs(mu_mat(n,t) - prev_mu[t]);
    }
    return difference;

}


double compute_elbo(int n, NumericMatrix mu_mat, NumericMatrix q_mat, NumericMatrix Feat, int T, NumericVector singleton, double** pairwise_mat, double sum_beta_feat) {
    double elbo = 0;
    for (int t = 0; t < T; t++) {
        elbo = elbo + sum_beta_feat*mu_mat(n,t) + singleton(t)*mu_mat(n,t);
        elbo = elbo - q_mat(n,t)*log(q_mat(n,t)) - (1.0-q_mat(n,t))*log(1.0-q_mat(n,t));
        for (int t_2 = 0; t_2 < T; t_2++) {
            if (t_2 > t) {
                elbo = elbo + pairwise_mat[t][t_2]*mu_mat(n,t)*mu_mat(n,t_2);
            }

        }
    }
    return elbo;
}

// [[Rcpp::export]]
NumericMatrix variational_inference_mu_updates(NumericMatrix Feat, NumericMatrix q_mat, NumericMatrix mu_mat, NumericVector singleton, NumericVector beta, NumericVector pairwise, double vi_damping_factor) {
    int N = mu_mat.nrow();
    int T = mu_mat.ncol();
    int num_feat = Feat.ncol();

    Rcout << std::fixed;
    Rcout << std::setprecision(20);




    double** pairwise_mat = new double*[T];
    int counter = 0;
    for(int i = 0; i < T; ++i) {
        pairwise_mat[i] = new double[T];
    }
    for (int t_1 = 0; t_1 < T; t_1++) {
        pairwise_mat[t_1][t_1] = 0;
    }
    for (int t_1 = 0; t_1 < (T-1); t_1++) {
        for (int t_2 = (t_1+1); t_2 < T; t_2++) {
            pairwise_mat[t_1][t_2] = pairwise[counter];
            pairwise_mat[t_2][t_1] = pairwise[counter];
            counter = counter + 1;
        }
    }


    double temp = 0.0;
    double elbo = 0.0;
    double elbo_old = -100000000000;
    for (int n = 0; n < N; n++) {
        bool converge = false;
        NumericVector prev_mu(T);
        elbo_old = -1000000000000000000;
        int itera = 0;
        while ( converge == false ) {
            double sum_beta_feat = 0;
            for (int d = 0; d < num_feat; d++) {
                sum_beta_feat += Feat(n,d)*beta[d];
            }
            for (int t = 0; t < T; t++) {
                prev_mu[t] = mu_mat(n,t);
                double sum_neighbor = 0;

                for (int t_2 = 0; t_2 < T; t_2++) {
                    sum_neighbor += pairwise_mat[t][t_2]*mu_mat(n,t_2);

                }
                temp = sigmoid_fxn(2.0*(sum_neighbor + sum_beta_feat + singleton[t]));
                mu_mat(n,t) = (2.0*temp - 1.0)*(1.0-vi_damping_factor) + (vi_damping_factor)*prev_mu[t];
                q_mat(n,t) = (mu_mat(n,t) + 1.0)/2.0;
            }
            itera = itera + 1;

            if (itera > 10000000) {
                converge = true;
                Rcout << "miss: " << n << std::endl;
            }
  
            //elbo = compute_elbo(n,mu_mat,q_mat,Feat,T,singleton,pairwise_mat,sum_beta_feat);
            //if (elbo - elbo_old < -0.000000000001) {
             //   printf("ERROR %d\n", n);
             //   Rcout << "ELBO_OLD: " << elbo_old << std::endl;
              //  Rcout << "ELBO_CUR: " << elbo << std::endl;
            //}
            //elbo_old = elbo;
            if (abs_difference(mu_mat, n, prev_mu, T) < .0000005) {
                converge = true;
            }
        }

    }
    return q_mat;
}


// [[Rcpp::export]]
NumericMatrix pseudo_mu_updates(NumericMatrix Feat, NumericMatrix q_mat, NumericMatrix mu_mat, NumericVector singleton, NumericVector beta, NumericVector pairwise, NumericMatrix posterior_mu) {
    int N = mu_mat.nrow();
    int T = mu_mat.ncol();
    int num_feat = Feat.ncol();

    Rcout << std::fixed;
    Rcout << std::setprecision(20);




    double** pairwise_mat = new double*[T];
    int counter = 0;
    for(int i = 0; i < T; ++i) {
        pairwise_mat[i] = new double[T];
    }
    for (int t_1 = 0; t_1 < T; t_1++) {
        pairwise_mat[t_1][t_1] = 0;
    }
    for (int t_1 = 0; t_1 < (T-1); t_1++) {
        for (int t_2 = (t_1+1); t_2 < T; t_2++) {
            pairwise_mat[t_1][t_2] = pairwise[counter];
            pairwise_mat[t_2][t_1] = pairwise[counter];
            counter = counter + 1;
        }
    }


    double temp = 0.0;
    for (int n = 0; n < N; n++) {
        double sum_beta_feat = 0;
        for (int d = 0; d < num_feat; d++) {
            sum_beta_feat += Feat(n,d)*beta[d];
        }
        for (int t = 0; t < T; t++) {
            double sum_neighbor = 0;

            for (int t_2 = 0; t_2 < T; t_2++) {
                sum_neighbor += pairwise_mat[t][t_2]*posterior_mu(n,t_2);

            }
            temp = sigmoid_fxn(2.0*(sum_neighbor + sum_beta_feat + singleton[t]));
            mu_mat(n,t) = (2.0*temp - 1.0);
            q_mat(n,t) = (mu_mat(n,t) + 1.0)/2.0;
        }
    }
    return q_mat;
}


// [[Rcpp::export]]
NumericMatrix variational_inference_posterior_updates(NumericMatrix Feat, NumericMatrix q_mat, NumericMatrix mu_mat, NumericVector singleton, NumericVector beta, NumericVector pairwise, NumericMatrix ll_plus, NumericMatrix ll_minus, double vi_damping_factor) {
    int N = mu_mat.nrow();
    int T = mu_mat.ncol();
    int num_feat = Feat.ncol();

    Rcout << std::fixed;
    Rcout << std::setprecision(20);




    double** pairwise_mat = new double*[T];
    int counter = 0;
    for(int i = 0; i < T; ++i) {
        pairwise_mat[i] = new double[T];
    }
    for (int t_1 = 0; t_1 < T; t_1++) {
        pairwise_mat[t_1][t_1] = 0;
    }
    for (int t_1 = 0; t_1 < (T-1); t_1++) {
        for (int t_2 = (t_1+1); t_2 < T; t_2++) {
            pairwise_mat[t_1][t_2] = pairwise[counter];
            pairwise_mat[t_2][t_1] = pairwise[counter];
            counter = counter + 1;
        }
    }

    int itera = 0.0;
    double temp = 0.0;
    double elbo = 0.0;
    double elbo_old = -100000000000;
    for (int n = 0; n < N; n++) {
        bool converge = false;
        NumericVector prev_mu(T);
        elbo_old = -1000000000000000000;
        int itera = 0;
        while ( converge == false ) {
            double sum_beta_feat = 0;
            for (int d = 0; d < num_feat; d++) {
                sum_beta_feat += Feat(n,d)*beta[d];
            }
            for (int t = 0; t < T; t++) {
                prev_mu[t] = mu_mat(n,t);
                double sum_neighbor = 0;

                for (int t_2 = 0; t_2 < T; t_2++) {
                    sum_neighbor += pairwise_mat[t][t_2]*mu_mat(n,t_2);

                }
                temp = sigmoid_fxn(ll_plus(n,t) - ll_minus(n,t) + 2.0*(sum_neighbor + sum_beta_feat + singleton[t]));
                mu_mat(n,t) = (2.0*temp - 1.0)*(1.0-vi_damping_factor) + (vi_damping_factor)*prev_mu[t];
                q_mat(n,t) = (mu_mat(n,t) + 1.0)/2.0;
            }

            itera = itera + 1;
            if (itera > 10000000) {
                converge = true;
                Rcout << "missy: " << n << std::endl;
            }

            //elbo = compute_elbo(n,mu_mat,q_mat,Feat,T,singleton,pairwise_mat,sum_beta_feat);
            //if (elbo - elbo_old < -0.000000000001) {
             //   printf("ERROR %d\n", n);
             //   Rcout << "ELBO_OLD: " << elbo_old << std::endl;
              //  Rcout << "ELBO_CUR: " << elbo << std::endl;
            //}
            //elbo_old = elbo;
     
            if (abs_difference(mu_mat, n, prev_mu, T) < .0000005) {
                converge = true;
            }
        }

    }
    return q_mat;
}




// [[Rcpp::export]]
NumericVector compute_singleton_gradient_cpp(NumericVector singleton_vec, NumericMatrix posterior_mu, NumericMatrix mu) {
    int N = mu.nrow();
    int T = mu.ncol();

    NumericVector gradient(T);

    for (int n = 0; n < N; n++) {
        for (int t = 0; t < T; t++) {
            gradient[t] = gradient[t] + posterior_mu(n,t) - mu(n,t);
        }
    }

    for (int t = 0; t < T; t++) {
        gradient[t] = gradient[t]/(double)N;
    }
    return gradient;
}

// [[Rcpp::export]]
NumericVector compute_beta_gradient_cpp(NumericVector beta_vec, NumericMatrix posterior_mu, NumericMatrix mu, NumericMatrix feat, double lambda) {
    int N = mu.nrow();
    int T = mu.ncol();
    int num_feat = feat.ncol();

    NumericVector gradient(num_feat);

    for (int n = 0; n < N; n++) {
        for (int t = 0; t < T; t++) {
            for (int d = 0; d < num_feat; d++) {
                gradient[d] = gradient[d] + posterior_mu(n,t)*feat(n,d) - mu(n,t)*feat(n,d);
            }
        }
    }

    for (int d = 0; d < num_feat; d++) {
        gradient[d] = (gradient[d]/(double)N) - lambda*beta_vec[d];
    }
    return gradient;
}

// [[Rcpp::export]]
NumericVector compute_pair_gradient_cpp(NumericVector pair_vec, NumericMatrix posterior_mu, NumericMatrix mu, double lambda) {
    int N = mu.nrow();
    int T = mu.ncol();
    int num_pairs = (T*T - T)/2;

    NumericVector gradient(num_pairs);

    for (int n = 0; n < N; n++) {
        int counter = 0;
        for (int t_1 = 0; t_1 < (T-1); t_1++) {
            for (int t_2 = (t_1+1); t_2 < T; t_2++) {
                gradient[counter] = gradient[counter] + posterior_mu(n,t_1)*posterior_mu(n,t_2) - mu(n,t_1)*mu(n,t_2);
                counter = counter + 1;
            }
        }
    }

    for (int d = 0; d < num_pairs; d++) {
        gradient[d] = (gradient[d]/(double)N) - lambda*pair_vec[d];
    }
    return gradient;
}






// [[Rcpp::export]]
NumericVector compute_singleton_gradient_cpp_2(NumericVector singleton_vec, NumericMatrix posterior_mu, NumericMatrix mu) {
    int N = mu.nrow();
    int T = mu.ncol();

    NumericVector gradient(T);

    for (int n = 0; n < N; n++) {
        for (int t = 0; t < T; t++) {
            gradient[t] = gradient[t] + posterior_mu(n,t) - mu(n,t);
        }
    }

    for (int t = 0; t < T; t++) {
        gradient[t] = gradient[t];
    }
    return gradient;
}

// [[Rcpp::export]]
NumericVector compute_singleton_gradient_pseudo_cpp(NumericVector singleton_vec, NumericMatrix posterior_mu, NumericMatrix mu) {
    int N = mu.nrow();
    int T = mu.ncol();

    NumericVector gradient(T);

    for (int n = 0; n < N; n++) {
        for (int t = 0; t < T; t++) {
            gradient[t] = gradient[t] + posterior_mu(n,t) - mu(n,t);
        }
    }
    return gradient;
}


// [[Rcpp::export]]
NumericVector compute_beta_gradient_cpp_2(NumericVector beta_vec, NumericMatrix posterior_mu, NumericMatrix mu, NumericMatrix feat, double lambda) {
    int N = mu.nrow();
    int T = mu.ncol();
    int num_feat = feat.ncol();

    NumericVector gradient(num_feat);

    for (int n = 0; n < N; n++) {
        for (int t = 0; t < T; t++) {
            for (int d = 0; d < num_feat; d++) {
                gradient[d] = gradient[d] + posterior_mu(n,t)*feat(n,d) - mu(n,t)*feat(n,d);
            }
        }
    }

    for (int d = 0; d < num_feat; d++) {
        gradient[d] = gradient[d] - (double)N*lambda*beta_vec[d];
    }
    return gradient;
}

// [[Rcpp::export]]
NumericVector compute_beta_gradient_pseudo_cpp(NumericVector beta_vec, NumericMatrix posterior_mu, NumericMatrix mu, NumericMatrix feat, double lambda) {
    int N = mu.nrow();
    int T = mu.ncol();
    int num_feat = feat.ncol();

    NumericVector gradient(num_feat);

    for (int n = 0; n < N; n++) {
        for (int t = 0; t < T; t++) {
            for (int d = 0; d < num_feat; d++) {
                gradient[d] = gradient[d] + posterior_mu(n,t)*feat(n,d) - mu(n,t)*feat(n,d);
            }
        }
    }

    for (int d = 0; d < num_feat; d++) {
        gradient[d] = gradient[d] - (double)N*lambda*beta_vec[d];
    }
    return gradient;
}


// [[Rcpp::export]]
NumericVector compute_pair_gradient_cpp_2(NumericVector pair_vec, NumericMatrix posterior_mu, NumericMatrix mu, double lambda) {
    int N = mu.nrow();
    int T = mu.ncol();
    int num_pairs = (T*T - T)/2;

    NumericVector gradient(num_pairs);

    for (int n = 0; n < N; n++) {
        int counter = 0;
        for (int t_1 = 0; t_1 < (T-1); t_1++) {
            for (int t_2 = (t_1+1); t_2 < T; t_2++) {
                gradient[counter] = gradient[counter] + posterior_mu(n,t_1)*posterior_mu(n,t_2) - mu(n,t_1)*mu(n,t_2);
                counter = counter + 1;
            }
        }
    }

    for (int d = 0; d < num_pairs; d++) {
        gradient[d] = gradient[d] - (double)N*lambda*pair_vec[d];
    }
    return gradient;
}

// [[Rcpp::export]]
NumericVector compute_pair_gradient_pseudo_cpp(NumericVector pair_vec, NumericMatrix posterior_mu, NumericMatrix mu, double lambda) {
    int N = mu.nrow();
    int T = mu.ncol();
    int num_pairs = (T*T - T)/2;

    NumericVector gradient(num_pairs);

    for (int n = 0; n < N; n++) {
        int counter = 0;
        for (int t_1 = 0; t_1 < (T-1); t_1++) {
            for (int t_2 = (t_1+1); t_2 < T; t_2++) {
                gradient[counter] = gradient[counter] + 2.0*posterior_mu(n,t_1)*posterior_mu(n,t_2) - posterior_mu(n,t_1)*mu(n,t_2) - mu(n,t_1)*posterior_mu(n,t_2);
                counter = counter + 1;
            }
        }
    }

    for (int d = 0; d < num_pairs; d++) {
        gradient[d] = gradient[d] - (double)N*lambda*pair_vec[d];
    }
    return gradient;
}


// [[Rcpp::export]]
double crf_log_likelihood_helper_cpp(NumericVector intercept, NumericVector pair_vec, NumericVector beta, NumericMatrix mu, NumericMatrix q, NumericMatrix posterior_mu, NumericMatrix posterior, NumericMatrix Feat, int T, double lambda, double lambda_theta_pair) {
    int N = mu.nrow();
    int num_feat = Feat.ncol();
    double beta_term = 0;
    double pair_term = 0;
    double beta_reg_term = 0;
    double pair_reg_term = 0;
    double singleton_term = 0;
    double non_differentiable_term = 0;
    for (int n = 0; n < N; n++) {
        for (int t=0; t < T; t++) {
            singleton_term = singleton_term + posterior_mu(n,t)*intercept[t] - mu(n,t)*intercept[t];
            double beta_g = 0;
            for (int d = 0; d < num_feat; d++) {
                beta_g = beta_g + beta[d]*Feat(n,d);
            }
            beta_term = beta_term + posterior_mu(n,t)*beta_g - mu(n,t)*beta_g;
            if (q(n,t) != 0.0) {
                non_differentiable_term = non_differentiable_term + q(n,t)*log(q(n,t));
            }
            if (mu(n,t) != 1.0) {
                non_differentiable_term = non_differentiable_term + (1.0-q(n,t))*log(1.0 - q(n,t));
            }
        }
        int counter = 0;
        for (int t_1 = 0; t_1 < (T-1); t_1++) {
            for (int t_2 = (t_1+1); t_2 < T; t_2++) {
                pair_term = pair_term + pair_vec[counter]*posterior_mu(n,t_1)*posterior_mu(n,t_2) - pair_vec[counter]*mu(n,t_1)*mu(n,t_2);
                counter = counter + 1;
            }
        }
    }

    double dot_prod = 0;
    for (int d=0; d < num_feat;d++) {
        dot_prod = dot_prod + beta[d]*beta[d];
    }
    beta_reg_term = beta_reg_term + (lambda/2.0)*dot_prod;
    int counter = 0;
    for (int t_1 = 0; t_1 < (T-1); t_1++) {
        for (int t_2 = (t_1+1); t_2 < T; t_2++) {
            pair_reg_term = pair_reg_term + pair_vec[counter]*pair_vec[counter];
            counter = counter + 1;
        }
    }
    pair_reg_term = pair_reg_term*(lambda_theta_pair/2.0);



    double ll =  singleton_term + beta_term + pair_term - (pair_reg_term)*(double)N - (beta_reg_term)*(double)N + non_differentiable_term;
    return ll;
}





// [[Rcpp::export]]
double crf_log_likelihood_helper_pseudo_cpp(NumericVector intercept, NumericVector pair_vec, NumericVector beta, NumericMatrix mu, NumericMatrix q, NumericMatrix posterior_mu, NumericMatrix posterior, NumericMatrix Feat, int T, double lambda, double lambda_theta_pair) {
    int N = mu.nrow();
    int num_feat = Feat.ncol();
    double beta_term = 0;
    double pair_term = 0;
    double beta_reg_term = 0;
    double pair_reg_term = 0;
    double singleton_term = 0;
    double normalization_constant = 0;


    double** pairwise_mat = new double*[T];
    int countery = 0;
    for(int i = 0; i < T; ++i) {
        pairwise_mat[i] = new double[T];
    }
    for (int t_1 = 0; t_1 < T; t_1++) {
        pairwise_mat[t_1][t_1] = 0;
    }
    for (int t_1 = 0; t_1 < (T-1); t_1++) {
        for (int t_2 = (t_1+1); t_2 < T; t_2++) {
            pairwise_mat[t_1][t_2] = pair_vec[countery];
            pairwise_mat[t_2][t_1] = pair_vec[countery];
            countery = countery + 1;
        }
    }




    for (int n = 0; n < N; n++) {
        for (int t=0; t < T; t++) {
            double normalization_temp_value = 0;
            normalization_temp_value = normalization_temp_value + intercept[t];
            singleton_term = singleton_term + posterior_mu(n,t)*intercept[t];
            double beta_g = 0;
            for (int d = 0; d < num_feat; d++) {
                beta_g = beta_g + beta[d]*Feat(n,d);
            }
            beta_term = beta_term + posterior_mu(n,t)*beta_g;
            normalization_temp_value = normalization_temp_value + beta_g;
            for (int t2=0;t2 < T;t2++) {
                normalization_temp_value = normalization_temp_value + posterior_mu(n,t2)*pairwise_mat[t][t2];
            }
            normalization_constant = normalization_constant + log(exp(-normalization_temp_value) + exp(normalization_temp_value));

        }
        int counter = 0;
        for (int t_1 = 0; t_1 < (T-1); t_1++) {
            for (int t_2 = (t_1+1); t_2 < T; t_2++) {
                pair_term = pair_term + 2.0*pair_vec[counter]*posterior_mu(n,t_1)*posterior_mu(n,t_2);
                counter = counter + 1;
            }
        }
    }

    double dot_prod = 0;
    for (int d=0; d < num_feat;d++) {
        dot_prod = dot_prod + beta[d]*beta[d];
    }
    beta_reg_term = beta_reg_term + (lambda/2.0)*dot_prod;
    int counter = 0;
    for (int t_1 = 0; t_1 < (T-1); t_1++) {
        for (int t_2 = (t_1+1); t_2 < T; t_2++) {
            pair_reg_term = pair_reg_term + pair_vec[counter]*pair_vec[counter];
            counter = counter + 1;
        }
    }
    pair_reg_term = pair_reg_term*(lambda_theta_pair/2.0);



    double ll =  singleton_term + beta_term + pair_term - (pair_reg_term)*(double)N - (beta_reg_term)*(double)N - normalization_constant;
    return ll;
}

