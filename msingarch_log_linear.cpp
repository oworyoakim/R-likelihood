#include <TMB.hpp>
#include "utils.cpp"  // Load functions needed (Gamma_w2n, Stat_dist)

// Likelihood for a Markov-switching INGARCH model.
template<class Type>
Type objective_function<Type>::operator() ()
{

  // Data
  DATA_VECTOR(y);             // timeseries vector
  DATA_INTEGER(m);            // Number of regimes m
  DATA_VECTOR(lambda0);       // Initial value of lambda0
  DATA_SCALAR(y0);            // Initial value of y0
  DATA_MATRIX(x);             // n x p design matrix

  // Parameters
  PARAMETER_VECTOR(a);       // m  values of a
  PARAMETER_VECTOR(b);       // m  values of b
  PARAMETER_VECTOR(d);       // m  values of d
  PARAMETER_VECTOR(tgamma);   // m(m-1) working parameters of TPM
  PARAMETER_MATRIX(Beta);     // p x m matrix of covariate effects

  // Uncomment only when initial distribution of the Markov-Chain should be estimated
  //PARAMETER_VECTOR(tdelta);    // transformed stationary distribution,

  // Transform working parameters to natural parameters:
  matrix<Type> gamma = Gamma_w2n(m, tgamma);

  // Construct stationary distribution
  vector<Type> delta = Stat_dist(m, gamma);

  // Uncomment only when initial distribution of the Markov-Chain should be estimated
  //vector<Type> delta = Delta_w2n(m, tdelta);


  // Compute fixed effects
  matrix<Type> x_Beta = x*Beta;

  // Get number of timesteps
  int n = y.size();

  // Evaluate conditional distribution: Put conditional
  // probabilities of observed y in n times m matrix
  // (one column for each regime, one row for each datapoint):

  matrix<Type> lambda_mat(n,m);
  matrix<Type> eta_mat(n,m);
  matrix<Type> emission_probs(n, m);
  matrix<Type> row1vec(1, m);
  row1vec.setOnes();

  // for t = 1;
  if (y[0] != y[0]) { // f != f returns true if and only if f is NaN.
    // Replace missing values (NA in R, NaN in C++) with 1
    emission_probs.row(0) = row1vec;
  }
    else {
      vector<Type> x_Beta_row0 = x_Beta.row(0);
      vector<Type> eta0_vec = d + x_Beta_row0 + a*lambda0.log() + b*log(y0 + 1);
      vector<Type> lambda0_vec = eta0_vec.exp();
      lambda_mat.row(0) = lambda0_vec;
      emission_probs.row(0) = dpois(y(0), lambda0_vec, false);
      eta_mat.row(0) = eta0_vec;
    }

  // for t = 1,2,...,T
  for (int i = 1; i < n; i++) {
    if (y[i] != y[i]) { // f != f returns true if and only if f is NaN.
      // Replace missing values (NA in R, NaN in C++) with 1
      emission_probs.row(i) = row1vec;
    }
    else {
      vector<Type> x_Beta_row_i = x_Beta.row(i);
      vector<Type> eta_t_minus_one = eta_mat.row(i - 1);
      vector<Type> eta_vec = d + x_Beta_row_i + a*eta_t_minus_one + b*log(y(i - 1) + 1);
      vector<Type> lambda_vec = eta_vec.exp();
      lambda_mat.row(i) = lambda_vec;
      emission_probs.row(i) = dpois(y(i), lambda_vec, false);
      eta_mat.row(i) = eta_vec;
    }
  }

  // Forward-algorithm. Corresponds to Zucchini (book) page 333
  matrix<Type> foo, P;
  Type mllk, sumfoo, lscale;

  foo = (delta * vector<Type>(emission_probs.row(0))).matrix();
  sumfoo = foo.sum();
  lscale = log(sumfoo);
  foo.transposeInPlace();
  foo /= sumfoo;
  for (int i = 2; i <= n; i++) {
    P = emission_probs.row(i - 1);
    foo = ((foo * gamma).array() * P.array()).matrix();
    sumfoo = foo.sum();
    lscale += log(sumfoo);
    foo /= sumfoo;
  }
  mllk = -lscale;

  // Use adreport on variables for which we want standard errors
  ADREPORT(a);
  ADREPORT(b);
  ADREPORT(d);
  ADREPORT(gamma);
  ADREPORT(delta);
  ADREPORT(Beta);

  // Variables we need for local decoding and in a convenient format
  REPORT(a);
  REPORT(b);
  REPORT(d);
  REPORT(gamma);
  REPORT(delta);
  REPORT(Beta);
  REPORT(n);
  REPORT(emission_probs);
  REPORT(mllk);
  REPORT(lambda_mat);

  return mllk;
}
