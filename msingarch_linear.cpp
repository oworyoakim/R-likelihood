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
  
  // Parameters
  PARAMETER_VECTOR(ta);       // m working values of a
  PARAMETER_VECTOR(tb);       // m working values of b
  PARAMETER_VECTOR(td);       // m working values of d 
  PARAMETER_VECTOR(tgamma);   // m(m-1) working parameters of TPM
  
  // Uncomment only when initial distribution of the Markov-Chain should be estimated
  //PARAMETER_VECTOR(tdelta);    // transformed stationary distribution,

  // Transform working parameters to natural parameters:
  vector<Type> a = ta.exp();
  vector<Type> b = tb.exp();
  vector<Type> d = td.exp();
  matrix<Type> gamma = Gamma_w2n(m, tgamma);
  
  // Construct stationary distribution
  vector<Type> delta = Stat_dist(m, gamma);
  
  // Uncomment only when initial distribution of the Markov-Chain should be estimated
  //vector<Type> delta = Delta_w2n(m, tdelta);
  
  // Get number of timesteps 
  int n = y.size();
  
  // Evaluate conditional distribution: Put conditional
  // probabilities of observed y in n times m matrix
  // (one column for each regime, one row for each datapoint):
  
  matrix<Type> lambda_mat(n,m);
  matrix<Type> emission_probs(n, m);
  matrix<Type> row1vec(1, m);
  row1vec.setOnes();
  
  // for t = 1;
  if (y[0] != y[0]) { // f != f returns true if and only if f is NaN.
    // Replace missing values (NA in R, NaN in C++) with 1
    emission_probs.row(0) = row1vec;
  }
    else {
      vector<Type> lambda0_vec = d + a*lambda0 + b*y0; 
      lambda_mat.row(0) = lambda0_vec;
      emission_probs.row(0) = dpois(y(0), lambda0_vec, false);
    }
  
  // for t = 1,2,...,T  
  for (int i = 1; i < n; i++) {
    if (y[i] != y[i]) { // f != f returns true if and only if f is NaN.
      // Replace missing values (NA in R, NaN in C++) with 1
      emission_probs.row(i) = row1vec;
    }
    else {
      vector<Type> lambda_t_minus_one = lambda_mat.row(i - 1);
      vector<Type> lambda_vec = d + a*lambda_t_minus_one + b*y(i - 1);
      lambda_mat.row(i) = lambda_vec;
      emission_probs.row(i) = dpois(y(i), lambda_vec, false);
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

  // Variables we need for local decoding and in a convenient format
  REPORT(a);
  REPORT(b);
  REPORT(d);
  REPORT(gamma);
  REPORT(delta);
  REPORT(n);
  REPORT(emission_probs);
  REPORT(mllk);
  REPORT(lambda_mat);
  
  return mllk;
}
