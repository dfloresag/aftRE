// AFT Model with random slope and intercept.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_FACTOR(group_1);     // already expanded
  DATA_FACTOR(group_2);     // already expanded
  DATA_FACTOR(group_12);
  DATA_MATRIX(X);
  DATA_VECTOR(w);           // weights just in case
  DATA_VECTOR(t);           // durations 
  DATA_VECTOR(d);           // censoring indicator
  DATA_INTEGER(cond_resp);  // 1 = ln ; 2 = ll ; 3 = wb ;
  
  PARAMETER_VECTOR(u1);     // non-expanded
  PARAMETER_VECTOR(u2);     // non-expanded
  PARAMETER_VECTOR(betas);   
  PARAMETER(log_sigma_0);   // Scale
  PARAMETER(log_sigma_1);   // RE 1
  PARAMETER(log_sigma_2);   // RE 2
  
  using namespace density;
  
  int n = t.size(); int n_1 = u1.size(); int n_2 = u2.size();
  
  Type res(0.0); 
  
  Type sigma_1(exp(log_sigma_1));
  Type sigma_2(exp(log_sigma_2));
  Type sigma_0(exp(log_sigma_0));
  
  Type sigma2_1(exp(2*log_sigma_1));
  Type sigma2_2(exp(2*log_sigma_2));
  Type sigma2_0(exp(2*log_sigma_0));
  
  ADREPORT(sigma_0);
  ADREPORT(sigma_1);
  ADREPORT(sigma_2);
  
  ADREPORT(sigma2_0);
  ADREPORT(sigma2_1);
  ADREPORT(sigma2_2);
  
  
  /* Prior: intercept~N(mu,sd) */
  
  for(int i=0;i<n_1;i++){
    res-=w[i]*dnorm(u1[i], Type(0), sigma_1 , 1);
  }
  
  for(int ij=0; ij<n_2; ij++){
    int i = group_12[ij];
    res-= w[i]*dnorm(u2[ij], Type(0), sigma_2, 1);
  }
  
  vector<Type> log_t  = log(t);
  vector<Type> e(n);
  vector<Type> mu = X*betas;  
  
  /* Observations: T|u ~ logNormal(e,sigma_0) */
  for(int ijk=0; ijk<n; ijk++){
    
    int i =group_1[ijk];
    int ij=group_2[ijk];
    
    e[ijk]= (log_t[ijk]-mu[ijk]-u1[i]-u2[ij])/exp(log_sigma_0);
    
    if (cond_resp==1){
      // log-normal
      res-= w[i]*(d[ijk]*(-log_t[ijk]-log_sigma_0 + dnorm(e[ijk], Type(0), Type(1), 1) - log(pnorm(-e[ijk], Type(0), Type(1)))) + log(1 - pnorm(e[ijk],Type(0), Type(1))));
    } else if(cond_resp==2){
      // log-logistic
      res-= w[i]*(d[ijk]*(-log_t[ijk]-log_sigma_0 -log(1+exp(-e[ijk]))) - log(1+exp(e[ijk])));
    } else if(cond_resp==3){
      // Weibull
      res-= w[i]*(d[ijk]*(-log_t[ijk]-log_sigma_0 + e[ijk])-exp(e[ijk]));  
    }
  }
  return res;
}
