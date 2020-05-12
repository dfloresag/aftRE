// AFT Model with random slope and intercept.
#include <TMB.hpp>
#include<cmath>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_FACTOR(group_1);       // already expanded
  DATA_FACTOR(group_2);       // already expanded
  DATA_FACTOR(group_12);
  DATA_MATRIX(X);
  DATA_VECTOR(w);             // weights just in case
  DATA_VECTOR(t);             // durations
  DATA_VECTOR(d);             // censoring indicator
  DATA_VECTOR(rand_strc_1);   // 1 = ri ; 2 = ri rs
  DATA_VECTOR(rand_strc_2);   // 0 = no ; 1 = ri ; 2 = ri rs

  DATA_INTEGER(cond_resp);    // 1 = ln ; 2 = ll ; 3 = wb

  DATA_VECTOR(s1);             // random slope
  DATA_VECTOR(s2);             // random slope

  PARAMETER_VECTOR(u11);     // non-expanded
  PARAMETER_VECTOR(u12);     // non-expanded
  PARAMETER_VECTOR(u21);     // non-expanded
  PARAMETER_VECTOR(u22);     // non-expanded
  PARAMETER_VECTOR(betas);
  PARAMETER(log_sigma_0);    // Scale
  PARAMETER(log_sigma_11);   // RE 1
  PARAMETER(log_sigma_12);   // RE 2
  PARAMETER(log_sigma_21);   // RE 1
  PARAMETER(log_sigma_22);   // RE 2

  PARAMETER(trans_rho_1);   // RE 1
  PARAMETER(trans_rho_2);   // RE 2

  int n = t.size(); int n_1 = u11.size(); int n_2 = u21.size();

  Type res(0.0);
  Type P_I(3.141592653589793238463);


  Type sigma_11(exp(log_sigma_11));
  Type sigma_12(exp(log_sigma_12));
  Type sigma_21(exp(log_sigma_21));
  Type sigma_22(exp(log_sigma_22));
  Type sigma_0(exp(log_sigma_0));

  Type sigma2_11(exp(2*log_sigma_11));
  Type sigma2_12(exp(2*log_sigma_12));
  Type sigma2_21(exp(2*log_sigma_21));
  Type sigma2_22(exp(2*log_sigma_22));
  Type sigma2_0(exp(2*log_sigma_0));

  ADREPORT(sigma_0);
  ADREPORT(sigma_11);
  ADREPORT(sigma_12);
  ADREPORT(sigma_21);
  ADREPORT(sigma_22);

  ADREPORT(sigma2_0);
  ADREPORT(sigma2_11);
  ADREPORT(sigma2_12);
  ADREPORT(sigma2_21);
  ADREPORT(sigma2_22);

  /* Prior: intercept~N(mu,sd) */

  if(rand_strc_1 == 1){
    for(int i=0;i<n_1;i++){
      res-=w[i]*dnorm(u11[i], Type(0), sigma_11 , 1);
    }
  } else if(rand_strc_1 == 2){
    Type rho_1(2/P_I*atan(trans_rho_1));
    ADREPORT(rho_1);

    matrix<Type> Sigma_1(2,2);
    Sigma_1.fill(rho_1);
    Sigma_1.diagonal() *= 1.0/rho_1;

    vector<Type> u(2);
    for(int i=0;i<n_1;i++){
      u[0] = u12[i];
      u[1] = u12[i];
      res += w[i]*MVNORM(Sigma_1)(u);
    }
  }

  /* Prior: intercept~N(mu,sd) */

  if(rand_strc_2 == 1){
    for(int ij=0; ij<n_2; ij++){
      int i = group_12[ij];
      res-= w[i]*dnorm(u21[ij], Type(0), sigma_21, 1);
    }
  } else if(rand_strc_2 == 2){
    Type rho_2(2/P_I*atan(trans_rho_2));
    ADREPORT(rho_2);

    matrix<Type> Sigma_2(2,2);
    Sigma_2.fill(rho_2);
    Sigma_2.diagonal() *= 1.0/rho_2;

    vector<Type> u(2);
    for(int ij=0; ij<n_2; ij++){
      int i = group_12[ij];
      u[0] = u21[ij];
      u[1] = u22[ij];
      res += w[ij]*MVNORM(Sigma_2)(u);
    }
  } else {
    res = res;
  }

  vector<Type> log_t  = log(t);
  vector<Type> xtb = X*betas;
  vector<Type> eta(n);
  vector<Type> ztu(n);

  /* Observations: T|u ~ logNormal(e,sigma_0) */

  for(int ijk=0; ijk<n; ijk++){

    if(rand_strc_1 == 0){
      ztu[ijk] = 0.0;
    } else {
      int i =group_1[ijk];
      if(rand_strc_1 == 1){
        ztu[ijk] = sigma_11*u11[i];
        if(rand_strc_2 == 0){
          ztu[ijk] += 0.0;
        } else {
          int ij =group_2[ijk];
          if (rand_strc_2 == 1){
            ztu[ijk] += sigma_21*u21[ij];
          } else if (rand_strc_2 == 2){
            ztu[ijk] += sigma_21*u21[ij]-sigma_22*s2*u22[ij];
          }
        }
      } else if(rand_strc_1 == 2){
          ztu[ijk] = sigma_11*u11[i]+;
        } else {
          int ij =group_2[ijk];
          if (rand_strc_2 == 1){
            ztu[ijk] = sigma_11*u11[i]+sigma_21*u21[ij];
          } else if (rand_strc_2 == 2){
            ztu[ijk] = (log_t[ijk]-xtb[ijk]-sigma_11*u11[i]-sigma_21*u21[ij]-sigma_22*s2*u22[ij])/exp(log_sigma_0);
          }
        }

    }


    }
    eta[ijk] = (log_t[ijk]-xtb[ijk]-sigma_11*u11[i])/exp(log_sigma_0);
  }
    int ij=group_2[ijk];
    if (rand_strc_1 ==1){
      if(rand_strc_2 ==0){
        eta[ijk] = (log_t[ijk]-xtb[ijk]-u11[i]-u21[ij])/exp(log_sigma_0);

      }

      }

  }

  for(int ijk=0; ijk<n; ijk++){

    int i =group_1[ijk];
    int ij=group_2[ijk];

    if (cond_resp==1){
      // log-normal
      res-= w[i]*(d[ijk]*(-log_t[ijk]-log_sigma_0 + dnorm(eta[ijk], Type(0), Type(1), 1) - log(pnorm(-eta[ijk], Type(0), Type(1)))) + log(1 - pnorm(eta[ijk],Type(0), Type(1))));
    } else if(cond_resp==2){
      // log-logistic
      res-= w[i]*(d[ijk]*(-log_t[ijk]-log_sigma_0 -log(1+exp(-eta[ijk]))) - log(1+exp(eta[ijk])));
    } else if(cond_resp==3){
      // Weibull
      res-= w[i]*(d[ijk]*(-log_t[ijk]-log_sigma_0 + eta[ijk])-exp(eta[ijk]));
    }
  }
  return res;
}
