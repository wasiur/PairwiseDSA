
functions{
  real[] pairwiseSIR(real t, real[] y, real[] parms, real[] x_r, int[] x_i){
    real R0 = parms[1];
    real n = parms[2];
    real gamma = parms[3];
    real rho = parms[4];
    
    real xi = 1.0*(n-1)/n;
    real tau = R0*gamma/(n*xi - R0);
    real b = n*tau; 
    real a = tau + gamma;
    
    real dydt[1];
    dydt[1] = - a*y[1]*(1-pow(y[1], xi-1))*inv(1-xi) - b*(1- pow(y[1], xi))*pow(y[1], xi) - b*rho*pow(y[1], xi);
    return dydt;
  }
  real[] recovery_time_density(real t, real[] y, real[] parms, real[] x_r, int[] x_i){
    real R0 = parms[1];
    real n = parms[2];
    real gamma = parms[3];
    real rho = parms[4];
    
    real xi = 1.0*(n-1)/n;
    real tau = R0*gamma/(n*xi - R0);
    real b = n*tau; 
    real a = tau + gamma;
    
    real dydt[5];
    dydt[1] = - a*y[1]*(1-pow(y[1], xi-1))*inv(1-xi) - b*(1- pow(y[1], xi))*pow(y[1], xi) - b*rho*pow(y[1], xi); // S
    dydt[2] = -gamma*y[2]; //h_1
    dydt[3] = -dydt[1]*exp(gamma*t); //h_2
    dydt[4] = dydt[2]*y[3]+ y[2]*dydt[3] ; //g
    dydt[5] = y[4] ; //G
    //dydt[3] = (- a*y[1]*(1-pow(y[1], xi-1))*inv(1-xi) - b*(1- pow(y[1], xi))*pow(y[1], xi) - b*rho*pow(y[1], xi))*exp(gamma*t) ;
    //dydt[4] = (-gamma*y[2])*y[3] + y[2]*(- a*y[1]*(1-pow(y[1], xi-1))*inv(1-xi) - b*(1- pow(y[1], xi))*pow(y[1], xi) - b*rho*pow(y[1], xi))*exp(gamma*t) ;
    dydt[5] = y[4];
    return dydt;
  }
}

data{
  int<lower=0> N; 
  real<lower=0.0> infection_times[N];
  int<lower=0> M;
  real<lower=0.0> recovery_times[M];
}

transformed data {
    real x_r[0];
    int x_i[0];
}

parameters {
    real<lower=0.50, upper=10.0> R0;
    real<lower=R0+1.0, upper=12.0> n;
    real<lower=0.03, upper=0.3> gamma;
    real<lower=0.0, upper=0.3> rho;
}

transformed parameters{
  real<lower=2.0> deg = round(n);
  real<lower=0.0, upper=1.0> xi = 1.0*(n-1)/n;
  //real gamma = gamma_0;
  real<lower=0.0> tau = R0*gamma/(n*xi - R0);
  real<lower=0.0> b = n*tau; 
  real<lower=0.0> a = tau + gamma;
  //real rho=1.0/10000;
}

model{
  real parms[4];
  real ic[1];
  real s[N,1];
  real r[M, 5];
  real r_ic[5];
  real t0;
  real Smax;
  real factor;
  int nSample;
  real y[3,1];
  real temp_sum;

  parms[1] = R0;
  parms[2] = n;
  parms[3] = gamma;
  parms[4] = rho;
  
  ic[1] = 1.0;
  
  t0 = 0.0; 
  
  s = integrate_ode_rk45(pairwiseSIR,ic,t0,infection_times,parms,x_r,x_i);
  Smax = s[N,1];
  factor = 1.0 - Smax;
  
  for (i in 1:N){
        target += log((a*s[i,1]*(1- pow(s[i,1], xi-1))*inv(1-xi) + b*(1-pow(s[i,1], xi))*pow(s[i,1], xi) + b*rho*pow(s[i,1],xi))/factor);
    }
  
  r_ic[1] = 1.0; r_ic[2] = gamma; r_ic[3] = n*tau*rho; r_ic[4] = n*tau*gamma*rho; r_ic[5] = 0.0; 
  r = integrate_ode_rk45(recovery_time_density,r_ic, t0, recovery_times, parms, x_r,x_i);
  Smax = r[M, 5];
    
  for (i in 1:M){
    target += log(r[i, 4]/Smax);
  }
  
  // target += exponential_lpdf(gamma | gamma_0);

  //target += gamma_lpdf(tau|2,2) + uniform_lpdf(R0 | 1.0, 5.0)+ uniform_lpdf(gamma|0.03,0.5) + uniform_lpdf(n|2.0,20.0);
  
}





