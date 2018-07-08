data {
  int<lower=1> K; // number of hidden states
  int<lower=1> span; //number of years
  int<lower=1> nt; //number of trees
  int start[nt]; //start year of tree[nt]
  int end[nt]; //end year of tree[nt]
  row_vector [span] age[nt]; //nt x age matrix
  int<lower=1> N_mis;//417
  matrix [nt,span] y; //logy
  vector[span-N_mis] x_obs; //observed temp
  vector[K] alpha; //c(1,1,1)
  int b[nt,span];
  int len[span];
}

transformed data{}

parameters {
  real<lower=0> sigmay; //sd of y
  real<lower=0> sigmax; //sd of mux
  real<lower=0> sigmaeta; //sd of eta
  real<lower=0> sigbeta0; //sd of beta0
  //real<lower=0> sigbeta1; sd of beta1
  real<lower=0> sigbeta2; //sd of beta2
  
  real<lower=-4,upper=4> mu_x[K]; //ordered mu
  //vector[nt] beta1;
  vector[nt] beta0;
  real <lower=0>beta2; 
  
  real betamu0;
  real betamu1;
  real betamu2;
  
  vector<lower=-4> [N_mis] x_mis; //no bound
  simplex[K] A[K];
  vector[span] eta;
  
  vector[42] x_pred1; //prediction of temperature from 1913-1954
  vector[42] eta_pred1; //prediction of temperature effect eta
  vector[41] x_pred2; //prediction of temperature from 1955-1995
  vector[41] eta_pred2; //prediction of temperature efect eta
  
}

transformed parameters{
  real<lower=-4> mux[K]; //lower bounded by 4,mux[1]..
  mux=sort_asc(mu_x); //ascending order
}

model {
  //transition matrix
  for (u in 1:K) A[u] ~ dirichlet(alpha);
  
  //piror for missing x, x ~N(0,100^2)
  target+= normal_lpdf(x_mis|0,100);
  
  //priror of mux, mux[i]~N(0,100^2)
  for(i in 1:K){
    target+= normal_lpdf(mux[i]|0,100);
  }
  
  //for eta with observed temperature x[i], eta[i]~N(beta2*x[i], sigmaeta^2)
  for (i in (N_mis+1):span){
    target+=normal_lpdf(eta[i]|beta2*x_obs[i-N_mis],sigmaeta);
  }
  
  //for eta with unobserved temperature x[i]
  for (i in 1:N_mis){
    target+=normal_lpdf(eta[i]|beta2*x_mis[i],sigmaeta);
  }
  
  // log y[j,t]~N(beta0[j]+age[j,t]*betamu1+eta[t],sigmay^2)
  for (j in 1:nt){
    for(t in start[j]:end[j]){
      target += normal_lpdf(y[j,t] | beta0[j]+betamu1*age[j,t]+eta[t],sigmay);
    }
  }
  
  //beta prior
  beta0~normal(betamu0,sigbeta0);
  //beta1~normal(betamu1,sigbeta1);
  beta2~normal(betamu2,sigbeta2);
  
  //mu prior
  betamu0~normal(0,100);
  betamu1~normal(0,100);
  betamu2~normal(0,100);

  //sigma prior
  sigmay ~ student_t(3,0,1);
  sigmax ~ student_t(3,0,1);
  sigmaeta ~ student_t(3,0,1);
  sigbeta0~student_t(3,0,1);
  //sigbeta1~student_t(3,0,1);
  sigbeta2~student_t(3,0,1);

  {//forward algorithm
    real acc[K];
    real GAMMA[span,K];
    real maxacc[span];
    real GAMMA_pred1[84,K];//GAMMA_pred1 share first 417 rows as GAMMA..
    real GAMMA_pred2[42,K];//GAMMA_pred2 share first 
    
    #initialise
    for (k in 1:K){
      GAMMA[1,k]= normal_lpdf(x_mis[1]|mux[k],sigmax);
    }
    maxacc[1]= 0;
    
    for (t in 2:N_mis){
      for (k in 1:K){
        for (j in 1:K){
          acc[j]=GAMMA[t-1,j]+log(A[j,k])+normal_lpdf(x_mis[t]|mux[k],sigmax);
        }
        maxacc[t]=max(acc);
        GAMMA[t,k]=maxacc[t]+log_sum_exp(to_vector(acc)-rep_vector(maxacc[t],3));
      }
    }
    
    for (t in (N_mis+1):span){
      for (k in 1:K){
        for (j in 1:K){
          acc[j]=GAMMA[t-1,j]+log(A[j,k])+normal_lpdf(x_obs[t-N_mis]|mux[k],sigmax);
        }
        maxacc[t]=max(acc);
        GAMMA[t,k]=maxacc[t]+log_sum_exp(to_vector(acc)-rep_vector(maxacc[t],3));
      }
    }
    
    target+= log_sum_exp(GAMMA[span]);
    
    #prediction part
    for (k in 1:K){
      GAMMA_pred1[1,k]=GAMMA[N_mis,k];
      GAMMA_pred2[1,k]=GAMMA[N_mis+42,k];
    }
    #first half temperature
    for (t in (N_mis+1):(N_mis+42)){
      for (k in 1:K){
        for (j in 1:K){
          acc[j]=GAMMA_pred1[t-N_mis,j]+log(A[j,k])+normal_lpdf(x_pred1[t-N_mis]|mux[k],sigmax);
        }
        maxacc[t]=max(acc);
        GAMMA_pred1[t-N_mis+1,k]=maxacc[t]+log_sum_exp(to_vector(acc)-rep_vector(maxacc[t],3));
      }
    }
    
    for (t in (N_mis+43):(N_mis+83)){
      for (k in 1:K){
        for (j in 1:K){
          acc[j]=GAMMA_pred1[t-N_mis,j]+log(A[j,k])+normal_lpdf(x_obs[t-N_mis]|mux[k],sigmax);
        }
        maxacc[t]=max(acc);
        GAMMA_pred1[t-N_mis+1,k]=maxacc[t]+log_sum_exp(to_vector(acc)-rep_vector(maxacc[t],3));
      }
    }
    target+=log_sum_exp(GAMMA_pred1[84]);
    
    #second half temperature
    for (t in (N_mis+43):(N_mis+83)){
      for (k in 1:K){
        for (j in 1:K){
          acc[j]=GAMMA_pred2[t-N_mis-42,j]+log(A[j,k])+normal_lpdf(x_pred2[t-N_mis-42]|mux[k],sigmax);
        }
        maxacc[t]=max(acc);
        GAMMA_pred2[t-N_mis-42+1,k]=maxacc[t]+log_sum_exp(to_vector(acc)-rep_vector(maxacc[t],3));
      }
    }
    target+=log_sum_exp(GAMMA_pred2[42]);
  }
  
  #prediction
  for (i in (N_mis+1):(N_mis+42))
    target+=normal_lpdf(eta_pred1[i-N_mis]|beta2*x_pred1[i-N_mis],sigmaeta);
  
  
  for (i in (N_mis+43):(N_mis+83)){
    target+=normal_lpdf(eta_pred2[i-N_mis-42]|beta2*x_pred2[i-N_mis-42],sigmaeta);
  }
  
  #get the tree index between 418:459
  
  
  for (t in (N_mis+1):(N_mis+42)){
    for (l in 1:len[t]){
      int ii=b[l,t];
      target+=normal_lpdf(y[ii,t]|beta0[ii]+betamu1*age[ii,t]+eta_pred1[t-N_mis],sigmay);
    }
  }
  
  
  for (t in (N_mis+43):(N_mis+83)){
    for (l in 1:len[t]){
      int iii=b[l,t];
      target+=normal_lpdf(y[iii,t]|beta0[iii]+betamu1*age[iii,t]+eta_pred2[t-N_mis-42],sigmay);
    }
  }
  
  #now forward algorithm again..include code in previous chunk
  
  
}


generated quantities{
  int<lower=1,upper=K> x_star[span]; //predictions here is the most likely state sequenc
  real log_p_x_star;
  {
    int back_ptr[span, K]; //best
    real best_logp[span,K]; //best of log likelihood of pi
    real best_total_logp; //best of the log likelihood (total)
    for (k in 1:K)
      best_logp[1,k] = normal_lpdf(x_mis[1]|mux[k],sigmax);
    for (t in 2:N_mis) {
      for (k in 1:K) {
        best_logp[t, k] = negative_infinity();
        for (j in 1:K) {
          real logp;
          logp= best_logp[t-1,j]+ log(A[j,k])+normal_lpdf(x_mis[t]|mux[k],sigmax);
          if (logp > best_logp[t, k]) {
            back_ptr[t, k] = j;
            best_logp[t, k] = logp;
          }
        }
      }
    }
    
    for (t in (N_mis+1):span) {
      for(k in 1:K){
        best_logp[t, k] = negative_infinity();
        for (j in 1:K) {
          real logp1;
          logp1 = best_logp[t-1, j]+ log(A[j, k]) + normal_lpdf(x_obs[t-N_mis]|mux[k],sigmax);
          if (logp1 > best_logp[t, k]) {
            back_ptr[t, k] = j;
            best_logp[t, k] = logp1;
          }
        }
      }
    }
    
    log_p_x_star = max(best_logp[span]);
    for (k in 1:K){
      if (best_logp[span, k] == log_p_x_star)
        x_star[span] = k;
    }
    for (t in 1:(span - 1))
      x_star[span - t] = back_ptr[span - t + 1,x_star[span - t + 1]];
  }
}
