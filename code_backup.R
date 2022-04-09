# generate data
generate.data = function(N, alpha, beta1, beta2, omega){
  alpha1 = alpha[1]
  alpha2 = alpha[2]
  beta10 = beta1[1]
  beta11 = beta1[2]
  beta12 = beta1[3]
  beta20 = beta2[1]
  beta21 = beta2[2]
  beta22 = beta2[3]
  
  # Generate theta1 and theta2
  t1 = c()
  t2 = c()
  while (length(t1) < N){
    u1 = runif(1, 0, 1)
    z = runif(1, 0, 1)
    theta1 = qgamma(u1, alpha1, rate = alpha1)
    C = function(u2){
      f1 = exp(-1 * theta1) - (alpha1 / (alpha1 + 1))^alpha1
      f2 = (alpha2 / (alpha2 + 1))^alpha2 * (pgamma(qgamma(u2, alpha2, rate = alpha2), alpha2, rate = (alpha2 + 1)) - u2)
      return (u2 + omega * f1 * f2 - z)
    }
    r = unlist(unname(uniroot(C, lower = 0, upper = 1)[1]))
    theta2 = qgamma(r, alpha2, rate = alpha2)       
    if (theta1 != Inf && theta2 != Inf){
      t1 = c(t1, theta1)
      t2 = c(t2, theta2)
    }else{
      next
    }
  }
  
  # Generate covariate
  z1 = rnorm(N)
  z2 = sample(c(0, 1), N, p = c(.5, .5), replace = T)
  
  # Calculate the Lambda_i
  t = c(0, 1, 2, 3, 4, 5)
  lambda1_vec = c()
  lambda2_vec = c(0.5, 1, 1.5, 2, 2.5)
  for (i in 1:(length(t)-1)){
    v1 = integrate(f = function(x) 0.5 * x, lower = t[i], upper = t[i+1])$val
    #v2 = integrate(f = function(x) 0.5 * x^2, lower = t[i], upper = t[i+1])$val
    lambda1_vec = c(lambda1_vec, v1)
    #lambda2_vec = c(lambda2_vec, v2)
  }
  #print(lambda1_vec)
  #print(lambda2_vec)
  
  ind = length(lambda1_vec)
  n1_5 = c()
  n2_5 = c()    
  for (i in 1:ind){
    lambda1 = t1 * lambda1_vec[i] * exp(beta10 + beta11 * z1 + beta12 * z2)
    n1 = rpois(N, lambda = lambda1)
    n1_5 = cbind(n1_5, n1) 
  }
  
  for (i in 1:ind){
    lambda2 = t2 * lambda2_vec[i] * exp(beta20 + beta21 * z1 + beta22 * z2)
    n2 = rpois(N, lambda = lambda2)
    n2_5 = cbind(n2_5, n2)
  }
  
  data = cbind(n1_5, n2_5, z1, z2)
  colnames(data) = c("N11", "N21", "N31", "N41", "N51", "N12", "N22", "N32", "N42", "N52", "Z1", "Z2")
  return (data)
}


# Loglikelihood function with exponential reparameteric
logL5 = function(par, data){
  alpha1 = par[1]
  beta10 = par[2]
  beta11 = par[3]
  beta12 = par[4]
  dl11 = par[5]
  dl12 = par[6]
  dl13 = par[7]
  dl14 = par[8]
  dl15 = par[9]
  alpha2 = par[10]
  beta20 = par[11]
  beta21 = par[12]
  beta22 = par[13]
  dl21 = par[14]
  dl22 = par[15]
  dl23 = par[16]
  dl24 = par[17]
  dl25 = par[18]
  omega = par[19]
  
  yi11 = data[, 1]
  yi12 = data[, 2]
  yi13 = data[, 3]
  yi14 = data[, 4]
  yi15 = data[, 5]
  yi21 = data[, 6]
  yi22 = data[, 7]
  yi23 = data[, 8]
  yi24 = data[, 9]
  yi25 = data[, 10]
  z1 = data[, 11]
  z2 = data[, 12]
  
  sum_k_yi1k = yi11 + yi12 + yi13 + yi14 + yi15
  sum_k_dl1k = exp(dl11) + exp(dl12) + exp(dl13) + exp(dl14) + exp(dl15)
  
  sum_k_yi2k = yi21 + yi22 + yi23 + yi24 + yi25
  sum_k_dl2k = exp(dl21) + exp(dl22) + exp(dl23) + exp(dl24) + exp(dl25)
  
  exp1 = exp(exp(beta10) + exp(beta11) * z1 + exp(beta12) * z2)
  exp2 = exp(exp(beta20) + exp(beta21) * z1 + exp(beta22) * z2)
  
  C1 = sum(yi11 * dl11 + yi12 * dl12 + yi13 * dl13 + yi14 * dl14 + dl15 * dl15)
  C2 = lgamma(exp(alpha1) + sum_k_yi1k)
  C3 = lgamma(exp(alpha1))
  C4 = lgamma(sum_k_yi1k + 1)
  C5 = alpha1 * exp(alpha1)
  C6 = (exp(beta10) + exp(beta11) * z1 + exp(beta12) * z2) * sum_k_yi1k
  C7 = (exp(alpha1) + sum_k_yi1k) * log(exp(alpha1) + exp1 * sum_k_dl1k)
  
  C8 = sum(yi21 * dl21 + yi22 * dl22 + yi23 * dl23 + yi24 * dl24 + yi25 * dl25)
  C9 = lgamma(exp(alpha2) + sum_k_yi2k)
  C10 = lgamma(exp(alpha2))
  C11 = lgamma(sum_k_yi2k + 1)
  C12 = alpha2 * exp(alpha2)
  C13 = (exp(beta20) + exp(beta21) * z1 + exp(beta22) * z2) * sum_k_yi2k
  C14 = (exp(alpha2) + sum_k_yi2k) * log(exp(alpha2) + exp2 * sum_k_dl2k)
  
  C151 = ((exp(alpha1) + exp1 * sum_k_dl1k)/(exp(alpha1) + exp1 * sum_k_dl1k + 1))^(exp(alpha1) + sum_k_yi1k)-(exp(alpha1)/(exp(alpha1) + 1))^(exp(alpha1))
  C153 = ((exp(alpha2) + exp2 * sum_k_dl2k)/(exp(alpha2) + exp2 * sum_k_dl2k + 1))^(exp(alpha2) + sum_k_yi2k)-(exp(alpha2)/(exp(alpha2) + 1))^(exp(alpha2))
  C15 = log(1 + exp(omega) * C151 * C153)
  
  return ((-1)*sum(C1 + C2 - C3 - C4 + C5 + C6 - C7 + C8 + C9 - C10 - C11 + C12 + C13 - C14 + C15))
}

g5 = function(par, data){
  alpha1 = par[1]
  beta10 = par[2]
  beta11 = par[3]
  beta12 = par[4]
  dl11 = par[5]
  dl12 = par[6]
  dl13 = par[7]
  dl14 = par[8]
  dl15 = par[9]
  alpha2 = par[10]
  beta20 = par[11]
  beta21 = par[12]
  beta22 = par[13]
  dl21 = par[14]
  dl22 = par[15]
  dl23 = par[16]
  dl24 = par[17]
  dl25 = par[18]
  omega = par[19]
  
  yi11 = data[, 1]
  yi12 = data[, 2]
  yi13 = data[, 3]
  yi14 = data[, 4]
  yi15 = data[, 5]
  yi21 = data[, 6]
  yi22 = data[, 7]
  yi23 = data[, 8]
  yi24 = data[, 9]
  yi25 = data[, 10]
  z1 = data[, 11]
  z2 = data[, 12]
  
  sum_k_yi1k = yi11 + yi12 + yi13 + yi14 + yi15
  sum_k_dl1k = exp(dl11) + exp(dl12) + exp(dl13) + exp(dl14) + exp(dl15)
  sum_k_yi2k = yi21 + yi22 + yi23 + yi24 + yi25
  sum_k_dl2k = exp(dl21) + exp(dl22) + exp(dl23) + exp(dl24) + exp(dl25)
  exp1 = exp(exp(beta10) + exp(beta11) * z1 + exp(beta12) * z2)
  exp2 = exp(exp(beta20) + exp(beta21) * z1 + exp(beta22) * z2)
  exp_alpha1_exp_exp_beta1z_sum_k_dl1k = exp(alpha1) + exp1 * sum_k_dl1k + 1
  exp_alpha2_exp_exp_beta2z_sum_k_dl2k = exp(alpha2) + exp2 * sum_k_dl2k + 1
  C151 = ((exp(alpha1) + exp1 * sum_k_dl1k)/(exp(alpha1) + exp1 * sum_k_dl1k + 1))^(exp(alpha1) + sum_k_yi1k)-(exp(alpha1)/(exp(alpha1) + 1))^(exp(alpha1))
  C153 = ((exp(alpha2) + exp2 * sum_k_dl2k)/(exp(alpha2) + exp2 * sum_k_dl2k + 1))^(exp(alpha2) + sum_k_yi2k)-(exp(alpha2)/(exp(alpha2) + 1))^(exp(alpha2))
  C15 = log(1 + exp(omega) * C151 * C153)
  
  p_alpha1 = sum(exp(alpha1)*(digamma(exp(alpha1)+sum_k_yi1k)-digamma(exp(alpha1))+1+alpha1-log(exp(alpha1)+exp1*sum_k_dl1k)) + 
                   (exp(omega)*C153/(1+exp(omega)*C151*C153))*(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))^(exp(alpha1)+sum_k_yi1k) * 
                   (exp(alpha1)*log(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))+(exp(alpha1)*(exp(alpha1)+sum_k_yi1k)*(exp_alpha1_exp_exp_beta1z_sum_k_dl1k)^(-2))/(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))))
  p_alpha2 = sum(exp(alpha2)*(digamma(exp(alpha2)+sum_k_yi2k)-digamma(exp(alpha2))+1+alpha2-log(exp(alpha2)+exp2*sum_k_dl2k)) + 
                   (exp(omega)*C151/(1+exp(omega)*C151*C153))*(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))^(exp(alpha2)+sum_k_yi2k) * 
                   (exp(alpha2)*log(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))+(exp(alpha2)*(exp(alpha2)+sum_k_yi2k)*(exp_alpha2_exp_exp_beta2z_sum_k_dl2k)^(-2))/(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))))
  
  p_beta10 = sum(exp(beta10)*1*sum_k_yi1k - (exp(alpha1)+sum_k_yi1k)*exp(beta10)*1*exp1*sum_k_dl1k/(exp_alpha1_exp_exp_beta1z_sum_k_dl1k-1) + 
                   (exp(omega)*C153/(1+exp(omega)*C151*C153))*(exp(alpha1)+sum_k_yi1k)*(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))^(exp(alpha1)+sum_k_yi1k-1) * 
                   (1*exp(beta10)*exp1*sum_k_dl1k*exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-2)))
  p_beta11 = sum(exp(beta11)*z1*sum_k_yi1k - (exp(alpha1)+sum_k_yi1k)*exp(beta11)*z1*exp1*sum_k_dl1k/(exp_alpha1_exp_exp_beta1z_sum_k_dl1k-1) + 
                   (exp(omega)*C153/(1+exp(omega)*C151*C153))*(exp(alpha1)+sum_k_yi1k)*(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))^(exp(alpha1)+sum_k_yi1k-1) * 
                   (z1*exp(beta11)*exp1*sum_k_dl1k*exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-2)))
  p_beta12 = sum(exp(beta12)*z2*sum_k_yi1k - (exp(alpha1)+sum_k_yi1k)*exp(beta12)*z2*exp1*sum_k_dl1k/(exp_alpha1_exp_exp_beta1z_sum_k_dl1k-1) + 
                   (exp(omega)*C153/(1+exp(omega)*C151*C153))*(exp(alpha1)+sum_k_yi1k)*(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))^(exp(alpha1)+sum_k_yi1k-1) * 
                   (z2*exp(beta12)*exp1*sum_k_dl1k*exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-2)))
  p_beta20 = sum(exp(beta20)*1*sum_k_yi2k - (exp(alpha2)+sum_k_yi2k)*exp(beta20)*1*exp2*sum_k_dl2k/(exp_alpha2_exp_exp_beta2z_sum_k_dl2k-1) + 
                   (exp(omega)*C151/(1+exp(omega)*C151*C153))*(exp(alpha2)+sum_k_yi2k)*(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))^(exp(alpha2)+sum_k_yi2k-1) * 
                   (1*exp(beta20)*exp2*sum_k_dl2k*exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-2)))
  p_beta21 = sum(exp(beta21)*z1*sum_k_yi2k - (exp(alpha2)+sum_k_yi2k)*exp(beta21)*z1*exp2*sum_k_dl2k/(exp_alpha2_exp_exp_beta2z_sum_k_dl2k-1) + 
                   (exp(omega)*C151/(1+exp(omega)*C151*C153))*(exp(alpha2)+sum_k_yi2k)*(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))^(exp(alpha2)+sum_k_yi2k-1) * 
                   (z1*exp(beta21)*exp2*sum_k_dl2k*exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-2)))
  p_beta22 = sum(exp(beta22)*z2*sum_k_yi2k - (exp(alpha2)+sum_k_yi2k)*exp(beta22)*z2*exp2*sum_k_dl2k/(exp_alpha2_exp_exp_beta2z_sum_k_dl2k-1) + 
                   (exp(omega)*C151/(1+exp(omega)*C151*C153))*(exp(alpha2)+sum_k_yi2k)*(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))^(exp(alpha2)+sum_k_yi2k-1) * 
                   (z2*exp(beta22)*exp2*sum_k_dl2k*exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-2)))
  
  p_dl11 = sum(yi11 - (exp(alpha1)+sum_k_yi1k)*exp(dl11)*exp1/(exp_alpha1_exp_exp_beta1z_sum_k_dl1k-1) +
                 (exp(omega)*C153/(1+exp(omega)*C151*C153)) * (exp(dl11)*exp1*(exp(alpha1)+sum_k_yi1k)*exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-2)*(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))^(exp(alpha1)+sum_k_yi1k-1)))
  p_dl12 = sum(yi12 - (exp(alpha1)+sum_k_yi1k)*exp(dl12)*exp1/(exp_alpha1_exp_exp_beta1z_sum_k_dl1k-1) +
                 (exp(omega)*C153/(1+exp(omega)*C151*C153)) * (exp(dl12)*exp1*(exp(alpha1)+sum_k_yi1k)*exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-2)*(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))^(exp(alpha1)+sum_k_yi1k-1)))
  p_dl13 = sum(yi13 - (exp(alpha1)+sum_k_yi1k)*exp(dl13)*exp1/(exp_alpha1_exp_exp_beta1z_sum_k_dl1k-1) +
                 (exp(omega)*C153/(1+exp(omega)*C151*C153)) * (exp(dl13)*exp1*(exp(alpha1)+sum_k_yi1k)*exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-2)*(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))^(exp(alpha1)+sum_k_yi1k-1)))
  p_dl14 = sum(yi14 - (exp(alpha1)+sum_k_yi1k)*exp(dl14)*exp1/(exp_alpha1_exp_exp_beta1z_sum_k_dl1k-1) +
                 (exp(omega)*C153/(1+exp(omega)*C151*C153)) * (exp(dl14)*exp1*(exp(alpha1)+sum_k_yi1k)*exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-2)*(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))^(exp(alpha1)+sum_k_yi1k-1)))
  p_dl15 = sum(yi15 - (exp(alpha1)+sum_k_yi1k)*exp(dl15)*exp1/(exp_alpha1_exp_exp_beta1z_sum_k_dl1k-1) +
                 (exp(omega)*C153/(1+exp(omega)*C151*C153)) * (exp(dl15)*exp1*(exp(alpha1)+sum_k_yi1k)*exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-2)*(1-exp_alpha1_exp_exp_beta1z_sum_k_dl1k^(-1))^(exp(alpha1)+sum_k_yi1k-1)))
  p_dl21 = sum(yi21 - (exp(alpha2)+sum_k_yi2k)*exp(dl21)*exp2/(exp_alpha2_exp_exp_beta2z_sum_k_dl2k-1) +
                 (exp(omega)*C151/(1+exp(omega)*C151*C153)) * (exp(dl21)*exp2*(exp(alpha2)+sum_k_yi2k)*exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-2)*(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))^(exp(alpha2)+sum_k_yi2k-1)))
  p_dl22 = sum(yi22 - (exp(alpha2)+sum_k_yi2k)*exp(dl22)*exp2/(exp_alpha2_exp_exp_beta2z_sum_k_dl2k-1) +
                 (exp(omega)*C151/(1+exp(omega)*C151*C153)) * (exp(dl22)*exp2*(exp(alpha2)+sum_k_yi2k)*exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-2)*(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))^(exp(alpha2)+sum_k_yi2k-1)))
  p_dl23 = sum(yi23 - (exp(alpha2)+sum_k_yi2k)*exp(dl23)*exp2/(exp_alpha2_exp_exp_beta2z_sum_k_dl2k-1) +
                 (exp(omega)*C151/(1+exp(omega)*C151*C153)) * (exp(dl23)*exp2*(exp(alpha2)+sum_k_yi2k)*exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-2)*(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))^(exp(alpha2)+sum_k_yi2k-1)))
  p_dl24 = sum(yi24 - (exp(alpha2)+sum_k_yi2k)*exp(dl24)*exp2/(exp_alpha2_exp_exp_beta2z_sum_k_dl2k-1) +
                 (exp(omega)*C151/(1+exp(omega)*C151*C153)) * (exp(dl24)*exp2*(exp(alpha2)+sum_k_yi2k)*exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-2)*(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))^(exp(alpha2)+sum_k_yi2k-1)))
  p_dl25 = sum(yi25 - (exp(alpha2)+sum_k_yi2k)*exp(dl25)*exp2/(exp_alpha2_exp_exp_beta2z_sum_k_dl2k-1) +
                 (exp(omega)*C151/(1+exp(omega)*C151*C153)) * (exp(dl25)*exp2*(exp(alpha2)+sum_k_yi2k)*exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-2)*(1-exp_alpha2_exp_exp_beta2z_sum_k_dl2k^(-1))^(exp(alpha2)+sum_k_yi2k-1)))
  
  p_omega = sum(exp(omega)*C151*C153/(1 + exp(omega) * C151 * C153))
  
  return (c(p_alpha1, p_beta10, p_beta11, p_beta12, p_dl11, p_dl12, p_dl13, p_dl14, p_dl15, p_alpha2, p_beta20, p_beta21, p_beta22, p_dl21, p_dl22, p_dl23, p_dl24, p_dl25, p_omega))
}


# Loglikelihood function normal mode
logL = function(par, data){
  alpha1 = par[1]
  beta10 = par[2]; beta11 = par[3]; beta12 = par[4]
  dl11 = par[5]; dl12 = par[6]; dl13 = par[7]; dl14 = par[8]; dl15 = par[9]
  alpha2 = par[10]
  beta20 = par[11]; beta21 = par[12]; beta22 = par[13]
  dl21 = par[14]; dl22 = par[15]; dl23 = par[16]; dl24 = par[17]; dl25 = par[18]
  omega = par[19]
  
  yi11 = data[, 1]; yi12 = data[, 2]; yi13 = data[, 3]; yi14 = data[, 4]; yi15 = data[, 5]
  yi21 = data[, 6]; yi22 = data[, 7]; yi23 = data[, 8]; yi24 = data[, 9]; yi25 = data[, 10]
  z1 = data[, 11]; z2 = data[, 12]
  
  sum_k_yi1k = yi11 + yi12 + yi13 + yi14 + yi15
  sum_k_yi2k = yi21 + yi22 + yi23 + yi24 + yi25
  sum_k_dl1k = dl11 + dl12 + dl13 + dl14 + dl15
  sum_k_dl2k = dl21 + dl22 + dl23 + dl24 + dl25
  
  beta1_z = beta10 + beta11 * z1 + beta12 * z2
  beta2_z = beta20 + beta21 * z1 + beta22 * z2
  exp_beta1_z = exp(beta10 + beta11 * z1 + beta12 * z2)
  exp_beta2_z = exp(beta20 + beta21 * z1 + beta22 * z2)
  
  c1 = yi11*log(dl11) + yi12*log(dl12) + yi13*log(dl13) + yi14*log(dl14) + yi15*log(dl15)
  c2 = lgamma(alpha1 + sum_k_yi1k)
  c3 = lgamma(alpha1)
  c4 = lgamma(sum_k_yi1k + 1)
  c5 = alpha1 * log(alpha1)
  c6 = beta1_z * sum_k_yi1k
  c7 = (alpha1 + sum_k_yi1k) * log(alpha1 + exp_beta1_z * sum_k_dl1k)
  
  c8 = yi21*log(dl21) + yi22*log(dl22) + yi23*log(dl23) + yi24*log(dl24) + yi25*log(dl25)
  c9 = lgamma(alpha2 + sum_k_yi2k)
  c10 = lgamma(alpha2)
  c11 = lgamma(sum_k_yi2k + 1)
  c12 = alpha2 * log(alpha2)
  c13 = beta2_z * sum_k_yi2k
  c14 = (alpha2 + sum_k_yi2k) * log(alpha2 + exp_beta2_z * sum_k_dl2k)
  
  c151 = ((alpha1 + exp_beta1_z * sum_k_dl1k)/(alpha1 + exp_beta1_z * sum_k_dl1k + 1))^(alpha1 + sum_k_yi1k) - (alpha1/(alpha1 + 1))^(alpha1)
  c152 = ((alpha2 + exp_beta2_z * sum_k_dl2k)/(alpha2 + exp_beta2_z * sum_k_dl2k + 1))^(alpha2 + sum_k_yi2k) - (alpha2/(alpha2 + 1))^(alpha2)
  c15 = log(1 + omega * c151 * c152)
  
  return (-sum(c1 + c2 - c3 - c4 + c5 + c6 - c7 + 
                 c8 + c9 - c10 - c11 + c12 + c13 - c14 + c15))
}

score = function(par, data){
  alpha1 = par[1]
  beta10 = par[2]; beta11 = par[3]; beta12 = par[4]
  dl11 = par[5]; dl12 = par[6]; dl13 = par[7]; dl14 = par[8]; dl15 = par[9]
  alpha2 = par[10]
  beta20 = par[11]; beta21 = par[12]; beta22 = par[13]
  dl21 = par[14]; dl22 = par[15]; dl23 = par[16]; dl24 = par[17]; dl25 = par[18]
  omega = par[19]
  
  yi11 = data[, 1]; yi12 = data[, 2]; yi13 = data[, 3]; yi14 = data[, 4]; yi15 = data[, 5]
  yi21 = data[, 6]; yi22 = data[, 7]; yi23 = data[, 8]; yi24 = data[, 9]; yi25 = data[, 10]
  z1 = data[, 11]; z2 = data[, 12]
  
  sum_k_yi1k = yi11 + yi12 + yi13 + yi14 + yi15
  sum_k_yi2k = yi21 + yi22 + yi23 + yi24 + yi25
  sum_k_dl1k = dl11 + dl12 + dl13 + dl14 + dl15
  sum_k_dl2k = dl21 + dl22 + dl23 + dl24 + dl25
  beta1_z = beta10 + beta11 * z1 + beta12 * z2
  beta2_z = beta20 + beta21 * z1 + beta22 * z2
  exp_beta1_z = exp(beta10 + beta11 * z1 + beta12 * z2)
  exp_beta2_z = exp(beta20 + beta21 * z1 + beta22 * z2)
  denominator1 = alpha1 + exp_beta1_z * sum_k_dl1k + 1
  denominator2 = alpha2 + exp_beta2_z * sum_k_dl2k + 1
  c151 = ((alpha1 + exp_beta1_z * sum_k_dl1k)/(alpha1 + exp_beta1_z * sum_k_dl1k + 1))^(alpha1 + sum_k_yi1k) - (alpha1/(alpha1 + 1))^(alpha1)
  c152 = ((alpha2 + exp_beta2_z * sum_k_dl2k)/(alpha2 + exp_beta2_z * sum_k_dl2k + 1))^(alpha2 + sum_k_yi2k) - (alpha2/(alpha2 + 1))^(alpha2)
  
  p_alpha1 = sum(digamma(alpha1 + sum_k_yi1k) - digamma(alpha1) + 1 + log(alpha1) - log(alpha1 + exp_beta1_z * sum_k_dl1k) - (alpha1 + sum_k_yi1k)/(alpha1 + exp_beta1_z * sum_k_dl1k) + 
                   (omega*c152/(1+omega*c151*c152)) * (1 - denominator1^(-1))^(alpha1 + sum_k_yi1k) * (log(1 - denominator1^(-1)) + (alpha1 + sum_k_yi1k)*denominator1^(-2)/(1-denominator1^(-1))))
  p_alpha2 = sum(digamma(alpha2 + sum_k_yi2k) - digamma(alpha2) + 1 + log(alpha2) - log(alpha2 + exp_beta2_z * sum_k_dl2k) - (alpha2 + sum_k_yi2k)/(alpha2 + exp_beta2_z * sum_k_dl2k) + 
                   (omega*c151/(1+omega*c151*c152)) * (1 - denominator2^(-1))^(alpha2 + sum_k_yi2k) * (log(1 - denominator2^(-1)) + (alpha2 + sum_k_yi2k)*denominator2^(-2)/(1-denominator2^(-1))))
  
  p_beta10 = sum(1 * sum_k_yi1k - (alpha1 + sum_k_yi1k) * 1 * exp_beta1_z * sum_k_dl1k / (alpha1 + exp_beta1_z * sum_k_dl1k) + 
                   (omega*c152/(1+omega*c151*c152)) * 1 * exp_beta1_z * (alpha1 + sum_k_yi1k) * denominator1^(-2) * (1 - (denominator1)^(-1))^(alpha1 + sum_k_yi1k - 1))
  p_beta11 = sum(z1 * sum_k_yi1k - (alpha1 + sum_k_yi1k) * z1 * exp_beta1_z * sum_k_dl1k / (alpha1 + exp_beta1_z * sum_k_dl1k) + 
                   (omega*c152/(1+omega*c151*c152)) * z1 * exp_beta1_z * (alpha1 + sum_k_yi1k) * denominator1^(-2) * (1 - (denominator1)^(-1))^(alpha1 + sum_k_yi1k - 1))
  p_beta12 = sum(z2 * sum_k_yi1k - (alpha1 + sum_k_yi1k) * z2 * exp_beta1_z * sum_k_dl1k / (alpha1 + exp_beta1_z * sum_k_dl1k) + 
                   (omega*c152/(1+omega*c151*c152)) * z2 * exp_beta1_z * (alpha1 + sum_k_yi1k) * denominator1^(-2) * (1 - (denominator1)^(-1))^(alpha1 + sum_k_yi1k - 1))
  p_beta20 = sum(1 * sum_k_yi2k - (alpha2 + sum_k_yi2k) * 1 * exp_beta2_z * sum_k_dl2k / (alpha2 + exp_beta2_z * sum_k_dl2k) + 
                   (omega*c151/(1+omega*c151*c152)) * 1 * exp_beta2_z * (alpha2 + sum_k_yi2k) * denominator2^(-2) * (1 - (denominator2)^(-1))^(alpha2 + sum_k_yi2k - 1))
  p_beta21 = sum(z1 * sum_k_yi2k - (alpha2 + sum_k_yi2k) * z1 * exp_beta2_z * sum_k_dl2k / (alpha2 + exp_beta2_z * sum_k_dl2k) + 
                   (omega*c151/(1+omega*c151*c152)) * z1 * exp_beta2_z * (alpha2 + sum_k_yi2k) * denominator2^(-2) * (1 - (denominator2)^(-1))^(alpha2 + sum_k_yi2k - 1))
  p_beta22 = sum(z2 * sum_k_yi2k - (alpha2 + sum_k_yi2k) * z2 * exp_beta2_z * sum_k_dl2k / (alpha2 + exp_beta2_z * sum_k_dl2k) + 
                   (omega*c151/(1+omega*c151*c152)) * z2 * exp_beta2_z * (alpha2 + sum_k_yi2k) * denominator2^(-2) * (1 - (denominator2)^(-1))^(alpha2 + sum_k_yi2k - 1))
  
  p_dl11 = sum((yi11/dl11) - (alpha1 + sum_k_yi1k) * exp_beta1_z / (denominator1 - 1) + 
                 (omega*c152/(1+omega*c151*c152)) * exp_beta1_z * (alpha1 + sum_k_yi1k) * denominator1^(-2) * (1 - denominator1^(-1))^(alpha1 + sum_k_yi1k - 1))
  p_dl12 = sum((yi12/dl12) - (alpha1 + sum_k_yi1k) * exp_beta1_z / (denominator1 - 1) + 
                 (omega*c152/(1+omega*c151*c152)) * exp_beta1_z * (alpha1 + sum_k_yi1k) * denominator1^(-2) * (1 - denominator1^(-1))^(alpha1 + sum_k_yi1k - 1))
  p_dl13 = sum((yi13/dl13) - (alpha1 + sum_k_yi1k) * exp_beta1_z / (denominator1 - 1) + 
                 (omega*c152/(1+omega*c151*c152)) * exp_beta1_z * (alpha1 + sum_k_yi1k) * denominator1^(-2) * (1 - denominator1^(-1))^(alpha1 + sum_k_yi1k - 1))
  p_dl14 = sum((yi14/dl14) - (alpha1 + sum_k_yi1k) * exp_beta1_z / (denominator1 - 1) + 
                 (omega*c152/(1+omega*c151*c152)) * exp_beta1_z * (alpha1 + sum_k_yi1k) * denominator1^(-2) * (1 - denominator1^(-1))^(alpha1 + sum_k_yi1k - 1))
  p_dl15 = sum((yi15/dl15) - (alpha1 + sum_k_yi1k) * exp_beta1_z / (denominator1 - 1) + 
                 (omega*c152/(1+omega*c151*c152)) * exp_beta1_z * (alpha1 + sum_k_yi1k) * denominator1^(-2) * (1 - denominator1^(-1))^(alpha1 + sum_k_yi1k - 1))
  p_dl21 = sum((yi21/dl21) - (alpha2 + sum_k_yi2k) * exp_beta2_z / (denominator2 - 1) + 
                 (omega*c151/(1+omega*c151*c152)) * exp_beta2_z / (alpha2 + sum_k_yi2k) * denominator2^(-2) * (1 - denominator2^(-2))^(alpha2 + sum_k_yi2k - 1))
  p_dl22 = sum((yi22/dl22) - (alpha2 + sum_k_yi2k) * exp_beta2_z / (denominator2 - 1) + 
                 (omega*c151/(1+omega*c151*c152)) * exp_beta2_z / (alpha2 + sum_k_yi2k) * denominator2^(-2) * (1 - denominator2^(-2))^(alpha2 + sum_k_yi2k - 1))
  p_dl23 = sum((yi23/dl23) - (alpha2 + sum_k_yi2k) * exp_beta2_z / (denominator2 - 1) + 
                 (omega*c151/(1+omega*c151*c152)) * exp_beta2_z / (alpha2 + sum_k_yi2k) * denominator2^(-2) * (1 - denominator2^(-2))^(alpha2 + sum_k_yi2k - 1))
  p_dl24 = sum((yi24/dl24) - (alpha2 + sum_k_yi2k) * exp_beta2_z / (denominator2 - 1) + 
                 (omega*c151/(1+omega*c151*c152)) * exp_beta2_z / (alpha2 + sum_k_yi2k) * denominator2^(-2) * (1 - denominator2^(-2))^(alpha2 + sum_k_yi2k - 1))
  p_dl25 = sum((yi25/dl25) - (alpha2 + sum_k_yi2k) * exp_beta2_z / (denominator2 - 1) + 
                 (omega*c151/(1+omega*c151*c152)) * exp_beta2_z / (alpha2 + sum_k_yi2k) * denominator2^(-2) * (1 - denominator2^(-2))^(alpha2 + sum_k_yi2k - 1))
  
  p_omega = sum(c151*c152/(1+omega*c151*c152))
  
  return (c(p_alpha1, p_alpha2, 
            p_beta10, p_beta11, p_beta12, 
            p_beta20, p_beta21, p_beta22,
            p_dl11, p_dl12, p_dl13, p_dl14, p_dl15,
            p_dl21, p_dl22, p_dl23, p_dl24, p_dl25,
            p_omega))
}