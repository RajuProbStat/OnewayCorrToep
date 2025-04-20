################################################################################
##The LRT T_2 and asymptotic LRT T_2A are developed for testing the equality of mean treatment effects in a onw-way ANOVA with Toeplitz error structure
##function to evaluate the observed, critical and P-values of the LRT test T_2 and the asymptotic LRT test T_2A
#' Function to obtain the observed, critical and P-values of the likelihood ratio and asymptotic likelihood ratio tests for testing homogeneity of mean effects in a one-way correlated ANOVA with a special Toeplitz error structure
#'
#' @param data A numeric data frame consists of cell data sets in the one-way ANOVA
#' @param size The level of significance
#' @param t The number of time (space) points
#' @param k The number of groups
#' @param n A numeric vector consists of cell sample sizes in the one-way ANOVA
#'
#' @return Numeric values consist of observed, critical and P-values of the likelihood ratio and asymptotic likelihood ratio tests
#' @export
#'
#' @examples
#' data<-read.csv("https://raw.githubusercontent.com/RajuProbStat/OnewayCorrToep/refs/heads/main/DBP_data.csv")
#' T_2_T_2_A(data,0.05,3,3,c(6,6,8))
T_2_T_2_A<-function(data,size,t,k,n){
  library(foreach)
  library(doParallel)
  library(MASS)
  library(matlib)
  library(stats)
  data<-data #data consists of cell data sets
  size<-size #the level of significance
  k<-k #number of treatments/groups
  t<-t #number of time (space) points
  n<-n #vector consists of cell sample sizes
  x<-array(NA,dim=c(t,max(n),k))  #three dimensional array to store all the cell data sets
  n_1<-c(0,n) #require to define sequence of lengths of the individual cell data sets
  n_limit<-c(0,n[1]) #initialization
  for(i in 1:k)
  {
    sum<-0
    for(j in 1:i)
    {
      sum<-sum+n_1[j]
    }
    n_limit<-sum+c(1,n[i])
    x[,1:n[i],i]<-data.matrix(data[,n_limit[1]:n_limit[2]]) #store the individual cell data sets
  }
  ##define the orthogonal matrix P
  P<-array(NA,dim=c(t,t))
  for(i in 1:t)
  {
    for(j in 1:t)
    {
      P[i,j]<-sqrt(2/(t+1))*sin((pi*i*j)/(t+1))
    }
  }
  ##define the inverse of the diagonal matrix
  D_inv_fun<-function(rho){
    D_inv<-array(NA,dim=c(t,t))
    for(i in 1:t)
    {
      for(j in 1:t)
      {
        if(i==j)
        {
          D_inv[i,j]<-(1/(1+2*rho*cos((i*pi)/(t+1))))
        }
        else
        {
          D_inv[i,j]<-0
        }
      }
    }
    return(D_inv)
  }
  ##define the derivative of the inverse of the diagonal matrix
  D_diff_inv_fun<-function(rho){
    D_diff_inv<-array(NA,dim=c(t,t))
    for(i in 1:t)
    {
      for(j in 1:t)
      {
        if(i==j)
        {
          D_diff_inv[i,j]<- -(2*cos((i*pi)/(t+1)))/(1+2*rho*cos((i*pi)/(t+1)))^2
        }
        else
        {
          D_diff_inv[i,j]<-0
        }
      }
    }
    return(D_diff_inv)
  }
  ##constraint matrix for constrOptim function under the null space
  cons_mat_null<-array(NA,dim=c(3*k,(2*k+1)))
  for(i in 1:k)
  {
    for(j in 1:(2*k+1))
    {
      if(j==(i+1))
      {
        cons_mat_null[i,j]<-1
      }
      else
      {
        cons_mat_null[i,j]<-0
      }
    }
  }
  for(i in seq((k+1),(3*k-1),2))
  {
    for(j in 1:(2*k+1))
    {
      if(j==((i+k+3)/2))
      {
        cons_mat_null[i,j]<-1
      }
      else
      {
        cons_mat_null[i,j]<-0
      }
    }
  }
  for(i in seq((k+2),(3*k),2))
  {
    cons_mat_null[i,]<--cons_mat_null[i-1,]
  }
  ##constraint matrix for constrOptim function under the full space
  cons_mat_full<-array(NA,dim=c(3*k,(2*k+k)))
  for(i in 1:k)
  {
    for(j in 1:(2*k+k))
    {
      if(j==(i+k))
      {
        cons_mat_full[i,j]<-1
      }
      else
      {
        cons_mat_full[i,j]<-0
      }
    }
  }
  for(i in seq((k+1),(3*k-1),2))
  {
    for(j in 1:(2*k+k))
    {
      if(j==((i+3*k+1)/2))
      {
        cons_mat_full[i,j]<-1
      }
      else
      {
        cons_mat_full[i,j]<-0
      }
    }
  }
  for(i in seq((k+2),(3*k),2))
  {
    cons_mat_full[i,]<--cons_mat_full[i-1,]
  }
  ##vector for constrOptim function
  cons_vec<-rep(NA,3*k)
  for(i in 1:k)
  {
    cons_vec[i]<-0
  }
  for(i in (k+1):(3*k))
  {
    cons_vec[i]<-(-1/2)
  }
  ##orthogonal transformation y=Px, transform data y
  y<-array(NA,dim=c(t,max(n),k))
  for(l in 1:k)
  {
    y[,1:n[l],l]<-P%*%x[,1:n[l],l]
  }
  w<-P%*%rep(1,t)
  #### MLEs of the parameters under the null parameter space
  ##the log-likelihood function under the null parameter space
  frr_null<-function(x){
    mu<-x[1]
    v<-x[2:(k+1)]
    rho<-x[(k+2):(2*k+1)]
    sum1<-0
    for(i in 1:k)
    {
      for(j in 1:n[i])
      {
        sum1<-sum1+(1/(2*v[i]))*(t(y[,j,i]-mu*w)%*%D_inv_fun(rho[i])%*%(y[,j,i]-mu*w))
      }
    }
    sum2<-0
    for(i in 1:k)
    {
      sum2<-sum2+((n[i]*t)/2)*log(v[i])-(n[i]/2)*log(det(D_inv_fun(rho[i])))
    }

    ((sum(n)*t)/2)*log(2*pi)+sum2+sum1
  }
  ##gradient of the log-likelihood function under the null parameter space
  grr_null<-function(x){
    mu<-x[1]
    v<-x[2:(k+1)]
    rho<-x[(k+2):(2*k+1)]

    der<-rep(0,(2*k+1))

    for(i in 1:k)
    {
      der[1]<-der[1]-(n[i]/v[i])*((t((apply(y[,1:n[i],i],1,mean))-mu*w))%*%D_inv_fun(rho[i])%*%w)
    }
    for(r in 2:(k+1))
    {
      der[r]<-((n[r-1]*t)/(2*v[r-1]))
      for(j in 1:n[r-1])
      {
        der[r]<-der[r]-(1/(2*(v[r-1])^2))*(t(y[,j,r-1]-mu*w)%*%D_inv_fun(rho[r-1])%*%(y[,j,r-1]-mu*w))
      }
    }
    for(r in (k+2):(2*k+1))
    {
      for(s in 1:t)
      {
        der[r]<-der[r]+(n[r-k-1]/2)*((2*cos((s*pi)/(t+1)))/(1+2*rho[r-k-1]*cos((s*pi)/(t+1))))
      }
      for(j in 1:n[r-k-1])
      {
        der[r]<-der[r]+(1/(2*(v[r-k-1])))*(t(y[,j,r-k-1]-mu*w)%*%D_diff_inv_fun(rho[r-k-1])%*%(y[,j,r-k-1]-mu*w))
      }
    }
    der
  }
  ##initialization of the parameters under the null parameter space
  cor_vec<-rep(NA,k)
  for(i in 1:k)
  {
    cor_vec[i]<-cor(x[,1,i],x[,2,i])
    repeat
    {
      if(abs(cor_vec[i])<0.5)
      {
        break
      }
      cor_vec[i]<-(0.8)*cor_vec[i]
    }
  }
  var_vec<-rep(NA,k)
  for(i in 1:k)
  {
    var_vec[i]<-var(x[,1,i])
  }

  mu_vec<-rep(NA,k)
  for(i in 1:k)
  {
    mu_vec[i]<-mean(x[,1:n[i],i])
  }
  par_initial_null<-c(mean(mu_vec),var_vec,cor_vec)
  par_null<-constrOptim(par_initial_null,frr_null,grr_null,ui=cons_mat_null,ci=cons_vec, method="BFGS")$par
  value_null<-frr_null(par_null) #value of the log-likelihood function under the null parameter space
  #### MLEs of the parameters under the full parameter space
  ##the log-likelihood function under the full parameter space
  frr_full<-function(x){
    mu<-x[1:k]
    v<-x[(k+1):(2*k)]
    rho<-x[(2*k+1):(3*k)]
    sum1<-0
    for(i in 1:k)
    {
      for(j in 1:n[i])
      {
        sum1<-sum1+(1/(2*v[i]))*(t(y[,j,i]-mu[i]*w)%*%D_inv_fun(rho[i])%*%(y[,j,i]-mu[i]*w))
      }
    }
    sum2<-0
    for(i in 1:k)
    {
      sum2<-sum2+((n[i]*t)/2)*log(v[i])-(n[i]/2)*log(det(D_inv_fun(rho[i])))
    }
    ((sum(n)*t)/2)*log(2*pi)+sum2+sum1
  }
  ##gradient of the log-likelihood function under the full parameter space
  grr_full<-function(x){
    mu<-x[1:k]
    v<-x[(k+1):(2*k)]
    rho<-x[(2*k+1):(3*k)]
    der<-rep(0,(3*k))
    for(r in 1:k)
    {
      der[r]<-der[r]-(n[r]/v[r])*((t((apply(y[,1:n[r],r],1,mean))-mu[r]*w))%*%D_inv_fun(rho[r])%*%w)
    }
    for(r in (k+1):(2*k))
    {
      der[r]<-((n[r-k]*t)/(2*v[r-k]))
      for(j in 1:n[r-k])
      {
        der[r]<-der[r]-(1/(2*(v[r-k])^2))*(t(y[,j,r-k]-mu[r-k]*w)%*%D_inv_fun(rho[r-k])%*%(y[,j,r-k]-mu[r-k]*w))
      }
    }
    for(r in (2*k+1):(3*k))
    {
      for(s in 1:t)
      {
        der[r]<-der[r]+(n[r-2*k]/2)*((2*cos((s*pi)/(t+1)))/(1+2*rho[r-2*k]*cos((s*pi)/(t+1))))
      }
      for(j in 1:n[r-2*k])
      {
        der[r]<-der[r]+(1/(2*(v[r-2*k])))*(t(y[,j,r-2*k]-mu[r-2*k]*w)%*%D_diff_inv_fun(rho[r-2*k])%*%(y[,j,r-2*k]-mu[r-2*k]*w))
      }
    }
    der
  }
  ##initialization of the parameters under the full parameter space
  par_initial_full<-c(mu_vec,var_vec,cor_vec)
  par_full<-constrOptim(par_initial_full,frr_full,grr_full,ui=cons_mat_full,ci=cons_vec, method="BFGS")$par #the optimal parameters
  value_full<-frr_full(par_full) #the value of the log-likelihood function at the MLEs under the full parameter space
  LRT_value<-exp(value_full-value_null) #the likelihood ratio test statistic value
  values<-c(LRT_value,par_full[(k+1):(3*k)]) #contains the LRT values and the MLEs of the parameters under the full space
  LRT_value<-values[1] #LRT value

  ##bootstrap step
  cl<-makeCluster(4)    #number of cores used in the parallelization
  registerDoParallel(cl)
  B<-1000               #the number of bootstrap replications
  start_time <- Sys.time()
  LRT_value_boot<-foreach(b=1:B, .combine=rbind)%dopar%{
    library(foreach)
    library(doParallel)
    library(MASS)
    library(matlib)
    library(stats)
  ##function to evaluate the LRT value for parametric bootstrap samples
  corr_LRT<-function(k,t_data,n_data,mu_data,variance_data,rho_data){
    k<-k
    t<-t_data
    n<-n_data
    mu<-mu_data
    var<-variance_data
    rho<-rho_data
    ##function for tridiagonal banded covariance matrix
    cov_fun<-function(t,rho,v)
    {
      Sigma<-array(NA,dim=c(t,t))
      for(i in 1:t)
      {
        for(j in 1:t)
        {
          if(i==j)
          {
            Sigma[i,j]<-v
          }
          else
          {
            if(abs(i-j)==1)
            {
              Sigma[i,j]<-rho*v
            }
            else
            {
              Sigma[i,j]<-0
            }
          }
        }
      }
      return(Sigma)
    }
    ##data generation from the multivariate normal distributions
    x<-array(NA,dim=c(t,max(n),k))
    for(l in 1:k)
    {
      for(j in 1:n[l])
      {
        x[,j,l]<-mvrnorm(1,mu[l]*rep(1,t),cov_fun(t,rho[l],var[l]))
      }
    }
    ##transform data y=Px
    y<-array(NA,dim=c(t,max(n),k))
    for(l in 1:k)
    {
      y[,1:n[l],l]<-P%*%x[,1:n[l],l]
    }
    w<-P%*%rep(1,t)
    ####MLEs of the parameters under the null parameter space
    ##the log-likelihood function under the null parameter space
    frr_null<-function(x){
      mu<-x[1]
      v<-x[2:(k+1)]
      rho<-x[(k+2):(2*k+1)]
      sum1<-0
      for(i in 1:k)
      {
        for(j in 1:n[i])
        {
          sum1<-sum1+(1/(2*v[i]))*(t(y[,j,i]-mu*w)%*%D_inv_fun(rho[i])%*%(y[,j,i]-mu*w))
        }
      }
      sum2<-0
      for(i in 1:k)
      {
        sum2<-sum2+((n[i]*t)/2)*log(v[i])-(n[i]/2)*log(det(D_inv_fun(rho[i])))
      }
      ((sum(n)*t)/2)*log(2*pi)+sum2+sum1
    }
    ##gradient of the log-likelihood function under the null parameter space
    grr_null<-function(x){
      mu<-x[1]
      v<-x[2:(k+1)]
      rho<-x[(k+2):(2*k+1)]
      der<-rep(0,(2*k+1))
      for(i in 1:k)
      {
        der[1]<-der[1]-(n[i]/v[i])*((t((apply(y[,1:n[i],i],1,mean))-mu*w))%*%D_inv_fun(rho[i])%*%w)
      }
      for(r in 2:(k+1))
      {
        der[r]<-((n[r-1]*t)/(2*v[r-1]))
        for(j in 1:n[r-1])
        {
          der[r]<-der[r]-(1/(2*(v[r-1])^2))*(t(y[,j,r-1]-mu*w)%*%D_inv_fun(rho[r-1])%*%(y[,j,r-1]-mu*w))
        }
      }
      for(r in (k+2):(2*k+1))
      {
        for(s in 1:t)
        {
          der[r]<-der[r]+(n[r-k-1]/2)*((2*cos((s*pi)/(t+1)))/(1+2*rho[r-k-1]*cos((s*pi)/(t+1))))
        }
        for(j in 1:n[r-k-1])
        {
          der[r]<-der[r]+(1/(2*(v[r-k-1])))*(t(y[,j,r-k-1]-mu*w)%*%D_diff_inv_fun(rho[r-k-1])%*%(y[,j,r-k-1]-mu*w))
        }
      }
      der
    }
    ##initialization of the parameters under the null parameter space
    cor_vec<-rep(NA,k)
    for(i in 1:k)
    {
      cor_vec[i]<-cor(x[,1,i],x[,2,i])
      repeat
      {
        if(abs(cor_vec[i])<0.5)
        {
          break
        }
        cor_vec[i]<-(0.8)*cor_vec[i]
      }
    }
    var_vec<-rep(NA,k)
    for(i in 1:k)
    {
      var_vec[i]<-var(x[,1,i])
    }
    mu_vec<-rep(NA,k)
    for(i in 1:k)
    {
      mu_vec[i]<-mean(x[,1:n[i],i])
    }
    par_initial_null<-c(mean(mu_vec),var_vec,cor_vec)
    par_null<-constrOptim(par_initial_null,frr_null,grr_null,ui=cons_mat_null,ci=cons_vec, method="BFGS")$par
    value_null<-frr_null(par_null)
    #### MLEs of the parameters under the full parameter space
    ##the log-likelihood function under the full parameter space
    frr_full<-function(x){
      mu<-x[1:k]
      v<-x[(k+1):(2*k)]
      rho<-x[(2*k+1):(3*k)]
      sum1<-0
      for(i in 1:k)
      {
        for(j in 1:n[i])
        {
          sum1<-sum1+(1/(2*v[i]))*(t(y[,j,i]-mu[i]*w)%*%D_inv_fun(rho[i])%*%(y[,j,i]-mu[i]*w))
        }
      }
      sum2<-0
      for(i in 1:k)
      {
        sum2<-sum2+((n[i]*t)/2)*log(v[i])-(n[i]/2)*log(det(D_inv_fun(rho[i])))
      }
      ((sum(n)*t)/2)*log(2*pi)+sum2+sum1
    }
    ##gradient of the log-likelihood function under the full parameter space
    grr_full<-function(x){
      mu<-x[1:k]
      v<-x[(k+1):(2*k)]
      rho<-x[(2*k+1):(3*k)]
      der<-rep(0,(3*k))
      for(r in 1:k)
      {
        der[r]<-der[r]-(n[r]/v[r])*((t((apply(y[,1:n[r],r],1,mean))-mu[r]*w))%*%D_inv_fun(rho[r])%*%w)
      }
      for(r in (k+1):(2*k))
      {
        der[r]<-((n[r-k]*t)/(2*v[r-k]))
        for(j in 1:n[r-k])
        {
          der[r]<-der[r]-(1/(2*(v[r-k])^2))*(t(y[,j,r-k]-mu[r-k]*w)%*%D_inv_fun(rho[r-k])%*%(y[,j,r-k]-mu[r-k]*w))
        }
      }
      for(r in (2*k+1):(3*k))
      {
        for(s in 1:t)
        {
          der[r]<-der[r]+(n[r-2*k]/2)*((2*cos((s*pi)/(t+1)))/(1+2*rho[r-2*k]*cos((s*pi)/(t+1))))
        }
        for(j in 1:n[r-2*k])
        {
          der[r]<-der[r]+(1/(2*(v[r-2*k])))*(t(y[,j,r-2*k]-mu[r-2*k]*w)%*%D_diff_inv_fun(rho[r-2*k])%*%(y[,j,r-2*k]-mu[r-2*k]*w))
        }
      }
      der
    }
    par_initial_full<-c(mu_vec,var_vec,cor_vec)
    par_full<-constrOptim(par_initial_full,frr_full,grr_full,ui=cons_mat_full,ci=cons_vec, method="BFGS")$par
    value_full<-frr_full(par_full)
    LRT_value<-exp(value_full-value_null)
    return(c(LRT_value,par_full[(k+1):(3*k)])) # returns the LRT value as well as the parameters under full space cell variances and correlations
  }
  ##for bootstrap critical point
  k<-k
  t_data<-t
  n_data<-n
  set.seed(17*b)
  corr_LRT(k,t_data,n_data,rep(0,k),values[2:(k+1)],values[(k+2):(2*k+1)])[1]
  }
  critical_boot<-quantile(LRT_value_boot,size,names=FALSE) #bootstrap critical value of the LRT
  cat("The observed value, critical value and the P-value of the LRT test statistic are\n",LRT_value,critical_boot,mean(LRT_value>LRT_value_boot))
  if(LRT_value<critical_boot)
  {
    cat("\n The LRT test rejects the null hypothesis")
  }
  else
  {
    cat("\n The LRT test do not rejects the null hypothesis")
  }
  critical_value_asym<-qchisq(p =1-size ,df =(k-1) ,ncp = ,lower.tail = TRUE) #critical value of the asymptotic LRT T_2A
  cat("\n The observed value, critical value and the P-value of the asymptotic LRT test statistic are\n",-2*log(LRT_value),critical_value_asym,1-pchisq(q = -2*log(LRT_value),df = (k-1),ncp = 0))
  if(-2*log(LRT_value)>critical_value_asym)
  {
    cat("\n The asymptotic LRT test rejects the null hypothesis")
  }
  else
  {
    cat("\n The asymptotic LRT test do not rejects the null hypothesis")
  }
}
################################################################################





