################################################################################
##Simulatneous confidence intervals of the pairwise differences of the mean effects
#' Function to evaluate the simultaneous confidence intervals of the pairwise differences of the mean effects
#'
#' @param data A numeric data frame consists of cell data sets in the one-way ANOVA
#' @param size The level of significance
#' @param t The number of time (space) points
#' @param k The number of groups
#' @param n A numeric vector consists of cell sample sizes in the one-way ANOVA
#'
#' @return Simultaneous confidence intervals of the pairwise differences of the mean treatment effects
#' @export
#'
#' @examples
#' data<-read.csv("https://raw.githubusercontent.com/RajuProbStat/OnewayCorrToep/refs/heads/main/DBP_data.csv")
#' SCI_ANOVA_Toep(data,0.05,3,3,c(6,6,8))
SCI_ANOVA_Toep<-function(data,size,t,k,n){
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
  values<-c(par_full[1:(3*k)]) #contains the LRT values and the MLEs of the parameters under the full space
  #LRT_value<-values[1] #LRT value
  mu_hat_full<-values[1:k]
  var_hat_full<-values[(k+1):(2*k)]
  rho_hat_full<-values[(2*k+1):(3*k)]

  ##for the value of psi(r,s)
  psi<-array(0,dim=c(k-1,k))
  for(r in 1:(k-1))
  {
    for(s in (r+1):k)
    {
      psi[r,s]<-(1/(n[r]*t))*(t(rep(1,t))%*%(apply(x[,1:n[r],r],1,sum)))-(1/(n[s]*t))*(t(rep(1,t))%*%(apply(x[,1:n[s],s],1,sum)))
    }
  }
  var_psi<-array(0,dim=c(k-1,k))
  for(r in 1:(k-1))
  {
    for(s in (r+1):k)
    {
      var_psi[r,s]<-(1/t^2)*((var_hat_full[r]/n[r])*(t+2*(t-1)*rho_hat_full[r])+(var_hat_full[s]/n[s])*(t+2*(t-1)*rho_hat_full[s]))
    }
  }
  ## parametric bootstrap step
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
      ##for the initial values
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
      values_boot<-c(par_full[1:(3*k)]) #contains the LRT values and the MLEs of the parameters under the full space
      #LRT_value<-values[1] #LRT value
      mu_hat_full_boot<-values_boot[1:k]
      var_hat_full_boot<-values_boot[(k+1):(2*k)]
      rho_hat_full_boot<-values_boot[(2*k+1):(3*k)]
      ####
      psi_boot<-array(0,dim=c(k-1,k))
      for(r in 1:(k-1))
      {
        for(s in (r+1):k)
        {
          psi_boot[r,s]<-(1/(n[r]*t))*(t(rep(1,t))%*%(apply(x[,1:n[r],r],1,sum)))-(1/(n[s]*t))*(t(rep(1,t))%*%(apply(x[,1:n[s],s],1,sum)))
        }
      }
      var_psi_boot<-array(NA,dim=c(k-1,k))
      for(r in 1:(k-1))
      {
        for(s in (r+1):k)
        {
          var_psi_boot[r,s]<-(1/t^2)*((var_hat_full_boot[r]/n[r])*(t+2*(t-1)*rho_hat_full_boot[r])+(var_hat_full_boot[s]/n[s])*(t+2*(t-1)*rho_hat_full_boot[s]))
        }
      }
      U<-array(0,dim=c((k-1),k))
      for(r in 1:(k-1))
      {
        for(s in (r+1):k)
        {
          U[r,s]<-abs((psi_boot[r,s]-psi[r,s])/sqrt(var_psi_boot[r,s]))
        }
      }
      return(max(U)) # returns the LRT value as well as the parameters under full space cell variances and correlations
    }
    ##for bootstrap critical point
    k<-k
    t_data<-t
    n_data<-n
    set.seed(17*b)
    corr_LRT(k,t_data,n_data,values[1:k],values[(k+1):(2*k)],values[(2*k+1):(3*k)])
  }
  critical_boot<-quantile(LRT_value_boot,1-size,names=FALSE) #bootstrap critical value of the LRT
  tau_value<-critical_boot
  CIs_matrix_lower<-array(0,dim=c((k-1),k))
  for(r in 1:(k-1))
  {
    for(s in (r+1):k)
    {
      CIs_matrix_lower[r,s]<-psi[r,s]-(tau_value*sqrt(var_psi[r,s]))
    }
  }

  CIs_matrix_lower<-array(0,dim=c((k-1),k))
  for(r in 1:(k-1))
  {
    for(s in (r+1):k)
    {
      CIs_matrix_lower[r,s]<-psi[r,s]-(tau_value*sqrt(var_psi[r,s]))
    }
  }
  CIs_matrix_upper<-array(0,dim=c((k-1),k))
  for(r in 1:(k-1))
  {
    for(s in (r+1):k)
    {
      CIs_matrix_upper[r,s]<-psi[r,s]+(tau_value*sqrt(var_psi[r,s]))
    }
  }
  for(r in 1:(k-1))
  {
    for(s in (r+1):k)
    {
      cat(paste("lower and upper limits of the SCI for Psi",r,s, "are",CIs_matrix_lower[r,s], "and", CIs_matrix_upper[r,s], "respectively\n"))
    }
  }
}
################################################################################



