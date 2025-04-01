################################################################################
##The LRT T_1 and the asymptotic LRT T_1A are developed for testing the structure of the covariance matrix
##function to evaluate the observed, critical and P-values of the LRT test statistic T_1 and the asymptotic LRT T_1A

#' Function to evaluate the observed, critical and P-values of likelihood ratio and asymptotic likelihood ratio tests for testing covariance structure
#'
#' @param data A data frame of dimension t by n consists of n number of t-dimensional cell data
#' @param size The level of significance, say 0.05
#' @param t The number of time (space) points, say t
#' @param n Sample size, say n
#'
#' @return Numeric values consist of observed, critical and P-values of the likelihood ratio and asymptotic likelihood ratio tests
#' @export
#'
#' @examples
#' data<-read.csv("https://raw.githubusercontent.com/RajuProbStat/OnewayCorrToep/refs/heads/main/DBP_data.csv")
#' T_1_T_1_A(data[1:6],0.05,3,6)
#' @examples
#' data<-read.csv("https://raw.githubusercontent.com/RajuProbStat/OnewayCorrToep/refs/heads/main/DBP_data.csv")
#' T_1_T_1_A(data[7:12],0.05,3,6)
#' @examples
#' data<-read.csv("https://raw.githubusercontent.com/RajuProbStat/OnewayCorrToep/refs/heads/main/DBP_data.csv")
#' T_1_T_1_A(data[13:20],0.05,3,8)
T_1_T_1_A<-function(data,size,t,n){
  library(foreach)
  library(doParallel)
  library(MASS)
  library(matlib)
  library(stats)
  data<-data #the data set
  size<-size #the level of significance
  t<-t #the number of time (space) points
  n<-n #sample size
  x<-data.matrix(data)
  ##speical Toeplitz matrix
  sig<-function(t,variance,rho){
    sig_m<-array(NA,dim=c(t,t))
    for(i in 1:t)
    {
      for(j in 1:t)
      {
        if(i==j)
        {
          sig_m[i,j]<-variance
        }
        else
        {
          if(j==(i-1)|j==(i+1))
          {
            sig_m[i,j]<-variance*rho
          }
          else
          {
            sig_m[i,j]<-0
          }
        }
      }
    }
    return(sig_m)
  }
  ##the orthogonal matrix P
  P<-array(NA,dim=c(t,t))
  for(i in 1:t)
  {
    for(j in 1:t)
    {
      P[i,j]<-sqrt(2/(t+1))*sin((pi*i*j)/(t+1))
    }
  }
  ##the inverse of the diagonal matrix D
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
  ##the diagonal matrix having the diagonal elements as the derivatives of the diagonal elements of the diagonal matrix D
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
  x<-as.matrix(x)
  mean_x<-apply(x,1,mean)
  ## initialization of the value of rho
  bvn_data_length<-n*floor(t/2)
  bvn_data<-array(NA,dim=c(2,bvn_data_length))
  for(i in 1:n)
  {
    for(k in ((i-1)*floor(t/2)+1):(i*floor(t/2)))
    {
      bvn_data[,k]<-x[,i][(2*(k-((i-1)*floor(t/2)))-1):(2*(k-((i-1)*floor(t/2))))]
    }
  }
  cor_value<-cor(bvn_data[1,],bvn_data[2,])
  repeat
  {
    if(abs(cor_value)<0.5)
    {
      break
    }
    cor_value<-0.8*cor_value
  }
  ##MLEs and the log-likelihood function value under the full parameter space
  mu_0<-(1/t)*(rep(1,t)%*%mean_x)
  sigma_0<-array(0,dim=c(t,t))
  for(j in 1:n)
  {
    sigma_0<-sigma_0+(1/n)*(x[,j]-rep(mu_0,t))%*%t((x[,j]-rep(mu_0,t)))
  }
  repeat
  {
    mu_1<-(t(mean_x)%*%solve(sigma_0)%*%rep(1,t))/(t(rep(1,t))%*%solve(sigma_0)%*%rep(1,t))
    sigma_1<-array(0,dim=c(t,t))
    for(j in 1:n)
    {
      sigma_1<-sigma_1+(1/n)*(x[,j]-rep(mu_1,t))%*%t((x[,j]-rep(mu_1,t)))
    }
    ##
    sum1<-0
    for(j in 1:n)
    {
      sum1<-sum1+t(x[,j]-rep(mu_1,t))%*%solve(sigma_1)%*%(x[,j]-rep(mu_1,t))
    }
    loglike_value_full<- ((n*t)/2)*log(2*pi)+(n/2)*log(det(sigma_1))+(1/2)*sum1
    ##
    if(abs(mu_1-mu_0)<0.000001)
    {
      break
    }
    mu_0<-mu_1
    sigma_0<-sigma_1
  }
  ##MLEs and the log-likelihood function value under the null parameter space
  y<-P%*%x
  w<-P%*%rep(1,t)
  det_sigma_full<-det(sigma_0)
  fr<-function(z){
    mu_fun<-z[1]
    rho_fun<-z[2]
    var_fun<-z[3]
    sum1<-0
    for(j in 1:n)
    {
      sum1<-sum1+t(y[,j]-mu_fun*w)%*%D_inv_fun(rho_fun)%*%(y[,j]-mu_fun*w)
    }
    return(((n*t)/2)*log(2*pi*var_fun)-(n/2)*log(det(D_inv_fun(rho_fun)))+(1/(2*var_fun))*sum1)
  }
  grr<-function(z){
    mu_fun<-z[1]
    rho_fun<-z[2]
    var_fun<-z[3]
    sum1<-0
    for(j in 1:n)
    {
      sum1<-sum1+t(y[,j]-mu_fun*w)%*%D_inv_fun(rho_fun)%*%w
    }
    sum2<-0
    for(j in 1:n)
    {
      sum2<-sum2+t(y[,j]-mu_fun*w)%*%D_inv_fun(rho_fun)%*%(y[,j]-mu_fun*w)
    }
    sum3<-0
    for(i in 1:t)
    {
      sum3<-sum3+(2*cos((i*pi)/(t+1)))/(1+2*rho_fun*cos((i*pi)/(t+1)))
    }
    sum4<-0
    for(j in 1:n)
    {
      sum4<-sum4+t(y[,j]-mu_fun*w)%*%D_diff_inv_fun(rho_fun)%*%(y[,j]-mu_fun*w)
    }
    c(-(1/var_fun)*sum1,(n/2)*sum3+(1/(2*var_fun))*sum4,(n*t)/(2*var_fun)-(1/(2*var_fun^2))*sum2)
  }
  par_null<-constrOptim(c(mean(x),cor_value,mean(apply(x,1,var))),fr,grr,ui=rbind(c(0,1,0),c(0,-1,0),c(0,0,1)),ci=c(-0.5,-0.5,0), method="BFGS")$par
  loglike_value_null<-fr(c(par_null[1],par_null[2],par_null[3]))
  LRT_value<-exp(loglike_value_full-loglike_value_null) #the value of the likelihood ratio

  ##for bootstrap step
  cl<-makeCluster(4) #number of cores used in the parallelization
  registerDoParallel(cl)
  B<-1000 #the number of bootstrap replications
  start_time <- Sys.time()
  LRT_value_boot<-foreach(b=1:B, .combine=rbind)%dopar%{
    library(MASS)
    library(foreach)
    library(doParallel)
    library(matlib)
    library(stats)
    set.seed(17*b)
    sigma_boot<-sig(t,par_null[3],par_null[2])
    x_boot<-array(NA,dim=c(t,n))
    for(j in 1:n)
    {
      x_boot[,j]<-mvrnorm(1,rep(0,t),sigma_boot)
    }
    mean_x_boot<-apply(x_boot,1,mean)
    ##
    bvn_data_length_boot<-n*floor(t/2)
    bvn_data_boot<-array(NA,dim=c(2,bvn_data_length_boot))
    for(i in 1:n)
    {
      for(k in ((i-1)*floor(t/2)+1):(i*floor(t/2)))
      {
        bvn_data_boot[,k]<-x_boot[,i][(2*(k-((i-1)*floor(t/2)))-1):(2*(k-((i-1)*floor(t/2))))]
      }
    }
    cor_value_boot<-cor(bvn_data_boot[1,],bvn_data_boot[2,])
    repeat
    {
      if(abs(cor_value_boot)<0.5)
      {
        break
      }
      cor_value_boot<-0.8*cor_value_boot
    }
    ##
    mu_0_boot<-(1/t)*(rep(1,t)%*%mean_x_boot)
    sigma_0_boot<-array(0,dim=c(t,t))
    for(j in 1:n)
    {
      sigma_0_boot<-sigma_0_boot+(1/n)*(x_boot[,j]-rep(mu_0_boot,t))%*%t((x_boot[,j]-rep(mu_0_boot,t)))
    }
    repeat
    {
      mu_1_boot<-(t(mean_x_boot)%*%solve(sigma_0_boot)%*%rep(1,t))/(t(rep(1,t))%*%solve(sigma_0_boot)%*%rep(1,t))
      sigma_1_boot<-array(0,dim=c(t,t))
      for(j in 1:n)
      {
        sigma_1_boot<-sigma_1_boot+(1/n)*(x_boot[,j]-rep(mu_1_boot,t))%*%t((x_boot[,j]-rep(mu_1_boot,t)))
      }
      ##
      sum1<-0
      for(j in 1:n)
      {
        sum1<-sum1+t(x_boot[,j]-rep(mu_1_boot,t))%*%solve(sigma_1_boot)%*%(x_boot[,j]-rep(mu_1_boot,t))
      }
      loglike_value_full_boot<- ((n*t)/2)*log(2*pi)+(n/2)*log(det(sigma_1_boot))+(1/2)*sum1
      ##
      if(abs(mu_1_boot-mu_0_boot)<0.000001)
      {
        break
      }
      mu_0_boot<-mu_1_boot
      sigma_0_boot<-sigma_1_boot
    }
    y_boot<-P%*%x_boot
    fr_boot<-function(z){
      mu_fun_boot<-z[1]
      rho_fun_boot<-z[2]
      var_fun_boot<-z[3]
      sum1<-0
      for(j in 1:n)
      {
        sum1<-sum1+t(y_boot[,j]-mu_fun_boot*w)%*%D_inv_fun(rho_fun_boot)%*%(y_boot[,j]-mu_fun_boot*w)
      }
      return(((n*t)/2)*log(2*pi*var_fun_boot)-(n/2)*log(det(D_inv_fun(rho_fun_boot)))+(1/(2*var_fun_boot))*sum1)
    }
    grr_boot<-function(z){
      mu_fun_boot<-z[1]
      rho_fun_boot<-z[2]
      var_fun_boot<-z[3]
      sum1<-0
      for(j in 1:n)
      {
        sum1<-sum1+t(y_boot[,j]-mu_fun_boot*w)%*%D_inv_fun(rho_fun_boot)%*%w
      }
      sum2<-0
      for(j in 1:n)
      {
        sum2<-sum2+t(y_boot[,j]-mu_fun_boot*w)%*%D_inv_fun(rho_fun_boot)%*%(y_boot[,j]-mu_fun_boot*w)
      }
      sum3<-0
      for(i in 1:t)
      {
        sum3<-sum3+(2*cos((i*pi)/(t+1)))/(1+2*rho_fun_boot*cos((i*pi)/(t+1)))
      }
      sum4<-0
      for(j in 1:n)
      {
        sum4<-sum4+t(y_boot[,j]-mu_fun_boot*w)%*%D_diff_inv_fun(rho_fun_boot)%*%(y_boot[,j]-mu_fun_boot*w)
      }
      c(-(1/var_fun_boot)*sum1,(n/2)*sum3+(1/(2*var_fun_boot))*sum4,(n*t)/(2*var_fun_boot)-(1/(2*var_fun_boot^2))*sum2)
    }
    par_null_boot<-constrOptim(c(mean(x_boot),cor_value_boot,mean(apply(x_boot,1,var))),fr_boot,grr_boot,ui=rbind(c(0,1,0),c(0,-1,0),c(0,0,1)),ci=c(-0.5,-0.5,0), method="BFGS")$par
    loglike_value_null_boot<-fr_boot(c(par_null_boot[1],par_null_boot[2],par_null_boot[3]))
    exp(loglike_value_full_boot-loglike_value_null_boot) #the value of likelihood ratio based on the bootstrap samples
  }
  critical_boot<-quantile(LRT_value_boot,size,names=F) #the bootstrap critical value of the LRT T_1
  cat("The observed value, critical value and the P-value of the LRT statistic are\n",LRT_value,critical_boot,mean(as.numeric(LRT_value)>LRT_value_boot))
  if(LRT_value<critical_boot)
  {
    cat("\n The LRT test rejects the null hypothesis")
  }
  else
  {
    cat("\n THE LRT test do not reject the null hypothesis")
  }
  critical_boot_asym<-qchisq(p =1-size ,df =((t*(t+1))/2)-2 ,ncp =0 ,lower.tail =TRUE) #the critical value of the asymptotic LRT T_1A
  cat("\n The observed value, critical value and the P-value of the aymptotic LRT statistic are\n",-2*log(LRT_value),critical_boot_asym,1-pchisq(q =-2*log(LRT_value) ,df =((t*(t+1))/2)-2 ,ncp =0 ,lower.tail = TRUE))
  if(-2*log(LRT_value)>critical_boot_asym)
  {
    cat("\n The asymptotic LRT rejects the null hypothesis")
  }
  else
  {
    cat("\n The asymptotic LRT do not reject the null hypothesis")
  }
}
################################################################################






