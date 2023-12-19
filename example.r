##
#
# Test the correlation estimation on a simulated dataset
#
##

source("./setup.r")

library("MASS")

#Set mean, variance, correlation, and covariance
test_mu<-c(3,1)
test_sig<-c(1,0.25)
test_corr<-matrix(c(1,0.7,0.7,1),nrow=2)
test_cov<-diag(sqrt(test_sig))%*%test_corr%*%diag(sqrt(test_sig))

#Create dataset from a multivariate normal distribution
test_data<-mvrnorm(n = 25920, mu=test_mu, test_cov)
test_data<-cbind(test_data,rep(1:1080, each = 24)) #add a group variable to aggregate

#Create interval-valued dataset, aggregating by maximum and minimum
test_max<-aggregate(test_data[,1:2], by=list(test_data[,3]), FUN="max")
test_min<-aggregate(test_data[,1:2], by=list(test_data[,3]), FUN="min")

#Estimate the mean vector, covariance matrix, and correlation matrix
n<-24 #number of observations per group
alpha<-0.375

test_mu<-mean_est(test_min[,-1],test_max[,-1])
test_sig<-variance_est(test_min[,-1],test_max[,-1],n,alpha)

test_gau<-gaussian_corr(test_min[,-1],test_max[,-1],n,alpha)
test_ken_clayton<-kendall_corr(test_min[,-1],sim_max[,-1],n,alpha,"clayton")
test_ken_gumbel<-kendall_corr(test_min[,-1],sim_max[,-1],n,alpha,"gumbel")