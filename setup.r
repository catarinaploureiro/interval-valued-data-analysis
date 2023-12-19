##
#
# Define functions based on the method proposed by Lin et al.(2022).
#
##

library("copula")

#Mean estimation as defined on Theorem 2.1
mean_est <- function(Ymin,Ymax){
    mu<-list()
    for (j in seq_len(ncol(Ymax))){
        mu[[j]]<-0
        for (i in seq_len(nrow(Ymax))){
            mu[[j]]<-mu[[j]]+(Ymin[i,j]+Ymax[i,j])
        }
        mu[[j]]<-mu[[j]]/(2*nrow(Ymax))
    }
    return(mu)
}

#Variance estimation as defined on Theorem 2.1
variance_est <- function(Ymin,Ymax,n,alpha){
    sigma<-list()
    x<-nrow(Ymax)*(qnorm((n-alpha)/(n-2*alpha+1))-qnorm((1-alpha)/(n-2*alpha+1)))
    for (j in seq_len(ncol(Ymax))){
        sigma[[j]]<-0
        for (i in seq_len(nrow(Ymax))){
            sigma[[j]]<-sigma[[j]]+(Ymax[i,j]-Ymin[i,j])
        }
        sigma[[j]]<-(sigma[[j]]/x)^2
    }
    return(sigma)
}

#Clayton copula density function
clayton_cop_df <- function(u,v,rho){
    c1<-(1/rho+1)*(u*v)^(-1/rho-1)
    c2<-(u^(-1/rho)+v^(-1/rho)-1)^(-rho-2)
    return(c1*c2)
}

#Clayton copula joint density function as defined in (2.8)
clayton_joint_df <- function(x,y,mu_x,sig_x,mu_y,sig_y,rho){
    xs<-(x-mu_x)/sqrt(sig_x)
    ys<-(y-mu_y)/sqrt(sig_y)
    c1<-clayton_cop_df(pnorm(xs),pnorm(ys),rho)
    c2<-dnorm(xs)*dnorm(ys)
    return(c1*c2)
}

#Clayton copula joint cumulative distribution function
clayton_joint_cdf <- function(xl,xu,yl,yu,mu_x,sig_x,mu_y,sig_y,rho){
    cl<-mvdc(claytonCopula(param = 1/rho, dim = 2), c("norm", "norm"), list(list(mean = 0, sd =1), list(mean = 0, sd =1)))
    
    xls<-(xl-mu_x)/sqrt(sig_x)
    yls<-(yl-mu_y)/sqrt(sig_y)
    xus<-(xu-mu_x)/sqrt(sig_x)
    yus<-(yu-mu_y)/sqrt(sig_y)

    return(pMvdc(c(xus,yus),cl)-pMvdc(c(xls,yus),cl)-pMvdc(c(xus,yls),cl)+pMvdc(c(xls,yls),cl))
}

#Gumbel copula density function
gumbel_cop_df <- function(u,v,rho){
    c1<-exp(-((-log(u))^rho+(-log(v))^rho)^(1/rho))
    c2<-((log(u)*log(v))^(rho-1))/(u*v)
    c3<-((-log(u))^rho+(-log(v))^rho)^(2/rho-2)
    c4<-(rho-1)*((-log(u))^rho+(-log(v))^rho)^(1/rho-2)
    return(c1*c2*(c3+c4))
}

#Gumbel copula joint density function as defined in (2.9)
gumbel_joint_df <- function(x,y,mu_x,sig_x,mu_y,sig_y,rho){
    xs<-(x-mu_x)/sqrt(sig_x)
    ys<-(y-mu_y)/sqrt(sig_y)
    c1<-gumbel_cop_df(pnorm(xs),pnorm(ys),rho)
    c2<-dnorm(xs)*dnorm(ys)
    return(c1*c2)
}

#Gumbel copula joint cumulative distribution function
gumbel_joint_cdf <- function(xl,xu,yl,yu,mu_x,sig_x,mu_y,sig_y,rho){
    cl<-mvdc(gumbelCopula(param = rho, dim = 2), c("norm", "norm"), list(list(mean = 0, sd =1), list(mean = 0, sd =1)))
    
    xls<-(xl-mu_x)/sqrt(sig_x)
    yls<-(yl-mu_y)/sqrt(sig_y)
    xus<-(xu-mu_x)/sqrt(sig_x)
    yus<-(yu-mu_y)/sqrt(sig_y)

    return(pMvdc(c(xus,yus),cl)-pMvdc(c(xls,yus),cl)-pMvdc(c(xus,yls),cl)+pMvdc(c(xls,yls),cl))
}

#Gaussian copula joint density function as defined in (2.7)
gaussian_joint_df<-function(x,y,mu_x,sig_x,mu_y,sig_y,rho){
    xs<-(x-mu_x)^2/sig_x
    ys<-(y-mu_y)^2/sig_y
    frac<-(2*pi*sqrt(sig_x)*sqrt(sig_y)*sqrt(1-rho^2))^(-1)
    e<-exp(-(2*(1-rho^2))^(-1)*(xs+ys-(2*rho*(x-mu_x)*(y-mu_y))/(sqrt(sig_x)*sqrt(sig_y))))
    return(frac*e)
}

#Gaussian copula joint cumulative distribution function
gaussian_joint_cdf <- function(xl,xu,yl,yu,mu_x,sig_x,mu_y,sig_y,rho){
    cl<-mvdc(normalCopula(param = rho, dim = 2), c("norm", "norm"), list(list(mean = 0, sd =1), list(mean = 0, sd =1)))
    
    xls<-(xl-mu_x)/sqrt(sig_x)
    yls<-(yl-mu_y)/sqrt(sig_y)
    xus<-(xu-mu_x)/sqrt(sig_x)
    yus<-(yu-mu_y)/sqrt(sig_y)

    return(pMvdc(c(xus,yus),cl)-pMvdc(c(xls,yus),cl)-pMvdc(c(xus,yls),cl)+pMvdc(c(xls,yls),cl))
}

#Relationship between Kendall's τ and the Clayton copula parameter ρ
kendall_tau_clayton <- function(rho){
    return(1/(1+2*rho))
}

#Relationship between Kendall’s τ and the Gumbel copula parameter ρ
kendall_tau_gumbel<-function(rho){
    return((rho-1)/rho)
}

#Bivariate g function used in the likelihood function as defined in Theorem 2.2 (2.5)
g_function <- function(xl,xu,yl,yu,n,f,cdf){
    
    I<-cdf(xl,xu,yl,yu)
    Ix<-function(b){sapply(b,function(b){integrate(function(x){f(x,b)},xl,xu)$value})}
    Iy<-function(a){sapply(a,function(a){integrate(function(y){f(a,y)},yl,yu)$value})}
    
    a1<-f(xu,yu)*f(xl,yl)+f(xu,yl)*f(xl,yu)
    a2<-f(xu,yu)*Ix(yl)*Iy(xl)+f(xu,yl)*Ix(yu)*Iy(xl)+f(xl,yu)*Ix(yl)*Iy(xu)+f(xl,yl)*Ix(yu)*Iy(xu)
    a3<-Ix(yu)*Ix(yl)*Iy(xu)*Iy(xl)

    g1<-n*(n-1)*I^(n-2)*a1
    g2<-n*(n-1)*(n-2)*I^(n-3)*a2
    g3<-n*(n-1)*(n-2)*(n-3)*I^(n-4)*a3

    return(g1+g2+g3)
}

#Log of the bivariate likelihood function defined in Theorem 2.2 (2.4)
likelihood_log<-function(g,Xmin,Xmax,Ymin,Ymax){
    l<-list()
    for (i in seq_along(Xmin)){
        if (Xmin[[i]]<Xmax[[i]] && Ymin[[i]]<Ymax[[i]]) {
            g_fun<-g(Xmin[[i]],Xmax[[i]],Ymin[[i]],Ymax[[i]])
            if (g_fun!=0){
                l<-append(l,g_fun)
            }
        }
    }
    return(sum(log(unlist(l))))
}

#Estimate the Gaussian copula parameter ρ using MLE, correlation, and covariance
gaussian_corr<-function(Mmin,Mmax,n,alpha){
    mu<-mean_est(Mmin,Mmax)
    sig<-variance_est(Mmin,Mmax,n,alpha)
    corr<-diag(ncol(Mmax))

    for (i in seq_len(ncol(Mmax))){
        mu_x<-mu[[i]]
        sig_x<-sig[[i]]

        for (j in (1:i)[-i]){
            mu_y<-mu[[j]]
            sig_y<-sig[[j]]

            likelihood_rho<-function(rho){
                f<-function(x,y){gaussian_joint_df(x,y,mu_x,sig_x,mu_y,sig_y,rho)}
                cdf<-function(x_l,x_u,y_l,y_u){gaussian_joint_cdf(x_l,x_u,y_l,y_u,mu_x,sig_x,mu_y,sig_y,rho)}
                g<-function(x_l,x_u,y_l,y_u){g_function(x_l,x_u,y_l,y_u,n,f,cdf)}
                return(likelihood_log(g,Mmin[,i],Mmax[,i],Mmin[,j],Mmax[,j]))
            }

            corr[i,j]<-optimize(likelihood_rho,c(-0.99999,0.99999),maximum = TRUE)$maximum
            corr[j,i]<-corr[i,j]
        }
    }
    cov_est<-diag(sqrt(as.numeric(sig)))%*%corr%*%diag(sqrt(as.numeric(sig)))
    result<-list("corr_est"=corr, "cov_est"=cov_est)
    return(result)
}

#Estimate the copula parameter ρ using MLE, Kendall's correlation τ, and covariance 
kendall_corr <- function(Mmin, Mmax, n, alpha, copula) {

    if (!(copula %in% c("clayton", "gumbel"))){
        stop("Invalid copula type. Supported types: 'clayton', 'gumbel'")
    }

    mu <- mean_est(Mmin, Mmax)
    sig <- variance_est(Mmin, Mmax, n, alpha)
    tau <- diag(ncol(Mmax))

    for (i in seq_len(ncol(Mmax))) {
        mu_x <- mu[[i]]
        sig_x <- sig[[i]]

        for (j in (1:i)[-i]) {
            mu_y <- mu[[j]]
            sig_y <- sig[[j]]

            likelihood_rho <- function(rho) {
                if (copula == "clayton") {
                    f <- function(x, y){clayton_joint_df(x, y, mu_x, sig_x, mu_y, sig_y, rho)}
                    cdf <- function(x_l, x_u, y_l, y_u){clayton_joint_cdf(x_l, x_u, y_l, y_u, mu_x, sig_x, mu_y, sig_y, rho)}
                } else if (copula == "gumbel") {
                    f <- function(x, y){gumbel_joint_df(x, y, mu_x, sig_x, mu_y, sig_y, rho)}
                    cdf <- function(x_l, x_u, y_l, y_u){gumbel_joint_cdf(x_l, x_u, y_l, y_u, mu_x, sig_x, mu_y, sig_y, rho)}
                }

                g <- function(x_l, x_u, y_l, y_u){g_function(x_l, x_u, y_l, y_u, n, f, cdf)}
                return(likelihood_log(g, Mmin[, i], Mmax[, i], Mmin[, j], Mmax[, j]))
            }

            if (copula == "clayton") {
                tau[i, j] <- kendall_tau_clayton(optimize(likelihood_rho, c(0.00001, 20), maximum = TRUE)$maximum)
            } else if (copula == "gumbel") {
                tau[i, j] <- kendall_tau_gumbel(optimize(likelihood_rho, c(1, 20), maximum = TRUE)$maximum)
            }
            tau[j, i] <- tau[i, j]
        }
    }
    cov_est <- diag(sqrt(as.numeric(sig)))%*%tau%*%diag(sqrt(as.numeric(sig)))
    result <- list("tau" = tau, "cov_est" = cov_est)
    return(result)
}

#Minimum and Maximum principal component scores as in (4.2)
pca_scores_symbolic<-function(Mmin,Mmax,Vpca){
    S_L<-list()
    S_U<-list()
    for (i in seq_len(nrow(Mmin))){
        sl<-0
        su<-0
        for (j in seq_len(ncol(Mmin))){
            if(Vpca[j]<0){
                sl<-sl+Mmax[i,j]*Vpca[j]
                su<-su+Mmin[i,j]*Vpca[j]
            }else{
                sl<-sl+Mmin[i,j]*Vpca[j]
                su<-su+Mmax[i,j]*Vpca[j]
            }
        }
        S_L[i]<-sl
        S_U[i]<-su
    }
    result<-list("MinimumPC"=S_L, "MaximumPC"=S_U)
    return(result)
}