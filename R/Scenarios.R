#' Scenario-4:  Location shift in single component of additive logistics normal distribution
#'
#'
#' Simulate data from two classes using scenario-4. Data are sampled from the additive logistics normal distribution (see Aitichinson (1986))
#'
#' @param n1 number of sample in class 1
#' @param n2 number of samples in class 2
#' @param dms_ number of features
#' @param seed random seed
#'
#' @examples
#' scenario4(n1=40,n2=40,dms_ = 75)
#'
#' @references
#' Aitchinson, J 1986
#'
#' @return A data frame with simulated data from two classes
#'
#' @export
scenario4 <-
  function( n1 = 30 , n2 = 30, dms_ = 75,seed = 08272008){


    set.seed(seed)
    set.seed(seed)
    mu_ = rep(0,dms_)
    #t25 = round(dms_*.25)
    mu2 = rep(0,dms_)
    mu2[1] = log(dms_)/3
    sigma_ = diag(dms_)


    U = matrix(runif(dms_*dms_,0,32/dms_^2),nrow = dms_)
    U_ = sigma_ + U
    U_ = Matrix::nearPD(U_)$mat
    eg = min(eigen(as.matrix(sigma_))$values)
    sig = min(eigen(as.matrix(U_))$values)
    dd = min(eg,sig)+0.05
    sig1 = Matrix::nearPD( sigma_ + dd*diag(dms_) )$mat
    sig2 = Matrix::nearPD( sigma_+ U + dd*diag(dms_) )$mat

    s1 = sampleAddLogisticNormal(n = n1,dims_ = dms_,mu = mu_ ,sigma = sig1,sigmaScale = 1,sampleName = "S1")
    s2 = sampleAddLogisticNormal(n = n2,dims_ = dms_,mu = mu2,sigma = sig2,sigmaScale = 1,sampleName = "S2")
    df = rbind(s1,s2)
    return(df)
  }


#' Sample from Additive Logisitcs Normal
#'
#'
#' Simulate data from additive logistics normal distribution (see Aitichinson (1986))
#'
#' @param n number of samples
#' @param dims_ number of features/dimensions
#' @param sigma dims_ x dims_ covaraince matrix
#' @param sigmaScale scale factor for covaraince matrix
#' @param sampleName class names for data
#'
#' @examples
#' sampleAddLogisticNormal(n = 40,dims_ = 75,mu = rep(0,75) ,sigma = diag(1,75),sigmaScale = 1,sampleName = "S1")
#'
#' @references
#' Aitchinson, J 1986
#'
#' @return A data frame with simulated data
#'
#' @export
#'
sampleAddLogisticNormal <-
  function(n, dims_,mu,sigma,sigmaScale=1,sampleName){
    sigma = sigma*sigmaScale
    y1 = MASS::mvrnorm(n , mu = mu, Sigma = sigma, tol = 1e-6, empirical = FALSE)
    y1 = compositions::alrInv(y1)
    s1  = data.frame(Status = sampleName,y1)
    return(s1)
  }
