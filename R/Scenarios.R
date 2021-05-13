#' Scenario-1:     Large Location shift and COvariance shift using Dirichilet Distribution
#'
#'
#' Simulate data for two classes using scenario-1(Hinton (2021)). Data are sampled from the Dirichilet with different alpha vectors
#'
#' @param n1 number of sample in class 1
#' @param n2 number of samples in class 2
#' @param dms_ number of features
#' @param seed random seed
#'
#' @examples
#' \dontrun{
#' scenario1(n1=40,n2=40,dms_ = 75)
#'}
#'
#' @references
#' Aitchinson, J 1986
#'
#' @return A data frame with simulated data from two classes
#'
#' @export
scenario1 <-
  function( n1 = 30 ,n2 = 30, dms_ = 75,seed){

    a = rep(1,dms_)
    a = a
    expectedScale = dms_
    scale = sum(a)
    scaleRatio = scale/(dms_)
    set.seed(seed)
    s1 = sampleDirichlet(n1 = n1,dims = dms_,sampleName = "S1",a1 = a*3)$Sample
    s2 = sampleDirichlet(n1 = n2,dims = dms_,sampleName = "S2",a1 = a*(log(dms_)/5))$Sample
    df = rbind(s1,s2)

    return(df)

  }


#' Scenario-2:  Sparse Location shift in single component of Dirichilet Distribution
#'
#'
#' Simulate data for two classes using scenario-2(Hinton (2021)). Data are sampled from the Dirichilet and converted to counts where library size is sample from a negative binomial dsitribution
#'
#' @param n1 number of sample in class 1
#' @param n2 number of samples in class 2
#' @param dms_ number of features
#' @param seed random seed
#' @param avgReads avg number reads; simulates libray size using a negative binomial distribution
#' @param size_ dispersion parameter for negative binomial distribution
#'
#' @examples
#' \dontrun{
#' scenario2(n1=40,n2=40,dms_ = 75)
#'}
#'
#' @references
#' Aitchinson, J 1986
#'
#' @return A data frame with simulated data from two classes
#'
#' @export
scenario2 <-
  function( n1 = 30 , n2 = 30,dms_ = 75,seed,avgReads = 1e7,size_ = 1){

    set.seed(seed)
    factScale = 0.01

    libSize = stats::rnbinom(n = n1,size = size_,mu = avgReads)
    s1 = data.frame()
    for(n in 1:n1){
      t25 = 10
      a = rep(stats::runif(n = dms_,min = 1,max = 5))
      a[1] = stats::runif(1,3000,5000)
      #a[2] = runif(1,9000,12500)
      for (i in 2:t25) {
        a[i] = stats::runif(1,500,1500)
      }
      #Random Sparsty
      a[t25:dms_] = sample(a[t25:dms_],replace = F)
      ph = sampleDirichlet(n1 = 2,dims = dms_,sampleName = "S1",a1 = a*factScale)$Sample[1,]
      ph = round(as.numeric(ph[1,-1])*libSize[n])
      ph = data.frame(Status = "S1",t(ph))
      s1 = rbind(s1,ph)
    }


    libSize = stats::rnbinom(n = n2,size = size_,mu = avgReads)
    s2 = data.frame()
    for(n in 1:n2){
      t25 = 10
      a = rep( stats::runif(n = dms_,min = 1,max = 5))
      a[1] = stats::runif(1,12500,17500)
      for (i in 2:t25) {
        a[i] = stats::runif(1,500,1500)
      }
      #Random Sparsty
      a[t25:dms_] = sample(a[t25:dms_],replace = F)
      ph = sampleDirichlet(n1 = 2,dims = dms_,sampleName = "S2",a1 = a*factScale)$Sample[1,]
      ph = round(as.numeric(ph[1,-1])*libSize[n])
      ph = data.frame(Status = "S2",t(ph))
      s2 = rbind(s2,ph)
    }

    df = rbind(s1,s2)

    return(df)
  }



#' Scenario-3:  Large location shift, and small covariance difference from an additive logistics normal distribtion
#'
#'
#' Simulate data for two classes using scenario-2(Hinton (2021)). Data are sampled from the additive logistics normal distribution (see Aitichinson (1986))
#'
#' @param n1 number of sample in class 1
#' @param n2 number of samples in class 2
#' @param dms_ number of features
#' @param seed random seed
#'
#' @examples
#' \dontrun{
#' scenario3(n1=40,n2=40,dms_ = 75)
#'}
#'
#'
#' @references
#' Aitchinson, J 1986
#'
#' @return A data frame with simulated data from two classes
#'
#' @export
scenario3 <-
  function( n1 = 30 , n2 = 30, dms_ = 75,seed){

    set.seed(seed)
    mu_ = rep(0,dms_)
    t25 = round(dms_*.25)
    mu2 = rep(0,dms_)
    mu2[1:t25] = 1/sqrt(dms_)
    sigma_ = diag(dms_)

    diag(sigma_[-1,])<-rep(.2,ncol(sigma_)-1)
    sigma_ = t(sigma_)
    diag(sigma_[-1,])<-rep(.2,ncol(sigma_)-1)
    sigma_ = Matrix::nearPD(sigma_)$mat

    U = matrix(stats::runif(dms_*dms_,0,32/(dms_^2)),nrow = dms_)
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
  }





#' Scenario-4:  Location shift in single component of additive logistics normal distribution
#'
#'
#' Simulate data from two classes using scenario-4 (Hinton (2021)). Data are sampled from the additive logistics normal distribution (see Aitichinson (1986))
#'
#' @param n1 number of sample in class 1
#' @param n2 number of samples in class 2
#' @param dms_ number of features
#' @param seed random seed
#'
#' @examples
#' \dontrun{
#' scenario4(n1=40,n2=40,dms_ = 75)
#'}
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


    U = matrix(stats::runif(dms_*dms_,0,32/dms_^2),nrow = dms_)
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




#' Scenario-5:  No Location shift with large covariance difference using additive logistics normal distribution
#'
#'
#' Simulate data from two classes using scenario-5 (Hinton (2021)). Data are sampled from the additive logistics normal distribution (see Aitichinson (1986))
#'
#' @param n1 number of sample in class 1
#' @param n2 number of samples in class 2
#' @param dms_ number of features
#' @param seed random seed
#'
#' @examples
#' \dontrun{
#' scenario5(n1=40,n2=40,dms_ = 75)
#'}
#'
#' @references
#' Aitchinson, J 1986
#'
#' @return A data frame with simulated data from two classes
#'
#' @export
scenario5 <-
  function( n1 = 30 , n2 = 30, dms_ = 75,seed){

    set.seed(seed)
    mu_ = rep(1,dms_)
    #t25 = round(dms_*.25)
    mu2 = rep(1/sqrt(n1+n2),dms_)
    # mu2[1:t25] = 0
    sigma_ = diag(dms_)

    diag(sigma_[-1,])<-rep(.2,ncol(sigma_)-1)
    sigma_ = t(sigma_)
    diag(sigma_[-1,])<-rep(.2,ncol(sigma_)-1)
    sigma_ = Matrix::nearPD(sigma_)$mat

    U = matrix(stats::runif(dms_*dms_,0,32),nrow = dms_)
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
  }



#' Sample from Additive Logisitcs Normal
#'
#'
#' Simulate data from additive logistics normal distribution (see Aitichinson (1986))
#'
#' @param n1 number of samples
#' @param dims number of features/dimensions
#' @param a1 alpha vector of dirichiliet distribution
#' @param sampleName class names for data
#'
#' @examples
#' \dontrun{
#' sampleDirichlet()
#'}
#'
#' @references
#' Aitchinson, J 1986
#'
#' @return A data frame with simulated data
#'
#' @export
#'
sampleDirichlet <-
  function(n1,dims,a1,sampleName){

    mean1 = sapply(1:dims,function(x) a1[x]/sum(a1))
    variance1 = sapply(1:dims,function(x) {
      a0 = sum(a1)
      ai = a1[x]/a0
      (ai*(1-ai)) / (a0+1)
    } )

    a0 = sum(a1)
    entropy1 = log(prod(gamma(a1))/gamma(sum(a1))) + (a0-length(a1))*digamma(a0) - sum(sapply(1:dims, function(x) (a1[x]-1)*digamma(a1[x])))

    s1 = data.frame(Status = sampleName,compositions::rDirichlet.acomp(n1,a1))
    return(list(Sample = s1,mean = mean1,variance = variance1,entropy = entropy1))

  }



#' Sample from Additive Logisitcs Normal
#'
#'
#' Simulate data from additive logistics normal distribution (see Aitichinson (1986))
#'
#' @param n number of samples
#' @param dims_ number of features/dimensions
#' @param mu mean vector of length 'dims_'
#' @param sigma dims_ x dims_ covaraince matrix
#' @param sigmaScale scale factor for covaraince matrix
#' @param sampleName class names for data
#'
#' @examples
#' \dontrun{
#' sampleAddLogisticNormal(n = 40,dims_ = 75,mu = rep(0,75) ,
#'   sigma = diag(1,75),sigmaScale = 1,sampleName = "S1")
#'}
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



#' Simulate large covaraince shift from empirical data
#'
#'
#' Simulate  covariance shift in data assuming  count data are distributied as an additive logistics normal distribution on the simplex (see Aitichinson (1986))
#'
#' @param raMatrix empirical relative abundance matrix
#' @param n1 number of samples class 1
#' @param n2 number of samples class 2
#' @param maxCov controls the dispersion of \eqn{U} in the formation of the covariance matrix see Hinton(2021)
#' @param perFixedFeatures percent of lowest variance features to shift mean vector
#' @param shiftPercent percent to shift mean vector of the 'perFixedFeatures' number of features
#' @param seed random seed
#'
#' @examples
#' \dontrun{
#' simFromExpData.covarianceShift()
#' }
#'
#' @references
#' Aitchinson, J 1986
#'
#' @return A data frame with simulated data
#'
#' @export
#'
simFromExpData.covarianceShift <-
  function( raMatrix ,n1 = 30 , n2 = 30,maxCov = 32,seed= 08272008,perFixedFeatures = .8,shiftPercent = 1){

    set.seed(seed)

    var = NULL
    y.alr = compositions::alr( raMatrix )
    alr.cov = stats::cov(y.alr)
    var.alr = data.frame(ratio = names(diag(alr.cov)),var = diag(alr.cov)) %>%
      dplyr::arrange(-var)
    y.alr = subset(y.alr,select = var.alr$ratio)
    alr.mean = colMeans(y.alr)
    sigma_ = stats::cov(y.alr)
    nfeats = length(alr.mean)


    # #Adjust cov matrices
    U = matrix(stats::runif(nfeats*nfeats,0,maxCov ),nrow = nfeats)
    U_ = sigma_ + U
    U_ = Matrix::nearPD(U_)$mat
    eg = min(eigen(as.matrix(sigma_))$values)
    sig = min(eigen(as.matrix(U_))$values)
    dd = min(eg,sig)+0.05
    sig1 = Matrix::nearPD( sigma_ + dd*diag(nfeats) )$mat
    sig2 = Matrix::nearPD( sigma_+ U + dd*diag(nfeats) )$mat

    #Sample1
    alr.mean = colMeans(y.alr)
    alr.sim = MASS::mvrnorm(n = n1,mu = alr.mean,Sigma = sig1)
    alr.mean2 = colMeans(alr.sim)
    alr.ra = compositions::alrInv(alr.sim); colnames(alr.ra) = colnames(raMatrix)
    s1 = data.frame(Status = "S1",alr.ra)


    #Sample2
    hold = round(length(alr.mean2)*perFixedFeatures)
    alr.mean2[(hold+1):nfeats] = shiftPercent * alr.mean[(hold+1):nfeats]
    alr.sim = MASS::mvrnorm(n = n2,mu = alr.mean2,Sigma = sig2)
    alr.ra = compositions::alrInv(alr.sim); colnames(alr.ra) = colnames(raMatrix)
    s2 = data.frame(Status = "S2",alr.ra)


    rbind(s1,s2)
  }




#' Simulate location shift in mean vector of Additive logistics normal distribution
#'
#'
#' Simulate  location shift in data assuming  count data are distributied as an additive logistics normal distribution on the simplex (see Aitichinson (1986))
#'
#' @param raMatrix empirical relative abundance matrix
#' @param n1 number of samples class 1
#' @param n2 number of samples class 2
#' @param featureShiftPercent percent to shift mean vector of the 'perFixedFeatures' number of features
#' @param seed random seed
#' @param perFixedFeatures percent of lowest variance features to shift mean vector
#'
#' @examples
#' \dontrun{
#' simFromExpData.largeMeanShft()
#' }
#'
#' @references
#' Aitchinson, J 1986
#'
#' @return A data frame with simulated data
#'
#' @export
#'
simFromExpData.largeMeanShft <-
  function( raMatrix ,n1 = 30 , n2 = 30,perFixedFeatures = 0.95,featureShiftPercent = 1.25,seed = 08272008){
    mu = NULL
    var = NULL
    set.seed(seed)
    y.alr =  compositions::alr( raMatrix )
    alr.cov = stats::cov(y.alr)
    var.alr = data.frame(ratio = names(diag(alr.cov)),var = diag(alr.cov)) %>%
      dplyr::arrange(-var)
    y.alr = subset(y.alr,select = var.alr$ratio)
    alr.mean = colMeans(y.alr)
    alr.cov = stats::cov(y.alr)
    cv = mean(abs( stats::sd(var.alr$var)/alr.mean))


    alr.sim = MASS::mvrnorm(n = n1,mu = alr.mean,Sigma = alr.cov)
    alr.mean2 = colMeans(alr.sim)
    alr.cov = stats::cov(alr.sim)
    alr.ra =  compositions::alrInv(alr.sim); colnames(alr.ra) = colnames(raMatrix)
    s1 = data.frame(Status = "S1",alr.ra)

    ####
    nfeats = length(alr.mean)
    hold = round(length(alr.mean2)*perFixedFeatures)
    alr.mean2[(hold+1):nfeats] = featureShiftPercent * alr.mean[(hold+1):nfeats]
    alr.sim = MASS::mvrnorm(n = n2,mu = alr.mean2,Sigma = alr.cov)
    alr.ra =  compositions::alrInv(alr.sim); colnames(alr.ra) = colnames(raMatrix)
    s2 = data.frame(Status = "S2",alr.ra)
    ###


    rbind(s1,s2)
  }


