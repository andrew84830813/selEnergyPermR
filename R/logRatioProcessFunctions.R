#' Calculate all pairwise log-ratios
#'
#'
#' Compute all pairwise log-ratios froma  data.frame where rows=samples and columns=features/parts
#'
#' @param df a data.frame where the first column contains class labels and the remaining p columns are features
#'
#' @examples
#' calcLogRatio()
#'
#' @references
#'
#' @return A data frame with all choose(p,2) log ratios
#'
#' @export
calcLogRatio <-
  function(df){
    df1 = stageData(featureMatrix = df[,-1],labels = df[,1],
                    permuteLabel = F,permuteFeatures = F)
    df2 = df1$allData
    fn = length(colnames(df2[,-1]))
    mat = matrix(rep(0,fn*fn),nrow = fn)
    colnames(mat) = colnames(df2[,-1])
    rownames(mat) =  colnames(df2[,-1])
    mat[lower.tri(mat,diag = T)]=-1000000000000000000999
    mat = subset(reshape2::melt(mat),value!=-1000000000000000000999)
    mat = tidyr::unite(data = mat,col = Ratio,sep = "___",c("Var1","Var2"))
    keyRats = dplyr::separate(data.frame(mat),1,into = c("Num","Denom"),sep = "___",remove = F)
    el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
    df1 = data.frame(Status = df[,1], getLogRatios(z = df2,ratioList  = el_))

    return(df1)
  }



#' Impute zeroes
#'
#'
#' Impute zeroes using a lightly modified multiplicative replacement strategy (Martin-Fernandez (2003))
#'
#' @param df.cdata2 count or relative abundance matrix where rows are samples and columns are features
#' @param impFactor delta imputed value
#'
#' @examples
#' fastImputeZeroes()
#'
#' @references
#' Martin-Fernandez, J. A., Barcelo-Vidal, C., & Pawlowsky-Glahn, V. (2003). Dealing with Zeros and Missing Values in Compositional Data Sets Using Nonparametric Imputation. Mathematical Geology, 26.
#'
#' @return A matrix with multiplicatively imputed zeroes
#'
#' @export
fastImputeZeroes <-
  function(df.cdata2, impFactor = 1e-11){

    require(doParallel)

    cn = colnames(df.cdata2)
    rn = rownames(df.cdata2)
    df.cdata2 = as.matrix(compositions::clo(df.cdata2))

    if(is.null(impFactor)){
      impFactor =min(df.cdata2[df.cdata2>0])/10
    }

    df_imputed = foreach::foreach(i  = 1:nrow(df.cdata2),.combine = rbind) %dopar% {
      sampleData = df.cdata2[i,]
      nz = sum(sampleData==0)
      sampleData[sampleData==0] = impFactor
      sampleData[sampleData!=0] =  sampleData[sampleData!=0]*(1-impFactor*nz)
    }

    colnames(df_imputed) = cn
    rownames(df_imputed) = rn
    return(df_imputed)
  }


#' Get non-sparse feature
#'
#'
#' Computes names of features present above a thresold
#'
#' @param compmat a matrix/data.frame with counts (sample x features)
#' @param mxSparsePercent thresold for percent of samples a feature must be present with at least 1 count
#'
#' @examples
#' getFeats_bySparsity()
#'
#' @references
#'
#'
#' @return A vector of colnames of non sparse features to keep
#'
#' @export
getFeats_bySparsity <-
  function(compmat,mxSparsePercent = .5,ncountRequired = 1){
    zeroFeats = compmat==0
    zeroFeats = colSums(zeroFeats)
    sparseThreshold = round(mxSparsePercent*nrow(compmat))
    keep = zeroFeats[zeroFeats<=sparseThreshold]
    return(names(keep))
  }

#' Calculate select log-ratios
#'
#'
#' Calculate specific log-ratios defined in a matrix where rows are ratio and columns describe the numerator, denominator and log-ratio names
#'
#' @param z a matrix/data.frame with imputed relative abundace data (sample x features)
#' @param ratioList a data.frame with columns eqaul to the numerator of log-ratio, denominator of log-ratio, and the log-ratio names
#'
#' @examples
#' getLogRatios()
#'
#' @references
#'
#'
#' @return A data frame with simulated data from two classes
#'
#' @export
getLogRatios <-
  function(z,ratioList){
    #Check To make sure status variable is removed
    if(sum(colnames(z)=="Status")>0){
      c = which(colnames(z)=="Status")
      z = z[,-c]
    }

    #Create feature to location library
    features = data.frame(loc = 1:ncol(z),f = colnames(z))

    #Map Num and denom feature locations
    num = data.frame(f = ratioList[,1])
  }


#' Wrapper function for removing spasre features from count data prior to log ratio analysis
#'
#'
#' Removes sparse features based on a provided threshold,computes \eqn{\delta} impute factor for zero impuation and returns the minority and majority class for binary group data
#'
#' @param tbl a n-samples x p-features data.frame with count or relative abundance data
#' @param minPrevalence class labels for data
#' @param permuteLabel should lass labels be permuted
#'
#' @examples
#' processCompData()
#'
#' @references
#'
#'
#' @return A list containing:\tabular{ll}{
#'    \code{processedData} \tab a data.frame where the first column contains class labels and the remaining p columns are features \cr
#'    \tab \cr
#'    \code{impFactor} \tab A vector of class lables \cr
#'    \tab \cr
#'    \code{minClss} \tab A vector of class lables \cr
#'    \tab \cr
#'    \code{majClass} \tab A vector of class lables \cr
#' }
#' @export
processCompData <-
  function(tbl,minPrevalence = 0.9,permuteLabel = F){
    print(table(tbl[,1]))
    minorityClass = names(which.min(table(tbl[,1])))
    majorityClass = as.character(unique(tbl[tbl[,1]!=minorityClass,1]))
    keep = getFeats_bySparsity(tbl[,-1],mxSparsePercent = minPrevalence)
    fcol = colnames(tbl)[1]
    tbl = subset(tbl,select = c(fcol,keep))
    tbl = tbl[rowSums(tbl[,-1])>0,]
    factor = 1
    ph = compositions::clo(tbl[,-1])
    impFact = min(ph[ph>0]) / factor
    tbl = data.frame(Status = factor(tbl[,1]), tbl[,-1]   )
    #permuteLabels
    permuteLabel = F
    if(permuteLabel==T){
      tbl$Status = sample(tbl$Status)
    }
    return(list(processedData = tbl,impFactor = impFact,minClss = minorityClass,majClass = majorityClass))
  }



#' Wrapper function for processing count/relative data for log-ratio analysis
#'
#'
#' Close data and Impute zeroes. Allow permuting labels or columns for permutation distribution analysis
#'
#' @param featureMatrix a n-samples x p-features data.frame with count or relative abundance data
#' @param labels class labels for data
#' @param permuteLabel should lass labels be permuted
#' @param permuteFeatures should features be permuted
#'
#' @examples
#' stageData()
#'
#' @references
#'
#'
#' @return A list containing:\tabular{ll}{
#'    \code{allData} \tab a data.frame where the first column contains class labels and the remaining p columns are features \cr
#'    \tab \cr
#'    \code{y_train} \tab A vector of class lables \cr
#' }
#' @export
stageData <-
  function(featureMatrix,labels,permuteLabel=F,permuteFeatures=F){


    #Data Preprocess - Closure/Normalization ###########################
    xtrain = compositions::clo(featureMatrix)

    ytrain = labels
    if(permuteLabel==T){
      ytrain = sample(ytrain)
    }else{
      ytrain = ytrain
    }
    ############################################################################*
    ###  Permute Labels to CHeck for statistically significant pattern ###
    ############################################################################*
    classes = as.character(unique(ytrain))
    ytrain = factor(as.vector(as.matrix(ytrain)),levels = c(classes[1],classes[2]))
    if(permuteFeatures==T){
      c.names = colnames(xtrain)
      ncol=length(c.names)
      posCases = xtrain[ytrain==classes[1],]
      negCases = xtrain[ytrain==classes[2],]
      posCases = posCases[,sample(1:ncol,ncol,replace = F)]
      negCases = negCases[,sample(1:ncol,ncol,replace = F)]
      colnames(posCases) = c.names
      colnames(negCases)=c.names
      xtrain=rbind(posCases,negCases)
    }
    ytrain1 = as.factor(if_else(ytrain==classes[1],1,0))



    ############################################################################*
    ####    Data Preprocessing - Zero Imputation ####
    ############################################################################*
    xtrain1_imputed = fastImputeZeroes(compositions::clo(xtrain))
    allTrain.method = data.frame(Status = ytrain1,xtrain1_imputed)
    ############# ZERO IMPUTATION (END) ###########################

    return(list(allData = allTrain.method,
                y_train = ytrain)
    )

  }



