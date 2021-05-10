#' Selection Energy Permutation
#'
#'
#' A sparse non parametric association test for compositional/omics data  (Hinton (2021))
#'
#' @importFrom magrittr %>%
#' @param inputData xxx
#' @param optimizationMetric xxx
#' @param dcv_useInfoGain xxx
#' @param dcv_nfold xxx
#' @param dcv_numRepeats xxx
#' @param nreps_energy xxx
#' @param eps xxx
#' @param targetFeats xxx
#' @param patience xxx
#' @param alpha_ xxx
#'
#' @examples
#' selectionEnergy.scaled()
#'
#' @references
#' Hinton, A.L., Mucha, P.J., (2021). Simultaneous variable selection and group association testing in sparse high dimensional omics data. XXXX.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{lrs} \tab a p-row log-ratio scoring matrix aggregated across nfolds and numRepeats  \cr
#'    \tab \cr
#'    \code{rawDCV} \tab The raw log-ratio scoring matrix with scores from each metric (Hinton(2021)) unaggreagted \cr
#' }
#' use_pipe(export = TRUE)
#'
#' @export
selectionEnergy.scaled = function(inputData , optimizationMetric = NULL,
                                  dcv_useInfoGain = T,dcv_nfold = 5,dcv_numRepeats = 1,
                                  nreps_energy = 1e4,
                                  eps = 1e-8,
                                  targetFeats = NULL,
                                  patience = 25,alpha_ = 0.05){


  ## Pre Process Input Data
  inputData$i  = 1:nrow(inputData)
  inputData = data.frame(Status = inputData[,1], inputData[,-1])
  inputData = dplyr::arrange(.data = inputData,Status)
  ids = inputData$i
  inputData =  dplyr::select(inputData,-i)
  inputData$Status = factor(inputData$Status)


  ## compute Logratio and Run DCV
  lrs.train = calcLogRatio(inputData)
  lrs.scaled = scale(lrs.train[,-1])
  lrs.train = data.frame(Status = lrs.train[,1],lrs.scaled)

  dcv  = diffCompVarRcpp::finalDCV(logRatioMatrix = lrs.train,includeInfoGain = dcv_useInfoGain,nfolds = dcv_nfold,numRepeats = dcv_numRepeats )
  dcv_ = dcv$rawDCV
  keep = dcv$lrs
  end = which.max(keep$nDistinct)
  keep = keep[1:end,]
  # select top N tsuch that every part is included at least once; subset that explains 100% of total variace
  trainData2 = subset(lrs.train,select = keep$Ratio)
  # compute max span tree for indepenet set
  allFeats = diffCompVarRcpp::mstAll(featMatrix = trainData2,dcvRanking = keep)
  labels = inputData[,1]
  message("Forward Selection ranking of logratio sets..")
  rollmean.df = data.frame()
  maxRLM = 0
  impTime = 0


  ###-----------------------------------------*
  ## All features
  ##-----------------------------------------*
  cn = colnames(allFeats)
  metrics = featureSlectionPerformance(data.frame(Status = labels,allFeats),plot_ = F,trainModel = F)
  allFeatureMetrics = metrics$performance
  d.compdat = parallelDist::parDist(as.matrix(allFeats),method = "euclidean")
  a.df = data.frame(Type = labels)
  mod = vegan::betadisper(d.compdat,group = labels)
  bd1 = vegan::permutest(mod,permutations = nreps_energy)



  ## automatic optimization mode selection
  ## if dispersion effects are detected the the search optimizes combinedF
  ## if dispersion not detected then the search optimizes the scaled F
  if(is.null(optimizationMetric)){
    optimizationMetric = if_else(bd1$tab$`Pr(>F)`[1]<alpha_,"combinedF","scale")
    message("dispersion p = ",round(bd1$tab$`Pr(>F)`[1],3)," ;  optmization is:  ", optimizationMetric)
  }

  ###----------------------------*
  ## Stepwise Forward Selection
  ###------------------------------*
  baseSet = 1:3
  cn = colnames(allFeats)
  ph = subset(allFeats,select =  cn[baseSet])

  ###----------------------------*
  ### Compute estat
  ###----------------------------*
  lbs = labels
  if(optimizationMetric == "scale"){

    d = parallelDist::parallelDist(as.matrix(ph))
    enf = energy::disco(x = d,factors = lbs,distance = T,R = 2,method = "discoF")
    estat = enf$statistic
    e_stat = estat
    newF = estat

  }else if (optimizationMetric == "norm"){

    newF = normalizedEnergy(lrMat = ph,labels = lbs)$H
    e_stat = newF

  }else if (optimizationMetric == "combinedF"){

    a.df = data.frame(Type = lbs)
    d.compdat = parallelDist::parDist(as.matrix(ph),method = "euclidean")
    pmv1 = vegan::adonis(d.compdat~Type,data = a.df,permutations = 2)
    f1 = pmv1$aov.tab$F.Model[1]
    mod = vegan::betadisper(d.compdat,group = lbs)
    bd1 = vegan::permutest(mod,permutations = 1)
    f2 = bd1$statistic
    newF = f1+f2
    e_stat = newF

  }else{

    metrics = featureSlectionPerformance(data.frame(Status = lbs,ph),plot_ = F,trainModel = F)
    e_stat = metrics$performance[optimizationMetric]
    estat = e_stat
    newF = metrics$performance[optimizationMetric]

  }
  ## set max value
  maxF = newF

  ## Run Forward Selection
  astat = data.frame()
  baseSet.list = list()
  astat = rbind(astat,data.frame(Value = maxF,optF = maxF,numFeats = length(baseSet)))
  baseSet.list[[1]] = baseSet
  bs = 2
  ss = 1

  end  = ncol(allFeats)

  if(end>4){
    for(x in 4:end){
      message(x, " of ", ncol(allFeats) )
      newSet = c(baseSet,x)
      cn = colnames(allFeats[newSet])
      ph = subset(allFeats,select =  cn)

      ###----------------------------*
      ### Compute estat
      ###----------------------------*
      if(optimizationMetric == "scale"){

        d = parallelDist::parallelDist(as.matrix(ph))
        enf = energy::disco(x = d,factors = lbs,distance = T,R = 2,method = "discoF")
        estat = enf$statistic
        e_stat = estat
        newF = estat

      }else if (optimizationMetric == "norm"){

        newF = normalizedEnergy(lrMat = ph,labels = lbs)$H
        e_stat = newF

      }else if (optimizationMetric == "combinedF"){

        a.df = data.frame(Type = lbs)
        d.compdat = parallelDist::parDist(as.matrix(ph),method = "euclidean")
        pmv1 = vegan::adonis(d.compdat~Type,data = a.df,permutations = 1)
        f1 = pmv1$aov.tab$F.Model[1]
        mod = vegan::betadisper(d.compdat,group = lbs)
        bd1 = vegan::permutest(mod,permutations = 1)
        f2 = bd1$statistic
        newF = f1+f2
        e_stat = newF

      }else{

        metrics = featureSlectionPerformance(data.frame(Status = lbs,ph),plot_ = F,trainModel = F)
        e_stat = metrics$performance[optimizationMetric]
        estat = e_stat
        newF = metrics$performance[optimizationMetric]

      }

      diff_ = newF-maxF
      if(diff_>eps){
        baseSet = c(baseSet,x)
        maxF = newF
        e_stat = newF
      }

      baseSet.list[[bs]] = baseSet
      bs = bs + 1


      ###---------------------------------------------------*
      ### Energy
      ###--------------------------------------------------*
      astat = rbind(astat,data.frame(Value = newF,optF = maxF,numFeats = length(baseSet)))

      #Breaking criteria - if target features are reached
      if(!is.null(targetFeats)){
        if(targetFeats<=length(baseSet)){
          break
        }
      }

      message("#Features = ",length(baseSet),", maxF = ",round(maxF,4)," , impTime = ",impTime)
      #stop if perm test and test statistic has been exceeded


      if(diff_>eps){
        impTime = 0
      }else{
        impTime = impTime+1
      }
      if(impTime>patience){
        break
      }
    }
  }



  ## get final subset
  avec = c(astat$Value)
  mx = which.max(astat$Value)
  pp = baseSet.list[[mx]]
  feats = allFeats[,pp]




  return(list(finalSubset = data.frame(Status = labels,feats),
              optimization_Metric = optimizationMetric,
              optimResult = astat,
              rawDCV = dcv_,retainedDCV = keep,
              mstSubset = allFeats))

}
