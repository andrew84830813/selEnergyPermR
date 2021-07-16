
#' Mulitvariate Statistics of a set of log-ratios
#'
#'
#' Computes various multivariate statistics to be used in permutation testing
#'
#'
#' @importFrom magrittr %>%
#'
#' @param tbl sorted matrix with n-samples
#' @param plot_ class sizes
#' @param nreps_energy is a distance matrix
#' @param Method_Name method
#' @param scaleData a data.frame where the first column contains class labels and the remaining p columns are features
#' @param scenario a data.frame where the first column contains class labels and the remaining p columns are features
#' @param combinedF_reps a data.frame where the first column contains class labels and the remaining p columns are features
#'
#' @examples
#' \dontrun{
#' featureSlectionPerformance()
#' }
#'
#' @references
#' Szekely, G. J. and Rizzo, M. L. (2004) Testing for Equal Distributions in High Dimension, InterStat, November (5).
#'
#' M. L. Rizzo and G. J. Szekely (2010). DISCO Analysis: A Nonparametric Extension of Analysis of Variance, Annals of Applied Statistics, Vol. 4, No. 2, 1034-1055. http://dx.doi.org/10.1214/09-AOAS245
#'
#' Szekely, G. J. (2000) Technical Report 03-05: E-statistics: Energy of Statistical Samples, Department of Mathematics and Statistics, Bowling Green State University.
#'
#' @return A data frame with all choose(p,2) log ratios
#'
#' @export
featureSlectionPerformance = function(tbl,
                                      plot_=F,
                                      nreps_energy = 1e5,
                                      Method_Name = "Feat_Selection",
                                      scaleData = T,
                                      scenario = "Scenario_X",
                                      combinedF_reps = 5000){

  ##-------------------------*
  ## Scaled E
  ##-------------------------*
  Status = NULL
  tbl[,1] = as.character(tbl[,1])
  mat = dplyr::arrange(.data = data.frame(tbl),dplyr::desc(Status))

  if(scaleData){
    mat[,-1] = scale(mat[,-1])
  }


  ## Energy disco
  mat[,1] = factor(mat[,1])
  levs = data.frame(Labels = unique(mat[,1]))
  labels = mat[,1]
  tb = data.frame((table(mat[,1] )))
  colnames(tb)[1] = "Labels"
  sz = dplyr::left_join(levs,tb)
  d = stats::dist(mat[,-1])
  enf = energy::disco(x = d,factors = mat[,1],distance = T,R = 2)




  ##z test combined beta and perma
  a.df = data.frame(Type = mat[,1])
  pmv1 = vegan::adonis(d~Type,data = a.df,permutations = combinedF_reps)
  f1 = (pmv1$aov.tab$F.Model[1]-mean(pmv1$f.perms)) / stats::sd(pmv1$f.perms)
  mod = vegan::betadisper( stats::as.dist(d),group = mat[,1])
  bd1 = vegan::permutest(mod,permutations = combinedF_reps)
  f2 = as.numeric((bd1$statistic-mean(bd1$perm)) / stats::sd(bd1$perm))



  return(list(
    performance =  data.frame(Method = Method_Name,
                              Scenario = scenario,
                              NumRatios = ncol(tbl)-1,
                              EnergyF = enf$statistic,
                              combinedF = as.numeric(f2+f1))
  )
  )

}






#' Normalized Energy Statistic
#'
#'
#' Normalizes the enrgy statistic as described in Rizzo (2000)
#'
#' @importFrom magrittr %>%
#'
#' @param lrMat sorted matrix with n-samples
#' @param labels class sizes
#'
#' @examples
#' \dontrun{
#' normalizedEnergy()
#' }
#'
#' @references
#' Szekely, G. J. and Rizzo, M. L. (2004) Testing for Equal Distributions in High Dimension, InterStat, November (5).
#'
#' M. L. Rizzo and G. J. Szekely (2010). DISCO Analysis: A Nonparametric Extension of Analysis of Variance, Annals of Applied Statistics, Vol. 4, No. 2, 1034-1055. http://dx.doi.org/10.1214/09-AOAS245
#'
#' Szekely, G. J. (2000) Technical Report 03-05: E-statistics: Energy of Statistical Samples, Department of Mathematics and Statistics, Bowling Green State University.
#'
#' @return A data frame with all choose(p,2) log ratios
#'
#' @export
normalizedEnergy <-
  function(lrMat,labels){
    Labels = NULL
    mat = dplyr::arrange(data.frame(Labels = labels,lrMat),dplyr::desc(Labels))

    classes = unique(labels)
    cc = combinat::combn2(as.character(classes))
    W = data.frame(Count = (table(labels)),Weight = compositions::clo(table(labels)))
    N = sum(W$Count.Freq)

    Estat = c()
    H = c()
    for(c in 1:nrow(cc)){
      l = c(cc[c,1],cc[c,2])
      mat.ph = mat[mat$Labels %in% l,]
      mat.ph$Labels = factor(mat.ph$Labels)
      levs = data.frame(Labels = unique(mat.ph$Labels))
      labels = mat.ph$Labels
      tb = data.frame((table(mat.ph$Labels )))
      colnames(tb)[1] = "Labels"
      sz = dplyr::left_join(levs,tb)
      ## Compute Distance
      d = parallelDist::parDist(as.matrix(mat.ph[,-1]),method = "minkowski",p=1.99) ## p here is eqv to the alpha in the (rizzo/szekely - energy distanmce advanced review feb 2016 review paper)
      d = as.matrix(d)
      #energy::eqdist.e(x = d,distance = T,sizes = sz$Freq,method = "discoB")

      ## Compute A
      a1  = 1:sz[1,2]
      a2  = (sz[1,2]+1):nrow(mat.ph)
      A = d[a1,a2]
      A = as.vector(A)
      A = mean(A)
      ## Compute B
      B = d[a1,a1]
      B= as.vector(B)
      B = sum(B) / length(a1)^2
      ## Compute C
      C = d[a2,a2]
      C = as.vector(C)
      C = mean(C)
      ## Compute and Norm Energy
      E = 2*A-B-C
      T_ = (length(a1)*length(a2)) / ((length(a1)+length(a2)))
      estat =E*T_
      h = E/(2*A)
      ## Compute weights
      w = sum(sz$Freq)/(2*N)
      ## Weight stats
      estat = w*estat
      h = h*w

      H[c] = h
      Estat[c] = estat

    }

    h = sum(H)
    e = sum(Estat)


    return(list(H = h, Estat = e))

  }
