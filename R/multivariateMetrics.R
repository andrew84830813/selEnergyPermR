#' Wrapper around the Energy package Multisample E-statistic (Energy) Test of Equal Distributions
#'
#'
#' Modifies the energy::eqdist.etest to return permuation distribution
#'
#' @useDynLib selEnergyPermR
#'
#' @param x sorted matrix with n-samples
#' @param sizes class sizes
#' @param distance is a distance matrix
#' @param method method
#' @param R number of permuattions
#'
#' @examples
#' etest2()
#'
#' @references
#' Szekely, G. J. and Rizzo, M. L. (2004) Testing for Equal Distributions in High Dimension, InterStat, November (5).
#' M. L. Rizzo and G. J. Szekely (2010). DISCO Analysis: A Nonparametric Extension of Analysis of Variance, Annals of Applied Statistics, Vol. 4, No. 2, 1034-1055. http://dx.doi.org/10.1214/09-AOAS245
#' Szekely, G. J. (2000) Technical Report 03-05: E-statistics: Energy of Statistical Samples, Department of Mathematics and Statistics, Bowling Green State University.
#'
#' @return A data frame with all choose(p,2) log ratios
#'
#' @seealso [energy::eqdist.etest()]
#' @export
etest2 <-
  function (x, sizes, distance = FALSE, method = c("original",
                                                   "discoB", "discoF"), R)
  {
    method <- match.arg(method)
    if (method == "discoB" || method == "discoF") {
      g <- as.factor(rep(1:length(sizes), sizes))
      return(energy::disco(x, factors = g, distance = distance, index = 1,
                   R = R, method = method))
    }
    nsamples <- length(sizes)
    if (nsamples < 2)
      return(NA)
    if (min(sizes) < 1)
      return(NA)
    if (!is.null(attr(x, "Size")))
      distance <- TRUE
    x <- as.matrix(x)
    if (NROW(x) != sum(sizes))
      stop("nrow(x) should equal sum(sizes)")
    if (distance == FALSE && nrow(x) == ncol(x))
      warning("square data matrix with distance==FALSE")
    d <- NCOL(x)
    if (distance == TRUE)
      d <- 0
    str <- "Multivariate "
    if (d == 1)
      str <- "Univariate "
    if (d == 0)
      str <- ""
    e0 <- 0
    repl <- rep(0, R)
    pval <- 1
    b <- .C("ksample.e", x = as.double(t(x)), byrow = as.integer(1),
            nsamples = as.integer(nsamples), sizes = as.integer(sizes),
            dim = as.integer(d), R = as.integer(R), e0 = as.double(e0),
            e = as.double(repl), pval = as.double(pval), PACKAGE = "energy")
    names(b$e0) <- "E-statistic"
    sz <- paste(sizes, collapse = " ", sep = "")
    methodname <- paste(str, length(sizes), "-sample E-test of equal distributions",
                        sep = "")
    dataname <- paste("sample sizes ", sz, ", replicates ",
                      R, sep = "")
    e <- list(call = match.call(), method = methodname, statistic = b$e0, perms = b$e,
              p.value = b$pval, data.name = dataname)
    class(e) <- "htest"
    e
  }


#' Mulitvariate Statistics of a set of log-ratios
#'
#'
#' Computes various multivariate statistics to be used in permutation testing
#'
#' @useDynLib selEnergyPermR
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
#' featureSlectionPerformance()
#'
#' @references
#' Szekely, G. J. and Rizzo, M. L. (2004) Testing for Equal Distributions in High Dimension, InterStat, November (5).
#' M. L. Rizzo and G. J. Szekely (2010). DISCO Analysis: A Nonparametric Extension of Analysis of Variance, Annals of Applied Statistics, Vol. 4, No. 2, 1034-1055. http://dx.doi.org/10.1214/09-AOAS245
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

  tbl[,1] = as.character(tbl[,1])
  mat = dplyr::arrange(.data = data.frame(tbl),desc(Status))

  if(scaleData){
    mat[,-1] = scale(mat[,-1])
  }

  classes = c(mat[1,1],mat[nrow(mat),1])
  cc = combinat::combn2(as.character(classes))
  labels = tbl$Status
  c.df = data.frame(Count.labels = classes)
  W = data.frame(Count = (table(labels)),Weight = compositions::clo(table(labels)))
  W = dplyr::left_join(c.df,W)
  N = sum(W$Count.Freq)
  d = parallelDist::parDist(as.matrix(mat[,-1])) ## p here is eqv to the alpha in the (rizzo/szekely - energy distanmce advanced review feb 2016 review paper)
  energy_aitch = energy::eqdist.e(x = d,sizes = W$Count.Freq,distance = T)
  # estat = energy_aitch$statistic
  # ii = energy_aitch$perms
  allFeatures_eqdist = energy_aitch#as.numeric((estat-mean(ii))/sd(ii)) ##



  #dispersion test
  labels = mat[,1]
  mod = vegan::betadisper(d,group = labels)
  bd = vegan::permutest(mod,permutations  = 2)
  #permanova
  a.df = data.frame(Type = labels)
  pmv = vegan::adonis2(d~Type,data = a.df,permutations = 0)
  ## ANOSIM
  ano = vegan::anosim(x = d,grouping = labels,permutations = 0)

  ###----------------------------
  ## mds projection
  ###----------------------------
  mds = cmdscale(d,k = 2 )
  mds_dist = dist((mds))
  normE_mds = normalizedEnergy(mds,labels)
  mod = vegan::betadisper(mds_dist,group = labels)
  bd_mds = vegan::permutest(mod)
  pmv_mds = vegan::adonis2(mds_dist~Type,data = a.df,permutations = 0)

  #scaled mds e stat
  energy_aitch = energy::eqdist.e(x = mds_dist,sizes = W$Count.Freq,distance = T)
  estat = energy_aitch$statistic
  # ii = energy_aitch$perms
  EnergyE_mds = estat#as.numeric((estat-mean(ii))/sd(ii)) ##



  ## distance cov and cor test for indpedence
  tt = as.matrix(as.integer(as.factor(mat[,1]))-1)
  bcd = energy::bcdcor(as.matrix(mat[,-1]),as.matrix(tt))
  dcr = energy::dcor(as.matrix(mat[,-1]),as.matrix(tt))


  ## Energy disco
  mat.ph = data.frame(tbl) %>%
    dplyr::arrange(desc(Status))
  mat.ph[,1] = factor(mat.ph[,1])
  levs = data.frame(Labels = unique(mat.ph[,1]))
  labels = mat.ph[,1]
  tb = data.frame((table(mat.ph[,1] )))
  colnames(tb)[1] = "Labels"
  sz = dplyr::left_join(levs,tb)
  d1 = dist(mat.ph[,-1])
  #energy::eqdist.etest(d,sizes = sz$Freq,distance = T,R = 100,method = "discoB")
  enf = energy::disco(x = d1,factors = mat.ph$Status,distance = T,R = 2)



  ns = NA
  g = NA
  netStr = NA
  if(plot_){

    ##-----------------------------------------*
    ## LR Network ####
    ##-----------------------------------------*
    feature.df = tbl
    #### Kruskall Test
    krus.test_df = data.frame()
    # tbl =  mst.empirical$features[[1]]
    lrs_ = data.frame(feature.df[,-1])
    colnames(lrs_) = colnames(feature.df)[-1]
    cnames = colnames(lrs_)
    for(i in 1:ncol(lrs_)){
      ph = kruskal.test(x = lrs_[,i],g  = factor(feature.df[,1]))
      ph = data.frame(Ratio =cnames[i],pval = ph$p.value,Statistic = ph$statistic )
      krus.test_df = rbind(krus.test_df,ph)
    }
    #network construction
    krus.test_df$p.adjust = p.adjust(krus.test_df$pval,method = "BH")
    pval_level = 0.05
    fdrLevel =  max(krus.test_df$p.adjust[krus.test_df$pval <= pval_level])
    krus.test_df$signf = dplyr::if_else(krus.test_df$p.adjust>0.05,F,T)
    #Importance df
    imp.df = data.frame(Ratio = krus.test_df$Ratio,Imp = krus.test_df$Statistic)
    keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
    el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
    g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
    g = igraph::simplify(g, remove.loops = TRUE,
                         edge.attr.comb = igraph_opt("edge.attr.comb"))
    E(g)$weight = dplyr::if_else((imp.df$Imp)<0,0,(imp.df$Imp))

    ns = igraph::diameter(g)/length(igraph::strength(g))
    netStr = length(igraph::strength(g))


    plot(g,layout = igraph::layout_with_fr,vertex.size = log(1/compositions::clo(igraph::strength(g)))+1,
         vertex.label.cex = .75,edge.curved = .2,edge.width = igraph::E(g)$weight*.15,
         edge.arrow.size = .25,edge.arrow.width = 1)
    plot(mod,label = F)
  }


  ##z test combined beta and perma
  a.df = data.frame(Type = mat[,1])
  pmv1 = vegan::adonis(d~Type,data = a.df,permutations = combinedF_reps)
  f1 = (pmv1$aov.tab$F.Model[1]-mean(pmv1$f.perms)) / sd(pmv1$f.perms)
  mod = vegan::betadisper(as.dist(d),group = mat[,1])
  bd1 = vegan::permutest(mod,permutations = combinedF_reps)
  f2 = as.numeric((bd1$statistic-mean(bd1$perm)) / sd(bd1$perm))



  return(list(
    performance =  data.frame(Method = Method_Name,
                              Scenario = scenario,
                              NumRatios = ncol(tbl)-1,
                              NumParts = netStr,raw_energyStat = estat,
                              mds_normEnergy = normE_mds$H,
                              mds_rawEnergy = normE_mds$Estat,
                              scaleEStat_mds = EnergyE_mds,
                              PermanovaF = pmv$F[1],
                              PermanovaF_mds = pmv_mds$F[1],
                              betaDispF = bd$tab$F[1],
                              betaDispF_mds = bd_mds$`F value`[1],
                              AnosimR = ano$statistic,
                              EnergyF = enf$statistic,
                              energyCor = dcr,
                              energyBiasedCorrecCor = bcd,
                              combinedF = as.numeric(f2+f1),
                              eqdist = allFeatures_eqdist,
                              netDiamByStrength = ns),
    graph = g
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
#' normalizedEnergy()
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
    mat = arrange(data.frame(Labels = labels,lrMat),desc(Labels))

    classes = unique(labels)
    cc = combinat::combn2(as.character(classes))
    W = data.frame(Count = (table(labels)),Weight = compositions::clo(table(labels)))
    N = sum(W$Count.Freq)

    Estat = c()
    H = c()
    for(c in 1:nrow(cc)){
      l = c(cc[c,1],cc[c,2])
      mat.ph = mat[mat$Labels%in%l,]
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
