require(jsonlite)
require(futile.logger)
require(dplyr)

## helper functions
rows <- function(x) lapply(seq_len(nrow(x)), function(i) lapply(x,"[",i))
btest <- function (x, p, n){binom.test(x,n,p, alternative='greater')$p.value} 

#### functions imported here from BLMA package
averageCumulative <- function (x) {
    n <- length(x)
    a <- mean(x)
    1/factorial(n) * sum(sapply(0:floor(n * a), function(k) (-1)^k * 
        choose(n, k) * (n * a - k)^(n)))
}
addCLT <- function (x) {
    if (sum(is.na(x)) > 0) 
        NA
    else {
        n <- length(x)
        if (n <= 20) {
            averageCumulative(x)
        }
        else {
            pnorm(mean(x), 1/2, sqrt(1/(12 * n)), lower=TRUE)
        }
    }
}
fisherMethod <- function (x) {
    if (sum(is.na(x)) > 0) 
        NA
    else pchisq(-2 * sum(log(x)), df=2 * length(x), lower=FALSE)
}
stoufferMethod <- function (x) {
    if (sum(is.na(x)) > 0) 
        NA
    else pnorm(sum(qnorm(x))/sqrt(length(x)))
}

minP <- function(x) {
    if (sum(is.na(x)) > 0) 
      NA
    else min(x)
}

# Combining p-value method used in SPIA paper
pGMethod <- function (x) { 
    if (sum(is.na(x)) > 0) 
      NA
    else { prod(x) - prod(x)*log(prod(x)) }
}

# different options for combining p-value as provided by BLMA package
combFuncMap <- list (
                add=function(x, y, z, t){ifelse(sign(z) == sign(t) && abs(z)>=abs(t), addCLT(c(x,y)), y)},
                fisher=function(x, y, z, t){ifelse(sign(z) == sign(t) && abs(z)>=abs(t), fisherMethod(c(x,y)) ,y)},
                stouffer=function(x, y, z, t){ifelse(sign(z) == sign(t) && abs(z)>=abs(t), stoufferMethod(c(x,y)), y)},
                minP=function(x, y, z, t){ifelse(sign(z) == sign(t) && abs(z)>=abs(t), minP(c(x,y)), y)},
                pGMethod=function(x, y, z, t){ifelse(sign(z) == sign(t) && abs(z)>=abs(t), pGMethod(c(x,y)), y)}
                )
#####
# network analysis
#####
###
# computes upstream regulators predicted as activated and predicted as inhibited
# combines p-value for z-score with ora p-value if |zscore|>=2 with fisher method (alternative methods for combining pvalues are provided from BLMA)
# 
# count_de - number of de targets
# count_all - number of measured targets
# 
# pv - naive approach for computing p-value just by looking at number of DE targets vs number of all measured targets  (ignore edges)
# pv_zscore - ASSUME z-score ~ N(0,1) for large enough number of targets
# pv_up - p-value - predicted as activated using phyper
# pv_down - p-value - predicted as inhibted using phyper
# pv_binom_up - p-value - predicted as activated using binomial test 
# pv_binom_down - p-value - predicted as inhibted using binomial test
# pv_comb_up - combined pv_up and pv_zscore when z-score>=2 or pv_up otherwise
# pv_comb_down - combined pv_down and pv_zcore when z-score<=-2 or pv_down otherwise
upstreamAnalysis <- function(inputData, allData, network, measuredTargets, deTargets, thr=list(edgeConfidence=400, significantZscore=2), combFunc=combFuncMap$fisher, overlapFunc, signFunc, tgtGeneIdColumn){
  
  if (nrow(inputData) == 0) {
    flog.error('no input data. Size= %s. Stop.', nrow(inputData))
    return()
  }
  flog.info('size of input data: %s', nrow(inputData))

  deGenes <- inputData$entrez
  deTargets <- backgroundOverlap(inputData, network, thr)

  a <- subset(network, maxScore>=thr$edgeConfidence & target_gene_id %in% deGenes)
  if (nrow(a) == 0) {
    flog.error('no upstream regulators. Size= %s. Stop.', nrow(a))
    return()
  }
  p <- aggregate(a$target_gene_id, by=list(entrez=a$src_gene_id, symbol=a$src_gene_symbol), FUN=function(a){length(unique(a))})
  colnames(p) <- c('entrez','symbol','count_de')

  flog.info("collected upstream regulators")
  q <- subset(network, maxScore>=thr$edgeConfidence & src_gene_id %in% p$entrez)
  q <- q[,c('src_gene_id','target_gene_id','effect')]

  flog.info("collecting background for each upstream regulator")
  qAllData <- merge(x=q, y=allData, by.x="target_gene_id", by.y='entrez', all.x=FALSE)
  
  countAllData <- aggregate(
      x=qAllData$target_gene_id, by=list(src_gene_id=qAllData$src_gene_id),
      FUN=function(x){length(unique(x))})
  p <- merge(x=p, y=countAllData, by.x="entrez", by.y='src_gene_id', all.x=TRUE)
  p$count_all <- p$x

  qAllData$sign <- ifelse(qAllData$effect == 'activation',1,-1) * sign(qAllData$logfc)
  z_p <- NA
  if (length(qAllData[qAllData$sign == 1,'target_gene_id'])>0) {
    countP <- aggregate(
      x=qAllData[qAllData$sign == 1,'target_gene_id'], by=list(src_gene_id=qAllData[qAllData$sign == 1,'src_gene_id']),
      FUN=function(x){length(unique(x))})
    z_p=merge(x=p, y=countP, by.x="entrez", by.y='src_gene_id', all.x=TRUE) 
  }
  p$count_p <- ifelse(is.na(z_p)||is.na(z_p$x.y), 0, z_p$x.y)
  
  z_n <- NA
  if (length(qAllData[qAllData$sign == -1,'target_gene_id'])>0) {
    countN <- aggregate(
      x=qAllData[qAllData$sign == -1,'target_gene_id'], by=list(src_gene_id=qAllData[qAllData$sign == -1,'src_gene_id']),
      FUN=function(x){length(unique(x))})
  
    z_n <- merge(x=p, y=countN, by.x="entrez", by.y='src_gene_id', all.x=TRUE)
  }
  p$count_n <- ifelse(is.na(z_n)||is.na(z_n$x.y), 0, z_n$x.y)

  flog.info("computing upstream regulators p-values - naive approach")
  p["pv"] <- phyper(p$count_de-1, deTargets$count, measuredTargets$count-deTargets$count, p$count_all, lower.tail=FALSE)
  p["pv_fdr"] <- p.adjust(p$pv, method="fdr")
  p["pv_bonferroni"] <- p.adjust(p$pv, method="bonferroni")
  p["pv_binom"] <- mapply(btest, p$count_de, deTargets$count/measuredTargets$count, p$count_all)

  resultNEO4J <- p[,c("entrez","symbol","count_de","count_all","pv","pv_fdr","pv_bonferroni","pv_binom","count_n","count_p")]
  flog.info("computing upstream regulators zscore")
  
  resultNEO4J["zscore"] <- zscoreUpstream(network, resultNEO4J$entrez, inputData, overlapFunc, signFunc, tgtGeneIdColumn, thr)
  resultNEO4J["pv_zscore"] <- mapply(function(z){pnorm(-abs(z))}, resultNEO4J$zscore)
  resultNEO4J["pv_zscore_fdr"] <- p.adjust(resultNEO4J$pv_zscore, method="fdr")
  resultNEO4J["pv_zscore_bonferroni"] <- p.adjust(resultNEO4J$pv_zscore, method="bonferroni")

  flog.info("computing upstream regulators upstream p-values")
  r <- pvUpstream(resultNEO4J, inputData, network, measuredTargets, deTargets, thr)
  resultNEO4J <- merge(x=resultNEO4J, y=r, by="entrez", by.y='entrez', all.x=TRUE)
  
  resultNEO4J["pv_comb_up"] <- mapply(combFunc, resultNEO4J$pv_zscore, resultNEO4J$pv_up, resultNEO4J$zscore, rep(thr$significantZscore, length(resultNEO4J$pv_zscore)))
  
  resultNEO4J["pv_comb_down"] <- mapply(combFunc, resultNEO4J$pv_zscore, resultNEO4J$pv_down, resultNEO4J$zscore, rep(-1*thr$significantZscore, length(resultNEO4J$pv_zscore))) 
  
  resultNEO4J["pv_up_fdr"] <- p.adjust(resultNEO4J$pv_up, method="fdr")
  resultNEO4J["pv_up_bonferroni"] <- p.adjust(resultNEO4J$pv_up, method="bonferroni")
  
  resultNEO4J["pv_down_fdr"] <- p.adjust(resultNEO4J$pv_down, method="fdr")
  resultNEO4J["pv_down_bonferroni"] <- p.adjust(resultNEO4J$pv_down, method="bonferroni")
  
  resultNEO4J["pv_comb_up_fdr"] <- p.adjust(resultNEO4J$pv_comb_up, method="fdr")
  resultNEO4J["pv_comb_up_bonferroni"] <- p.adjust(resultNEO4J$pv_comb_up, method="bonferroni")
  
  resultNEO4J["pv_comb_down_fdr"] <- p.adjust(resultNEO4J$pv_comb_down, method="fdr")
  resultNEO4J["pv_comb_down_bonferroni"] <- p.adjust(resultNEO4J$pv_comb_down, method="bonferroni")
  resultNEO4J=resultNEO4J[order(resultNEO4J$pv_comb_up),]

  return(resultNEO4J)
}

getUpstreamType <- function(identifier, suffix="db_id.csv"){
  analysisDataFile <- paste(env[["ANALYSIS_TYPE"]], identifier, suffix, sep='_')
  file <- file.path(env[['DATA_FOLDER']], analysisDataFile)
  inputData <- read.csv(file=file, header=TRUE, stringsAsFactors=FALSE)
  return (inputData$type[1])
}

drugUpstreamAnalysis <- function (inputData, allData, network, measuredTargets, deTargets, thr=list(edgeConfidence=400, significantZscore=2), combFunc=combFuncMap$fisher, overlapFunc, signFunc, tgtGeneIdColumn) {
  
  if (nrow(inputData) == 0) {
    flog.error('no input data. Size= %s. Stop.', nrow(inputData))
    return()
  }
  flog.info('size of input data: %s', nrow(inputData))
  
  deGenes <- inputData$entrez
  
  DENet <- subset(network, GeneID %in% deGenes) # list of all edges that have DE genes as targets. 
  if (nrow(DENet) == 0) { 
    # If it is NULL, then no upstream regulator found
    flog.error('no upstream regulators. Size= %s. Stop.', nrow(DENet))
    return()
  }
  ChemDE <- aggregate(DENet$GeneID, by=list(DENet$ChemicalID), FUN=function(DENet){length(unique(DENet))}) # p: number of DE genes for each regulator
  # p sumarizes number of DE under each upregulator
  colnames(ChemDE) <- c('ChemicalID', 'count_de')
  
  flog.info("collected upstream regulators")
  ChemNet <- subset(network, ChemicalID %in% ChemDE$ChemicalID) # sub-network of all regulators
  ChemNet <- ChemNet[,c('ChemicalID', 'GeneID','edge_type')]
  
  flog.info("collecting background for each upstream regulator")
  ChemNetData <- merge(x=ChemNet, y=allData, by.x="GeneID", by.y='entrez', all.x=FALSE)
  
  ChemAllData <- aggregate(x=ChemNetData$GeneID,
      by=list(ChemicalID=ChemNetData$ChemicalID),
      FUN=function(x){length(unique(x))})
  colnames(ChemAllData) <- c('ChemicalID', 'count')
  
  ChemSummary <- merge(x=ChemDE, y=ChemAllData, by.x="ChemicalID", by.y='ChemicalID', all.x=TRUE)
  ChemDE$count_all <- ChemSummary$count
  
  ChemNetData$sign <- ChemNetData$edge_type * sign(ChemNetData$logfc) # sign(edge) * sign(fc)
  countP=aggregate( # number of "white balls" in Hypothesis 1 for each regulators
      x=ChemNetData[ChemNetData$sign == 1,'GeneID'], by=list(ChemicalID=ChemNetData[ChemNetData$sign == 1,'ChemicalID']),
      FUN=function(x){length(unique(x))})
  z <- merge(x=ChemDE, y=countP, by.x="ChemicalID", by.y='ChemicalID', all.x=TRUE) 
  ChemDE$count_p <- ifelse(is.na(z$x), 0, z$x)
  
  countN <- aggregate( # number of "white balls" in Hypothesis 2 for each regulators
      x=ChemNetData[ChemNetData$sign == -1,'GeneID'], by=list(ChemicalID=ChemNetData[ChemNetData$sign == -1, "ChemicalID"]),
      FUN=function(x){length(unique(x))})
  z <- merge(x=ChemDE, y=countN, by.x="ChemicalID", by.y="ChemicalID", all.x=TRUE) 
  ChemDE$count_n=ifelse(is.na(z$x), 0, z$x)
  
  flog.info("computing upstream regulators p-values - naive approach")
  ChemDE["pv"] <- phyper(ChemDE$count_de-1, deTargets$count, measuredTargets$count-deTargets$count, ChemDE$count_all, lower.tail=FALSE)
  ChemDE["pv_fdr"] <- p.adjust(ChemDE$pv, method="fdr")
  ChemDE["pv_bonferroni"] <- p.adjust(ChemDE$pv, method="bonferroni")
  ChemDE["pv_binom"] <- mapply(btest, ChemDE$count_de, deTargets$count/measuredTargets$count, ChemDE$count_all)
  
  resultNEO4J <- ChemDE[,c("ChemicalID", "count_de", "count_all" , "pv", "pv_fdr", "pv_bonferroni", "pv_binom", "count_n", "count_p")]
  flog.info("computing upstream regulators zscore") # First evidence (zscore and pvalue based on zscore)

  resultNEO4J["zscore"] <- zscoreUpstream(ChemNetData, resultNEO4J$ChemicalID, inputData, overlapFunc, signFunc, tgtGeneIdColumn, thr)
  resultNEO4J["pv_zscore"] <- mapply(function(z){pnorm(-abs(z))}, resultNEO4J$zscore)
  resultNEO4J["pv_zscore_fdr"] <- p.adjust(resultNEO4J$pv_zscore, method="fdr")
  resultNEO4J["pv_zscore_bonferroni"] <- p.adjust(resultNEO4J$pv_zscore, method="bonferroni")



  flog.info("computing upstream regulators upstream p-values") # Second evidence (pvalue upstream)
  r <- drugPvUpstream(resultNEO4J, inputData, network, measuredTargets, deTargets)
  resultNEO4J <- merge(x=resultNEO4J, y=r, by="ChemicalID", by.y='ChemicalID', all.x=TRUE)
  
  resultNEO4J["pv_up_fdr"] <- p.adjust(resultNEO4J$pv_up, method="fdr")
  resultNEO4J["pv_up_bonferroni"] <- p.adjust(resultNEO4J$pv_up, method="bonferroni")
  
  resultNEO4J["pv_down_fdr"] <- p.adjust(resultNEO4J$pv_down, method="fdr")
  resultNEO4J["pv_down_bonferroni"] <- p.adjust(resultNEO4J$pv_down, method="bonferroni")

  #Combining p-values using Fisher method
  resultNEO4J["pv_comb_up"] <- mapply(combFunc, resultNEO4J[,"pv_zscore"], resultNEO4J[,"pv_up"], resultNEO4J[,"zscore"], rep(thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  resultNEO4J["pv_comb_up_fdr"] <- p.adjust(resultNEO4J$pv_comb_up, method="fdr")
  resultNEO4J["pv_comb_up_bonferroni"] <- p.adjust(resultNEO4J$pv_comb_up, method="bonferroni")
  resultNEO4J["pv_comb_down"] <- mapply(combFunc, resultNEO4J[,"pv_zscore"], resultNEO4J[,"pv_down"], resultNEO4J[,"zscore"], rep(-1*thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  resultNEO4J["pv_comb_down_fdr"] <- p.adjust(resultNEO4J$pv_comb_down, method="fdr")
  resultNEO4J["pv_comb_down_bonferroni"] <- p.adjust(resultNEO4J$pv_comb_down, method="bonferroni")

  return(resultNEO4J)
}

# given a data frame of measured genes, 
# 
# # count - number of genes that overlap in the network
# # count_n - number of target genes that contradict the sign of an incomming edge
# # count_p - number of target genes that are in sync with the sign of an incomming edge 
# sign(gene)=sign(g$logfc), 
# sign (edge)=-1 for inhibition, or +1 for activation
backgroundOverlap <- function(allData, network, thr=list(edgeConfidence=400)) {

    if(upstream_type == "gene")
    {
      p <- subset(network, maxScore>=thr$edgeConfidence & target_gene_id %in% allData$entrez)
      overlap_data <- merge(x=p, y=allData, by.x="target_gene_id", by.y='entrez', all.x=FALSE)
      overlap_data$sign <- with(overlap_data, ifelse(effect == 'activation',1,ifelse(effect == 'inhibition',-1,0))*sign(logfc))
      r <- list( 
        count=length(unique(overlap_data$target_gene_id)),
        count_n=length(unique(subset(overlap_data, sign == -1)$target_gene_id)), 
        count_p=length(unique(subset(overlap_data, sign == 1)$target_gene_id)))
    } else {
      p <- subset(network, GeneID %in% allData$entrez)
      overlap_data <- merge(x=p, y=allData, by.x="GeneID", by.y='entrez', all.x=FALSE)
      overlap_data$sign <- with(overlap_data, sign(edge_type)*sign(logfc))
      r <- list ( 
        count=length(unique(overlap_data$GeneID)),
        count_n=length(unique(subset(overlap_data, sign == -1)$GeneID)), 
        count_p=length(unique(subset(overlap_data, sign == 1)$GeneID))
      ) 
    } 

    flog.info('size of background data: count=%s, count_n=%s, count_p=%s', r$count, r$count_n, r$count_p)
    return(r)
}

# computes pv_up, pv_down, pv_binom_up and pv_binom_down as defined above
# writes intermediate counts used to compute the p-values for each upstream regulator into a file 
pvUpstream <- function(upregulators, inputData, network, measuredTargets, deOverlap, thr=list(edgeConfidence=400)){
  analysisDatailsFile <- paste(env[["ANALYSIS_TYPE"]], identifier, upstream_type, "details.csv", sep='_')
  for (i in 1:nrow(upregulators)) {

    overlap <- subset(network, maxScore>=thr$edgeConfidence & target_gene_id %in% inputData$entrez & src_gene_id == upregulators[i,]$entrez)
    overlap <- overlap[,c('src_gene_id','target_gene_id','maxScore','effect')]
    if (!is.null(overlap)) {
      v<-merge(x=overlap, y=inputData, by.x="target_gene_id", by.y='entrez', all.x=TRUE)
      v['sign']<-ifelse(v$effect == 'activation',1,-1) * sign(v$logfc)
      upregulators[i,'count_up'] <- length(unique(subset(v, sign == 1)$target_gene_id))
      upregulators[i,'count_down'] <- length(unique(subset(v, sign == -1)$target_gene_id))
      upregulators[i,'pv_up'] <- phyper(upregulators[i,'count_up'] -1, deOverlap$count_p, measuredTargets$count-deOverlap$count_p, upregulators[i,]$count_all, lower.tail=FALSE)
      #upregulators[i,'pv_up'] <- phyper(upregulators[i,'count_up'] -1, upregulators[i,]$count_p, upregulators[i,]$count_all-upregulators[i,]$count_p, upregulators[i,]$count_de, lower.tail=FALSE)
      upregulators[i,'pv_down'] <- phyper(upregulators[i,'count_down'] -1, deOverlap$count_n, measuredTargets$count-deOverlap$count_n, upregulators[i,]$count_all, lower.tail=FALSE)
      #upregulators[i,'pv_down'] <- phyper(upregulators[i,'count_down'] -1, upregulators[i,]$count_n, upregulators[i,]$count_all-upregulators[i,]$count_n, upregulators[i,]$count_de, lower.tail=FALSE)
      upregulators[i, 'pv_binom_up'] <- binom.test(x=upregulators[i,'count_up'], n=upregulators[i,]$count_all, p=deOverlap$count_p/measuredTargets$count, alternative='greater')$p.value
      upregulators[i, 'pv_binom_down'] <- binom.test(x=upregulators[i,'count_down'], n=upregulators[i,]$count_all, p=deOverlap$count_n/measuredTargets$count, alternative='greater')$p.value
      
	  write.table(v, file=file.path(env[["DATA_FOLDER"]], analysisDatailsFile), row.names=FALSE, na ="", sep=",", col.names=(i == 1), append=(i>1))
    }
  }
  
  return(upregulators[,c('entrez', 'count_up','pv_up','count_down', 'pv_down', 'pv_binom_up', 'pv_binom_down')])
}

################ p-value of second evidence ############
# pvUpstream <- function(upregulators, inputData, network, measuredTargets, deOverlap, thr=list(edgeConfidence=400)){
drugPvUpstream <- function(resultNEO4J, inputData, network, measuredTargets, deTargets){
  # resultNE04J is the table of regulators' information
  analysisDatailsFile <- paste(env[["ANALYSIS_TYPE"]], identifier, upstream_type, "details.csv", sep='_')
  for (i in 1:nrow(resultNEO4J)) {
    overlap <- subset(network,  GeneID %in% inputData$entrez & ChemicalID == resultNEO4J[i,]$ChemicalID)
    
    if (!is.null(overlap)) {
      v <- merge(x=overlap, y=inputData, by.x="GeneID", by.y='entrez', all.x=TRUE)
      v['sign'] <- v$edge_type * sign(v$logfc)
      resultNEO4J[i,'count_up'] <- length(unique(subset(v, sign == 1)$GeneID))
      resultNEO4J[i,'count_down'] <- length(unique(subset(v, sign == -1)$GeneID))
      resultNEO4J[i,'pv_up'] <- phyper(resultNEO4J[i,'count_up'] -1, deTargets$count_p, measuredTargets$count-deTargets$count_p, resultNEO4J[i,]$count_all, lower.tail=FALSE)
      #upregulators[i,'pv_up'] <- phyper(upregulators[i,'count_up'] -1, upregulators[i,]$count_p, upregulators[i,]$count_all-upregulators[i,]$count_p, upregulators[i,]$count_de, lower.tail=FALSE)
      resultNEO4J[i,'pv_down'] <- phyper(resultNEO4J[i,'count_down'] -1, deTargets$count_n, measuredTargets$count-deTargets$count_n, resultNEO4J[i,]$count_all, lower.tail=FALSE)
      #upregulators[i,'pv_down'] <- phyper(upregulators[i,'count_down'] -1, upregulators[i,]$count_n, upregulators[i,]$count_all-upregulators[i,]$count_n, upregulators[i,]$count_de, lower.tail=FALSE)
      resultNEO4J[i, 'pv_binom_up'] <- binom.test(x=resultNEO4J[i,'count_up'], n=resultNEO4J[i,]$count_all, p=deTargets$count_p/measuredTargets$count, alternative='greater')$p.value
      resultNEO4J[i, 'pv_binom_down'] <- binom.test(x=resultNEO4J[i,'count_down'], n=resultNEO4J[i,]$count_all, p=deTargets$count_n/measuredTargets$count, alternative='greater')$p.value
      
      write.table(v[,c("id","GeneID","edge_type","logfc","adjpv","symbol","sign")], file=file.path(env[["DATA_FOLDER"]], analysisDatailsFile), row.names=FALSE, na ="", sep=",", col.names=(i == 1), append=(i>1))
    }
  }
  
  return(resultNEO4J[,c('ChemicalID', 'count_up','pv_up','count_down', 'pv_down', 'pv_binom_up', 'pv_binom_down')])
}


geneOverlapFunc <- function(upregulator, network, inputdata, thr=list(edgeConfidence=400)) {
	subset(network, maxScore>=thr$edgeConfidence & target_gene_id %in% inputData$entrez & src_gene_id == upregulator)
}

geneSignFunc <- function(overlap, j) {
	sR = 0;
    if(overlap[j,"effect"] == "activation") {
      sR = 1
    } else if(overlap[j,"effect"] == "inhibition") {
      sR = -1
    }
    return(sR)
}

drugOverlapFunc <- function(upregulator, network, inputdata, thr=list(edgeConfidence=400)) {
	subset(network, GeneID %in% inputData$entrez & ChemicalID == upregulator)
}

drugSignFunc <- function(overlap, j) {
	return(overlap[j,"edge_type"])
}

upstreamTypefuncMap = list(
	gene = list(
		overlapFunc = geneOverlapFunc,
		signFunc = geneSignFunc
		),
	drug=list(
		overlapFunc = drugOverlapFunc,
		signFunc = drugSignFunc
		)
	)

zscoreUpstream = function(network, upregulators, inputData, overlapFunc, signFunc, tgtGeneIdColumn="target_gene_id", thr=list(edgeConfidence=400)) {
  zscoreUP = c(rep(NA,length(upregulators)))
  for (i in 1:length(upregulators)) {
    
    overlap <- overlapFunc(upregulators[i], network, inputdata, thr)
    #flog.info(paste("p=",nrow(overlap), "p1=",nrow(p1)))
    # stopifnot(nrow(overlap)==nrow(p1))
    numerator <- 0
    denominator <- 0
    if(!is.null(overlap) && nrow(overlap)>0){
      for (j in 1:nrow(overlap)) {
        sR = signFunc(overlap, j)
        if (abs(sR)==1) {
          w = 1#overlap[j,"r.maxScore"]*.001  #we can use the evidence score if needed. decided to ignore it
          numerator <- numerator + w*sR*sign(subset(inputData, entrez==overlap[j,tgtGeneIdColumn])$logfc)
          denominator <- denominator + (w)^2
        }
      }
      denominator <- sqrt(denominator)
      zscoreUP[i] <- numerator/denominator
    }
  }
  
  return(zscoreUP)
}

# read de genes and background measured genes from file
getInputData <- function(identifier, suffix="ga_de.csv"){
  analysisDataFile <- paste(env[["ANALYSIS_TYPE"]], identifier, suffix, sep='_')
  file <- file.path(env[['DATA_FOLDER']], analysisDataFile)
  inputData <- read.csv(file=file, header=TRUE, stringsAsFactors=FALSE)
  return(inputData)
}

# save results to file
postAnalysisData <- function(identifier, daData=data.frame()) {
  analysisDataFile <- paste(env[["ANALYSIS_TYPE"]], identifier, upstream_type, "data.csv", sep='_')
  write.csv(daData, file=file.path(env[['DATA_FOLDER']], analysisDataFile), row.names=FALSE, na ="")
  return()
}

# save parameters to file
saveParameters <- function(identifier, parameters){
  paramsFile <- paste(env[["ANALYSIS_TYPE"]], identifier, "parameters.json", sep='_')
  write(jsonlite::toJSON(parameters, auto_unbox=TRUE), file=file.path(env[['DATA_FOLDER']], paramsFile))
  return()
}


# main logic starts here

env <- Sys.getenv(c("IDENTIFIER","ANALYSIS_TYPE", "CACHE_FOLDER", "DATA_FOLDER", "EKB_VERSION" ))
env

identifier <- env[['IDENTIFIER']]
upstream_type <- 'chemical'#getUpstreamType(identifier, suffix="db_id.csv")
inputData <- getInputData(identifier, suffix="ga_de.csv")
allData <- getInputData(identifier, suffix="ga_all.csv")
params <- list(
    thr=0.05, 
     edgeConfidence=400,
     significantZscore=2,
     corrections=data.frame(do.call(rbind, list(c("No correction", "", ""),
      c("Bonferroni", "FWER", "_bonferroni"),
      c("False discovery rate (FDR)", "FDR", "_fdr")))),
     upTypes=data.frame(do.call(rbind, list(
      c("predicted as activated", "expression patterns in the differentially expressed genes are consistent to the effect observed in the literature", "_p","+"),
      c("predicted as inhibited", "expression patterns in the differentially expressed genes are opposite to the effect observed in the literature", "_n","-")
     ))))

res <- data.frame()
if(upstream_type == "gene") {
  network <- readRDS(file=file.path(env[['CACHE_FOLDER']], env[['EKB_VERSION']], 'network-db', 'networkActions.RDS'))
  measuredTargets <- backgroundOverlap(allData, network, thr=params)
  deTargets <- backgroundOverlap(inputData, network, thr=params)
  res_columns <- c("entrez","symbol","count_de","count_all","pv","pv_fdr","pv_bonferroni","pv_binom","zscore", "pv_zscore", "pv_zscore_fdr", "pv_zscore_bonferroni","count_up","pv_up","pv_up_fdr", "pv_up_bonferroni","count_down","pv_down","pv_down_fdr", "pv_down_bonferroni","pv_binom_up", "pv_binom_down", "pv_comb_up", "pv_comb_up_fdr", "pv_comb_up_bonferroni", "pv_comb_down", "pv_comb_down_fdr", "pv_comb_down_bonferroni")
  
  res <- upstreamAnalysis(inputData, allData, network, measuredTargets, deTargets, thr=params, combFuncMap$fisher, upstreamTypefuncMap$gene$overlapFunc, upstreamTypefuncMap$gene$signFunc, 'target_gene_id')
  res <- merge(allData[,c('entrez','logfc','adjpv')], res[,res_columns], by="entrez", all.y=TRUE)
} else if(upstream_type == "chemical") {
  network <- readRDS(file=file.path(env[["CACHE_FOLDER"]], env[["EKB_VERSION"]], "ctd", "drugNetwork.RDS"))
  
  measuredTargets <- backgroundOverlap(allData, network, thr=params)
  deTargets <- backgroundOverlap(inputData, network, thr=params)
  res <- drugUpstreamAnalysis(inputData, allData, network, measuredTargets, deTargets, thr=params, combFuncMap$fisher, upstreamTypefuncMap$drug$overlapFunc, upstreamTypefuncMap$drug$signFunc, 'GeneID')
  res_columns <- c("ChemicalID","count_de","count_all","pv","pv_fdr","pv_bonferroni","zscore", "pv_zscore", "pv_zscore_fdr", "pv_zscore_bonferroni","count_up","pv_up","pv_up_fdr", "pv_up_bonferroni","count_down","pv_down","pv_down_fdr", "pv_down_bonferroni", "pv_comb_up", "pv_comb_up_fdr", "pv_comb_up_bonferroni", "pv_comb_down", "pv_comb_down_fdr", "pv_comb_down_bonferroni")
  
  params$upTypes=data.frame(do.call(rbind, list(
      c("predicted as present (or overly abundant)", "expression patterns in the differentially expressed genes are consistent to the effect observed in the literature", "_p","+"),
      c("predicted as absent (or insufficient)", "expression patterns in the differentially expressed genes are opposite to the effect observed in the literature", "_n","-")
     )))

  drugs <- network[,c('ChemicalID','id','ChemicalName')] 
  drugs <- distinct(drugs)  
  res <- merge(drugs, res[,res_columns], by="ChemicalID", all.y=TRUE) 
} else {
  flog.info("FAILURE! upstream_type must be either 'gene' or 'chemical'")
  quit()  
}
names(params$corrections) <- c("displayName", "shortName", "suffix")
names(params$upTypes) <- c("displayName", "description", "suffix", "sign")

saveParameters(identifier, parameters=params)
flog.info("results size %s", nrow(res))
if (is.null(res) || nrow(res) == 0) {
  flog.info("empty results. stop.")
} else {
	# Sorting result based on pv_comb_up
  res<-res[order(res$pv_comb_up),]
  postAnalysisData(identifier, res)
  flog.info("done")
}