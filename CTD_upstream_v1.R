require("dplyr")

dataPath = "~/upstream_analysis/data/"

ctdCsvPath = paste0(dataPath, "CTD_chem_gene_ixns_subset.csv")
deCsvPath = paste0(dataPath, "upstream_135_ga_de.csv")

loadNetwork <- function(OrganismID = 9606) {
  graph = read.csv(file = ctdCsvPath)
  graph <- graph[which(graph$OrganismID == OrganismID),]
  graph <- data.frame(lapply(graph, as.character), stringsAsFactors=FALSE)
  edge_type <- rep(NA, nrow(graph))
  edge_type[grep("decreases expression", graph$InteractionActions)] = -1
  edge_type[grep("increases expression", graph$InteractionActions)] = 1
  graph <- cbind(graph, edge_type)
  #graph <- na.omit(graph)
  graph <- graph[,c("ChemicalID", "GeneSymbol", "GeneID", "edge_type")]
  # head(sort(table(graph$ChemicalName), decreasing = TRUE), 20)
  graph <- distinct(graph)
  graph
}

btest <- function (x, p, n){binom.test(x,n,p, alternative = 'greater')$p.value} 
rows <- function(x) lapply(seq_len(nrow(x)), function(i) lapply(x,"[",i))

############## Combining methods #######################################
averageCumulative<-function (x) {
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
      pnorm(mean(x), 1/2, sqrt(1/(12 * n)), lower = TRUE)
    }
  }
}
fisherMethod <- function (x) {
  if (sum(is.na(x)) > 0) 
    NA
  else pchisq(-2 * sum(log(x)), df = 2 * length(x), lower = FALSE)
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
pGMethod <- function (x) { # Combining p-value method used in SPIA paper
  if (sum(is.na(x)) > 0) 
    NA
  else { prod(x) - prod(x)*log(prod(x)) }
}
combFuncMap = list(
  add = function(x, y, z, t){ifelse(sign(z)==sign(t) && abs(z)>=abs(t), addCLT(c(x,y)), y)},
  fisher = function(x, y, z, t){ifelse(sign(z)==sign(t) && abs(z)>=abs(t), fisherMethod(c(x,y)) ,y)},
  # stouffer = function(x, y, z, t){ifelse(sign(z)==sign(t) && abs(z)>=abs(t), stoufferMethod(c(x,y)), y)},
  minP = function(x, y, z, t){ifelse(sign(z)==sign(t) && abs(z)>=abs(t), minP(c(x,y)), y)},
  pGMethod = function(x, y, z, t){ifelse(sign(z)==sign(t) && abs(z)>=abs(t), pGMethod(c(x,y)), y)}
)


getInputData <- function(identifier, cohort, cacheFolder="./data"){
  analysisDataFile <- deCsvPath #paste("upstream", identifier, cohort, "ga_de.csv", sep='_')
  #file = file.workingPath(cacheFolder, analysisDataFile)
  inputData <- read.csv(file=analysisDataFile, header=TRUE,  stringsAsFactors = FALSE)
  return(inputData)
}

backgroundOverlap = function(data, network) {
  p <- subset(network, GeneID %in% data$entrez)
  overlap_data = merge(x = p, y = data, by.x = "GeneID", by.y='entrez', all.x = FALSE)
  overlap_data$sign = with(overlap_data, sign(edge_type)*sign(logfc))
  r = list( count = length(unique(overlap_data$GeneID)),
            count_n = length(unique(subset(overlap_data, sign == -1)$GeneID)), 
            count_p = length(unique(subset(overlap_data, sign == 1)$GeneID)))
  # flog.info('size of background data: count=%s, count_n=%s, count_p=%s', r$count, r$count_n, r$count_p)
  return(r)
}

upstreamAnalysis = function(DEgenes, # DE
                            allData, # all genes measured
                            network, # graph database (e.g. STRING DB)
                            measuredTargets, # genes measured in database
                            deTargets, # DE in database
                            thr = list(edgeConfidence=400, significantZscore=2), # threshold: 400 for edgeConf. and 2 for z-score
                            combFunc = combFuncMap){ # combining p-value functions
  
  if (nrow(DEgenes) == 0) {
    # flog.error('no input data. Size= %s. Stop.', nrow(DEgenes))
    return()
  }
  # flog.info('size of input data: %s', nrow(DEgenes))
  
  deGenes <- DEgenes$entrez
  deTargets = backgroundOverlap(DEgenes, network)
  
  DENet <- subset(network, GeneID %in% deGenes) # list of all edges that have DE genes as targets. 
  if (nrow(DENet) == 0) { # If it is NULL, then no upstream regulator found
    flog.error('no upstream regulators. Size= %s. Stop.', nrow(DENet))
    return()
  }
  ChemDE <- aggregate(DENet$GeneID, by = list(DENet$ChemicalID), FUN=function(DENet){length(unique(DENet))}) # p: number of DE genes for each regulator
  # p sumarizes number of DE under each upregulator
  colnames(ChemDE)<-c('ChemicalID','count_de')
  
  # flog.info("collected upstream regulators")
  ChemNet <- subset(network, ChemicalID %in% ChemDE$ChemicalID) # sub-network of all regulators
  ChemNet <- ChemNet[,c('ChemicalID','GeneID','edge_type')]
  
  # flog.info("collecting background for each upstream regulator")
  ChemNetData = merge(x = ChemNet, y = allData, by.x = "GeneID", by.y='entrez', all.x = FALSE)
  
  ChemAllData = aggregate(x = ChemNetData$GeneID,
                          by = list(ChemicalID=ChemNetData$ChemicalID),
                          FUN = function(x){length(unique(x))})
  colnames(ChemAllData)<-c('ChemicalID','count')
  
  ChemSummary = merge(x = ChemDE, y = ChemAllData, by.x = "ChemicalID", by.y='ChemicalID', all.x = TRUE)
  ChemDE$count_all = ChemSummary$count
  
  ChemNetData$sign <- ChemNetData$edge_type * sign(ChemNetData$logfc) # sign(edge) * sign(fc)
  countP = aggregate( # number of "white balls" in Hypothesis 1 for each regulators
    x = ChemNetData[ChemNetData$sign == 1,'GeneID'], by = list(ChemicalID = ChemNetData[ChemNetData$sign == 1,'ChemicalID']),
    FUN = function(x){length(unique(x))})
  z = merge(x = ChemDE, y = countP, by.x = "ChemicalID", by.y='ChemicalID', all.x = TRUE) 
  ChemDE$count_p = ifelse(is.na(z$x), 0, z$x)
  
  countN = aggregate( # number of "white balls" in Hypothesis 2 for each regulators
    x = ChemNetData[ChemNetData$sign == -1,'GeneID'], by = list(ChemicalID = ChemNetData[ChemNetData$sign == -1, "ChemicalID"]),
    FUN = function(x){length(unique(x))})
  z = merge(x = ChemDE, y = countN, by.x = "ChemicalID", by.y = "ChemicalID", all.x = TRUE) 
  ChemDE$count_n = ifelse(is.na(z$x), 0, z$x)
  
  # flog.info("computing upstream regulators p-values - naive approach")
  ChemDE["pv"] = phyper(ChemDE$count_de-1, deTargets$count, measuredTargets$count-deTargets$count, ChemDE$count_all, lower.tail = FALSE)
  ChemDE["pv_fdr"] = p.adjust(ChemDE$pv, method="fdr")
  ChemDE["pv_bonferroni"] = p.adjust(ChemDE$pv, method="bonferroni")
  ChemDE["pv_binom"] = mapply(btest, ChemDE$count_de, deTargets$count/measuredTargets$count, ChemDE$count_all)
  
  resultNEO4J = ChemDE[,c("ChemicalID", "count_de", "count_all" , "pv", "pv_fdr", "pv_bonferroni", "pv_binom", "count_n", "count_p")]
  # flog.info("computing upstream regulators zscore") # First evidence (zscore and pvalue based on zscore)
  resultNEO4J["zscore"] = zscoreUpstream(resultNEO4J$ChemicalID, DEgenes, ChemNetData)
  resultNEO4J["pv_zscore"]=mapply(function(z){pnorm(-abs(z))}, resultNEO4J["zscore"])
  
  # flog.info("computing upstream regulators upstream p-values") # Second evidence (pvalue upstream)
  r <- pvUpstream(resultNEO4J, DEgenes, network, measuredTargets, deTargets)
  resultNEO4J <- merge(x = resultNEO4J, y = r, by = "ChemicalID", by.y='ChemicalID', all.x = TRUE)
  
  #Combining p-values using Fisher method
  resultNEO4J["pv_Fisher_up"] = mapply(combFunc$fisher, 
                                       resultNEO4J[,"pv_zscore"],
                                       resultNEO4J[,"pv_up"], 
                                       resultNEO4J[,"zscore"], 
                                       rep(thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  resultNEO4J["pv_Fisher_up_fdr"] = p.adjust(resultNEO4J$pv_Fisher_up, method="fdr")
  resultNEO4J["pv_Fisher_up_bonferroni"] = p.adjust(resultNEO4J$pv_Fisher_up, method="bonferroni")
  resultNEO4J["pv_Fisher_down"] = mapply(combFunc$fisher, 
                                         resultNEO4J[,"pv_zscore"],
                                         resultNEO4J[,"pv_down"], 
                                         resultNEO4J[,"zscore"], 
                                         rep(-1*thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  resultNEO4J["pv_Fisher_down_fdr"] = p.adjust(resultNEO4J$pv_Fisher_down, method="fdr")
  resultNEO4J["pv_Fisher_down_bonferroni"] = p.adjust(resultNEO4J$pv_Fisher_down, method="bonferroni")
  
  #Combining p-values using addCLT method
  resultNEO4J["pv_addCLT_up"] = mapply(combFunc$add, 
                                       resultNEO4J[,"pv_zscore"],
                                       resultNEO4J[,"pv_up"], 
                                       resultNEO4J[,"zscore"], 
                                       rep(thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  resultNEO4J["pv_addCLT_up_fdr"] = p.adjust(resultNEO4J$pv_addCLT_up, method="fdr")
  resultNEO4J["pv_addCLT_up_bonferroni"] = p.adjust(resultNEO4J$pv_addCLT_up, method="bonferroni")
  resultNEO4J["pv_addCLT_down"] = mapply(combFunc$add, 
                                         resultNEO4J[,"pv_zscore"],
                                         resultNEO4J[,"pv_down"], 
                                         resultNEO4J[,"zscore"], 
                                         rep(-1*thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  resultNEO4J["pv_addCLT_down_fdr"] = p.adjust(resultNEO4J$pv_addCLT_down, method="fdr")
  resultNEO4J["pv_addCLT_down_bonferroni"] = p.adjust(resultNEO4J$pv_addCLT_down, method="bonferroni")
  
  #Combining p-values using stouffer method
  # resultNEO4J["pv_Stouffer_up"] = mapply(combFunc$stouffer, 
  #                                        resultNEO4J[,"pv_zscore"],
  #                                        resultNEO4J[,"pv_up"], 
  #                                        resultNEO4J[,"zscore"], 
  #                                        rep(thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  # resultNEO4J["pv_Stouffer_up_fdr"] = p.adjust(resultNEO4J$pv_Stouffer_up, method="fdr")
  # resultNEO4J["pv_Stouffer_up_bonferroni"] = p.adjust(resultNEO4J$pv_Stouffer_up, method="bonferroni")
  # resultNEO4J["pv_Stouffer_down"] = mapply(combFunc$stouffer, 
  #                                          resultNEO4J[,"pv_zscore"],
  #                                          resultNEO4J[,"pv_down"], 
  #                                          resultNEO4J[,"zscore"], 
  #                                          rep(-1*thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  # resultNEO4J["pv_Stouffer_down_fdr"] = p.adjust(resultNEO4J$pv_Stouffer_down, method="fdr")
  # resultNEO4J["pv_Stouffer_down_bonferroni"] = p.adjust(resultNEO4J$pv_Stouffer_down, method="bonferroni")
 
  # Combining p-values using minP method
  resultNEO4J["pv_minP_up"] = mapply(combFunc$minP,
                                         resultNEO4J[,"pv_zscore"],
                                         resultNEO4J[,"pv_up"],
                                         resultNEO4J[,"zscore"],
                                         rep(thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  resultNEO4J["pv_minP_up_fdr"] = p.adjust(resultNEO4J$pv_minP_up, method="fdr")
  resultNEO4J["pv_minP_up_bonferroni"] = p.adjust(resultNEO4J$pv_minP_up, method="bonferroni")
  resultNEO4J["pv_minP_down"] = mapply(combFunc$minP,
                                           resultNEO4J[,"pv_zscore"],
                                           resultNEO4J[,"pv_down"],
                                           resultNEO4J[,"zscore"],
                                           rep(-1*thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  resultNEO4J["pv_minP_down_fdr"] = p.adjust(resultNEO4J$pv_minP_down, method="fdr")
  resultNEO4J["pv_minP_down_bonferroni"] = p.adjust(resultNEO4J$pv_minP_down, method="bonferroni")
  
  #Combining p-values using pGMethod method
  resultNEO4J["pv_pGMethod_up"] = mapply(combFunc$pGMethod, 
                                         resultNEO4J[,"pv_zscore"],
                                         resultNEO4J[,"pv_up"], 
                                         resultNEO4J[,"zscore"], 
                                         rep(thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  resultNEO4J["pv_pGMethod_up_fdr"] = p.adjust(resultNEO4J$pv_pGMethod_up, method="fdr")
  resultNEO4J["pv_pGMethod_up_bonferroni"] = p.adjust(resultNEO4J$pv_pGMethod_up, method="bonferroni")
  resultNEO4J["pv_pGMethod_down"] = mapply(combFunc$pGMethod, 
                                           resultNEO4J[,"pv_zscore"],
                                           resultNEO4J[,"pv_down"], 
                                           resultNEO4J[,"zscore"], 
                                           rep(-1*thr$significantZscore, length(resultNEO4J[,"pv_zscore"])))
  resultNEO4J["pv_pGMethod_down_fdr"] = p.adjust(resultNEO4J$pv_pGMethod_down, method="fdr")
  resultNEO4J["pv_pGMethod_down_bonferroni"] = p.adjust(resultNEO4J$pv_pGMethod_down, method="bonferroni")
  
  # Sorting result based on pv_comb_up
  resultNEO4J = resultNEO4J[order(resultNEO4J$pv_Fisher_up),]
  return(resultNEO4J)
}

# ################ INPUT: PREPROCESSED DATA USING IPG ##################################
# # First step:
# # upload Cel files into iPG, set the thresholds: pvalue 0.05, fold-change 1
# # get the csv {gene, p-value, fold-change} of DE genes
# #data = read.csv("GSE11352_allSamples.csv", header = TRUE)
# dataset <- "GSE2639"
# data = read.csv(paste0("dataset/", dataset, "/", dataset,"_DEGenes.csv"), header = TRUE)
# rownames(data) = data$symbol
# data = data[,-1]
# DEGenesSYM = rownames(data)
# foldChange = data$logfc
# names(foldChange) = rownames(data)
# 
# measured_gene = intersect(rownames(data), graph$GeneSymbol) # m + n
# graph <- graph[graph$GeneSymbol %in% measured_gene,] # only consider measured genes

############# Z-Score (IPA) #############

zscoreUpstream = function(upregulators, DEgenes, ChemNetData){
  zscoreUP = c(rep(NA,length(upregulators)))
  names(zscoreUP) = upregulators
  for (i in 1:length(upregulators)) {
    overlap <- subset(ChemNetData, ChemNetData$GeneID %in% DEgenes$entrez & ChemNetData$ChemicalID == upregulators[i])
    numerator <- 0
    for (j in 1:nrow(overlap)) {
      numerator <- numerator + overlap[j, "sign"]
    }
    zscoreUP[i] <- numerator/sqrt(nrow(overlap))
  }
  # zscoreUP <- na.omit(zscoreUP)
  zscoreUP
}



################ p-value of second evidence ############

# pvUpstream = function(upregulators, inputData, network, measuredTargets, deOverlap, thr=list(edgeConfidence=400)){
pvUpstream = function(resultNEO4J, DEgenes, network, measuredTargets, deTargets){ # resultNE04J is the table of regulators' information
  # analysisDatailsFile <- paste(env[["ANALYSIS_TYPE"]], identifier, "details.csv", sep='_')
  for (i in 1:nrow(resultNEO4J)) {
    # print(i)
    overlap <- subset(network,  GeneID %in% DEgenes$entrez & ChemicalID == resultNEO4J[i,]$ChemicalID)
    
    if (!is.null(overlap)) {
      v <- merge(x = overlap, y = DEgenes, by.x = "GeneID", by.y='entrez', all.x = TRUE)
      v['sign'] <- v$edge_type * sign(v$logfc)
      resultNEO4J[i,'count_up'] = length(unique(subset(v, sign == 1)$GeneID))
      resultNEO4J[i,'count_down'] = length(unique(subset(v, sign == -1)$GeneID))
      resultNEO4J[i,'pv_up'] = phyper(resultNEO4J[i,'count_up'] -1, deTargets$count_p, measuredTargets$count-deTargets$count_p, resultNEO4J[i,]$count_all, lower.tail = FALSE)
      #upregulators[i,'pv_up'] = phyper(upregulators[i,'count_up'] -1, upregulators[i,]$count_p, upregulators[i,]$count_all-upregulators[i,]$count_p, upregulators[i,]$count_de, lower.tail = FALSE)
      resultNEO4J[i,'pv_down'] = phyper(resultNEO4J[i,'count_down'] -1, deTargets$count_n, measuredTargets$count-deTargets$count_n, resultNEO4J[i,]$count_all, lower.tail = FALSE)
      #upregulators[i,'pv_down'] = phyper(upregulators[i,'count_down'] -1, upregulators[i,]$count_n, upregulators[i,]$count_all-upregulators[i,]$count_n, upregulators[i,]$count_de, lower.tail = FALSE)
      resultNEO4J[i, 'pv_binom_up'] = binom.test(x=resultNEO4J[i,'count_up'], n=resultNEO4J[i,]$count_all, p=deTargets$count_p/measuredTargets$count, alternative = 'greater')$p.value
      resultNEO4J[i, 'pv_binom_down'] = binom.test(x=resultNEO4J[i,'count_down'], n=resultNEO4J[i,]$count_all, p=deTargets$count_n/measuredTargets$count, alternative = 'greater')$p.value
      
      # write.table(v, file=file.workingPath("./data", analysisDatailsFile), row.names=FALSE, na ="", sep = ",", col.names = (i == 1), append = (i>1))
    }
  }
  
  return(resultNEO4J[,c('ChemicalID', 'count_up','pv_up','count_down', 'pv_down', 'pv_binom_up', 'pv_binom_down')])
}

############ MAIN LOGIC STARTS HERE ################
network <- loadNetwork()
# network <- readRDS(paste(dataPath, "CTD_human.rds"))
identifier = "GSE11352" # https://ipathwayguide.advaitabio.com/intake/54168/contrasts
group <- "48h"
DEgenes <- getInputData(identifier, group, workingPath) #ga_de (list of DE genes, their logfc, and adj.pvalue)
allData <- getInputData(identifier, group, workingPath)

params <- list(thr = 0.05, 
               significantZscore=2,
               corrections = data.frame(do.call(rbind, list(c("No correction", "", ""),
                                                            c("Bonferroni", "FWER", "_bonferroni"),
                                                            c("False discovery rate (FDR)", "FDR", "_fdr")))),
               upTypes = data.frame(do.call(rbind, list(
                 c("predicted as activated", "expression patterns in the differentially expressed genes are consistent to the effect observed in the literature", "_p","+"),
                 c("predicted as inhibited", "expression patterns in the differentially expressed genes are opposite to the effect observed in the literature", "_n","-")
               ))))
names(params$corrections) <- c("displayName", "shortName", "suffix")
names(params$upTypes) <- c("displayName", "description", "suffix", "sign")

measuredTargets <- backgroundOverlap(allData, network)
deTargets <- backgroundOverlap(DEgenes, network)

res <- upstreamAnalysis(DEgenes, allData, network, measuredTargets, deTargets, thr=params)

write.csv(res, file = paste0(dataPath, "upstream_", identifier, "_", group, "_data_v1.csv"), row.names=FALSE, na ="")

