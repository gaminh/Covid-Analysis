require(dplyr)

cacheDrugNetwork <- function(file="CTD_chem_gene_ixns.csv", drug_info="chemicals_with_internal_ids.csv") {
    print("Processing CSV")
    graph<-read.csv(file=file, stringsAsFactors=FALSE)
    
    edge_type <- rep(NA, nrow(graph))
    edge_type[graph$InteractionActions == "decreases^expression"] <- -1
    edge_type[graph$InteractionActions == "increases^expression"] <- 1
    graph <- cbind(graph, edge_type)
    graph <- na.omit(graph)
    graph <- graph[,c("ChemicalName", "ChemicalID", "GeneSymbol", "GeneID", "edge_type")]
    graph <- distinct(graph)
    drugs<-read.csv(file=drug_info, stringsAsFactors=FALSE)
    colnames(drugs)<-c('id', "name", "ChemicalID")
    res <- merge(drugs, graph, by="ChemicalID", all.y=TRUE)
    print("Saving")
    print(colnames(res))
    saveRDS(res[,c("id","ChemicalName","ChemicalID", "GeneSymbol", "GeneID", "edge_type")], file="drugNetwork.RDS")
}

cacheDrugNetwork(file="~/refdata/cache/v1906/network-db/CTD_chem_gene_ixns_parsed.csv", drug_info="~/refdata/cache/v1906/ctd/chemicals_with_internal_ids.csv")  
