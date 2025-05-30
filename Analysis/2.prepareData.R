
OTU <- read.csv("data/OTU_table.csv", row.names = 1)
patientDat <- read.csv("data_filtered/meta_filtered.csv", row.names = 1 )
taxonomy <- read.csv("data/taxonomy_table.csv", row.names = 1 )

# Unify patient names in OTU count matrix and patient metadata
newNames <- sub(".*MB", "MB", colnames(OTU))
newNames <- sub(".*MF", "MF", newNames)
colnames(OTU) <- newNames
commonPatients <- intersect(newNames, rownames(patientDat))

query <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    counts = as.matrix(OTU[,commonPatients])
  ),
  colData = patientDat[commonPatients,]
) 

saveRDS(query, "data_filtered/sce_all.rds")

filteredOTU <- read.csv("data_filtered/OTU_filtered_clr.csv")
selectedOTU <- filteredOTU %>% select(!X) %>% colnames
saveRDS(selectedOTU, "data_filtered/selectedOTUs.rds")


