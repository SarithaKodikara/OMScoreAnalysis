
library(dplyr)

require_files <- function(paths) {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop(
      "Missing required input files:\n",
      paste0("  - ", missing_paths, collapse = "\n"),
      call. = FALSE
    )
  }
}

require_files(c(
  "data/OTU_table.csv",
  "data/taxonomy_table.csv",
  "data_filtered/meta_filtered.csv",
  "data_filtered/OTU_filtered_clr.csv"
))

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

