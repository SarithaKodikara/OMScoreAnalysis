 
library(readr)
library(magrittr)
library(stringr)
library(dplyr)
library(mixOmics)
library(tidyverse)
library(RColorBrewer)

#### Reading Data ####
meta_raw <- read.csv('data/patientdata_3obs_allresponses.csv', row.names=1) 
otu_data<- read.csv('data/OTU_table.csv', row.names=1)%>%t()
rownames(otu_data) <- gsub("Gibby.|Reyes.", "", rownames(otu_data))

########### Filtering meta data ############

# Filter the IDs where redcap_event_name is "follow_up_1" (time of treatment; this is considered at time 0)
ids <- meta_raw %>%
  filter(redcap_event_name == "follow_up_1")%>%
  dplyr::select(id) %>%
  unique()

#140 out of 145 patients have follow up 1

# Filter the metadata to only include rows with IDs in the ids dataframe
# Create a new column 'timepoint' which is the difference in days from the treatment date
meta_filtered <- meta_raw %>%
  rownames_to_column() %>%
  filter(id %in% ids$id) %>%      
  mutate(date = as.Date(date, format = "%m/%d/%y")) %>%  
  group_by(id) %>%
  mutate(followup1_date = date[redcap_event_name == "follow_up_1"],
         timepoint = as.numeric(date - followup1_date)) %>% 
  ungroup()%>% 
  column_to_rownames(var="rowname")

# Sometimes the timepoint is 0 for baseline (baseline and treatment happening on the same day), 
# convert theses into negative number to indicate baseline
meta_filtered$timepoint <- ifelse(meta_filtered$redcap_event_name == "baseline" & meta_filtered$timepoint == 0, 
                                  -1e-100, meta_filtered$timepoint)

# Save the filtered meta data
write.csv(meta_filtered, 'data_filtered/meta_filtered.csv')

########### Filtering microbiome data ############

# Filter the OTU data to only include rows with IDs in the meta_filtered dataframe
OTU_data <- otu_data[rownames(otu_data) %in% rownames(meta_filtered),]

# Remove OTUs with less than 0.01% relative abundance across all samples
OTUtoKeep<-t(apply(OTU_data,1, function(x) x/sum(x))) %>%
  data.frame() %>%
  rownames_to_column('sample_name')%>%
  pivot_longer(cols = -sample_name,
               names_to = 'sOTUs',
               values_to = 'RelativeAbundance') %>%
  mutate(RelativeAbundance = 100*RelativeAbundance) %>%
  left_join(meta_filtered %>%
              rownames_to_column('sample_name'), by = 'sample_name')  %>%
  group_by(sOTUs) %>%
  summarise(MeanRA = mean(RelativeAbundance)) %>%
  filter(MeanRA >= 0.01)%>%
  dplyr::select(sOTUs)%>%
  c()%>%unlist()

# Keep only the OTUs that have a mean relative abundance of at least 0.01%
OTU_filtered_counts<-OTU_data[,OTUtoKeep]

# Convert the OTU counts to clr transformed data
OTU_filtered_clr <- logratio.transfo(as.matrix(OTU_filtered_counts), 
                                     logratio = 'CLR', offset = 1)

# Save the filtered microbiome data
write.csv(OTU_filtered_counts, 'data_filtered/OTU_filtered_counts.csv')
write.csv(OTU_filtered_clr, 'data_filtered/OTU_filtered_clr.csv')
