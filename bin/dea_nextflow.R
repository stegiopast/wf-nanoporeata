# Library dependencies
library(data.table)
library(BiocParallel)
library(DESeq2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
input_metadata <- args[1]
input_gene_or_trans <- args[2] # = "gene" or "transcript"
input_counttable <- args[3]
input_threads <- args[4]
output_file <- args[5] # = "gene.RData" or "transcript.RData"

# read metadata file (fread)
# input_metadata = path/to/metadata_file, piped in by nextflow pipeline
metadata <- fread(input_metadata)

if(input_gene_or_trans == "gene"){
  # fread count table: input_counttable = path/to/counttable, piped in by nextflow pipeline
  counts <- fread(input_counttable, data.table = F)
  geneNames = counts$Geneid
  counts$Geneid <- NULL
  counts_numeric = sapply(counts, as.numeric)
  row.names(counts_numeric) = geneNames
  nulls = sapply(0,function(counts)rowSums(counts_numeric==counts))
  num_nulls = length(which(nulls[,1] == 0))
  
  # make counts global variable
  counts <- as.data.frame(counts_numeric)
  
} else if (input_gene_or_trans == "transcript"){
  counts = fread(input_counttable, data.table = T)
  geneNames = counts$Name
  counts = counts %>% dplyr::select(-1:-4)
  counts = as.data.frame(counts)
  counts_numeric = sapply(counts, as.numeric)
  row.names(counts_numeric) = geneNames
  nulls = sapply(0,function(counttable)rowSums(counts_numeric==counts))
  
  counts <- as.data.frame(counts_numeric)
}

# threads as input parameter!
register(MulticoreParam(input_threads))

createDDS2 <- function(counts, metadata){
  #flog.info("########## Create DDS object ###########")
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~ Conditions)
  
  dds <- DESeq(dds, parallel = T)
  
  return(dds)
}

run_preprocessing_dea <- function(metadata, counts, output_file){
  #flog.info("########## Differential Expression Analysis ###########")
  
  metadata <- as.data.frame(metadata)
  
  # Rename metadata columns
  colnames(metadata)[c(1,2)] <- c('Samples', 'Conditions') 
  
  row.names(metadata) <- metadata$Samples
  
  missingSampleInfos = colnames(counts)[-which(colnames(counts) %in% metadata$Samples)]
  if (length(missingSampleInfos) > 0){
    print(paste0("No metadata found for the following samples: ", paste(missingSampleInfos, collapse = ",")))
    print(paste0(">>>>> Counts will be excluded!"))
    
  } else {
    print("All required information is included!")
  }
  metadata = metadata[which(row.names(metadata) %in% colnames(counts)),]
  
  metadata = metadata[which(row.names(metadata) %in% intersect(row.names(metadata), colnames(counts))),]
  
  counts = counts[, match(row.names(metadata), colnames(counts))]
  
  dds = createDDS2(counts, metadata)
  rld = rlog(dds)
  
  # save objects 
  save(dds, rld, file=output_file)
  
}


run_preprocessing_dea(metadata, counts, output_file)