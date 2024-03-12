library("data.table")
library("dplyr")
library("DESeq2")
library("futile.logger")
library("BiocParallel")

args = commandArgs(trailingOnly=TRUE)
input_metadata <- args[1]
input_gene_or_trans <- args[2]
input_counttable <- args[3]
input_threads <- args[4]
output_file <- args[5]

# read metadata file (fread)
# input_metadata = path/to/metadata_file, piped in by nextflow pipeline
metadata <- fread(input_metadata)


if(input_gene_or_trans == "gene"){
  # fread count table: input_counttable = path/to/counttable, piped in by nextflow pipeline
  counttable <- fread(input_counttable, data.table = F)
  geneNames = counts$Geneid
  counttable$Geneid <- NULL
  counts_numeric = sapply(counttable, as.numeric)
  row.names(counts_numeric) = geneNames
  nulls = sapply(0,function(counts)rowSums(counts_numeric==counttable))
  num_nulls = length(which(nulls[,1] == 0))
  
  # make counts global variable
  counts$df <<- as.data.frame(counts_numeric)
  
  } else if (input_gene_or_trans == "transcript"){
    counttable = fread(input_counttable, data.table = T)
    geneNames = counttable$Name
    counttable = counttable %>% dplyr::select(-1:-4)
    counttable = as.data.frame(counttable)
    counts_numeric = sapply(counttable, as.numeric)
    row.names(counts_numeric) = geneNames
    nulls = sapply(0,function(counttable)rowSums(counts_numeric==counttable))
    
    countsfile$df_trans <<- as.data.frame(counts_numeric)
}

# threads as input parameter!
register(MulticoreParam(input_threads))

createDDS2 <- function(counts, metadata){
  flog.info("########## Create DDS object ###########")
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~ Conditions)
  
  dds <- DESeq(dds, parallel = T)
  
  return(dds)
}

run_preprocessing_dea <- function(metadata, counts, output_file){
  flog.info("########## Differential Expression Analysis ###########")
  
  # Rename metadata columns
  colnames(metadata)[c(1,2)] <- list('Samples', 'Conditions') 
  
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

run_preprocessing_dea(input_metadata, input_counttable, output_file)