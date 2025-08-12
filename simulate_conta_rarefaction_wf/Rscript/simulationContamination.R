library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(utils)
library(momr)

args = commandArgs(trailingOnly=TRUE)  
gene_counts.path = args[1]
name.contaminated_sample = args[2]
name.sink = args[3]
name.source = args[4]
contamination_rate = as.numeric(args[5])
sequencing.depth.sink = as.integer(args[6])
sequencing.depth.source = as.integer(args[7])
output = paste0(name.sink, "_", name.source, "_", contamination_rate, "_", sequencing.depth.sink, "_", sequencing.depth.source, ".tsv")


selMGS <- function(genebag, size=50, marker=F, MG=100){
  res <- genebag[as.numeric(summary(genebag)[,"Length"])>= size]
  if(marker ==T){
    f <- function(x){x <- x[1:min(MG,length(x))]}
    res <- lapply(res,f)
  }
  return(res)
}

computeFilteredVectors <- function (profile, type = "mean", filt = 10, debug = FALSE) 
{
  if (class(profile)== "list") {              # corrected elc 07/08/2017 bug Magali
    res <- matrix(data = NA, ncol = ncol(profile[[1]]), nrow = length(profile))
    for (i in 1:length(profile)) {
      if (debug) 
        if (i%%100 == 0) {
          print(i)
        }
      if (class(profile[[i]])[1] == "numeric"){
        res[i, ] <- profile[[i]]              # only one gene
      } else {
        if (type == "mean") {
          res[i, ] <- apply(filterMat(as.matrix(profile[[i]]), 
                                      filt = filt), 2, mean)
        }
        else {
          if (type == "median") {
            res[i, ] <- apply(filterMat(as.matrix(profile[[i]]), 
                                        filt = filt), 2, median)
          }
          else {
            res[i, ] <- apply(filterMat(as.matrix(profile[[i]]), 
                                        filt = filt), 2, sum)
          }
        } 
      }
    }
    rownames(res) <- names(profile)
    colnames(res) <- colnames(profile[[1]])
  }
  else {
    profile <- as.matrix(profile)              
    if (type == "mean") {
      res <- apply(filterMat(as.matrix(profile), filt = filt), 
                   2, mean)
    }
    else {
      if (type == "median") {
        res <- apply(filterMat(as.matrix(profile), filt = filt), 
                     2, median)
      }
      else {
        res <- apply(filterMat(as.matrix(profile), filt = filt), 
                     2, sum)
      }
    }
  }
  return(res)
}

## _________________

hs_10_4_igc2.id = read_tsv(gzcon(url('https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/D07SEU')), col_names = c('gene_id','gene_name','gene_length'))  %>%
mutate(gene_id = as.character(gene_id))

MSP_data = read_tsv(gzcon(url("https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/SOKOXS")))

mgs_gut = split(MSP_data$gene_id,MSP_data$msp_name)
mgs_gut = selMGS(mgs_gut, marker = T)

mgs_gut_tbl = MSP_data %>%
  select(gene_id) %>%
  distinct() %>%
  mutate(
    MG_id = as.integer(gene_id),
    MG_length = hs_10_4_igc2.id$gene_length[MG_id]
  )

gene_counts = read_tsv(gene_counts.path)

sample_source.gene_counts = gene_counts %>%
  select(all_of(name.source)) %>% as.matrix()
rownames(sample_source.gene_counts) = gene_counts$gene_id

sample_sink.gene_counts = gene_counts %>%
  select(all_of(name.sink)) %>% as.matrix()
rownames(sample_sink.gene_counts) = gene_counts$gene_id


sample_sink.available_HQ_reads = sum(sample_sink.gene_counts)
sample_source.available_HQ_reads = sum(sample_source.gene_counts)

final_sample_source.gene_counts = downsizeMatrix(sample_source.gene_counts,
                                                 level = sequencing.depth.source,
                                                 HQ_reads = sample_source.available_HQ_reads,
                                                 silent = F)[,1]

final_contaminated_sample.gene_counts.from_source = downsizeMatrix(sample_source.gene_counts,
                                                             level = sequencing.depth.sink * contamination_rate,
                                                             HQ_reads = sample_source.available_HQ_reads)[,1]
final_contaminated_sample.gene_counts.from_sink = downsizeMatrix(sample_sink.gene_counts,
                                                             level = sequencing.depth.sink * (1-contamination_rate),
                                                             HQ_reads = sample_sink.available_HQ_reads)[,1]

final_contaminated_sample.gene_counts = final_contaminated_sample.gene_counts.from_source + final_contaminated_sample.gene_counts.from_sink
rm(final_contaminated_sample.gene_counts.from_source, final_contaminated_sample.gene_counts.from_sink)

mgs_gut_tbl = mgs_gut_tbl %>%
  mutate(
    source_abundance = replace_na(final_sample_source.gene_counts[gene_id], 0),
    contaminated_sample_abundance = replace_na(final_contaminated_sample.gene_counts[gene_id], 0),
    source_abundance_cov = source_abundance / MG_length * 80,
    contaminated_sample_abundance_cov = contaminated_sample_abundance / MG_length * 80
  )


cur_contamination_case.gene_counts.cov = mgs_gut_tbl %>%
  select(gene_id, source_abundance_cov, contaminated_sample_abundance_cov) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()


prof = extractProfiles(mgs_gut,cur_contamination_case.gene_counts.cov)
cur_contamination_case.mgs_profiles.cov = computeFilteredVectors(prof)
rm(prof)

cur_contamination_case.mgs_profiles.cov.df = cur_contamination_case.mgs_profiles.cov %>%
  as_tibble() %>%
  mutate(id_msp = rownames(.), .before = 1) %>%
  rename(
    !!name.source := 2,
    !!name.contaminated_sample := 3
  )

write.table(cur_contamination_case.mgs_profiles.cov.df, output, sep='\t', row.names=FALSE)
