Sleuth <- function(dplyr, select, run_accession, condition, mutate, biomaRt, useMart, getBM, rename, ensembl_transcript_id, ensembl_gene_id, external_gene_name) {
  base_dir <- "~/Downloads/"

  sample_id <- dir(file.path(base_dir,"results"))

  kal_dirs <- sapply(sample_id, function(id) file.path(base_dir,"results",id,"kallisto"))

  s2c <- read.table(file.path(base_dir, "hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)
  s2c <- dplyr::select(s2c, sample = run_accession, condition)
  s2c <- dplyr::mutate(s2c, path = kal_dirs)
  print(s2c)
  so <- sleuth_prep(s2c, ~ condition)
  so <- sleuth_fit(so)
  so <- sleuth_fit(so, ~1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  models(so)
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "hsapiens_gene_ensembl",
                           host = 'ensembl.org')

  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)

  so <- sleuth_fit(so)

  so <- sleuth_fit(so, ~1, 'reduced')

  so <- sleuth_lrt(so, 'reduced', 'full')

  sleuth_live(so)
  results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

  #gene level
  so <- sleuth_prep(s2c, ~condition, target_mapping = t2g,
                    aggregation_column = 'ens_gene')
}


TxsDE <- function(dplyr, select, run_accession, condition, mutate, biomaRt, useMart, getBM, rename, ensembl_transcript_id, ensembl_gene_id, external_gene_name) {
  base_dir <- "~/Downloads/"

  sample_id <- dir(file.path(base_dir,"results"))

  kal_dirs <- sapply(sample_id, function(id) file.path(base_dir,"results",id,"kallisto"))

  s2c <- read.table(file.path(base_dir, "hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)
  s2c <- dplyr::select(s2c, sample = run_accession, condition)
  s2c <- dplyr::mutate(s2c, path = kal_dirs)
  print(s2c)

  #transcript based
  so <- sleuth_prep(s2c, ~ condition)
  so <- sleuth_fit(so)
  so <- sleuth_fit(so, ~1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  models(so)

  #Including gene names into transcript-level analysis
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "hsapiens_gene_ensembl",
                           host = 'ensembl.org')

  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


  so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)

  so <- sleuth_fit(so)

  so <- sleuth_fit(so, ~1, 'reduced')

  so <- sleuth_lrt(so, 'reduced', 'full')

  sleuth_live(so)
  results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

  #gene level based
  so <- sleuth_prep(s2c, ~condition, target_mapping = t2g,
                    aggregation_column = 'ens_gene')
  so <- sleuth_fit(so)
  so <- sleuth_fit(so, ~1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  models(so)

  #get results
  results_table_gene <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

}
