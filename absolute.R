args <- commandArgs(trailingOnly = TRUE)
  library(ABSOLUTE)
  genome <- "hg19"
  platform <- "Illumina_WES"
  meta1.disease <- "Lung"
  sigma.p <- 0
  max.sigma.h <- 0.02
  min.ploidy <- 0.95
  max.ploidy <- 6.05
  max.as.seg.count <- 1000
  max.non.clonal <- 0.2
  max.neg.genome <- 0
  copy_num_type <- "total"
  maf.fn <- paste(as.character(args[3]),"/",as.character(args[1]),sep="")
  seg.dat.fn <- paste(as.character(args[4]),"/",as.character(args[2]),sep="")
  sample.name <- strsplit(as.character(args[1]),"-VS")[[1]][1]
                          

  results.dir <- paste(as.character(args[5]),"/",sample.name,"/",sep="")
  
  RunAbsolute(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, meta1.disease, 
              platform, sample.name, results.dir, max.as.seg.count, max.non.clonal, 
              max.neg.genome, copy_num_type, maf.fn, min.mut.af=0.01, verbose=TRUE)
  
  absolute.files <- paste(results.dir,sample.name,".ABSOLUTE.RData",sep="")

  CreateReviewObject(sample.name, absolute.files, results.dir, "total", verbose=TRUE)
  
  calls.path = paste(results.dir,sample.name,".PP-calls_tab.txt",sep="") 
  modes.path = paste(results.dir,sample.name,".PP-modes.data.RData",sep="") 
  output.path = results.dir
  ExtractReviewedResults(calls.path, sample.name, modes.path, output.path, "absolute", "total")
