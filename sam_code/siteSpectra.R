# functions to calculate the number of a,c,g,t and ambiguities
# these are in coronavirus_2020_gisaid_R_pipeline.R but extracted here for other use
# S J Lycett
# 29 May 2020

numAmbs <- function(seq=seq, nucls=c("a","c","g","t","-","n","x")) {
  cseq <- as.character(seq)
  minds<- match(cseq, nucls)
  kk   <- which(!is.finite(minds))
  cseq[kk] <- "x"
  tbl <- table(factor(cseq,levels=nucls))
  return(tbl)
}

siteSpectra <- function(seqs=seqs, nucls=c("a","c","g","t","-","n","x")) {
  site_spec <- matrix(0, length(seqs[1,]), length(nucls))
  colnames(site_spec) <- nucls
  for (j in 1:length(seqs[1,])) {
    site_spec[j,] <- numAmbs(seqs[,j])
  }
  numValid      <- apply(as.matrix(site_spec[,1:4]), 1, sum)
  majority      <- apply(site_spec, 1, which.max)
  majorityVal   <- apply(site_spec, 1, max)
  fractMajority <- majorityVal/numValid
  return( list(site_spec=site_spec, numValid=numValid, 
               majority=majority, majorityVal=majorityVal,
               fractMajority=fractMajority, nucls=nucls))
}
