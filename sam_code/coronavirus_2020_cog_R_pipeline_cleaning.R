# script to process COG download 
# like GISAID pipeline - but modified
# S J Lycett
# 16 May 2020
# 23 May 2020

#####################################################
### coronavirus_2020_cog_R_pipeline_cleaning.R ###
#####################################################

##################################################################################################
# load libraries
library(ape)
library(seqinr)
#library(ips) - would use this if using mafft from within R

#Rpath <- "Rcode//Coronavirus_2020_R//"
source(paste(Rpath,"siteSpectra.R",sep=""))
source(paste(Rpath,"coronavirus_2020_gisaid_R_pipeline_cleaning.R",sep=""))

##################################################################################################
# load utility functions

# Sam HP workstation only - setwd
#setwd("Rcode//Coronavirus_2020_R//")

#Rpath <- ""
#source(paste(Rpath,"getEl.R",sep=""))
#source(paste(Rpath,"get_BEAST_cols.R",sep=""))
#source(paste(Rpath,"calcDecimalDate.R",sep=""))
#source(paste(Rpath,"get_udates.R",sep=""))
#source(paste(Rpath,"coronavirus_2020_gisaid_R_pipeline_cleaning.R",sep=""))

COG_clean_2 <- function( dataDate="2020-05-08",
                            rootPath="E://data//Coronavirus_2020//",
                            dataPath=paste(rootPath,"COG_downloads//",dataDate,"//",sep=""),
                            dataName=paste("cog_",dataDate,"_alignment.fasta",sep=""),
                            echo=TRUE, nthres=300
                            ) {
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"COG_clean_2_log.txt",sep="")
  write("# Log file for COG_clean_2",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  
  write(paste("dataName=",dataName),file=logName,append=TRUE)
  write(paste("# reading",dataName), file=logName,append=TRUE)
  seqs <- read.dna( paste(dataPath,dataName,sep=""), format="fasta", as.matrix=FALSE)
  
  taxa <- as.matrix(attributes(seqs)$names)

  if (echo) {
    print(paste("Number of sequences=",length(taxa)))
  }

  len      <- unlist(lapply(seqs,length))
  lenThres <- 28000
  if (min(len) < lenThres) {
    inds     <- which(len >= lenThres)
    seqs     <- seqs[inds]
   
    write(paste("Length threshold=",lenThres),file=logName,append=TRUE)
    write(paste("Number human sequences > length =",length(inds)), file=logName, append=TRUE)
  
    write(paste("Start file writing=",Sys.time()),file=logName,append=TRUE)
  
    # replace with R object beecause very slow
    #sname    <- paste(dataPath,"human_cov2020_len",lenThres,".fas",sep="")
    #write.dna(seqs, file=sname, format="fasta", nbcol=-1, colsep="")
  
    sname <- paste(dataPath,"cog_cov2020_len",lenThres,".seqs.Rdata",sep="")
    save(seqs, file=sname)
  
    write(paste("Finished file writing=",Sys.time()),file=logName,append=TRUE)
    #write(paste("Sequences file name=human_cov2020_",lenThres,".fas",sep=""),file=logName,append=TRUE)
    write(paste("Sequences file name=",sname,sep=""),file=logName,append=TRUE)
  }
  
  # read-read the sequences (not aligned)
  #seqs     <- read.dna(sname, format="fasta", as.matrix=FALSE)
  taxa     <- attributes(seqs)$names
  if (echo) {
    print(paste("There are",length(taxa),"cog sequences with length >=",lenThres))
  }
  write(paste("# There are",length(taxa),"cog sequences with length >=",lenThres),
        file=logName,append=TRUE)
  
  # quality
  write(paste("# Start quality=",Sys.time()),file=logName,append=TRUE)
  nuclsTbl <- matrix(0, length(taxa), 7)
  colnames(nuclsTbl) <- c("a","c","g","t","-","n","x")
  rownames(nuclsTbl) <- taxa   # added on 29 april 2020 11.30
  for (j in 1:length(taxa)) {
    nuclsTbl[j,] <- numAmbs( unlist(as.character(seqs[j])) )
  }
  write.csv(nuclsTbl,file=paste(dataPath,"cog_cov2020_",lenThres,"_nuclsTbl.csv",sep=""))
  write(paste("nuclsTbl=cog_cov2020_",lenThres,"_nuclsTbl.csv",sep=""),file=logName,append=TRUE)
  
  
  #nthres <- 300 #28000*0.005
  ok_inds <- which(nuclsTbl[,6] <= nthres)
  if (echo) {
    print(paste("of these, there are",length(ok_inds),"sequences with <=",nthres,"Ns"))
  }
  write(paste("nthres=",nthres),file=logName,append=TRUE)
  write(paste("Number of sequences <",nthres,"Ns =",length(ok_inds)), file=logName, append=TRUE)
  write(paste("# There are",length(ok_inds),"sequences with <=",nthres,"Ns"), file=logName, append=TRUE)
  write(paste("# Finished quality=",Sys.time()),file=logName,append=TRUE)
  
  
  seqs <- seqs[ok_inds]
  taxa <- attributes(seqs)$names
  
  write(paste("Start file writing=",Sys.time()),file=logName,append=TRUE)
  write.dna(seqs, file=paste(dataPath,"cog_cov2020_nthres",nthres,".fas",sep=""),
            format="fasta", nbcol=-1, colsep="")

  write(paste("Finished file writing=",Sys.time()),file=logName,append=TRUE)
  write(paste("Sequences file name= cog_cov2020_nthres",nthres,".fas",sep=""),
        file=logName,append=TRUE)
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
}


#use refPath = "refSeq//" or refPath = dataPath
COG_clean_3 <- function(   dataDate="2020-05-08",
                           rootPath="E://data//Coronavirus_2020//",
                           dataPath=paste(rootPath,"COG_downloads//",dataDate,"//",sep=""),
                           refPath ="refSeq//",
                           refName ="COVID-19_human_373_coding_nuclsConsensus.fas",
                           nthres=300,
                           refFragLen=30,
                           simpleCutPos=0,
                           degapFirst=FALSE, goodPart=1000,
                           echo=TRUE
                        ) {
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"COG_clean_3_log.txt",sep="")
  write("# Log file for COG_clean_3",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  write(paste("simpleCutPos=",simpleCutPos),file=logName,append=TRUE)
  write(paste("refFragLen=",refFragLen),file=logName,append=TRUE)
  
  if (echo) {
    print(paste("Reading reference sequence=",refName))
  }
  write(paste("refPath=",refPath),file=logName,append=TRUE)
  write(paste("refName=",refName),file=logName,append=TRUE)
  write("# reading reference sequence", file=logName,append=TRUE)
  
  # read in reference sequence - this is a consensus sequence of good sequences up until 18 Mar 2020
  y              <- read.dna(paste(refPath,refName,sep=""), format="fasta", as.matrix=FALSE)
  refNucls       <- unlist( as.character( y ) )
  refNuclsString <- paste(refNucls,collapse="")
  refStartFrag   <- refNucls[1:refFragLen]
  refStartString <- paste(refStartFrag,collapse="")
  
  # 21 april 2020 - also do 1 shorter because of Australian sequences
  refStartString2<- substring(refStartString, 2)
  
  dataName       <- paste("cog_cov2020_nthres",nthres,".fas",sep="")
  if (echo) {
    print(paste("Reading COG aligned sequences=",dataName))
  }
  write(paste("dataName=",dataName),file=logName,append=TRUE)
  write("# reading COG aligned sequences", file=logName,append=TRUE)
  
  x              <- read.dna(paste(dataPath,dataName,sep=""), format="fasta", as.matrix=FALSE)
  taxa           <- as.matrix(attributes(x)$names)
  nseqs          <- length(x)
  
  write("# trimed sequences to cog_cov2020_start_trim_with_ref.fas", file=logName, append=TRUE)
  write("# list of non-trimmed sequences to trim_start_ORF1ab_fail.txt", file=logName, append=TRUE)
  write("# non-trimmed sequences to cog_cov2020_bad_failed_trim_with_ref.fas", file=logName, append=TRUE)
  write(paste("# Trimming start",Sys.time()), file=logName, append=TRUE)
  
  outName        <- paste(dataPath,"cog_cov2020_start_trim_with_ref.fas",sep="")
  write(paste(">",attributes(y)$names,sep=""), file=outName, append=FALSE)
  write(refNuclsString,file=outName,append=TRUE)
  
  badOutName      <- paste(dataPath,"cog_cov2020_bad_failed_trim_with_ref.fas",sep="")
  write(paste(">",attributes(y)$names,sep=""), file=badOutName, append=FALSE)
  write(refNuclsString,file=badOutName,append=TRUE)
  
  errorLog <- paste(dataPath,"trim_start_ORF1ab_fail.txt")
  write("# sequences which dont have exact 30 nucls start or otherwise failed trimming", file=errorLog, append=FALSE)
  
  if (echo) {
    print("Trimmed sequences to cog_cov2020_start_trim_with_ref.fas")
    print("list of non-trimmed sequences to trim_start_ORF1ab_fail.txt")
  }
  
  if (simpleCutPos>0) {
    if (echo) print(paste("Using simple cut pos=",simpleCutPos,"and not matching with start of reference"))
    write(paste("# Using simple cut pos=",simpleCutPos,"and not matching with start of reference"),
          file=logName, append=TRUE)
  }
  
  badCount <- 0
  goodCount<- 0
  for (j in 1:nseqs) {
    x_nucls <- unlist( as.character(x[j]) )
    if (simpleCutPos>0) {
      pos     <- simpleCutPos
      x_string<- paste(x_nucls,collapse="")
    } else {
      if (degapFirst) {
        if (j==1) {
          write(paste("# degapping first 1-",goodPart,"nucls"), file=logName, append=TRUE)
        }
        part1_string <- gsub("-","",paste(x_nucls[1:goodPart],collapse=""))
        part2_string <- paste(x_nucls[(goodPart+1):length(x_nucls)],collapse="")
        x_string <- paste(part1_string,part2_string,sep="")
      } else {
        x_string<- paste(x_nucls,collapse="")
      }
      pos     <- gregexpr(refStartString,x_string)
    }
    if (pos>0) {
      x_string <- substring(x_string,pos)
      write(paste(">",taxa[j],sep=""), file=outName, append=TRUE)
      write(x_string,file=outName,append=TRUE)
      goodCount <- goodCount+1
    } else {
      pos2     <- gregexpr(refStartString2,x_string)
      if (pos2>0) {
        x_string <- paste("-",substring(x_string,pos),sep="")
        write(paste(">",taxa[j],sep=""), file=outName, append=TRUE)
        write(x_string,file=outName,append=TRUE)
        goodCount <- goodCount+1
      } else {
        if (echo) {
          print(paste("Sorry didnt do",j,taxa[j]))
        }
        write(paste(">",taxa[j],sep=""), file=badOutName, append=TRUE)
        write(x_string,file=badOutName,append=TRUE)
        write(paste(j,taxa[j],sep=","),file=errorLog,append=TRUE)
        badCount <- badCount+1
      }
    }
  }
  
  if (echo) {
    print(paste("Number bad=",badCount))
    print(paste("Number good=",goodCount))
  }
  
  write(paste("Number bad=",badCount), file=logName, append=TRUE)
  write(paste("Number good=",goodCount), file=logName, append=TRUE)
  write("trimName= cog_cov2020_start_trim_with_ref.fas", file=logName, append=TRUE)
  write("trimErrorLog= trim_start_ORF1ab_fail.txt", file=logName, append=TRUE)
  
  write(paste("# Trimming finished",Sys.time()), file=logName, append=TRUE)
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
  
}

COG_align <- function(dataDate="2020-05-08",
                      rootPath="E://data//Coronavirus_2020//",
                      dataPath=paste(rootPath,"COG_downloads//",dataDate,"//",sep=""),
                      refName="consensus",
                      echo=TRUE
                      ) {
  if (echo) {
    print("COG sequences are already aligned, so remove the reference sequence which is first in the file")
  }
  dataName <- "cog_cov2020_start_trim_with_ref.fas"
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"COG_align_log.txt",sep="")
  write("# Log file for COG_align",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  write(paste("dataName=",dataName),file=logName,append=TRUE)
  write("# COG sequences are already aligned, so remove the reference sequence which is first in the file",
        file=logName, append=TRUE)
  
  #sname <- paste(dataPath,dataName,sep="")
  #txtLines <- readLines(sname)
  #txtLines <- txtLines[3:length(txtLines)]
  
  #dataName <- "cog_cov2020_start_trim.fas"
  #write(txtLines, file=paste(dataPath,dataName,sep=""))
  #write(paste("output dataName=",dataName),file=logName,append=TRUE)
  
  if (echo) {
    print("replacing original ref sequence with consensus of the others")
    print("doing this because using already aligned downloaded data")
    print("this data has excess gaps and not all sequences are the same length")
  }
  
  write("# replacing original ref sequence with consensus of the others", file=logName,append=TRUE)
  write("# doing this because using already mafft downloaded data", file=logName,append=TRUE)
  
  seqs <- read.dna(paste(dataPath,dataName,sep=""), as.matrix=FALSE, format="fasta")
  taxa <- attributes(seqs)$names
  ref_i<- grep(refName,taxa)
  inds <- setdiff(1:length(taxa),ref_i)
  seqs <- seqs[inds]
  taxa <- attributes(seqs)$names
  seqlen <- unlist(lapply(seqs, length))
  len  <- max(seqlen)
  
  
  if (echo) {
    print("Finding new consensus")
  }
  write(paste("# start finding new consensus",Sys.time()), file=logName, append=TRUE)
  
  # hideous but necessary because the sequences are not the same length
  seqMatrix <- matrix("-", length(seqs), len)
  for (j in 1:length(seqs)) {
    cseq <- unlist(as.character(seqs[j]))
    seqMatrix[j,1:length(cseq)] <- cseq
  }
  consensusNucls <- apply(seqMatrix, 2, numAmbs)
  kk             <- apply(consensusNucls[1:6,], 2, which.max)
  consensusSeq   <- rownames(consensusNucls)[kk]
  
  
  write(paste("# finished finding new consensus",Sys.time()), file=logName, append=TRUE)
  
  temptaxa <- c(">consensus",paste(">",taxa,sep=""))
  tempseq  <- apply(seqMatrix, 1, paste, collapse="")
  tempseq  <- c( paste(consensusSeq, collapse=""), tempseq   )
  tempM    <- matrix( c(temptaxa,tempseq), length(temptaxa), 2)
  tempM    <- t(tempM)
  tempM    <- matrix(tempM, 1, length(temptaxa)*2)
  write(tempM, file=paste(dataPath,"cog_cov2020_start_trim_with_ref_mafft.fas",sep=""))
  
  write("# new psuedo alignment file= cog_cov2020_start_trim_with_ref_mafft.fas", file=logName, append=TRUE)
  
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
}



##########################################################################################

#maxRefLen = 41672
COG_tidy_1 <- function(dataDate="2020-05-08",
                       rootPath="E://data//Coronavirus_2020//",
                       dataPath=paste(rootPath,"COG_downloads//",dataDate,"//",sep=""),
                          refName="consensus",
                          maxFractGap=0.1,
                          maxDivergence=-1,
                          maxRefLen=-1,calcMaxRefLen=TRUE,
                          echo=TRUE
) {
  
  # old refName="COVID-19_human_373_coding_consensus"
  
  if (maxDivergence<=0) {
    dataMon <- getEl(dataDate,ind=2,sep="-")
    namedMonths <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    #mm <- match(tolower(dataMon), tolower(namedMonths))
    mm <- as.integer(dataMon)
    maxDivergence = 10*mm
  }
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"COG_tidy_1_log.txt",sep="")
  write("# Log file for COG_tidy_1",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  
  dataName <- "cog_cov2020_start_trim_with_ref_mafft.fas"
  fname    <- paste(dataPath,dataName,sep="")
  write(paste("dataName=",fname), file=logName, append=TRUE)
  
  if (echo) {
    print("Part A")
    print(paste("Reading aligned sequences=",fname))
  }
  
  write(paste("# Reading aligned sequences",dataName), file=logName, append=TRUE)
  
  seqs   <- read.dna(fname, format="fasta", as.matrix=TRUE)
  taxa   <- attributes(seqs)$dimnames[[1]]
  ref_i  <- grep(refName,taxa)
  refseq <- unlist(as.character(seqs[ref_i,]))
  
  if (calcMaxRefLen) {
    maxRefLen <- max(which(refseq !="n" & refseq !="-"))
  } else {
    if (maxRefLen>0) {
      if (length(refseq)<maxRefLen) {
        maxRefLen <- length(refseq)
      }
    } else {
      maxRefLen <- length(refseq)
    }
  }
  write(paste("maxRefLen=",maxRefLen),file=logName,append=TRUE)
  write(paste("# trimming sequences to maxRefLen=",maxRefLen),file=logName,append=TRUE)
  if (echo) {
    print(paste("trimming sequences to maxRefLen=",maxRefLen))
  }
  
  seqs <- seqs[,1:maxRefLen]
  refseq <- refseq[1:maxRefLen]
  
  write("# Calculating differences from reference sequence", file=logName, append=TRUE)
  write(paste("maxFractGap=",maxFractGap), file=logName, append=TRUE)
  write(paste("maxDivergence=",maxDivergence), file=logName, append=TRUE)
  if (echo) {
    print("Calculating differences from reference sequence")
  }
  
  site_res <- siteSpectra(seqs)
  ref_ambs <- numAmbs(seqs[ref_i,])
  ref_valid<- sum(ref_ambs[1:4])
  nseqs    <- length(seqs[,1])
  write(paste("ref_valid=",ref_valid), file=logName,append=TRUE)
  write(paste("ref_gaps=",ref_ambs[5]),file=logName,append=TRUE)
  write(paste("ref_ambs=",ref_ambs[6]+ref_ambs[7]),file=logName,append=TRUE)
  write(paste("ref_length_with_gaps=",length(refseq)),file=logName,append=TRUE)
  
  diff_to_ref <- seqs != t(matrix( rep(seqs[ref_i,],nseqs), length(refseq), nseqs))
  ndiff       <- apply(diff_to_ref, 1, sum)
  ndiff_valid <- apply(diff_to_ref[,which((refseq!="-") & (refseq!="n"))], 1, sum)
  
  mostly_gaps <- which(site_res$numValid<=maxFractGap*nseqs)
  write(paste("number of mostly_gaps=",length(mostly_gaps)), file=logName, append=TRUE)
  gap_causing <- c()
  for (i in 1:length(mostly_gaps)) {
    inds <- which(as.character(seqs[,mostly_gaps[i]]) != "-")
    gap_causing <- c(gap_causing,inds)
  }
  top_gap_causing_tbl <- sort(table(gap_causing),decreasing=TRUE)
  seqs_to_remove      <- as.integer(rownames(top_gap_causing_tbl)[which(top_gap_causing_tbl>3)])
  seqs_to_remove      <- intersect(seqs_to_remove, which(ndiff_valid > maxDivergence))
  
  if (echo) {
    print(paste("ref_valid=",ref_valid))
    print(paste("ref_gaps=",ref_ambs[5]))
    print(paste("ref_ambs=",ref_ambs[6]+ref_ambs[7]))
    print(paste("ref_length_with_gaps=",length(refseq)))
    print(paste("number of mostly_gaps=",length(mostly_gaps)))
    print("Remove sequences")
    print(paste(seqs_to_remove,taxa[seqs_to_remove],ndiff[seqs_to_remove],ndiff_valid[seqs_to_remove]))
  }
  errorLog <- paste(dataPath,"gap_causing_fail.txt",sep="")
  write("# sequences which caused exessive gaps and too different to reference",file=errorLog,append=FALSE)
  write(c("Index,Taxa,Ndiff,Ndiff_valid"),file=errorLog,append=TRUE)
  write(paste(seqs_to_remove,taxa[seqs_to_remove],ndiff[seqs_to_remove],ndiff_valid[seqs_to_remove],sep=","), file=errorLog, append=TRUE)
  
  write.dna(seqs[seqs_to_remove,], file=paste(dataPath,"gap_causing_fail.fas",sep=""),
            format="fasta", nbcol=-1, colsep="")
  
  seqs_to_keep <- setdiff(1:nseqs,seqs_to_remove)
  seqs <- seqs[seqs_to_keep,]
  
  write(paste("number sequences removed=",length(seqs_to_remove)), file=logName,append=TRUE)
  write(paste("number sequences keep=",length(seqs_to_keep)),file=logName,append=TRUE)
  
  #dataName <- "cog_cov2020_with_ref_al0.fas"
  dataName <- "cog_cov2020_with_ref_al0.seqs.Rdata"
  fname    <- paste(dataPath,dataName,sep="")
  
  if (echo) {
    print(paste("number sequences removed=",length(seqs_to_remove)))
    print(paste("number sequences keep=",length(seqs_to_keep)))
    print(paste("Writing to",fname))
  }
  
  #write.dna(seqs, file=fname, format="fasta", nbcol=-1, colsep="")
  save(seqs, file=fname)
  write(paste("alignment 0=",dataName), file=logName,append=TRUE)
  
  if (echo) {
    print("Part B - Degapping begin")
  }
  write("# Begin de-gapping", file=logName, append=TRUE)
  
  # give it a few s to make sure that file is written properly
  #Sys.sleep(5)
  
  #seqs     <- read.dna(file=fname,format="fasta",as.matrix=TRUE)
  site_res <- siteSpectra(seqs)
  all_gaps <- which(site_res$numValid <= 1)
  ok_pos   <- which(site_res$numValid > 1)
  seqs <- seqs[,ok_pos]
  
  #dataName <- "cog_cov2020_with_ref_al1.fas"
  dataName <- "cog_cov2020_with_ref_al1.seqs.Rdata"
  fname    <- paste(dataPath,dataName,sep="")
  
  if (echo) {
    print(paste("number gap sites removed=",length(all_gaps)))
    print(paste("number sites keep=",length(ok_pos)))
    print(paste("Writing to",fname))
  }
  
  #write.dna(seqs, file=fname, format="fasta", nbcol=-1, colsep="")
  save(seqs, file=fname)
  write(paste("alignment 1=",dataName), file=logName,append=TRUE)
  
  
  if (echo) {
    print("Part C - trim end and remove reference")
  }
  write("# trim end and remove reference")
  
  # give it a few s to make sure that file is written properly
  #Sys.sleep(5)
  
  #seqs     <- read.dna(file=fname,format="fasta",as.matrix=TRUE)
  taxa     <- attributes(seqs)$dimnames[[1]]
  ref_i    <- grep(refName,taxa)
  refseq   <- unlist(as.character(seqs[ref_i,]))
  
  pos <- gregexpr("tag",paste(refseq,collapse=""))
  endpos <- pos[[1]][length(pos[[1]])]
  
  pos <- gregexpr("atg",paste(refseq,collapse=""))
  startpos <- pos[[1]][1]
  
  seqs     <- seqs[,startpos:(endpos+2)]
  
  trimmed_refseq <- paste(unlist(as.character(seqs[ref_i,])),collapse="")
  write(paste(">",taxa[ref_i],sep=""), file=paste(dataPath,"trimmed_refSeq.fas",sep=""))
  write(trimmed_refseq, file=paste(dataPath,"trimmed_refSeq.fas",sep=""), append=TRUE)
  
  inc_inds <- setdiff(1:length(seqs[,1]),ref_i)
  seqs     <- seqs[inc_inds,]
  
  dataName <- "cog_cov2020_al2.seqs.Rdata"
  fname    <- paste(dataPath,dataName,sep="")
  save(seqs, file=fname)
  
  dataName <- "cog_cov2020_al2.fas"
  fname    <- paste(dataPath,dataName,sep="")
  
  if (echo) {
    print("Trimmed to start and end of ref")
    print("Removed ref")
    print(paste("number of sequences=",length(seqs[,1])))
    print(paste("number of sites=",length(seqs[1,])))
    print(paste("Writing to",fname))
  }
  
  write.dna(seqs, file=fname, format="fasta", nbcol=-1, colsep="")
  write(paste("alignment 2=",dataName), file=logName,append=TRUE)
  write(paste("number sequences=",length(seqs[,1])),file=logName,append=TRUE)
  write(paste("number sites=",length(seqs[1,])),file=logName,append=TRUE)
  
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
  
}

COG_tidy_2 <- function(dataDate="2020-05-08",
                       rootPath="E://data//Coronavirus_2020//",
                       dataPath=paste(rootPath,"COG_downloads//",dataDate,"//",sep=""),
                       refPath ="refSeq//",
                       refName ="COVID-19_human_373_coding_nuclsConsensus.fas",
                       echo=TRUE) {
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"COG_tidy_2_log.txt",sep="")
  write("# Log file for COG_tidy_2",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  
  dataName <- "cog_cov2020_al2.fas"
  fname    <- paste(dataPath,dataName,sep="")
  write(paste("dataName=",fname), file=logName, append=TRUE)
  write(paste("refPath=",refPath),file=logName, append=TRUE)
  write(paste("refName=",refName),file=logName, append=TRUE)
  
  if (echo) {
    print(paste("Reading sequences",fname))
  }
  
  seqs    <- read.dna(fname, format="fasta", as.matrix=TRUE, as.character=TRUE)
  taxa    <- as.matrix(attributes(seqs)$dimnames[[1]])
  siteSpec<- siteSpectra(seqs)
  consensusSeq <- siteSpec$nucls[siteSpec$majority]
  nseqs        <- length(seqs[,1])
  refNucls     <- read.dna(paste(refPath,refName,sep=""), format="fasta", as.matrix=TRUE, as.character=TRUE)
  added_gaps   <- setdiff(which(refNucls=="-"),which(consensusSeq=="-"))
  kk           <- added_gaps[1]
 
  if (echo) {
    print("Adding gaps for inframe coding")
  }
  write("# Adding gaps for inframe coding",file=logName,append=TRUE)
  while(length(kk)>0) {
    if (echo) {
      print(paste("Added gap",kk))
    }
    write(paste("added gap",kk),file=logName,append=TRUE)
    nlen    <- length(consensusSeq)
    newSeqs                   <- matrix("-", nseqs, nlen+1)
    newSeqs[,1:(kk-1)]        <- seqs[,1:(kk-1)]
    newSeqs[,(kk+1):(nlen+1)] <- seqs[,kk:nlen]
    seqs         <- newSeqs
    consensusSeq <- c(consensusSeq[1:(kk-1)],"-",consensusSeq[kk:nlen])
    added_gaps   <- setdiff(which(refNucls=="-"),which(consensusSeq=="-"))
    if (length(added_gaps)>0) {
      kk <- added_gaps[1]
    } else {
      kk <- c()
    }
  }
   
  consensusSeq <- consensusSeq[1:length(refNucls)]
  non_matching <- (which(refNucls!=consensusSeq))
  seqs         <- seqs[,1:length(refNucls)]
  
  write("#Positions where consensus not match reference", file=logName, append=TRUE)
  write(paste("different consensus to ref =",non_matching), file=logName, append=TRUE)
  
  seqsTxt      <- apply(seqs, 1, paste, collapse="")
  temp         <- matrix( c(paste(">",taxa,sep=""),seqsTxt), length(seqsTxt), 2)
  temp         <- matrix( t(temp), 2*length(seqsTxt), 1)
  
  
  dataName <- "cog_cov2020_al3.fas"
  if (echo) {
    print(paste("Writing to file",dataName))
  }
  fname    <- paste(dataPath,dataName,sep="")
  write(temp, file=fname)
  write(paste("alignment 3=",dataName),file=logName,append=TRUE)
  
  write(">coding_consensus",file=paste(dataPath,"cog_cov2020_coding_consensus.fas",sep=""))
  write(paste(consensusSeq,collapse=""),file=paste(dataPath,"cog_cov2020_coding_consensus.fas",sep=""),append=TRUE)
  
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
  
}

# formerly COG_tidy_2
COG_tree <- function(dataDate="2020-05-08",
                       rootPath="E://data//Coronavirus_2020//",
                       dataPath=paste(rootPath,"COG_downloads//",dataDate,"//",sep=""),
                       dataName="cog_cov2020_al2.fas",
                          echo=TRUE
            ) {
  
  
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"COG_tree_log.txt",sep="")
  write("# Log file for COG_tree",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  write(paste("dataName=",dataName), file=logName, append=TRUE)
  fname    <- paste(dataPath,dataName,sep="")
  
  if (echo) {
    print(paste("Reading sequences",fname))
  }
  
  seqs <- read.dna(fname, format="fasta", as.matrix=TRUE)
  
  write("# making initial NJ tree with TN93, gamma=0.5, pairwise.deletion=TRUE", file=logName, append=TRUE)
  write(paste("# tree start",Sys.time()), file=logName, append=TRUE)
  
  if (echo) {
    print("Making initial NJ tree with TN93, gamma=0.5, pairwise.deletion=TRUE")
  }
  
  dd <- dist.dna(seqs, model="TN93", gamma=0.5, pairwise.deletion=TRUE, as.matrix=TRUE)
  save(dd, file=paste(dataPath,dataName,"_genetic_dist.dd.Rdata",sep=""))
  tr <- nj(dd)
  tr <- ladderize(tr)
  save(tr, file=paste(dataPath,dataName,"_tree.tr.Rdata",sep=""))
  
  trName <- paste(dataPath,dataName,"_ape_TN93_nj.nwk",sep="")
  write.tree(tr, file=trName)
  
  write(paste("# tree end",Sys.time()), file=logName, append=TRUE)
  write(paste("treeName=",trName), file=logName, append=TRUE)
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
  
}


#not used - see elsewhere now

#get_COG_metadata <- function(dataDate="2020-05-08",
#                             rootPath="E://data//Coronavirus_2020//",
#                             dataPath=paste(rootPath,"COG_downloads//",dataDate,"//",sep="")) {
#  
#  fname    <- paste(dataPath,"cog_",dataDate,"_metadata.txt",sep="")
#  metadata <- read.csv(fname)
#  decDate  <- apply(as.matrix(metadata$sample_date), 1, calcDecimalDate_fromTxt, sep="-", dayFirst=FALSE)
#  country_state <- as.matrix( metadata$adm1 )
#  uadm1 <- c("UK-ENG","UK-WLS","UK-SCT","UK-NIR")
#  ustate<- c("England","Wales","Scotland","Northern Ireland")
#  for (i in 1:4) {
#    country_state <- gsub(uadm1[i],ustate[i],country_state)
#  }
#  country_state <- factor(country_state, levels=ustate)
#  metadata      <- cbind(metadata,country_state,decDate)
#  
#  return( metadata )
#}

COG_lineage_clusters <- function(dataDate="2020-05-08",
                                 rootPath="E://data//Coronavirus_2020//",
                                 dataPath=paste(rootPath,"COG_downloads//",dataDate,"//",sep="")) {
  
  ################################################
  metadata <- get_COG_metadata(dataDate=dataDate,dataPath=dataPath)
  
  max_ii    <- which.max(metadata$decDate)
  min_ii    <- which.min(metadata$decDate)
  min_date  <- as.integer(strsplit(paste(metadata$sample_date[min_ii]),"-")[[1]])
  max_date  <- as.integer(strsplit(paste(metadata$sample_date[max_ii]),"-")[[1]])
  
  if (min(decDate) > 2020) {
    min_date<- c(2020,1,1)
  }
  
  dates_res <- get_udates(min_date=min_date,max_date=max_date)
  
  # all the UK sequnences not just the good ones
  u_uk_states <- c("England","Wales","Scotland","Northern Ireland")
  uk_cols   <- get_BEAST_cols(length(u_uk_states), bright=0.8, sat=0.7)
  
  uk_dates  <- factor( paste(metadata$sample_date), levels=dates_res$udate)
  uk_states <- factor( paste(metadata$country_state), levels=c(u_uk_states) )
  
  imageName <- paste(dataPath,"cog_all_UK_summary_barplot_6x12.png",sep="")
  png(file=imageName, height=6*300, width=12*300, res=300)
  uk_date_tbl <- table(uk_dates,uk_states)
  barplot(t(uk_date_tbl),col=uk_cols,names=dates_res$xlabs,xlab="Collection date",ylab="Number of sequences")
  legend("topleft",paste(colnames(uk_date_tbl)," (",apply(uk_date_tbl, 2, sum),")",sep=""),
         pch=22,pt.bg=uk_cols,bty="n")
  title(paste("Number of UK Sequences",toupper(gsub("_"," ",dataDate))))
  dev.off()
  write.csv(uk_date_tbl, file=paste(dataPath,"cog_all_UK_summary_date_table.csv",sep=""))
  ###############################################################
  
  #for (j in 1:length(u_uk_states)) {
  #  inds <- which(metadata$country_state==u_uk_states[j])
  #  ltbl <- table(metadata$lineage[inds], metadata$epi_week[inds])
  #}
  
}

# mostly finished 29 May 2020
COG_clusters <- function(dataDate="2020-05-08",
                         rootPath="E://data//Coronavirus_2020//",
                         dataPath=paste(rootPath,"COG_downloads//",dataDate,"//",sep=""),
                         dataName="cog_cov2020_al2.fas",
                         metaPath=paste(rootPath,"COG_downloads//",dataDate,"//microreact//",sep=""),
                         echo=TRUE) {
  
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"COG_clusters_log.txt",sep="")
  write("# Log file for COG_clusters",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  
  fname    <- paste(dataPath,dataName,sep="")
  write(paste("dataName=",fname), file=logName, append=TRUE)
  
  if (echo) {
    print(paste("Reading sequences",fname))
  }
  
  seqs      <- read.dna(fname, format="fasta", as.matrix=FALSE)
  seqlen    <- length(seqs[[1]])
  numSeqs   <- length(seqs)
  taxa      <- attributes(seqs)$names
  #rm(seqs)
  write(paste("seqlen=",seqlen),file=logName, append=TRUE)
  write(paste("numSeqs=",numSeqs),file=logName, append=TRUE)
  if (echo) print(paste("seqlen=",seqlen))
  if (echo) print(paste("numSeqs=",numSeqs))
  
  load( paste(metaPath,"cog_global_",dataDate,"_meta_extended.metadata.Rdata",sep="") )
  maxEpiWeek     <- max(metadata$epi_week[which(is.finite(metadata$epi_week))])
  epiweek        <- 1:maxEpiWeek
  epidate        <- array(0,maxEpiWeek)
  for (i in 1:maxEpiWeek) {
    ii <- which(metadata$epi_week==i)
    jj <- ii[ which.min(metadata$decDates[ii])]
    epidate[i] <- paste(metadata$sample_date[jj])
  }
  epidate <- gsub("2020-","",epidate)
  epidate <- gsub("2019-","",epidate)
  epidate[1] <- "12-29"
  epidate[2] <- "01-05"
  epidate[3] <- "01-12"
  epidate[21]<- "05-17"
  
  load(paste(dataPath,dataName,"_genetic_dist.dd.Rdata",sep=""))
  load(paste(dataPath,dataName,"_tree.tr.Rdata",sep=""))
  ntips <- length(tr$tip.label)
  
  write(paste("num rows dd=",length(dd[,1])),file=logName, append=TRUE)
  write(paste("num cols dd=",length(dd[1,])),file=logName, append=TRUE)
  write(paste("num tips tr=",ntips), file=logName, append=TRUE)
  if (echo) paste("num rows dd=",length(dd[,1]))
  if (echo) paste("num cols dd=",length(dd[1,]))
  if (echo) paste("num tips tr=",ntips)
  
  
  ucont <- sort(unique(metadata$continent))
  u_uk_state <- c("England","Wales","Scotland","Northern Ireland")
  u_scot<- sort(unique(metadata$uk_county[which(metadata$country_state=="Scotland")]))
  u_scot<- setdiff(u_scot,"UNKNOWN")
  
  minds <- match(rownames(dd), gsub(" ","_",metadata$newTaxa))
  all(rownames(dd)==gsub(" ","_",metadata$newTaxa)[minds])
  metadata <- metadata[minds,]
  #save(metadata, file=paste(dataPath,dataName,"_meta_extended.metadata.Rdata",sep=""))
  
  # quality
  #write(paste("# Start quality=",Sys.time()),file=logName,append=TRUE)
  #nuclsTbl <- matrix(0, numSeqs, 7)
  #colnames(nuclsTbl) <- c("a","c","g","t","-","n","x")
  #rownames(nuclsTbl) <- taxa   # added on 29 april 2020 11.30
  #for (j in 1:length(taxa)) {
  #  nuclsTbl[j,] <- numAmbs( unlist(as.character(seqs[j])) )
  #}
  #write.csv(nuclsTbl,file=paste(dataPath,dataName,"_nuclsTbl.csv",sep=""))
  
  #nuclsTbl <- read.csv(paste(dataPath,"cog_cov2020_28000_nuclsTbl.csv",sep=""))
  #rownames(nuclsTbl) <- nuclsTbl[,1]
  #nuclsTbl <- nuclsTbl[,2:length(nuclsTbl[1,])]
  #dinds <- match(rownames(dd), rownames(nuclsTbl))
  #nuclsTbl <- nuclsTbl[dinds,]
  
  #ok_inds <- which((nuclsTbl[,5] <= 100) & nuclsTbl[,6] <= 100)
  ok_inds <- 1:numSeqs
  gtol    <- 0.1/seqlen
  GG      <- 1*(dd <= gtol)
  M       <- GG
  M       <- M[ok_inds,ok_inds]
  g       <- graph.adjacency(M)
  clusts  <- clusters(g)
  
  numclusts    <- clusts$no
  maxclust     <- which.max(clusts$csize)
  maxinds      <- which(clusts$membership==maxclust)
  maxClustSize <- clusts$csize[maxclust]
  if (maxClustSize>1000) {
    maxval <- maxClustSize
  } else {
    maxval <- 1000
  }
  breaks <- c(0,1,10,20,50,100,maxval)
  
  h<-hist(clusts$csize, breaks=breaks, plot=FALSE)
  clust_sizes_tbl <- t(as.matrix(h$counts))
  colnames(clust_sizes_tbl) <- paste(breaks[1:(length(breaks)-1)]+1,breaks[2:length(breaks)],sep="->")
  colnames(clust_sizes_tbl)[1] <- "1"
  
  write(paste("Num clusters=",clusts$no), file=logName, append=TRUE)
  write("# Cluster sizes table", file=logName, append=TRUE)
  write.table(clust_sizes_tbl, file=logName, append=TRUE, sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
  
  all_cluster_matrix  <- matrix(0, clusts$no, maxEpiWeek)
  uk_cluster_matrix   <- matrix(0, clusts$no, maxEpiWeek)
  scot_cluster_matrix <- matrix(0, clusts$no, maxEpiWeek)
  england_cluster_matrix <- matrix(0, clusts$no, maxEpiWeek)
  wales_cluster_matrix <- matrix(0, clusts$no, maxEpiWeek)
  northernireland_cluster_matrix <- matrix(0, clusts$no, maxEpiWeek)
  for (i in 1:clusts$no) {
    cinds <- which(clusts$membership==i)
    all_cluster_matrix[i,] <- table(factor(metadata$epi_week[cinds], levels=epiweek))

    jj    <- metadata$country[cinds]=="UK"
    uk_cluster_matrix[i,] <- table(factor(metadata$epi_week[cinds[jj]], levels=epiweek))
        
    jj    <- metadata$country_state[cinds]=="Scotland"
    scot_cluster_matrix[i,] <- table(factor(metadata$epi_week[cinds[jj]], levels=epiweek))
    
    jj    <- metadata$country_state[cinds]=="England"
    england_cluster_matrix[i,] <- table(factor(metadata$epi_week[cinds[jj]], levels=epiweek))
    
    jj    <- metadata$country_state[cinds]=="Wales"
    wales_cluster_matrix[i,] <- table(factor(metadata$epi_week[cinds[jj]], levels=epiweek))
    
    jj    <- metadata$country_state[cinds]=="Northern Ireland"
    northernireland_cluster_matrix[i,] <- table(factor(metadata$epi_week[cinds[jj]], levels=epiweek))
  }
  
  cluster_max    <- matrix(0, 6, maxEpiWeek)
  clustered_seqs <- matrix(0, 6, maxEpiWeek)
  rownames(cluster_max)    <- c("All","UK","England","Wales","Scotland","Northern Ireland")
  rownames(clustered_seqs) <- c("All","UK","England","Wales","Scotland","Northern Ireland")
  for (i in 1:6) {
    if (i==1) cluster_matrix <- all_cluster_matrix
    if (i==2) cluster_matrix <- uk_cluster_matrix
    if (i==3) cluster_matrix <- england_cluster_matrix
    if (i==4) cluster_matrix <- wales_cluster_matrix
    if (i==5) cluster_matrix <- scot_cluster_matrix
    if (i==6) cluster_matrix <- northernireland_cluster_matrix
  
    cluster_max[i,]    <- apply(cluster_matrix, 2, max)
    clustered_seqs[i,] <- apply(cluster_matrix, 2, sum)
    
    # counting clusters (number of clusters) of different size ranges per epi week
    num_singles    <- apply(cluster_matrix==1, 2, sum)
    num_2_10       <- apply( (cluster_matrix>=2) & (cluster_matrix<=10), 2, sum   )
    num_11_20      <- apply( (cluster_matrix>=11) & (cluster_matrix<=20), 2, sum   )
    num_21_50      <- apply( (cluster_matrix>=21) & (cluster_matrix<=50), 2, sum   )
    num_51_100      <- apply( (cluster_matrix>=51) & (cluster_matrix<=100), 2, sum   )
    num_over_100    <- apply( (cluster_matrix>=101) , 2, sum   )
    cluster_nums  <- t(cbind(num_singles, num_2_10, num_11_20, num_21_50, num_51_100, num_over_100))
  
    # counting number of sequences within clusters of different size ranges per epi week
    val_singles    <- apply( cluster_matrix*(cluster_matrix==1), 2, sum)
    val_2_10       <- apply( cluster_matrix*((cluster_matrix>=2) & (cluster_matrix<=10)), 2, sum   )
    val_11_20      <- apply( cluster_matrix*((cluster_matrix>=11) & (cluster_matrix<=20)), 2, sum   )
    val_21_50      <- apply( cluster_matrix*((cluster_matrix>=21) & (cluster_matrix<=50)), 2, sum   )
    val_51_100      <- apply( cluster_matrix*((cluster_matrix>=51) & (cluster_matrix<=100)), 2, sum   )
    val_over_100    <- apply( cluster_matrix*((cluster_matrix>=101)) , 2, sum   )
    val_cols <- get_BEAST_cols(6, sat=0.7, bright=0.7, transparency=1)[6:1]
    cluster_vals <- t(cbind(val_singles, val_2_10, val_11_20, val_21_50, val_51_100, val_over_100))
    colnames(cluster_vals) <- epiweek
    rownames(cluster_vals) <- c("1","2-10","11-20","21-50","51-100","100+")
    ltxt <- paste("Cluster size=",rownames(cluster_vals))[6:1]
    #ltxt <- paste(ltxt,"num seqs=",apply(cluster_vals, 1, sum)[6:1],"in num clusts=",apply(cluster_nums, 1, sum)[6:1])
    
    imageName <- paste(dataPath,gsub("\\.fas","",dataName),"__sequences_in_cluster_sizes__",rownames(cluster_max)[i],"_barplot.png",sep="")
    png(file=imageName, height=6*150, width=12*150, res=150)
    barplot(cluster_vals, col=val_cols, xlab="Epi Week", names=epidate)
    title(paste(rownames(cluster_max)[i],sum(cluster_matrix),"sequences in clusters"))
    legend("topleft",pch=22,ltxt,pt.bg=val_cols[6:1],bty="n")
    dev.off()
    
    
    cname <- paste(dataPath,gsub("\\.fas","",dataName),"__cluster_matrix__",rownames(cluster_max)[i],".csv",sep="")
    write(paste("# Cluster matrix to",cname), file=logName, append=TRUE)
    write.table(cluster_matrix,file=cname,quote=FALSE,row.names=FALSE,col.names=TRUE,sep=",")
  }
    
    #cluster_size_across_region <- apply(cluster_matrix, 1, sum)
    #tbl <- table(cluster_size_across_region)
    #tbl2<- tbl*as.integer(rownames(tbl))
    bigClusts <- which(clusts$csize>=10)
    cluster_report <- matrix(0, length(bigClusts), 1+1+length(epiweek)+length(ucont)+length(u_uk_state)+length(u_scot)+1)
    colnames(cluster_report) <- c("ClustNo","ClustSize",paste("EpiWeek",epiweek,sep="_"),ucont,u_uk_state,u_scot,"MainLineage")
    for (j in 1:length(bigClusts)) {
      cinds      <- which(clusts$membership==bigClusts[j])
      clust_meta <- metadata[cinds,]
      clust_epiweek <- table(factor(clust_meta$epi_week, levels=epiweek))
      clust_continent<- table(factor(clust_meta$continent, levels=ucont))
      clust_uk <- table(factor(clust_meta$country_state, levels=u_uk_state))
      clust_scot <- table(factor(clust_meta$uk_county, levels=u_scot))
      clust_lineage <- names(sort(table(paste(clust_meta$lineage)),decreasing=TRUE))[1]
      cluster_report[j,] <- c(bigClusts[j],length(cinds),clust_epiweek,clust_continent,clust_uk,clust_scot,clust_lineage)
    }
    cname <- paste(dataPath,gsub("\\.fas","",dataName),"__cluster_report.csv",sep="")
    write.table(cluster_report,file=cname,quote=FALSE,row.names=FALSE,col.names=TRUE,sep=",")
    write(paste("# Cluster matrix to",cname), file=logName, append=TRUE)
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
  
}


################################################################

do_UK_COG <- FALSE
if (do_UK_COG) {
  dataDate <- "2020-05-08"
  
  COG_clean_2(dataDate=dataDate)
  COG_clean_3(dataDate=dataDate)
  COG_align(dataDate=dataDate)
  COG_tidy_1(dataDate=dataDate)
  
  COG_tidy_2(dataDate=dataDate)
  
}
