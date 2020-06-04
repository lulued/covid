# script to process Sequences of COG metadata - extracting Scotland
# S J Lycett
# 21 May 2020
# 23 May 2020

library(ape)

# packages for tree dating - should all be under treedater, but might need to separarely install
# this is Erik Volz tree dating
#library(quadprog)
#library(limSolve)
#library(treedater)

#Rpath <- "~/Documents//GitLab//RiseFall-Scottish-COVID19//coronavirus_2020_r//"
#Rpath <- ""
source(paste(Rpath,"getEl.R",sep=""))
source(paste(Rpath,"calcDecimalDate.R",sep=""))
source(paste(Rpath,"get_BEAST_cols.R",sep=""))
source(paste(Rpath,"get_udates.R",sep=""))
source(paste(Rpath,"birdSpecies_and_continents.R",sep=""))



####################################################################################
# CUSTOM FUNCTIONS
# used in main functions

get_selection_from_sub_trees <- function(dataPath=dataPath, ext="_sub_tree_dated.nwk", 
                                         usePart=TRUE, ind=1, sep="\\|") {
  fnames    <- dir(dataPath)
  fnames    <- fnames[grep(ext,fnames,fixed=TRUE)]
  listNames <- gsub(ext,"",fnames)
  taxaList  <- vector("list",length(listNames))
  for (i in 1:length(listNames)) {
    tr <- read.tree(paste(dataPath,fnames[i],sep=""))
    trTaxa <- as.matrix(tr$tip.label)
    if (usePart) {
      trTaxa <- as.matrix(apply(trTaxa, 1, getEl, ind=ind, sep=sep))
    }
    taxaList[[i]] <- trTaxa
  }
  return( list(taxaList=taxaList, listNames=listNames)  )
}

find_selection_in_seqs <- function(seqs=seqs, selectedTaxa=selectedTaxa, 
                                   exactMatch=TRUE, expart="hCoV-19/",
                                   usePart=TRUE, ind=1, sep="\\|") {
  taxa <- as.matrix(attributes(seqs)$names)
  if (expart!="") {
    taxa         <- gsub(expart,"",taxa)
    selectedTaxa <- gsub(expart,"",selectedTaxa)
  }
  taxa         <- gsub(" ","_",taxa)
  selectedTaxa <- gsub(" ","_",selectedTaxa)
  if (exactMatch) {
    inds <- match(selectedTaxa, taxa)
    return(inds)
  } else {
    if (usePart) {
      partTaxa <- apply(taxa, 1, getEl, ind=ind, sep=sep)
      inds <- match(selectedTaxa, partTaxa)
      return(inds)
    } else {
      inds <- apply(as.matrix(selectedTaxa), 1, grep, taxa)
      return(inds)
    }
  }
}

find_seqs <- function(seqPath=seqPath,seqName=c("x"),ext=".fas",outPath=outPath,listNames=sel$listNames,taxaList=sel$taxaList,
                      exactMatch=TRUE, usePart=TRUE, ind=1, sep="\\|") {
  fnames <- dir(seqPath)
  if (seqName[1]=="x") {
    fnames <- fnames[grep(ext,fnames,fixed=TRUE)]
  } else {
    fnames <- seqName
  }
  
  outNames <- paste(outPath,listNames,"_selected",ext,sep="")
  apps     <- array(FALSE, length(taxaList))
  
  # for each of the potentially big alignments in the directory
  for (s in 1:length(fnames)) {
    seqs <- read.dna(paste(seqPath,fnames[s],sep=""), format="fasta", as.matrix=FALSE)
    
    # for each of the selections
    for (i in 1:length(taxaList)) {
      sinds <- find_selection_in_seqs(seqs=seqs, selectedTaxa=taxaList[[i]], 
                                      exactMatch=exactMatch, 
                                      usePart=usePart, ind=ind, sep=sep)
      jj    <- which(is.finite(sinds))
      sinds <- sinds[jj]
      if (length(sinds)>0) {
        write.dna(seqs[sinds],file=outNames[i],format="fasta",nbcol=-1,colsep="",append=apps[i])
        print(paste("Wrote",length(sinds),"found sequences to",outNames[i]))
        apps[i] <- TRUE
      } else {
        print(paste("Warning didnt find any for",listNames[i],"in",fnames[s]))
      }
    }
    print(paste("Finished searching in",fnames[s]))
  }
  
}

check_found_seqs <- function(ext=".fas",outPath=outPath,
                             listNames=sel$listNames,taxaList=sel$taxaList,
                             isCountry="Scotland") {
  
  outNames <- paste(outPath,listNames,"_selected",ext,sep="")
  res      <- matrix(0, length(outNames), 7)
  colnames(res) <- c("Number_in_List","Number_in_SeqFile","Number_Found","Number_Missing",
                     paste("Number_in_List",isCountry,sep="_"),
                     paste("Number_Found",isCountry,sep="_"),
                     paste("Number_Missing",isCountry,sep="_"))
  rownames(res) <- listNames
  fnames <- dir(outPath)
  fnames <- fnames[grep("_selected",fnames)]
  for (i in 1:length(outNames)) {
    pos <- grep(paste(listNames[i],"_selected",sep=""),fnames)
    if (length(pos)>0) {
      found_seqs <- read.dna(outNames[i],format="fasta",as.matrix=FALSE)
      found_inds <- find_selection_in_seqs(seqs=found_seqs, selectedTaxa=taxaList[[i]], 
                          exactMatch=exactMatch, 
                          usePart=usePart, ind=ind, sep=sep)
      jj         <- which(is.finite(found_inds))
      kk         <- which(!is.finite(found_inds))
    
      list_country <- apply(as.matrix(taxaList[[i]]), 1, getEl, ind=1, sep="/")
      num_country  <- length(which(list_country==isCountry))
      if (length(jj>0)) {
        num_found_country <- length(which(list_country[jj]==isCountry))
      } else {
        num_found_country <- 0
      }
      if (length(kk>0)) {
        num_not_found_country <- length(which(list_country[kk]==isCountry))
      } else {
        num_not_found_country <- 0
      }
      
    
      res[i,] <- c(length(taxaList[[i]]),length(found_seqs),length(jj),length(kk),
                 num_country,num_found_country,num_not_found_country)
    } else {
      list_country <- apply(as.matrix(taxaList[[i]]), 1, getEl, ind=1, sep="/")
      num_country  <- length(which(list_country==isCountry))
      num_found_country <- 0
      num_not_found_country <- num_country
      res[i,] <- c(length(taxaList[[i]]),0,0,length(taxaList[[i]]),
                   num_country,num_found_country,num_not_found_country)
    }
  }
  return(res)
}

####################################################################################

####################################################################################
# MAIN FUNCTIONS

# GLOBAL SET UP (not a function)
doSetUp <- TRUE
echo    <- TRUE
if (doSetUp) {
  if (echo) {
    print("PART 0")
    print("Setup the path names")
  }
  
  # SETUP
  procType  <- c("scotland_subtrees_and_related","scotland_region_only")
  procExt   <- c("_sub_tree_dated.nwk",".nwk")

  # COG METADATA and PROCESSED METADATA
  dataDate  <- "2020-05-15"
  #rootPath <- "Documents//data//Coronavirus_2020//"
  rootPath  <- "E://data//Coronavirus_2020//"

  # SEQUENCE DATA
  # may not be in the same path - could use GISAID data for example - especially the msa 
  seqType  <- "MSA" #"GISAID" #COG"
  seqDate  <- "2020_may_21_gisaid" #2020-05-08" 
  sp2      <- "alignments_and_proteins//msa_0520" #"United_Kingdom" # ""
  seqPath  <- paste(rootPath,seqDate,"//",sp2,"//",sep="")
  if (seqType=="COG") {
    exactMatch <- TRUE
    usePart    <- FALSE
  } else {
    exactMatch <- FALSE
    usePart    <- TRUE
  }
  ind <- 1
  sep <- "\\|"
  
  simpleCutPos <- 1000
  
  pp <- 1
  
  dataPath <- paste(rootPath,"COG_downloads//",dataDate,"_processed//",procType[pp],"//",sep="")
  outPath <- paste(rootPath,seqDate,"//",procType[pp],"//",sep="")
  
  # load extended metadata (see risefall_scotland_COG_metadata_extract.R)
  metaPath <- paste(rootPath,"COG_downloads//",dataDate,"//microreact//",sep="")
  metaName <- paste("cog_global_",dataDate,"_meta_extended",sep="")
  load( paste(metaPath,metaName,".metadata.Rdata",sep=""))
  
}

# MAIN FUNCTION - extract_1 - extracts all the relevant sequences from the msa (as default)
extract_1 <- function() {
  
  if (echo) {
    print("EXTRACT 1")
  }
  
  # step 1 get all the sequences from the msa
  #for (pp in 1:length(pp)) {
  pp <- 1
  
  dataPath <- paste(rootPath,"COG_downloads//",dataDate,"_processed//",procType[pp],"//",sep="")
  outPath <- paste(rootPath,seqDate,"//",procType[pp],"//",sep="")
  
  if (echo) {
    print(paste("Sequence data in path=",seqPath))
    print(paste("Meta and processed data in path=",dataPath))
    print(paste("Outputs to path=",outPath))
  }
  
  startTime <- Sys.time()
  logName <- paste(outPath,"risefall_scotland_COG_sequences_extract_1_log.txt",sep="")
  write("# Log file for risefall_scotland_COG_sequences_extract_1",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("seqPath=",seqPath), file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  write(paste("outPath=",outPath),file=logName,append=TRUE)
  
  if (echo) print("Retrieving lists of sequences to process from subtrees")
  write("# Retrieving lists of sequences to process from subtrees", file=logName,append=TRUE)
  
  sel     <- get_selection_from_sub_trees(dataPath=dataPath,ext=procExt[pp],usePart=usePart, ind=ind, sep=sep)
  if (echo) print(paste("There are",length(sel$listNames),"subtree files with sequences"))
  write(paste("# There are",length(sel$listNames),"subtree files with sequences"), file=logName,append=TRUE)
  
  all_taxa<- unique(unlist(sel$taxaList))
  if (echo) print(paste("There are",length(all_taxa),"unique sequences in all of these files"))
  write(paste("# There are",length(all_taxa),"unique sequences in all of these files"), file=logName,append=TRUE)
  
  taxaList<- vector("list",1)
  taxaList[[1]] <- all_taxa
  combined_name <- "cog_global_selected"
  listNames <- array(combined_name,1)
  sel$taxaList <- taxaList
  sel$listNames <- listNames
  
  if (echo) print("Extracting unique sequences from sequences file(s)")
  write("# Extracting unique sequences from sequences file(s)", file=logName, append=TRUE)
  
  find_seqs(seqPath=seqPath,ext=".fas",outPath=outPath,listNames=sel$listNames,taxaList=sel$taxaList,
            exactMatch=exactMatch, usePart=usePart, ind=ind, sep=sep)
  res     <- check_found_seqs(ext=".fas",outPath=outPath,listNames=sel$listNames,taxaList=sel$taxaList)
  write.csv(res, file=paste(outPath,"check_found_all_sequences.csv",sep=""))
  save(dataPath, outPath, seqPath, sel, res, file=paste(outPath,"check_found_all_sequences.objs.Rdata",sep=""))
  #}
  
  if (echo) {
    print(paste("unique sequences filename= ",combined_name,".fas",sep=""))
    print("check table filename= check_found_all_sequences.csv")
  }
  write(paste("unique sequences filename= ",combined_name,"_selected.fas",sep=""), file=logName, append=TRUE)
  write(paste("check table filename= check_found_all_sequences.csv",sep=""), file=logName, append=TRUE)
  write(paste("R objects filename= check_found_all_sequences.objs.Rdata",sep=""), file=logName, append=TRUE)
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
}

# MAIN FUNCTION - extract_2 - extracts the sequences per individual subtree files
#seqName="cog_cov2020_start_trim_with_ref_mafft.fas"
extract_2 <- function(seqName="cog_cov2020_al3.fas") {
  if (echo) {
    print("EXTRACT 2")
  }
  
  startTime <- Sys.time()
  logName <- paste(outPath,"risefall_scotland_COG_sequences_extract_2_log.txt",sep="")
  write("# Log file for risefall_scotland_COG_sequences_extract_2",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)

  
  for (pp in 1:length(pp)) {
    
    dataPath <- paste(rootPath,"COG_downloads//",dataDate,"_processed//",procType[pp],"//",sep="")
    outPath <- paste(rootPath,seqDate,"//",procType[pp],"//",sep="")
    seqPath <- paste(rootPath,seqDate,"//",procType[1],"//",sep="")
    
    if (echo) {
      print(paste("# Extracting",procType[pp],"from",seqName,"in",seqPath))
    }
    write(paste("# Extracting",procType[pp],"from",seqName,"in",seqPath), file=logName,append=TRUE)
    write(paste("seqPath=",seqPath), file=logName, append=TRUE)
    write(paste("seqName=",seqName), file=logName, append=TRUE)
    write(paste("dataPath=",dataPath),file=logName,append=TRUE)
    write(paste("outPath=",outPath),file=logName,  append=TRUE)
    
    sel     <- get_selection_from_sub_trees(dataPath=dataPath,ext=procExt[pp],usePart=usePart, ind=ind, sep=sep)
    find_seqs(seqPath=seqPath,seqName=seqName,ext=".fas",outPath=outPath,listNames=sel$listNames,taxaList=sel$taxaList,
              exactMatch=exactMatch, usePart=usePart, ind=ind, sep=sep)
    res     <- check_found_seqs(ext=".fas",outPath=outPath,listNames=sel$listNames,taxaList=sel$taxaList)
    write.csv(res, file=paste(outPath,"check_found_sequences.csv",sep=""))
    save(dataPath, outPath, seqPath, sel, res, file=paste(outPath,"check_found_sequences.objs.Rdata",sep=""))
    
    if (echo) print("check table filename= check_found_sequences.csv")
    write(paste("check table filename= ",".check_found_sequences.csv",sep=""), file=logName, append=TRUE)
    write(paste("R objects filename= ",".check_found_sequences.objs.Rdata",sep=""), file=logName, append=TRUE)
    
  }
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
}

extract_3 <- function() {
  
  if (echo) {
    print("EXTRACT 3")
  }
  
  startTime <- Sys.time()
  logName <- paste(outPath,"risefall_scotland_COG_sequences_extract_3_log.txt",sep="")
  write("# Log file for risefall_scotland_COG_sequences_extract_3",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  
  # get the long sequence names with admin-2 and decimal dates
  pp <- 1
  sel     <- get_selection_from_sub_trees(dataPath=dataPath,ext=procExt[pp],usePart=FALSE, ind=ind, sep=sep)
  if (echo) print(paste("There are",length(sel$listNames),"subtree files with sequences"))
  write(paste("# There are",length(sel$listNames),"subtree files with sequences"), file=logName,append=TRUE)
  
  all_taxa<- unique(unlist(sel$taxaList))
  all_ids <- apply(as.matrix(all_taxa), 1, getEl, ind=ind, sep=sep)
  all_country <- apply(as.matrix(all_ids), 1, getEl, ind=1, sep="/")
  if (echo) print(paste("There are",length(all_taxa),"unique sequences in all of these files"))
  write(paste("# There are",length(all_taxa),"unique sequences in all of these files"), file=logName,append=TRUE)
  
  #get the aligned sequences file
  seqPath <- paste(rootPath,seqDate,"//",procType[1],"//",sep="")
  seqName <- "cog_cov2020_al3.fas"
  
  write(paste("seqPath=",seqPath),file=logName,append=TRUE)
  write(paste("seqName=",seqName),file=logName,append=TRUE)
  
  seqs    <- read.dna(paste(seqPath,seqName,sep=""), format="fasta", as.matrix=FALSE)
  taxa    <- as.matrix( attributes(seqs)$names )
  s_ids   <- apply(taxa, 1, getEl, ind=ind, sep=sep)
  s_ids   <- gsub("hCoV-19/","",s_ids)
  s_ids   <- gsub(" ","_",s_ids)
  
  # rename the sequences
  minds   <- match(s_ids, all_ids)
  all(s_ids==all_ids[minds])
  newTaxa <- all_taxa[minds]
  attributes(seqs)$names <- newTaxa
  
  kinds <- match(newTaxa, gsub(" ","_",metadata$newTaxa))
  all(newTaxa==gsub(" ","_",metadata$newTaxa[kinds]))
  sel_meta <- metadata[kinds,]
  
  sel_country <- "Scotland"
  num_country <- length(which(all_country==sel_country))
  num_found_country <- length(which(sel_meta$country_state==sel_country))
  
  not_found <- match(all_ids, s_ids)
  not_found <- not_found[which(is.finite(not_found))]
  not_found <- setdiff(1:length(all_ids),not_found)
  num_not_found_country <- length(which(all_country[not_found]==sel_country))
  
  res <- matrix(0, 1, 7)
  colnames(res) <- c("Number_in_List","Number_in_SeqFile","Number_Found","Number_Missing",        
                     "Number_in_List_Scotland","Number_Found_Scotland","Number_Missing_Scotland")
  res[1,] <- c(length(all_taxa),length(s_ids),length(which(is.finite(minds))),length(all_taxa)-length(s_ids),
               num_country,num_found_country,num_not_found_country)
  
  if (echo) {
    print(t(res))
  }
  write("# sequences found and aligned table:",file=logName, append=TRUE)
  write.table(t(res), file=logName, append=TRUE, col.names=FALSE, row.names=TRUE, quote=FALSE, sep="=")

  write.csv(res, file=paste(outPath,combined_name,"_check_found_all_sequences.csv",sep=""), row.names=FALSE)
  write.table(sel_meta, file=paste(outPath,combined_name,"_extended_metadata.txt",sep=""), 
              col.names=TRUE, row.names=FALSE, quote=FALSE)
  
  save(seqs, file=paste(outPath,combined_name,".seqs.Rdata",sep=""))
  write.dna(seqs, file=paste(outPath,combined_name,".fas",sep=""),format="fasta",nbcol=-1, colsep="")
  
  write(paste("renamed taxa sequence file name=",combined_name,".fas",sep=""), file=logName, append=TRUE)

  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
  
}

#####################################################################
# PIPELINE

doPipeline <- FALSE

if (doPipeline) {
  # set up done in compile

  extract_1()
  
  # need to clean the selected msa because is not very good
  COG_clean_2(dataDate=dataDate,dataPath=outPath,dataName=paste(combined_name,"_selected.fas",sep=""))
  COG_clean_3(dataDate=dataDate,dataPath=outPath,simpleCutPos=simpleCutPos)
  
  #print("Do alignment with MAFFT outside of R")
  # can use COG_align instead - this just takes off the added reference sequence re-calcs consensus
  COG_align(dataDate=dataDate,dataPath=outPath)
  
  # tidy the alignment with de-dapping and remove consensus etc
  COG_tidy_1(dataDate=dataDate,dataPath=outPath)
  
  # add gaps to make coding
  COG_tidy_2(dataDate=dataDate,dataPath=outPath)
  
  # extract the cleaned sequences according to the sub_trees
  extract_2(seqName="cog_cov2020_al3.fas")
  
  # rename all of the sequences to beast names with admin-2 and also get extended metadata
  extract_3()
  
  # calculate quick nj tree and genetic distance - useful for clustering to do on all ?
  COG_tree(dataDate=dataDate,dataPath=outPath,dataName="cog_global_scotland_focus.fas")
  
  COG_clusters(dataPath=outPath,dataName="cog_global_scotland_focus.fas",dataDate=dataDate,rootPath=rootPath)
  
}