# function to process SARS-CoV-2 from GISAID as pipeline - backend - initial data clean
# S J Lycett
# 8  April 2020
# 27 April 2020
# 7 May 2020 - next strain metadata available (contains more details); need to download in batches
# 8 May 2020 - updated calcDecimalDates.R to account for leap year

#####################################################
### coronavirus_2020_gisaid_R_pipeline_cleaning.R ###
#####################################################

# see process_Wuhan_coronavirus_11Mar2020.R for development - some of the things are now extracted here

# Bioconductor packages needed for sequence alignment ?
# no actually not; use mafft it is better, faster and install is easier

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msa")

##################################################################################################
# load libraries
library(ape)
library(seqinr)
#library(ips) - would use this if using mafft from within R


##################################################################################################
# load utility functions

# Sam HP workstation only - setwd
# setwd("~//Rcode//Coronavirus_2020_R//")

Rpath <- "./lib/"
source(paste(Rpath,"getEl.R",sep=""))
source(paste(Rpath,"get_BEAST_cols.R",sep=""))
source(paste(Rpath,"calcDecimalDate.R",sep=""))
source(paste(Rpath,"get_udates.R",sep=""))


##################################################################################################
# define custom functions

# function to read strain details extract from the acknowledgement table
# this will read the file name as specified in rootname (readtable=TRUE)
# or will read all the strain details tables in the directory - useful if need to download in batches (readtable=FALSE)
get_strain_info <- function(dataPath, 
                            rootname="gisaid_cov2020_acknowledgement_table", 
                            ext="_strain_details.txt",
                            readtable=FALSE) {
  if (readtable) {
    info <- read.table(paste(dataPath,rootname,ext,sep=""), sep="\t", header=TRUE)
  } else {
    # get all strain_details tables
    fnames <- dir(dataPath)
    inds   <- grep("acknowledgement_table_strain_details.txt", fnames, fixed=TRUE)
    fnames <- fnames[inds]
    
    outname <- "all_acknowledgement_table_strain_details.txt"
    pos     <- grep(outname,fnames,fixed=TRUE)
    if (length(pos)>0) {
      info <- read.table(paste(dataPath,outname,sep=""), sep="\t", header=TRUE)
    } else {
      for (j in 1:length(fnames)) {
        part_info <- readLines(paste(dataPath,fnames[j],sep=""))
        part_info <- gsub("'","",part_info,fixed=TRUE)
        part_info <- apply(as.matrix(part_info), 1, getEl, sep="\t")
        cn        <- part_info[,1]
        part_info <- part_info[,2:length(part_info[1,])]
        part_info <- t(part_info)
        colnames(part_info) <- cn
        if (j==1) {
          write.table(part_info, file=paste(dataPath,outname,sep=""), sep="\t", 
                      col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)
        } else {
          write.table(part_info, file=paste(dataPath,outname,sep=""), sep="\t",
                      col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
        }
        info <- read.table(paste(dataPath,outname,sep=""), sep="\t", header=TRUE)
        
        # de-dup
        #uepi <- table(info[,1])
        #jj   <- which(uepi>1)
        #non_u_epi <- rownames(uepi)[jj]
        #ex_inds <- c()
        #if (length(non_u_epi)>0) {
        #  for (j in 1:length(non_u_epi)) {
        #    kk <- which(info[,1]==non_u_epi[j])
        #    nkk<- length(kk)
        #    ex_inds <- c(ex_inds,kk[1:(nkk-1)])
        #  }
        #}
        #inc_inds <- 1:length(info[,1])
        #inc_inds <- setdiff(inc_inds, ex_inds)
        #info <- info[inc_inds,]
        #write.table(info, file=paste(dataPath,outname,sep=""), sep="\t",
         #            col.names=FALSE, row.names=FALSE, quote=FALSE, append=FALSE)
        
        #info <- read.table(paste(dataPath,outname,sep=""), sep="\t", header=TRUE)
      }
    }
    
  }
  return(info)
}


# get metadata from bulk download this contains more info than the acknowledge table
# data from Next Strain - not one entry for all sequences though
get_metadata <- function(dataPath, rootname="metadata", ext=".tsv") {
  
  fnames <- dir(dataPath)
  pos    <- grep(".metadata.Rdata", fnames, fixed=TRUE)
  
  if (ext==".tsv") {
    if (length(pos)==0) {
      fname    <- paste(dataPath,rootname,ext,sep="")
  
      # read as text lines 
      txtLines <- readLines(fname)
  
      # remove the evil '
      txtLines <- gsub("'","",txtLines)
  
      # split to tab
      metadata <- apply(as.matrix(txtLines), 1, getEl, sep="\t")
      metadata <- t(metadata)
      colnames(metadata) <- metadata[1,]
      metadata <- metadata[2:length(metadata[,1]),]
      metadata[,15] <- gsub(",","",metadata[,15])
  
      # re-write to cleaned file
      #this causes errors on windows due to French acute e symbols
      #write.table(metadata, file=paste(dataPath,rootname,"_tab.txt",sep=""), sep="\t",
      #          col.names=TRUE, row.names=FALSE, quote=FALSE)
      save(metadata, file=paste(dataPath,rootname,".metadata.Rdata",sep=""))
    } else {
      load(paste(dataPath,rootname,".metadata.Rdata",sep=""))
    }
    
  } else {
    metadata <- read.table(paste(dataPath,rootname,ext,sep=""), sep="\t", header=TRUE)
  }
  return( metadata )
  
}


# get all taxa from the individual batch downloaded fasta files from GISAID (unprocessed data)
get_all_taxa <- function(dataPath=dataPath,n2="submission_gisaid_hcov-19",ext=".fasta", echo=FALSE) {
  fnames <- dir(dataPath)
  
  pos <- grep("all_taxa.info.Rdata",fnames,fixed=TRUE)
  if (length(pos)==0) {
  
    fnames <- fnames[grep(n2,fnames,fixed=TRUE)]
    fnames <- fnames[grep(ext,fnames,fixed=TRUE)]
    download_date    <- apply(as.matrix(fnames), 1, getEl, ind=c(4,3,2), fromEnd=TRUE, sep="_", reconstruct=TRUE)
    download_hour    <- apply(as.matrix(fnames), 1, getEl, ind=1, fromEnd=TRUE, sep="_")
    download_decDate <- apply(as.matrix(download_date), 1, calcDecimalDate_fromTxt, sep="_", dayFirst=FALSE, namedMonths=FALSE)
  
    download_order <- sort(download_decDate, index.return=TRUE)$ix
    for (j in 1:length(download_order)) {
      k    <- download_order[j]
      seqs <- read.dna(paste(dataPath,fnames[k],sep=""), format="fasta", as.matrix=FALSE)
      taxa <- attributes(seqs)$names
      taxa <- as.matrix(taxa)
      epi_isl <- apply(taxa, 1, getEl, ind=2, sep="\\|")
      latest_seq<- array(1, length(taxa))
    
      part_info <- cbind(epi_isl, taxa, array(fnames[k],length(taxa)), latest_seq)
      rm(seqs)
    
      if (j>1) {
        minds <- match(info[,1], epi_isl)
        jj    <- which(is.finite(minds))
        if (length(jj)>0) {
          info[jj,4] <- 0
          if (echo) {
            print( all( info[jj,1]==epi_isl[minds[jj]] ) )
            print("Sequence updates")
            print(info[jj,])
          }
        }
      
        info <- rbind(info,part_info)
      } else {
        info <- part_info
      }
    }
    save(info, file=paste(dataPath,"all_taxa.info.Rdata",sep=""))
  }
  
  load( paste(dataPath,"all_taxa.info.Rdata",sep="") )
  
  return(info)
}

get_orig_dataNames <- function(dataPath=dataPath) {
  load(paste(dataPath,"all_taxa.info.Rdata",sep=""))
  dataName <- unique(info[,3])
  return( dataName )
}

# load the sequences corresponding to the all_info
load_all_seqs <- function(dataPath=dataPath, all_dataName="combined_all_seqs.fas") {
  
  load(paste(dataPath,"all_taxa.info.Rdata",sep=""))
  dataName <- unique(info[,3])
  fnames   <- dir(dataPath)
  
  pos <- grep("all_seqs.all_seqs.Rdata",fnames,fixed=TRUE)
  if (length(pos)==0) {
    pos <- grep(all_dataName,fnames,fixed=TRUE)
    outname  <- paste(dataPath,all_dataName,sep="")
    if (length(pos)==0) {
      seqFiles <- paste(gsub("\\.fasta","",dataName),".seqs.Rdata",sep="")
      for (j in 1:length(dataName)) {
        pos     <- grep(seqFiles[j], fnames)
        if (length(pos)==0) {
          seqs    <- read.dna(paste(dataPath,dataName[j],sep=""), format="fasta", as.matrix=FALSE)
          ok_taxa <- info[which(info[,3]==dataName[j] & as.integer(info[,4])==1),2]
          ok_inds <- match(ok_taxa, attributes(seqs)$names)
          seqs    <- seqs[ok_inds]
          all( ok_taxa==attributes(seqs)$names )
          save(seqs, file=paste(dataPath,seqFiles[j],sep=""))
        } else {
          load( paste(dataPath,seqFiles[j],sep="") )
        }
        write.dna(seqs, file=outname, format="fasta", nbcol=-1, colsep="", append=(j>1))
      }
    }
    all_seqs <- read.dna(outname, format="fasta", as.matrix=FALSE)
    save(all_seqs, file=paste(dataPath,"all_seqs.all_seqs.Rdata",sep=""))
  } else {
    load( paste(dataPath,"all_seqs.all_seqs.Rdata",sep="") )
  }
  all_taxa    <- attributes(all_seqs)$names
  latest_inds <- which(as.integer(info[,4])==1)
  all( all_taxa == info[latest_inds,2])
  return( all_seqs )
}

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

# 27 April 2020 - dont need to define this now because have put into a separate file get_udates.R
doDef <- FALSE
if (doDef) {
  
get_udates <- function(min_date=min_date, max_date=max_date) {
  
  days_in_month <- c(31,29,31,30,31,30,31,31,30,31,30,31)
  month_names   <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
 # min_date <- as.integer(strsplit(paste(info$Collection.date[min_ii]),"-")[[1]])
#  max_date <- as.integer(strsplit(paste(info$Collection.date[max_ii]),"-")[[1]])
  
  udate <- c()
  seldates <- c()
  seldates2<- c()
  for (y in min_date[1]:max_date[1]) {
    min_mon <- 1
    if (y==min_date[1]) {
      min_mon <- min_date[2]
    }
    
    max_mon <- 12
    if (y==max_date[1]) {
      max_mon <- max_date[2]
    }
    
    for (m in min_mon:max_mon) {
      min_day <- 1
      if (y==min_date[1] & m==min_date[2]) {
        min_day <- min_date[3]
      }
      
      max_day <- days_in_month[m]
      if (m<10) {
        temp1 <- paste(y,"-0",m,"-",sep="")
      } else {
        temp1 <- paste(y,"-",m,"-",sep="")
      }
      
      if (min_day > 10) {
        temp2<- paste(temp1,min_day:max_day,sep="")
      } else {
        temp2 <- c(paste(temp1,"0",1:9,sep=""),paste(temp1,10:max_day,sep=""))
      }
      udate <- c(udate, temp2)
      
      if (y >= 2020) {
        seldates <- c(seldates, paste(temp1,"01",sep=""), paste(temp1,"15",sep=""))
        seldates2<- c(seldates2, paste(month_names[m],"-01",sep=""), paste(month_names[m],"-15",sep=""))
      }
    }
  }
  linds <- match(seldates, udate)
  xlabs <- array(NA,length(udate))
  xlabs[linds] <- seldates2
  
  return(list(udate=udate,seldates=seldates,seldates2=seldates2,xlabs=xlabs))
  
}

} # end doDef

##################################################################################################
# MAIN PIPELINE FUNCTIONS
##################################################################################################

GISAID_clean_1 <- function( dataDate="2020_apr_8",
                            rootPath="E://data//Coronavirus_2020//",
                            dataPath=paste(rootPath,dataDate,"_gisaid//",sep=""),
                            dataName="gisaid_cov2020_sequences.fasta", use_sub_parts=TRUE,n2="submission_gisaid_hcov-19",
                            echo=TRUE, incMetaData=TRUE
                  ) {
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"GISAID_clean_1_log.txt",sep="")
  write("# Log file for GISAID_clean_1",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  
  if (use_sub_parts) {
    write("# using several fasta files by submission date", file=logName,append=TRUE)
  } else {
    write(paste("dataName=",dataName),file=logName,append=TRUE)
  }
  
  # read and process strain information
  info     <- get_strain_info(dataPath,readtable=FALSE)
  uepi <- unique(as.matrix(info[,1]))
  minds<- match(uepi,info[,1])
  info <- info[minds,]
  
  if (echo) {
    print(paste("Number of rows of strain information=",length(info[,1])))
  }
  
  write("# Read strain information",file=logName,append=TRUE)
  write(paste("Number of rows of strain information table=",length(info[,1])),file=logName,append=TRUE)
  
  if (incMetaData) {
    #metadata also available
    metadata <- get_metadata(dataPath)
    if (echo) {
      print(paste("Number of rows of metadata=",length(metadata[,1])))
    }
    write(paste("Number of rows of metadata=",length(metadata[,1])),file=logName,append=TRUE)
  }

  if (use_sub_parts) {
  # 7 May 2020 - get all taxa from individual fasta files in the directory
    taxa_info <- get_all_taxa(dataPath=dataPath,n2=n2)
    minds <- which(taxa_info[,4]==1)
    taxa <- as.matrix(taxa_info[minds,2])
  } else {
    seqs <- read.dna( paste(dataPath,dataName,sep=""), format="fasta", as.matrix=FALSE)
    taxa <- as.matrix(attributes(seqs)$names)
    rm(seqs)
  }
  
  taxa_accn <- apply(taxa, 1, getEl, ind=2, sep="\\|")
  

  if (echo) {
    print(paste("Number of sequences=",length(taxa_accn)))
  }
  
  write("# Read sequences", file=logName, append=TRUE)
  write(paste("Number of sequences",length(taxa_accn)),file=logName,append=TRUE)
  
  diff_accn <- setdiff(taxa_accn, info[,1])
  if (echo) {
    print(paste("There are",length(diff_accn),"sequences which have no corresponding strain information table entry"))
  }
  write(paste("There are",length(diff_accn),"sequences which have no corresponding strain information table entry"),
        file=logName,append=TRUE )
  
  diff_accn <- setdiff(info[,1],taxa_accn)
  if (echo) {
    print(paste("There are",length(diff_accn),"strain information table entries which have no corresponding sequences"))
  }
  write(paste("There are",length(diff_accn),"strain information table entries which have no corresponding sequences"),
        file=logName,append=TRUE)
  
  # need metadata to get this right
  host     <- array("Human",length(info[,1]))
  host[grep("[Bb]at",info[,2])]      <- "Bat"
  host[grep("[Pp]angolin",info[,2])] <- "Pangolin"
  host[grep("[Ee]nv",info[,2])]      <- "Environment"
  host[grep("[Mm]ink",info[,2])]     <- "Mink"
  host[grep("[Tt]iger",info[,2])]    <- "Tiger"
  host[grep("[Cc]anine",info[,2])]   <- "Canine"
  host[grep("[Uu]nknown",info[,2])]  <- "Unknown"
  host[grep("[Ff]elis catus",info[,2])]  <- "Cat"
  host[grep("h[Cc]o[vV]-19/cat",info[,2])] <- "Cat"
  
  if (incMetaData) {
    cn <- colnames(metadata)
    metadata[,which(cn=="host")] <- gsub("human","Human",metadata[,which(cn=="host")])
    uinds <- grep("[Uu]nknown",metadata[,which(cn=="host")])
    kk <- match(metadata[uinds,which(cn=="gisaid_epi_isl")], info[,1])
    host[kk] <- "Unknown"
    non_human_inds <- which(metadata[,which(cn=="host")] != "Human" )
    non_human_epi_isl <- metadata[non_human_inds,which(cn=="gisaid_epi_isl")]
    kk <- match(non_human_epi_isl, info[,1])
  }

  
  if (echo) {
    print(table(host))
  }
  #  4 Apr 2020
  #  host
  #  Bat Environment       Human    Pangolin 
  #  1           5         3849           9   
  
  if (echo & incMetaData) {
    print( table( metadata[non_human_inds,15], host[kk]) )
  }
  
  write("# Host table",file=logName,append=TRUE)
  write.table(table(host),file=logName,append=TRUE,col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
  
  if (incMetaData) {
    write("# Host latin names from metadata", file=logName, append=TRUE)
    write.table( table( metadata[non_human_inds,15], host[kk]), 
               file=logName,append=TRUE,col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t" )
  }
  
  write("# Processing decimal date",file=logName,append=TRUE)
  decDate  <- as.numeric(apply(as.matrix(info$Collection.date), 1, calcDecimalDate_fromTxt, dayFirst=FALSE, namedMonth=FALSE, sep="-"))
  dateLen <- apply(as.matrix(info$Collection.date), 1, nchar)
  hasFullDate <- dateLen>=10
  hasMonthOnly<- dateLen==7
  hasYearOnly <- dateLen==4
  
  write("# Processing Location to continent, country, state",file=logName,append=TRUE)
  location <- as.matrix(info$Location)
  location <- gsub(" /","/",location)
  location <- gsub("/ ","/",location)
  location <- gsub(" /","/",location) # do this twice because some have double space
  location <- gsub("/ ","/",location)
  location <- gsub("Europe/England", "Europe/United Kingdom/England",location)
  location[which(location=="China/Wuhan")] <- "Asia/China/Hubei/Wuhan"
  location[which(location=="Malaysia")] <- "Asia/Malaysia"
  location[which(location=="USA/LA")] <- "North America/USA/LA"
  location <- gsub("Chongqinq","Chongqing",location)
  location <- gsub("Guandong","Guangdong",location)
  location <- gsub("Gen?ve","Geneva",location)
  location <- gsub("Compi?gne","Compiegne",location)
  location <- gsub("?","e",location)
  continent<- apply(location, 1, getEl, ind=1, sep="/")
  country  <- apply(location, 1, getEl, ind=2, sep="/")
  state    <- apply(location, 1, getEl, ind=3, sep="/")
  
  location <- gsub("Eurasia/Russia/Moscow ","Eurasia/Russia/Moscow",location)
  russia_inds <- which(location=="Eurasia/Russia/Moscow")
  continent[russia_inds] <- "Europe"
  
  if (echo) {
    print(table(continent))
  }
  #  11 Mar 2020
  #  continent
  # Africa  Asia Central America          Europe   North America         Oceania   South America 
  # 42      668               1            1772             949             384              48
  
  write("# Continent table",file=logName,append=TRUE)
  write.table(table(continent),file=logName,append=TRUE,col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
  
  if (echo) {
    pie(table(continent),col=colourize_BEAST_cols(continent)$ucols,main="Number of entries")
  }
  
  write("# making new sequence names",file=logName,append=TRUE)
  seqName  <- paste(info$Accession.ID,host,country,info$Virus.name,info$Collection.date,format(decDate,digits=7),sep="|")
  
  info <- cbind(seqName,info,location,continent,country,state,host,hasYearOnly,hasMonthOnly,hasFullDate,decDate)
  write.table(info, file=paste(dataPath,"all_expanded_strain_info.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  save(info, file=paste(dataPath,"all_expanded_strain_info.info.Rdata",sep=""))
  
  write("# writing info to file",file=logName,append=TRUE)
  write(paste("info file name=",paste(dataPath,"all_expanded_strain_info.info.Rdata",sep="")),
        file=logName,append=TRUE)
  
  
  ok_inds <- which(info$hasFullDate & info$host=="Human")
  
  write(paste("Number of human sequences with full date=",length(ok_inds)), 
        file=logName,append=TRUE)
  
  maxDate   <- max(decDate[ok_inds])
  max_ii    <- ok_inds[which.max(decDate[ok_inds])]
  min_ii    <- ok_inds[which.min(decDate[ok_inds])]
  min_date  <- as.integer(strsplit(paste(info$Collection.date[min_ii]),"-")[[1]])
  if (min(decDate[ok_inds]) > 2020) {
    min_date<- c(2020,1,1)
  }
  
  max_date  <- as.integer(strsplit(paste(info$Collection.date[max_ii]),"-")[[1]])
  dates_res <- get_udates(min_date=min_date,max_date=max_date)
  
  udate    <- dates_res$udate
  seldates <- dates_res$seldates
  seldates2<- dates_res$seldates2
  xlabs    <- dates_res$xlabs
  
  write(paste("Min ok date=",info$Collection.date[min_ii]), 
        file=logName,append=TRUE)
  write(paste("Max ok date=",info$Collection.date[max_ii]), 
        file=logName,append=TRUE)
  
  ccols <- colourize_BEAST_cols(paste(info$continent[ok_inds]))
  cont_tbl <- table(factor(info$continent[ok_inds],levels=ccols$ustates))
  
  imageName <- paste(dataPath,"all_expanded_strain_info_continent_pie.png",sep="")
  png(file=imageName, height=6*300, width=8*300, res=300)
  pie(cont_tbl, col=ccols$ucols)
  legend("bottomleft",paste(ccols$ustates,"(",cont_tbl,")"),pch=22,bty="n",pt.bg=ccols$ucols)
  title(paste("Number of sequences",toupper(gsub("_"," ",dataDate))))
  dev.off()
  
  imageName <- paste(dataPath,"all_expanded_strain_info_date_vs_continent_samples_6x12.png",sep="")
  png(file=imageName, height=6*300, width=12*300, res=300)
  barplot(  table( info$continent[ok_inds], factor(info$Collection.date[ok_inds],levels=udate ) ), 
            names=xlabs, col=ccols$ucols, xlab="Collection date", ylab="Number of sequences" )
  legend("topleft",paste(ccols$ustates,"(",cont_tbl,")"),pch=22,bty="n",pt.bg=ccols$ucols)
  title(paste("Number of sequences",toupper(gsub("_"," ",dataDate))))
  dev.off()
  
  write("# written all_expanded_strain_info_continent_pie.png", file=logName, append=TRUE)
  write("# written all_expanded_strain_info_date_vs_continent_samples_6x12.png", file=logName, append=TRUE)
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
  
}

GISAID_clean_2 <- function( dataDate="2020_apr_8",
                            rootPath="E://data//Coronavirus_2020//",
                            dataPath=paste(rootPath,dataDate,"_gisaid//",sep=""),
                            dataName="gisaid_cov2020_sequences.fasta",
                            use_sub_parts=TRUE,
                            echo=TRUE
                            ) {
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"GISAID_clean_2_log.txt",sep="")
  write("# Log file for GISAID_clean_2",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  
  
  write("# reading all_expanded_strain_info.info.Rdata", file=logName,append=TRUE)
  load(file=paste(dataPath,"all_expanded_strain_info.info.Rdata",sep=""))
  
  write("# starting to populate all_strains_list", file=logName,append=TRUE)
  all_strains_list <- matrix(-1, length(info[,1]), 7)
  colnames(all_strains_list)  <- c("Accn","HumanFullDate","Long","HighQuality","AlignOK","TempestOK","Include")
  all_strains_list[,1]        <- as.matrix(info$Accession.ID)
  
  ok_inds <- which(info$hasFullDate & info$host=="Human")
  all_strains_list[,2]        <- 0
  all_strains_list[ok_inds,2] <- 1
  
  info <- info[ok_inds,]
  
  if (use_sub_parts) {
    write(paste("# reading sequences from multiple part files"), file=logName, append=TRUE)
    dn   <- get_orig_dataNames(dataPath=dataPath)
    for (j in 1:length(dn)) {
      write(paste("dataName=",dn[j]),file=logName,append=TRUE)
    }
    seqs <- load_all_seqs(dataPath=dataPath)
  } else {
    write(paste("dataName=",dataName),file=logName,append=TRUE)
    write(paste("# reading",dataName), file=logName,append=TRUE)
    seqs <- read.dna( paste(dataPath,dataName,sep=""), format="fasta", as.matrix=FALSE)
  }
  
  
  taxa <- as.matrix(attributes(seqs)$names)
  taxa_accn <- apply(taxa, 1, getEl, ind=2, sep="\\|")
  
  
  inds <- match(info$Accession.ID, taxa_accn)
  if (echo) {
    print(paste("Number of info entries=",length(info[,1])))
    print(paste("Number of sequences=",length(taxa)))
  }
  jj <- which(is.finite(inds))
  all( info$Accession.ID[jj]==taxa_accn[inds[jj]])
  info <- info[jj,]
  inds <- match(info$Accession.ID, taxa_accn)
  if (echo) {
    print(paste("Removed some entries from info"))
    print(paste("Number of info entries=",length(info[,1])))
    print(paste("Number of sequences=",length(taxa)))
  }
  
  write(paste("# Matching sequence accns to human accns metadata"), file=logName, append=TRUE)
  ok <- all( taxa_accn[inds]==info$Accession.ID)
  write(paste("Matching ok=",ok),file=logName,append=TRUE)
  if (echo) {
    print(paste("Matching human sequence accns",ok))
  }
    
  seqs <- seqs[inds]
  taxa <- as.matrix(attributes(seqs)$names)
  taxa_accn <- apply(taxa, 1, getEl, ind=2, sep="\\|")
  
  ok <- all(taxa_accn==info$Accession.ID)
  write(paste("Check matching=",ok), file=logName,append=TRUE)
  if (echo) {
    print(paste("Check matching human sequence accns",ok))
  }
  
  attributes(seqs)$names <- paste(info$seqName)
  
  write("# re-namining human sequences",file=logName,append=TRUE)
  write(paste("Number of human sequences=",length(seqs)),file=logName,append=TRUE)
  write(paste("Start file writing=",Sys.time()),file=logName,append=TRUE)
  
  # 8 May 2020
  # this is very slow and not often use this file - do it in R object instead
  #write.dna(seqs, file=paste(dataPath,"human_cov2020.fas",sep=""), format="fasta", nbcol=-1, colsep="")
  save(seqs, file=paste(dataPath,"human_cov2020.seqs.Rdata",sep=""))
  
  write(paste("Finished file writing=",Sys.time()),file=logName,append=TRUE)
  #write("Sequences file name=human_cov2020.fas",file=logName,append=TRUE)
  write("Sequences file name=human_cov2020.seqs.Rdata",file=logName,append=TRUE)
  if (echo) {
    #print("Human sequences written to human_cov2020.fas")
    print("Human sequences written to human_cov2020.seqs.Rdata")
  }
  
  
  len      <- unlist(lapply(seqs,length))
  lenThres <- 28000
  inds     <- which(len >= lenThres)
  seqs     <- seqs[inds]
  taxa     <- as.matrix(attributes(seqs)$names)
  taxa_accn<- apply(taxa, 1, getEl, ind=1, sep="\\|")
  linds    <- match(taxa_accn, all_strains_list[,1])
  all( all_strains_list[linds,1]==taxa_accn)
  all_strains_list[,3]      <- 0
  all_strains_list[linds,3] <- 1
  write(paste("Length threshold=",lenThres),file=logName,append=TRUE)
  write(paste("Number human sequences > length =",length(linds)), file=logName, append=TRUE)
  
  write(paste("Start file writing=",Sys.time()),file=logName,append=TRUE)
  
  # replace with R object beecause very slow
  #sname    <- paste(dataPath,"human_cov2020_len",lenThres,".fas",sep="")
  #write.dna(seqs, file=sname, format="fasta", nbcol=-1, colsep="")
  
  sname <- paste(dataPath,"human_cov2020_len",lenThres,".seqs.Rdata",sep="")
  save(seqs, file=sname)
  
  write(paste("Finished file writing=",Sys.time()),file=logName,append=TRUE)
  #write(paste("Sequences file name=human_cov2020_",lenThres,".fas",sep=""),file=logName,append=TRUE)
  write(paste("Sequences file name=human_cov2020_",lenThres,".seqs.Rdata",sep=""),file=logName,append=TRUE)
  
  # read-read the sequences (not aligned)
  #seqs     <- read.dna(sname, format="fasta", as.matrix=FALSE)
  taxa     <- attributes(seqs)$names
  if (echo) {
    print(paste("There are",length(taxa),"human sequences with full date and length >=",lenThres))
  }
  write(paste("# There are",length(taxa),"human sequences with full date and length >=",lenThres),
        file=logName,append=TRUE)
  
  # quality
  write(paste("# Start quality=",Sys.time()),file=logName,append=TRUE)
  nuclsTbl <- matrix(0, length(taxa), 7)
  colnames(nuclsTbl) <- c("a","c","g","t","-","n","x")
  rownames(nuclsTbl) <- taxa   # added on 29 april 2020 11.30
  for (j in 1:length(taxa)) {
    nuclsTbl[j,] <- numAmbs( unlist(as.character(seqs[j])) )
  }
  write.csv(nuclsTbl,file=paste(dataPath,"human_cov2020_",lenThres,"_nuclsTbl.csv",sep=""))
  write(paste("nuclsTbl=human_cov2020_",lenThres,"_nuclsTbl.csv",sep=""),file=logName,append=TRUE)
  
  
  nthres <- 300 #28000*0.005
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
  taxa_accn <- apply(as.matrix(taxa), 1, getEl, ind=1, sep="\\|")
  inds <- match(taxa, info$seqName)
  all(taxa==info$seqname[inds])
  info <- info[inds,]
  all(taxa==info$seqName)
  
  
  write(paste("Start file writing=",Sys.time()),file=logName,append=TRUE)
  write.dna(seqs, file=paste(dataPath,"human_cov2020_nthres",nthres,"_unal.fas",sep=""),
            format="fasta", nbcol=-1, colsep="")
  save(info, file=paste(dataPath,"human_cov2020_nthres",nthres,".info.Rdata",sep=""))
  write.table(info,file=paste(dataPath,"human_cov2020_nthres",nthres,"_info.txt",sep=""), sep="\t", 
              col.names=TRUE, row.names=FALSE, quote=FALSE)
  
  linds <- match(taxa_accn, all_strains_list[,1])
  all( all_strains_list[linds,1]==taxa_accn )
  all_strains_list[,4]      <- 0
  all_strains_list[linds,4] <- 1
  
  write.csv(all_strains_list, row.names=FALSE, file=paste(dataPath,"all_strains_list.csv",sep=""))
  
  write(paste("Finished file writing=",Sys.time()),file=logName,append=TRUE)
  write(paste("Sequences file name= human_cov2020_nthres",nthres,"_unal.fas",sep=""),
        file=logName,append=TRUE)
  write(paste("Info file= human_cov2020_nthres",nthres,".info.Rdata",sep=""),
        file=logName,append=TRUE)
  write(paste("Info file= human_cov2020_nthres",nthres,"_info.txt",sep=""),
        file=logName,append=TRUE)
  write("All strains list file= all_strains_list.csv", 
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
GISAID_clean_3 <- function(dataDate="2020_apr_8",
                           rootPath="E://data//Coronavirus_2020//",
                           dataPath=paste(rootPath,dataDate,"_gisaid//",sep=""),
                           refPath ="refSeq//",
                           refName ="COVID-19_human_373_coding_nuclsConsensus.fas", 
                           refFragLen=30,
                           echo=TRUE
                  ) {
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"GISAID_clean_3_log.txt",sep="")
  write("# Log file for GISAID_clean_3",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
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
  
  # 29 May 2020 - also do 1 shifted fragment because of early orf1ab muts (genuine I think)
  refStartFrag_shift <- refNucls[(refFragLen+1):(2*refFragLen)]
  refStartString_shift <- paste(refStartFrag_shift,collapse="")
  
  dataName       <- "human_cov2020_nthres300_unal.fas"
  if (echo) {
    print(paste("Reading unaliged sequences=",dataName))
  }
  write(paste("dataName=",dataName),file=logName,append=TRUE)
  write("# reading unaligned sequences", file=logName,append=TRUE)
  
  x              <- read.dna(paste(dataPath,dataName,sep=""), format="fasta", as.matrix=FALSE)
  taxa           <- as.matrix(attributes(x)$names)
  nseqs          <- length(x)
  
  write("# trimed sequences to human_cov2020_start_trim_with_ref.fas", file=logName, append=TRUE)
  write("# list of non-trimmed sequences to trim_start_ORF1ab_fail.txt", file=logName, append=TRUE)
  write("# non-trimmed sequences to human_cov2020_bad_failed_trim_with_ref.fas", file=logName, append=TRUE)
  write(paste("# Trimming start",Sys.time()), file=logName, append=TRUE)
  
  outName        <- paste(dataPath,"human_cov2020_start_trim_with_ref.fas",sep="")
  write(paste(">",attributes(y)$names,sep=""), file=outName, append=FALSE)
  write(refNuclsString,file=outName,append=TRUE)
  
  badOutName      <- paste(dataPath,"human_cov2020_bad_failed_trim_with_ref.fas",sep="")
  write(paste(">",attributes(y)$names,sep=""), file=badOutName, append=FALSE)
  write(refNuclsString,file=badOutName,append=TRUE)
  
  errorLog <- paste(dataPath,"trim_start_ORF1ab_fail.txt")
  write("# sequences which dont have exact 30 nucls start", file=errorLog, append=FALSE)
  
  if (echo) {
    print("Trimmed sequences to human_cov2020_start_trim_with_ref.fas")
    print("list of non-trimmed sequences to trim_start_ORF1ab_fail.txt")
  }
  
  badCount <- 0
  goodCount<- 0
  for (j in 1:nseqs) {
    x_nucls <- unlist( as.character(x[j]) )
    x_string<- paste(x_nucls,collapse="")
    pos     <- gregexpr(refStartString,x_string)
    if (pos>0) {
      x_string <- substring(x_string,pos)
      write(paste(">",taxa[j],sep=""), file=outName, append=TRUE)
      write(x_string,file=outName,append=TRUE)
      goodCount <- goodCount+1
    } else {
      pos2     <- gregexpr(refStartString2,x_string)
      if (pos2>=1) {
        # 29 May 2020 - used pos2 here (formerly pos which was wrong)
        #x_string <- paste("-",substring(x_string,pos),sep="")
        x_string <- substring(x_string,pos2[[1]]-1)
        write(paste(">",taxa[j],sep=""), file=outName, append=TRUE)
        write(x_string,file=outName,append=TRUE)
        goodCount <- goodCount+1
      } else {
        pos3 <- gregexpr(refStartString_shift,x_string)
        if (pos3>=refFragLen) {
          x_string <- substring(x_string,pos3[[1]]-refFragLen+1)
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
  }
  
  if (echo) {
    print(paste("Number bad=",badCount))
    print(paste("Number good=",goodCount))
  }
  
  write(paste("Number bad=",badCount), file=logName, append=TRUE)
  write(paste("Number good=",goodCount), file=logName, append=TRUE)
  write("trimName= human_cov2020_start_trim_with_ref.fas", file=logName, append=TRUE)
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


##########################################################################################

GISAID_align <- function(dataDate="2020_apr_8",
                         rootPath="E://data//Coronavirus_2020//",
                         dataPath=paste(rootPath,dataDate,"_gisaid//",sep=""),
                         refName="consensus",
                         doConsensus=TRUE,
                         echo=TRUE
                         ) {
  
  fname <- paste(dataPath,"human_cov2020_start_trim_with_ref.fas",sep="")
  print(paste("Use MAFFT to align",fname))
  print("Not implemented on windows")
  print("Output format = 4; method = 1 (auto)")
  print("do this and save as human_cov2020_start_trim_with_ref_mafft.fas")
  logName <- ""
  
  if (doConsensus) {
    startTime <- Sys.time()
    logName <- paste(dataPath,"GISAID_align_log.txt",sep="")
    write("# Log file for GISAID_align",file=logName,append=FALSE)
    write("# START", file=logName,append=TRUE)
    write(paste("# process start",startTime), file=logName,append=TRUE)
    write(paste("dataDate=",dataDate),file=logName,append=TRUE)
    write(paste("rootPath=",rootPath),file=logName,append=TRUE)
    write(paste("dataPath=",dataPath),file=logName,append=TRUE)
    
    if (echo) {
      print("replacing original ref sequence with consensus of the others")
      print("doing this because using already mafft downloaded data")
      print("this data has excess gaps and not all sequences are the same length")
    }
    dataName <- "human_cov2020_start_trim_with_ref"
    
    write(paste("dataName=",dataName),file=logName,append=TRUE)
    
    write("# replacing original ref sequence with consensus of the others", file=logName,append=TRUE)
    write("# doing this because using already mafft downloaded data", file=logName,append=TRUE)
    write("# this data has excess gaps and not all sequences are the same length", file=logName,append=TRUE)
    
    seqs <- read.dna(paste(dataPath,dataName,".fas",sep=""), as.matrix=FALSE, format="fasta")
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
    write(tempM, file=paste(dataPath,"human_cov2020_start_trim_with_ref_mafft.fas",sep=""))
    
    write("# new psuedo alignment file= human_cov2020_start_trim_with_ref_mafft.fas", file=logName, append=TRUE)
    
    endTime <- Sys.time()
    write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
    write(paste("# process end",endTime), file=logName,append=TRUE)
    write(paste("# END"),file=logName,append=TRUE)
    
    if (echo) {
      print("Finshed")
    }
   
  }
  return( logName )
}



##########################################################################################

GISAID_tidy_1 <- function(dataDate="2020_apr_8",
                        rootPath="E://data//Coronavirus_2020//",
                        dataPath=paste(rootPath,dataDate,"_gisaid//",sep=""),
                        refName="consensus",
                        maxFractGap=0.1,
                        maxDivergence=-1,
                        echo=TRUE
                        ) {
  
  # old refName="COVID-19_human_373_coding_consensus"
  
  if (maxDivergence<=0) {
    dataMon <- getEl(dataDate,ind=2,sep="_")
    namedMonths <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    mm <- match(tolower(dataMon), tolower(namedMonths))
    maxDivergence = 10*mm
  }
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"GISAID_tidy_1_log.txt",sep="")
  write("# Log file for GISAID_tidy_1",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)

  dataName <- "human_cov2020_start_trim_with_ref_mafft.fas"
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
  
  diff_to_ref <- seqs != t(matrix( rep(seqs[ref_i,],nseqs), length(refseq), nseqs))
  ndiff       <- apply(diff_to_ref, 1, sum)
  ndiff_valid <- apply(diff_to_ref[,which(refseq!="-")], 1, sum)
  
  mostly_gaps <- which(site_res$numValid<=maxFractGap*nseqs)
  gap_causing <- c()
  for (i in 1:length(mostly_gaps)) {
    inds <- which(as.character(seqs[,mostly_gaps[i]]) != "-")
    gap_causing <- c(gap_causing,inds)
  }
  top_gap_causing_tbl <- sort(table(gap_causing),decreasing=TRUE)
  seqs_to_remove      <- as.integer(rownames(top_gap_causing_tbl)[which(top_gap_causing_tbl>3)])
  seqs_to_remove      <- intersect(seqs_to_remove, which(ndiff_valid > maxDivergence))
  
  if (echo) {
    print("Remove sequences")
    print(paste(seqs_to_remove,taxa[seqs_to_remove]))
  }
  errorLog <- paste(dataPath,"gap_causing_fail.txt",sep="")
  write("# sequences which caused exessive gaps and too different to reference",file=errorLog,append=FALSE)
  write(paste(seqs_to_remove,taxa[seqs_to_remove],sep=","), file=errorLog, append=TRUE)

  write.dna(seqs[seqs_to_remove,], file=paste(dataPath,"gap_causing_fail.fas",sep=""),
            format="fasta", nbcol=-1, colsep="")
  
  seqs_to_keep <- setdiff(1:nseqs,seqs_to_remove)
  seqs <- seqs[seqs_to_keep,]
  
  write(paste("number sequences removed=",length(seqs_to_remove)), file=logName,append=TRUE)
  write(paste("number sequences keep=",length(seqs_to_keep)),file=logName,append=TRUE)

  dataName <- "human_cov2020_with_ref_al0.fas"
  fname    <- paste(dataPath,dataName,sep="")
  
  if (echo) {
    print(paste("number sequences removed=",length(seqs_to_remove)))
    print(paste("number sequences keep=",length(seqs_to_keep)))
    print(paste("Writing to",fname))
  }
  
  write.dna(seqs, file=fname, format="fasta", nbcol=-1, colsep="")
  write(paste("alignment 0=",dataName), file=logName,append=TRUE)
  
  if (echo) {
    print("Part B - Degapping begin")
  }
  write("# Begin de-gapping", file=logName, append=TRUE)
  
  # give it a few s to make sure that file is written properly
  Sys.sleep(5)
  
  seqs     <- read.dna(file=fname,format="fasta",as.matrix=TRUE)
  site_res <- siteSpectra(seqs)
  all_gaps <- which(site_res$numValid <= 1)
  ok_pos   <- which(site_res$numValid > 1)
  seqs <- seqs[,ok_pos]
  
  dataName <- "human_cov2020_with_ref_al1.fas"
  fname    <- paste(dataPath,dataName,sep="")
  
  if (echo) {
    print(paste("number gap sites removed=",length(all_gaps)))
    print(paste("number sites keep=",length(ok_pos)))
    print(paste("Writing to",fname))
  }
  
  write.dna(seqs, file=fname, format="fasta", nbcol=-1, colsep="")
  write(paste("alignment 1=",dataName), file=logName,append=TRUE)
  

  if (echo) {
    print("Part C - trim end and remove reference")
  }
  write("# trim end and remove reference")
  
  # give it a few s to make sure that file is written properly
  Sys.sleep(5)
  
  seqs     <- read.dna(file=fname,format="fasta",as.matrix=TRUE)
  taxa     <- attributes(seqs)$dimnames[[1]]
  ref_i    <- grep(refName,taxa)
  refseq   <- unlist(as.character(seqs[ref_i,]))
  pos <- gregexpr("tag",paste(refseq,collapse=""))
  endpos <- pos[[1]][length(pos[[1]])]
  seqs <- seqs[,1:(endpos+2)]
  seqs <- seqs[2:length(seqs[,1]),]
  
  dataName <- "human_cov2020_al2.fas"
  fname    <- paste(dataPath,dataName,sep="")
  
  if (echo) {
    print("Trimmed to end of ref")
    print("Removed ref")
    print(paste("number of sequences=",length(seqs[,1])))
    print(paste("Writing to",fname))
  }
  
  write.dna(seqs, file=fname, format="fasta", nbcol=-1, colsep="")
  write(paste("alignment 2=",dataName), file=logName,append=TRUE)
  write(paste("number sequences=",length(seqs[,1])),file=logName,append=TRUE)
  
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
  
}

GISAID_tidy_2 <- function(dataDate="2020_apr_8",
                          rootPath="E://data//Coronavirus_2020//",
                          dataPath=paste(rootPath,dataDate,"_gisaid//",sep=""),
                          echo=TRUE
                          ) {
  

  
  startTime <- Sys.time()
  logName <- paste(dataPath,"GISAID_tidy_2_log.txt",sep="")
  write("# Log file for GISAID_tidy_2",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  
  dataName <- "human_cov2020_al2.fas"
  fname    <- paste(dataPath,dataName,sep="")
  write(paste("dataName=",fname), file=logName, append=TRUE)
  
  if (echo) {
    print(paste("Reading sequences",fname))
  }
  
  seqs <- read.dna(fname, format="fasta", as.matrix=FALSE)
  
  write("# making initial NJ tree with TN93, gamma=0.5, pairwise.deletion=TRUE", file=logName, append=TRUE)
  write(paste("# tree start",Sys.time()), file=logName, append=TRUE)
  
  if (echo) {
    print("Making initial NJ tree with TN93, gamma=0.5, pairwide.deletion=TRUE")
  }
  
  dd <- dist.dna(seqs, model="TN93", gamma=0.5, pairwise.deletion=TRUE)
  save(dd, file=paste(dataPath,"human_cov2020_al2_genetic_dist.dd.Rdata",sep=""))
  tr <- nj(dd)
  tr <- ladderize(tr)
  save(tr, file=paste(dataPath,"human_cov2020_al2_tree.tr.Rdata",sep=""))
  
  trName <- paste(dataPath,"human_cov2020_al_ape_TN93_nj.nwk",sep="")
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


GISAID_tidy_3 <- function(dataDate="2020_apr_8",
                          rootPath="E://data//Coronavirus_2020//",
                          dataPath=paste(rootPath,dataDate,"_gisaid//",sep=""),
                          longthres=5e-4,
                          echo=TRUE
                          ) {
  
  
  
  startTime <- Sys.time()
  logName <- paste(dataPath,"GISAID_tidy_3_log.txt",sep="")
  write("# Log file for GISAID_tidy_3",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  
  #All strains list file= all_strains_list.csv
  all_strains_list <- read.csv(paste(dataPath,"all_strains_list.csv",sep=""))
  
  dataName <- "human_cov2020_al2"
  write(paste("dataName=",dataName),file=logName,append=TRUE)
  write(paste("longthres=",longthres),file=logName,append=TRUE)
  
  write("# loading NJ tree", file=logName,append=TRUE)
  if (echo) {
    print("loading NJ tree")
  }
  load(paste(dataPath,dataName,"_tree.tr.Rdata",sep=""))
  
  write("# rooting NJ tree on oldest", file=logName,append=TRUE)
  if (echo) {
    print("rooting NJ tree oldest")
  }
  taxa <- as.matrix(tr$tip.label)
  taxa_accn<- apply(taxa, 1, getEl, ind=1, sep="\\|")
  decDate <- as.numeric(apply(taxa, 1, getEl, ind=1, fromEnd=TRUE, sep="\\|"))
  oldest  <- which.min(decDate)
  tr      <- root(tr, oldest)
  
  # update all strains list with all the sequences in the tree before drop tip
  all_strains_list$AlignOK <- 0
  ainds <- match(taxa_accn, all_strains_list$Accn)
  all( taxa_accn==all_strains_list$Accn[ainds])
  all_strains_list$AlignOK[ainds] <- 1
  
  ntips <- length(taxa)
  #longthres <- 5e-4
  longbranches <- which((tr$edge.length >= longthres) & (tr$edge[,2] <= ntips))
  longtips <- tr$edge[longbranches,2]
  
  errorLog <- paste(dataPath,"longbranch_tips_excluded.txt",sep="")
  write("# sequences excluded because tip branch length too long", file=errorLog, append=FALSE)
  write(paste(longtips,taxa[longtips],sep=","),file=errorLog,append=TRUE)
  
  write(paste("ntips=",ntips),file=logName,append=TRUE)
  write(paste("longbranches=",length(longbranches)),file=logName,append=TRUE)
  write(paste("# dropping tips of long branches"),file=logName,append=TRUE)
  write(paste("longtipErroLog= longbranch_tips_excluded.txt"), file=logName,append=TRUE)
  tr <- drop.tip(tr, longtips)
  
  taxa <- as.matrix(tr$tip.label)
  taxa_accn <- apply(taxa, 1, getEl, ind=1, sep="\\|")
  decDate <- as.numeric(apply(taxa, 1, getEl, ind=1, fromEnd=TRUE, sep="\\|"))
  oldest  <- which.min(decDate)
  tr      <- multi2di(tr)
  tr      <- ladderize(tr)
  tr      <- root(tr, oldest)
  tr      <- multi2di(tr)
  
  all_strains_list$TempestOK <- 0
  ainds <- match(taxa_accn, all_strains_list$Accn)
  all( taxa_accn==all_strains_list$Accn[ainds])
  all_strains_list$TempestOK[ainds] <- 1
  
  load(paste(dataPath,"human_cov2020_nthres300.info.Rdata",sep=""))
  tinds <- match(taxa, info[,1])
  all(taxa == info[tinds,1])
  info <- info[tinds,]
  
  tr$continent <- info$continent
  ccols <- colourize_BEAST_cols(info$continent, transparency = 0.7)
  ccols2<- colourize_BEAST_cols(info$continent, transparency = 0.2)
  
  save(tr, file=paste(dataPath,dataName,"_tree_droptips_with_continent.tr.Rdata",sep=""))
  write.tree(tr, paste(dataPath,dataName,"_tree_droptips.nwk",sep=""))
  
  save(info, file=paste(dataPath,"human_cov2020_droptips.info.Rdata",sep=""))
  write.table(info, file=paste(dataPath,"human_cov2020_droptips_info.txt",sep=""),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  
  write.table(all_strains_list, file=paste(dataPath,"all_strains_list.csv",sep=""), 
              col.names=TRUE,row.names=FALSE,quote=FALSE,sep=",")
  
  imageName <- paste(dataPath,dataName,"_tree_droptips_with_continent.png",sep="")
  png(file=imageName, height=8*300, width=8*300, res=300)
  
  plot(tr, show.tip=FALSE)
  tiplabels(pch=21, col=ccols$statecols, bg=ccols$statecols)
  ctbl   <- table(info$continent)
  legtxt <- paste(rownames(ctbl)," (",ctbl,")",sep="")
  legend("bottomleft", legtxt, pch=21, pt.bg=ccols2$ucols, col=ccols$ucols, bty="n")
  
  dev.off()
  
  write(paste("trName= ",dataName,"_tree_droptips_with_continent.tr.Rdata",sep=""),
        file=logName,append=TRUE)
  write(paste("trName= ",dataName,"_tree_droptips.nwk",sep=""),
        file=logName,append=TRUE)
  write(paste("treeImageName= ",dataName,"_tree_droptips_with_continent.png",sep=""),
        file=logName,append=TRUE)
  write("infoName= human_cov2020_droptips.info.Rdata",file=logName,append=TRUE)
  write("infoName= human_cov2020_droptips_info.txt",file=logName,append=TRUE)
  write("# all_strains_list.csv updated", file=logName,append=TRUE)
  
  if (echo) {
    print("Extracting good sequences")
  }
  
  write("# writing good sequences to file", file=logName, append=TRUE)
  seqs <- read.dna(paste(dataPath,dataName,".fas",sep=""), format="fasta", as.matrix=FALSE)
  taxa <- attributes(seqs)$names
  sinds<- match(info[,1],taxa)
  all( info[,1]==taxa[sinds])
  write(paste("number sequences=",length(sinds)), file=logName, append=TRUE)
  write.dna(seqs[sinds], file=paste(dataPath,"human_cov2020_al3.fas",sep=""),
            format="fasta", nbcol=-1, colsep="")
  
  write("dataName= human_cov2020_al3", file=logName, append=TRUE)
  
  
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
  
}
 
# 27 April 2020
# run this to summarise data after cleaning
# writes file GISAID_summary_log.txt and others
GISAID_summary <- function(dataDate="2020_apr_8",
                              rootPath="E://data//Coronavirus_2020//",
                              dataPath=paste(rootPath,dataDate,"_gisaid//",sep=""),
                              echo=TRUE,
                              selcol=6
                              ) {

  startTime <- Sys.time()
  logName <- paste(dataPath,"GISAID_summary_log.txt",sep="")
  write("# Log file for GISAID_summary",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  
  #All strains list file= all_strains_list.csv
  all_strains_list <- read.csv(paste(dataPath,"all_strains_list.csv",sep=""))
  
  #all expanded strains info
  load(paste(dataPath,"all_expanded_strain_info.info.Rdata",sep=""))
  
  minds <- match(all_strains_list$Accn, info$Accession.ID)
  all(all_strains_list$Accn==info$Accession.ID[minds])
  #jj <- which(all_strains_list$HumanFullDate==1)
  
  # dates processing
  jj        <- which(all_strains_list[,2]==1)
  ok_inds   <- minds[jj]
  max_ii    <- ok_inds[which.max(info$decDate[ok_inds])]
  min_ii    <- ok_inds[which.min(info$decDate[ok_inds])]
  min_date  <- as.integer(strsplit(paste(info$Collection.date[min_ii]),"-")[[1]])
  max_date  <- as.integer(strsplit(paste(info$Collection.date[max_ii]),"-")[[1]])
  
  if (min(decDate[ok_inds]) > 2020) {
    min_date<- c(2020,1,1)
  }
  
  dates_res <- get_udates(min_date=min_date,max_date=max_date)
 
  #selcol <- 6
  write(paste("# from all_strains_list.csv file, using col=",selcol), file=logName,append=TRUE)
  write(paste("Good colname=",colnames(all_strains_list)[selcol]), file=logName,append=TRUE)
  
  #get continent table
  jj   <- which(info$host[minds]=="Human")
  ctbl <- table(all_strains_list[jj,selcol],info$continent[minds[jj]])
  rownames(ctbl)<-paste("Good",rownames(ctbl),sep="=")
  total<- table(info$continent[minds[jj]])
  ctbl <- t(rbind(ctbl,total))
  write.csv(ctbl, file=paste(dataPath,"GISAID_Continent_summary_tbl.csv",sep=""))
  
  cont_dates    <- factor( paste(info$Collection.date[minds[jj]]), levels=dates_res$udate)
  cont_date_tbl <- table(cont_dates,info$continent[minds[jj]])
  write.csv(cont_date_tbl, file=paste(dataPath,"GISAID_all_Continent_summary_date_table.csv",sep=""))
  
  jj   <- which(info$host[minds]=="Human" & all_strains_list[,selcol]==1)
  cont_dates    <- factor( paste(info$Collection.date[minds[jj]]), levels=dates_res$udate)
  cont_date_tbl <- table(cont_dates,info$continent[minds[jj]])
  write.csv(cont_date_tbl, file=paste(dataPath,"GISAID_good_Continent_summary_date_table.csv",sep=""))
  
  #get numbers of european sequences
  jj     <- which(info$continent[minds]=="Europe" & info$host[minds]=="Human")
  eu_tbl <- table(all_strains_list[jj,selcol],paste(info$country[jj]))
  rownames(eu_tbl)<-paste("Good",rownames(eu_tbl),sep="=")
  total  <- table(paste(info$country[minds[jj]]))
  eu_tbl <- t(rbind(eu_tbl,total))
  write.csv(eu_tbl, file=paste(dataPath,"GISAID_Europe_summary_tbl.csv",sep=""))
  
  #get numbers of uk sequences, by England, Scotland, Wales, Northern Ireland
  jj     <- which(info$country[minds]=="United Kingdom" & info$host[minds]=="Human")
  uk_tbl <- table(all_strains_list[jj,selcol],paste(info$state[minds[jj]]))
  rownames(uk_tbl)<-paste("Good",rownames(uk_tbl),sep="=")
  total  <- table(paste(info$state[minds[jj]]))
  uk_tbl <- t(rbind(uk_tbl,total))
  write.csv(uk_tbl, file=paste(dataPath,"GISAID_UK_summary_tbl.csv",sep=""))
  
  # all the UK sequnences not just the good ones
  u_uk_states <- c("England","Wales","Scotland","Northern Ireland")
  uk_cols   <- get_BEAST_cols(length(u_uk_states), bright=0.8, sat=0.7)
  
  jj        <- which(info$country[minds]=="United Kingdom" & info$host[minds]=="Human")
  uk_dates  <- factor( paste(info$Collection.date[minds[jj]]), levels=dates_res$udate)
  uk_states <- factor( paste(info$state[minds[jj]]), levels=c(u_uk_states) )
  
  imageName <- paste(dataPath,"GISAID_all_UK_summary_barplot_6x12.png",sep="")
  png(file=imageName, height=6*300, width=12*300, res=300)
  uk_date_tbl <- table(uk_dates,uk_states)
  barplot(t(uk_date_tbl),col=uk_cols,names=dates_res$xlabs,xlab="Collection date",ylab="Number of sequences")
  legend("topleft",paste(colnames(uk_date_tbl)," (",apply(uk_date_tbl, 2, sum),")",sep=""),
         pch=22,pt.bg=uk_cols,bty="n")
  title(paste("Number of UK Sequences",toupper(gsub("_"," ",dataDate))))
  dev.off()
  write.csv(uk_date_tbl, file=paste(dataPath,"GISAID_all_UK_summary_date_table.csv",sep=""))
  
  # the good uk sequences only
  jj        <- which(info$country[minds]=="United Kingdom" & info$host[minds]=="Human" &
                       all_strains_list[,selcol]==1)
  uk_dates  <- factor( paste(info$Collection.date[minds[jj]]), levels=dates_res$udate)
  uk_states <- factor( paste(info$state[minds[jj]]), levels=c(u_uk_states) )
  
  imageName <- paste(dataPath,"GISAID_good_UK_summary_barplot_6x12.png",sep="")
  png(file=imageName, height=6*300, width=12*300, res=300)
  uk_date_tbl <- table(uk_dates,uk_states)
  barplot(t(uk_date_tbl),col=uk_cols,names=dates_res$xlabs,xlab="Collection date",ylab="Number of sequences")
  legend("topleft",paste(colnames(uk_date_tbl)," (",apply(uk_date_tbl, 2, sum),")",sep=""),
         pch=22,pt.bg=uk_cols,bty="n")
  title(paste("Number of UK Sequences",toupper(gsub("_"," ",dataDate))))
  dev.off()
  write.csv(uk_date_tbl, file=paste(dataPath,"GISAID_good_UK_summary_date_table.csv",sep=""))
  
  write("##############",file=logName,append=TRUE)
  write("Continent,Good=0,Good=1,Total", file=logName,append=TRUE)
  write.table(ctbl, file=logName, append=TRUE, sep=",", row.names=TRUE, col.names=FALSE)

  
  write("##############",file=logName,append=TRUE)
  write("Europe,Good=0,Good=1,Total", file=logName,append=TRUE)
  write.table(eu_tbl, file=logName, append=TRUE, sep=",", row.names=TRUE, col.names=FALSE)
  
  
  write("##############",file=logName,append=TRUE)
  write("UK,Good=0,Good=1,Total", file=logName,append=TRUE)
  write.table(uk_tbl, file=logName, append=TRUE, sep=",", row.names=TRUE, col.names=FALSE)
  
  
  write("##############",file=logName,append=TRUE)
  endTime <- Sys.time()
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  if (echo) {
    print("Finshed")
  }
  return( logName )
}


##################################################################################
# temp example run
doPipeline <- FALSE
if (doPipeline) {
  #dataDate <- "2020_apr_20"
  #dataDate <- "2020_may_7"  # re-run on 8 May with new calc decimal dates
  #dataDate <- "2020_may_15"
  dataDate <- "2020_may_27"
  GISAID_clean_1(dataDate=dataDate)
  
  GISAID_clean_2(dataDate=dataDate)
  GISAID_clean_3(dataDate=dataDate)
  
  stop()
  
  GISAID_align(dataDate=dataDate)
  
  stop()
  
  GISAID_tidy_1(dataDate=dataDate)
  GISAID_tidy_2(dataDate=dataDate)
  GISAID_tidy_3(dataDate=dataDate)
  
  GISAID_summary(dataDate=dataDate)
  
}

# example with MSA - ok but downloaded alignment is not very good anyway
doMSA <- FALSE
if (doMSA) {
  dataDate <-"2020_apr_28"
  dataPath <- "E://data//Coronavirus_2020//2020_apr_28_gisaid//msa_2020-04-28//"
  dataName <- "msa_2020-04-28.fasta"
  GISAID_clean_1(dataDate=dataDate,dataPath=dataPath,dataName=dataName)
  GISAID_clean_2(dataDate=dataDate,dataPath=dataPath,dataName=dataName)
  GISAID_clean_3(dataDate=dataDate,dataPath=dataPath)
  
  GISAID_align(dataDate=dataDate,dataPath=dataPath,doConsensus = TRUE)
  
  
  GISAID_tidy_1(dataDate=dataDate, dataPath=dataPath)
  #GISAID_tidy_2(dataDate=dataDate, dataPath=dataPath)
  #GISAID_tidy_3(dataDate=dataDate, dataPath=dataPath)
  
  #GISAID_summary(dataDate=dataDate, dataPath=dataPath)
  
}

# example with UK sequences only
doSingleCountry <- FALSE
if (doSingleCountry) {
  #dataDate <-"2020_apr_29"
  #dataDate <- "2020_may_4"
  #dataDate <- "2020_may_7"
  dataDate <- "2020_may_15"
  country  <- "United_Kingdom"
  rootPath <- "E://data//Coronavirus_2020//"
  dataPath <- paste(rootPath,dataDate,"_gisaid//",country,"//",sep="")
  fnames   <- dir(dataPath)
  inds     <- grep("gisaid_hcov-19_2020_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]\\.fasta",fnames)
  fnames   <- fnames[inds]
  dataName <- fnames[grep(country,fnames)]
  
  #GISAID_clean_1(dataDate=dataDate,dataPath=dataPath,dataName=dataName,use_sub_parts=FALSE,incMetaData=FALSE)
  #GISAID_clean_2(dataDate=dataDate,dataPath=dataPath,dataName=dataName)
  
  # 15 may 2020 - had to split download England, Scotland, Wales, Northern Ireland because > 10,000
  GISAID_clean_1(dataDate=dataDate,dataPath=dataPath,dataName=dataName,
                 use_sub_parts=TRUE,n2="_gisaid_hcov-19_",incMetaData=FALSE)
  GISAID_clean_2(dataDate=dataDate,dataPath=dataPath,dataName=dataName,use_sub_parts=TRUE)
  
  GISAID_clean_3(dataDate=dataDate,dataPath=dataPath)
  
  #do mafft on command line outside R
  #GISAID_align(dataDate=dataDate,dataPath=dataPath,doConsensus=FALSE)
  
  
  GISAID_tidy_1(dataDate=dataDate, dataPath=dataPath)
  GISAID_tidy_2(dataDate=dataDate, dataPath=dataPath)
  GISAID_tidy_3(dataDate=dataDate, dataPath=dataPath)
  
  GISAID_summary(dataDate=dataDate, dataPath=dataPath)
}
