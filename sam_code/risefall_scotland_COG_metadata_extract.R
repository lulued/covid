# script to process COG metadata - extracting Scotland
# S J Lycett
# 19 May 2020
# 20 May 2020
# 24 May 2020

# cog metadata now includes admin 2
# means that can get places within Scotland

#  path     <- "Documents//data//Coronavirus_2020/COG_downloads//2020-05-15//microreact//"
#  metadata <- read.csv(paste(path,"cog_global_2020-05-15_metadata_private.csv",sep=""))

library(ape)

# packages for tree dating - should all be under treedater, but might need to separarely install
# this is Erik Volz tree dating
library(quadprog)
library(limSolve)
library(treedater)

Rpath <- "~/Documents//GitLab//RiseFall-Scottish-COVID19//coronavirus_2020_r//"
#Rpath <- ""
source(paste(Rpath,"getEl.R",sep=""))
source(paste(Rpath,"calcDecimalDate.R",sep=""))
source(paste(Rpath,"get_BEAST_cols.R",sep=""))
source(paste(Rpath,"get_udates.R",sep=""))
source(paste(Rpath,"birdSpecies_and_continents.R",sep=""))
source(paste(Rpath,"birdSpecies_and_continents.R",sep=""))
source(paste(Rpath,"read_MCC_tree.R",sep=""))


####################################################################################
# GLOBAL SETUP

#dataDate <- "2020-05-15"
dataDate <- "2020-05-22"
rootPath <- "~//Documents//data//Coronavirus_2020/COG_downloads//"
dataPath <- paste(rootPath,dataDate,"_processed//",sep="")
dataName <- "cog_global_scotland"

# tree dater global parameters
clock           <- "strict"
etol            <- 1e-12
seqlen          <- 30000
omega0          <- 1e-3
minblen         <- etol/omega0
meanRateLimits  <- c(1e-5,5e-3)

# run chunk parts
doSetUp       <- TRUE
doGlobalPlots <- FALSE
doUKPlots     <- FALSE
doScotPlots   <- FALSE

doPlotGlobalTree<- FALSE

doScotTree    <- FALSE
doScotRegions <- FALSE
doScotLineages<- FALSE


####################################################################################
# CUSTOM FUNCTIONS
# define custom plot functions
do_continent_pie <- function(cont_tbl=cont_tbl, ccols=ccols, 
                             ctxt=ctxt, dataDate=dataDate, show.legend=TRUE) {
  op <- par(mar=c(0,0,2,0))
  pie(cont_tbl, col=ccols$ucols)
  title(paste("Number of sequences",dataDate))
  if (show.legend) {
    if (length(ctxt)<10) {
      legend("topleft",ctxt,pch=22,bty="n",pt.bg=ccols$ucols)
    } else {
      nn <- length(ctxt)
      nn2<- ceiling(nn/2)
      legend("topleft",ctxt[1:nn2],pch=22,bty="n",pt.bg=ccols$ucols[1:nn2])
      legend("topright",ctxt[(nn2+1):nn],pch=22,bty="n",pt.bg=ccols$ucols[(nn2+1):nn])
    }
  }
  par(op)
}

do_continent_bars <- function(continent, ctxt=ctxt, ccols=ccols, dates=dates, 
                              dataDate=dataDate, dates_res=dates_res) {
  udate    <- dates_res$udate
  seldates <- dates_res$seldates
  seldates2<- dates_res$seldates2
  xlabs    <- dates_res$xlabs
  
  barplot(  table( continent, factor(dates,levels=udate ) ), 
            names=xlabs, col=ccols$ucols, xlab="Sample date", ylab="Number of sequences" )
  legend("topleft",ctxt,pch=22,bty="n",pt.bg=ccols$ucols)
  title(paste("Number of sequences",dataDate))
}

do_uk_bars <- function(country_state, uk_inds=uk_inds, uk_txt=uk_txt, uk_cols=uk_cols, dates=dates, 
                       dataDate=dataDate, dates_res=dates_res, show.legend=TRUE) {
  udate    <- dates_res$udate
  seldates <- dates_res$seldates
  seldates2<- dates_res$seldates2
  xlabs    <- dates_res$xlabs
  
  barplot(  table( country_state[uk_inds], factor(dates[uk_inds],levels=udate ) ), 
            names=xlabs, col=uk_cols$ucols, xlab="Sample date", ylab="Number of sequences" )
  
  if (show.legend) {
    if (length(uk_txt)<10) {
      legend("topleft",uk_txt,pch=22,bty="n",pt.bg=uk_cols$ucols)
    } else {
      nn <- length(uk_txt)
      nn2<- ceiling(nn/2)
      legend("topleft",uk_txt[1:nn2],pch=22,bty="n",pt.bg=uk_cols$ucols)
      legend("topright",uk_txt[(nn2+1):nn],pch=22,bty="n",pt.bg=uk_cols$ucols)
    }
  }
  title(paste("Number of sequences",dataDate))
}

top_pie_lineages <- function(ltbl, lincols=lincols, show.legend=TRUE, legpos="topleft", titleTxt="") {
  sorted_ltbl <- sort(ltbl,decreasing=TRUE, index.return=TRUE)
  incs <- which(sorted_ltbl$x>0)
  sorted_labels <- attributes(sorted_ltbl$x)$names
  pie(sorted_ltbl$x[incs],col=lincols$ucols[sorted_ltbl$ix[incs]])
  
  if (show.legend) {
    legtxt <- paste(sorted_labels[incs]," (",sorted_ltbl$x[incs],")",sep="")
    legend(legpos,legtxt,pch=22,pt.bg=lincols$ucols[sorted_ltbl$ix[incs]],bty="n")
  }
  if (titleTxt!="") {
    title(titleTxt)
  }
}

get_metadata <- function(metaPath=metaPath, metaName=metaName) {
  metadata <- read.csv(paste(metaPath,metaName,sep=""))
  return( metadata )
}

repair_meta <- function(rootPath="~/Documents//data//Coronavirus_2020/COG_downloads//", dataDate1="2020-05-15", dataDate2="2020-05-22", 
                          path2="microreact", ext="_metadata_private.csv", ext2="_metadata_private_repaired.csv") {
  
  metaPath1 <- paste(rootPath,dataDate1,"//",path2,"//",sep="")
  metaName1 <- paste("cog_global_",dataDate1,ext,sep="")
  meta1 <- get_metadata(metaPath1, metaName1)
  
  metaPath2 <- paste(rootPath,dataDate2,"//",path2,"//",sep="")
  metaName2 <- paste("cog_global_",dataDate2,ext,sep="")
  meta2 <- get_metadata(metaPath2, metaName2)
  
  metaName3 <- paste("cog_global_",dataDate2,ext2,sep="")
  
  minds12 <- match(meta1[,1],meta2[,1])
  jj      <- which(is.finite(minds12))
  ncols   <- length(meta1[1,])
  ok      <- matrix(0, ncols, 2)
  rownames(ok) <- colnames(meta1)
  colnames(ok) <- c("Num match values","Num mismatch")
  
  print(paste("Num meta1=",length(meta1[,1])))
  print(paste("Num meta2=",length(meta2[,1])))
  print(paste("Num match names=",length(jj)))

  
  for (k in 1:ncols) {
    ok[k,1] <-length(which(paste(meta1[jj,k])==paste(meta2[minds12[jj],k])))
    ok[k,2] <-length(which(paste(meta1[jj,k])!=paste(meta2[minds12[jj],k])))
  }
  
  print(ok)
  
  minds21 <- match(meta2[,1],meta1[,1])
  kk      <- which(is.finite(minds21))
  all( paste(meta2[kk,1])==paste(meta1[minds21[kk],1]))
  kk_i    <- kk[which(meta2$lineage!="")]
  
  kk_l    <- kk[ which(meta2$lineage[kk]=="")]
  m2_lineage <- paste(meta2$lineage)
  m2_lineage[kk_l] <- paste(meta1$lineage[minds21[kk_l]])
  print("Which lineage = (none)")
  print(which(m2_lineage==""))
  
  meta2$lineage <- factor(m2_lineage)
  write.csv(meta2, file=paste(metaPath2,metaName3,sep=""))
  print(paste("Repaired lineages written to ",metaPath2,metaName3,sep=""))
}

####################################################################################
# MAIN SCRIPT
####################################################################################

# SPECIFIC SET UP and first process of metadata

if (doSetUp) {
  # set up paths
  #rootPath <- "~/Documents//data//Coronavirus_2020/COG_downloads//"
  path2    <- "microreact"
  #ext      <-"_metadata_private.csv"
  ext      <-"_metadata_private_repaired.csv"
  metaPath <- paste(rootPath,dataDate,"//",path2,"//",sep="")
  metaName <- paste("cog_global_",dataDate,ext,sep="")
  treeName <- paste("cog_global_",dataDate,"_tree.newick",sep="")
  
  metadata <- get_metadata(metaPath=metaPath,metaName=metaName)

  # configure the metadata information
  # basic continents
  region<- unGeoScheme(metadata$country)
  continent <- region
  # re-configure so is like the GISAID
  continent[grep("Oceania", continent)] <- "Oceania"
  continent[grep("Asia", continent)]    <- "Asia"
  continent[grep("Europe", continent)]  <- "Europe"
  continent[grep("Africa", continent)]  <- "Africa"
  continent[grep("NorthernAmerica", continent)] <- "North America"
  continent[grep("SouthernAmerica", continent)] <- "South America"
  continent[grep("CentralAmerica", continent)] <- "Central America"
  ccols <- colourize_BEAST_cols(continent)
  cont_tbl <- table(continent)
  ctxt  <- paste(rownames(cont_tbl)," (",cont_tbl,")",sep="")

  # country_state and UK
  country       <- metadata$country
  country_state <- paste(metadata$adm1)
  u_uk_states   <- c("England","Wales","Scotland","Northern Ireland")
  country_state[grep("UK-ENG", country_state)] <- "England"
  country_state[grep("UK-SCT", country_state)] <- "Scotland"
  country_state[grep("UK-WLS", country_state)] <- "Wales"
  country_state[grep("UK-NIR", country_state)] <- "Northern Ireland"
  country_state[which(country_state=="")] <- "-"
  country_state <- factor(country_state, levels=c(u_uk_states,"-"))
  uk_cols     <- colourize_BEAST_cols(country_state, ustates=u_uk_states, 
                                    ucols=c(get_BEAST_cols(4,bright=0.8, sat=0.7),"grey90"))
  uk_inds     <- which(country=="UK")
  uk_tbl <- table(country_state[uk_inds])
  uk_txt <- paste(rownames(uk_tbl)," (",uk_tbl,")",sep="")
  uk_txt <- uk_txt[1:4]

  # uk_county
  uk_county <- paste(metadata$adm2)
  uk_county[grep("UNKNOWN",uk_county)] <- "UNKNOWN"
  uk_county[which(uk_county=="")] <- "UNKNOWN"
  u_uk_county <- setdiff(sort(unique(uk_county)),"UNKNOWN")
  uk_county <- factor(uk_county, levels=c(u_uk_county,"UNKNOWN"))

  # scotland
  scot_inds <- which(country_state=="Scotland")
  u_scot    <- c(setdiff(sort(unique(uk_county[scot_inds])),"UNKNOWN"),"UNKNOWN")
  scot_cols <- colourize_BEAST_cols(uk_county, ustates=u_scot, 
                                  ucols=c(get_BEAST_cols(length(u_scot)-1,bright=0.8, sat=0.7),"grey90"))
  scot_tbl  <- table(factor(paste(uk_county[scot_inds]), levels=u_scot))
  scot_txt  <- paste(rownames(scot_tbl)," (",scot_tbl,")",sep="")

  # decimal dates
  # correction
  metadata$sample_date[which(metadata$sample_date=="2022-04-05")] <- "2020-04-05"
  decDates <- apply(as.matrix(metadata$sample_date), 1, calcDecimalDate_fromTxt, sep="-", dayFirst=FALSE)

  maxDate   <- max(decDates)
  max_ii    <- which.max(decDates)
  min_ii    <- which.min(decDates)
  min_date  <- as.integer(strsplit(paste(metadata$sample_date[min_ii]),"-")[[1]])
  if (min(decDates) > 2020) {
    min_date<- c(2020,1,1)
  }

  max_date  <- as.integer(strsplit(paste(metadata$sample_date[max_ii]),"-")[[1]])
  dates_res <- get_udates(min_date=min_date,max_date=max_date)
  
  newTaxa   <- paste(metadata$sequence_name, uk_county, format(decDates,digits=7), sep="|")
  
  metadata$decDates   <- decDates
  metadata$newTaxa    <- newTaxa
  metadata$continent  <- continent
  metadata$country_state <- country_state
  metadata$uk_county  <- uk_county
  
  #write.table(metadata, row.names=FALSE, col.names=TRUE, quote=FALSE,
  #            file=paste(metaPath,"cog_global_",dataDate,"_meta_extended.txt",sep=""), sep="\t")
  #save(metadata, file=paste(metaPath,"cog_global_",dataDate,"_meta_extended.metadata.Rdata",sep=""))
}
  

########################################################################
# START OF SUMMARY FIGURES

# GLOBAL
if (doGlobalPlots) {
  imageName <- paste(dataPath,"cog_global_continent_pie_6x12.png",sep="")
  png(file=imageName, height=6*300, width=12*300, res=300)
  do_continent_pie(cont_tbl=cont_tbl, ccols=ccols, ctxt=ctxt, dataDate=dataDate)
  dev.off()
  
  
  imageName <- paste(dataPath,"cog_global_continent_bars_6x12.png",sep="")
  png(file=imageName, height=6*300, width=12*300, res=300)
  do_continent_bars(continent=continent, ctxt=ctxt, ccols=ccols, 
                    dates=metadata$sample_date, dataDate=dataDate,
                    dates_res=dates_res)
  dev.off()
}

# UK
if (doUKPlots) {
  imageName <- paste(dataPath,"cog_global_uk_bars_6x12.png",sep="")
  png(file=imageName, height=6*300, width=12*300, res=300)
  do_uk_bars(country_state, uk_inds=uk_inds, uk_txt=uk_txt,
             uk_cols=uk_cols, dates=metadata$sample_date, 
             dataDate=dataDate, dates_res=dates_res)
  dev.off()
  
  # "UK Sequences by Region"
  write.csv( table(uk_county[uk_inds],factor(country_state[uk_inds],levels=u_uk_states)),
             file=paste(dataPath,"cog_global_uk_sequences_by_region.csv",sep=""))
}

# SCOTLAND
if (doScotPlots) {
  imageName <- paste(dataPath,"cog_global_scotland_pie_6x12.png",sep="")
  png(file=imageName, height=6*300, width=12*300, res=300)
  do_continent_pie(cont_tbl=scot_tbl, ccols=scot_cols, ctxt=scot_txt, dataDate=dataDate)
  dev.off()
  
  imageName <- paste(dataPath,"cog_global_scotland_pie_no_legend_6x12.png",sep="")
  png(file=imageName, height=6*300, width=12*300, res=300)
  do_continent_pie(cont_tbl=scot_tbl, ccols=scot_cols, ctxt=scot_txt, dataDate=dataDate, show.legend=FALSE)
  dev.off()
  
  imageName <- paste(dataPath,"cog_global_scotland_bars_6x12.png",sep="")
  png(file=imageName, height=6*300, width=12*300, res=300)
  do_uk_bars(uk_county, uk_inds=scot_inds, uk_txt=scot_txt,
             uk_cols=scot_cols, dates=metadata$sample_date, 
             dataDate=dataDate, dates_res=dates_res, show.legend=FALSE)
  dev.off()
  
  # global pie, but with England, Wales, Scotland, Northern Ireland as well
  cont2 <- paste(metadata$continent)
  cont2[uk_inds] <- paste(metadata$country_state[uk_inds])
  uc2    <- c(ccols$ustates, uk_cols$ustates)
  c2order<- c(1:4,8:11,5:7)
  uc2    <- uc2[c2order]
  cont2 <- factor(cont2, levels=uc2)
  cont2_tbl <- table(cont2)
  ccols2 <- colourize_BEAST_cols(cont2, ustates=uc2, ucols=c(ccols$ucols,uk_cols$ucols)[c2order])
  imageName <- paste(dataPath,"cog_global_and_uk_pie_6x8.png",sep="")
  png(file=imageName, height=6*300, width=8*300, res=300)
  do_continent_pie(cont2_tbl, dataDate=dataDate, show.legend=FALSE, ccols=ccols2 )
  dev.off()
  
  # "Scottish sequences by epi-week"
  write.csv( table(  factor(metadata$adm2[scot_inds], levels=u_scot),metadata$epi_week[scot_inds]),
             file=paste(dataPath,"cog_global_scotland_by_epiweek.csv",sep=""))
  
  # "Scottish sequences by lineage"
  write.csv( table(  factor(metadata$adm2[scot_inds], levels=u_scot),paste(metadata$lineage[scot_inds])),
             file=paste(dataPath,"cog_global_scotland_by_lineage.csv",sep=""))
}

# LINEAGES
if (doLineagePlots) {
  europe_inds           <- which(metadata$continent=="Europe")
  europe_not_uk_inds    <- setdiff(europe_inds,uk_inds)
  uk_not_scotland_inds  <- setdiff(uk_inds,scot_inds)
  
  lintbl <- table(metadata$lineage)
  lintbl <- rbind(lintbl,table(metadata$lineage[europe_not_uk_inds]))
  lintbl <- rbind(lintbl,table(metadata$lineage[uk_not_scotland_inds]))
  lintbl <- rbind(lintbl,table(metadata$lineage[scot_inds]))
  rownames(lintbl) <- c("Global","Rest of Europe","Rest of UK","Scotland")
  lincols <- colourize_BEAST_cols(metadata$lineage)
  totseqs <- apply(as.matrix(lintbl), 1, sum)
  fractLinTbl <- lintbl/totseqs
  
  # for lineages in scotland only
  sorted_ltbl <- sort(lintbl[4,],decreasing=TRUE, index.return=TRUE)
  incs <- sorted_ltbl$ix[which(sorted_ltbl$x>0)]
  
  for (j in 1:4) {
    piePath <- paste(dataPath,"scotland_region_only//",dataName,"_",rownames(lintbl)[j],"_lineage_pie.png",sep="")
    png(file=piePath,height=6*150,width=6*150,res=150)
    top_pie_lineages(lintbl[j,], show.legend=FALSE, titleTxt=rownames(lintbl)[j], lincols=lincols)
    dev.off()
  }
  piePath2 <- paste(dataPath,"scotland_region_only//",dataName,"_lineage_pies.csv",sep="")
  write.csv( lintbl, file=piePath2)
  
  
  pieImPath <- paste(dataPath,"scotland_region_only//",dataName,"_global_lineages_x4_pie.png",sep="")
  png(file=pieImPath, height=4*150, width=12*150, res=150)
  op <- par(mfrow=c(1,4), mar=c(0.5,0.5,0.5,0.5)) 
  for (j in 1:4) {
    top_pie_lineages(lintbl[j,], show.legend=FALSE, titleTxt="", lincols=lincols)
    text(0,-1.2,rownames(lintbl)[j])
  }
  par(op)
  dev.off()
  
  print("Lineages in Scotland and not in Rest of UK")
  print(which(lintbl[4,]>0 & lintbl[3,]<=0))
  
  for (j in 1:length(u_scot)) {
    focusRegion <- u_scot[j]
    
    region_inds <- which(metadata$uk_county==focusRegion)
    rlintbl <- rbind(lintbl,table(metadata$lineage[region_inds]))
    rownames(rlintbl) <- c("Global","Rest of Europe","Rest of UK","Scotland",focusRegion)
    
    piePath <- paste(dataPath,"scotland_region_only//",dataName,"_",focusRegion,"_lineage_pie.png",sep="")
    png(file=piePath,height=6*150,width=6*150,res=150)
    top_pie_lineages(rlintbl[5,], show.legend=TRUE, titleTxt=rownames(rlintbl)[5], lincols=lincols)
    dev.off()
    
    piePath2 <- paste(dataPath,"scotland_region_only//",dataName,"_",focusRegion,"_lineage_pie.csv",sep="")
    write.csv( rlintbl[4:5,incs], file=piePath2)
  }
  
}


# END OF SUMMARY FIGURES
########################################################################

########################################################################
# START OF TREE PROCESSING

# PLOT GLOBAL TREE WITH LINEAGES
if (doPlotGlobalTree) {
  tr        <- read.tree(paste(metaPath,treeName,sep=""))
  ntips     <- length(tr$tip.label)
  bad_edges <- which((tr$edge.length >= 1/1000) & tr$edge[,2]<=ntips)
  todrop    <- tr$edge[bad_edges,2]
  tr        <- drop.tip(tr, todrop)

  einds     <- which(tr$edge.length < etol)
  tr$edge.length[einds] <- etol
  tr        <- multi2di(tr)
  tr        <- ladderize(tr)
  
  # repeat
  einds     <- which(tr$edge.length < etol)
  tr$edge.length[einds] <- etol
  tr        <- multi2di(tr)
  tr        <- ladderize(tr)
  tr        <- multi2di(tr)
  
  ntips     <- length(tr$tip.label)
  tinds <- match(tr$tip.label, metadata$sequence_name)
  all(tr$tip.label==metadata$sequence_name[tinds])
  tr$lineage <- metadata$lineage[tinds]
  tr$continent <- metadata$continent[tinds]
  l2 <- paste(tr$lineage)
  l2[which(l2=="")] <- "X"
  l2[grep("A.1",l2)] <- "A.1"
  l2[grep("B.1",l2)] <- "B.1"
  l2[grep("B.2",l2)] <- "B.2"
  tr$l2 <- l2
  l2tbl <- table(tr$l2)
  
  lintip_cols <- colourize_BEAST_cols(tr$lineage, transparency=0.5)
  lintxt      <- paste(colnames(lintbl)," (",lintbl[1,],")",sep="")
  #l2_ucols    <- lintip_cols$ucols[match(names(l2tbl), lintip_cols$ustates)]
  #l2_cols     <- lintip_cols$ucols[match(tr$l2, lintip_cols$ustates)]
  l2txt       <- paste( names(l2tbl)," (",l2tbl,")",sep="" )
  ace_l2 <- ace(l2, tr, type="discrete", model="ER")
  
  ancs   <- colnames(ace_l2$lik.anc)[apply(ace_l2$lik.anc, 1, which.max)]
  state  <- c(paste(l2),ancs)
  tr$state <- state
  estates  <- tr$state[tr$edge[,1]]
  l2cols   <- colourize_BEAST_cols(tr$state)
  ecols    <- l2cols$ucols[match(estates,l2cols$ustates)]
  
  trImPath <- paste(dataPath,"scotland_region_only//",dataName,"_global_lineage_tree_image.png",sep="")
  png(file=trImPath, height=12*150, width=8*150, res=150)
  plot(tr, show.tip.label=FALSE, edge.color=ecols)
  #tiplabels(pch=21,col=lintip_cols$statecols,bg=lintip_cols$statecols)
  #legend("bottomleft",lintxt,pch=21,col=lintip_cols$ucols,pt.bg=lintip_cols$ucols,bty="n")
  tiplabels(pch=21,col=l2cols$statecols,bg=l2cols$statecols)
  legend("bottomleft",l2txt,pch=21,col=l2cols$ucols,pt.bg=l2cols$ucols,bty="n")
  title("Global Tree with Lineages")
  add.scale.bar(2e-4,0,1e-4)
  dev.off()
  
  pieImPath <- paste(dataPath,"scotland_region_only//",dataName,"_global_lineage_pies_image.png",sep="")
  png(file=pieImPath, height=4*150, width=8*150, res=150)
  op <- par(mfrow=c(2,4), mar=c(0,0,0,0))
  for (j in 1:length(ccols$ustates)) {
    cinds <- which(tr$continent==ccols$ustates[j])
    pie(table(factor(tr$l2[cinds], levels=l2cols$ustates)),col=l2cols$ucols)
    text(0,-1,ccols$ustates[j])
  }
  par(op)
  dev.off()
  
}

# plot Scotland on tree
# also do quick discrete traits on tree so get crude estimate of imports
if (doImports) {
  
  startTime <- Sys.time()
  outPath   <- dataPath
  
  if (!doPlotGlobalTree) {
    tr        <- read.tree(paste(metaPath,treeName,sep=""))
    ntips     <- length(tr$tip.label)
    bad_edges <- which((tr$edge.length >= 1/1000) & tr$edge[,2]<=ntips)
    todrop    <- tr$edge[bad_edges,2]
    tr        <- drop.tip(tr, todrop)
    
    einds     <- which(tr$edge.length < etol)
    tr$edge.length[einds] <- etol
    tr        <- multi2di(tr)
    tr        <- ladderize(tr)
    
    # repeat
    einds     <- which(tr$edge.length < etol)
    tr$edge.length[einds] <- etol
    tr        <- multi2di(tr)
    tr        <- ladderize(tr)
    tr        <- multi2di(tr)
  }
  
  tinds <- match(tr$tip.label, metadata$sequence_name)
  all(tr$tip.label==metadata$sequence_name[tinds])
  tr$continent <- metadata$continent[tinds]
  tr$country   <- metadata$country[tinds]
  tr$country_state <- metadata$country_state[tinds]
  tr$uk_county <- metadata$uk_county[tinds]
  tr$lineages  <- metadata$lineage[tinds]
  
  place <- paste(tr$continent)
  place[which(tr$country=="UK")] <- paste(tr$country_state[which(tr$country=="UK")])
  tr$place <- place
  
  ace_place_ER <- ace(place, tr, type="discrete", model="ER")
  
  ancs     <- colnames(ace_place_ER$lik.anc)[apply(ace_place_ER$lik.anc, 1, which.max)]
  state    <- c(paste(place),ancs)
  tr$state <- state
  
  fromToTbl <- table(  paste("From",tr$state[tr$edge[,1]]), paste("To",tr$state[tr$edge[,2]]) )
  
  save(tr, place, ancs, state,
       ace_place_ER, fromToTbl,
       file=paste(outPath,"cog_global_tr_and_ace_place.Robjs.Rdata",sep=""))
  
  write.csv(fromToTbl, file=paste(outPath,"cog_global_ace_place_fromToTbl.csv",sep=""))
  
  #load(paste(dataPath,"cog_global_tr_and_ace_place.Robjs.Rdata",sep=""))
  
  place_cols <- colourize_BEAST_cols(tr$state)
  place_tbl  <- table(place)
  tinds      <- 1:length(tr$tip.label)
  ninds      <- (length(tr$tip.label)+1):(max(tr$edge))
  ecols      <- place_cols$statecols[tr$edge[,1]]
  ptxt <- paste(rownames(place_tbl)," (",place_tbl,")",sep="")
  scot_inds <- which(place=="Scotland")
  tcols     <- place_cols$statecols[tinds]
  #tcols     <- array(NA,length(place))
  tcols[scot_inds] <- "black"
  #tcols2    <- array(NA,length(place))
  tcols2    <- place_cols$statecols[tinds]
  #tcols2[scot_inds]<- place_cols$ucols[which(place_cols$ustates=="Scotland")]
  ppch      <- array(NA,length(place))
  ppch[scot_inds]<-21
  
  trImPath <- paste(outPath,dataName,"_place_tree_image.png",sep="")
  png(file=trImPath, height=12*150, width=8*150, res=150)
  plot(tr, show.tip=FALSE, edge.color=ecols)
  tiplabels(pch=ppch,col=tcols,bg=tcols2)
  upch <- array(NA,length(ptxt))
  upch[which(rownames(place_tbl)=="Scotland")] <- 21
  legend("bottomleft",ptxt,pch=upch,lty=1,col=place_cols$ucols,pt.bg=place_cols$ucols,bty="n")
  title("Global Tree with Place")
  add.scale.bar(2.5e-4,0,1e-4)
  dev.off()
  
  for (k in 1:length(u_uk_states)) {
    uk_part <- u_uk_states[k]
  
    ss <- which(rownames(place_tbl)==uk_part)
    notss <- setdiff(1:length(rownames(place_tbl)),ss)
  
    imPath <- paste(outPath,dataName,"_place_tree_imports_",uk_part,".png",sep="")
    png(file=imPath, height=4*150, width=12*150, res=150)
    barplot(fromToTbl[notss,ss],col=place_cols$ucols[notss],names=rownames(place_tbl)[notss],cex.names=0.75)
    title(paste("Number of import events to",uk_part,"=",sum(fromToTbl[notss,ss])))
    dev.off()
    
    write.table(t(fromToTbl[notss,ss]), row.names=FALSE, col.names=TRUE, sep=",",
                file=paste(outPath,dataName,"_place_tree_imports_exports_",uk_part,".csv",sep=""))
    
    write.table(t(fromToTbl[ss,notss]), row.names=FALSE, col.names=TRUE, sep=",",
                file=paste(outPath,dataName,"_place_tree_imports_exports_",uk_part,".csv",sep=""), append=TRUE)
    
    imPath <- paste(outPath,dataName,"_place_tree_exports_",uk_part,".png",sep="")
    png(file=imPath, height=4*150, width=12*150, res=150)
    barplot(-fromToTbl[ss,notss],col=place_cols$ucols[notss],names=rownames(place_tbl)[notss],cex.names=0.75)
    title(paste("Number of export events from",uk_part,"=",sum(fromToTbl[ss,notss])))
    dev.off()
  }
  
  logName <- paste(dataPath,"risefall_scotland_COG_metadata_IMPORTS_1.txt",sep="")
  #startTime <- Sys.time()
  write("# Log file for risefall_scotland_COG_metadata_IMPORTS_1",file=logName,append=FALSE)
  write("# START", file=logName,append=TRUE)
  write(paste("# process start",startTime), file=logName,append=TRUE)
  write(paste("dataDate=",dataDate),file=logName,append=TRUE)
  write(paste("rootPath=",rootPath),file=logName,append=TRUE)
  write(paste("metaPath=",metaPath),file=logName,append=TRUE)
  write(paste("treeName=",treeName),file=logName,append=TRUE)
  write(paste("dataPath=",dataPath),file=logName,append=TRUE)
  write(paste("outPath=",outPath),file=logName,append=TRUE)
  write("###########################",file=logName,append=TRUE)
  write(paste("# SUMMARY OF IMPORTS CALCULATION"), file=logName,append=TRUE)
  write(paste("# Using",treeName,"with place"),file=logName,append=TRUE)
  write(paste("numPlaces=",length(unique(place))), file=logName,append=TRUE)
  write(paste("numTips=",length(tr$tip.label)), file=logName,append=TRUE)
  write(paste("numNodes=",length(tr$Nnode)), file=logName,append=TRUE)
  write(paste("numEdges=",length(tr$edge.length)), file=logName,append=TRUE)
  write(paste("model=",ace_place_ER$call),file=logName, append=TRUE)
  write("# model details: ACE discrete traits with ER model on place", file=logName, append=TRUE)
  write("# Places",file=logname,append=TRUE)
  write.table(table(place), file=logName, append=TRUE, row.names=TRUE, col.names=TRUE, sep=",")
  write("# From-to table",file=logName,append=TRUE)
  write(paste(c("fromToTbl",colnames(fromToTbl)),collapse=","),file=logName,append=TRUE)
  write.table(fromToTbl, col.names=FALSE, row.names=TRUE, file=logName,append=TRUE, sep=",")
  write("###########################",file=logName,append=TRUE)
  write(paste("fromToTbl filename= cog_global_ace_place_fromToTbl.csv"),file=logName,append=TRUE)
  write(paste("objects filename= cog_global_tr_and_ace_place.Robjs.Rdata"),file=logName,append=TRUE)
  
  endTime <- Sys.time()
  write(paste("# END SUMMARY",file=logName, append=TRUE))
  write(paste("# process run time",endTime-startTime),file=logName, append=TRUE)
  write(paste("# process end",endTime), file=logName,append=TRUE)
  write(paste("# END"),file=logName,append=TRUE)
  
  #if (echo) {
  #  print("Finshed")
  #}
  #return( logName )
  
}

# EXTRACT SCOTLAND TIPS FROM GLOBAL TREE
if (doScotTree) {
  tr <- read.tree(paste(metaPath,treeName,sep=""))
  save(tr, file=paste(dataPath,"cog_global.tr.Rdata",sep=""))
  
  # extract the scottish sequences fron the tree
  scot_inds <- which(metadata$country_state=="Scotland")
  scot_tips <- match( metadata$sequence_name[scot_inds], tr$tip.label )
  all( metadata$sequence_name[scot_inds]==tr$tip.label[scot_tips])
  ex_tips   <- setdiff(1:length(tr$tip.label), scot_tips)
  scot_tr   <- drop.tip(tr, ex_tips)
  scot_tr   <- ladderize(scot_tr)

  # needed for treedater
  scot_tr <- multi2di(scot_tr)
  
  # get rid of negative branch lengths - they will cause grief later
  #etol  <- 1e-12
  einds <- which(scot_tr$edge.length < etol)
  scot_tr$edge.length[einds] <- etol
  scot_tr <- multi2di(scot_tr)
  scot_tr <- ladderize(scot_tr)
  
  # by region
  tinds <- match(scot_tr$tip.label, metadata$sequence_name)
  all( scot_tr$tip.label == metadata$sequence_name[tinds] )
  scot_tr$tip.label <- newTaxa[tinds]
  
  scot_tr$uk_county <- factor(paste(uk_county[tinds]), levels=u_scot)
  scot_tr$decDate   <- decDates[tinds]
  
  sts               <- scot_tr$decDate
  names(sts)        <- scot_tr$tip.label
  #seqlen            <- 30000
  #omega0            <- 1e-3
  #minblen           <- etol/omega0
  #meanRateLimits    <- c(1e-5,5e-3)
  
  # compose object for tree dater ?
  tr_res <- list(tr=scot_tr, sts=sts, seqlen=seqlen, 
                 omega0=omega0, 
                 meanRateLimits=meanRateLimits,
                 minblen=minblen,
                 info=metadata[tinds,])
  
  save(tr_res,  file=paste(dataPath,"cog_global_scot_tree.tr_res.Rdata",sep=""))
  save(scot_tr, file=paste(dataPath,"cog_global_scot_tree.scot_tr.Rdata",sep=""))
  write.tree(scot_tr, file=paste(dataPath,"cog_global_scot_tree.nwk",sep=""))
  
  ace_res <- ace(scot_tr$uk_county, scot_tr, type="discrete", model="ER")
  ancs    <- colnames(ace_res$lik.anc)[apply(ace_res$lik.anc, 1, which.max)]
  #ancs    <- factor(ancs, levels=u_scot)
  scot_tr$state <- factor( c(paste(scot_tr$uk_county), ancs), levels=u_scot)
  scot_tr$ancs  <- ancs
  save(scot_tr, ace_res, file=paste(dataPath,"cog_global_scot_tree_with_ace_ER.scot_tr.Rdata",sep=""))
  
  
  tr_cols <- colourize_BEAST_cols(scot_tr$state, ustates=scot_cols$ustates,
                                  ucols=scot_cols$ucols)
  tipinds <- 1:length(scot_tr$tip.label)
  nodeinds<- (length(scot_tr$tip.label)+1):(max(scot_tr$edge[,1]))
  ecols   <- tr_cols$statecols[scot_tr$edge[,1]]
  
  imageName <- paste(dataPath,"cog_global_scot_tree_with_region.png",sep="")
  png(file=imageName, height=12*300, width=12*300, res=300)
  op <- par(mar=c(0,0,2,0))
  plot(scot_tr, show.tip.label=FALSE, edge.color = ecols)
  tiplabels(pch=21, bg=tr_cols$statecols[tipinds])
  nodelabels(pch=1, col=tr_cols$statecols[nodeinds])
  legend("bottomleft",scot_txt,pch=21,pt.bg=tr_cols$ucols,bty="n")
  title("Scotland Sequences")
  par(op)
  dev.off()
  
}

# ADD SCOTLAND REGIONS TO SCOTLAND TREE AND EXTRACT SUB-TREE PER REGION ONLY
if (doScotRegions) {
  if (!doScotTree) {
    load(paste(dataPath,"cog_global_scot_tree_with_ace_ER.scot_tr.Rdata",sep=""))
  }
  
  dp2 <- "scotland_region_only//"
  
  nthres     <- 10
  bigRegions <- setdiff(rownames(scot_tbl)[which(scot_tbl >= nthres)],"UNKNOWN")
  nbig       <- length(bigRegions)
  for (i in 1:nbig) {
    inds      <- which(scot_tr$uk_county==bigRegions[i])
    ex_inds   <- setdiff(1:length(scot_tr$tip.label), inds)
    region_tr <- drop.tip(scot_tr, ex_inds)
    write.tree(region_tr, file=paste(dataPath,dp2,"cog_global_scotland_",gsub(" ","",bigRegions[i]),".nwk",sep=""))
    save(region_tr, file=paste(dataPath,dp2,"cog_global_scotland_",gsub(" ","",bigRegions[i]),".tr.Rdata",sep=""))
    region_taxa <- apply(as.matrix(region_tr$tip.label), 1, getEl, ind=1, sep="\\|")
    write(region_taxa, file=paste(dataPath,dp2,"cog_global_scotland_",gsub(" ","",bigRegions[i]),"_taxaList.txt",sep=""))
  }
}

# EXTRACT SUBTREES FROM SCOTLAND-REGION - THESE HAVE SCOTLAND SEQUENCES AND OTHER RELATED
if (doScotLineages) {
  if (!doScotTree) {
    load(paste(dataPath,"cog_global.tr.Rdata",sep="") )
    load(paste(dataPath,"cog_global_scot_tree_with_ace_ER.scot_tr.Rdata",sep=""))
  }
  nthres     <- 10
  bigRegions <- setdiff(rownames(scot_tbl)[which(scot_tbl >= nthres)],"UNKNOWN")
  nbig       <- length(bigRegions)
  
  bars <- apply(as.matrix(dates_res$seldates), 1, calcDecimalDate_fromTxt, sep="-", dayFirst=FALSE)
  barsTxt <- dates_res$seldates2
  
  dp2 <- "scotland_subtrees_and_related//"
  
  for (i in 1:nbig) {
    inds        <- which(uk_county==bigRegions[i])
    region_meta <- metadata[inds,]
    ulin        <- unique(paste(region_meta$lineage))
    nulin       <- length(ulin)
    for (k in 1:nulin) {
  
      sub_trName <- paste("cog_global_scotland_",bigRegions[i],"_",gsub("\\.","-",ulin[k]),"_sub_tree.dtr.Rdata",sep="")
      fnames <- dir(paste(dataPath,dp2,sep=""))
      pos    <- grep(sub_trName,fnames,fixed=TRUE)
      dname  <- paste(dataPath,dp2,sub_trName,sep="")
      
      linds <- which(region_meta$lineage==ulin[k])
      
    if (length(pos)==0) {
      ltaxa <- region_meta$sequence_name[linds]
      tinds <- match(ltaxa, tr$tip.label)
      if (length(tinds)>=2) {
        anc   <- getMRCA(tr, tinds)
      } else {
        anc   <- tr$edge[which(tr$edge[,2]==tinds),1]
      }
      sub_tr <- extract.clade(tr, anc)
      minds  <- match(sub_tr$tip.label, metadata$sequence_name)
      sub_tr$tip.label <- metadata$newTaxa[minds]
      
      if (length(sub_tr$tip.label) > 200) {
        # really need to cut out some things
        sub_meta  <- metadata[minds,]
        ok_inds   <- which(sub_meta$uk_county==bigRegions[i])
        numRegion <- length(ok_inds)
        
        mper <- 20
        other_scot_inds <- which(sub_meta$country_state=="Scotland" & sub_meta$uk_county!=bigRegions[i])
        other_eng_inds  <- which(sub_meta$country_state=="England")
        other_wls_inds  <- which(sub_meta$country_state=="Wales")
        other_nir_inds  <- which(sub_meta$country_state=="Northern Ireland")
        other_europe_inds<- which(sub_meta$continent=="Europe" & sub_meta$country!="UK")
        other_world_inds<- which(sub_meta$continent!="Europe")
        
        if (length(other_scot_inds)>mper)   other_scot_inds   <- sample(other_scot_inds,mper)
        if (length(other_eng_inds)>mper)    other_eng_inds    <- sample(other_eng_inds,mper)
        if (length(other_wls_inds)>mper)    other_wls_inds    <- sample(other_wls_inds,mper)
        if (length(other_nir_inds)>mper)    other_nir_inds    <- sample(other_nir_inds,mper)
        if (length(other_europe_inds)>mper) other_europe_inds <- sample(other_europe_inds,mper)
        if (length(other_world_inds)>mper)  other_world_inds  <- sample(other_world_inds,mper)
        
        inc_inds <- c(ok_inds,other_scot_inds,other_eng_inds,other_wls_inds,other_nir_inds,other_europe_inds,other_world_inds)
        keep_taxa<- sub_meta$newTaxa[inc_inds]
        inc_inds <- match(keep_taxa, sub_tr$tip.label)
        ex_inds  <- setdiff(1:length(sub_tr$tip.label), inc_inds)
        sub_tr   <- drop.tip(sub_tr, ex_inds)
        minds    <- match(sub_tr$tip.label, metadata$newTaxa)
      }
      
      # get rid of negative branch lengths - they will cause grief later
      sub_tr<- multi2di(sub_tr)
      einds <- which(sub_tr$edge.length < etol)
      sub_tr$edge.length[einds] <- etol
      sub_tr <- multi2di(sub_tr)
      sub_tr <- ladderize(sub_tr)
      
      # needed for treedater
      sts <- metadata$decDates[minds]
      names(sts) <- sub_tr$tip.label
      
      # do tree dating
      dtr <- dater(sub_tr, sts=sts, s=seqlen, omega0=omega0, minblen=minblen,
                   clock=clock, numStartConditions=0,
                   searchRoot=0, meanRateLimits=meanRateLimits)
      dtr <- nodeTimes(dtr, youngestTip=max(sts))
      rootnode <- length(dtr$tip.label)+1
      nodeinds <- rootnode:max(dtr$edge[,1])
      
      # add location details to object
      tipstate <- array("Other",length(sub_tr$tip.label))
      tipcols  <- array("grey80", length(sub_tr$tip.label))
      tipcols2 <- array(NA, length(sub_tr$tip.label))
      minds    <- match(sub_tr$tip.label, metadata$newTaxa)
      sub_meta <- metadata[minds,]
      for (j in 1:length(u_uk_states)) {
        tipcols[which(sub_meta$country_state==u_uk_states[j])] <- uk_cols$ucols[j]
        tipstate[which(sub_meta$country_state==u_uk_states[j])]<- u_uk_states[j]
      }
      for (j in 1:(length(u_scot)-1)) {
        tipcols2[which(sub_meta$uk_county==u_scot[j])] <- scot_cols$ucols[j]
        tipstate[which(sub_meta$uk_county==u_scot[j])] <- u_scot[j]
      }
      dtr$tipstate <- tipstate
      dtr$tipcols  <- tipcols
      dtr$tipcols2 <- tipcols2
      save(dtr, file=dname)
    } else {
      load(dname)
      tipstate <- dtr$tipstate
      tipcols  <- dtr$tipcols
      tipcols2 <- dtr$tipcols2
      rootnode <- length(dtr$tip.label)+1
      nodeinds <- rootnode:max(dtr$edge[,1])
    }
      
      #ace_res <- ace(tipstate, dtr, type="discrete", model="SYM")
      #ancs    <- colnames(ace_res$lik.anc)[apply(ace_res$lik.anc, 1, which.max)]
      #state   <- c(tipstate, ancs)
      #dtr$state <- state
      
      sub_u   <- sort(unique(tipstate))
      if (any(sub_u=="Other")) {
        sub_u <- c(setdiff(sub_u, "Other"),"Other")
      }
      sub_tbl <- table(factor(tipstate, levels=sub_u))
      sub_txt <- paste(sub_u," (",sub_tbl,")",sep="")
      sub_ucol<- array(0,length(sub_u))
      sub_ucol2<- array(0, length(sub_u))
      for (j in 1:length(sub_u)) {
        sub_ucol[j] <- tipcols[which(tipstate==sub_u[j])[1]]
        sub_ucol2[j] <- tipcols2[which(tipstate==sub_u[j])[1]]
      }
      
      #rootDate <- invertDecimalDate(dtr$timeOfMRCA,formatAsTxt=TRUE,ddmmyy=TRUE)
      
      imageName <- paste(dataPath,dp2,"cog_global_scotland_",bigRegions[i],"_",gsub("\\.","-",ulin[k]),"_sub_tree.png",sep="")
      png(file=imageName, height=8*150, width=8*150, res=150)
        op <- par(mar=c(1,0,3,0))
        plot(dtr, tip.color=tipcols)
        for (j in 1:length(bars)) {
          lines(  c(bars[j],bars[j])-dtr$timeOfMRCA,  c(-1,length(dtr$tip.label)+1), lty=2, col=hsv(0,0,0.8,0.5))
        }
        text(bars-dtr$timeOfMRCA, array(1,length(bars)), barsTxt, col=hsv(0,0,0.8,0.5), pos=1)
        text(bars-dtr$timeOfMRCA, array(length(dtr$tip.label),length(bars)), barsTxt, col=hsv(0,0,0.8,0.5), pos=3)
        tiplabels(pch=21,col=tipcols,bg=tipcols2)
        nodelabels(format(dtr$timeOfMRCA,digits=6), rootnode, pos=4,bg=NULL, frame="none")
        legend("bottomleft",sub_txt,pch=21,col=sub_ucol,pt.bg=sub_ucol2,bty="n")
        title( paste(bigRegions[i],ulin[k],"only n=",length(linds)))
        par(op)
      dev.off()
      
      out_tr_name <- paste(dataPath,dp2,"cog_global_scotland_",bigRegions[i],"_",gsub("\\.","-",ulin[k]),"_sub_tree_dated.nwk",sep="")
      write.tree(dtr, file=out_tr_name)
      
    }
  }
  
}

if (doScotLineages) {
  
  # PART 2 of Scot Lineages
  if (!doScotTree) {
    load(paste(dataPath,"cog_global.tr.Rdata",sep="") )
    load(paste(dataPath,"cog_global_scot_tree_with_ace_ER.scot_tr.Rdata",sep=""))
  }
  nthres     <- 10
  bigRegions <- setdiff(rownames(scot_tbl)[which(scot_tbl >= nthres)],"UNKNOWN")
  nbig       <- length(bigRegions)
  
  bars <- apply(as.matrix(dates_res$seldates), 1, calcDecimalDate_fromTxt, sep="-", dayFirst=FALSE)
  barsTxt <- dates_res$seldates2
  
  dp2 <- "scotland_subtrees_and_related//"
  
  for (i in 1:nbig) {
    inds        <- which(uk_county==bigRegions[i])
    region_meta <- metadata[inds,]
    ulin        <- unique(paste(region_meta$lineage))
    nulin       <- length(ulin)
    
    # recompose list of all sequences in this region and related
    region_rel_taxa <- c()
    for (k in 1:nulin) {
      sub_trName <- paste("cog_global_scotland_",bigRegions[i],"_",gsub("\\.","-",ulin[k]),"_sub_tree.dtr.Rdata",sep="")
      fnames <- dir(paste(dataPath,dp2,sep=""))
      pos    <- grep(sub_trName,fnames,fixed=TRUE)
      dname  <- paste(dataPath,dp2,sub_trName,sep="")
      linds  <- which(region_meta$lineage==ulin[k])
      
      load(dname)
      region_rel_taxa <- c(region_rel_taxa, dtr$tip.label)
    }
    region_rel_taxa <- unique(region_rel_taxa)
    region_rel_ids  <- metadata$sequence_name[match(region_rel_taxa, metadata$newTaxa)]
    
    # re-extract from global tree and do dating
    todrop <- setdiff(1:length(tr$tip.label),match(region_rel_ids,tr$tip.label))
    sub_tr <- drop.tip(tr, todrop)
    minds  <- match(sub_tr$tip.label, metadata$sequence_name)
    #all(sub_tr$tip.label==metadata$sequence_name[minds])
    sub_tr$tip.label <- metadata$newTaxa[minds]
    minds  <- match(sub_tr$tip.label, metadata$newTaxa)
  
    # get rid of negative branch lengths - they will cause grief later
    sub_tr<- multi2di(sub_tr)
    einds <- which(sub_tr$edge.length < etol)
    sub_tr$edge.length[einds] <- etol
    sub_tr <- multi2di(sub_tr)
    sub_tr <- ladderize(sub_tr)
  
    # needed for treedater
    sts <- metadata$decDates[minds]
    names(sts) <- sub_tr$tip.label
  
    # do tree dating
    dtr <- dater(sub_tr, sts=sts, s=seqlen, omega0=omega0, minblen=minblen,
               clock=clock, numStartConditions=0,
               searchRoot=0, meanRateLimits=meanRateLimits)
    dtr <- nodeTimes(dtr, youngestTip=max(sts))
    rootnode <- length(dtr$tip.label)+1
    nodeinds <- rootnode:max(dtr$edge[,1])
  
    # add location details to object
    tipstate <- array("Other",length(sub_tr$tip.label))
    tipcols  <- array("grey80", length(sub_tr$tip.label))
    tipcols2 <- array(NA, length(sub_tr$tip.label))
    minds    <- match(sub_tr$tip.label, metadata$newTaxa)
    sub_meta <- metadata[minds,]
    for (j in 1:length(u_uk_states)) {
      tipcols[which(sub_meta$country_state==u_uk_states[j])] <- uk_cols$ucols[j]
      tipstate[which(sub_meta$country_state==u_uk_states[j])]<- u_uk_states[j]
    }
    for (j in 1:(length(u_scot)-1)) {
      tipcols2[which(sub_meta$uk_county==u_scot[j])] <- scot_cols$ucols[j]
      tipstate[which(sub_meta$uk_county==u_scot[j])] <- u_scot[j]
    }
    dtr$tipstate <- tipstate
    dtr$tipcols  <- tipcols
    dtr$tipcols2 <- tipcols2
    
    tt <- metadata$sequence_name[match(dtr$tip.label,metadata$newTaxa)]
    dtr$sequence_name <- tt
    dtr$lineage <- metadata$lineage[match(dtr$tip.label,metadata$newTaxa)]
    
    ace_ER <- ace(tipstate, dtr, type="discrete", model="ER")
    ancs   <- colnames(ace_ER$lik.anc)[apply(ace_ER$lik.anc, 1, which.max)]
    state  <- c(paste(tipstate),ancs)
    dtr$state  <- state
    dtr$ace_ER <- ace_ER
    dtr$ancs   <- ancs
    
    fromToTbl <- table(paste("From",dtr$state[dtr$edge[,1]]),paste("To",dtr$state[dtr$edge[,2]]))
    write.csv(fromToTbl, file=paste("cog_global_scotland_",bigRegions[i],"_all_sub_tree_fromToTbl.csv",sep=""))
    
    dtr$fromToTbl <- fromToTbl
    
    sub_trName  <- paste("cog_global_scotland_",bigRegions[i],"_all_sub_tree.dtr.Rdata",sep="")
    dname2      <- paste(dataPath,dp2,sub_trName,sep="")
    save(dtr, file=dname2)
    
    sub_u   <- sort(unique(tipstate))
    if (any(sub_u=="Other")) {
      sub_u <- c(setdiff(sub_u, "Other"),"Other")
    }
    sub_tbl <- table(factor(tipstate, levels=sub_u))
    sub_txt <- paste(sub_u," (",sub_tbl,")",sep="")
    sub_ucol<- array(0,length(sub_u))
    sub_ucol2<- array(0, length(sub_u))
    for (j in 1:length(sub_u)) {
      sub_ucol[j] <- tipcols[which(tipstate==sub_u[j])[1]]
      sub_ucol2[j] <- tipcols2[which(tipstate==sub_u[j])[1]]
    }
    
    imageName <- paste(dataPath,dp2,"cog_global_scotland_",bigRegions[i],"_all_sub_tree.png",sep="")
    png(file=imageName, height=8*150, width=8*150, res=150)
    op <- par(mar=c(1,0,3,0))
    #plot(dtr, tip.color=tipcols)
    plot(dtr, show.tip.label=FALSE)
    for (j in 1:length(bars)) {
      lines(  c(bars[j],bars[j])-dtr$timeOfMRCA,  c(-1,length(dtr$tip.label)+1), lty=2, col=hsv(0,0,0.8,0.5))
    }
    text(bars-dtr$timeOfMRCA, array(1,length(bars)), barsTxt, col=hsv(0,0,0.8,0.5), pos=1)
    text(bars-dtr$timeOfMRCA, array(length(dtr$tip.label),length(bars)), barsTxt, col=hsv(0,0,0.8,0.5), pos=3)
    tiplabels(pch=21,col=tipcols,bg=tipcols2)
    tiplabels(tt,col=tipcols,cex=0.5,frame="none",pos=4)
    nodelabels(format(dtr$timeOfMRCA,digits=6), rootnode, pos=4,bg=NULL, frame="none")
    legend("bottomleft",sub_txt,pch=21,col=sub_ucol,pt.bg=sub_ucol2,bty="n")
    title( paste(bigRegions[i]) )
    par(op)
    dev.off()
    
    out_tr_name <- paste(dataPath,dp2,"cog_global_scotland_",bigRegions[i],"_all_sub_tree_dated.nwk",sep="")
    write.tree(dtr, file=out_tr_name)
      
  }
  
  # reload objects for further inspection
  for (i in 1:length(bigRegions)) {
    
    sub_trName  <- paste("cog_global_scotland_",bigRegions[i],"_all_sub_tree.dtr.Rdata",sep="")
    dname2      <- paste(dataPath,dp2,sub_trName,sep="")
    load(dname2)
    place_tbl <- table(dtr$tipstate)
    uplace    <- rownames(place_tbl)

    uplace_cols <- array("grey70",length(uplace)) 
    for (j in 1:length(u_uk_states)) {
      uplace_cols[which(uplace==u_uk_states[j])] <- uk_cols$ucols[j]
    }
    for (j in 1:length(u_scot)) {
      uplace_cols[which(uplace==u_scot[j])]      <- scot_cols$ucols[j]
    }
    
    out_inds  <- which(dtr$state[dtr$edge[,1]]==bigRegions[i])
    exports   <- table(factor(dtr$state[dtr$edge[out_inds,2]],levels=uplace))
    in_inds   <- which(dtr$state[dtr$edge[,2]]==bigRegions[i])
    imports   <- table(factor(dtr$state[dtr$edge[in_inds,1]],levels=uplace))
    ietbl   <- cbind(place_tbl,imports,exports)
    colnames(ietbl) <- c("Number of Sequences",paste("Import to",bigRegions[i]),paste("Export from",bigRegions[i]))
    
    tot_self   <- length(which(dtr$state[dtr$edge[,1]]==bigRegions[i] & dtr$state[dtr$edge[,2]]==bigRegions[i]))
    tot_import <- length(which(dtr$state[dtr$edge[,1]]!=bigRegions[i] & dtr$state[dtr$edge[,2]]==bigRegions[i]))
    tot_export <- length(which(dtr$state[dtr$edge[,1]]==bigRegions[i] & dtr$state[dtr$edge[,2]]!=bigRegions[i]))
    
    self_sum <- matrix(c(tot_self,tot_import,tot_export), 3, 1)
    rownames(self_sum) <- c("Within","Imports","Exports")
    
    
    ppName      <- paste("cog_global_scotland_",bigRegions[i],"_all_sub_tree_imports_exports_2x2pies.png",sep="")
    piePath     <- paste(dataPath,dp2,ppName,sep="")
    png(file=piePath, height=10*150, width=10*150, res=150)
    
    op <- par(mfrow=c(2,2))
    sscols <- c("white","skyblue","pink")
    pie( self_sum, labels=rownames(self_sum), col=sscols )
    title(paste("Summary of",bigRegions[i]))
    legend("topleft",paste(rownames(self_sum)," (",self_sum,")",sep=""),
           pch=22,pt.bg=sscols,bty="n")
    for (k in 1:3) {
      if (k>1) {
        #notss <- which(uplace != bigRegions[i])
        sort_ie <- sort(ietbl[,k],index.return=TRUE,decreasing=TRUE)
        sort_ie <- sort_ie$ix[which(sort_ie$x>0)]
        notss   <- setdiff(sort_ie,which(uplace==bigRegions[i]))
      } else {
        sort_ie <- sort(ietbl[,k],index.return=TRUE,decreasing=TRUE)
        sort_ie <- sort_ie$ix[which(sort_ie$x>0)]
        notss   <- c(which(uplace==bigRegions[i]),setdiff(sort_ie,which(uplace==bigRegions[i])))
      }
      if (length(notss)==0) {
        pie(ietbl[which(uplace==bigRegions[i]),k],col="grey90",labels="(none)")
        title(colnames(ietbl)[k])
        legend("topleft","(none)",pch=22,pt.bg="grey90",bty="n")
      } else {
        pie(ietbl[notss,k], col=uplace_cols[notss])
        title(colnames(ietbl)[k])
        if (k>1) {
          legend("topleft",paste(rownames(ietbl)[notss]," (",ietbl[notss,k],")",sep=""),
             pch=22,pt.bg=uplace_cols[notss],bty="n")
        }
      }
    }
    par(op)
    
    dev.off()
    
  }
  
}

