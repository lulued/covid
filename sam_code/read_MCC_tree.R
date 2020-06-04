# functions to read traits mcc trees
# S J Lycett
# 15 Feb 2016
# 21 July 2016
# 28 Dec 2018 - modified for use in R shiny
# 29 Dec 2018 - read_multi_mcc_tr correction for extra ] in discrete traits

# based on process_genbank_zika.R and global_consortium_h5nx_usa_h5n2_h5n1.R


#RcodePath <- "//ris-fas1a.roslin.ed.ac.uk//slycett2//Rcode//"
#RcodePath <- "//cmvm.datastore.ed.ac.uk//cmvm//eb//users//slycett2//Rcode//"

#source( paste(RcodePath,"birdSpecies_and_continents.R",sep="") )
#source( paste(RcodePath,"getEl.R",sep="") )
#source( paste(RcodePath,"calcDecimalDate.R",sep="") )
#source( paste(RcodePath,"figTreeSelectionAsR.R",sep="") )
#source( paste(RcodePath,"read_BEAST_MultiStateTrees.R",sep=""))
#source( paste(RcodePath,"read_BEAST_tree_latlon.R",sep=""))
#source( paste(RcodePath,"get_BEAST_cols.R",sep=""))
library(ape)

####################################################################################################
# from population_fitness_features.R

distFromRoot	<- function( tr ) {

	rootNode	<- length(tr$tip.label)+1
	nodeDists	<- array(0, max(tr$edge))
	
	toProcess	<- c(rootNode)

	while ( length(toProcess) > 0 ) {
		einds 			<- which(tr$edge[,1]==toProcess[1])

		if (length(einds) > 0) {
			children 			<- tr$edge[einds,2]
			nodeDists[children] 	<- nodeDists[toProcess[1]] + tr$edge.length[einds]

			toProcess			<- c(toProcess, children)
		}
		toProcess			<- setdiff(toProcess, toProcess[1])
	}

	tr$nodeDists <- nodeDists

	return ( tr )
}

nodeTimes	<- function(tr, youngestTip=2011.027) {

	if ( !any(attributes(tr)$names == "nodeDists") ) {
		tr <- distFromRoot(tr)
	}

	tr$rootHeight<- max(tr$nodeDists)
	tr$nodeTimes <- youngestTip - tr$rootHeight + tr$nodeDists

	return ( tr )
}

####################################################################################################

# used to get the numbers -> names translation for the taxa
getTranslation	<- function( lines ) {
	ts 		<- grep("Translate", lines)
	te 		<- grep(";", lines)
	te 		<- te[which( te > ts )][1]
	taxaLines 	<- lines[ (ts+1) : (te-1) ]

	taxaLines   <- gsub("'", "", taxaLines)

	#taxa		<- apply(as.matrix(taxaLines), 1, getEl, sep="'", ind=2)
	#tipNumbers	<- apply(as.matrix(taxaLines), 1, getEl, sep="'", ind=1)

	taxa		<- apply(as.matrix(taxaLines), 1, getEl, sep=" ", ex=1, reconstruct=TRUE)
	tipNumbers  <- apply(as.matrix(taxaLines), 1, getEl, sep=" ", ind=1)

	taxa		<- gsub(",", "", taxa)
	tipNumbers  <- gsub("\t", "", tipNumbers)
	tipNumbers  <- gsub(" ", "", tipNumbers)

	translateTbl <- cbind(tipNumbers, taxa)

	return( translateTbl )
}

getTranslation_BEAST2 <- function( lines ) {

	ts 		<- grep("Translate", lines)
	te 		<- grep(";", lines)
	te 		<- te[which( te > ts )][1]
	taxaLines 	<- lines[ (ts+1) : (te-1) ]

	taxaLines   <- gsub("'", "", taxaLines)
	pos		<- gregexpr(" [0-9]+ ",taxaLines)
	taxaLines	<- substring(taxaLines, pos)
	taxaLines	<- substring(taxaLines, 2)

	taxa		<- apply(as.matrix(taxaLines), 1, getEl, sep=" ", ex=1, reconstruct=TRUE)
	tipNumbers  <- apply(as.matrix(taxaLines), 1, getEl, sep=" ", ind=1)

	taxa		<- gsub(",", "", taxa)
	tipNumbers  <- gsub("\t", "", tipNumbers)
	tipNumbers  <- gsub(" ", "", tipNumbers)

	translateTbl <- cbind(tipNumbers, taxa)

	return( translateTbl )

}

# mod 24 april 2019 - include removal of last ] if present
addPosterior <- function( tr ) {
	posterior	 <- apply(as.matrix(tr$details), 1, getEl, ind=2, sep="posterior=")
	posterior	 <- apply(as.matrix(posterior), 1, getEl, ind=1, sep=",")
	posterior	 <- gsub("\\]","",posterior)
	posterior	 <- as.numeric(posterior)
	tr$posterior <- posterior
	return( tr )
}

addHeights 	<- function(tr) {
	h95		<- apply(as.matrix(tr$details), 1, getEl, ind=2, sep="height_95%_HPD=\\{")
	h95		<- apply(as.matrix(h95), 1, getEl, ind=1, sep="\\}")
	h95_1		<- as.numeric(apply(as.matrix(h95), 1, getEl, ind=1, sep=","))
	h95_2		<- as.numeric(apply(as.matrix(h95), 1, getEl, ind=2, sep=","))
	hmedian	<- apply(as.matrix(tr$details), 1, getEl, ind=2, sep="height_median=")
	hmedian	<- as.numeric(apply(as.matrix(hmedian), 1, getEl, ind=1, sep=","))
	hmean		<- apply(as.matrix(tr$details), 1, getEl, ind=2, sep="height=")
	hmean		<- as.numeric(apply(as.matrix(hmean), 1, getEl, ind=1, sep="\\]"))
	tr$height.median <- hmedian
	tr$height.mean   <- hmean
	tr$height.95_1  <- h95_1
	tr$height.95_2  <- h95_2

	tr$youngestTip <- max(tr$nodeTimes)
	tr$height.lower95  <- tr$youngestTip - h95_2
	tr$height.upper95  <- tr$youngestTip - h95_1
		
	return( tr ) 
}

#19 feb 2020
addRates <- function(tr) {
  rate_range      <- apply(as.matrix(tr$details), 1, getEl, ind=2, sep="rate_range=\\{")
  rate_range		  <- apply(as.matrix(rate_range), 1, getEl, ind=1, sep="\\}")
  rate_range_1		<- as.numeric(apply(as.matrix(rate_range), 1, getEl, ind=1, sep=","))
  rate_range_2		<- as.numeric(apply(as.matrix(rate_range), 1, getEl, ind=2, sep=","))
  rate <- apply(as.matrix(tr$details), 1, getEl, ind=2, sep="rate=")
  rate <- as.numeric(apply(as.matrix(rate), 1, getEl, ind=1, sep=","))
  rate_median <- apply(as.matrix(tr$details), 1, getEl, ind=2, sep="rate_median=")
  rate_median <- as.numeric(apply(as.matrix(rate_median), 1, getEl, ind=1, sep=","))
  tr$rate <- rate
  tr$rate_median <- rate_median
  tr$rate_range_1 <- rate_range_1
  tr$rate_range_2 <- rate_range_2
  
  return(tr)
}

# 24 april 2019
addMascot <- function( tr ) {
	ca_pos <- gregexpr("CA[a-zA-Z0-9_\\%\\-]+=[\\{]?[a-zA-Z0-9_\\%\\-\\,\\.]+[\\}]?",tr$details)
	firstpos<- array(-1, length(ca_pos))
	lastpos <- array(-1, length(ca_pos))
	for (j in 1:length(ca_pos)) {
		lp 		<- ca_pos[[j]]
		lastpos[j] 	<- (lp + attributes(lp)$match.length)[length(lp)]
		firstpos[j] <- lp[1]
	}
	ca_details <- substring(tr$details, firstpos, lastpos-1)
	details    <- substring(tr$details, lastpos+1)
	kk <- which(firstpos==-1)
	if (length(kk) > 0) {
		details[kk] <- gsub("\\[\\&","",details[kk])
	}

	hpos <- gregexpr("height",details)
	firstpos<- array(-1, length(hpos))
	for (j in 1:length(hpos)) {
		firstpos[j] <- hpos[[j]][1]
	}
	h_details 	   	<- substring(details, firstpos)
	struct_details 	<- substring(details, 1, firstpos-1)

	tr$details 		<- h_details
	tr$ca_details 	<- ca_details
	tr$struct_details <- struct_details

	tr <- addHeights(tr)
	tr <- addDiscreteTraits(tr)
	tr <- correctUprops(tr=tr)
	
	return( tr )

}

# 4 march 2019
addLatLonHPD <- function(tr, regtype=3) {
  # get the lat-lon 80% HPD (not done in original)
  res1			<- gregexpr("latlon1_80\\%HPD_1=\\{[0-9eE-\\.,]+}",tr$details)
  res2			<- gregexpr("latlon2_80\\%HPD_1=\\{[0-9eE-\\.,]+}",tr$details)
  
  pos1			<- unlist(res1)
  mat1			<- as.integer(matrix(unlist(lapply(res1, attributes)), regtype, length(res1))[1,])
  inds1			<- which(pos1 > 0)
  ninds1		<- which(pos1 < 0)
  latlon1_80		<- array(0, length(res1))
  latlon1_80[inds1] <- substring(tr$details[inds1], pos1[inds1]+18, pos1[inds1]+mat1[inds1]-2)
  latlon1_80[ninds1]<- tr$latlon[ninds1,1]
  
  pos2			<- unlist(res2)
  mat2			<- as.integer(matrix(unlist(lapply(res2, attributes)), regtype, length(res2))[1,])
  inds2			<- which(pos2 > 0)
  ninds2		<- which(pos2 < 0)
  latlon2_80		<- array(0, length(res2))
  latlon2_80[inds2] <- substring(tr$details[inds2], pos2[inds2]+18, pos2[inds2]+mat2[inds2]-2)
  latlon2_80[ninds2]<- tr$latlon[ninds2,2]
  
  lat80 <- lapply(strsplit(latlon1_80,","), as.numeric)
  lon80 <- lapply(strsplit(latlon2_80,","), as.numeric)
  tr$lat80 <- lat80
  tr$lon80 <- lon80
  
  return( tr )
}

#11 mar 2019
addDiscreteTraits <- function(tr) {
# find the discrete traits
	discreteRegex <- "\\,[A-Za-z0-9_\\-]+\\.prob"
	pos		  <- gregexpr(discreteRegex,tr$details[1])[[1]]
	is		  <- pos+1
	ie		  <- is + attributes(pos)$match.length - 7
	traitNames	  <- substring(tr$details[1], is, ie)
	numTraits	  <- length(traitNames)
	
	# get discrete values
	res		  <- apply(as.matrix(paste(traitNames,"=",sep="")), 1, gregexpr, tr$details)
	props		  <- matrix(0, length(tr$details), numTraits)
	colnames(props)<- traitNames
	uprops	  <- vector("list",numTraits)
	for (j in 1:numTraits) {
		is 		<- unlist(res[[j]]) + nchar(traitNames[j]) + 1
		nodeVals 	<- substring(tr$details, is)
		nodeVals	<- apply(as.matrix(nodeVals), 1, getEl, ind=1, sep=",")
		nodeVals	<- gsub("\"", "", nodeVals)
		nodeVals  <- gsub("\\]","",nodeVals)
		props[,j]	<- nodeVals
		uprops[[j]] <- sort(unique(nodeVals))
	}
	tr$props 	<- props
	tr$uprops	<- uprops
	tr$propNames<- traitNames

	# get discrete probabilities
	res			<- apply(as.matrix(paste(traitNames,"\\.prob=",sep="")), 1, gregexpr, tr$details)
	propProb		<- matrix(0, length(tr$details), numTraits)
	colnames(propProb)<- traitNames
	for (j in 1:numTraits) {
		is 		<- unlist(res[[j]]) + nchar(traitNames[j]) + 6
		nodeVals 	<- substring(tr$details, is)
		nodeVals	<- apply(as.matrix(nodeVals), 1, getEl, ind=1, sep=",")
		nodeVals	<- gsub("\"", "", nodeVals)
		propProb[,j]<- as.numeric(nodeVals)
	}
	tr$propProb <- propProb
	return (tr)
}

# 16 arch 2019 - added utraits
getTraitSet <- function( detail, traitName, orderTraits=TRUE, utraits=c(-1) ) {
	r1 <- paste(traitName,"\\.set",sep="")
	r2 <- paste(traitName,"\\.set\\.prob",sep="")
	pos2<- gregexpr(r2,detail)[[1]][1]
	pos1<- setdiff(gregexpr(r1,detail)[[1]],pos2)
	trait_set <- substring(detail,pos1)
	trait_set <- strsplit(trait_set,"\\}")[[1]][1]
	trait_set <- strsplit(trait_set,"\\{")[[1]][2]
	trait_set <- gsub("\"","",trait_set)
	trait_set <- strsplit(trait_set,",")[[1]]
	trait_prob<- substring(detail,pos2)
	trait_prob<- strsplit(trait_prob,"\\}")[[1]][1]
	trait_prob<- strsplit(trait_prob,"\\{")[[1]][2]
	trait_prob<- strsplit(trait_prob,",")[[1]]
	trait_prob<- as.numeric(trait_prob)
	trait_prob<- as.matrix(trait_prob)
	rownames(trait_prob) <- trait_set
	if (orderTraits) {
		torder     		   <- sort(trait_set, index.return=TRUE)$ix
		trait_prob 		   <- trait_prob[torder,]
	} else {
		if (utraits[1] != -1) {
			torder <- match(trait_set, utraits)
			tp	 <- matrix(0, length(utraits), 1)
			tp[torder] <- trait_prob
			rownames(tp) <- utraits
			trait_prob <- tp
		}
	}
	trait_prob <- t(trait_prob)
	
	return( trait_prob )
}

correctUprops   <- function(tr=tr) {
	for (i in 1:length(tr$uprops)) {
		uprops <- as.matrix(tr$uprops[[i]])
		jj 	 <- grep("\\+",uprops)
		if (length(jj) > 0) {
			uprops <- setdiff(uprops, uprops[jj])
			tr$uprops[[i]] <- uprops
			print(paste("Corrected",tr$propNames[i]))
		}
	}
	return(tr)
}

getFullTraitSet <- function(tr=tr, propIndex=1) {
	utraits 	<- tr$uprops[[propIndex]]
	traitName	<- tr$propNames[[propIndex]]
	ntips 	<- length(tr$tip.label)
	ancs		<- (ntips+1):(ntips+tr$Nnode)
	res		<- apply(as.matrix(tr$details[ancs]), 1, getTraitSet, traitName, orderTraits=FALSE, utraits=utraits)
	rownames(res) <- utraits

	tres		<- matrix(0, length(utraits), ntips)
	for (j in 1:length(utraits)) {
		inds <- which(tr$props[1:ntips,propIndex]==utraits[j])
		if (length(inds)>0) {
			tres[j,inds] <- 1
		}
	}
	rownames(tres) <- utraits
	res <- t( cbind(tres,res) )
	return( res )
}

######################################################################

# for plain mcc tree (no traits)
read_mcc_tr <- function( trName, sep="\\|", BEAST2=FALSE ) {

	trTxt <- readLines(trName)
	trLine<- trTxt[ length(trTxt) - 1]
	if (BEAST2) {
		trLine<- strsplit(trLine, "tree TREE1 = ")[[1]][2]
	} else {
		trLine<- strsplit(trLine, "\\[\\&R\\] ")[[1]][2]
	}
	pos1	<- gregexpr("\\[\\&",trLine)
	pos2  <- gregexpr("\\]",trLine)
	nnodes<- length(pos1[[1]])
	nodeDetails <- array(0, nnodes)
	for (j in 1:nnodes) {
		nodeDetails[j] <- substring(trLine, pos1[[1]][j], pos2[[1]][j])
	}
	for (j in 1:nnodes) {
		trLine <- gsub(nodeDetails[j], paste("_nodeDetails_",j,sep=""), trLine, fixed=TRUE)
	}
	tr <- read.tree(text=trLine)
	tr$tip.details <- as.integer(apply(as.matrix(tr$tip.label), 1, getEl, ind=1, fromEnd=TRUE, sep="_"))
	tr$node.details<- as.integer(apply(as.matrix(tr$node.label), 1, getEl, ind=1, fromEnd=TRUE, sep="_"))
	tr$tip.details <- nodeDetails[tr$tip.details]
	tr$node.details<- nodeDetails[tr$node.details]
	tr$tip.label   <- apply(as.matrix(tr$tip.label), 1, getEl, ind=1, sep="_")
	
	if (BEAST2) {
		translateTbl <- getTranslation_BEAST2(trTxt)
	} else {
		translateTbl <- getTranslation(trTxt)
	}
	tr$tip.label   <- translateTbl[match(tr$tip.label, translateTbl[,1]),2]
	tr$node.label  <- NULL

	#tr2 <- read.nexus(trName)
	#dist.topo(tr2,tr)

	tr$details	   <- c(tr$tip.details, tr$node.details)
	
	decDates	   <- as.numeric(apply(as.matrix(tr$tip.label), 1, getEl, ind=1, sep=sep, fromEnd=TRUE))
	youngestTip	   <- max(decDates)
	tr		   <- nodeTimes(tr, youngestTip=youngestTip)

	#all(tr$nodeTimes[1:length(tr$tip.label)]==decDates)
	all( abs(tr$nodeTimes[1:length(tr$tip.label)]-decDates) < 1e-5 )

	tr		   <- addPosterior(tr)
	if (!BEAST2) {
		tr	   <- addHeights(tr)
	}

	return( tr )
}

read_latlon_mcc_tr <- function( trName ) {

	trTxt <- readLines(trName)
	trLine<- trTxt[ length(trTxt) - 1]
	trLine<- strsplit(trLine, "\\[\\&R\\] ")[[1]][2]
	pos1	<- gregexpr("\\[\\&",trLine)
	pos2  <- gregexpr("\\]",trLine)
	nnodes<- length(pos1[[1]])
	nodeDetails <- array(0, nnodes)
	for (j in 1:nnodes) {
		nodeDetails[j] <- substring(trLine, pos1[[1]][j], pos2[[1]][j])
	}
	for (j in 1:nnodes) {
		trLine <- gsub(nodeDetails[j], paste("_nodeDetails_",j,sep=""), trLine, fixed=TRUE)
	}
	tr <- read.tree(text=trLine)
	tr$tip.details <- as.integer(apply(as.matrix(tr$tip.label), 1, getEl, ind=1, fromEnd=TRUE, sep="_"))
	tr$node.details<- as.integer(apply(as.matrix(tr$node.label), 1, getEl, ind=1, fromEnd=TRUE, sep="_"))
	tr$tip.details <- nodeDetails[tr$tip.details]
	tr$node.details<- nodeDetails[tr$node.details]
	tr$tip.label   <- apply(as.matrix(tr$tip.label), 1, getEl, ind=1, sep="_")
	
	translateTbl <- getTranslation(trTxt)
	tr$tip.label   <- translateTbl[match(tr$tip.label, translateTbl[,1]),2]
	tr$node.label  <- NULL

	#tr2 <- read.nexus(trName)
	#dist.topo(tr2,tr)

	tr$details	   <- c(tr$tip.details, tr$node.details)
	
	decDates	   <- as.numeric(apply(as.matrix(tr$tip.label), 1, getEl, ind=1, sep="\\|", fromEnd=TRUE))
	youngestTip	   <- max(decDates)
	tr		   <- nodeTimes(tr, youngestTip=youngestTip)

	#all(tr$nodeTimes[1:length(tr$tip.label)]==decDates)
	all( abs(tr$nodeTimes[1:length(tr$tip.label)]-decDates) < 1e-5 )
	
	latlon1	   <- apply(as.matrix(tr$details), 1, getEl, ind=2, sep="latlon1=")
	latlon1	   <- as.numeric(apply(as.matrix(latlon1), 1, getEl, ind=1, sep=","))

	latlon2	   <- apply(as.matrix(tr$details), 1, getEl, ind=2, sep="latlon2=")
	latlon2	   <- as.numeric(apply(as.matrix(latlon2), 1, getEl, ind=1, sep=","))

	tr$latlon	   <- cbind(latlon1,latlon2)
	colnames(tr$latlon) <- c("lat","lon")

	return( tr ) 

}

read_discrete_mcc_tr	<- function( trName, trait="Subtype" ) {

	trTxt <- readLines(trName)
	trLine<- trTxt[ length(trTxt) - 1]
	trLine<- strsplit(trLine, "\\[\\&R\\] ")[[1]][2]
	pos1	<- gregexpr("\\[\\&",trLine)
	pos2  <- gregexpr("\\]",trLine)
	nnodes<- length(pos1[[1]])
	nodeDetails <- array(0, nnodes)
	for (j in 1:nnodes) {
		nodeDetails[j] <- substring(trLine, pos1[[1]][j], pos2[[1]][j])
	}
	for (j in 1:nnodes) {
		trLine <- gsub(nodeDetails[j], paste("_nodeDetails_",j,sep=""), trLine, fixed=TRUE)
	}
	tr <- read.tree(text=trLine)
	tr$tip.details <- as.integer(apply(as.matrix(tr$tip.label), 1, getEl, ind=1, fromEnd=TRUE, sep="_"))
	tr$node.details<- as.integer(apply(as.matrix(tr$node.label), 1, getEl, ind=1, fromEnd=TRUE, sep="_"))
	tr$tip.details <- nodeDetails[tr$tip.details]
	tr$node.details<- nodeDetails[tr$node.details]
	tr$tip.label   <- apply(as.matrix(tr$tip.label), 1, getEl, ind=1, sep="_")
	
	translateTbl <- getTranslation(trTxt)
	tr$tip.label   <- translateTbl[match(tr$tip.label, translateTbl[,1]),2]
	tr$node.label  <- NULL

	tr2 <- read.nexus(trName)
	dist.topo(tr2,tr)

	tr$details	   <- c(tr$tip.details, tr$node.details)
	
	decDates	   <- as.numeric(apply(as.matrix(tr$tip.label), 1, getEl, ind=1, sep="\\|", fromEnd=TRUE))
	youngestTip	   <- max(decDates)
	tr		   <- nodeTimes(tr, youngestTip=youngestTip)

	#all(tr$nodeTimes[1:length(tr$tip.label)]==decDates)
	all( abs(tr$nodeTimes[1:length(tr$tip.label)]-decDates) < 1e-5 )
	
	traitVals	<- apply(as.matrix(tr$details), 1, getEl, ind=2, sep=paste(trait,"=",sep=""))
	traitVals	<- gsub("'","",traitVals)
	traitVals	<- apply(as.matrix(traitVals), 1, getEl, ind=1, sep=",")
	traitVals	<- gsub("\"","",traitVals)	   

	traitProb	<- apply(as.matrix(tr$details), 1, getEl, ind=2, sep=paste(trait,".prob=",sep=""))
	traitProb	<- as.numeric(apply(as.matrix(traitProb), 1, getEl, ind=1, sep=","))
	
	tr$state		<- traitVals
	tr$state.prob 	<- traitProb

	traitSet		<- apply(as.matrix(tr$details), 1, getEl, ind=2, sep=paste(trait,".set=\\{",sep=""))
	traitSet		<- gsub("\"","",traitSet)
	traitSet		<- apply(as.matrix(traitSet), 1, getEl, ind=1, sep="\\}")
	traitSet		<- apply(as.matrix(traitSet), 1, getEl, sep=",")
	traitSetProb	<- apply(as.matrix(tr$details), 1, getEl, ind=2, sep=paste(trait,".set.prob=\\{",sep=""))
	traitSetProb	<- apply(as.matrix(traitSetProb), 1, getEl, ind=1, sep="\\}")
	traitSetProb	<- apply(as.matrix(traitSetProb), 1, getEl, sep=",")

	ustates		<- sort(unique(tr$state))
	temp			<- matrix(0, length(traitSet), length(ustates))
	minds			<- lapply(traitSet, match, ustates)
	for (k in 1:length(traitSet)) {
		temp[k,unlist(minds[k])] <- as.numeric(unlist(traitSetProb[k]))
	}
	colnames(temp) 	<- ustates
	tr$state.set.prob <- temp

	return( tr ) 

}


# 21 July 2016
# 29 Dec 2018 - regtype=3, but formerly 2 (suspect version issue on gregexpr)
read_multi_mcc_tr <- function( trName = trName, regtype=3 ) {

	# first read the lat-lon tree
	tr	<- read_latlon_mcc_tr( trName = trName )

	# find the discrete traits
	discreteRegex <- "\\,[A-Za-z0-9_\\-]+\\.prob"
	pos		  <- gregexpr(discreteRegex,tr$details[1])[[1]]
	is		  <- pos+1
	ie		  <- is + attributes(pos)$match.length - 7
	traitNames	  <- substring(tr$details[1], is, ie)
	numTraits	  <- length(traitNames)
	
	# get discrete values
	res		  <- apply(as.matrix(paste(traitNames,"=",sep="")), 1, gregexpr, tr$details)
	props		  <- matrix(0, length(tr$details), numTraits)
	colnames(props)<- traitNames
	uprops	  <- vector("list",numTraits)
	for (j in 1:numTraits) {
		is 		<- unlist(res[[j]]) + nchar(traitNames[j]) + 1
		nodeVals 	<- substring(tr$details, is)
		nodeVals	<- apply(as.matrix(nodeVals), 1, getEl, ind=1, sep=",")
		nodeVals	<- gsub("\"", "", nodeVals)
		nodeVals  <- gsub("\\]","",nodeVals)
		props[,j]	<- nodeVals
		uprops[[j]] <- sort(unique(nodeVals))
	}
	tr$props 	<- props
	tr$uprops	<- uprops
	tr$propNames<- traitNames

	# get discrete probabilities
	res			<- apply(as.matrix(paste(traitNames,"\\.prob=",sep="")), 1, gregexpr, tr$details)
	propProb		<- matrix(0, length(tr$details), numTraits)
	colnames(propProb)<- traitNames
	for (j in 1:numTraits) {
		is 		<- unlist(res[[j]]) + nchar(traitNames[j]) + 6
		nodeVals 	<- substring(tr$details, is)
		nodeVals	<- apply(as.matrix(nodeVals), 1, getEl, ind=1, sep=",")
		nodeVals	<- gsub("\"", "", nodeVals)
		propProb[,j]<- as.numeric(nodeVals)
	}
	tr$propProb <- propProb

	# get the lat-lon 80% HPD (not done in original)
	res1			<- gregexpr("latlon1_80\\%HPD_1=\\{[0-9eE-\\.,]+}",tr$details)
	res2			<- gregexpr("latlon2_80\\%HPD_1=\\{[0-9eE-\\.,]+}",tr$details)

	pos1			<- unlist(res1)
	mat1			<- as.integer(matrix(unlist(lapply(res1, attributes)), regtype, length(res1))[1,])
	inds1			<- which(pos1 > 0)
	ninds1		<- which(pos1 < 0)
	latlon1_80		<- array(0, length(res1))
	latlon1_80[inds1] <- substring(tr$details[inds1], pos1[inds1]+18, pos1[inds1]+mat1[inds1]-2)
	latlon1_80[ninds1]<- tr$latlon[ninds1,1]

	pos2			<- unlist(res2)
	mat2			<- as.integer(matrix(unlist(lapply(res2, attributes)), regtype, length(res2))[1,])
	inds2			<- which(pos2 > 0)
	ninds2		<- which(pos2 < 0)
	latlon2_80		<- array(0, length(res2))
	latlon2_80[inds2] <- substring(tr$details[inds2], pos2[inds2]+18, pos2[inds2]+mat2[inds2]-2)
	latlon2_80[ninds2]<- tr$latlon[ninds2,2]

	lat80 <- lapply(strsplit(latlon1_80,","), as.numeric)
	lon80 <- lapply(strsplit(latlon2_80,","), as.numeric)
	tr$lat80 <- lat80
	tr$lon80 <- lon80

	return( tr )
}

#############################################################################################

plot_discrete_tree <- function(tr, propIndex=1, 
                               legpos="bottomleft",
                               show.tip.label=TRUE, 
                               tpch=23, npch=21, 
                               tcex=1.2, ncex=1,
                               edgeCols=TRUE) {
  states <- tr$props[,propIndex]
  bcols  <- colourize_BEAST_cols(states)
  ntips  <- length(tr$tip.label)
  nnodes <- length(states)
  ncols  <- bcols$statecols[(ntips+1):nnodes]
  tcols  <- bcols$statecols[1:ntips]
  if (edgeCols) {
    ecols <- bcols$statecols[tr$edge[,2]]
    plot(tr, show.tip.label=show.tip.label, edge.color=ecols)
  } else {
    plot(tr, show.tip.label=show.tip.label)
  }
  tiplabels(pch=tpch, col=tcols, bg=tcols, cex=tcex)
  nodelabels(pch=npch, col=ncols, bg=ncols, cex=ncex)
  if (legpos != "x") {
    legend(legpos,bcols$ustates,pch=tpch,col=bcols$ucols,pt.bg=bcols$ucols,bty="n")
  }
  title(tr$propNames[propIndex])
}

