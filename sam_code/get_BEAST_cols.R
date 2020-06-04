# functions to get BEAST discrete trait colours
# S J Lycett
# 22 Dec 2015
# see H5NX work
# 11 March 2017 - added transparency
# 5 Sept 2018 - add method to get mapping
# 20 Aug 2019 - added function for continous colours

# these are very slightly darker than the defaults below
get_Host_cols <- function() {
	# beast colours from figTree
	hCols		<- c(rgb(204/255,82/255,82/255),  rgb(143/255,204/255,82/255), rgb(82/255,204/255,204/255), rgb(143/255,82/255,204/255) )
	return( hCols )
}

get_BEAST_cols <- function( nstates, sat=0.7, bright=0.9, transparency=1 ) {
	if (transparency < 1) {
		ucols		<- hsv( (0:(nstates-1))/nstates, sat, bright, transparency)
	} else {
		ucols		<- hsv( (0:(nstates-1))/nstates, sat, bright)
	}
	return( ucols )
}

colourize_BEAST_cols <- function( states, ustates = sort(unique(states)), 
				nstates=length(ustates), 
				sat=0.8, bright=0.8, transparency=0.8, default=0,
				ucols = get_BEAST_cols(nstates, sat=sat, bright=bright, transparency=transparency)  ) {
	
	#ucols <- get_BEAST_cols(nstates, sat=sat, bright=bright, transparency=transparency)
	if (default==0) {
		defaultCol = hsv(0,0,bright,transparency)
	} else {
		defaultCol = ucols[default]
	}

	minds <- match(states,ustates)
	scols <- ucols[minds]
	jj	<- which(!is.finite(minds))
	if (length(jj) > 0) {
		scols[jj] <- defaultCol
	}

	return( list(ustates=ustates, ucols=ucols, statecols=scols) )

}

get_continous_cols <- function( vals, normalise=TRUE, minVal=0, maxVal=1, transparency=0.7,
                                minCol=c(0.6,0.5,0.5,transparency), maxCol=c(0,1,1,transparency) ) {
  hramp <- (maxCol[1]-minCol[1])/(maxVal-minVal)
  sramp <- (maxCol[2]-minCol[2])/(maxVal-minVal)
  vramp <- (maxCol[3]-minCol[3])/(maxVal-minVal)
  tramp <- (maxCol[4]-minCol[4])/(maxVal-minVal)
  
  if (normalise) {
    vals <- vals-min(vals)
    vals <- vals/max(vals)
  }
  
  inds <- which(vals < minVal)
  if (length(inds)>0) {
    vals[inds] <- minVal
  }
  
  inds <- which(vals > maxVal)
  if (length(inds)>0) {
    vals[inds] <- maxVal
  }
  
  h_contCols <- vals*hramp+minCol[1]
  s_contCols <- vals*sramp+minCol[2]
  v_contCols <- vals*vramp+minCol[3]
  t_contCols <- vals*tramp+minCol[4]
  contCols   <- hsv(h_contCols,s_contCols,v_contCols,t_contCols)
  return( contCols )
}
