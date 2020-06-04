# function to get a set of dates year-month-day for making plots etc
# set up for COVID019 and year 2020 only
# S J Lycett
# 27 APril 2020


get_udates <- function(min_date=min_date, max_date=max_date,leap_year=TRUE) {
  
  days_in_month <- c(31,29,31,30,31,30,31,31,30,31,30,31)
  if (!leap_year) {
    days_in_month[2] <- 28
  }
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
