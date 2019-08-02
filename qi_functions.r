# Some everyday R functions that I have implemented to speed up workflow


RiskCalc <- function(df,risk_type,const,prior,inst,examp){
     # Computes the Bayes Average (instances per example)
     # i.e., calculate average with a prior
     # Args:
     #   df: dataframe
     #   group_arg: name to call bayes risk
     #   const: constant proportional to typical group size
     #   prior: prior for the group, i.e., tot_inst/tot_examp
     #   inst: number of instances in group
     #   examp: number of examples, i.e., number in group
     #
   df[[risk_type]] <- (const * prior + inst)/(const + examp)
   df

}


MissingToAvg <- function(df,col,values=c(NA),avg='mean'){
    # Replaces NA with average of column
    # Args:
    #   avg: 'mean' or 'median'

 if (is.numeric(df[[col]]) == TRUE){
   switch(avg,
          'mean' = {df[[col]][df[[col]] %in% values] <- mean(df[[col]][!df[[col]] %in% values],na.rm=TRUE)},
          'median' = {df[[col]][df[[col]] %in% values] <- median(df[[col]][!df[[col]] %in% values],na.rm=TRUE)},
           message("unknown average type")
         )
 } else { message("yet to acount for non-numeric data types") }

 df
}


LowerLim <- function(df,column,lim){
   # Replaces all values below lim with the value lim
   # Args:
   #   df: the dataframe
   #   lim:
 df[[column]][df[[column]] < lim ] <- lim
 df
}



TruePredConversion <- function(df,uw.prob,uw.ir,t.ip){
   # convert the model output prediction to a true probability when
   # model has been run on oversampled data
   # Args:
   #   uwpred: up-weight probability of inc
   #   uw.ir: up-weighted incident ratio
   #   t.ip: true incident probability

   prob <- (df[[uw.prob]]*(1-uw.ir)*t.ip/
                      (uw.ir*(1-t.ip)*(1-df[[uw.prob]])+df[[uw.prob]]*(1-uw.ir)*t.ip))
   prob
}



SqlAppend <- function(channel,tablename,data){
 # Function to write SQL data using ODBC connection
 #
 # Args:
 #  channel: channel to write to
 #  tablename: tablename
 #  data: data as data frame

 ch <- odbcDriverConnect(channel)

 start_time <- Sys.time()
 sqlSave(ch,data,tablename,append=TRUE,rownames = FALSE)
 end_time <- Sys.time()


 message(paste0("Succesfully written '",deparse(substitute(data)),"' to '",tablename,"'"))
 message(paste0('Write time: ',end_time-start_time))
 odbcClose(ch)
 message('Closing channel')

}


SqlSource <- function(channel,query.file,...){
 # Function to get SQL data from RightShip DW
 #
 # Args:
 #   query.file: Query file
 #   ...: SQL as character, adds 'AND (...)' to query.file
 sql.queries = sql.queries

 if (is.character(query.file)== FALSE){
   stop("File name must be a string")
 }

 query <- read_file(paste0(sql.queries, query.file))

 ch <- odbcDriverConnect(channel)

 if(missing(...)){
   query <- query
   args <- NULL

 } else{
   args <- paste("\nAND (",list(...),")",sep = "")
   args <- paste(args,sep="",collapse=" ")
   query <- paste(query,args)
 }

 data <- sqlQuery(ch,query)
 odbcClose(ch)
 message(paste0('Query file : ',query.file,' from ',sql.queries))
 message(paste0('With conditions: ',args))

 data

}
