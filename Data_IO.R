#'***********************************************************************
#'*This file is part of RP-2000.                                        *
#'*Copyright (c) 2017 Nikolay Kondakov (k2n@yandex.ru).                 *
#'*                                                                     *
#'*RP-2000 is free software: you can redistribute it and/or modify      *
#'*it under the terms of the GNU General Public License as published by *
#'*the Free Software Foundation, either version 3 of the License, or    *
#'*any later version.                                                   *
#'*                                                                     *
#'*RP-2000 is distributed in the hope that it will be useful,           *
#'*but WITHOUT ANY WARRANTY; without even the implied warranty of       *
#'*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
#'*GNU General Public License for more details.                         *
#'*                                                                     *
#'*You should have received a copy of the GNU General Public License    *
#'*along with RP-2000. If not, see <http://www.gnu.org/licenses/>.      *
#'***********************************************************************

#'This module is responsible for importing data from text files.

add_lab_code <- function(name, weigth, lab_code){
  source("SQL_functions.R")
  add_compound(name, weigth, lab_code)
}

import_pol_descr <- function(route, filename, concentration, path, temp, type, lab_code, solvent, shift){
  checkpoint <- get_pol_descr_id(filename)  
  if(length(checkpoint) == 0){
    #Read file
    list <- readLines(route)
    
    #Find and read number of datapoints
    sizestr <- pmatch("NPOINTS,",list)
    size <- strsplit(list[sizestr], ',')
    npoint <- as.numeric(gsub(" ", "", size[[1]][2], fixed = TRUE))
    
    #Find and read file creation time.
    timeline <- pmatch("TIME,",list)
    timestring <- strsplit(list[timeline], ',')
    tim <- gsub(" ", "", timestring[[1]][2], fixed = TRUE)
    dateline <- pmatch("DATE,",list)
    datestring <- strsplit(list[dateline], ',')
    dt <- gsub(" ", "", datestring[[1]][2], fixed = TRUE)
    datet <- paste(dt, tim, sep = " ")
    datetime <- as.POSIXct(strptime(datet, "%y/%m/%d %H:%M:%S"))
    datetime <- datetime - npoint
  
    #Find Comment and splits it apart.
    com <- pmatch("Comment,",list)
    commen <- strsplit(gsub(" ", "", list[com]), ',', fixed = FALSE, perl = FALSE, useBytes = FALSE)
    
    #Read concentration data from comment.
    if(missing(concentration)) {
      concentration <- as.numeric((strsplit(commen[[1]][2], " "))[[1]][1])
      if(is.na(concentration)){
        warning("The information about sample concentration wasn't given. Now it is set to 0")
        concentration = 0
      }
      if(concentration > 5){
        warning(paste("Are you sure that concentration is", concentration, sep = " "))
      }
    }  
    
    #Read temperature from comment.
    if(missing(temp)) {
      temp <- as.numeric(strsplit(commen[[1]][3], " ")[[1]][1])
      #'!!!!!!!!!!!!!!Add conditional check if temp == "NA" than get average value from POL_Measurements table
      if(is.na(temp)){
        warning("The information about desired experiment temperature wasn't given. Later I will set it based on date stored in measurements table.")
      }
    } 
    
    #Read experiment type from comment.
    if(missing(type)) {
      type <- (strsplit(commen[[1]][4], " "))[[1]][1]
      if(is.na(type)){
        warning("The information about experiment type wasn't given (it wasn't set manualy or stored in the file). Experiment type automatically set as 'Static'")
        type = "static"}
    } 
  
    #Find "Sample name" and splits it apart.
    sam <- pmatch("Sample name,",list)
    samp <- strsplit(gsub(" ", "", list[sam]), ',', fixed = FALSE, perl = FALSE, useBytes = FALSE)
    
    #Add labcode information
    if(missing(lab_code)) {
      lab_code <- strsplit(samp[[1]][2], " ")[[1]][1]
      if(is.na(lab_code)){
        warning("The information about sample labcode wasn't given. It is set to AAAAAA")
        lab_code <- "AAAAAA"
        }
    }
    
    #Add solvent information
    if(missing(solvent)) {
      solvent <- strsplit(samp[[1]][3], " ")[[1]][1]
      if(is.na(solvent)){
        warning("The information about solvent wasn't given. It is set to AAAAAA")
        solvent <- "AAAAAA"
      }
    }
  
    #Add time-shift information
    if(missing(shift)) {
      shift <- strsplit(samp[[1]][4], " ")[[1]][1]
      if(is.na(shift)){
        warning("The information about time-shift wasn't given. It is set to 0")
        shift <- 0
      }
    }
    
    #Get aperture value
    apert <- pmatch("Aperture(S),", list)
    ap <- strsplit(list[apert], ',')[[1]][2]
    #regular expression is used to get values before mm
    aperture <- as.numeric(sub('mm.*', '', ap))
    
    #Get path Length
    if(missing(path)){
      pat <- pmatch("  Path Length,", list)
      pa <-  strsplit(list[pat], ',')[[1]][2]
      path <- as.numeric(sub(' mm.*', '', pa))
    }
    #Get device name
    dc <- pmatch("Instrument Name,", list)
    device <- strsplit(list[dc], ',')[[1]][2]
    
    #Get interval
    inter <- pmatch("D.I.T.,", list)
    interval <- strsplit(list[inter], ',')[[1]][2]
    source("SQL_functions.R")
    add_pol_descr(filename, datetime, npoint, temp, concentration, interval,
                   lab_code, solvent, path, aperture, type, device, shift)
    return(print("Experiments desription is read."))
  }
}

import_pol_mes <- function(route, filename){
  source("SQL_functions.R")
  sam_id <- get_pol_descr_id(filename)
  if(get_num_mes(filename, sam_id) == 0){
    size <- get_npoints(filename)
    #Read CSV file
    df <- read.csv(route, sep = ",", dec =".", skip = 21, nrows = size, header = FALSE, colClasses = rep("numeric",4), col.names = c("time","alpha", "temperature",""))[,1:3]
    df <- cbind("pol_descr_id" = sam_id, df)
    invisible(dbWriteTable(db, "POL_Measurments", df, append = TRUE));
  }
}

import_pol <- function(folder, type, concentration, path, temper, solvent, lab_code, shift){
  source("SQL_functions.R")
  track <- paste("./",folder,"/",sep = "")
  temp = list.files(path = track, pattern="*.csv")
  ds <- length(temp)
  if(!missing(concentration)){
    if(length(concentration) == 1){
      concentration = rep(concentration, times = ds)
    }
    if(length(concentration) != ds){stop("Length of concentration vector isn't equal to the number of files.")}
  }
  if(!missing(shift)){
    if(length(shift) == 1){
      shift = rep(shift, times = ds)
    }
    if(length(shift) != ds){stop("Length of time shift vector isn't equal to the number of files.")}
  }
    
  for(i in 1:length(temp)){
    route <- paste(track, temp[i], sep = "")
    print(temp[i])
    if(missing(concentration)){
      import_pol_descr(route, temp[i], type = type, path = path, temp = temper, solvent = solvent, lab_code = lab_code)
      import_pol_mes(route, temp[i])
    }
    else{
      import_pol_descr(route, temp[i], type = type, path = path, concentration = concentration[i], shift = shift[i] ,temp = temper, solvent = solvent, lab_code = lab_code)
      import_pol_mes(route, temp[i])
    }
  }
}

pol_entry <- function(foldername, code, type){
  source("SQL_functions.R")
  if(missing(code)){
    import_pol(foldername)
    track <- paste("./",foldername,"/",sep = "")
    temp = list.files(path = track, pattern="*.csv")
    code <- substring(temp[1], first = 1, last = 6)
  }
  if(missing(type)){type = 'static'}
  files <- c()
  for(i in 1:length(code)){
    print(code[i])
    print(type)
    files <- c(files, get_filenames(code[i], type))
  }
  files <- files[!is.na(files)]
  return(files)
}

pol_get_one <- function(filename, pass, size, aver, steps, mshift){
  df <- get_pol_data(filename)

  if(missing(pass)){pass = 0}
  if((mshift != 0)&(pass == 0)){pass = 900} #If we want to ,take in account time required for sample preparation and doesn't set pass, then pass value will be set to 15 minutes.
  if(mshift != 0){pass = pass - (get_shift(filename) - mshift) * 60}
  if(missing(size)){sz = nrow(df) - pass}
  else{sz = size}
  
  if((pass + sz) > nrow(df)){
    sz = nrow(df) - pass
    }
  df <- df[pass:(pass + sz),]
  
  if(missing(steps)){steps = 1}
  if(steps != 1){
    tdf <- c()
    last = steps
    while(last < nrow(df)){
      tdf <- rbind(tdf, c(df[last, 1], df[last, 2], df[last, 3], df[last, 4]))
      last = last + steps
    }
    tdf <- rbind(tdf, c(df[nrow(df), 1], mean(df[last:nrow(df), 2]), df[nrow(df), 3], df[nrow(df), 4]))
    tdf <- tdf[complete.cases(tdf),]
    colnames(tdf) <- c("time","alpha","temperature","concentration")
    df <- data.frame(tdf)
    df$time <- as.POSIXct(df$time, origin="1970-01-01")
  }
  
  if(missing(aver)){aver = 1}
  if(aver != 1){
    tdf <- c()
    first = 1
    last = aver
    while(last < nrow(df)){
      tdf <- rbind(tdf, c(df[last, 1], mean(df[first:last, 2]), df[last, 3], df[last, 4], sd(df[first:last, 2])))
      first = last
      last = last + aver
    }
    tdf <- rbind(tdf, c(df[nrow(df), 1], mean(df[last:nrow(df), 2]), df[nrow(df), 3], df[nrow(df), 4], sd(df[last:nrow(df), 2])))
    tdf <- tdf[complete.cases(tdf),]
    colnames(tdf) <- c("time","alpha","temperature", "concentration", "SD")
    df <- data.frame(tdf)
    df$time <- as.POSIXct(df$time, origin="1970-01-01")
  }
  else{df$SD <- 0}
  df$filename <- filename
  return(df)
}

pol_combine_files <- function(temp, code, pass, size, aver, steps, shift){
  #'This function is really complex. It imports data from the folder into the database, then it
  #'get values from database. Then it combines them according to settings.
  #'Settings are:
  #'code - name of the experiment. Example OA0902 Note that there is no simbols after digits.
  #'pass - the number of points from the begining of single dataframe that should be passed.
  #'size - the number of points of the single dataframe that should be imported.
  #'aver - the number of points that should be averaged.
  #'steps - step in seconds between two following points.
  #'Then in performs Shapiro-Wilk of normality for that dataframe.
  #'At the end this function returns dataframe with following columns "time, alpha, temperature"
  source("SQL_functions.R")

  tshift <- c()
  if(shift == TRUE){
    for(n in 1:length(temp)){
      tshift <- cbind(tshift, as.numeric(get_shift(temp[n])))
    }
    k <- which.min(tshift)
    min_shift <- tshift[k]
  }
  else(min_shift = 0)

  newdf <- c()
  for(i in 1:length(temp)){
    df <- pol_get_one(filename = temp[i], size = size, pass = pass, aver = aver, mshift = min_shift)
    newdf <- rbind(newdf, df)
  }
  newdf$time <- as.integer(newdf$time - min(newdf[,1]))
  newdf$alpha <- round(newdf[,2],5)
  factor(newdf$filename)
  return(newdf)
}

#'Change log
#'25.01.17 Edited pol_entry function. Added
#'20.01.17 Changed import_pol function. Hopefully it now will be able to accept vector of concentrations
#'and time shifts.
#'19.01.17 Added shift paramenter to import_pol_descr. This parameter describes time between
#'sample preparation and start of analisys. 