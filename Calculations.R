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

#' This file contains main function which are used in calculation of specific rotation
#' and building it dependencies on time, concentration and temperature. 

add_comp <- function(name, mw, lab_code){
  #' This function adds new compound to the database
  #' name - name of the compound
  #' mw - molecular weight of the compound
  #' lab_code - code of the compound in laboratory notebook
  #' return - write data to the database
  
  source("SQL_functions.R")
  add_compound(name = name, weigth = mw, lab_code = lab_code)
}

ud <- function(ds, filename){
  #' This function calculates specific rotation and it's standart deviation.
  #' ds - dataset in format: time, alpha, temperature
  #' filename - name of the file with corresponding dataset
  #' return - specific rotation of this dataset
  
  conc <- get_concentration(filename) # read concentration data from the database
  mw <- get_molweight(filename) # read substance molecular weight
  pt <- get_path(filename) # read optical path lenght data
  mult <- pt/(conc*mw/10) # calculates multiplicator for futher calculation
  alphaud <- ds*mult # calculates specific rotation
  return(alphaud)
}

pol_alpha_conc <- function(foldername, alphaud, code, pass, size, aver, steps, type, error, shift, exclude, y_lim){
  #' This function build rotation-concentration dependency 
  #' foldername - path to the folder with experimental files (csv) 
  #' alphaud - bool, if True then function calculates specific rotation else it use experimental data
  #' code - name of the experiment. Example OA0902 Note that there is no simbols after digits.
  #' pass - the number of points from the begining of single dataframe that should be passed.
  #' size - the number of points of the single dataframe that should be imported.
  #' aver - the number of points that should be averaged.
  #' steps - step in seconds between two following points.
  #' type - the type of an experiment ('static' for measurement with constant temperature,
  #'        'dynamic' for cases experiment with variable temperature, 'reference' - for references, 
  #'        'wash' - for washing of sapmle cell, 'empty' - for empty sample container)
  #' error - bool, show or hide error on the graphs
  #' shift - time in minutes between making the sample and start of measurements.
  #' exclude - which experiments should be excluded from analysis
  #' y_lim - y axis bounds
  
  # import required files
  source("SQL_functions.R")
  source("Data_IO.R")
  source("Output.R")
  
  # procesing missing values
  if(missing(alphaud)){alphaud = FALSE}
  if(missing(error)){error = TRUE}
  if(missing(shift)){shift = FALSE}
  if(missing(aver)){aver = 1}
  if(missing(exclude)){ex_flag = FALSE}
  else{ex_flag = TRUE}
  flag = TRUE # have no idea what it is doing, however when I wrote thi function it looked reasonable
  temp <- pol_entry(foldername, code = code, type = type) # read data from the folder into database
  if(ex_flag == TRUE){temp <- setdiff(temp, exclude)}
  
  # read data from the database and combine it into one dataframe
  df <- pol_combine_files(temp, pass = pass, size = size, aver = aver, steps = steps, shift = shift)
  print(df)
  
  route <- substring(temp[1], first = 1, last = 6)
  result <- c()
  
  if(alphaud == FALSE){
    for(i in 1:length(temp)){
      newdf <- subset(df, filename == temp[i])
      res <- c(get_concentration(temp[i]), round(mean(newdf[,2]), 5), round(sd(newdf[,2]), 5), round(mean(newdf[,3]), 2))
      result <- rbind(result, res)
    }
    rownames(result) <- temp
    colnames(result) <- c("Concentration","Alpha","SD","Temperature")
    result <- data.frame(result)
    ending = "_alpha_conc"
    plot_alpha_conc(df = result, code = route, ending = ending, error = error, flag = flag, y_lim = y_lim)
  }
  else{
    for(i in 1:length(temp)){
      newdf <- subset(df, filename == temp[i])
      alud <- ud(newdf$alpha, temp[i])
      res <- c(get_concentration(temp[i]), round(mean(alud), 5), round(sd(alud), 5), round(mean(newdf[,3]), 2))
      result <- rbind(result, res)
    }
    rownames(result) <- temp
    colnames(result) <- c("Concentration","AlphaUD","SD","Temperature")
    result <- data.frame(result)
    ending = "_alphaud_conc"
    plot_alphaud_conc(df = result, code = route, ending = ending, error = error, flag = flag, y_lim = y_lim)
  }
  write_file(route = route, result = result, ending = ending)
  return(result)
}

pol_alpha_time <- function(foldername, alphaud, code, pass, size, aver, steps, type, error, shift, exclude, y_lim){
  source("SQL_functions.R")
  source("Data_IO.R")
  source("Output.R")
  
  if(missing(alphaud)){alphaud = FALSE}
  if(missing(error)){error = TRUE}
  if(missing(aver)){aver = 1}
  if(missing(shift)){shift = FALSE}
  if(aver != 1){flag = TRUE}
  else{flag = FALSE}
  if(missing(exclude)){ex_flag = FALSE}
  else{ex_flag = TRUE}
  temp <- pol_entry(foldername, code = code, type = type)
  if(ex_flag == TRUE){temp <- setdiff(temp, exclude)}
  df <- pol_combine_files(temp = temp, pass = pass, size = size, aver = aver, steps = steps, shift = shift)
  route <- substring(temp[1], first = 1, last = 6)
  result <- c()
  if(alphaud == FALSE){
    df$time <- round(df[,1]/60, 2)
    df$SD <- round(df$SD, 5)
    sect <- c("time", "alpha", "SD", "temperature")
    df <- df[sect]
    ending = "_alpha_time"
    plot_alpha_time(df = df, code = route, ending = ending, error = error, flag = flag, y_lim = y_lim)
    if(flag){
      sec <- c("time", "alpha", "temperature")
      df <- df[sec]
    }
  }
  else{
    df$time <- round(df[,1]/60, 2)
    #print(df)
    for(i in 1:length(temp)){
      newdf <- subset(df, filename == temp[i])
      newdf$alphaud <- round(ud(newdf$alpha, temp[i]), 5)
      newdf$SD <- round(ud(newdf$SD, temp[i]), 5)
      sect <- c("time", "alphaud", "SD", "temperature")
      newdf <- newdf[sect]
      result <- rbind(result, newdf)
    }
    df <- result
    ending = "_alphaud_time"
    plot_alphaud_time(df = df, code = route, ending = ending, error = error, flag = flag, y_lim = y_lim)
    if(flag){
      sec <- c("time", "alphaud", "temperature")
      df <- df[sec]
    }
  }
  
  write_file(route = route, result = df, ending = ending)
  return(df)
}

pol_alpha_temp <- function(foldername, alphaud, code, pass, size, aver, steps, type, error, shift, exclude, y_lim){
  source("SQL_functions.R")
  source("Data_IO.R")
  source("Output.R")
  
  if(missing(alphaud)){alphaud = FALSE}
  if(missing(error)){error = TRUE}
  if(missing(aver)){aver = 1}
  if(missing(shift)){shift = FALSE}
  if(aver != 1){flag = TRUE}
  else{flag = FALSE}
  if(missing(exclude)){ex_flag = FALSE}
  else{ex_flag = TRUE}
  
  temp <- pol_entry(foldername, code = code, type = type)
  if(ex_flag == TRUE){temp <- setdiff(temp, exclude)}
  df <- pol_combine_files(temp = temp, pass = pass, size = size, aver = aver, steps = steps, shift = shift)
  route <- substring(temp[1], first = 1, last = 6)
  result <- c()
  
  if(alphaud == FALSE){
    df$time <- round(df[,1]/60, 2)
    df$SD <- round(df$SD, 5)
    sect <- c("time", "alpha", "SD", "temperature")
    df <- df[sect]
    ending = "_alpha_temp"
    plot_alpha_temp(df = df, code = route, ending = ending, error = error, flag = flag, y_lim = y_lim)
    if(flag){
      sec <- c("time", "alpha", "temperature")
      df <- df[sec]
    }
  }
  else{
    df$time <- round(df[,1]/60, 2)
    df$alphaud <- round(ud(df$alpha, temp[1]), 5)
    df$SD <- round(ud(df$SD, temp[1]), 5)
    sect <- c("time", "alphaud", "SD", "temperature")
    df <- df[sect]
    ending = "_alphaud_temp"
    plot_alphaud_temp(df = df, code = route, ending = ending, error = error, flag = flag, y_lim = y_lim)
    if(flag){
      sec <- c("time", "alphaud", "temperature")
      df <- df[sec]
    }
  }
  write_file(route = route, result = df, ending = ending)
  return(df)
}

pol_temp_time <- function(foldername, code, pass, size, aver, steps, type, shift, exclude){
  source("SQL_functions.R")
  source("Data_IO.R")
  source("Output.R")
  
  if(missing(shift)){shift = FALSE}
  if(missing(exclude)){ex_flag = FALSE}
  else{ex_flag = TRUE}

  temp <- pol_entry(foldername, code = code, type = type)
  if(ex_flag == TRUE){temp <- setdiff(temp, exclude)}
  df <- pol_combine_files(temp = temp, pass = pass, size = size, aver = aver, steps = steps, shift = shift)
  route <- substring(temp[1], first = 1, last = 6)
  df$time <- round(df[,1]/60, 2)
  sect <- c("time", "temperature")
  df <- df[sect]
  ending = "_temp_time"
  plot_temp_time(df = df, code = route, ending = ending)
  write_file(route = route, result = df, ending = ending)
  return(df)
}