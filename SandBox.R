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


#'This file is used for testing idea of new functions and options.

source("SQL_functions.R")
source("Data_IO.R")
source("Output.R")

ud <- function(ds, filename){
  #This function calculates specific rotation and it's standart deviation.
  conc <- get_concentration(filename)
  mw <- get_molweight(filename)
  pt <- get_path(filename)
  mult <- pt/(conc*mw/10)
  alphaud <- ds*mult
  return(alphaud)
}

pol_alpha_time_2 <- function(foldername, alphaud, code, pass, size, aver, steps, type, error, shift, exclude, y_lim){
  source("SQL_functions.R")
  source("Data_IO.R")
  source("Output.R")
  source("Calculations.R")
  
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
    dualAxis(df = df, code = route, ending = ending, error = error, flag = flag, y_lim = y_lim)
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
    ending = "_alphaud_time"
    dualAxis(df = df, code = route, ending = ending, error = error, flag = flag, y_lim = y_lim)
    if(flag){
      sec <- c("time", "alphaud", "temperature")
      df <- df[sec]
    }
  }
  write_file(route = route, result = df, ending = ending)
  return(df)
}

dualAxis <- function(df, code, ending, error, flag, y_lim){
  name <- paste(code, ending, sep = "")
  setEPS()
  postscript(paste("./Output/",name,".eps",sep=""))
  if(missing(y_lim)){y_lim=range(c(min(df[,2])-max(df[,3]), max(df[,2])+max(df[,3])))}
  if(error&flag){
    postscript(paste("./Output/",name,".eps",sep=""))
    plot(df[,1], df[,2], ylim=y_lim, pch=16, bty = "l", xlab="t, min", ylab = expression("["*alpha*"]"[D]))
    arrows(df[,1], df[,2]-df[,3], df[,1], df[,2]+df[,3], length=0.05, angle=90, code=3)
    par(new = T)
    plot(df[,1], df[,4], type="l", col="red3",yaxt="n", xlab="", ylab="")
    axis(side = 4, bty = 'l')
    #mtext(side = 4, line = 3, 'Temperature')
}
  else{
    plot(df[,1], df[,2], pch=19, bty = "l", xlab="t, min", ylab = expression("["*alpha*"]"[D]))
  }
  dev.off()
  dev.copy2pdf(file = paste("./Output/",name,".pdf",sep = ""))
}

pol_exp_aver <- function(foldername, alphaud, code, pass, size, aver, steps, type, error, shift, exclude, level){
  source("SQL_functions.R")
  source("Data_IO.R")
  source("Output.R")
  
  if(missing(alphaud)){alphaud = FALSE}
  if(missing(level)){level = 0.95}
  if(missing(error)){error = TRUE}
  if(missing(shift)){shift = FALSE}
  if(missing(aver)){aver = 1}
  if(missing(exclude)){ex_flag = FALSE}
  else{ex_flag = TRUE}
  flag = TRUE
  temp <- pol_entry(foldername, code = code, type = type)
  if(ex_flag == TRUE){temp <- setdiff(temp, exclude)}
  
  df <- pol_combine_files(temp, pass = pass, size = size, aver = aver, steps = steps, shift = shift)
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
  }
  stud <- t.test(result[,2], conf.level = level)
  print(stud$estimate)
  print(stud$conf.int)
  # write_file(route = route, result = result, ending = ending)
}