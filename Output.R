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


write_file <- function(route, result, ending){
    dir.create(file.path("./Output/"), showWarnings = FALSE)
    name <- paste(route, ending, sep = "")
    write.csv(result, file = paste("./Output/", name, ".csv",sep = ""))
    wd <- getwd()
    print(paste("File was saved to ", wd, " ./Output/", name, ".csv",sep = ""))
}

plot_alpha_conc <- function(df, code, ending, error, flag, y_lim){
  name <- paste(code, ending, sep = "")
#  setEPS()
#  postscript(paste("./Output/",name,".eps",sep=""))
  if(missing(y_lim)){y_lim=range(c(df[,2]-df[,3], df[,2]+df[,3]))}
  if(error&flag){
  plot(df[,1], df[,2], ylim=y_lim, pch=19, bty = "l", xlab="c, mol/L", ylab = expression(alpha[D]))
  arrows(df[,1], df[,2]-df[,3], df[,1], df[,2]+df[,3], length=0.05, angle=90, code=3)}
  else{
    plot(df[,1], df[,2], pch=19, bty = "l", xlab="c, mol/L", ylab = expression(alpha[D]))
  }
#  dev.off()
  dev.copy2pdf(file = paste("./Output/",name,".pdf",sep = ""))
}

plot_alphaud_conc <- function(df, code, ending, error, flag, y_lim){
  name <- paste(code, ending, sep = "")
#  setEPS()
#  postscript(paste("./Output/",name,".eps",sep=""))
  if(missing(y_lim)){y_lim=range(c(df[,2]-df[,3], df[,2]+df[,3]))}
  if(error&flag){
    plot(df[,1], df[,2], ylim=y_lim , pch=19, bty = "l", xlab="c, mol/L", ylab = expression("["*alpha*"]"[D]))
    arrows(df[,1], df[,2]-df[,3], df[,1], df[,2]+df[,3], length=0.05, angle=90, code=3)}
  else{
    plot(df[,1], df[,2], pch=19, bty = "l", xlab="c, mol/L", ylab = expression("["*alpha*"]"[D]))
  }
  dev.copy2pdf(file = paste("./Output/",name,".pdf",sep = ""))
#  dev.off()
}

plot_alpha_time <- function(df, code, ending, error, flag, y_lim){
  name <- paste(code, ending, sep = "")
#  setEPS()
#  postscript(paste("./Output/",name,".eps",sep=""))
  if(missing(y_lim)){y_lim=range(c(min(df[,2])-max(df[,3]), max(df[,2])+max(df[,3])))}
  if(error&flag){
    plot(df[,1], df[,2], ylim=y_lim, pch=19, bty = "l", xlab="t, min", ylab = expression(alpha[D]))
    arrows(df[,1], df[,2]-df[,3], df[,1], df[,2]+df[,3], length=0.05, angle=90, code=3)}
  else{
    plot(df[,1], df[,2], pch=19, bty = "l", xlab="t, min", ylab = expression(alpha[D]))
  }
#  dev.off()
  dev.copy2pdf(file = paste("./Output/",name,".pdf",sep = ""))
}

plot_alphaud_time <- function(df, code, ending, error, flag, y_lim){
  name <- paste(code, ending, sep = "")
#  setEPS()
#  postscript(paste("./Output/",name,".eps",sep=""))
  if(missing(y_lim)){y_lim=range(c(min(df[,2])-max(df[,3]), max(df[,2])+max(df[,3])))}
  if(error&flag){
    plot(df[,1], df[,2], ylim=y_lim, pch=19, bty = "l", xlab="t, min", ylab = expression("["*alpha*"]"[D]))
    arrows(df[,1], df[,2]-df[,3], df[,1], df[,2]+df[,3], length=0.05, angle=90, code=3)}
  else{
    plot(df[,1], df[,2], pch=19, bty = "l", xlab="t, min", ylab = expression("["*alpha*"]"[D]))
  }
#  dev.off()
  dev.copy2pdf(file = paste("./Output/",name,".pdf",sep = ""))
}

plot_alpha_temp <- function(df, code, ending, error, flag, y_lim){
  name <- paste(code, ending, sep = "")
#  setEPS()
#  postscript(paste("./Output/",name,".eps",sep=""))
  if(missing(y_lim)){y_lim=range(c(min(df[,2])-max(df[,3]), max(df[,2])+max(df[,3])))}
  if(error&flag){
    plot(df[,4], df[,2], ylim=y_lim, pch=19, bty = "l", xlab="t, min", ylab = expression(alpha[D]))
    arrows(df[,4], df[,2]-df[,3], df[,4], df[,2]+df[,3], length=0.05, angle=90, code=3)}
  else{
    plot(df[,4], df[,2], pch=19, bty = "l", xlab="t, min", ylab = expression(alpha[D]))
  }
#  dev.off()
  dev.copy2pdf(file = paste("./Output/",name,".pdf",sep = ""))
}

plot_alphaud_temp <- function(df, code, ending, error, flag, y_lim){
  name <- paste(code, ending, sep = "")
#  setEPS()
#  postscript(paste("./Output/",name,".eps",sep=""))
  if(missing(y_lim)){y_lim = range(c(min(df[,2])-max(df[,3]), max(df[,2])+max(df[,3])))}
  if(error&flag){
    plot(df[,4], df[,2], ylim=y_lim, pch=19, bty = "l", xlab=expression(t^o*C), ylab = expression("["*alpha*"]"[D]))
    arrows(df[,4], df[,2]-df[,3], df[,4], df[,2]+df[,3], length=0.05, angle=90, code=3)}
  else{
    plot(df[,4], df[,2], pch=19, bty = "l", xlab=expression(t^o*C), ylab = expression("["*alpha*"]"[D]))
  }
#  dev.off()
  dev.copy2pdf(file = paste("./Output/",name,".pdf",sep = ""))
}

plot_temp_time <- function(df, code, ending){
  name <- paste(code, ending, sep = "")
#  setEPS()
#  postscript(paste("./Output/",name,".eps",sep=""))
  plot(df[,1], df[,2], pch=19, bty = "l", xlab="t, min", ylab = expression(t^o*C))
#  dev.off()
  dev.copy2pdf(file = paste("./Output/",name,".pdf",sep = ""))
}

#'Change log
#'19.01.17 Fixed "flag" issue in plotting concentration-dependent graphs