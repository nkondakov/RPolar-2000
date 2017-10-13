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


#' SQL input-output and supporting functions are stored here.

#Connecting to existing database, or creation an empty database
library("RSQLite")
db <- dbConnect(SQLite(), dbname="supramers.sqlite")

#Creation of database structure and dependencies
create_database <- function(){
  dbSendQuery(conn = db, "PRAGMA encoding='UTF-8';")
  
  dbSendQuery(conn = db,
              "CREATE TABLE Solvents(id INTEGER PRIMARY KEY, solvent TEXT UNIQUE COLLATE NOCASE);"
  )
  dbSendQuery(conn = db,
              "CREATE TABLE Samples(id INTEGER PRIMARY KEY, name TEXT UNIQUE COLLATE NOCASE, molweight NUMERIC);"
  )
  dbSendQuery(conn = db,
              "CREATE TABLE DTypes(id INTEGER PRIMARY KEY, name TEXT UNIQUE COLLATE NOCASE);"
  )
  dbSendQuery(conn = db,
              "CREATE TABLE Compounds(id INTEGER PRIMARY KEY, name_id INTEGER,
              lab_code TEXT UNIQUE COLLATE NOCASE, comment TEXT COLLATE NOCASE, FOREIGN KEY(name_id) REFERENCES Samples(id));"
              )
  dbSendQuery(conn = db,
              "CREATE TABLE Devices(id INTEGER PRIMARY KEY, name TEXT UNIQUE COLLATE NOCASE, dtype_id INTEGER,
              FOREIGN KEY(dtype_id) REFERENCES DTypes(id));"
              )
  dbSendQuery(conn = db,
              "CREATE TABLE POL_Apertures(id INTEGER PRIMARY KEY, aperture NUMERIC UNIQUE);"
  )
  dbSendQuery(conn = db,
              "CREATE TABLE POL_Intervals(id INTEGER PRIMARY KEY, interval INTEGER UNIQUE);"
  )
  dbSendQuery(conn = db,
              "CREATE TABLE POL_Types(id INTEGER PRIMARY KEY, type TEXT UNIQUE COLLATE NOCASE);"
  )
  dbSendQuery(conn = db,
              "CREATE TABLE POL_Paths(id INTEGER PRIMARY KEY, path INTEGER UNIQUE);"
  )
  dbSendQuery(conn = db,
              "CREATE TABLE Experiments(id INTEGER PRIMARY KEY, experiment TEXT UNIQUE COLLATE NOCASE);"
  )
  dbSendQuery(conn = db,
              "CREATE TABLE POL_Description(id INTEGER PRIMARY KEY, filename TEXT UNIQUE COLLATE NOCASE, datetime TEXT,
              npoint INTEGER, temperature NUMERIC, concentration NUMERIC, interval_id INTEGER,
              compound_id INTEGER, solvent_id INTEGER, path_id INTEGER, aperture_id INTEGER, type_id INTEGER,
              exp_id INTEGER, device_id INTEGER, shift INTEGER,
              FOREIGN KEY(compound_id) REFERENCES Compounds(id),
              FOREIGN KEY(interval_id) REFERENCES POL_Intervals(id),
              FOREIGN KEY(path_id) REFERENCES POL_Paths(id),
              FOREIGN KEY(aperture_id) REFERENCES POL_Apertures(id),
              FOREIGN KEY(solvent_id) REFERENCES Solvents(id),
              FOREIGN KEY(type_id) REFERENCES POL_Types(id),
              FOREIGN KEY(exp_id) REFERENCES Experiments(id),
              FOREIGN KEY(device_id) REFERENCES Devices(id));"
  )
  dbSendQuery(conn = db,
              "CREATE TABLE POL_Measurments(pol_descr_id INTEGER, time INTEGER, alpha NUMERIC, temperature NUMERIC,
              FOREIGN KEY(pol_descr_id) REFERENCES POL_Description(id));"
              )
}

#Work with universal database tables

get_solvent_id <- function(solvent){
  trans <- dbSendQuery(db, "INSERT OR IGNORE INTO Solvents (solvent) VALUES (?)", list(solvent))
  dbClearResult(trans)
  sol_id <- dbSendQuery(db, "SELECT id FROM Solvents WHERE solvent = ?")
  dbBind(sol_id, list(solvent))
  sol <- dbFetch(sol_id, n = 1)
  dbClearResult(sol_id)
  return(sol$id)
}



get_exp_id <- function(experiment){
  trans <- dbSendQuery(db, "INSERT OR IGNORE INTO Experiments (experiment) VALUES (?)", list(experiment))
  dbClearResult(trans)
  i_id <- dbSendQuery(db, "SELECT id FROM Experiments WHERE experiment = ?")
  dbBind(i_id, list(experiment))
  in_id <- dbFetch(i_id, n = 1)
  dbClearResult(i_id)
  return(in_id$id)
}

get_compound_id <- function(lab_code){
  c_id <- dbSendQuery(db, "SELECT id FROM Compounds WHERE lab_code = ?")
  dbBind(c_id, list(lab_code))
  co_id <- dbFetch(c_id, n = 1)
  dbClearResult(c_id)
  return(co_id$id)
}

get_device_id <- function(name){
  c_id <- dbSendQuery(db, "SELECT id FROM Devices WHERE name = ?")
  dbBind(c_id, list(name))
  co_id <- dbFetch(c_id, n = 1)
  dbClearResult(c_id)
  return(co_id$id)
}

add_compound <- function(name, weigth, lab_code, comments){
  Encoding(name) <- "UTF-8"
  trans <- dbSendQuery(db, "INSERT OR IGNORE INTO Samples (name, molweight) VALUES (?, ?)", list(name, weigth))
  dbClearResult(trans)
  i_id <- dbSendQuery(db, "SELECT id FROM Samples WHERE name = ?")
  dbBind(i_id, list(name))
  in_id <- dbFetch(i_id, n = 1)
  dbClearResult(i_id)
  sam_id <- in_id$id
  if(missing(comments)){comments <- "NA"}
  trans_2 <- dbSendQuery(db, "INSERT OR IGNORE INTO Compounds (name_id, lab_code, comment) VALUES (?, ?, ?)"
                         , list(sam_id, lab_code, comments))
  dbClearResult(trans_2)
  return(get_compound_id(lab_code))
}

add_device <- function(name, type){
  trans <- dbSendQuery(db, "INSERT OR IGNORE INTO DTypes (name) VALUES (?)", list(type))
  dbClearResult(trans)
  i_id <- dbSendQuery(db, "SELECT id FROM DTypes WHERE name = ?")
  dbBind(i_id, list(type))
  in_id <- dbFetch(i_id, n = 1)
  dbClearResult(i_id)
  dtype_id <- in_id$id
  trans_2 <- dbSendQuery(db, "INSERT OR IGNORE INTO Devices (name, dtype_id) VALUES (?, ?)"
                         , list(name, dtype_id))
  dbClearResult(trans_2)
  return(get_device_id(name))
}

#Works with polarimeter-only tables

get_path_id <- function(path){
  trans <- dbSendQuery(db, "INSERT OR IGNORE INTO POL_Paths (path) VALUES (?)", list(path))
  dbClearResult(trans)
  p_id <- dbSendQuery(db, "SELECT id FROM POL_Paths WHERE path = ?")
  dbBind(p_id, list(path))
  path_id <- dbFetch(p_id, n = 1)
  dbClearResult(p_id)
  return(path_id$id)
}

get_aperture_id <- function(aperture){
  trans <- dbSendQuery(db, "INSERT OR IGNORE INTO POL_Apertures (aperture) VALUES (?)", list(aperture))
  dbClearResult(trans)
  a_id <- dbSendQuery(db, "SELECT id FROM POL_Apertures WHERE aperture = ?")
  dbBind(a_id, list(aperture))
  ap_id <- dbFetch(a_id, n = 1)
  dbClearResult(a_id)
  return(ap_id$id)
}

get_type_id <- function(type){
  trans <- dbSendQuery(db, "INSERT OR IGNORE INTO POL_Types (type) VALUES (?)", list(type))
  dbClearResult(trans)
  t_id <- dbSendQuery(db, "SELECT id FROM POL_Types WHERE type = ?")
  dbBind(t_id, list(type))
  ty_id <- dbFetch(t_id, n = 1)
  dbClearResult(t_id)
  return(ty_id$id)
}

get_interval_id <- function(interval){
  trans <- dbSendQuery(db, "INSERT OR IGNORE INTO POL_Intervals (interval) VALUES (?)", list(interval))
  dbClearResult(trans)
  i_id <- dbSendQuery(db, "SELECT id FROM POL_Intervals WHERE interval = ?")
  dbBind(i_id, list(interval))
  in_id <- dbFetch(i_id, n = 1)
  dbClearResult(i_id)
  return(in_id$id)
}

get_pol_descr_id <- function(filename){
  sam_id <- dbSendQuery(db, "SELECT id FROM POL_Description WHERE filename = ?")
  dbBind(sam_id, list(filename))
  sa_id <- dbFetch(sam_id, n = 1)
  dbClearResult(sam_id)
  return(sa_id$id)
}

get_npoints <- function(filename){
  c_id <- dbSendQuery(db, "SELECT npoint FROM POL_Description WHERE filename = ?")
  dbBind(c_id, list(filename))
  co_id <- dbFetch(c_id, n = 1)
  dbClearResult(c_id)
  return(co_id$npoint)
}

add_pol_descr <- function(filename, datetime, npoint, temp, concentration, interval,
                            lab_code, solvent, path, aperture, type, device, shift){
  interval <- get_interval_id(interval)
  compound <- get_compound_id(lab_code)
  solvent <- get_solvent_id(solvent)
  path <- get_path_id(path)
  aperture <- get_aperture_id(aperture)
  type <- get_type_id(type)
  exper <- get_exp_id(substring(filename, first = 1, last = 6)) #!!! this is how you can extract substring!!!
  device <- get_device_id(device)
  trans <- dbSendQuery(db, "INSERT OR IGNORE INTO POL_Description (filename, datetime, npoint, temperature,
                       concentration, interval_id, compound_id, solvent_id, path_id, aperture_id, type_id,
                       exp_id, device_id, shift) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
                         , list(filename, datetime, npoint, temp, concentration, interval, compound, 
                                solvent, path, aperture, type, exper, device, shift))
  dbClearResult(trans)
}

add_pol_measurement <- function(filename){
  sam_id <- get_pol_descr_id(filename)
  ts <- dbSendQuery(db, "SELECT id FROM POL_Measurements WHERE pol_descr_id = ?")
  dbBind(ts, list(sam_id))
  co_id <- dbFetch(ts, n = 1)
  dbClearResult(ts)
  sz = co_id$pol_descr_id
  if(sz == 0){
    size <- get_npoints(filename)
    #Read CSV file
    df <- read.csv(filename, sep = ",", dec =".", skip = 21, nrows = size, header = FALSE, colClasses = rep("numeric",4), col.names = c("time","alpha", "temperature",""))[,1:3]
    df <- cbind("sample_id" = sam_id, df)
    sink("/dev/null");
    invisible(dbWriteTable(db, "Measurments", df, append = TRUE));
    return(data)
  }
  else{
    warning("This sample already stored in Database. You are trying to store it for second time.")
  }
}

get_pol_data <- function(filename){
  sam_id <- get_pol_descr_id(filename)
  ts <- dbSendQuery(db, "SELECT time, alpha, temperature FROM POL_Measurments WHERE pol_descr_id = ?")
  dbBind(ts, list(sam_id))
  df <- dbFetch(ts)
  dbClearResult(ts)
  
  c_id <- dbSendQuery(db, "SELECT datetime FROM POL_Description WHERE filename = ?")
  dbBind(c_id, list(filename))
  co_id <- dbFetch(c_id, n = 1)
  dbClearResult(c_id)
  
  df$time <- as.integer(df$time) + as.integer(co_id$datetime)
  df$time <- as.POSIXct(df$time, origin="1970-01-01")
  df$concentration <- as.numeric(get_concentration(filename))
  return(df)
}

get_num_mes <- function(filename, sam_id){
  tsa <- dbSendQuery(db, "SELECT pol_descr_id FROM POL_Measurments WHERE pol_descr_id = ?")
  dbBind(tsa, list(sam_id))
  ts <- dbFetch(tsa, n = 1)
  dbClearResult(tsa)
  return(length(ts$pol_descr_id))
}

get_molweight <- function(filename){
  #Get compound_id based on filename
  c_id <- dbSendQuery(db, "SELECT compound_id FROM POL_Description WHERE filename = ?")
  dbBind(c_id, list(filename))
  comp_id <- dbFetch(c_id, n = 1)
  dbClearResult(c_id)
  
  #Get name_id based on compound_id
  n_id <- dbSendQuery(db, "SELECT name_id FROM Compounds WHERE id = ?")
  dbBind(n_id, list(comp_id$compound_id))
  name_id <- dbFetch(n_id, n = 1)
  dbClearResult(n_id)

  #Get molecular weight based on name_id
  w <- dbSendQuery(db, "SELECT molweight FROM Samples WHERE id = ?")
  dbBind(w, name_id$name_id)
  mw <- dbFetch(w, n = 1)
  dbClearResult(w)
  return(mw$molweight)
}

get_concentration <- function(filename){
  c_id <- dbSendQuery(db, "SELECT concentration FROM POL_Description WHERE filename = ?")
  dbBind(c_id, list(filename))
  conc <- dbFetch(c_id, n = 1)
  dbClearResult(c_id)
  return(conc$concentration)
}

get_shift <- function(filename){
  c_id <- dbSendQuery(db, "SELECT shift FROM POL_Description WHERE filename = ?")
  dbBind(c_id, list(filename))
  conc <- dbFetch(c_id, n = 1)
  dbClearResult(c_id)
  return(as.numeric(conc$shift))
}

get_path <- function(filename){
  p_id <- dbSendQuery(db, "SELECT path_id FROM POL_Description WHERE filename = ?")
  dbBind(p_id, list(filename))
  path_id <- dbFetch(p_id, n = 1)
  dbClearResult(p_id)
  
  #Get path based on path_id
  pa_id <- dbSendQuery(db, "SELECT path FROM POL_Paths WHERE id = ?")
  dbBind(pa_id, list(path_id$path_id))
  path <- dbFetch(pa_id, n = 1)
  dbClearResult(pa_id)
  return(path$path)
}

get_filenames <- function(experiment, type){
  exp_id <- dbSendQuery(db, "SELECT id FROM Experiments WHERE experiment = ?")
  dbBind(exp_id, list(experiment))
  e_id <- dbFetch(exp_id, n = 1)
  dbClearResult(exp_id)
  print(exp_id)
  t_id <- get_type_id(type)

  #Get path based on path_id
  pa_id <- dbSendQuery(db, "SELECT filename FROM POL_Description WHERE exp_id = ? AND type_id = ?")
  dbBind(pa_id, list(e_id$id, t_id))
  filenames <- dbFetch(pa_id)
  dbClearResult(pa_id)
  
  name <- filenames$filename[1]

  for(i in 2:length(filenames$filename)){
    temp_name <- filenames$filename[i]
    name <- c(name, filenames$filename[i])
  }
  return(name)
}

get_subs_name <- function(lab_code){
  #Get name_id based on compound_id
  n_id <- dbSendQuery(db, "SELECT name_id FROM Compounds WHERE lab_code = ?")
  dbBind(n_id, list(lab_code))
  name_id <- dbFetch(n_id, n = 1)
  dbClearResult(n_id)
  
  #Get molecular weight based on name_id
  w <- dbSendQuery(db, "SELECT name FROM Samples WHERE id = ?")
  dbBind(w, name_id$name_id)
  mw <- dbFetch(w, n = 1)
  dbClearResult(w)
  return(mw$name)
}

#'Added "shift" field to create_database, add_pol_descr functions and get_shift. This parameter describes time between sample preparation 
#'and start of analisys. This functions create database, writes shift to database and reads it from there.