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

#' This block is needed for initialisation of emtpy database
install.packages("RSQLite")

source("SQL_functions.R")
create_database()
add_device("P-2000", "Polarimeter")
add_compound("Missing data", 1, "AAAAAA")

solvents <- c("NA", "water", "EtOH", "MeOH", "tBuOH", "iPrOH", "Et2O", "THF", "DCM", "CHl3", "DCE", "DMF", "DMSO", "CH3CN", "acetone", "AAAAAA")

dbSendQuery(conn = db, "ALTER TABLE POL_Description ADD COLUMN shift")

for(i in 1:length(solvents)){
  source("SQL_functions.R")
  get_solvent_id(solvents[i])  
}

types <- c('static', 'dynamic', 'reference', 'wash', 'empty')
for(n in 1:length(types)){
  source("SQL_functions.R")
  get_type_id(types[n])
}

add_compound("Solvent", 1, "Solvent")
add_compound("Air", 1, "Zero")